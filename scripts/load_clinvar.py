# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load ClinVar GRCh38 variant_summary into Parquet.

alleles.parquet is one row per GRCh38 record. Header columns are kept whole,
including PhenotypeIDS and PhenotypeList. allele_dedup.parquet is one row per
variant_key, or per AlleleID when the record has no VCF allele.
"""

from __future__ import annotations

import csv
import gzip
import os
import shutil
import sys
import time
from collections import Counter
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coords import canon_chrom, norm_rsid, variant_key  # noqa: E402

SRC = ROOT / "raw" / "clinvar" / "variant_summary.txt.gz"
OUT_DIR = ROOT / "warehouse" / "clinvar"
ALLELES_PATH = OUT_DIR / "alleles.parquet"
DEDUP_PATH = OUT_DIR / "allele_dedup.parquet"
QC_PATH = ROOT / "warehouse" / "qc" / "clinvar.md"
TRUTH_PATH = ROOT / "warehouse" / "qc" / "truth_set.tsv"
TMP_DIR = OUT_DIR / "_duckdb_tmp"

BATCH = 50_000
SEP_ID = "|"
SEP_TEXT = " | "
SEP_LONG = " || "

REQUIRED = [
    "AlleleID",
    "Type",
    "Name",
    "GeneID",
    "GeneSymbol",
    "HGNC_ID",
    "ClinicalSignificance",
    "ClinSigSimple",
    "LastEvaluated",
    "RS# (dbSNP)",
    "nsv/esv (dbVar)",
    "RCVaccession",
    "PhenotypeIDS",
    "PhenotypeList",
    "Origin",
    "OriginSimple",
    "Assembly",
    "ChromosomeAccession",
    "Chromosome",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "Cytogenetic",
    "ReviewStatus",
    "NumberSubmitters",
    "Guidelines",
    "TestedInGTR",
    "OtherIDs",
    "SubmitterCategories",
    "VariationID",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "SomaticClinicalImpact",
    "SomaticClinicalImpactLastEvaluated",
    "ReviewStatusClinicalImpact",
    "Oncogenicity",
    "OncogenicityLastEvaluated",
    "ReviewStatusOncogenicity",
    "SCVsForAggregateGermlineClassification",
    "SCVsForAggregateSomaticClinicalImpact",
    "SCVsForAggregateOncogenicityClassification",
]

DERIVED = [
    "chrom_raw",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "rsid",
    "review_stars",
]

# ClinVar uses these strings when a field is blank. They are not allele sequences
# and they are not oncogenicity calls.
_MISSING = {"", ".", "-", "na", "n/a", "none", "not provided"}

DEDUP_SCHEMA = pa.schema(
    [
        ("dedup_key", pa.string()),
        ("variant_key", pa.string()),
        ("n_rows", pa.int64()),
        ("AlleleID", pa.string()),
        ("VariationID", pa.string()),
        ("GeneID", pa.string()),
        ("GeneSymbol", pa.string()),
        ("HGNC_ID", pa.string()),
        ("Type", pa.string()),
        ("chrom_raw", pa.string()),
        ("chrom", pa.string()),
        ("pos", pa.int64()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("Start", pa.string()),
        ("Stop", pa.string()),
        ("rsid", pa.string()),
        ("ClinicalSignificance", pa.string()),
        ("ReviewStatus", pa.string()),
        ("review_stars", pa.int8()),
        ("PhenotypeIDS", pa.string()),
        ("PhenotypeList", pa.string()),
        ("OriginSimple", pa.string()),
        ("germline_somatic", pa.string()),
        ("SomaticClinicalImpact", pa.string()),
        ("Oncogenicity", pa.string()),
        ("germline_pathogenic", pa.bool_()),
    ]
)


def review_stars(status: str | None) -> int | None:
    """Map a ClinVar review phrase to 4..0. Unmatched phrases stay null."""
    if status is None:
        return 0
    text = status.strip().lower()
    if not text:
        return 0
    if "practice guideline" in text:
        return 4
    if "expert panel" in text:
        return 3
    if "multiple submitters" in text and "no conflicts" in text:
        return 2
    if "single submitter" in text or "conflicting classifications" in text:
        return 1
    if "no assertion" in text or "no interpretation" in text:
        return 0
    return None


def usable_allele(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    if text.lower() in _MISSING:
        return None
    return text.upper()


def usable_pos(value: str | None) -> int | None:
    if value is None:
        return None
    text = value.strip()
    if not text.isdigit():
        return None
    number = int(text)
    if number < 1:
        return None
    return number


def _pieces(value: str | None, sep: str):
    if not value:
        return
    for part in value.split(sep):
        text = part.strip()
        if text:
            yield text


def origin_class(value: str | None) -> str:
    """germline, somatic, germline/somatic, or other. Mixed rows stay mixed."""
    has_g = False
    has_s = False
    for token in _pieces(value, SEP_ID):
        low = token.lower()
        if "germline" in low:
            has_g = True
        if "somatic" in low:
            has_s = True
    if has_g and has_s:
        return "germline/somatic"
    if has_s:
        return "somatic"
    if has_g:
        return "germline"
    return "other"


def onco_has_call(value: str | None) -> bool:
    for token in _pieces(value, SEP_TEXT):
        low = token.strip().lower()
        if low in _MISSING:
            continue
        if low.startswith("no classification") or low.startswith("no assertion"):
            continue
        if low.startswith("no interpretation"):
            continue
        return True
    return False


def sig_flags(value: str | None) -> tuple[bool, bool]:
    """Pathogenic piece, and conflicting piece. Conflict is not pathogenic."""
    pathogenic = False
    conflict = False
    for phrase in _pieces(value, SEP_TEXT):
        for piece in phrase.split(";"):
            text = piece.strip().lower()
            if not text or text in _MISSING:
                continue
            if "conflict" in text:
                conflict = True
            elif "pathogenic" in text:
                pathogenic = True
    return pathogenic, conflict


def germline_pathogenic(sig: str | None, origin: str | None, onco: str | None) -> bool:
    """Germline-only pathogenic or likely pathogenic, with no oncogenicity call.

    Origin class germline/somatic is not germline. A conflicting piece blocks
    the flag even when another piece says pathogenic.
    """
    if origin_class(origin) != "germline":
        return False
    if onco_has_call(onco):
        return False
    pathogenic, conflict = sig_flags(sig)
    return pathogenic and not conflict


def _check_rules() -> None:
    assert review_stars("practice guideline") == 4
    assert review_stars("Reviewed by expert panel") == 3
    assert review_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert review_stars("criteria provided, single submitter") == 1
    assert review_stars("criteria provided, conflicting classifications") == 1
    assert review_stars("no assertion criteria provided") == 0
    assert review_stars("no interpretation for the single variant") == 0
    assert review_stars("") == 0
    assert review_stars(None) == 0
    assert review_stars("-") is None
    assert review_stars("no classification provided") is None
    assert review_stars("no classification for the single variant") is None
    assert review_stars("no classifications from unflagged records") is None
    assert usable_allele("na") is None
    assert usable_allele("NA") is None
    assert usable_allele("atct") == "ATCT"
    assert usable_pos("-1") is None
    assert usable_pos("117559590") == 117559590
    assert variant_key("15", 84799209, "G", "A") == "chr15:84799209:G:A"
    assert variant_key("15", 84799209, usable_allele("na"), usable_allele("na")) is None
    assert norm_rsid("267606829") == "rs267606829"
    assert norm_rsid("-1") is None
    assert canon_chrom("1") == "chr1"
    assert origin_class("germline") == "germline"
    assert origin_class("somatic") == "somatic"
    assert origin_class("germline/somatic") == "germline/somatic"
    assert origin_class("germline|somatic") == "germline/somatic"
    assert origin_class("unknown") == "other"
    assert origin_class("not applicable") == "other"
    assert germline_pathogenic("Pathogenic", "germline", "-") is True
    assert germline_pathogenic("Likely pathogenic", "germline", "-") is True
    assert germline_pathogenic("Pathogenic/Likely pathogenic", "germline", "-") is True
    assert germline_pathogenic("Pathogenic, low penetrance", "germline", "-") is True
    assert germline_pathogenic("Pathogenic; drug response", "germline", "-") is True
    assert germline_pathogenic("Pathogenic", "germline/somatic", "Likely oncogenic") is False
    assert germline_pathogenic("Pathogenic", "germline", "Likely oncogenic") is False
    assert germline_pathogenic("Pathogenic", "somatic", "-") is False
    assert germline_pathogenic("Pathogenic", "unknown", "-") is False
    assert (
        germline_pathogenic(
            "Conflicting classifications of pathogenicity", "germline", "-"
        )
        is False
    )
    assert (
        germline_pathogenic(
            "Conflicting classifications of pathogenicity | Pathogenic",
            "germline",
            "-",
        )
        is False
    )
    assert germline_pathogenic("Uncertain significance", "germline", "-") is False
    assert germline_pathogenic("drug response", "germline", "-") is False
    assert (
        germline_pathogenic(
            "Pathogenic", "germline", "no classification for the single variant"
        )
        is True
    )
    assert germline_pathogenic("Pathogenic", "germline", "Uncertain significance") is False
    assert onco_has_call("- | Likely oncogenic") is True
    assert onco_has_call("-") is False


def alleles_schema(header: list[str]) -> pa.Schema:
    fields = [(name, pa.string()) for name in header]
    fields.extend(
        [
            ("chrom_raw", pa.string()),
            ("chrom", pa.string()),
            ("pos", pa.int64()),
            ("ref", pa.string()),
            ("alt", pa.string()),
            ("variant_key", pa.string()),
            ("rsid", pa.string()),
            ("review_stars", pa.int8()),
        ]
    )
    return pa.schema(fields)


def dedup_sql(relation: str) -> str:
    return f"""
    SELECT
      grp AS dedup_key,
      CASE WHEN grp LIKE 'allele:%' THEN NULL ELSE grp END AS variant_key,
      CAST(count(*) AS BIGINT) AS n_rows,
      string_agg(
        DISTINCT "AlleleID", '{SEP_ID}'
        ORDER BY try_cast("AlleleID" AS BIGINT), "AlleleID"
      ) AS "AlleleID",
      string_agg(
        DISTINCT "VariationID", '{SEP_ID}'
        ORDER BY try_cast("VariationID" AS BIGINT), "VariationID"
      ) AS "VariationID",
      string_agg(DISTINCT "GeneID", '{SEP_ID}' ORDER BY "GeneID") AS "GeneID",
      string_agg(DISTINCT "GeneSymbol", '{SEP_ID}' ORDER BY "GeneSymbol") AS "GeneSymbol",
      string_agg(DISTINCT "HGNC_ID", '{SEP_ID}' ORDER BY "HGNC_ID") AS "HGNC_ID",
      string_agg(DISTINCT "Type", '{SEP_ID}' ORDER BY "Type") AS "Type",
      string_agg(DISTINCT chrom_raw, '{SEP_ID}' ORDER BY chrom_raw) AS chrom_raw,
      string_agg(DISTINCT chrom, '{SEP_ID}' ORDER BY chrom) AS chrom,
      CASE
        WHEN count(DISTINCT pos) = 1 THEN CAST(min(pos) AS BIGINT)
        ELSE CAST(NULL AS BIGINT)
      END AS pos,
      string_agg(DISTINCT ref, '{SEP_ID}' ORDER BY ref) AS ref,
      string_agg(DISTINCT alt, '{SEP_ID}' ORDER BY alt) AS alt,
      string_agg(
        DISTINCT "Start", '{SEP_ID}'
        ORDER BY try_cast("Start" AS BIGINT), "Start"
      ) AS "Start",
      string_agg(
        DISTINCT "Stop", '{SEP_ID}'
        ORDER BY try_cast("Stop" AS BIGINT), "Stop"
      ) AS "Stop",
      string_agg(DISTINCT rsid, '{SEP_ID}' ORDER BY rsid) AS rsid,
      string_agg(
        DISTINCT "ClinicalSignificance", '{SEP_TEXT}'
        ORDER BY "ClinicalSignificance"
      ) AS "ClinicalSignificance",
      string_agg(
        DISTINCT "ReviewStatus", '{SEP_TEXT}' ORDER BY "ReviewStatus"
      ) AS "ReviewStatus",
      CAST(max(review_stars) AS TINYINT) AS review_stars,
      string_agg(
        DISTINCT "PhenotypeIDS", '{SEP_LONG}' ORDER BY "PhenotypeIDS"
      ) AS "PhenotypeIDS",
      string_agg(
        DISTINCT "PhenotypeList", '{SEP_LONG}' ORDER BY "PhenotypeList"
      ) AS "PhenotypeList",
      string_agg(
        DISTINCT "OriginSimple", '{SEP_ID}' ORDER BY "OriginSimple"
      ) AS "OriginSimple",
      string_agg(
        DISTINCT "SomaticClinicalImpact", '{SEP_TEXT}'
        ORDER BY "SomaticClinicalImpact"
      ) AS "SomaticClinicalImpact",
      string_agg(
        DISTINCT "Oncogenicity", '{SEP_TEXT}' ORDER BY "Oncogenicity"
      ) AS "Oncogenicity"
    FROM (
      SELECT
        variant_key,
        "AlleleID",
        "VariationID",
        "GeneID",
        "GeneSymbol",
        "HGNC_ID",
        "Type",
        chrom_raw,
        chrom,
        pos,
        ref,
        alt,
        "Start",
        "Stop",
        rsid,
        "ClinicalSignificance",
        "ReviewStatus",
        review_stars,
        "PhenotypeIDS",
        "PhenotypeList",
        "OriginSimple",
        "SomaticClinicalImpact",
        "Oncogenicity",
        CASE
          WHEN variant_key IS NULL THEN 'allele:' || "AlleleID"
          ELSE variant_key
        END AS grp
      FROM {relation}
    )
    GROUP BY grp
    ORDER BY dedup_key
    """


def apply_flags(batch: pa.RecordBatch) -> tuple[pa.RecordBatch, dict[str, int]]:
    origins = batch.column("OriginSimple").to_pylist()
    sigs = batch.column("ClinicalSignificance").to_pylist()
    oncos = batch.column("Oncogenicity").to_pylist()
    classes: list[str] = []
    flags: list[bool] = []
    counts = Counter()
    for origin, sig, onco in zip(origins, sigs, oncos):
        klass = origin_class(origin)
        flagged = germline_pathogenic(sig, origin, onco)
        classes.append(klass)
        flags.append(flagged)
        pathogenic, conflict = sig_flags(sig)
        call = onco_has_call(onco)
        counts["rows"] += 1
        if flagged:
            counts["germline_pathogenic"] += 1
        elif pathogenic and not conflict:
            if klass != "germline":
                counts["excluded_origin"] += 1
            if call:
                counts["excluded_onco"] += 1
            if klass != "germline" and call:
                counts["excluded_both"] += 1
            if klass == "germline" and call:
                counts["excluded_onco_only"] += 1
    batch = batch.append_column("germline_somatic", pa.array(classes, pa.string()))
    batch = batch.append_column("germline_pathogenic", pa.array(flags, pa.bool_()))
    return batch, counts


def align(batch: pa.RecordBatch, schema: pa.Schema) -> pa.RecordBatch:
    arrays = []
    for field in schema:
        col = batch.column(field.name)
        if col.type != field.type:
            col = col.cast(field.type, safe=True)
        arrays.append(col)
    return pa.record_batch(arrays, schema=schema)


def _sample_table() -> pa.Table:
    rows = [
        {
            "variant_key": "chr11:126275389:C:T",
            "AlleleID": "15044",
            "VariationID": "5",
            "GeneID": "1",
            "GeneSymbol": "FOXRED1",
            "HGNC_ID": "HGNC:1",
            "Type": "single nucleotide variant",
            "chrom_raw": "11",
            "chrom": "chr11",
            "pos": 126275389,
            "ref": "C",
            "alt": "T",
            "Start": "126275389",
            "Stop": "126275389",
            "rsid": "rs267606829",
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
            "review_stars": 2,
            "PhenotypeIDS": "MONDO:MONDO:1,MedGen:C1",
            "PhenotypeList": "disease A|disease B",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": "chr11:126275389:C:T",
            "AlleleID": "15044",
            "VariationID": "5",
            "GeneID": "1",
            "GeneSymbol": "FOXRED1",
            "HGNC_ID": "HGNC:1",
            "Type": "single nucleotide variant",
            "chrom_raw": "11",
            "chrom": "chr11",
            "pos": 126275389,
            "ref": "C",
            "alt": "T",
            "Start": "126275389",
            "Stop": "126275389",
            "rsid": "rs267606829",
            "ClinicalSignificance": "Likely pathogenic",
            "ReviewStatus": "criteria provided, single submitter",
            "review_stars": 1,
            "PhenotypeIDS": "MedGen:C9",
            "PhenotypeList": "disease C",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": "chr17:31334927:C:T",
            "AlleleID": "15382",
            "VariationID": "343",
            "GeneID": "4763",
            "GeneSymbol": "NF1",
            "HGNC_ID": "HGNC:7765",
            "Type": "single nucleotide variant",
            "chrom_raw": "17",
            "chrom": "chr17",
            "pos": 31334927,
            "ref": "C",
            "alt": "T",
            "Start": "31334927",
            "Stop": "31334927",
            "rsid": "rs137854552",
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
            "review_stars": 2,
            "PhenotypeIDS": "MedGen:NF1",
            "PhenotypeList": "Neurofibromatosis",
            "OriginSimple": "germline/somatic",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "Likely oncogenic",
        },
        {
            "variant_key": None,
            "AlleleID": "17413",
            "VariationID": "2374",
            "GeneID": "-",
            "GeneSymbol": "-",
            "HGNC_ID": "-",
            "Type": "Indel",
            "chrom_raw": "9",
            "chrom": "chr9",
            "pos": None,
            "ref": None,
            "alt": None,
            "Start": "12891379",
            "Stop": "12891380",
            "rsid": None,
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "criteria provided, single submitter",
            "review_stars": 1,
            "PhenotypeIDS": "MedGen:C1",
            "PhenotypeList": "phenotype one",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": None,
            "AlleleID": "17413",
            "VariationID": "2374",
            "GeneID": "-",
            "GeneSymbol": "-",
            "HGNC_ID": "-",
            "Type": "Indel",
            "chrom_raw": "20",
            "chrom": "chr20",
            "pos": None,
            "ref": None,
            "alt": None,
            "Start": "4720489",
            "Stop": "4720490",
            "rsid": None,
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "criteria provided, single submitter",
            "review_stars": 1,
            "PhenotypeIDS": "MedGen:C1",
            "PhenotypeList": "phenotype one",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": "chrY:641037:C:T",
            "AlleleID": "24911",
            "VariationID": "9872",
            "GeneID": "6473",
            "GeneSymbol": "SHOX",
            "HGNC_ID": "HGNC:10853",
            "Type": "single nucleotide variant",
            "chrom_raw": "Y",
            "chrom": "chrY",
            "pos": 641037,
            "ref": "C",
            "alt": "T",
            "Start": "641037",
            "Stop": "641037",
            "rsid": "rs1",
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "practice guideline",
            "review_stars": 4,
            "PhenotypeIDS": "MedGen:SHOX",
            "PhenotypeList": "short stature",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": "chrX:641037:C:T",
            "AlleleID": "24911",
            "VariationID": "9872",
            "GeneID": "6473",
            "GeneSymbol": "SHOX",
            "HGNC_ID": "HGNC:10853",
            "Type": "single nucleotide variant",
            "chrom_raw": "X",
            "chrom": "chrX",
            "pos": 641037,
            "ref": "C",
            "alt": "T",
            "Start": "641037",
            "Stop": "641037",
            "rsid": "rs1",
            "ClinicalSignificance": "Pathogenic",
            "ReviewStatus": "practice guideline",
            "review_stars": 4,
            "PhenotypeIDS": "MedGen:SHOX",
            "PhenotypeList": "short stature",
            "OriginSimple": "germline",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
        {
            "variant_key": "chr1:10:A:G",
            "AlleleID": "9",
            "VariationID": "9",
            "GeneID": "9",
            "GeneSymbol": "DEMO",
            "HGNC_ID": "-",
            "Type": "single nucleotide variant",
            "chrom_raw": "1",
            "chrom": "chr1",
            "pos": 10,
            "ref": "A",
            "alt": "G",
            "Start": "10",
            "Stop": "10",
            "rsid": None,
            "ClinicalSignificance": "-",
            "ReviewStatus": "-",
            "review_stars": None,
            "PhenotypeIDS": "-",
            "PhenotypeList": "-",
            "OriginSimple": "unknown",
            "SomaticClinicalImpact": "-",
            "Oncogenicity": "-",
        },
    ]
    columns = {name: [] for name in rows[0]}
    for row in rows:
        for name, value in row.items():
            columns[name].append(value)
    schema = pa.schema(
        [
            ("variant_key", pa.string()),
            ("AlleleID", pa.string()),
            ("VariationID", pa.string()),
            ("GeneID", pa.string()),
            ("GeneSymbol", pa.string()),
            ("HGNC_ID", pa.string()),
            ("Type", pa.string()),
            ("chrom_raw", pa.string()),
            ("chrom", pa.string()),
            ("pos", pa.int64()),
            ("ref", pa.string()),
            ("alt", pa.string()),
            ("Start", pa.string()),
            ("Stop", pa.string()),
            ("rsid", pa.string()),
            ("ClinicalSignificance", pa.string()),
            ("ReviewStatus", pa.string()),
            ("review_stars", pa.int8()),
            ("PhenotypeIDS", pa.string()),
            ("PhenotypeList", pa.string()),
            ("OriginSimple", pa.string()),
            ("SomaticClinicalImpact", pa.string()),
            ("Oncogenicity", pa.string()),
        ]
    )
    arrays = [pa.array(columns[field.name], type=field.type) for field in schema]
    return pa.Table.from_arrays(arrays, schema=schema)


def _check_dedup() -> None:
    con = duckdb.connect()
    con.register("sample_alleles", _sample_table())
    reader = con.sql(dedup_sql("sample_alleles")).to_arrow_reader(batch_size=20)
    batches = [align(apply_flags(batch)[0], DEDUP_SCHEMA) for batch in reader]
    table = pa.Table.from_batches(batches, schema=DEDUP_SCHEMA)
    by_key = {row["dedup_key"]: row for row in table.to_pylist()}
    assert len(by_key) == 6, len(by_key)
    fox = by_key["chr11:126275389:C:T"]
    assert fox["n_rows"] == 2
    assert fox["review_stars"] == 2
    assert fox["germline_pathogenic"] is True
    assert fox["germline_somatic"] == "germline"
    assert "Pathogenic" in fox["ClinicalSignificance"]
    assert "Likely pathogenic" in fox["ClinicalSignificance"]
    assert "MedGen:C9" in fox["PhenotypeIDS"]
    assert "MONDO:MONDO:1,MedGen:C1" in fox["PhenotypeIDS"]
    assert "disease A|disease B" in fox["PhenotypeList"]
    assert "disease C" in fox["PhenotypeList"]
    nf1 = by_key["chr17:31334927:C:T"]
    assert nf1["germline_pathogenic"] is False
    assert nf1["germline_somatic"] == "germline/somatic"
    assert nf1["Oncogenicity"] == "Likely oncogenic"
    assert nf1["ClinicalSignificance"] == "Pathogenic"
    assert nf1["SomaticClinicalImpact"] == "-"
    bare = by_key["allele:17413"]
    assert bare["variant_key"] is None
    assert bare["n_rows"] == 2
    assert bare["chrom_raw"] == "20|9"
    assert "chr20" in bare["chrom"] and "chr9" in bare["chrom"]
    assert bare["Start"] == "4720489|12891379"
    assert bare["germline_pathogenic"] is True
    assert bare["pos"] is None
    assert "chrX:641037:C:T" in by_key
    assert "chrY:641037:C:T" in by_key
    assert by_key["chrX:641037:C:T"]["review_stars"] == 4
    assert by_key["chr1:10:A:G"]["review_stars"] is None
    assert by_key["chr1:10:A:G"]["germline_pathogenic"] is False
    con.close()


def _flush(writer: pq.ParquetWriter, schema: pa.Schema, cols: dict[str, list]) -> int:
    n = len(cols["AlleleID"])
    if n == 0:
        return 0
    arrays = [pa.array(cols[field.name], type=field.type) for field in schema]
    writer.write_batch(pa.record_batch(arrays, schema=schema))
    for values in cols.values():
        values.clear()
    return n


def write_alleles() -> dict:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if ALLELES_PATH.exists():
        ALLELES_PATH.unlink()
    stats = {
        "assemblies": Counter(),
        "bad_len": 0,
        "grch38": 0,
        "written": 0,
        "max_pheno": 0,
        "max_pheno_allele": None,
        "max_pheno_chrom": None,
        "max_pheno_text": "",
        "max_ids": 0,
        "max_ids_allele": None,
        "max_ids_chrom": None,
        "max_ids_text": "",
        "sep_hits": [],
        "empty_allele": 0,
        "null_key": 0,
        "seconds": 0.0,
    }
    started = time.perf_counter()
    writer = None
    try:
        with gzip.open(SRC, "rt", encoding="utf-8", newline="") as fh:
            header = fh.readline().lstrip("\ufeff#").rstrip("\r\n").split("\t")
            missing = [name for name in REQUIRED if name not in header]
            if missing:
                raise SystemExit(f"variant_summary missing columns: {missing}")
            overlap = [name for name in DERIVED if name in header]
            if overlap:
                raise SystemExit(f"derived column names collide with header: {overlap}")
            schema = alleles_schema(header)
            idx = {name: i for i, name in enumerate(header)}
            writer = pq.ParquetWriter(
                ALLELES_PATH,
                schema,
                compression="zstd",
                compression_level=3,
                use_dictionary=False,
            )
            cols = {name: [] for name in schema.names}
            reader = csv.reader(fh, delimiter="\t")
            for rec in reader:
                if len(rec) != len(header):
                    stats["bad_len"] += 1
                    continue
                assembly = rec[idx["Assembly"]]
                stats["assemblies"][assembly] += 1
                if assembly != "GRCh38":
                    continue
                stats["grch38"] += 1
                chrom_raw = rec[idx["Chromosome"]]
                pos = usable_pos(rec[idx["PositionVCF"]])
                ref = usable_allele(rec[idx["ReferenceAlleleVCF"]])
                alt = usable_allele(rec[idx["AlternateAlleleVCF"]])
                key = variant_key(chrom_raw, pos, ref, alt)
                allele_id = rec[idx["AlleleID"]]
                if not allele_id:
                    stats["empty_allele"] += 1
                if key is None:
                    stats["null_key"] += 1
                phenotype = rec[idx["PhenotypeList"]]
                phenotype_ids = rec[idx["PhenotypeIDS"]]
                if len(phenotype) > stats["max_pheno"]:
                    stats["max_pheno"] = len(phenotype)
                    stats["max_pheno_allele"] = allele_id
                    stats["max_pheno_chrom"] = chrom_raw
                    stats["max_pheno_text"] = phenotype
                if len(phenotype_ids) > stats["max_ids"]:
                    stats["max_ids"] = len(phenotype_ids)
                    stats["max_ids_allele"] = allele_id
                    stats["max_ids_chrom"] = chrom_raw
                    stats["max_ids_text"] = phenotype_ids
                for name, sep in (
                    ("ClinicalSignificance", SEP_TEXT),
                    ("ReviewStatus", SEP_TEXT),
                    ("Oncogenicity", SEP_TEXT),
                    ("SomaticClinicalImpact", SEP_TEXT),
                    ("PhenotypeIDS", SEP_LONG),
                    ("PhenotypeList", SEP_LONG),
                    ("OriginSimple", SEP_ID),
                    ("GeneSymbol", SEP_ID),
                    ("ReferenceAlleleVCF", SEP_ID),
                    ("AlternateAlleleVCF", SEP_ID),
                ):
                    if sep in rec[idx[name]]:
                        if len(stats["sep_hits"]) < 8:
                            stats["sep_hits"].append((name, sep, rec[idx[name]][:80]))
                for name in header:
                    cols[name].append(rec[idx[name]])
                cols["chrom_raw"].append(chrom_raw)
                cols["chrom"].append(canon_chrom(chrom_raw))
                cols["pos"].append(pos)
                cols["ref"].append(ref)
                cols["alt"].append(alt)
                cols["variant_key"].append(key)
                cols["rsid"].append(norm_rsid(rec[idx["RS# (dbSNP)"]]))
                cols["review_stars"].append(review_stars(rec[idx["ReviewStatus"]]))
                if len(cols["AlleleID"]) >= BATCH:
                    stats["written"] += _flush(writer, schema, cols)
                    if stats["grch38"] % 500_000 < BATCH:
                        print(f"alleles {stats['grch38']}", flush=True)
            stats["written"] += _flush(writer, schema, cols)
    finally:
        if writer is not None:
            writer.close()
    stats["seconds"] = round(time.perf_counter() - started, 1)
    print(
        f"alleles written {stats['written']} in {stats['seconds']}s",
        flush=True,
    )
    return stats


def write_dedup() -> dict[str, int]:
    if DEDUP_PATH.exists():
        DEDUP_PATH.unlink()
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)
    counts: Counter = Counter()
    written = 0
    con = duckdb.connect()
    writer = None
    try:
        con.execute(f"SET temp_directory='{TMP_DIR.as_posix()}'")
        con.execute("SET memory_limit='6GB'")
        con.execute("SET threads=4")
        con.execute("SET preserve_insertion_order=false")
        con.execute("SET enable_progress_bar=false")
        path = ALLELES_PATH.as_posix().replace("'", "''")
        reader = con.sql(dedup_sql(f"read_parquet('{path}')")).to_arrow_reader(
            batch_size=100_000
        )
        writer = pq.ParquetWriter(
            DEDUP_PATH,
            DEDUP_SCHEMA,
            compression="zstd",
            compression_level=3,
            use_dictionary=False,
        )
        for batch in reader:
            flagged, batch_counts = apply_flags(batch)
            counts.update(batch_counts)
            out = align(flagged, DEDUP_SCHEMA)
            writer.write_batch(out)
            written += out.num_rows
            if written % 500_000 < out.num_rows:
                print(f"dedup {written}", flush=True)
    finally:
        if writer is not None:
            writer.close()
        con.close()
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    counts["written"] = written
    print(f"dedup written {written}", flush=True)
    return dict(counts)


def _connect() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("SET memory_limit='6GB'")
    con.execute("SET threads=4")
    con.execute("SET enable_progress_bar=false")
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)
    con.execute(f"SET temp_directory='{TMP_DIR.as_posix()}'")
    return con


def _q(con: duckdb.DuckDBPyConnection, sql: str):
    return con.execute(sql).fetchall()


def _one(con: duckdb.DuckDBPyConnection, sql: str):
    return con.execute(sql).fetchone()


def _rows(con: duckdb.DuckDBPyConnection, sql: str, params: list) -> list[dict]:
    cur = con.execute(sql, params)
    names = [item[0] for item in cur.description]
    return [dict(zip(names, row)) for row in cur.fetchall()]


def load_truth() -> dict[str, dict]:
    with TRUTH_PATH.open(encoding="utf-8", newline="") as fh:
        return {row["case_id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def _pass(rows: list[dict], ok: bool, detail: str) -> None:
    rows.append({"ok": ok, "detail": detail})


def gather(stats: dict, dedup_counts: dict) -> dict:
    con = _connect()
    try:
        alleles = ALLELES_PATH.as_posix().replace("'", "''")
        dedup = DEDUP_PATH.as_posix().replace("'", "''")
        a = f"read_parquet('{alleles}')"
        d = f"read_parquet('{dedup}')"
        allele_n, distinct_var, distinct_allele = _one(
            con,
            f"""
            SELECT count(*), count(DISTINCT "VariationID"), count(DISTINCT "AlleleID")
            FROM {a}
            """,
        )
        assembly_left = _q(con, f"SELECT \"Assembly\", count(*) FROM {a} GROUP BY 1")
        star_hist = _q(
            con,
            f"""
            SELECT review_stars, count(*)
            FROM {a}
            GROUP BY 1
            ORDER BY review_stars NULLS LAST
            """,
        )
        unmatched = _q(
            con,
            f"""
            SELECT "ReviewStatus", count(*)
            FROM {a}
            WHERE review_stars IS NULL
            GROUP BY 1
            ORDER BY 2 DESC, 1
            """,
        )
        matched = _q(
            con,
            f"""
            SELECT review_stars, "ReviewStatus", count(*)
            FROM {a}
            WHERE review_stars IS NOT NULL
            GROUP BY 1, 2
            ORDER BY 1, 3 DESC
            """,
        )
        null_keys, null_key_alleles = _one(
            con,
            f"""
            SELECT count(*), count(DISTINCT "AlleleID")
            FROM {a}
            WHERE variant_key IS NULL
            """,
        )
        distinct_keys = _one(
            con,
            f"SELECT count(DISTINCT variant_key) FROM {a} WHERE variant_key IS NOT NULL",
        )[0]
        dedup_n, dedup_keyed, dedup_bare, germ_n, sum_rows = _one(
            con,
            f"""
            SELECT
              count(*),
              count(*) FILTER (WHERE variant_key IS NOT NULL),
              count(*) FILTER (WHERE variant_key IS NULL),
              count(*) FILTER (WHERE germline_pathogenic),
              sum(n_rows)
            FROM {d}
            """,
        )
        dup_keys = _one(
            con,
            f"SELECT count(*) - count(DISTINCT dedup_key) FROM {d}",
        )[0]
        multi_key = _one(
            con,
            f"""
            SELECT count(*) FROM (
              SELECT variant_key FROM {a}
              WHERE variant_key IS NOT NULL
              GROUP BY variant_key
              HAVING count(*) > 1
            )
            """,
        )[0]
        multi_allele_same_key = _one(
            con,
            f"""
            SELECT count(*) FROM (
              SELECT variant_key FROM {a}
              WHERE variant_key IS NOT NULL
              GROUP BY variant_key
              HAVING count(DISTINCT "AlleleID") > 1
            )
            """,
        )[0]
        allele_dups = _one(
            con,
            f"""
            SELECT
              count(*),
              count(*) FILTER (WHERE n_chrom > 1),
              count(*) FILTER (WHERE n_sig > 1),
              count(*) FILTER (WHERE n_pheno > 1),
              count(*) FILTER (WHERE n_key > 1),
              count(*) FILTER (WHERE n_key = 0)
            FROM (
              SELECT
                "AlleleID",
                count(DISTINCT chrom_raw) AS n_chrom,
                count(DISTINCT "ClinicalSignificance") AS n_sig,
                count(DISTINCT "PhenotypeList") AS n_pheno,
                count(DISTINCT variant_key) AS n_key
              FROM {a}
              GROUP BY "AlleleID"
              HAVING count(*) > 1
            )
            """,
        )
        pheno_max, ids_max, pheno_over_500 = _one(
            con,
            f"""
            SELECT
              max(length("PhenotypeList")),
              max(length("PhenotypeIDS")),
              count(*) FILTER (WHERE length("PhenotypeList") > 500)
            FROM {a}
            """,
        )
        bad_chrom = _q(
            con,
            f"""
            SELECT coalesce(chrom_raw, '<null>'), count(*)
            FROM {a}
            WHERE chrom IS NULL
            GROUP BY 1
            ORDER BY 2 DESC
            """,
        )
        sigs = _q(
            con,
            f"""
            SELECT "ClinicalSignificance", count(*)
            FROM {a}
            GROUP BY 1
            ORDER BY 2 DESC, 1
            """,
        )
        dedup_stars = _q(
            con,
            f"""
            SELECT review_stars, count(*)
            FROM {d}
            GROUP BY 1
            ORDER BY review_stars NULLS LAST
            """,
        )
        non_acgt = _q(
            con,
            f"""
            SELECT ref, alt, count(*) AS n
            FROM {a}
            WHERE (ref IS NOT NULL AND regexp_matches(ref, '[^ACGTN]'))
               OR (alt IS NOT NULL AND regexp_matches(alt, '[^ACGTN]'))
            GROUP BY 1, 2
            ORDER BY n DESC
            LIMIT 15
            """,
        )
        truth = load_truth()
        checks = []
        _truth_checks(con, a, d, truth, checks)
        pheno_row = _rows(
            con,
            f"""
            SELECT "PhenotypeList", "PhenotypeIDS"
            FROM {a}
            WHERE "AlleleID" = ? AND chrom_raw = ?
            """,
            [stats["max_pheno_allele"], stats["max_pheno_chrom"]],
        )
        pheno_ok = any(row["PhenotypeList"] == stats["max_pheno_text"] for row in pheno_row)
        ids_row = _rows(
            con,
            f"""
            SELECT "PhenotypeIDS"
            FROM {a}
            WHERE "AlleleID" = ? AND chrom_raw = ?
            """,
            [stats["max_ids_allele"], stats["max_ids_chrom"]],
        )
        ids_ok = any(row["PhenotypeIDS"] == stats["max_ids_text"] for row in ids_row)
        checks.append(
            {
                "case_id": "phenotype_text",
                "ok": pheno_ok and ids_ok and pheno_max == stats["max_pheno"] and ids_max == stats["max_ids"],
                "detail": (
                    f"PhenotypeList max {pheno_max} source {stats['max_pheno']}, "
                    f"stored intact {pheno_ok}; PhenotypeIDS max {ids_max} "
                    f"source {stats['max_ids']}, stored intact {ids_ok}; "
                    f"rows longer than 500 characters: {pheno_over_500}"
                ),
            }
        )
    finally:
        con.close()
        shutil.rmtree(TMP_DIR, ignore_errors=True)

    integrity = {
        "allele_rows_match": allele_n == stats["grch38"] == stats["written"],
        "assembly_only_grch38": assembly_left == [("GRCh38", allele_n)],
        "bad_len": stats["bad_len"] == 0,
        "no_sep_collision": not stats["sep_hits"],
        "empty_allele": stats["empty_allele"] == 0,
        "dedup_covers_rows": sum_rows == allele_n,
        "dedup_keys_unique": dup_keys == 0,
        "dedup_keyed_match": dedup_keyed == distinct_keys,
        "dedup_bare_match": dedup_bare == null_key_alleles,
        "dedup_written": dedup_counts["written"] == dedup_n,
        "flag_count": dedup_counts["germline_pathogenic"] == germ_n,
    }
    return {
        "stats": stats,
        "allele_n": allele_n,
        "distinct_var": distinct_var,
        "distinct_allele": distinct_allele,
        "star_hist": star_hist,
        "unmatched": unmatched,
        "matched": matched,
        "null_keys": null_keys,
        "null_key_alleles": null_key_alleles,
        "distinct_keys": distinct_keys,
        "dedup_n": dedup_n,
        "dedup_keyed": dedup_keyed,
        "dedup_bare": dedup_bare,
        "germ_n": germ_n,
        "sum_rows": sum_rows,
        "multi_key": multi_key,
        "multi_allele_same_key": multi_allele_same_key,
        "allele_dups": allele_dups,
        "pheno_max": pheno_max,
        "ids_max": ids_max,
        "pheno_over_500": pheno_over_500,
        "bad_chrom": bad_chrom,
        "sigs": sigs,
        "dedup_stars": dedup_stars,
        "non_acgt": non_acgt,
        "dedup_counts": dedup_counts,
        "checks": checks,
        "integrity": integrity,
    }


def _truth_checks(con, a: str, d: str, truth: dict, checks: list) -> None:
    expected_keys = {
        "coding_snp": "chr11:126275389:C:T",
        "cftr_indel": "chr7:117559590:ATCT:A",
        "mismatched_ref": "chr15:84799209:G:A",
        "chrom_crosswalk": "chr1:243344276:AG:A",
        "somatic_or_oncogenicity": "chr17:31334927:C:T",
    }
    for case_id, key in expected_keys.items():
        got = truth.get(case_id, {}).get("variant_key", "")
        if got != key:
            checks.append(
                {
                    "case_id": case_id,
                    "ok": False,
                    "detail": f"truth_set.tsv key {got!r} != required {key}",
                }
            )

    coding = _rows(
        con,
        f"""
        SELECT "GeneSymbol", rsid, "ClinicalSignificance", "ReviewStatus",
               review_stars, ref, alt, chrom, chrom_raw, pos, variant_key
        FROM {a} WHERE variant_key = ?
        """,
        [expected_keys["coding_snp"]],
    )
    coding_ok = bool(coding) and all(
        row["GeneSymbol"] == "FOXRED1"
        and row["rsid"] == "rs267606829"
        and row["ClinicalSignificance"] == "Pathogenic"
        and row["review_stars"] == 2
        and "multiple submitters, no conflicts" in row["ReviewStatus"]
        and row["ref"] == "C"
        and row["alt"] == "T"
        and row["chrom"] == "chr11"
        and row["chrom_raw"] == "11"
        and row["pos"] == 126275389
        for row in coding
    )
    coding_dedup = _rows(
        con,
        f"""
        SELECT review_stars, germline_pathogenic, germline_somatic,
               "ClinicalSignificance", "Oncogenicity"
        FROM {d} WHERE variant_key = ?
        """,
        [expected_keys["coding_snp"]],
    )
    coding_flag = (
        len(coding_dedup) == 1
        and coding_dedup[0]["review_stars"] == 2
        and coding_dedup[0]["germline_pathogenic"] is True
        and coding_dedup[0]["germline_somatic"] == "germline"
        and coding_dedup[0]["Oncogenicity"] == "-"
    )
    checks.append(
        {
            "case_id": "coding_snp",
            "ok": coding_ok and coding_flag,
            "detail": (
                f"{len(coding)} allele row(s); dedup {coding_dedup}; "
                "FOXRED1 rs267606829 Pathogenic stars 2, germline_pathogenic true"
            ),
        }
    )

    cftr = _rows(
        con,
        f"""
        SELECT "GeneSymbol", "ClinicalSignificance", "ReviewStatus", review_stars,
               ref, alt, variant_key, length("PhenotypeList") AS pheno_len
        FROM {a} WHERE variant_key = ?
        """,
        [expected_keys["cftr_indel"]],
    )
    cftr_ok = bool(cftr) and all(
        row["ClinicalSignificance"] == "Pathogenic"
        and row["review_stars"] == 4
        and "practice guideline" in row["ReviewStatus"].lower()
        and row["ref"] == "ATCT"
        and row["alt"] == "A"
        and row["variant_key"] == expected_keys["cftr_indel"]
        and row["pheno_len"] > 500
        for row in cftr
    )
    checks.append(
        {
            "case_id": "cftr_indel",
            "ok": cftr_ok,
            "detail": (
                f"{cftr}; Pathogenic, stars 4, review contains practice guideline, "
                "alleles ATCT/A, PhenotypeList longer than 500"
            ),
        }
    )

    mismatch = _rows(
        con,
        f"""
        SELECT ref, alt, variant_key, "ReferenceAllele", "AlternateAllele",
               "ReferenceAlleleVCF", "AlternateAlleleVCF", "Name", "GeneSymbol"
        FROM {a} WHERE variant_key = ?
        """,
        [expected_keys["mismatched_ref"]],
    )
    mismatch_ok = bool(mismatch) and all(
        row["ref"] == "G"
        and row["alt"] == "A"
        and row["variant_key"] == "chr15:84799209:G:A"
        and row["ReferenceAllele"].strip().lower() == "na"
        and row["AlternateAllele"].strip().lower() == "na"
        and row["ReferenceAlleleVCF"] == "G"
        and row["AlternateAlleleVCF"] == "A"
        and row["Name"] != row["variant_key"]
        and "NM_" in row["Name"]
        and "chr15:84799209:G:A" not in row["Name"]
        for row in mismatch
    )
    checks.append(
        {
            "case_id": "mismatched_ref",
            "ok": mismatch_ok,
            "detail": (
                f"{[{k: row[k] for k in ('ref', 'alt', 'variant_key', 'ReferenceAllele', 'Name')} for row in mismatch]}; "
                "key uses VCF G/A, not display na and not the HGVS name"
            ),
        }
    )

    cross = _rows(
        con,
        f"""
        SELECT chrom_raw, chrom, variant_key, ref, alt, "GeneSymbol"
        FROM {a} WHERE variant_key = ?
        """,
        [expected_keys["chrom_crosswalk"]],
    )
    cross_counts = _one(
        con,
        f"""
        SELECT
          count(*) FILTER (WHERE chrom_raw = '1' AND chrom = 'chr1'),
          count(*) FILTER (WHERE chrom_raw = '1' AND chrom IS DISTINCT FROM 'chr1')
        FROM {a}
        """,
    )
    cross_ok = (
        bool(cross)
        and all(
            row["chrom_raw"] == "1" and row["chrom"] == "chr1" and row["ref"] == "AG"
            for row in cross
        )
        and cross_counts[0] > 0
        and cross_counts[1] == 0
    )
    checks.append(
        {
            "case_id": "chrom_crosswalk",
            "ok": cross_ok,
            "detail": (
                f"example {cross}; chrom_raw 1 -> chr1 on {cross_counts[0]} rows; "
                f"mismatches {cross_counts[1]}"
            ),
        }
    )

    somatic = _rows(
        con,
        f"""
        SELECT "GeneSymbol", "OriginSimple", "ClinicalSignificance", "Oncogenicity",
               "SomaticClinicalImpact", review_stars, ref, alt
        FROM {a} WHERE variant_key = ?
        """,
        [expected_keys["somatic_or_oncogenicity"]],
    )
    somatic_dedup = _rows(
        con,
        f"""
        SELECT germline_pathogenic, germline_somatic, "ClinicalSignificance",
               "Oncogenicity", "SomaticClinicalImpact", "OriginSimple", review_stars
        FROM {d} WHERE variant_key = ?
        """,
        [expected_keys["somatic_or_oncogenicity"]],
    )
    somatic_ok = (
        bool(somatic)
        and all(
            row["GeneSymbol"] == "NF1"
            and row["OriginSimple"] == "germline/somatic"
            and row["ClinicalSignificance"] == "Pathogenic"
            and row["Oncogenicity"] == "Likely oncogenic"
            and row["ref"] == "C"
            and row["alt"] == "T"
            for row in somatic
        )
        and len(somatic_dedup) == 1
        and somatic_dedup[0]["germline_pathogenic"] is False
        and somatic_dedup[0]["germline_somatic"] == "germline/somatic"
        and somatic_dedup[0]["ClinicalSignificance"] == "Pathogenic"
        and somatic_dedup[0]["Oncogenicity"] == "Likely oncogenic"
        and "SomaticClinicalImpact" in somatic_dedup[0]
    )
    checks.append(
        {
            "case_id": "somatic_or_oncogenicity",
            "ok": somatic_ok,
            "detail": (
                f"allele {somatic}; dedup {somatic_dedup}; "
                "not counted as germline pathogenic"
            ),
        }
    )


def _num(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def _table(headers: list[str], rows: list[tuple]) -> str:
    head = "| " + " | ".join(headers) + " |"
    rule = "|" + "|".join("---:" if i else "---" for i in range(len(headers)))
    # first column is text, the rest are numeric when they look like counts
    rule = "| " + " | ".join(["---"] + ["---:"] * (len(headers) - 1)) + " |"
    body = []
    for row in rows:
        cells = []
        for i, value in enumerate(row):
            if i == 0:
                cells.append(f"`{value}`" if value is not None else "`null`")
            else:
                cells.append(_num(value))
        body.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, rule, *body])


def render(report: dict) -> str:
    stats = report["stats"]
    assemblies = stats["assemblies"]
    asm_rows = []
    for name, count in sorted(assemblies.items(), key=lambda item: -item[1]):
        loaded = "loaded" if name == "GRCh38" else "excluded"
        asm_rows.append((name, count, loaded))
    star_rows = []
    for stars, count in report["star_hist"]:
        label = "null" if stars is None else str(stars)
        star_rows.append((label, count))
    unmatched_rows = [(status, count) for status, count in report["unmatched"]]
    matched_rows = [
        (str(stars), status, count) for stars, status, count in report["matched"]
    ]
    check_rows = []
    for check in report["checks"]:
        check_rows.append(
            (
                check["case_id"],
                "PASS" if check["ok"] else "FAIL",
                check["detail"].replace("|", "\\|"),
            )
        )
    dup = report["allele_dups"]
    integrity = report["integrity"]
    failed = [name for name, ok in integrity.items() if not ok]
    failed_checks = [c["case_id"] for c in report["checks"] if not c["ok"]]
    lines = [
        "# ClinVar",
        "",
        f"Source: `{SRC}` ({SRC.stat().st_size:,} bytes).",
        f"Alleles: `{ALLELES_PATH}` ({ALLELES_PATH.stat().st_size:,} bytes).",
        f"Dedup: `{DEDUP_PATH}` ({DEDUP_PATH.stat().st_size:,} bytes).",
        "GRCh38 only. GRCh37 is excluded. NCBI36 and `na` assemblies are excluded too.",
        "",
        "## Assemblies",
        "",
        _table(["Assembly", "Rows", "Disposition"], asm_rows),
        "",
        f"Bad-length records skipped: {stats['bad_len']:,}.",
        f"Allele scan and write: {stats['seconds']} seconds.",
        "",
        "## Alleles",
        "",
        f"- GRCh38 rows: {_num(report['allele_n'])}.",
        f"- Distinct VariationID: {_num(report['distinct_var'])}.",
        f"- Distinct AlleleID: {_num(report['distinct_allele'])}.",
        f"- Rows with a variant_key: {_num(report['allele_n'] - report['null_keys'])}.",
        f"- Rows with no VCF allele (null variant_key): {_num(report['null_keys'])}, "
        f"across {_num(report['null_key_alleles'])} AlleleIDs.",
        f"- Distinct variant_key: {_num(report['distinct_keys'])}.",
        f"- PhenotypeList max length: {_num(report['pheno_max'])} "
        f"(source {_num(stats['max_pheno'])}). Not truncated.",
        f"- PhenotypeIDS max length: {_num(report['ids_max'])} "
        f"(source {_num(stats['max_ids'])}). Kept.",
        f"- PhenotypeList longer than 500 characters: {_num(report['pheno_over_500'])}.",
        "",
        "Distinct VariationID is only modestly smaller than the GRCh38 row count: "
        f"{_num(report['allele_n'] - report['distinct_var'])} extra rows. "
        "Phenotypes are already packed into one row with `|`. The extra rows are the "
        "same AlleleID placed on more than one chromosome, not one row per phenotype.",
        "",
        "AlleleIDs with more than one GRCh38 row: "
        f"{_num(dup[0])}. Of those, {_num(dup[1])} differ by chromosome, "
        f"{_num(dup[2])} differ by clinical significance, "
        f"{_num(dup[3])} differ by PhenotypeList, "
        f"{_num(dup[4])} have more than one variant_key, "
        f"and {_num(dup[5])} have no variant_key.",
        f"variant_keys shared by more than one source row: {_num(report['multi_key'])}.",
        "variant_keys shared by more than one AlleleID: "
        f"{_num(report['multi_allele_same_key'])}.",
        "",
        "Chromosome values that do not canonicalize:",
        "",
        _table(
            ["chrom_raw", "Rows"],
            report["bad_chrom"] or [("(none)", 0)],
        ),
        "",
        "The variant key uses PositionVCF, ReferenceAlleleVCF, and AlternateAlleleVCF. "
        "A VCF allele of `na` is missing, not a sequence. Start, Stop, ReferenceAllele, "
        "and AlternateAllele stay in the source columns. `rsid` is `norm_rsid` of "
        "`RS# (dbSNP)`.",
        "",
        "## Review stars",
        "",
        "Match is case-insensitive on the full ReviewStatus phrase. "
        "4 practice guideline; 3 expert panel; 2 multiple submitters and no conflicts; "
        "1 single submitter or conflicting classifications; "
        "0 no assertion, no interpretation, or a blank status. "
        "Any other phrase, including `-`, is null.",
        "",
        _table(["review_stars", "Rows"], star_rows),
        "",
        "Unmatched phrases (review_stars null), by raw ReviewStatus:",
        "",
        _table(["ReviewStatus", "Rows"], unmatched_rows or [("(none)", 0)]),
        "",
        "Matched phrases:",
        "",
        _table(["review_stars", "ReviewStatus", "Rows"], matched_rows),
        "",
        "Dedup uses the max star in the group:",
        "",
        _table(
            ["review_stars", "Rows"],
            [
                ("null" if stars is None else str(stars), count)
                for stars, count in report["dedup_stars"]
            ],
        ),
        "",
        "## Dedup",
        "",
        f"- Rows: {_num(report['dedup_n'])}.",
        f"- Keyed by variant_key: {_num(report['dedup_keyed'])}.",
        f"- Keyed by AlleleID because there is no VCF allele: {_num(report['dedup_bare'])}. "
        "dedup_key is `allele:` plus AlleleID.",
        f"- Source rows covered by n_rows: {_num(report['sum_rows'])}.",
        f"- germline_pathogenic true: {_num(report['germ_n'])}.",
        "",
        "germline_pathogenic is true only when the aggregated OriginSimple class is "
        "germline (not somatic, not germline/somatic, not unknown), Oncogenicity has "
        "no call, and ClinicalSignificance has a pathogenic or likely pathogenic piece "
        "without a conflicting piece. Low-penetrance pathogenic labels count. "
        "Drug response alone does not. Oncogenicity and SomaticClinicalImpact stay "
        "in their own columns and are not copied into ClinicalSignificance.",
        "",
        "Pathogenic pieces that were not flagged: "
        f"{_num(report['dedup_counts'].get('excluded_origin', 0))} because origin was not germline-only, "
        f"{_num(report['dedup_counts'].get('excluded_onco', 0))} because of an oncogenicity call, "
        f"{_num(report['dedup_counts'].get('excluded_both', 0))} hit both, "
        f"{_num(report['dedup_counts'].get('excluded_onco_only', 0))} were germline origin blocked only by oncogenicity.",
        "",
        "## Clinical significance",
        "",
        _table(["ClinicalSignificance", "Rows"], report["sigs"]),
        "",
        "## Non-ACGTN VCF alleles",
        "",
        _table(
            ["ref", "alt", "Rows"],
            report["non_acgt"] or [("(none)", "", 0)],
        ),
        "",
        "## Truth set",
        "",
        "Checks read `warehouse/qc/truth_set.tsv` and the Parquet files.",
        "",
        _table(["Case", "Result", "Detail"], check_rows),
        "",
        "## Integrity",
        "",
    ]
    for name, ok in integrity.items():
        lines.append(f"- {name}: {'PASS' if ok else 'FAIL'}")
    lines.append("")
    if failed or failed_checks:
        lines.append(
            "FAILED: " + ", ".join(failed + failed_checks)
        )
    else:
        lines.append("All integrity checks and truth-set cases passed.")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    _check_rules()
    _check_dedup()
    print("rule checks passed", flush=True)
    stats = write_alleles()
    if stats["sep_hits"]:
        raise SystemExit(f"aggregation separator collides with source text: {stats['sep_hits']}")
    if stats["bad_len"]:
        raise SystemExit(f"bad length rows: {stats['bad_len']}")
    if stats["written"] != stats["grch38"]:
        raise SystemExit(f"wrote {stats['written']} but counted {stats['grch38']}")
    dedup_counts = write_dedup()
    report = gather(stats, dedup_counts)
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)
    QC_PATH.write_text(render(report), encoding="utf-8", newline="\n")
    print(f"wrote {QC_PATH}", flush=True)
    failed = [name for name, ok in report["integrity"].items() if not ok]
    failed += [c["case_id"] for c in report["checks"] if not c["ok"]]
    if failed:
        raise SystemExit("QC failed: " + ", ".join(failed))
    print(
        f"GRCh38 {report['allele_n']} VariationID {report['distinct_var']} "
        f"AlleleID {report['distinct_allele']} dedup {report['dedup_n']} "
        f"germline_pathogenic {report['germ_n']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
