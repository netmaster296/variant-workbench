# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Run the Wave 0 truth set against finished Wave 1 tables.

Writes warehouse/qc/truth_set_report.md. A table that is not built is
waiting, not a failure. This script does not modify any Parquet table.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import coords  # noqa: E402
import vcf_query  # noqa: E402

TRUTH = ROOT / "warehouse" / "qc" / "truth_set.tsv"
REPORT = ROOT / "warehouse" / "qc" / "truth_set_report.md"

ALLELES = ROOT / "warehouse" / "clinvar" / "alleles.parquet"
DEDUP = ROOT / "warehouse" / "clinvar" / "allele_dedup.parquet"
ASSOC = ROOT / "warehouse" / "gwas" / "associations.parquet"
LOCI = ROOT / "warehouse" / "gwas" / "loci.parquet"
GENE_CONSTRAINT = ROOT / "warehouse" / "constraint" / "gene_constraint.parquet"
CDS = ROOT / "warehouse" / "intervals" / "cds.parquet"
EGENES = ROOT / "warehouse" / "gtex" / "egenes.parquet"
SCORED = ROOT / "warehouse" / "mavedb" / "mapped_variants_scored.parquet"
SCORES = ROOT / "warehouse" / "mavedb" / "scores.parquet"
EXPERIMENTS = ROOT / "warehouse" / "mavedb" / "experiments.parquet"
MAPPED = ROOT / "warehouse" / "mavedb" / "mapped_variants.parquet"
CROSSWALK = ROOT / "warehouse" / "crosswalk" / "genes.parquet"
GENE_MART = ROOT / "warehouse" / "marts" / "gene_mart.parquet"
TIER_A = ROOT / "warehouse" / "variants" / "tier_a.parquet"
TIER_A_GENOME = ROOT / "warehouse" / "variants" / "tier_a_genome.parquet"
TIER_A_EXOME = ROOT / "warehouse" / "variants" / "tier_a_exome.parquet"
LEDGER = ROOT / "warehouse" / "ledger" / "evidence.parquet"

LOEUF_TOLERANCE = 0.001
AF_TOLERANCE = 1e-8
MISSING_TOKENS = {"", ".", "-", "na", "n/a", "none", "not provided"}
URN_RE = re.compile(r"urn:mavedb:[A-Za-z0-9._-]+(?:#[A-Za-z0-9._-]+)?")
ENSG_RE = re.compile(r"ENSG\d+")
SPAN_RE = re.compile(r"(chr(?:[0-9]+|X|Y|M)):(\d+)-(\d+)")

STATUSES = ("pass", "fail", "waiting")


class Check:
    def __init__(self, name: str, status: str, observed: str, file: Path | str):
        if status not in STATUSES:
            raise ValueError(status)
        self.name = name
        self.status = status
        self.observed = observed
        self.file = str(file)


def shown(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return "''" if value == "" else value
    if isinstance(value, float):
        return num_text(value)
    return str(value)


def num_text(value) -> str:
    number = float(value)
    if number.is_integer() and abs(number) < 1e15:
        return str(int(number))
    return f"{number:.6g}"


def close_enough(left, right, tolerance: float) -> bool:
    if left is None or right is None:
        return False
    return abs(float(left) - float(right)) <= tolerance


def md_cell(value: str) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def split_notes(text: str) -> list[str]:
    return [part.strip() for part in (text or "").split(";") if part.strip()]


def note_value(notes: str, prefix: str) -> str | None:
    for part in split_notes(notes):
        if part.startswith(prefix):
            return part[len(prefix) :].strip()
    return None


def constraint_spec(row: dict) -> dict:
    spec = {"gene": row["gene"], "transcript": None, "loeuf": None, "decile": None, "ensg": None}
    for part in split_notes(row["notes"]) + split_notes(row["expected"]):
        if part.startswith("LOEUF "):
            spec["loeuf"] = float(part.split()[1])
        elif part.startswith("decile "):
            spec["decile"] = int(float(part.split()[1]))
        elif part.startswith("transcript "):
            spec["transcript"] = part.split()[1]
        elif part.startswith("ENST"):
            spec["transcript"] = part.split()[0]
        elif part.startswith("ENSG"):
            spec["ensg"] = part.split()[0].split(".")[0]
    expected = (row["expected"] or "").strip()
    if expected.startswith("ENST"):
        spec["transcript"] = expected.split()[0]
    ensg = ENSG_RE.search(row["notes"] or "")
    if ensg and spec["ensg"] is None:
        spec["ensg"] = ensg.group(0)
    return spec


def locator_line(locator: str) -> int | None:
    parts = (locator or "").split()
    if len(parts) >= 2 and parts[0] == "line" and parts[1].isdigit():
        return int(parts[1])
    return None


def vcf_path(dataset: str, chrom: str) -> Path:
    if dataset == "genomes":
        name = f"gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz"
    elif dataset == "exomes":
        name = f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
    else:
        raise ValueError(dataset)
    return ROOT / "raw" / "gnomad" / dataset / name


def parse_b38(variant_id: str | None) -> str | None:
    if not variant_id:
        return None
    parts = variant_id.split("_")
    if len(parts) != 5 or parts[4] != "b38" or not parts[1].isdigit():
        return None
    return coords.variant_key(parts[0], parts[1], parts[2], parts[3])


def fetch_all(con: duckdb.DuckDBPyConnection, sql: str, params: list):
    cursor = con.execute(sql, params)
    names = [item[0] for item in cursor.description]
    return names, cursor.fetchall()


def fetch_dicts(con: duckdb.DuckDBPyConnection, sql: str, params: list) -> list[dict]:
    names, rows = fetch_all(con, sql, params)
    return [dict(zip(names, row)) for row in rows]


def one_int(con: duckdb.DuckDBPyConnection, sql: str, params: list) -> int:
    return int(con.execute(sql, params).fetchone()[0])


def table_ready(path: Path) -> Check | None:
    if path.exists():
        return None
    return Check(path.name, "waiting", "not built", path)


def allele_rows(con: duckdb.DuckDBPyConnection, variant_key: str) -> list[dict]:
    return fetch_dicts(
        con,
        """
        SELECT
          "GeneSymbol", "ReviewStatus", review_stars, ref, alt, variant_key, rsid,
          chrom, chrom_raw, pos,
          "ReferenceAllele", "AlternateAllele",
          "ReferenceAlleleVCF", "AlternateAlleleVCF",
          "Name", "ClinicalSignificance", "Oncogenicity", "OriginSimple", "Chromosome"
        FROM read_parquet(?)
        WHERE variant_key = ?
        """,
        [ALLELES.as_posix(), variant_key],
    )


def gnomad_hits(dataset: str, chrom: str, pos: int, ref: str, alt: str):
    path = vcf_path(dataset, chrom)
    if not path.exists():
        return path, None
    return path, vcf_query.query(path, chrom, pos, ref, alt)


def describe_hit(hit: dict, alt: str) -> str:
    alts = ",".join(hit["alts"])
    return f"ref {hit['ref']} alts {alts} id {shown(hit['id'])} contains {alt} {str(alt in hit['alts']).lower()}"


def check_vcf_alleles(dataset: str, chrom: str, pos: int, ref: str, alt: str) -> Check:
    path, hits = gnomad_hits(dataset, chrom, pos, ref, alt)
    label = f"gnomAD {dataset}"
    if hits is None:
        return Check(label, "fail", "miss: file not present", path)
    if not hits:
        return Check(label, "fail", f"miss: no record with ref {ref} alt {alt}", path)
    matched = [
        hit
        for hit in hits
        if hit["ref"] == ref and alt in hit["alts"]
    ]
    observed = "; ".join(describe_hit(hit, alt) for hit in hits)
    if len(matched) == len(hits) and matched:
        return Check(label, "pass", observed, path)
    return Check(label, "fail", observed, path)


def af_piece(raw: str | None, alts: list[str], alt: str) -> str | None:
    if raw is None or raw == ".":
        return None
    parts = raw.split(",")
    if len(parts) == 1:
        value = parts[0]
    elif alt in alts and len(parts) == len(alts):
        value = parts[alts.index(alt)]
    else:
        return None
    if value in {".", ""}:
        return None
    return value


def tier_a_waiting(purpose: str) -> Check:
    related = [TIER_A, TIER_A_GENOME, TIER_A_EXOME]
    present = [path.name for path in related if path.exists()]
    if present:
        extra = "also on disk: " + ", ".join(present)
    else:
        extra = "tier_a_genome.parquet and tier_a_exome.parquet are not built"
    return Check("Tier A", "waiting", f"not built; {extra}. {purpose}", TIER_A)


def check_coding_snp(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    key = row["variant_key"]
    checks = []
    if not ALLELES.exists():
        checks.append(Check("ClinVar alleles", "waiting", "not built", ALLELES))
    else:
        found = allele_rows(con, key)
        if len(found) != 1:
            observed = f"{len(found)} rows for {key}"
            checks.append(Check("ClinVar alleles", "fail", observed, ALLELES))
        else:
            item = found[0]
            problems = []
            if item["review_stars"] is None or int(item["review_stars"]) != 2:
                problems.append(f"review_stars {shown(item['review_stars'])}")
            if item["ref"] != "C" or item["alt"] != "T":
                problems.append(f"ref {shown(item['ref'])} alt {shown(item['alt'])}")
            if item["rsid"] != row["rsid"]:
                problems.append(f"rsid {shown(item['rsid'])}")
            if item["GeneSymbol"] != row["gene"]:
                problems.append(f"GeneSymbol {shown(item['GeneSymbol'])}")
            if item["ClinicalSignificance"] != "Pathogenic":
                problems.append(f"ClinicalSignificance {shown(item['ClinicalSignificance'])}")
            if item["variant_key"] != key:
                problems.append(f"variant_key {shown(item['variant_key'])}")
            observed = (
                f"variant_key {item['variant_key']}; GeneSymbol {shown(item['GeneSymbol'])}; "
                f"rsid {shown(item['rsid'])}; ref {shown(item['ref'])}; alt {shown(item['alt'])}; "
                f"review_stars {shown(item['review_stars'])}; "
                f"ReviewStatus {shown(item['ReviewStatus'])}; "
                f"ClinicalSignificance {shown(item['ClinicalSignificance'])}"
            )
            checks.append(Check("ClinVar alleles", "fail" if problems else "pass", observed, ALLELES))
    if not ASSOC.exists():
        checks.append(Check("GWAS rsid", "waiting", "not built", ASSOC))
    else:
        raw_id = row["rsid"][2:] if row["rsid"].startswith("rs") else row["rsid"]
        n_rsid, n_raw, n_prefixed = con.execute(
            """
            SELECT
              count(*) FILTER (WHERE rsid = ?),
              count(*) FILTER (WHERE "SNP_ID_CURRENT" = ?),
              count(*) FILTER (WHERE "SNP_ID_CURRENT" = ?)
            FROM read_parquet(?)
            """,
            [row["rsid"], raw_id, row["rsid"], ASSOC.as_posix()],
        ).fetchone()
        observed = (
            f"rsid {row['rsid']} rows {n_rsid}; "
            f"SNP_ID_CURRENT {raw_id} rows {n_raw}; "
            f"SNP_ID_CURRENT {row['rsid']} rows {n_prefixed}"
        )
        status = "pass" if n_rsid == 0 and n_raw == 0 and n_prefixed == 0 else "fail"
        checks.append(Check("GWAS rsid", status, observed, ASSOC))
    chrom = key.split(":", 1)[0]
    pos = int(row["pos"])
    print(f"vcf genomes {chrom}:{pos}", flush=True)
    checks.append(check_vcf_alleles("genomes", chrom, pos, row["ref"], row["alt"]))
    print(f"vcf exomes {chrom}:{pos}", flush=True)
    checks.append(check_vcf_alleles("exomes", chrom, pos, row["ref"], row["alt"]))
    if TIER_A.exists():
        checks.append(check_tier_a_key(con, key, row["ref"], row["alt"]))
    else:
        checks.append(
            tier_a_waiting("Genome and exome ref/alt above are from the sites VCFs, not from Tier A.")
        )
    return checks


def check_tier_a_key(con: duckdb.DuckDBPyConnection, key: str, ref: str, alt: str) -> Check:
    columns = pq.read_schema(TIER_A).names
    if "variant_key" not in columns:
        return Check("Tier A", "fail", "no variant_key column", TIER_A)
    names, rows = fetch_all(
        con,
        "SELECT * FROM read_parquet(?) WHERE variant_key = ? LIMIT 5",
        [TIER_A.as_posix(), key],
    )
    if not rows:
        return Check("Tier A", "fail", f"no row for {key}", TIER_A)
    item = dict(zip(names, rows[0]))
    ref_cols = [name for name in names if name.lower() in {"ref", "ref_genome", "genome_ref"}]
    alt_cols = [name for name in names if name.lower() in {"alt", "alt_genome", "genome_alt"}]
    bits = [f"rows {len(rows)}"]
    ok = True
    for name in ref_cols:
        bits.append(f"{name} {shown(item[name])}")
        if item[name] not in {None, ref}:
            ok = False
    for name in alt_cols:
        bits.append(f"{name} {shown(item[name])}")
        if item[name] not in {None, alt}:
            ok = False
    if not ref_cols and not alt_cols:
        bits.append("no ref/alt column")
        ok = False
    return Check("Tier A", "pass" if ok else "fail", "; ".join(bits), TIER_A)


def check_cftr(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    key = row["variant_key"]
    checks = []
    if not ALLELES.exists():
        checks.append(Check("ClinVar stars", "waiting", "not built", ALLELES))
    else:
        found = allele_rows(con, key)
        if len(found) != 1:
            checks.append(Check("ClinVar stars", "fail", f"{len(found)} rows for {key}", ALLELES))
        else:
            item = found[0]
            status_text = shown(item["ReviewStatus"])
            stars = None if item["review_stars"] is None else int(item["review_stars"])
            guideline = "practice guideline" in status_text.lower()
            problems = []
            if stars != 4:
                problems.append(f"review_stars {shown(stars)}")
            if not guideline:
                problems.append("ReviewStatus has no practice guideline")
            if item["ref"] != row["ref"] or item["alt"] != row["alt"]:
                problems.append(f"ref {shown(item['ref'])} alt {shown(item['alt'])}")
            if item["ClinicalSignificance"] != "Pathogenic":
                problems.append(f"ClinicalSignificance {shown(item['ClinicalSignificance'])}")
            observed = (
                f"variant_key {item['variant_key']}; GeneSymbol {shown(item['GeneSymbol'])}; "
                f"ref {shown(item['ref'])}; alt {shown(item['alt'])}; "
                f"review_stars {shown(stars)}; ReviewStatus {status_text}; "
                f"ClinicalSignificance {shown(item['ClinicalSignificance'])}"
            )
            checks.append(Check("ClinVar stars", "fail" if problems else "pass", observed, ALLELES))
    chrom = key.split(":", 1)[0]
    pos = int(row["pos"])
    print(f"vcf genomes {chrom}:{pos}", flush=True)
    checks.append(check_vcf_alleles("genomes", chrom, pos, row["ref"], row["alt"]))
    print(f"vcf exomes {chrom}:{pos}", flush=True)
    checks.append(check_vcf_alleles("exomes", chrom, pos, row["ref"], row["alt"]))
    return checks


def association_at_line(line_no: int) -> dict:
    index = line_no - 2
    columns = [
        "STUDY ACCESSION",
        "PUBMEDID",
        "SNPS",
        "CHR_POS",
        "CHR_ID",
        "MAPPED_TRAIT_URI",
        "MAPPED_TRAIT",
        "CONTEXT",
        "INTERGENIC",
        "SNP_ID_CURRENT",
        "MAPPED_GENE",
        "chrom",
        "pos",
        "rsid",
        "locus_key",
        "variant_key",
    ]
    batch = next(pq.ParquetFile(ASSOC).iter_batches(batch_size=index + 1, columns=columns))
    rows = batch.to_pylist()
    if index < 0 or index >= len(rows):
        raise RuntimeError(f"associations.parquet has no data index for source line {line_no}")
    return rows[index]


def check_gwas_intergenic(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    checks = []
    line_no = locator_line(row["source_locator"])
    if not ASSOC.exists():
        checks.append(Check("associations", "waiting", "not built", ASSOC))
        assoc = None
    elif line_no is None:
        checks.append(Check("associations", "fail", f"no source line in {row['source_locator']!r}", ASSOC))
        assoc = None
    else:
        assoc = association_at_line(line_no)
        problems = []
        if assoc["chrom"] != "chr22" or assoc["pos"] != int(row["pos"]):
            problems.append(f"chrom {shown(assoc['chrom'])} pos {shown(assoc['pos'])}")
        if assoc["rsid"] != row["rsid"]:
            problems.append(f"rsid {shown(assoc['rsid'])}")
        uri = assoc["MAPPED_TRAIT_URI"]
        if (uri or "").strip() != row["expected"].strip():
            problems.append(f"MAPPED_TRAIT_URI {shown(uri)}")
        observed = (
            f"source line {line_no}; chrom {shown(assoc['chrom'])}; pos {shown(assoc['pos'])}; "
            f"rsid {shown(assoc['rsid'])}; MAPPED_TRAIT_URI {shown(uri)}; "
            f"CONTEXT {shown(assoc['CONTEXT'])}; INTERGENIC {shown(assoc['INTERGENIC'])}; "
            f"locus_key {shown(assoc['locus_key'])}"
        )
        checks.append(Check("associations", "fail" if problems else "pass", observed, ASSOC))
    if not CDS.exists():
        checks.append(Check("CDS overlap", "waiting", "not built", CDS))
    else:
        chrom = "chr22"
        pos = int(row["pos"])
        if assoc is not None and assoc["chrom"] and assoc["pos"]:
            chrom = assoc["chrom"]
            pos = int(assoc["pos"])
        count = one_int(
            con,
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE chrom = ? AND start <= ? AND "end" >= ?
            """,
            [CDS.as_posix(), chrom, pos, pos],
        )
        observed = f"CDS intervals covering {chrom}:{pos}: {count}"
        if count:
            sample = fetch_dicts(
                con,
                """
                SELECT ensg, enst, start, "end", feature
                FROM read_parquet(?)
                WHERE chrom = ? AND start <= ? AND "end" >= ?
                LIMIT 3
                """,
                [CDS.as_posix(), chrom, pos, pos],
            )
            bits = [
                f"{item['ensg']} {item['enst']} {item['start']}-{item['end']} {item['feature']}"
                for item in sample
            ]
            observed += "; " + "; ".join(bits)
        checks.append(Check("CDS overlap", "fail" if count else "pass", observed, CDS))
    return checks


def check_gwas_empty(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    line_no = locator_line(row["source_locator"])
    if not ASSOC.exists():
        return [Check("association row", "waiting", "not built", ASSOC)]
    if line_no is None:
        return [Check("association row", "fail", f"locator {row['source_locator']!r} has no line", ASSOC)]
    assoc = association_at_line(line_no)
    study = note_value(row["notes"], "study")
    pubmed = note_value(row["notes"], "pubmed")
    uri = note_value(row["notes"], "trait uri")
    snps = note_value(row["notes"], "snps")
    problems = []
    if assoc["STUDY ACCESSION"] != study:
        problems.append(f"STUDY ACCESSION {shown(assoc['STUDY ACCESSION'])}")
    if assoc["PUBMEDID"] != pubmed:
        problems.append(f"PUBMEDID {shown(assoc['PUBMEDID'])}")
    if assoc["SNPS"] != snps:
        problems.append(f"SNPS {shown(assoc['SNPS'])}")
    if (assoc["MAPPED_TRAIT_URI"] or "").strip() != uri:
        problems.append(f"MAPPED_TRAIT_URI {shown(assoc['MAPPED_TRAIT_URI'])}")
    if assoc["pos"] is not None:
        problems.append(f"pos {shown(assoc['pos'])}")
    if not LOCI.exists():
        loci_observed = "loci.parquet is not built"
        loci_status = "waiting"
        loci_file = LOCI
    else:
        null_leads = one_int(
            con,
            "SELECT count(*) FROM read_parquet(?) WHERE lead_pos IS NULL OR chrom IS NULL",
            [LOCI.as_posix()],
        )
        if assoc["rsid"]:
            same_rsid = one_int(
                con,
                """
                SELECT count(*)
                FROM read_parquet(?)
                WHERE lead_rsid = ? AND trait_uri = ?
                """,
                [LOCI.as_posix(), assoc["rsid"], assoc["MAPPED_TRAIT_URI"]],
            )
            rsid_note = f"loci with lead_rsid {assoc['rsid']} and this trait {same_rsid}"
        else:
            rsid_note = "row rsid is null, so it cannot be a locus lead_rsid"
        loci_problems = []
        if assoc["locus_key"] is not None:
            loci_problems.append(f"locus_key {shown(assoc['locus_key'])}")
        if null_leads != 0:
            loci_problems.append(f"loci with null lead_pos or chrom {null_leads}")
        loci_observed = (
            f"no study column on loci; this row pos {shown(assoc['pos'])}; "
            f"locus_key {shown(assoc['locus_key'])}; "
            f"loci with null lead_pos or chrom {null_leads}; "
            f"{rsid_note}"
        )
        loci_status = "fail" if loci_problems else "pass"
        loci_file = LOCI
    observed = (
        f"source line {line_no}; STUDY ACCESSION {shown(assoc['STUDY ACCESSION'])}; "
        f"PUBMEDID {shown(assoc['PUBMEDID'])}; SNPS {shown(assoc['SNPS'])}; "
        f"CHR_POS {shown(assoc['CHR_POS'])}; CHR_ID {shown(assoc['CHR_ID'])}; "
        f"chrom {shown(assoc['chrom'])}; pos {shown(assoc['pos'])}; "
        f"rsid {shown(assoc['rsid'])}; locus_key {shown(assoc['locus_key'])}; "
        f"MAPPED_TRAIT_URI {shown(assoc['MAPPED_TRAIT_URI'])}"
    )
    return [
        Check("association row", "fail" if problems else "pass", observed, ASSOC),
        Check("absent from loci", loci_status, loci_observed, loci_file),
    ]


def check_mismatched(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    if not ALLELES.exists():
        return [Check("ClinVar key", "waiting", "not built", ALLELES)]
    found = allele_rows(con, row["variant_key"])
    if len(found) != 1:
        return [Check("ClinVar key", "fail", f"{len(found)} rows for {row['variant_key']}", ALLELES)]
    item = found[0]
    display_ref = shown(item["ReferenceAllele"])
    display_alt = shown(item["AlternateAllele"])
    name = item["Name"]
    expected_name = note_value(row["notes"], "name")
    problems = []
    if item["variant_key"] != row["variant_key"]:
        problems.append(f"variant_key {shown(item['variant_key'])}")
    if item["ref"] != row["ref"] or item["alt"] != row["alt"]:
        problems.append(f"ref {shown(item['ref'])} alt {shown(item['alt'])}")
    if str(item["ReferenceAlleleVCF"] or "").upper() != row["ref"]:
        problems.append(f"ReferenceAlleleVCF {shown(item['ReferenceAlleleVCF'])}")
    if str(item["AlternateAlleleVCF"] or "").upper() != row["alt"]:
        problems.append(f"AlternateAlleleVCF {shown(item['AlternateAlleleVCF'])}")
    if str(item["ReferenceAllele"] or "").strip().lower() != "na":
        problems.append(f"ReferenceAllele {display_ref}")
    if str(item["AlternateAllele"] or "").strip().lower() != "na":
        problems.append(f"AlternateAllele {display_alt}")
    if name == item["variant_key"]:
        problems.append("variant_key equals Name")
    if expected_name and name != expected_name:
        problems.append(f"Name {shown(name)}")
    observed = (
        f"variant_key {shown(item['variant_key'])}; ref {shown(item['ref'])}; alt {shown(item['alt'])}; "
        f"ReferenceAlleleVCF {shown(item['ReferenceAlleleVCF'])}; "
        f"AlternateAlleleVCF {shown(item['AlternateAlleleVCF'])}; "
        f"ReferenceAllele {display_ref}; AlternateAllele {display_alt}; "
        f"Name {shown(name)}"
    )
    return [Check("ClinVar key", "fail" if problems else "pass", observed, ALLELES)]


def check_af(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    chrom = row["variant_key"].split(":", 1)[0]
    pos = int(row["pos"])
    print(f"vcf genomes AF {chrom}:{pos}", flush=True)
    path, hits = gnomad_hits("genomes", chrom, pos, row["ref"], row["alt"])
    expected = {
        "AF_nfe": note_value(row["notes"], "AF_nfe"),
        "AF_afr": note_value(row["notes"], "AF_afr"),
    }
    if hits is None:
        vcf_check = Check("AF_nfe and AF_afr", "fail", "miss: file not present", path)
    elif not hits:
        vcf_check = Check(
            "AF_nfe and AF_afr",
            "fail",
            f"miss: no record with ref {row['ref']} alt {row['alt']}",
            path,
        )
    else:
        pieces = []
        problems = []
        for hit in hits:
            nfe = af_piece(hit["info"].get("AF_nfe"), hit["alts"], row["alt"])
            afr = af_piece(hit["info"].get("AF_afr"), hit["alts"], row["alt"])
            pieces.append(
                f"ref {hit['ref']} alts {','.join(hit['alts'])} "
                f"AF_nfe {shown(nfe)} AF_afr {shown(afr)} "
                f"raw AF_nfe {shown(hit['info'].get('AF_nfe'))} "
                f"raw AF_afr {shown(hit['info'].get('AF_afr'))}"
            )
            if nfe is None or afr is None:
                problems.append("a frequency is missing")
                continue
            try:
                differ = float(nfe) != float(afr)
            except ValueError:
                differ = nfe != afr
            if not differ:
                problems.append(f"AF_nfe {nfe} equals AF_afr {afr}")
            for field, want in expected.items():
                got = nfe if field == "AF_nfe" else afr
                if want is not None and not close_enough(got, want, AF_TOLERANCE):
                    problems.append(f"{field} {got} != truth-set {want}")
        vcf_check = Check(
            "AF_nfe and AF_afr",
            "fail" if problems else "pass",
            "; ".join(pieces),
            path,
        )
    checks = [vcf_check]
    if TIER_A.exists():
        checks.append(check_tier_a_af(con, row["variant_key"]))
    else:
        checks.append(
            tier_a_waiting("AF_nfe and AF_afr above are from the genome sites VCF, not stored on Tier A.")
        )
    return checks


def check_tier_a_af(con: duckdb.DuckDBPyConnection, key: str) -> Check:
    columns = pq.read_schema(TIER_A).names
    if "variant_key" not in columns:
        return Check("Tier A frequencies", "fail", "no variant_key column", TIER_A)
    nfe_cols = [name for name in columns if "AF_nfe" in name]
    afr_cols = [name for name in columns if "AF_afr" in name]
    if not nfe_cols or not afr_cols:
        shown_cols = ", ".join(columns[:40])
        return Check(
            "Tier A frequencies",
            "fail",
            f"AF_nfe columns {nfe_cols or 'none'}; AF_afr columns {afr_cols or 'none'}; columns {shown_cols}",
            TIER_A,
        )
    select = ", ".join(f'"{name}"' for name in nfe_cols + afr_cols)
    rows = fetch_dicts(
        con,
        f"SELECT {select} FROM read_parquet(?) WHERE variant_key = ? LIMIT 5",
        [TIER_A.as_posix(), key],
    )
    if not rows:
        return Check("Tier A frequencies", "fail", f"no row for {key}", TIER_A)
    item = rows[0]
    problems = []
    bits = []
    # Genome and exome frequencies are separate columns. A site absent from one
    # VCF keeps a null on that side. The case passes when each ancestry has a
    # stored value and no stored AF_nfe equals a stored AF_afr.
    for name in nfe_cols + afr_cols:
        bits.append(f"{name} {shown(item[name])}")
    if not any(item[name] is not None for name in nfe_cols):
        problems.append("AF_nfe missing")
    if not any(item[name] is not None for name in afr_cols):
        problems.append("AF_afr missing")
    for nfe_name in nfe_cols:
        for afr_name in afr_cols:
            if item[nfe_name] is None or item[afr_name] is None:
                continue
            try:
                same = float(item[nfe_name]) == float(item[afr_name])
            except (TypeError, ValueError):
                same = item[nfe_name] == item[afr_name]
            if same:
                problems.append(f"{nfe_name} equals {afr_name}")
    return Check("Tier A frequencies", "fail" if problems else "pass", "; ".join(bits), TIER_A)


def check_constraint(con: duckdb.DuckDBPyConnection, row: dict, mane_case: bool) -> list[Check]:
    spec = constraint_spec(row)
    checks = []
    if not GENE_CONSTRAINT.exists():
        checks.append(Check("gene_constraint", "waiting", "not built", GENE_CONSTRAINT))
    else:
        found = fetch_dicts(
            con,
            """
            SELECT gene, gene_id, transcript, bare_ensg, bare_enst, loeuf, loeuf_decile,
                   mane_select, "lof.oe_ci.upper" AS loeuf_upper
            FROM read_parquet(?)
            WHERE gene = ?
            """,
            [GENE_CONSTRAINT.as_posix(), spec["gene"]],
        )
        if len(found) != 1:
            checks.append(
                Check("gene_constraint", "fail", f"{len(found)} rows for gene {spec['gene']}", GENE_CONSTRAINT)
            )
        else:
            item = found[0]
            problems = []
            transcript = item["transcript"]
            if spec["transcript"] and transcript != spec["transcript"] and item["bare_enst"] != spec["transcript"]:
                problems.append(f"transcript {shown(transcript)}")
            if spec["loeuf"] is not None and not close_enough(item["loeuf"], spec["loeuf"], LOEUF_TOLERANCE):
                problems.append(f"LOEUF {shown(item['loeuf'])}")
            if spec["decile"] is not None and (
                item["loeuf_decile"] is None or int(float(item["loeuf_decile"])) != spec["decile"]
            ):
                problems.append(f"decile {shown(item['loeuf_decile'])}")
            if item["loeuf"] is not None and item["loeuf_upper"] is not None:
                if not close_enough(item["loeuf"], item["loeuf_upper"], LOEUF_TOLERANCE):
                    problems.append(f"loeuf {shown(item['loeuf'])} != lof.oe_ci.upper {shown(item['loeuf_upper'])}")
            if spec["ensg"] and item["bare_ensg"] != spec["ensg"] and item["gene_id"] != spec["ensg"]:
                problems.append(f"gene_id {shown(item['gene_id'])}")
            if mane_case and item["mane_select"] is not True:
                problems.append(f"mane_select {shown(item['mane_select'])}")
            observed = (
                f"gene {shown(item['gene'])}; gene_id {shown(item['gene_id'])}; "
                f"transcript {shown(transcript)}; bare_ensg {shown(item['bare_ensg'])}; "
                f"LOEUF {shown(item['loeuf'])}; decile {shown(item['loeuf_decile'])}; "
                f"mane_select {shown(item['mane_select'])}"
            )
            checks.append(Check("gene_constraint", "fail" if problems else "pass", observed, GENE_CONSTRAINT))
    checks.append(gene_mart_check(con, spec))
    return checks


def gene_mart_check(con: duckdb.DuckDBPyConnection, spec: dict) -> Check:
    if not GENE_MART.exists():
        return Check(
            "gene mart",
            "waiting",
            "not built; decile was not compared to gene_constraint",
            GENE_MART,
        )
    columns = pq.read_schema(GENE_MART).names
    if "symbol" not in columns:
        return Check("gene mart", "fail", "no symbol column", GENE_MART)
    decile_col = next((name for name in ("loeuf_decile", "LOEUF_decile", "loeuf_bin") if name in columns), None)
    loeuf_col = next((name for name in ("loeuf", "LOEUF") if name in columns), None)
    transcript_col = next(
        (name for name in ("transcript", "mane_enst", "mane_transcript", "constraint_transcript") if name in columns),
        None,
    )
    wanted = [name for name in (decile_col, loeuf_col, transcript_col) if name]
    if decile_col is None:
        return Check("gene mart", "fail", "no LOEUF decile column", GENE_MART)
    select = ", ".join(f'"{name}"' for name in ["symbol", *wanted])
    found = fetch_dicts(
        con,
        f"SELECT {select} FROM read_parquet(?) WHERE symbol = ?",
        [GENE_MART.as_posix(), spec["gene"]],
    )
    if len(found) != 1:
        return Check("gene mart", "fail", f"{len(found)} rows for {spec['gene']}", GENE_MART)
    item = found[0]
    problems = []
    if spec["decile"] is not None and (
        item[decile_col] is None or int(float(item[decile_col])) != spec["decile"]
    ):
        problems.append(f"{decile_col} {shown(item[decile_col])}")
    if loeuf_col and spec["loeuf"] is not None and not close_enough(item[loeuf_col], spec["loeuf"], LOEUF_TOLERANCE):
        problems.append(f"{loeuf_col} {shown(item[loeuf_col])}")
    if transcript_col and spec["transcript"] and item[transcript_col] not in {spec["transcript"], None}:
        if item[transcript_col] != spec["transcript"]:
            problems.append(f"{transcript_col} {shown(item[transcript_col])}")
    observed = "; ".join(f"{name} {shown(item[name])}" for name in ["symbol", *wanted])
    return Check("gene mart", "fail" if problems else "pass", observed, GENE_MART)


def check_gtex(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    if not EGENES.exists():
        return [Check("eGene keys", "waiting", "not built", EGENES)]
    tissue_match = re.search(r"/([^/]+)\.v8\.egenes", row["source_locator"])
    tissue = tissue_match.group(1) if tissue_match else None
    variant_id = note_value(row["notes"], "variant_id")
    if tissue is None or variant_id is None:
        return [Check("eGene keys", "fail", f"locator {row['source_locator']!r} notes {row['notes']!r}", EGENES)]
    found = fetch_dicts(
        con,
        """
        SELECT tissue, gene_id, ensg, variant_id, "chr", variant_pos, ref, alt, variant_key
        FROM read_parquet(?)
        WHERE tissue = ? AND gene_id = ? AND variant_id = ?
        """,
        [EGENES.as_posix(), tissue, row["gene"], variant_id],
    )
    if not found:
        return [
            Check(
                "eGene keys",
                "fail",
                f"no row tissue {tissue} gene_id {row['gene']} variant_id {variant_id}",
                EGENES,
            )
        ]
    problems = []
    bits = []
    for item in found:
        split_key = coords.variant_key(item["chr"], item["variant_pos"], item["ref"], item["alt"])
        id_key = parse_b38(item["variant_id"])
        bits.append(
            f"tissue {item['tissue']}; gene_id {shown(item['gene_id'])}; "
            f"chr {shown(item['chr'])}; variant_pos {shown(item['variant_pos'])}; "
            f"ref {shown(item['ref'])}; alt {shown(item['alt'])}; "
            f"split_key {shown(split_key)}; id_key {shown(id_key)}; "
            f"variant_key {shown(item['variant_key'])}"
        )
        if split_key != id_key:
            problems.append(f"split_key {split_key} != id_key {id_key}")
        if split_key != row["variant_key"] or id_key != row["variant_key"]:
            problems.append(f"keys are not {row['variant_key']}")
        if item["variant_key"] != row["variant_key"]:
            problems.append(f"stored variant_key {shown(item['variant_key'])}")
    observed = f"{len(found)} rows; " + " | ".join(bits)
    return [Check("eGene keys", "fail" if problems else "pass", observed, EGENES)]


def check_mavedb_mapped(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    if not SCORED.exists():
        return [Check("mapped score", "waiting", "not built", SCORED)]
    urns = URN_RE.findall(row["notes"] or "")
    variant_urns = [urn for urn in urns if "#" in urn]
    if len(variant_urns) != 1:
        return [Check("mapped score", "fail", f"notes urns {urns}", SCORED)]
    urn = variant_urns[0]
    found = fetch_dicts(
        con,
        """
        SELECT urn, score_set_urn, target_gene, chrom, pos, ref, alt, variant_key, score, score_text
        FROM read_parquet(?)
        WHERE urn = ?
        """,
        [SCORED.as_posix(), urn],
    )
    if not found:
        return [Check("mapped score", "fail", f"no row for {urn}", SCORED)]
    problems = []
    bits = []
    for item in found:
        bits.append(
            f"urn {item['urn']}; variant_key {shown(item['variant_key'])}; "
            f"chrom {shown(item['chrom'])}; pos {shown(item['pos'])}; "
            f"ref {shown(item['ref'])}; alt {shown(item['alt'])}; "
            f"score {shown(item['score'])}; score_text {shown(item['score_text'])}; "
            f"target_gene {shown(item['target_gene'])}"
        )
        if item["variant_key"] != row["variant_key"]:
            problems.append(f"variant_key {shown(item['variant_key'])}")
        if item["score"] is None:
            problems.append("score null")
    observed = f"{len(found)} rows; " + " | ".join(bits)
    return [Check("mapped score", "fail" if problems else "pass", observed, SCORED)]


def check_mavedb_unmapped(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    urns = URN_RE.findall(row["notes"] or "")
    score_urns = [urn for urn in urns if "#" not in urn]
    if len(score_urns) != 1:
        return [Check("unmapped score", "fail", f"notes urns {urns}", EXPERIMENTS)]
    urn = score_urns[0]
    checks = []
    if not EXPERIMENTS.exists():
        checks.append(Check("experiment", "waiting", "not built", EXPERIMENTS))
    else:
        found = fetch_dicts(
            con,
            """
            SELECT target_name, mapped_file_exists, mapping_state, score_set_title
            FROM read_parquet(?)
            WHERE score_set_urn = ?
            """,
            [EXPERIMENTS.as_posix(), urn],
        )
        names = sorted({shown(item["target_name"]) for item in found})
        flags = sorted({shown(item["mapped_file_exists"]) for item in found})
        titles = sorted({shown(item["score_set_title"]) for item in found})
        states = sorted({shown(item["mapping_state"]) for item in found})
        problems = []
        if not found:
            problems.append("no experiment row")
        if names != [row["gene"]]:
            problems.append(f"target_name {names}")
        if flags != ["false"]:
            problems.append(f"mapped_file_exists {flags}")
        observed = (
            f"score_set_urn {urn}; rows {len(found)}; target_name {names}; "
            f"mapped_file_exists {flags}; mapping_state {states}; score_set_title {titles}"
        )
        checks.append(Check("experiment", "fail" if problems else "pass", observed, EXPERIMENTS))
    if not SCORES.exists():
        checks.append(Check("scores", "waiting", "not built", SCORES))
    else:
        n_rows, n_score, n_text = con.execute(
            """
            SELECT count(*), count(score), count(score_text)
            FROM read_parquet(?)
            WHERE score_set_urn = ?
            """,
            [SCORES.as_posix(), urn],
        ).fetchone()
        observed = f"score rows {n_rows}; non-null score {n_score}; non-null score_text {n_text}"
        status = "pass" if n_rows > 0 and n_score > 0 else "fail"
        checks.append(Check("scores", status, observed, SCORES))
    if not MAPPED.exists():
        checks.append(Check("coordinates", "waiting", "not built", MAPPED))
    else:
        n_rows, n_chrom, n_pos, n_key = con.execute(
            """
            SELECT count(*), count(chrom), count(pos), count(variant_key)
            FROM read_parquet(?)
            WHERE score_set_urn = ?
            """,
            [MAPPED.as_posix(), urn],
        ).fetchone()
        observed = (
            f"mapped rows {n_rows}; non-null chrom {n_chrom}; "
            f"non-null pos {n_pos}; non-null variant_key {n_key}"
        )
        status = "pass" if n_chrom == 0 and n_pos == 0 and n_key == 0 else "fail"
        checks.append(Check("coordinates", status, observed, MAPPED))
    return checks


def check_crosswalk(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    checks = []
    key = row["variant_key"]
    ensg_match = ENSG_RE.search(row["notes"] or "")
    ensg = ensg_match.group(0) if ensg_match else None
    span = SPAN_RE.search(row["notes"] or "")
    if not ALLELES.exists():
        checks.append(Check("ClinVar chromosome", "waiting", "not built", ALLELES))
        allele = None
    else:
        found = allele_rows(con, key)
        if len(found) != 1:
            checks.append(Check("ClinVar chromosome", "fail", f"{len(found)} rows for {key}", ALLELES))
            allele = None
        else:
            allele = found[0]
            problems = []
            if str(allele["chrom_raw"]) != "1":
                problems.append(f"chrom_raw {shown(allele['chrom_raw'])}")
            if allele["chrom"] != "chr1":
                problems.append(f"chrom {shown(allele['chrom'])}")
            if allele["Chromosome"] != "1":
                problems.append(f"Chromosome {shown(allele['Chromosome'])}")
            if allele["GeneSymbol"] != row["gene"]:
                problems.append(f"GeneSymbol {shown(allele['GeneSymbol'])}")
            if allele["variant_key"] != key:
                problems.append(f"variant_key {shown(allele['variant_key'])}")
            observed = (
                f"variant_key {shown(allele['variant_key'])}; Chromosome {shown(allele['Chromosome'])}; "
                f"chrom_raw {shown(allele['chrom_raw'])}; chrom {shown(allele['chrom'])}; "
                f"GeneSymbol {shown(allele['GeneSymbol'])}; pos {shown(allele['pos'])}; "
                f"ref {shown(allele['ref'])}; alt {shown(allele['alt'])}"
            )
            checks.append(Check("ClinVar chromosome", "fail" if problems else "pass", observed, ALLELES))
    if not CROSSWALK.exists():
        checks.append(Check("crosswalk gene", "waiting", "not built", CROSSWALK))
    elif ensg is None:
        checks.append(Check("crosswalk gene", "fail", "truth-set notes have no ENSG", CROSSWALK))
    else:
        found = fetch_dicts(
            con,
            """
            SELECT ensg, ensg_version, symbol, chrom, start, "end", biotype
            FROM read_parquet(?)
            WHERE ensg = ?
            """,
            [CROSSWALK.as_posix(), ensg],
        )
        if len(found) != 1:
            checks.append(Check("crosswalk gene", "fail", f"{len(found)} rows for {ensg}", CROSSWALK))
        else:
            item = found[0]
            pos = int(row["pos"]) if row["pos"] else None
            if allele is not None and allele["pos"] is not None:
                pos = int(allele["pos"])
            inside = (
                pos is not None
                and item["start"] is not None
                and item["end"] is not None
                and int(item["start"]) <= pos <= int(item["end"])
            )
            problems = []
            if item["chrom"] != "chr1":
                problems.append(f"chrom {shown(item['chrom'])}")
            if item["symbol"] != row["gene"]:
                problems.append(f"symbol {shown(item['symbol'])}")
            if item["ensg"] != ensg:
                problems.append(f"ensg {shown(item['ensg'])}")
            if not inside:
                problems.append(f"pos {shown(pos)} outside {shown(item['start'])}-{shown(item['end'])}")
            if span and (int(item["start"]) != int(span.group(2)) or int(item["end"]) != int(span.group(3))):
                problems.append(
                    f"span {item['start']}-{item['end']} != note {span.group(2)}-{span.group(3)}"
                )
            if allele is not None and allele["chrom"] != item["chrom"]:
                problems.append(f"allele chrom {shown(allele['chrom'])} != gene chrom {shown(item['chrom'])}")
            observed = (
                f"ensg {shown(item['ensg'])}; ensg_version {shown(item['ensg_version'])}; "
                f"symbol {shown(item['symbol'])}; chrom {shown(item['chrom'])}; "
                f"start {shown(item['start'])}; end {shown(item['end'])}; "
                f"biotype {shown(item['biotype'])}; pos inside span {str(inside).lower()}"
            )
            checks.append(Check("crosswalk gene", "fail" if problems else "pass", observed, CROSSWALK))
    return checks


def is_missing_token(value) -> bool:
    if value is None:
        return True
    return str(value).strip().lower() in MISSING_TOKENS


def check_somatic(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    if not DEDUP.exists():
        return [Check("allele_dedup", "waiting", "not built", DEDUP)]
    columns = pq.read_schema(DEDUP).names
    if "Oncogenicity" not in columns or "ClinicalSignificance" not in columns:
        return [Check("allele_dedup", "fail", f"columns {', '.join(columns)}", DEDUP)]
    if "germline_pathogenic" not in columns:
        return [Check("allele_dedup", "fail", "germline_pathogenic column absent", DEDUP)]
    found = fetch_dicts(
        con,
        """
        SELECT variant_key, "GeneSymbol", "ClinicalSignificance", "Oncogenicity",
               "OriginSimple", germline_somatic, germline_pathogenic, review_stars
        FROM read_parquet(?)
        WHERE variant_key = ?
        """,
        [DEDUP.as_posix(), row["variant_key"]],
    )
    if len(found) != 1:
        return [Check("allele_dedup", "fail", f"{len(found)} rows for {row['variant_key']}", DEDUP)]
    item = found[0]
    problems = []
    if item["germline_pathogenic"] is not False:
        problems.append(f"germline_pathogenic {shown(item['germline_pathogenic'])}")
    if is_missing_token(item["Oncogenicity"]):
        problems.append(f"Oncogenicity {shown(item['Oncogenicity'])}")
    if item["Oncogenicity"] == item["ClinicalSignificance"]:
        problems.append("Oncogenicity equals ClinicalSignificance")
    if item["ClinicalSignificance"] != "Pathogenic":
        problems.append(f"ClinicalSignificance {shown(item['ClinicalSignificance'])}")
    if item["Oncogenicity"] != "Likely oncogenic":
        problems.append(f"Oncogenicity {shown(item['Oncogenicity'])}")
    if item["GeneSymbol"] != row["gene"]:
        problems.append(f"GeneSymbol {shown(item['GeneSymbol'])}")
    observed = (
        f"variant_key {shown(item['variant_key'])}; GeneSymbol {shown(item['GeneSymbol'])}; "
        f"germline_pathogenic {shown(item['germline_pathogenic'])}; "
        f"Oncogenicity {shown(item['Oncogenicity'])}; "
        f"ClinicalSignificance {shown(item['ClinicalSignificance'])}; "
        f"OriginSimple {shown(item['OriginSimple'])}; "
        f"germline_somatic {shown(item['germline_somatic'])}"
    )
    return [Check("allele_dedup", "fail" if problems else "pass", observed, DEDUP)]


def check_mane(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    return check_constraint(con, row, mane_case=True)


def check_loeuf(con: duckdb.DuckDBPyConnection, row: dict) -> list[Check]:
    return check_constraint(con, row, mane_case=False)


HANDLERS = {
    "coding_snp": check_coding_snp,
    "cftr_indel": check_cftr,
    "gwas_intergenic": check_gwas_intergenic,
    "gwas_empty_pos": check_gwas_empty,
    "mismatched_ref": check_mismatched,
    "af_pop_diff": check_af,
    "mane_select": check_mane,
    "loeuf_low_1": check_loeuf,
    "loeuf_low_2": check_loeuf,
    "loeuf_low_3": check_loeuf,
    "loeuf_near_1_1": check_loeuf,
    "loeuf_near_1_2": check_loeuf,
    "loeuf_near_1_3": check_loeuf,
    "gtex_egene": check_gtex,
    "mavedb_mapped": check_mavedb_mapped,
    "mavedb_unmapped": check_mavedb_unmapped,
    "chrom_crosswalk": check_crosswalk,
    "somatic_or_oncogenicity": check_somatic,
}


def case_status(checks: list[Check]) -> str:
    if any(item.status == "fail" for item in checks):
        return "fail"
    if any(item.status == "waiting" for item in checks):
        return "waiting"
    return "pass"


def load_cases() -> list[dict]:
    with TRUTH.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def absent_lines() -> list[str]:
    lines = []
    for path in (GENE_MART, TIER_A, TIER_A_GENOME, TIER_A_EXOME, LEDGER):
        state = "absent" if not path.exists() else "present"
        lines.append(f"- `{path}` — {state}")
    return lines


def render(cases: list[dict], results: list[tuple[dict, list[Check]]]) -> str:
    statuses = [case_status(checks) for _row, checks in results]
    counts = {name: statuses.count(name) for name in STATUSES}
    lines = [
        "# Truth-set report",
        "",
        "| status | cases |",
        "| --- | ---: |",
        f"| pass | {counts['pass']} |",
        f"| fail | {counts['fail']} |",
        f"| waiting | {counts['waiting']} |",
        "",
        f"{len(results)} cases from `{TRUTH}`.",
        "A case is fail when a file that exists does not match.",
        "A case is waiting when every check against an existing file matched and a table that case still needs is not built.",
        "Waiting is not a failure. No Parquet table was edited.",
        "",
    ]
    failed = [(row, checks) for row, checks in results if case_status(checks) == "fail"]
    if not failed:
        lines.append("Failing cases: none.")
    else:
        lines.append("Failing cases:")
        lines.append("")
        for row, checks in failed:
            for item in checks:
                if item.status == "fail":
                    lines.append(f"- `{row['case_id']}`. `{Path(item.file).name}` observed {item.observed}.")
    lines.append("")
    waiting = [row["case_id"] for row, checks in results if case_status(checks) == "waiting"]
    if waiting:
        lines.append("Waiting cases: " + ", ".join(f"`{case_id}`" for case_id in waiting) + ".")
    else:
        lines.append("Waiting cases: none.")
    lines.append("")
    lines.append("## Not built")
    lines.append("")
    lines.extend(absent_lines())
    lines.append("")
    lines.append(
        f"`{LEDGER.name}` is not a case result. The somatic case is read from `allele_dedup.parquet`."
    )
    lines.append("")
    lines.append("## Cases")
    lines.append("")
    lines.append("| case | status |")
    lines.append("| --- | --- |")
    for row, checks in results:
        lines.append(f"| `{row['case_id']}` | {case_status(checks)} |")
    lines.append("")
    for row, checks in results:
        status = case_status(checks)
        lines.append(f"### {row['case_id']}: {status}")
        lines.append("")
        lines.append(f"Requirement: {row['requirement']}")
        lines.append("")
        lines.append(f"Truth-set expected: {row['expected'] or '—'}")
        lines.append("")
        if status == "fail":
            for item in checks:
                if item.status == "fail":
                    lines.append(f"Fail. `{Path(item.file).name}` observed {item.observed}.")
                    lines.append("")
        lines.append("| check | status | observed | file |")
        lines.append("| --- | --- | --- | --- |")
        for item in checks:
            lines.append(
                f"| {md_cell(item.name)} | {item.status} | {md_cell(item.observed)} | `{item.file}` |"
            )
        lines.append("")
    text = "\n".join(lines)
    if not text.endswith("\n"):
        text += "\n"
    return text


def main() -> None:
    cases = load_cases()
    known = set(HANDLERS)
    found_ids = [row["case_id"] for row in cases]
    missing = [case_id for case_id in HANDLERS if case_id not in found_ids]
    extra = [case_id for case_id in found_ids if case_id not in known]
    if missing or extra or len(cases) != 18:
        print(f"truth set rows {len(cases)} missing {missing} extra {extra}", flush=True)
    con = duckdb.connect()
    con.execute("SET threads = 4")
    results = []
    try:
        for row in cases:
            print(f"case {row['case_id']}", flush=True)
            handler = HANDLERS.get(row["case_id"])
            if handler is None:
                checks = [Check("harness", "fail", "no check implemented", TRUTH)]
            else:
                try:
                    checks = handler(con, row)
                except Exception as exc:
                    checks = [Check("error", "fail", f"{type(exc).__name__}: {exc}", TRUTH)]
            print(f"  {case_status(checks)}", flush=True)
            results.append((row, checks))
    finally:
        con.close()
    report = render(cases, results)
    REPORT.write_text(report, encoding="utf-8", newline="\n")
    counts = {name: 0 for name in STATUSES}
    for _row, checks in results:
        counts[case_status(checks)] += 1
    print(
        f"wrote {REPORT} pass {counts['pass']} fail {counts['fail']} waiting {counts['waiting']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
