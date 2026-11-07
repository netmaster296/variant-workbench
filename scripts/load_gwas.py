# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load the GWAS Catalog association index into Parquet.

The catalog TSV is a literature index of reported associations, not
genome-wide summary statistics. This script does not compute a score.
`variant_key` stays null: the file has no ref or alt, and the risk-allele
token is not a genotype.
"""

from __future__ import annotations

import bisect
import csv
import gc
import io
import math
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from coords import canon_chrom, norm_rsid

ZIP_PATH = ROOT / "raw" / "gwas" / "gwas-catalog-associations_ontology-annotated-full.zip"
MEMBER = "gwas-catalog-download-associations-alt-full.tsv"
GWAS_DIR = ROOT / "warehouse" / "gwas"
QC_PATH = ROOT / "warehouse" / "qc" / "gwas.md"

EXPECTED_ROWS = 1_205_644
WINDOW_BP = 1_000_000
BATCH_ROWS = 100_000

TRUTH_URI = "http://www.ebi.ac.uk/efo/EFO_0004348"
TRUTH_CHROM = "chr22"
TRUTH_POS = 45_604_058
TRUTH_RSID = "rs4823308"
EMPTY_STUDY = "GCST004604"
EMPTY_PUBMED = "27863252"
EMPTY_SNPS = "chr22:37462936"

COLUMNS = [
    "DATE ADDED TO CATALOG",
    "PUBMEDID",
    "FIRST AUTHOR",
    "DATE",
    "JOURNAL",
    "LINK",
    "STUDY",
    "DISEASE/TRAIT",
    "INITIAL SAMPLE SIZE",
    "REPLICATION SAMPLE SIZE",
    "REGION",
    "CHR_ID",
    "CHR_POS",
    "REPORTED GENE(S)",
    "MAPPED_GENE",
    "UPSTREAM_GENE_ID",
    "DOWNSTREAM_GENE_ID",
    "SNP_GENE_IDS",
    "UPSTREAM_GENE_DISTANCE",
    "DOWNSTREAM_GENE_DISTANCE",
    "STRONGEST SNP-RISK ALLELE",
    "SNPS",
    "MERGED",
    "SNP_ID_CURRENT",
    "CONTEXT",
    "INTERGENIC",
    "RISK ALLELE FREQUENCY",
    "P-VALUE",
    "PVALUE_MLOG",
    "P-VALUE (TEXT)",
    "OR or BETA",
    "95% CI (TEXT)",
    "PLATFORM [SNPS PASSING QC]",
    "CNV",
    "MAPPED_TRAIT",
    "MAPPED_TRAIT_URI",
    "STUDY ACCESSION",
    "GENOTYPING TECHNOLOGY",
]

ADDED = [
    "chrom",
    "chrom_raw",
    "pos",
    "rsid",
    "p_value_float",
    "effect_float",
    "risk_allele",
    "variant_key",
    "locus_key",
]

# Clump tuple layout.
POS, PVAL, RSID, GENE, LINE = range(5)

_PLAIN_NUM = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
_RISK_SPLIT = re.compile(r"\s*;\s*|\s+x\s+", re.IGNORECASE)
_MISSING_ALLELE = {"?", ".", "-", "NA", "na", "N/A", "n/a"}

MAPPED_GENES_TYPE = pa.struct(
    [
        pa.field("raw", pa.string()),
        pa.field("tokens", pa.list_(pa.string())),
    ]
)

ASSOC_SCHEMA = pa.schema(
    [pa.field(name, pa.string()) for name in COLUMNS]
    + [
        pa.field("chrom", pa.string()),
        pa.field("chrom_raw", pa.string()),
        pa.field("pos", pa.int64()),
        pa.field("rsid", pa.string()),
        pa.field("p_value_float", pa.float64()),
        pa.field("effect_float", pa.float64()),
        pa.field("risk_allele", pa.string()),
        pa.field("variant_key", pa.string()),
        pa.field("locus_key", pa.string()),
    ]
)

TRAIT_SCHEMA = pa.schema(
    [
        pa.field("trait_uri", pa.string()),
        pa.field("mapped_trait", pa.string()),
        pa.field("n_associations", pa.int64()),
        pa.field("n_author_traits", pa.int64()),
    ]
)

LOCI_SCHEMA = pa.schema(
    [
        pa.field("trait_uri", pa.string()),
        pa.field("chrom", pa.string()),
        pa.field("lead_pos", pa.int64()),
        pa.field("lead_rsid", pa.string()),
        pa.field("lead_p", pa.float64()),
        pa.field("n_associations", pa.int64()),
        pa.field("start", pa.int64()),
        pa.field("end", pa.int64()),
        pa.field("mapped_genes", MAPPED_GENES_TYPE),
    ]
)


class TraitAcc:
    __slots__ = ("n", "names", "authors")

    def __init__(self) -> None:
        self.n = 0
        self.names: Counter[str] = Counter()
        self.authors: set[str] = set()


def _is_literal_zero(text: str) -> bool:
    body = text[1:] if text[:1] in "+-" else text
    coeff = body.lower().split("e", 1)[0].replace(".", "")
    return coeff != "" and coeff.strip("0") == ""


def _parse_finite(text: str) -> float | None:
    if not text or not _PLAIN_NUM.match(text):
        return None
    value = float(text)
    if not math.isfinite(value):
        return None
    return value


def parse_p_value(cell: str) -> tuple[float | None, str | None]:
    """Return (p, flag). Flag is 'underflow', 'negative', or None.

    A p-value that underflows float64 stays the most significant rank.
    It is stored as the smallest positive float, not as zero.
    """
    text = cell.strip()
    value = _parse_finite(text)
    if value is None:
        return None, None
    if value < 0:
        return None, "negative"
    if value == 0.0 and not _is_literal_zero(text):
        return sys.float_info.min, "underflow"
    return value, None


def parse_effect(cell: str) -> tuple[float | None, str | None]:
    """Float only when the whole cell is a plain finite number."""
    text = cell.strip()
    value = _parse_finite(text)
    if value is None:
        return None, None
    if value == 0.0 and not _is_literal_zero(text):
        return None, "underflow"
    return value, None


def parse_pos(cell: str) -> int | None:
    text = cell.strip()
    if not text or not text.isdigit():
        return None
    value = int(text)
    if value < 1:
        return None
    return value


def parse_risk_allele(cell: str) -> str | None:
    """Allele token from `rs128494-C`. Not a REF/ALT pair.

    Multiple SNP tokens (`rs1-A; rs2-G` or `rs1-A x rs2-G`) stay in order,
    joined by ';'. `?` and other missing markers become null.
    """
    text = cell.strip()
    if not text:
        return None
    alleles: list[str] = []
    for part in _RISK_SPLIT.split(text):
        part = part.strip()
        if "-" not in part:
            continue
        allele = part.rsplit("-", 1)[1].strip()
        if not allele or allele in _MISSING_ALLELE:
            continue
        alleles.append(allele)
    if not alleles:
        return None
    return ";".join(alleles)


def split_uris(cell: str) -> list[str]:
    if not cell or not cell.strip():
        return []
    return [part.strip() for part in cell.split(",") if part.strip()]


def trait_name_pairs(mapped_trait: str, uris: list[str]) -> list[tuple[str | None, str | None]]:
    """Pair each URI with a mapped name. Names that contain commas are not split.

    A single URI keeps the whole MAPPED_TRAIT string. Several URIs are paired
    only when a comma split yields the same count. Otherwise the name is left
    unset for this row and filled from rows that are unambiguous.
    """
    if not uris:
        return [(None, None)]
    if len(uris) == 1:
        name = mapped_trait.strip() or None
        return [(uris[0], name)]
    parts = [part.strip() for part in mapped_trait.split(",")]
    if len(parts) == len(uris):
        return [(uri, part or None) for uri, part in zip(uris, parts)]
    return [(uri, None) for uri in uris]


def gene_tokens(raw: str | None) -> list[str]:
    """Split MAPPED_GENE on commas only. Drop empty tokens and a bare '-'."""
    if not raw:
        return []
    tokens: list[str] = []
    for part in raw.split(","):
        token = part.strip()
        if token and token != "-":
            tokens.append(token)
    return tokens


def is_cnv(cell: str) -> bool:
    return cell.strip().upper() == "Y"


def chrom_rank(chrom: str) -> int:
    body = chrom[3:] if chrom.startswith("chr") else chrom
    if body.isdigit():
        return int(body)
    if body == "X":
        return 23
    if body == "Y":
        return 24
    if body == "M":
        return 25
    return 99


def _p_sort(value: float | None) -> float:
    return math.inf if value is None else value


def clump_hits(hits: list[tuple]) -> list[list[int]]:
    """Lead-centered clumps. Every member is within 1 Mb of the min-p lead.

    Sorting by position only locates the window. The lead is the unassigned
    association with the smallest p-value, not the leftmost association, and
    not a chain of neighbors each within 1 Mb of the previous one.
    """
    n = len(hits)
    if n == 0:
        return []
    hits.sort(key=lambda hit: (hit[POS], hit[LINE]))
    positions = [hit[POS] for hit in hits]
    order = sorted(
        range(n),
        key=lambda i: (
            _p_sort(hits[i][PVAL]),
            hits[i][POS],
            hits[i][RSID] or "",
            hits[i][LINE],
        ),
    )
    used = bytearray(n)
    clumps: list[list[int]] = []
    for lead_i in order:
        if used[lead_i]:
            continue
        lead_pos = positions[lead_i]
        lo = _bisect_left(positions, lead_pos - WINDOW_BP)
        hi = _bisect_right(positions, lead_pos + WINDOW_BP)
        members: list[int] = []
        for j in range(lo, hi):
            if used[j]:
                continue
            distance = positions[j] - lead_pos
            if distance < 0:
                distance = -distance
            if distance > WINDOW_BP:
                continue
            used[j] = 1
            members.append(j)
        if not members:
            raise RuntimeError("clump lead was not within its own window")
        clumps.append(members)
    if any(flag == 0 for flag in used):
        raise RuntimeError("an eligible association was left outside every clump")
    lead_positions: list[int] = []
    for members in clumps:
        lead_i = _lead_index(hits, members)
        lead_positions.append(positions[lead_i])
        for j in members:
            distance = positions[j] - positions[lead_i]
            if distance < 0:
                distance = -distance
            if distance > WINDOW_BP:
                raise RuntimeError("clump member is more than 1 Mb from its lead")
    lead_positions.sort()
    for left, right in zip(lead_positions, lead_positions[1:]):
        if right - left <= WINDOW_BP:
            raise RuntimeError("two leads on one chromosome are within 1 Mb")
    return clumps


def _lead_index(hits: list[tuple], members: list[int]) -> int:
    return min(
        members,
        key=lambda i: (
            _p_sort(hits[i][PVAL]),
            hits[i][POS],
            hits[i][RSID] or "",
            hits[i][LINE],
        ),
    )


def _bisect_left(positions: list[int], target: int) -> int:
    return bisect.bisect_left(positions, target)


def _bisect_right(positions: list[int], target: int) -> int:
    return bisect.bisect_right(positions, target)


def _check() -> None:
    assert len(COLUMNS) == 38
    assert len(ASSOC_SCHEMA) == 47
    assert canon_chrom("22") == "chr22"
    assert canon_chrom("") is None
    assert canon_chrom("6 x 6") is None
    assert norm_rsid("4823308") == "rs4823308"
    assert norm_rsid("7658266b") is None
    assert norm_rsid("") is None
    assert parse_pos("") is None
    assert parse_pos("45604058") == 45604058
    assert parse_pos("70610804;70611059") is None
    assert parse_pos("0") is None
    assert parse_p_value("2E-16")[0] == 2e-16
    assert parse_p_value("")[0] is None
    assert parse_p_value("-1")[0] is None
    assert parse_effect("0.03481525")[0] == 0.03481525
    assert parse_effect("Infinity")[0] is None
    assert parse_effect("")[0] is None
    assert parse_effect("1.5 (approx)")[0] is None
    assert parse_risk_allele("rs128494-C") == "C"
    assert parse_risk_allele("rs145611121-GC") == "GC"
    assert parse_risk_allele("rs4823308-A") == "A"
    assert parse_risk_allele("chr22:37462936-G") == "G"
    assert parse_risk_allele("rs12140153-?") is None
    assert parse_risk_allele("rs147903261-G x rs6493454-C") == "G;C"
    assert parse_risk_allele("rs17005943-T; rs17005947-G") == "T;G"
    assert parse_risk_allele("rs1-?; rs2-A") == "A"
    assert parse_risk_allele("") is None
    assert gene_tokens("CLDN14, LNCTSI") == ["CLDN14", "LNCTSI"]
    assert gene_tokens("FBLN1 - LINC01589") == ["FBLN1 - LINC01589"]
    assert gene_tokens("BRCA1, -, NBR2") == ["BRCA1", "NBR2"]
    assert gene_tokens("HLA-DQB1") == ["HLA-DQB1"]
    assert gene_tokens("-") == []
    assert gene_tokens("") == []
    assert is_cnv("Y") and is_cnv(" y ")
    assert not is_cnv("N") and not is_cnv("")
    pairs = trait_name_pairs(
        "osteoarthritis, hand",
        ["http://purl.obolibrary.org/obo/MONDO_0006632"],
    )
    assert pairs == [("http://purl.obolibrary.org/obo/MONDO_0006632", "osteoarthritis, hand")]
    aligned = trait_name_pairs(
        "neutrophil count, eosinophil count",
        [
            "http://www.ebi.ac.uk/efo/EFO_0004833",
            "http://www.ebi.ac.uk/efo/EFO_0004842",
        ],
    )
    assert aligned[0][1] == "neutrophil count"
    assert aligned[1][1] == "eosinophil count"
    unaligned = trait_name_pairs(
        "osteoarthritis, hip, osteoarthritis, knee",
        [
            "http://purl.obolibrary.org/obo/MONDO_0006629",
            "http://purl.obolibrary.org/obo/MONDO_0005416",
        ],
    )
    assert unaligned[0][1] is None and unaligned[1][1] is None

    # pos 0 is outside 1 Mb of the strong lead; pos 100 is inside it.
    hits = [
        (0, 0.4, None, None, 1),
        (100, 0.3, None, None, 2),
        (1_000_100, 1e-30, None, "BRCA1, NBR2", 3),
    ]
    clumps = clump_hits(hits)
    assert len(clumps) == 2
    groups = [sorted(hits[i][POS] for i in members) for members in clumps]
    assert [0] in groups
    assert [100, 1_000_100] in groups

    row = {
        "trait_uri": TRUTH_URI,
        "chrom": "chr22",
        "lead_pos": 10,
        "lead_rsid": "rs1",
        "lead_p": 1e-8,
        "n_associations": 1,
        "start": 10,
        "end": 10,
        "mapped_genes": {"raw": "FBLN1 - LINC01589", "tokens": ["FBLN1 - LINC01589"]},
    }
    table = pa.Table.from_pylist([row], schema=LOCI_SCHEMA)
    got = table.column("mapped_genes").to_pylist()[0]
    assert got["raw"] == "FBLN1 - LINC01589"
    assert got["tokens"] == ["FBLN1 - LINC01589"]
    print("load_gwas.py checks passed", flush=True)


def _flush(buffers: list[list], writer: pq.ParquetWriter) -> int:
    n = len(buffers[0])
    if n == 0:
        return 0
    arrays = [
        pa.array(buffers[i], type=ASSOC_SCHEMA.field(i).type)
        for i in range(len(ASSOC_SCHEMA))
    ]
    writer.write_table(pa.Table.from_arrays(arrays, schema=ASSOC_SCHEMA))
    for column in buffers:
        column.clear()
    return n


def _choose_name(names: Counter[str]) -> str | None:
    if not names:
        return None
    return sorted(names.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _locus_row(
    trait_uri: str | None,
    chrom: str,
    hits: list[tuple],
    members: list[int],
    token_cache: dict[str, tuple[str, ...]],
) -> dict:
    lead_i = _lead_index(hits, members)
    lead = hits[lead_i]
    positions = [hits[i][POS] for i in members]
    start = min(positions)
    end = max(positions)
    if end - start > 2 * WINDOW_BP:
        raise RuntimeError("clump span is wider than 2 Mb")
    seen: set[str] = set()
    tokens: list[str] = []
    ordered = [lead_i] + [i for i in members if i != lead_i]
    for i in ordered:
        gene = hits[i][GENE]
        if gene is None:
            continue
        cached = token_cache.get(gene)
        if cached is None:
            cached = tuple(gene_tokens(gene))
            token_cache[gene] = cached
        for token in cached:
            if token not in seen:
                seen.add(token)
                tokens.append(token)
    return {
        "trait_uri": trait_uri,
        "chrom": chrom,
        "lead_pos": lead[POS],
        "lead_rsid": lead[RSID],
        "lead_p": lead[PVAL],
        "n_associations": len(members),
        "start": start,
        "end": end,
        "mapped_genes": {"raw": lead[GENE], "tokens": tokens},
    }


def _write_pylist(path: Path, schema: pa.Schema, rows: list[dict]) -> None:
    tmp = path.with_name(path.name + ".partial")
    if tmp.exists():
        tmp.unlink()
    writer = pq.ParquetWriter(tmp, schema, compression="zstd")
    try:
        for start in range(0, len(rows), BATCH_ROWS):
            batch = rows[start : start + BATCH_ROWS]
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
        writer.close()
    except Exception:
        writer.close()
        if tmp.exists():
            tmp.unlink()
        raise
    tmp.replace(path)


def _fmt_int(value: int) -> str:
    return f"{value:,}"


def _fmt_frac(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return f"{numerator:,} / 0"
    return f"{numerator:,} / {denominator:,} ({numerator / denominator:.6f})"


def _render_qc(stats: dict) -> str:
    truth_intergenic = "PASS" if stats["truth_intergenic_pass"] else "FAIL"
    truth_empty = "PASS" if stats["truth_empty_pass"] else "FAIL"
    row_note = (
        "equal to the 1,205,644 baseline"
        if stats["n_rows"] == EXPECTED_ROWS
        else f"differs from the 1,205,644 baseline by {stats['n_rows'] - EXPECTED_ROWS:+,}"
    )
    multi_name = stats["n_uris_multi_name"]
    unnamed = stats["n_uris_unnamed"]
    return f"""# GWAS Catalog

This table is a literature index of reported associations. It is not a file of genome-wide summary statistics, and it does not contain individual genotypes. Do not compute a polygenic score, or any other score of a person, from these rows. Lead SNPs are the associations an author reported, not an effect-size vector that can be summed.

## Source

- Zip: `raw/gwas/gwas-catalog-associations_ontology-annotated-full.zip`
- Member: `gwas-catalog-download-associations-alt-full.tsv`
- Compressed bytes: {_fmt_int(stats["zip_compressed"])}
- Uncompressed bytes: {_fmt_int(stats["zip_uncompressed"])}
- Catalog columns kept: all 38, in file order, as text. Empty cells stay empty strings.
- `DATE ADDED TO CATALOG` ranges from `{stats["date_min"]}` to `{stats["date_max"]}` (the values are `YYYY-MM-DD`).
- There is no assembly column. Current `CHR_POS` values are GRCh38. An empty or multi-valued `CHR_POS` is an unplaced association (haplotype, interaction, or similar), not a second assembly. Those cells are not split and are not lifted over.

## Columns added

| Column | Rule |
|---|---|
| `chrom` | Canonical chromosome from `CHR_ID` via `canon_chrom`. Null when `CHR_ID` is empty or is not one chromosome (`6 x 6`, `2;2`, and the same kind of cell). |
| `chrom_raw` | `CHR_ID` with surrounding whitespace removed. Null when that is empty. Multi-chromosome text is kept here and does not become `chrom`. |
| `pos` | `CHR_POS` as an integer >= 1. Null when the cell is empty, zero, or not a single integer. Multi-position cells stay null. |
| `rsid` | `norm_rsid(SNP_ID_CURRENT)`. The file stores a bare integer (`4823308` becomes `rs4823308`). Empty ids and suffixed junk such as `7658266b` stay null. The `SNPS` column is not a fallback. |
| `p_value_float` | `P-VALUE` when the whole cell is a plain finite number >= 0. The raw text stays in `P-VALUE`. |
| `effect_float` | `OR or BETA` when the whole cell is a plain finite number. `Infinity`, blanks, and prose stay null. The raw text stays in `OR or BETA`. |
| `risk_allele` | Allele token from `STRONGEST SNP-RISK ALLELE`. `rs128494-C` stores `C`. `rs145611121-GC` stores `GC`. This is not a REF/ALT pair. `?` is null. Several SNP tokens stay in order, joined by `;`. |
| `variant_key` | Null on every row. This file has no ref or alt. Alleles are not guessed from `risk_allele`. A later rsID match to ClinVar or gnomAD is how a row can gain a `chrom:pos:ref:alt` key. |
| `locus_key` | `chrom:pos` when both `chrom` and `pos` exist. Otherwise null. |

Sample size and ancestry stay as prose in `INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE`. They are not parsed.

P-values that underflow float64: {_fmt_int(stats["n_p_underflow"])}. Those ranks are stored as the smallest positive float so they remain more significant than any representable p-value. Negative p-values left null: {_fmt_int(stats["n_p_negative"])}. Non-representable effects left null: {_fmt_int(stats["n_effect_underflow"])}.

## Counts

| Measure | Value |
|---|---|
| Association rows | {_fmt_int(stats["n_rows"])} ({row_note}) |
| Null `pos` | {_fmt_frac(stats["n_null_pos"], stats["n_rows"])} |
| Null `chrom` | {_fmt_frac(stats["n_null_chrom"], stats["n_rows"])} |
| `chrom` null but `chrom_raw` set (multi-chromosome or other non-canonical `CHR_ID`) | {_fmt_int(stats["n_chrom_raw_unmapped"])} |
| `CONTEXT` contains `intergenic` | {_fmt_frac(stats["n_context_intergenic"], stats["n_rows"])} |
| `INTERGENIC` flag equal to `1` | {_fmt_int(stats["n_intergenic_flag"])} |
| Distinct non-empty `MAPPED_TRAIT_URI` | {_fmt_int(stats["n_traits"])} |
| Trait-table rows, including one empty-URI bucket | {_fmt_int(stats["n_trait_rows"])} |
| Loci (clumps) | {_fmt_int(stats["n_loci"])} |
| Raw rows per locus | {stats["rows_per_locus"]:.2f} |
| Eligible rows | {_fmt_int(stats["n_eligible"])} |
| Eligible rows per locus | {stats["eligible_per_locus"]:.2f} |
| Sum of locus `n_associations` | {_fmt_int(stats["sum_memberships"])} |
| Singleton loci | {_fmt_frac(stats["n_singletons"], stats["n_loci"])} |
| `CNV` equal to `Y` (kept in associations, excluded from clumps) | {_fmt_int(stats["n_cnv"])} |
| Rows with a `locus_key` | {_fmt_int(stats["n_locus_key"])} |
| Rows with an `rsid` | {_fmt_int(stats["n_rsid"])} |
| Rows with a parsed `effect_float` | {_fmt_int(stats["n_effect"])} |
| Rows with a `risk_allele` | {_fmt_int(stats["n_risk"])} |
| Rows whose `risk_allele` has more than one token | {_fmt_int(stats["n_risk_multi"])} |
| Rows with more than one mapped trait URI | {_fmt_int(stats["n_multi_uri"])} |
| Rows with an empty `MAPPED_TRAIT_URI` | {_fmt_int(stats["n_empty_uri"])} |
| Multi-URI rows whose trait-name list could not be aligned | {_fmt_int(stats["n_unaligned_names"])} |
| URIs with more than one mapped name | {_fmt_int(multi_name)} |
| Non-empty URIs with no mapped name | {_fmt_int(unnamed)} |

`CONTEXT` and the `INTERGENIC` flag are not the same claim. The intergenic fraction above is the `CONTEXT` substring. Some rows have `INTERGENIC` of `1` and a `CONTEXT` that does not contain `intergenic`.

The sum of locus `n_associations` is larger than the eligible row count when a row lists more than one `MAPPED_TRAIT_URI`: that row is clumped once inside each URI. It is smaller than the raw row count because null positions never enter a locus. File check of that sum: {_fmt_int(stats["sum_memberships_file"])}.

Parquet row counts read back from the files: associations {_fmt_int(stats["file_rows"])}, traits {_fmt_int(stats["file_trait_rows"])}, loci {_fmt_int(stats["file_loci"])}. Null `pos` read back from associations: {_fmt_int(stats["file_null_pos"])}. `CONTEXT` containing `intergenic`, read back: {_fmt_int(stats["file_context_intergenic"])}. `variant_key` nulls read back: {_fmt_int(stats["file_variant_key_nulls"])}.

## Clumping

Within one `MAPPED_TRAIT_URI`, and within one chromosome, associations with a canonical `chrom`, a `pos` >= 1, and `CNV` not equal to `Y` are clumped as follows.

1. The lead is the still-unassigned association with the smallest `p_value_float`. Ties break by smaller position, then `rsid`, then file line.
2. Every other still-unassigned association on that chromosome within 1,000,000 bp of the lead joins the clump.
3. Repeat until none remain.

A missing p-value sorts after every real p-value. This is not single-linkage: a SNP within 1 Mb of a member but more than 1 Mb from the lead stays out and can seed its own clump. Leads on the same trait and chromosome are therefore more than 1,000,000 bp apart, and every member is within 1,000,000 bp of its own lead.

`start` and `end` are the minimum and maximum member positions (1-based, inclusive). They are the span of the collapsed associations, not the ±1 Mb search window. A one-SNP locus has `start` = `end` = `lead_pos`.

Rows with a null position stay in `associations.parquet` and stay out of `loci.parquet`. `CNV` of `Y` would do the same; this release has {_fmt_int(stats["n_cnv"])} such rows. Associations with an empty `MAPPED_TRAIT_URI` are clumped together under a null `trait_uri` when they have a position, and they are the empty-URI row in `traits.parquet`.

## Traits and mapped genes

`traits.parquet` has one row per `MAPPED_TRAIT_URI`:

- `trait_uri`
- `mapped_trait`, the ontology label. A URI that is alone on a row keeps the whole `MAPPED_TRAIT` text, commas included (`osteoarthritis, hand`). Several URIs on one row are paired with comma-separated labels only when the counts match. The stored name is the most common label, then the alphabetical label on a tie.
- `n_associations`, catalog rows that list that URI. A row with two URIs counts in both.
- `n_author_traits`, distinct non-empty `DISEASE/TRAIT` strings on those rows.

`loci.parquet` column `mapped_genes` is a struct, not a plain string:

- `mapped_genes.raw` is the lead association's `MAPPED_GENE` text, unchanged. Null when the lead gene is blank.
- `mapped_genes.tokens` is the comma-split token list from every member of the clump, lead first, then the other members in position order, with duplicates removed. Tokens are trimmed. Empty tokens and a token that is exactly `-` are removed. `HLA-DQB1` stays one token.

The intergenic catalog form `FBLN1 - LINC01589` has no comma, so it stays one token. A symbol matches a token, not a substring: `FBLN1` does not match that token. The raw string is kept so a later step can see the unsplit text. `BRCA1` does match a token of `BRCA1, NBR2`.

## Truth set

`warehouse/qc/truth_set.tsv`.

### gwas_intergenic ({truth_intergenic})

File line 8 (header is line 1). Read back from `associations.parquet` in file order: chrom `{stats["line8_chrom"]}`, pos `{stats["line8_pos"]}`, rsid `{stats["line8_rsid"]}`, locus_key `{stats["line8_locus_key"]}`, variant_key `{stats["line8_variant_key"]}`, CONTEXT `{stats["line8_context"]}`, trait URI `{stats["line8_uri"]}`, CNV `{stats["line8_cnv"]}`, risk_allele `{stats["line8_risk"]}`.

The position is not null, so the row is outside the null-position set. `CNV` is not `Y`, and both chrom and pos exist, so the row is eligible for a clump. The clump that contains it is trait `{stats["truth_locus_trait"]}`, chrom `{stats["truth_locus_chrom"]}`, lead `{stats["truth_locus_rsid"]}` at `{stats["truth_locus_lead_pos"]}` (p `{stats["truth_locus_lead_p"]}`), span `{stats["truth_locus_start"]}`–`{stats["truth_locus_end"]}`, members `{stats["truth_locus_n"]}`. That locus row is present in `loci.parquet`.

### gwas_empty_pos ({truth_empty})

File line 5. Read back from `associations.parquet`: study `{stats["line5_study"]}`, pubmed `{stats["line5_pubmed"]}`, SNPS `{stats["line5_snps"]}`, CHR_POS `{stats["line5_chr_pos"]!r}`, chrom `{stats["line5_chrom"]}`, pos `{stats["line5_pos"]}`, locus_key `{stats["line5_locus_key"]}`, variant_key `{stats["line5_variant_key"]}`, risk_allele `{stats["line5_risk"]}`, trait URI `{stats["line5_uri"]}`.

The row is in the association table. `pos` and `locus_key` are null. `STRONGEST SNP-RISK ALLELE` on this row looks like a coordinate plus an allele (`chr22:37462936-G`); only the allele token is stored, and it is not written into `pos`. The clumper never sees this row. Every locus has a non-null `lead_pos` and `chrom` (null lead positions in the locus file: {_fmt_int(stats["loci_null_lead_pos"])}), so this association is absent from `loci.parquet`.
"""


def _capture_line(rec: list[str], parsed: dict) -> dict:
    out = {
        "CHR_ID": rec[COLUMNS.index("CHR_ID")],
        "CHR_POS": rec[COLUMNS.index("CHR_POS")],
        "CONTEXT": rec[COLUMNS.index("CONTEXT")],
        "CNV": rec[COLUMNS.index("CNV")],
        "MAPPED_GENE": rec[COLUMNS.index("MAPPED_GENE")],
        "MAPPED_TRAIT_URI": rec[COLUMNS.index("MAPPED_TRAIT_URI")],
        "STUDY ACCESSION": rec[COLUMNS.index("STUDY ACCESSION")],
        "PUBMEDID": rec[COLUMNS.index("PUBMEDID")],
        "SNPS": rec[COLUMNS.index("SNPS")],
        "STRONGEST SNP-RISK ALLELE": rec[COLUMNS.index("STRONGEST SNP-RISK ALLELE")],
    }
    out.update(parsed)
    return out


def load() -> dict:
    if not ZIP_PATH.exists():
        raise SystemExit(f"missing {ZIP_PATH}")
    GWAS_DIR.mkdir(parents=True, exist_ok=True)
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)

    stats: dict = {
        "n_rows": 0,
        "n_null_pos": 0,
        "n_null_chrom": 0,
        "n_chrom_raw_unmapped": 0,
        "n_context_intergenic": 0,
        "n_intergenic_flag": 0,
        "n_p_underflow": 0,
        "n_p_negative": 0,
        "n_p_null": 0,
        "n_effect": 0,
        "n_effect_underflow": 0,
        "n_risk": 0,
        "n_risk_multi": 0,
        "n_rsid": 0,
        "n_locus_key": 0,
        "n_cnv": 0,
        "n_eligible": 0,
        "n_multi_uri": 0,
        "n_empty_uri": 0,
        "n_unaligned_names": 0,
        "expected_memberships": 0,
        "date_min": None,
        "date_max": None,
    }
    traits: dict[str | None, TraitAcc] = defaultdict(TraitAcc)
    buckets: dict[tuple[str | None, str], list[tuple]] = defaultdict(list)
    line5 = None
    line8 = None
    truth_member = None

    assoc_tmp = GWAS_DIR / "associations.parquet.partial"
    if assoc_tmp.exists():
        assoc_tmp.unlink()
    writer = pq.ParquetWriter(assoc_tmp, ASSOC_SCHEMA, compression="zstd")
    buffers: list[list] = [[] for _ in ASSOC_SCHEMA.names]
    n_schema = len(ASSOC_SCHEMA)

    try:
        csv.field_size_limit(2**31 - 1)
        with zipfile.ZipFile(ZIP_PATH) as zf:
            info = zf.getinfo(MEMBER)
            stats["zip_compressed"] = info.compress_size
            stats["zip_uncompressed"] = info.file_size
            with zf.open(info, "r") as binary:
                text = io.TextIOWrapper(binary, encoding="utf-8", newline="")
                reader = csv.reader(text, delimiter="\t")
                header = next(reader)
                if header and header[0].startswith("\ufeff"):
                    header[0] = header[0].lstrip("\ufeff")
                if header != COLUMNS:
                    raise SystemExit("GWAS header does not match the 38 catalog columns")
                for line_no, rec in enumerate(reader, start=2):
                    if len(rec) != 38:
                        raise SystemExit(f"line {line_no} has {len(rec)} columns, expected 38")
                    stats["n_rows"] += 1
                    chrom_raw_cell = rec[11].strip()
                    chrom_raw = chrom_raw_cell or None
                    chrom = canon_chrom(chrom_raw_cell)
                    if chrom is not None:
                        chrom = sys.intern(chrom)
                    pos = parse_pos(rec[12])
                    rsid = norm_rsid(rec[23])
                    p_value, p_flag = parse_p_value(rec[27])
                    effect, effect_flag = parse_effect(rec[30])
                    risk = parse_risk_allele(rec[20])
                    if risk is not None:
                        risk = sys.intern(risk)
                    cnv_yes = is_cnv(rec[33])
                    if chrom_raw is not None:
                        chrom_raw = sys.intern(chrom_raw)
                    if rsid is not None:
                        rsid = sys.intern(rsid)
                    gene_cell = rec[14]
                    gene = None if not gene_cell.strip() else sys.intern(gene_cell)
                    locus_key = None
                    if chrom is not None and pos is not None:
                        locus_key = f"{chrom}:{pos}"

                    if pos is None:
                        stats["n_null_pos"] += 1
                    if chrom is None:
                        stats["n_null_chrom"] += 1
                        if chrom_raw is not None:
                            stats["n_chrom_raw_unmapped"] += 1
                    context = rec[24]
                    if "intergenic" in context.lower():
                        stats["n_context_intergenic"] += 1
                    if rec[25].strip() == "1":
                        stats["n_intergenic_flag"] += 1
                    if p_value is None:
                        stats["n_p_null"] += 1
                    if p_flag == "underflow":
                        stats["n_p_underflow"] += 1
                    elif p_flag == "negative":
                        stats["n_p_negative"] += 1
                    if effect is not None:
                        stats["n_effect"] += 1
                    if effect_flag == "underflow":
                        stats["n_effect_underflow"] += 1
                    if risk is not None:
                        stats["n_risk"] += 1
                        if ";" in risk:
                            stats["n_risk_multi"] += 1
                    if rsid is not None:
                        stats["n_rsid"] += 1
                    if locus_key is not None:
                        stats["n_locus_key"] += 1
                    if cnv_yes:
                        stats["n_cnv"] += 1

                    added_on = rec[0].strip()
                    if added_on:
                        if stats["date_min"] is None or added_on < stats["date_min"]:
                            stats["date_min"] = added_on
                        if stats["date_max"] is None or added_on > stats["date_max"]:
                            stats["date_max"] = added_on

                    uris = split_uris(rec[35])
                    if len(uris) > 1:
                        stats["n_multi_uri"] += 1
                        if len([part.strip() for part in rec[34].split(",")]) != len(uris):
                            stats["n_unaligned_names"] += 1
                    if not uris:
                        stats["n_empty_uri"] += 1
                    author = rec[7].strip()
                    if author:
                        author = sys.intern(author)
                    for uri, name in trait_name_pairs(rec[34], uris):
                        if uri is not None:
                            uri = sys.intern(uri)
                        acc = traits[uri]
                        acc.n += 1
                        if name:
                            acc.names[name] += 1
                        if author:
                            acc.authors.add(author)

                    if chrom is not None and pos is not None and not cnv_yes:
                        hit = (pos, p_value, rsid, gene, line_no)
                        targets = uris if uris else [None]
                        for uri in targets:
                            buckets[(uri, chrom)].append(hit)
                        stats["n_eligible"] += 1
                        stats["expected_memberships"] += len(targets)

                    parsed = {
                        "chrom": chrom,
                        "pos": pos,
                        "rsid": rsid,
                        "locus_key": locus_key,
                        "variant_key": None,
                        "risk_allele": risk,
                        "p_value_float": p_value,
                    }
                    if line_no == 5:
                        line5 = _capture_line(rec, parsed)
                    elif line_no == 8:
                        line8 = _capture_line(rec, parsed)

                    row = rec + [
                        chrom,
                        chrom_raw,
                        pos,
                        rsid,
                        p_value,
                        effect,
                        risk,
                        None,
                        locus_key,
                    ]
                    if len(row) != n_schema:
                        raise RuntimeError("association row width does not match the schema")
                    for i, value in enumerate(row):
                        buffers[i].append(value)
                    if len(buffers[0]) >= BATCH_ROWS:
                        _flush(buffers, writer)
                        print(f"associations {stats['n_rows']:,}", flush=True)
                _flush(buffers, writer)
        writer.close()
    except Exception:
        writer.close()
        if assoc_tmp.exists():
            assoc_tmp.unlink()
        raise

    if line5 is None or line8 is None:
        assoc_tmp.unlink()
        raise SystemExit("truth-set lines 5 and 8 were not in the TSV")

    print(
        f"scanned {stats['n_rows']:,} rows; clumping {len(buckets):,} trait-chromosome buckets",
        flush=True,
    )
    gc.collect()

    token_cache: dict[str, tuple[str, ...]] = {}
    loci_rows: list[dict] = []
    sum_memberships = 0
    n_singletons = 0
    keys = sorted(
        buckets,
        key=lambda key: (key[0] is None, key[0] or "", chrom_rank(key[1]), key[1]),
    )
    for key in keys:
        uri, chrom = key
        hits = buckets.pop(key)
        for members in clump_hits(hits):
            locus = _locus_row(uri, chrom, hits, members, token_cache)
            sum_memberships += locus["n_associations"]
            if locus["n_associations"] == 1:
                n_singletons += 1
            if (
                uri == TRUTH_URI
                and chrom == TRUTH_CHROM
                and any(
                    hits[i][POS] == TRUTH_POS and hits[i][RSID] == TRUTH_RSID and hits[i][LINE] == 8
                    for i in members
                )
            ):
                truth_member = locus
            loci_rows.append(locus)
    if sum_memberships != stats["expected_memberships"]:
        assoc_tmp.unlink()
        raise RuntimeError(
            f"clump memberships {sum_memberships} != eligible pairs {stats['expected_memberships']}"
        )
    del buckets
    gc.collect()

    trait_rows = []
    n_uris_multi_name = 0
    n_uris_unnamed = 0
    n_traits = 0
    for uri, acc in traits.items():
        if uri is not None:
            n_traits += 1
            if not acc.names:
                n_uris_unnamed += 1
        if len(acc.names) > 1:
            n_uris_multi_name += 1
        trait_rows.append(
            {
                "trait_uri": uri,
                "mapped_trait": _choose_name(acc.names),
                "n_associations": acc.n,
                "n_author_traits": len(acc.authors),
            }
        )
    trait_rows.sort(key=lambda row: (row["trait_uri"] is None, row["trait_uri"] or ""))
    del traits

    assoc_path = GWAS_DIR / "associations.parquet"
    trait_path = GWAS_DIR / "traits.parquet"
    loci_path = GWAS_DIR / "loci.parquet"
    print(f"writing {len(trait_rows):,} traits and {len(loci_rows):,} loci", flush=True)
    _write_pylist(trait_path, TRAIT_SCHEMA, trait_rows)
    _write_pylist(loci_path, LOCI_SCHEMA, loci_rows)
    assoc_tmp.replace(assoc_path)

    file_stats = _verify_files(assoc_path, trait_path, loci_path, line5, line8, truth_member)
    stats.update(file_stats)
    stats["n_traits"] = n_traits
    stats["n_trait_rows"] = len(trait_rows)
    stats["n_loci"] = len(loci_rows)
    stats["sum_memberships"] = sum_memberships
    stats["n_singletons"] = n_singletons
    stats["n_uris_multi_name"] = n_uris_multi_name
    stats["n_uris_unnamed"] = n_uris_unnamed
    stats["rows_per_locus"] = stats["n_rows"] / stats["n_loci"] if stats["n_loci"] else 0.0
    stats["eligible_per_locus"] = (
        stats["n_eligible"] / stats["n_loci"] if stats["n_loci"] else 0.0
    )
    stats["line5_study"] = line5["STUDY ACCESSION"]
    stats["line5_pubmed"] = line5["PUBMEDID"]
    stats["line5_snps"] = line5["SNPS"]
    stats["line5_chr_pos"] = line5["CHR_POS"]
    stats["line5_chrom"] = line5["chrom"]
    stats["line5_pos"] = line5["pos"]
    stats["line5_locus_key"] = line5["locus_key"]
    stats["line5_variant_key"] = line5["variant_key"]
    stats["line5_risk"] = line5["risk_allele"]
    stats["line5_uri"] = line5["MAPPED_TRAIT_URI"]
    stats["line8_chrom"] = line8["chrom"]
    stats["line8_pos"] = line8["pos"]
    stats["line8_rsid"] = line8["rsid"]
    stats["line8_locus_key"] = line8["locus_key"]
    stats["line8_variant_key"] = line8["variant_key"]
    stats["line8_context"] = line8["CONTEXT"]
    stats["line8_uri"] = line8["MAPPED_TRAIT_URI"]
    stats["line8_cnv"] = line8["CNV"]
    stats["line8_risk"] = line8["risk_allele"]
    if truth_member is None:
        stats["truth_locus_trait"] = None
        stats["truth_locus_chrom"] = None
        stats["truth_locus_rsid"] = None
        stats["truth_locus_lead_pos"] = None
        stats["truth_locus_lead_p"] = None
        stats["truth_locus_start"] = None
        stats["truth_locus_end"] = None
        stats["truth_locus_n"] = None
    else:
        stats["truth_locus_trait"] = truth_member["trait_uri"]
        stats["truth_locus_chrom"] = truth_member["chrom"]
        stats["truth_locus_rsid"] = truth_member["lead_rsid"]
        stats["truth_locus_lead_pos"] = truth_member["lead_pos"]
        stats["truth_locus_lead_p"] = truth_member["lead_p"]
        stats["truth_locus_start"] = truth_member["start"]
        stats["truth_locus_end"] = truth_member["end"]
        stats["truth_locus_n"] = truth_member["n_associations"]

    stats["truth_intergenic_pass"] = _truth_intergenic_ok(stats)
    stats["truth_empty_pass"] = _truth_empty_ok(stats)
    QC_PATH.write_text(_render_qc(stats), encoding="utf-8", newline="\n")
    return stats


def _truth_intergenic_ok(stats: dict) -> bool:
    return (
        stats["line8_chrom"] == TRUTH_CHROM
        and stats["line8_pos"] == TRUTH_POS
        and stats["line8_rsid"] == TRUTH_RSID
        and stats["line8_locus_key"] == f"{TRUTH_CHROM}:{TRUTH_POS}"
        and stats["line8_variant_key"] is None
        and stats["line8_context"] == "intergenic_variant"
        and stats["line8_uri"] == TRUTH_URI
        and stats["line8_cnv"] != "Y"
        and stats["line8_readback_ok"] is True
        and stats["truth_locus_in_file"] is True
        and stats["truth_locus_trait"] == TRUTH_URI
        and stats["truth_locus_chrom"] == TRUTH_CHROM
        and stats["truth_locus_start"] is not None
        and stats["truth_locus_start"] <= TRUTH_POS <= stats["truth_locus_end"]
        and abs(stats["truth_locus_lead_pos"] - TRUTH_POS) <= WINDOW_BP
    )


def _truth_empty_ok(stats: dict) -> bool:
    return (
        stats["line5_study"] == EMPTY_STUDY
        and stats["line5_pubmed"] == EMPTY_PUBMED
        and stats["line5_snps"] == EMPTY_SNPS
        and stats["line5_chr_pos"] == ""
        and stats["line5_uri"] == TRUTH_URI
        and stats["line5_chrom"] is None
        and stats["line5_pos"] is None
        and stats["line5_locus_key"] is None
        and stats["line5_variant_key"] is None
        and stats["line5_risk"] == "G"
        and stats["line5_readback_ok"] is True
        and stats["loci_null_lead_pos"] == 0
        and stats["loci_null_chrom"] == 0
    )


def _verify_files(
    assoc_path: Path,
    trait_path: Path,
    loci_path: Path,
    line5: dict,
    line8: dict,
    truth_member: dict | None,
) -> dict:
    assoc_meta = pq.read_metadata(assoc_path)
    trait_meta = pq.read_metadata(trait_path)
    loci_meta = pq.read_metadata(loci_path)
    pos = pq.read_table(assoc_path, columns=["pos"]).column("pos")
    variant_key = pq.read_table(assoc_path, columns=["variant_key"]).column("variant_key")
    context = pq.read_table(assoc_path, columns=["CONTEXT"]).column("CONTEXT")
    context_hits = pc.sum(
        pc.match_substring(pc.utf8_lower(context), "intergenic", ignore_case=False)
    ).as_py()
    early = next(
        pq.ParquetFile(assoc_path).iter_batches(
            batch_size=10,
            columns=[
                "STUDY ACCESSION",
                "PUBMEDID",
                "SNPS",
                "CHR_POS",
                "CHR_ID",
                "CONTEXT",
                "CNV",
                "MAPPED_TRAIT_URI",
                "chrom",
                "pos",
                "rsid",
                "locus_key",
                "variant_key",
                "risk_allele",
            ],
        )
    )
    early_rows = early.to_pylist()
    # File line 5 is data index 3; file line 8 is data index 6.
    got5 = early_rows[3]
    got8 = early_rows[6]
    line5_ok = (
        got5["STUDY ACCESSION"] == line5["STUDY ACCESSION"]
        and got5["PUBMEDID"] == line5["PUBMEDID"]
        and got5["SNPS"] == line5["SNPS"]
        and got5["CHR_POS"] == line5["CHR_POS"]
        and got5["pos"] == line5["pos"]
        and got5["chrom"] == line5["chrom"]
        and got5["locus_key"] == line5["locus_key"]
        and got5["variant_key"] is None
        and got5["risk_allele"] == line5["risk_allele"]
        and got5["MAPPED_TRAIT_URI"] == line5["MAPPED_TRAIT_URI"]
    )
    line8_ok = (
        got8["chrom"] == line8["chrom"]
        and got8["pos"] == line8["pos"]
        and got8["rsid"] == line8["rsid"]
        and got8["locus_key"] == line8["locus_key"]
        and got8["variant_key"] is None
        and got8["CONTEXT"] == line8["CONTEXT"]
        and got8["CNV"] == line8["CNV"]
        and got8["MAPPED_TRAIT_URI"] == line8["MAPPED_TRAIT_URI"]
        and got8["risk_allele"] == line8["risk_allele"]
    )

    loci = pq.read_table(
        loci_path,
        columns=["trait_uri", "chrom", "lead_pos", "lead_rsid", "lead_p", "n_associations", "start", "end"],
    )
    sum_memberships_file = pc.sum(loci.column("n_associations")).as_py()
    loci_null_lead = pc.sum(pc.is_null(loci.column("lead_pos")).cast(pa.int64())).as_py()
    loci_null_chrom = pc.sum(pc.is_null(loci.column("chrom")).cast(pa.int64())).as_py()
    truth_locus_in_file = False
    if truth_member is not None:
        uri = loci.column("trait_uri")
        chrom = loci.column("chrom")
        lead_pos = loci.column("lead_pos")
        mask = pc.and_(
            pc.and_(pc.equal(uri, truth_member["trait_uri"]), pc.equal(chrom, truth_member["chrom"])),
            pc.equal(lead_pos, truth_member["lead_pos"]),
        )
        picked = loci.filter(mask)
        for row in picked.to_pylist():
            if (
                row["lead_rsid"] == truth_member["lead_rsid"]
                and row["start"] == truth_member["start"]
                and row["end"] == truth_member["end"]
                and row["n_associations"] == truth_member["n_associations"]
            ):
                truth_locus_in_file = True
                break

    return {
        "file_rows": assoc_meta.num_rows,
        "file_trait_rows": trait_meta.num_rows,
        "file_loci": loci_meta.num_rows,
        "file_null_pos": pos.null_count,
        "file_variant_key_nulls": variant_key.null_count,
        "file_context_intergenic": int(context_hits),
        "sum_memberships_file": int(sum_memberships_file),
        "loci_null_lead_pos": int(loci_null_lead),
        "loci_null_chrom": int(loci_null_chrom),
        "line5_readback_ok": line5_ok,
        "line8_readback_ok": line8_ok,
        "truth_locus_in_file": truth_locus_in_file,
    }


def main() -> None:
    _check()
    stats = load()
    ok = stats["truth_intergenic_pass"] and stats["truth_empty_pass"]
    ok = ok and stats["n_rows"] == stats["file_rows"] == EXPECTED_ROWS
    ok = ok and stats["file_null_pos"] == stats["n_null_pos"]
    ok = ok and stats["file_context_intergenic"] == stats["n_context_intergenic"]
    ok = ok and stats["file_variant_key_nulls"] == stats["n_rows"]
    ok = ok and stats["sum_memberships_file"] == stats["sum_memberships"]
    ok = ok and stats["file_loci"] == stats["n_loci"]
    ok = ok and stats["file_trait_rows"] == stats["n_trait_rows"]
    print(
        f"rows {stats['n_rows']:,} null_pos {stats['n_null_pos']:,} "
        f"intergenic {stats['n_context_intergenic']:,} traits {stats['n_traits']:,} "
        f"loci {stats['n_loci']:,}",
        flush=True,
    )
    print(f"truth intergenic {'PASS' if stats['truth_intergenic_pass'] else 'FAIL'}", flush=True)
    print(f"truth empty_pos {'PASS' if stats['truth_empty_pass'] else 'FAIL'}", flush=True)
    print(f"wrote {QC_PATH}", flush=True)
    if not ok:
        raise SystemExit("GWAS QC checks failed")


if __name__ == "__main__":
    main()
