# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load GENCODE v49 gene models into interval Parquet tables.

Coordinates stay 1-based and inclusive, as they are in the GTF.
A VCF position P falls in a feature when start <= P <= end.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from collections import Counter
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from coords import bare_ensg, bare_enst, canon_chrom

GTF = ROOT / "raw" / "gencode" / "gencode.v49.annotation.gtf.gz"
OUT = ROOT / "warehouse" / "intervals"
QC_PATH = ROOT / "warehouse" / "qc" / "intervals.md"
TRUTH = ROOT / "warehouse" / "qc" / "truth_set.tsv"

KEEP = frozenset({"gene", "transcript", "exon", "CDS", "UTR", "stop_codon"})
WANTED = frozenset(
    {
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_id",
        "transcript_type",
        "tag",
        "hgnc_id",
        "exon_number",
    }
)
GENE_BASELINE = 78691
PROTEIN_BASELINE = 20097
CANON = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

_ATTR = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
# GENCODE v49 leaves exon_number unquoted (`exon_number 1;`). Quoted form is accepted too.
_EXON_NUM = re.compile(r'\bexon_number\s+"?(\d+)"?')
_SPAN = re.compile(r"(ENSG\d+\.\d+)\s+(chr[0-9XYM]+):(\d+)-(\d+)\s+(\S+)")

GENE_SCHEMA = pa.schema(
    [
        pa.field("ensg", pa.string(), nullable=False),
        pa.field("ensg_version", pa.string(), nullable=False),
        pa.field("symbol", pa.string()),
        pa.field("chrom", pa.string(), nullable=False),
        pa.field("start", pa.int64(), nullable=False),
        pa.field("end", pa.int64(), nullable=False),
        pa.field("strand", pa.string(), nullable=False),
        pa.field("biotype", pa.string()),
        pa.field("mane_select", pa.bool_(), nullable=False),
        pa.field("hgnc_id", pa.string()),
    ]
)
TRANSCRIPT_SCHEMA = pa.schema(
    list(GENE_SCHEMA)
    + [
        pa.field("enst", pa.string(), nullable=False),
        pa.field("enst_version", pa.string(), nullable=False),
        pa.field("transcript_type", pa.string()),
        pa.field("tags", pa.list_(pa.string()), nullable=False),
    ]
)
INTERVAL_FIELDS = [
    pa.field("ensg", pa.string(), nullable=False),
    pa.field("enst", pa.string(), nullable=False),
    pa.field("chrom", pa.string(), nullable=False),
    pa.field("start", pa.int64(), nullable=False),
    pa.field("end", pa.int64(), nullable=False),
    pa.field("strand", pa.string(), nullable=False),
    pa.field("exon_number", pa.int64()),
]
EXON_SCHEMA = pa.schema(INTERVAL_FIELDS)
CDS_SCHEMA = pa.schema(INTERVAL_FIELDS + [pa.field("feature", pa.string(), nullable=False)])
UTR_SCHEMA = EXON_SCHEMA

PATHS = {
    "genes": OUT / "genes.parquet",
    "transcripts": OUT / "transcripts.parquet",
    "exons": OUT / "exons.parquet",
    "cds": OUT / "cds.parquet",
    "utrs": OUT / "utrs.parquet",
}


def log(message: str) -> None:
    print(message, flush=True)


def parse_attrs(field: str) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for key, value in _ATTR.findall(field):
        if key not in WANTED:
            continue
        found.setdefault(key, []).append(value)
    return found


def first(attrs: dict[str, list[str]], key: str) -> str | None:
    values = attrs.get(key)
    if not values:
        return None
    text = values[0]
    if text == "":
        return None
    return text


def opt_intern(value: str | None) -> str | None:
    if value is None or value == "":
        return None
    return sys.intern(value)


def columns(names: list[str]) -> dict[str, list]:
    return {name: [] for name in names}


def append_row(dest: dict[str, list], row: dict) -> None:
    for key, value in row.items():
        dest[key].append(value)


def exon_number_of(field: str, bad: Counter) -> int | None:
    match = _EXON_NUM.search(field)
    if match:
        return int(match.group(1))
    if "exon_number" in field:
        bad["exon_number"] += 1
    return None


def _check_parser() -> None:
    sample = (
        'gene_id "ENSG00000186092.7"; transcript_id "ENST00000641515.2"; '
        'gene_type "protein_coding"; gene_name "OR4F5"; '
        'transcript_type "protein_coding"; transcript_name "OR4F5-201"; '
        'hgnc_id "HGNC:14825"; tag "basic"; tag "MANE_Select"; exon_number 1;'
    )
    attrs = parse_attrs(sample)
    assert attrs["gene_id"] == ["ENSG00000186092.7"]
    assert attrs["tag"] == ["basic", "MANE_Select"]
    assert "transcript_name" not in attrs
    assert bare_ensg(attrs["gene_id"][0]) == "ENSG00000186092"
    assert bare_enst(attrs["transcript_id"][0]) == "ENST00000641515"
    assert bare_ensg("ENST00000641515.2") is None
    assert exon_number_of(sample, Counter()) == 1
    assert exon_number_of('exon_number "2";', Counter()) == 2
    assert canon_chrom("chr1") == "chr1"
    assert canon_chrom("MT") == "chrM"
    assert canon_chrom("chr1_KI270706v1_random") is None


def _check_duckdb() -> None:
    con = duckdb.connect(database=":memory:")
    row = con.execute(
        """
        SELECT list_contains(['MANE_Select', 'basic'], 'MANE_Select')
           AND NOT list_contains(['MANE_Plus_Clinical'], 'MANE_Select')
           AND length(['a', 'b']) = 2
        """
    ).fetchone()
    if row is None or row[0] is not True:
        raise RuntimeError("duckdb list_contains check failed")
    con.close()


def load() -> dict:
    if not GTF.is_file():
        raise FileNotFoundError(GTF)
    OUT.mkdir(parents=True, exist_ok=True)

    genes = columns(
        [
            "ensg",
            "ensg_version",
            "symbol",
            "chrom",
            "start",
            "end",
            "strand",
            "biotype",
            "hgnc_id",
        ]
    )
    transcripts = columns(
        [
            "ensg",
            "ensg_version",
            "symbol",
            "chrom",
            "start",
            "end",
            "strand",
            "biotype",
            "mane_select",
            "hgnc_id",
            "enst",
            "enst_version",
            "transcript_type",
            "tags",
        ]
    )
    exons = columns(
        ["ensg", "enst", "chrom", "start", "end", "strand", "exon_number"]
    )
    cds = columns(
        ["ensg", "enst", "chrom", "start", "end", "strand", "exon_number", "feature"]
    )
    utrs = columns(
        ["ensg", "enst", "chrom", "start", "end", "strand", "exon_number"]
    )

    header: list[str] = []
    skipped: Counter[str] = Counter()
    dropped: Counter[str] = Counter()
    failed_chrom: Counter[str] = Counter()
    rewritten: Counter[str] = Counter()
    written: Counter[str] = Counter()
    bad = Counter()
    comments = 0
    data_rows = 0
    mane_on_gene_feature = 0
    mane_tx_versions: set[str] = set()
    gene_level_mane: set[str] = set()

    log(f"reading {GTF}")
    with gzip.open(GTF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                comments += 1
                if len(header) < 40:
                    header.append(line.rstrip("\n\r"))
                continue
            data_rows += 1
            if data_rows % 1_000_000 == 0:
                log(f"read {data_rows:,} data rows")
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 9:
                dropped["short_line"] += 1
                continue
            if len(parts) > 9:
                parts = parts[:8] + ["\t".join(parts[8:])]
            feature = parts[2]
            if feature not in KEEP:
                skipped[feature] += 1
                continue
            raw_chrom = parts[0]
            chrom = canon_chrom(raw_chrom)
            if chrom is None:
                dropped["chrom"] += 1
                failed_chrom[raw_chrom] += 1
                continue
            if chrom != raw_chrom:
                rewritten[f"{raw_chrom} -> {chrom}"] += 1
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                dropped["coord"] += 1
                continue
            attrs = parse_attrs(parts[8])
            raw_gene = first(attrs, "gene_id")
            ensg = bare_ensg(raw_gene)
            if ensg is None or raw_gene is None:
                dropped["gene_id"] += 1
                continue
            ensg = sys.intern(ensg)
            ensg_version = sys.intern(raw_gene)
            chrom = sys.intern(chrom)
            strand = sys.intern(parts[6])
            tags = [sys.intern(tag) for tag in attrs.get("tag", [])]
            has_mane = "MANE_Select" in tags
            if feature == "gene":
                if has_mane:
                    mane_on_gene_feature += 1
                    gene_level_mane.add(ensg_version)
                append_row(
                    genes,
                    {
                        "ensg": ensg,
                        "ensg_version": ensg_version,
                        "symbol": opt_intern(first(attrs, "gene_name")),
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "biotype": opt_intern(first(attrs, "gene_type")),
                        "hgnc_id": opt_intern(first(attrs, "hgnc_id")),
                    },
                )
                written["gene"] += 1
                continue
            raw_tx = first(attrs, "transcript_id")
            enst = bare_enst(raw_tx)
            if enst is None or raw_tx is None:
                dropped["transcript_id"] += 1
                continue
            enst = sys.intern(enst)
            enst_version = sys.intern(raw_tx)
            if feature == "transcript":
                if has_mane:
                    mane_tx_versions.add(ensg_version)
                append_row(
                    transcripts,
                    {
                        "ensg": ensg,
                        "ensg_version": ensg_version,
                        "symbol": opt_intern(first(attrs, "gene_name")),
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "biotype": opt_intern(first(attrs, "gene_type")),
                        "mane_select": has_mane,
                        "hgnc_id": opt_intern(first(attrs, "hgnc_id")),
                        "enst": enst,
                        "enst_version": enst_version,
                        "transcript_type": opt_intern(first(attrs, "transcript_type")),
                        "tags": tags,
                    },
                )
                written["transcript"] += 1
                continue
            interval = {
                "ensg": ensg,
                "enst": enst,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "exon_number": exon_number_of(parts[8], bad),
            }
            if feature == "exon":
                append_row(exons, interval)
                written["exon"] += 1
            elif feature == "UTR":
                append_row(utrs, interval)
                written["UTR"] += 1
            else:
                interval["feature"] = feature
                append_row(cds, interval)
                written[feature] += 1

    # Gene rows in this GTF do not carry MANE_Select. The flag is true only
    # when a kept tag on that versioned gene is exactly MANE_Select.
    genes["mane_select"] = [
        version in mane_tx_versions or version in gene_level_mane
        for version in genes["ensg_version"]
    ]

    kept_total = sum(written.values())
    dropped_total = sum(dropped.values())
    skipped_total = sum(skipped.values())
    if kept_total + dropped_total + skipped_total != data_rows:
        raise RuntimeError(
            "row accounting failed: "
            f"kept {kept_total} + dropped {dropped_total} + skipped {skipped_total} "
            f"!= data {data_rows}"
        )

    log("writing parquet")
    written_rows = {
        "genes": write_table(PATHS["genes"], genes, GENE_SCHEMA),
        "transcripts": write_table(PATHS["transcripts"], transcripts, TRANSCRIPT_SCHEMA),
        "exons": write_table(PATHS["exons"], exons, EXON_SCHEMA),
        "cds": write_table(PATHS["cds"], cds, CDS_SCHEMA),
        "utrs": write_table(PATHS["utrs"], utrs, UTR_SCHEMA),
    }
    if written_rows["genes"] != written["gene"]:
        raise RuntimeError("gene parquet row count does not match the parse")
    if written_rows["transcripts"] != written["transcript"]:
        raise RuntimeError("transcript parquet row count does not match the parse")
    if written_rows["exons"] != written["exon"]:
        raise RuntimeError("exon parquet row count does not match the parse")
    if written_rows["cds"] != written["CDS"] + written["stop_codon"]:
        raise RuntimeError("cds parquet row count does not match the parse")
    if written_rows["utrs"] != written["UTR"]:
        raise RuntimeError("utr parquet row count does not match the parse")

    return {
        "header": header,
        "comments": comments,
        "data_rows": data_rows,
        "skipped": skipped,
        "dropped": dropped,
        "failed_chrom": failed_chrom,
        "rewritten": rewritten,
        "written": written,
        "bad_exon_number": bad["exon_number"],
        "mane_on_gene_feature": mane_on_gene_feature,
        "input_bytes": GTF.stat().st_size,
    }


def write_table(path: Path, data: dict[str, list], schema: pa.Schema) -> int:
    table = pa.Table.from_pydict(data, schema=schema)
    pq.write_table(table, path, compression="zstd", compression_level=3)
    log(f"wrote {path.name} rows={table.num_rows:,}")
    return table.num_rows


def sql_path(path: Path) -> str:
    return path.as_posix()


def md_table(headers: list[str], rows: list[list[object]]) -> str:
    head = "| " + " | ".join(headers) + " |"
    sep = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join("" if cell is None else str(cell) for cell in row) + " |" for row in rows]
    return "\n".join([head, sep, *body])


def chrom_sort(name: str) -> tuple:
    body = name[3:] if name.startswith("chr") else name
    if body.isdigit():
        return (0, int(body), name)
    order = {"X": 23, "Y": 24, "M": 25}
    return (1, order.get(body, 99), name)


def null_report(con: duckdb.DuckDBPyConnection, view: str, fields: list[str]) -> list[list[object]]:
    pieces = [f'SUM("{name}" IS NULL)::BIGINT AS "{name}"' for name in fields]
    sql = f'SELECT COUNT(*)::BIGINT AS n, {", ".join(pieces)} FROM {view}'
    row = con.execute(sql).fetchone()
    if row is None:
        raise RuntimeError(f"null query failed for {view}")
    total = int(row[0])
    lines = []
    for index, name in enumerate(fields, start=1):
        nulls = int(row[index] or 0)
        if total == 0:
            rate = "n/a"
        elif nulls == 0:
            rate = "0"
        else:
            rate = f"{nulls / total:.6f}"
        lines.append([name, f"{nulls:,}", rate])
    return lines


def open_views(con: duckdb.DuckDBPyConnection) -> None:
    for name, path in PATHS.items():
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{sql_path(path)}')")


def build_qc(stats: dict) -> str:
    con = duckdb.connect(database=":memory:")
    open_views(con)

    def one(sql: str, params: list | None = None):
        row = con.execute(sql, params or []).fetchone()
        if row is None:
            raise RuntimeError(sql)
        return row

    gene_n = int(one("SELECT COUNT(*) FROM genes")[0])
    tx_n = int(one("SELECT COUNT(*) FROM transcripts")[0])
    exon_n = int(one("SELECT COUNT(*) FROM exons")[0])
    cds_n = int(one("SELECT COUNT(*) FROM cds")[0])
    utr_n = int(one("SELECT COUNT(*) FROM utrs")[0])
    protein_n = int(one("SELECT COUNT(*) FROM genes WHERE biotype = 'protein_coding'")[0])
    null_ensg = int(
        one(
            """
            SELECT
              (SELECT COUNT(*) FROM genes WHERE ensg IS NULL)
            + (SELECT COUNT(*) FROM transcripts WHERE ensg IS NULL)
            + (SELECT COUNT(*) FROM exons WHERE ensg IS NULL)
            + (SELECT COUNT(*) FROM cds WHERE ensg IS NULL)
            + (SELECT COUNT(*) FROM utrs WHERE ensg IS NULL)
            """
        )[0]
    )
    bad_coord = int(
        one(
            """
            SELECT
              (SELECT COUNT(*) FROM genes WHERE start < 1 OR "end" < start)
            + (SELECT COUNT(*) FROM transcripts WHERE start < 1 OR "end" < start)
            + (SELECT COUNT(*) FROM exons WHERE start < 1 OR "end" < start)
            + (SELECT COUNT(*) FROM cds WHERE start < 1 OR "end" < start)
            + (SELECT COUNT(*) FROM utrs WHERE start < 1 OR "end" < start)
            """
        )[0]
    )
    gene_dup_version = int(one("SELECT COUNT(*) - COUNT(DISTINCT ensg_version) FROM genes")[0])
    gene_dup_bare = int(one("SELECT COUNT(*) - COUNT(DISTINCT ensg) FROM genes")[0])
    tx_dup_version = int(one("SELECT COUNT(*) - COUNT(DISTINCT enst_version) FROM transcripts")[0])
    tx_dup_bare = int(one("SELECT COUNT(*) - COUNT(DISTINCT enst) FROM transcripts")[0])
    mane_tx = int(one("SELECT COUNT(*) FROM transcripts WHERE list_contains(tags, 'MANE_Select')")[0])
    mane_flag_mismatch = int(
        one(
            """
            SELECT COUNT(*) FROM transcripts
            WHERE mane_select <> list_contains(tags, 'MANE_Select')
            """
        )[0]
    )
    mane_gene_rows = int(one("SELECT COUNT(*) FROM genes WHERE mane_select")[0])
    mane_gene_distinct = int(
        one(
            """
            SELECT COUNT(DISTINCT ensg_version) FROM transcripts
            WHERE list_contains(tags, 'MANE_Select')
            """
        )[0]
    )
    flag_without_tx = int(
        one(
            """
            SELECT COUNT(*) FROM genes g
            WHERE g.mane_select
              AND g.ensg_version NOT IN (
                SELECT ensg_version FROM transcripts
                WHERE list_contains(tags, 'MANE_Select')
              )
            """
        )[0]
    )
    tx_without_flag = int(
        one(
            """
            WITH mane_genes AS (
              SELECT DISTINCT ensg_version FROM transcripts
              WHERE list_contains(tags, 'MANE_Select')
            )
            SELECT COUNT(*) FROM mane_genes m
            WHERE m.ensg_version NOT IN (SELECT ensg_version FROM genes WHERE mane_select)
            """
        )[0]
    )
    empty_tags = int(one("SELECT COUNT(*) FROM transcripts WHERE length(tags) = 0")[0])
    odd_hgnc = int(
        one(
            """
            SELECT COUNT(*) FROM genes
            WHERE hgnc_id IS NOT NULL AND NOT starts_with(hgnc_id, 'HGNC:')
            """
        )[0]
    )
    odd_strand = con.execute(
        """
        SELECT strand, COUNT(*)::BIGINT
        FROM genes
        WHERE strand NOT IN ('+', '-')
        GROUP BY strand
        ORDER BY strand
        """
    ).fetchall()
    biotypes = con.execute(
        """
        SELECT COALESCE(biotype, '(null)'), COUNT(*)::BIGINT
        FROM genes
        GROUP BY biotype
        ORDER BY COUNT(*) DESC, biotype
        """
    ).fetchall()
    protein_like = [row for row in biotypes if "protein" in str(row[0])]
    feature_counts = con.execute(
        """
        SELECT feature, COUNT(*)::BIGINT
        FROM cds
        GROUP BY feature
        ORDER BY feature
        """
    ).fetchall()
    tag_counts = con.execute(
        """
        SELECT tag, COUNT(*)::BIGINT AS n
        FROM (SELECT UNNEST(tags) AS tag FROM transcripts)
        WHERE tag IN ('MANE_Select', 'MANE_Plus_Clinical', 'Ensembl_canonical')
           OR starts_with(tag, 'appris_')
        GROUP BY tag
        ORDER BY tag
        """
    ).fetchall()
    mane_not_coding = int(
        one(
            """
            SELECT COUNT(*) FROM transcripts
            WHERE list_contains(tags, 'MANE_Select')
              AND transcript_type IS DISTINCT FROM 'protein_coding'
            """
        )[0]
    )
    tx_missing_gene = int(
        one(
            """
            SELECT COUNT(*) FROM transcripts t
            WHERE t.ensg_version NOT IN (SELECT ensg_version FROM genes)
            """
        )[0]
    )
    genes_without_tx = int(
        one(
            """
            SELECT COUNT(*) FROM genes g
            WHERE g.ensg_version NOT IN (SELECT ensg_version FROM transcripts)
            """
        )[0]
    )
    child_miss = {}
    ensg_mismatch = {}
    for view in ("exons", "cds", "utrs"):
        child_miss[view] = int(
            one(
                f"""
                SELECT COUNT(*) FROM {view} c
                WHERE NOT EXISTS (
                  SELECT 1 FROM transcripts t
                  WHERE t.enst = c.enst AND t.ensg = c.ensg
                )
                """
            )[0]
        )
        ensg_mismatch[view] = int(
            one(
                f"""
                SELECT COUNT(*) FROM {view} c
                JOIN transcripts t ON c.enst = t.enst
                WHERE c.ensg <> t.ensg
                """
            )[0]
        )
    dup_symbol_n = int(
        one(
            """
            SELECT COUNT(*) FROM (
              SELECT symbol FROM genes
              WHERE symbol IS NOT NULL
              GROUP BY symbol
              HAVING COUNT(DISTINCT ensg) > 1
            )
            """
        )[0]
    )
    dup_examples = con.execute(
        """
        SELECT symbol, COUNT(DISTINCT ensg)::BIGINT AS genes, COUNT(DISTINCT chrom)::BIGINT AS chroms
        FROM genes
        WHERE symbol IS NOT NULL
        GROUP BY symbol
        HAVING COUNT(DISTINCT ensg) > 1
        ORDER BY genes DESC, symbol
        LIMIT 10
        """
    ).fetchall()
    stored_chroms: dict[str, dict[str, int]] = {}
    noncanon_out: dict[str, list[str]] = {}
    for name in PATHS:
        rows = con.execute(f"SELECT chrom, COUNT(*)::BIGINT FROM {name} GROUP BY chrom").fetchall()
        for chrom, count in rows:
            stored_chroms.setdefault(chrom, {})[name] = int(count)
            if chrom not in CANON:
                noncanon_out.setdefault(chrom, []).append(name)
    no_version = int(
        one(
            """
            SELECT
              (SELECT COUNT(*) FROM genes WHERE ensg_version NOT LIKE '%.%')
            + (SELECT COUNT(*) FROM transcripts WHERE enst_version NOT LIKE '%.%')
            """
        )[0]
    )
    truth_lines = truth_checks(con)

    parts: list[str] = []
    add = parts.append
    add("# Interval QC")
    add("")
    add(
        f"GENCODE v49 intervals from `{GTF.relative_to(ROOT).as_posix()}` "
        f"({stats['input_bytes']:,} bytes gzip). "
        f"Gene rows: {gene_n:,} (baseline {GENE_BASELINE:,}, delta {gene_n - GENE_BASELINE:+,}). "
        f"Protein-coding (`biotype = protein_coding`): {protein_n:,} "
        f"(baseline {PROTEIN_BASELINE:,}, delta {protein_n - PROTEIN_BASELINE:+,}). "
        f"Transcript rows with tag `MANE_Select`: {mane_tx:,}. "
        f"Genes with `mane_select` true: {mane_gene_rows:,}."
    )
    add("")
    if stats["failed_chrom"]:
        failed_list = ", ".join(
            f"`{name}` ({count:,} rows)" for name, count in sorted(stats["failed_chrom"].items())
        )
        add(f"Chromosomes that failed canonicalization: {failed_list}. Those rows were not written.")
    else:
        add(
            "Chromosomes that failed canonicalization: none. "
            "Every kept GTF chromosome passed `canon_chrom`, and every `chrom` value "
            "in the five Parquet files is one of chr1–chr22, chrX, chrY, or chrM."
        )
    add("")
    add("## Source and rules")
    add("")
    add("Header lines from the GTF:")
    add("")
    add("```text")
    for line in stats["header"]:
        add(line)
    add("```")
    add("")
    add(
        "GTF `start` and `end` are stored as 1-based inclusive integers. "
        "They were not converted to 0-based half-open coordinates. "
        "`ensg` and `enst` are version-stripped. `ensg_version` and `enst_version` "
        "are the original ids, including the version suffix when the GTF has one. "
        f"Rows whose version field has no `.` suffix: {no_version:,}."
    )
    add("")
    add(
        "Feature types kept: gene, transcript, exon, CDS, UTR, stop_codon. "
        "`stop_codon` rows are in `cds.parquet` with `feature = stop_codon`. "
        "CDS rows have `feature = CDS`. Exon and UTR files do not repeat those features. "
        "`start_codon` and `Selenocysteine` were not kept."
    )
    add("")
    add(
        "Gene feature rows in this file do not carry `MANE_Select` "
        f"({stats['mane_on_gene_feature']:,} gene rows had that tag). "
        "`genes.mane_select` is true only when a transcript of that versioned gene id "
        "has a tag whose value is exactly `MANE_Select`, or when the gene row itself has that tag. "
        "`transcripts.mane_select` is true only when that transcript row's own tags contain "
        "`MANE_Select`. `MANE_Plus_Clinical` is not treated as MANE Select. Tags were not invented."
    )
    add("")
    add("Counts below are DuckDB queries against the Parquet files just written, except the read-time drop table.")
    add("")
    add("## Row counts")
    add("")
    add(
        md_table(
            ["file", "rows", "bytes"],
            [
                [path.name, f"{con.execute(f'SELECT COUNT(*) FROM {name}').fetchone()[0]:,}", f"{path.stat().st_size:,}"]
                for name, path in PATHS.items()
            ],
        )
    )
    add("")
    add("CDS file by feature:")
    add("")
    add(md_table(["feature", "rows"], [[name, f"{count:,}"] for name, count in feature_counts]))
    add("")
    add(
        f"Parse versus file: gene {stats['written']['gene']:,}, "
        f"transcript {stats['written']['transcript']:,}, "
        f"exon {stats['written']['exon']:,}, "
        f"CDS {stats['written']['CDS']:,}, "
        f"stop_codon {stats['written']['stop_codon']:,}, "
        f"UTR {stats['written']['UTR']:,}. "
        f"These match the Parquet row counts "
        f"({gene_n:,}, {tx_n:,}, {exon_n:,}, {cds_n:,}, {utr_n:,}; "
        f"CDS plus stop_codon = {stats['written']['CDS'] + stats['written']['stop_codon']:,})."
    )
    add("")
    add("## Null rates")
    add("")
    add(f"Null `ensg` across all five files: {null_ensg:,}.")
    add("")
    for view, schema in (
        ("genes", GENE_SCHEMA),
        ("transcripts", TRANSCRIPT_SCHEMA),
        ("exons", EXON_SCHEMA),
        ("cds", CDS_SCHEMA),
        ("utrs", UTR_SCHEMA),
    ):
        add(f"### {view}")
        add("")
        add(md_table(["column", "nulls", "null_rate"], null_report(con, view, [field.name for field in schema])))
        add("")
    add(
        f"Transcript rows with an empty tag list: {empty_tags:,}. "
        f"`hgnc_id` values that do not start with `HGNC:`: {odd_hgnc:,}. "
        f"Non-integer `exon_number` values stored as null: {stats['bad_exon_number']:,}."
    )
    add("")
    if odd_strand:
        add("Gene strands other than `+` or `-`:")
        add("")
        add(md_table(["strand", "rows"], [[strand, f"{count:,}"] for strand, count in odd_strand]))
    else:
        add("Every gene strand is `+` or `-`.")
    add("")
    add(f"Rows with `start < 1` or `end < start`, all files: {bad_coord:,}.")
    add("")
    add("## Gene count versus 78691")
    add("")
    add(
        f"`genes.parquet` has {gene_n:,} rows. The SQLite gene-table baseline in the plan is "
        f"{GENE_BASELINE:,}. Delta: {gene_n - GENE_BASELINE:+,}."
    )
    add("")
    add(
        f"Duplicate versioned gene ids: {gene_dup_version:,}. "
        f"Duplicate bare `ensg`: {gene_dup_bare:,}. "
        f"Duplicate versioned transcript ids: {tx_dup_version:,}. "
        f"Duplicate bare `enst`: {tx_dup_bare:,}."
    )
    add("")
    add("Gene biotypes:")
    add("")
    add(md_table(["biotype", "genes"], [[name, f"{count:,}"] for name, count in biotypes]))
    add("")
    add("## Protein-coding versus 20097")
    add("")
    add(
        f"Genes with `biotype = protein_coding`: {protein_n:,}. "
        f"Baseline: {PROTEIN_BASELINE:,}. Delta: {protein_n - PROTEIN_BASELINE:+,}."
    )
    add("")
    add("Biotypes whose name contains `protein`:")
    add("")
    add(md_table(["biotype", "genes"], [[name, f"{count:,}"] for name, count in protein_like]))
    add("")
    add(
        "The baseline comparison uses the exact biotype `protein_coding` only. "
        "Other protein-like biotypes are not folded into that count."
    )
    add("")
    add("## MANE tag count")
    add("")
    if mane_tx == 0 and stats["mane_on_gene_feature"] == 0:
        add(
            "The GTF has no `MANE_Select` tag on a kept feature. "
            "No tag was invented. The constraint file remains the MANE source."
        )
    else:
        add(
            f"The GTF does contain `MANE_Select`. Transcript rows with that exact tag: {mane_tx:,}. "
            f"Distinct versioned gene ids on those transcripts: {mane_gene_distinct:,}. "
            f"`genes.mane_select` true: {mane_gene_rows:,}. "
            f"Gene-feature rows with the tag: {stats['mane_on_gene_feature']:,}."
        )
    add("")
    add(
        f"Transcripts where `mane_select` disagrees with `list_contains(tags, 'MANE_Select')`: "
        f"{mane_flag_mismatch:,}. "
        f"Gene flags with no MANE transcript: {flag_without_tx:,}. "
        f"MANE transcripts whose gene flag is false: {tx_without_flag:,}. "
        f"MANE transcripts whose `transcript_type` is not `protein_coding`: {mane_not_coding:,}."
    )
    add("")
    add("Tag occurrences on transcript rows (unnested; one transcript contributes one if the tag is present once):")
    add("")
    add(md_table(["tag", "occurrences"], [[tag, f"{count:,}"] for tag, count in tag_counts]))
    add("")
    add("Predicate used: `list_contains(tags, 'MANE_Select')`.")
    add("")
    add("## Chromosomes that failed canonicalization")
    add("")
    if stats["failed_chrom"] or noncanon_out:
        if stats["failed_chrom"]:
            add("Rejected while reading kept features:")
            add("")
            add(
                md_table(
                    ["raw_chrom", "rows_dropped"],
                    [[name, f"{count:,}"] for name, count in sorted(stats["failed_chrom"].items())],
                )
            )
            add("")
        if noncanon_out:
            add("Non-canonical `chrom` values found in the Parquet files:")
            add("")
            add(
                md_table(
                    ["chrom", "files"],
                    [[name, ", ".join(files)] for name, files in sorted(noncanon_out.items())],
                )
            )
            add("")
    else:
        add("None. No kept row was dropped for chromosome name, and no Parquet `chrom` is outside the canonical set.")
    add("")
    if stats["rewritten"]:
        add("Chromosomes rewritten by `canon_chrom` (not failures):")
        add("")
        add(
            md_table(
                ["rewrite", "rows"],
                [[name, f"{count:,}"] for name, count in sorted(stats["rewritten"].items())],
            )
        )
        add("")
    else:
        add("No chromosome name was rewritten. Raw GTF names and stored `chrom` values are the same.")
    add("")
    add("Stored chromosome counts:")
    add("")
    chrom_rows = []
    for chrom in sorted(stored_chroms, key=chrom_sort):
        counts = stored_chroms[chrom]
        chrom_rows.append(
            [
                chrom,
                f"{counts.get('genes', 0):,}",
                f"{counts.get('transcripts', 0):,}",
                f"{counts.get('exons', 0):,}",
                f"{counts.get('cds', 0):,}",
                f"{counts.get('utrs', 0):,}",
            ]
        )
    add(md_table(["chrom", "genes", "transcripts", "exons", "cds", "utrs"], chrom_rows))
    add("")
    add("## Read-time drops")
    add("")
    add(
        f"Comment lines: {stats['comments']:,}. Data rows: {stats['data_rows']:,}. "
        "Skipped features are not drops."
    )
    add("")
    add(
        md_table(
            ["feature", "rows_skipped"],
            [[name, f"{count:,}"] for name, count in sorted(stats["skipped"].items())],
        )
    )
    add("")
    if stats["dropped"]:
        add(
            md_table(
                ["reason", "rows_dropped"],
                [[name, f"{count:,}"] for name, count in sorted(stats["dropped"].items())],
            )
        )
    else:
        add("Rows dropped for chromosome, id, or coordinate: 0.")
    add("")
    add("## Join checks inside this lane")
    add("")
    add(
        f"Transcript rows whose `ensg_version` is absent from genes: {tx_missing_gene:,}. "
        f"Genes with no transcript: {genes_without_tx:,}."
    )
    add("")
    add(
        md_table(
            ["file", "rows_without_matching_transcript", "rows_whose_ensg_disagrees_with_transcript"],
            [
                [view, f"{child_miss[view]:,}", f"{ensg_mismatch[view]:,}"]
                for view in ("exons", "cds", "utrs")
            ],
        )
    )
    add("")
    add(
        f"Symbols shared by more than one bare `ensg`: {dup_symbol_n:,}. "
        "Symbol is not the primary key. Later joins should use `ensg`."
    )
    add("")
    if dup_examples:
        add(md_table(["symbol", "genes", "chroms"], [[sym, f"{n:,}", f"{c:,}"] for sym, n, c in dup_examples]))
        add("")
    add("## Truth-set interval checks")
    add("")
    add(
        "These are the truth-set cases this lane can test from coordinates and tags. "
        "Constraint, ClinVar, GWAS, and gnomAD comparisons belong to other lanes. "
        "Positions are 1-based and match a feature when `start <= pos <= end`."
    )
    add("")
    add(md_table(["case", "result", "detail"], truth_lines))
    add("")
    con.close()
    return "\n".join(parts) + "\n"


def truth_checks(con: duckdb.DuckDBPyConnection) -> list[list[object]]:
    if not TRUTH.is_file():
        return [["(missing)", "not run", str(TRUTH)]]
    with TRUTH.open(encoding="utf-8", newline="") as handle:
        cases = {row["case_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    lines: list[list[object]] = []
    lines.append(check_coding_overlap(con, cases.get("coding_snp"), "coding_snp"))
    lines.append(check_coding_overlap(con, cases.get("cftr_indel"), "cftr_indel"))
    lines.append(check_intergenic(con, cases.get("gwas_intergenic")))
    lines.append(check_crosswalk(con, cases.get("chrom_crosswalk")))
    lines.append(check_mane(con, cases.get("mane_select")))
    return lines


def check_coding_overlap(con: duckdb.DuckDBPyConnection, row: dict | None, case_id: str) -> list[object]:
    if row is None or not row.get("pos"):
        return [case_id, "not run", "case or position missing from truth_set.tsv"]
    chrom = canon_chrom(row["chrom_raw"])
    pos = int(row["pos"])
    symbol = row["gene"]
    hits = con.execute(
        """
        SELECT COUNT(*) FROM cds
        WHERE feature = 'CDS' AND chrom = ? AND start <= ? AND "end" >= ?
          AND ensg IN (SELECT ensg FROM genes WHERE symbol = ?)
        """,
        [chrom, pos, pos, symbol],
    ).fetchone()
    any_cds = con.execute(
        """
        SELECT COUNT(*) FROM cds
        WHERE feature = 'CDS' AND chrom = ? AND start <= ? AND "end" >= ?
        """,
        [chrom, pos, pos],
    ).fetchone()
    gene_n = con.execute("SELECT COUNT(*) FROM genes WHERE symbol = ?", [symbol]).fetchone()
    covered = int(hits[0]) if hits else 0
    any_n = int(any_cds[0]) if any_cds else 0
    genes = int(gene_n[0]) if gene_n else 0
    result = "PASS" if covered >= 1 else "FAIL"
    detail = (
        f"{symbol} genes {genes}; CDS rows of that symbol covering {chrom}:{pos}: {covered}; "
        f"CDS rows of any gene covering it: {any_n}"
    )
    return [case_id, result, detail]


def check_intergenic(con: duckdb.DuckDBPyConnection, row: dict | None) -> list[object]:
    if row is None or not row.get("pos"):
        return ["gwas_intergenic", "not run", "case or position missing from truth_set.tsv"]
    chrom = canon_chrom(row["chrom_raw"])
    pos = int(row["pos"])
    cds_n = scalar_overlap(con, "cds", chrom, pos, "feature = 'CDS'")
    stop_n = scalar_overlap(con, "cds", chrom, pos, "feature = 'stop_codon'")
    exon_n = scalar_overlap(con, "exons", chrom, pos, None)
    utr_n = scalar_overlap(con, "utrs", chrom, pos, None)
    gene_n = scalar_overlap(con, "genes", chrom, pos, None)
    symbols = con.execute(
        """
        SELECT symbol FROM genes
        WHERE chrom = ? AND start <= ? AND "end" >= ?
        ORDER BY symbol
        """,
        [chrom, pos, pos],
    ).fetchall()
    symbol_text = ", ".join(item[0] or "(null)" for item in symbols) or "(none)"
    result = "PASS" if cds_n == 0 else "FAIL"
    detail = (
        f"{chrom}:{pos} CDS {cds_n}, stop_codon {stop_n}, exon {exon_n}, "
        f"UTR {utr_n}, gene span {gene_n} ({symbol_text})"
    )
    return ["gwas_intergenic", result, detail]


def scalar_overlap(con: duckdb.DuckDBPyConnection, view: str, chrom: str, pos: int, extra: str | None) -> int:
    where = 'chrom = ? AND start <= ? AND "end" >= ?'
    if extra:
        where = f"{where} AND {extra}"
    row = con.execute(f"SELECT COUNT(*) FROM {view} WHERE {where}", [chrom, pos, pos]).fetchone()
    return int(row[0]) if row else 0


def check_crosswalk(con: duckdb.DuckDBPyConnection, row: dict | None) -> list[object]:
    if row is None:
        return ["chrom_crosswalk", "not run", "case missing from truth_set.tsv"]
    symbol = row["gene"]
    pos = int(row["pos"]) if row.get("pos") else None
    match = _SPAN.search(row.get("notes") or "")
    found = con.execute(
        """
        SELECT ensg_version, chrom, start, "end", biotype, mane_select
        FROM genes
        WHERE symbol = ?
        ORDER BY ensg_version
        """,
        [symbol],
    ).fetchall()
    if not found:
        return ["chrom_crosswalk", "FAIL", f"no gene with symbol {symbol}"]
    described = "; ".join(
        f"{ensg} {chrom}:{start}-{end} {biotype} mane_select={mane}"
        for ensg, chrom, start, end, biotype, mane in found
    )
    if match is None:
        inside = ""
        if pos is not None and found:
            chrom, start, end = found[0][1], found[0][2], found[0][3]
            inside = f"; pos {pos} inside first row: {start <= pos <= end}"
        return ["chrom_crosswalk", "FAIL", f"could not parse the truth note; found {described}{inside}"]
    ensg_version, chrom, start_s, end_s, biotype = match.groups()
    start = int(start_s)
    end = int(end_s)
    ok = [
        item
        for item in found
        if item[0] == ensg_version
        and item[1] == chrom
        and item[2] == start
        and item[3] == end
        and item[4] == biotype
    ]
    pos_ok = False
    if ok and pos is not None:
        pos_ok = ok[0][2] <= pos <= ok[0][3]
    result = "PASS" if ok and (pos is None or pos_ok) else "FAIL"
    detail = (
        f"expected {ensg_version} {chrom}:{start}-{end} {biotype}; "
        f"pos {pos} inside span: {pos_ok}; found {described}"
    )
    return ["chrom_crosswalk", result, detail]


def check_mane(con: duckdb.DuckDBPyConnection, row: dict | None) -> list[object]:
    if row is None:
        return ["mane_select", "not run", "case missing from truth_set.tsv"]
    symbol = row["gene"]
    bare = bare_enst(row["expected"]) or row["expected"]
    ensg_note = None
    note_match = re.search(r"ENSG\d+", row.get("notes") or "")
    if note_match:
        ensg_note = note_match.group(0)
    found = con.execute(
        """
        SELECT ensg, ensg_version, enst, enst_version, symbol, chrom,
               mane_select, list_contains(tags, 'MANE_Select')
        FROM transcripts
        WHERE enst = ?
        ORDER BY enst_version
        """,
        [bare],
    ).fetchall()
    if not found:
        return ["mane_select", "FAIL", f"no transcript {bare}"]
    described = "; ".join(
        f"{ensg_version} {enst_version} {sym} {chrom} mane_select={flag} tag={tag}"
        for ensg, ensg_version, enst, enst_version, sym, chrom, flag, tag in found
    )
    ok = [
        item
        for item in found
        if item[4] == symbol and item[6] is True and item[7] is True and (ensg_note is None or item[0] == ensg_note)
    ]
    gene_flag = None
    if ok:
        flag_row = con.execute(
            "SELECT mane_select FROM genes WHERE ensg_version = ?",
            [ok[0][1]],
        ).fetchone()
        gene_flag = None if flag_row is None else bool(flag_row[0])
    result = "PASS" if ok and gene_flag is True else "FAIL"
    detail = (
        f"expected {symbol} {bare} tag MANE_Select"
        + (f" gene {ensg_note}" if ensg_note else "")
        + f"; gene.mane_select={gene_flag}; found {described}. "
        "The constraint-file half of this case is not this lane."
    )
    return ["mane_select", result, detail]


def main() -> int:
    _check_parser()
    _check_duckdb()
    stats = load()
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)
    QC_PATH.write_text(build_qc(stats), encoding="utf-8", newline="\n")
    log(f"wrote {QC_PATH}")
    con = duckdb.connect(database=":memory:")
    missing = con.execute(
        """
        SELECT
          (SELECT COUNT(*) FROM read_parquet(?) WHERE exon_number IS NULL),
          (SELECT COUNT(*) FROM read_parquet(?)),
          (SELECT COUNT(*) FROM read_parquet(?) WHERE exon_number IS NULL),
          (SELECT COUNT(*) FROM read_parquet(?)),
          (SELECT COUNT(*) FROM read_parquet(?) WHERE exon_number IS NULL),
          (SELECT COUNT(*) FROM read_parquet(?))
        """,
        [
            sql_path(PATHS["exons"]),
            sql_path(PATHS["exons"]),
            sql_path(PATHS["cds"]),
            sql_path(PATHS["cds"]),
            sql_path(PATHS["utrs"]),
            sql_path(PATHS["utrs"]),
        ],
    ).fetchone()
    con.close()
    if missing is None:
        raise RuntimeError("exon_number check failed")
    pairs = (("exons", missing[0], missing[1]), ("cds", missing[2], missing[3]), ("utrs", missing[4], missing[5]))
    for name, nulls, rows in pairs:
        if rows and nulls == rows:
            raise RuntimeError(f"exon_number is null on every {name} row")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
