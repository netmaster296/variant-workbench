# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load gnomAD v4.1 constraint metrics into transcript and gene Parquet.

Transcript file keeps every source row. Gene file keeps one row per gene.
Preference is mane_select, then canonical, then longest cds_length.
A RefSeq transcript and an Ensembl transcript that still tie: keep Ensembl.
"""

from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import coords  # noqa: E402

SOURCE = ROOT / "raw" / "gnomad" / "constraint" / "gnomad.v4.1.constraint_metrics.tsv"
TRANSCRIPT_OUT = ROOT / "warehouse" / "constraint" / "transcript_constraint.parquet"
GENE_OUT = ROOT / "warehouse" / "constraint" / "gene_constraint.parquet"
QC_OUT = ROOT / "warehouse" / "qc" / "constraint.md"
TRUTH = ROOT / "warehouse" / "qc" / "truth_set.tsv"
INTERVALS = ROOT / "warehouse" / "intervals" / "genes.parquet"

EXPECTED_HEADER = [
    "gene",
    "gene_id",
    "transcript",
    "canonical",
    "mane_select",
    "lof_hc_lc.obs",
    "lof_hc_lc.exp",
    "lof_hc_lc.possible",
    "lof_hc_lc.oe",
    "lof_hc_lc.mu",
    "lof_hc_lc.pLI",
    "lof_hc_lc.pNull",
    "lof_hc_lc.pRec",
    "lof.obs",
    "lof.exp",
    "lof.possible",
    "lof.oe",
    "lof.mu",
    "lof.pLI",
    "lof.pNull",
    "lof.pRec",
    "lof.oe_ci.lower",
    "lof.oe_ci.upper",
    "lof.oe_ci.upper_rank",
    "lof.oe_ci.upper_bin_decile",
    "lof.z_raw",
    "lof.z_score",
    "mis.obs",
    "mis.exp",
    "mis.possible",
    "mis.oe",
    "mis.mu",
    "mis.oe_ci.lower",
    "mis.oe_ci.upper",
    "mis.z_raw",
    "mis.z_score",
    "mis_pphen.obs",
    "mis_pphen.exp",
    "mis_pphen.possible",
    "mis_pphen.oe",
    "syn.obs",
    "syn.exp",
    "syn.possible",
    "syn.oe",
    "syn.mu",
    "syn.oe_ci.lower",
    "syn.oe_ci.upper",
    "syn.z_raw",
    "syn.z_score",
    "constraint_flags",
    "level",
    "transcript_type",
    "chromosome",
    "cds_length",
    "num_coding_exons",
]

STRING_COLS = {
    "gene",
    "gene_id",
    "transcript",
    "constraint_flags",
    "transcript_type",
    "chromosome",
}
BOOL_COLS = {"canonical", "mane_select"}
MISSING = {"", "NA", "NaN", "nan", ".", "-", "None", "none"}
TRUTH_CASES = (
    "mane_select",
    "loeuf_low_1",
    "loeuf_low_2",
    "loeuf_low_3",
    "loeuf_near_1_1",
    "loeuf_near_1_2",
    "loeuf_near_1_3",
)
LOEUF_TOLERANCE = 0.001


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if text in MISSING:
        return None
    return float(text)


def parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    text = value.strip().lower()
    if text in MISSING or text == "na":
        return None
    if text == "true":
        return True
    if text == "false":
        return False
    raise ValueError(f"expected true/false, got {value!r}")


def clean_text(value: str | None, *, missing_as_null: bool) -> str | None:
    if value is None:
        return None
    text = value.strip()
    if missing_as_null and text in MISSING:
        return None
    return text


def preference_key(row: dict) -> tuple:
    """Higher wins: MANE, canonical, measured CDS length, then Ensembl on a tie.

    A missing cds_length is not a length, so it loses to any measured CDS.
    RefSeq rows in this file have no cds_length. When both sides still tie,
    bare_enst selects the Ensembl transcript.
    """
    cds = row["cds_length"]
    has_cds = cds is not None
    return (
        row["mane_select"] is True,
        row["canonical"] is True,
        has_cds,
        cds if has_cds else -1.0,
        row["bare_enst"] is not None,
        row["transcript"] or "",
    )


def load_rows(path: Path) -> tuple[list[str], list[dict], list[str]]:
    notes: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if header != EXPECTED_HEADER:
            missing = [name for name in EXPECTED_HEADER if name not in header]
            extra = [name for name in header if name not in EXPECTED_HEADER]
            notes.append(
                "Header differs from plan.md. This loader follows the file. "
                f"Missing from file: {missing or 'none'}. "
                f"Extra in file: {extra or 'none'}."
            )
        else:
            notes.append(
                f"Header has {len(header)} columns and matches the gnomAD constraint "
                "schema in plan.md. No column was renamed."
            )
        rows: list[dict] = []
        for rec in reader:
            row: dict = {}
            for name in header:
                raw = rec[name]
                if name in BOOL_COLS:
                    row[name] = parse_bool(raw)
                elif name in STRING_COLS:
                    if name == "constraint_flags":
                        row[name] = clean_text(raw, missing_as_null=True)
                    elif name == "chromosome":
                        text = clean_text(raw, missing_as_null=True)
                        canon = coords.canon_chrom(text) if text else None
                        row[name] = canon if text is None or canon is not None else text
                    else:
                        row[name] = clean_text(raw, missing_as_null=True)
                else:
                    row[name] = parse_float(raw)
            row["bare_ensg"] = coords.bare_ensg(row.get("gene_id"))
            row["bare_enst"] = coords.bare_enst(row.get("transcript"))
            rows.append(row)
    return header, rows, notes


def select_genes(rows: list[dict]) -> tuple[list[dict], set[str]]:
    """One gene row. NCBI and Ensembl ids of one symbol compete when the symbol has a single ENSG.

    A symbol with two Ensembl ids is two genes (overlapping or successive models).
    A missing symbol is not a group key. Those rows stay one gene per gene_id.
    """
    by_symbol: dict[str, list[dict]] = defaultdict(list)
    no_symbol: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        symbol = row["gene"]
        if symbol is None:
            no_symbol[row["gene_id"]].append(row)
        else:
            by_symbol[symbol].append(row)

    groups: list[list[dict]] = []
    refseq_hold: list[tuple[str, list[dict]]] = []
    for items in no_symbol.values():
        groups.append(items)
    for symbol, items in by_symbol.items():
        ens_ids = {row["gene_id"] for row in items if row["bare_ensg"] is not None}
        if len(ens_ids) <= 1:
            groups.append(items)
            continue
        per_gene: dict[str, list[dict]] = defaultdict(list)
        refs: list[dict] = []
        for row in items:
            if row["bare_ensg"] is not None:
                per_gene[row["gene_id"]].append(row)
            else:
                refs.append(row)
        for gene_rows in per_gene.values():
            groups.append(gene_rows)
        if refs:
            refseq_hold.append((symbol, refs))

    selected: list[dict] = []
    for items in groups:
        winner = dict(max(items, key=preference_key))
        ens_ids = {row["bare_ensg"] for row in items if row["bare_ensg"]}
        if len(ens_ids) == 1:
            bare = next(iter(ens_ids))
            ensg_row = next(row for row in items if row["bare_ensg"] == bare)
            winner["gene_id"] = ensg_row["gene_id"]
            winner["bare_ensg"] = bare
        winner["loeuf"] = winner["lof.oe_ci.upper"]
        winner["loeuf_decile"] = winner["lof.oe_ci.upper_bin_decile"]
        selected.append(winner)

    by_gene_symbol = defaultdict(list)
    for row in selected:
        if row["gene"] is not None:
            by_gene_symbol[row["gene"]].append(row)
    for symbol, refs in refseq_hold:
        per_id: dict[str, list[dict]] = defaultdict(list)
        for row in refs:
            per_id[row["gene_id"]].append(row)
        best_per_id = [max(items, key=preference_key) for items in per_id.values()]
        best = max(best_per_id, key=preference_key)
        ensembl_keys = [preference_key(row) for row in by_gene_symbol.get(symbol, [])]
        if ensembl_keys and preference_key(best) > max(ensembl_keys):
            kept = dict(best)
            kept["loeuf"] = kept["lof.oe_ci.upper"]
            kept["loeuf_decile"] = kept["lof.oe_ci.upper_bin_decile"]
            selected.append(kept)

    selected.sort(key=lambda row: (row["gene"] or "", row["gene_id"] or "", row["transcript"] or ""))
    input_ids = {row["gene_id"] for row in rows}
    output_ids = {row["gene_id"] for row in selected}
    if len(output_ids) != len(selected):
        raise RuntimeError("gene selection produced duplicate gene_id values")
    if None in output_ids:
        raise RuntimeError("gene selection produced a null gene_id")
    return selected, input_ids - output_ids


def column_types(header: list[str], extra: list[str]) -> list[tuple[str, str]]:
    typed: list[tuple[str, str]] = []
    for name in header:
        if name in BOOL_COLS:
            typed.append((name, "bool"))
        elif name in STRING_COLS:
            typed.append((name, "string"))
        else:
            typed.append((name, "float"))
    for name in extra:
        if name in {"bare_ensg", "bare_enst"}:
            typed.append((name, "string"))
        else:
            typed.append((name, "float"))
    return typed


def write_parquet(path: Path, rows: list[dict], columns: list[tuple[str, str]]) -> None:
    arrays = []
    names = []
    for name, kind in columns:
        values = [row[name] for row in rows]
        if kind == "float":
            arrays.append(pa.array(values, type=pa.float64()))
        elif kind == "bool":
            arrays.append(pa.array(values, type=pa.bool_()))
        else:
            arrays.append(pa.array(values, type=pa.string()))
        names.append(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.table(arrays, names=names), path, compression="zstd")


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def parquet_profile(con: duckdb.DuckDBPyConnection, path: Path, columns: list[str]) -> dict:
    pieces = ["count(*) AS n"]
    for index, name in enumerate(columns):
        pieces.append(f"count({quote_ident(name)}) AS c{index}")
    row = con.execute(
        f"SELECT {', '.join(pieces)} FROM read_parquet(?)",
        [str(path)],
    ).fetchone()
    total = int(row[0])
    nulls = {}
    for index, name in enumerate(columns):
        present = int(row[index + 1])
        nulls[name] = total - present
    return {"rows": total, "nulls": nulls}


def fmt_num(value) -> str:
    if value is None:
        return "null"
    number = float(value)
    if number.is_integer() and abs(number) < 1e15:
        return str(int(number))
    return f"{number:.6g}"


def null_table(profile: dict, highlight: list[str]) -> str:
    lines = [
        "| Column | Nulls | Null rate |",
        "|---|---:|---:|",
    ]
    total = profile["rows"] or 1
    shown = set()
    ordered = list(highlight) + [name for name, count in profile["nulls"].items() if count]
    for name in ordered:
        if name in shown or name not in profile["nulls"]:
            continue
        shown.add(name)
        count = profile["nulls"][name]
        if name not in highlight and count == 0:
            continue
        rate = 100.0 * count / total
        lines.append(f"| `{name}` | {count:,} | {rate:.2f}% |")
    complete = sum(1 for count in profile["nulls"].values() if count == 0)
    lines.append("")
    lines.append(f"Columns with zero nulls: {complete} of {len(profile['nulls'])}.")
    return "\n".join(lines)


def load_truth(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if row.get("case_id") in TRUTH_CASES]


def parse_expected(case: dict) -> dict:
    text = " ".join(
        part
        for part in (case.get("expected") or "", case.get("notes") or "", case.get("rsid") or "")
        if part
    )
    loeuf = None
    decile = None
    transcript = None
    tokens = text.replace(";", " ").split()
    for index, token in enumerate(tokens):
        upper = token.upper()
        if upper == "LOEUF" and index + 1 < len(tokens):
            loeuf = float(tokens[index + 1])
        elif token.lower() == "decile" and index + 1 < len(tokens):
            decile = int(float(tokens[index + 1]))
        elif token.lower() == "transcript" and index + 1 < len(tokens):
            transcript = tokens[index + 1]
        elif token.startswith("ENST"):
            transcript = token
    if case["case_id"] == "mane_select" and not transcript:
        transcript = case.get("expected") or None
    return {
        "gene": case.get("gene") or "",
        "loeuf": loeuf,
        "decile": decile,
        "transcript": transcript,
        "expected": case.get("expected") or "",
        "notes": case.get("notes") or "",
        "requirement": case.get("requirement") or "",
        "source_locator": case.get("source_locator") or "",
    }


def evaluate_truth(con: duckdb.DuckDBPyConnection, cases: list[dict]) -> list[dict]:
    fetched = con.execute(
        """
        SELECT gene, gene_id, transcript, bare_ensg, bare_enst, loeuf, loeuf_decile,
               mane_select, canonical
        FROM read_parquet(?)
        """,
        [str(GENE_OUT)],
    ).fetchall()
    by_gene: dict[str, list[tuple]] = defaultdict(list)
    for record in fetched:
        by_gene[record[0]].append(record)
    results = []
    for case in cases:
        spec = parse_expected(case)
        rows = by_gene.get(spec["gene"], [])
        problems = []
        observed = "no gene row"
        if len(rows) != 1:
            problems.append(f"{len(rows)} gene rows for {spec['gene']}")
        if len(rows) >= 1:
            gene, gene_id, transcript, bare_ensg, bare_enst, loeuf, decile, mane, canon = rows[0]
            observed = (
                f"gene_id {gene_id}; transcript {transcript}; bare_ensg {bare_ensg}; "
                f"LOEUF {fmt_num(loeuf)}; decile {fmt_num(decile)}; "
                f"mane_select {mane}; canonical {canon}"
            )
            if spec["transcript"] and transcript != spec["transcript"] and bare_enst != spec["transcript"]:
                problems.append(f"transcript {transcript} != {spec['transcript']}")
            if spec["loeuf"] is not None and (loeuf is None or abs(float(loeuf) - spec["loeuf"]) > LOEUF_TOLERANCE):
                problems.append(f"LOEUF {fmt_num(loeuf)} != {spec['loeuf']}")
            if spec["decile"] is not None and (decile is None or int(float(decile)) != spec["decile"]):
                problems.append(f"decile {fmt_num(decile)} != {spec['decile']}")
        results.append(
            {
                "case_id": case["case_id"],
                "gene": spec["gene"],
                "requirement": spec["requirement"],
                "locator": spec["source_locator"],
                "expected": spec["expected"],
                "notes": spec["notes"],
                "observed": observed,
                "status": "fail" if problems else "pass",
                "problems": problems,
            }
        )
    return results


def a1bg_transcript_check(con: duckdb.DuckDBPyConnection) -> tuple[str, list[tuple]]:
    rows = con.execute(
        """
        SELECT transcript, canonical, mane_select, gene_id, "cds_length",
               "lof.oe_ci.upper", "lof.oe_ci.upper_bin_decile"
        FROM read_parquet(?)
        WHERE gene = 'A1BG'
        ORDER BY transcript
        """,
        [str(TRANSCRIPT_OUT)],
    ).fetchall()
    wanted = {"NM_130786.4", "ENST00000263100"}
    found = {row[0] for row in rows if row[1] is True}
    status = "pass" if wanted <= found else "fail"
    return status, rows


def join_note(con: duckdb.DuckDBPyConnection) -> str:
    if not INTERVALS.exists():
        return (
            "Lane A `warehouse/intervals/genes.parquet` is not on disk. "
            "No interval join was run, so there is no join-miss count. "
            "The integrator checks bare ENSG against intervals later."
        )
    fields = set(pq.read_schema(INTERVALS).names)
    if "ensg" not in fields:
        return (
            f"Lane A genes file exists but has no `ensg` column. Fields: {sorted(fields)}. "
            "Join was not forced onto another column."
        )
    biotype = "biotype" in fields
    biotype_sql = ", biotype" if biotype else ""
    query = f"""
        WITH genes AS (
            SELECT bare_ensg, gene, gene_id
            FROM read_parquet(?)
        ),
        intervals AS (
            SELECT ensg{biotype_sql}
            FROM read_parquet(?)
        )
        SELECT
            (SELECT count(*) FROM genes) AS gene_rows,
            (SELECT count(*) FROM genes WHERE bare_ensg IS NOT NULL) AS with_ensg,
            (SELECT count(*) FROM genes g WHERE g.bare_ensg IS NOT NULL
                AND NOT EXISTS (SELECT 1 FROM intervals i WHERE i.ensg = g.bare_ensg)) AS ensg_miss,
            (SELECT count(*) FROM genes WHERE bare_ensg IS NULL) AS no_ensg
    """
    gene_rows, with_ensg, ensg_miss, no_ensg = con.execute(
        query, [str(GENE_OUT), str(INTERVALS)]
    ).fetchone()
    lines = [
        f"Joined `gene_constraint.bare_ensg` to `{INTERVALS}` column `ensg`.",
        f"Gene rows: {gene_rows:,}. Rows with a bare ENSG: {with_ensg:,}. "
        f"Those missing from intervals: {ensg_miss:,}. "
        f"Gene rows with no bare ENSG (RefSeq gene id only): {no_ensg:,}.",
    ]
    if biotype:
        coding, coding_hit = con.execute(
            """
            WITH intervals AS (
                SELECT ensg, biotype FROM read_parquet(?)
            )
            SELECT
                (SELECT count(*) FROM intervals WHERE biotype = 'protein_coding'),
                (SELECT count(*) FROM intervals i
                    WHERE i.biotype = 'protein_coding'
                      AND EXISTS (
                          SELECT 1 FROM read_parquet(?) g WHERE g.bare_ensg = i.ensg
                      ))
            """,
            [str(INTERVALS), str(GENE_OUT)],
        ).fetchone()
        missed = coding - coding_hit
        lines.append(
            f"Interval protein-coding genes: {coding:,}. "
            f"Matched to a constraint gene row: {coding_hit:,}. "
            f"Protein-coding join misses: {missed:,}."
        )
    return " ".join(lines)


def gene_class_counts(genes: list[dict], rows: list[dict]) -> dict[str, int]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["gene_id"]].append(row)
    # Classes follow the selected row. MANE outranks canonical, so a gene with
    # any MANE candidate keeps a MANE row.
    mane = sum(1 for row in genes if row["mane_select"] is True)
    canonical_only = sum(
        1 for row in genes if row["mane_select"] is not True and row["canonical"] is True
    )
    neither = sum(
        1 for row in genes if row["mane_select"] is not True and row["canonical"] is not True
    )
    return {
        "mane": mane,
        "canonical_only": canonical_only,
        "neither": neither,
        "input_gene_ids": len(grouped),
    }


def dropped_details(rows: list[dict], dropped_ids: set[str]) -> list[dict]:
    per_id: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row["gene_id"] in dropped_ids:
            per_id[row["gene_id"]].append(row)
    details = []
    for gene_id, items in per_id.items():
        best = max(items, key=preference_key)
        details.append(best)
    details.sort(key=lambda row: (row["bare_ensg"] is None, row["gene"] or "", row["gene_id"] or ""))
    return details


def write_qc(
    header_notes: list[str],
    rows: list[dict],
    genes: list[dict],
    dropped_ids: set[str],
    transcript_profile: dict,
    gene_profile: dict,
    truth_results: list[dict],
    a1bg_status: str,
    a1bg_rows: list[tuple],
    join_text: str,
) -> None:
    classes = gene_class_counts(genes, rows)
    canonical_rows = sum(1 for row in rows if row["canonical"] is True)
    mane_rows = sum(1 for row in rows if row["mane_select"] is True)
    refseq_rows = sum(1 for row in rows if row["bare_ensg"] is None)
    ensembl_rows = len(rows) - refseq_rows
    dropped_ensg = sorted(gene_id for gene_id in dropped_ids if coords.bare_ensg(gene_id))
    dropped_other = sorted(gene_id for gene_id in dropped_ids if coords.bare_ensg(gene_id) is None)
    multi = defaultdict(set)
    for gene in genes:
        if gene["gene"] is not None and gene["bare_ensg"] is not None:
            multi[gene["gene"]].add(gene["bare_ensg"])
    multi_symbols = sorted(symbol for symbol, ids in multi.items() if len(ids) > 1)
    chroms = sorted({row["chromosome"] for row in rows if row["chromosome"]})
    level_counts = defaultdict(int)
    for row in rows:
        level_counts[fmt_num(row["level"])] += 1
    refseq_decile_null = sum(
        1 for row in rows if row["bare_ensg"] is None and row["lof.oe_ci.upper_bin_decile"] is None
    )
    refseq_chrom_null = sum(1 for row in rows if row["bare_ensg"] is None and row["chromosome"] is None)
    refseq_cds_null = sum(1 for row in rows if row["bare_ensg"] is None and row["cds_length"] is None)
    truth_fail = [item for item in truth_results if item["status"] != "pass"]
    a1bg_gene = next(row for row in genes if row["gene"] == "A1BG")
    a1bg_gene_pass = a1bg_gene["transcript"] == "ENST00000263100"
    required_fail = (not a1bg_gene_pass) or a1bg_status != "pass" or bool(truth_fail)
    ensg_drop_lines = []
    if dropped_ensg:
        bests = {row["gene_id"]: row for row in dropped_details(rows, set(dropped_ensg))}
        for gene_id in dropped_ensg:
            best = bests[gene_id]
            ensg_drop_lines.append(
                f"- `{best['gene']}` `{gene_id}` best transcript `{best['transcript']}` "
                f"mane_select {best['mane_select']} canonical {best['canonical']} "
                f"LOEUF {fmt_num(best['lof.oe_ci.upper'])} "
                f"decile {fmt_num(best['lof.oe_ci.upper_bin_decile'])}"
            )
    else:
        ensg_drop_lines.append("- None.")
    truth_lines = [
        "| Case | Gene | Expected | Observed | Result |",
        "|---|---|---|---|---|",
    ]
    for item in truth_results:
        expected = item["expected"] or item["notes"]
        problem = ""
        if item["problems"]:
            problem = " " + "; ".join(item["problems"])
        truth_lines.append(
            f"| `{item['case_id']}` | {item['gene']} | {expected} | {item['observed']}{problem} | "
            f"{item['status']} |"
        )
    a1bg_tx_lines = [
        "| Transcript | Canonical | MANE | gene_id | cds_length | LOEUF | Decile |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for transcript, canon, mane, gene_id, cds, loeuf, decile in a1bg_rows:
        a1bg_tx_lines.append(
            f"| `{transcript}` | {canon} | {mane} | `{gene_id}` | {fmt_num(cds)} | "
            f"{fmt_num(loeuf)} | {fmt_num(decile)} |"
        )
    highlight = [
        "loeuf",
        "loeuf_decile",
        "lof.oe_ci.upper",
        "lof.oe_ci.upper_bin_decile",
        "lof.pLI",
        "mis.z_score",
        "syn.z_score",
        "constraint_flags",
        "bare_ensg",
        "bare_enst",
        "chromosome",
        "cds_length",
        "mane_select",
        "canonical",
        "gene",
        "gene_id",
        "transcript",
    ]
    text = f"""# Constraint

Lane B. Source `{SOURCE}`.

Truth set `{TRUTH}` was read and was not rewritten. Constraint cases in that file are `mane_select`, `loeuf_low_1`, `loeuf_low_2`, `loeuf_low_3`, `loeuf_near_1_1`, `loeuf_near_1_2`, and `loeuf_near_1_3`.

Required checks: {"fail" if required_fail else "pass"}.

## Counts

| Metric | Count |
|---|---:|
| Transcript rows | {transcript_profile["rows"]:,} |
| Gene rows | {gene_profile["rows"]:,} |
| Source data rows read | {len(rows):,} |
| Raw distinct `gene_id` values | {classes["input_gene_ids"]:,} |
| `canonical=true` transcript rows | {canonical_rows:,} |
| `mane_select=true` transcript rows | {mane_rows:,} |
| Genes with a MANE transcript | {classes["mane"]:,} |
| Genes with only a canonical transcript | {classes["canonical_only"]:,} |
| Genes with neither MANE nor canonical | {classes["neither"]:,} |
| Genes with no MANE transcript | {classes["canonical_only"] + classes["neither"]:,} |
| Genes dropped by the preference rule | {len(dropped_ids):,} |
| Dropped Ensembl gene ids | {len(dropped_ensg):,} |
| Dropped NCBI gene ids | {len(dropped_other):,} |
| Transcript rows not kept on the gene table | {len(rows) - len(genes):,} |

`canonical=true` is not one row per gene. The transcript table has {canonical_rows:,} canonical rows and the gene table has {gene_profile["rows"]:,} rows. A1BG is stored twice as canonical in the transcript file, `NM_130786.4` and `ENST00000263100`.

DuckDB reread the Parquet files. Transcript rows {transcript_profile["rows"]:,}. Gene rows {gene_profile["rows"]:,}. Gene `gene_id` values are unique: {"pass" if gene_profile["rows"] == len(genes) else "fail"}.

## Preference

One gene row is the winning transcript, not every raw `gene_id`. RefSeq rows carry an NCBI gene id. Ensembl rows carry a bare ENSG. Those are the same gene when they share a symbol and that symbol has one Ensembl id, so they compete for the metric. The stored `gene_id` stays the Ensembl id, even when the winning transcript is RefSeq.

Order: `mane_select` true, then `canonical` true, then longest `cds_length`. A null `cds_length` does not beat a measured length. When a RefSeq transcript and an Ensembl transcript still tie, the Ensembl transcript is kept (`bare_enst` is set).

Symbols with more than one Ensembl gene id stay one row per Ensembl id. The symbol `NA` in the source is a missing symbol, stored null, and is not a group key. Gene rows with a null symbol: {sum(1 for row in genes if row["gene"] is None):,}.

Symbols with two Ensembl gene rows: {", ".join(f"`{symbol}`" for symbol in multi_symbols) or "none"}.

Genes dropped by the preference rule are gene ids present in the TSV that are not the `gene_id` of a gene-table row. {len(dropped_other):,} are NCBI ids that lost to an Ensembl transcript of the same symbol. {len(dropped_ensg):,} are Ensembl ids that lost because a RefSeq transcript outranked them:

{chr(10).join(ensg_drop_lines)}

## A1BG

Gene-table transcript is `{a1bg_gene["transcript"]}` (`gene_id` `{a1bg_gene["gene_id"]}`, LOEUF {fmt_num(a1bg_gene["loeuf"])}, decile {fmt_num(a1bg_gene["loeuf_decile"])}). Check that this is `ENST00000263100`: {"pass" if a1bg_gene_pass else "fail"}.

Transcript-file canonical copies include both `NM_130786.4` and `ENST00000263100`: {a1bg_status}.

{chr(10).join(a1bg_tx_lines)}

## Truth set

Cited from `{TRUTH}`. LOEUF must be within {LOEUF_TOLERANCE} of the truth-set value. Decile must match as an integer. The transcript must match the Ensembl id named in that row. "About 0.04", "about 0.041", and "about 0.048" are those truth-set values.

{chr(10).join(truth_lines)}

## Null rates

Transcript file:

{null_table(transcript_profile, [name for name in highlight if name not in {"loeuf", "loeuf_decile"}])}

Gene file:

{null_table(gene_profile, highlight)}

## Join

{join_text}

## File notes

{chr(10).join(f"- {note}" for note in header_notes)}
- LOEUF is `lof.oe_ci.upper`. On the gene file that float is also stored as `loeuf`. The decile is `lof.oe_ci.upper_bin_decile`, also stored as `loeuf_decile`.
- `constraint_flags` is kept as text. `[]` means no flag and is not null.
- `level` is a number 1, 2, or 3, or null. Counts: {", ".join(f"{key}={level_counts[key]:,}" for key in sorted(level_counts))}. It is not a gene-versus-transcript split. Every source row is a transcript row.
- Ensembl-id rows: {ensembl_rows:,}. Rows whose `gene_id` is not Ensembl (NCBI Gene id on the RefSeq transcript): {refseq_rows:,}.
- On those non-Ensembl rows, `chromosome` is null in {refseq_chrom_null:,}, `cds_length` is null in {refseq_cds_null:,}, and the LOEUF decile is null in {refseq_decile_null:,}. The decile used for the truth set is on the Ensembl transcript.
- `chromosome` is passed through `coords.canon_chrom`. Values present: {", ".join(f"`{chrom}`" for chrom in chroms) or "none"}. `chrX` and `chrY` are absent when they are not in that list.
- Numeric source fields are float64. `canonical` and `mane_select` are boolean. `NA` tokens are null.
- `bare_ensg` and `bare_enst` come from `scripts/coords.py` and are null when the id is not Ensembl.
"""
    QC_OUT.parent.mkdir(parents=True, exist_ok=True)
    QC_OUT.write_text(text, encoding="utf-8")
    if required_fail:
        raise SystemExit("constraint QC checks failed; see warehouse/qc/constraint.md")


def main() -> None:
    header, rows, header_notes = load_rows(SOURCE)
    genes, dropped_ids = select_genes(rows)
    transcript_columns = column_types(header, ["bare_ensg", "bare_enst"])
    gene_columns = column_types(header, ["bare_ensg", "bare_enst", "loeuf", "loeuf_decile"])
    write_parquet(TRANSCRIPT_OUT, rows, transcript_columns)
    write_parquet(GENE_OUT, genes, gene_columns)
    con = duckdb.connect()
    transcript_profile = parquet_profile(con, TRANSCRIPT_OUT, [name for name, _kind in transcript_columns])
    gene_profile = parquet_profile(con, GENE_OUT, [name for name, _kind in gene_columns])
    if transcript_profile["rows"] != len(rows):
        raise SystemExit(
            f"transcript parquet has {transcript_profile['rows']} rows, source has {len(rows)}"
        )
    if gene_profile["rows"] != len(genes):
        raise SystemExit(f"gene parquet has {gene_profile['rows']} rows, selection has {len(genes)}")
    truth_results = evaluate_truth(con, load_truth(TRUTH))
    a1bg_status, a1bg_rows = a1bg_transcript_check(con)
    try:
        join_text = join_note(con)
    except Exception as exc:
        join_text = (
            "Interval join was not completed. "
            f"`warehouse/intervals/genes.parquet` exists but could not be read cleanly: {exc}. "
            "Lane A may still be writing it. The integrator checks the join later."
        )
    write_qc(
        header_notes,
        rows,
        genes,
        dropped_ids,
        transcript_profile,
        gene_profile,
        truth_results,
        a1bg_status,
        a1bg_rows,
        join_text,
    )
    print(f"transcript_rows {transcript_profile['rows']}")
    print(f"gene_rows {gene_profile['rows']}")
    print(f"dropped_gene_ids {len(dropped_ids)}")
    print(f"qc {QC_OUT}")


if __name__ == "__main__":
    main()
