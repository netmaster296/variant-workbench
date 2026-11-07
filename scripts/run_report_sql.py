# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Run one report's DuckDB SQL against the warehouse under ROOT.

The report SQL files in warehouse/reports/sql/ are published exactly as the
build wrote them, so they hard-code the original project folder
(D:/variant-workbench). This runner reads a file, replaces that prefix with
the configured ROOT at run time, and runs the statements one by one in an
in-memory DuckDB connection. The .sql file on disk is never changed.

ROOT is the VW_ROOT environment variable when it is set, otherwise the folder
above scripts/.

Usage:
    python scripts/run_report_sql.py ancestry_audit
    python scripts/run_report_sql.py warehouse/reports/sql/splicing.sql --max-rows 50
    python scripts/run_report_sql.py splicing --dry-run      # print the rewritten SQL only

Notes:
- Some report SQL reads predictor columns (CADD, REVEL, SpliceAI, Pangolin,
  phyloP, SIFT, PolyPhen) that are not in the published Tier A tables. Those
  statements fail until you regenerate the columns locally (see README).
- The SQL sets DuckDB memory and thread limits and spills to
  ROOT/warehouse/tmp, which this runner creates if it is missing.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import duckdb

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
SQL_DIR = ROOT / "warehouse" / "reports" / "sql"
REPO_SQL_DIR = Path(__file__).resolve().parents[1] / "warehouse" / "reports" / "sql"
# The prefix the published report SQL was written with, either slash style.
ORIGINAL_PREFIX = re.compile(r"D:[/\\]variant-workbench", re.IGNORECASE)
RESULT_TYPES = {"SELECT", "EXPLAIN", "PRAGMA", "CALL"}


def find_sql(name: str) -> Path:
    candidate = Path(name)
    if candidate.is_file():
        return candidate
    stem = name[:-4] if name.lower().endswith(".sql") else name
    for folder in (SQL_DIR, REPO_SQL_DIR):
        path = folder / f"{stem}.sql"
        if path.is_file():
            return path
    raise SystemExit(f"no report SQL named {name!r} in {SQL_DIR} or {REPO_SQL_DIR}")


def rewrite(sql: str) -> tuple[str, int]:
    # Every hard-coded path sits inside a single-quoted SQL literal, so a quote
    # in ROOT is doubled to keep the literal valid.
    target = ROOT.as_posix().rstrip("/").replace("'", "''")
    return ORIGINAL_PREFIX.subn(lambda _m: target, sql)


def show(con: duckdb.DuckDBPyConnection, max_rows: int) -> None:
    columns = [d[0] for d in con.description or []]
    rows = con.fetchall()
    print("| " + " | ".join(columns) + " |")
    print("|" + "---|" * len(columns))
    shown = rows if max_rows <= 0 else rows[:max_rows]
    for row in shown:
        print("| " + " | ".join("" if v is None else str(v) for v in row) + " |")
    if len(shown) < len(rows):
        print(f"... {len(rows) - len(shown)} more rows (use --max-rows 0 to print all)")
    print(f"({len(rows)} rows)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("report", help="report name (e.g. ancestry_audit) or path to a .sql file")
    parser.add_argument("--max-rows", type=int, default=200, help="rows to print per result; 0 prints all")
    parser.add_argument("--dry-run", action="store_true", help="print the rewritten SQL and stop")
    args = parser.parse_args()

    path = find_sql(args.report)
    sql, count = rewrite(path.read_text(encoding="utf-8-sig"))
    print(f"-- {path.name}: ROOT = {ROOT}; replaced {count} path prefix(es)", file=sys.stderr)
    if args.dry_run:
        print(sql)
        return 0

    (ROOT / "warehouse" / "tmp").mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    statements = duckdb.extract_statements(sql)
    for index, statement in enumerate(statements, 1):
        kind = getattr(statement.type, "name", str(statement.type))
        try:
            con.execute(statement.query)
        except duckdb.Error as exc:
            print(f"\n-- statement {index}/{len(statements)} ({kind}) failed: {exc}", file=sys.stderr)
            return 1
        if kind in RESULT_TYPES and con.description:
            first = " ".join(statement.query.strip().split())[:100]
            print(f"\n-- result {index}/{len(statements)}: {first}")
            show(con, args.max_rows)
    con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
