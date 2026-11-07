# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Build the GENCODE gene crosswalk.

One row per versioned GENCODE gene. The gene join is bare ENSG.
Symbol attaches ClinVar GeneID and HGNC_ID only. Symbols collide, so a
symbol match never creates or drops a gene row.

MANE Ensembl id: GENCODE transcript tag MANE_Select when that tag exists,
otherwise a constraint transcript with mane_select and an ENST id.
MANE RefSeq: a constraint transcript with mane_select whose id starts with
NM_ or NR_, paired to the same bare ENSG inside the constraint file.
RefSeq rows have no ENSG. They are not symbol-joined to GENCODE.
A constraint gene row with a null bare ENSG stays unmatched on id.
"""

from __future__ import annotations

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

GENES = ROOT / "warehouse" / "intervals" / "genes.parquet"
TRANSCRIPTS = ROOT / "warehouse" / "intervals" / "transcripts.parquet"
GENE_CONSTRAINT = ROOT / "warehouse" / "constraint" / "gene_constraint.parquet"
TRANSCRIPT_CONSTRAINT = ROOT / "warehouse" / "constraint" / "transcript_constraint.parquet"
ALLELES = ROOT / "warehouse" / "clinvar" / "alleles.parquet"
GENE_OUT = ROOT / "warehouse" / "crosswalk" / "genes.parquet"
UNMATCHED_OUT = ROOT / "warehouse" / "crosswalk" / "unmatched.parquet"
QC_OUT = ROOT / "warehouse" / "qc" / "crosswalk.md"

EXPECTED_GENES = 78691
HGNC_MISSING = {"", "-", ".", "na", "NA", "n/a", "N/A"}

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
        pa.field("ncbi_gene_id", pa.string()),
        pa.field("ncbi_gene_id_primary", pa.string()),
        pa.field("ncbi_gene_id_n", pa.int64(), nullable=False),
        pa.field("ncbi_gene_id_rule", pa.string()),
        pa.field("ncbi_clinvar_symbols", pa.string()),
        pa.field("hgnc_id", pa.string()),
        pa.field("hgnc_id_gtf", pa.string()),
        pa.field("hgnc_id_clinvar", pa.string()),
        pa.field("hgnc_disagree", pa.bool_(), nullable=False),
        pa.field("mane_enst", pa.string()),
        pa.field("mane_enst_version", pa.string()),
        pa.field("mane_enst_source", pa.string()),
        pa.field("mane_enst_gencode", pa.string()),
        pa.field("mane_enst_constraint", pa.string()),
        pa.field("mane_enst_disagree", pa.bool_(), nullable=False),
        pa.field("mane_refseq", pa.string()),
        pa.field("constraint_ensg_matched", pa.bool_(), nullable=False),
        pa.field("symbol_n_gencode", pa.int64(), nullable=False),
    ]
)

UNMATCHED_SCHEMA = pa.schema(
    [
        pa.field("source", pa.string(), nullable=False),
        pa.field("key", pa.string(), nullable=False),
        pa.field("reason", pa.string(), nullable=False),
        pa.field("n", pa.int64(), nullable=False),
        pa.field("n_single", pa.int64()),
        pa.field("n_multi", pa.int64()),
        pa.field("symbol", pa.string()),
        pa.field("bare_ensg", pa.string()),
        pa.field("gene_id", pa.string()),
        pa.field("transcript", pa.string()),
        pa.field("mane_refseq", pa.string()),
        pa.field("mane_select", pa.bool_()),
        pa.field("note", pa.string(), nullable=False),
    ]
)


def connect() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("SET memory_limit='8GB'")
    # SPECIAL: bare_ensg/bare_enst return None for NCBI and RefSeq ids.
    con.create_function(
        "coords_bare_ensg",
        lambda value: coords.bare_ensg(value),
        ["VARCHAR"],
        "VARCHAR",
        null_handling="special",
    )
    con.create_function(
        "coords_bare_enst",
        lambda value: coords.bare_enst(value),
        ["VARCHAR"],
        "VARCHAR",
        null_handling="special",
    )
    return con


def records(con: duckdb.DuckDBPyConnection, sql: str, params: list | None = None) -> list[dict]:
    cursor = con.execute(sql, params or [])
    names = [col[0] for col in cursor.description]
    return [dict(zip(names, row)) for row in cursor.fetchall()]


def one(con: duckdb.DuckDBPyConnection, sql: str, params: list | None = None):
    row = con.execute(sql, params or []).fetchone()
    if row is None:
        raise RuntimeError("query returned no row")
    return row[0]


def pq_path(path: Path) -> str:
    return path.as_posix()


def norm_hgnc(value: str | None) -> str | None:
    """Canonical HGNC id. Numeric ids and leading zeros compare equal."""
    if value is None:
        return None
    text = str(value).strip()
    if text in HGNC_MISSING:
        return None
    prefix, sep, rest = text.partition(":")
    if sep and prefix.upper() == "HGNC" and rest.isdigit():
        return "HGNC:" + str(int(rest))
    if text.isdigit():
        return "HGNC:" + str(int(text))
    return text


def delimited(values, key=None) -> str | None:
    seen: set[str] = set()
    items: list[str] = []
    for value in values:
        if value is None:
            continue
        text = str(value)
        if text in seen:
            continue
        seen.add(text)
        items.append(text)
    if not items:
        return None
    items.sort(key=key)
    return ";".join(items)


def ncbi_sort_key(value: str):
    if value.isdigit():
        return (0, int(value))
    return (1, value)


def choose_ncbi(id_info: dict[str, dict]) -> tuple[str, str]:
    """Most allele rows wins. A tie keeps the smaller numeric GeneID."""
    ranking = sorted(
        id_info,
        key=lambda ncbi: (
            -id_info[ncbi]["n"],
            int(ncbi) if ncbi.isdigit() else 10**30,
            ncbi,
        ),
    )
    primary = ranking[0]
    if len(ranking) == 1:
        return primary, "sole"
    top_n = id_info[primary]["n"]
    tied = [ncbi for ncbi in ranking if id_info[ncbi]["n"] == top_n]
    if len(tied) > 1:
        return primary, "smallest_id"
    return primary, "most_alleles"


def check_identifiers(con: duckdb.DuckDBPyConnection) -> list[str]:
    """Bare ids used here are coords.bare_ensg and coords.bare_enst."""
    assert coords.bare_ensg("ENSG00000121410.15") == "ENSG00000121410"
    assert coords.bare_ensg("1") is None
    assert coords.bare_enst("ENST00000263100.8") == "ENST00000263100"
    assert coords.bare_enst("NM_130786.4") is None
    notes = []
    gene_bad = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE coords_bare_ensg(ensg_version) IS DISTINCT FROM ensg
        """,
        [pq_path(GENES)],
    )
    tx_bad = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE list_contains(tags, 'MANE_Select')
          AND coords_bare_enst(enst_version) IS DISTINCT FROM enst
        """,
        [pq_path(TRANSCRIPTS)],
    )
    constraint_gene_bad = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE coords_bare_ensg(gene_id) IS DISTINCT FROM bare_ensg
        """,
        [pq_path(GENE_CONSTRAINT)],
    )
    constraint_tx_bad = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE coords_bare_ensg(gene_id) IS DISTINCT FROM bare_ensg
           OR coords_bare_enst(transcript) IS DISTINCT FROM bare_enst
        """,
        [pq_path(TRANSCRIPT_CONSTRAINT)],
    )
    notes.append(
        "Identifier check against coords.bare_ensg and coords.bare_enst: "
        f"GENCODE gene rows differing {gene_bad:,}; "
        f"MANE transcripts differing {tx_bad:,}; "
        f"gene-constraint rows differing {constraint_gene_bad:,}; "
        f"transcript-constraint rows differing {constraint_tx_bad:,}."
    )
    if any((gene_bad, tx_bad, constraint_gene_bad, constraint_tx_bad)):
        raise SystemExit(notes[-1])
    return notes


def load_genes(con: duckdb.DuckDBPyConnection) -> list[dict]:
    return records(
        con,
        """
        SELECT ensg, ensg_version, symbol, chrom, start, "end" AS end_pos,
               strand, biotype, hgnc_id, mane_select
        FROM read_parquet(?)
        ORDER BY ensg
        """,
        [pq_path(GENES)],
    )


def load_gencode_mane(con: duckdb.DuckDBPyConnection) -> tuple[dict[str, set[str]], dict[str, str]]:
    rows = records(
        con,
        """
        SELECT ensg, enst, enst_version
        FROM read_parquet(?)
        WHERE list_contains(tags, 'MANE_Select')
        """,
        [pq_path(TRANSCRIPTS)],
    )
    by_ensg: dict[str, set[str]] = defaultdict(set)
    version_of: dict[str, str] = {}
    for row in rows:
        bare = coords.bare_enst(row["enst_version"])
        if bare is None:
            raise SystemExit(f"MANE transcript has no ENST id: {row}")
        by_ensg[row["ensg"]].add(bare)
        previous = version_of.get(bare)
        if previous is not None and previous != row["enst_version"]:
            raise SystemExit(f"bare ENST {bare} has versions {previous} and {row['enst_version']}")
        version_of[bare] = row["enst_version"]
    return by_ensg, version_of


def load_constraint(con: duckdb.DuckDBPyConnection) -> dict:
    matched = {
        row["bare_ensg"]
        for row in records(
            con,
            """
            SELECT bare_ensg
            FROM read_parquet(?)
            WHERE bare_ensg IS NOT NULL
            """,
            [pq_path(GENE_CONSTRAINT)],
        )
    }
    tx_counts = {
        row["gene_id"]: int(row["n"])
        for row in records(
            con,
            """
            SELECT gene_id, count(*)::BIGINT AS n
            FROM read_parquet(?)
            GROUP BY 1
            """,
            [pq_path(TRANSCRIPT_CONSTRAINT)],
        )
    }
    unmatched_rows = records(
        con,
        """
        SELECT gene_id, gene AS symbol, bare_ensg, transcript, mane_select
        FROM read_parquet(?)
        WHERE bare_ensg IS NULL
           OR bare_ensg NOT IN (SELECT ensg FROM read_parquet(?))
        ORDER BY gene_id
        """,
        [pq_path(GENE_CONSTRAINT), pq_path(GENES)],
    )
    mane_enst_rows = records(
        con,
        """
        SELECT coords_bare_ensg(gene_id) AS bare_ensg,
               coords_bare_enst(transcript) AS bare_enst,
               gene AS symbol,
               transcript
        FROM read_parquet(?)
        WHERE mane_select
          AND coords_bare_enst(transcript) IS NOT NULL
        """,
        [pq_path(TRANSCRIPT_CONSTRAINT)],
    )
    refseq_rows = records(
        con,
        """
        SELECT gene AS symbol, transcript, gene_id,
               coords_bare_ensg(gene_id) AS bare_ensg
        FROM read_parquet(?)
        WHERE mane_select
          AND (starts_with(transcript, 'NM_') OR starts_with(transcript, 'NR_'))
        """,
        [pq_path(TRANSCRIPT_CONSTRAINT)],
    )
    symbol_ensg = records(
        con,
        """
        SELECT gene AS symbol, coords_bare_ensg(gene_id) AS bare_ensg
        FROM read_parquet(?)
        WHERE gene IS NOT NULL
          AND coords_bare_ensg(gene_id) IS NOT NULL
        GROUP BY 1, 2
        """,
        [pq_path(TRANSCRIPT_CONSTRAINT)],
    )
    return {
        "matched": matched,
        "tx_counts": tx_counts,
        "unmatched_rows": unmatched_rows,
        "mane_enst_rows": mane_enst_rows,
        "refseq_rows": refseq_rows,
        "symbol_ensg": symbol_ensg,
    }


def pair_refseq(constraint: dict) -> tuple[dict[str, set[str]], list[dict], dict[str, set[str]]]:
    """Pair NM_/NR_ MANE transcripts to one bare ENSG inside the constraint file.

    A symbol with one MANE ENST gives that transcript's ENSG.
    A symbol with no MANE ENST and exactly one ENSG on any transcript uses that ENSG.
    Anything else is unattached. It is not matched to GENCODE by symbol.
    """
    mane_enst_by_ensg: dict[str, set[str]] = defaultdict(set)
    mane_ensg_by_symbol: dict[str, set[str]] = defaultdict(set)
    for row in constraint["mane_enst_rows"]:
        if row["bare_ensg"] is None or row["bare_enst"] is None:
            continue
        mane_enst_by_ensg[row["bare_ensg"]].add(row["bare_enst"])
        if row["symbol"] is not None:
            mane_ensg_by_symbol[row["symbol"]].add(row["bare_ensg"])

    ensgs_by_symbol: dict[str, set[str]] = defaultdict(set)
    for row in constraint["symbol_ensg"]:
        ensgs_by_symbol[row["symbol"]].add(row["bare_ensg"])

    refseq_by_ensg: dict[str, set[str]] = defaultdict(set)
    unattached: list[dict] = []
    for row in constraint["refseq_rows"]:
        symbol = row["symbol"]
        target = None
        if symbol is not None and len(mane_ensg_by_symbol.get(symbol, ())) == 1:
            target = next(iter(mane_ensg_by_symbol[symbol]))
        elif (
            symbol is not None
            and symbol not in mane_ensg_by_symbol
            and len(ensgs_by_symbol.get(symbol, ())) == 1
        ):
            target = next(iter(ensgs_by_symbol[symbol]))
        if target is None:
            unattached.append(row)
            continue
        refseq_by_ensg[target].add(row["transcript"])
    return mane_enst_by_ensg, refseq_by_ensg, unattached


def load_clinvar(con: duckdb.DuckDBPyConnection) -> dict:
    singles = records(
        con,
        """
        SELECT GeneSymbol AS symbol,
               GeneID AS ncbi,
               CASE
                 WHEN HGNC_ID IS NULL OR HGNC_ID IN ('-', '.', '', 'na', 'NA', 'n/a', 'N/A')
                   THEN NULL
                 ELSE HGNC_ID
               END AS hgnc,
               count(*)::BIGINT AS n_alleles
        FROM read_parquet(?)
        WHERE regexp_full_match(GeneID, '^[0-9]+$')
          AND GeneSymbol IS NOT NULL
          AND GeneSymbol NOT IN ('-', '.', '')
          AND position(';' IN GeneSymbol) = 0
          AND position(' ' IN GeneSymbol) = 0
          AND position(',' IN GeneSymbol) = 0
        GROUP BY 1, 2, 3
        """,
        [pq_path(ALLELES)],
    )
    multi = records(
        con,
        """
        SELECT trim(tok) AS symbol, count(*)::BIGINT AS n
        FROM read_parquet(?),
             unnest(string_split(GeneSymbol, ';')) AS u(tok)
        WHERE GeneID = '-1'
          AND position(';' IN GeneSymbol) > 0
        GROUP BY 1
        """,
        [pq_path(ALLELES)],
    )
    other = records(
        con,
        """
        SELECT
          CASE WHEN GeneSymbol = '-' THEN 'missing' ELSE 'prose' END AS kind,
          GeneSymbol AS key,
          count(*)::BIGINT AS n
        FROM read_parquet(?)
        WHERE GeneID = '-1'
          AND position(';' IN GeneSymbol) = 0
        GROUP BY 1, 2
        """,
        [pq_path(ALLELES)],
    )
    rejected = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE regexp_full_match(GeneID, '^[0-9]+$')
          AND (
            GeneSymbol IS NULL
            OR GeneSymbol IN ('-', '.', '')
            OR position(';' IN GeneSymbol) > 0
            OR position(' ' IN GeneSymbol) > 0
            OR position(',' IN GeneSymbol) > 0
          )
        """,
        [pq_path(ALLELES)],
    )
    bad_tokens = one(
        con,
        """
        SELECT count(*)
        FROM read_parquet(?)
        WHERE GeneID = '-1'
          AND position(';' IN GeneSymbol) > 0
          AND EXISTS (
            SELECT 1
            FROM unnest(string_split(GeneSymbol, ';')) AS u(tok)
            WHERE trim(tok) = ''
               OR position(' ' IN tok) > 0
               OR position(',' IN tok) > 0
          )
        """,
        [pq_path(ALLELES)],
    )
    return {
        "singles": singles,
        "multi": multi,
        "other": other,
        "rejected_single_rows": int(rejected),
        "bad_multigene_rows": int(bad_tokens),
    }


def clinvar_by_symbol(singles: list[dict]) -> tuple[dict[str, dict], dict[str, set[str]]]:
    by_symbol: dict[str, dict] = {}
    for row in singles:
        slot = by_symbol.get(row["symbol"])
        if slot is None:
            slot = {"ids": {}}
            by_symbol[row["symbol"]] = slot
        info = slot["ids"].get(row["ncbi"])
        if info is None:
            info = {"n": 0, "hgncs": set()}
            slot["ids"][row["ncbi"]] = info
        info["n"] += int(row["n_alleles"])
        hgnc = norm_hgnc(row["hgnc"])
        if hgnc is not None:
            info["hgncs"].add(hgnc)
    symbols_by_ncbi: dict[str, set[str]] = defaultdict(set)
    for symbol, slot in by_symbol.items():
        for ncbi in slot["ids"]:
            symbols_by_ncbi[ncbi].add(symbol)
    return by_symbol, symbols_by_ncbi


def build_gene_rows(
    genes: list[dict],
    clinvar: dict[str, dict],
    symbols_by_ncbi: dict[str, set[str]],
    gencode_mane: dict[str, set[str]],
    mane_versions: dict[str, str],
    mane_enst_by_ensg: dict[str, set[str]],
    refseq_by_ensg: dict[str, set[str]],
    matched_ensg: set[str],
) -> list[dict]:
    symbol_n: dict[str, int] = defaultdict(int)
    for gene in genes:
        if gene["symbol"] is not None:
            symbol_n[gene["symbol"]] += 1
    rows = []
    for gene in genes:
        symbol = gene["symbol"]
        slot = clinvar.get(symbol) if symbol is not None else None
        if slot is None:
            ncbi_all = None
            primary = None
            n_ids = 0
            rule = None
            clinvar_symbols = None
            hgnc_clinvar_values: set[str] = set()
        else:
            primary, rule = choose_ncbi(slot["ids"])
            ncbi_all = delimited(slot["ids"], key=ncbi_sort_key)
            n_ids = len(slot["ids"])
            clinvar_symbols = delimited(symbols_by_ncbi.get(primary, ()))
            hgnc_clinvar_values = set()
            for info in slot["ids"].values():
                hgnc_clinvar_values.update(info["hgncs"])
        hgnc_gtf = norm_hgnc(gene["hgnc_id"])
        hgnc_clinvar = delimited(hgnc_clinvar_values)
        if hgnc_gtf is not None:
            hgnc = hgnc_gtf
        elif len(hgnc_clinvar_values) == 1:
            hgnc = next(iter(hgnc_clinvar_values))
        else:
            # Several ClinVar HGNC ids and no GTF id: keep the list, do not pick one.
            hgnc = None
        hgnc_disagree = bool(
            hgnc_gtf is not None
            and hgnc_clinvar_values
            and hgnc_clinvar_values != {hgnc_gtf}
        )
        gencode_ids = gencode_mane.get(gene["ensg"], set())
        constraint_ids = mane_enst_by_ensg.get(gene["ensg"], set())
        if gencode_ids:
            chosen = min(gencode_ids)
            source = "gencode"
            version = mane_versions.get(chosen)
        elif constraint_ids:
            chosen = min(constraint_ids)
            source = "constraint"
            version = None
        else:
            chosen = None
            source = None
            version = None
        refseqs = refseq_by_ensg.get(gene["ensg"], set())
        rows.append(
            {
                "ensg": gene["ensg"],
                "ensg_version": gene["ensg_version"],
                "symbol": symbol,
                "chrom": gene["chrom"],
                "start": int(gene["start"]),
                "end": int(gene["end_pos"]),
                "strand": gene["strand"],
                "biotype": gene["biotype"],
                "ncbi_gene_id": ncbi_all,
                "ncbi_gene_id_primary": primary,
                "ncbi_gene_id_n": n_ids,
                "ncbi_gene_id_rule": rule,
                "ncbi_clinvar_symbols": clinvar_symbols,
                "hgnc_id": hgnc,
                "hgnc_id_gtf": hgnc_gtf,
                "hgnc_id_clinvar": hgnc_clinvar,
                "hgnc_disagree": hgnc_disagree,
                "mane_enst": chosen,
                "mane_enst_version": version,
                "mane_enst_source": source,
                "mane_enst_gencode": delimited(gencode_ids),
                "mane_enst_constraint": delimited(constraint_ids),
                "mane_enst_disagree": bool(
                    gencode_ids and constraint_ids and gencode_ids != constraint_ids
                ),
                "mane_refseq": delimited(refseqs),
                "constraint_ensg_matched": gene["ensg"] in matched_ensg,
                "symbol_n_gencode": symbol_n[symbol] if symbol is not None else 1,
            }
        )
    rows.sort(key=lambda row: row["ensg"])
    return rows


def build_unmatched(
    genes: list[dict],
    clinvar_raw: dict,
    constraint: dict,
    refseq_by_ensg: dict[str, set[str]],
) -> list[dict]:
    gencode_symbols = {gene["symbol"] for gene in genes if gene["symbol"] is not None}
    single_n: dict[str, int] = defaultdict(int)
    for row in clinvar_raw["singles"]:
        single_n[row["symbol"]] += int(row["n_alleles"])
    multi_n: dict[str, int] = {}
    for row in clinvar_raw["multi"]:
        symbol = row["symbol"]
        if symbol is None or symbol == "":
            continue
        multi_n[symbol] = multi_n.get(symbol, 0) + int(row["n"])

    unmatched = []
    symbols = set(single_n) | set(multi_n)
    for symbol in symbols:
        if symbol in gencode_symbols:
            continue
        n_single = int(single_n.get(symbol, 0))
        n_multi = int(multi_n.get(symbol, 0))
        unmatched.append(
            {
                "source": "clinvar",
                "key": symbol,
                "reason": "no_gencode_symbol",
                "n": n_single + n_multi,
                "n_single": n_single,
                "n_multi": n_multi,
                "symbol": symbol,
                "bare_ensg": None,
                "gene_id": None,
                "transcript": None,
                "mane_refseq": None,
                "mane_select": None,
                "note": (
                    "ClinVar symbol has no GENCODE symbol. "
                    "n sums single-gene allele rows and multigene rows that list this symbol. "
                    "Multigene rows have GeneID -1 and do not assign an NCBI GeneID."
                ),
            }
        )
    for row in clinvar_raw["other"]:
        if row["kind"] == "missing":
            reason = "missing_symbol"
            note = "ClinVar GeneSymbol is '-'. Missing symbol, not a failed symbol join."
            symbol = None
        else:
            reason = "not_a_gene_symbol"
            note = "ClinVar GeneID -1 label is not a gene symbol. Not joined."
            symbol = None
        unmatched.append(
            {
                "source": "clinvar",
                "key": row["key"],
                "reason": reason,
                "n": int(row["n"]),
                "n_single": None,
                "n_multi": None,
                "symbol": symbol,
                "bare_ensg": None,
                "gene_id": None,
                "transcript": None,
                "mane_refseq": None,
                "mane_select": None,
                "note": note,
            }
        )
    for row in constraint["unmatched_rows"]:
        bare = row["bare_ensg"]
        if bare is None:
            reason = "bare_ensg_null"
            note = "Constraint gene row has no bare ENSG. Unmatched on id, not a failed symbol join."
        else:
            reason = "bare_ensg_not_in_gencode"
            note = "Constraint bare ENSG is not a GENCODE gene."
        n_tx = int(constraint["tx_counts"].get(row["gene_id"], 0))
        transcript = row["transcript"]
        if bare is not None:
            refseq = delimited(refseq_by_ensg.get(bare, ()))
        elif row["mane_select"] and transcript and (
            transcript.startswith("NM_") or transcript.startswith("NR_")
        ):
            refseq = transcript
        else:
            refseq = None
        unmatched.append(
            {
                "source": "constraint",
                "key": row["gene_id"],
                "reason": reason,
                "n": n_tx,
                "n_single": None,
                "n_multi": None,
                "symbol": row["symbol"],
                "bare_ensg": bare,
                "gene_id": row["gene_id"],
                "transcript": transcript,
                "mane_refseq": refseq,
                "mane_select": bool(row["mane_select"]),
                "note": note,
            }
        )
    unmatched.sort(key=lambda row: (row["source"], row["reason"], -row["n"], row["key"]))
    return unmatched


def table_from_rows(rows: list[dict], schema: pa.Schema) -> pa.Table:
    arrays = []
    for field in schema:
        arrays.append(pa.array([row[field.name] for row in rows], type=field.type))
    return pa.Table.from_arrays(arrays, schema=schema)


def write_outputs(genes: list[dict], unmatched: list[dict]) -> None:
    GENE_OUT.parent.mkdir(parents=True, exist_ok=True)
    QC_OUT.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table_from_rows(genes, GENE_SCHEMA), GENE_OUT, compression="zstd")
    pq.write_table(table_from_rows(unmatched, UNMATCHED_SCHEMA), UNMATCHED_OUT, compression="zstd")


def md_table(headers: list[str], rows: list[list]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        cells = []
        for value in row:
            if isinstance(value, int):
                cells.append(f"{value:,}")
            elif value is None:
                cells.append("")
            else:
                cells.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def pct(part: int, whole: int) -> str:
    if whole == 0:
        return "n/a"
    return f"{part / whole:.2%}"


def write_qc(con: duckdb.DuckDBPyConnection, id_notes: list[str], build_notes: list[str]) -> None:
    genes = pq_path(GENE_OUT)
    unmatched = pq_path(UNMATCHED_OUT)
    intervals = pq_path(GENES)
    transcripts = pq_path(TRANSCRIPTS)
    gene_constraint = pq_path(GENE_CONSTRAINT)
    transcript_constraint = pq_path(TRANSCRIPT_CONSTRAINT)

    row_n = int(one(con, "SELECT count(*) FROM read_parquet(?)", [genes]))
    distinct_version = int(one(con, "SELECT count(DISTINCT ensg_version) FROM read_parquet(?)", [genes]))
    distinct_ensg = int(one(con, "SELECT count(DISTINCT ensg) FROM read_parquet(?)", [genes]))
    interval_n = int(one(con, "SELECT count(*) FROM read_parquet(?)", [intervals]))
    pc = records(
        con,
        """
        SELECT
          count(*) FILTER (WHERE biotype = 'protein_coding') AS n,
          count(*) FILTER (
            WHERE biotype = 'protein_coding' AND constraint_ensg_matched
          ) AS matched,
          count(*) FILTER (
            WHERE biotype = 'protein_coding' AND NOT constraint_ensg_matched
          ) AS unmatched
        FROM read_parquet(?)
        """,
        [genes],
    )[0]
    matched_all = int(
        one(con, "SELECT count(*) FROM read_parquet(?) WHERE constraint_ensg_matched", [genes])
    )
    dup = records(
        con,
        """
        WITH d AS (
          SELECT symbol, count(*) AS n
          FROM read_parquet(?)
          GROUP BY symbol
          HAVING count(*) > 1
        )
        SELECT count(*) AS symbols,
               coalesce(sum(n), 0) AS genes,
               count(*) FILTER (
                 WHERE symbol IN (
                   SELECT symbol FROM read_parquet(?) WHERE biotype = 'protein_coding'
                 )
               ) AS symbols_touching_protein
        FROM d
        """,
        [genes, genes],
    )[0]
    pc_dup_symbols = int(
        one(
            con,
            """
            WITH pc AS (
              SELECT symbol
              FROM read_parquet(?)
              WHERE biotype = 'protein_coding' AND symbol IS NOT NULL
              GROUP BY symbol
              HAVING count(*) > 1
            )
            SELECT count(*) FROM pc
            """,
            [genes],
        )
    )
    pc_dup_genes = int(
        one(
            con,
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE biotype = 'protein_coding'
              AND symbol IN (
                SELECT symbol
                FROM read_parquet(?)
                WHERE biotype = 'protein_coding'
                GROUP BY symbol
                HAVING count(*) > 1
              )
            """,
            [genes, genes],
        )
    )
    multiplicity = records(
        con,
        """
        SELECT n_genes, count(*) AS n_symbols
        FROM (
          SELECT symbol, count(*) AS n_genes
          FROM read_parquet(?)
          GROUP BY symbol
          HAVING count(*) > 1
        )
        GROUP BY 1
        ORDER BY 1
        """,
        [genes],
    )
    largest_dups = records(
        con,
        """
        SELECT symbol, count(*) AS n_genes,
               count(*) FILTER (WHERE biotype = 'protein_coding') AS protein_coding,
               count(DISTINCT biotype) AS n_biotype,
               string_agg(DISTINCT biotype, ', ' ORDER BY biotype) AS biotypes
        FROM read_parquet(?)
        GROUP BY symbol
        HAVING count(*) > 1
        ORDER BY n_genes DESC, symbol
        LIMIT 15
        """,
        [genes],
    )
    pc_dup_list = records(
        con,
        """
        SELECT symbol, count(*) AS n_genes,
               string_agg(ensg, ', ' ORDER BY ensg) AS ensgs
        FROM read_parquet(?)
        WHERE biotype = 'protein_coding'
        GROUP BY symbol
        HAVING count(*) > 1
        ORDER BY n_genes DESC, symbol
        """,
        [genes],
    )
    clinvar_symbols = records(
        con,
        """
        SELECT reason, count(*) AS symbols, sum(n)::BIGINT AS alleles,
               sum(n_single)::BIGINT AS single_alleles,
               sum(n_multi)::BIGINT AS multi_alleles
        FROM read_parquet(?)
        WHERE source = 'clinvar'
        GROUP BY reason
        ORDER BY reason
        """,
        [unmatched],
    )
    constraint_unmatched = records(
        con,
        """
        SELECT reason, count(*) AS gene_rows, sum(n)::BIGINT AS transcript_rows,
               count(*) FILTER (WHERE mane_select) AS mane_gene_rows
        FROM read_parquet(?)
        WHERE source = 'constraint'
        GROUP BY reason
        ORDER BY reason
        """,
        [unmatched],
    )
    hgnc = records(
        con,
        """
        SELECT
          count(*) FILTER (WHERE hgnc_disagree) AS disagree,
          count(*) FILTER (WHERE hgnc_id_gtf IS NOT NULL AND hgnc_id_clinvar IS NOT NULL) AS both,
          count(*) FILTER (WHERE hgnc_id_gtf IS NOT NULL AND hgnc_id_clinvar IS NULL) AS gtf_only,
          count(*) FILTER (WHERE hgnc_id_gtf IS NULL AND hgnc_id_clinvar IS NOT NULL) AS clinvar_only,
          count(*) FILTER (WHERE hgnc_id IS NULL) AS resolved_null,
          count(*) FILTER (WHERE hgnc_id_clinvar LIKE '%;%') AS multi_clinvar_hgnc
        FROM read_parquet(?)
        """,
        [genes],
    )[0]
    hgnc_examples = records(
        con,
        """
        SELECT symbol, ensg, hgnc_id_gtf, hgnc_id_clinvar, hgnc_id
        FROM read_parquet(?)
        WHERE hgnc_disagree
        ORDER BY symbol, ensg
        LIMIT 20
        """,
        [genes],
    )
    ncbi = records(
        con,
        """
        SELECT
          count(*) FILTER (WHERE ncbi_gene_id IS NOT NULL) AS genes_with_ncbi,
          count(*) FILTER (WHERE ncbi_gene_id IS NULL) AS genes_without_ncbi,
          count(*) FILTER (WHERE ncbi_gene_id_n > 1) AS genes_multi_ncbi,
          count(*) FILTER (WHERE ncbi_gene_id_rule = 'sole') AS rule_sole,
          count(*) FILTER (WHERE ncbi_gene_id_rule = 'most_alleles') AS rule_most,
          count(*) FILTER (WHERE ncbi_gene_id_rule = 'smallest_id') AS rule_smallest,
          count(*) FILTER (
            WHERE ncbi_clinvar_symbols IS NOT NULL AND position(';' IN ncbi_clinvar_symbols) > 0
          ) AS genes_shared_ncbi_symbol
        FROM read_parquet(?)
        """,
        [genes],
    )[0]
    shared_ncbi = records(
        con,
        """
        SELECT ncbi_gene_id_primary AS ncbi, ncbi_clinvar_symbols AS symbols, count(*) AS gencode_genes
        FROM read_parquet(?)
        WHERE ncbi_clinvar_symbols IS NOT NULL
          AND position(';' IN ncbi_clinvar_symbols) > 0
        GROUP BY 1, 2
        ORDER BY gencode_genes DESC, ncbi
        LIMIT 20
        """,
        [genes],
    )
    shared_ncbi_n = int(
        one(
            con,
            """
            SELECT count(*) FROM (
              SELECT ncbi_gene_id_primary, ncbi_clinvar_symbols
              FROM read_parquet(?)
              WHERE ncbi_clinvar_symbols IS NOT NULL
                AND position(';' IN ncbi_clinvar_symbols) > 0
              GROUP BY 1, 2
            )
            """,
            [genes],
        )
    )
    mane = records(
        con,
        """
        WITH gtf AS (
          SELECT ensg FROM read_parquet(?) WHERE mane_select
        ),
        con_m AS (
          SELECT gene_id, gene AS symbol, bare_ensg, transcript
          FROM read_parquet(?)
          WHERE mane_select
        )
        SELECT
          (SELECT count(*) FROM gtf) AS gtf_genes,
          (SELECT count(*) FROM con_m) AS constraint_genes,
          (SELECT count(*) FROM con_m c JOIN gtf g ON c.bare_ensg = g.ensg) AS overlap,
          (SELECT count(*) FROM gtf g
            WHERE g.ensg NOT IN (SELECT bare_ensg FROM con_m WHERE bare_ensg IS NOT NULL)
          ) AS gtf_only,
          (SELECT count(*) FROM con_m
            WHERE bare_ensg IS NOT NULL
              AND bare_ensg NOT IN (SELECT ensg FROM gtf)
          ) AS constraint_only_with_ensg,
          (SELECT count(*) FROM con_m WHERE bare_ensg IS NULL) AS constraint_only_null_ensg,
          (SELECT count(*) FROM con_m
            WHERE bare_ensg IS NOT NULL
              AND bare_ensg NOT IN (SELECT ensg FROM read_parquet(?))
          ) AS constraint_ensg_absent_from_gencode
        """,
        [intervals, gene_constraint, intervals],
    )[0]
    gtf_tag_genes = int(
        one(
            con,
            """
            SELECT count(DISTINCT ensg)
            FROM read_parquet(?)
            WHERE list_contains(tags, 'MANE_Select')
            """,
            [transcripts],
        )
    )
    gtf_tag_tx = int(
        one(
            con,
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE list_contains(tags, 'MANE_Select')
            """,
            [transcripts],
        )
    )
    cross_mane = records(
        con,
        """
        SELECT
          count(*) FILTER (WHERE mane_enst_gencode IS NOT NULL) AS gencode_rows,
          count(*) FILTER (WHERE mane_enst_constraint IS NOT NULL) AS constraint_enst_rows,
          count(*) FILTER (WHERE mane_enst_disagree) AS enst_disagree,
          count(*) FILTER (
            WHERE mane_enst_gencode IS NOT NULL AND mane_enst_constraint IS NOT NULL
              AND NOT mane_enst_disagree
          ) AS enst_agree,
          count(*) FILTER (WHERE mane_enst_source = 'gencode') AS source_gencode,
          count(*) FILTER (WHERE mane_enst_source = 'constraint') AS source_constraint,
          count(*) FILTER (WHERE mane_enst IS NULL) AS mane_enst_null,
          count(*) FILTER (WHERE mane_refseq IS NOT NULL) AS with_refseq,
          count(*) FILTER (WHERE mane_refseq LIKE '%;%') AS multi_refseq,
          count(*) FILTER (
            WHERE position(';' IN coalesce(mane_enst_gencode, '')) > 0
          ) AS multi_gencode_enst,
          count(*) FILTER (
            WHERE position(';' IN coalesce(mane_enst_constraint, '')) > 0
          ) AS multi_constraint_enst
        FROM read_parquet(?)
        """,
        [genes],
    )[0]
    disagree_rows = records(
        con,
        """
        SELECT symbol, ensg, mane_enst_gencode, mane_enst_constraint, mane_enst, mane_refseq
        FROM read_parquet(?)
        WHERE mane_enst_disagree
        ORDER BY symbol, ensg
        """,
        [genes],
    )
    constraint_only_examples = records(
        con,
        """
        SELECT c.gene AS symbol, c.bare_ensg, c.gene_id, c.transcript,
               g.ensg IS NOT NULL AS in_gencode,
               g.biotype, g.mane_enst_gencode, g.mane_enst_constraint,
               coalesce(g.mane_refseq, u.mane_refseq) AS mane_refseq
        FROM read_parquet(?) c
        LEFT JOIN read_parquet(?) g ON g.ensg = c.bare_ensg
        LEFT JOIN read_parquet(?) u
          ON u.source = 'constraint' AND u.gene_id = c.gene_id
        WHERE c.mane_select
          AND (c.bare_ensg IS NULL OR c.bare_ensg NOT IN (
                SELECT ensg FROM read_parquet(?) WHERE mane_select
              ))
        ORDER BY c.bare_ensg NULLS FIRST, c.gene
        """,
        [gene_constraint, genes, unmatched, intervals],
    )
    moved_mane = records(
        con,
        """
        WITH gtf AS (
          SELECT enst, ensg, symbol
          FROM read_parquet(?)
          WHERE list_contains(tags, 'MANE_Select')
        ),
        con_tx AS (
          SELECT bare_enst AS enst, bare_ensg AS ensg, gene AS symbol
          FROM read_parquet(?)
          WHERE mane_select AND bare_enst IS NOT NULL
        )
        SELECT g.symbol AS gtf_symbol, g.ensg AS gtf_ensg,
               c.symbol AS constraint_symbol, c.ensg AS constraint_ensg, g.enst,
               c.ensg IN (SELECT ensg FROM read_parquet(?)) AS constraint_ensg_in_gencode
        FROM gtf g
        JOIN con_tx c ON g.enst = c.enst
        WHERE g.ensg IS DISTINCT FROM c.ensg
        ORDER BY g.enst
        """,
        [transcripts, transcript_constraint, intervals],
    )
    pc_other_ensg = records(
        con,
        """
        WITH pc AS (
          SELECT ensg, symbol
          FROM read_parquet(?)
          WHERE biotype = 'protein_coding' AND NOT constraint_ensg_matched
        )
        SELECT p.symbol, p.ensg AS gencode_ensg, c.bare_ensg AS constraint_ensg, c.transcript
        FROM pc p
        JOIN read_parquet(?) c ON c.gene = p.symbol
        WHERE c.bare_ensg IS NOT NULL AND c.bare_ensg <> p.ensg
        ORDER BY p.symbol, p.ensg
        """,
        [genes, gene_constraint],
    )
    pc_null_ensg = records(
        con,
        """
        WITH pc AS (
          SELECT ensg, symbol
          FROM read_parquet(?)
          WHERE biotype = 'protein_coding' AND NOT constraint_ensg_matched
        )
        SELECT p.symbol, p.ensg AS gencode_ensg, c.gene_id, c.transcript
        FROM pc p
        JOIN read_parquet(?) c ON c.gene = p.symbol
        WHERE c.bare_ensg IS NULL
          AND NOT EXISTS (
            SELECT 1 FROM read_parquet(?) x
            WHERE x.gene = p.symbol AND x.bare_ensg IS NOT NULL
          )
        ORDER BY p.symbol
        """,
        [genes, gene_constraint, gene_constraint],
    )
    refseq_on_unmatched = int(
        one(
            con,
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE source = 'constraint' AND mane_refseq IS NOT NULL
            """,
            [unmatched],
        )
    )
    source_refseq_n = int(
        one(
            con,
            """
            SELECT count(*)
            FROM read_parquet(?)
            WHERE mane_select
              AND (starts_with(transcript, 'NM_') OR starts_with(transcript, 'NR_'))
            """,
            [transcript_constraint],
        )
    )
    null_ensg_constraint = int(
        one(
            con,
            "SELECT count(*) FROM read_parquet(?) WHERE bare_ensg IS NULL",
            [gene_constraint],
        )
    )
    pc_symbol_gap = records(
        con,
        """
        WITH pc AS (
          SELECT ensg, symbol
          FROM read_parquet(?)
          WHERE biotype = 'protein_coding' AND NOT constraint_ensg_matched
        ),
        con AS (
          SELECT gene AS symbol, bare_ensg
          FROM read_parquet(?)
          WHERE gene IS NOT NULL
        )
        SELECT
          (SELECT count(*) FROM pc) AS unmatched_pc,
          (SELECT count(*) FROM pc p
            WHERE EXISTS (
              SELECT 1 FROM con c
              WHERE c.symbol = p.symbol AND c.bare_ensg IS NOT NULL AND c.bare_ensg <> p.ensg
            )
          ) AS symbol_hits_other_ensg,
          (SELECT count(*) FROM pc p
            WHERE EXISTS (
              SELECT 1 FROM con c WHERE c.symbol = p.symbol AND c.bare_ensg IS NULL
            )
            AND NOT EXISTS (
              SELECT 1 FROM con c
              WHERE c.symbol = p.symbol AND c.bare_ensg IS NOT NULL
            )
          ) AS symbol_hits_null_ensg_only,
          (SELECT count(*) FROM pc p
            WHERE NOT EXISTS (SELECT 1 FROM con c WHERE c.symbol = p.symbol)
          ) AS no_constraint_symbol
        """,
        [genes, gene_constraint],
    )[0]
    status = "pass" if row_n == EXPECTED_GENES == interval_n == distinct_version == distinct_ensg else "FAIL"
    overlap = int(mane["overlap"])
    gtf_genes = int(mane["gtf_genes"])
    constraint_genes = int(mane["constraint_genes"])
    gtf_only = int(mane["gtf_only"])
    constraint_only = int(mane["constraint_only_with_ensg"]) + int(mane["constraint_only_null_ensg"])
    mane_identity = (
        overlap + gtf_only == gtf_genes and overlap + constraint_only == constraint_genes
    )
    lines = [
        "# Crosswalk",
        "",
        "Lane G. One row per versioned GENCODE gene id. "
        "Counts below were queried from the Parquet files after they were written.",
        "",
        f"Genes: `{GENE_OUT}` ({GENE_OUT.stat().st_size:,} bytes).",
        f"Unmatched: `{UNMATCHED_OUT}` ({UNMATCHED_OUT.stat().st_size:,} bytes).",
        "",
        "## Row count",
        "",
        f"Status: {status}.",
        "",
        md_table(
            ["check", "count"],
            [
                ["genes.parquet rows", row_n],
                ["distinct ensg_version", distinct_version],
                ["distinct ensg", distinct_ensg],
                ["interval genes.parquet rows", interval_n],
                ["expected", EXPECTED_GENES],
            ],
        ),
        "",
        "The crosswalk does not add a row for an NCBI id or a RefSeq id that failed to join. "
        "Those ids stay on the gene they matched, or in `unmatched.parquet`.",
        "",
        "## Join rules",
        "",
        "Primary key: bare ENSG (`coords.bare_ensg`). Symbol is not a gene key.",
        "",
        "ClinVar `GeneID` and `HGNC_ID` attach by exact `GeneSymbol` on single-gene rows "
        "(`GeneID` is digits, symbol has no space, comma, or semicolon). "
        "`GeneID` `-1` is not an NCBI gene id. Semicolon-separated symbols on those rows are split "
        "only to see which symbols failed to match. They do not assign an NCBI id.",
        "",
        "If one symbol has several NCBI GeneIDs, `ncbi_gene_id` keeps every id, semicolon-separated, "
        "sorted numerically. `ncbi_gene_id_primary` is the chosen id: the most allele rows, "
        "then the smaller numeric id on a tie. `ncbi_gene_id_rule` is `sole`, `most_alleles`, "
        "or `smallest_id`. The gene is not dropped. Every GENCODE gene with that symbol gets the same ClinVar ids. "
        "`symbol_n_gencode` is how many GENCODE genes share the symbol. "
        "`ncbi_clinvar_symbols` lists every ClinVar symbol that uses the primary NCBI id.",
        "",
        "`hgnc_id` is the GTF `hgnc_id` when the gene has one. Otherwise it is the ClinVar HGNC id "
        "when that symbol has exactly one. If ClinVar has several and the GTF has none, `hgnc_id` stays null "
        "and `hgnc_id_clinvar` keeps the list. No id is dropped. "
        "`hgnc_disagree` is true when both sources are present and the sets differ. "
        "Ids are compared as `HGNC:` plus the integer, so a bare number or leading zeros would not count as a disagreement.",
        "",
        "`mane_enst` is the bare ENST. The GENCODE tag `MANE_Select` wins when the gene has one. "
        "Otherwise the id comes from a constraint transcript with `mane_select` and an ENST id on that bare ENSG. "
        "`mane_enst_gencode` and `mane_enst_constraint` both stay populated. "
        "`mane_enst_disagree` is true when both are present and the sets differ. "
        "If a source has several ENST ids, the stored `mane_enst` is the lexicographically smallest and the full set stays in the source column. "
        "`mane_enst_version` is the GENCODE versioned id when the chosen source is GENCODE. Constraint ENST ids in this file have no version, so the version is null when the source is constraint.",
        "",
        "`mane_refseq` is the constraint transcript with `mane_select` whose id starts with `NM_` or `NR_`. "
        "RefSeq rows have a null bare ENSG. Pairing stays inside the constraint file: "
        "if the symbol has exactly one MANE ENST, the RefSeq id attaches to that transcript's bare ENSG; "
        "if it has no MANE ENST and exactly one bare ENSG on any constraint transcript, it attaches to that ENSG; "
        "otherwise it is unattached. It is not joined to GENCODE by symbol. Null when none attach.",
        "",
        "`constraint_ensg_matched` is true when the gene's bare ENSG equals `gene_constraint.bare_ensg`. "
        "A constraint row whose `bare_ensg` is null is unmatched on id, not a failed symbol join.",
        "",
        "`MANE_Plus_Clinical` is not MANE Select.",
        "",
        id_notes[0],
        "",
    ]
    lines.extend(build_notes)
    lines.extend(
        [
            "",
            "## Protein-coding constraint match",
            "",
            "Match is `biotype = protein_coding` and `constraint_ensg_matched`. "
            "An earlier partial join saw 18,567 of 20,097 protein-coding genes and about 1,530 misses. "
            "That figure is not a target. These counts are recomputed from `genes.parquet`.",
            "",
            md_table(
                ["set", "genes", "rate of protein-coding"],
                [
                    ["protein-coding", int(pc["n"]), "100%"],
                    ["matched to a constraint bare ENSG", int(pc["matched"]), pct(int(pc["matched"]), int(pc["n"]))],
                    ["unmatched", int(pc["unmatched"]), pct(int(pc["unmatched"]), int(pc["n"]))],
                    ["all biotypes matched", matched_all, ""],
                    ["all biotypes unmatched", row_n - matched_all, ""],
                ],
            ),
            "",
            "Protein-coding genes that missed on bare ENSG were not rescued by symbol. "
            "The next table only describes what a symbol comparison would have seen. It is not the join.",
            "",
            md_table(
                ["protein-coding unmatched on ENSG", "genes"],
                [
                    ["total", int(pc_symbol_gap["unmatched_pc"])],
                    ["symbol hits a constraint row with a different bare ENSG", int(pc_symbol_gap["symbol_hits_other_ensg"])],
                    ["symbol hits only a constraint row whose bare ENSG is null", int(pc_symbol_gap["symbol_hits_null_ensg_only"])],
                    ["no constraint row with that symbol", int(pc_symbol_gap["no_constraint_symbol"])],
                ],
            ),
            "",
            "These symbol hits were not used as a join.",
            "",
        ]
    )
    if pc_other_ensg:
        lines.extend(
            [
                md_table(
                    ["symbol", "GENCODE ENSG", "constraint ENSG", "constraint transcript"],
                    [
                        [row["symbol"], row["gencode_ensg"], row["constraint_ensg"], row["transcript"]]
                        for row in pc_other_ensg
                    ],
                ),
                "",
            ]
        )
    if pc_null_ensg:
        lines.extend(
            [
                "Protein-coding symbols that meet only a constraint row with a null bare ENSG. "
                "That row is unmatched on id.",
                "",
                md_table(
                    ["symbol", "GENCODE ENSG", "constraint gene_id", "constraint transcript"],
                    [
                        [row["symbol"], row["gencode_ensg"], row["gene_id"], row["transcript"]]
                        for row in pc_null_ensg
                    ],
                ),
                "",
            ]
        )
    refseq_accounted = int(cross_mane["with_refseq"]) + refseq_on_unmatched
    refseq_identity = refseq_accounted == source_refseq_n
    lines.extend(
        [
            "## Duplicate symbols",
            "",
            f"Symbols on more than one GENCODE gene: {int(dup['symbols']):,}, covering {int(dup['genes']):,} genes. "
            f"Symbols that also occur on a protein-coding gene: {int(dup['symbols_touching_protein']):,}. "
            f"Symbols shared by more than one protein-coding gene: {pc_dup_symbols:,} "
            f"({pc_dup_genes:,} protein-coding genes).",
            "",
            "Every gene is kept. ClinVar ids for that symbol are copied onto each gene. "
            "`symbol_n_gencode` records the fan-out.",
            "",
            md_table(
                ["genes with the symbol", "symbols"],
                [[int(row["n_genes"]), int(row["n_symbols"])] for row in multiplicity],
            ),
            "",
            "Largest duplicate symbols:",
            "",
            md_table(
                ["symbol", "genes", "protein_coding", "biotypes"],
                [
                    [row["symbol"], int(row["n_genes"]), int(row["protein_coding"]), row["biotypes"]]
                    for row in largest_dups
                ],
            ),
            "",
        ]
    )
    if pc_dup_list:
        lines.extend(
            [
                "Protein-coding symbols that are not unique:",
                "",
                md_table(
                    ["symbol", "genes", "ensgs"],
                    [[row["symbol"], int(row["n_genes"]), row["ensgs"]] for row in pc_dup_list],
                ),
                "",
            ]
        )
    else:
        lines.extend(["No protein-coding symbol is shared by two protein-coding genes.", ""])
    lines.extend(
        [
            "## ClinVar symbols with no GENCODE symbol",
            "",
            "Rows are in `unmatched.parquet` with `source = clinvar`. "
            "`n` is allele rows. For `no_gencode_symbol`, `n_single` is single-gene rows and "
            "`n_multi` is multigene rows whose semicolon-separated symbol list contains the symbol. "
            "A symbol that matches any GENCODE gene is absent here, including when several genes share it.",
            "",
            md_table(
                ["reason", "keys", "allele rows", "single-gene rows", "multigene rows"],
                [
                    [
                        row["reason"],
                        int(row["symbols"]),
                        int(row["alleles"]),
                        "" if row["single_alleles"] is None else int(row["single_alleles"]),
                        "" if row["multi_alleles"] is None else int(row["multi_alleles"]),
                    ]
                    for row in clinvar_symbols
                ],
            ),
            "",
            "`no_gencode_symbol` is the ClinVar symbol miss. "
            "`not_a_gene_symbol` is prose such as a dosage-coverage label. "
            "`missing_symbol` is GeneSymbol `-`. Those two are not symbol joins.",
            "",
            "## Constraint ids with no GENCODE gene",
            "",
            f"`gene_constraint` rows with null `bare_ensg`: {null_ensg_constraint:,}. "
            "A null bare ENSG is unmatched on id. It is not counted as a symbol miss.",
            "",
            md_table(
                ["reason", "gene rows", "transcript rows with that gene_id", "mane_select gene rows"],
                [
                    [
                        row["reason"],
                        int(row["gene_rows"]),
                        int(row["transcript_rows"]),
                        int(row["mane_gene_rows"]),
                    ]
                    for row in constraint_unmatched
                ],
            ),
            "",
            "`n` on a constraint unmatched row is the number of `transcript_constraint` rows with that `gene_id`. "
            "NCBI ids that the constraint loader collapsed onto an Ensembl gene are not listed again here. "
            "The gene id that remains is the bare ENSG, and it is unmatched only when that ENSG is absent from GENCODE.",
            "",
            "## NCBI and HGNC",
            "",
            md_table(
                ["metric", "count"],
                [
                    ["genes with an NCBI GeneID", int(ncbi["genes_with_ncbi"])],
                    ["genes with no NCBI GeneID", int(ncbi["genes_without_ncbi"])],
                    ["genes whose symbol has more than one NCBI GeneID", int(ncbi["genes_multi_ncbi"])],
                    ["primary rule sole", int(ncbi["rule_sole"])],
                    ["primary rule most_alleles", int(ncbi["rule_most"])],
                    ["primary rule smallest_id", int(ncbi["rule_smallest"])],
                    ["genes whose primary NCBI id is shared by more than one ClinVar symbol", int(ncbi["genes_shared_ncbi_symbol"])],
                    ["distinct primary NCBI ids shared by more than one ClinVar symbol", shared_ncbi_n],
                    ["genes with GTF and ClinVar HGNC", int(hgnc["both"])],
                    ["genes with GTF HGNC only", int(hgnc["gtf_only"])],
                    ["genes with ClinVar HGNC only", int(hgnc["clinvar_only"])],
                    ["genes whose resolved hgnc_id is null", int(hgnc["resolved_null"])],
                    ["genes whose ClinVar HGNC list has more than one id", int(hgnc["multi_clinvar_hgnc"])],
                    ["HGNC disagreements", int(hgnc["disagree"])],
                ],
            ),
            "",
        ]
    )
    if shared_ncbi:
        lines.extend(
            [
                "NCBI ids used by more than one ClinVar symbol (up to 20):",
                "",
                md_table(
                    ["ncbi_gene_id", "clinvar symbols", "gencode genes carrying that id"],
                    [
                        [row["ncbi"], row["symbols"], int(row["gencode_genes"])]
                        for row in shared_ncbi
                    ],
                ),
                "",
            ]
        )
    if int(hgnc["disagree"]) == 0:
        lines.append("No gene has both a GTF HGNC id and a different ClinVar HGNC id.")
        lines.append("")
    else:
        lines.extend(
            [
                f"HGNC disagreements: {int(hgnc['disagree']):,}. "
                "Both ids are kept. `hgnc_id` stays the GTF id. First 20:",
                "",
                md_table(
                    ["symbol", "ensg", "hgnc_id_gtf", "hgnc_id_clinvar", "hgnc_id"],
                    [
                        [
                            row["symbol"],
                            row["ensg"],
                            row["hgnc_id_gtf"],
                            row["hgnc_id_clinvar"],
                            row["hgnc_id"],
                        ]
                        for row in hgnc_examples
                    ],
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## MANE",
            "",
            "GENCODE MANE Select is the transcript tag `MANE_Select`, not `MANE_Plus_Clinical`. "
            "Constraint MANE is `gene_constraint.mane_select`. "
            "Overlap is a constraint MANE row whose `bare_ensg` is a GENCODE MANE gene. "
            "The two difference sets are GENCODE-only and constraint-only. Neither is dropped from the note.",
            "",
            md_table(
                ["set", "genes"],
                [
                    ["GENCODE genes with mane_select", gtf_genes],
                    ["GENCODE genes with a MANE_Select transcript tag", gtf_tag_genes],
                    ["GENCODE MANE_Select transcripts", gtf_tag_tx],
                    ["constraint gene rows with mane_select", constraint_genes],
                    ["overlap on bare ENSG", overlap],
                    ["GENCODE only", gtf_only],
                    ["constraint only", constraint_only],
                    ["constraint only, bare ENSG set, not a GENCODE MANE gene", int(mane["constraint_only_with_ensg"])],
                    ["constraint only, bare ENSG null (unmatched on id)", int(mane["constraint_only_null_ensg"])],
                    ["constraint MANE bare ENSG absent from GENCODE entirely", int(mane["constraint_ensg_absent_from_gencode"])],
                ],
            ),
            "",
            f"Set identity (overlap + GENCODE only = GENCODE; overlap + constraint only = constraint): "
            f"{'pass' if mane_identity else 'FAIL'}.",
            "",
            "On the crosswalk row, a GENCODE MANE transcript is preferred for `mane_enst` and the constraint id is still stored. "
            "A constraint-only MANE gene that joins on bare ENSG keeps the constraint ENST when it has one, and the RefSeq id when the pair attached.",
            "",
            md_table(
                ["crosswalk count", "genes"],
                [
                    ["mane_enst_gencode set", int(cross_mane["gencode_rows"])],
                    ["mane_enst_constraint set", int(cross_mane["constraint_enst_rows"])],
                    ["both ENST ids present and equal", int(cross_mane["enst_agree"])],
                    ["both ENST ids present and different", int(cross_mane["enst_disagree"])],
                    ["mane_enst chosen from GENCODE", int(cross_mane["source_gencode"])],
                    ["mane_enst chosen from constraint", int(cross_mane["source_constraint"])],
                    ["mane_enst null", int(cross_mane["mane_enst_null"])],
                    ["mane_refseq set", int(cross_mane["with_refseq"])],
                    ["more than one mane_refseq", int(cross_mane["multi_refseq"])],
                    ["more than one GENCODE MANE ENST", int(cross_mane["multi_gencode_enst"])],
                    ["more than one constraint MANE ENST", int(cross_mane["multi_constraint_enst"])],
                    ["constraint NM_/NR_ mane transcripts in the transcript file", source_refseq_n],
                    ["those RefSeq ids stored on unmatched constraint rows", refseq_on_unmatched],
                    ["crosswalk plus unmatched RefSeq ids", refseq_accounted],
                ],
            ),
            "",
            "Every constraint `NM_` or `NR_` MANE transcript is either `mane_refseq` on the GENCODE gene "
            "whose bare ENSG the constraint file paired it to, or `mane_refseq` on the unmatched constraint row "
            "when that ENSG is missing or null. "
            f"Accounting: {'pass' if refseq_identity else 'FAIL'} "
            f"({int(cross_mane['with_refseq']):,} + {refseq_on_unmatched:,} = {refseq_accounted:,}; "
            f"transcript file {source_refseq_n:,}).",
            "",
        ]
    )
    if not disagree_rows:
        lines.append(
            "No gene has a GENCODE MANE ENST and a different constraint MANE ENST. "
            "The gene-set differences above are genes present on only one side, not conflicting transcript ids."
        )
        lines.append("")
    else:
        shown = disagree_rows if len(disagree_rows) <= 100 else disagree_rows[:20]
        lines.append(
            f"Transcript id disagreements: {len(disagree_rows):,}. "
            "`mane_enst` uses the GENCODE id. The constraint id stays in `mane_enst_constraint`."
        )
        if len(disagree_rows) > 100:
            lines.append("First 20:")
        lines.append("")
        lines.append(
            md_table(
                ["symbol", "ensg", "gencode", "constraint", "chosen", "mane_refseq"],
                [
                    [
                        row["symbol"],
                        row["ensg"],
                        row["mane_enst_gencode"],
                        row["mane_enst_constraint"],
                        row["mane_enst"],
                        row["mane_refseq"] or "",
                    ]
                    for row in shown
                ],
            )
        )
        lines.append("")
    lines.append(
        f"Constraint MANE rows that are not GENCODE MANE genes: {len(constraint_only_examples):,}. "
        "Null `bare_ensg` rows are unmatched on id. "
        "`in_gencode` false means that bare ENSG is not a GENCODE gene, so the row is not in `genes.parquet`. "
        "A same symbol in GENCODE was not given this constraint id."
    )
    lines.append("")
    if len(constraint_only_examples) <= 40:
        lines.append(
            md_table(
                ["symbol", "bare_ensg", "in_gencode", "gene_id", "transcript", "biotype", "gencode ENST", "constraint ENST", "refseq"],
                [
                    [
                        row["symbol"] or "",
                        row["bare_ensg"] or "",
                        "yes" if row["in_gencode"] else "no",
                        row["gene_id"],
                        row["transcript"],
                        row["biotype"] or "",
                        row["mane_enst_gencode"] or "",
                        row["mane_enst_constraint"] or "",
                        row["mane_refseq"] or "",
                    ]
                    for row in constraint_only_examples
                ],
            )
        )
        lines.append("")
    else:
        lines.append("The constraint-only set is larger than 40 and is not listed row by row. The split counts are above.")
        lines.append("")
    if moved_mane:
        lines.extend(
            [
                f"The same MANE ENST is on different ENSG ids in GENCODE and in constraint: {len(moved_mane):,}. "
                "Each source stays on its own gene. The transcript was not copied across.",
                "",
                md_table(
                    ["ENST", "GENCODE symbol", "GENCODE ENSG", "constraint symbol", "constraint ENSG", "constraint ENSG in GENCODE"],
                    [
                        [
                            row["enst"],
                            row["gtf_symbol"],
                            row["gtf_ensg"],
                            row["constraint_symbol"] or "",
                            row["constraint_ensg"] or "",
                            "yes" if row["constraint_ensg_in_gencode"] else "no",
                        ]
                        for row in moved_mane
                    ],
                ),
                "",
            ]
        )
    else:
        lines.append("No MANE ENST is attached to a different ENSG in GENCODE than in the constraint file.")
        lines.append("")
    lines.extend(
        [
            "## Column notes",
            "",
            "`ncbi_gene_id` is every ClinVar GeneID for the symbol. "
            "`ncbi_gene_id_primary` is the chosen one. They are equal when the symbol has one GeneID.",
            "",
            "`hgnc_id` is the resolved id. `hgnc_id_gtf` and `hgnc_id_clinvar` are the sources, including when they agree.",
            "",
            "`mane_enst` is bare. `mane_enst_version` has the GENCODE version suffix when the source is GENCODE.",
            "",
        ]
    )
    QC_OUT.write_text("\n".join(lines), encoding="utf-8")
    if status != "pass" or not mane_identity or not refseq_identity:
        raise SystemExit(
            f"crosswalk QC failed: row status {status}; "
            f"mane identity {mane_identity}; refseq identity {refseq_identity}"
        )


def main() -> None:
    for path in (GENES, TRANSCRIPTS, GENE_CONSTRAINT, TRANSCRIPT_CONSTRAINT, ALLELES):
        if not path.is_file():
            raise SystemExit(f"missing input: {path}")
    con = connect()
    id_notes = check_identifiers(con)
    genes = load_genes(con)
    if len(genes) != EXPECTED_GENES:
        raise SystemExit(f"interval genes are {len(genes)}, expected {EXPECTED_GENES}")
    gencode_mane, mane_versions = load_gencode_mane(con)
    constraint = load_constraint(con)
    mane_enst_by_ensg, refseq_by_ensg, unattached = pair_refseq(constraint)
    clinvar_raw = load_clinvar(con)
    clinvar, symbols_by_ncbi = clinvar_by_symbol(clinvar_raw["singles"])
    gene_rows = build_gene_rows(
        genes,
        clinvar,
        symbols_by_ncbi,
        gencode_mane,
        mane_versions,
        mane_enst_by_ensg,
        refseq_by_ensg,
        constraint["matched"],
    )
    unmatched_rows = build_unmatched(genes, clinvar_raw, constraint, refseq_by_ensg)
    if len(gene_rows) != EXPECTED_GENES:
        raise SystemExit(f"crosswalk rows {len(gene_rows)}, expected {EXPECTED_GENES}")
    versions = [row["ensg_version"] for row in gene_rows]
    if len(set(versions)) != len(versions):
        raise SystemExit("duplicate versioned gene id in the crosswalk")
    write_outputs(gene_rows, unmatched_rows)
    multi_ncbi_symbols = sum(1 for slot in clinvar.values() if len(slot["ids"]) > 1)
    build_notes = [
        f"ClinVar single-gene symbols: {len(clinvar):,}. "
        f"Symbols with more than one NCBI GeneID: {multi_ncbi_symbols:,}.",
        f"Single-gene ClinVar rows rejected as non-symbols: {clinvar_raw['rejected_single_rows']:,}. "
        f"Multigene rows with an empty or prose token: {clinvar_raw['bad_multigene_rows']:,}.",
        f"MANE RefSeq transcripts left unattached: {len(unattached):,}. "
        "They have no unique bare ENSG inside the constraint file, so they are not symbol-joined to GENCODE.",
    ]
    if unattached:
        shown = ", ".join(
            f"{row['symbol'] or '(no symbol)'} {row['transcript']}" for row in unattached[:15]
        )
        build_notes.append(f"Unattached MANE RefSeq transcripts (up to 15): {shown}.")
    write_qc(con, id_notes, build_notes)
    print(f"wrote {GENE_OUT}")
    print(f"wrote {UNMATCHED_OUT}")
    print(f"wrote {QC_OUT}")
    print(f"genes {len(gene_rows):,} unmatched {len(unmatched_rows):,}")


if __name__ == "__main__":
    main()
