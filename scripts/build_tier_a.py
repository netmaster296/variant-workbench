# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Build the Tier A variant-key list.

Tier A is the union of ClinVar keys, MaveDB GRCh38 keys, and lead eQTL and
sQTL keys. A GWAS row has no ref or alt. It adds a key only when its rsID
matches a ClinVar rsID that already has a key. Unmatched GWAS rsIDs that
still have a position are written aside. This pass does not scan gnomAD for
those IDs.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import duckdb

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coords import canon_chrom, variant_key  # noqa: E402

DEDUP = ROOT / "warehouse" / "clinvar" / "allele_dedup.parquet"
ASSOC = ROOT / "warehouse" / "gwas" / "associations.parquet"
SCORED = ROOT / "warehouse" / "mavedb" / "mapped_variants_scored.parquet"
EGENES = ROOT / "warehouse" / "gtex" / "egenes.parquet"
SGENES = ROOT / "warehouse" / "gtex" / "sgenes.parquet"
OUT_DIR = ROOT / "warehouse" / "variants"
KEYS = OUT_DIR / "tier_a_keys.tsv"
UNRESOLVED = OUT_DIR / "tier_a_unresolved.tsv"
TMP = ROOT / "warehouse" / "tmp"

TRUTH_KEYS = (
    "chr11:126275389:C:T",
    "chr7:117559590:ATCT:A",
    "chr22:10510183:A:G",
)
REASON = (
    "no ClinVar variant_key for this rsid; "
    "genome-wide gnomAD ID scan deferred because bcftools is not available "
    "and this pass does not stream every sites VCF to recover alleles"
)


def sql_path(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET temp_directory = '{sql_path(TMP)}'")
    con.execute("SET memory_limit = '64GB'")
    con.execute("SET threads = 8")
    con.execute("SET preserve_insertion_order = false")

    dedup = sql_path(DEDUP)
    assoc = sql_path(ASSOC)
    scored = sql_path(SCORED)
    egenes = sql_path(EGENES)
    sgenes = sql_path(SGENES)

    print("collect catalog keys", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE src AS
        SELECT variant_key, 'clinvar' AS src
        FROM read_parquet('{dedup}')
        WHERE variant_key IS NOT NULL
        UNION ALL
        SELECT variant_key, 'mavedb'
        FROM read_parquet('{scored}')
        WHERE variant_key IS NOT NULL
        UNION ALL
        SELECT variant_key, 'eqtl'
        FROM read_parquet('{egenes}')
        WHERE variant_key IS NOT NULL
        UNION ALL
        SELECT variant_key, 'sqtl'
        FROM read_parquet('{sgenes}')
        WHERE variant_key IS NOT NULL
        """
    )
    for name, count in con.execute(
        "SELECT src, count(DISTINCT variant_key) FROM src GROUP BY 1 ORDER BY 1"
    ).fetchall():
        print(f"QC keys_{name}={count}", flush=True)

    print("match GWAS rsids to ClinVar keys", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE clin_rs AS
        SELECT rsid, variant_key, chrom, pos
        FROM read_parquet('{dedup}')
        WHERE rsid IS NOT NULL AND variant_key IS NOT NULL
        """
    )
    con.execute(
        f"""
        INSERT INTO src
        SELECT DISTINCT c.variant_key, 'gwas'
        FROM read_parquet('{assoc}') a
        JOIN clin_rs c ON a.rsid = c.rsid
        WHERE a.rsid IS NOT NULL
        """
    )
    # The frequency truth-set site is a gnomAD sample. It stays in the annotated
    # table even when no catalog already carried that key.
    present = set()
    for key in TRUTH_KEYS:
        hit = con.execute("SELECT count(*) FROM src WHERE variant_key = ?", [key]).fetchone()[0]
        if hit:
            present.add(key)
        print(f"QC truth_key_in_catalogs {key}={hit > 0}", flush=True)
    for key in TRUTH_KEYS:
        if key in present:
            continue
        con.execute("INSERT INTO src VALUES (?, 'truth_set')", [key])
        print(f"QC truth_key_added={key}", flush=True)

    con.execute(
        """
        CREATE TEMP TABLE keyed AS
        SELECT
          variant_key,
          concat_ws(
            ',',
            max(CASE WHEN src = 'clinvar' THEN 'clinvar' END),
            max(CASE WHEN src = 'gwas' THEN 'gwas' END),
            max(CASE WHEN src = 'mavedb' THEN 'mavedb' END),
            max(CASE WHEN src = 'eqtl' THEN 'eqtl' END),
            max(CASE WHEN src = 'sqtl' THEN 'sqtl' END),
            max(CASE WHEN src = 'truth_set' THEN 'truth_set' END)
          ) AS source
        FROM src
        GROUP BY variant_key
        """
    )
    bad = con.execute(
        """
        SELECT count(*)
        FROM keyed
        WHERE variant_key <> concat_ws(
          ':',
          split_part(variant_key, ':', 1),
          split_part(variant_key, ':', 2),
          split_part(variant_key, ':', 3),
          split_part(variant_key, ':', 4)
        )
        OR split_part(variant_key, ':', 5) <> ''
        OR try_cast(split_part(variant_key, ':', 2) AS BIGINT) IS NULL
        """
    ).fetchone()[0]
    print(f"QC keys_not_four_fields={bad}", flush=True)
    if bad:
        raise SystemExit(f"{bad} variant keys are not chrom:pos:ref:alt")

    n_keys = con.execute("SELECT count(*) FROM keyed").fetchone()[0]
    print(f"QC tier_a_keys={n_keys}", flush=True)
    for key in TRUTH_KEYS:
        source = con.execute("SELECT source FROM keyed WHERE variant_key = ?", [key]).fetchone()
        if source is None:
            raise SystemExit(f"missing truth key {key}")
        chrom, pos, ref, alt = key.split(":")
        rebuilt = variant_key(canon_chrom(chrom), pos, ref, alt)
        if rebuilt != key:
            raise SystemExit(f"coords.variant_key rebuilt {rebuilt} from {key}")
        print(f"QC truth_source {key}={source[0]}", flush=True)

    print("write keys", flush=True)
    tab = "\t"
    con.execute(
        f"""
        COPY (
          SELECT
            split_part(variant_key, ':', 1) AS chrom,
            try_cast(split_part(variant_key, ':', 2) AS BIGINT) AS pos,
            split_part(variant_key, ':', 3) AS ref,
            split_part(variant_key, ':', 4) AS alt,
            variant_key,
            source
          FROM keyed
          ORDER BY
            CASE split_part(variant_key, ':', 1)
              WHEN 'chrX' THEN 23
              WHEN 'chrY' THEN 24
              WHEN 'chrM' THEN 25
              ELSE try_cast(regexp_extract(split_part(variant_key, ':', 1), '[0-9]+') AS INTEGER)
            END,
            2, 3, 4
        ) TO '{sql_path(KEYS)}' (HEADER, DELIMITER '{tab}')
        """
    )

    print("write unresolved GWAS rsids", flush=True)
    con.execute(
        f"""
        COPY (
          SELECT
            a.rsid,
            a.chrom,
            a.pos,
            count(*)::BIGINT AS n_associations,
            '{REASON}' AS reason
          FROM read_parquet('{assoc}') a
          WHERE a.rsid IS NOT NULL
            AND a.chrom IS NOT NULL
            AND a.pos IS NOT NULL
            AND NOT EXISTS (
              SELECT 1 FROM clin_rs c WHERE c.rsid = a.rsid
            )
          GROUP BY a.rsid, a.chrom, a.pos
          ORDER BY
            CASE a.chrom
              WHEN 'chrX' THEN 23
              WHEN 'chrY' THEN 24
              WHEN 'chrM' THEN 25
              ELSE try_cast(regexp_extract(a.chrom, '[0-9]+') AS INTEGER)
            END,
            a.pos,
            a.rsid
        ) TO '{sql_path(UNRESOLVED)}' (HEADER, DELIMITER '{tab}')
        """
    )

    gwas_stats = con.execute(
        f"""
        SELECT
          count(*) FILTER (WHERE rsid IS NOT NULL),
          count(*) FILTER (WHERE rsid IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL),
          count(*) FILTER (WHERE rsid IS NOT NULL AND (chrom IS NULL OR pos IS NULL)),
          count(DISTINCT rsid) FILTER (WHERE rsid IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL)
        FROM read_parquet('{assoc}')
        """
    ).fetchone()
    print(f"QC gwas_rows_with_rsid={gwas_stats[0]}", flush=True)
    print(f"QC gwas_rows_with_rsid_and_pos={gwas_stats[1]}", flush=True)
    print(f"QC gwas_rows_with_rsid_without_pos={gwas_stats[2]}", flush=True)
    print(f"QC gwas_distinct_rsid_with_pos={gwas_stats[3]}", flush=True)
    unresolved_rows, unresolved_rs = con.execute(
        f"""
        SELECT count(*), count(DISTINCT rsid)
        FROM read_csv('{sql_path(UNRESOLVED)}', header = true, delim = '{tab}')
        """
    ).fetchone()
    print(f"QC unresolved_rows={unresolved_rows}", flush=True)
    print(f"QC unresolved_rsids={unresolved_rs}", flush=True)
    matched_rs = con.execute(
        f"""
        SELECT count(DISTINCT a.rsid)
        FROM read_parquet('{assoc}') a
        SEMI JOIN clin_rs c ON a.rsid = c.rsid
        WHERE a.rsid IS NOT NULL AND a.chrom IS NOT NULL AND a.pos IS NOT NULL
        """
    ).fetchone()[0]
    print(f"QC gwas_rsids_matched_clinvar={matched_rs}", flush=True)
    multi = con.execute(
        """
        SELECT count(*) FROM (
          SELECT rsid FROM clin_rs GROUP BY rsid HAVING count(DISTINCT variant_key) > 1
        )
        """
    ).fetchone()[0]
    print(f"QC clinvar_rsids_with_multiple_keys={multi}", flush=True)
    disagree = con.execute(
        f"""
        WITH gwas_pos AS (
          SELECT DISTINCT rsid, chrom, pos
          FROM read_parquet('{assoc}')
          WHERE rsid IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL
        ),
        matched AS (
          SELECT g.rsid, g.chrom, g.pos
          FROM gwas_pos g
          SEMI JOIN clin_rs c ON g.rsid = c.rsid
        )
        SELECT count(*)
        FROM matched m
        WHERE NOT EXISTS (
          SELECT 1 FROM clin_rs c
          WHERE c.rsid = m.rsid AND c.chrom = m.chrom AND c.pos = m.pos
        )
        """
    ).fetchone()[0]
    print(f"QC gwas_positions_whose_clinvar_keys_differ={disagree}", flush=True)
    for label, count in con.execute(
        """
        SELECT label, count(*)
        FROM (
          SELECT unnest(string_split(source, ',')) AS label FROM keyed
        )
        GROUP BY 1
        ORDER BY 1
        """
    ).fetchall():
        print(f"QC source_{label}={count}", flush=True)
    con.close()
    print(f"wrote {KEYS}", flush=True)
    print(f"wrote {UNRESOLVED}", flush=True)


if __name__ == "__main__":
    build()
