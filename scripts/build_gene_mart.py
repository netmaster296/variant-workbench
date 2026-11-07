# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Build the gene mart from Wave 1 Parquet. One row per GENCODE gene.

Constraint joins on bare ENSG only. The protein-coding misses are not
filled in by symbol. ClinVar counts join on a gene-symbol token, which is
a secondary key and collides when several GENCODE genes share a symbol.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coords import bare_ensg, canon_chrom  # noqa: E402

GENES = ROOT / "warehouse" / "intervals" / "genes.parquet"
CROSSWALK = ROOT / "warehouse" / "crosswalk" / "genes.parquet"
CONSTRAINT = ROOT / "warehouse" / "constraint" / "gene_constraint.parquet"
EXPRESSION = ROOT / "warehouse" / "gtex" / "gene_expression.parquet"
DEDUP = ROOT / "warehouse" / "clinvar" / "allele_dedup.parquet"
LOCI = ROOT / "warehouse" / "gwas" / "loci.parquet"
COVERAGE = ROOT / "warehouse" / "mavedb" / "gene_coverage.parquet"
SCORED = ROOT / "warehouse" / "mavedb" / "mapped_variants_scored.parquet"
EGENES = ROOT / "warehouse" / "gtex" / "egenes.parquet"
SGENES = ROOT / "warehouse" / "gtex" / "sgenes.parquet"
CDS = ROOT / "warehouse" / "intervals" / "cds.parquet"
OUT_DIR = ROOT / "warehouse" / "marts"
MART = OUT_DIR / "gene_mart.parquet"
CODING = OUT_DIR / "gene_mart_protein_coding.parquet"
TMP = ROOT / "warehouse" / "tmp"

EXPECTED_GENES = 78_691
EXPECTED_CODING = 20_097

CLASSES = (
    "germline_pathogenic",
    "likely_pathogenic",
    "vus",
    "benign",
    "likely_benign",
    "conflicting",
    "drug_response",
)
STAR_LABELS = ("0", "1", "2", "3", "4", "null")
MISSING = {"", ".", "-", "na", "n/a", "none", "not provided"}

TRUTH_GENES = (
    ("A1BG", "ENST00000263100", 1.34, 7),
    ("XPO6", "ENST00000304658", 0.04, 0),
    ("MED13L", "ENST00000281928", 0.041, 0),
    ("EIF4G1", "ENST00000346169", 0.048, 0),
    ("CABLES2", "ENST00000279101", 1.0, 5),
    ("CFAP69", "ENST00000389297", 1.0, 5),
    ("CHI3L1", "ENST00000255409", 1.0, 5),
)


def clinvar_columns() -> list[str]:
    names = ["clinvar_n"]
    for name in CLASSES:
        names.append(f"clinvar_{name}_n")
        for star in STAR_LABELS:
            names.append(f"clinvar_{name}_star_{star}")
    return names


def pieces(value: str | None) -> list[str]:
    if not value:
        return []
    out = []
    for chunk in value.split("|"):
        for part in chunk.split(";"):
            for bit in part.split("/"):
                text = bit.strip().lower()
                if text and text not in MISSING:
                    out.append(text)
    return out


def is_likely_pathogenic(text: str) -> bool:
    return text.startswith("likely pathogenic")


def is_strong_pathogenic(text: str) -> bool:
    return "pathogenic" in text and not is_likely_pathogenic(text)


def is_likely_benign(text: str) -> bool:
    return text.startswith("likely benign")


def is_strong_benign(text: str) -> bool:
    return "benign" in text and not is_likely_benign(text)


def is_vus(text: str) -> bool:
    return "uncertain significance" in text or text.startswith("vus")


def star_index(stars) -> int:
    """0..4 match the star label. 5 is the null bucket."""
    if stars is None:
        return 5
    try:
        number = int(stars)
    except (TypeError, ValueError):
        return 5
    if 0 <= number <= 4:
        return number
    return 5


def symbol_tokens(value: str | None) -> list[str]:
    """ClinVar dedup joins distinct symbols with ';'. A token is not a substring."""
    if not value:
        return []
    tokens = []
    for part in str(value).split(";"):
        text = part.strip()
        if text and text not in {".", "-"}:
            tokens.append(text)
    return tokens


def classify(sig: str | None, flagged: bool) -> set[str]:
    """Star classes. Pathogenic buckets follow the dedup germline_pathogenic flag.

    That flag is already false for somatic origin, mixed origin, oncogenicity,
    and conflicting text. Likely pathogenic is split out of that flag. Drug
    response is counted as well as any other class. Benign, VUS, and
    conflicting are not limited to germline origin.
    """
    parts = pieces(sig)
    conflict = any("conflict" in part for part in parts)
    drug = any("drug response" in part for part in parts)
    strong_p = any(is_strong_pathogenic(part) for part in parts)
    likely_p = any(is_likely_pathogenic(part) for part in parts)
    strong_b = any(is_strong_benign(part) for part in parts)
    likely_b = any(is_likely_benign(part) for part in parts)
    vus = any(is_vus(part) for part in parts)
    found = set()
    if flagged:
        if strong_p or not likely_p:
            found.add("germline_pathogenic")
        else:
            found.add("likely_pathogenic")
    if conflict and not flagged:
        found.add("conflicting")
    if strong_b and not conflict:
        found.add("benign")
    elif likely_b and not conflict:
        found.add("likely_benign")
    if vus and not conflict and not flagged and not strong_b and not likely_b:
        found.add("vus")
    if drug:
        found.add("drug_response")
    return found


def _check_classify() -> None:
    assert classify("Pathogenic", True) == {"germline_pathogenic"}
    assert classify("Likely pathogenic", True) == {"likely_pathogenic"}
    assert classify("Pathogenic/Likely pathogenic", True) == {"germline_pathogenic"}
    assert classify("Pathogenic, low penetrance", True) == {"germline_pathogenic"}
    assert classify("Pathogenic; drug response", True) == {"germline_pathogenic", "drug_response"}
    assert "germline_pathogenic" not in classify("Pathogenic", False)
    assert classify("Conflicting classifications of pathogenicity", False) == {"conflicting"}
    assert classify("Uncertain significance", False) == {"vus"}
    assert classify("Benign/Likely benign", False) == {"benign"}
    assert classify("Likely benign", False) == {"likely_benign"}
    assert classify("drug response", False) == {"drug_response"}
    assert star_index(None) == 5
    assert star_index(4) == 4
    assert star_index(0) == 0
    assert symbol_tokens("KCNQ1;KCNQ1OT1") == ["KCNQ1", "KCNQ1OT1"]
    assert symbol_tokens("MIR4268 - EPHA4") == ["MIR4268 - EPHA4"]
    assert bare_ensg("ENSG00000121410.15") == "ENSG00000121410"
    assert canon_chrom("1") == "chr1"


def aggregate_clinvar(symbols: set[str]) -> tuple[pa.Table, dict]:
    width = len(CLASSES) * len(STAR_LABELS)
    counts: dict[str, list[int]] = {}
    stats = {
        "alleles": 0,
        "alleles_with_token": 0,
        "alleles_matched": 0,
        "alleles_unmatched": 0,
        "alleles_multi_token": 0,
        "flag_without_bucket": 0,
        "semicolon_symbols": 0,
    }

    def row_for(symbol: str) -> list[int]:
        found = counts.get(symbol)
        if found is None:
            found = [0] * (1 + width)
            counts[symbol] = found
        return found

    parquet = pq.ParquetFile(DEDUP)
    columns = ["GeneSymbol", "ClinicalSignificance", "review_stars", "germline_pathogenic"]
    for batch in parquet.iter_batches(batch_size=200_000, columns=columns):
        data = batch.to_pydict()
        for symbol, sig, stars, flagged in zip(
            data["GeneSymbol"],
            data["ClinicalSignificance"],
            data["review_stars"],
            data["germline_pathogenic"],
        ):
            stats["alleles"] += 1
            tokens = symbol_tokens(symbol)
            if len(tokens) > 1:
                stats["alleles_multi_token"] += 1
                stats["semicolon_symbols"] += 1
            if not tokens:
                stats["alleles_unmatched"] += 1
                continue
            stats["alleles_with_token"] += 1
            classes = classify(sig, bool(flagged))
            if flagged and not (classes & {"germline_pathogenic", "likely_pathogenic"}):
                stats["flag_without_bucket"] += 1
            star = star_index(stars)
            matched = False
            seen = set()
            for token in tokens:
                if token in seen:
                    continue
                seen.add(token)
                if token not in symbols:
                    continue
                matched = True
                bucket = row_for(token)
                bucket[0] += 1
                for name in classes:
                    class_index = CLASSES.index(name)
                    bucket[1 + class_index * len(STAR_LABELS) + star] += 1
            if matched:
                stats["alleles_matched"] += 1
            else:
                stats["alleles_unmatched"] += 1
    names = clinvar_columns()
    columns_out = {name: [] for name in ["symbol", *names]}
    for symbol, bucket in counts.items():
        columns_out["symbol"].append(symbol)
        columns_out["clinvar_n"].append(bucket[0])
        cursor = 1
        for name in CLASSES:
            stars = bucket[cursor : cursor + len(STAR_LABELS)]
            cursor += len(STAR_LABELS)
            columns_out[f"clinvar_{name}_n"].append(sum(stars))
            for label, value in zip(STAR_LABELS, stars):
                columns_out[f"clinvar_{name}_star_{label}"].append(value)
    schema_fields = [pa.field("symbol", pa.string())]
    schema_fields.extend(pa.field(name, pa.int64()) for name in names)
    table = pa.table(columns_out, schema=pa.schema(schema_fields))
    stats["symbols_with_alleles"] = table.num_rows
    return table, stats


def sql_path(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def build() -> None:
    _check_classify()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    gene_symbols = pq.read_table(GENES, columns=["symbol"]).column("symbol").to_pylist()
    symbols = {symbol for symbol in gene_symbols if symbol}
    bad_chrom = [
        chrom
        for chrom in pq.read_table(GENES, columns=["chrom"]).column("chrom").to_pylist()
        if canon_chrom(chrom) != chrom
    ]
    if bad_chrom:
        raise SystemExit(f"genes.parquet has non-canonical chromosomes: {sorted(set(bad_chrom))[:5]}")
    print(f"clinvar aggregate over {len(symbols)} symbols", flush=True)
    clinvar, clinvar_stats = aggregate_clinvar(symbols)
    for key, value in clinvar_stats.items():
        print(f"QC clinvar_{key}={value}", flush=True)
    if clinvar_stats["flag_without_bucket"]:
        raise SystemExit("germline_pathogenic rows did not land in a pathogenic bucket")

    con = duckdb.connect()
    con.execute(f"SET temp_directory = '{sql_path(TMP)}'")
    con.execute("SET memory_limit = '64GB'")
    con.execute("SET threads = 8")
    con.execute("SET preserve_insertion_order = false")
    con.register("clinvar_symbol", clinvar)

    g = sql_path(GENES)
    x = sql_path(CROSSWALK)
    c = sql_path(CONSTRAINT)
    e = sql_path(EXPRESSION)
    loci = sql_path(LOCI)
    coverage = sql_path(COVERAGE)
    scored = sql_path(SCORED)
    egenes = sql_path(EGENES)
    sgenes = sql_path(SGENES)

    dup_constraint = con.execute(
        f"""
        SELECT count(*) FROM (
          SELECT bare_ensg
          FROM read_parquet('{c}')
          WHERE bare_ensg IS NOT NULL
          GROUP BY bare_ensg
          HAVING count(*) > 1
        )
        """
    ).fetchone()[0]
    if dup_constraint:
        raise SystemExit(f"gene_constraint bare_ensg is not unique: {dup_constraint}")
    dup_expr = con.execute(
        f"""
        SELECT count(*) FROM (
          SELECT ensg FROM read_parquet('{e}')
          WHERE ensg IS NOT NULL
          GROUP BY ensg
          HAVING count(*) > 1
        )
        """
    ).fetchone()[0]
    if dup_expr:
        raise SystemExit(f"gene_expression ensg is not unique: {dup_expr}")

    print("gwas overlap", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE gwas_overlap AS
        SELECT
          g.ensg,
          count(DISTINCT concat_ws(
            '|',
            coalesce(l.trait_uri, ''),
            l.chrom,
            cast(l.lead_pos AS VARCHAR),
            coalesce(l.lead_rsid, '')
          )) AS gwas_loci_by_overlap,
          count(DISTINCT l.trait_uri) AS gwas_traits_by_overlap
        FROM read_parquet('{g}') g
        JOIN read_parquet('{loci}') l
          ON g.chrom = l.chrom
         AND (
              (l.start <= g."end" AND l."end" >= g.start)
              OR (l.lead_pos BETWEEN g.start AND g."end")
         )
        GROUP BY g.ensg
        """
    )
    print("gwas mapped symbol", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE gwas_mapped AS
        SELECT
          g.ensg,
          count(DISTINCT concat_ws(
            '|',
            coalesce(t.trait_uri, ''),
            coalesce(t.chrom, ''),
            cast(t.lead_pos AS VARCHAR),
            coalesce(t.lead_rsid, '')
          )) AS gwas_loci_by_mapped_gene,
          count(DISTINCT t.trait_uri) AS gwas_traits_by_mapped_gene
        FROM read_parquet('{g}') g
        JOIN (
          SELECT trait_uri, chrom, lead_pos, lead_rsid, unnest(mapped_genes.tokens) AS token
          FROM read_parquet('{loci}')
        ) t
          ON t.token = g.symbol
        GROUP BY g.ensg
        """
    )
    print("mavedb span", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE mave_span AS
        SELECT
          g.ensg,
          count(DISTINCT v.variant_key) AS mavedb_variant_keys_in_span,
          count(DISTINCT v.score_set_urn) AS mavedb_score_sets_in_span
        FROM read_parquet('{g}') g
        JOIN (
          SELECT variant_key, score_set_urn, chrom, pos
          FROM read_parquet('{scored}')
          WHERE variant_key IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL
        ) v
          ON v.chrom = g.chrom
         AND v.pos BETWEEN g.start AND g."end"
        GROUP BY g.ensg
        """
    )

    clinvar_select = ",\n          ".join(
        f"coalesce(cv.{name}, 0)::BIGINT AS {name}" for name in clinvar_columns()
    )
    order = """
        CASE i.chrom
          WHEN 'chrX' THEN 23
          WHEN 'chrY' THEN 24
          WHEN 'chrM' THEN 25
          ELSE try_cast(regexp_extract(i.chrom, '[0-9]+') AS INTEGER)
        END,
        i.start,
        i.ensg
    """
    select = f"""
        SELECT
          i.ensg,
          i.ensg_version,
          i.symbol,
          i.chrom,
          i.start,
          i."end" AS "end",
          i.strand,
          i.biotype,
          x.mane_enst,
          x.mane_enst_gencode,
          x.mane_enst_constraint,
          x.symbol_n_gencode,
          k.loeuf,
          CAST(k.loeuf_decile AS INTEGER) AS loeuf_decile,
          k."lof.pLI" AS "lof.pLI",
          k."mis.z_score" AS "mis.z_score",
          k."syn.z_score" AS "syn.z_score",
          k.constraint_flags,
          expr.max_tissue,
          expr.max_median_tpm,
          expr.n_tissues_tpm_ge_1,
          expr.n_tissues_tpm_ge_10,
          expr.tau,
          {clinvar_select},
          coalesce(go.gwas_loci_by_overlap, 0)::BIGINT AS gwas_loci_by_overlap,
          coalesce(go.gwas_traits_by_overlap, 0)::BIGINT AS gwas_traits_by_overlap,
          coalesce(gm.gwas_loci_by_mapped_gene, 0)::BIGINT AS gwas_loci_by_mapped_gene,
          coalesce(gm.gwas_traits_by_mapped_gene, 0)::BIGINT AS gwas_traits_by_mapped_gene,
          coalesce(mc.mavedb_experiments_by_symbol, 0)::BIGINT AS mavedb_experiments_by_symbol,
          coalesce(mc.mavedb_mapped_variants_by_symbol, 0)::BIGINT AS mavedb_mapped_variants_by_symbol,
          coalesce(ms.mavedb_variant_keys_in_span, 0)::BIGINT AS mavedb_variant_keys_in_span,
          coalesce(ms.mavedb_score_sets_in_span, 0)::BIGINT AS mavedb_score_sets_in_span,
          eq.eqtl_tissues,
          sq.sqtl_tissues
        FROM read_parquet('{g}') i
        LEFT JOIN read_parquet('{x}') x ON i.ensg = x.ensg
        LEFT JOIN (
          SELECT
            bare_ensg AS ensg,
            loeuf,
            loeuf_decile,
            "lof.pLI",
            "mis.z_score",
            "syn.z_score",
            constraint_flags
          FROM read_parquet('{c}')
          WHERE bare_ensg IS NOT NULL
        ) k ON i.ensg = k.ensg
        LEFT JOIN (
          SELECT ensg, max_tissue, max_median_tpm, n_tissues_tpm_ge_1, n_tissues_tpm_ge_10, tau
          FROM read_parquet('{e}')
          WHERE ensg IS NOT NULL
        ) expr ON i.ensg = expr.ensg
        LEFT JOIN clinvar_symbol cv ON i.symbol = cv.symbol
        LEFT JOIN gwas_overlap go ON i.ensg = go.ensg
        LEFT JOIN gwas_mapped gm ON i.ensg = gm.ensg
        LEFT JOIN (
          SELECT target_gene AS symbol,
                 sum(n_score_sets)::BIGINT AS mavedb_experiments_by_symbol,
                 sum(n_mapped_variants)::BIGINT AS mavedb_mapped_variants_by_symbol
          FROM read_parquet('{coverage}')
          GROUP BY target_gene
        ) mc ON i.symbol = mc.symbol
        LEFT JOIN mave_span ms ON i.ensg = ms.ensg
        LEFT JOIN (
          SELECT ensg, count(DISTINCT tissue) FILTER (WHERE qval_le_0_05)::INTEGER AS eqtl_tissues
          FROM read_parquet('{egenes}')
          WHERE ensg IS NOT NULL
          GROUP BY ensg
        ) eq ON i.ensg = eq.ensg
        LEFT JOIN (
          SELECT ensg, count(DISTINCT tissue) FILTER (WHERE qval_le_0_05)::INTEGER AS sqtl_tissues
          FROM read_parquet('{sgenes}')
          WHERE ensg IS NOT NULL
          GROUP BY ensg
        ) sq ON i.ensg = sq.ensg
        ORDER BY {order}
    """
    print("write gene mart", flush=True)
    con.execute(
        f"COPY ({select}) TO '{sql_path(MART)}' (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
    con.execute(
        f"""
        COPY (
          SELECT * FROM read_parquet('{sql_path(MART)}')
          WHERE biotype = 'protein_coding'
        ) TO '{sql_path(CODING)}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )
    rows, coding, other = con.execute(
        f"""
        SELECT
          count(*),
          count(*) FILTER (WHERE biotype = 'protein_coding'),
          count(*) FILTER (WHERE biotype = 'protein_coding')
            - (SELECT count(*) FROM read_parquet('{sql_path(CODING)}'))
        FROM read_parquet('{sql_path(MART)}')
        """
    ).fetchone()
    coding_file = con.execute(f"SELECT count(*) FROM read_parquet('{sql_path(CODING)}')").fetchone()[0]
    print(f"QC gene_mart_rows={rows}", flush=True)
    print(f"QC protein_coding_rows={coding}", flush=True)
    print(f"QC protein_coding_file_rows={coding_file}", flush=True)
    if rows != EXPECTED_GENES or coding != EXPECTED_CODING or coding_file != EXPECTED_CODING or other:
        raise SystemExit(f"row count mismatch mart {rows} coding {coding} file {coding_file} delta {other}")

    stats = con.execute(
        f"""
        SELECT
          count(*) FILTER (WHERE loeuf IS NULL),
          count(*) FILTER (WHERE biotype = 'protein_coding' AND loeuf IS NULL),
          count(*) FILTER (WHERE max_tissue IS NULL AND max_median_tpm IS NULL AND tau IS NULL),
          count(*) FILTER (WHERE eqtl_tissues IS NULL),
          count(*) FILTER (WHERE sqtl_tissues IS NULL),
          count(*) FILTER (WHERE symbol_n_gencode > 1),
          count(*) FILTER (WHERE symbol_n_gencode > 1 AND clinvar_n > 0),
          count(*) FILTER (WHERE symbol_n_gencode > 1 AND mavedb_experiments_by_symbol > 0),
          count(*) FILTER (WHERE mavedb_experiments_by_symbol > 0),
          count(*) FILTER (WHERE mavedb_variant_keys_in_span > 0),
          count(*) FILTER (WHERE gwas_loci_by_overlap > 0),
          count(*) FILTER (WHERE gwas_loci_by_mapped_gene > 0),
          count(*) FILTER (WHERE clinvar_n > 0),
          count(*) FILTER (WHERE mane_enst IS NOT NULL),
          count(*) FILTER (WHERE mane_enst_constraint IS NOT NULL),
          sum(clinvar_germline_pathogenic_n),
          sum(clinvar_likely_pathogenic_n),
          sum(clinvar_drug_response_n),
          sum(clinvar_conflicting_n),
          sum(clinvar_vus_n)
        FROM read_parquet('{sql_path(MART)}')
        """
    ).fetchone()
    labels = (
        "constraint_miss",
        "protein_coding_constraint_miss",
        "expression_miss",
        "eqtl_miss",
        "sqtl_miss",
        "symbol_collision_genes",
        "clinvar_symbol_collision_genes",
        "mavedb_symbol_collision_genes",
        "mavedb_symbol_genes",
        "mavedb_span_genes",
        "gwas_overlap_genes",
        "gwas_mapped_genes",
        "clinvar_genes",
        "mane_enst_present",
        "mane_enst_constraint_present",
        "sum_germline_pathogenic",
        "sum_likely_pathogenic",
        "sum_drug_response",
        "sum_conflicting",
        "sum_vus",
    )
    for label, value in zip(labels, stats):
        print(f"QC {label}={value}", flush=True)

    expr_orphan = con.execute(
        f"""
        SELECT count(*)
        FROM read_parquet('{e}') e
        ANTI JOIN read_parquet('{g}') g ON e.ensg = g.ensg
        WHERE e.ensg IS NOT NULL
        """
    ).fetchone()[0]
    expr_null_id = con.execute(
        f"SELECT count(*) FROM read_parquet('{e}') WHERE ensg IS NULL"
    ).fetchone()[0]
    print(f"QC expression_ensg_not_in_gencode={expr_orphan}", flush=True)
    print(f"QC expression_null_ensg={expr_null_id}", flush=True)
    coord_mismatch = con.execute(
        f"""
        SELECT count(*)
        FROM read_parquet('{g}') i
        JOIN read_parquet('{x}') x ON i.ensg = x.ensg
        WHERE i.chrom IS DISTINCT FROM x.chrom
           OR i.start IS DISTINCT FROM x.start
           OR i."end" IS DISTINCT FROM x."end"
           OR i.symbol IS DISTINCT FROM x.symbol
           OR i.biotype IS DISTINCT FROM x.biotype
        """
    ).fetchone()[0]
    print(f"QC identity_coord_mismatch={coord_mismatch}", flush=True)

    for symbol, transcript, loeuf, decile in TRUTH_GENES:
        found = con.execute(
            f"""
            SELECT count(*),
                   any_value(mane_enst),
                   any_value(loeuf),
                   any_value(loeuf_decile),
                   any_value(mane_enst_constraint)
            FROM read_parquet('{sql_path(MART)}')
            WHERE symbol = ?
            """,
            [symbol],
        ).fetchone()
        if found[0] != 1:
            raise SystemExit(f"{symbol} rows {found[0]}")
        if found[1] != transcript:
            raise SystemExit(f"{symbol} mane_enst {found[1]} != {transcript}")
        if found[2] is None or abs(float(found[2]) - loeuf) > 0.001:
            raise SystemExit(f"{symbol} loeuf {found[2]} != {loeuf}")
        if found[3] is None or int(found[3]) != decile:
            raise SystemExit(f"{symbol} decile {found[3]} != {decile}")
        print(
            f"QC truth {symbol} mane_enst={found[1]} loeuf={found[2]} decile={found[3]} "
            f"mane_enst_constraint={found[4]}",
            flush=True,
        )
    con.close()
    print(f"wrote {MART}", flush=True)
    print(f"wrote {CODING}", flush=True)


if __name__ == "__main__":
    build()
