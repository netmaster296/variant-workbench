# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Evidence ledger and the dark-genome base count.

Six classes stay separate: clinvar, gwas, eqtl, sqtl, mavedb, constraint.
Significant-pair tables are not read. Genome and exome allele counts are not
used here. Chromosome lengths come from a gnomAD genome contig header.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf_stream  # noqa: E402
from coords import bare_ensg, canon_chrom  # noqa: E402

DEDUP = ROOT / "warehouse" / "clinvar" / "allele_dedup.parquet"
LOCI = ROOT / "warehouse" / "gwas" / "loci.parquet"
ASSOC = ROOT / "warehouse" / "gwas" / "associations.parquet"
EGENES = ROOT / "warehouse" / "gtex" / "egenes.parquet"
SGENES = ROOT / "warehouse" / "gtex" / "sgenes.parquet"
SCORED = ROOT / "warehouse" / "mavedb" / "mapped_variants_scored.parquet"
CONSTRAINT = ROOT / "warehouse" / "constraint" / "gene_constraint.parquet"
GENES = ROOT / "warehouse" / "intervals" / "genes.parquet"
CDS = ROOT / "warehouse" / "intervals" / "cds.parquet"
KEYS = ROOT / "warehouse" / "variants" / "tier_a_keys.tsv"
LEDGER = ROOT / "warehouse" / "ledger" / "evidence.parquet"
REPORT = ROOT / "warehouse" / "reports" / "dark_genome.md"
LENGTH_VCF = ROOT / "raw" / "gnomad" / "genomes" / "gnomad.genomes.v4.1.sites.chr21.vcf.bgz"
TMP = ROOT / "warehouse" / "tmp" / "ledger"

SCHEMA = pa.schema(
    [
        ("evidence_id", pa.string()),
        ("subject_type", pa.string()),
        ("subject_id", pa.string()),
        ("object_type", pa.string()),
        ("object_id", pa.string()),
        ("evidence_class", pa.string()),
        ("strength_value", pa.float64()),
        ("strength_unit", pa.string()),
        ("source_row", pa.string()),
        ("ancestry_note", pa.string()),
        ("n_note", pa.string()),
    ]
)
MISSING = {"", ".", "-", "na", "n/a", "none", "not provided"}
CLASSES = ("clinvar", "gwas", "eqtl", "sqtl", "mavedb", "constraint")


def sql_path(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def evidence_id(natural: str) -> str:
    return hashlib.sha1(natural.encode("utf-8")).hexdigest()


def blank(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in MISSING:
        return None
    return text


class Ledger:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.buf = {field.name: [] for field in SCHEMA}
        self.writer = None
        self.n = 0
        self.by_class = {name: 0 for name in CLASSES}

    def add(self, row: dict) -> None:
        self.by_class[row["evidence_class"]] += 1
        for field in SCHEMA:
            self.buf[field.name].append(row.get(field.name))
        if len(self.buf["evidence_id"]) >= 100_000:
            self.flush()

    def flush(self) -> None:
        if not self.buf["evidence_id"]:
            return
        table = pa.table(self.buf, schema=SCHEMA)
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.path, SCHEMA, compression="zstd")
        self.writer.write_table(table)
        self.n += table.num_rows
        for column in self.buf.values():
            column.clear()
        print(f"  ledger rows {self.n}", flush=True)

    def close(self) -> None:
        self.flush()
        if self.writer is not None:
            self.writer.close()


def unique_symbols() -> dict[str, str]:
    table = pq.read_table(GENES, columns=["ensg", "symbol"])
    grouped: dict[str, list[str]] = {}
    for ensg, symbol in zip(table.column("ensg").to_pylist(), table.column("symbol").to_pylist()):
        if not symbol:
            continue
        if bare_ensg(ensg) != ensg:
            raise SystemExit(f"gene id is not bare: {ensg}")
        grouped.setdefault(symbol, []).append(ensg)
    return {symbol: ensgs[0] for symbol, ensgs in grouped.items() if len(ensgs) == 1}


def add_clinvar(ledger: Ledger, symbols: dict[str, str]) -> dict:
    stats = {"rows": 0, "variant_subject": 0, "gene_subject": 0, "allele_id_subject": 0}
    columns = [
        "dedup_key",
        "variant_key",
        "VariationID",
        "GeneSymbol",
        "review_stars",
        "PhenotypeIDS",
    ]
    for batch in pq.ParquetFile(DEDUP).iter_batches(batch_size=100_000, columns=columns):
        data = batch.to_pydict()
        for dedup_key, key, variation, symbol, stars, pheno in zip(
            data["dedup_key"],
            data["variant_key"],
            data["VariationID"],
            data["GeneSymbol"],
            data["review_stars"],
            data["PhenotypeIDS"],
        ):
            stats["rows"] += 1
            if key:
                subject_type = "variant"
                subject_id = key
                stats["variant_subject"] += 1
            else:
                token = symbol.strip() if isinstance(symbol, str) else ""
                ensg = symbols.get(token)
                if ensg:
                    subject_type = "gene"
                    subject_id = ensg
                    stats["gene_subject"] += 1
                else:
                    subject_type = "variant"
                    subject_id = dedup_key
                    stats["allele_id_subject"] += 1
            phenotype = blank(pheno)
            strength = None if stars is None else float(stars)
            natural = f"clinvar|{dedup_key}|{phenotype or ''}"
            ledger.add(
                {
                    "evidence_id": evidence_id(natural),
                    "subject_type": subject_type,
                    "subject_id": subject_id,
                    "object_type": "phenotype",
                    "object_id": phenotype,
                    "evidence_class": "clinvar",
                    "strength_value": strength,
                    "strength_unit": "stars",
                    "source_row": f"warehouse/clinvar/allele_dedup.parquet#VariationID={variation}",
                    "ancestry_note": None,
                    "n_note": None,
                }
            )
    return stats


def add_gwas(ledger: Ledger, con: duckdb.DuckDBPyConnection) -> dict:
    print("gwas leads", flush=True)
    con.execute(
        f"""
        CREATE TEMP TABLE lead_sample AS
        SELECT trait_uri, chrom, lead_pos, lead_rsid, study, initial_n, rep_n
        FROM (
          SELECT
            l.trait_uri,
            l.chrom,
            l.lead_pos,
            l.lead_rsid,
            a."STUDY ACCESSION" AS study,
            a."INITIAL SAMPLE SIZE" AS initial_n,
            a."REPLICATION SAMPLE SIZE" AS rep_n,
            row_number() OVER (
              PARTITION BY l.trait_uri, l.chrom, l.lead_pos, l.lead_rsid
              ORDER BY a.p_value_float NULLS LAST, a."STUDY ACCESSION"
            ) AS rn
          FROM read_parquet('{sql_path(LOCI)}') l
          LEFT JOIN read_parquet('{sql_path(ASSOC)}') a
            ON a.chrom = l.chrom
           AND a.pos = l.lead_pos
           AND a.rsid IS NOT DISTINCT FROM l.lead_rsid
           AND list_contains(
             list_transform(
               string_split(replace(a."MAPPED_TRAIT_URI", ', ', ','), ','),
               lambda x: trim(x)
             ),
             l.trait_uri
           )
        )
        WHERE rn = 1
        """
    )
    clin = con.execute(
        f"""
        SELECT c.rsid, c.variant_key, c.chrom, c.pos
        FROM read_parquet('{sql_path(ROOT / "warehouse" / "clinvar" / "allele_dedup.parquet")}') c
        SEMI JOIN (
          SELECT DISTINCT lead_rsid AS rsid
          FROM read_parquet('{sql_path(LOCI)}')
          WHERE lead_rsid IS NOT NULL
        ) l ON c.rsid = l.rsid
        WHERE c.variant_key IS NOT NULL AND c.rsid IS NOT NULL
        """
    ).fetchall()
    by_rsid: dict[str, list[tuple]] = {}
    for rsid, key, chrom, pos in clin:
        by_rsid.setdefault(rsid, []).append((key, chrom, int(pos)))
    rows = con.execute(
        f"""
        SELECT
          l.trait_uri, l.chrom, l.lead_pos, l.lead_rsid, l.lead_p,
          s.study, s.initial_n, s.rep_n
        FROM read_parquet('{sql_path(LOCI)}') l
        LEFT JOIN lead_sample s
          ON s.trait_uri IS NOT DISTINCT FROM l.trait_uri
         AND s.chrom IS NOT DISTINCT FROM l.chrom
         AND s.lead_pos IS NOT DISTINCT FROM l.lead_pos
         AND s.lead_rsid IS NOT DISTINCT FROM l.lead_rsid
        """
    ).fetchall()
    stats = {
        "loci": 0,
        "rows": 0,
        "variant_keys": 0,
        "locus_keys": 0,
        "multi_alt": 0,
        "rsid_match_other_position": 0,
        "with_sample_size": 0,
    }
    for trait, chrom, pos, rsid, lead_p, study, initial, rep in rows:
        stats["loci"] += 1
        if chrom is None or pos is None:
            continue
        pos = int(pos)
        matches = by_rsid.get(rsid or "", [])
        same = [item for item in matches if item[1] == chrom and item[2] == pos]
        if same:
            subjects = [item[0] for item in same]
            stats["variant_keys"] += 1
            if len(same) > 1:
                stats["multi_alt"] += 1
        else:
            subjects = [f"{chrom}:{pos}"]
            stats["locus_keys"] += 1
            if matches:
                stats["rsid_match_other_position"] += 1
        notes = []
        if blank(initial):
            notes.append(blank(initial))
        if blank(rep):
            notes.append("replication: " + blank(rep))
        n_note = " | ".join(notes) if notes else None
        if n_note:
            stats["with_sample_size"] += 1
        source = (
            "warehouse/gwas/loci.parquet#"
            f"{chrom}:{pos}:{rsid or ''}:{trait or ''}"
        )
        if study:
            source += f"#study={study}"
        for subject in subjects:
            natural = f"gwas|{trait or ''}|{chrom}|{pos}|{rsid or ''}|{subject}"
            stats["rows"] += 1
            ledger.add(
                {
                    "evidence_id": evidence_id(natural),
                    "subject_type": "variant",
                    "subject_id": subject,
                    "object_type": "trait",
                    "object_id": trait,
                    "evidence_class": "gwas",
                    "strength_value": None if lead_p is None else float(lead_p),
                    "strength_unit": "p_value",
                    "source_row": source,
                    "ancestry_note": None,
                    "n_note": n_note,
                }
            )
    return stats


def add_qtl(ledger: Ledger, con: duckdb.DuckDBPyConnection, kind: str) -> dict:
    if kind == "eqtl":
        path = EGENES
        object_type = "gene_expression"
        evidence_class = "eqtl"
        phenotype = "NULL"
        count_col = "NULL"
    else:
        path = SGENES
        object_type = "splicing"
        evidence_class = "sqtl"
        phenotype = "phenotype_id"
        count_col = "ma_count"
    print(kind, flush=True)
    cursor = con.execute(
        f"""
        SELECT tissue, ensg, variant_key, pval_nominal, variant_id, {phenotype}
        FROM read_parquet('{sql_path(path)}')
        WHERE qval_le_0_05 AND variant_key IS NOT NULL AND ensg IS NOT NULL
        """
    )
    stats = {"rows": 0}
    while True:
        batch = cursor.fetchmany(100_000)
        if not batch:
            break
        for tissue, ensg, key, pvalue, variant_id, phenotype_id in batch:
            if bare_ensg(ensg) != ensg:
                ensg = bare_ensg(ensg) or ensg
            if kind == "sqtl":
                natural = f"sqtl|{phenotype_id}|{tissue}|{key}"
                source = (
                    "warehouse/gtex/sgenes.parquet#"
                    f"{tissue}:{phenotype_id}:{variant_id}"
                )
            else:
                natural = f"eqtl|{ensg}|{tissue}|{key}"
                source = f"warehouse/gtex/egenes.parquet#{tissue}:{variant_id}"
            stats["rows"] += 1
            ledger.add(
                {
                    "evidence_id": evidence_id(natural),
                    "subject_type": "variant",
                    "subject_id": key,
                    "object_type": object_type,
                    "object_id": f"{ensg}|{tissue}",
                    "evidence_class": evidence_class,
                    "strength_value": None if pvalue is None else float(pvalue),
                    "strength_unit": "p_value",
                    "source_row": source,
                    "ancestry_note": None,
                    "n_note": None,
                }
            )
    if count_col:
        pass
    return stats


def add_mavedb(ledger: Ledger, symbols: dict[str, str]) -> dict:
    print("mavedb", flush=True)
    stats = {
        "rows_seen": 0,
        "written": 0,
        "no_score": 0,
        "variant_subject": 0,
        "gene_subject": 0,
        "no_subject": 0,
    }
    columns = ["urn", "score_set_urn", "target_gene", "variant_key", "member_index", "score", "score_text"]
    for batch in pq.ParquetFile(SCORED).iter_batches(batch_size=100_000, columns=columns):
        data = batch.to_pydict()
        for urn, score_set, gene, key, member, score, score_text in zip(
            data["urn"],
            data["score_set_urn"],
            data["target_gene"],
            data["variant_key"],
            data["member_index"],
            data["score"],
            data["score_text"],
        ):
            stats["rows_seen"] += 1
            if score is None and not blank(score_text):
                stats["no_score"] += 1
                continue
            if key:
                subject_type = "variant"
                subject_id = key
                stats["variant_subject"] += 1
            else:
                ensg = symbols.get(gene) if isinstance(gene, str) else None
                if not ensg:
                    stats["no_subject"] += 1
                    continue
                subject_type = "gene"
                subject_id = ensg
                stats["gene_subject"] += 1
            member_text = "" if member is None else str(member)
            natural = f"mavedb|{score_set}|{urn}|{member_text}"
            stats["written"] += 1
            ledger.add(
                {
                    "evidence_id": evidence_id(natural),
                    "subject_type": subject_type,
                    "subject_id": subject_id,
                    "object_type": "assay",
                    "object_id": score_set,
                    "evidence_class": "mavedb",
                    "strength_value": None if score is None else float(score),
                    "strength_unit": "score",
                    "source_row": f"warehouse/mavedb/mapped_variants_scored.parquet#{urn}#{member_text}",
                    "ancestry_note": None,
                    "n_note": None,
                }
            )
    return stats


def add_constraint(ledger: Ledger, con: duckdb.DuckDBPyConnection) -> dict:
    rows = con.execute(
        f"""
        SELECT bare_ensg, gene_id, transcript, loeuf
        FROM read_parquet('{sql_path(CONSTRAINT)}')
        WHERE bare_ensg IS NOT NULL
        """
    ).fetchall()
    dropped = con.execute(
        f"SELECT count(*) FROM read_parquet('{sql_path(CONSTRAINT)}') WHERE bare_ensg IS NULL"
    ).fetchone()[0]
    for ensg, gene_id, transcript, loeuf in rows:
        if bare_ensg(ensg) != ensg:
            raise SystemExit(f"constraint id is not bare: {ensg}")
        natural = f"constraint|{ensg}|loeuf"
        ledger.add(
            {
                "evidence_id": evidence_id(natural),
                "subject_type": "gene",
                "subject_id": ensg,
                "object_type": "constraint",
                "object_id": "loeuf",
                "evidence_class": "constraint",
                "strength_value": None if loeuf is None else float(loeuf),
                "strength_unit": "loeuf",
                "source_row": f"warehouse/constraint/gene_constraint.parquet#{gene_id}:{transcript}",
                "ancestry_note": None,
                "n_note": None,
            }
        )
    return {"rows": len(rows), "null_bare_ensg": dropped}


def chrom_rank(name: str) -> int:
    if name == "chrX":
        return 23
    if name == "chrY":
        return 24
    if name == "chrM":
        return 25
    if name.startswith("chr") and name[3:].isdigit():
        return int(name[3:])
    return 99


def load_lengths() -> dict[str, int]:
    if canon_chrom("21") != "chr21":
        raise SystemExit("coords.canon_chrom failed")
    lengths = vcf_stream.contig_lengths(LENGTH_VCF)
    wanted = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]
    missing = [name for name in wanted if name not in lengths]
    if missing:
        raise SystemExit(f"contig header missing {missing}")
    print(f"QC contig_lengths={len(lengths)} source={LENGTH_VCF.name}", flush=True)
    return {name: lengths[name] for name in wanted}


def covered(intervals: list[tuple[int, int]], length: int) -> int:
    clipped = []
    for start, end in intervals:
        if start > end:
            start, end = end, start
        if end < 1 or start > length:
            continue
        clipped.append((max(1, start), min(length, end)))
    if not clipped:
        return 0
    clipped.sort()
    total = 0
    cs, ce = clipped[0]
    for start, end in clipped[1:]:
        if start <= ce + 1:
            if end > ce:
                ce = end
        else:
            total += ce - cs + 1
            cs, ce = start, end
    total += ce - cs + 1
    return total


def dark_genome(lengths: dict[str, int], con: duckdb.DuckDBPyConnection) -> None:
    print("dark genome", flush=True)
    cds_rows = con.execute(
        f"""
        SELECT c.chrom, c.start, c."end"
        FROM read_parquet('{sql_path(CDS)}') c
        JOIN read_parquet('{sql_path(GENES)}') g ON c.ensg = g.ensg
        WHERE g.biotype = 'protein_coding' AND c.feature = 'CDS'
        """
    ).fetchall()
    gene_rows = con.execute(
        f"""SELECT chrom, start, "end" FROM read_parquet('{sql_path(GENES)}')"""
    ).fetchall()
    loci_rows = con.execute(
        f"""
        SELECT chrom, start, "end"
        FROM read_parquet('{sql_path(LOCI)}')
        WHERE chrom IS NOT NULL AND start IS NOT NULL AND "end" IS NOT NULL
        """
    ).fetchall()
    singletons = con.execute(
        f"""
        SELECT count(*) FROM read_parquet('{sql_path(LOCI)}')
        WHERE start = "end"
        """
    ).fetchone()[0]
    tier: dict[str, list[tuple[int, int]]] = {name: [] for name in lengths}
    tier_outside = 0
    with KEYS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row["chrom"]
            pos = int(row["pos"])
            if chrom not in lengths:
                tier_outside += 1
                continue
            tier[chrom].append((pos - 1000, pos + 1000))
    buckets = {
        "cds": {name: [] for name in lengths},
        "gene": {name: [] for name in lengths},
        "gwas": {name: [] for name in lengths},
    }
    excluded = {"cds": 0, "gene": 0, "gwas": 0}
    for chrom, start, end in cds_rows:
        if chrom not in lengths:
            excluded["cds"] += 1
            continue
        buckets["cds"][chrom].append((int(start), int(end)))
    for chrom, start, end in gene_rows:
        if chrom not in lengths:
            excluded["gene"] += 1
            continue
        buckets["gene"][chrom].append((int(start), int(end)))
    for chrom, start, end in loci_rows:
        if chrom not in lengths:
            excluded["gwas"] += 1
            continue
        buckets["gwas"][chrom].append((int(start), int(end)))

    lines = []
    totals = {name: 0 for name in ("length", "cds", "gene", "tier", "gwas", "any", "none")}
    for chrom in sorted(lengths, key=chrom_rank):
        length = lengths[chrom]
        cds_n = covered(buckets["cds"][chrom], length)
        gene_n = covered(buckets["gene"][chrom], length)
        tier_n = covered(tier[chrom], length)
        gwas_n = covered(buckets["gwas"][chrom], length)
        any_n = covered(
            buckets["cds"][chrom] + buckets["gene"][chrom] + tier[chrom] + buckets["gwas"][chrom],
            length,
        )
        none_n = length - any_n
        if any_n > length or none_n < 0:
            raise SystemExit(f"{chrom} coverage {any_n} exceeds length {length}")
        totals["length"] += length
        totals["cds"] += cds_n
        totals["gene"] += gene_n
        totals["tier"] += tier_n
        totals["gwas"] += gwas_n
        totals["any"] += any_n
        totals["none"] += none_n
        lines.append((chrom, length, cds_n, gene_n, tier_n, gwas_n, any_n, none_n))
    lines.append(
        (
            "total",
            totals["length"],
            totals["cds"],
            totals["gene"],
            totals["tier"],
            totals["gwas"],
            totals["any"],
            totals["none"],
        )
    )

    def fraction(part: int) -> str:
        return f"{part / totals['length']:.6f}"

    report = [
        "# Dark genome",
        "",
        "This is a count of bases, not a claim that the uncounted bases have no function.",
        "",
        f"Contig lengths were read from the `##contig` header of `{LENGTH_VCF}` by `scripts/vcf_stream.py`.",
        "The table is not hard-coded. The header lists chr1–chr22, chrX, and chrY. It has no chrM contig.",
        f"Denominator: {totals['length']} bases.",
        "",
        "Each category is the union of its own intervals, so overlap inside a category is counted once.",
        "A base that falls in two categories is counted in both. The categories are not an exclusive partition.",
        "The remainder is the bases that are in none of the four categories.",
        "",
        "Definitions:",
        "",
        "- Protein-coding CDS: `warehouse/intervals/cds.parquet` rows with `feature = CDS` on a gene whose biotype is `protein_coding`. `stop_codon` rows are not included.",
        "- Gene span: `start`–`end` of every GENCODE gene on these contigs, any biotype.",
        "- Tier A neighborhood: inclusive `[pos - 1000, pos + 1000]` around each key in `warehouse/variants/tier_a_keys.tsv`, clipped to the contig.",
        "- GWAS clump: `start`–`end` in `warehouse/gwas/loci.parquet`. Those columns are the minimum and maximum member positions. They are not expanded to the ±1,000,000 bp search radius used to build the clump.",
        "",
        f"Loci whose start equals end (single association): {singletons}.",
        f"CDS intervals outside the contig set, excluded: {excluded['cds']}.",
        f"Gene spans outside the contig set, excluded: {excluded['gene']}.",
        f"GWAS intervals outside the contig set, excluded: {excluded['gwas']}.",
        f"Tier A keys outside the contig set, excluded: {tier_outside}.",
        "",
        "| chrom | length | protein_coding_cds | gene_span | tier_a_1kb | gwas_clump | in_any | in_none |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for chrom, length, cds_n, gene_n, tier_n, gwas_n, any_n, none_n in lines:
        report.append(
            f"| {chrom} | {length} | {cds_n} | {gene_n} | {tier_n} | {gwas_n} | {any_n} | {none_n} |"
        )
    report.extend(
        [
            "",
            "| set | bases | fraction of denominator |",
            "| --- | ---: | ---: |",
            f"| protein_coding_cds | {totals['cds']} | {fraction(totals['cds'])} |",
            f"| gene_span | {totals['gene']} | {fraction(totals['gene'])} |",
            f"| tier_a_1kb | {totals['tier']} | {fraction(totals['tier'])} |",
            f"| gwas_clump | {totals['gwas']} | {fraction(totals['gwas'])} |",
            f"| in_any | {totals['any']} | {fraction(totals['any'])} |",
            f"| in_none | {totals['none']} | {fraction(totals['none'])} |",
            "",
            "The four category fractions are not exclusive, so they do not sum to the in_any fraction.",
            "in_none is the exclusive remainder.",
            "",
        ]
    )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(report), encoding="utf-8", newline="\n")
    print(f"wrote {REPORT}", flush=True)
    for name in ("length", "cds", "gene", "tier", "gwas", "any", "none"):
        print(f"QC dark_{name}={totals[name]}", flush=True)


def main() -> None:
    if not KEYS.exists():
        raise SystemExit(f"missing {KEYS}")
    TMP.mkdir(parents=True, exist_ok=True)
    lengths = load_lengths()
    symbols = unique_symbols()
    ledger = Ledger(LEDGER)
    con = duckdb.connect()
    con.execute(f"SET temp_directory = '{sql_path(TMP)}'")
    con.execute("SET memory_limit = '64GB'")
    con.execute("SET threads = 8")
    con.execute("SET preserve_insertion_order = false")
    try:
        clin = add_clinvar(ledger, symbols)
        print("QC clinvar " + " ".join(f"{key}={value}" for key, value in clin.items()), flush=True)
        gwas = add_gwas(ledger, con)
        print("QC gwas " + " ".join(f"{key}={value}" for key, value in gwas.items()), flush=True)
        eqtl = add_qtl(ledger, con, "eqtl")
        sqtl = add_qtl(ledger, con, "sqtl")
        print(f"QC eqtl_rows={eqtl['rows']}", flush=True)
        print(f"QC sqtl_rows={sqtl['rows']}", flush=True)
        mave = add_mavedb(ledger, symbols)
        print("QC mavedb " + " ".join(f"{key}={value}" for key, value in mave.items()), flush=True)
        constraint = add_constraint(ledger, con)
        print("QC constraint " + " ".join(f"{key}={value}" for key, value in constraint.items()), flush=True)
    finally:
        ledger.close()
        con.close()
    for name in CLASSES:
        print(f"QC ledger_{name}={ledger.by_class[name]}", flush=True)
        if ledger.by_class[name] <= 0:
            raise SystemExit(f"ledger has no {name} rows")
    con = duckdb.connect()
    dup, n_rows, n_id = con.execute(
        f"""
        SELECT count(*) - count(DISTINCT evidence_id), count(*), count(DISTINCT evidence_id)
        FROM read_parquet('{sql_path(LEDGER)}')
        """
    ).fetchone()
    print(f"QC ledger_rows={n_rows} distinct_id={n_id} duplicate_id={dup}", flush=True)
    if dup:
        raise SystemExit(f"duplicate evidence_id rows: {dup}")
    classes = con.execute(
        f"""
        SELECT evidence_class, count(*)
        FROM read_parquet('{sql_path(LEDGER)}')
        GROUP BY 1
        ORDER BY 1
        """
    ).fetchall()
    found = {name for name, _count in classes}
    if found != set(CLASSES):
        raise SystemExit(f"ledger classes {sorted(found)}")
    for name, count in classes:
        print(f"QC ledger_file_{name}={count}", flush=True)
    con.close()
    con = duckdb.connect()
    con.execute(f"SET temp_directory = '{sql_path(TMP)}'")
    con.execute("SET threads = 8")
    dark_genome(lengths, con)
    con.close()
    print(f"wrote {LEDGER}", flush=True)


if __name__ == "__main__":
    main()
