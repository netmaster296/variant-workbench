# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load the small GRCh38 tables into one SQLite catalog.

Big VCFs and QTL archives stay on disk. This file records where they are
and loads the tables that can be joined: genes, ClinVar, GWAS hits,
gnomAD constraint.
"""

from __future__ import annotations

import csv
import gzip
import os
import sqlite3
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
DB = ROOT / "db" / "variants.sqlite"
RAW = ROOT / "raw"


def connect() -> sqlite3.Connection:
    DB.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS source_file (
            name TEXT PRIMARY KEY,
            path TEXT NOT NULL,
            bytes INTEGER,
            note TEXT
        );
        CREATE TABLE IF NOT EXISTS gene (
            gene_id TEXT PRIMARY KEY,
            symbol TEXT,
            chrom TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            biotype TEXT
        );
        CREATE TABLE IF NOT EXISTS clinvar (
            allele_id TEXT,
            variation_id TEXT,
            gene_symbol TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            rsid TEXT,
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            phenotype TEXT
        );
        CREATE TABLE IF NOT EXISTS gwas_hit (
            pubmed_id TEXT,
            trait TEXT,
            mapped_gene TEXT,
            chrom TEXT,
            pos INTEGER,
            snp TEXT,
            risk_allele TEXT,
            p_value TEXT,
            rsid TEXT
        );
        CREATE TABLE IF NOT EXISTS gene_constraint (
            gene TEXT,
            transcript TEXT,
            canonical TEXT,
            lof_oe TEXT,
            lof_pLI TEXT,
            lof_oe_ci_upper TEXT,
            mis_z TEXT
        );
        """
    )


def record(conn: sqlite3.Connection, name: str, path: Path, note: str) -> None:
    if not path.exists():
        return
    conn.execute(
        "INSERT OR REPLACE INTO source_file(name, path, bytes, note) VALUES (?,?,?,?)",
        (name, str(path), path.stat().st_size, note),
    )


def attr(field: str, key: str) -> str:
    token = key + ' "'
    start = field.find(token)
    if start < 0:
        return ""
    start += len(token)
    end = field.find('"', start)
    return field[start:end]


def load_genes(conn: sqlite3.Connection) -> int:
    path = RAW / "gencode" / "gencode.v49.annotation.gtf.gz"
    if not path.exists():
        return 0
    conn.execute("DELETE FROM gene")
    rows = []
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            rows.append(
                (
                    attr(parts[8], "gene_id"),
                    attr(parts[8], "gene_name"),
                    parts[0],
                    int(parts[3]),
                    int(parts[4]),
                    parts[6],
                    attr(parts[8], "gene_type"),
                )
            )
            if len(rows) >= 5000:
                conn.executemany("INSERT OR REPLACE INTO gene VALUES (?,?,?,?,?,?,?)", rows)
                n += len(rows)
                rows.clear()
    if rows:
        conn.executemany("INSERT OR REPLACE INTO gene VALUES (?,?,?,?,?,?,?)", rows)
        n += len(rows)
    conn.execute("CREATE INDEX IF NOT EXISTS gene_symbol ON gene(symbol)")
    conn.execute("CREATE INDEX IF NOT EXISTS gene_locus ON gene(chrom, start, end)")
    record(conn, "gencode_v49", path, "GRCh38 gene models, gene rows only")
    return n


def load_clinvar(conn: sqlite3.Connection) -> int:
    path = RAW / "clinvar" / "variant_summary.txt.gz"
    if not path.exists():
        return 0
    conn.execute("DELETE FROM clinvar")
    n = 0
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        reader = csv.reader(fh, delimiter="\t")
        for rec in reader:
            if len(rec) <= idx.get("Assembly", 0):
                continue
            if rec[idx["Assembly"]] != "GRCh38":
                continue
            pos_raw = rec[idx["PositionVCF"]] or rec[idx["Start"]]
            try:
                pos = int(pos_raw)
            except ValueError:
                continue
            ref = rec[idx["ReferenceAlleleVCF"]] or rec[idx["ReferenceAllele"]]
            alt = rec[idx["AlternateAlleleVCF"]] or rec[idx["AlternateAllele"]]
            rows.append(
                (
                    rec[idx["AlleleID"]],
                    rec[idx["VariationID"]],
                    rec[idx["GeneSymbol"]],
                    rec[idx["ClinicalSignificance"]],
                    rec[idx["ReviewStatus"]],
                    rec[idx["RS# (dbSNP)"]],
                    rec[idx["Chromosome"]],
                    pos,
                    ref,
                    alt,
                    rec[idx["PhenotypeList"]][:500],
                )
            )
            if len(rows) >= 5000:
                conn.executemany("INSERT INTO clinvar VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
                n += len(rows)
                rows.clear()
    if rows:
        conn.executemany("INSERT INTO clinvar VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
        n += len(rows)
    conn.execute("CREATE INDEX IF NOT EXISTS clinvar_locus ON clinvar(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS clinvar_gene ON clinvar(gene_symbol)")
    record(conn, "clinvar_variant_summary", path, "GRCh38 rows from the weekly summary")
    return n


def load_gwas(conn: sqlite3.Connection) -> int:
    zpath = RAW / "gwas" / "gwas-catalog-associations_ontology-annotated-full.zip"
    if not zpath.exists():
        return 0
    conn.execute("DELETE FROM gwas_hit")
    n = 0
    rows = []
    with zipfile.ZipFile(zpath) as zf:
        name = next(x for x in zf.namelist() if x.endswith(".tsv"))
        with zf.open(name) as raw:
            text = (line.decode("utf-8", "replace") for line in raw)
            reader = csv.DictReader(text, delimiter="\t")
            for rec in reader:
                rows.append(
                    (
                        rec.get("PUBMEDID", ""),
                        rec.get("DISEASE/TRAIT", ""),
                        rec.get("MAPPED_GENE", ""),
                        rec.get("CHR_ID", ""),
                        int(rec["CHR_POS"]) if (rec.get("CHR_POS") or "").isdigit() else None,
                        rec.get("SNPS", ""),
                        rec.get("STRONGEST SNP-RISK ALLELE", ""),
                        rec.get("P-VALUE", ""),
                        rec.get("SNP_ID_CURRENT", ""),
                    )
                )
                if len(rows) >= 5000:
                    conn.executemany("INSERT INTO gwas_hit VALUES (?,?,?,?,?,?,?,?,?)", rows)
                    n += len(rows)
                    rows.clear()
    if rows:
        conn.executemany("INSERT INTO gwas_hit VALUES (?,?,?,?,?,?,?,?,?)", rows)
        n += len(rows)
    conn.execute("CREATE INDEX IF NOT EXISTS gwas_locus ON gwas_hit(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS gwas_gene ON gwas_hit(mapped_gene)")
    record(conn, "gwas_associations", zpath, "GWAS Catalog associations, ontology-annotated, 2026-11-02 release")
    return n


def pick(header: list[str], *names: str) -> str | None:
    have = {h.lower(): h for h in header}
    for name in names:
        if name.lower() in have:
            return have[name.lower()]
    return None


def load_constraint(conn: sqlite3.Connection) -> int:
    path = RAW / "gnomad" / "constraint" / "gnomad.v4.1.constraint_metrics.tsv"
    if not path.exists():
        return 0
    conn.execute("DELETE FROM gene_constraint")
    n = 0
    rows = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        gene = pick(header, "gene", "gene_symbol")
        transcript = pick(header, "transcript", "transcript_id")
        canonical = pick(header, "canonical")
        oe = pick(header, "lof.oe", "lof_oe")
        pli = pick(header, "lof.pLI", "lof_pli")
        ci = pick(header, "lof.oe_ci.upper", "lof_oe_ci_upper")
        mis = pick(header, "mis.z_score", "mis_z_score")
        if not gene:
            record(conn, "gnomad_constraint", path, "downloaded; header did not match expected columns")
            return 0
        for rec in reader:
            rows.append(
                (
                    rec.get(gene, ""),
                    rec.get(transcript, "") if transcript else "",
                    rec.get(canonical, "") if canonical else "",
                    rec.get(oe, "") if oe else "",
                    rec.get(pli, "") if pli else "",
                    rec.get(ci, "") if ci else "",
                    rec.get(mis, "") if mis else "",
                )
            )
            if len(rows) >= 5000:
                conn.executemany("INSERT INTO gene_constraint VALUES (?,?,?,?,?,?,?)", rows)
                n += len(rows)
                rows.clear()
    if rows:
        conn.executemany("INSERT INTO gene_constraint VALUES (?,?,?,?,?,?,?)", rows)
        n += len(rows)
    conn.execute("CREATE INDEX IF NOT EXISTS gene_constraint_gene ON gene_constraint(gene)")
    record(conn, "gnomad_v4_1_constraint", path, "gene constraint, not variant sites")
    return n


def catalog_archives(conn: sqlite3.Connection) -> None:
    notes = {
        RAW / "gtex" / "GTEx_Analysis_v8_eQTL.tar": "significant cis-eQTL pairs, 49 tissues, kept as archive",
        RAW / "gtex" / "GTEx_Analysis_v8_sQTL.tar": "significant cis-sQTL pairs, kept as archive",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz": "gene TPM by sample",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz": "median TPM by tissue",
        RAW / "mavedb" / "mavedb-dump.2026062418131.tar.gz": "MaveDB CC0 bulk, measured variant effects",
    }
    for path, note in notes.items():
        record(conn, path.name, path, note)
    for kind in ("genomes", "exomes"):
        folder = RAW / "gnomad" / kind
        if not folder.exists():
            continue
        for path in sorted(folder.glob("*.vcf.bgz")):
            record(conn, path.name, path, f"gnomAD v4.1 {kind} sites VCF, GRCh38, not loaded into SQLite")


def main() -> None:
    conn = connect()
    schema(conn)
    counts = {
        "genes": load_genes(conn),
        "clinvar": load_clinvar(conn),
        "gwas": load_gwas(conn),
        "constraint": load_constraint(conn),
    }
    catalog_archives(conn)
    conn.commit()
    print("LOADED", counts)
    for table in ("gene", "clinvar", "gwas_hit", "gene_constraint", "source_file"):
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"{table}\t{n}")


if __name__ == "__main__":
    main()
