# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Freeze the on-disk raw tree and rebuild the SQLite source catalog.

Wave 0. Does not download, and does not read VCF rows.
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
RAW = ROOT / "raw"
DB = ROOT / "db" / "variants.sqlite"
LOG = ROOT / "logs" / "fetch.log"
STATUS = ROOT / "logs" / "fetch.status"
OUT_DIR = ROOT / "warehouse" / "manifest"
# Label for the disk lines in disk.md: the drive letter on Windows ("D:").
DRIVE = ROOT.drive or "volume"
FILES = OUT_DIR / "files.tsv"
DISK = OUT_DIR / "disk.md"

GENOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]


def parse_log(path: Path) -> dict[str, dict]:
    """Map a fetch name to its destination, reported size, and OK/SKIP state."""
    found: dict[str, dict] = {}
    if not path.exists():
        return found
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("GET "):
            # GET name size -> dest
            try:
                head, dest = line.split(" -> ", 1)
            except ValueError:
                continue
            parts = head.split()
            if len(parts) < 3:
                continue
            name, size = parts[1], parts[2]
            row = found.setdefault(name, {})
            row["dest"] = dest.strip()
            try:
                row["reported"] = int(size)
            except ValueError:
                row["reported"] = None
        elif line.startswith("OK "):
            parts = line.split()
            if len(parts) >= 3:
                row = found.setdefault(parts[1], {})
                try:
                    row["ok_bytes"] = int(parts[2])
                except ValueError:
                    row["ok_bytes"] = None
                row["state"] = "OK"
        elif line.startswith("SKIP "):
            parts = line.split()
            if len(parts) >= 2:
                row = found.setdefault(parts[1], {})
                row["state"] = "SKIP"
                for token in parts:
                    if token.endswith("bytes)"):
                        token = token[:-1]
                    if token.isdigit():
                        row["ok_bytes"] = int(token)
    return found


def expected_notes() -> dict[str, str]:
    notes = {
        "clinvar": "ClinVar weekly summary, GRCh37 and GRCh38 mixed",
        "gencode": "GENCODE v49 gene models",
        "gwas": "GWAS Catalog ontology-annotated associations",
        "constraint": "gnomAD v4.1 gene and transcript constraint",
        "gene_median_tpm": "GTEx v8 median TPM by tissue",
        "gene_tpm": "GTEx v8 sample-level TPM, deferred",
        "eQTL": "GTEx v8 single-tissue cis-eQTL tar",
        "sQTL": "GTEx v8 single-tissue cis-sQTL tar",
        "mavedb": "MaveDB CC0 bulk dump",
    }
    return notes


def note_for(path: Path) -> str:
    text = str(path).replace("/", "\\").lower()
    if "clinvar" in text:
        return "ClinVar weekly summary, GRCh37 and GRCh38 mixed"
    if "gencode" in text:
        return "GENCODE v49 annotation"
    if "gwas" in text:
        return "GWAS Catalog associations"
    if "constraint" in text:
        return "gnomAD v4.1 constraint metrics"
    if "gene_median_tpm" in text:
        return "GTEx v8 median TPM"
    if "gene_tpm" in text:
        return "GTEx v8 sample-level TPM, deferred"
    if "eqtl" in text:
        return "GTEx v8 eQTL tar"
    if "sqtl" in text:
        return "GTEx v8 sQTL tar"
    if "mavedb" in text:
        return "MaveDB bulk dump"
    if "gnomad\\exomes" in text or "gnomad/exomes" in text:
        kind = "index" if path.name.endswith(".tbi") else "sites VCF"
        return f"gnomAD v4.1 exome {kind}, not loaded as rows"
    if "gnomad\\genomes" in text or "gnomad/genomes" in text:
        kind = "index" if path.name.endswith(".tbi") else "sites VCF"
        return f"gnomAD v4.1 genome {kind}, not loaded as rows"
    return "raw file"


def chromosome_checklist(label: str, folder: Path) -> list[str]:
    lines = []
    missing = []
    for chrom in GENOMES:
        vcf = folder / f"gnomad.{label}.v4.1.sites.chr{chrom}.vcf.bgz"
        tbi = Path(str(vcf) + ".tbi")
        if vcf.exists() and tbi.exists() and vcf.stat().st_size > 0 and tbi.stat().st_size > 0:
            lines.append(f"- {label} chr{chrom}: OK {vcf.stat().st_size} bytes, tbi {tbi.stat().st_size} bytes")
        else:
            missing.append(chrom)
            lines.append(f"- {label} chr{chrom}: MISSING vcf={vcf.exists()} tbi={tbi.exists()}")
    lines.append(f"- {label} missing: {', '.join(missing) if missing else 'none'}")
    return lines


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_rows = parse_log(LOG)
    by_dest = {}
    for name, row in log_rows.items():
        dest = row.get("dest")
        if dest:
            by_dest[str(Path(dest)).lower()] = (name, row)

    files = sorted(p for p in RAW.rglob("*") if p.is_file())
    header = [
        "path",
        "bytes",
        "fetch_name",
        "fetch_state",
        "fetch_bytes",
        "size_matches_log",
        "note",
    ]
    matched = 0
    with FILES.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(header) + "\n")
        for path in files:
            key = str(path).lower()
            name, row = by_dest.get(key, ("", {}))
            size = path.stat().st_size
            state = row.get("state", "")
            ok_bytes = row.get("ok_bytes")
            size_ok = bool(state in {"OK", "SKIP"} and ok_bytes is not None and ok_bytes == size)
            if size_ok:
                matched += 1
            fh.write(
                "\t".join(
                    [
                        str(path),
                        str(size),
                        name,
                        state,
                        "" if ok_bytes is None else str(ok_bytes),
                        "yes" if size_ok else "no",
                        note_for(path),
                    ]
                )
                + "\n"
            )

    import shutil

    usage = shutil.disk_usage(ROOT)
    status = STATUS.read_text(encoding="utf-8").strip() if STATUS.exists() else "MISSING"
    disk_lines = [
        "# Disk",
        "",
        f"- fetch.status: `{status}`",
        f"- {DRIVE} total bytes: {usage.total}",
        f"- {DRIVE} used bytes: {usage.used}",
        f"- {DRIVE} free bytes: {usage.free}",
        f"- Raw files: {len(files)}",
        f"- Raw files whose size matches an OK or SKIP log line: {matched}",
        "",
        "## Chromosomes",
        "",
    ]
    disk_lines.extend(chromosome_checklist("genomes", RAW / "gnomad" / "genomes"))
    disk_lines.extend(chromosome_checklist("exomes", RAW / "gnomad" / "exomes"))
    core = [
        RAW / "clinvar" / "variant_summary.txt.gz",
        RAW / "gencode" / "gencode.v49.annotation.gtf.gz",
        RAW / "gwas" / "gwas-catalog-associations_ontology-annotated-full.zip",
        RAW / "gnomad" / "constraint" / "gnomad.v4.1.constraint_metrics.tsv",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
        RAW / "gtex" / "GTEx_Analysis_v8_eQTL.tar",
        RAW / "gtex" / "GTEx_Analysis_v8_sQTL.tar",
        RAW / "mavedb" / "mavedb-dump.2026062418131.tar.gz",
    ]
    disk_lines.extend(["", "## Core files", ""])
    for path in core:
        disk_lines.append(f"- {'OK' if path.exists() and path.stat().st_size > 0 else 'MISSING'} `{path}`")
    DISK.write_text("\n".join(disk_lines) + "\n", encoding="utf-8")

    conn = sqlite3.connect(DB)
    conn.execute("DELETE FROM source_file")
    conn.executemany(
        "INSERT INTO source_file(name, path, bytes, note) VALUES (?,?,?,?)",
        [
            (
                path.name,
                str(path),
                path.stat().st_size,
                note_for(path),
            )
            for path in files
        ],
    )
    conn.commit()
    count = conn.execute("SELECT COUNT(*) FROM source_file").fetchone()[0]
    conn.close()
    print(f"files {len(files)} log_matches {matched} source_file {count} status {status} free {usage.free}")
    print(f"wrote {FILES}")
    print(f"wrote {DISK}")
    return 0 if status == "DONE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
