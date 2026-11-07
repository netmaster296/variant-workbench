# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Download the public human-variant sources that fit on this disk.

Assembly for every coordinate file is GRCh38. Huge VCFs stay as bgzip files.
They are not loaded row-by-row. The SQLite catalog is built by load_core.py.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
RAW = ROOT / "raw"
LOG = ROOT / "logs" / "fetch.log"
STATUS = ROOT / "logs" / "fetch.status"
RESERVE = 80 * 1024**3
UA = "variant-workbench/0.1 (research; local catalog)"

CORE = [
    (
        "clinvar_variant_summary",
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        RAW / "clinvar" / "variant_summary.txt.gz",
    ),
    (
        "gencode_v49_gtf",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
        RAW / "gencode" / "gencode.v49.annotation.gtf.gz",
    ),
    (
        "gwas_associations",
        "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
        RAW / "gwas" / "gwas-catalog-associations_ontology-annotated-full.zip",
    ),
    (
        "gnomad_v4_1_constraint",
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv",
        RAW / "gnomad" / "constraint" / "gnomad.v4.1.constraint_metrics.tsv",
    ),
    (
        "gtex_v8_eqtl",
        "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar",
        RAW / "gtex" / "GTEx_Analysis_v8_eQTL.tar",
    ),
    (
        "gtex_v8_sqtl",
        "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL.tar",
        RAW / "gtex" / "GTEx_Analysis_v8_sQTL.tar",
    ),
    (
        "gtex_v8_gene_tpm",
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
    ),
    (
        "gtex_v8_gene_median_tpm",
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        RAW / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    ),
    (
        "mavedb_bulk",
        "https://zenodo.org/api/records/20840937/files/mavedb-dump.2026062418131.tar.gz/content",
        RAW / "mavedb" / "mavedb-dump.2026062418131.tar.gz",
    ),
]

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def log(msg: str, milestone: bool = False) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    line = msg.rstrip()
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    if milestone:
        print(line, flush=True)
        STATUS.write_text(line + "\n", encoding="utf-8")


def remote_size(url: str) -> int | None:
    proc = subprocess.run(
        ["curl.exe", "-sI", "-L", "-A", UA, "--max-time", "60", url],
        capture_output=True,
        text=True,
        check=False,
    )
    size = None
    for line in proc.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            try:
                size = int(line.split(":", 1)[1].strip())
            except ValueError:
                pass
    return size


def free_bytes() -> int:
    return shutil.disk_usage(ROOT).free


def fetch(name: str, url: str, dest: Path) -> str:
    """Return ok, fail, or stop."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    size = remote_size(url)
    if dest.exists() and size and dest.stat().st_size == size:
        log(f"SKIP {name} already complete ({dest.stat().st_size} bytes)")
        return "ok"
    need = size or 0
    if free_bytes() < need + RESERVE:
        log(
            f"STOP {name} needs {need} bytes plus {RESERVE} reserve, free {free_bytes()}",
            milestone=True,
        )
        return "stop"
    log(f"GET {name} {size} -> {dest}")
    proc = subprocess.run(
        [
            "curl.exe",
            "-L",
            "--fail",
            "--retry",
            "8",
            "--retry-all-errors",
            "--retry-delay",
            "5",
            "-C",
            "-",
            "-A",
            UA,
            "-o",
            str(dest),
            url,
        ],
        check=False,
    )
    ok = proc.returncode == 0 and dest.exists() and (not size or dest.stat().st_size == size)
    if not ok:
        log(f"FILE_FAIL {name} exit={proc.returncode} path={dest}", milestone=True)
        return "fail"
    log(f"OK {name} {dest.stat().st_size}")
    return "ok"


def gnomad_urls(kind: str) -> list[tuple[str, str, Path]]:
    items = []
    folder = RAW / "gnomad" / kind
    for chrom in CHROMS:
        base = f"gnomad.{kind}.v4.1.sites.chr{chrom}.vcf.bgz"
        url = (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
            f"{kind}/{base}"
        )
        items.append((f"gnomad_{kind}_chr{chrom}", url, folder / base))
        items.append((f"gnomad_{kind}_chr{chrom}_tbi", url + ".tbi", folder / (base + ".tbi")))
    return items


def main() -> int:
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    failed = False
    if phase in {"core", "all"}:
        for item in CORE:
            if fetch(*item) == "fail":
                failed = True
        log("PHASE1_DONE", milestone=True)
    if phase in {"gnomad", "all"}:
        for kind in ("genomes", "exomes"):
            for item in gnomad_urls(kind):
                result = fetch(*item)
                if result == "stop":
                    failed = True
                    break
                if result == "fail":
                    failed = True
            else:
                continue
            break
    if failed:
        log("FAILED fetch incomplete", milestone=True)
        return 1
    log("DONE", milestone=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
