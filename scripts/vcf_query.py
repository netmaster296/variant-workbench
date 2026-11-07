# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Region queries for bgzip+tabix VCFs.

conda-forge and bioconda do not publish bcftools for win-64, and pysam does
not ship a Windows wheel. This module reads the htslib TBI index and seeks
the BGZF file. It is the Wave 0 substitute for `bcftools view -H -r`.
"""

from __future__ import annotations

import gzip
import os
import struct
from pathlib import Path

from Bio import bgzf

WINDOW = 1 << 14


class TabixIndex:
    def __init__(self, path: Path):
        data = gzip.decompress(Path(path).read_bytes())
        if data[:4] != b"TBI\x01":
            raise ValueError(f"not a TBI index: {path}")
        self.n_ref = struct.unpack_from("<i", data, 4)[0]
        off = 8
        self.conf = struct.unpack_from("<6i", data, off)
        off += 24
        self.names: list[str] = []
        for _ in range(self.n_ref):
            length = struct.unpack_from("<i", data, off)[0]
            off += 4
            raw = data[off : off + length]
            off += length
            self.names.append(raw.split(b"\x00", 1)[0].decode("ascii"))
        self.linear: dict[str, tuple[int, ...]] = {}
        for name in self.names:
            n_bin = struct.unpack_from("<i", data, off)[0]
            off += 4
            for _ in range(n_bin):
                _bin, n_chunk = struct.unpack_from("<II", data, off)
                off += 8 + n_chunk * 16
            n_intv = struct.unpack_from("<i", data, off)[0]
            off += 4
            self.linear[name] = struct.unpack_from("<" + "Q" * n_intv, data, off)
            off += n_intv * 8
        self.consumed = off
        self.size = len(data)

    def virtual_offset(self, chrom: str, pos: int) -> int | None:
        offsets = self.linear.get(chrom)
        if not offsets:
            return None
        window = max(0, int(pos)) >> 14
        if window >= len(offsets):
            window = len(offsets) - 1
        return offsets[window]


def _info(text: str) -> dict[str, str]:
    out = {}
    for field in text.split(";"):
        if not field or field == ".":
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value
        else:
            out[field] = "1"
    return out


def _usable(line: str) -> bool:
    return line.startswith("#") or line.startswith("chr")


def query(vcf_path, chrom: str, pos: int, ref: str | None = None, alt: str | None = None, limit: int = 20000):
    """Return VCF records at this 1-based position. ref and alt are optional filters."""
    vcf_path = Path(vcf_path)
    index = TabixIndex(Path(str(vcf_path) + ".tbi"))
    virtual = index.virtual_offset(chrom, pos)
    if virtual is None:
        return []
    found = []
    with bgzf.BgzfReader(str(vcf_path)) as fh:
        fh.seek(virtual)
        line = fh.readline()
        if line and not _usable(line):
            line = fh.readline()
        seen = 0
        while line and seen < limit:
            seen += 1
            if line.startswith("#"):
                line = fh.readline()
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                line = fh.readline()
                continue
            row_pos = int(cols[1])
            if row_pos > pos:
                break
            if row_pos == pos and cols[0] == chrom:
                alts = cols[4].split(",")
                if ref is not None and cols[3] != ref:
                    line = fh.readline()
                    continue
                if alt is not None and alt not in alts:
                    line = fh.readline()
                    continue
                found.append(
                    {
                        "chrom": cols[0],
                        "pos": row_pos,
                        "id": cols[2],
                        "ref": cols[3],
                        "alts": alts,
                        "filter": cols[6],
                        "info": _info(cols[7]),
                        "source": str(vcf_path),
                    }
                )
            line = fh.readline()
    return found


def find_af_disagreement(vcf_path, chrom: str, samples: int = 24, per_window: int = 1500):
    """Find one site where AF_nfe and AF_afr are both present and differ."""
    vcf_path = Path(vcf_path)
    index = TabixIndex(Path(str(vcf_path) + ".tbi"))
    offsets = index.linear.get(chrom) or ()
    if not offsets:
        return None
    step = max(1, len(offsets) // samples)
    with bgzf.BgzfReader(str(vcf_path)) as fh:
        for window in range(0, len(offsets), step):
            fh.seek(offsets[window])
            line = fh.readline()
            if line and not _usable(line):
                line = fh.readline()
            for _ in range(per_window):
                if not line or line.startswith("#"):
                    line = fh.readline()
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8 or cols[0] != chrom:
                    line = fh.readline()
                    continue
                info = _info(cols[7])
                nfe = info.get("AF_nfe")
                afr = info.get("AF_afr")
                if nfe not in (None, ".") and afr not in (None, ".") and nfe != afr:
                    return {
                        "chrom": cols[0],
                        "pos": int(cols[1]),
                        "ref": cols[3],
                        "alt": cols[4].split(",")[0],
                        "AF_nfe": nfe,
                        "AF_afr": afr,
                        "source": str(vcf_path),
                    }
                line = fh.readline()
    return None


if __name__ == "__main__":
    root = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
    vcf = root / "raw" / "gnomad" / "genomes" / "gnomad.genomes.v4.1.sites.chrY.vcf.bgz"
    index = TabixIndex(Path(str(vcf) + ".tbi"))
    print("chrY contigs", index.names, "windows", len(index.linear["chrY"]))
    hit = find_af_disagreement(vcf, "chrY", samples=8, per_window=400)
    print("af disagreement", hit)
