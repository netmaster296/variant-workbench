# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Chromosome, variant, and identifier rules for the GRCh38 warehouse.

Every later script imports this module. Coordinate policy lives here once.
"""

from __future__ import annotations

import re

_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
_ENST = re.compile(r"^(ENS[TG]\d+)(?:\.\d+)?$")
_RS = re.compile(r"^rs(\d+)$", re.IGNORECASE)
_DIGITS = re.compile(r"^(\d+)$")


def canon_chrom(raw) -> str | None:
    """Return chr1..chr22, chrX, chrY, chrM, or None."""
    if raw is None:
        return None
    text = str(raw).strip()
    if not text or text in {".", "-", "NA", "na", "None", "nan"}:
        return None
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.strip().upper()
    if text in {"M", "MT"}:
        return "chrM"
    if text in _CHROMS:
        return "chr" + text
    return None


def variant_key(chrom, pos, ref, alt) -> str | None:
    """1-based key chrom:pos:ref:alt with uppercase alleles."""
    name = canon_chrom(chrom)
    if name is None or pos is None or ref is None or alt is None:
        return None
    try:
        position = int(pos)
    except (TypeError, ValueError):
        return None
    if position < 1:
        return None
    ref_s = str(ref).strip().upper()
    alt_s = str(alt).strip().upper()
    if not ref_s or not alt_s or ref_s in {".", "-"} or alt_s in {".", "-"}:
        return None
    return f"{name}:{position}:{ref_s}:{alt_s}"


def bare_ensg(gene_id) -> str | None:
    return _bare_ensembl(gene_id, "ENSG")


def bare_enst(transcript_id) -> str | None:
    return _bare_ensembl(transcript_id, "ENST")


def _bare_ensembl(value, prefix: str) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {".", "-"}:
        return None
    match = _ENST.match(text)
    if not match:
        return None
    bare = match.group(1)
    if not bare.startswith(prefix):
        return None
    return bare


def norm_rsid(value) -> str | None:
    """Accept rs123, RS123, or a bare integer. Empty, '.', and -1 become None."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {".", "-", "NA", "na", "None", "nan"}:
        return None
    match = _RS.match(text)
    if match:
        number = match.group(1)
        if number.lstrip("0") == "" or int(number) < 1:
            return None
        return "rs" + str(int(number))
    match = _DIGITS.match(text)
    if match:
        number = int(match.group(1))
        if number < 1:
            return None
        return "rs" + str(number)
    return None


def _check() -> None:
    assert canon_chrom("1") == "chr1"
    assert canon_chrom("chr1") == "chr1"
    assert canon_chrom("CHR1") == "chr1"
    assert canon_chrom("X") == "chrX"
    assert canon_chrom("chrX") == "chrX"
    assert canon_chrom("MT") == "chrM"
    assert canon_chrom("M") == "chrM"
    assert canon_chrom("chrM") == "chrM"
    assert canon_chrom("Un") is None
    assert canon_chrom("") is None
    assert variant_key("1", 12345, "a", "g") == "chr1:12345:A:G"
    assert variant_key("chr1", "12345", "A", "G") == "chr1:12345:A:G"
    assert variant_key("1", 10, ".", "A") is None
    assert bare_ensg("ENSG00000123456.15") == "ENSG00000123456"
    assert bare_ensg("ENSG00000123456") == "ENSG00000123456"
    assert bare_enst("ENST00000263100.8") == "ENST00000263100"
    assert bare_ensg("ENST00000263100") is None
    assert norm_rsid("rs123") == "rs123"
    assert norm_rsid("RS123") == "rs123"
    assert norm_rsid("123") == "rs123"
    assert norm_rsid("-1") is None
    assert norm_rsid(".") is None
    assert norm_rsid("") is None
    assert norm_rsid("junk") is None
    print("coords.py checks passed")


if __name__ == "__main__":
    _check()
