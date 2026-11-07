# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Chromosome map of lead SNPs for four cognitive and educational traits.

Colors mark 2 Mb bins that contain a clumped GWAS lead. The clump span is
measured and reported. It is not painted, because those spans are the
catalog's neighborhood around a lead, not a segment that causes the trait.
"""

from __future__ import annotations

import os
from pathlib import Path

import duckdb
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
LOCI = str(ROOT / "warehouse" / "gwas" / "loci.parquet")
OUT = str(ROOT / "warehouse" / "reports" / "figures" / "cognitive_loci_map.png")

LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
TRAITS = [
    ("Educational attainment", "http://www.ebi.ac.uk/efo/EFO_0011015", (0, 114, 178)),
    ("Self-reported educational attainment", "http://www.ebi.ac.uk/efo/EFO_0004784", (0, 158, 115)),
    ("Intelligence", "http://www.ebi.ac.uk/efo/EFO_0004337", (168, 78, 138)),
    ("Cognitive function measurement", "http://www.ebi.ac.uk/efo/EFO_0008354", (214, 137, 16)),
]
BIN = 2_000_000


def font(size: int, bold: bool = False):
    path = r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf"
    try:
        return ImageFont.truetype(path, size)
    except OSError:
        pass
    # Arial is missing (not Windows, or no Fonts folder): use PIL's default font.
    # Pillow 10.1+ can scale it; older Pillow only has the fixed bitmap font.
    try:
        return ImageFont.load_default(size=size)
    except (TypeError, OSError, ImportError):
        return ImageFont.load_default()


def load_rows():
    uris = [uri for _, uri, _ in TRAITS]
    con = duckdb.connect()
    rows = con.execute(
        """
        SELECT trait_uri, chrom, lead_pos, "start", "end"
        FROM read_parquet(?)
        WHERE trait_uri IN (SELECT * FROM unnest(?))
          AND chrom IS NOT NULL
          AND lead_pos IS NOT NULL
          AND lead_pos > 0
        """,
        [LOCI, uris],
    ).fetchall()
    return rows


def shade(rgb, count: int):
    strength = min(1.0, 0.35 + 0.2 * count)
    return tuple(int(255 - (255 - channel) * strength) for channel in rgb)


def main() -> None:
    rows = load_rows()
    by_trait: dict[str, list[tuple[str, int, int, int]]] = {uri: [] for _, uri, _ in TRAITS}
    for uri, chrom, lead, start, end in rows:
        if chrom not in LENGTHS:
            continue
        by_trait[uri].append((chrom, int(lead), int(start or lead), int(end or lead)))

    max_len = LENGTHS["chr1"]
    left, right = 118, 36
    plot_w = 1680
    scale = plot_w / max_len
    track_h, track_gap, chrom_gap = 11, 2, 16
    block_h = 4 * track_h + 3 * track_gap
    top = 168
    height = top + len(ORDER) * (block_h + chrom_gap) + 78
    width = left + plot_w + right
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_f = font(28, True)
    sub_f = font(16)
    label_f = font(15, True)
    small_f = font(14)
    tiny_f = font(13)

    draw.text((left, 22), "Where the catalog places lead variants for cognition and schooling", fill=(20, 20, 20), font=title_f)
    draw.text(
        (left, 62),
        "Each colored cell is a 2 Mb stretch that contains one or more lead SNPs. Darker means more leads in that stretch.",
        fill=(50, 50, 50),
        font=sub_f,
    )
    draw.text(
        (left, 86),
        "A lead is the variant a paper highlighted. The cell is not a gene, and it is not a segment that causes intelligence.",
        fill=(50, 50, 50),
        font=sub_f,
    )

    legend_x = left
    for name, _, rgb in TRAITS:
        draw.rectangle((legend_x, 122, legend_x + 22, 138), fill=shade(rgb, 3))
        draw.text((legend_x + 28, 120), name, fill=(30, 30, 30), font=small_f)
        legend_x += 360 if "Self" not in name else 430

    covered = {uri: 0 for _, uri, _ in TRAITS}
    any_covered = 0
    leads_per = {uri: 0 for _, uri, _ in TRAITS}

    for index, chrom in enumerate(ORDER):
        y0 = top + index * (block_h + chrom_gap)
        length = LENGTHS[chrom]
        bar_w = max(1, int(length * scale))
        draw.text((12, y0 + block_h // 2 - 8), chrom, fill=(20, 20, 20), font=label_f)
        for track, (_, uri, rgb) in enumerate(TRAITS):
            y = y0 + track * (track_h + track_gap)
            draw.rectangle((left, y, left + bar_w, y + track_h), fill=(236, 236, 236))
            bins: dict[int, int] = {}
            intervals = []
            for c, lead, start, end in by_trait[uri]:
                if c != chrom:
                    continue
                leads_per[uri] += 1
                bins[lead // BIN] = bins.get(lead // BIN, 0) + 1
                a, b = sorted((max(1, start), min(length, end)))
                intervals.append((a, b))
            intervals.sort()
            merged = []
            for a, b in intervals:
                if not merged or a > merged[-1][1] + 1:
                    merged.append([a, b])
                else:
                    merged[-1][1] = max(merged[-1][1], b)
            covered[uri] += sum(b - a + 1 for a, b in merged)
            for bin_index, count in bins.items():
                x1 = left + int(bin_index * BIN * scale)
                x2 = left + int(min(length, (bin_index + 1) * BIN) * scale)
                x2 = max(x2, x1 + 1)
                draw.rectangle((x1, y, min(left + bar_w, x2), y + track_h), fill=shade(rgb, count))
        # Union of the four traits' clump spans, for the printed summary only.
        all_intervals = []
        for _, uri, _ in TRAITS:
            for c, _, start, end in by_trait[uri]:
                if c != chrom:
                    continue
                a, b = sorted((max(1, start), min(length, end)))
                all_intervals.append((a, b))
        all_intervals.sort()
        union = []
        for a, b in all_intervals:
            if not union or a > union[-1][1] + 1:
                union.append([a, b])
            else:
                union[-1][1] = max(union[-1][1], b)
        any_covered += sum(b - a + 1 for a, b in union)

    genome = sum(LENGTHS.values())
    footer_y = height - 62
    draw.text(
        (left, footer_y),
        "Gray means no lead SNP from that trait in the 2 Mb cell. Chromosome lengths are the GRCh38 lengths used in the dark-genome note.",
        fill=(60, 60, 60),
        font=tiny_f,
    )
    parts = []
    for name, uri, _ in TRAITS:
        parts.append(f"{name}: {leads_per[uri]} leads, clump neighborhoods cover {covered[uri] / genome:.1%} of the assembly")
    draw.text((left, footer_y + 22), "  |  ".join(parts[:2]), fill=(60, 60, 60), font=tiny_f)
    draw.text((left, footer_y + 40), "  |  ".join(parts[2:]), fill=(60, 60, 60), font=tiny_f)

    image.save(OUT, "PNG")
    print("wrote", OUT)
    print("union_fraction", round(any_covered / genome, 4))
    for name, uri, _ in TRAITS:
        print(name, leads_per[uri], round(covered[uri] / genome, 4))


if __name__ == "__main__":
    main()
