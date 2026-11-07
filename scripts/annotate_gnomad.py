# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Annotate Tier A keys from the gnomAD sites VCFs.

Each chromosome is streamed once per dataset. Genome and exome stay in
separate files and separate columns. A key the VCF does not contain is kept
with in_genomes or in_exomes false. Allele counts are not added together.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

import vcf_stream  # noqa: E402
from coords import variant_key  # noqa: E402

KEYS = ROOT / "warehouse" / "variants" / "tier_a_keys.tsv"
OUT_DIR = ROOT / "warehouse" / "variants"
PARTS = OUT_DIR / "_parts"
GENOME_OUT = OUT_DIR / "tier_a_genome.parquet"
EXOME_OUT = OUT_DIR / "tier_a_exome.parquet"
JOIN_OUT = OUT_DIR / "tier_a.parquet"
FIELD_MAP = OUT_DIR / "field_map.tsv"
TMP = ROOT / "warehouse" / "tmp" / "annotate"

ORDER = (
    ["chrY", "chr22"]
    + [f"chr{n}" for n in range(21, 10, -1)]
    + [f"chr{n}" for n in range(10, 2, -1)]
    + ["chrX", "chr2", "chr1"]
)

INT_FIELDS = ("AC", "AN", "nhomalt")
FLOAT_FIELDS = (
    "AF",
    "AF_afr",
    "AF_amr",
    "AF_asj",
    "AF_eas",
    "AF_fin",
    "AF_mid",
    "AF_nfe",
    "AF_sas",
    "AF_ami",
    "AF_remaining",
    "AF_grpmax",
    "cadd_phred",
    "revel_max",
    "spliceai_ds_max",
    "pangolin_largest_ds",
    "phylop",
    "sift_max",
    "polyphen_max",
)
STR_FIELDS = ("grpmax", "allele_type")
VEP_COLUMNS = (
    ("Consequence", "vep_consequence"),
    ("IMPACT", "vep_impact"),
    ("SYMBOL", "vep_symbol"),
    ("Feature", "vep_enst"),
    ("HGVSc", "vep_hgvsc"),
    ("HGVSp", "vep_hgvsp"),
)
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}

# Acceptance sites. Values were read from these VCFs before the stream.
CHECKS = {
    "chr22:10510183:A:G": {
        "genomes": {"in_genomes": True, "ref": "A", "alt": "G", "AC": 2, "AF_nfe": 0.0, "AF_afr": 0.00229358},
        "exomes": {"in_exomes": False},
    },
    "chr11:126275389:C:T": {
        "genomes": {"in_genomes": True, "ref": "C", "alt": "T", "AC": 5},
        "exomes": {"in_exomes": True, "AC_exome": 53},
    },
    "chr7:117559590:ATCT:A": {
        "genomes": {"in_genomes": True, "ref": "ATCT", "alt": "A", "AC": 1200},
        "exomes": {"in_exomes": True, "AC_exome": 18037},
    },
}


def vcf_path(dataset: str, chrom: str) -> Path:
    if dataset == "genomes":
        name = f"gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz"
    elif dataset == "exomes":
        name = f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
    else:
        raise ValueError(dataset)
    return ROOT / "raw" / "gnomad" / dataset / name


def sql_path(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def schema_for(dataset: str) -> pa.Schema:
    suffix = "" if dataset == "genomes" else "_exome"
    flag = "in_genomes" if dataset == "genomes" else "in_exomes"
    fields = [
        pa.field("variant_key", pa.string()),
        pa.field("chrom", pa.string()),
        pa.field("pos", pa.int64()),
        pa.field("ref", pa.string()),
        pa.field("alt", pa.string()),
        pa.field(flag, pa.bool_()),
        pa.field("filter" + suffix, pa.string()),
    ]
    for name in INT_FIELDS:
        fields.append(pa.field(name + suffix, pa.int64()))
    for name in FLOAT_FIELDS:
        fields.append(pa.field(name + suffix, pa.float64()))
    for name in STR_FIELDS:
        fields.append(pa.field(name + suffix, pa.string()))
    for _src, column in VEP_COLUMNS:
        fields.append(pa.field(column + suffix, pa.string()))
    fields.append(pa.field("vep_on_mane" + suffix, pa.bool_()))
    return pa.schema(fields)


def load_keys() -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    with KEYS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["pos"] = int(row["pos"])
            row["ref_b"] = row["ref"].encode("ascii")
            row["alt_b"] = row["alt"].encode("ascii")
            grouped.setdefault(row["chrom"], []).append(row)
    return grouped


def as_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        try:
            number = float(value)
        except ValueError:
            return None
    if not number.is_integer():
        return None
    return int(number)


def as_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def choose_vep(raw: str | None, names: list[str]) -> dict:
    empty = {column: None for _src, column in VEP_COLUMNS}
    empty["vep_on_mane"] = None
    if not raw or not names:
        return empty
    index = {name: i for i, name in enumerate(names)}
    blocks = []
    for block in raw.split(","):
        parts = block.split("|")
        if len(parts) < len(names):
            parts.extend([""] * (len(names) - len(parts)))
        blocks.append(parts)
    if not blocks:
        return empty

    def rank(parts: list[str]) -> int:
        if "IMPACT" not in index:
            return 0
        return IMPACT_RANK.get(parts[index["IMPACT"]], 0)

    def pool(predicate) -> list[list[str]]:
        return [parts for parts in blocks if predicate(parts)]

    mane_idx = index.get("MANE_SELECT")
    canon_idx = index.get("CANONICAL")
    chosen_pool = blocks
    if mane_idx is not None:
        mane = pool(lambda parts: bool(parts[mane_idx].strip()))
        if mane:
            chosen_pool = mane
    if chosen_pool is blocks and canon_idx is not None:
        canon = pool(lambda parts: parts[canon_idx].strip().upper() == "YES")
        if canon:
            chosen_pool = canon
    chosen = max(enumerate(chosen_pool), key=lambda item: (rank(item[1]), -item[0]))[1]

    def field(name: str) -> str | None:
        slot = index.get(name)
        if slot is None:
            return None
        text = chosen[slot].strip()
        return text or None

    out = {column: field(src) for src, column in VEP_COLUMNS}
    if mane_idx is None:
        out["vep_on_mane"] = None
    else:
        out["vep_on_mane"] = bool(chosen[mane_idx].strip())
    return out


def extract(info: dict[str, str], index: int, numbers: dict[str, str], vep_names: list[str]) -> dict:
    out = {}
    for name in INT_FIELDS:
        raw = vcf_stream.allele_value(info.get(name), index, numbers.get(name, "A"))
        out[name] = as_int(raw)
    for name in FLOAT_FIELDS:
        raw = vcf_stream.allele_value(info.get(name), index, numbers.get(name, "A"))
        out[name] = as_float(raw)
    for name in STR_FIELDS:
        out[name] = vcf_stream.allele_value(info.get(name), index, numbers.get(name, "A"))
    out.update(choose_vep(info.get("vep"), vep_names))
    return out


def rename(record: dict, dataset: str) -> dict:
    if dataset == "genomes":
        return record
    kept = {"variant_key", "chrom", "pos", "ref", "alt", "in_genomes", "in_exomes"}
    out = {}
    for name, value in record.items():
        if name in kept:
            out[name] = value
        else:
            out[name + "_exome"] = value
    return out


def table_from_rows(rows: list[dict], dataset: str) -> pa.Table:
    layout = schema_for(dataset)
    columns = {field.name: [] for field in layout}
    for row in rows:
        for field in layout:
            columns[field.name].append(row.get(field.name))
    arrays = [pa.array(columns[field.name], type=field.type) for field in layout]
    return pa.table(arrays, schema=layout)


def gate(path: Path, chrom: str) -> None:
    header = vcf_stream.info_header(path)
    if "AF_nfe" not in header or "AF_afr" not in header:
        raise SystemExit(f"{path} header does not list AF_nfe and AF_afr as distinct INFO fields")
    prefix = (chrom + "\t").encode("ascii")
    for line in vcf_stream.iter_lines(path):
        if not line.startswith(prefix):
            continue
        info = vcf_stream.parse_info(line.split(b"\t", 7)[7].decode("ascii", errors="replace"))
        if "AF_nfe" in info and "AF_afr" in info:
            print(
                f"QC gate {path.name} AF_nfe={info['AF_nfe']} AF_afr={info['AF_afr']}",
                flush=True,
            )
            return
    raise SystemExit(f"{path.name}: streamer could not read AF_nfe and AF_afr as distinct fields")


def write_field_map(genome_header: dict, exome_header: dict, vep_names: list[str]) -> None:
    kept = ("AC", "AN", "AF", *FLOAT_FIELDS[1:], *STR_FIELDS, "vep")
    # FLOAT_FIELDS already includes AF. Don't list it twice in the file.
    wanted = []
    seen = set()
    for name in ("AC", "AN", "AF", "nhomalt", *FLOAT_FIELDS, *STR_FIELDS, "vep"):
        if name not in seen:
            wanted.append(name)
            seen.add(name)
    with FIELD_MAP.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["dataset", "info_id", "number", "type", "column", "present"])
        for dataset, header, suffix in (
            ("genomes", genome_header, ""),
            ("exomes", exome_header, "_exome"),
        ):
            for name in wanted:
                meta = header.get(name)
                column = "vep_*" if name == "vep" else name + suffix
                if meta is None:
                    writer.writerow([dataset, name, "", "", column, "false"])
                else:
                    writer.writerow(
                        [dataset, name, meta.get("Number", ""), meta.get("Type", ""), column, "true"]
                    )
        for src, column in VEP_COLUMNS:
            writer.writerow(["vep", src, "1", "String", column, "true" if src in vep_names else "false"])
        writer.writerow(["vep", "MANE_SELECT", "1", "String", "vep_on_mane", "true" if "MANE_SELECT" in vep_names else "false"])
    print(f"wrote {FIELD_MAP}", flush=True)
    print("QC vep_format=" + "|".join(vep_names), flush=True)


def annotate_chrom(dataset: str, chrom: str, keys: list[dict]) -> Path:
    dest = PARTS / dataset / f"{chrom}.parquet"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        rows = pq.read_metadata(dest).num_rows
        if rows == len(keys):
            print(f"reuse {dest.name} {dataset} {chrom} rows {rows}", flush=True)
            return dest
        dest.unlink()
    path = vcf_path(dataset, chrom)
    index = Path(str(path) + ".tbi")
    flag = "in_genomes" if dataset == "genomes" else "in_exomes"
    matches: dict[str, dict] = {}
    if not path.exists() or not index.exists() or path.stat().st_size <= 0:
        print(f"missing {path}", flush=True)
    else:
        header = vcf_stream.info_header(path)
        numbers = {name: meta.get("Number", "") for name, meta in header.items()}
        vep_names = vcf_stream.vep_format(header.get("vep", {}).get("Description", ""))
        by_pos: dict[int, list[dict]] = {}
        for key in keys:
            by_pos.setdefault(key["pos"], []).append(key)
        positions = sorted(by_pos)
        print(
            f"stream {dataset} {chrom} keys {len(keys)} positions {len(positions)} bytes {path.stat().st_size}",
            flush=True,
        )
        for pos, line in vcf_stream.stream_positions(path, chrom, positions):
            parts = line.split(b"\t")
            if len(parts) < 8:
                continue
            ref = parts[3]
            alts = parts[4].split(b",")
            hits = []
            for cand in by_pos[pos]:
                if cand["ref_b"] == ref and cand["alt_b"] in alts:
                    hits.append((cand, alts.index(cand["alt_b"])))
            if not hits:
                continue
            info = vcf_stream.parse_info(parts[7].decode("ascii", errors="replace"))
            filt = parts[6].decode("ascii", errors="replace")
            for cand, index in hits:
                record = extract(info, index, numbers, vep_names)
                record["filter"] = filt
                record[flag] = True
                previous = matches.get(cand["variant_key"])
                if previous is None or (previous.get("filter") != "PASS" and filt == "PASS"):
                    matches[cand["variant_key"]] = record
    rows_out = []
    for key in keys:
        hit = matches.get(key["variant_key"])
        base = {
            "variant_key": key["variant_key"],
            "chrom": key["chrom"],
            "pos": key["pos"],
            "ref": key["ref"],
            "alt": key["alt"],
            flag: hit is not None,
        }
        if hit is not None:
            base.update(hit)
            base[flag] = True
        renamed = rename(base, dataset)
        rows_out.append(renamed)
    table = table_from_rows(rows_out, dataset)
    temporary = dest.with_suffix(".parquet.tmp")
    pq.write_table(table, temporary, compression="zstd")
    temporary.replace(dest)
    print(
        f"wrote {dest} rows {table.num_rows} matched {len(matches)}",
        flush=True,
    )
    return dest


def check_part(dataset: str, key: str, expect: dict) -> None:
    chrom = key.split(":", 1)[0]
    path = PARTS / dataset / f"{chrom}.parquet"
    if not path.exists():
        raise SystemExit(f"missing part {path}")
    columns = list(expect)
    table = pq.read_table(path, columns=["variant_key", *columns], filters=[("variant_key", "=", key)])
    if table.num_rows != 1:
        raise SystemExit(f"{path.name} rows for {key}: {table.num_rows}")
    row = {name: table.column(name)[0].as_py() for name in columns}
    for name, want in expect.items():
        got = row[name]
        if isinstance(want, float):
            if got is None or abs(float(got) - want) > 1e-8:
                raise SystemExit(f"{key} {dataset} {name} {got} != {want}")
        elif got != want:
            raise SystemExit(f"{key} {dataset} {name} {got} != {want}")
    if dataset == "genomes" and "AF_nfe" in expect and "AF_afr" in expect:
        if float(row["AF_nfe"]) == float(row["AF_afr"]):
            raise SystemExit(f"{key} AF_nfe equals AF_afr")
    print(f"QC check {dataset} {key} ok", flush=True)


def combine(dataset: str, grouped: dict[str, list[dict]], dest: Path) -> None:
    writer = None
    total = 0
    try:
        for chrom in ORDER:
            part = PARTS / dataset / f"{chrom}.parquet"
            if chrom not in grouped:
                if part.exists():
                    raise SystemExit(f"unexpected part {part}")
                continue
            if not part.exists():
                raise SystemExit(f"missing {part}")
            table = pq.read_table(part)
            if table.num_rows != len(grouped[chrom]):
                raise SystemExit(f"{part} rows {table.num_rows} != keys {len(grouped[chrom])}")
            if writer is None:
                writer = pq.ParquetWriter(dest, table.schema, compression="zstd")
            writer.write_table(table)
            total += table.num_rows
        extra = [chrom for chrom in grouped if chrom not in ORDER]
        for chrom in extra:
            part = PARTS / dataset / f"{chrom}.parquet"
            table = pq.read_table(part)
            writer.write_table(table)
            total += table.num_rows
    finally:
        if writer is not None:
            writer.close()
    print(f"wrote {dest} rows {total}", flush=True)


def join(grouped: dict[str, list[dict]]) -> None:
    n_keys = sum(len(rows) for rows in grouped.values())
    con = duckdb.connect()
    con.execute(f"SET temp_directory = '{sql_path(TMP)}'")
    con.execute("SET memory_limit = '64GB'")
    con.execute("SET threads = 8")
    tab = "\t"
    con.execute(
        f"""
        COPY (
          SELECT k.source, g.*, e.* EXCLUDE (variant_key, chrom, pos, ref, alt)
          FROM read_parquet('{sql_path(GENOME_OUT)}') g
          JOIN read_parquet('{sql_path(EXOME_OUT)}') e USING (variant_key)
          JOIN read_csv(
            '{sql_path(KEYS)}',
            header = true,
            delim = '{tab}',
            columns = {{
              'chrom': 'VARCHAR',
              'pos': 'BIGINT',
              'ref': 'VARCHAR',
              'alt': 'VARCHAR',
              'variant_key': 'VARCHAR',
              'source': 'VARCHAR'
            }}
          ) k USING (variant_key)
        ) TO '{sql_path(JOIN_OUT)}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )
    rows = con.execute(f"SELECT count(*), count(DISTINCT variant_key) FROM read_parquet('{sql_path(JOIN_OUT)}')").fetchone()
    print(f"QC tier_a_rows={rows[0]} distinct={rows[1]}", flush=True)
    if rows[0] != n_keys or rows[1] != n_keys:
        raise SystemExit(f"joined rows {rows} != keys {n_keys}")
    both = con.execute(
        f"""
        SELECT
          count(*) FILTER (WHERE in_genomes),
          count(*) FILTER (WHERE in_exomes),
          count(*) FILTER (WHERE in_genomes AND in_exomes)
        FROM read_parquet('{sql_path(JOIN_OUT)}')
        """
    ).fetchone()
    print(f"QC in_genomes={both[0]} in_exomes={both[1]} in_both={both[2]}", flush=True)
    # The two count columns stay independent. This sum is only a sanity check
    # that the join did not replace one with the other.
    swapped = con.execute(
        f"""
        SELECT count(*)
        FROM read_parquet('{sql_path(JOIN_OUT)}')
        WHERE in_genomes AND in_exomes AND AC IS NOT NULL AND AC_exome IS NOT NULL AND AC = AC_exome
          AND variant_key = 'chr11:126275389:C:T'
        """
    ).fetchone()[0]
    if swapped:
        raise SystemExit("coding SNP genome AC equals exome AC; the columns may have been copied")
    con.close()


def main() -> None:
    if not KEYS.exists():
        raise SystemExit(f"missing {KEYS}")
    TMP.mkdir(parents=True, exist_ok=True)
    grouped = load_keys()
    print(f"QC loaded_keys={sum(len(v) for v in grouped.values())} chroms={len(grouped)}", flush=True)
    gate(vcf_path("genomes", "chrY"), "chrY")
    gate(vcf_path("exomes", "chrY"), "chrY")
    genome_header = vcf_stream.info_header(vcf_path("genomes", "chrY"))
    exome_header = vcf_stream.info_header(vcf_path("exomes", "chrY"))
    vep_names = vcf_stream.vep_format(genome_header.get("vep", {}).get("Description", ""))
    write_field_map(genome_header, exome_header, vep_names)
    missing = [name for name in ("AF_afr", "AF_nfe", "AF_eas", "AF_sas", "AF_amr") if name not in genome_header]
    if missing:
        raise SystemExit(f"genome header missing {missing}")
    for chrom in ORDER:
        keys = grouped.get(chrom, [])
        if not keys:
            print(f"no keys {chrom}", flush=True)
            continue
        for dataset in ("genomes", "exomes"):
            annotate_chrom(dataset, chrom, keys)
        for key, expect in CHECKS.items():
            if key.split(":", 1)[0] != chrom:
                continue
            rebuilt = variant_key(*key.split(":"))
            if rebuilt != key:
                raise SystemExit(f"truth key {key} failed coords.variant_key")
            for dataset, fields in expect.items():
                check_part(dataset, key, fields)
    for chrom, keys in grouped.items():
        if chrom in ORDER:
            continue
        for dataset in ("genomes", "exomes"):
            annotate_chrom(dataset, chrom, keys)
    combine("genomes", grouped, GENOME_OUT)
    combine("exomes", grouped, EXOME_OUT)
    join(grouped)
    for key, expect in CHECKS.items():
        check_joined(key, expect)
    print(f"wrote {JOIN_OUT}", flush=True)


def check_joined(key: str, expect: dict) -> None:
    wanted = ["variant_key", "ref", "alt", "in_genomes", "in_exomes", "AC", "AC_exome", "AF_nfe", "AF_afr", "AF_nfe_exome", "AF_afr_exome"]
    table = pq.read_table(JOIN_OUT, columns=wanted, filters=[("variant_key", "=", key)])
    if table.num_rows != 1:
        raise SystemExit(f"tier_a rows for {key}: {table.num_rows}")
    row = {name: table.column(name)[0].as_py() for name in wanted}
    if "ref" in expect.get("genomes", {}):
        if row["ref"] != expect["genomes"]["ref"] or row["alt"] != expect["genomes"]["alt"]:
            raise SystemExit(f"{key} alleles {row['ref']}:{row['alt']}")
    if key == "chr22:10510183:A:G":
        if row["AF_nfe"] is None or row["AF_afr"] is None or float(row["AF_nfe"]) == float(row["AF_afr"]):
            raise SystemExit(f"joined AF fields {row['AF_nfe']} {row['AF_afr']}")
        if row["AF_nfe_exome"] is not None or row["AF_afr_exome"] is not None or row["in_exomes"]:
            raise SystemExit("exome frequencies were filled for a site absent from the exome VCF")
    print(f"QC joined {key} ref={row['ref']} alt={row['alt']} AC={row['AC']} AC_exome={row['AC_exome']}", flush=True)


if __name__ == "__main__":
    main()
