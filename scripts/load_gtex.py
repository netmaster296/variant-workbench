# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load GTEx v8 median expression and QTL summaries into Parquet.

Streams tar members. Does not extract the archives and does not read the
sample-level gene TPM matrix.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import re
import shutil
import sys
import tarfile
import time
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from coords import bare_ensg, norm_rsid, variant_key

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
RAW = ROOT / "raw" / "gtex"
GENCODE = ROOT / "raw" / "gencode" / "gencode.v49.annotation.gtf.gz"
TRUTH = ROOT / "warehouse" / "qc" / "truth_set.tsv"
OUT = ROOT / "warehouse" / "gtex"
QC_PATH = ROOT / "warehouse" / "qc" / "gtex.md"

MEDIAN_GZ = RAW / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
# Sample-level TPM. Left unused on purpose.
SAMPLE_TPM = RAW / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
EQTL_TAR = RAW / "GTEx_Analysis_v8_eQTL.tar"
SQTL_TAR = RAW / "GTEx_Analysis_v8_sQTL.tar"

MEDIAN_OUT = OUT / "median_tpm.parquet"
EXPR_OUT = OUT / "gene_expression.parquet"
EGENE_OUT = OUT / "egenes.parquet"
SGENE_OUT = OUT / "sgenes.parquet"
EQTL_PAIR_OUT = OUT / "eqtl_pairs.parquet"
SQTL_PAIR_OUT = OUT / "sqtl_pairs.parquet"

EQTL_EGENE_SUFFIX = ".egenes.txt.gz"
EQTL_PAIR_SUFFIX = ".signif_variant_gene_pairs.txt.gz"
SQTL_GENE_SUFFIX = ".sgenes.txt.gz"
SQTL_PAIR_SUFFIX = ".sqtl_signifpairs.txt.gz"

EXPECTED_TISSUES = 54
EXPECTED_QTL_FILES = 49
FREE_FLOOR = 120 * 1024**3
QVAL_MAX = 0.05

EGENE_HEADER = [
    "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand",
    "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df",
    "variant_id", "tss_distance", "chr", "variant_pos", "ref", "alt",
    "num_alt_per_site", "rs_id_dbSNP151_GRCh38p7", "minor_allele_samples",
    "minor_allele_count", "maf", "ref_factor", "pval_nominal", "slope",
    "slope_se", "pval_perm", "pval_beta", "qval", "pval_nominal_threshold",
    "log2_aFC", "log2_aFC_lower", "log2_aFC_upper",
]
SGENE_HEADER = [
    "phenotype_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand",
    "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df",
    "variant_id", "tss_distance", "chr", "variant_pos", "ref", "alt",
    "num_alt_per_site", "rs_id_dbSNP151_GRCh38p7", "ma_samples", "ma_count",
    "maf", "ref_factor", "pval_nominal", "slope", "slope_se", "pval_perm",
    "pval_beta", "gene_id", "group_size", "qval", "pval_nominal_threshold",
]
EQTL_PAIR_HEADER = [
    "variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf",
    "pval_nominal", "slope", "slope_se", "pval_nominal_threshold",
    "min_pval_nominal", "pval_beta",
]
SQTL_PAIR_HEADER = [
    "variant_id", "phenotype_id", "tss_distance", "ma_samples", "ma_count", "maf",
    "pval_nominal", "slope", "slope_se", "pval_nominal_threshold",
    "min_pval_nominal", "pval_beta",
]

RS_COLUMN = "rs_id_dbSNP151_GRCh38p7"
_ACGTN = frozenset("ACGTN")
# Empty cells are missing log2 allelic fold changes. '.' is a missing rs id.
_NULL_STRINGS = ["", "NA", "na", "NaN", "nan", "."]
_VARIANT_ID_RE = re.compile(r"variant_id\s+(\S+)")


class Tally:
    def __init__(self, kind: str):
        self.kind = kind
        self.rows = 0
        self.files = 0
        self.tissues: list[str] = []
        self.per_tissue: list[tuple[str, int, int | None]] = []
        self.b38_fail: Counter[str] = Counter()
        self.b38_examples: list[str] = []
        self.disagree = 0
        self.disagree_examples: list[str] = []
        self.fallback = 0
        self.key_null = 0
        self.ensg_null = 0
        self.gene_id_null = 0
        self.rsid_null = 0
        self.rsid_unparsed = 0
        self.rsid_examples: list[str] = []
        self.q_le = 0
        self.q_null = 0
        self.q_border = 0
        self.slope_null = 0
        self.pval_nominal_null = 0
        self.pval_beta_null = 0
        self.dup_extra = 0
        self.phenotype_mismatch = 0
        self.phenotype_examples: list[str] = []
        self.phenotype_parse_fail = 0
        self.genes: set[str] = set()
        self.sig_genes: set[str] = set()
        self.gzip_bytes = 0
        self.text_bytes = 0


def log(msg: str) -> None:
    print(msg, flush=True)


def remember(bucket: list[str], text: str, limit: int = 12) -> None:
    if len(bucket) < limit:
        bucket.append(text)


def free_bytes() -> int:
    # ROOT's drive on Windows; the filesystem holding ROOT elsewhere.
    return shutil.disk_usage(ROOT.anchor if ROOT.drive else ROOT).free


def fmt_int(value: int | None) -> str:
    if value is None:
        return "NA"
    return f"{value:,}"


def rate(nulls: int | None, rows: int) -> str:
    if nulls is None or rows == 0:
        return "NA"
    return f"{nulls:,} / {rows:,} ({nulls / rows:.4%})"


def tissue_from_name(name: str) -> str:
    base = Path(name).name
    marker = ".v8."
    if marker not in base:
        raise SystemExit(f"member name has no .v8. marker: {name}")
    return base.split(marker)[0]


def parse_b38(variant_id: str | None) -> tuple[str | None, str | None]:
    """Return (variant_key, failure reason). Reason is None when the id parses."""
    if variant_id is None or variant_id == "":
        return None, "null"
    parts = variant_id.split("_")
    if len(parts) != 5 or parts[4] != "b38" or not parts[1].isdigit():
        return None, "shape"
    ref = parts[2].upper()
    alt = parts[3].upper()
    if (
        not ref
        or not alt
        or any(base not in _ACGTN for base in ref)
        or any(base not in _ACGTN for base in alt)
    ):
        return None, "allele"
    key = variant_key(parts[0], parts[1], ref, alt)
    if key is None:
        return None, "rejected"
    return key, None


def gene_from_phenotype(phenotype_id: str | None) -> str | None:
    if not phenotype_id:
        return None
    tail = phenotype_id.rsplit(":", 1)[-1]
    if bare_ensg(tail) is None:
        return None
    return tail


def column_types(header: list[str], ints: set[str], floats: set[str]) -> dict:
    types = {}
    for name in header:
        if name in ints:
            types[name] = pa.int64()
        elif name in floats:
            types[name] = pa.float64()
        else:
            types[name] = pa.string()
    return types


def field_type(name: str, ints: set[str], floats: set[str]) -> pa.DataType:
    if name in ints:
        return pa.int64()
    if name in floats:
        return pa.float64()
    return pa.string()


def list_tar(path: Path) -> list[tarfile.TarInfo]:
    with tarfile.open(path, "r") as tar:
        return [member for member in tar.getmembers() if member.isfile()]


def partition(members: list[tarfile.TarInfo], suffixes: tuple[str, ...]):
    groups = {suffix: [] for suffix in suffixes}
    ignored = []
    for member in members:
        matched = None
        for suffix in suffixes:
            if member.name.endswith(suffix):
                matched = suffix
                break
        if matched is None:
            ignored.append(member.name)
        else:
            groups[matched].append(member)
    for suffix in suffixes:
        groups[suffix].sort(key=lambda member: member.name)
    return groups, ignored


def compression_ratio(tar_path: Path, member: tarfile.TarInfo) -> float:
    with tarfile.open(tar_path, "r") as tar:
        raw = tar.extractfile(member).read()
    if not raw:
        return 1.0
    return len(gzip.decompress(raw)) / len(raw)


def read_member_table(tar: tarfile.TarFile, member: tarfile.TarInfo, types: dict) -> tuple[pa.Table, int, int]:
    raw = tar.extractfile(member).read()
    text = gzip.decompress(raw)
    table = pacsv.read_csv(
        io.BytesIO(text),
        read_options=pacsv.ReadOptions(use_threads=True),
        parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False),
        convert_options=pacsv.ConvertOptions(
            column_types=types,
            null_values=_NULL_STRINGS,
            strings_can_be_null=True,
        ),
    )
    return table, len(raw), len(text)


def write_parquet(path: Path, schema: pa.Schema, tables) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(path) + ".partial")
    if partial.exists():
        partial.unlink()
    rows = 0
    writer = pq.ParquetWriter(
        str(partial),
        schema,
        compression="zstd",
        compression_level=3,
    )
    try:
        for table in tables:
            writer.write_table(table, row_group_size=500_000)
            rows += table.num_rows
    except BaseException:
        writer.close()
        if partial.exists():
            partial.unlink()
        raise
    writer.close()
    partial.replace(path)
    return rows


def repeated_string(value: str, n: int) -> pa.Array:
    return pa.array(np.full(n, value, dtype=object), type=pa.string())


def map_ensg(gene_ids: list[str | None]) -> list[str | None]:
    mapping = {gene_id: bare_ensg(gene_id) for gene_id in set(gene_ids)}
    return [mapping[gene_id] for gene_id in gene_ids]


def self_check() -> None:
    if variant_key("chr1", 64764, "C", "T") != "chr1:64764:C:T":
        raise SystemExit("variant_key self-check failed")
    if bare_ensg("ENSG00000227232.5") != "ENSG00000227232":
        raise SystemExit("bare_ensg self-check failed")
    key, reason = parse_b38("chr1_64764_C_T_b38")
    if key != "chr1:64764:C:T" or reason is not None:
        raise SystemExit(f"b38 self-check failed: {key} {reason}")
    indel, reason = parse_b38("chr1_739465_TTTTG_T_b38")
    if indel != "chr1:739465:TTTTG:T" or reason is not None:
        raise SystemExit(f"indel self-check failed: {indel} {reason}")
    bad, reason = parse_b38("not_a_variant")
    if bad is not None or reason != "shape":
        raise SystemExit("malformed id self-check failed")
    parsed = gene_from_phenotype("chr1:497299:498399:clu_51878:ENSG00000237094.11")
    if parsed != "ENSG00000237094.11" or bare_ensg(parsed) != "ENSG00000237094":
        raise SystemExit("phenotype self-check failed")


def load_truth_case() -> dict:
    with TRUTH.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["case_id"] == "gtex_egene":
                break
        else:
            raise SystemExit("truth set has no gtex_egene row")
    match = _VARIANT_ID_RE.search(row["notes"])
    if not match:
        raise SystemExit("truth set gtex_egene notes have no variant_id")
    variant_id = match.group(1).rstrip(";,")
    tissue = tissue_from_name(row["source_locator"].split()[0])
    expected = variant_key(row["chrom_raw"], row["pos"], row["ref"], row["alt"])
    if expected != row["variant_key"] or expected != row["expected"]:
        raise SystemExit(
            f"truth set gtex_egene is inconsistent: {expected} {row['variant_key']} {row['expected']}"
        )
    row["parsed_variant_id"] = variant_id
    row["parsed_tissue"] = tissue
    row["parsed_key"] = expected
    return row


def load_gencode() -> tuple[set[str], set[str]]:
    versioned: set[str] = set()
    bare: set[str] = set()
    with gzip.open(GENCODE, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            tab1 = line.find("\t")
            if tab1 < 0:
                continue
            tab2 = line.find("\t", tab1 + 1)
            tab3 = line.find("\t", tab2 + 1)
            if tab2 < 0 or tab3 < 0 or line[tab2 + 1:tab3] != "gene":
                continue
            attr = line[line.rfind("\t") + 1:]
            marker = 'gene_id "'
            start = attr.find(marker)
            if start < 0:
                continue
            start += len(marker)
            end = attr.find('"', start)
            gene_id = attr[start:end]
            versioned.add(gene_id)
            ensg = bare_ensg(gene_id)
            if ensg:
                bare.add(ensg)
    return versioned, bare


def load_median() -> dict:
    with gzip.open(MEDIAN_GZ, "rt", encoding="utf-8", newline="") as handle:
        version_line = handle.readline().rstrip("\r\n")
        dim_line = handle.readline().rstrip("\r\n")
        header = handle.readline().rstrip("\r\n").split("\t")
        if not version_line.startswith("#"):
            raise SystemExit(f"median TPM is not GCT: {version_line!r}")
        if header[0] != "Name" or header[1] != "Description":
            raise SystemExit(f"unexpected GCT columns: {header[:2]}")
        dim_genes, dim_tissues = (int(part) for part in dim_line.split("\t"))
        tissues = header[2:]
        if len(tissues) != dim_tissues:
            raise SystemExit(
                f"GCT dimension says {dim_tissues} tissues but the header has {len(tissues)}"
            )
        values = np.empty((dim_genes, dim_tissues), dtype=np.float64)
        versions: list[str] = []
        ensgs: list[str | None] = []
        seen: set[str] = set()
        duplicate_version = 0
        for index, line in enumerate(handle):
            if index >= dim_genes:
                raise SystemExit("median TPM has more gene rows than the dimension line")
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != dim_tissues + 2:
                raise SystemExit(f"{cols[0]} has {len(cols)} columns, expected {dim_tissues + 2}")
            name = cols[0].strip()
            if name in seen:
                duplicate_version += 1
            seen.add(name)
            versions.append(name)
            ensgs.append(bare_ensg(name))
            for tissue_index, raw in enumerate(cols[2:]):
                try:
                    values[index, tissue_index] = float(raw)
                except ValueError as exc:
                    raise SystemExit(f"non-numeric TPM for {name} / {tissues[tissue_index]}: {raw!r}") from exc
        n_genes = index + 1 if dim_genes else 0
    if n_genes != dim_genes:
        raise SystemExit(f"median TPM rows {n_genes} != dimension {dim_genes}")
    if not np.isfinite(values).all():
        raise SystemExit("median TPM has non-finite values")
    if (values < 0).any():
        negative = int((values < 0).sum())
    else:
        negative = 0

    n_tissues = dim_tissues
    xmax = values.max(axis=1)
    max_index = values.argmax(axis=1)
    ties = ((values == xmax[:, None]).sum(axis=1) > 1) & (xmax > 0)
    ge1 = (values >= 1).sum(axis=1).astype(np.int32)
    ge10 = (values >= 10).sum(axis=1).astype(np.int32)
    tau = np.full(n_genes, np.nan, dtype=np.float64)
    positive = xmax > 0
    if n_tissues > 1 and positive.any():
        ratio = values[positive] / xmax[positive, None]
        tau[positive] = (1.0 - ratio).sum(axis=1) / (n_tissues - 1)
    finite_tau = tau[np.isfinite(tau)]
    max_tissues = [
        None if xmax[i] <= 0 else tissues[int(max_index[i])]
        for i in range(n_genes)
    ]
    bare_counts: Counter[str] = Counter(ensg for ensg in ensgs if ensg)
    duplicate_bare = sum(count - 1 for count in bare_counts.values() if count > 1)
    non_ensg = [version for version, ensg in zip(versions, ensgs) if ensg is None]

    long_n = n_genes * n_tissues
    median = pa.table(
        {
            "ensg_version": pa.array(np.repeat(np.array(versions, dtype=object), n_tissues), type=pa.string()),
            "ensg": pa.array(np.repeat(np.array(ensgs, dtype=object), n_tissues), type=pa.string()),
            "tissue": pa.array(np.tile(np.array(tissues, dtype=object), n_genes), type=pa.string()),
            "median_tpm": pa.array(values.reshape(-1), type=pa.float64()),
        },
        schema=pa.schema([
            ("ensg_version", pa.string()),
            ("ensg", pa.string()),
            ("tissue", pa.string()),
            ("median_tpm", pa.float64()),
        ]),
    )
    tau_mask = ~np.isfinite(tau)
    tau_filled = np.where(tau_mask, 0.0, tau)
    expression = pa.table(
        {
            "ensg": pa.array(ensgs, type=pa.string()),
            "ensg_version": pa.array(versions, type=pa.string()),
            "max_tissue": pa.array(max_tissues, type=pa.string()),
            "max_median_tpm": pa.array(xmax, type=pa.float64()),
            "n_tissues_tpm_ge_1": pa.array(ge1, type=pa.int32()),
            "n_tissues_tpm_ge_10": pa.array(ge10, type=pa.int32()),
            "tau": pa.array(tau_filled, mask=tau_mask, type=pa.float64()),
            "n_tissues": pa.array(np.full(n_genes, n_tissues, dtype=np.int32), type=pa.int32()),
        },
        schema=pa.schema([
            ("ensg", pa.string()),
            ("ensg_version", pa.string()),
            ("max_tissue", pa.string()),
            ("max_median_tpm", pa.float64()),
            ("n_tissues_tpm_ge_1", pa.int32()),
            ("n_tissues_tpm_ge_10", pa.int32()),
            ("tau", pa.float64()),
            ("n_tissues", pa.int32()),
        ]),
    )
    write_parquet(MEDIAN_OUT, median.schema, [median])
    write_parquet(EXPR_OUT, expression.schema, [expression])
    return {
        "version": version_line,
        "dim_line": dim_line,
        "n_genes": n_genes,
        "n_tissues": n_tissues,
        "tissues": tissues,
        "long_rows": long_n,
        "ensg_null": len(non_ensg),
        "non_ensg": non_ensg,
        "duplicate_version": duplicate_version,
        "duplicate_bare": duplicate_bare,
        "negative_cells": negative,
        "zero_max": int((xmax <= 0).sum()),
        "ties": int(ties.sum()),
        "tau_null": int(tau_mask.sum()),
        "tau_min": None if finite_tau.size == 0 else float(finite_tau.min()),
        "tau_max": None if finite_tau.size == 0 else float(finite_tau.max()),
        "genes": set(versions),
        "ge1_genes": int((ge1 > 0).sum()),
        "ge10_genes": int((ge10 > 0).sum()),
    }


def gene_schema(header: list[str], ints: set[str], floats: set[str]) -> pa.Schema:
    fields = [
        pa.field("tissue", pa.string()),
        pa.field("ensg", pa.string()),
        pa.field("variant_key", pa.string()),
        pa.field("rsid", pa.string()),
        pa.field("qval_le_0_05", pa.bool_()),
    ]
    for name in header:
        fields.append(pa.field(name, field_type(name, ints, floats)))
    return pa.schema(fields)


def transform_tested(table: pa.Table, tissue: str, tally: Tally, truth: dict, phenotype: bool) -> pa.Table:
    n = table.num_rows
    gene_ids = table.column("gene_id").to_pylist()
    variant_ids = table.column("variant_id").to_pylist()
    chroms = table.column("chr").to_pylist()
    positions = table.column("variant_pos").to_pylist()
    refs = table.column("ref").to_pylist()
    alts = table.column("alt").to_pylist()
    rs_raw = table.column(RS_COLUMN).to_pylist()
    phenotypes = table.column("phenotype_id").to_pylist() if phenotype else None
    ensg = map_ensg(gene_ids)
    keys: list[str | None] = []
    rsids: list[str | None] = []
    track = (not phenotype) and tissue == truth["parsed_tissue"]
    id_cache: dict[str | None, tuple[str | None, str | None]] = {}

    for i in range(n):
        column_key = variant_key(chroms[i], positions[i], refs[i], alts[i])
        variant_id = variant_ids[i]
        if variant_id not in id_cache:
            id_cache[variant_id] = parse_b38(variant_id)
        id_key, reason = id_cache[variant_id]
        if reason is not None:
            tally.b38_fail[reason] += 1
            remember(tally.b38_examples, f"{tissue} {variant_id!r} ({reason})")
        if column_key is None and id_key is not None:
            tally.fallback += 1
            final = id_key
        else:
            final = column_key
        if column_key is not None and id_key is not None and column_key != id_key:
            tally.disagree += 1
            remember(
                tally.disagree_examples,
                f"{tissue} {variant_id} columns={column_key} id={id_key}",
            )
        if final is None:
            tally.key_null += 1
        keys.append(final)
        rsid = norm_rsid(rs_raw[i])
        if rsid is None:
            tally.rsid_null += 1
            if rs_raw[i] not in (None, "", ".", "NA", "na"):
                tally.rsid_unparsed += 1
                if rs_raw[i] not in tally.rsid_examples and len(tally.rsid_examples) < 8:
                    tally.rsid_examples.append(str(rs_raw[i]))
        rsids.append(rsid)
        if track and gene_ids[i] == truth["gene"] and variant_id == truth["parsed_variant_id"]:
            truth["seen_during_load"] = True
            truth["load_index"] = i
            truth["load_is_first"] = i == 0
            truth["load_column_key"] = column_key
            truth["load_id_key"] = id_key
            truth["load_final_key"] = final
            truth["load_chr"] = chroms[i]
            truth["load_pos"] = positions[i]
            truth["load_ref"] = refs[i]
            truth["load_alt"] = alts[i]

    if phenotype:
        mismatch = 0
        for phenotype_id, gene_id in zip(phenotypes, gene_ids):
            tail = phenotype_id.rsplit(":", 1)[-1] if phenotype_id else None
            if tail != gene_id:
                mismatch += 1
                remember(
                    tally.phenotype_examples,
                    f"{tissue} phenotype={phenotype_id} gene_id={gene_id}",
                )
        tally.phenotype_mismatch += mismatch
        tally.dup_extra += n - len(set(phenotypes))
    else:
        tally.dup_extra += n - len(set(gene_ids))

    qvals = table.column("qval").to_numpy(zero_copy_only=False)
    finite = np.isfinite(qvals)
    flag = finite & (qvals <= QVAL_MAX)
    tally.q_le += int(flag.sum())
    tally.q_null += int((~finite).sum())
    tally.q_border += int(np.sum((qvals > QVAL_MAX) & (qvals <= QVAL_MAX + 1e-8)))
    tally.ensg_null += sum(value is None for value in ensg)
    tally.gene_id_null += sum(value is None for value in gene_ids)
    tally.slope_null += table.column("slope").null_count
    tally.pval_nominal_null += table.column("pval_nominal").null_count
    tally.pval_beta_null += table.column("pval_beta").null_count
    tally.genes.update(gene_id for gene_id in gene_ids if gene_id)
    tally.sig_genes.update(gene_ids[i] for i in range(n) if flag[i] and gene_ids[i])
    tally.rows += n
    tally.files += 1
    tally.tissues.append(tissue)
    tally.per_tissue.append((tissue, n, int(flag.sum())))

    columns = {
        "tissue": repeated_string(tissue, n),
        "ensg": pa.array(ensg, type=pa.string()),
        "variant_key": pa.array(keys, type=pa.string()),
        "rsid": pa.array(rsids, type=pa.string()),
        "qval_le_0_05": pa.array(flag, type=pa.bool_()),
    }
    for name in table.column_names:
        columns[name] = table.column(name)
    return pa.Table.from_arrays([columns[field.name] for field in tally.schema], schema=tally.schema)


def load_tested(tar_path: Path, members: list[tarfile.TarInfo], dest: Path, header: list[str], ints: set[str], floats: set[str], tally: Tally, truth: dict, phenotype: bool) -> None:
    types = column_types(header, ints, floats)
    tally.schema = gene_schema(header, ints, floats)
    seen: set[str] = set()

    def batches():
        with tarfile.open(tar_path, "r") as tar:
            for member in members:
                tissue = tissue_from_name(member.name)
                if tissue in seen:
                    raise SystemExit(f"duplicate tissue {tissue} in {tar_path.name}")
                seen.add(tissue)
                table, gzip_n, text_n = read_member_table(tar, member, types)
                if table.column_names != header:
                    raise SystemExit(f"{member.name} header {table.column_names} != expected")
                tally.gzip_bytes += gzip_n
                tally.text_bytes += text_n
                out = transform_tested(table, tissue, tally, truth, phenotype)
                log(f"{tally.kind} {tissue} rows={table.num_rows:,} cumulative={tally.rows:,}")
                yield out

    written = write_parquet(dest, tally.schema, batches())
    if written != tally.rows:
        raise SystemExit(f"{dest.name} wrote {written} rows, tally {tally.rows}")


def pair_schema(phenotype: bool) -> pa.Schema:
    names = [
        "tissue", "variant_id",
        "phenotype_id", "gene_id", "ensg", "variant_key",
        "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal",
        "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta",
    ]
    if not phenotype:
        names.remove("phenotype_id")
    types = {
        "tissue": pa.string(),
        "variant_id": pa.string(),
        "phenotype_id": pa.string(),
        "gene_id": pa.string(),
        "ensg": pa.string(),
        "variant_key": pa.string(),
        "tss_distance": pa.int64(),
        "ma_samples": pa.int64(),
        "ma_count": pa.int64(),
        "maf": pa.float64(),
        "pval_nominal": pa.float64(),
        "slope": pa.float64(),
        "slope_se": pa.float64(),
        "pval_nominal_threshold": pa.float64(),
        "min_pval_nominal": pa.float64(),
        "pval_beta": pa.float64(),
    }
    return pa.schema([pa.field(name, types[name]) for name in names])


def transform_pairs(table: pa.Table, tissue: str, tally: Tally, schema: pa.Schema, phenotype: bool) -> pa.Table:
    n = table.num_rows
    variant_ids = table.column("variant_id").to_pylist()
    if phenotype:
        phenotypes = table.column("phenotype_id").to_pylist()
        gene_ids = []
        for phenotype_id in phenotypes:
            gene_id = gene_from_phenotype(phenotype_id)
            if gene_id is None:
                tally.phenotype_parse_fail += 1
                remember(tally.phenotype_examples, f"{tissue} {phenotype_id}")
            gene_ids.append(gene_id)
        gene_col = pa.array(gene_ids, type=pa.string())
    else:
        gene_ids = table.column("gene_id").to_pylist()
        gene_col = table.column("gene_id")
    ensg = map_ensg(gene_ids)
    keys: list[str | None] = []
    cache: dict[str | None, tuple[str | None, str | None]] = {}
    for variant_id in variant_ids:
        if variant_id not in cache:
            cache[variant_id] = parse_b38(variant_id)
        key, reason = cache[variant_id]
        if reason is not None:
            tally.b38_fail[reason] += 1
            remember(tally.b38_examples, f"{tissue} {variant_id!r} ({reason})")
        if key is None:
            tally.key_null += 1
        keys.append(key)
    tally.ensg_null += sum(value is None for value in ensg)
    tally.gene_id_null += sum(value is None for value in gene_ids)
    tally.slope_null += table.column("slope").null_count
    tally.pval_nominal_null += table.column("pval_nominal").null_count
    tally.pval_beta_null += table.column("pval_beta").null_count
    tally.genes.update(gene_id for gene_id in gene_ids if gene_id)
    tally.rows += n
    tally.files += 1
    tally.tissues.append(tissue)
    tally.per_tissue.append((tissue, n, None))

    columns = {
        "tissue": repeated_string(tissue, n),
        "variant_id": table.column("variant_id"),
        "gene_id": gene_col,
        "ensg": pa.array(ensg, type=pa.string()),
        "variant_key": pa.array(keys, type=pa.string()),
        "tss_distance": table.column("tss_distance"),
        "ma_samples": table.column("ma_samples"),
        "ma_count": table.column("ma_count"),
        "maf": table.column("maf"),
        "pval_nominal": table.column("pval_nominal"),
        "slope": table.column("slope"),
        "slope_se": table.column("slope_se"),
        "pval_nominal_threshold": table.column("pval_nominal_threshold"),
        "min_pval_nominal": table.column("min_pval_nominal"),
        "pval_beta": table.column("pval_beta"),
    }
    if phenotype:
        columns["phenotype_id"] = table.column("phenotype_id")
    return pa.Table.from_arrays([columns[field.name] for field in schema], schema=schema)


def load_pairs(tar_path: Path, members: list[tarfile.TarInfo], dest: Path, header: list[str], ints: set[str], floats: set[str], tally: Tally, phenotype: bool) -> None:
    types = column_types(header, ints, floats)
    schema = pair_schema(phenotype)
    tally.schema = schema
    seen: set[str] = set()

    def batches():
        with tarfile.open(tar_path, "r") as tar:
            for member in members:
                tissue = tissue_from_name(member.name)
                if tissue in seen:
                    raise SystemExit(f"duplicate tissue {tissue} in {tar_path.name}")
                seen.add(tissue)
                table, gzip_n, text_n = read_member_table(tar, member, types)
                if table.column_names != header:
                    raise SystemExit(f"{member.name} header {table.column_names} != expected")
                tally.gzip_bytes += gzip_n
                tally.text_bytes += text_n
                out = transform_pairs(table, tissue, tally, schema, phenotype)
                log(f"{tally.kind} {tissue} rows={table.num_rows:,} cumulative={tally.rows:,}")
                yield out

    written = write_parquet(dest, schema, batches())
    if written != tally.rows:
        raise SystemExit(f"{dest.name} wrote {written} rows, tally {tally.rows}")


def column_nulls(path: Path, name: str) -> int | None:
    parquet = pq.ParquetFile(path)
    index = parquet.schema_arrow.get_field_index(name)
    if index < 0:
        return None
    total = 0
    for group in range(parquet.metadata.num_row_groups):
        stats = parquet.metadata.row_group(group).column(index).statistics
        if stats is None or stats.null_count is None:
            return None
        total += stats.null_count
    return total


def file_facts(path: Path) -> dict:
    parquet = pq.ParquetFile(path)
    return {
        "rows": parquet.metadata.num_rows,
        "bytes": path.stat().st_size,
        "row_groups": parquet.metadata.num_row_groups,
        "columns": parquet.schema_arrow.names,
    }


def classify_genes(ids: set[str], versioned: set[str], bare: set[str]) -> dict:
    version_only = []
    gone = []
    not_ensg = []
    bare_ids = set()
    for gene_id in ids:
        ensg = bare_ensg(gene_id)
        if ensg:
            bare_ids.add(ensg)
        if gene_id in versioned:
            continue
        if ensg is None:
            not_ensg.append(gene_id)
        elif ensg in bare:
            version_only.append(gene_id)
        else:
            gone.append(gene_id)
    bare_miss = sorted(bare_ids - bare)
    return {
        "versioned": len(ids),
        "bare": len(bare_ids),
        "version_only": sorted(version_only),
        "gone": sorted(gone),
        "not_ensg": sorted(not_ensg),
        "bare_miss": bare_miss,
    }


def confirm_truth(truth: dict) -> dict:
    dataset = ds.dataset(str(EGENE_OUT))
    table = dataset.to_table(
        filter=(
            (pc.field("tissue") == truth["parsed_tissue"])
            & (pc.field("gene_id") == truth["gene"])
            & (pc.field("variant_id") == truth["parsed_variant_id"])
        )
    )
    result = {
        "rows": table.num_rows,
        "pass": False,
        "detail": "",
    }
    if table.num_rows != 1:
        result["detail"] = f"expected 1 row on disk, found {table.num_rows}"
        return result
    row = table.to_pylist()[0]
    id_key, reason = parse_b38(row["variant_id"])
    column_key = variant_key(row["chr"], row["variant_pos"], row["ref"], row["alt"])
    checks = {
        "tissue": row["tissue"] == truth["parsed_tissue"],
        "gene_id": row["gene_id"] == truth["gene"],
        "ensg": row["ensg"] == bare_ensg(truth["gene"]),
        "variant_id": row["variant_id"] == truth["parsed_variant_id"],
        "variant_key": row["variant_key"] == truth["parsed_key"],
        "column_key": column_key == truth["parsed_key"],
        "id_key": id_key == truth["parsed_key"] and reason is None,
        "keys_agree": column_key == id_key == row["variant_key"],
        "chr": row["chr"] == truth["chrom_raw"],
        "pos": row["variant_pos"] == int(truth["pos"]),
        "ref": row["ref"] == truth["ref"],
        "alt": row["alt"] == truth["alt"],
        "seen_during_load": truth.get("seen_during_load") is True,
        "first_data_row": truth.get("load_is_first") is True,
    }
    result["row"] = row
    result["column_key"] = column_key
    result["id_key"] = id_key
    result["checks"] = checks
    result["pass"] = all(checks.values())
    failed = [name for name, ok in checks.items() if not ok]
    result["detail"] = "all checks passed" if result["pass"] else "failed: " + ", ".join(failed)
    return result


def par_note(names: list[str]) -> str:
    if names and all(name.endswith("_PAR_Y") for name in names):
        return (
            f"All {len(names)} are Names ending in `_PAR_Y`. `bare_ensg` does not accept that suffix, "
            "so `ensg` is null and `ensg_version` keeps the Name. No row was dropped."
        )
    return f"Examples: {examples(names)}."


def examples(values: list[str], limit: int = 8) -> str:
    if not values:
        return "none"
    shown = values[:limit]
    suffix = "" if len(values) <= limit else f" (+{len(values) - limit:,} more)"
    return ", ".join(shown) + suffix


def tissue_table(rows: list[tuple[str, int, int | None]], sig: bool) -> str:
    if sig:
        lines = ["| Tissue | Rows | qval <= 0.05 |", "|---|---:|---:|"]
        for tissue, n, n_sig in rows:
            lines.append(f"| {tissue} | {n:,} | {n_sig:,} |")
    else:
        lines = ["| Tissue | Rows |", "|---|---:|"]
        for tissue, n, _sig in rows:
            lines.append(f"| {tissue} | {n:,} |")
    return "\n".join(lines)


def b38_line(tally: Tally) -> str:
    total = sum(tally.b38_fail.values())
    if total == 0:
        return "0"
    parts = ", ".join(f"{reason} {count:,}" for reason, count in sorted(tally.b38_fail.items()))
    return f"{total:,} ({parts})"


def render(ctx: dict) -> str:
    median = ctx["median"]
    truth = ctx["truth_result"]
    free_before = ctx["free_before"]
    free_after = ctx["free_after"]
    lines = [
        "# GTEx lane E",
        "",
        f"Loaded {date.today().isoformat()}. Assembly GRCh38. Variant keys and bare ENSG ids come from `scripts/coords.py`.",
        "",
        "## Free space",
        "",
        f"- Before: {fmt_int(free_before)} bytes.",
        f"- After: {fmt_int(free_after)} bytes.",
        f"- Floor that would have blocked pair output: {fmt_int(FREE_FLOOR)} bytes (120 GiB).",
        f"- Pair text estimate, used only for that floor check: {fmt_int(ctx['pair_estimate'])} bytes.",
        f"- Estimate method: gzip size of both pair-member sets times the larger measured uncompressed/gzip ratio ({ctx['pair_ratio']:.2f}), times 1.5. The ratio came from `{ctx['ratio_member']}`.",
        f"- Decision: {ctx['pair_decision']}",
        "",
        "## Inputs",
        "",
        f"- Median TPM: `{MEDIAN_GZ.name}`, {fmt_int(MEDIAN_GZ.stat().st_size)} bytes. GCT `{median['version']}`, dimension line `{median['dim_line']}`.",
        f"- eQTL tar: `{EQTL_TAR.name}`, {fmt_int(EQTL_TAR.stat().st_size)} bytes, {fmt_int(ctx['eqtl_members'])} file members. Ignored members: {ctx['eqtl_ignored'] or 'none'}.",
        f"- sQTL tar: `{SQTL_TAR.name}`, {fmt_int(SQTL_TAR.stat().st_size)} bytes, {fmt_int(ctx['sqtl_members'])} file members. Ignored members: {ctx['sqtl_ignored'] or 'none'}.",
        "- Members were read from the tars and parsed in memory. The tars were not extracted to disk.",
        "- Cells that are empty, `NA`, `NaN`, or `.` are null. A `.` rs id therefore has a null raw value and a null `rsid`.",
        "",
        "## Sample-level TPM",
        "",
        f"`{SAMPLE_TPM.name}` ({fmt_int(ctx['sample_bytes'])} bytes) was not opened. It is the sample-level gene TPM matrix. This lane uses the median TPM file instead.",
        "",
        "## Median TPM",
        "",
        f"- Tissue columns: {median['n_tissues']}. The file's dimension line says {median['n_tissues']}. GTEx v8 is described as 54 tissues. Match: {median['n_tissues'] == EXPECTED_TISSUES}.",
        f"- Genes: {fmt_int(median['n_genes'])}. Long rows written: {fmt_int(median['long_rows'])} ({fmt_int(median['n_genes'])} x {median['n_tissues']}).",
        f"- `ensg` null (Name was not an ENSG id): {rate(median['ensg_null'], median['n_genes'])}. {par_note(median['non_ensg'])}",
        f"- Duplicate versioned ids: {fmt_int(median['duplicate_version'])}. Extra rows sharing a bare ENSG: {fmt_int(median['duplicate_bare'])}.",
        f"- Negative TPM cells: {fmt_int(median['negative_cells'])}. Non-finite cells: 0.",
        "- Tissue labels are the GCT column names, in file order. QTL tissues use a different label (`Adipose_Subcutaneous`, not `Adipose - Subcutaneous`). This lane does not rename either set.",
        "",
        "GCT tissues:",
        "",
        ", ".join(median["tissues"]),
        "",
        "## Gene expression",
        "",
        "One row per GCT gene. `ensg` is `bare_ensg(ensg_version)`.",
        "",
        "Tau, over the n tissue columns in the median file, only when x_max > 0:",
        "",
        "```text",
        "tau = sum_i(1 - x_i / x_max) / (n - 1)",
        "```",
        "",
        "x_i is median TPM in tissue i. x_max is the maximum of those values. n is the tissue count in the file, including tissues where the gene is 0. No pseudocount is added. Tau is null when x_max is 0, because the ratio is undefined. Under this definition a single-tissue gene has tau 1 and a gene with the same median in every tissue has tau 0.",
        "",
        f"- n = {median['n_tissues']}.",
        f"- Tau null: {fmt_int(median['tau_null'])}. Genes whose maximum is 0: {fmt_int(median['zero_max'])}. For those, `max_tissue` is null and `max_median_tpm` is 0.",
        f"- Finite tau range: {median['tau_min']} to {median['tau_max']}.",
        f"- Genes with median TPM >= 1 in at least one tissue: {fmt_int(median['ge1_genes'])}. Genes with median TPM >= 10 in at least one tissue: {fmt_int(median['ge10_genes'])}.",
        f"- Ties for the maximum, counted only when x_max > 0: {fmt_int(median['ties'])}. `max_tissue` is the leftmost GCT column among the tied tissues.",
        "",
        "## eGenes",
        "",
        f"- Files found: {ctx['egene'].files}. Expected: {EXPECTED_QTL_FILES}. Match: {ctx['egene'].files == EXPECTED_QTL_FILES}.",
        "- Source header is the same in every `*.egenes.txt.gz` member. The effect column is `slope`. `slope_se` is its standard error. There is no effect column named `beta`. `beta_shape1` and `beta_shape2` are beta-distribution parameters for the permutation approximation, not effect sizes.",
        "- `pval_nominal` is the lead variant's nominal p-value. `pval_perm` is the permutation p-value. `pval_beta` is the beta-approximated p-value. `qval` is the gene-level FDR.",
        f"- The file has one row per tested gene per tissue, including genes that are not significant. `qval_le_0_05` is true when `qval` <= {QVAL_MAX}. That is the eGene call. Rows with a larger q-value are kept.",
        f"- Rows: {fmt_int(ctx['egene'].rows)}. qval <= 0.05: {fmt_int(ctx['egene'].q_le)}. qval null: {fmt_int(ctx['egene'].q_null)}. Values in (0.05, 0.05+1e-8]: {fmt_int(ctx['egene'].q_border)}.",
        f"- Distinct versioned genes: {fmt_int(len(ctx['egene'].genes))}. Distinct versioned genes with qval <= 0.05 in at least one tissue: {fmt_int(len(ctx['egene'].sig_genes))}.",
        f"- Duplicate gene ids within a tissue (extra rows): {fmt_int(ctx['egene'].dup_extra)}.",
        f"- `variant_key` null: {rate(ctx['nulls']['egenes']['variant_key'], ctx['egene'].rows)}. Loader tally of null keys: {fmt_int(ctx['egene'].key_null)}.",
        f"- `ensg` null: {rate(ctx['nulls']['egenes']['ensg'], ctx['egene'].rows)}.",
        f"- `rsid` null: {rate(ctx['nulls']['egenes']['rsid'], ctx['egene'].rows)}. Raw rs text present but not a single id `norm_rsid` accepts: {fmt_int(ctx['egene'].rsid_unparsed)}. Examples: {examples(ctx['egene'].rsid_examples)}. A comma-separated pair stays in the raw column; `rsid` stays null.",
        f"- `slope` null: {rate(ctx['nulls']['egenes']['slope'], ctx['egene'].rows)}. `pval_nominal` null: {rate(ctx['nulls']['egenes']['pval_nominal'], ctx['egene'].rows)}. `qval` null: {rate(ctx['nulls']['egenes']['qval'], ctx['egene'].rows)}.",
        f"- `log2_aFC` null: {rate(ctx['nulls']['egenes']['log2_aFC'], ctx['egene'].rows)}. `log2_aFC_lower` null: {rate(ctx['nulls']['egenes']['log2_aFC_lower'], ctx['egene'].rows)}. `log2_aFC_upper` null: {rate(ctx['nulls']['egenes']['log2_aFC_upper'], ctx['egene'].rows)}. Empty allelic-fold-change cells in the source are null here, not zeros.",
        f"- Split-column key missing and id key used instead: {fmt_int(ctx['egene'].fallback)}.",
        f"- Non-null disagreements between the split-column key and the b38 id key: {fmt_int(ctx['egene'].disagree)}. Examples: {examples(ctx['egene'].disagree_examples)}.",
        f"- Variant ids that failed the b38 parse: {b38_line(ctx['egene'])}. Examples: {examples(ctx['egene'].b38_examples)}.",
        "- The stored `variant_key` prefers `chr`, `variant_pos`, `ref`, and `alt`. The id is the check, and the fallback only when those columns do not make a key.",
        "",
        tissue_table(ctx["egene"].per_tissue, sig=True),
        "",
        "## sGenes",
        "",
        f"- Files found: {ctx['sgene'].files}. Expected: {EXPECTED_QTL_FILES}. Match: {ctx['sgene'].files == EXPECTED_QTL_FILES}.",
        "- Source header is the same in every `*.sgenes.txt.gz` member. The effect column is again `slope`, not `beta`.",
        "- Grain is one tested splicing phenotype (intron cluster) per tissue, not one row per gene. `phenotype_id` looks like `chrom:start:end:cluster:ENSG`. `gene_id` is a source column.",
        f"- A gene is an sGene in a tissue when any phenotype row for that gene and tissue has qval <= {QVAL_MAX}. `qval_le_0_05` is that phenotype-level flag. It is not collapsed to one row per gene.",
        f"- Rows: {fmt_int(ctx['sgene'].rows)}. qval <= 0.05: {fmt_int(ctx['sgene'].q_le)}. qval null: {fmt_int(ctx['sgene'].q_null)}. Values in (0.05, 0.05+1e-8]: {fmt_int(ctx['sgene'].q_border)}.",
        f"- Distinct versioned genes: {fmt_int(len(ctx['sgene'].genes))}. Distinct versioned genes with at least one phenotype at qval <= 0.05: {fmt_int(len(ctx['sgene'].sig_genes))}.",
        f"- Duplicate phenotype ids within a tissue (extra rows): {fmt_int(ctx['sgene'].dup_extra)}.",
        f"- Phenotype tail differing from the `gene_id` column: {fmt_int(ctx['sgene'].phenotype_mismatch)}. Examples: {examples(ctx['sgene'].phenotype_examples)}.",
        f"- `variant_key` null: {rate(ctx['nulls']['sgenes']['variant_key'], ctx['sgene'].rows)}. Loader tally: {fmt_int(ctx['sgene'].key_null)}.",
        f"- `ensg` null: {rate(ctx['nulls']['sgenes']['ensg'], ctx['sgene'].rows)}.",
        f"- `rsid` null: {rate(ctx['nulls']['sgenes']['rsid'], ctx['sgene'].rows)}. Raw rs text present but not a single id `norm_rsid` accepts: {fmt_int(ctx['sgene'].rsid_unparsed)}. Examples: {examples(ctx['sgene'].rsid_examples)}. A comma-separated pair stays in the raw column; `rsid` stays null.",
        f"- `slope` null: {rate(ctx['nulls']['sgenes']['slope'], ctx['sgene'].rows)}. `pval_nominal` null: {rate(ctx['nulls']['sgenes']['pval_nominal'], ctx['sgene'].rows)}. `qval` null: {rate(ctx['nulls']['sgenes']['qval'], ctx['sgene'].rows)}.",
        f"- Column-key fallback to the id: {fmt_int(ctx['sgene'].fallback)}. Non-null key disagreements: {fmt_int(ctx['sgene'].disagree)}. Examples: {examples(ctx['sgene'].disagree_examples)}.",
        f"- Variant ids that failed the b38 parse: {b38_line(ctx['sgene'])}. Examples: {examples(ctx['sgene'].b38_examples)}.",
        "- Minor-allele columns in this file are `ma_samples` and `ma_count`. The eGene file calls them `minor_allele_samples` and `minor_allele_count`. Both names are kept as published.",
        "",
        tissue_table(ctx["sgene"].per_tissue, sig=True),
        "",
        "## eQTL pairs",
        "",
        f"- Files found: {ctx['eqtl_pairs'].files}. Expected: {EXPECTED_QTL_FILES}. Match: {ctx['eqtl_pairs'].files == EXPECTED_QTL_FILES}.",
        "- Member suffix `*.signif_variant_gene_pairs.txt.gz`. These rows are the significant pairs GTEx already filtered. This loader does not filter them again.",
        "- The pair header has no `chr`, `variant_pos`, `ref`, or `alt`. `variant_key` is parsed from `variant_id` (`chr_pos_ref_alt_b38`).",
        f"- Rows: {fmt_int(ctx['eqtl_pairs'].rows)}. Distinct versioned genes: {fmt_int(len(ctx['eqtl_pairs'].genes))}.",
        f"- `variant_key` null: {rate(ctx['nulls']['eqtl_pairs']['variant_key'], ctx['eqtl_pairs'].rows)}. Loader tally: {fmt_int(ctx['eqtl_pairs'].key_null)}.",
        f"- `ensg` null: {rate(ctx['nulls']['eqtl_pairs']['ensg'], ctx['eqtl_pairs'].rows)}. `slope` null: {rate(ctx['nulls']['eqtl_pairs']['slope'], ctx['eqtl_pairs'].rows)}. `pval_nominal` null: {rate(ctx['nulls']['eqtl_pairs']['pval_nominal'], ctx['eqtl_pairs'].rows)}. `pval_beta` null: {rate(ctx['nulls']['eqtl_pairs']['pval_beta'], ctx['eqtl_pairs'].rows)}.",
        f"- Variant ids that failed the b38 parse: {b38_line(ctx['eqtl_pairs'])}. Examples: {examples(ctx['eqtl_pairs'].b38_examples)}.",
        f"- Versioned pair genes absent from the eGene file: {fmt_int(len(ctx['eqtl_pairs'].genes - ctx['egene'].genes))}.",
        "",
        tissue_table(ctx["eqtl_pairs"].per_tissue, sig=False),
        "",
        "## sQTL pairs",
        "",
        f"- Files found: {ctx['sqtl_pairs'].files}. Expected: {EXPECTED_QTL_FILES}. Match: {ctx['sqtl_pairs'].files == EXPECTED_QTL_FILES}.",
        "- Member suffix `*.sqtl_signifpairs.txt.gz`. This is not the eQTL pair name. The header uses `phenotype_id` instead of `gene_id` and has no split allele columns.",
        "- `gene_id` in the Parquet is the last colon-separated field of `phenotype_id` when that field is an ENSG id. `ensg` is the bare id. Rows whose phenotype does not end in an ENSG are kept, with those two columns null.",
        f"- Phenotype ids that did not yield an ENSG: {fmt_int(ctx['sqtl_pairs'].phenotype_parse_fail)}. Examples: {examples(ctx['sqtl_pairs'].phenotype_examples)}.",
        f"- Rows: {fmt_int(ctx['sqtl_pairs'].rows)}. Distinct versioned genes: {fmt_int(len(ctx['sqtl_pairs'].genes))}.",
        f"- `variant_key` null: {rate(ctx['nulls']['sqtl_pairs']['variant_key'], ctx['sqtl_pairs'].rows)}. Loader tally: {fmt_int(ctx['sqtl_pairs'].key_null)}.",
        f"- `ensg` null: {rate(ctx['nulls']['sqtl_pairs']['ensg'], ctx['sqtl_pairs'].rows)}. `gene_id` null: {rate(ctx['nulls']['sqtl_pairs']['gene_id'], ctx['sqtl_pairs'].rows)}. `slope` null: {rate(ctx['nulls']['sqtl_pairs']['slope'], ctx['sqtl_pairs'].rows)}. `pval_nominal` null: {rate(ctx['nulls']['sqtl_pairs']['pval_nominal'], ctx['sqtl_pairs'].rows)}. `pval_beta` null: {rate(ctx['nulls']['sqtl_pairs']['pval_beta'], ctx['sqtl_pairs'].rows)}.",
        f"- Variant ids that failed the b38 parse: {b38_line(ctx['sqtl_pairs'])}. Examples: {examples(ctx['sqtl_pairs'].b38_examples)}.",
        f"- Versioned pair genes absent from the sGene file: {fmt_int(len(ctx['sqtl_pairs'].genes - ctx['sgene'].genes))}.",
        "",
        tissue_table(ctx["sqtl_pairs"].per_tissue, sig=False),
        "",
        "## Variant id parse, all QTL tables",
        "",
        "Accepted form: five underscore fields, `chrom`, `pos`, `ref`, `alt`, `b38`, with pos an integer and alleles A/C/G/T/N. The key is then `variant_key` from `scripts/coords.py`, so the chromosome is canonical and the alleles are uppercase.",
        "",
        f"- eGenes: {b38_line(ctx['egene'])}.",
        f"- sGenes: {b38_line(ctx['sgene'])}.",
        f"- eQTL pairs: {b38_line(ctx['eqtl_pairs'])}.",
        f"- sQTL pairs: {b38_line(ctx['sqtl_pairs'])}.",
        f"- Total failed ids: {fmt_int(ctx['b38_total'])}.",
        "",
        "## GENCODE v49",
        "",
        f"Gene features in `{GENCODE.name}`: {fmt_int(ctx['gencode_versioned'])} versioned ids, {fmt_int(ctx['gencode_bare'])} bare ENSG ids.",
        "",
        "Join rule is the bare ENSG. A version-only difference still joins. A bare id that is absent from GENCODE v49 does not. No GTEx row was dropped for either reason.",
        "",
        "| Source | Versioned ids | Bare ids | Version-only difference | Bare id absent | Not an ENSG |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, status in ctx["joins"]:
        lines.append(
            f"| {label} | {len(status['ids']):,} | {status['class']['bare']:,} | {len(status['class']['version_only']):,} | {len(status['class']['bare_miss']):,} | {len(status['class']['not_ensg']):,} |"
        )
    lines.extend([
        "",
        "Examples of bare ENSG ids absent from GENCODE v49 (these will not join):",
        "",
    ])
    for label, status in ctx["joins"]:
        lines.append(f"- {label}: {examples(status['class']['bare_miss'])}.")
    lines.extend([
        "",
        "Examples of versioned ids whose bare id is in GENCODE v49 under a different version (these still join on the bare id):",
        "",
    ])
    for label, status in ctx["joins"]:
        lines.append(f"- {label}: {examples(status['class']['version_only'])}.")
    lines.extend([
        "",
        "## Truth set `gtex_egene`",
        "",
        f"- Source locator: {ctx['truth']['source_locator']}.",
        f"- Tissue `{ctx['truth']['parsed_tissue']}`, gene `{ctx['truth']['gene']}`, variant id `{ctx['truth']['parsed_variant_id']}`.",
        f"- Expected key `{ctx['truth']['parsed_key']}` from both the truth-set columns and `variant_key`.",
        f"- Seen while reading the eGene member: {ctx['truth'].get('seen_during_load')}. First data row of that member: {ctx['truth'].get('load_is_first')}. Load index: {ctx['truth'].get('load_index')}.",
        f"- Key from the split columns at load: {ctx['truth'].get('load_column_key')}. Key from the variant id at load: {ctx['truth'].get('load_id_key')}. Stored key: {ctx['truth'].get('load_final_key')}.",
        f"- Re-read from `egenes.parquet`: {truth['rows']} row(s). {truth['detail']}.",
        f"- Result: {'PASS' if truth['pass'] else 'FAIL'}.",
        "",
        "## Output files",
        "",
        "| File | Rows | Bytes | Row groups |",
        "|---|---:|---:|---:|",
    ])
    for path in (MEDIAN_OUT, EXPR_OUT, EGENE_OUT, SGENE_OUT, EQTL_PAIR_OUT, SQTL_PAIR_OUT):
        facts = ctx["facts"][path.name]
        lines.append(f"| `{path.name}` | {facts['rows']:,} | {facts['bytes']:,} | {facts['row_groups']:,} |")
    lines.extend([
        "",
        "Columns:",
        "",
    ])
    for path in (MEDIAN_OUT, EXPR_OUT, EGENE_OUT, SGENE_OUT, EQTL_PAIR_OUT, SQTL_PAIR_OUT):
        facts = ctx["facts"][path.name]
        lines.append(f"- `{path.name}`: {', '.join(facts['columns'])}.")
    lines.extend([
        "",
        f"Elapsed: {ctx['elapsed']:.1f} seconds. Python {sys.version.split()[0]}. pyarrow {pa.__version__}. bcftools was not used.",
        "",
    ])
    return "\n".join(lines)


def null_map(path: Path, names: list[str], rows: int) -> dict[str, int | None]:
    return {name: column_nulls(path, name) for name in names}


def require_inputs() -> None:
    missing = [str(path) for path in (MEDIAN_GZ, EQTL_TAR, SQTL_TAR, GENCODE, TRUTH) if not path.exists()]
    if missing:
        raise SystemExit("missing inputs: " + ", ".join(missing))


def main() -> int:
    started = time.time()
    self_check()
    require_inputs()
    truth = load_truth_case()
    free_before = free_bytes()
    log(f"free bytes {free_before:,}")

    eqtl_members = list_tar(EQTL_TAR)
    sqtl_members = list_tar(SQTL_TAR)
    eqtl_groups, eqtl_ignored = partition(eqtl_members, (EQTL_EGENE_SUFFIX, EQTL_PAIR_SUFFIX))
    sqtl_groups, sqtl_ignored = partition(sqtl_members, (SQTL_GENE_SUFFIX, SQTL_PAIR_SUFFIX))
    log(
        f"eqtl egenes {len(eqtl_groups[EQTL_EGENE_SUFFIX])} pairs {len(eqtl_groups[EQTL_PAIR_SUFFIX])}; "
        f"sqtl sgenes {len(sqtl_groups[SQTL_GENE_SUFFIX])} pairs {len(sqtl_groups[SQTL_PAIR_SUFFIX])}"
    )

    ratio_member = min(eqtl_groups[EQTL_PAIR_SUFFIX], key=lambda member: member.size)
    pair_ratio = compression_ratio(EQTL_TAR, ratio_member)
    sqtl_small = min(sqtl_groups[SQTL_PAIR_SUFFIX], key=lambda member: member.size)
    pair_ratio = max(pair_ratio, compression_ratio(SQTL_TAR, sqtl_small))
    pair_gzip = sum(member.size for member in eqtl_groups[EQTL_PAIR_SUFFIX])
    pair_gzip += sum(member.size for member in sqtl_groups[SQTL_PAIR_SUFFIX])
    pair_estimate = int(pair_gzip * pair_ratio * 1.5)
    free_now = free_bytes()
    if free_now - pair_estimate < FREE_FLOOR:
        raise SystemExit(
            f"pair estimate {pair_estimate} would leave free space under {FREE_FLOOR}; not writing pairs"
        )
    pair_decision = (
        f"wrote both pair files. Free space {free_now:,} minus estimate {pair_estimate:,} "
        f"stays above {FREE_FLOOR:,}."
    )
    log(pair_decision)

    log("reading GENCODE gene ids")
    gencode_versioned, gencode_bare = load_gencode()
    log(f"GENCODE genes {len(gencode_versioned):,} bare {len(gencode_bare):,}")

    log("median TPM")
    median = load_median()
    log(f"median genes {median['n_genes']:,} tissues {median['n_tissues']}")

    gene_ints = {
        "gene_start", "gene_end", "num_var", "tss_distance", "variant_pos",
        "num_alt_per_site", "minor_allele_samples", "minor_allele_count",
        "ma_samples", "ma_count", "ref_factor", "group_size",
    }
    gene_floats = {
        "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "maf",
        "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval",
        "pval_nominal_threshold", "log2_aFC", "log2_aFC_lower", "log2_aFC_upper",
    }
    pair_ints = {"tss_distance", "ma_samples", "ma_count"}
    pair_floats = {
        "maf", "pval_nominal", "slope", "slope_se", "pval_nominal_threshold",
        "min_pval_nominal", "pval_beta",
    }

    egene = Tally("egene")
    sgene = Tally("sgene")
    eqtl_pairs = Tally("eqtl_pair")
    sqtl_pairs = Tally("sqtl_pair")

    load_tested(
        EQTL_TAR, eqtl_groups[EQTL_EGENE_SUFFIX], EGENE_OUT, EGENE_HEADER,
        gene_ints, gene_floats, egene, truth, phenotype=False,
    )
    load_tested(
        SQTL_TAR, sqtl_groups[SQTL_GENE_SUFFIX], SGENE_OUT, SGENE_HEADER,
        gene_ints, gene_floats, sgene, truth, phenotype=True,
    )
    load_pairs(
        EQTL_TAR, eqtl_groups[EQTL_PAIR_SUFFIX], EQTL_PAIR_OUT, EQTL_PAIR_HEADER,
        pair_ints, pair_floats, eqtl_pairs, phenotype=False,
    )
    load_pairs(
        SQTL_TAR, sqtl_groups[SQTL_PAIR_SUFFIX], SQTL_PAIR_OUT, SQTL_PAIR_HEADER,
        pair_ints, pair_floats, sqtl_pairs, phenotype=True,
    )

    log("confirming truth row on disk")
    truth_result = confirm_truth(truth)
    log(f"truth {truth_result['detail']}")

    joins = []
    for label, ids in (
        ("median TPM", median["genes"]),
        ("eGenes", egene.genes),
        ("sGenes", sgene.genes),
        ("eQTL pairs", eqtl_pairs.genes),
        ("sQTL pairs", sqtl_pairs.genes),
    ):
        joins.append((label, {"ids": ids, "class": classify_genes(ids, gencode_versioned, gencode_bare)}))

    facts = {path.name: file_facts(path) for path in (
        MEDIAN_OUT, EXPR_OUT, EGENE_OUT, SGENE_OUT, EQTL_PAIR_OUT, SQTL_PAIR_OUT,
    )}
    nulls = {
        "egenes": null_map(
            EGENE_OUT,
            ["variant_key", "ensg", "rsid", "slope", "pval_nominal", "qval", "log2_aFC", "log2_aFC_lower", "log2_aFC_upper"],
            egene.rows,
        ),
        "sgenes": null_map(SGENE_OUT, ["variant_key", "ensg", "rsid", "slope", "pval_nominal", "qval"], sgene.rows),
        "eqtl_pairs": null_map(EQTL_PAIR_OUT, ["variant_key", "ensg", "slope", "pval_nominal", "pval_beta"], eqtl_pairs.rows),
        "sqtl_pairs": null_map(SQTL_PAIR_OUT, ["variant_key", "ensg", "gene_id", "slope", "pval_nominal", "pval_beta"], sqtl_pairs.rows),
    }
    for label, tally, key in (
        ("egenes", egene, "variant_key"),
        ("sgenes", sgene, "variant_key"),
        ("eqtl_pairs", eqtl_pairs, "variant_key"),
        ("sqtl_pairs", sqtl_pairs, "variant_key"),
    ):
        recorded = nulls[label][key]
        if recorded is not None and recorded != tally.key_null:
            log(f"WARNING {label} key nulls file={recorded} tally={tally.key_null}")

    b38_total = sum(sum(tally.b38_fail.values()) for tally in (egene, sgene, eqtl_pairs, sqtl_pairs))
    sample_bytes = SAMPLE_TPM.stat().st_size if SAMPLE_TPM.exists() else None
    ctx = {
        "median": median,
        "egene": egene,
        "sgene": sgene,
        "eqtl_pairs": eqtl_pairs,
        "sqtl_pairs": sqtl_pairs,
        "truth": truth,
        "truth_result": truth_result,
        "free_before": free_before,
        "free_after": free_bytes(),
        "pair_estimate": pair_estimate,
        "pair_ratio": pair_ratio,
        "ratio_member": ratio_member.name,
        "pair_decision": pair_decision,
        "eqtl_members": len(eqtl_members),
        "sqtl_members": len(sqtl_members),
        "eqtl_ignored": ", ".join(eqtl_ignored),
        "sqtl_ignored": ", ".join(sqtl_ignored),
        "sample_bytes": sample_bytes,
        "gencode_versioned": len(gencode_versioned),
        "gencode_bare": len(gencode_bare),
        "joins": joins,
        "facts": facts,
        "nulls": nulls,
        "b38_total": b38_total,
        "elapsed": time.time() - started,
    }
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)
    QC_PATH.write_text(render(ctx), encoding="utf-8", newline="\n")
    log(f"wrote {QC_PATH}")
    if not truth_result["pass"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
