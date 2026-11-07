# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Load the MaveDB CC0 dump into the Lane F warehouse tables.

Streams the tar. Does not extract it. Mapped JSON has no scores.
`python load_mavedb.py --scores-only` reads csv/*.scores.csv and joins on
the accession string, which is the mapped variant urn.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import re
import sqlite3
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coords import bare_ensg, canon_chrom, variant_key  # noqa: E402

TAR = ROOT / "raw" / "mavedb" / "mavedb-dump.2026062418131.tar.gz"
CLINVAR = ROOT / "raw" / "clinvar" / "variant_summary.txt.gz"
GENE_DB = ROOT / "db" / "variants.sqlite"
TRUTH = ROOT / "warehouse" / "qc" / "truth_set.tsv"
OUT_DIR = ROOT / "warehouse" / "mavedb"
QC_PATH = ROOT / "warehouse" / "qc" / "mavedb.md"

EXPERIMENTS_PATH = OUT_DIR / "experiments.parquet"
VARIANTS_PATH = OUT_DIR / "mapped_variants.parquet"
COVERAGE_PATH = OUT_DIR / "gene_coverage.parquet"
SCORES_PATH = OUT_DIR / "scores.parquet"
SCORED_VARIANTS_PATH = OUT_DIR / "mapped_variants_scored.parquet"

_PLAIN_NUMBER = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")

BATCH = 200_000

EXPERIMENT_SCHEMA = pa.schema(
    [
        ("experiment_set_urn", pa.string()),
        ("experiment_urn", pa.string()),
        ("experiment_title", pa.string()),
        ("score_set_urn", pa.string()),
        ("score_set_title", pa.string()),
        ("target_name", pa.string()),
        ("target_category", pa.string()),
        ("ensembl_id", pa.string()),
        ("ensg", pa.string()),
        ("organism_name", pa.string()),
        ("assay_text", pa.string()),
        ("assay_text_source", pa.string()),
        ("mapping_state", pa.string()),
        ("num_variants", pa.int64()),
        ("mapped_file_exists", pa.bool_()),
    ]
)

VARIANT_SCHEMA = pa.schema(
    [
        ("urn", pa.string()),
        ("score_set_urn", pa.string()),
        ("target_gene", pa.string()),
        ("protein_change", pa.string()),
        ("sequence_label", pa.string()),
        ("chrom_raw", pa.string()),
        ("chrom", pa.string()),
        ("pos", pa.int64()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("variant_key", pa.string()),
        ("post_mapped_type", pa.string()),
        ("member_index", pa.int32()),
        ("alignment_level", pa.string()),
        ("error_message", pa.string()),
    ]
)

COVERAGE_SCHEMA = pa.schema(
    [
        ("target_gene", pa.string()),
        ("n_score_sets", pa.int64()),
        ("n_mapped_variants", pa.int64()),
        ("has_genomic_mapping", pa.bool_()),
    ]
)


def is_appledouble(name: str) -> bool:
    base = name.rsplit("/", 1)[-1]
    return "/._" in name or base.startswith("._")


def as_int(value):
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value


def as_text(value):
    if value is None:
        return None
    if isinstance(value, str):
        return value
    return str(value)


def join_unique(values: list[str]) -> str | None:
    seen = set()
    out = []
    for value in values:
        if value is None or value in seen:
            continue
        seen.add(value)
        out.append(value)
    if not out:
        return None
    return out[0] if len(out) == 1 else "|".join(out)


def load_accession_map(path: Path):
    """GRCh38 ChromosomeAccession -> Chromosome, taken only from ClinVar.

    One accession is listed as both X and Y. The chromosome with more GRCh38
    rows is kept. A tie is left unmapped. Nothing is hard-coded.
    """
    counts = defaultdict(Counter)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        i_asm = header.index("Assembly")
        i_acc = header.index("ChromosomeAccession")
        i_chr = header.index("Chromosome")
        width = len(header)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < width or parts[i_asm] != "GRCh38":
                continue
            accession = parts[i_acc].strip()
            chrom = parts[i_chr].strip()
            if not accession or accession in {"-", ".", "na", "NA"}:
                continue
            counts[accession][chrom] += 1

    chosen = {}
    disagreements = []
    ties = []
    for accession, counter in counts.items():
        best_chrom, best_n = counter.most_common(1)[0]
        leaders = [chrom for chrom, n in counter.items() if n == best_n]
        if len(leaders) > 1:
            ties.append((accession, dict(counter)))
            continue
        if len(counter) > 1:
            disagreements.append((accession, dict(counter), best_chrom))
        chosen[accession] = best_chrom
    return chosen, disagreements, ties


def load_gene_symbols(path: Path):
    uri = "file:" + path.as_posix() + "?mode=ro"
    con = sqlite3.connect(uri, uri=True)
    try:
        rows = con.execute("select gene_id, symbol, biotype from gene").fetchall()
    finally:
        con.close()
    protein = defaultdict(set)
    any_biotype = defaultdict(set)
    biotype_counts = Counter()
    for gene_id, symbol, biotype in rows:
        biotype_counts[biotype] += 1
        if symbol:
            any_biotype[symbol].add((gene_id, biotype))
            if biotype == "protein_coding":
                protein[symbol].add(gene_id)
    return protein, any_biotype, biotype_counts


def load_truth(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {row["case_id"]: row for row in rows}
    if "mavedb_mapped" not in by_id or "mavedb_unmapped" not in by_id:
        raise SystemExit("truth_set.tsv is missing a MaveDB case")
    return by_id["mavedb_mapped"], by_id["mavedb_unmapped"]


def urn_from_mapped_filename(base: str) -> str:
    suffix = ".mapped-variants.json"
    if not base.endswith(suffix):
        raise ValueError(base)
    stem = base[: -len(suffix)]
    prefix = "urn-mavedb-"
    if not stem.startswith(prefix):
        raise ValueError(base)
    return "urn:mavedb:" + stem[len(prefix) :]


def ensembl_from_target(target: dict) -> list[str]:
    found = []
    for item in target.get("externalIdentifiers") or []:
        if not isinstance(item, dict):
            continue
        ident = item.get("identifier")
        if not isinstance(ident, dict):
            continue
        if ident.get("dbName") == "Ensembl" and ident.get("identifier"):
            found.append(str(ident["identifier"]))
    return found


def organism_of(target: dict) -> str | None:
    sequence = target.get("targetSequence") or {}
    if not isinstance(sequence, dict):
        return None
    taxonomy = sequence.get("taxonomy") or {}
    if not isinstance(taxonomy, dict):
        return None
    name = taxonomy.get("organismName")
    return name if isinstance(name, str) and name else None


def assay_from(score_set: dict, experiment: dict):
    """MaveDB has no assay-type field. methodText is the assay writeup."""
    score_method = score_set.get("methodText")
    if isinstance(score_method, str) and score_method.strip():
        return score_method, "scoreSets.methodText"
    experiment_method = experiment.get("methodText")
    if isinstance(experiment_method, str) and experiment_method.strip():
        return experiment_method, "experiments.methodText"
    return None, None


def local_hgvs_p(obj) -> list[str]:
    if not isinstance(obj, dict):
        return []
    values = []
    for expr in obj.get("expressions") or []:
        if not isinstance(expr, dict):
            continue
        if expr.get("syntax") == "hgvs.p" and isinstance(expr.get("value"), str):
            values.append(expr["value"])
    return values


def iter_alleles(obj):
    if not isinstance(obj, dict):
        return
    kind = obj.get("type")
    if kind == "Allele":
        yield obj
        return
    if kind in {"Haplotype", "CisPhasedBlock"}:
        for member in obj.get("members") or []:
            yield from iter_alleles(member)


def allele_fields(allele: dict):
    refs = []
    n_ext = 0
    for ext in allele.get("extensions") or []:
        if not isinstance(ext, dict) or ext.get("name") != "vrs_ref_allele_seq":
            continue
        n_ext += 1
        if isinstance(ext.get("value"), str):
            refs.append(ext["value"])
    distinct = set(refs)
    ref = refs[0] if len(distinct) == 1 else None
    ambiguous_ref = len(distinct) > 1

    state = allele.get("state")
    alt = None
    if isinstance(state, str):
        alt = state
    elif isinstance(state, dict) and isinstance(state.get("sequence"), str):
        alt = state["sequence"]

    loc = allele.get("location") if isinstance(allele.get("location"), dict) else {}
    start = as_int(loc.get("start"))
    end = as_int(loc.get("end"))
    reference = loc.get("sequenceReference")
    label = None
    if isinstance(reference, dict) and isinstance(reference.get("label"), str):
        text = reference["label"].strip()
        label = text or None
    return label, start, end, ref, alt, ambiguous_ref, n_ext


def protein_change_for(post_allele, pre_obj, member_index: int, post_is_allele: bool):
    values = local_hgvs_p(post_allele)
    if values:
        return join_unique(values)
    if not isinstance(pre_obj, dict):
        return None
    if post_is_allele and pre_obj.get("type") == "Allele":
        return join_unique(local_hgvs_p(pre_obj))
    if not post_is_allele:
        pre_alleles = list(iter_alleles(pre_obj))
        if member_index < len(pre_alleles):
            return join_unique(local_hgvs_p(pre_alleles[member_index]))
    return None


class Writer:
    def __init__(self, path: Path, schema: pa.Schema):
        self.path = path
        self.schema = schema
        self.writer = pq.ParquetWriter(path, schema, compression="snappy")
        self.cols = {name: [] for name in schema.names}
        self.pending = 0
        self.rows = 0

    def add(self, row: dict):
        for name in self.schema.names:
            self.cols[name].append(row.get(name))
        self.pending += 1
        self.rows += 1
        if self.pending >= BATCH:
            self.flush()

    def flush(self):
        if self.pending == 0:
            return
        self.writer.write_table(pa.table(self.cols, schema=self.schema))
        for name in self.cols:
            self.cols[name].clear()
        self.pending = 0

    def close(self):
        self.flush()
        self.writer.close()


def target_summary(targets: list[dict]):
    names = []
    categories = []
    organisms = []
    ensembl = []
    bare = []
    for target in targets:
        name = target.get("name")
        if isinstance(name, str) and name:
            names.append(name)
        category = target.get("category")
        if isinstance(category, str) and category:
            categories.append(category)
        organism = organism_of(target)
        if organism:
            organisms.append(organism)
        for ident in ensembl_from_target(target):
            ensembl.append(ident)
            bare_id = bare_ensg(ident)
            if bare_id:
                bare.append(bare_id)
    distinct_names = list(dict.fromkeys(names))
    return {
        "names": names,
        "distinct_names": distinct_names,
        "target_name": join_unique(distinct_names),
        "target_category": join_unique(categories),
        "organism_name": join_unique(organisms),
        "ensembl_id": join_unique(ensembl),
        "ensg": join_unique(bare),
    }


def build(acc_map: dict[str, str], stats: dict):
    experiments = []
    coverage = {}
    mapped_urns = set()
    score_header = None
    score_header_name = None
    main = None
    pending_mapped = []

    def ensure_gene(name: str):
        row = coverage.get(name)
        if row is None:
            row = {
                "target_gene": name,
                "n_score_sets": 0,
                "n_mapped_variants": 0,
                "has_genomic_mapping": False,
                "has_variant_key": False,
                "organisms": set(),
                "categories": set(),
            }
            coverage[name] = row
        return row

    def consume_main(payload: dict):
        nonlocal main
        main = payload
        for experiment_set in payload.get("experimentSets") or []:
            set_urn = experiment_set.get("urn")
            for experiment in experiment_set.get("experiments") or []:
                for score_set in experiment.get("scoreSets") or []:
                    summary = target_summary(score_set.get("targetGenes") or [])
                    for name in summary["distinct_names"]:
                        gene = ensure_gene(name)
                        gene["n_score_sets"] += 1
                    for target in score_set.get("targetGenes") or []:
                        name = target.get("name")
                        if not isinstance(name, str) or not name:
                            stats["unnamed_targets"] += 1
                            continue
                        gene = ensure_gene(name)
                        organism = organism_of(target)
                        if organism:
                            gene["organisms"].add(organism)
                        category = target.get("category")
                        if isinstance(category, str) and category:
                            gene["categories"].add(category)
                    text, source = assay_from(score_set, experiment)
                    stats["assay_source"][source or "missing"] += 1
                    num = score_set.get("numVariants")
                    experiments.append(
                        {
                            "experiment_set_urn": as_text(set_urn),
                            "experiment_urn": as_text(experiment.get("urn")),
                            "experiment_title": as_text(experiment.get("title")),
                            "score_set_urn": as_text(score_set.get("urn")),
                            "score_set_title": as_text(score_set.get("title")),
                            "target_name": summary["target_name"],
                            "target_category": summary["target_category"],
                            "ensembl_id": summary["ensembl_id"],
                            "ensg": summary["ensg"],
                            "organism_name": summary["organism_name"],
                            "assay_text": text,
                            "assay_text_source": source,
                            "mapping_state": as_text(score_set.get("mappingState")),
                            "num_variants": num if isinstance(num, int) and not isinstance(num, bool) else None,
                            "mapped_file_exists": False,
                            "_names": summary["distinct_names"],
                        }
                    )
                    columns = ((score_set.get("datasetColumns") or {}).get("scoreColumns") or [])
                    if isinstance(columns, list) and "score" in columns:
                        stats["score_column_listed"] += 1
                    stats["score_sets"] += 1

    def consume_mapped(base: str, payload):
        score_urn = urn_from_mapped_filename(base)
        mapped_urns.add(score_urn)
        if not isinstance(payload, list):
            stats["mapped_not_list"] += 1
            return
        meta = stats["by_urn"].get(score_urn)
        names = meta["_names"] if meta else []
        single = len(names) == 1
        target_gene = meta["target_name"] if meta else None
        if meta is None:
            stats["mapped_without_metadata"] += 1
        seen = {}
        kept_records = 0
        genomic_records = 0
        has_key_records = 0
        for record in payload:
            stats["raw_records"] += 1
            if not isinstance(record, dict):
                stats["bad_record"] += 1
                continue
            if record.get("current") is not True:
                stats["non_current"] += 1
                continue
            stats["current_records"] += 1
            variant_urn = record.get("variantUrn")
            if isinstance(variant_urn, str) and "#" in variant_urn:
                prefix = variant_urn.split("#", 1)[0]
                if prefix != score_urn:
                    stats["urn_filename_mismatch"] += 1
                    if len(stats["urn_mismatch_examples"]) < 8:
                        stats["urn_mismatch_examples"].append((score_urn, variant_urn))
            else:
                stats["urn_filename_mismatch"] += 1
                if len(stats["urn_mismatch_examples"]) < 8:
                    stats["urn_mismatch_examples"].append((score_urn, variant_urn))
            post = record.get("postMapped")
            allele_sig = []
            if isinstance(post, dict):
                for allele in iter_alleles(post):
                    label, start, end, ref, alt, _ambiguous, _n_ext = allele_fields(allele)
                    allele_sig.append((label, start, end, ref, alt))
            signature = (
                post.get("type") if isinstance(post, dict) else None,
                tuple(allele_sig),
                record.get("alignmentLevel"),
            )
            previous = seen.get(variant_urn)
            if previous == signature:
                stats["exact_duplicate_records"] += 1
                continue
            if previous is not None:
                stats["conflicting_duplicate_records"] += 1
            seen[variant_urn] = signature
            kept_records += 1
            rows, has_key, on_grch38, labels = rows_for_record(
                record, score_urn, target_gene, acc_map, stats
            )
            if post is None:
                stats["null_post_records"] += 1
            elif on_grch38:
                stats["grch38_records"] += 1
                genomic_records += 1
                if has_key:
                    stats["grch38_with_key"] += 1
                    has_key_records += 1
                else:
                    stats["grch38_without_key"] += 1
            else:
                stats["non_grch38_records"] += 1
                if labels and all(label.startswith("NP_") for label in labels):
                    stats["non_grch38_protein"] += 1
                elif not labels:
                    stats["non_grch38_unlabeled"] += 1
                else:
                    stats["non_grch38_other"] += 1
            for row in rows:
                stats["writer"].add(row)
                if row["variant_key"] is None:
                    stats["null_variant_key_rows"] += 1
                if row["protein_change"] is None:
                    stats["null_protein_rows"] += 1
                if row["chrom"] is None:
                    stats["null_chrom_rows"] += 1
                if row["pos"] is None:
                    stats["null_pos_rows"] += 1
                if row["ref"] is None:
                    stats["null_ref_rows"] += 1
                if row["alt"] is None:
                    stats["null_alt_rows"] += 1
                if row["sequence_label"] is None:
                    stats["null_label_rows"] += 1
                if row["target_gene"] is None:
                    stats["null_target_rows"] += 1
        if single:
            gene = ensure_gene(names[0])
            gene["n_mapped_variants"] += kept_records
            if genomic_records:
                gene["has_genomic_mapping"] = True
            if has_key_records:
                gene["has_variant_key"] = True
        elif names:
            stats["multi_target_mapped_records"] += kept_records
            stats["multi_target_mapped_sets"] += 1
        stats["kept_records"] += kept_records

    stats["by_urn"] = {}
    with tarfile.open(TAR, "r:gz") as archive:
        for member in archive:
            name = member.name
            if is_appledouble(name):
                stats["appledouble"] += 1
                continue
            base = name.rsplit("/", 1)[-1]
            if base == "main.json" and name.endswith("/main.json"):
                consume_main(json.loads(archive.extractfile(member).read()))
                stats["by_urn"] = {row["score_set_urn"]: row for row in experiments}
                continue
            if base.endswith(".scores.csv") and score_header is None:
                line = archive.extractfile(member).readline()
                score_header = line.decode("utf-8", errors="replace").rstrip("\r\n")
                score_header_name = name
                continue
            if not (name.endswith(".mapped-variants.json") and "/mapped/" in name):
                continue
            stats["mapped_files"] += 1
            if stats["mapped_files"] % 100 == 0:
                print(f"mapped files {stats['mapped_files']}", file=sys.stderr, flush=True)
            payload = json.loads(archive.extractfile(member).read())
            if main is None:
                pending_mapped.append((base, payload))
            else:
                consume_mapped(base, payload)
        for base, payload in pending_mapped:
            consume_mapped(base, payload)

    if main is None:
        raise SystemExit("main.json was not in the tar")
    for row in experiments:
        row["mapped_file_exists"] = row["score_set_urn"] in mapped_urns
        state = row["mapping_state"]
        flag = "file" if row["mapped_file_exists"] else "no_file"
        stats["state_vs_file"][f"{state}|{flag}"] += 1
    stats["score_header"] = score_header
    stats["score_header_name"] = score_header_name
    stats["dump_title"] = main.get("title")
    stats["dump_as_of"] = main.get("asOf")
    stats["experiment_sets"] = len(main.get("experimentSets") or [])
    return experiments, coverage


def rows_for_record(record, score_urn, target_gene, acc_map, stats):
    post = record.get("postMapped")
    pre = record.get("preMapped")
    post_type = post.get("type") if isinstance(post, dict) else None
    if post is not None and not isinstance(post, dict):
        stats["post_not_object"] += 1
    alleles = list(iter_alleles(post)) if isinstance(post, dict) else []
    if isinstance(post, dict) and post_type not in {"Allele", "Haplotype", "CisPhasedBlock", None}:
        stats["other_post_types"][post_type] += 1
    if not alleles:
        stats["allele_rows"] += 1
        return (
            [
                {
                    "urn": as_text(record.get("variantUrn")),
                    "score_set_urn": score_urn,
                    "target_gene": target_gene,
                    "protein_change": None,
                    "sequence_label": None,
                    "chrom_raw": None,
                    "chrom": None,
                    "pos": None,
                    "ref": None,
                    "alt": None,
                    "variant_key": None,
                    "post_mapped_type": as_text(post_type),
                    "member_index": None,
                    "alignment_level": as_text(record.get("alignmentLevel")),
                    "error_message": as_text(record.get("errorMessage")),
                }
            ],
            False,
            False,
            [],
        )

    rows = []
    labels = []
    resolved = False
    on_grch38 = False
    post_is_allele = post_type == "Allele"
    for index, allele in enumerate(alleles):
        label, start, end, ref, alt, ambiguous_ref, n_ext = allele_fields(allele)
        if n_ext > 1:
            stats["multi_ref_extension"] += 1
        if ambiguous_ref:
            stats["ambiguous_ref"] += 1
        if label:
            labels.append(label)
            if label in acc_map:
                stats["label_in_grch38_map"] += 1
            elif label.startswith("NP_"):
                stats["label_np"] += 1
            else:
                stats["label_other"] += 1
                stats["other_labels"][label] += 1
        else:
            stats["label_missing"] += 1

        chrom_raw = None
        chrom = None
        pos = None
        key = None
        # VRS locations are 0-based and half-open. A 1-based VCF position is
        # start+1 only when the span equals the reference length, and only on
        # a GRCh38 chromosome accession from ClinVar. Protein labels stay
        # unplaced. Empty alleles do not become a variant key.
        if label in acc_map:
            chrom_raw = acc_map[label]
            chrom = canon_chrom(chrom_raw)
            if chrom is None:
                stats["accession_not_canonical"] += 1
            elif ref is not None and start is not None and end is not None:
                if (end - start) == len(ref):
                    pos = start + 1
                    key = variant_key(chrom, pos, ref, alt)
                    if key is None:
                        stats["grch38_no_key"] += 1
                        if ref == "":
                            stats["empty_ref_on_grch38"] += 1
                        if alt == "":
                            stats["empty_alt_on_grch38"] += 1
                else:
                    stats["span_mismatch_grch38"] += 1
            if chrom is not None:
                on_grch38 = True
            if key is not None:
                resolved = True
        if record.get("variantUrn") == stats["truth_variant_urn"]:
            stats["truth_hits"].append(
                {
                    "key": key,
                    "label": label,
                    "chrom_raw": chrom_raw,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "start": start,
                    "end": end,
                    "target": target_gene,
                }
            )
        rows.append(
            {
                "urn": as_text(record.get("variantUrn")),
                "score_set_urn": score_urn,
                "target_gene": target_gene,
                "protein_change": protein_change_for(allele, pre, index, post_is_allele),
                "sequence_label": label,
                "chrom_raw": chrom_raw,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "variant_key": key,
                "post_mapped_type": as_text(post_type),
                "member_index": index,
                "alignment_level": as_text(record.get("alignmentLevel")),
                "error_message": as_text(record.get("errorMessage")),
            }
        )
        stats["allele_rows"] += 1
    return rows, resolved, on_grch38, labels


def write_experiments(rows: list[dict]):
    columns = {name: [] for name in EXPERIMENT_SCHEMA.names}
    for row in rows:
        for name in EXPERIMENT_SCHEMA.names:
            columns[name].append(row.get(name))
    table = pa.table(columns, schema=EXPERIMENT_SCHEMA)
    pq.write_table(table, EXPERIMENTS_PATH, compression="snappy")
    return table.num_rows


def write_coverage(coverage: dict):
    names = sorted(coverage)
    columns = {
        "target_gene": names,
        "n_score_sets": [coverage[name]["n_score_sets"] for name in names],
        "n_mapped_variants": [coverage[name]["n_mapped_variants"] for name in names],
        "has_genomic_mapping": [coverage[name]["has_genomic_mapping"] for name in names],
    }
    table = pa.table(columns, schema=COVERAGE_SCHEMA)
    pq.write_table(table, COVERAGE_PATH, compression="snappy")
    return table.num_rows


def pct(part, whole) -> str:
    if not whole:
        return "0/0"
    return f"{part}/{whole} ({part / whole:.4%})"


def markdown(stats, experiments, coverage, protein, any_biotype, biotype_counts, disagreements, ties, truth_mapped, truth_unmapped, checks):
    score_sets = [row for row in experiments]
    with_file = [row for row in score_sets if row["mapped_file_exists"]]
    without_file = [row for row in score_sets if not row["mapped_file_exists"]]
    experiments_with = set()
    experiments_all = set()
    for row in score_sets:
        experiments_all.add(row["experiment_urn"])
        if row["mapped_file_exists"]:
            experiments_with.add(row["experiment_urn"])

    protein_names = []
    category_protein = []
    human_names = []
    human_protein_names = []
    phrase_unmatched = []
    other_biotype = []
    unmatched = []
    matched_gene_ids = set()
    for name, gene in coverage.items():
        ids = protein.get(name)
        if ids:
            protein_names.append(name)
            matched_gene_ids.update(ids)
            if "Homo sapiens" in gene["organisms"]:
                human_protein_names.append(name)
        else:
            unmatched.append(name)
            if any(ch in name for ch in " ()/"):
                phrase_unmatched.append(name)
            hits = any_biotype.get(name)
            if hits:
                other_biotype.append((name, sorted({item[1] for item in hits})))
        if "protein_coding" in gene["categories"]:
            category_protein.append(name)
        if "Homo sapiens" in gene["organisms"]:
            human_names.append(name)

    n_rows = stats["writer"].rows
    lines = [
        "# MaveDB QC",
        "",
        f"Dump title: {stats['dump_title']}. asOf: {stats['dump_as_of']}.",
        f"Archive: `{TAR}`.",
        "The tar was streamed. It was not extracted.",
        "",
        "## AppleDouble",
        "",
        f"Skipped {stats['appledouble']} members whose names contain `/._` or whose filename starts with `._`. Their bytes were not read.",
        "",
        "## Experiments and score sets",
        "",
        f"Experiment sets: {stats['experiment_sets']}. Experiments: {len(experiments_all)}. Score sets: {len(score_sets)}.",
        f"Mapped score sets (a `mapped/*.mapped-variants.json` file is in the tar): {len(with_file)}.",
        f"Unmapped score sets (no mapped file): {len(without_file)}.",
        f"Experiments with at least one mapped score set: {len(experiments_with)}.",
        f"Experiments with no mapped score set: {len(experiments_all) - len(experiments_with)}.",
        f"Mapped JSON files read: {stats['mapped_files']}.",
        "",
        "MaveDB has no assay-type field. `assay_text` is `methodText`.",
        "",
        "Assay text source:",
        "",
    ]
    for key, count in sorted(stats["assay_source"].items(), key=lambda item: (-item[1], str(item[0]))):
        lines.append(f"- {key}: {count}")
    lines.extend(["", "mappingState versus mapped file:", ""])
    for key, count in sorted(stats["state_vs_file"].items()):
        lines.append(f"- {key}: {count}")
    lines.extend(
        [
            "",
            f"Score sets whose `datasetColumns.scoreColumns` lists `score`: {stats['score_column_listed']}/{stats['score_sets']}.",
            "That name is metadata on the score set. It is not a value in the mapped JSON.",
            "",
            "## Scores",
            "",
            "Mapped records do not contain a score. A mapped record's keys are variant identity, preMapped, postMapped, mapping status, and ClinGen id. `mapped_variants.parquet` therefore has no score column and no author-normalized score. Inventing either from a column name would be wrong.",
            "",
            "Scores live in `csv/{urn}.scores.csv`. The first scores file reached in the tar was",
            f"`{stats['score_header_name']}` with header:",
            "",
            "```",
            stats["score_header"] or "",
            "```",
            "",
            "The primary CSV column is `scores.score`. Extra `scores.*` columns differ by score set and were not loaded. README says score values are not normalized across score sets.",
            "",
            "## Mapped records",
            "",
            f"Raw records in mapped JSON files: {stats['raw_records']}.",
            f"`current` is not true, excluded: {stats['non_current']}.",
            f"Exact duplicate records of the same variant urn in one file, excluded: {stats['exact_duplicate_records']}.",
            f"Same variant urn with a different postMapped payload, kept: {stats['conflicting_duplicate_records']}.",
            f"Records kept: {stats['kept_records']}.",
            f"Allele rows written: {n_rows}. A haplotype is one record and one row per member allele.",
            f"Records with postMapped null: {stats['null_post_records']}. These stay in the table with no coordinate.",
            f"Records placed on a GRCh38 chromosome: {stats['grch38_records']}. Of those, {stats['grch38_with_key']} have a variant_key and {stats['grch38_without_key']} have a chromosome but no key.",
            f"Mapped records that are non-GRCh38 (postMapped present, and no allele label is a GRCh38 chromosome accession): {stats['non_grch38_records']}.",
            f"- protein labels (`NP_`): {stats['non_grch38_protein']}",
            f"- no sequence label: {stats['non_grch38_unlabeled']}",
            f"- other labels: {stats['non_grch38_other']}",
            "",
            "Allele labels:",
            "",
            f"- on a GRCh38 ChromosomeAccession: {stats['label_in_grch38_map']}",
            f"- `NP_` protein, variant_key left null: {stats['label_np']}",
            f"- missing label: {stats['label_missing']}",
            f"- any other label, variant_key left null: {stats['label_other']}",
            "",
            "A chromosome is built only when the label is a GRCh38 `ChromosomeAccession` in ClinVar `variant_summary.txt.gz` (`Assembly == GRCh38`). The VCF position is `start + 1` when `end - start` equals the `vrs_ref_allele_seq` length. `variant_key` comes from `coords.variant_key`, so an empty reference or alternate stays null. No anchor base was added.",
            "",
            f"GRCh38 alleles whose span matched but `variant_key` was still null: {stats['grch38_no_key']} (empty ref {stats['empty_ref_on_grch38']}, empty alt {stats['empty_alt_on_grch38']}).",
            f"GRCh38 alleles whose span did not equal the reference length, position left null: {stats['span_mismatch_grch38']}.",
            f"Accessions whose ClinVar chromosome did not canonicalize: {stats['accession_not_canonical']}.",
            f"Ambiguous `vrs_ref_allele_seq` values: {stats['ambiguous_ref']}. Extra copies of that extension: {stats['multi_ref_extension']}.",
            f"Variant urns that do not start with the score-set urn from the filename: {stats['urn_filename_mismatch']}. These are ids stored on the mapped record. The examples are temporary ids, not a filename parse failure.",
            "Examples:",
            "",
        ]
    )
    if stats["urn_mismatch_examples"]:
        for file_urn, variant_urn in stats["urn_mismatch_examples"]:
            lines.append(f"- `{file_urn}` vs `{variant_urn}`")
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            f"Mapped files with no main.json score set: {stats['mapped_without_metadata']}.",
            f"Mapped JSON members that were not arrays: {stats['mapped_not_list']}.",
            f"Mapped records that were not objects: {stats['bad_record']}.",
            "",
        ]
    )
    if stats["other_labels"]:
        lines.append("Other labels, up to 15:")
        lines.append("")
        for label, count in stats["other_labels"].most_common(15):
            lines.append(f"- {label}: {count}")
        lines.append("")
    if disagreements:
        lines.append("ClinVar GRCh38 accessions with more than one chromosome value. The majority value was used:")
        lines.append("")
        for accession, counter, chosen in disagreements:
            lines.append(f"- {accession}: {counter}; used {chosen}")
        lines.append("")
    if ties:
        lines.append("ClinVar GRCh38 accessions with a tied chromosome count were not used:")
        lines.append("")
        for accession, counter in ties:
            lines.append(f"- {accession}: {counter}")
        lines.append("")
    lines.extend(
        [
            "Multi-target score sets are not copied onto every gene. Their mapped records stay on the score-set row with the joined target names.",
            f"Multi-target score sets that had a mapped file: {stats['multi_target_mapped_sets']}. Records in those files: {stats['multi_target_mapped_records']}.",
            "",
            "## Null rates",
            "",
            f"mapped_variants rows: {n_rows}.",
            f"- variant_key null: {pct(stats['null_variant_key_rows'], n_rows)}",
            f"- chrom null: {pct(stats['null_chrom_rows'], n_rows)}",
            f"- pos null: {pct(stats['null_pos_rows'], n_rows)}",
            f"- ref null: {pct(stats['null_ref_rows'], n_rows)}",
            f"- alt null: {pct(stats['null_alt_rows'], n_rows)}",
            f"- protein_change null: {pct(stats['null_protein_rows'], n_rows)}",
            f"- sequence_label null: {pct(stats['null_label_rows'], n_rows)}",
            f"- target_gene null: {pct(stats['null_target_rows'], n_rows)}",
            "",
        ]
    )
    n_exp = len(score_sets)
    null_ensembl = sum(1 for row in score_sets if not row["ensembl_id"])
    null_assay = sum(1 for row in score_sets if not row["assay_text"])
    null_target = sum(1 for row in score_sets if not row["target_name"])
    null_state = sum(1 for row in score_sets if not row["mapping_state"])
    lines.extend(
        [
            f"experiments rows: {n_exp}.",
            f"- ensembl_id null: {pct(null_ensembl, n_exp)}",
            f"- assay_text null: {pct(null_assay, n_exp)}",
            f"- target_name null: {pct(null_target, n_exp)}",
            f"- mapping_state null: {pct(null_state, n_exp)}",
            "",
            "## Gene coverage",
            "",
            f"Distinct raw target names: {len(coverage)}.",
            f"Names MaveDB categories as `protein_coding`: {len(category_protein)}. This category includes phrases and non-human targets. It is not a GENCODE biotype.",
            f"Names that are Homo sapiens: {len(human_names)}.",
            f"Target names that exactly match a GENCODE protein-coding symbol: {len(protein_names)}.",
            f"Of those, Homo sapiens: {len(human_protein_names)}.",
            f"GENCODE protein-coding genes in `db/variants.sqlite` (`biotype == protein_coding`): {biotype_counts['protein_coding']}.",
            f"Protein-coding gene ids touched by an exact symbol match: {len(matched_gene_ids)}.",
            f"Protein-coding gene ids with no MaveDB score set under that rule: {biotype_counts['protein_coding'] - len(matched_gene_ids)}.",
            f"Target names with no exact protein-coding symbol match: {len(unmatched)}.",
            f"- names that contain a space or punctuation and do not match: {len(phrase_unmatched)}",
            f"- names that match a GENCODE symbol of another biotype: {len(other_biotype)}",
            "",
            "The unmatched human protein-coding gene list is the rest of the catalog. It is not written out. MaveDB covers a small slice of proteins.",
            "",
            "Sample of target names that are not protein-coding symbols:",
            "",
        ]
    )
    for name in unmatched[:20]:
        lines.append(f"- {name}")
    if other_biotype:
        lines.extend(["", "Sample of names that match a non-protein-coding GENCODE symbol:", ""])
        for name, biotypes in other_biotype[:10]:
            lines.append(f"- {name}: {', '.join(biotypes)}")
    genomic_genes = sum(1 for gene in coverage.values() if gene["has_genomic_mapping"])
    keyed_genes = sum(1 for gene in coverage.values() if gene["has_variant_key"])
    mapped_genes = sum(1 for gene in coverage.values() if gene["n_mapped_variants"])
    lines.extend(
        [
            "",
            f"Target genes with at least one mapped variant record: {mapped_genes}.",
            f"Target genes with any GRCh38 genomic mapping: {genomic_genes}. This flag is true when a chromosome was built, including alleles that have no variant_key.",
            f"Target genes with any variant_key: {keyed_genes}.",
            f"Sum of per-gene mapped-variant counts: {sum(gene['n_mapped_variants'] for gene in coverage.values())}. This counts kept variant records for score sets with one target name, not allele rows.",
            "",
            "## Truth set",
            "",
            f"Mapped case `{truth_mapped['case_id']}` expects variant `{stats['truth_variant_urn']}` to produce `{truth_mapped['variant_key']}`.",
            f"Hits written for that urn: {len(stats['truth_hits'])}.",
        ]
    )
    for hit in stats["truth_hits"]:
        lines.append(
            f"- label {hit['label']}, chrom_raw {hit['chrom_raw']}, chrom {hit['chrom']}, pos {hit['pos']}, ref {hit['ref']}, alt {hit['alt']}, key {hit['key']}, target {hit['target']}"
        )
    lines.extend(
        [
            "",
            f"Unmapped case `{truth_unmapped['case_id']}` expects score set `{stats['truth_unmapped_urn']}` gene `{truth_unmapped['gene']}` to stay without a coordinate.",
            f"On-disk checks: {checks['text']}",
            "",
            "## Result",
            "",
            checks["result"],
            "",
        ]
    )
    return "\n".join(lines)


def disk_checks(truth_mapped, truth_unmapped, stats):
    import pyarrow.dataset as ds

    variant_urn = stats["truth_variant_urn"]
    unmapped_urn = stats["truth_unmapped_urn"]
    variants = ds.dataset(str(VARIANTS_PATH))
    hit = variants.to_table(
        filter=ds.field("urn") == variant_urn,
        columns=["urn", "target_gene", "sequence_label", "chrom", "pos", "ref", "alt", "variant_key"],
    )
    keys = hit.column("variant_key").to_pylist()
    labels = hit.column("sequence_label").to_pylist()
    targets = hit.column("target_gene").to_pylist()
    expected = truth_mapped["variant_key"]
    mapped_ok = bool(keys) and all(key == expected for key in keys) and any(
        label == truth_mapped["chrom_raw"] for label in labels
    )
    gene_ok = bool(targets) and all(target == truth_mapped["gene"] for target in targets)

    unmapped_rows = variants.to_table(
        filter=ds.field("score_set_urn") == unmapped_urn,
        columns=["score_set_urn"],
    )
    experiments = ds.dataset(str(EXPERIMENTS_PATH)).to_table(
        filter=ds.field("score_set_urn") == unmapped_urn,
        columns=["score_set_urn", "target_name", "mapped_file_exists"],
    )
    exp_urns = experiments.column("score_set_urn").to_pylist()
    exp_genes = experiments.column("target_name").to_pylist()
    exp_files = experiments.column("mapped_file_exists").to_pylist()
    coverage_table = pq.read_table(COVERAGE_PATH)
    gene_rows = [
        {
            "target_gene": gene,
            "n_mapped_variants": n,
            "has_genomic_mapping": flag,
        }
        for gene, n, flag in zip(
            coverage_table.column("target_gene").to_pylist(),
            coverage_table.column("n_mapped_variants").to_pylist(),
            coverage_table.column("has_genomic_mapping").to_pylist(),
        )
        if gene == truth_unmapped["gene"]
    ]
    unmapped_ok = (
        unmapped_rows.num_rows == 0
        and exp_urns == [unmapped_urn]
        and exp_genes == [truth_unmapped["gene"]]
        and exp_files == [False]
        and len(gene_rows) == 1
        and gene_rows[0]["n_mapped_variants"] == 0
        and gene_rows[0]["has_genomic_mapping"] is False
    )
    variant_rows = pq.ParquetFile(VARIANTS_PATH).metadata.num_rows
    experiment_rows = pq.ParquetFile(EXPERIMENTS_PATH).metadata.num_rows
    coverage_rows = pq.ParquetFile(COVERAGE_PATH).metadata.num_rows
    counts_ok = variant_rows == stats["writer"].rows and experiment_rows == stats["score_sets"] and coverage_rows > 0
    passed = mapped_ok and gene_ok and unmapped_ok and counts_ok
    text = (
        f"mapped key ok {mapped_ok}; target name ok {gene_ok}; "
        f"unmapped score set absent from mapped_variants ({unmapped_rows.num_rows} rows) "
        f"and experiments say no file and gene {truth_unmapped['gene']} coverage {gene_rows}; "
        f"parquet rows variants {variant_rows}, experiments {experiment_rows}, genes {coverage_rows}."
    )
    result = "PASS" if passed else "FAIL"
    return {"passed": passed, "text": text, "result": result, "gene_rows": gene_rows}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)
    truth_mapped, truth_unmapped = load_truth(TRUTH)
    print("loading ClinVar accessions", file=sys.stderr, flush=True)
    acc_map, disagreements, ties = load_accession_map(CLINVAR)
    protein, any_biotype, biotype_counts = load_gene_symbols(GENE_DB)
    stats = {
        "appledouble": 0,
        "mapped_files": 0,
        "score_sets": 0,
        "score_column_listed": 0,
        "raw_records": 0,
        "current_records": 0,
        "non_current": 0,
        "exact_duplicate_records": 0,
        "conflicting_duplicate_records": 0,
        "kept_records": 0,
        "allele_rows": 0,
        "null_post_records": 0,
        "grch38_records": 0,
        "grch38_with_key": 0,
        "grch38_without_key": 0,
        "non_grch38_records": 0,
        "non_grch38_protein": 0,
        "non_grch38_unlabeled": 0,
        "non_grch38_other": 0,
        "label_in_grch38_map": 0,
        "label_np": 0,
        "label_missing": 0,
        "label_other": 0,
        "other_labels": Counter(),
        "grch38_no_key": 0,
        "empty_ref_on_grch38": 0,
        "empty_alt_on_grch38": 0,
        "span_mismatch_grch38": 0,
        "accession_not_canonical": 0,
        "ambiguous_ref": 0,
        "multi_ref_extension": 0,
        "urn_filename_mismatch": 0,
        "urn_mismatch_examples": [],
        "mapped_without_metadata": 0,
        "mapped_not_list": 0,
        "post_not_object": 0,
        "bad_record": 0,
        "other_post_types": Counter(),
        "unnamed_targets": 0,
        "multi_target_mapped_records": 0,
        "multi_target_mapped_sets": 0,
        "null_variant_key_rows": 0,
        "null_protein_rows": 0,
        "null_chrom_rows": 0,
        "null_pos_rows": 0,
        "null_ref_rows": 0,
        "null_alt_rows": 0,
        "null_label_rows": 0,
        "null_target_rows": 0,
        "assay_source": Counter(),
        "state_vs_file": Counter(),
        "truth_variant_urn": "urn:mavedb:00000352-0-1#464",
        "truth_unmapped_urn": "urn:mavedb:00000004-a-1",
        "truth_hits": [],
        "writer": None,
    }
    notes = truth_mapped.get("notes") or ""
    if "urn:mavedb:00000352-0-1#464" not in notes:
        raise SystemExit("truth mapped notes did not contain the expected variant urn")
    if stats["truth_unmapped_urn"] not in (truth_unmapped.get("notes") or ""):
        raise SystemExit("truth unmapped notes did not contain the expected score-set urn")

    writer = Writer(VARIANTS_PATH, VARIANT_SCHEMA)
    stats["writer"] = writer
    try:
        experiments, coverage = build(acc_map, stats)
        write_experiments(experiments)
        write_coverage(coverage)
    finally:
        writer.close()

    checks = disk_checks(truth_mapped, truth_unmapped, stats)
    text = markdown(
        stats,
        experiments,
        coverage,
        protein,
        any_biotype,
        biotype_counts,
        disagreements,
        ties,
        truth_mapped,
        truth_unmapped,
        checks,
    )
    QC_PATH.write_text(text, encoding="utf-8")
    print(checks["text"])
    print(checks["result"])
    if not checks["passed"]:
        raise SystemExit(1)


SCORES_SCHEMA = pa.schema(
    [
        ("score_set_urn", pa.string()),
        ("accession", pa.string()),
        ("source_header", pa.string()),
        ("score", pa.float64()),
        ("score_text", pa.string()),
        ("sigma", pa.float64()),
        ("sigma_text", pa.string()),
        ("sigma_column", pa.string()),
        ("raw_score", pa.float64()),
        ("raw_score_text", pa.string()),
        ("raw_score_column", pa.string()),
        ("raw_sigma", pa.float64()),
        ("raw_sigma_text", pa.string()),
        ("raw_sigma_column", pa.string()),
        ("extra_scores", pa.string()),
    ]
)


def plain_float(text):
    """A float only when the whole cell is a plain number. NA stays null."""
    if text is None:
        return None
    stripped = text.strip()
    if not stripped or not _PLAIN_NUMBER.fullmatch(stripped):
        return None
    value = float(stripped)
    if value != value or value in {float("inf"), float("-inf")}:
        return None
    return value


def urn_from_scores_filename(base: str) -> str:
    suffix = ".scores.csv"
    if not base.endswith(suffix):
        raise ValueError(base)
    stem = base[: -len(suffix)]
    prefix = "urn-mavedb-"
    if not stem.startswith(prefix):
        raise ValueError(base)
    return "urn:mavedb:" + stem[len(prefix) :]


def _cell(record, index):
    if index is None or index >= len(record):
        return None
    return record[index]


def load_scores():
    """Read score CSVs. Join key is accession, the same string as mapped urn."""
    import duckdb

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    writer = Writer(SCORES_PATH, SCORES_SCHEMA)
    stats = Counter()
    truth_score = None
    try:
        with tarfile.open(TAR, "r:gz") as archive:
            for member in archive:
                name = member.name
                base = name.rsplit("/", 1)[-1]
                if is_appledouble(name):
                    stats["appledouble"] += 1
                    continue
                if not base.endswith(".scores.csv") or not member.isfile():
                    continue
                stats["score_files"] += 1
                score_urn = urn_from_scores_filename(base)
                text = io.TextIOWrapper(archive.extractfile(member), encoding="utf-8", newline="")
                reader = csv.reader(text)
                try:
                    header = next(reader)
                except StopIteration:
                    stats["empty_file"] += 1
                    continue
                source_header = ",".join(header)
                index = {column: i for i, column in enumerate(header)}
                if "accession" not in index or "scores.score" not in index:
                    stats["bad_header"] += 1
                    for _ in reader:
                        stats["skipped_rows"] += 1
                    continue
                if "scores.sigma" in index:
                    sigma_at, sigma_column = index["scores.sigma"], "scores.sigma"
                elif "scores.Sigma" in index:
                    sigma_at, sigma_column = index["scores.Sigma"], "scores.Sigma"
                else:
                    sigma_at, sigma_column = None, None
                if "scores.raw_score" in index:
                    raw_at, raw_column = index["scores.raw_score"], "scores.raw_score"
                elif "scores.score_raw" in index:
                    raw_at, raw_column = index["scores.score_raw"], "scores.score_raw"
                else:
                    raw_at, raw_column = None, None
                raw_sigma_at = index.get("scores.raw_sigma")
                raw_sigma_column = "scores.raw_sigma" if raw_sigma_at is not None else None
                lifted = {sigma_column, raw_column, raw_sigma_column, "scores.score"}
                extra_indexes = [
                    (column, i)
                    for column, i in index.items()
                    if column.startswith("scores.") and column not in lifted
                ]
                if sigma_column:
                    stats["files_with_sigma"] += 1
                if raw_column:
                    stats["files_with_raw_score"] += 1
                if raw_sigma_column:
                    stats["files_with_raw_sigma"] += 1
                if extra_indexes:
                    stats["files_with_extra"] += 1
                for record in reader:
                    stats["score_rows"] += 1
                    accession = _cell(record, index["accession"])
                    if not accession or accession == "NA":
                        stats["empty_accession"] += 1
                    score_text = _cell(record, index["scores.score"])
                    score = plain_float(score_text)
                    if score is None:
                        stats["score_not_numeric"] += 1
                    else:
                        stats["score_numeric"] += 1
                    sigma_text = _cell(record, sigma_at)
                    raw_text = _cell(record, raw_at)
                    raw_sigma_text = _cell(record, raw_sigma_at)
                    extra = {}
                    for column, i in extra_indexes:
                        extra[column] = _cell(record, i)
                    extra_text = json.dumps(extra, ensure_ascii=False, separators=(",", ":")) if extra else None
                    if accession == "urn:mavedb:00000352-0-1#464":
                        truth_score = score_text
                    writer.add(
                        {
                            "score_set_urn": score_urn,
                            "accession": accession,
                            "source_header": source_header,
                            "score": score,
                            "score_text": score_text,
                            "sigma": plain_float(sigma_text),
                            "sigma_text": sigma_text,
                            "sigma_column": sigma_column,
                            "raw_score": plain_float(raw_text),
                            "raw_score_text": raw_text,
                            "raw_score_column": raw_column,
                            "raw_sigma": plain_float(raw_sigma_text),
                            "raw_sigma_text": raw_sigma_text,
                            "raw_sigma_column": raw_sigma_column,
                            "extra_scores": extra_text,
                        }
                    )
                if stats["score_files"] % 400 == 0:
                    print(f"score files {stats['score_files']} rows {stats['score_rows']}", file=sys.stderr, flush=True)
    finally:
        writer.close()

    if writer.rows != stats["score_rows"]:
        raise SystemExit(f"score writer rows {writer.rows} != {stats['score_rows']}")

    scores = SCORES_PATH.as_posix()
    mapped = VARIANTS_PATH.as_posix()
    experiments = EXPERIMENTS_PATH.as_posix()
    scored = SCORED_VARIANTS_PATH.as_posix()
    con = duckdb.connect()
    con.execute(
        f"""
        CREATE TEMP TABLE attach AS
        SELECT
          score_set_urn,
          accession,
          any_value(score) AS score,
          any_value(score_text) AS score_text,
          any_value(sigma) AS sigma,
          any_value(sigma_text) AS sigma_text,
          any_value(sigma_column) AS sigma_column,
          any_value(raw_score) AS raw_score,
          any_value(raw_score_text) AS raw_score_text,
          any_value(raw_score_column) AS raw_score_column,
          any_value(raw_sigma) AS raw_sigma,
          any_value(raw_sigma_text) AS raw_sigma_text,
          any_value(raw_sigma_column) AS raw_sigma_column,
          any_value(extra_scores) AS extra_scores,
          any_value(source_header) AS source_header
        FROM read_parquet('{scores}')
        GROUP BY score_set_urn, accession
        HAVING count(*) = 1
            OR (
              count(DISTINCT coalesce(score_text, '')) = 1
              AND count(DISTINCT coalesce(sigma_text, '')) = 1
              AND count(DISTINCT coalesce(raw_score_text, '')) = 1
              AND count(DISTINCT coalesce(raw_sigma_text, '')) = 1
              AND count(DISTINCT coalesce(extra_scores, '')) = 1
            )
        """
    )
    duplicate_keys, conflicting_keys = con.execute(
        f"""
        SELECT
          count(*) FILTER (WHERE n > 1),
          count(*) FILTER (WHERE n > 1 AND n_score > 1)
        FROM (
          SELECT accession, score_set_urn, count(*) AS n,
                 count(DISTINCT coalesce(score_text, '')) AS n_score
          FROM read_parquet('{scores}')
          GROUP BY accession, score_set_urn
        )
        """
    ).fetchone()
    attached_keys = con.execute("SELECT count(*) FROM attach").fetchone()[0]
    print("writing mapped_variants_scored.parquet", file=sys.stderr, flush=True)
    con.execute(
        f"""
        COPY (
          SELECT
            m.*,
            s.score,
            s.score_text,
            s.sigma,
            s.sigma_text,
            s.sigma_column,
            s.raw_score,
            s.raw_score_text,
            s.raw_score_column,
            s.raw_sigma,
            s.raw_sigma_text,
            s.raw_sigma_column,
            s.extra_scores,
            s.source_header
          FROM read_parquet('{mapped}') m
          LEFT JOIN attach s
            ON m.urn = s.accession
           AND m.score_set_urn = s.score_set_urn
        ) TO '{scored}' (FORMAT PARQUET, COMPRESSION SNAPPY)
        """
    )
    original_rows = pq.ParquetFile(VARIANTS_PATH).metadata.num_rows
    scored_rows = pq.ParquetFile(SCORED_VARIANTS_PATH).metadata.num_rows
    old_names = list(pq.read_schema(VARIANTS_PATH).names)
    new_names = list(pq.read_schema(SCORED_VARIANTS_PATH).names)
    if scored_rows != original_rows or scored_rows != 6350898 or new_names[: len(old_names)] != old_names:
        SCORED_VARIANTS_PATH.unlink(missing_ok=True)
        raise SystemExit(
            f"scored file rejected: rows {scored_rows} vs mapped {original_rows}; schema {new_names[:len(old_names)]}"
        )

    (
        score_rows,
        distinct_accessions,
        score_rows_matched,
        distinct_accessions_matched,
    ) = con.execute(
        f"""
        WITH mapped_keys AS (
          SELECT DISTINCT urn, score_set_urn FROM read_parquet('{mapped}')
        )
        SELECT
          (SELECT count(*) FROM read_parquet('{scores}')),
          (SELECT count(*) FROM (SELECT DISTINCT accession, score_set_urn FROM read_parquet('{scores}'))),
          (SELECT count(*) FROM read_parquet('{scores}') s
             JOIN mapped_keys m ON s.accession = m.urn AND s.score_set_urn = m.score_set_urn),
          (SELECT count(*) FROM (
             SELECT DISTINCT s.accession, s.score_set_urn
             FROM read_parquet('{scores}') s
             JOIN mapped_keys m ON s.accession = m.urn AND s.score_set_urn = m.score_set_urn
          ))
        """
    ).fetchone()
    mapped_rows, mapped_urns = con.execute(
        f"""
        SELECT count(*), count(DISTINCT urn)
        FROM read_parquet('{mapped}')
        """
    ).fetchone()
    mapped_rows_hit, mapped_urns_hit, mapped_numeric, mapped_text_only = con.execute(
        f"""
        SELECT
          count(*) FILTER (WHERE score IS NOT NULL OR score_text IS NOT NULL),
          count(DISTINCT urn) FILTER (WHERE score IS NOT NULL OR score_text IS NOT NULL),
          count(score),
          count(*) FILTER (WHERE score IS NULL AND score_text IS NOT NULL)
        FROM read_parquet('{scored}')
        """
    ).fetchone()
    mapped_set_scores, mapped_set_scores_hit, mapped_set_scores_miss = con.execute(
        f"""
        WITH mapped_sets AS (
          SELECT score_set_urn FROM read_parquet('{experiments}') WHERE mapped_file_exists
        ),
        mapped_keys AS (
          SELECT DISTINCT urn, score_set_urn FROM read_parquet('{mapped}')
        )
        SELECT
          (SELECT count(*) FROM read_parquet('{scores}') s JOIN mapped_sets e ON s.score_set_urn = e.score_set_urn),
          (SELECT count(*) FROM read_parquet('{scores}') s
             JOIN mapped_keys m ON s.accession = m.urn AND s.score_set_urn = m.score_set_urn),
          (SELECT count(*) FROM read_parquet('{scores}') s
             JOIN mapped_sets e ON s.score_set_urn = e.score_set_urn
             ANTI JOIN mapped_keys m ON s.accession = m.urn AND s.score_set_urn = m.score_set_urn)
        """
    ).fetchone()
    truth = con.execute(
        f"""
        SELECT variant_key, score, score_text, count(*)
        FROM read_parquet('{scored}')
        WHERE urn = 'urn:mavedb:00000352-0-1#464'
        GROUP BY variant_key, score, score_text
        """
    ).fetchall()
    e4b_mapped = con.execute(
        f"""
        SELECT count(*) FROM read_parquet('{scored}')
        WHERE score_set_urn = 'urn:mavedb:00000004-a-1'
        """
    ).fetchone()[0]
    e4b_scores = con.execute(
        f"""
        SELECT count(*), count(*) FILTER (WHERE score IS NOT NULL)
        FROM read_parquet('{scores}')
        WHERE score_set_urn = 'urn:mavedb:00000004-a-1'
        """
    ).fetchone()
    tmp_set = con.execute(
        f"""
        WITH keys AS (
          SELECT DISTINCT urn FROM read_parquet('{mapped}')
          WHERE score_set_urn = 'urn:mavedb:00001239-a-6'
        )
        SELECT
          (SELECT count(*) FROM keys),
          (SELECT count(*) FROM read_parquet('{scores}') WHERE score_set_urn = 'urn:mavedb:00001239-a-6'),
          (SELECT count(*) FROM read_parquet('{scores}') s
             JOIN keys k ON s.accession = k.urn
             WHERE s.score_set_urn = 'urn:mavedb:00001239-a-6')
        """
    ).fetchone()

    lines = [
        "",
        "## Score join",
        "",
        "Mapped JSON still has no score. Scores were read from `csv/*.scores.csv` without extracting the tar.",
        f"AppleDouble members skipped in this pass: {stats['appledouble']}. Score files read: {stats['score_files']}. Files missing `accession` or `scores.score`: {stats['bad_header']}.",
        "",
        "The column both sides share is `accession` in the CSV and `urn` in `mapped_variants.parquet` (the mapped `variantUrn`). It is not HGVS. For `urn:mavedb:00000352-0-1#464` the CSV `hgvs_pro` is `p.Asn22Ile` and the mapped `protein_change` is null, while `accession` equals the mapped urn.",
        "Score set `urn:mavedb:00001239-a-6` uses `tmp:` ids on both sides.",
        f"That set has {tmp_set[0]} mapped urns, {tmp_set[1]} score rows, and {tmp_set[2]} accessions equal to a mapped urn.",
        "README says score values are not normalized across score sets. No normalized score was invented.",
        "",
        f"`scores.score` is a float only when the cell is a plain number. Numeric: {stats['score_numeric']}. Not numeric, float left null, raw text kept: {stats['score_not_numeric']}.",
        f"Files with a sigma column (`scores.sigma` or `scores.Sigma`): {stats['files_with_sigma']}.",
        f"Files with a raw score column (`scores.raw_score` or `scores.score_raw`): {stats['files_with_raw_score']}.",
        f"Files with `scores.raw_sigma`: {stats['files_with_raw_sigma']}.",
        f"Files with other `scores.*` columns, stored as JSON in `extra_scores`: {stats['files_with_extra']}.",
        "Those other columns are not given their own fields. Empty accession cells: {empty}.".format(empty=stats["empty_accession"]),
        "",
        f"Score rows: {score_rows}. Distinct `(score_set_urn, accession)` keys: {distinct_accessions}.",
        f"Duplicate keys: {duplicate_keys}. Duplicate keys whose `scores.score` text disagrees: {conflicting_keys}.",
        f"Keys used for the attach (unique, or duplicates with identical score text): {attached_keys}. Conflicting keys were not attached.",
        "",
        f"Score rows whose accession equals a mapped urn in the same score set: {score_rows_matched}.",
        f"Score rows with no mapped urn: {score_rows - score_rows_matched}.",
        f"Distinct score keys that match a mapped urn: {distinct_accessions_matched}.",
        f"Distinct score keys with no mapped urn: {distinct_accessions - distinct_accessions_matched}.",
        "",
        f"Score rows in score sets that have a mapped file: {mapped_set_scores}.",
        f"Of those, accession matched a mapped urn: {mapped_set_scores_hit}.",
        f"Of those, accession missed: {mapped_set_scores_miss}.",
        "",
        f"Mapped rows: {mapped_rows}. Rows with a score cell attached: {mapped_rows_hit}. Rows with no score attached: {mapped_rows - mapped_rows_hit}.",
        f"Of those attached cells, {mapped_numeric} are plain numbers and {mapped_text_only} are not, so `score` is null and `score_text` keeps the cell.",
        f"Distinct mapped urns: {mapped_urns}. Urns with a score attached: {mapped_urns_hit}. Urns with no score: {mapped_urns - mapped_urns_hit}.",
        "A miss leaves `score` and `score_text` null. It is not filled from another column.",
        "",
        f"`mapped_variants.parquet` was not rewritten ({original_rows} rows). `mapped_variants_scored.parquet` has the same rows and the same columns, then the score columns. Scored rows: {scored_rows}.",
        "",
        f"Truth variant `urn:mavedb:00000352-0-1#464` in the scored file: {truth}.",
        f"CSV score text seen while reading: {truth_score}.",
        f"Unmapped score set `urn:mavedb:00000004-a-1` (E4B) has {e4b_scores[0]} score rows ({e4b_scores[1]} numeric) and {e4b_mapped} mapped rows. It still has no coordinate.",
        "",
    ]
    qc = QC_PATH.read_text(encoding="utf-8")
    pointer = (
        "## Scores\n\n"
        "Mapped JSON has no score field. Scores were loaded from `csv/*.scores.csv`. "
        "The join key both sides share is the CSV `accession`, the same string as the mapped variant urn. "
        "Scores are not normalized across score sets. Counts are in Score join.\n\n"
    )
    start = qc.find("## Scores\n")
    end = qc.find("\n## Mapped records\n")
    if start == -1 or end == -1:
        raise SystemExit("QC note is missing the Scores section")
    qc = qc[:start] + pointer + qc[end + 1 :]
    if "## Score join\n" in qc:
        qc = qc[: qc.find("\n## Score join\n")]
    if not qc.endswith("\n"):
        qc += "\n"
    qc += "\n".join(lines)
    QC_PATH.write_text(qc, encoding="utf-8")
    print(
        f"scores {score_rows} matched_rows {score_rows_matched} "
        f"mapped_hit {mapped_rows_hit}/{mapped_rows} truth {truth} csv {truth_score}"
    )
    if not truth or truth[0][0] != "chr5:146469691:A:T" or truth[0][2] != truth_score:
        raise SystemExit("truth variant score or key failed")
    if e4b_mapped != 0:
        raise SystemExit("E4B gained a mapped row")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--scores-only":
        load_scores()
    else:
        main()

