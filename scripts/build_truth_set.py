# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Choose Wave 0 truth-set cases by reading the source files.

Positions come from those files and from tabix lookups. Nothing is typed in
from memory.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import sqlite3
import sys
import tarfile
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("VW_ROOT") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parent))

import coords  # noqa: E402
import vcf_query  # noqa: E402

RAW = ROOT / "raw"
OUT = ROOT / "warehouse" / "qc" / "truth_set.tsv"
CLINVAR = RAW / "clinvar" / "variant_summary.txt.gz"
GWAS_ZIP = RAW / "gwas" / "gwas-catalog-associations_ontology-annotated-full.zip"
CONSTRAINT = RAW / "gnomad" / "constraint" / "gnomad.v4.1.constraint_metrics.tsv"
GTF = RAW / "gencode" / "gencode.v49.annotation.gtf.gz"
EQTL = RAW / "gtex" / "GTEx_Analysis_v8_eQTL.tar"
MAVE = RAW / "mavedb" / "mavedb-dump.2026062418131.tar.gz"
DB = ROOT / "db" / "variants.sqlite"
GENOMES = RAW / "gnomad" / "genomes"
EXOMES = RAW / "gnomad" / "exomes"

COLUMNS = [
    "case_id",
    "requirement",
    "source_file",
    "source_locator",
    "chrom_raw",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "gene",
    "rsid",
    "expected",
    "notes",
]


def stars(status: str):
    text = (status or "").lower()
    if "practice guideline" in text:
        return 4
    if "expert panel" in text:
        return 3
    if "multiple submitters" in text and "no conflicts" in text:
        return 2
    if "single submitter" in text or "conflicting" in text:
        return 1
    if not text.strip() or "no assertion" in text or "no interpretation" in text:
        return 0
    return None


def blank(value: str) -> bool:
    return value is None or value.strip() in {"", ".", "-", "na", "NA"}


def row_out(**kwargs) -> dict:
    return {key: "" if kwargs.get(key) is None else str(kwargs.get(key)) for key in COLUMNS}


def scan_clinvar():
    coding = None
    indel = None
    mismatch = None
    somatic = None
    chrom1 = None
    accessions: dict[str, str] = {}
    with gzip.open(CLINVAR, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line_no, rec in enumerate(csv.reader(fh, delimiter="\t"), start=2):
            if len(rec) <= idx["Assembly"] or rec[idx["Assembly"]] != "GRCh38":
                continue
            chrom = rec[idx["Chromosome"]]
            accession = rec[idx["ChromosomeAccession"]]
            if accession and chrom and accession not in accessions:
                accessions[accession] = chrom
            gene = rec[idx["GeneSymbol"]]
            sig = rec[idx["ClinicalSignificance"]]
            origin = rec[idx["OriginSimple"]]
            onco = rec[idx["Oncogenicity"]]
            review = rec[idx["ReviewStatus"]]
            name = rec[idx["Name"]]
            ref = rec[idx["ReferenceAlleleVCF"]]
            alt = rec[idx["AlternateAlleleVCF"]]
            display_ref = rec[idx["ReferenceAllele"]]
            rsid = coords.norm_rsid(rec[idx["RS# (dbSNP)"]])
            key = coords.variant_key(chrom, rec[idx["PositionVCF"]], ref, alt)
            locator = f"line {line_no} AlleleID {rec[idx['AlleleID']]} VariationID {rec[idx['VariationID']]}"

            if (
                coding is None
                and rec[idx["Type"]] == "single nucleotide variant"
                and sig == "Pathogenic"
                and origin.lower() == "germline"
                and rsid
                and key
                and len(ref) == 1
                and len(alt) == 1
                and stars(review) is not None
                and stars(review) >= 2
            ):
                coding = {
                    "line": line_no,
                    "locator": locator,
                    "chrom": chrom,
                    "pos": rec[idx["PositionVCF"]],
                    "ref": ref,
                    "alt": alt,
                    "key": key,
                    "gene": gene,
                    "rsid": rsid,
                    "sig": sig,
                    "review": review,
                    "stars": stars(review),
                    "name": name,
                }
            upper_name = name.upper().replace(" ", "")
            if (
                indel is None
                and gene == "CFTR"
                and ("F508DEL" in upper_name or "PHE508DEL" in upper_name or "DELTAF508" in upper_name)
                and "Pathogenic" in sig
                and key
                and (len(ref) != 1 or len(alt) != 1)
            ):
                indel = {
                    "locator": locator,
                    "chrom": chrom,
                    "pos": rec[idx["PositionVCF"]],
                    "ref": ref,
                    "alt": alt,
                    "key": key,
                    "gene": gene,
                    "rsid": rsid or "",
                    "sig": sig,
                    "review": review,
                    "stars": stars(review),
                    "name": name,
                }
            display_missing = display_ref.strip().lower() in {"", ".", "-", "na"}
            if (
                mismatch is None
                and key
                and not blank(ref)
                and (display_missing or display_ref.upper() != ref.upper())
                and len(ref) == 1
                and len(alt) == 1
            ):
                mismatch = {
                    "locator": locator,
                    "chrom": chrom,
                    "pos": rec[idx["PositionVCF"]],
                    "ref": ref,
                    "alt": alt,
                    "display_ref": display_ref,
                    "display_alt": rec[idx["AlternateAllele"]],
                    "key": key,
                    "gene": gene,
                    "name": name,
                }
            if somatic is None and (
                origin.lower() == "somatic" or (onco and onco not in {"-", "not provided"})
            ):
                somatic = {
                    "locator": locator,
                    "chrom": chrom,
                    "pos": rec[idx["PositionVCF"]],
                    "ref": ref,
                    "alt": alt,
                    "key": key or "",
                    "gene": gene,
                    "sig": sig,
                    "origin": origin,
                    "onco": onco,
                    "review": review,
                }
            if chrom1 is None and chrom == "1" and gene and gene != "-" and key:
                chrom1 = {
                    "locator": locator,
                    "chrom": chrom,
                    "pos": rec[idx["PositionVCF"]],
                    "ref": ref,
                    "alt": alt,
                    "key": key,
                    "gene": gene,
                    "rsid": rsid or "",
                }
            if coding and indel and mismatch and somatic and chrom1 and len(accessions) >= 24:
                # Accession coverage can settle early. Keep scanning until the
                # five case rows exist; 24 accessions is the autosomes plus X/Y/MT.
                pass
            if coding and indel and mismatch and somatic and chrom1 and len(accessions) >= 24:
                break
    return coding, indel, mismatch, somatic, chrom1, accessions


def scan_gwas(rsid: str):
    hit = None
    intergenic = None
    empty = None
    with zipfile.ZipFile(GWAS_ZIP) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as raw:
            text = (line.decode("utf-8", "replace") for line in raw)
            header = next(text).rstrip("\n").split("\t")
            idx = {col: i for i, col in enumerate(header)}
            for line_no, line in enumerate(text, start=2):
                rec = line.rstrip("\n").split("\t")
                if len(rec) < len(header):
                    continue
                if hit is None and rsid and coords.norm_rsid(rec[idx["SNP_ID_CURRENT"]]) == rsid:
                    hit = (line_no, rec, idx)
                flag = rec[idx["INTERGENIC"]].strip().lower()
                pos = rec[idx["CHR_POS"]].strip()
                trait_uri = rec[idx["MAPPED_TRAIT_URI"]].strip()
                context = rec[idx["CONTEXT"]].strip().lower()
                if (
                    intergenic is None
                    and "intergenic" in context
                    and pos.isdigit()
                    and trait_uri
                    and rec[idx["CHR_ID"]].strip()
                ):
                    intergenic = (line_no, rec, idx)
                if empty is None and not pos:
                    empty = (line_no, rec, idx)
                if hit and intergenic and empty:
                    break
    return hit, intergenic, empty


def cds_contains(chrom: str, pos: int) -> bool:
    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5 or parts[2] != "CDS" or parts[0] != chrom:
                continue
            if int(parts[3]) <= pos <= int(parts[4]):
                return True
    return False


def gene_on_chrom(symbol: str, chrom: str):
    conn = sqlite3.connect(DB)
    row = conn.execute(
        "SELECT gene_id, chrom, start, end, biotype FROM gene WHERE symbol=? AND chrom=? ORDER BY CASE biotype WHEN 'protein_coding' THEN 0 ELSE 1 END LIMIT 1",
        (symbol, chrom),
    ).fetchone()
    conn.close()
    return row


def scan_constraint():
    mane = []
    with CONSTRAINT.open(encoding="utf-8", errors="replace", newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line_no, line in enumerate(fh, start=2):
            rec = line.rstrip("\n").split("\t")
            if rec[idx["mane_select"]].lower() != "true":
                continue
            gene_id = rec[idx["gene_id"]]
            if not gene_id.startswith("ENSG"):
                continue
            try:
                loeuf = float(rec[idx["lof.oe_ci.upper"]])
            except ValueError:
                continue
            mane.append(
                {
                    "line": line_no,
                    "gene": rec[idx["gene"]],
                    "gene_id": gene_id,
                    "transcript": rec[idx["transcript"]],
                    "loeuf": loeuf,
                    "decile": rec[idx["lof.oe_ci.upper_bin_decile"]],
                    "chrom": rec[idx["chromosome"]],
                }
            )
    by_gene = {}
    for row in mane:
        current = by_gene.get(row["gene_id"])
        if current is None or (row["transcript"].startswith("ENST") and not current["transcript"].startswith("ENST")):
            by_gene[row["gene_id"]] = row
    rows = list(by_gene.values())
    low = sorted(rows, key=lambda item: item["loeuf"])[:3]
    near = sorted(rows, key=lambda item: abs(item["loeuf"] - 1.0))[:3]
    return rows[0] if rows else None, low, near


def scan_gtex():
    with tarfile.open(EQTL, "r:") as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith(".egenes.txt.gz"))
        fh = gzip.GzipFile(fileobj=tf.extractfile(member))
        header = fh.readline().decode("utf-8").rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line_no, raw in enumerate(fh, start=2):
            rec = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            variant_id = rec[idx["variant_id"]]
            parts = variant_id.split("_")
            if len(parts) != 5 or parts[4] != "b38":
                continue
            chrom, pos, ref, alt, _build = parts
            split_key = coords.variant_key(rec[idx["chr"]], rec[idx["variant_pos"]], rec[idx["ref"]], rec[idx["alt"]])
            id_key = coords.variant_key(chrom, pos, ref, alt)
            return {
                "member": member.name,
                "line": line_no,
                "variant_id": variant_id,
                "gene": rec[idx["gene_id"]],
                "key": split_key,
                "id_key": id_key,
                "chrom": rec[idx["chr"]],
                "pos": rec[idx["variant_pos"]],
                "ref": rec[idx["ref"]],
                "alt": rec[idx["alt"]],
            }
    return None


def mapped_urn(name: str) -> str:
    base = Path(name).name
    if base.startswith("._") or "/._" in name.replace("\\", "/"):
        return ""
    marker = "urn-mavedb-"
    if marker not in base:
        return ""
    rest = base.split(marker, 1)[1].removesuffix(".mapped-variants.json")
    return "urn:mavedb:" + rest


def scan_mavedb(accessions: dict[str, str]):
    wanted = {accession: coords.canon_chrom(chrom) for accession, chrom in accessions.items()}
    with tarfile.open(MAVE, "r:gz") as tf:
        mapped_members = [
            m
            for m in tf.getmembers()
            if m.name.endswith(".mapped-variants.json") and "/._" not in m.name and not Path(m.name).name.startswith("._")
        ]
        urns = {mapped_urn(m.name): m for m in mapped_members}
        main = next(m for m in tf.getmembers() if m.name.endswith("/main.json") and not Path(m.name).name.startswith("._"))
        meta = json.load(tf.extractfile(main))
        unmapped = None
        gene_by_urn = {}
        for exp_set in meta["experimentSets"]:
            for exp in exp_set.get("experiments") or []:
                for score in exp.get("scoreSets") or []:
                    urn = score.get("urn") or ""
                    genes = score.get("targetGenes") or []
                    symbol = genes[0].get("name") if genes else ""
                    if urn:
                        gene_by_urn[urn] = symbol
                    if unmapped is None and urn and urn not in urns:
                        unmapped = {"urn": urn, "gene": symbol, "title": score.get("title") or ""}
        mapped_record = None
        for member_info in mapped_members:
            if member_info.size > 25_000_000:
                continue
            blob = tf.extractfile(member_info).read()
            if b'"label": "NC_' not in blob:
                continue
            payload = json.loads(blob)
            urn = mapped_urn(member_info.name)
            for record in payload:
                post = record.get("postMapped")
                if not isinstance(post, dict):
                    continue
                if post.get("type") == "Allele":
                    nodes = [post]
                elif post.get("type") == "Haplotype" and len(post.get("members") or []) == 1:
                    nodes = post["members"]
                else:
                    continue
                member = nodes[0]
                location = member.get("location") or {}
                label = (location.get("sequenceReference") or {}).get("label") or ""
                chrom = wanted.get(label)
                if chrom is None:
                    continue
                start = location.get("start")
                end = location.get("end")
                alt = (member.get("state") or {}).get("sequence") or ""
                ref = ""
                for ext in member.get("extensions") or []:
                    if ext.get("name") == "vrs_ref_allele_seq":
                        ref = ext.get("value") or ""
                if start is None or end is None or not ref or not alt:
                    continue
                if end - start != len(ref):
                    continue
                if set((ref + alt).upper()) - set("ACGTN"):
                    continue
                pos = int(start) + 1
                key = coords.variant_key(chrom, pos, ref, alt)
                if not key:
                    continue
                mapped_record = {
                    "urn": urn,
                    "gene": gene_by_urn.get(urn, ""),
                    "file": member_info.name,
                    "variant_urn": record.get("variantUrn") or "",
                    "accession": label,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "key": key,
                }
                break
            if mapped_record:
                break
    return mapped_record, unmapped


def gnomad_hit(kind: str, chrom: str, pos: int, ref: str, alt: str):
    folder = GENOMES if kind == "genomes" else EXOMES
    path = folder / f"gnomad.{kind}.v4.1.sites.{chrom}.vcf.bgz"
    if not path.exists():
        return None
    rows = vcf_query.query(path, chrom, pos, ref, alt)
    return rows[0] if rows else None


def main() -> int:
    print("scanning ClinVar", flush=True)
    coding, indel, mismatch, somatic, chrom1, accessions = scan_clinvar()
    print("clinvar cases", bool(coding), bool(indel), bool(mismatch), bool(somatic), bool(chrom1), "accessions", len(accessions), flush=True)
    print("scanning GWAS", flush=True)
    gwas_hit, intergenic, empty = scan_gwas(coding["rsid"] if coding else "")
    print("gwas", bool(gwas_hit), bool(intergenic), bool(empty), flush=True)
    print("scanning constraint, GTEx, MaveDB", flush=True)
    mane_example, low, near = scan_constraint()
    gtex = scan_gtex()
    mapped, unmapped = scan_mavedb(accessions)
    print("mapped", bool(mapped), "unmapped", bool(unmapped), "mane", bool(mane_example), flush=True)

    rows = []
    if coding:
        genome = gnomad_hit("genomes", coords.canon_chrom(coding["chrom"]), int(coding["pos"]), coding["ref"], coding["alt"])
        exome = gnomad_hit("exomes", coords.canon_chrom(coding["chrom"]), int(coding["pos"]), coding["ref"], coding["alt"])
        gwas_note = "no GWAS row whose SNP_ID_CURRENT is this rsID" if gwas_hit is None else (
            f"GWAS line {gwas_hit[0]} SNP_ID_CURRENT {gwas_hit[1][gwas_hit[2]['SNP_ID_CURRENT']]} trait {gwas_hit[1][gwas_hit[2]['DISEASE/TRAIT']]}"
        )
        genome_note = "absent from genome sites" if genome is None else f"genome id {genome['id']} ref {genome['ref']} alts {','.join(genome['alts'])}"
        exome_note = "absent from exome sites" if exome is None else f"exome id {exome['id']} ref {exome['ref']} alts {','.join(exome['alts'])}"
        rows.append(
            row_out(
                case_id="coding_snp",
                requirement="ClinVar, GWAS if present, and gnomAD agree on variant_key",
                source_file=str(CLINVAR),
                source_locator=coding["locator"],
                chrom_raw=coding["chrom"],
                pos=coding["pos"],
                ref=coding["ref"],
                alt=coding["alt"],
                variant_key=coding["key"],
                gene=coding["gene"],
                rsid=coding["rsid"],
                expected=coding["key"],
                notes=f"ClinVar {coding['sig']}; review {coding['review']}; stars {coding['stars']}; {gwas_note}; {genome_note}; {exome_note}",
            )
        )
    if indel:
        canon = coords.canon_chrom(indel["chrom"])
        genome = gnomad_hit("genomes", canon, int(indel["pos"]), indel["ref"], indel["alt"])
        note = "gnomAD genome miss" if genome is None else f"gnomAD genome ref {genome['ref']} alts {','.join(genome['alts'])}"
        rows.append(
            row_out(
                case_id="cftr_indel",
                requirement="Indel alleles match gnomAD when gnomAD has that allele, or the miss is recorded",
                source_file=str(CLINVAR),
                source_locator=indel["locator"],
                chrom_raw=indel["chrom"],
                pos=indel["pos"],
                ref=indel["ref"],
                alt=indel["alt"],
                variant_key=indel["key"],
                gene=indel["gene"],
                rsid=indel["rsid"],
                expected=f"stars {indel['stars']}; significance {indel['sig']}",
                notes=f"{indel['name']}; review {indel['review']}; {note}",
            )
        )
    if intergenic:
        rec, idx = intergenic[1], intergenic[2]
        chrom = coords.canon_chrom(rec[idx["CHR_ID"]])
        pos = int(rec[idx["CHR_POS"]])
        print(f"checking CDS overlap {chrom}:{pos}", flush=True)
        inside = cds_contains(chrom, pos) if chrom else True
        rows.append(
            row_out(
                case_id="gwas_intergenic",
                requirement="Lands outside coding exons and still has a trait URI",
                source_file=str(GWAS_ZIP),
                source_locator=f"line {intergenic[0]}",
                chrom_raw=rec[idx["CHR_ID"]],
                pos=pos,
                variant_key="",
                gene=rec[idx["MAPPED_GENE"]],
                rsid=coords.norm_rsid(rec[idx["SNP_ID_CURRENT"]]) or "",
                expected=rec[idx["MAPPED_TRAIT_URI"]],
                notes=f"INTERGENIC {rec[idx['INTERGENIC']]}; CONTEXT {rec[idx['CONTEXT']]}; CDS overlap {inside}; trait {rec[idx['DISEASE/TRAIT']]}",
            )
        )
    if empty:
        rec, idx = empty[1], empty[2]
        rows.append(
            row_out(
                case_id="gwas_empty_pos",
                requirement="Stays in the association table and stays out of the clump",
                source_file=str(GWAS_ZIP),
                source_locator=f"line {empty[0]}",
                chrom_raw=rec[idx["CHR_ID"]],
                pos="",
                gene=rec[idx["MAPPED_GENE"]],
                rsid=coords.norm_rsid(rec[idx["SNP_ID_CURRENT"]]) or "",
                expected="in associations; absent from loci",
                notes=f"study {rec[idx['STUDY ACCESSION']]}; pubmed {rec[idx['PUBMEDID']]}; trait uri {rec[idx['MAPPED_TRAIT_URI']]}; snps {rec[idx['SNPS']]}",
            )
        )
    if mismatch:
        rows.append(
            row_out(
                case_id="mismatched_ref",
                requirement="The key uses the VCF ref and alt, not the ClinVar display allele",
                source_file=str(CLINVAR),
                source_locator=mismatch["locator"],
                chrom_raw=mismatch["chrom"],
                pos=mismatch["pos"],
                ref=mismatch["ref"],
                alt=mismatch["alt"],
                variant_key=mismatch["key"],
                gene=mismatch["gene"],
                expected=mismatch["key"],
                notes=f"display ref {mismatch['display_ref']} alt {mismatch['display_alt']}; name {mismatch['name']}",
            )
        )
    print("looking for AF_nfe != AF_afr", flush=True)
    af = vcf_query.find_af_disagreement(GENOMES / "gnomad.genomes.v4.1.sites.chr22.vcf.bgz", "chr22")
    if af:
        rows.append(
            row_out(
                case_id="af_pop_diff",
                requirement="AF_nfe and AF_afr are both stored and neither overwrites the other",
                source_file=af["source"],
                source_locator=f"{af['chrom']}:{af['pos']}",
                chrom_raw=af["chrom"],
                pos=af["pos"],
                ref=af["ref"],
                alt=af["alt"],
                variant_key=coords.variant_key(af["chrom"], af["pos"], af["ref"], af["alt"]),
                expected=f"AF_nfe {af['AF_nfe']}; AF_afr {af['AF_afr']}",
                notes="Sampled chr22 genome sites windows until the two fields differed",
            )
        )
    if mane_example:
        rows.append(
            row_out(
                case_id="mane_select",
                requirement="The gene-constraint row uses the MANE Select transcript",
                source_file=str(CONSTRAINT),
                source_locator=f"line {mane_example['line']}",
                chrom_raw=mane_example["chrom"],
                gene=mane_example["gene"],
                variant_key="",
                expected=mane_example["transcript"],
                notes=f"{mane_example['gene_id']}; LOEUF {mane_example['loeuf']}; decile {mane_example['decile']}",
            )
        )
    for label, group in (("loeuf_low", low), ("loeuf_near_1", near)):
        for n, item in enumerate(group, start=1):
            rows.append(
                row_out(
                    case_id=f"{label}_{n}",
                    requirement="Gene mart carries this LOEUF decile",
                    source_file=str(CONSTRAINT),
                    source_locator=f"line {item['line']}",
                    chrom_raw=item["chrom"],
                    gene=item["gene"],
                    expected=f"LOEUF {item['loeuf']}; decile {item['decile']}; transcript {item['transcript']}",
                    notes=item["gene_id"],
                )
            )
    if gtex:
        rows.append(
            row_out(
                case_id="gtex_egene",
                requirement="chr_pos_ref_alt_b38 becomes the same variant key as the split columns",
                source_file=str(EQTL),
                source_locator=f"{gtex['member']} line {gtex['line']}",
                chrom_raw=gtex["chrom"],
                pos=gtex["pos"],
                ref=gtex["ref"],
                alt=gtex["alt"],
                variant_key=gtex["key"],
                gene=gtex["gene"],
                expected=gtex["id_key"],
                notes=f"variant_id {gtex['variant_id']}; keys agree {gtex['key'] == gtex['id_key']}",
            )
        )
    if mapped:
        rows.append(
            row_out(
                case_id="mavedb_mapped",
                requirement="A mapped GRCh38 variant gets a variant key",
                source_file=str(MAVE),
                source_locator=mapped["file"],
                chrom_raw=mapped["accession"],
                pos=mapped["pos"],
                ref=mapped["ref"],
                alt=mapped["alt"],
                variant_key=mapped["key"],
                gene=mapped["gene"],
                expected=mapped["key"],
                notes=f"score set {mapped['urn']}; variant {mapped['variant_urn']}; VRS start converted by adding 1",
            )
        )
    if unmapped:
        rows.append(
            row_out(
                case_id="mavedb_unmapped",
                requirement="An unmapped score stays at gene level and does not get a coordinate",
                source_file=str(MAVE),
                source_locator="main.json",
                gene=unmapped["gene"],
                expected="no variant_key",
                notes=f"score set {unmapped['urn']}; title {unmapped['title']}",
            )
        )
    if chrom1:
        canon = coords.canon_chrom(chrom1["chrom"])
        gene = gene_on_chrom(chrom1["gene"], canon)
        rows.append(
            row_out(
                case_id="chrom_crosswalk",
                requirement="ClinVar chromosome 1 joins a GENCODE gene on chr1",
                source_file=str(CLINVAR),
                source_locator=chrom1["locator"],
                chrom_raw=chrom1["chrom"],
                pos=chrom1["pos"],
                ref=chrom1["ref"],
                alt=chrom1["alt"],
                variant_key=chrom1["key"],
                gene=chrom1["gene"],
                rsid=chrom1["rsid"],
                expected=canon,
                notes="no GENCODE gene" if gene is None else f"GENCODE {gene[0]} {gene[1]}:{gene[2]}-{gene[3]} {gene[4]}",
            )
        )
    if somatic:
        rows.append(
            row_out(
                case_id="somatic_or_oncogenicity",
                requirement="Does not count as germline pathogenic",
                source_file=str(CLINVAR),
                source_locator=somatic["locator"],
                chrom_raw=somatic["chrom"],
                pos=somatic["pos"],
                ref=somatic["ref"],
                alt=somatic["alt"],
                variant_key=somatic["key"],
                gene=somatic["gene"],
                expected="exclude from germline pathogenic counts",
                notes=f"origin {somatic['origin']}; significance {somatic['sig']}; oncogenicity {somatic['onco']}; review {somatic['review']}",
            )
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} cases to {OUT}")
    for row in rows:
        print(row["case_id"], row["variant_key"] or row["gene"], row["expected"][:120])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
