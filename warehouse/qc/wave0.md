# Wave 0

Checked 2026-11-03 after `logs/fetch.status` was `DONE`.

## Download

- Status file: `DONE`. The log ends in `DONE`.
- Genome chr1–22, X, and Y: every `.vcf.bgz` has a `.tbi`. None missing.
- Exome chr1–22, X, and Y: every `.vcf.bgz` has a `.tbi`, including chrY. None missing.
- Core files present: ClinVar, GENCODE v49, GWAS zip, gnomAD constraint, both GTEx TPM files, both QTL tars, MaveDB tar.
- Raw files: 105. All 105 match an `OK` or `SKIP` size in `logs/fetch.log`.
- Free space on D: at the manifest: 676,664,020,992 bytes. Total 2,097,150,947,328. Used 1,420,486,926,336.
- SQLite `source_file` rebuilt from those 105 files. The old 10-row table, which had recorded the chr1 genome VCF at 16,484,257,792 bytes, is gone. The file on disk is 44,078,103,205 bytes.

Manifest: `warehouse/manifest/files.tsv`. Disk note: `warehouse/manifest/disk.md`.

## Tools

- Python: `%USERPROFILE%\miniconda3\python.exe` [user name redacted for publication]
- duckdb 2.0.0
- pyarrow 26.0.1
- biopython 1.88
- `coords.py` checks passed for `1`, `chr1`, `CHR1`, `X`, `MT`, and `M`.
- bcftools and htslib were not installed. conda-forge has no win-64 bcftools package. pysam has no Windows wheel and did not build.
- `scripts/vcf_query.py` is the region reader. A chrY seek to the middle of the linear index returned `chrY 26250338 T C`, and seeking back to that position found the same record. The truth-set gnomAD rows below were read with the same function.

## Truth set

`warehouse/qc/truth_set.tsv` has 18 cases. Each row cites the file it came from. Positions were read from those files or from a tabix lookup. Pass/fail against warehouse tables waits for Lane H. The cases are:

| Case | Resolved value |
|---|---|
| coding_snp | `chr11:126275389:C:T`, FOXRED1, rs267606829. ClinVar pathogenic, two stars. Same alleles in the genome VCF and the exome VCF. No GWAS row with this current SNP id. |
| cftr_indel | `chr7:117559590:ATCT:A`, CFTR p.Phe508del. Practice guideline, four stars. Genome VCF has the same alleles. |
| gwas_intergenic | chr22:45,604,058, rs4823308, context `intergenic_variant`, trait URI `EFO_0004348`. No GENCODE CDS covers the position. |
| gwas_empty_pos | GWAS line 5, empty `CHR_POS`, study GCST004604. It must stay in associations and out of the clump. |
| mismatched_ref | `chr15:84799209:G:A`. ClinVar display alleles are `na`. The HGVS name is not the key. |
| af_pop_diff | `chr22:10510183:A:G`. `AF_nfe` 0.00000 and `AF_afr` 0.00229358. |
| mane_select | A1BG, `ENST00000263100`, LOEUF 1.34, decile 7. |
| loeuf_low_1..3 | XPO6 0.04, MED13L 0.041, EIF4G1 0.048. All decile 0. |
| loeuf_near_1_1..3 | CABLES2, CFAP69, CHI3L1. LOEUF 1.0, decile 5. |
| gtex_egene | `chr1_64764_C_T_b38` and the split columns both make `chr1:64764:C:T`. |
| mavedb_mapped | `chr5:146469691:A:T` from accession `NC_000005.10`, score set `urn:mavedb:00000352-0-1`, gene TCERG1. VRS start was converted by adding 1. |
| mavedb_unmapped | `urn:mavedb:00000004-a-1`, gene E4B. No coordinate. |
| chrom_crosswalk | ClinVar chromosome `1` at `chr1:243344276:AG:A`, SDCCAG8, joins GENCODE `ENSG00000054282.17` on chr1. |
| somatic_or_oncogenicity | NF1 `chr17:31334927:C:T`. Origin `germline/somatic`, significance Pathogenic, oncogenicity Likely oncogenic. It must not be counted as germline pathogenic. |

## Acceptance

- Manifest covers every file in `raw/`.
- A sites VCF can be queried by chromosome and position.
- `coords.py` canonicalizes chromosomes.
- The truth set exists and each row cites its source file.

Wave 1 may start. Phase 2 downloads stay off.
