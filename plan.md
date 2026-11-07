# Variant workbench plan

This is the handoff for a fresh agent. The goal is a queryable GRCh38 warehouse of genes, studied variants, and the evidence that ties them to traits, molecular effects, and population frequencies. Build that from the files already under `D:\variant-workbench`. The download has finished. Start at [After the download](#after-the-download). The wave sections are the specifications for those steps. Re-verify `logs/fetch.status` and free space before trusting the completion note.

The owner wants to understand traits, including polygenic ones, from public genomic catalogs. The deliverable is the warehouse plus written answers to the question tracks. It is a model of public evidence. It is not a complete causal map of the genome, and it is not a score for any person.

## Read this first

Work only under `D:\variant-workbench`. C: and D: are partitions of one 4 TB SSD [model redacted for publication]. The Grok process may be started in a Windows system folder [exact folder redacted for publication]. Do not write project files there.

Python is `%USERPROFILE%\miniconda3\python.exe` (3.13) [user name redacted for publication]. The `py` launcher also exists. pandas and numpy import. As of Wave 0, duckdb 2.0.0 and pyarrow 26.0.1 import, and biopython 1.88 is installed for BGZF seeks.

`conda install -c conda-forge bcftools htslib` failed on this machine. bcftools is not published for win-64 on conda-forge, and the bioconda build is linux and macOS only. `pip install pysam` has no Windows wheel, and building it from source failed. Do not spend a lane reinstalling bcftools. Region queries go through `scripts/vcf_query.py`, which reads the TBI index and seeks the bgzip VCF. Wave 0 used it for the chrY round trip and for the truth-set gnomAD lookups.

Existing code:

- `scripts/fetch_sources.py` downloads public GRCh38 sources. It keeps an 80 GB free-space reserve and refuses a file that would break that reserve. Log is `logs/fetch.log`. Milestone status is `logs/fetch.status`.
- `scripts/load_core.py` loads a thin SQLite catalog at `db/variants.sqlite`. It already did its job. Extend it or replace the thin tables with the Parquet layouts below. Do not treat the current SQLite columns as the full source schemas.

Machine context that matters for this job: a 12-core desktop CPU, about 127 GB RAM, and a desktop GPU [hardware models redacted for publication]. Bioinformatics here is CPU, disk, and RAM. Do not plan on GPU acceleration.

## Rules

1. Assembly is GRCh38 for every coordinate. If a row has no assembly, prove the coordinate on the truth set before joining it.
2. Canonical chromosome names are `chr1` through `chr22`, `chrX`, `chrY`. Keep the raw value in `chrom_raw`. GENCODE and gnomAD already use the `chr` prefix. ClinVar and the GWAS Catalog do not.
3. The variant key is `chrom:pos:ref:alt` with `pos` as a 1-based VCF position and alleles uppercase. rsIDs are a second key, never the only key. Indels stay in VCF normalization. Do not invent a left-align step unless a truth-set indel fails to match.
4. Ensembl gene ids are joined with the version stripped (`ENSG00000123456.15` and `ENSG00000123456` match). Keep the versioned id beside the bare id.
5. gnomAD sites VCFs stay bgzip+tabix. Query them. Do not insert site rows into SQLite or materialize an uncompressed copy.
6. Genome allele counts and exome allele counts stay in separate columns. Adding them together is wrong. The sample sets differ.
7. Evidence classes stay distinct. A ClinVar label, a GWAS association, an eQTL, an sQTL, a MaveDB score, and a constraint metric are different claims. The ledger records the class on every row.
8. Wave 1 lanes write their own directories. A lane does not edit another lane's files. One integrator builds marts in Wave 2.
9. Every lane runs the truth set and writes a QC note with row counts, null rates, and join misses.
10. Do not start Phase 2 downloads in this run. The exome fetch has finished. Phase 2 still waits until Wave 2 is accepted and free space is measured again.
11. This warehouse describes public biology. Do not build a tool that ranks people, scores a person who is not in these files, or treats a behavioral polygenic association as a reading of someone's nature. Behavioral and cognitive GWAS hits stay in the catalog as weak observational claims. Population-frequency work reports coverage and missingness. It does not rank groups.

## What this disk can and cannot support

These files can support:

- A gene-level model: identity, canonical transcript, constraint, expression across GTEx tissues, ClinVar burden, published GWAS overlap, MaveDB coverage, eGene and sGene status.
- A variant-level model of variants someone has already studied or measured.
- Population frequencies and in-silico scores for those variants, pulled out of the gnomAD sites VCFs. Those VCFs already contain the scores. This is extraction, not a new annotation run.
- An evidence ledger that can answer which genes sit at the intersection of rare-disease evidence, common-variant associations, and a molecular effect.

These files cannot support:

- An explanation of every base. Most of GRCh38 will fall outside gene exons, lead QTLs, and GWAS loci. Measure that fraction and publish it.
- Linkage disequilibrium or a real polygenic score. Sites VCFs have no individual genotypes. The GWAS Catalog stores reported significant associations, not genome-wide summary statistics.
- A personal genome report. No individual genotypes are in this directory.
- Colocalization in the statistical sense. Overlap between a GWAS mapped gene and an eGene in a relevant tissue is a hypothesis. Colocalization needs full summary statistics.
- A complete functional map. MaveDB covers a small slice of proteins. Say so whenever a result depends on it.

## Snapshot (re-verify)

Observed 2026-11-03 while exome downloads were still running. The SQLite row counts below were rechecked the same evening and still match. Disk sizes, free space, and exome progress do not. Use the triage remeasure for those.

## Triage remeasure

Checked later on 2026-11-03. `fetch_sources.py all` was still running (a python process under a parent PowerShell process; process IDs redacted for publication). Do not kill it. `logs/fetch.status` was still `PHASE1_DONE`. In this script that line is written when the small core catalogs finish, before the gnomAD loop. `DONE` is written only after genomes and exomes both finish.

D: free space was 632,980,602,880 bytes (`Get-PSDrive`) and 632,832,360,448 bytes (`Get-Volume SizeRemaining`) a moment apart. About 589 GiB. The older 322 GiB figure does not match this volume reading. The rest of D: was not inventoried, so the increase is unexplained here. The fetch reserve is still 80 GiB. This free-space reading leaves room for the remaining exome chromosomes.

| Group | Files | Bytes |
|---|---:|---:|
| `raw/gnomad/genomes` | 48 | 563,052,329,190 |
| `raw/gnomad/exomes` | 37 | 167,887,914,439 |
| `raw/gtex` | 4 | 3,764,758,198 |
| `raw/mavedb` | 1 | 1,907,736,956 |
| `db/variants.sqlite` | 1 | 1,013,321,728 |
| `raw/clinvar` | 1 | 457,862,838 |
| `raw/gnomad/constraint` | 1 | 95,546,041 |
| `raw/gencode` | 1 | 93,374,019 |
| `raw/gwas` | 1 | 72,739,430 |

Genome sites chr1–22, X, and Y are each present as `.vcf.bgz` plus `.tbi`. Exome sites chr1–18 are present with indexes. The log's last line was `GET gnomad_exomes_chr19` for 12,620,525,657 bytes. At the directory listing, `gnomad.exomes.v4.1.sites.chr19.vcf.bgz` was 0 bytes and had no `.tbi`. Exome chr20, chr21, chr22, chrX, and chrY had not been requested yet. The fetch list includes exome chrY. Whether that URL exists was not checked.

SQLite counts rechecked: `gene` 78,691; `clinvar` 4,639,270; `gwas_hit` 1,205,644; `gene_constraint` 211,523; `source_file` 10. `source_file` is still stale. It records the chr1 genome VCF at 16,484,257,792 bytes. The file on disk is 44,078,103,205 bytes. The other genome VCFs are still unregistered.

Protein-coding genes: 20,097. By exact symbol against the thin tables: 17,895 have a constraint row, 18,868 have a ClinVar row, 14,876 have a GWAS row whose `mapped_gene` is exactly that symbol. Exact symbol undercounts GWAS. `BRCA1` has 26 exact rows and 31 rows whose mapped-gene text contains `BRCA1` (`BRCA1, NBR2` and `NBR2, BRCA1`). 211,802 GWAS rows have an empty mapped gene. 216,880 have `pos` null or 0. Distinct author trait strings: 54,630. Constraint is 18,204 distinct gene symbols and 211,523 transcript rows, with 36,315 marked `canonical=true`, because RefSeq and Ensembl copies are both stored. Canonical A1BG appears twice (`NM_130786.4` and `ENST00000263100`).

Header checks, not full loads:

- Constraint TSV has 55 columns, including `gene_id`, `mane_select`, and `level`.
- GTEx median TPM GCT is version 1.2, dimensions `56200` by `54` tissues.
- eQTL tar: 98 members, 49 `*.egenes.txt.gz`, 49 `*.signif_variant_gene_pairs.txt.gz`. An eGene header has `variant_id` and also separate `chr`, `variant_pos`, `ref`, `alt`, and `rs_id_dbSNP151_GRCh38p7`. Prefer those columns over parsing the id alone.
- sQTL tar: 98 members, 49 `*.sgenes.txt.gz`, 49 `*.sqtl_signifpairs.txt.gz`. The pair names do not match the eQTL pair pattern. Follow the filenames in the tar.
- MaveDB tar: 13,887 members, `main.json` present, 2,340 paths with `mapped` in the name. AppleDouble members (`._*`) are present and stay skipped.
- GWAS zip contains one TSV, 38 header columns, 749,786,313 bytes uncompressed.

No Parquet warehouse, `coords.py`, DuckDB, PyArrow, or `bcftools` was created in this pass. The readable layer is still the design below. Wave 0 has not started.

### Download completion check

Rechecked after the triage listing. `logs/fetch.status` is `DONE`. The log ends in `DONE`. Genome sites chr1–22, X, and Y each have a `.vcf.bgz` and a `.tbi` (24 VCFs). Exome sites chr1–22, X, and Y each have a `.vcf.bgz` and a `.tbi`. Exome chrY is on disk at 113,561,604 bytes, with a 5,364-byte index. The log's `OK` lines for exome chr19 through chrY match those file sizes, including chr19 at 12,620,525,657 bytes. Free space was 676,664,074,240 bytes. The Python fetch process was gone. A PowerShell process whose command line still matched `fetch_sources.py` was still listed (process ID redacted for publication). It was not restarted and it was not killed. Treat the files and the `DONE` line as the completion evidence. Wave 0 still rechecks both before it writes the manifest.

SQLite `db/variants.sqlite`, loaded by `load_core.py`:

| Table | Rows | What was actually loaded |
|---|---:|---|
| `gene` | 78,691 | GENCODE v49 `gene` features only. 20,097 protein-coding, 34,880 lncRNA, plus pseudogenes and small RNAs |
| `clinvar` | 4,639,270 | GRCh38 rows from the weekly summary. Phenotype text truncated to 500 characters. Many source columns omitted |
| `gwas_hit` | 1,205,644 | 9 of 38 GWAS Catalog columns. p-values stored as text |
| `gene_constraint` | 211,523 | Transcript-level. RefSeq and Ensembl rows are both present, so this is not one row per gene |
| `source_file` | 10 | Stale. Only the chr1 genome VCF was registered. The other genome VCFs were already on disk |

ClinVar `clinical_significance` counts at load time: Uncertain significance 2,397,869; Likely benign 1,144,428; empty/`-` 244,654; Pathogenic 197,557; Benign 195,919; Conflicting 183,629; Likely pathogenic 132,766; Benign/Likely benign 83,396; Pathogenic/Likely pathogenic 45,648; not provided 7,065; drug response 1,889; other 1,519. These rows are not unique variants. The same allele repeats across phenotypes. Count distinct `VariationID` and `AlleleID`.

On-disk raw data:

| Path | Role |
|---|---|
| `raw/clinvar/variant_summary.txt.gz` | ClinVar weekly summary, GRCh37 and GRCh38 mixed. Filter `Assembly == GRCh38` |
| `raw/gencode/gencode.v49.annotation.gtf.gz` | GENCODE v49 gene models. Current loader keeps `gene` rows only |
| `raw/gwas/gwas-catalog-associations_ontology-annotated-full.zip` | One TSV: `gwas-catalog-download-associations-alt-full.tsv`. Release noted by the loader as 2026-11-02 |
| `raw/gnomad/constraint/gnomad.v4.1.constraint_metrics.tsv` | Gene and transcript constraint |
| `raw/gnomad/genomes/gnomad.genomes.v4.1.sites.chr*.vcf.bgz` | gnomAD v4.1 genome sites. Log showed chr1–22, X, and Y complete, each with a `.tbi` |
| `raw/gnomad/exomes/gnomad.exomes.v4.1.sites.chr*.vcf.bgz` | gnomAD v4.1 exome sites. At snapshot, chr1–8 were `OK` in the log and chr9 was mid-`GET` |
| `raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz` | Median TPM by tissue. Small. Load this |
| `raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` | Sample-level TPM, about 1.6 GB gzip. Defer |
| `raw/gtex/GTEx_Analysis_v8_eQTL.tar` | Per-tissue `*.egenes.txt.gz` and `*.signif_variant_gene_pairs.txt.gz` |
| `raw/gtex/GTEx_Analysis_v8_sQTL.tar` | Splicing equivalent of the eQTL tar |
| `raw/mavedb/mavedb-dump.2026062418131.tar.gz` | MaveDB CC0 bulk. `main.json` is about 21 MB of experiment metadata. Mapped genomic JSON lives under `mapped/` |

Approximate sizes when the gnomAD folder was about 581 GB: genomes on the order of 525 GB, exomes still growing, GTEx about 3.5 GB, MaveDB about 1.8 GB, ClinVar about 0.4 GB. D: had about 322 GB free (`345,652,195,328` bytes). The fetch script will stop rather than eat the last 80 GB.

A gnomAD genome sites header (chrY) contained 240 `INFO` fields and no sample columns. Fields to extract for the variant-of-interest set, not for every site:

- Frequency: `AF`, `AC`, `AN`, `nhomalt`, `AF_<pop>`, `AC_<pop>`, `AN_<pop>` for `afr`, `amr`, `asj`, `eas`, `fin`, `mid`, `nfe`, `sas`, `ami`, `remaining`. Also `faf95` and `AF_grpmax` / `grpmax` when present.
- Sex-stratified `AF_XX` and `AF_XY` only if a question asks. Skip them in the first extract.
- Quality and context: `FILTER` column, `allele_type`, `variant_type`, `lcr`, `segdup`, `non_par`.
- Predictors already in the VCF: `cadd_phred`, `revel_max`, `spliceai_ds_max`, `pangolin_largest_ds`, `phylop`, `sift_max`, `polyphen_max`.
- Consequence: the `vep` INFO field. Parse the worst consequence on the MANE transcript when that transcript is in the annotation. Keep the raw `vep` string only on Tier A, and only if it fits. Prefer parsed columns.
- Identity: `VRS_Allele_IDs` when present. The VCF `ID` column may be an rsID or `.`.

Skip histogram fields (`age_hist_*`, `gq_hist_*`, `dp_hist_*`, `ab_hist_*`) and allele-specific site quality (`AS_*`) in the first extract.

## Coordinate and identifier contract

Write this once, in `scripts/coords.py`, during Wave 0. Every later script imports it.

```text
canon_chrom(raw) -> chr1 | chr2 | ... | chrX | chrY | chrM | None
variant_key(chrom, pos, ref, alt) -> "chr1:12345:A:G"
bare_ensg(gene_id) -> "ENSG00000123456"
bare_enst(transcript_id) -> "ENST00000123456"
norm_rsid(value) -> "rs123" | None
```

`norm_rsid` accepts `rs123`, `RS123`, and a bare integer from ClinVar `RS# (dbSNP)` or GWAS `SNP_ID_CURRENT`. Empty, `.`, `-1`, and non-numeric junk become null.

Gene crosswalk columns, one row per versioned GENCODE gene id, plus rows that attach NCBI and RefSeq ids:

| Column | Source |
|---|---|
| `ensg`, `ensg_version` | GENCODE `gene_id` |
| `symbol` | GENCODE `gene_name` |
| `ncbi_gene_id` | ClinVar `GeneID`, many genes will have none |
| `hgnc_id` | ClinVar `HGNC_ID` and GTF `hgnc_id` when present |
| `mane_enst`, `mane_refseq` | gnomAD constraint `mane_select` and GTF tags |
| `chrom`, `start`, `end`, `strand`, `biotype` | GENCODE gene feature |

GTEx v8 gene ids are older than GENCODE v49. Join on bare ENSG and report genes that fail to match. Do not force a symbol join as the primary key. Symbols change and collide.

## Source schemas the current loader dropped

### ClinVar (`variant_summary.txt.gz`)

Header columns:

`AlleleID`, `Type`, `Name`, `GeneID`, `GeneSymbol`, `HGNC_ID`, `ClinicalSignificance`, `ClinSigSimple`, `LastEvaluated`, `RS# (dbSNP)`, `nsv/esv (dbVar)`, `RCVaccession`, `PhenotypeIDS`, `PhenotypeList`, `Origin`, `OriginSimple`, `Assembly`, `ChromosomeAccession`, `Chromosome`, `Start`, `Stop`, `ReferenceAllele`, `AlternateAllele`, `Cytogenetic`, `ReviewStatus`, `NumberSubmitters`, `Guidelines`, `TestedInGTR`, `OtherIDs`, `SubmitterCategories`, `VariationID`, `PositionVCF`, `ReferenceAlleleVCF`, `AlternateAlleleVCF`, `SomaticClinicalImpact`, `SomaticClinicalImpactLastEvaluated`, `ReviewStatusClinicalImpact`, `Oncogenicity`, `OncogenicityLastEvaluated`, `ReviewStatusOncogenicity`, `SCVsForAggregateGermlineClassification`, `SCVsForAggregateSomaticClinicalImpact`, `SCVsForAggregateOncogenicityClassification`.

Load GRCh38 only. Prefer `PositionVCF`, `ReferenceAlleleVCF`, and `AlternateAlleleVCF` for the variant key. Keep `Start`/`Stop` and the non-VCF alleles. Keep `PhenotypeIDS` intact. That column links to MedGen, OMIM, HPO, and Orphanet. Do not truncate `PhenotypeList`.

Star buckets from `ReviewStatus`, stored as an integer `review_stars`:

| Stars | Review status text contains |
|---:|---|
| 4 | practice guideline |
| 3 | reviewed by expert panel |
| 2 | multiple submitters, no conflicts |
| 1 | single submitter or conflicting classifications |
| 0 | no assertion, no interpretation, or status missing |

Match case-insensitively against the full ClinVar phrase. Write the raw `ReviewStatus` next to the bucket. If a phrase does not match this table, leave `review_stars` null and count those phrases in the QC note.

`OriginSimple` separates germline from somatic. Oncogenicity and somatic clinical impact are their own evidence, not germline pathogenicity.

### GWAS Catalog

38 columns in `gwas-catalog-download-associations-alt-full.tsv`:

`DATE ADDED TO CATALOG`, `PUBMEDID`, `FIRST AUTHOR`, `DATE`, `JOURNAL`, `LINK`, `STUDY`, `DISEASE/TRAIT`, `INITIAL SAMPLE SIZE`, `REPLICATION SAMPLE SIZE`, `REGION`, `CHR_ID`, `CHR_POS`, `REPORTED GENE(S)`, `MAPPED_GENE`, `UPSTREAM_GENE_ID`, `DOWNSTREAM_GENE_ID`, `SNP_GENE_IDS`, `UPSTREAM_GENE_DISTANCE`, `DOWNSTREAM_GENE_DISTANCE`, `STRONGEST SNP-RISK ALLELE`, `SNPS`, `MERGED`, `SNP_ID_CURRENT`, `CONTEXT`, `INTERGENIC`, `RISK ALLELE FREQUENCY`, `P-VALUE`, `PVALUE_MLOG`, `P-VALUE (TEXT)`, `OR or BETA`, `95% CI (TEXT)`, `PLATFORM [SNPS PASSING QC]`, `CNV`, `MAPPED_TRAIT`, `MAPPED_TRAIT_URI`, `STUDY ACCESSION`, `GENOTYPING TECHNOLOGY`.

There is no assembly column. Current Catalog `CHR_POS` values are GRCh38, and some rows are empty because the association is a haplotype, a CNV, or an interaction. Parse `P-VALUE` to float. Keep the raw text. `OR or BETA` is one mixed column. Keep it as text and also store a float when the cell is a plain number. Sample size and ancestry live in `INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE` as prose. Keep the prose. A later parse can pull ancestry tokens. Do not drop rows that fail the parse.

`MAPPED_TRAIT_URI` is the EFO id used to roll traits up. `DISEASE/TRAIT` is the author string and is messier.

Distance clump: within one `MAPPED_TRAIT_URI`, sort by chromosome and position and collapse associations within 1,000,000 bp of the current lead (smallest p-value in the window). Rows with null positions stay in the association table and stay out of the clump. CNV rows (`CNV == Y`) stay out of the SNP clump.

This file is a literature index of reported associations. It is the wrong file for polygenic scores. Say that in the GWAS QC note so later agents do not compute a fake score from lead SNPs.

### gnomAD constraint

Header:

`gene`, `gene_id`, `transcript`, `canonical`, `mane_select`, `lof_hc_lc.obs`, `lof_hc_lc.exp`, `lof_hc_lc.possible`, `lof_hc_lc.oe`, `lof_hc_lc.mu`, `lof_hc_lc.pLI`, `lof_hc_lc.pNull`, `lof_hc_lc.pRec`, `lof.obs`, `lof.exp`, `lof.possible`, `lof.oe`, `lof.mu`, `lof.pLI`, `lof.pNull`, `lof.pRec`, `lof.oe_ci.lower`, `lof.oe_ci.upper`, `lof.oe_ci.upper_rank`, `lof.oe_ci.upper_bin_decile`, `lof.z_raw`, `lof.z_score`, `mis.obs`, `mis.exp`, `mis.possible`, `mis.oe`, `mis.mu`, `mis.oe_ci.lower`, `mis.oe_ci.upper`, `mis.z_raw`, `mis.z_score`, `mis_pphen.obs`, `mis_pphen.exp`, `mis_pphen.possible`, `mis_pphen.oe`, `syn.obs`, `syn.exp`, `syn.possible`, `syn.oe`, `syn.mu`, `syn.oe_ci.lower`, `syn.oe_ci.upper`, `syn.z_raw`, `syn.z_score`, `constraint_flags`, `level`, `transcript_type`, `chromosome`, `cds_length`, `num_coding_exons`.

LOEUF is `lof.oe_ci.upper`. The decile is `lof.oe_ci.upper_bin_decile`. The current loader kept a thin subset and inserted both the RefSeq transcript and the Ensembl transcript when both were marked canonical. Emit:

- `warehouse/constraint/transcript_constraint.parquet`, all rows, numeric columns parsed to float.
- `warehouse/constraint/gene_constraint.parquet`, one row per `gene_id` (bare ENSG when `gene_id` is Ensembl). Preference order: `mane_select` true, then `canonical` true, then the transcript with the longest `cds_length`. QC the number of genes with no MANE row.

### GENCODE GTF

Feature types to keep: `gene`, `transcript`, `exon`, `CDS`, `UTR`, `stop_codon`. Attributes to parse when present: `gene_id`, `gene_name`, `gene_type`, `transcript_id`, `transcript_type`, `transcript_name`, `tag` (look for `MANE_Select`, `Ensembl_canonical`, `appris_principal`, and related `appris_*` tags), `hgnc_id`, `exon_number`. GTF coordinates are 1-based inclusive. Store them as `start` and `end` in that system. A VCF position `P` falls inside a feature when `start <= P <= end`.

### GTEx

Median TPM is a GCT. Two header lines, then a gene id column, a description column, then one column per tissue. Load gene × tissue as long Parquet: `ensg_version`, `ensg`, `tissue`, `median_tpm`.

eQTL tar members look like `GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt.gz` and `...signif_variant_gene_pairs.txt.gz`. sQTL summaries are `*.sgenes.txt.gz`. sQTL pair files are `*.sqtl_signifpairs.txt.gz`, not the eQTL pair name. Confirmed from the tar listing: 49 tissues each, 98 members in each tar. An eGene header already splits `chr`, `variant_pos`, `ref`, `alt`, and `rs_id_dbSNP151_GRCh38p7` beside `variant_id`. Use the split columns for the variant key. The id form `chr1_13550_G_A_b38` is a check, not the only parser. If a variant id does not match that pattern and the split columns are also unusable, keep the raw id and leave the key null.

Load all eGene and sGene files in Wave 1. Stream significant pairs to Parquet in the same wave if disk allows. Pairs are Tier B. They are not all joined to gnomAD in Wave 2.

The full sample TPM gzip waits until the gene mart exists and free space is still comfortable.

### MaveDB

The tar has `main.json`, `README.md`, `LICENSE.txt`, `mapped/*.mapped-variants.json`, plus `csv/`, `va/`, and other score payloads. Skip AppleDouble members whose names contain `/._` or start with `._`.

Wave 1 reads `README.md` and `main.json` first and writes an experiment inventory: score-set urn, target gene, assay type, genome build, and whether a mapped file exists. Then it extracts mapped variants. Unmapped scores stay attached to the gene and the urn. They do not get a fake genomic coordinate.

Many experiments are protein-coordinate only. Coverage of the proteome is sparse. The QC note states how many GENCODE protein-coding genes have at least one score set.

## Warehouse layout

```text
D:\variant-workbench\
  plan.md
  db\variants.sqlite
  raw\...
  scripts\
    coords.py
    build_manifest.py
    load_intervals.py
    load_constraint.py
    load_clinvar.py
    load_gwas.py
    load_gtex.py
    load_mavedb.py
    build_crosswalk.py
    build_gene_mart.py
    build_tier_a.py
    annotate_gnomad.py
    build_ledger.py
    truth_set.py
  warehouse\
    manifest\
    qc\
    crosswalk\
    intervals\
    constraint\
    clinvar\
    gwas\
    gtex\
    mavedb\
    variants\
    marts\
    ledger\
    reports\
```

Parquet is the analytical store. DuckDB is the query engine. SQLite remains the small catalog the existing script knows how to rebuild. After Wave 1, refresh `source_file` so every file on disk is registered. Do not pour the new wide tables back into SQLite unless a query truly needs it. DuckDB can attach the existing SQLite file for comparisons.

## Truth set

Build `warehouse/qc/truth_set.tsv` in Wave 0 and run it after every lane. Resolve coordinates from ClinVar and from a gnomAD tabix lookup. Do not hard-code positions from memory. Required cases:

| Case | What it must prove |
|---|---|
| A coding SNP with an rsID and a ClinVar interpretation | ClinVar, GWAS if present, and gnomAD agree on `variant_key` |
| A famous pathogenic coding allele | Review stars and clinical significance survive the reload |
| An indel, such as a well-known CFTR deletion if it is in ClinVar | Indel alleles match gnomAD when gnomAD has that allele, or the miss is recorded as a real miss |
| An intergenic GWAS association with a position | It lands outside coding exons and still has a trait URI |
| A GWAS row with an empty `CHR_POS` | It stays in the association table and stays out of the clump |
| A multi-allelic or mismatched ref site | The key uses the VCF ref and alt, not the ClinVar display name |
| A variant whose `AF_nfe` and `AF_afr` both exist and differ | Both frequencies are stored, and neither is overwritten by the other |
| A gene with MANE Select in the constraint file | The gene-constraint row uses that transcript |
| A GTEx eGene lead variant id of the form `chr_pos_ref_alt_b38` | It becomes the same variant key |
| A MaveDB mapped variant and a MaveDB unmapped score | The mapped one gets a key. The unmapped one stays gene-level |
| A ClinVar row with `Chromosome` `1` and a GENCODE gene on `chr1` | The crosswalk joins |
| A somatic or oncogenicity ClinVar row | It does not count as germline pathogenic |

Also pick three genes with very low LOEUF and three with LOEUF near 1 and confirm the gene mart carries the decile. Write pass/fail per case to `warehouse/qc/truth_set_report.md`. A lane is not done with a failing case that it was supposed to satisfy.

## After the download

Follow this section in order. The sections below it are the specifications. Do not start at a later wave because an earlier wave's files already exist in `raw/`. `raw/` is the download. `warehouse/` is the readable layer, and it does not exist yet.

The readable dataset has three layers:

1. Raw files stay where they are. gnomAD stays bgzip plus tabix. Do not insert site rows into SQLite, Parquet, or any other table. Do not write a second uncompressed copy.
2. Parquet tables under `warehouse/`, queried with DuckDB. One table family per evidence class, then one gene mart, one Tier A variant table, and one evidence ledger.
3. `db/variants.sqlite` stays the small catalog. Rebuild `source_file` from the manifest. Do not pour the wide tables back into it.

Four join keys, applied in every lane through `scripts/coords.py`:

| Key | Rule |
|---|---|
| Chromosome | `chr1`–`chr22`, `chrX`, `chrY`, `chrM`. Keep `chrom_raw`. ClinVar and GWAS use `1`. GENCODE and gnomAD use `chr1`. |
| Gene id | Strip the Ensembl version. Keep the versioned id beside the bare id. This is the primary gene join. |
| Symbol | Secondary. Split GWAS `MAPPED_GENE` on commas, trim tokens, and drop empty tokens and `-`. A symbol matches a token. It does not match a substring. `BRCA1` matches `BRCA1, NBR2`. It does not match a different gene that merely contains those letters. |
| Variant | `chrom:pos:ref:alt`, position 1-based, alleles uppercase. rsIDs are a second key. GWAS rows with no ref/alt keep a null `variant_key` and a `locus_key` of `chrom:pos`. Do not invent alleles from `STRONGEST SNP-RISK ALLELE`. |

Evidence classes stay in their own columns and their own ledger rows: ClinVar, GWAS, eQTL, sQTL, MaveDB, constraint. Genome allele counts and exome allele counts stay in separate columns.

### Step 0 — confirm the download

Stop if `logs/fetch.status` is not exactly `DONE`. Do not start a second fetch. Do not kill a fetch that is still running.

Then confirm all of these, and write the result at the top of `warehouse/qc/wave0.md`:

- Genome `chr1`–`chr22`, `chrX`, `chrY`: each `.vcf.bgz` has a `.tbi`, and the byte size matches an `OK` line in `logs/fetch.log`.
- Exome `chr1`–`chr22`, `chrX`, `chrY`: the same check. Exome chrY is expected. If a later recheck finds it missing, record the log line and continue with the chromosomes that exist. Do not block the warehouse on one missing index.
- Core files exist: ClinVar summary, GENCODE v49 GTF, GWAS zip, gnomAD constraint TSV, both GTEx TPM gzips, both GTEx QTL tars, MaveDB tar.
- Free space on D: is recorded in bytes. Leave the 80 GB fetch reserve untouched. Do not extract a second full copy of any gnomAD VCF. Stream GTEx and MaveDB members. If a load would drop free space below about 120 GB, stream it instead of extracting the archive.

The completion check above already passed once, on 2026-11-03. Step 0 still runs, because that is the handoff. The historical triage numbers are not a substitute for it.

### Step 1 — Wave 0, one agent

Specification: [Wave 0](#wave-0--one-agent-finish-before-fan-out). Finish before any lane starts.

1. Run Step 0.
2. Write `warehouse/manifest/disk.md` with free space.
3. Write `warehouse/manifest/files.parquet` or a TSV: path, byte size, and whether the log has a matching `OK` line. Cover every file in `raw/`.
4. Rebuild SQLite `source_file` from that manifest. The current table has 10 rows and records the chr1 genome VCF at 16,484,257,792 bytes. The file on disk at triage was 44,078,103,205 bytes. The manifest's size wins.
5. Install `bcftools` and `htslib` from conda-forge into the existing miniconda, plus Python packages `duckdb` and `pyarrow`. Confirm `bcftools --version` and one variant query against the chrY genome sites VCF.
6. Write `scripts/coords.py` and test `1`, `chr1`, `CHR1`, `X`, and `MT`/`M`.
7. Write `warehouse/qc/truth_set.tsv` by looking variants up in ClinVar and with tabix. Do not hard-code positions from memory. Cases are listed in the truth-set section. Each row cites the file it came from.
8. Write `warehouse/qc/wave0.md` with free space, the chromosome checklist, tool versions, and the truth-set size.

Wave 1 does not start until this note exists and the acceptance lines in the Wave 0 section pass.

### Step 2 — Wave 1, parallel lanes

Specification: [Wave 1](#wave-1--parallel-lanes). Each lane writes only its own `warehouse/<lane>/` directory and `warehouse/qc/<lane>.md`. Shared reads of `raw/` and `scripts/coords.py` are fine. Shared writes to SQLite are not.

| Lane | Script | Writes |
|---|---|---|
| A intervals | `scripts/load_intervals.py` | `genes`, `transcripts`, `exons`, `cds`, `utrs` Parquet. Features `gene`, `transcript`, `exon`, `CDS`, `UTR`, `stop_codon`. |
| B constraint | `scripts/load_constraint.py` | All transcript rows, then one row per gene. Preference: `mane_select`, then `canonical`, then longest `cds_length`. Canonical A1BG is stored twice in the raw file (`NM_130786.4` and `ENST00000263100`). Do not count `canonical=true` rows as genes. |
| C ClinVar | `scripts/load_clinvar.py` | GRCh38 alleles, full phenotype text and `PhenotypeIDS`, review stars, then one deduplicated row per `variant_key`. |
| D GWAS | `scripts/load_gwas.py` | All 38 columns, traits by EFO URI, and loci clumped within 1,000,000 bp inside one trait. Split mapped-gene tokens as specified above. Rows with no position stay in the association table and stay out of the clump. |
| E GTEx | `scripts/load_gtex.py` | Median TPM long table, per-gene tau, eGenes, sGenes, significant eQTL pairs, significant sQTL pairs. eQTL pairs are `*.signif_variant_gene_pairs.txt.gz`. sQTL pairs are `*.sqtl_signifpairs.txt.gz`. Build the variant key from `chr`, `variant_pos`, `ref`, and `alt` when those columns exist. Defer sample-level `gene_tpm.gct.gz`. |
| F MaveDB | `scripts/load_mavedb.py` | Experiment inventory from `main.json`, mapped variants, gene coverage. Skip AppleDouble members. Unmapped scores stay on the gene and the urn. |
| G crosswalk | `scripts/build_crosswalk.py` | One gene-id crosswalk. Start after Lane A, and rerun after Lanes A, B, and C exist. |
| H truth set | `scripts/truth_set.py` | `warehouse/qc/truth_set_report.md`. Unfinished lanes are skipped, not failed. Rerun when the other lanes have finished. |

Lane QC recomputes join coverage. The thin-SQLite baselines are a hint, not a pass line: 20,097 protein-coding genes, 17,895 with a constraint symbol match, 18,868 with a ClinVar symbol match, 14,876 with a GWAS mapped gene equal to the symbol. Better keys should change those rates. Report the new rates and the miss list.

### Step 3 — Wave 2, one integrator

Specification: [Wave 2](#wave-2--integrator-after-wave-1-outputs-exist). Start only after Lane H's finished report passes the cases those lanes were supposed to satisfy.

1. `warehouse/marts/gene_mart.parquet`, one row per GENCODE gene, plus the protein-coding copy. Columns cover identity, constraint, expression, ClinVar star buckets, both GWAS counts (`gwas_loci_by_overlap` and `gwas_loci_by_mapped_gene`), MaveDB, and eGene/sGene tissue counts.
2. Tier A keys: ClinVar variants with a key, GWAS variants that gained a key through an rsID match, MaveDB variants with a GRCh38 key, and lead eQTL/sQTL variants. Significant pairs that are not the lead stay in the pair files.
3. Annotate Tier A from the gnomAD sites VCFs. Run chrY first and check the truth set. Then one process per remaining chromosome. Write genome scores and exome scores to separate columns. A missing site stays in the table with `in_genomes = false`.
4. `warehouse/ledger/evidence.parquet` with all six evidence classes. Tier B QTL pairs stay out of the first ledger.
5. `warehouse/reports/dark_genome.md`. Chromosome lengths come from the gnomAD contig header. Report bases in protein-coding CDS, bases in any gene span, bases within 1 kb of a Tier A variant, bases in a clumped GWAS window, and the exclusive remainder in none of those sets. Overlaps are not exclusive.
6. `warehouse/qc/wave2.md` starts with a fresh free-space measurement and the exome checklist.

### Step 4 — Wave 3, sixteen read-only reports

Specification: [Wave 3](#wave-3--question-tracks). Each track writes `warehouse/reports/<track>.md` and the SQL or script that produced the numbers. No track rewrites a mart. A broken join goes to `warehouse/qc/issues.md` and that claim stops.

1. `rare_disease.md` — LOEUF, germline pathogenic alleles at two or more stars, expression breadth. Single-submitter and conflicting calls are split out.
2. `common_disease.md` — clumped loci, mapped genes, and eGene status in a tissue named for that trait. Call it a shared-gene hypothesis.
3. `coding_versus_noncoding.md` — missense, intronic, UTR, and intergenic shares of clumped GWAS loci.
4. `pleiotropy.md` — genes and Tier A variants linked to many traits. Show shared genes and shared lead variants separately.
5. `conflicts.md` — pathogenic alleles with high frequency, common GWAS risk alleles, high missense-z genes with mild MaveDB scores, conflicting ClinVar labels with high REVEL or SpliceAI.
6. `constraint_versus_function.md` — LOEUF decile against pathogenic fraction and against MaveDB effect sizes. State how few genes have an assay.
7. `splicing.md` — ClinVar splice-region variants, SpliceAI or Pangolin, and sQTLs for the same gene, in separate columns.
8. `tissue_programs.md` — expression of mapped genes across GTEx tissues, split by tau.
9. `ancestry_audit.md` — coverage of `AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, and `AF_amr` on Tier A GWAS variants. Report missingness. Do not rank groups.
10. `clinvar_versus_frequency.md` — clinical significance against genome AF and against exome AF. Check exceptions back to the source files.
11. `noncoding_genes.md` — lncRNA and other noncoding genes. Keep these conclusions out of the protein-coding rankings.
12. `drug_response.md` — the ClinVar drug-response set with constraint and expression. Say how small it is.
13. `somatic_versus_germline.md` — somatic origin and oncogenicity, held apart from the germline ranking.
14. `study_redundancy.md` — how far the association rows shrink after clumping.
15. `predictor_disagreement.md` — CADD, REVEL, SpliceAI, MaveDB, and ClinVar stars. Report missingness.
16. `genes_with_no_evidence.md` — genes with no expression, no constraint call, no ClinVar submission, no GWAS locus by either rule, and no MaveDB assay.

Every report starts with the question, the filters, the denominator, the result, and the caveat.

### Step 5 — stop

This run ends at the definition of done. Phase 2 is a later decision. It is not the next automatic download. Do not fetch ENCODE, HPO, Gene Ontology, Reactome, GWAS summary statistics, a linkage panel, gnomAD structural variants, OMIM full text, or extra pathogenicity scores while building Steps 1–4.

## Wave 0 — one agent, finish before fan-out

Goal: a frozen input manifest, the coordinate module, the tools, and the truth-set file.

Steps:

1. Run the Step 0 checklist in [After the download](#after-the-download). Re-read `logs/fetch.log` and list `raw/gnomad/exomes` and `raw/gnomad/genomes`. A chromosome is complete only when the log has an `OK` line and a `.tbi` sits beside the `.vcf.bgz`. If `fetch.status` is not `DONE`, stop. Do not kill a running fetch. Do not start a second fetch.
2. Measure free space on D:. Record it in `warehouse/manifest/disk.md`.
3. Write `warehouse/manifest/files.parquet` or a TSV with path, byte size, and whether the fetch log has a matching `OK <name> <bytes>` line. Include ClinVar, GENCODE, GWAS, constraint, GTEx, MaveDB, and every gnomAD VCF and index.
4. Rebuild the SQLite `source_file` table from that manifest. The current 10-row table is incomplete.
5. Install tooling into the existing miniconda, not into a new prefix unless conda requires it. Need `bcftools` and `htslib` (conda-forge) and Python packages `duckdb` and `pyarrow`. Confirm `bcftools --version` and a one-variant query against chrY, which is the small genome VCF (`raw/gnomad/genomes/gnomad.genomes.v4.1.sites.chrY.vcf.bgz`).
6. Write `scripts/coords.py` as specified above, with tests for `1`, `chr1`, `CHR1`, `X`, and `MT`/`M`.
7. Write the truth-set TSV by looking variants up, not by guessing coordinates. The first gnomAD probe can be a `bcftools view -H -r` on chrY or on a known ClinVar site once the key exists.
8. Write `warehouse/qc/wave0.md` with free space, exome progress, tool versions, and truth-set size.

Acceptance: manifest covers every file in `raw/`. `bcftools` can read a sites VCF. `coords.py` canonicalizes chromosomes. The truth set exists and each row cites the file it was resolved from. Wave 1 does not start without this note.

Disk budget: leave the fetch's 80 GB reserve alone. Do not extract a second full copy of any gnomAD VCF. GTEx tars and the MaveDB tar may be streamed, or extracted under `warehouse/` only for the pieces a loader needs. If extraction would drop free space below about 120 GB, stream instead of extracting.

## Wave 1 — parallel lanes

Each lane writes only its own `warehouse/<lane>/` directory plus `warehouse/qc/<lane>.md`. Shared reads from `raw/` and from `scripts/coords.py` are fine. Shared writes to SQLite are not, except the Wave 0 manifest agent already finished `source_file`.

### Lane A — intervals

Script: `scripts/load_intervals.py`.

Input: `raw/gencode/gencode.v49.annotation.gtf.gz`.

Output:

- `warehouse/intervals/genes.parquet`
- `warehouse/intervals/transcripts.parquet`
- `warehouse/intervals/exons.parquet`
- `warehouse/intervals/cds.parquet`
- `warehouse/intervals/utrs.parquet`

Required gene columns: `ensg`, `ensg_version`, `symbol`, `chrom`, `start`, `end`, `strand`, `biotype`, `mane_select` (bool), `hgnc_id`. Transcript columns add `enst`, `enst_version`, `transcript_type`, `tags`. Exon, CDS, and UTR rows carry `ensg`, `enst`, `chrom`, `start`, `end`, `strand`, `exon_number` when present.

QC: gene count near 78,691, protein-coding count near 20,097, every `chrom` uses the canonical prefix, no null gene id. Report MANE tag counts. If the GTF has no `MANE_Select` tag, say so and let the constraint file remain the MANE source. Do not invent tags.

### Lane B — constraint

Script: `scripts/load_constraint.py`.

Input: `raw/gnomad/constraint/gnomad.v4.1.constraint_metrics.tsv`.

Output: `warehouse/constraint/transcript_constraint.parquet` and `warehouse/constraint/gene_constraint.parquet` as specified above. Parse the numeric fields with a float parser. Keep `constraint_flags`.

QC: one row per gene in the gene file. `canonical=true` is not one row per gene. The raw file stores RefSeq and Ensembl copies, so canonical A1BG appears as both `NM_130786.4` and `ENST00000263100`. Count genes with a MANE transcript, genes with only a canonical transcript, and genes dropped by the preference rule. Join bare `gene_id` to Lane A when Lane A exists. This lane may finish first. The integrator checks the join later.

### Lane C — ClinVar

Script: `scripts/load_clinvar.py`.

Input: `raw/clinvar/variant_summary.txt.gz`.

Output:

- `warehouse/clinvar/alleles.parquet`, one row per GRCh38 record, wide columns listed in the schema section. Add `chrom`, `pos`, `ref`, `alt`, `variant_key`, `rsid`, `review_stars`.
- `warehouse/clinvar/allele_dedup.parquet`, one row per `variant_key` plus a row for records that have no VCF allele, keyed by `AlleleID`. Aggregate clinical significance, max review stars, phenotype ids, and a germline-versus-somatic flag.

QC: GRCh37 rows excluded. Distinct `VariationID` much smaller than 4.6 million raw rows. Star-bucket nulls listed by raw status. Truth-set pathogenic allele present.

### Lane D — GWAS

Script: `scripts/load_gwas.py`.

Input: the TSV inside the GWAS zip.

Output:

- `warehouse/gwas/associations.parquet`, all 38 columns plus `chrom`, `pos`, `rsid`, `p_value_float`, `effect_float` when parseable, `variant_key` when both chrom and pos exist. The Catalog often has no ref/alt. When ref and alt are absent, `variant_key` may be null and `locus_key` = `chrom:pos` is the join to later annotation. Record that limitation in the QC note. rsID plus a gnomAD or ClinVar lookup is how those rows gain alleles. Do not guess the alleles from the risk-allele field alone. `STRONGEST SNP-RISK ALLELE` looks like `rs128494-C`. Store `risk_allele` from that field. It is not a full REF/ALT pair.
- `warehouse/gwas/traits.parquet`, one row per `MAPPED_TRAIT_URI` with the mapped trait name, the number of associations, and the number of distinct author trait strings.
- `warehouse/gwas/loci.parquet`, distance-clumped within trait as specified above. Columns: `trait_uri`, `chrom`, `lead_pos`, `lead_rsid`, `lead_p`, `n_associations`, `start`, `end`, `mapped_genes`. `mapped_genes` is the raw field plus a token list split on commas, trimmed, with empty tokens and `-` removed. Gene matches use those tokens. A substring match is wrong.

QC: row count near 1,205,644. Fraction with null position. Fraction intergenic. Number of traits. Number of clumps versus number of raw rows, so the redundancy is visible. Truth-set intergenic hit and empty-position row behave as specified.

### Lane E — GTEx expression and QTL summaries

Script: `scripts/load_gtex.py`.

Input: median TPM gzip, both QTL tars.

Output:

- `warehouse/gtex/median_tpm.parquet` long gene × tissue.
- `warehouse/gtex/gene_expression.parquet`, one row per gene: `ensg`, max tissue, max median TPM, number of tissues with median TPM `>= 1` and `>= 10`, and a simple specificity statistic. Specificity can be tau: `sum(1 - x_i / x_max) / (n - 1)` over tissues with `x_max > 0`. Document the formula in the QC note.
- `warehouse/gtex/egenes.parquet` and `warehouse/gtex/sgenes.parquet` with tissue, gene, lead variant id, parsed variant key, p-value, and slope or beta column using whatever the file actually names. Read the header. Do not assume the slope column name.
- `warehouse/gtex/eqtl_pairs.parquet` and `warehouse/gtex/sqtl_pairs.parquet` for significant pairs. eQTL pair members match `*.signif_variant_gene_pairs.txt.gz`. sQTL pair members match `*.sqtl_signifpairs.txt.gz`. If a pair file's header differs from the eGene header, follow the pair header. Build `variant_key` from `chr`, `variant_pos`, `ref`, and `alt` when those columns exist. Include tissue in every row.

QC: tissue count for median TPM (GTEx v8 has 54 tissues; if the file disagrees, report the file). eGene files found versus tissues. Variant ids that failed the `chr_pos_ref_alt_b38` parse. Gene ids that will not match GENCODE v49 because the version or the gene is gone. Do not drop those genes.

Defer `gene_tpm.gct.gz` (sample-level). Mention it as unused in the QC note.

### Lane F — MaveDB

Script: `scripts/load_mavedb.py`.

Input: `raw/mavedb/mavedb-dump.2026062418131.tar.gz`.

Output:

- `warehouse/mavedb/experiments.parquet` from `main.json` after reading `README.md`.
- `warehouse/mavedb/mapped_variants.parquet` with urn, target gene, protein change, genomic chrom/pos/ref/alt when the mapped JSON provides them, `variant_key` when they are GRCh38, score, and any author-normalized score the dump actually contains. Read a small mapped JSON before writing the schema. The plan does not invent MaveDB field names beyond what a first file shows.
- `warehouse/mavedb/gene_coverage.parquet`, one row per target gene: number of score sets, number of mapped variants, whether any genomic mapping exists.

QC: AppleDouble files skipped. Count of mapped versus unmapped experiments. Count of protein-coding genes covered. Build named in the mapped records, and a count of non-GRCh38 coordinates left unjoined.

### Lane G — crosswalk

Script: `scripts/build_crosswalk.py`.

This lane should start as soon as Lane A has written genes, and it should be rerun if it started too early. Inputs: interval genes, constraint transcripts, ClinVar gene ids.

Output: `warehouse/crosswalk/genes.parquet` with the columns in the contract section. A second file `warehouse/crosswalk/unmatched.parquet` lists ClinVar symbols and constraint gene ids that did not join, with counts.

QC: protein-coding GENCODE genes matched to a constraint row, unmatched rate, and duplicate symbols.

### Lane H — truth-set harness

Script: `scripts/truth_set.py`.

Run the Wave 0 truth set against whatever Wave 1 outputs exist. Write `warehouse/qc/truth_set_report.md`. Lanes that have not finished show as skipped, not as failures. Rerun at the end of Wave 1 so the report is complete.

## Wave 2 — integrator, after Wave 1 outputs exist

Three products. Chromosome annotation can be parallel. The mart and the ledger are single writers.

### Gene mart

Script: `scripts/build_gene_mart.py`.

Output: `warehouse/marts/gene_mart.parquet`, one row per GENCODE gene, and `warehouse/marts/gene_mart_protein_coding.parquet` as a view or a filtered copy.

Columns:

- Identity from intervals: `ensg`, `symbol`, `chrom`, `start`, `end`, `strand`, `biotype`, MANE transcript id.
- Constraint: LOEUF, LOEUF decile, pLI, `mis.z_score`, `syn.z_score`, `constraint_flags`. Null when the gene has no constraint row. Noncoding genes usually will not.
- Expression: max tissue, max median TPM, tissues at or above 1 TPM, tau.
- ClinVar: counts of alleles at each star bucket split by germline pathogenic / likely pathogenic, VUS, benign / likely benign, conflicting, and drug response. Use the deduplicated allele table.
- GWAS: number of clumped loci whose `mapped_genes` contain the symbol or whose lead variant falls inside the gene span, and the number of distinct trait URIs. Document which rule fired. Span overlap and mapped-gene name will disagree. Keep both counts: `gwas_loci_by_overlap` and `gwas_loci_by_mapped_gene`.
- MaveDB: number of experiments, number of mapped variants.
- QTL: number of tissues where the gene is an eGene, number where it is an sGene.

QC: row count equals the gene interval count. Protein-coding file contains the protein-coding biotype only. Truth-set genes carry the expected constraint and a ClinVar count.

### Tier A variants

Script: `scripts/build_tier_a.py` builds the key list. `scripts/annotate_gnomad.py` fills frequencies and scores.

Tier A is the union of:

- ClinVar rows with a VCF `variant_key`.
- GWAS associations that gained a `variant_key` through rsID match to ClinVar or through a successful gnomAD rsID lookup. A genome-wide rsID lookup against every sites VCF is expensive. Do it per chromosome for the GWAS rsIDs that lack alleles, using `bcftools` with the ID column, and cache hits. Unresolved GWAS rsIDs stay in `warehouse/variants/tier_a_unresolved.tsv` and stay out of the annotated table.
- MaveDB mapped variants with a GRCh38 key.
- Lead eQTL and lead sQTL variants from the eGene and sGene tables, not every significant pair.

Output:

- `warehouse/variants/tier_a_keys.tsv` sorted by chromosome and position.
- `warehouse/variants/tier_a_genome.parquet` with the extract fields listed in the snapshot section, plus `in_genomes` bool.
- `warehouse/variants/tier_a_exome.parquet` with the same fields prefixed or column-suffixed `_exome`, and only for chromosomes whose exome VCF and index are complete. Skip a chromosome that is still downloading.
- `warehouse/variants/tier_a.parquet` joined on `variant_key`.

A missing gnomAD record is `in_genomes = false`. That is expected for many rare ClinVar alleles. Count it. Do not drop the variant.

Annotate one chromosome per process. chr1 and chr2 are the large ones. chrY is the test chromosome. Run chrY first and check the truth set before launching the rest.

Consequence parsing: from `vep`, keep `consequence`, `impact`, `symbol`, `enst`, `hgvsc`, `hgvsp`, and a boolean `on_mane` when the transcript matches the gene's MANE transcript. If the VEP version in the header uses different subfield names, follow the header and document the mapping.

### Evidence ledger

Script: `scripts/build_ledger.py`.

Output: `warehouse/ledger/evidence.parquet`.

| Column | Meaning |
|---|---|
| `evidence_id` | Stable hash of the natural key |
| `subject_type` | `variant` or `gene` |
| `subject_id` | `variant_key` or bare `ensg` |
| `object_type` | `trait`, `phenotype`, `gene_expression`, `splicing`, `assay`, `constraint` |
| `object_id` | EFO URI, phenotype id, `ensg|tissue`, score-set urn, or metric name |
| `evidence_class` | `clinvar`, `gwas`, `eqtl`, `sqtl`, `mavedb`, `constraint` |
| `strength` | Stars, p-value, absolute score, or LOEUF. One column `strength_value` plus `strength_unit` |
| `source_row` | File and a source id (`VariationID`, study accession, urn) |
| `ancestry_note` | Free text when the source has it, else null |
| `n_note` | Sample size prose or allele number, else null |

Populate it from the deduplicated ClinVar alleles, clumped GWAS loci plus the lead association, lead eQTL and sQTL rows, MaveDB gene coverage and mapped-variant scores, and gene constraint. Tier B QTL pairs stay out of the first ledger. They remain in the pair Parquet files for a question that needs them.

### Dark-genome inventory

From the interval CDS and exon tables, the clumped GWAS loci, and the lead QTL variant positions, compute fractions of the GRCh38 primary assembly. Use chromosome lengths from the gnomAD VCF `contig` header rather than a hard-coded length table. Write `warehouse/reports/dark_genome.md` with bases in protein-coding CDS, bases in any gene span, bases within 1 kb of a Tier A variant, and bases in a clumped GWAS window. Overlaps are not exclusive. Also report the exclusive remainder that is in none of those sets. This note is an input to Wave 3 so later writers do not talk as if the annotated subset were the whole genome.

Wave 2 acceptance: gene mart row count matches genes. Tier A chrY annotation matches the truth set. Ledger has all six evidence classes. Dark-genome note exists. Free space and exome completion are remeasured and written at the top of `warehouse/qc/wave2.md`.

## Wave 3 — question tracks

Read-only against the marts. Each track writes `warehouse/reports/<track>.md` and the DuckDB SQL, or the script, that produced the numbers. No track rewrites a mart. If a track finds a broken join, it files the problem in `warehouse/qc/issues.md` and stops that claim. It does not quietly patch the mart.

1. **Rare-disease ranking.** Protein-coding genes ordered by LOEUF, count of germline pathogenic or likely pathogenic alleles at two or more stars, and expression breadth. Split out genes whose pathogenic calls are single-submitter or conflicting. Name the filters in the report.
2. **Common-disease gene hypotheses.** For each EFO trait with at least a stated number of clumped loci, list loci, catalog mapped genes, and whether that gene is an eGene in a tissue that matches the trait. Use an explicit tissue map for the obvious cases: blood traits with whole blood, brain traits with the brain tissues present in GTEx v8, liver enzymes with liver, heart traits with the heart tissues. Traits with no tissue mapping stay in an unmapped table. Call this a shared-gene hypothesis. Do not call it colocalization.
3. **Coding versus noncoding.** Using GWAS `CONTEXT` and `INTERGENIC`, report the share of clumped loci that are missense, intronic, UTR, or intergenic, overall and for the largest trait groups. This is the part of the catalog that ClinVar and constraint barely describe.
4. **Pleiotropy.** Genes and Tier A variants linked to many EFO traits. Build a trait–trait table of shared genes from clumped loci. A shared gene inside a huge window is weaker than a shared lead variant. Show both.
5. **Conflict mining.** Germline ClinVar pathogenic or likely pathogenic alleles with high genome AF. State the AF threshold in the report and show more than one threshold. GWAS risk alleles that are the common allele. High missense-z genes whose MaveDB scores are mostly mild, with the number of assayed variants visible. Conflicting ClinVar labels on variants with high REVEL or high SpliceAI. These rows audit the warehouse as much as they describe biology.
6. **Constraint versus function.** LOEUF decile against the fraction of ClinVar alleles that are pathogenic at two or more stars, and against MaveDB effect sizes where an assay exists. State how few protein-coding genes have a MaveDB assay. A genome-wide claim from MaveDB alone is out of bounds.
7. **Splicing.** ClinVar splice-region variants that have SpliceAI or Pangolin scores, overlapped with sQTLs for the same gene. This is one of the few mechanism-shaped results in these files. Keep clinical, predictor, and sQTL evidence in separate columns.
8. **Tissue programs.** For trait groups, summarize median expression of mapped genes across GTEx tissues. Separate broadly expressed genes from tissue-private genes using the tau already on the gene mart.
9. **Ancestry audit.** For Tier A GWAS variants with enough `AN` in each group to be meaningful, report `AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, and `AF_amr`, plus missingness. Most catalog discoveries were made in European-ancestry samples. A frequency in the discovery group is not a frequency for everyone. Report coverage. Do not rank groups. Do not infer disease risk for a population from a lead SNP.
10. **ClinVar against frequency.** Pathogenic and likely pathogenic alleles should mostly be rare. Plot or bin clinical significance against genome AF and against exome AF. Investigate the common exceptions. Expected causes include allele mismatches, founder alleles, risk factors, drug-response rows mislabeled by a filter, and bad joins. Write what the exceptions actually were after checking a sample of them back to the source files.
11. **Noncoding genes.** A separate report on lncRNA and other non-protein-coding genes: expression, GWAS overlap, nearby protein-coding genes. Keep these conclusions out of the protein-coding rankings.
12. **Drug response.** The ClinVar drug-response set joined to constraint and expression. Say how small it is. A pharmacogenetics catalog needs CPIC and PharmGKB, which are not on disk. Do not imply this subset is that catalog.
13. **Somatic versus germline.** Counts and top genes for somatic origin and oncogenicity, held apart from the germline rare-disease ranking.
14. **Study redundancy.** How far 1.2 million association rows shrink after trait-and-distance clumping, and which traits are dominated by repeated studies of the same locus.
15. **Predictor disagreement.** On Tier A coding variants, compare CADD, REVEL, SpliceAI, a MaveDB score when one exists, and ClinVar stars. The interesting set is disagreement. Agreement is the set to trust first. Report missingness of each predictor. gnomAD does not score every site with every tool.
16. **Genes with no evidence.** Low or no expression, no constraint call, no ClinVar submission, no GWAS locus by either rule, no MaveDB assay. Publish the count and a table. This is the boundary of what the warehouse knows.

Each report starts with the question, the filters, the denominator, the result, and the caveat that changes how to read it. A number without a denominator is incomplete.

## Phase 2 — only after Wave 2 and a new free-space measurement

Do not download these in the first swarm.

Small, and the right next annotation layer for intergenic GWAS hits:

- ENCODE candidate cis-regulatory elements, GRCh38.
- Human Phenotype Ontology, to structure ClinVar `PhenotypeIDS`.
- Gene Ontology or Reactome, for gene-set summaries of the rankings.
- A MANE summary file only if the GTF and the constraint file disagree or lack tags.

Large, and only for a trait list the owner chooses:

- GWAS Catalog full summary statistics for those traits. Many studies have none. Some are GRCh37 and need a liftover policy. Some are access-controlled and are out of scope.
- A linkage-disequilibrium panel, such as high-coverage 1000 Genomes, if the owner wants real polygenic scores or clumping by LD instead of by distance. This is hundreds of gigabytes. It waits on free space.
- gnomAD structural-variant VCFs, if the owner wants CNVs. ClinVar's `nsv/esv (dbVar)` column is only an identifier until those calls exist locally. Another large download.

OMIM full text is licensed. Use the phenotype identifiers already in ClinVar.

AlphaMissense and extra pathogenicity scores are optional. CADD, REVEL, SpliceAI, Pangolin, phyloP, SIFT, and PolyPhen are already in the sites VCFs for the sites gnomAD annotated.

## Swarm briefing

Copy this preamble onto every agent. Then give the agent exactly one lane.

```text
You are working on D:\variant-workbench. Read plan.md and follow it.
GRCh38 only. Import scripts/coords.py for chromosome names, variant keys, and Ensembl ids.
Write only the output paths in your lane. Do not insert gnomAD sites into SQLite.
Do not add exome allele counts to genome allele counts.
Do not download new datasets.
Do not kill or restart scripts/fetch_sources.py.
Run the truth-set cases your lane is responsible for.
Write a QC note with row counts, null rates, and join misses.
If a source header disagrees with plan.md, follow the file and record the difference in the QC note.
```

Order:

1. One agent runs Step 0 and Wave 0 to acceptance. The download is already complete. Confirm that. Do not download anything.
2. Launch Lanes A through H together. Lane G reruns after A, B, and C exist. Lane H reruns after the others finish.
3. One integrator builds the gene mart, Tier A keys, the ledger, and the dark-genome note. Chromosome annotation of Tier A can be one agent per chromosome after chrY passes.
4. Launch the sixteen question tracks together, read-only.

The first product that makes new questions possible is the gene mart plus a passing truth-set report. The Tier A table is the second product. Question tracks interpret those two. They do not replace them.

## Definition of done

- Wave 0 note, manifest, `coords.py`, and a passing tool check.
- Wave 1 Parquet datasets and QC notes for intervals, constraint, ClinVar, GWAS, GTEx, MaveDB, and the crosswalk.
- Gene mart, Tier A annotation for every complete gnomAD chromosome, evidence ledger, dark-genome note.
- Sixteen reports, each with a denominator and a caveat.
- `warehouse/qc/issues.md` either empty or explicit about what remains broken.
- No second copy of the gnomAD VCFs, and free space still above the fetch reserve.
