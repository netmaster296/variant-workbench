# Interval QC

GENCODE v49 intervals from `raw/gencode/gencode.v49.annotation.gtf.gz` (93,374,019 bytes gzip). Gene rows: 78,691 (baseline 78,691, delta +0). Protein-coding (`biotype = protein_coding`): 20,097 (baseline 20,097, delta +0). Transcript rows with tag `MANE_Select`: 19,276. Genes with `mane_select` true: 19,276.

Chromosomes that failed canonicalization: none. Every kept GTF chromosome passed `canon_chrom`, and every `chrom` value in the five Parquet files is one of chr1–chr22, chrX, chrY, or chrM.

## Source and rules

Header lines from the GTF:

```text
##description: evidence-based annotation of the human genome (GRCh38), version 49 (Ensembl 115)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2025-07-08
```

GTF `start` and `end` are stored as 1-based inclusive integers. They were not converted to 0-based half-open coordinates. `ensg` and `enst` are version-stripped. `ensg_version` and `enst_version` are the original ids, including the version suffix when the GTF has one. Rows whose version field has no `.` suffix: 0.

Feature types kept: gene, transcript, exon, CDS, UTR, stop_codon. `stop_codon` rows are in `cds.parquet` with `feature = stop_codon`. CDS rows have `feature = CDS`. Exon and UTR files do not repeat those features. `start_codon` and `Selenocysteine` were not kept.

Gene feature rows in this file do not carry `MANE_Select` (0 gene rows had that tag). `genes.mane_select` is true only when a transcript of that versioned gene id has a tag whose value is exactly `MANE_Select`, or when the gene row itself has that tag. `transcripts.mane_select` is true only when that transcript row's own tags contain `MANE_Select`. `MANE_Plus_Clinical` is not treated as MANE Select. Tags were not invented.

Counts below are DuckDB queries against the Parquet files just written, except the read-time drop table.

## Row counts

| file | rows | bytes |
| --- | --- | --- |
| genes.parquet | 78,691 | 2,364,509 |
| transcripts.parquet | 507,365 | 9,305,175 |
| exons.parquet | 3,673,949 | 14,772,240 |
| cds.parquet | 2,498,156 | 8,146,448 |
| utrs.parquet | 770,688 | 4,635,827 |

CDS file by feature:

| feature | rows |
| --- | --- |
| CDS | 2,283,012 |
| stop_codon | 215,144 |

Parse versus file: gene 78,691, transcript 507,365, exon 3,673,949, CDS 2,283,012, stop_codon 215,144, UTR 770,688. These match the Parquet row counts (78,691, 507,365, 3,673,949, 2,498,156, 770,688; CDS plus stop_codon = 2,498,156).

## Null rates

Null `ensg` across all five files: 0.

### genes

| column | nulls | null_rate |
| --- | --- | --- |
| ensg | 0 | 0 |
| ensg_version | 0 | 0 |
| symbol | 0 | 0 |
| chrom | 0 | 0 |
| start | 0 | 0 |
| end | 0 | 0 |
| strand | 0 | 0 |
| biotype | 0 | 0 |
| mane_select | 0 | 0 |
| hgnc_id | 37,114 | 0.471642 |

### transcripts

| column | nulls | null_rate |
| --- | --- | --- |
| ensg | 0 | 0 |
| ensg_version | 0 | 0 |
| symbol | 0 | 0 |
| chrom | 0 | 0 |
| start | 0 | 0 |
| end | 0 | 0 |
| strand | 0 | 0 |
| biotype | 0 | 0 |
| mane_select | 0 | 0 |
| hgnc_id | 115,824 | 0.228285 |
| enst | 0 | 0 |
| enst_version | 0 | 0 |
| transcript_type | 0 | 0 |
| tags | 0 | 0 |

### exons

| column | nulls | null_rate |
| --- | --- | --- |
| ensg | 0 | 0 |
| enst | 0 | 0 |
| chrom | 0 | 0 |
| start | 0 | 0 |
| end | 0 | 0 |
| strand | 0 | 0 |
| exon_number | 0 | 0 |

### cds

| column | nulls | null_rate |
| --- | --- | --- |
| ensg | 0 | 0 |
| enst | 0 | 0 |
| chrom | 0 | 0 |
| start | 0 | 0 |
| end | 0 | 0 |
| strand | 0 | 0 |
| exon_number | 0 | 0 |
| feature | 0 | 0 |

### utrs

| column | nulls | null_rate |
| --- | --- | --- |
| ensg | 0 | 0 |
| enst | 0 | 0 |
| chrom | 0 | 0 |
| start | 0 | 0 |
| end | 0 | 0 |
| strand | 0 | 0 |
| exon_number | 0 | 0 |

Transcript rows with an empty tag list: 62,029. `hgnc_id` values that do not start with `HGNC:`: 0. Non-integer `exon_number` values stored as null: 0.

Every gene strand is `+` or `-`.

Rows with `start < 1` or `end < start`, all files: 0.

## Gene count versus 78691

`genes.parquet` has 78,691 rows. The SQLite gene-table baseline in the plan is 78,691. Delta: +0.

Duplicate versioned gene ids: 0. Duplicate bare `ensg`: 0. Duplicate versioned transcript ids: 0. Duplicate bare `enst`: 0.

Gene biotypes:

| biotype | genes |
| --- | --- |
| lncRNA | 34,880 |
| protein_coding | 20,097 |
| processed_pseudogene | 9,487 |
| misc_RNA | 2,207 |
| unprocessed_pseudogene | 1,949 |
| snRNA | 1,901 |
| miRNA | 1,879 |
| transcribed_unprocessed_pseudogene | 1,587 |
| transcribed_processed_pseudogene | 1,149 |
| TEC | 1,019 |
| snoRNA | 942 |
| rRNA_pseudogene | 497 |
| transcribed_unitary_pseudogene | 201 |
| IG_V_pseudogene | 187 |
| IG_V_gene | 146 |
| TR_V_gene | 107 |
| unitary_pseudogene | 89 |
| TR_J_gene | 79 |
| scaRNA | 49 |
| rRNA | 47 |
| IG_D_gene | 37 |
| TR_V_pseudogene | 33 |
| Mt_tRNA | 22 |
| artifact | 19 |
| IG_J_gene | 18 |
| IG_C_gene | 14 |
| IG_C_pseudogene | 9 |
| ribozyme | 8 |
| TR_C_gene | 6 |
| TR_D_gene | 5 |
| sRNA | 5 |
| TR_J_pseudogene | 4 |
| vault_RNA | 4 |
| IG_J_pseudogene | 3 |
| Mt_rRNA | 2 |
| translated_processed_pseudogene | 2 |
| IG_pseudogene | 1 |

## Protein-coding versus 20097

Genes with `biotype = protein_coding`: 20,097. Baseline: 20,097. Delta: +0.

Biotypes whose name contains `protein`:

| biotype | genes |
| --- | --- |
| protein_coding | 20,097 |

The baseline comparison uses the exact biotype `protein_coding` only. Other protein-like biotypes are not folded into that count.

## MANE tag count

The GTF does contain `MANE_Select`. Transcript rows with that exact tag: 19,276. Distinct versioned gene ids on those transcripts: 19,276. `genes.mane_select` true: 19,276. Gene-feature rows with the tag: 0.

Transcripts where `mane_select` disagrees with `list_contains(tags, 'MANE_Select')`: 0. Gene flags with no MANE transcript: 0. MANE transcripts whose gene flag is false: 0. MANE transcripts whose `transcript_type` is not `protein_coding`: 50.

Tag occurrences on transcript rows (unnested; one transcript contributes one if the tag is present once):

| tag | occurrences |
| --- | --- |
| Ensembl_canonical | 78,691 |
| MANE_Plus_Clinical | 66 |
| MANE_Select | 19,276 |
| appris_alternative_1 | 12,898 |
| appris_alternative_2 | 12,368 |
| appris_principal_1 | 40,915 |
| appris_principal_2 | 5,790 |
| appris_principal_3 | 16,481 |
| appris_principal_4 | 6,772 |
| appris_principal_5 | 228 |

Predicate used: `list_contains(tags, 'MANE_Select')`.

## Chromosomes that failed canonicalization

None. No kept row was dropped for chromosome name, and no Parquet `chrom` is outside the canonical set.

No chromosome name was rewritten. Raw GTF names and stored `chrom` values are the same.

Stored chromosome counts:

| chrom | genes | transcripts | exons | cds | utrs |
| --- | --- | --- | --- | --- | --- |
| chr1 | 7,090 | 48,010 | 357,629 | 248,494 | 75,349 |
| chr2 | 5,686 | 38,255 | 288,688 | 196,037 | 55,930 |
| chr3 | 4,158 | 31,373 | 246,383 | 169,571 | 54,489 |
| chr4 | 3,577 | 21,415 | 152,218 | 100,948 | 29,291 |
| chr5 | 3,946 | 24,027 | 162,413 | 100,322 | 31,053 |
| chr6 | 4,230 | 25,648 | 172,707 | 111,187 | 35,896 |
| chr7 | 3,957 | 26,951 | 195,011 | 124,871 | 36,990 |
| chr8 | 3,265 | 20,905 | 134,738 | 82,936 | 28,698 |
| chr9 | 3,188 | 21,361 | 156,322 | 100,820 | 28,967 |
| chr10 | 3,268 | 20,905 | 163,994 | 117,621 | 33,272 |
| chr11 | 4,166 | 27,524 | 205,031 | 149,264 | 47,944 |
| chr12 | 3,860 | 25,993 | 189,985 | 132,250 | 40,541 |
| chr13 | 2,001 | 10,011 | 67,995 | 42,488 | 11,102 |
| chr14 | 2,797 | 17,051 | 124,402 | 86,602 | 28,479 |
| chr15 | 2,818 | 17,938 | 127,127 | 78,485 | 24,460 |
| chr16 | 3,188 | 20,639 | 149,064 | 104,404 | 32,627 |
| chr17 | 3,780 | 25,996 | 195,809 | 145,300 | 43,232 |
| chr18 | 1,653 | 9,482 | 61,418 | 38,751 | 11,822 |
| chr19 | 3,545 | 25,328 | 175,475 | 134,130 | 46,998 |
| chr20 | 1,970 | 12,664 | 88,176 | 60,481 | 18,929 |
| chr21 | 1,134 | 6,806 | 43,732 | 24,561 | 7,217 |
| chr22 | 1,747 | 11,615 | 90,130 | 64,263 | 18,226 |
| chrX | 2,958 | 15,421 | 114,237 | 79,827 | 27,993 |
| chrY | 672 | 2,010 | 11,228 | 4,525 | 1,178 |
| chrM | 37 | 37 | 37 | 18 | 5 |

## Read-time drops

Comment lines: 5. Data rows: 7,750,154. Skipped features are not drops.

| feature | rows_skipped |
| --- | --- |
| Selenocysteine | 130 |
| start_codon | 221,175 |

Rows dropped for chromosome, id, or coordinate: 0.

## Join checks inside this lane

Transcript rows whose `ensg_version` is absent from genes: 0. Genes with no transcript: 0.

| file | rows_without_matching_transcript | rows_whose_ensg_disagrees_with_transcript |
| --- | --- | --- |
| exons | 0 | 0 |
| cds | 0 | 0 |
| utrs | 0 | 0 |

Symbols shared by more than one bare `ensg`: 483. Symbol is not the primary key. Later joins should use `ensg`.

| symbol | genes | chroms |
| --- | --- | --- |
| Y_RNA | 756 | 24 |
| Metazoa_SRP | 170 | 23 |
| U3 | 50 | 18 |
| U6 | 33 | 17 |
| SNORA70 | 27 | 16 |
| U8 | 22 | 12 |
| U2 | 19 | 8 |
| 5S_rRNA | 9 | 8 |
| SNORA72 | 8 | 4 |
| 7SK | 7 | 5 |

## Truth-set interval checks

These are the truth-set cases this lane can test from coordinates and tags. Constraint, ClinVar, GWAS, and gnomAD comparisons belong to other lanes. Positions are 1-based and match a feature when `start <= pos <= end`.

| case | result | detail |
| --- | --- | --- |
| coding_snp | PASS | FOXRED1 genes 1; CDS rows of that symbol covering chr11:126275389: 24; CDS rows of any gene covering it: 24 |
| cftr_indel | PASS | CFTR genes 1; CDS rows of that symbol covering chr7:117559590: 22; CDS rows of any gene covering it: 22 |
| gwas_intergenic | PASS | chr22:45604058 CDS 0, stop_codon 0, exon 0, UTR 0, gene span 0 ((none)) |
| chrom_crosswalk | PASS | expected ENSG00000054282.17 chr1:243255349-243500095 protein_coding; pos 243344276 inside span: True; found ENSG00000054282.17 chr1:243255349-243500095 protein_coding mane_select=True |
| mane_select | PASS | expected A1BG ENST00000263100 tag MANE_Select gene ENSG00000121410; gene.mane_select=True; found ENSG00000121410.13 ENST00000263100.8 A1BG chr19 mane_select=True tag=True. The constraint-file half of this case is not this lane. |

