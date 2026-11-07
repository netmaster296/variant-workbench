# ClinVar

Source: `D:\variant-workbench\raw\clinvar\variant_summary.txt.gz` (457,862,838 bytes).
Alleles: `D:\variant-workbench\warehouse\clinvar\alleles.parquet` (391,492,660 bytes).
Dedup: `D:\variant-workbench\warehouse\clinvar\allele_dedup.parquet` (176,020,865 bytes).
GRCh38 only. GRCh37 is excluded. NCBI36 and `na` assemblies are excluded too.

## Assemblies

| Assembly | Rows | Disposition |
| --- | ---: | ---: |
| `GRCh37` | 4,692,274 | excluded |
| `GRCh38` | 4,639,270 | loaded |
| `na` | 9,552 | excluded |
| `NCBI36` | 4,771 | excluded |

Bad-length records skipped: 0.
Allele scan and write: 142.3 seconds.

## Alleles

- GRCh38 rows: 4,639,270.
- Distinct VariationID: 4,636,853.
- Distinct AlleleID: 4,636,853.
- Rows with a variant_key: 4,618,645.
- Rows with no VCF allele (null variant_key): 20,625, across 20,608 AlleleIDs.
- Distinct variant_key: 4,618,645.
- PhenotypeList max length: 1,478 (source 1,478). Not truncated.
- PhenotypeIDS max length: 4,741 (source 4,741). Kept.
- PhenotypeList longer than 500 characters: 701.

Distinct VariationID is only modestly smaller than the GRCh38 row count: 2,417 extra rows. Phenotypes are already packed into one row with `|`. The extra rows are the same AlleleID placed on more than one chromosome, not one row per phenotype.

AlleleIDs with more than one GRCh38 row: 2,417. Of those, 2,415 differ by chromosome, 0 differ by clinical significance, 0 differ by PhenotypeList, 2,399 have more than one variant_key, and 17 have no variant_key.
variant_keys shared by more than one source row: 0.
variant_keys shared by more than one AlleleID: 0.

Chromosome values that do not canonicalize:

| chrom_raw | Rows |
| --- | ---: |
| `Un` | 3 |

The variant key uses PositionVCF, ReferenceAlleleVCF, and AlternateAlleleVCF. A VCF allele of `na` is missing, not a sequence. Start, Stop, ReferenceAllele, and AlternateAllele stay in the source columns. `rsid` is `norm_rsid` of `RS# (dbSNP)`.

## Review stars

Match is case-insensitive on the full ReviewStatus phrase. 4 practice guideline; 3 expert panel; 2 multiple submitters and no conflicts; 1 single submitter or conflicting classifications; 0 no assertion, no interpretation, or a blank status. Any other phrase, including `-`, is null.

| review_stars | Rows |
| --- | ---: |
| `0` | 111,292 |
| `1` | 3,486,444 |
| `2` | 766,045 |
| `3` | 22,853 |
| `4` | 54 |
| `null` | 252,582 |

Unmatched phrases (review_stars null), by raw ReviewStatus:

| ReviewStatus | Rows |
| --- | ---: |
| `-` | 244,654 |
| `no classification provided` | 7,064 |
| `no classification for the single variant` | 674 |
| `no classifications from unflagged records` | 190 |

Matched phrases:

| review_stars | ReviewStatus | Rows |
| --- | ---: | ---: |
| `0` | no assertion criteria provided | 111,292 |
| `1` | criteria provided, single submitter | 3,303,135 |
| `1` | criteria provided, conflicting classifications | 183,309 |
| `2` | criteria provided, multiple submitters, no conflicts | 766,045 |
| `3` | reviewed by expert panel | 22,853 |
| `4` | practice guideline | 54 |

Dedup uses the max star in the group:

| review_stars | Rows |
| --- | ---: |
| `0` | 111,281 |
| `1` | 3,486,441 |
| `2` | 766,045 |
| `3` | 22,853 |
| `4` | 54 |
| `null` | 252,579 |

## Dedup

- Rows: 4,639,253.
- Keyed by variant_key: 4,618,645.
- Keyed by AlleleID because there is no VCF allele: 20,608. dedup_key is `allele:` plus AlleleID.
- Source rows covered by n_rows: 4,639,270.
- germline_pathogenic true: 338,587.

germline_pathogenic is true only when the aggregated OriginSimple class is germline (not somatic, not germline/somatic, not unknown), Oncogenicity has no call, and ClinicalSignificance has a pathogenic or likely pathogenic piece without a conflicting piece. Low-penetrance pathogenic labels count. Drug response alone does not. Oncogenicity and SomaticClinicalImpact stay in their own columns and are not copied into ClinicalSignificance.

Pathogenic pieces that were not flagged: 37,763 because origin was not germline-only, 836 because of an oncogenicity call, 836 hit both, 0 were germline origin blocked only by oncogenicity.

## Clinical significance

| ClinicalSignificance | Rows |
| --- | ---: |
| `Uncertain significance` | 2,397,869 |
| `Likely benign` | 1,144,428 |
| `-` | 244,654 |
| `Pathogenic` | 197,557 |
| `Benign` | 195,919 |
| `Conflicting classifications of pathogenicity` | 183,629 |
| `Likely pathogenic` | 132,766 |
| `Benign/Likely benign` | 83,396 |
| `Pathogenic/Likely pathogenic` | 45,648 |
| `not provided` | 7,065 |
| `drug response` | 1,889 |
| `other` | 1,519 |
| `no classification for the single variant` | 674 |
| `risk factor` | 352 |
| `association` | 334 |
| `no classifications from unflagged records` | 190 |
| `Uncertain significance/Uncertain risk allele` | 139 |
| `conflicting data from submitters` | 132 |
| `Affects` | 130 |
| `Likely risk allele` | 87 |
| `Pathogenic, low penetrance` | 80 |
| `VUS-mid` | 73 |
| `Pathogenic; other` | 68 |
| `VUS-high` | 65 |
| `Likely pathogenic, low penetrance` | 59 |
| `Pathogenic; drug response` | 40 |
| `protective` | 40 |
| `Uncertain risk allele` | 36 |
| `Uncertain significance/VUS-high` | 29 |
| `VUS-low` | 28 |
| `Pathogenic; risk factor` | 22 |
| `Conflicting classifications of pathogenicity; other` | 19 |
| `Benign; risk factor` | 18 |
| `Likely pathogenic; drug response` | 17 |
| `Conflicting classifications of pathogenicity; risk factor` | 16 |
| `Likely pathogenic/Likely risk allele` | 14 |
| `Uncertain significance/VUS-mid` | 14 |
| `Benign; other` | 12 |
| `Benign/Likely benign; other` | 11 |
| `Likely benign; other` | 11 |
| `Uncertain significance; other` | 11 |
| `confers sensitivity` | 10 |
| `Benign; association` | 9 |
| `Uncertain significance/VUS-low` | 9 |
| `Likely pathogenic/Pathogenic, low penetrance` | 8 |
| `Pathogenic/Likely pathogenic/Pathogenic, low penetrance` | 8 |
| `Pathogenic/Likely pathogenic; other` | 8 |
| `Pathogenic/Likely risk allele` | 8 |
| `Pathogenic; Affects` | 8 |
| `Conflicting classifications of pathogenicity; association` | 7 |
| `Uncertain significance; drug response` | 7 |
| `Uncertain significance; risk factor` | 7 |
| `Uncertain significance; association` | 6 |
| `Likely pathogenic; association` | 5 |
| `Pathogenic/Likely pathogenic/Likely risk allele` | 5 |
| `Pathogenic/Likely pathogenic; risk factor` | 5 |
| `Pathogenic/Pathogenic, low penetrance` | 5 |
| `Benign; Affects` | 4 |
| `Benign; protective` | 4 |
| `Likely benign; drug response; other` | 4 |
| `Likely pathogenic/Likely pathogenic, low penetrance` | 4 |
| `Likely pathogenic; risk factor` | 4 |
| `Pathogenic; association` | 4 |
| `association not found` | 4 |
| `protective; risk factor` | 4 |
| `Benign/Likely benign; association` | 3 |
| `Benign/Likely benign; drug response` | 3 |
| `Benign/Likely benign; drug response; other` | 3 |
| `Benign; drug response` | 3 |
| `Conflicting classifications of pathogenicity; drug response` | 3 |
| `Likely benign; drug response` | 3 |
| `Pathogenic/Likely pathogenic, low penetrance` | 3 |
| `Conflicting classifications of pathogenicity; association; risk factor` | 2 |
| `Conflicting classifications of pathogenicity; other; risk factor` | 2 |
| `Conflicting classifications of pathogenicity; protective` | 2 |
| `Established risk allele` | 2 |
| `Likely benign; risk factor` | 2 |
| `Pathogenic/Likely pathogenic; drug response` | 2 |
| `Pathogenic; protective` | 2 |
| `drug response; risk factor` | 2 |
| `Benign/Likely benign; other; risk factor` | 1 |
| `Benign/Likely benign; risk factor` | 1 |
| `Benign; Affects; association; other` | 1 |
| `Benign; confers sensitivity` | 1 |
| `Conflicting classifications of pathogenicity; Affects` | 1 |
| `Conflicting classifications of pathogenicity; drug response; other` | 1 |
| `Established risk allele; association` | 1 |
| `Likely benign; association` | 1 |
| `Likely pathogenic; Affects` | 1 |
| `Pathogenic/Likely pathogenic/Established risk allele` | 1 |
| `Pathogenic/Likely pathogenic/Likely pathogenic, low penetrance` | 1 |
| `Pathogenic/Likely pathogenic/Pathogenic, low penetrance/Established risk allele; risk factor` | 1 |
| `Pathogenic/Likely pathogenic/Pathogenic, low penetrance; other` | 1 |
| `Pathogenic/Likely pathogenic; association` | 1 |
| `Pathogenic/Pathogenic, low penetrance; other` | 1 |
| `Pathogenic/Pathogenic, low penetrance; risk factor` | 1 |
| `Pathogenic; confers sensitivity` | 1 |
| `Uncertain risk allele; risk factor` | 1 |
| `Uncertain significance; Affects` | 1 |
| `Uncertain significance; Pathogenic/Likely pathogenic` | 1 |
| `association; drug response; risk factor` | 1 |
| `association; protective` | 1 |
| `association; risk factor` | 1 |
| `confers sensitivity; other` | 1 |
| `drug response; other` | 1 |
| `other; risk factor` | 1 |

## Non-ACGTN VCF alleles

| ref | alt | Rows |
| --- | ---: | ---: |
| `C` | R | 10 |
| `G` | Y | 9 |
| `A` | Y | 3 |
| `T` | Y | 2 |
| `C` | D | 1 |
| `T` | R | 1 |
| `C` | CAAAACAARATCGAGGACCTCACAGAGCTGGAAGGCTGGGAAGAGCTGGAAGATAATGAGAAGGAACAGATAACCCAGCACATTGCAGGT | 1 |
| `AG` | YT | 1 |
| `G` | H | 1 |

## Truth set

Checks read `warehouse/qc/truth_set.tsv` and the Parquet files.

| Case | Result | Detail |
| --- | ---: | ---: |
| `coding_snp` | PASS | 1 allele row(s); dedup [{'review_stars': 2, 'germline_pathogenic': True, 'germline_somatic': 'germline', 'ClinicalSignificance': 'Pathogenic', 'Oncogenicity': '-'}]; FOXRED1 rs267606829 Pathogenic stars 2, germline_pathogenic true |
| `cftr_indel` | PASS | [{'GeneSymbol': 'CFTR', 'ClinicalSignificance': 'Pathogenic', 'ReviewStatus': 'practice guideline', 'review_stars': 4, 'ref': 'ATCT', 'alt': 'A', 'variant_key': 'chr7:117559590:ATCT:A', 'pheno_len': 783}]; Pathogenic, stars 4, review contains practice guideline, alleles ATCT/A, PhenotypeList longer than 500 |
| `mismatched_ref` | PASS | [{'ref': 'G', 'alt': 'A', 'variant_key': 'chr15:84799209:G:A', 'ReferenceAllele': 'na', 'Name': 'NM_014630.3(ZNF592):c.3136G>A (p.Gly1046Arg)'}]; key uses VCF G/A, not display na and not the HGVS name |
| `chrom_crosswalk` | PASS | example [{'chrom_raw': '1', 'chrom': 'chr1', 'variant_key': 'chr1:243344276:AG:A', 'ref': 'AG', 'alt': 'A', 'GeneSymbol': 'SDCCAG8'}]; chrom_raw 1 -> chr1 on 417899 rows; mismatches 0 |
| `somatic_or_oncogenicity` | PASS | allele [{'GeneSymbol': 'NF1', 'OriginSimple': 'germline/somatic', 'ClinicalSignificance': 'Pathogenic', 'Oncogenicity': 'Likely oncogenic', 'SomaticClinicalImpact': '-', 'review_stars': 2, 'ref': 'C', 'alt': 'T'}]; dedup [{'germline_pathogenic': False, 'germline_somatic': 'germline/somatic', 'ClinicalSignificance': 'Pathogenic', 'Oncogenicity': 'Likely oncogenic', 'SomaticClinicalImpact': '-', 'OriginSimple': 'germline/somatic', 'review_stars': 2}]; not counted as germline pathogenic |
| `phenotype_text` | PASS | PhenotypeList max 1478 source 1478, stored intact True; PhenotypeIDS max 4741 source 4741, stored intact True; rows longer than 500 characters: 701 |

## Integrity

- allele_rows_match: PASS
- assembly_only_grch38: PASS
- bad_len: PASS
- no_sep_collision: PASS
- empty_allele: PASS
- dedup_covers_rows: PASS
- dedup_keys_unique: PASS
- dedup_keyed_match: PASS
- dedup_bare_match: PASS
- dedup_written: PASS
- flag_count: PASS

All integrity checks and truth-set cases passed.
