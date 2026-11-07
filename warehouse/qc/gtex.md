# GTEx lane E

Loaded 2026-11-03. Assembly GRCh38. Variant keys and bare ENSG ids come from `scripts/coords.py`.

## Free space

- Before: 676,485,353,472 bytes.
- After: 669,209,763,840 bytes.
- Floor that would have blocked pair output: 128,849,018,880 bytes (120 GiB).
- Pair text estimate, used only for that floor check: 32,048,621,909 bytes.
- Estimate method: gzip size of both pair-member sets times the larger measured uncompressed/gzip ratio (11.09), times 1.5. The ratio came from `GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.txt.gz`.
- Decision: wrote both pair files. Free space 676,485,353,472 minus estimate 32,048,621,909 stays above 128,849,018,880.

## Inputs

- Median TPM: `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`, 6,952,331 bytes. GCT `#1.2`, dimension line `56200	54`.
- eQTL tar: `GTEx_Analysis_v8_eQTL.tar`, 1,562,828,800 bytes, 98 file members. Ignored members: none.
- sQTL tar: `GTEx_Analysis_v8_sQTL.tar`, 569,169,920 bytes, 98 file members. Ignored members: none.
- Members were read from the tars and parsed in memory. The tars were not extracted to disk.
- Cells that are empty, `NA`, `NaN`, or `.` are null. A `.` rs id therefore has a null raw value and a null `rsid`.

## Sample-level TPM

`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` (1,625,807,147 bytes) was not opened. It is the sample-level gene TPM matrix. This lane uses the median TPM file instead.

## Median TPM

- Tissue columns: 54. The file's dimension line says 54. GTEx v8 is described as 54 tissues. Match: True.
- Genes: 56,200. Long rows written: 3,034,800 (56,200 x 54).
- `ensg` null (Name was not an ENSG id): 44 / 56,200 (0.0783%). All 44 are Names ending in `_PAR_Y`. `bare_ensg` does not accept that suffix, so `ensg` is null and `ensg_version` keeps the Name. No row was dropped.
- Duplicate versioned ids: 0. Extra rows sharing a bare ENSG: 0.
- Negative TPM cells: 0. Non-finite cells: 0.
- Tissue labels are the GCT column names, in file order. QTL tissues use a different label (`Adipose_Subcutaneous`, not `Adipose - Subcutaneous`). This lane does not rename either set.

GCT tissues:

Adipose - Subcutaneous, Adipose - Visceral (Omentum), Adrenal Gland, Artery - Aorta, Artery - Coronary, Artery - Tibial, Bladder, Brain - Amygdala, Brain - Anterior cingulate cortex (BA24), Brain - Caudate (basal ganglia), Brain - Cerebellar Hemisphere, Brain - Cerebellum, Brain - Cortex, Brain - Frontal Cortex (BA9), Brain - Hippocampus, Brain - Hypothalamus, Brain - Nucleus accumbens (basal ganglia), Brain - Putamen (basal ganglia), Brain - Spinal cord (cervical c-1), Brain - Substantia nigra, Breast - Mammary Tissue, Cells - Cultured fibroblasts, Cells - EBV-transformed lymphocytes, Cervix - Ectocervix, Cervix - Endocervix, Colon - Sigmoid, Colon - Transverse, Esophagus - Gastroesophageal Junction, Esophagus - Mucosa, Esophagus - Muscularis, Fallopian Tube, Heart - Atrial Appendage, Heart - Left Ventricle, Kidney - Cortex, Kidney - Medulla, Liver, Lung, Minor Salivary Gland, Muscle - Skeletal, Nerve - Tibial, Ovary, Pancreas, Pituitary, Prostate, Skin - Not Sun Exposed (Suprapubic), Skin - Sun Exposed (Lower leg), Small Intestine - Terminal Ileum, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina, Whole Blood

## Gene expression

One row per GCT gene. `ensg` is `bare_ensg(ensg_version)`.

Tau, over the n tissue columns in the median file, only when x_max > 0:

```text
tau = sum_i(1 - x_i / x_max) / (n - 1)
```

x_i is median TPM in tissue i. x_max is the maximum of those values. n is the tissue count in the file, including tissues where the gene is 0. No pseudocount is added. Tau is null when x_max is 0, because the ratio is undefined. Under this definition a single-tissue gene has tau 1 and a gene with the same median in every tissue has tau 0.

- n = 54.
- Tau null: 11,981. Genes whose maximum is 0: 11,981. For those, `max_tissue` is null and `max_median_tpm` is 0.
- Finite tau range: 0.2898120532992479 to 1.0.
- Genes with median TPM >= 1 in at least one tissue: 31,910. Genes with median TPM >= 10 in at least one tissue: 19,201.
- Ties for the maximum, counted only when x_max > 0: 1. `max_tissue` is the leftmost GCT column among the tied tissues.

## eGenes

- Files found: 49. Expected: 49. Match: True.
- Source header is the same in every `*.egenes.txt.gz` member. The effect column is `slope`. `slope_se` is its standard error. There is no effect column named `beta`. `beta_shape1` and `beta_shape2` are beta-distribution parameters for the permutation approximation, not effect sizes.
- `pval_nominal` is the lead variant's nominal p-value. `pval_perm` is the permutation p-value. `pval_beta` is the beta-approximated p-value. `qval` is the gene-level FDR.
- The file has one row per tested gene per tissue, including genes that are not significant. `qval_le_0_05` is true when `qval` <= 0.05. That is the eGene call. Rows with a larger q-value are kept.
- Rows: 1,207,976. qval <= 0.05: 475,829. qval null: 0. Values in (0.05, 0.05+1e-8]: 0.
- Distinct versioned genes: 39,832. Distinct versioned genes with qval <= 0.05 in at least one tissue: 34,548.
- Duplicate gene ids within a tissue (extra rows): 0.
- `variant_key` null: 0 / 1,207,976 (0.0000%). Loader tally of null keys: 0.
- `ensg` null: 0 / 1,207,976 (0.0000%).
- `rsid` null: 2,379 / 1,207,976 (0.1969%). 2,297 of those source cells were `.` and were read as null. The other 82 are comma-separated rs id pairs (49 distinct strings, for example `rs4909,rs12012747`). `norm_rsid` does not choose one id from a pair. The raw column keeps the text and `rsid` is null.
- `slope` null: 0 / 1,207,976 (0.0000%). `pval_nominal` null: 0 / 1,207,976 (0.0000%). `qval` null: 0 / 1,207,976 (0.0000%).
- `log2_aFC` null: 0 / 1,207,976 (0.0000%). `log2_aFC_lower` null: 0 / 1,207,976 (0.0000%). `log2_aFC_upper` null: 43,520 / 1,207,976 (3.6027%). Empty allelic-fold-change cells in the source are null here, not zeros.
- Split-column key missing and id key used instead: 0.
- Non-null disagreements between the split-column key and the b38 id key: 0. Examples: none.
- Variant ids that failed the b38 parse: 0. Examples: none.
- The stored `variant_key` prefers `chr`, `variant_pos`, `ref`, and `alt`. The id is the check, and the fallback only when those columns do not make a key.

| Tissue | Rows | qval <= 0.05 |
|---|---:|---:|
| Adipose_Subcutaneous | 24,665 | 15,607 |
| Adipose_Visceral_Omentum | 24,724 | 12,482 |
| Adrenal_Gland | 23,820 | 8,123 |
| Artery_Aorta | 23,959 | 12,493 |
| Artery_Coronary | 24,529 | 6,296 |
| Artery_Tibial | 23,304 | 15,008 |
| Brain_Amygdala | 24,069 | 3,726 |
| Brain_Anterior_cingulate_cortex_BA24 | 24,342 | 5,640 |
| Brain_Caudate_basal_ganglia | 24,718 | 8,362 |
| Brain_Cerebellar_Hemisphere | 25,144 | 10,027 |
| Brain_Cerebellum | 25,461 | 11,240 |
| Brain_Cortex | 24,849 | 9,082 |
| Brain_Frontal_Cortex_BA9 | 24,676 | 7,335 |
| Brain_Hippocampus | 24,420 | 5,517 |
| Brain_Hypothalamus | 25,096 | 5,499 |
| Brain_Nucleus_accumbens_basal_ganglia | 24,890 | 8,198 |
| Brain_Putamen_basal_ganglia | 23,804 | 6,902 |
| Brain_Spinal_cord_cervical_c-1 | 24,669 | 4,483 |
| Brain_Substantia_nigra | 24,044 | 3,301 |
| Breast_Mammary_Tissue | 25,849 | 10,872 |
| Cells_Cultured_fibroblasts | 22,050 | 15,219 |
| Cells_EBV-transformed_lymphocytes | 22,759 | 4,934 |
| Colon_Sigmoid | 24,483 | 10,550 |
| Colon_Transverse | 25,379 | 11,686 |
| Esophagus_Gastroesophageal_Junction | 24,168 | 10,534 |
| Esophagus_Mucosa | 23,949 | 14,675 |
| Esophagus_Muscularis | 23,871 | 14,110 |
| Heart_Atrial_Appendage | 23,194 | 10,991 |
| Heart_Left_Ventricle | 21,353 | 9,642 |
| Kidney_Cortex | 24,807 | 1,260 |
| Liver | 22,262 | 5,734 |
| Lung | 26,095 | 14,113 |
| Minor_Salivary_Gland | 25,579 | 4,836 |
| Muscle_Skeletal | 21,031 | 13,532 |
| Nerve_Tibial | 25,873 | 17,285 |
| Ovary | 25,325 | 5,542 |
| Pancreas | 22,615 | 9,660 |
| Pituitary | 26,854 | 9,146 |
| Prostate | 26,529 | 7,356 |
| Skin_Not_Sun_Exposed_Suprapubic | 25,279 | 15,483 |
| Skin_Sun_Exposed_Lower_leg | 25,196 | 16,967 |
| Small_Intestine_Terminal_Ileum | 26,182 | 6,681 |
| Spleen | 25,479 | 10,783 |
| Stomach | 24,290 | 8,771 |
| Testis | 35,007 | 18,795 |
| Thyroid | 26,054 | 17,684 |
| Uterus | 25,188 | 3,577 |
| Vagina | 25,778 | 3,730 |
| Whole_Blood | 20,315 | 12,360 |

## sGenes

- Files found: 49. Expected: 49. Match: True.
- Source header is the same in every `*.sgenes.txt.gz` member. The effect column is again `slope`, not `beta`.
- Grain is one tested splicing phenotype (intron cluster) per tissue, not one row per gene. `phenotype_id` looks like `chrom:start:end:cluster:ENSG`. `gene_id` is a source column.
- A gene is an sGene in a tissue when any phenotype row for that gene and tissue has qval <= 0.05. `qval_le_0_05` is that phenotype-level flag. It is not collapsed to one row per gene.
- Rows: 618,932. qval <= 0.05: 142,320. qval null: 0. Values in (0.05, 0.05+1e-8]: 0.
- Distinct versioned genes: 23,879. Distinct versioned genes with at least one phenotype at qval <= 0.05: 16,895.
- Duplicate phenotype ids within a tissue (extra rows): 0.
- Phenotype tail differing from the `gene_id` column: 0. Examples: none.
- `variant_key` null: 0 / 618,932 (0.0000%). Loader tally: 0.
- `ensg` null: 0 / 618,932 (0.0000%).
- `rsid` null: 1,508 / 618,932 (0.2436%). 1,470 of those source cells were `.` and were read as null. The other 38 are comma-separated rs id pairs (25 distinct strings, for example `rs199730337,rs864309627`). `norm_rsid` does not choose one id from a pair. The raw column keeps the text and `rsid` is null.
- `slope` null: 0 / 618,932 (0.0000%). `pval_nominal` null: 0 / 618,932 (0.0000%). `qval` null: 0 / 618,932 (0.0000%).
- Column-key fallback to the id: 0. Non-null key disagreements: 0. Examples: none.
- Variant ids that failed the b38 parse: 0. Examples: none.
- Minor-allele columns in this file are `ma_samples` and `ma_count`. The eGene file calls them `minor_allele_samples` and `minor_allele_count`. Both names are kept as published.

| Tissue | Rows | qval <= 0.05 |
|---|---:|---:|
| Adipose_Subcutaneous | 12,392 | 5,113 |
| Adipose_Visceral_Omentum | 12,691 | 4,210 |
| Adrenal_Gland | 12,413 | 2,369 |
| Artery_Aorta | 12,220 | 3,740 |
| Artery_Coronary | 12,834 | 2,140 |
| Artery_Tibial | 11,847 | 4,791 |
| Brain_Amygdala | 12,583 | 892 |
| Brain_Anterior_cingulate_cortex_BA24 | 12,975 | 1,238 |
| Brain_Caudate_basal_ganglia | 13,062 | 1,809 |
| Brain_Cerebellar_Hemisphere | 12,647 | 2,397 |
| Brain_Cerebellum | 12,548 | 2,786 |
| Brain_Cortex | 13,013 | 2,048 |
| Brain_Frontal_Cortex_BA9 | 13,079 | 1,684 |
| Brain_Hippocampus | 12,671 | 1,185 |
| Brain_Hypothalamus | 13,293 | 1,414 |
| Brain_Nucleus_accumbens_basal_ganglia | 13,177 | 1,881 |
| Brain_Putamen_basal_ganglia | 12,427 | 1,350 |
| Brain_Spinal_cord_cervical_c-1 | 12,772 | 1,144 |
| Brain_Substantia_nigra | 12,474 | 814 |
| Breast_Mammary_Tissue | 13,303 | 4,124 |
| Cells_Cultured_fibroblasts | 11,338 | 4,586 |
| Cells_EBV-transformed_lymphocytes | 11,529 | 2,360 |
| Colon_Sigmoid | 12,550 | 3,269 |
| Colon_Transverse | 13,235 | 3,459 |
| Esophagus_Gastroesophageal_Junction | 12,281 | 3,286 |
| Esophagus_Mucosa | 12,014 | 3,986 |
| Esophagus_Muscularis | 12,063 | 4,088 |
| Heart_Atrial_Appendage | 12,173 | 3,055 |
| Heart_Left_Ventricle | 11,164 | 2,357 |
| Kidney_Cortex | 12,188 | 547 |
| Liver | 11,049 | 1,485 |
| Lung | 13,294 | 4,774 |
| Minor_Salivary_Gland | 13,346 | 1,648 |
| Muscle_Skeletal | 10,488 | 4,056 |
| Nerve_Tibial | 12,754 | 5,297 |
| Ovary | 12,713 | 1,998 |
| Pancreas | 11,279 | 2,250 |
| Pituitary | 13,691 | 2,901 |
| Prostate | 13,510 | 2,463 |
| Skin_Not_Sun_Exposed_Suprapubic | 12,629 | 4,652 |
| Skin_Sun_Exposed_Lower_leg | 12,584 | 5,134 |
| Small_Intestine_Terminal_Ileum | 13,721 | 2,084 |
| Spleen | 12,767 | 2,837 |
| Stomach | 12,646 | 2,638 |
| Testis | 19,027 | 8,626 |
| Thyroid | 12,801 | 5,358 |
| Uterus | 12,800 | 1,524 |
| Vagina | 13,292 | 1,460 |
| Whole_Blood | 9,585 | 3,013 |

## eQTL pairs

- Files found: 49. Expected: 49. Match: True.
- Member suffix `*.signif_variant_gene_pairs.txt.gz`. These rows are the significant pairs GTEx already filtered. This loader does not filter them again.
- The pair header has no `chr`, `variant_pos`, `ref`, or `alt`. `variant_key` is parsed from `variant_id` (`chr_pos_ref_alt_b38`).
- Rows: 71,478,479. Distinct versioned genes: 34,548.
- `variant_key` null: 0 / 71,478,479 (0.0000%). Loader tally: 0.
- `ensg` null: 0 / 71,478,479 (0.0000%). `slope` null: 0 / 71,478,479 (0.0000%). `pval_nominal` null: 0 / 71,478,479 (0.0000%). `pval_beta` null: 0 / 71,478,479 (0.0000%).
- Variant ids that failed the b38 parse: 0. Examples: none.
- Versioned pair genes absent from the eGene file: 0.

| Tissue | Rows |
|---|---:|
| Adipose_Subcutaneous | 2,930,370 |
| Adipose_Visceral_Omentum | 2,031,091 |
| Adrenal_Gland | 959,376 |
| Artery_Aorta | 2,038,599 |
| Artery_Coronary | 721,637 |
| Artery_Tibial | 2,958,254 |
| Brain_Amygdala | 347,359 |
| Brain_Anterior_cingulate_cortex_BA24 | 531,518 |
| Brain_Caudate_basal_ganglia | 906,716 |
| Brain_Cerebellar_Hemisphere | 1,183,391 |
| Brain_Cerebellum | 1,505,827 |
| Brain_Cortex | 1,062,868 |
| Brain_Frontal_Cortex_BA9 | 770,794 |
| Brain_Hippocampus | 541,493 |
| Brain_Hypothalamus | 568,169 |
| Brain_Nucleus_accumbens_basal_ganglia | 908,866 |
| Brain_Putamen_basal_ganglia | 695,347 |
| Brain_Spinal_cord_cervical_c-1 | 422,663 |
| Brain_Substantia_nigra | 283,774 |
| Breast_Mammary_Tissue | 1,616,203 |
| Cells_Cultured_fibroblasts | 2,890,486 |
| Cells_EBV-transformed_lymphocytes | 458,161 |
| Colon_Sigmoid | 1,474,957 |
| Colon_Transverse | 1,672,709 |
| Esophagus_Gastroesophageal_Junction | 1,556,679 |
| Esophagus_Mucosa | 2,616,828 |
| Esophagus_Muscularis | 2,529,273 |
| Heart_Atrial_Appendage | 1,620,487 |
| Heart_Left_Ventricle | 1,409,424 |
| Kidney_Cortex | 111,169 |
| Liver | 629,559 |
| Lung | 2,392,204 |
| Minor_Salivary_Gland | 422,971 |
| Muscle_Skeletal | 2,601,575 |
| Nerve_Tibial | 3,460,273 |
| Ovary | 570,916 |
| Pancreas | 1,340,448 |
| Pituitary | 1,226,883 |
| Prostate | 823,690 |
| Skin_Not_Sun_Exposed_Suprapubic | 2,697,337 |
| Skin_Sun_Exposed_Lower_leg | 3,244,342 |
| Small_Intestine_Terminal_Ileum | 644,572 |
| Spleen | 1,262,024 |
| Stomach | 1,156,158 |
| Testis | 2,887,027 |
| Thyroid | 3,712,240 |
| Uterus | 328,520 |
| Vagina | 338,599 |
| Whole_Blood | 2,414,653 |

## sQTL pairs

- Files found: 49. Expected: 49. Match: True.
- Member suffix `*.sqtl_signifpairs.txt.gz`. This is not the eQTL pair name. The header uses `phenotype_id` instead of `gene_id` and has no split allele columns.
- `gene_id` in the Parquet is the last colon-separated field of `phenotype_id` when that field is an ENSG id. `ensg` is the bare id. Rows whose phenotype does not end in an ENSG are kept, with those two columns null.
- Phenotype ids that did not yield an ENSG: 0. Examples: none.
- Rows: 24,445,206. Distinct versioned genes: 16,895.
- `variant_key` null: 0 / 24,445,206 (0.0000%). Loader tally: 0.
- `ensg` null: 0 / 24,445,206 (0.0000%). `gene_id` null: 0 / 24,445,206 (0.0000%). `slope` null: 0 / 24,445,206 (0.0000%). `pval_nominal` null: 0 / 24,445,206 (0.0000%). `pval_beta` null: 0 / 24,445,206 (0.0000%).
- Variant ids that failed the b38 parse: 0. Examples: none.
- Versioned pair genes absent from the sGene file: 0.

| Tissue | Rows |
|---|---:|
| Adipose_Subcutaneous | 1,019,648 |
| Adipose_Visceral_Omentum | 751,211 |
| Adrenal_Gland | 372,731 |
| Artery_Aorta | 667,895 |
| Artery_Coronary | 308,815 |
| Artery_Tibial | 955,738 |
| Brain_Amygdala | 112,753 |
| Brain_Anterior_cingulate_cortex_BA24 | 163,953 |
| Brain_Caudate_basal_ganglia | 261,977 |
| Brain_Cerebellar_Hemisphere | 344,489 |
| Brain_Cerebellum | 401,291 |
| Brain_Cortex | 287,562 |
| Brain_Frontal_Cortex_BA9 | 226,576 |
| Brain_Hippocampus | 172,761 |
| Brain_Hypothalamus | 202,412 |
| Brain_Nucleus_accumbens_basal_ganglia | 269,372 |
| Brain_Putamen_basal_ganglia | 186,304 |
| Brain_Spinal_cord_cervical_c-1 | 145,545 |
| Brain_Substantia_nigra | 104,747 |
| Breast_Mammary_Tissue | 716,593 |
| Cells_Cultured_fibroblasts | 889,990 |
| Cells_EBV-transformed_lymphocytes | 288,529 |
| Colon_Sigmoid | 527,835 |
| Colon_Transverse | 605,299 |
| Esophagus_Gastroesophageal_Junction | 548,916 |
| Esophagus_Mucosa | 759,227 |
| Esophagus_Muscularis | 786,214 |
| Heart_Atrial_Appendage | 526,642 |
| Heart_Left_Ventricle | 421,367 |
| Kidney_Cortex | 45,566 |
| Liver | 221,828 |
| Lung | 886,768 |
| Minor_Salivary_Gland | 189,537 |
| Muscle_Skeletal | 894,894 |
| Nerve_Tibial | 1,043,635 |
| Ovary | 260,517 |
| Pancreas | 355,195 |
| Pituitary | 482,318 |
| Prostate | 368,898 |
| Skin_Not_Sun_Exposed_Suprapubic | 902,558 |
| Skin_Sun_Exposed_Lower_leg | 1,037,195 |
| Small_Intestine_Terminal_Ileum | 265,416 |
| Spleen | 395,315 |
| Stomach | 423,319 |
| Testis | 1,537,363 |
| Thyroid | 1,119,675 |
| Uterus | 182,070 |
| Vagina | 185,912 |
| Whole_Blood | 620,835 |

## Variant id parse, all QTL tables

Accepted form: five underscore fields, `chrom`, `pos`, `ref`, `alt`, `b38`, with pos an integer and alleles A/C/G/T/N. The key is then `variant_key` from `scripts/coords.py`, so the chromosome is canonical and the alleles are uppercase.

- eGenes: 0.
- sGenes: 0.
- eQTL pairs: 0.
- sQTL pairs: 0.
- Total failed ids: 0.

## GENCODE v49

Gene features in `gencode.v49.annotation.gtf.gz`: 78,691 versioned ids, 78,691 bare ENSG ids.

Join rule is the bare ENSG. A version-only difference still joins. A bare id that is absent from GENCODE v49 does not. No GTEx row was dropped for either reason.

| Source | Versioned ids | Bare ids | Version-only difference | Bare id absent | Not an ENSG |
|---|---:|---:|---:|---:|---:|
| median TPM | 56,200 | 56,156 | 27,485 | 1,530 | 44 |
| eGenes | 39,832 | 39,832 | 26,468 | 1,289 | 0 |
| sGenes | 23,879 | 23,879 | 21,467 | 752 | 0 |
| eQTL pairs | 34,548 | 34,548 | 24,779 | 1,086 | 0 |
| sQTL pairs | 16,895 | 16,895 | 15,650 | 414 | 0 |

Examples of bare ENSG ids absent from GENCODE v49 (these will not join):

- median TPM: ENSG00000002079, ENSG00000069712, ENSG00000112096, ENSG00000116883, ENSG00000116957, ENSG00000130489, ENSG00000130723, ENSG00000131484 (+1,522 more).
- eGenes: ENSG00000002079, ENSG00000069712, ENSG00000112096, ENSG00000116883, ENSG00000116957, ENSG00000130489, ENSG00000130723, ENSG00000131484 (+1,281 more).
- sGenes: ENSG00000002079, ENSG00000069712, ENSG00000112096, ENSG00000116883, ENSG00000116957, ENSG00000130723, ENSG00000133808, ENSG00000139656 (+744 more).
- eQTL pairs: ENSG00000069712, ENSG00000112096, ENSG00000116883, ENSG00000116957, ENSG00000130489, ENSG00000130723, ENSG00000131484, ENSG00000133808 (+1,078 more).
- sQTL pairs: ENSG00000069712, ENSG00000112096, ENSG00000116883, ENSG00000116957, ENSG00000130723, ENSG00000133808, ENSG00000148362, ENSG00000170165 (+406 more).

Examples of versioned ids whose bare id is in GENCODE v49 under a different version (these still join on the bare id):

- median TPM: ENSG00000000003.14, ENSG00000000005.5, ENSG00000000419.12, ENSG00000000457.13, ENSG00000000460.16, ENSG00000000938.12, ENSG00000000971.15, ENSG00000001036.13 (+27,477 more).
- eGenes: ENSG00000000003.14, ENSG00000000005.5, ENSG00000000419.12, ENSG00000000457.13, ENSG00000000460.16, ENSG00000000938.12, ENSG00000000971.15, ENSG00000001036.13 (+26,460 more).
- sGenes: ENSG00000000419.12, ENSG00000000457.13, ENSG00000000460.16, ENSG00000000938.12, ENSG00000000971.15, ENSG00000001036.13, ENSG00000001084.10, ENSG00000001167.14 (+21,459 more).
- eQTL pairs: ENSG00000000003.14, ENSG00000000005.5, ENSG00000000419.12, ENSG00000000457.13, ENSG00000000460.16, ENSG00000000938.12, ENSG00000000971.15, ENSG00000001036.13 (+24,771 more).
- sQTL pairs: ENSG00000000419.12, ENSG00000000457.13, ENSG00000000460.16, ENSG00000000938.12, ENSG00000000971.15, ENSG00000001036.13, ENSG00000001084.10, ENSG00000001167.14 (+15,642 more).

## Truth set `gtex_egene`

- Source locator: GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt.gz line 2.
- Tissue `Adipose_Subcutaneous`, gene `ENSG00000227232.5`, variant id `chr1_64764_C_T_b38`.
- Expected key `chr1:64764:C:T` from both the truth-set columns and `variant_key`.
- Seen while reading the eGene member: True. First data row of that member: True. Load index: 0.
- Key from the split columns at load: chr1:64764:C:T. Key from the variant id at load: chr1:64764:C:T. Stored key: chr1:64764:C:T.
- Re-read from `egenes.parquet`: 1 row(s). all checks passed.
- Result: PASS.

## Output files

| File | Rows | Bytes | Row groups |
|---|---:|---:|---:|
| `median_tpm.parquet` | 3,034,800 | 12,305,410 | 7 |
| `gene_expression.parquet` | 56,200 | 1,464,827 | 1 |
| `egenes.parquet` | 1,207,976 | 212,890,014 | 49 |
| `sgenes.parquet` | 618,932 | 111,185,290 | 49 |
| `eqtl_pairs.parquet` | 71,478,479 | 1,991,325,385 | 168 |
| `sqtl_pairs.parquet` | 24,445,206 | 672,752,482 | 75 |

Columns:

- `median_tpm.parquet`: ensg_version, ensg, tissue, median_tpm.
- `gene_expression.parquet`: ensg, ensg_version, max_tissue, max_median_tpm, n_tissues_tpm_ge_1, n_tissues_tpm_ge_10, tau, n_tissues.
- `egenes.parquet`: tissue, ensg, variant_key, rsid, qval_le_0_05, gene_id, gene_name, gene_chr, gene_start, gene_end, strand, num_var, beta_shape1, beta_shape2, true_df, pval_true_df, variant_id, tss_distance, chr, variant_pos, ref, alt, num_alt_per_site, rs_id_dbSNP151_GRCh38p7, minor_allele_samples, minor_allele_count, maf, ref_factor, pval_nominal, slope, slope_se, pval_perm, pval_beta, qval, pval_nominal_threshold, log2_aFC, log2_aFC_lower, log2_aFC_upper.
- `sgenes.parquet`: tissue, ensg, variant_key, rsid, qval_le_0_05, phenotype_id, gene_name, gene_chr, gene_start, gene_end, strand, num_var, beta_shape1, beta_shape2, true_df, pval_true_df, variant_id, tss_distance, chr, variant_pos, ref, alt, num_alt_per_site, rs_id_dbSNP151_GRCh38p7, ma_samples, ma_count, maf, ref_factor, pval_nominal, slope, slope_se, pval_perm, pval_beta, gene_id, group_size, qval, pval_nominal_threshold.
- `eqtl_pairs.parquet`: tissue, variant_id, gene_id, ensg, variant_key, tss_distance, ma_samples, ma_count, maf, pval_nominal, slope, slope_se, pval_nominal_threshold, min_pval_nominal, pval_beta.
- `sqtl_pairs.parquet`: tissue, variant_id, phenotype_id, gene_id, ensg, variant_key, tss_distance, ma_samples, ma_count, maf, pval_nominal, slope, slope_se, pval_nominal_threshold, min_pval_nominal, pval_beta.

Elapsed: 454.7 seconds. Python 3.13.11. pyarrow 26.0.1. bcftools was not used.
