# Tissue programs

Question: for the trait groups that can be tied to a GTEx tissue, where are the catalog-mapped genes expressed, and which of those genes are broadly expressed versus tissue-private by tau?

The tau values are `tau` on `warehouse/marts/gene_mart.parquet`. They were not recomputed. The formula, from `warehouse/qc/gtex.md`, is `sum(1 - x_i / x_max) / (n - 1)` over the 54 median-TPM tissues, and it is null when `x_max` is 0. Expression locations use `warehouse/gtex/median_tpm.parquet`. Those tissue names look like `Whole Blood`. QTL tables, whose names look like `Whole_Blood`, were not read and were not treated as the same labels.

The counts come from `warehouse/reports/sql/tissue_programs.sql`.

## Filters

Trait groups are name rules on `mapped_trait`, applied to traits that have at least one row in `warehouse/gwas/loci.parquet`. A trait is kept only when it matches exactly one rule. The exact patterns are the `group_rule` values in the SQL. In short:

- Blood cell: counts and indices of platelets, red cells, white cells, or hemoglobin. Hemoglobin A1c and glycated hemoglobin are excluded. So are `amount of ...` protein and metabolite measurements that merely mention a blood cell.
- Brain: brain, brainstem, cerebral, intracerebral, Alzheimer, Parkinson, schizophrenia, hippocampus, cerebellum or cerebellar, amygdala, hypothalamus, substantia nigra. `cerebellin` does not match. `precerebral` does not match.
- Liver: liver, hepatic, ALT, AST, bilirubin, or gamma-glutamyl transferase. Alkaline phosphatase is excluded because it is also a bone enzyme. Protein measurements are excluded.
- Heart: heart, cardiac, atrial, ventricular, myocardial, cardiomyopathy. Drug-use measurements and protein measurements are excluded.
- Lung: lung, pulmonary, asthma, COPD, FEV, FVC, forced expiratory volume, forced vital capacity. `fever` does not match `FEV`.
- Kidney: kidney, renal, glomerular, serum creatinine, creatinine clearance, albumin-to-creatinine ratio, glomerular filtration rate.
- Adipose: adipose or body fat. Body-mass index is not in the rule.
- Thyroid: thyroid, including hypothyroidism, plus thyroxine, triiodothyronine, thyrotropin, or Hashimoto. Parathyroid is excluded.
- Skin: skin, psoriasis, dermatitis, eczema, or cutaneous. Uveal and choroidal melanoma are excluded.
- Pancreas: pancreatitis, pancreatic carcinoma, and other pancreas words, excluding `amount of` protein measurements.
- Breast: breast, excluding breastfeeding.
- Prostate: prostate carcinoma, prostate cancer, and prostate-specific antigen. Other prostate protein measurements are excluded.
- Ovary: ovary or ovarian.
- Gut: colon, colonic, colorectal, Crohn, ulcerative colitis, or inflammatory bowel. Colony-stimulating factor does not match.
- Artery: artery, arterial, or atherosclerosis.
- Skeletal muscle: skeletal muscle or thigh muscle volume. Fat-infiltration traits are excluded.

A trait whose name matches two rules is ambiguous and is in neither gene set. Six traits, 33 locus rows:

| mapped trait | locus rows | rules |
| --- | ---: | --- |
| cerebral atherosclerosis | 14 | artery, brain |
| pulmonary arterial hypertension | 11 | artery, lung |
| chronic pulmonary heart disease | 2 | heart, lung |
| pulmonary artery enlargement | 2 | artery, lung |
| idiopathic pulmonary arterial hypertension | 2 | artery, lung |
| cerebral artery occlusion | 2 | artery, brain |

Lipids, bone density, blood pressure, and depression were not given a tissue. GTEx v8 has no bone column. Blood pressure is not whole blood; mean arterial pressure is the exception that the artery rule does take. Depression does not contain a brain word, so it stays out. Cholesterol and triglycerides are not assigned to liver.

The index tissues are median-TPM names only:

| group | index tissues |
| --- | --- |
| blood cell | Whole Blood |
| brain | the 13 columns that begin `Brain - `, including `Brain - Spinal cord (cervical c-1)` |
| liver | Liver |
| heart | Heart - Atrial Appendage, Heart - Left Ventricle |
| lung | Lung |
| kidney | Kidney - Cortex, Kidney - Medulla |
| adipose | Adipose - Subcutaneous, Adipose - Visceral (Omentum) |
| thyroid | Thyroid |
| skin | Skin - Not Sun Exposed (Suprapubic), Skin - Sun Exposed (Lower leg) |
| pancreas | Pancreas |
| breast | Breast - Mammary Tissue |
| prostate | Prostate |
| ovary | Ovary |
| gut | Colon - Sigmoid, Colon - Transverse, Small Intestine - Terminal Ileum |
| artery | Artery - Aorta, Artery - Coronary, Artery - Tibial |
| skeletal muscle | Muscle - Skeletal |

`Cells - EBV-transformed lymphocytes` is not used as whole blood. The two cultured-cell columns are not index tissues for any group.

Mapped genes are exact symbol matches from `mapped_genes.tokens` on those loci to `gene_mart.symbol`. This is the catalog token, not a gene-span overlap. A token is kept only when it matches exactly one GENCODE gene. Tokens of the form `GENE1 - GENE2` stay one token, as in the GWAS load, and do not match a symbol. Symbols that match more than one ENSG are left out. Case-only mismatches: 0. Other unmatched tokens: 0.

Tau classes, using the mart value:

- Broad: `tau < 0.5`. Under the formula above, that means the other 53 tissues average more than half the peak.
- Tissue-private: `tau >= 0.95`. The other 53 tissues average at most 5 percent of the peak.
- Intermediate: `0.5 <= tau < 0.95`. These genes are neither class.
- Tau undefined: a GTEx row whose maximum median TPM is 0. Not broad and not private.
- No GTEx row: `max_median_tpm` is null. Not broad and not private.

The minimum tau on the mart is 0.289812. A cutoff of `tau < 0.3` would label 1 of 78,691 genes. It is not used. On an on/off pattern, `tau = (54 - k) / 53`. Then `tau < 0.5` needs expression in at least 28 tissues, and `tau >= 0.95` is at least as concentrated as equal expression in 3 tissues (`k = 4` gives 0.943).

Detection in a tissue means median TPM `>= 1`. That is the same threshold the mart uses for `n_tissues_tpm_ge_1`. The location tables use protein-coding genes that have a GTEx row. The background for those tables is every protein-coding mart gene with a GTEx row: 19,209 genes. Excess is the group's detection fraction minus that background fraction. It is a difference of fractions, not a test. Absolute rank 1 is the highest detection fraction. Background rank 1 is testis.

No minimum locus count was applied. Skeletal muscle is reported in the counts and is not interpreted as a tissue ranking.

## Denominator

The mart has 78,691 GENCODE genes. 24,065 have no GTEx expression row. That is the same count as an anti-join from the mart to `median_tpm` on `ensg`. Of the 24,065, 888 are protein-coding. The other 23,177 are not. GENCODE v49 is newer than the GTEx v8 gene set; the missing rows are mostly noncoding genes GTEx did not quantify. Every mart gene that does have a row has 54 median-TPM rows, and mart `tau` matches `gene_expression.tau` on every shared `ensg`.

Locus rows in `loci.parquet`: 311,160, from 15,109 traits. Assigned to one group: 52,697 locus rows (569 traits). Ambiguous: 33 locus rows (6 traits). In neither: 258,430 locus rows.

Left out on purpose, and not in the ambiguous set:

| left out | traits | locus rows |
| --- | ---: | ---: |
| cholesterol or triglyceride | 94 | 14,290 |
| bone mineral, bone tissue density, heel bone, or bone fracture | 14 | 2,964 |
| blood pressure | 7 | 2,305 |
| depression | 8 | 1,703 |

These four rows are examples, not a partition of the 258,430.

Across the 16 groups, the union is 14,177 distinct genes. 242 of them have no GTEx row. 11,210 are protein-coding, and 131 of those have no GTEx row. 7,845 genes (6,460 protein-coding) fall in two or more groups, so the group rows below are not an exclusive partition. A gene is counted once inside a group no matter how many loci name it.

| group | traits | locus rows | mapped genes | no GTEx row | protein-coding | protein-coding with no GTEx row | protein-coding with a GTEx row |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| blood cell | 96 | 17,470 | 8,128 | 133 | 6,726 | 80 | 6,646 |
| brain | 72 | 7,364 | 4,795 | 77 | 3,849 | 46 | 3,803 |
| lung | 50 | 4,892 | 3,411 | 59 | 2,704 | 26 | 2,678 |
| heart | 84 | 3,570 | 2,134 | 29 | 1,741 | 18 | 1,723 |
| liver | 37 | 3,388 | 2,656 | 43 | 2,225 | 30 | 2,195 |
| kidney | 41 | 2,973 | 2,515 | 36 | 2,088 | 22 | 2,066 |
| gut | 19 | 2,555 | 1,800 | 34 | 1,440 | 17 | 1,423 |
| artery | 25 | 2,083 | 1,633 | 24 | 1,349 | 15 | 1,334 |
| thyroid | 21 | 1,694 | 1,216 | 23 | 966 | 14 | 952 |
| skin | 45 | 1,643 | 1,342 | 28 | 1,086 | 15 | 1,071 |
| breast | 25 | 1,549 | 1,207 | 19 | 951 | 7 | 944 |
| adipose | 13 | 1,503 | 1,169 | 28 | 954 | 17 | 937 |
| prostate | 4 | 1,089 | 1,003 | 21 | 791 | 11 | 780 |
| ovary | 16 | 603 | 390 | 5 | 294 | 3 | 291 |
| pancreas | 18 | 302 | 247 | 3 | 189 | 1 | 188 |
| skeletal muscle | 3 | 19 | 7 | 0 | 5 | 0 | 5 |

The no-GTEx share inside these mapped sets is much smaller than 24,065 / 78,691. The highest is 28 of 1,169 adipose genes. The blood-cell share is 133 of 8,128. The catalog tokens are mostly genes GTEx already quantified.

Unsplit `GENE1 - GENE2` tokens are the unmatched tokens. There are 4,504 distinct ones in the blood-cell set, out of 12,711 distinct tokens, and 9 of 16 in skeletal muscle. Ambiguous symbols, excluded, are 79 distinct symbols in the blood-cell set and 0 in skeletal muscle. After that exclusion, each remaining token is one gene.

LncRNAs are in the mapped-gene counts and not in the location table. The blood-cell set has 1,193 lncRNAs of 8,128 genes. The pancreas set has 51 of 247.

Some groups are mostly lab measurements. Locus rows whose names match the lab phrase:

| group | lab phrase | traits | locus rows | share of the group's locus rows |
| --- | --- | ---: | ---: | ---: |
| liver | aminotransferase, bilirubin, or gamma-glutamyl transferase | 9 | 2,577 | 2,577 / 3,388 |
| kidney | filtration or creatinine | 4 | 1,782 | 1,782 / 2,973 |
| thyroid | TSH, thyroxine, triiodothyronine, or thyrotropin | 5 | 496 | 496 / 1,694 |
| prostate | prostate-specific antigen | 1 | 330 | 330 / 1,089 |
| artery | arterial pressure | 1 | 361 | 361 / 2,083 |
| heart | troponin | 3 | 255 | 255 / 3,570 |

## Result

Most mapped genes are neither broad nor private. Among protein-coding genes in a group, the intermediate class is 76 to 83 percent, except skeletal muscle (3 of 5). Broad genes are 2 to 7 percent of protein-coding genes in those larger groups, and 0 of 5 in skeletal muscle. Private genes are 11 to 20 percent there, and 2 of 5 in skeletal muscle.

All protein-coding genes on the mart, for scale: 1,196 / 20,097 broad (tau < 0.5), 13,358 intermediate, 4,313 private (tau >= 0.95), 888 with no GTEx row, and 342 with a row but tau undefined. The mapped protein-coding sets are not enriched for private genes relative to that 4,313 / 20,097. Their broad shares are close to 1,196 / 20,097.

| group | broad | intermediate | private | tau undefined | no GTEx row | pc broad | pc intermediate | pc private |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| blood cell | 482 | 5,875 | 1,578 | 60 | 133 | 476 | 5,115 | 1,045 |
| brain | 210 | 3,550 | 932 | 26 | 77 | 206 | 3,085 | 509 |
| lung | 147 | 2,501 | 690 | 14 | 59 | 143 | 2,146 | 389 |
| heart | 84 | 1,604 | 405 | 12 | 29 | 84 | 1,403 | 234 |
| liver | 132 | 2,006 | 464 | 11 | 43 | 130 | 1,765 | 297 |
| kidney | 110 | 1,887 | 471 | 11 | 36 | 108 | 1,641 | 316 |
| gut | 64 | 1,306 | 387 | 9 | 34 | 64 | 1,121 | 237 |
| artery | 47 | 1,244 | 312 | 6 | 24 | 47 | 1,092 | 193 |
| thyroid | 65 | 878 | 246 | 4 | 23 | 64 | 745 | 143 |
| skin | 55 | 969 | 289 | 1 | 28 | 53 | 831 | 187 |
| breast | 38 | 924 | 219 | 7 | 19 | 37 | 790 | 116 |
| adipose | 54 | 899 | 183 | 5 | 28 | 53 | 777 | 107 |
| prostate | 41 | 752 | 187 | 2 | 21 | 41 | 624 | 115 |
| ovary | 10 | 289 | 83 | 3 | 5 | 10 | 239 | 41 |
| pancreas | 4 | 173 | 66 | 1 | 3 | 4 | 147 | 37 |
| skeletal muscle | 0 | 4 | 2 | 1 | 0 | 0 | 3 | 2 |

Within each group the five all-gene classes sum to the mapped-gene count, and the three protein-coding classes plus protein-coding no-row and protein-coding undefined sum to the protein-coding count. The undefined protein-coding counts are 10, 3, 0, 2, 3, 1, 1, 2, 0, 0, 1, 0, 0, 1, 0, and 0 in the same group order.

Broad protein-coding genes are the ubiquitous set. Their median number of tissues at or above 1 TPM is 54 in every group that has one. At least 44 of 47 broad protein-coding artery genes, the lowest such share, are at or above 1 TPM in all 54 tissues. In the blood-cell set it is 454 of 476. A broad gene is expressed in the index tissue because it is expressed in every tissue. That is not a tissue program.

Private protein-coding genes do have a real peak. The median of their maximum median TPM is between 44 and 90 in every group except skeletal muscle (2 genes, median peak 137). In the blood-cell set, 869 of 1,045 private protein-coding genes peak at 10 TPM or more. The peak is usually not the index tissue.

| group | private protein-coding | max tissue in the index set | max tissue testis | max tissue liver | max tissue a cultured cell |
| --- | ---: | ---: | ---: | ---: | ---: |
| blood cell | 1,045 | 111 | 338 | 112 | 111 |
| brain | 509 | 59 | 184 | 47 | 58 |
| lung | 389 | 11 | 124 | 34 | 44 |
| kidney | 316 | 24 | 77 | 44 | 22 |
| liver | 297 | 86 | 64 | 86 | 25 |
| gut | 237 | 11 | 75 | 23 | 29 |
| heart | 234 | 11 | 72 | 29 | 14 |
| artery | 193 | 2 | 54 | 42 | 20 |
| skin | 187 | 18 | 49 | 15 | 28 |
| thyroid | 143 | 7 | 36 | 19 | 29 |
| breast | 116 | 0 | 40 | 8 | 13 |
| prostate | 115 | 10 | 30 | 9 | 12 |
| adipose | 107 | 1 | 26 | 17 | 6 |
| ovary | 41 | 2 | 11 | 1 | 7 |
| pancreas | 37 | 7 | 8 | 6 | 1 |
| skeletal muscle | 2 | 1 | 0 | 0 | 0 |

Liver is the only large group whose private genes most often peak in the index tissue: 86 of 297 peak in liver, against 64 in testis. In the blood-cell private set the order is testis 338, liver 112, then whole blood 111. Breast has no private protein-coding gene whose maximum tissue is breast. Cultured fibroblasts and EBV-transformed lymphocytes are GTEx columns, not organs; together they account for 111 of 1,045 blood-cell private genes.

Where the genes are detected is a different cut from where the private genes peak. Among protein-coding genes with a GTEx row, the tissue with the highest share at median TPM >= 1 is testis in every group with at least 188 such genes. The same is true if lncRNAs and other biotypes are put back in. Testis is also rank 1 in the background: 15,923 / 19,209 protein-coding genes with a row (0.829). Whole blood is background rank 54, the lowest, at 9,898 / 19,209 (0.515).

Testis fractions in the mapped protein-coding sets are 0.852 to 0.914 (ovary 248 / 291, prostate 713 / 780). That is above 0.829, and it is still the top tissue. It is the GTEx detection baseline, not a claim that these traits are testis programs.

Index tissues, protein-coding genes with a GTEx row. Median TPM is the median across those genes of that tissue's median TPM. Skeletal muscle is omitted here and stated after the table.

| group | index tissue | >= 1 TPM | fraction | median TPM | absolute rank | background fraction | background rank | excess | excess rank |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| blood cell | Whole Blood | 4,235 / 6,646 | 0.637 | 2.93 | 54 | 0.515 | 54 | +0.122 | 2 |
| brain | Brain - Hypothalamus | 3,104 / 3,803 | 0.816 | 5.88 | 3 | 0.704 | 16 | +0.113 | 3 |
| brain | Brain - Cortex | 3,093 / 3,803 | 0.813 | 6.80 | 4 | 0.698 | 21 | +0.115 | 1 |
| brain | Brain - Frontal Cortex (BA9) | 3,087 / 3,803 | 0.812 | 7.28 | 5 | 0.698 | 20 | +0.113 | 2 |
| brain | Brain - Cerebellum | 3,031 / 3,803 | 0.797 | 10.27 | 12 | 0.703 | 17 | +0.094 | 13 |
| brain | Brain - Nucleus accumbens (basal ganglia) | 3,027 / 3,803 | 0.796 | 5.81 | 13 | 0.688 | 29 | +0.108 | 5 |
| brain | Brain - Caudate (basal ganglia) | 3,019 / 3,803 | 0.794 | 5.59 | 15 | 0.687 | 32 | +0.107 | 6 |
| brain | Brain - Anterior cingulate cortex (BA24) | 3,015 / 3,803 | 0.793 | 5.63 | 16 | 0.682 | 35 | +0.110 | 4 |
| brain | Brain - Spinal cord (cervical c-1) | 2,984 / 3,803 | 0.785 | 5.84 | 18 | 0.689 | 27 | +0.096 | 12 |
| brain | Brain - Cerebellar Hemisphere | 2,972 / 3,803 | 0.781 | 10.49 | 20 | 0.692 | 25 | +0.089 | 14 |
| brain | Brain - Hippocampus | 2,964 / 3,803 | 0.779 | 4.83 | 22 | 0.676 | 40 | +0.103 | 9 |
| brain | Brain - Amygdala | 2,956 / 3,803 | 0.777 | 4.68 | 24 | 0.672 | 42 | +0.105 | 7 |
| brain | Brain - Substantia nigra | 2,940 / 3,803 | 0.773 | 4.77 | 27 | 0.674 | 41 | +0.099 | 11 |
| brain | Brain - Putamen (basal ganglia) | 2,933 / 3,803 | 0.771 | 4.61 | 29 | 0.667 | 44 | +0.104 | 8 |
| lung | Lung | 2,152 / 2,678 | 0.804 | 10.61 | 4 | 0.728 | 3 | +0.076 | 26 |
| heart | Heart - Atrial Appendage | 1,341 / 1,723 | 0.778 | 5.34 | 35 | 0.647 | 47 | +0.131 | 1 |
| heart | Heart - Left Ventricle | 1,218 / 1,723 | 0.707 | 3.28 | 48 | 0.593 | 51 | +0.114 | 2 |
| liver | Liver | 1,619 / 2,195 | 0.738 | 4.36 | 48 | 0.601 | 50 | +0.136 | 1 |
| kidney | Kidney - Medulla | 1,676 / 2,066 | 0.811 | 7.33 | 10 | 0.696 | 23 | +0.115 | 1 |
| kidney | Kidney - Cortex | 1,609 / 2,066 | 0.779 | 5.55 | 30 | 0.670 | 43 | +0.109 | 2 |
| gut | Small Intestine - Terminal Ileum | 1,107 / 1,423 | 0.778 | 7.67 | 5 | 0.719 | 7 | +0.059 | 13 |
| gut | Colon - Transverse | 1,087 / 1,423 | 0.764 | 6.91 | 15 | 0.708 | 11 | +0.056 | 19 |
| gut | Colon - Sigmoid | 1,057 / 1,423 | 0.743 | 7.49 | 30 | 0.688 | 30 | +0.055 | 22 |
| artery | Artery - Coronary | 1,039 / 1,334 | 0.779 | 9.77 | 24 | 0.689 | 26 | +0.090 | 20 |
| artery | Artery - Aorta | 1,019 / 1,334 | 0.764 | 9.34 | 36 | 0.679 | 38 | +0.085 | 25 |
| artery | Artery - Tibial | 1,012 / 1,334 | 0.759 | 9.64 | 42 | 0.664 | 45 | +0.095 | 15 |
| thyroid | Thyroid | 807 / 952 | 0.848 | 13.60 | 5 | 0.718 | 8 | +0.129 | 8 |
| skin | Skin - Sun Exposed (Lower leg) | 815 / 1,071 | 0.761 | 8.85 | 16 | 0.706 | 12 | +0.055 | 26 |
| skin | Skin - Not Sun Exposed (Suprapubic) | 811 / 1,071 | 0.757 | 8.52 | 21 | 0.702 | 18 | +0.055 | 25 |
| breast | Breast - Mammary Tissue | 752 / 944 | 0.797 | 8.44 | 15 | 0.705 | 13 | +0.092 | 22 |
| adipose | Adipose - Subcutaneous | 728 / 937 | 0.777 | 10.07 | 25 | 0.683 | 34 | +0.094 | 15 |
| adipose | Adipose - Visceral (Omentum) | 722 / 937 | 0.771 | 8.45 | 33 | 0.688 | 31 | +0.083 | 32 |
| prostate | Prostate | 672 / 780 | 0.862 | 12.23 | 2 | 0.728 | 4 | +0.134 | 1 |
| ovary | Ovary | 220 / 291 | 0.756 | 6.60 | 9 | 0.689 | 28 | +0.067 | 1 |
| pancreas | Pancreas | 125 / 188 | 0.665 | 2.78 | 46 | 0.618 | 49 | +0.047 | 8 |

Skeletal muscle: 3 of 5 protein-coding genes with a row are at or above 1 TPM in skeletal muscle. The same count of 3 occurs in other tissues, and 4 of 5 are at or above 1 TPM in both minor salivary gland and testis. Those ranks are not a muscle program.

Reading the index table:

- Testis is still where the most genes clear 1 TPM. For blood-cell genes the testis fraction is 6,014 / 6,646, and the median TPM there is 13.54, against 2.93 in whole blood. Liver is the same shape: median TPM 4.36 in liver against 11.47 in testis. Heart atrial appendage is 5.34 against 10.09 in testis.
- Prostate, thyroid, and lung are the index tissues where the typical level is at least as high as testis. Prostate median TPM is 12.23 against 10.44 in testis, with detection still higher in testis (0.914 against 0.862). Thyroid is 13.60 against 10.37. Lung is 10.61 against 9.67.
- An index tissue is the largest excess over the protein-coding baseline for brain, heart, kidney, liver, ovary, and prostate. Brain cortex is excess rank 1. The thirteen brain columns occupy thirteen of the fourteen largest excesses; pituitary is tenth and is not an index tissue. Heart atrial appendage and left ventricle are excess ranks 1 and 2. Kidney medulla and cortex are 1 and 2. Liver is 1, by +0.136. Prostate is 1, by +0.134, and absolute rank 2.
- Whole blood does not rise in the absolute order. It is rank 54 of 54 for blood-cell genes, the same place it occupies in the background. The fraction is higher than the baseline (0.637 against 0.515, excess +0.122, excess rank 2). The largest blood-cell excess is spleen (+0.127, 5,451 / 6,646), which is not in the tissue map.
- Adipose, artery, breast, gut, and lung have their largest excess in a brain region, not in the index tissue. Skin's largest excess is lung. Thyroid's largest excess is small intestine (+0.147) even though thyroid itself is +0.129 and absolute rank 5. Pancreas's largest excess is small intestine; the pancreas excess is about 9 genes on a denominator of 188.
- Ovary's excess rank is 1, on 291 genes, and the next tissues are cervix. That is a small set.

The all-biotype top tissue is also testis in every group, including skeletal muscle (5 of 7 genes). Restricting the location table to protein-coding genes raises the fractions and does not change which tissue is first, except the five-gene muscle tie noted above.

## Caveat

This is a description of catalog-mapped genes. It is not colocalization, not an eQTL result, and not a statement that the index tissue is the mechanism of the trait. The mapped-gene field includes neighboring genes the catalog listed on the clump. More than half of the union (7,845 / 14,177) sits in more than one group, so the tissue profiles are not independent samples.

The name rules miss true tissue traits that do not use the tissue word, and they still include imaging, family history, and circulating biomarkers. Liver enzymes and bilirubin are 2,577 of 3,388 liver locus rows. Creatinine and filtration are 1,782 of 2,973 kidney locus rows. Prostate-specific antigen is 330 of 1,089 prostate locus rows. Cardiac troponin is a blood measurement kept in the heart group because the name says cardiac. Mean arterial pressure is in the artery group. Coronary artery disease is in the artery group, not the heart group.

GTEx median TPM is bulk tissue, 54 columns, not cell types. Whole blood is one column. Many genes cross 1 TPM in testis, pituitary, lung, and prostate whether or not the trait names those organs. A private call can be a cell-culture column. Tau >= 0.95 is a concentration cut, not a claim that the gene is absent from every other tissue at a level that matters, and it is not proof of a single-tissue program: the modal peak is testis for every large group except liver.

The 0.5 and 0.95 cuts are choices. They are not a published housekeeping standard. The observed tau minimum is 0.289812, so a 0.3 housekeeping cut is empty here. Genes between the cuts are the majority and are reported as intermediate rather than forced into either class.

The 24,065 GENCODE genes with no expression row remain the genome-wide denominator. Inside these mapped sets the missing-row count is 242 of 14,177 genes, 131 of them protein-coding. Quoting only the mapped-set rate would hide the version gap. Quoting only 24,065 would hide that the catalog symbols mostly join.

Skeletal muscle is 19 locus rows and 5 protein-coding genes with a GTEx row. Pancreas is 188 such genes. Fractions on those sets move when a handful of genes move. They are counts, not programs.
