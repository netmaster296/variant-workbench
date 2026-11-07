# Constraint versus function

## Question

Does LOEUF decile track the fraction of ClinVar alleles that are pathogenic at two or more stars, and does it track MaveDB effect sizes where an assay exists?

## Filters

Protein-coding genes are the rows in `D:\variant-workbench\warehouse\marts\gene_mart_protein_coding.parquet`.

LOEUF decile is `loeuf_decile`, the gnomAD bin `lof.oe_ci.upper_bin_decile`. Decile 0 is the most constrained bin. In this file the observed LOEUF runs from 0.040 to 0.400 in decile 0 and from 1.801 to 1.998 in decile 9. A gene with no decile is not placed in a bin.

No constraint row means `constraint_flags` is null. A matched constraint row stores `[]` when it has no flag, so a null flag is the miss, not an empty flag list.

ClinVar counts are the mart's symbol-token counts from `allele_dedup.parquet`. The numerator is `clinvar_germline_pathogenic_star_2` + `star_3` + `star_4`. That bucket is germline pathogenic, including a `Pathogenic/Likely pathogenic` call. A call that is only `Likely pathogenic` is not in it. Somatic origin, oncogenicity, and conflicting text are not in it. Stars below 2, and a null star, are not in it. The denominator is `clinvar_n`, every ClinVar allele matched to the symbol, not only the pathogenic ones. The decile fraction is the sum of the numerator divided by the sum of `clinvar_n`. Genes with no ClinVar allele do not move that ratio. The unweighted mean uses genes with `clinvar_n > 0`.

Likely pathogenic at two or more stars is reported beside the pathogenic fraction. It is not added into the primary numerator.

MaveDB scores are not pooled across score sets. An assay counted here is a score set whose `target_name` equals the protein-coding symbol and whose `organism_name` is `Homo sapiens` or null. The score is numeric `score` on `mapped_variants_scored.parquet`, one value per `(score_set_urn, urn)`. Gallus gallus and Escherichia coli score sets are not treated as assays of the human gene. A position inside a gene span is not treated as that gene's assay. That join is in `D:\variant-workbench\warehouse\qc\issues_constraint_versus_function.md`. The effect-size claim that would use it is stopped.

## Denominator

20097 protein-coding genes.

1530 have no constraint row. 727 more have a constraint row and a null decile: 469 of those have a null LOEUF, and 258 have a LOEUF and no decile. 2257 genes have no decile. 17840 have a decile and are the rows in the table below.

19041 of the 20097 genes have at least one matched ClinVar allele. 3140 have at least one germline pathogenic allele at two or more stars.

Two symbols, PINX1 and POLR2J3, are two protein-coding genes each and fall in two deciles (PINX1 deciles 8 and 9, 107 alleles; POLR2J3 deciles 2 and 5, 9 alleles). Both have zero pathogenic alleles at two or more stars, so they do not change a numerator. No symbol is repeated inside one decile.

## Result

LOEUF decile does not order the pathogenic fraction. The median gene's fraction is 0 in every decile. The allele-weighted fraction is 0.016037 in decile 0 and 0.010789 in decile 9, and it is highest in decile 3 (0.028357), not in the most constrained bin. Adding likely pathogenic alleles does not move the peak.

The rank correlation is weak. Spearman correlation of decile with the gene-level fraction is -0.1334 on 17712 genes that have a decile and at least one ClinVar allele, and -0.1249 on 17118 genes with at least 10 alleles. Ties use average ranks.

What does fall toward the unconstrained end is how many genes have any such allele: 506 of 1851 genes in decile 0, and 73 of 1327 in decile 9. That indicator's Spearman correlation with decile is -0.1451 on all 17840 genes with a decile. The decline is not a straight line. Deciles 3 and 4 sit above deciles 1 and 2.

The decile 3 count is carried by a few genes. BRCA2 (decile 3, LOEUF 0.803) has 3605 of that decile's 17118 pathogenic alleles, out of 21849 ClinVar alleles on the gene. ATM (1554) and MSH6 (1065) are next in the same decile. BRCA1 is decile 4 (2755 of 15639). The largest decile 0 numerators are APC (989), FBN1 (874), and the rest of the top of decile 0 are smaller than BRCA2 alone.

MaveDB does not support an effect-size comparison, and it is not genome-wide.

134 target names in `gene_coverage.parquet` exactly match a protein-coding symbol. That is 134 genes. The other 19963 protein-coding genes do not. 18 of the 134 names have no score set of their own. They are only tokens on a multi-target score set, so `mavedb_experiments_by_symbol` on the mart is not "this gene was assayed." Those names are ABCD1, BAP1, BCL10, BRIP1, CCM2, FANCC, FARS2, FH, FLCN, NF2, PTPN11, RAF1, RET, RPS19, SDHD, SHOC2, SOS2, and ZC4H2. Four names match only Gallus gallus score sets: EN2, FYN, SRC, and VIL1.

112 genes have a sole-target score set with organism Homo sapiens or with no organism recorded (441 score sets). Six of those genes have only an unrecorded organism: BARD1, CTCF, PALB2, RAD51D, SFPQ, and XRCC2. Their score-set titles name those genes. They are not a named non-human organism. 102 of the 112 genes have a numeric mapped score, on 388 score sets and 1932514 variant urns. Ten genes have a qualifying score set and no numeric mapped score.

| decile | genes with a numeric score | score sets | scored urns |
| --- | ---: | ---: | ---: |
| 0 | 24 | 98 | 398188 |
| 1 | 18 | 69 | 486641 |
| 2 | 9 | 44 | 356230 |
| 3 | 6 | 25 | 142496 |
| 4 | 9 | 26 | 73905 |
| 5 | 10 | 48 | 140273 |
| 6 | 6 | 18 | 41063 |
| 7 | 5 | 10 | 38869 |
| 8 | 7 | 22 | 66849 |
| 9 | 1 | 3 | 9126 |
| no decile | 7 | 25 | 178874 |

Those 102 genes are not a draw of the 20097. One gene in decile 9 has a numeric score.

Scores are not on one scale. Across the 388 score sets, the score-set medians run from -29.3403863 to 127.0408005. The median of those medians is 0.3743261682276327. The raw scores run from -71.39392555 to 3462.809917355372. 141 sets have a negative median. 82 have a median above 1. 20 lie entirely inside 0 to 1. 70 have a width above 10. 300 cross zero. 86 genes have two or more score sets, and 45 of those have score-set medians more than 1 apart.

KCNQ2 is decile 0 and has 17 score sets. The lowest median is -29.3403863 (`urn:mavedb:00000664-a-14`, 57 urns). The highest is 127.0408005 (`urn:mavedb:00000664-a-13`, 57 urns). PTEN, decile 2, score set `urn:mavedb:00001244-b-1`, has 1228 urns with scores from 0.0891674398324491 to 0.8492537790420092 (median 0.5431423224794725). KCNH2, decile 1, score set `urn:mavedb:00001231-a-1`, has 22824 urns with scores from 0 to 352.187353715064 (median 85.057314389175). CYP2C9, decile 8, score set `urn:mavedb:00000062-a-1`, has 109 urns inside 0.25 to 0.74859721 (median 0.64589733). A mean of these scores is not an effect size, and it is not a LOEUF comparison.

Separately, 124 genes of any biotype contain a GRCh38 MaveDB key in the gene span: 98 protein-coding and 26 lncRNA. That count is not an assay count. The span join attributes keys to genes the assay did not name. It is not used for the paragraphs above.

| decile | genes | genes with a ClinVar allele | ClinVar alleles | pathogenic at 2+ stars | likely pathogenic at 2+ stars | allele fraction pathogenic | allele fraction pathogenic or likely pathogenic | mean gene fraction | genes with any pathogenic allele at 2+ stars | median LOEUF |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 1851 | 1849 | 847453 | 13591 | 2023 | 0.016037 | 0.018425 | 0.005000 | 506 | 0.288 |
| 1 | 1844 | 1841 | 638372 | 7074 | 1974 | 0.011081 | 0.014174 | 0.002588 | 331 | 0.497 |
| 2 | 1842 | 1838 | 581553 | 9242 | 1947 | 0.015892 | 0.019240 | 0.003089 | 333 | 0.664 |
| 3 | 1834 | 1825 | 603670 | 17118 | 2296 | 0.028357 | 0.032160 | 0.004600 | 381 | 0.804 |
| 4 | 1834 | 1824 | 499566 | 12854 | 2109 | 0.025730 | 0.029952 | 0.005293 | 383 | 0.924 |
| 5 | 1838 | 1829 | 414182 | 9899 | 2246 | 0.023900 | 0.029323 | 0.005145 | 315 | 1.038 |
| 6 | 1835 | 1821 | 291520 | 4935 | 814 | 0.016929 | 0.019721 | 0.003891 | 287 | 1.162 |
| 7 | 1825 | 1802 | 241630 | 4368 | 1148 | 0.018077 | 0.022828 | 0.003248 | 199 | 1.329 |
| 8 | 1810 | 1785 | 175136 | 2109 | 415 | 0.012042 | 0.014412 | 0.002322 | 141 | 1.607 |
| 9 | 1327 | 1298 | 105477 | 1138 | 222 | 0.010789 | 0.012894 | 0.001561 | 73 | 1.892 |

Fractions and the mean are `round(..., 6)`. The median gene fraction is 0 in every row. The largest numerators:

| symbol | decile | LOEUF | ClinVar alleles | pathogenic at 2+ stars | fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| BRCA2 | 3 | 0.803 | 21849 | 3605 | 0.1650 |
| BRCA1 | 4 | 0.885 | 15639 | 2755 | 0.1762 |
| NF1 | 1 | 0.449 | 17644 | 1706 | 0.0967 |
| ATM | 3 | 0.847 | 19656 | 1554 | 0.0791 |
| MSH6 | 3 | 0.841 | 10624 | 1065 | 0.1002 |
| APC | 0 | 0.321 | 17255 | 989 | 0.0573 |
| MSH2 | 2 | 0.672 | 8268 | 981 | 0.1187 |
| MLH1 | 2 | 0.639 | 6373 | 901 | 0.1414 |
| FBN1 | 0 | 0.105 | 9655 | 874 | 0.0905 |
| LDLR | 5 | 1.091 | 4659 | 826 | 0.1773 |

## Caveat

ClinVar alleles are the alleles that were submitted, not a sample of possible alleles. A gene that has been sequenced in many patients contributes more alleles, and more pathogenic labels, than a gene that has not. BRCA1 and BRCA2 are not in decile 0. Their LOEUF values sit in deciles 4 and 3, and they dominate those numerators. The fraction can be high because of that attention rather than because the decile is a measure of pathogenicity. A one-star or conflicting call is outside the numerator on purpose. The mart can still attach one symbol's alleles to two genes when the symbol is shared. Inside deciles 0 through 9 that happens only for PINX1 and POLR2J3, and neither has a pathogenic allele at two or more stars.

A MaveDB score is not comparable across score sets, including score sets of the same gene. No genome-wide functional claim follows from 102 assayed genes, or from 134 name matches, out of 20097. Names such as `BRCA1 RING domain` or `DDX3X Exon 12` do not exact-match a symbol, so the 134 is not every assay that mentions a gene. Keys in a gene span are a different and unreliable map, recorded in the issues note, and they are not the effect sizes above.

The numbers are the result of `D:\variant-workbench\warehouse\reports\sql\constraint_versus_function.sql`. The marts were not rewritten.
