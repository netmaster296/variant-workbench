# Noncoding genes

## Question

What do expression, GWAS overlap, and nearby protein-coding genes look like for lncRNA and other non-protein-coding genes?

These counts are not a protein-coding ranking. They should not be added to one.

## Filters

The gene set is `warehouse/marts/gene_mart.parquet`, one row per GENCODE gene (78,691). The protein-coding copy of the mart was not read. Counts below were produced by `warehouse/reports/sql/noncoding_genes.sql`.

Included biotypes are every biotype except `protein_coding` and the immunoglobulin and T-cell receptor gene segments. Those segments encode polypeptides, so they are not noncoding genes. `translated_processed_pseudogene` stays with the pseudogenes. `TEC` is "to be experimentally confirmed," not lncRNA. `artifact` is the GENCODE artifact biotype. Both are included and kept in their own rows.

Included, 58,182 genes:

| class | biotypes | genes |
| --- | --- | ---: |
| lncRNA | `lncRNA` | 34880 |
| other ncRNA | `misc_RNA`, `snRNA`, `miRNA`, `snoRNA`, `scaRNA`, `rRNA`, `Mt_tRNA`, `ribozyme`, `sRNA`, `vault_RNA`, `Mt_rRNA` | 7066 |
| pseudogene | `processed_pseudogene`, `unprocessed_pseudogene`, `transcribed_unprocessed_pseudogene`, `transcribed_processed_pseudogene`, `rRNA_pseudogene`, `transcribed_unitary_pseudogene`, `IG_V_pseudogene`, `unitary_pseudogene`, `TR_V_pseudogene`, `IG_C_pseudogene`, `TR_J_pseudogene`, `IG_J_pseudogene`, `translated_processed_pseudogene`, `IG_pseudogene` | 15198 |
| TEC | `TEC` | 1019 |
| artifact | `artifact` | 19 |

Excluded: `protein_coding` 20,097. Immunoglobulin and T-cell receptor segments 412: `IG_V_gene` 146, `TR_V_gene` 107, `TR_J_gene` 79, `IG_D_gene` 37, `IG_J_gene` 18, `IG_C_gene` 14, `TR_C_gene` 6, `TR_D_gene` 5. No other biotype is in the mart.

Expression is the GTEx v8 median-TPM summary already on the mart (`max_median_tpm`, `n_tissues_tpm_ge_1`, `n_tissues_tpm_ge_10`, `tau`, `max_tissue`). Tau is `sum(1 - x_i / x_max) / (n - 1)` over the 54 tissue columns, and it is null when `x_max` is 0. A gene has expression when `max_median_tpm > 0`. A null TPM means the bare ENSG is not in the GTEx v8 median table. A zero means it is in the table and the maximum is 0. Those are not the same fact. GTEx v8 is older than GENCODE v49.

GWAS overlap is `gwas_loci_by_overlap > 0`. That column counts clumped loci on the same chromosome whose member interval meets the gene span, or whose lead position falls inside the span. The interval is `loci.start`–`loci.end`, the minimum and maximum member positions, not the ±1,000,000 bp radius used to build the clump. The locus key is trait URI, chromosome, lead position, and lead rsID.

A GWAS mapped-gene token is `gwas_loci_by_mapped_gene > 0`. A token equals the gene symbol. The match does not require the same chromosome and does not use a substring. A hyphenated catalog string such as `FBLN1 - LINC01589` is one token.

Nearby means the same `chrom` value and a distance of at most 10,000 bp between gene spans. Spans are 1-based and inclusive. Distance is 0 when the spans share a base, including a noncoding gene inside a coding intron. Otherwise it is the number of bases strictly between the spans (`right_start - left_end - 1`). The neighbor has biotype `protein_coding` only. A 100,000 bp window is reported only as a sensitivity check. The mart chromosomes are `chr1`–`chr22`, `chrX`, `chrY`, and `chrM`. Distance is not measured from the transcription start site.

The same SQL recomputed both GWAS columns from `warehouse/gwas/loci.parquet`. Gene-level mismatches against the mart were 0 for overlap and 0 for mapped tokens, on all 78,691 genes. Mart coordinates, symbol, strand, and biotype matched `warehouse/intervals/genes.parquet` for all 78,691 genes. Expression columns matched `warehouse/gtex/gene_expression.parquet` for all 54,626 genes in that file (0 value mismatches). Those joins are intact.

## Denominator

| set | genes |
| --- | ---: |
| included biotypes | 58182 |
| has expression (`max_median_tpm > 0`) | 23621 |
| GWAS overlap | 55073 |
| GWAS mapped-gene token | 7185 |

The mapped-token count copies one token onto every GENCODE gene with that symbol. 5,627 of the 7,185 have a symbol that belongs to one GENCODE gene. 1,558 share a symbol. Those two figures are splits of the 7,185, not extra denominators.

Class denominators:

| class | genes | has expression | GWAS overlap | mapped token | unique-symbol token |
| --- | ---: | ---: | ---: | ---: | ---: |
| lncRNA | 34880 | 11968 | 33741 | 4508 | 4246 |
| other ncRNA | 7066 | 1404 | 6715 | 1177 | 129 |
| pseudogene | 15198 | 9426 | 13634 | 1494 | 1247 |
| TEC | 1019 | 811 | 983 | 6 | 5 |
| artifact | 19 | 12 | 0 | 0 | 0 |
| included | 58182 | 23621 | 55073 | 7185 | 5627 |

## Result

Most lncRNAs are missing from the GTEx table rather than measured at zero. Other ncRNA genes are usually in the table and usually undetected. A GWAS span overlap is common for every class that has one, including protein-coding genes, and it usually does not put the lead variant inside the noncoding gene. A mapped-gene token is uncommon, except where one RNA name is copied onto many genes. A bit more than half of included genes have a protein-coding gene within 10 kb, and about 40% share bases with one.

Expression, split into absence from the table, a measured zero, and detection:

| class | genes | not in GTEx | in GTEx, max TPM 0 | max TPM > 0 | TPM >= 1 in a tissue | TPM >= 10 in a tissue |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lncRNA | 34880 | 22044 | 868 | 11968 | 8168 | 2073 |
| other ncRNA | 7066 | 656 | 5006 | 1404 | 547 | 61 |
| pseudogene | 15198 | 418 | 5354 | 9426 | 3473 | 746 |
| TEC | 1019 | 57 | 151 | 811 | 447 | 99 |
| artifact | 19 | 0 | 7 | 12 | 6 | 2 |
| included | 58182 | 23175 | 11386 | 23621 | 12641 | 2981 |

Of 12,836 lncRNAs that are in the GTEx table, 11,968 have max TPM > 0 and 8,168 reach TPM >= 1 in at least one tissue. Of 6,410 other ncRNA genes in the table, 1,404 are detected and 547 reach TPM >= 1. Of 14,780 pseudogenes in the table, 9,426 are detected and 3,473 reach TPM >= 1.

Among genes with max TPM > 0, the signal is low and concentrated. Medians are DuckDB `median`. Tau quartiles are `quantile_cont`. Tissue counts are tissues with median TPM >= 1.

| class | expressed | median max TPM | tau p25 | tau p50 | tau p75 | tissues >= 1, p25 | p50 | p75 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| lncRNA | 11968 | 2.30358 | 0.837 | 0.949 | 0.995 | 0 | 1 | 8 |
| other ncRNA | 1404 | 0.650693 | 0.837 | 0.960 | 0.997 | 0 | 0 | 2 |
| pseudogene | 9426 | 0.529995 | 0.832 | 0.962 | 0.998 | 0 | 0 | 1 |
| TEC | 811 | 1.33281 | 0.760 | 0.894 | 0.977 | 0 | 1 | 10.5 |
| artifact | 12 | 1.3968825 | 0.551 | 0.828 | 0.921 | 0 | 2 | 47.25 |
| included | 23621 | 1.20867 | 0.832 | 0.953 | 0.996 | 0 | 1 | 4 |

Tau >= 0.85 for 17,060 of 23,621 expressed included genes, and for 8,717 of 11,968 expressed lncRNAs. Tau < 0.5 for 74 included genes and 39 lncRNAs. 10,980 expressed included genes, including 3,800 lncRNAs, never reach TPM >= 1 in any tissue. The tissue with the highest median TPM is testis for 5,301 of 11,968 expressed lncRNAs, then brain cerebellum 552, brain cerebellar hemisphere 549, pituitary 470, and thyroid 387. Testis is also the most common maximum for expressed pseudogenes (3,778 of 9,426), other ncRNA genes (266 of 1,404), and TEC genes (223 of 811). `max_tissue` is the highest median, not a specificity call by itself. The GTEx lane broke one genome-wide tie by taking the leftmost tissue column.

GWAS and distance to the nearest `protein_coding` gene:

| class | GWAS overlap | lead inside the gene | overlap, lead outside | mapped token | coding span overlap | within 10 kb | 10 kb flank only | within 100 kb | median nearest bp | median nearest bp if no span overlap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| lncRNA | 33741 | 15291 | 18450 | 4508 | 14392 | 20084 | 5692 | 28075 | 2579.5 | 41220.5 |
| other ncRNA | 6715 | 75 | 6640 | 1177 | 3764 | 4633 | 869 | 6118 | 0 | 37570.5 |
| pseudogene | 13634 | 1078 | 12556 | 1494 | 4418 | 6647 | 2229 | 11840 | 15980.5 | 39233 |
| TEC | 983 | 103 | 880 | 6 | 517 | 674 | 157 | 853 | 0 | 37133 |
| artifact | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4430145 | 4430145 |
| included | 55073 | 16547 | 38526 | 7185 | 23091 | 32038 | 8947 | 46886 | 4998.5 | 40075 |

Every included gene has a `protein_coding` gene on its chromosome. The nearest distance runs from 0 to 5,504,407 bp. Quartiles of that distance for all included genes are 0, 4,998.5, and 60,931.75 bp. Thirty genes abut a protein-coding gene without sharing a base. They are in the 10 kb flank, not in the span-overlap count. The other-ncRNA median of 0 bp is the span overlaps: 3,764 of 7,066 other ncRNA genes share bases with a protein-coding gene. `Mt_tRNA` (22) and `Mt_rRNA` (2) are all within 10 kb of a protein-coding gene, all have max TPM > 0, and none has a GWAS overlap or a mapped token. chrM is short, so that 10 kb window is not the same neighborhood as on a nuclear chromosome.

Of 55,073 included genes with a GWAS overlap, 31,364 are within 10 kb of a protein-coding gene, 22,727 share bases with one, and 23,709 are not within 10 kb. Only 7,130 of those 55,073 are themselves a mapped-gene token. The median overlapping gene meets 14 loci and 14 traits. For lncRNA the median is 15 loci and 15 traits. A lead position falls inside the gene for 16,547 included genes (15,291 lncRNAs, 1,078 pseudogenes, 103 TEC genes, 75 other ncRNA genes). The other 38,526 overlaps do not contain the lead.

Of 7,185 included genes with a mapped token, 4,590 are within 10 kb of a protein-coding gene, 3,465 share bases with one, and 2,595 are not within 10 kb. For lncRNA, 3,032 of 4,508 token matches are within 10 kb and 2,466 share bases with a protein-coding gene. Among token matches whose symbol is unique in GENCODE, the median is 4 loci and 4 traits overall, and 5 loci and 5 traits for lncRNA. It is 1 locus and 1 trait for other ncRNA, pseudogenes, and TEC. 55 included genes have a mapped token and no span overlap. 22,728 genes with expression also have a GWAS overlap. 5,569 genes with expression also have a mapped token, and 5,563 of those have both a token and a span overlap.

The mapped-token count is not 7,185 independent catalog assignments. Five repeated names account for 1,036 included genes, and each copy receives the same locus count:

| symbol | biotype | genes | mapped loci on each gene |
| --- | --- | ---: | ---: |
| Y_RNA | misc_RNA | 756 | 18 |
| Metazoa_SRP | misc_RNA | 170 | 14 |
| U3 | snoRNA | 50 | 8 |
| U6 | snRNA | 33 | 4 |
| SNORA70 | snoRNA | 27 | 1 |

That is why other ncRNA genes show 1,177 token matches but only 129 genes whose symbol is unique. `misc_RNA` alone is 967 token matches, 41 of them unique. lncRNA token matches are mostly unique symbols (4,246 of 4,508).

Nine lncRNAs use a symbol that also belongs to a `protein_coding` gene. The token is attached to both genes. Eight of the nine have a mapped-locus count. Those counts are not a noncoding assignment.

| symbol | ensg | mapped loci | overlap loci |
| --- | --- | ---: | ---: |
| CAST | ENSG00000153113 | 119 | 95 |
| ELFN2 | ENSG00000243902 | 4 | 23 |
| GOLGA8M | ENSG00000261480 | 1 | 7 |
| KBTBD11-OT1 | ENSG00000253696 | 43 | 12 |
| LINC02203 | ENSG00000284988 | 0 | 0 |
| MALL | ENSG00000293619 | 1 | 2 |
| SFTA3 | ENSG00000229415 | 8 | 17 |
| SIGLEC5 | ENSG00000105501 | 33 | 34 |
| SPATA13 | ENSG00000228741 | 48 | 15 |

LOEUF is non-null for 23 included genes (8 lncRNA, 15 pseudogene) and null for 58,159. The nulls are the expected state for this noncoding set. They are not treated as a missing join, and the 23 values are not listed and are not a ranking.

Same three rates for `protein_coding`, as a scale only. This row is not part of the 58,182 and is not a ranking.

| genes | not in GTEx | max TPM 0 | max TPM > 0 | TPM >= 1 | TPM >= 10 | GWAS overlap | mapped token | unique-symbol token | LOEUF non-null | median max TPM | median tau | median tissues >= 1 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20097 | 888 | 342 | 18867 | 17985 | 15830 | 19713 | 16499 | 16457 | 18098 | 44.5677 | 0.810 | 52 |

Biotype detail for the included set:

| biotype | genes | not in GTEx | max TPM 0 | max TPM > 0 | TPM >= 1 | TPM >= 10 | GWAS overlap | mapped token | unique-symbol token | coding span overlap | within 10 kb |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| lncRNA | 34880 | 22044 | 868 | 11968 | 8168 | 2073 | 33741 | 4508 | 4246 | 14392 | 20084 |
| misc_RNA | 2207 | 204 | 1368 | 635 | 207 | 10 | 2119 | 967 | 41 | 1119 | 1465 |
| snRNA | 1901 | 37 | 1552 | 312 | 103 | 6 | 1790 | 75 | 35 | 923 | 1139 |
| miRNA | 1879 | 306 | 1407 | 166 | 73 | 6 | 1793 | 43 | 43 | 1208 | 1372 |
| snoRNA | 942 | 96 | 600 | 246 | 130 | 26 | 915 | 85 | 6 | 473 | 574 |
| scaRNA | 49 | 9 | 25 | 15 | 7 | 0 | 49 | 2 | 2 | 33 | 43 |
| rRNA | 47 | 1 | 42 | 4 | 2 | 0 | 33 | 4 | 1 | 5 | 8 |
| Mt_tRNA | 22 | 0 | 0 | 22 | 22 | 11 | 0 | 0 | 0 | 0 | 22 |
| ribozyme | 8 | 2 | 6 | 0 | 0 | 0 | 8 | 0 | 0 | 2 | 6 |
| sRNA | 5 | 1 | 4 | 0 | 0 | 0 | 4 | 0 | 0 | 1 | 1 |
| vault_RNA | 4 | 0 | 2 | 2 | 1 | 0 | 4 | 1 | 1 | 0 | 1 |
| Mt_rRNA | 2 | 0 | 0 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 2 |
| processed_pseudogene | 9487 | 118 | 3617 | 5752 | 1566 | 187 | 8751 | 576 | 576 | 3293 | 4428 |
| unprocessed_pseudogene | 1949 | 88 | 940 | 921 | 209 | 31 | 1570 | 180 | 180 | 332 | 718 |
| transcribed_unprocessed_pseudogene | 1587 | 106 | 169 | 1312 | 889 | 327 | 1311 | 463 | 316 | 192 | 544 |
| transcribed_processed_pseudogene | 1149 | 58 | 139 | 952 | 572 | 126 | 1025 | 110 | 71 | 288 | 481 |
| rRNA_pseudogene | 497 | 16 | 404 | 77 | 24 | 0 | 474 | 15 | 9 | 241 | 291 |
| transcribed_unitary_pseudogene | 201 | 18 | 3 | 180 | 133 | 65 | 194 | 114 | 59 | 45 | 134 |
| IG_V_pseudogene | 187 | 0 | 51 | 136 | 54 | 7 | 171 | 9 | 9 | 0 | 0 |
| unitary_pseudogene | 89 | 13 | 17 | 59 | 10 | 1 | 86 | 21 | 21 | 22 | 45 |
| TR_V_pseudogene | 33 | 0 | 9 | 24 | 9 | 0 | 33 | 2 | 2 | 1 | 2 |
| IG_C_pseudogene | 9 | 0 | 1 | 8 | 6 | 2 | 9 | 4 | 4 | 3 | 3 |
| TR_J_pseudogene | 4 | 0 | 3 | 1 | 0 | 0 | 4 | 0 | 0 | 0 | 0 |
| IG_J_pseudogene | 3 | 0 | 1 | 2 | 0 | 0 | 3 | 0 | 0 | 0 | 0 |
| translated_processed_pseudogene | 2 | 1 | 0 | 1 | 0 | 0 | 2 | 0 | 0 | 1 | 1 |
| IG_pseudogene | 1 | 0 | 0 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| TEC | 1019 | 57 | 151 | 811 | 447 | 99 | 983 | 6 | 5 | 517 | 674 |
| artifact | 19 | 0 | 7 | 12 | 6 | 2 | 0 | 0 | 0 | 0 | 0 |

## Caveat

A GWAS overlap here is not evidence that the noncoding gene is the trait gene. Clump member spans cover most of the primary assembly (`warehouse/reports/dark_genome.md` records 0.886 of those bases). Protein-coding genes show the same pattern: 19,713 of 20,097 have an overlap, against 55,073 of 58,182 included genes. For 38,526 included genes the overlap does not place the lead variant inside the gene. The median overlapping gene meets 14 loci. The mapped-gene token is the catalog's own gene field, and it is the stricter of the two counts. It is still not colocalization. It has no chromosome filter, and a shared symbol is copied onto every gene that uses it. `Y_RNA` alone is 756 of the 1,177 other-ncRNA token matches.

A protein-coding gene within 10 kb, and especially a shared span, is the competing annotation. 22,727 genes with a GWAS overlap share bases with a protein-coding gene. 2,466 lncRNAs that are themselves a mapped token also share bases with one. Distance 0 includes an intronic or antisense gene whose transcription start site may be far from the coding start site. 238 included genes are within 10 kb of an immunoglobulin or T-cell receptor segment and are not within 10 kb of a `protein_coding` gene. They are not counted as nearby.

Not in the GTEx table is not a measured zero. 22,044 of 34,880 lncRNAs are absent because GTEx v8 uses an older gene set than GENCODE v49. Among lncRNAs that are in the table, detection is common and the level is still low: median max TPM 2.30, and the median gene reaches TPM >= 1 in one tissue. The GTEx v8 median-TPM file is poly(A) RNA-seq, not a small-RNA assay. Zeros for miRNA, snRNA, snoRNA, and the other short classes should be read that way. Testis as the maximum tissue is common in this file for several classes. It is not, by itself, a tissue-specific function.

Pseudogenes, TEC, and artifacts are in the denominator and are not lncRNAs. The artifact class is 19 genes, with no GWAS overlap and a median distance of 4,430,145 bp to the nearest protein-coding gene. Dropping them would not change the totals in a meaningful way. chrM contributes 24 genes that all fall inside the 10 kb window for a geometric reason.

LOEUF is null for 58,159 of 58,182 included genes. That absence is expected for noncoding genes. It is not a broken constraint join, and this report does not rank genes by constraint. eQTL, sQTL, ClinVar, and MaveDB columns were not used.
