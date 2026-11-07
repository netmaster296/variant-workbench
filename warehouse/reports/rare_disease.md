# Rare disease

Which protein-coding genes combine loss-of-function intolerance, high-confidence germline pathogenic alleles, and broad expression?

## Filters

Every filter used for the lists is named here. The counts come from `warehouse/marts/gene_mart_protein_coding.parquet`. Two audits read `warehouse/constraint/gene_constraint.parquet` and `warehouse/gtex/gene_expression.parquet` only to check the joins. No mart was rewritten, and a null LOEUF was not filled by gene symbol.

1. Universe. One row per protein-coding GENCODE gene in `gene_mart_protein_coding.parquet`. The file has no other biotype.

2. Loss-of-function intolerance. `loeuf` is not null and `loeuf < 0.6`. On this mart, `loeuf` is gnomAD v4.1 `lof.oe_ci.upper` from the constraint row joined on bare ENSG. The cutoff is the one gnomAD stated for v4 constraint: the 8 March 2024 v4.0 note recommended LOEUF < 0.6 in place of the v2 cutoff 0.35, and the October 2024 discussion of the v4.1 file restated `< 0.6`. A 30 March 2026 v4.1.1 note recommended `< 0.45` for a later file that is not on disk. That later cutoff is a sensitivity count in the result, not the filter. Equality at 0.6 fails. A null LOEUF fails. LOEUF decile, pLI, and `constraint_flags` are not filters.

3. High-confidence germline pathogenic alleles. `plp_stars_ge_2` is at least 1. That column, defined in the SQL, is the sum of `clinvar_germline_pathogenic_star_2`, `star_3`, and `star_4` and `clinvar_likely_pathogenic_star_2`, `star_3`, and `star_4`. Review stars are the mart's map of ClinVar review status: 4 practice guideline, 3 expert panel, 2 multiple submitters and no conflicts, 1 single submitter or conflicting classifications, 0 no assertion or no interpretation. A status that matched none of those phrases is the star-null bucket. Star 0, star 1, and star null are not in `plp_stars_ge_2`. The two classes already follow the deduplicated `germline_pathogenic` flag, so somatic origin, mixed germline/somatic origin, an oncogenicity call, and conflicting clinical-significance text are not in them. `Pathogenic/Likely pathogenic` is counted once, in the pathogenic class, not again as likely pathogenic.

4. Broad expression. `n_tissues_tpm_ge_1 >= 27`. The column is the number of GTEx v8 tissues with median TPM at least 1. Every expression row that joins has `n_tissues` 54, so 27 is at least half of those tissues. A null breadth fails. Tau is not a filter.

5. Order, inside the genes that pass filters 2, 3, and 4. `loeuf` ascending, then `plp_stars_ge_2` descending, then `n_tissues_tpm_ge_1` descending, then symbol, then ENSG. Rank 1 is the first row of that order. The rank is not a severity score.

6. Split-out. Same filters 1, 2, and 4. `plp_stars_ge_2` is 0. `plp_star_1` is at least 1. `plp_star_0` is 0 and `plp_star_null` is 0. These are the genes whose germline pathogenic or likely pathogenic alleles are only at star 1. Star 1 is one bucket for two review phrases, "criteria provided, single submitter" and "criteria provided, conflicting classifications." The mart does not store which phrase it was. Order is `loeuf` ascending, `plp_star_1` descending, breadth descending, symbol, ENSG.

Not filters: ancestry, allele frequency, GWAS, MaveDB, eQTL or sQTL, chromosome, drug response, VUS, benign calls, and the conflicting clinical-significance class. That conflicting class is not a pathogenic call on this mart.

## Denominator

The denominator is 20,097 protein-coding genes, 20,097 distinct ENSG ids, and 0 rows of another biotype.

18,098 of 20,097 have a non-null LOEUF. 1,999 of 20,097 have a null LOEUF. Of those 1,999, 1,530 of 20,097 have no constraint row on bare ENSG, and 469 of 20,097 have a constraint row whose LOEUF is null. The 1,530 were not filled by symbol.

The constraint file has 0 rows with chromosome chrX, chrY, or chrM. Of the 1,530 protein-coding genes with no constraint row, 867 are chrX, 61 are chrY, 13 are chrM, and 589 are on another chromosome. LOEUF on the mart matches the constraint file wherever both exist: 0 disagreements, and 0 mart LOEUF values without a constraint row. `bare_ensg` is unique in the constraint file (0 duplicates).

19,209 of 20,097 have an expression row, and 888 of 20,097 do not. Matched rows have 54 tissues, with no other value. Breadth and tau on the mart match the expression file (0 breadth disagreements, 0 tau disagreements). Among the 19,209 genes with a breadth value, the middle value of `n_tissues_tpm_ge_1` is 52 of 54. 7,807 of 20,097 are at 54, and 13,133 of 20,097 are at least 27.

Star buckets sum to each ClinVar class total on all 20,097 rows (0 mismatches). 63 of 20,097 protein-coding genes share a symbol with another GENCODE gene.

3,895 of 18,098 genes with a non-null LOEUF are below 0.6, which is 3,895 of 20,097. 2,354 of 18,098 are below 0.45. 18 of 20,097 have LOEUF equal to 0.6. 258 of 18,098 have a LOEUF and a null decile, including 27 of 18,098 below 0.6. 3,340 of 20,097 have at least one germline pathogenic or likely pathogenic allele at stars 2, 3, or 4.

## Result

755 of 20,097 protein-coding genes pass all three filters. That is 755 of 18,098 genes with a non-null LOEUF, and 755 of the 918 that pass the LOEUF and high-confidence filters before expression is applied.

| set | genes | denominator |
| --- | ---: | --- |
| protein-coding genes | 20,097 | 20,097 |
| non-null LOEUF | 18,098 | 20,097 |
| LOEUF < 0.6 | 3,895 | 18,098 with non-null LOEUF |
| breadth at least 27 of 54 | 13,133 | 20,097 |
| germline P or LP at stars >= 2 | 3,340 | 20,097 |
| LOEUF < 0.6 and stars >= 2 | 918 | 20,097 |
| all three filters | 755 | 918 that had passed the first two, and 20,097 |

Among the 755, LOEUF runs from 0.041 to 0.599, and the middle value is 0.343. The high-confidence allele count runs from 1 to 1,778, and the middle value is 6. The 755 genes carry 19,178 such alleles: 16,865 pathogenic and 2,313 likely pathogenic. 42 of the 755 have no pathogenic allele at stars >= 2 and are in the list on likely pathogenic alleles only. The largest count is 1,778, on NF1 (ENSG00000196712, LOEUF 0.449). Breadth runs from 27 to 54 of 54, and the middle value is 54.

LOEUF decile, which was not a filter: 456 of the 755 are decile 0, 274 of the 755 are decile 1, 25 of the 755 are decile 2, and 0 of the 755 have a null decile. None of the 27 genes with LOEUF < 0.6 and a null decile has a high-confidence allele, so a missing decile did not by itself remove anyone from the 755.

Review strength inside the high-confidence count: 726 of the 755 have no three-star or four-star allele, so every allele in their count is a two-star multiple-submitter call. 28 of the 755 have at least one three-star allele and no four-star allele. 1 of the 755 has a four-star allele: HTT (ENSG00000197386), and that gene's high-confidence count is the one pathogenic practice-guideline allele. These three groups do not overlap and sum to 755.

726 of the 755 also have at least one star-1 germline pathogenic or likely pathogenic allele. 697 of the 755 are in both that group and the two-star-only group. 29 of the 755 have a star-1 allele and also a three-star or four-star allele. 29 of the 755 are two-star-only and have no star-1 allele. 720 of the 755 also have at least one conflicting-classification allele. Those star-1 and conflicting alleles are not part of `plp_stars_ge_2`.

The same LOEUF and breadth filters, with no star >= 2, split out 202 of 20,097 genes whose germline pathogenic or likely pathogenic alleles are only at star 1. They carry 558 such alleles, from 1 to 31 per gene. LOEUF in that split runs from 0.082 to 0.592. Decile 0 holds 122 of the 202, decile 1 holds 73 of the 202, and decile 2 holds 7 of the 202. 96 of the 202 also have a conflicting-classification allele. 0 of the 202 share a symbol with another GENCODE gene.

Genes that miss one of the three filters, and are in neither list:

| set | genes | denominator | why it is outside both lists |
| --- | ---: | --- | --- |
| LOEUF < 0.6, stars >= 2, breadth under 27 | 162 | 918 | fails broad expression |
| LOEUF < 0.6, stars >= 2, breadth null | 1 | 918 | no GTEx v8 row for the ENSG |
| stars >= 2, breadth >= 27, LOEUF null | 164 | 20,097 | fails the LOEUF filter; 162 have no constraint row and 2 have a null metric |
| LOEUF = 0.6, stars >= 2, breadth >= 27 | 4 | 18 with LOEUF equal to 0.6 | the cutoff is strict `<` |
| LOEUF < 0.6, breadth >= 27, star 1 plus star 0 or star null | 208 | 20,097 | pathogenic calls are not only star 1 |
| LOEUF < 0.6, breadth >= 27, only star 0 or star null | 165 | 20,097 | review status is weaker than star 1 |
| LOEUF < 0.6, breadth >= 27, conflicting class, no germline P or LP | 197 | 20,097 | conflicting clinical significance is not a pathogenic call |

Using LOEUF < 0.45 instead of < 0.6, with the same ClinVar and breadth filters, would leave 525 of 20,097. 230 of the 755 have LOEUF from 0.45 up to but not including 0.6.

XPO6 (ENSG00000169180, LOEUF 0.040, breadth 54) and EIF4G1 (ENSG00000114867, LOEUF 0.048, breadth 54) have no germline pathogenic or likely pathogenic allele at any star, so they are not in either list. MED13L (ENSG00000123066, LOEUF 0.041, 36 high-confidence alleles, breadth 54) is rank 1 of the 755.

The ordered lists are below, after the caveat. `warehouse/reports/sql/rare_disease.sql` is the query that produces the counts and both lists.

## Caveat

This is the intersection of three filters. It is not a list of genes that cause rare disease, and the order is not a ranking of severity or of populations. LOEUF measures depletion of predicted loss-of-function variants in gnomAD v4.1. It does not say that a counted ClinVar allele is loss-of-function. The allele count is deduplicated ClinVar alleles, not patients. One two-star allele is enough to enter, and 726 of the 755 have no expert-panel or practice-guideline allele in the count. 720 of the 755 also carry a conflicting classification, so the high-confidence count is not the whole ClinVar record for the gene.

Broad expression, as defined here, is the usual state of a protein-coding gene in this file, not a rare property. The middle breadth among the 19,209 genes with a value is 52 of 54, and the middle breadth inside the 755 is 54 of 54. The expression filter removes 163 of the 918 genes that already passed LOEUF and ClinVar. Tau was not used. GTEx v8 is older than GENCODE v49. A null breadth means the bare ENSG is absent from the GTEx table, not that a measured count of tissues is 0.

SOX10 (ENSG00000100146) has LOEUF 0.248 and 23 high-confidence alleles, and its breadth is null. That ENSG is not in the GTEx v8 median-TPM file. It was not given a breadth by symbol, and it is not one of the 755.

A null LOEUF is not evidence that the gene tolerates loss of function. 1,530 of 20,097 have no constraint row, and 469 of 20,097 have a row with a null LOEUF. chrX, chrY, and chrM are absent from the constraint file, so none of the 867 chrX, 61 chrY, or 13 chrM protein-coding genes can enter, whatever ClinVar says. They were not matched by symbol. FANCF (ENSG00000183161) and SDHAF1 (ENSG00000205138) do have constraint rows, both flagged `no_exp_lof`, with null `lof.exp` and null LOEUF. FANCF has 17 high-confidence alleles and breadth 53. SDHAF1 has 1 high-confidence allele and breadth 54. They stay out. Another gene's LOEUF was not copied onto them by symbol. The other 162 of the 164 genes with high-confidence alleles, broad expression, and null LOEUF have no constraint row at all.

Four genes sit on the cutoff and otherwise pass. They are excluded because the filter is `loeuf < 0.6`, not `<=`. They are STK4 (ENSG00000101109, 4 alleles, breadth 54), RIPK1 (ENSG00000137275, 3 alleles, breadth 54), SLC9A3 (ENSG00000066230, 3 alleles, breadth 37), SMOC1 (ENSG00000198732, 3 alleles, breadth 44).

10 of the 755 have a non-empty `constraint_flags` value. None of the ten is `outlier_lof`. The flags were not a filter, so these genes stay in the list.

| symbol | ensg | loeuf | plp_stars_ge_2 | n_tissues_tpm_ge_1 | constraint_flags |
| --- | --- | ---: | ---: | ---: | --- |
| ANKRD11 | ENSG00000167522 | 0.132 | 163 | 54 | ["outlier_syn"] |
| PRR12 | ENSG00000126464 | 0.180 | 8 | 54 | ["outlier_syn"] |
| SRRM2 | ENSG00000167978 | 0.183 | 8 | 54 | ["outlier_mis","outlier_syn"] |
| TSC2 | ENSG00000103197 | 0.203 | 390 | 54 | ["outlier_syn"] |
| PTBP1 | ENSG00000011304 | 0.249 | 1 | 54 | ["outlier_syn"] |
| PKD1 | ENSG00000008710 | 0.364 | 636 | 54 | ["outlier_mis","outlier_syn"] |
| TUBB4B | ENSG00000188229 | 0.377 | 3 | 54 | ["outlier_syn"] |
| SBF1 | ENSG00000100241 | 0.381 | 3 | 54 | ["outlier_syn"] |
| SHANK3 | ENSG00000251322 | 0.446 | 27 | 52 | ["outlier_syn"] |
| IRF2BPL | ENSG00000119669 | 0.544 | 22 | 54 | ["outlier_syn"] |

Star 1 cannot be read as single-submitter alone. The same bucket is "criteria provided, conflicting classifications." Conflicting clinical significance is a different column. Under the LOEUF and breadth filters, 197 of 20,097 genes have only that class and no germline pathogenic or likely pathogenic allele. They are not listed as pathogenic.

ClinVar is attached by a gene-symbol token. An allele whose `GeneSymbol` lists more than one symbol is counted on each token that matches a GENCODE symbol. 63 of 20,097 protein-coding genes share a symbol with another GENCODE gene. None of those 63 are in the 755 or the 202. Thirteen symbols were recomputed from `allele_dedup.parquet` with the mart's own classifier: MED13L, TCF20, SON, WDR5, NIPBL, NF1, FBN1, APC, PSMD11, XPO1, WAPL, MALL, and C4orf36. All 13 matched the mart, including MALL and C4orf36, which share a symbol with another gene. That is a sample, not a second pass over every allele.

Allele frequency and ancestry were not filters. This report does not compare ancestry groups and does not estimate disease risk for a population.

## High-confidence genes

The 755 genes that pass filters 2, 3, and 4. `plp_stars_ge_2` is germline pathogenic plus likely pathogenic alleles at review stars 2, 3, or 4. `n_tissues_tpm_ge_1` is tissues with median TPM at least 1, out of 54. Order is filter 5.

| rank | symbol | ensg | loeuf | loeuf_decile | plp_stars_ge_2 | n_tissues_tpm_ge_1 |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 1 | MED13L | ENSG00000123066 | 0.041 | 0 | 36 | 54 |
| 2 | TCF20 | ENSG00000100207 | 0.050 | 0 | 21 | 54 |
| 3 | SON | ENSG00000159140 | 0.051 | 0 | 17 | 54 |
| 4 | WDR5 | ENSG00000196363 | 0.053 | 0 | 1 | 54 |
| 5 | NIPBL | ENSG00000164190 | 0.057 | 0 | 90 | 54 |
| 6 | CNOT3 | ENSG00000088038 | 0.057 | 0 | 10 | 54 |
| 7 | KDM2A | ENSG00000173120 | 0.057 | 0 | 2 | 54 |
| 8 | CNOT1 | ENSG00000125107 | 0.064 | 0 | 3 | 54 |
| 9 | KMT2A | ENSG00000118058 | 0.065 | 0 | 77 | 54 |
| 10 | UBR5 | ENSG00000104517 | 0.073 | 0 | 1 | 54 |
| 11 | SATB2 | ENSG00000119042 | 0.077 | 0 | 37 | 41 |
| 12 | DHX30 | ENSG00000132153 | 0.080 | 0 | 4 | 54 |
| 13 | ZNF462 | ENSG00000148143 | 0.082 | 0 | 8 | 49 |
| 14 | KAT6A | ENSG00000083168 | 0.083 | 0 | 29 | 54 |
| 15 | USP7 | ENSG00000187555 | 0.085 | 0 | 5 | 54 |
| 16 | PRPF3 | ENSG00000117360 | 0.087 | 0 | 4 | 54 |
| 17 | COL5A1 | ENSG00000130635 | 0.089 | 0 | 55 | 49 |
| 18 | BPTF | ENSG00000171634 | 0.090 | 0 | 7 | 54 |
| 19 | NFIX | ENSG00000008441 | 0.092 | 0 | 30 | 53 |
| 20 | NR4A2 | ENSG00000153234 | 0.093 | 0 | 2 | 54 |
| 21 | SPEN | ENSG00000065526 | 0.094 | 0 | 18 | 54 |
| 22 | TCF4 | ENSG00000196628 | 0.096 | 0 | 59 | 52 |
| 23 | RPL5 | ENSG00000122406 | 0.096 | 0 | 16 | 54 |
| 24 | KAT6B | ENSG00000156650 | 0.098 | 0 | 26 | 53 |
| 25 | COPA | ENSG00000122218 | 0.098 | 0 | 2 | 54 |
| 26 | STXBP1 | ENSG00000136854 | 0.099 | 0 | 89 | 52 |
| 27 | PAFAH1B1 | ENSG00000007168 | 0.100 | 0 | 32 | 54 |
| 28 | CREBBP | ENSG00000005339 | 0.101 | 0 | 77 | 54 |
| 29 | AGO2 | ENSG00000123908 | 0.102 | 0 | 5 | 54 |
| 30 | TLK2 | ENSG00000146872 | 0.103 | 0 | 21 | 54 |
| 31 | ASH1L | ENSG00000116539 | 0.104 | 0 | 15 | 54 |
| 32 | FUS | ENSG00000089280 | 0.104 | 0 | 14 | 54 |
| 33 | CLTC | ENSG00000141367 | 0.104 | 0 | 13 | 54 |
| 34 | FBN1 | ENSG00000166147 | 0.105 | 0 | 1,068 | 51 |
| 35 | GNAO1 | ENSG00000087258 | 0.105 | 0 | 27 | 44 |
| 36 | MAP1B | ENSG00000131711 | 0.105 | 0 | 7 | 52 |
| 37 | ADNP | ENSG00000101126 | 0.106 | 0 | 32 | 54 |
| 38 | COL5A2 | ENSG00000204262 | 0.106 | 0 | 6 | 50 |
| 39 | CHD7 | ENSG00000171316 | 0.109 | 0 | 158 | 54 |
| 40 | STAT3 | ENSG00000168610 | 0.109 | 0 | 27 | 54 |
| 41 | SETD1B | ENSG00000139718 | 0.109 | 0 | 13 | 54 |
| 42 | HNRNPU | ENSG00000153187 | 0.111 | 0 | 26 | 54 |
| 43 | AP1G1 | ENSG00000166747 | 0.113 | 0 | 2 | 54 |
| 44 | HNRNPH1 | ENSG00000169045 | 0.113 | 0 | 1 | 54 |
| 45 | ZEB2 | ENSG00000169554 | 0.114 | 0 | 67 | 52 |
| 46 | TRPS1 | ENSG00000104447 | 0.114 | 0 | 24 | 51 |
| 47 | HIVEP2 | ENSG00000010818 | 0.114 | 0 | 6 | 54 |
| 48 | AFF4 | ENSG00000072364 | 0.114 | 0 | 4 | 54 |
| 49 | TRRAP | ENSG00000196367 | 0.115 | 0 | 6 | 54 |
| 50 | DYNC1H1 | ENSG00000197102 | 0.117 | 0 | 41 | 54 |
| 51 | TAB2 | ENSG00000055208 | 0.117 | 0 | 13 | 54 |
| 52 | KDM3B | ENSG00000120733 | 0.117 | 0 | 8 | 54 |
| 53 | ANKRD17 | ENSG00000132466 | 0.118 | 0 | 3 | 54 |
| 54 | MATR3 | ENSG00000015479 | 0.122 | 0 | 1 | 54 |
| 55 | NTRK2 | ENSG00000148053 | 0.123 | 0 | 2 | 45 |
| 56 | BCL11A | ENSG00000119866 | 0.125 | 0 | 10 | 29 |
| 57 | UPF1 | ENSG00000005007 | 0.125 | 0 | 1 | 54 |
| 58 | SYNGAP1 | ENSG00000197283 | 0.127 | 0 | 91 | 54 |
| 59 | KMT2B | ENSG00000272333 | 0.127 | 0 | 26 | 54 |
| 60 | ZBTB18 | ENSG00000179456 | 0.128 | 0 | 16 | 54 |
| 61 | DDX6 | ENSG00000110367 | 0.128 | 0 | 1 | 54 |
| 62 | RAI1 | ENSG00000108557 | 0.129 | 0 | 17 | 54 |
| 63 | MED13 | ENSG00000108510 | 0.130 | 0 | 5 | 54 |
| 64 | ANKRD11 | ENSG00000167522 | 0.132 | 0 | 163 | 54 |
| 65 | SMC3 | ENSG00000108055 | 0.132 | 0 | 8 | 54 |
| 66 | NR2F1 | ENSG00000175745 | 0.133 | 0 | 21 | 50 |
| 67 | UBTF | ENSG00000108312 | 0.133 | 0 | 1 | 54 |
| 68 | VCP | ENSG00000165280 | 0.134 | 0 | 13 | 54 |
| 69 | NBEA | ENSG00000172915 | 0.137 | 0 | 8 | 51 |
| 70 | PURA | ENSG00000185129 | 0.138 | 0 | 43 | 53 |
| 71 | SLC2A1 | ENSG00000117394 | 0.139 | 0 | 86 | 54 |
| 72 | UBE3A | ENSG00000114062 | 0.139 | 0 | 48 | 54 |
| 73 | CTNNB1 | ENSG00000168036 | 0.139 | 0 | 35 | 54 |
| 74 | NSD2 | ENSG00000109685 | 0.139 | 0 | 7 | 54 |
| 75 | SIN3A | ENSG00000169375 | 0.140 | 0 | 13 | 54 |
| 76 | NCKAP1 | ENSG00000061676 | 0.140 | 0 | 1 | 53 |
| 77 | KMT2D | ENSG00000167548 | 0.141 | 0 | 159 | 54 |
| 78 | ARID1A | ENSG00000117713 | 0.141 | 0 | 16 | 54 |
| 79 | SETBP1 | ENSG00000152217 | 0.143 | 0 | 21 | 52 |
| 80 | NSD1 | ENSG00000165671 | 0.144 | 0 | 140 | 54 |
| 81 | WDR26 | ENSG00000162923 | 0.144 | 0 | 2 | 54 |
| 82 | EFTUD2 | ENSG00000108883 | 0.145 | 0 | 30 | 54 |
| 83 | TRIP12 | ENSG00000153827 | 0.145 | 0 | 11 | 54 |
| 84 | ARCN1 | ENSG00000095139 | 0.145 | 0 | 5 | 54 |
| 85 | PTCH1 | ENSG00000185920 | 0.146 | 0 | 149 | 50 |
| 86 | ZSWIM6 | ENSG00000130449 | 0.146 | 0 | 2 | 54 |
| 87 | PHF21A | ENSG00000135365 | 0.149 | 0 | 9 | 54 |
| 88 | DNMT1 | ENSG00000130816 | 0.150 | 0 | 3 | 54 |
| 89 | DLG4 | ENSG00000132535 | 0.151 | 0 | 50 | 54 |
| 90 | NOTCH2 | ENSG00000134250 | 0.151 | 0 | 11 | 54 |
| 91 | SLC6A1 | ENSG00000157103 | 0.153 | 0 | 38 | 44 |
| 92 | CHD8 | ENSG00000100888 | 0.154 | 0 | 33 | 54 |
| 93 | SOX9 | ENSG00000125398 | 0.154 | 0 | 19 | 49 |
| 94 | PTPN11 | ENSG00000179295 | 0.156 | 0 | 88 | 54 |
| 95 | PCBP1 | ENSG00000169564 | 0.156 | 0 | 1 | 54 |
| 96 | PRKAR1A | ENSG00000108946 | 0.158 | 0 | 25 | 54 |
| 97 | PPP1R12A | ENSG00000058272 | 0.158 | 0 | 3 | 54 |
| 98 | COL1A1 | ENSG00000108821 | 0.159 | 0 | 398 | 48 |
| 99 | EHMT1 | ENSG00000181090 | 0.159 | 0 | 48 | 54 |
| 100 | MYH9 | ENSG00000100345 | 0.159 | 0 | 23 | 54 |
| 101 | STAG1 | ENSG00000118007 | 0.159 | 0 | 8 | 54 |
| 102 | SPTBN1 | ENSG00000115306 | 0.159 | 0 | 6 | 54 |
| 103 | UBAP2L | ENSG00000143569 | 0.159 | 0 | 2 | 54 |
| 104 | CHAMP1 | ENSG00000198824 | 0.160 | 0 | 6 | 54 |
| 105 | ANK3 | ENSG00000151150 | 0.160 | 0 | 1 | 46 |
| 106 | ZNF292 | ENSG00000188994 | 0.161 | 0 | 18 | 54 |
| 107 | MBD5 | ENSG00000204406 | 0.161 | 0 | 13 | 45 |
| 108 | MECOM | ENSG00000085276 | 0.161 | 0 | 3 | 48 |
| 109 | EP300 | ENSG00000100393 | 0.164 | 0 | 27 | 54 |
| 110 | AGO1 | ENSG00000092847 | 0.164 | 0 | 4 | 54 |
| 111 | SOX6 | ENSG00000110693 | 0.165 | 0 | 5 | 40 |
| 112 | PUM1 | ENSG00000134644 | 0.165 | 0 | 3 | 54 |
| 113 | ZMYND11 | ENSG00000015171 | 0.166 | 0 | 17 | 54 |
| 114 | SMARCC2 | ENSG00000139613 | 0.167 | 0 | 13 | 54 |
| 115 | SET | ENSG00000119335 | 0.167 | 0 | 5 | 54 |
| 116 | AHDC1 | ENSG00000126705 | 0.168 | 0 | 19 | 54 |
| 117 | MSL2 | ENSG00000174579 | 0.169 | 0 | 1 | 54 |
| 118 | SRSF1 | ENSG00000136450 | 0.170 | 0 | 2 | 54 |
| 119 | TRIM28 | ENSG00000130726 | 0.170 | 0 | 2 | 54 |
| 120 | TANC2 | ENSG00000170921 | 0.170 | 0 | 2 | 51 |
| 121 | RPS17 | ENSG00000182774 | 0.170 | 0 | 1 | 54 |
| 122 | DYRK1A | ENSG00000157540 | 0.172 | 0 | 60 | 54 |
| 123 | ATP2B1 | ENSG00000070961 | 0.172 | 0 | 2 | 54 |
| 124 | CHD2 | ENSG00000173575 | 0.174 | 0 | 61 | 54 |
| 125 | TBL1XR1 | ENSG00000177565 | 0.174 | 0 | 16 | 54 |
| 126 | RNF220 | ENSG00000187147 | 0.174 | 0 | 1 | 54 |
| 127 | NFIA | ENSG00000162599 | 0.175 | 0 | 7 | 52 |
| 128 | SETD2 | ENSG00000181555 | 0.176 | 0 | 7 | 54 |
| 129 | NFKB1 | ENSG00000109320 | 0.178 | 0 | 20 | 54 |
| 130 | TUBB | ENSG00000196230 | 0.178 | 0 | 9 | 54 |
| 131 | PRR12 | ENSG00000126464 | 0.180 | 0 | 8 | 54 |
| 132 | FBXW7 | ENSG00000109670 | 0.181 | 0 | 2 | 54 |
| 133 | NUS1 | ENSG00000153989 | 0.182 | 0 | 9 | 54 |
| 134 | RPS24 | ENSG00000138326 | 0.182 | 0 | 2 | 54 |
| 135 | AP2M1 | ENSG00000161203 | 0.182 | 0 | 1 | 54 |
| 136 | SRRM2 | ENSG00000167978 | 0.183 | 0 | 8 | 54 |
| 137 | HSPD1 | ENSG00000144381 | 0.185 | 0 | 1 | 54 |
| 138 | CACNA1C | ENSG00000151067 | 0.186 | 0 | 25 | 48 |
| 139 | NF2 | ENSG00000186575 | 0.187 | 0 | 55 | 54 |
| 140 | CAMK2A | ENSG00000070808 | 0.187 | 0 | 4 | 28 |
| 141 | CTNND1 | ENSG00000198561 | 0.187 | 0 | 2 | 54 |
| 142 | ATP1A1 | ENSG00000163399 | 0.188 | 0 | 8 | 54 |
| 143 | ZFHX3 | ENSG00000140836 | 0.189 | 0 | 2 | 50 |
| 144 | KMT5B | ENSG00000110066 | 0.190 | 0 | 14 | 54 |
| 145 | ABL1 | ENSG00000097007 | 0.190 | 0 | 10 | 54 |
| 146 | BRD4 | ENSG00000141867 | 0.191 | 0 | 2 | 54 |
| 147 | CHD4 | ENSG00000111642 | 0.192 | 0 | 15 | 54 |
| 148 | HDAC4 | ENSG00000068024 | 0.192 | 0 | 2 | 54 |
| 149 | ACTB | ENSG00000075624 | 0.193 | 0 | 31 | 54 |
| 150 | NR6A1 | ENSG00000148200 | 0.193 | 0 | 1 | 37 |
| 151 | BMPR1A | ENSG00000107779 | 0.195 | 0 | 99 | 53 |
| 152 | DDB1 | ENSG00000167986 | 0.195 | 0 | 1 | 54 |
| 153 | NR2F2 | ENSG00000185551 | 0.197 | 0 | 3 | 53 |
| 154 | DLL4 | ENSG00000128917 | 0.197 | 0 | 2 | 51 |
| 155 | RAC2 | ENSG00000128340 | 0.197 | 0 | 1 | 52 |
| 156 | POGZ | ENSG00000143442 | 0.198 | 0 | 35 | 54 |
| 157 | PPP2R5C | ENSG00000078304 | 0.198 | 0 | 1 | 54 |
| 158 | WDFY3 | ENSG00000163625 | 0.200 | 0 | 16 | 53 |
| 159 | YY1 | ENSG00000100811 | 0.200 | 0 | 6 | 54 |
| 160 | SHOC2 | ENSG00000108061 | 0.201 | 0 | 3 | 54 |
| 161 | KMT2C | ENSG00000055609 | 0.202 | 0 | 23 | 54 |
| 162 | POLR2A | ENSG00000181222 | 0.202 | 0 | 6 | 54 |
| 163 | RPL35A | ENSG00000182899 | 0.202 | 0 | 1 | 54 |
| 164 | TSC2 | ENSG00000103197 | 0.203 | 0 | 390 | 54 |
| 165 | SOX5 | ENSG00000134532 | 0.204 | 0 | 14 | 40 |
| 166 | YWHAG | ENSG00000170027 | 0.204 | 0 | 7 | 54 |
| 167 | ARID1B | ENSG00000049618 | 0.205 | 0 | 119 | 54 |
| 168 | CTCF | ENSG00000102974 | 0.205 | 0 | 17 | 54 |
| 169 | BRPF1 | ENSG00000156983 | 0.205 | 0 | 12 | 54 |
| 170 | ATP2A2 | ENSG00000174437 | 0.205 | 0 | 7 | 54 |
| 171 | TAF4 | ENSG00000130699 | 0.205 | 0 | 4 | 54 |
| 172 | DNAJC5 | ENSG00000101152 | 0.205 | 0 | 2 | 54 |
| 173 | MEN1 | ENSG00000133895 | 0.206 | 0 | 266 | 54 |
| 174 | RORA | ENSG00000069667 | 0.207 | 0 | 2 | 54 |
| 175 | ZNF148 | ENSG00000163848 | 0.208 | 0 | 4 | 54 |
| 176 | TNFAIP3 | ENSG00000118503 | 0.209 | 0 | 8 | 52 |
| 177 | SATB1 | ENSG00000182568 | 0.209 | 0 | 7 | 53 |
| 178 | U2AF2 | ENSG00000063244 | 0.209 | 0 | 5 | 54 |
| 179 | JAG1 | ENSG00000101384 | 0.211 | 0 | 95 | 53 |
| 180 | KDM6B | ENSG00000132510 | 0.211 | 0 | 13 | 54 |
| 181 | MEF2C | ENSG00000081189 | 0.212 | 0 | 22 | 52 |
| 182 | PBX1 | ENSG00000185630 | 0.212 | 0 | 10 | 52 |
| 183 | BMP2 | ENSG00000125845 | 0.213 | 0 | 8 | 48 |
| 184 | TET3 | ENSG00000187605 | 0.213 | 0 | 3 | 52 |
| 185 | SUFU | ENSG00000107882 | 0.214 | 0 | 24 | 54 |
| 186 | IRF6 | ENSG00000117595 | 0.214 | 0 | 21 | 38 |
| 187 | FBXO11 | ENSG00000138081 | 0.215 | 0 | 20 | 54 |
| 188 | MYRF | ENSG00000124920 | 0.215 | 0 | 4 | 40 |
| 189 | TAOK1 | ENSG00000160551 | 0.216 | 0 | 9 | 54 |
| 190 | TRIO | ENSG00000038382 | 0.217 | 0 | 26 | 54 |
| 191 | KMT2E | ENSG00000005483 | 0.217 | 0 | 21 | 54 |
| 192 | NOTCH1 | ENSG00000148400 | 0.217 | 0 | 17 | 54 |
| 193 | KLF7 | ENSG00000118263 | 0.219 | 0 | 1 | 54 |
| 194 | ZMIZ1 | ENSG00000108175 | 0.220 | 0 | 9 | 54 |
| 195 | TCP1 | ENSG00000120438 | 0.220 | 0 | 1 | 54 |
| 196 | SMARCB1 | ENSG00000099956 | 0.221 | 0 | 25 | 54 |
| 197 | PIK3CA | ENSG00000121879 | 0.221 | 0 | 8 | 54 |
| 198 | SPTAN1 | ENSG00000197694 | 0.222 | 0 | 21 | 54 |
| 199 | GNB1 | ENSG00000078369 | 0.222 | 0 | 18 | 54 |
| 200 | SMARCA5 | ENSG00000153147 | 0.222 | 0 | 3 | 54 |
| 201 | HNRNPD | ENSG00000138668 | 0.222 | 0 | 1 | 54 |
| 202 | SRCAP | ENSG00000080603 | 0.223 | 0 | 15 | 54 |
| 203 | VCAN | ENSG00000038427 | 0.223 | 0 | 5 | 51 |
| 204 | ARID2 | ENSG00000189079 | 0.224 | 0 | 8 | 54 |
| 205 | SMARCA4 | ENSG00000127616 | 0.226 | 0 | 64 | 54 |
| 206 | TSC1 | ENSG00000165699 | 0.227 | 0 | 198 | 54 |
| 207 | SMAD4 | ENSG00000141646 | 0.227 | 0 | 86 | 54 |
| 208 | CAMTA1 | ENSG00000171735 | 0.227 | 0 | 17 | 54 |
| 209 | WNK1 | ENSG00000060237 | 0.227 | 0 | 8 | 54 |
| 210 | SPTB | ENSG00000070182 | 0.228 | 0 | 101 | 31 |
| 211 | SALL1 | ENSG00000103449 | 0.228 | 0 | 15 | 31 |
| 212 | GNAI2 | ENSG00000114353 | 0.228 | 0 | 1 | 54 |
| 213 | TARDBP | ENSG00000120948 | 0.229 | 0 | 10 | 54 |
| 214 | PIK3CD | ENSG00000171608 | 0.230 | 0 | 12 | 53 |
| 215 | RANBP2 | ENSG00000153201 | 0.231 | 0 | 1 | 54 |
| 216 | AUTS2 | ENSG00000158321 | 0.232 | 0 | 19 | 53 |
| 217 | PSMD12 | ENSG00000197170 | 0.232 | 0 | 5 | 54 |
| 218 | PPP2R5D | ENSG00000112640 | 0.234 | 0 | 13 | 54 |
| 219 | CDC42BPB | ENSG00000198752 | 0.235 | 0 | 9 | 54 |
| 220 | GATAD2B | ENSG00000143614 | 0.236 | 0 | 24 | 54 |
| 221 | ARFGEF1 | ENSG00000066777 | 0.236 | 0 | 8 | 54 |
| 222 | DLL1 | ENSG00000198719 | 0.236 | 0 | 5 | 53 |
| 223 | BRAF | ENSG00000157764 | 0.237 | 0 | 63 | 54 |
| 224 | CSNK2B | ENSG00000204435 | 0.237 | 0 | 41 | 54 |
| 225 | CLASP1 | ENSG00000074054 | 0.237 | 0 | 8 | 54 |
| 226 | COL4A1 | ENSG00000187498 | 0.238 | 0 | 72 | 52 |
| 227 | ACTN4 | ENSG00000130402 | 0.238 | 0 | 1 | 54 |
| 228 | MACF1 | ENSG00000127603 | 0.238 | 0 | 1 | 54 |
| 229 | RB1 | ENSG00000139687 | 0.239 | 0 | 265 | 54 |
| 230 | RPS19 | ENSG00000105372 | 0.239 | 0 | 24 | 54 |
| 231 | DNM2 | ENSG00000079805 | 0.239 | 0 | 20 | 54 |
| 232 | GPHN | ENSG00000171723 | 0.240 | 0 | 4 | 53 |
| 233 | KIT | ENSG00000157404 | 0.240 | 0 | 4 | 46 |
| 234 | HNRNPK | ENSG00000165119 | 0.242 | 0 | 11 | 54 |
| 235 | RAB7A | ENSG00000075785 | 0.242 | 0 | 3 | 54 |
| 236 | RRM1 | ENSG00000167325 | 0.242 | 0 | 1 | 54 |
| 237 | PPP3CA | ENSG00000138814 | 0.243 | 0 | 12 | 54 |
| 238 | HERC2 | ENSG00000128731 | 0.243 | 0 | 3 | 54 |
| 239 | PHIP | ENSG00000146247 | 0.246 | 0 | 34 | 54 |
| 240 | PPP1CB | ENSG00000213639 | 0.247 | 0 | 7 | 54 |
| 241 | ARF1 | ENSG00000143761 | 0.247 | 0 | 5 | 54 |
| 242 | WAC | ENSG00000095787 | 0.248 | 0 | 20 | 54 |
| 243 | TGFB3 | ENSG00000119699 | 0.249 | 0 | 22 | 53 |
| 244 | CIC | ENSG00000079432 | 0.249 | 0 | 7 | 54 |
| 245 | SMAD2 | ENSG00000175387 | 0.249 | 0 | 7 | 54 |
| 246 | KIF5C | ENSG00000168280 | 0.249 | 0 | 5 | 42 |
| 247 | PTBP1 | ENSG00000011304 | 0.249 | 0 | 1 | 54 |
| 248 | SYNCRIP | ENSG00000135316 | 0.249 | 0 | 1 | 54 |
| 249 | CDH11 | ENSG00000140937 | 0.249 | 0 | 1 | 48 |
| 250 | RAC1 | ENSG00000136238 | 0.250 | 0 | 8 | 54 |
| 251 | RPL13 | ENSG00000167526 | 0.250 | 0 | 1 | 54 |
| 252 | SEPHS1 | ENSG00000086475 | 0.254 | 0 | 1 | 54 |
| 253 | MEIS2 | ENSG00000134138 | 0.255 | 0 | 4 | 51 |
| 254 | CHD3 | ENSG00000170004 | 0.257 | 0 | 25 | 54 |
| 255 | SMARCAD1 | ENSG00000163104 | 0.257 | 0 | 1 | 54 |
| 256 | UBA2 | ENSG00000126261 | 0.258 | 0 | 3 | 54 |
| 257 | RPL26 | ENSG00000161970 | 0.258 | 0 | 1 | 54 |
| 258 | BMP4 | ENSG00000125378 | 0.258 | 0 | 1 | 42 |
| 259 | SF3B4 | ENSG00000143368 | 0.259 | 0 | 6 | 54 |
| 260 | PRPF8 | ENSG00000174231 | 0.260 | 0 | 23 | 54 |
| 261 | QRICH1 | ENSG00000198218 | 0.261 | 0 | 17 | 54 |
| 262 | CALM2 | ENSG00000143933 | 0.262 | 0 | 4 | 54 |
| 263 | KRAS | ENSG00000133703 | 0.264 | 0 | 24 | 54 |
| 264 | PPFIA3 | ENSG00000177380 | 0.264 | 0 | 3 | 51 |
| 265 | FLT4 | ENSG00000037280 | 0.265 | 0 | 10 | 50 |
| 266 | RPS10 | ENSG00000124614 | 0.265 | 0 | 2 | 54 |
| 267 | KDM2B | ENSG00000089094 | 0.266 | 0 | 3 | 54 |
| 268 | FGF10 | ENSG00000070193 | 0.267 | 0 | 3 | 27 |
| 269 | RBPJ | ENSG00000168214 | 0.267 | 0 | 1 | 54 |
| 270 | NR3C2 | ENSG00000151623 | 0.268 | 0 | 10 | 51 |
| 271 | MTOR | ENSG00000198793 | 0.268 | 0 | 5 | 54 |
| 272 | LRP6 | ENSG00000070018 | 0.270 | 0 | 7 | 53 |
| 273 | KDM5A | ENSG00000073614 | 0.270 | 0 | 1 | 54 |
| 274 | STAT1 | ENSG00000115415 | 0.271 | 0 | 15 | 54 |
| 275 | ANK2 | ENSG00000145362 | 0.271 | 0 | 11 | 51 |
| 276 | EZH2 | ENSG00000106462 | 0.272 | 0 | 7 | 43 |
| 277 | FGFR1 | ENSG00000077782 | 0.274 | 0 | 44 | 54 |
| 278 | CSNK2A1 | ENSG00000101266 | 0.274 | 0 | 23 | 54 |
| 279 | CYFIP2 | ENSG00000055163 | 0.274 | 0 | 4 | 54 |
| 280 | ZFPM2 | ENSG00000169946 | 0.274 | 0 | 2 | 45 |
| 281 | DICER1 | ENSG00000100697 | 0.275 | 0 | 271 | 54 |
| 282 | RELN | ENSG00000189056 | 0.276 | 0 | 8 | 35 |
| 283 | CALM1 | ENSG00000198668 | 0.276 | 0 | 4 | 54 |
| 284 | CUX1 | ENSG00000257923 | 0.277 | 0 | 8 | 54 |
| 285 | STAT5B | ENSG00000173757 | 0.277 | 0 | 3 | 54 |
| 286 | FGFR2 | ENSG00000066468 | 0.278 | 0 | 61 | 47 |
| 287 | SERPING1 | ENSG00000149131 | 0.278 | 0 | 46 | 54 |
| 288 | EIF5A | ENSG00000132507 | 0.278 | 0 | 2 | 54 |
| 289 | FOXP1 | ENSG00000114861 | 0.279 | 0 | 37 | 53 |
| 290 | SETD1A | ENSG00000099381 | 0.279 | 0 | 10 | 54 |
| 291 | PRICKLE1 | ENSG00000139174 | 0.279 | 0 | 1 | 44 |
| 292 | SCAF4 | ENSG00000156304 | 0.285 | 0 | 8 | 54 |
| 293 | UBAP1 | ENSG00000165006 | 0.285 | 0 | 1 | 54 |
| 294 | DPF2 | ENSG00000133884 | 0.287 | 0 | 4 | 54 |
| 295 | MARK2 | ENSG00000072518 | 0.287 | 0 | 2 | 54 |
| 296 | TCOF1 | ENSG00000070814 | 0.288 | 0 | 37 | 54 |
| 297 | SOS1 | ENSG00000115904 | 0.288 | 0 | 34 | 54 |
| 298 | MAX | ENSG00000125952 | 0.288 | 0 | 11 | 54 |
| 299 | NPTN | ENSG00000156642 | 0.288 | 0 | 2 | 54 |
| 300 | IGSF3 | ENSG00000143061 | 0.289 | 0 | 1 | 47 |
| 301 | SNAP25 | ENSG00000132639 | 0.290 | 0 | 9 | 41 |
| 302 | COL3A1 | ENSG00000168542 | 0.291 | 0 | 219 | 43 |
| 303 | TRIM8 | ENSG00000171206 | 0.291 | 0 | 4 | 54 |
| 304 | CELF2 | ENSG00000048740 | 0.294 | 0 | 3 | 54 |
| 305 | PRPF4 | ENSG00000136875 | 0.294 | 0 | 1 | 54 |
| 306 | PUF60 | ENSG00000179950 | 0.295 | 0 | 24 | 54 |
| 307 | DST | ENSG00000151914 | 0.297 | 0 | 18 | 53 |
| 308 | GLI3 | ENSG00000106571 | 0.298 | 0 | 20 | 44 |
| 309 | DAG1 | ENSG00000173402 | 0.298 | 0 | 2 | 54 |
| 310 | EIF4A2 | ENSG00000156976 | 0.298 | 0 | 2 | 54 |
| 311 | SIAH1 | ENSG00000196470 | 0.298 | 0 | 2 | 54 |
| 312 | DOCK3 | ENSG00000088538 | 0.298 | 0 | 1 | 43 |
| 313 | RASA1 | ENSG00000145715 | 0.300 | 0 | 43 | 54 |
| 314 | MORC2 | ENSG00000133422 | 0.300 | 0 | 11 | 54 |
| 315 | NFIB | ENSG00000147862 | 0.300 | 0 | 10 | 52 |
| 316 | RFX7 | ENSG00000181827 | 0.300 | 0 | 3 | 53 |
| 317 | PACS1 | ENSG00000175115 | 0.300 | 0 | 1 | 54 |
| 318 | PHF3 | ENSG00000118482 | 0.301 | 0 | 1 | 54 |
| 319 | WASF1 | ENSG00000112290 | 0.301 | 0 | 1 | 51 |
| 320 | GNAS | ENSG00000087460 | 0.302 | 0 | 62 | 54 |
| 321 | SOX2 | ENSG00000181449 | 0.302 | 0 | 14 | 31 |
| 322 | KCND3 | ENSG00000171385 | 0.303 | 0 | 6 | 50 |
| 323 | KIF1A | ENSG00000130294 | 0.304 | 0 | 43 | 27 |
| 324 | SRP54 | ENSG00000100883 | 0.304 | 0 | 3 | 54 |
| 325 | PRDM16 | ENSG00000142611 | 0.304 | 0 | 2 | 40 |
| 326 | ITSN1 | ENSG00000205726 | 0.304 | 0 | 1 | 51 |
| 327 | STX1B | ENSG00000099365 | 0.308 | 0 | 6 | 49 |
| 328 | NFKB2 | ENSG00000077150 | 0.308 | 0 | 5 | 54 |
| 329 | MAP4K4 | ENSG00000071054 | 0.308 | 0 | 1 | 54 |
| 330 | RPS29 | ENSG00000213741 | 0.308 | 0 | 1 | 54 |
| 331 | KCNB1 | ENSG00000158445 | 0.309 | 0 | 38 | 31 |
| 332 | SNRNP200 | ENSG00000144028 | 0.312 | 0 | 7 | 54 |
| 333 | CNNM2 | ENSG00000148842 | 0.312 | 0 | 6 | 50 |
| 334 | EEF1A2 | ENSG00000101210 | 0.314 | 0 | 13 | 48 |
| 335 | TMEM63B | ENSG00000137216 | 0.314 | 0 | 4 | 54 |
| 336 | KIF5A | ENSG00000155980 | 0.315 | 0 | 17 | 34 |
| 337 | PPP2CA | ENSG00000113575 | 0.315 | 0 | 7 | 54 |
| 338 | SMARCA2 | ENSG00000080503 | 0.317 | 0 | 28 | 54 |
| 339 | FAT4 | ENSG00000196159 | 0.317 | 0 | 5 | 37 |
| 340 | RPS7 | ENSG00000171863 | 0.317 | 0 | 2 | 54 |
| 341 | SMARCE1 | ENSG00000073584 | 0.318 | 0 | 17 | 54 |
| 342 | AFF3 | ENSG00000144218 | 0.318 | 0 | 2 | 43 |
| 343 | ATXN2 | ENSG00000204842 | 0.319 | 0 | 1 | 54 |
| 344 | RPS26 | ENSG00000197728 | 0.320 | 0 | 11 | 54 |
| 345 | IKZF1 | ENSG00000185811 | 0.320 | 0 | 4 | 31 |
| 346 | VPS35 | ENSG00000069329 | 0.320 | 0 | 1 | 54 |
| 347 | APC | ENSG00000134982 | 0.321 | 0 | 1,055 | 54 |
| 348 | WBP11 | ENSG00000084463 | 0.321 | 0 | 3 | 54 |
| 349 | TRAF3 | ENSG00000131323 | 0.321 | 0 | 1 | 54 |
| 350 | BAP1 | ENSG00000163930 | 0.324 | 0 | 155 | 54 |
| 351 | PFN1 | ENSG00000108518 | 0.324 | 0 | 1 | 54 |
| 352 | SNRPB | ENSG00000125835 | 0.324 | 0 | 1 | 54 |
| 353 | ITPR1 | ENSG00000150995 | 0.325 | 0 | 29 | 54 |
| 354 | KANSL1 | ENSG00000120071 | 0.325 | 0 | 21 | 54 |
| 355 | CYLD | ENSG00000083799 | 0.325 | 0 | 13 | 54 |
| 356 | ADCY5 | ENSG00000173175 | 0.326 | 0 | 10 | 51 |
| 357 | TTBK2 | ENSG00000128881 | 0.326 | 0 | 1 | 54 |
| 358 | NACC1 | ENSG00000160877 | 0.329 | 0 | 1 | 54 |
| 359 | HERC1 | ENSG00000103657 | 0.330 | 0 | 3 | 54 |
| 360 | ACTN1 | ENSG00000072110 | 0.331 | 0 | 6 | 54 |
| 361 | KCNMA1 | ENSG00000156113 | 0.331 | 0 | 5 | 44 |
| 362 | BMPR1B | ENSG00000138696 | 0.331 | 0 | 2 | 35 |
| 363 | TNRC6B | ENSG00000100354 | 0.332 | 0 | 8 | 54 |
| 364 | SYT1 | ENSG00000067715 | 0.332 | 0 | 7 | 29 |
| 365 | GATA2 | ENSG00000179348 | 0.335 | 0 | 59 | 50 |
| 366 | CYCS | ENSG00000172115 | 0.335 | 0 | 1 | 54 |
| 367 | CDH2 | ENSG00000170558 | 0.336 | 0 | 4 | 41 |
| 368 | TUBB2B | ENSG00000137285 | 0.338 | 0 | 13 | 47 |
| 369 | PRKACA | ENSG00000072062 | 0.338 | 0 | 1 | 54 |
| 370 | COL1A2 | ENSG00000164692 | 0.339 | 0 | 238 | 53 |
| 371 | EXT1 | ENSG00000182197 | 0.339 | 0 | 122 | 54 |
| 372 | ZEB1 | ENSG00000148516 | 0.339 | 0 | 5 | 54 |
| 373 | RAD21 | ENSG00000164754 | 0.341 | 0 | 8 | 54 |
| 374 | CLCN3 | ENSG00000109572 | 0.341 | 0 | 3 | 54 |
| 375 | SMCHD1 | ENSG00000101596 | 0.342 | 0 | 19 | 54 |
| 376 | RPL11 | ENSG00000142676 | 0.342 | 0 | 10 | 54 |
| 377 | PTPRC | ENSG00000081237 | 0.343 | 0 | 2 | 48 |
| 378 | ITCH | ENSG00000078747 | 0.343 | 0 | 1 | 54 |
| 379 | SUZ12 | ENSG00000178691 | 0.344 | 0 | 5 | 54 |
| 380 | UCHL1 | ENSG00000154277 | 0.344 | 0 | 3 | 52 |
| 381 | PPP2R1A | ENSG00000105568 | 0.345 | 0 | 12 | 54 |
| 382 | ANTXR1 | ENSG00000169604 | 0.346 | 0 | 1 | 52 |
| 383 | PRPF31 | ENSG00000105618 | 0.347 | 0 | 78 | 54 |
| 384 | DENND5B | ENSG00000170456 | 0.351 | 0 | 1 | 52 |
| 385 | RUNX2 | ENSG00000124813 | 0.352 | 0 | 25 | 31 |
| 386 | ATN1 | ENSG00000111676 | 0.352 | 0 | 4 | 54 |
| 387 | NUP85 | ENSG00000125450 | 0.352 | 0 | 1 | 54 |
| 388 | PSMC3 | ENSG00000165916 | 0.352 | 0 | 1 | 54 |
| 389 | RAF1 | ENSG00000132155 | 0.353 | 0 | 33 | 54 |
| 390 | EBF3 | ENSG00000108001 | 0.353 | 0 | 24 | 35 |
| 391 | CIT | ENSG00000122966 | 0.353 | 0 | 2 | 45 |
| 392 | HTT | ENSG00000197386 | 0.353 | 0 | 1 | 54 |
| 393 | FBXW11 | ENSG00000072803 | 0.354 | 0 | 5 | 54 |
| 394 | MAP3K1 | ENSG00000095015 | 0.355 | 0 | 1 | 54 |
| 395 | PLCB1 | ENSG00000182621 | 0.357 | 0 | 2 | 41 |
| 396 | DMXL2 | ENSG00000104093 | 0.360 | 0 | 3 | 53 |
| 397 | BRSK2 | ENSG00000174672 | 0.361 | 0 | 2 | 36 |
| 398 | PIK3R1 | ENSG00000145675 | 0.363 | 0 | 22 | 54 |
| 399 | PKD1 | ENSG00000008710 | 0.364 | 0 | 636 | 54 |
| 400 | STK11 | ENSG00000118046 | 0.364 | 0 | 108 | 54 |
| 401 | GLI2 | ENSG00000074047 | 0.364 | 0 | 10 | 38 |
| 402 | CARD11 | ENSG00000198286 | 0.364 | 0 | 4 | 38 |
| 403 | CUL3 | ENSG00000036257 | 0.365 | 0 | 20 | 54 |
| 404 | EPAS1 | ENSG00000116016 | 0.366 | 0 | 2 | 54 |
| 405 | DIAPH1 | ENSG00000131504 | 0.367 | 0 | 10 | 54 |
| 406 | CTR9 | ENSG00000198730 | 0.368 | 0 | 3 | 54 |
| 407 | TEK | ENSG00000120156 | 0.368 | 0 | 2 | 52 |
| 408 | PABPN1 | ENSG00000100836 | 0.369 | 0 | 6 | 54 |
| 409 | SMARCC1 | ENSG00000173473 | 0.369 | 0 | 2 | 54 |
| 410 | TUBA1A | ENSG00000167552 | 0.370 | 0 | 82 | 54 |
| 411 | KDM4B | ENSG00000127663 | 0.371 | 0 | 3 | 54 |
| 412 | CYP26B1 | ENSG00000003137 | 0.371 | 0 | 1 | 46 |
| 413 | FAS | ENSG00000026103 | 0.372 | 0 | 22 | 52 |
| 414 | GNB2 | ENSG00000172354 | 0.372 | 0 | 3 | 54 |
| 415 | CAPRIN1 | ENSG00000135387 | 0.372 | 0 | 2 | 54 |
| 416 | IL6ST | ENSG00000134352 | 0.372 | 0 | 2 | 54 |
| 417 | SMAD3 | ENSG00000166949 | 0.373 | 0 | 74 | 54 |
| 418 | EFEMP1 | ENSG00000115380 | 0.373 | 0 | 2 | 52 |
| 419 | PLAG1 | ENSG00000181690 | 0.373 | 0 | 1 | 35 |
| 420 | ENG | ENSG00000106991 | 0.374 | 0 | 242 | 54 |
| 421 | PMP22 | ENSG00000109099 | 0.374 | 0 | 32 | 53 |
| 422 | COL12A1 | ENSG00000111799 | 0.375 | 0 | 25 | 43 |
| 423 | DDR2 | ENSG00000162733 | 0.375 | 0 | 4 | 53 |
| 424 | BMPR2 | ENSG00000204217 | 0.376 | 0 | 83 | 54 |
| 425 | SETD5 | ENSG00000168137 | 0.376 | 0 | 44 | 54 |
| 426 | MAPK1 | ENSG00000100030 | 0.377 | 0 | 6 | 54 |
| 427 | MYH11 | ENSG00000133392 | 0.377 | 0 | 6 | 54 |
| 428 | KLHL3 | ENSG00000146021 | 0.377 | 0 | 6 | 52 |
| 429 | SPOP | ENSG00000121067 | 0.377 | 0 | 3 | 54 |
| 430 | TUBB4B | ENSG00000188229 | 0.377 | 0 | 3 | 54 |
| 431 | GLUL | ENSG00000135821 | 0.377 | 0 | 1 | 54 |
| 432 | GABBR1 | ENSG00000204681 | 0.378 | 0 | 1 | 54 |
| 433 | KCNN3 | ENSG00000143603 | 0.379 | 0 | 1 | 34 |
| 434 | GNA11 | ENSG00000088256 | 0.380 | 0 | 2 | 54 |
| 435 | SBF1 | ENSG00000100241 | 0.381 | 0 | 3 | 54 |
| 436 | CADM3 | ENSG00000162706 | 0.381 | 0 | 1 | 47 |
| 437 | TNFRSF1A | ENSG00000067182 | 0.382 | 0 | 10 | 54 |
| 438 | SRC | ENSG00000197122 | 0.383 | 0 | 1 | 54 |
| 439 | RUSC2 | ENSG00000198853 | 0.383 | 0 | 1 | 53 |
| 440 | GPC6 | ENSG00000183098 | 0.384 | 0 | 1 | 38 |
| 441 | CCND2 | ENSG00000118971 | 0.385 | 0 | 9 | 54 |
| 442 | BAZ2B | ENSG00000123636 | 0.386 | 0 | 1 | 54 |
| 443 | PSMC5 | ENSG00000087191 | 0.386 | 0 | 1 | 54 |
| 444 | HK1 | ENSG00000156515 | 0.387 | 0 | 7 | 54 |
| 445 | SEPTIN9 | ENSG00000184640 | 0.390 | 0 | 1 | 54 |
| 446 | GANAB | ENSG00000089597 | 0.391 | 0 | 16 | 54 |
| 447 | ATP6V1A | ENSG00000114573 | 0.391 | 0 | 5 | 54 |
| 448 | IGF1R | ENSG00000140443 | 0.393 | 0 | 12 | 52 |
| 449 | ADD3 | ENSG00000148700 | 0.395 | 0 | 1 | 54 |
| 450 | DNMT3B | ENSG00000088305 | 0.396 | 0 | 10 | 35 |
| 451 | FTH1 | ENSG00000167996 | 0.396 | 0 | 7 | 54 |
| 452 | CACNA2D2 | ENSG00000007402 | 0.396 | 0 | 4 | 47 |
| 453 | AXIN2 | ENSG00000168646 | 0.397 | 0 | 60 | 52 |
| 454 | MAP3K7 | ENSG00000135341 | 0.397 | 0 | 4 | 54 |
| 455 | LMNB2 | ENSG00000176619 | 0.399 | 0 | 1 | 54 |
| 456 | RELA | ENSG00000173039 | 0.400 | 0 | 6 | 54 |
| 457 | EGR2 | ENSG00000122877 | 0.401 | 1 | 9 | 37 |
| 458 | APBB1 | ENSG00000166313 | 0.401 | 1 | 2 | 54 |
| 459 | ATP2C1 | ENSG00000017260 | 0.402 | 1 | 10 | 54 |
| 460 | BCL2L2-PABPN1 | ENSG00000258643 | 0.402 | 1 | 1 | 47 |
| 461 | DEPDC5 | ENSG00000100150 | 0.403 | 1 | 72 | 54 |
| 462 | PDE4D | ENSG00000113448 | 0.403 | 1 | 2 | 51 |
| 463 | MINK1 | ENSG00000141503 | 0.403 | 1 | 1 | 54 |
| 464 | FN1 | ENSG00000115414 | 0.404 | 1 | 5 | 54 |
| 465 | CAMK2G | ENSG00000148660 | 0.404 | 1 | 1 | 54 |
| 466 | SEC61A1 | ENSG00000058262 | 0.404 | 1 | 1 | 54 |
| 467 | MAPK8IP3 | ENSG00000138834 | 0.407 | 1 | 8 | 54 |
| 468 | GOLGA2 | ENSG00000167110 | 0.407 | 1 | 2 | 54 |
| 469 | AGPS | ENSG00000018510 | 0.408 | 1 | 3 | 54 |
| 470 | CAD | ENSG00000084774 | 0.409 | 1 | 9 | 54 |
| 471 | APP | ENSG00000142192 | 0.410 | 1 | 10 | 54 |
| 472 | EED | ENSG00000074266 | 0.410 | 1 | 2 | 54 |
| 473 | RAB5C | ENSG00000108774 | 0.410 | 1 | 1 | 54 |
| 474 | LMNA | ENSG00000160789 | 0.412 | 1 | 189 | 54 |
| 475 | SOS2 | ENSG00000100485 | 0.415 | 1 | 6 | 54 |
| 476 | GRIP1 | ENSG00000155974 | 0.415 | 1 | 6 | 28 |
| 477 | PAX8 | ENSG00000125618 | 0.415 | 1 | 3 | 33 |
| 478 | DNAJC6 | ENSG00000116675 | 0.416 | 1 | 1 | 40 |
| 479 | DNM1 | ENSG00000106976 | 0.417 | 1 | 16 | 52 |
| 480 | SHANK2 | ENSG00000162105 | 0.417 | 1 | 9 | 29 |
| 481 | TOPORS | ENSG00000197579 | 0.418 | 1 | 9 | 54 |
| 482 | CDC73 | ENSG00000134371 | 0.420 | 1 | 30 | 54 |
| 483 | PCNA | ENSG00000132646 | 0.420 | 1 | 1 | 54 |
| 484 | COX6B1 | ENSG00000126267 | 0.421 | 1 | 1 | 54 |
| 485 | WDR37 | ENSG00000047056 | 0.423 | 1 | 6 | 54 |
| 486 | CTPS1 | ENSG00000171793 | 0.423 | 1 | 1 | 54 |
| 487 | FLNC | ENSG00000128591 | 0.424 | 1 | 133 | 49 |
| 488 | RYR2 | ENSG00000198626 | 0.424 | 1 | 65 | 37 |
| 489 | PSEN1 | ENSG00000080815 | 0.425 | 1 | 50 | 54 |
| 490 | NSF | ENSG00000073969 | 0.425 | 1 | 1 | 54 |
| 491 | VPS13D | ENSG00000048707 | 0.426 | 1 | 9 | 54 |
| 492 | PDE10A | ENSG00000112541 | 0.427 | 1 | 4 | 35 |
| 493 | UBE3C | ENSG00000009335 | 0.427 | 1 | 1 | 54 |
| 494 | ATP13A3 | ENSG00000133657 | 0.428 | 1 | 2 | 54 |
| 495 | CFH | ENSG00000000971 | 0.430 | 1 | 31 | 52 |
| 496 | THRB | ENSG00000151090 | 0.430 | 1 | 31 | 51 |
| 497 | CNTN1 | ENSG00000018236 | 0.431 | 1 | 1 | 46 |
| 498 | ACTN2 | ENSG00000077522 | 0.432 | 1 | 3 | 34 |
| 499 | KIF5B | ENSG00000170759 | 0.433 | 1 | 3 | 54 |
| 500 | CDK13 | ENSG00000065883 | 0.434 | 1 | 21 | 54 |
| 501 | SPECC1L | ENSG00000100014 | 0.434 | 1 | 1 | 54 |
| 502 | TEAD3 | ENSG00000007866 | 0.435 | 1 | 1 | 52 |
| 503 | ABCA2 | ENSG00000107331 | 0.436 | 1 | 2 | 53 |
| 504 | FBXO28 | ENSG00000143756 | 0.437 | 1 | 1 | 54 |
| 505 | WNT5A | ENSG00000114251 | 0.437 | 1 | 1 | 41 |
| 506 | NAA15 | ENSG00000164134 | 0.438 | 1 | 25 | 54 |
| 507 | THRA | ENSG00000126351 | 0.439 | 1 | 3 | 54 |
| 508 | EZH1 | ENSG00000108799 | 0.439 | 1 | 2 | 54 |
| 509 | ETV6 | ENSG00000139083 | 0.440 | 1 | 10 | 53 |
| 510 | SETX | ENSG00000107290 | 0.442 | 1 | 33 | 54 |
| 511 | DSP | ENSG00000096696 | 0.443 | 1 | 271 | 44 |
| 512 | ERF | ENSG00000105722 | 0.443 | 1 | 16 | 54 |
| 513 | PCGF2 | ENSG00000277258 | 0.443 | 1 | 2 | 53 |
| 514 | MRAS | ENSG00000158186 | 0.444 | 1 | 4 | 52 |
| 515 | SLC40A1 | ENSG00000138449 | 0.445 | 1 | 9 | 53 |
| 516 | RSPRY1 | ENSG00000159579 | 0.445 | 1 | 1 | 54 |
| 517 | SHANK3 | ENSG00000251322 | 0.446 | 1 | 27 | 52 |
| 518 | MET | ENSG00000105976 | 0.446 | 1 | 4 | 41 |
| 519 | HMGCR | ENSG00000113161 | 0.447 | 1 | 1 | 54 |
| 520 | DCTN1 | ENSG00000204843 | 0.448 | 1 | 4 | 54 |
| 521 | LONP1 | ENSG00000196365 | 0.448 | 1 | 1 | 54 |
| 522 | NF1 | ENSG00000196712 | 0.449 | 1 | 1,778 | 54 |
| 523 | TP53 | ENSG00000141510 | 0.449 | 1 | 221 | 54 |
| 524 | ADAR | ENSG00000160710 | 0.449 | 1 | 14 | 54 |
| 525 | TUBB2A | ENSG00000137267 | 0.449 | 1 | 10 | 54 |
| 526 | SPTBN2 | ENSG00000173898 | 0.452 | 1 | 8 | 43 |
| 527 | KIDINS220 | ENSG00000134313 | 0.452 | 1 | 5 | 54 |
| 528 | WARS1 | ENSG00000140105 | 0.456 | 1 | 1 | 54 |
| 529 | ATP6V1B2 | ENSG00000147416 | 0.457 | 1 | 6 | 54 |
| 530 | NEDD4L | ENSG00000049759 | 0.457 | 1 | 3 | 52 |
| 531 | NUP214 | ENSG00000126883 | 0.457 | 1 | 2 | 54 |
| 532 | FGFR3 | ENSG00000068078 | 0.458 | 1 | 30 | 49 |
| 533 | RERE | ENSG00000142599 | 0.459 | 1 | 14 | 54 |
| 534 | SERPINI1 | ENSG00000163536 | 0.459 | 1 | 1 | 52 |
| 535 | EIF3F | ENSG00000175390 | 0.460 | 1 | 1 | 54 |
| 536 | NCKAP1L | ENSG00000123338 | 0.460 | 1 | 1 | 43 |
| 537 | PLCG2 | ENSG00000197943 | 0.461 | 1 | 1 | 45 |
| 538 | TGFB2 | ENSG00000092969 | 0.463 | 1 | 34 | 46 |
| 539 | ANKH | ENSG00000154122 | 0.463 | 1 | 3 | 54 |
| 540 | MAST3 | ENSG00000099308 | 0.466 | 1 | 2 | 54 |
| 541 | ADAMTS10 | ENSG00000142303 | 0.467 | 1 | 2 | 53 |
| 542 | OSBPL2 | ENSG00000130703 | 0.467 | 1 | 1 | 54 |
| 543 | DIS3L2 | ENSG00000144535 | 0.468 | 1 | 14 | 54 |
| 544 | TBX3 | ENSG00000135111 | 0.468 | 1 | 1 | 52 |
| 545 | MAFB | ENSG00000204103 | 0.469 | 1 | 7 | 53 |
| 546 | UBR1 | ENSG00000159459 | 0.470 | 1 | 7 | 54 |
| 547 | IQSEC1 | ENSG00000144711 | 0.471 | 1 | 1 | 54 |
| 548 | HGF | ENSG00000019991 | 0.471 | 1 | 1 | 36 |
| 549 | ACTA2 | ENSG00000107796 | 0.472 | 1 | 19 | 54 |
| 550 | IKBKB | ENSG00000104365 | 0.472 | 1 | 3 | 54 |
| 551 | UBR2 | ENSG00000024048 | 0.472 | 1 | 1 | 54 |
| 552 | PDGFB | ENSG00000100311 | 0.473 | 1 | 2 | 52 |
| 553 | TCF7L2 | ENSG00000148737 | 0.474 | 1 | 2 | 53 |
| 554 | RAB3A | ENSG00000105649 | 0.474 | 1 | 1 | 54 |
| 555 | STIM1 | ENSG00000167323 | 0.475 | 1 | 5 | 54 |
| 556 | WLS | ENSG00000116729 | 0.475 | 1 | 2 | 53 |
| 557 | EGFR | ENSG00000146648 | 0.475 | 1 | 1 | 52 |
| 558 | RARB | ENSG00000077092 | 0.477 | 1 | 7 | 50 |
| 559 | TGIF1 | ENSG00000177426 | 0.477 | 1 | 1 | 52 |
| 560 | XYLT1 | ENSG00000103489 | 0.478 | 1 | 4 | 50 |
| 561 | AP2S1 | ENSG00000042753 | 0.478 | 1 | 3 | 54 |
| 562 | LRP5 | ENSG00000162337 | 0.482 | 1 | 51 | 53 |
| 563 | MITF | ENSG00000187098 | 0.482 | 1 | 29 | 49 |
| 564 | MYO6 | ENSG00000196586 | 0.482 | 1 | 22 | 53 |
| 565 | ZMYM2 | ENSG00000121741 | 0.482 | 1 | 11 | 54 |
| 566 | AKT3 | ENSG00000117020 | 0.482 | 1 | 2 | 54 |
| 567 | MN1 | ENSG00000169184 | 0.482 | 1 | 2 | 48 |
| 568 | MAPKBP1 | ENSG00000137802 | 0.484 | 1 | 1 | 54 |
| 569 | NDST1 | ENSG00000070614 | 0.484 | 1 | 1 | 54 |
| 570 | CSF1R | ENSG00000182578 | 0.485 | 1 | 15 | 52 |
| 571 | H3-3B | ENSG00000132475 | 0.485 | 1 | 4 | 54 |
| 572 | ARF3 | ENSG00000134287 | 0.485 | 1 | 2 | 54 |
| 573 | FLCN | ENSG00000154803 | 0.487 | 1 | 194 | 54 |
| 574 | REEP2 | ENSG00000132563 | 0.488 | 1 | 1 | 51 |
| 575 | ACO2 | ENSG00000100412 | 0.489 | 1 | 13 | 54 |
| 576 | CTNNA1 | ENSG00000044115 | 0.491 | 1 | 79 | 54 |
| 577 | DCHS1 | ENSG00000166341 | 0.491 | 1 | 4 | 51 |
| 578 | EHMT2 | ENSG00000204371 | 0.492 | 1 | 1 | 54 |
| 579 | SPTLC2 | ENSG00000100596 | 0.494 | 1 | 3 | 54 |
| 580 | MAPRE2 | ENSG00000166974 | 0.494 | 1 | 1 | 54 |
| 581 | JAG2 | ENSG00000184916 | 0.495 | 1 | 1 | 51 |
| 582 | SPRED1 | ENSG00000166068 | 0.496 | 1 | 35 | 53 |
| 583 | IGF2 | ENSG00000167244 | 0.496 | 1 | 5 | 52 |
| 584 | FOXP4 | ENSG00000137166 | 0.496 | 1 | 1 | 54 |
| 585 | PTPN4 | ENSG00000088179 | 0.496 | 1 | 1 | 54 |
| 586 | TMEM240 | ENSG00000205090 | 0.496 | 1 | 1 | 50 |
| 587 | STAT2 | ENSG00000170581 | 0.497 | 1 | 4 | 54 |
| 588 | HECW2 | ENSG00000138411 | 0.498 | 1 | 10 | 39 |
| 589 | ARNT2 | ENSG00000172379 | 0.498 | 1 | 1 | 51 |
| 590 | LBR | ENSG00000143815 | 0.499 | 1 | 6 | 54 |
| 591 | NFKBIA | ENSG00000100906 | 0.500 | 1 | 2 | 54 |
| 592 | LRP4 | ENSG00000134569 | 0.501 | 1 | 4 | 50 |
| 593 | OPA1 | ENSG00000198836 | 0.503 | 1 | 101 | 54 |
| 594 | COL11A2 | ENSG00000204248 | 0.503 | 1 | 43 | 43 |
| 595 | LOX | ENSG00000113083 | 0.503 | 1 | 20 | 42 |
| 596 | FLNB | ENSG00000136068 | 0.504 | 1 | 21 | 54 |
| 597 | PDGFRB | ENSG00000113721 | 0.504 | 1 | 4 | 54 |
| 598 | PSD | ENSG00000059915 | 0.504 | 1 | 1 | 44 |
| 599 | LYST | ENSG00000143669 | 0.505 | 1 | 56 | 54 |
| 600 | TBK1 | ENSG00000183735 | 0.505 | 1 | 19 | 54 |
| 601 | ATAD1 | ENSG00000138138 | 0.505 | 1 | 2 | 54 |
| 602 | DNM1L | ENSG00000087470 | 0.506 | 1 | 15 | 54 |
| 603 | PSAP | ENSG00000197746 | 0.507 | 1 | 20 | 54 |
| 604 | BICD2 | ENSG00000185963 | 0.507 | 1 | 10 | 54 |
| 605 | TGFBR1 | ENSG00000106799 | 0.509 | 1 | 24 | 54 |
| 606 | PLCB4 | ENSG00000101333 | 0.509 | 1 | 3 | 48 |
| 607 | MTHFD1 | ENSG00000100714 | 0.510 | 1 | 6 | 54 |
| 608 | STT3A | ENSG00000134910 | 0.510 | 1 | 3 | 54 |
| 609 | TPM3 | ENSG00000143549 | 0.511 | 1 | 8 | 54 |
| 610 | ACTC1 | ENSG00000159251 | 0.511 | 1 | 7 | 42 |
| 611 | EXTL3 | ENSG00000012232 | 0.511 | 1 | 1 | 54 |
| 612 | SPEG | ENSG00000072195 | 0.514 | 1 | 7 | 49 |
| 613 | CXCR4 | ENSG00000121966 | 0.514 | 1 | 5 | 53 |
| 614 | LRBA | ENSG00000198589 | 0.515 | 1 | 34 | 54 |
| 615 | MAP2K1 | ENSG00000169032 | 0.515 | 1 | 16 | 54 |
| 616 | GALNT2 | ENSG00000143641 | 0.515 | 1 | 1 | 54 |
| 617 | NPR2 | ENSG00000159899 | 0.517 | 1 | 18 | 53 |
| 618 | TIA1 | ENSG00000116001 | 0.517 | 1 | 1 | 54 |
| 619 | BLNK | ENSG00000095585 | 0.517 | 1 | 1 | 47 |
| 620 | ATR | ENSG00000175054 | 0.518 | 1 | 18 | 54 |
| 621 | FGD4 | ENSG00000139132 | 0.518 | 1 | 10 | 53 |
| 622 | DHX37 | ENSG00000150990 | 0.518 | 1 | 2 | 54 |
| 623 | AP3B2 | ENSG00000103723 | 0.520 | 1 | 8 | 33 |
| 624 | ITGA6 | ENSG00000091409 | 0.520 | 1 | 5 | 53 |
| 625 | SNX27 | ENSG00000143376 | 0.520 | 1 | 1 | 54 |
| 626 | OXCT1 | ENSG00000083720 | 0.521 | 1 | 5 | 53 |
| 627 | GJA1 | ENSG00000152661 | 0.522 | 1 | 11 | 52 |
| 628 | EIF2B3 | ENSG00000070785 | 0.522 | 1 | 5 | 54 |
| 629 | RIT1 | ENSG00000143622 | 0.523 | 1 | 26 | 54 |
| 630 | CNOT9 | ENSG00000144580 | 0.523 | 1 | 2 | 54 |
| 631 | PIKFYVE | ENSG00000115020 | 0.523 | 1 | 2 | 54 |
| 632 | EXOC8 | ENSG00000116903 | 0.523 | 1 | 1 | 54 |
| 633 | SLC25A12 | ENSG00000115840 | 0.525 | 1 | 3 | 54 |
| 634 | INF2 | ENSG00000203485 | 0.526 | 1 | 19 | 54 |
| 635 | GATA6 | ENSG00000141448 | 0.527 | 1 | 4 | 34 |
| 636 | SPTBN4 | ENSG00000160460 | 0.529 | 1 | 1 | 39 |
| 637 | PROS1 | ENSG00000184500 | 0.530 | 1 | 27 | 52 |
| 638 | TGFBR2 | ENSG00000163513 | 0.531 | 1 | 35 | 54 |
| 639 | CCNF | ENSG00000162063 | 0.531 | 1 | 1 | 46 |
| 640 | SEC63 | ENSG00000025796 | 0.532 | 1 | 11 | 54 |
| 641 | SLC25A3 | ENSG00000075415 | 0.532 | 1 | 1 | 54 |
| 642 | CTBP1 | ENSG00000159692 | 0.533 | 1 | 1 | 54 |
| 643 | JUP | ENSG00000173801 | 0.533 | 1 | 1 | 54 |
| 644 | PRKG1 | ENSG00000185532 | 0.534 | 1 | 1 | 49 |
| 645 | RAB3GAP2 | ENSG00000118873 | 0.535 | 1 | 3 | 54 |
| 646 | DHDDS | ENSG00000117682 | 0.536 | 1 | 15 | 54 |
| 647 | DVL3 | ENSG00000161202 | 0.536 | 1 | 1 | 54 |
| 648 | ANLN | ENSG00000011426 | 0.537 | 1 | 1 | 33 |
| 649 | BAG3 | ENSG00000151929 | 0.538 | 1 | 49 | 54 |
| 650 | KRT5 | ENSG00000186081 | 0.538 | 1 | 30 | 46 |
| 651 | RASGRP1 | ENSG00000172575 | 0.538 | 1 | 4 | 27 |
| 652 | KLHL20 | ENSG00000076321 | 0.538 | 1 | 1 | 54 |
| 653 | PC | ENSG00000173599 | 0.539 | 1 | 30 | 53 |
| 654 | SYNJ1 | ENSG00000159082 | 0.539 | 1 | 12 | 54 |
| 655 | COL27A1 | ENSG00000196739 | 0.539 | 1 | 12 | 52 |
| 656 | NUP188 | ENSG00000095319 | 0.539 | 1 | 1 | 54 |
| 657 | RALA | ENSG00000006451 | 0.540 | 1 | 2 | 54 |
| 658 | C3 | ENSG00000125730 | 0.541 | 1 | 13 | 53 |
| 659 | RDX | ENSG00000137710 | 0.542 | 1 | 3 | 54 |
| 660 | ROR2 | ENSG00000169071 | 0.543 | 1 | 6 | 33 |
| 661 | TNFRSF11B | ENSG00000164761 | 0.543 | 1 | 1 | 29 |
| 662 | EIF2B5 | ENSG00000145191 | 0.544 | 1 | 43 | 54 |
| 663 | IRF2BPL | ENSG00000119669 | 0.544 | 1 | 22 | 54 |
| 664 | CAMSAP1 | ENSG00000130559 | 0.544 | 1 | 2 | 54 |
| 665 | SLC6A9 | ENSG00000196517 | 0.544 | 1 | 2 | 51 |
| 666 | RHBDF2 | ENSG00000129667 | 0.545 | 1 | 1 | 54 |
| 667 | HMBS | ENSG00000256269 | 0.546 | 1 | 36 | 54 |
| 668 | RBM8A | ENSG00000265241 | 0.546 | 1 | 3 | 54 |
| 669 | PI4KA | ENSG00000241973 | 0.548 | 1 | 7 | 54 |
| 670 | DYRK1B | ENSG00000105204 | 0.548 | 1 | 1 | 54 |
| 671 | EPRS1 | ENSG00000136628 | 0.549 | 1 | 2 | 54 |
| 672 | NRAS | ENSG00000213281 | 0.550 | 1 | 6 | 54 |
| 673 | LMNB1 | ENSG00000113368 | 0.550 | 1 | 3 | 51 |
| 674 | DOCK7 | ENSG00000116641 | 0.551 | 1 | 9 | 53 |
| 675 | INSR | ENSG00000171105 | 0.556 | 1 | 8 | 54 |
| 676 | CALM3 | ENSG00000160014 | 0.556 | 1 | 3 | 54 |
| 677 | FERMT3 | ENSG00000149781 | 0.556 | 1 | 1 | 54 |
| 678 | ABCC9 | ENSG00000069431 | 0.557 | 1 | 10 | 49 |
| 679 | MYO1E | ENSG00000157483 | 0.557 | 1 | 2 | 52 |
| 680 | ARL3 | ENSG00000138175 | 0.557 | 1 | 1 | 54 |
| 681 | PTHLH | ENSG00000087494 | 0.558 | 1 | 4 | 27 |
| 682 | EFEMP2 | ENSG00000172638 | 0.560 | 1 | 7 | 54 |
| 683 | P4HB | ENSG00000185624 | 0.561 | 1 | 1 | 54 |
| 684 | PACS2 | ENSG00000179364 | 0.561 | 1 | 1 | 54 |
| 685 | ABCA1 | ENSG00000165029 | 0.563 | 1 | 9 | 53 |
| 686 | DSE | ENSG00000111817 | 0.563 | 1 | 1 | 53 |
| 687 | SYNE1 | ENSG00000131018 | 0.564 | 1 | 86 | 54 |
| 688 | COL6A3 | ENSG00000163359 | 0.564 | 1 | 55 | 39 |
| 689 | PLCE1 | ENSG00000138193 | 0.564 | 1 | 14 | 51 |
| 690 | NPRL2 | ENSG00000114388 | 0.564 | 1 | 5 | 54 |
| 691 | CPE | ENSG00000109472 | 0.565 | 1 | 2 | 52 |
| 692 | ADAM17 | ENSG00000151694 | 0.565 | 1 | 1 | 54 |
| 693 | NOTCH3 | ENSG00000074181 | 0.566 | 1 | 125 | 52 |
| 694 | KCNH2 | ENSG00000055118 | 0.568 | 1 | 268 | 38 |
| 695 | ADAMTS2 | ENSG00000087116 | 0.568 | 1 | 26 | 47 |
| 696 | C19orf12 | ENSG00000131943 | 0.568 | 1 | 12 | 54 |
| 697 | LCAT | ENSG00000213398 | 0.568 | 1 | 12 | 54 |
| 698 | RBM20 | ENSG00000203867 | 0.568 | 1 | 12 | 33 |
| 699 | ACVRL1 | ENSG00000139567 | 0.570 | 1 | 211 | 52 |
| 700 | ATL1 | ENSG00000198513 | 0.570 | 1 | 21 | 49 |
| 701 | NUP133 | ENSG00000069248 | 0.570 | 1 | 1 | 54 |
| 702 | THAP1 | ENSG00000131931 | 0.571 | 1 | 13 | 54 |
| 703 | GNPAT | ENSG00000116906 | 0.572 | 1 | 12 | 54 |
| 704 | ST3GAL3 | ENSG00000126091 | 0.573 | 1 | 2 | 54 |
| 705 | DENND5A | ENSG00000184014 | 0.574 | 1 | 3 | 54 |
| 706 | IPO8 | ENSG00000133704 | 0.575 | 1 | 5 | 54 |
| 707 | ADAM9 | ENSG00000168615 | 0.575 | 1 | 3 | 54 |
| 708 | HNRNPDL | ENSG00000152795 | 0.575 | 1 | 2 | 54 |
| 709 | RBP4 | ENSG00000138207 | 0.575 | 1 | 1 | 47 |
| 710 | AP4E1 | ENSG00000081014 | 0.576 | 1 | 4 | 54 |
| 711 | ATP5F1A | ENSG00000152234 | 0.576 | 1 | 4 | 54 |
| 712 | PTPN23 | ENSG00000076201 | 0.577 | 1 | 4 | 54 |
| 713 | ALB | ENSG00000163631 | 0.577 | 1 | 2 | 42 |
| 714 | RAB11A | ENSG00000103769 | 0.577 | 1 | 1 | 54 |
| 715 | CNTNAP1 | ENSG00000108797 | 0.578 | 1 | 8 | 53 |
| 716 | MYSM1 | ENSG00000162601 | 0.578 | 1 | 2 | 54 |
| 717 | NALCN | ENSG00000102452 | 0.579 | 1 | 40 | 27 |
| 718 | BLTP1 | ENSG00000138688 | 0.579 | 1 | 1 | 54 |
| 719 | EPG5 | ENSG00000152223 | 0.580 | 1 | 46 | 54 |
| 720 | MALT1 | ENSG00000172175 | 0.580 | 1 | 4 | 54 |
| 721 | BCAS3 | ENSG00000141376 | 0.580 | 1 | 1 | 54 |
| 722 | FBLN5 | ENSG00000140092 | 0.580 | 1 | 1 | 52 |
| 723 | SDHD | ENSG00000204370 | 0.581 | 1 | 78 | 54 |
| 724 | PPOX | ENSG00000143224 | 0.581 | 1 | 18 | 54 |
| 725 | AMOTL1 | ENSG00000166025 | 0.581 | 1 | 3 | 52 |
| 726 | SZT2 | ENSG00000198198 | 0.582 | 1 | 35 | 54 |
| 727 | LPIN2 | ENSG00000101577 | 0.582 | 1 | 12 | 54 |
| 728 | PRKAG2 | ENSG00000106617 | 0.582 | 1 | 7 | 54 |
| 729 | TUBG1 | ENSG00000131462 | 0.583 | 1 | 2 | 54 |
| 730 | LEMD3 | ENSG00000174106 | 0.584 | 2 | 9 | 54 |
| 731 | KIAA1549 | ENSG00000122778 | 0.584 | 1 | 1 | 40 |
| 732 | EDNRB | ENSG00000136160 | 0.585 | 2 | 9 | 52 |
| 733 | VPS4A | ENSG00000132612 | 0.585 | 2 | 2 | 54 |
| 734 | PDHB | ENSG00000168291 | 0.586 | 2 | 8 | 54 |
| 735 | SLC30A9 | ENSG00000014824 | 0.587 | 2 | 2 | 54 |
| 736 | NNT | ENSG00000112992 | 0.588 | 2 | 2 | 54 |
| 737 | LEPR | ENSG00000116678 | 0.589 | 2 | 5 | 48 |
| 738 | GRN | ENSG00000030582 | 0.590 | 2 | 38 | 54 |
| 739 | AP3B1 | ENSG00000132842 | 0.590 | 2 | 9 | 54 |
| 740 | SBF2 | ENSG00000133812 | 0.591 | 2 | 8 | 54 |
| 741 | SLC16A1 | ENSG00000155380 | 0.591 | 2 | 3 | 54 |
| 742 | TFG | ENSG00000114354 | 0.591 | 2 | 3 | 54 |
| 743 | GOT2 | ENSG00000125166 | 0.591 | 2 | 2 | 54 |
| 744 | SLC39A14 | ENSG00000104635 | 0.592 | 2 | 1 | 53 |
| 745 | LTBP2 | ENSG00000119681 | 0.595 | 2 | 12 | 53 |
| 746 | FAR1 | ENSG00000197601 | 0.595 | 2 | 4 | 54 |
| 747 | FHOD3 | ENSG00000134775 | 0.595 | 2 | 2 | 50 |
| 748 | GATM | ENSG00000171766 | 0.597 | 2 | 13 | 52 |
| 749 | DHCR24 | ENSG00000116133 | 0.597 | 2 | 2 | 54 |
| 750 | OSBPL8 | ENSG00000091039 | 0.597 | 2 | 1 | 54 |
| 751 | PDXK | ENSG00000160209 | 0.597 | 2 | 1 | 54 |
| 752 | VPS11 | ENSG00000160695 | 0.597 | 2 | 1 | 54 |
| 753 | ATP1A2 | ENSG00000018625 | 0.598 | 2 | 39 | 49 |
| 754 | FH | ENSG00000091483 | 0.599 | 2 | 229 | 54 |
| 755 | FAM171A2 | ENSG00000161682 | 0.599 | 2 | 2 | 31 |

## Star-1 pathogenic calls only

The 202 genes in the split-out. Same LOEUF and breadth filters, no allele at stars >= 2, at least one germline pathogenic or likely pathogenic allele at star 1, and none at star 0 or star null. `plp_star_1` is that star-1 count. Star 1 mixes single-submitter and conflicting review status.

| rank | symbol | ensg | loeuf | loeuf_decile | plp_star_1 | n_tissues_tpm_ge_1 |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 1 | PSMD11 | ENSG00000108671 | 0.082 | 0 | 2 | 54 |
| 2 | XPO1 | ENSG00000082898 | 0.083 | 0 | 18 | 54 |
| 3 | WAPL | ENSG00000062650 | 0.083 | 0 | 1 | 54 |
| 4 | PRPF19 | ENSG00000110107 | 0.099 | 0 | 1 | 54 |
| 5 | SF3B1 | ENSG00000115524 | 0.100 | 0 | 2 | 54 |
| 6 | EIF3B | ENSG00000106263 | 0.107 | 0 | 2 | 54 |
| 7 | RBBP5 | ENSG00000117222 | 0.144 | 0 | 2 | 54 |
| 8 | GIGYF2 | ENSG00000204120 | 0.150 | 0 | 2 | 54 |
| 9 | MACROH2A1 | ENSG00000113648 | 0.154 | 0 | 1 | 54 |
| 10 | JMJD1C | ENSG00000171988 | 0.157 | 0 | 4 | 54 |
| 11 | MAT2A | ENSG00000168906 | 0.158 | 0 | 1 | 54 |
| 12 | NUP62 | ENSG00000213024 | 0.158 | 0 | 1 | 54 |
| 13 | UPF2 | ENSG00000151461 | 0.159 | 0 | 1 | 54 |
| 14 | HELZ | ENSG00000198265 | 0.160 | 0 | 1 | 54 |
| 15 | INO80 | ENSG00000128908 | 0.165 | 0 | 2 | 54 |
| 16 | JAK1 | ENSG00000162434 | 0.168 | 0 | 1 | 54 |
| 17 | MAU2 | ENSG00000129933 | 0.175 | 0 | 11 | 54 |
| 18 | HMGCS1 | ENSG00000112972 | 0.176 | 0 | 6 | 54 |
| 19 | RPL6 | ENSG00000089009 | 0.176 | 0 | 1 | 54 |
| 20 | EIF3A | ENSG00000107581 | 0.180 | 0 | 1 | 54 |
| 21 | RUNX1T1 | ENSG00000079102 | 0.186 | 0 | 1 | 49 |
| 22 | CMIP | ENSG00000153815 | 0.189 | 0 | 1 | 54 |
| 23 | UBN1 | ENSG00000118900 | 0.193 | 0 | 1 | 54 |
| 24 | GTF2I | ENSG00000263001 | 0.195 | 0 | 1 | 54 |
| 25 | RALGAPB | ENSG00000170471 | 0.202 | 0 | 2 | 54 |
| 26 | RBM5 | ENSG00000003756 | 0.202 | 0 | 1 | 54 |
| 27 | ARHGAP35 | ENSG00000160007 | 0.204 | 0 | 14 | 54 |
| 28 | PHF5A | ENSG00000100410 | 0.205 | 0 | 3 | 54 |
| 29 | EIF4G3 | ENSG00000075151 | 0.213 | 0 | 1 | 54 |
| 30 | SNW1 | ENSG00000100603 | 0.216 | 0 | 1 | 54 |
| 31 | NCAM1 | ENSG00000149294 | 0.217 | 0 | 2 | 48 |
| 32 | RPL31 | ENSG00000071082 | 0.217 | 0 | 1 | 54 |
| 33 | DDX39B | ENSG00000198563 | 0.219 | 0 | 3 | 54 |
| 34 | SF1 | ENSG00000168066 | 0.220 | 0 | 2 | 54 |
| 35 | PIP4K2B | ENSG00000276293 | 0.224 | 0 | 1 | 54 |
| 36 | ATXN7L3 | ENSG00000087152 | 0.227 | 0 | 1 | 54 |
| 37 | DIP2C | ENSG00000151240 | 0.235 | 0 | 23 | 52 |
| 38 | CCT5 | ENSG00000150753 | 0.235 | 0 | 1 | 54 |
| 39 | SUPT5H | ENSG00000196235 | 0.235 | 0 | 1 | 54 |
| 40 | TLN1 | ENSG00000137076 | 0.236 | 0 | 1 | 54 |
| 41 | LDB1 | ENSG00000198728 | 0.240 | 0 | 12 | 54 |
| 42 | RABGAP1 | ENSG00000011454 | 0.241 | 0 | 2 | 54 |
| 43 | MYCBP2 | ENSG00000005810 | 0.243 | 0 | 16 | 54 |
| 44 | SAP130 | ENSG00000136715 | 0.243 | 0 | 1 | 54 |
| 45 | ZFHX4 | ENSG00000091656 | 0.244 | 0 | 21 | 36 |
| 46 | UBE2H | ENSG00000186591 | 0.244 | 0 | 1 | 54 |
| 47 | SIN3B | ENSG00000127511 | 0.245 | 0 | 2 | 54 |
| 48 | ZNF536 | ENSG00000198597 | 0.245 | 0 | 2 | 27 |
| 49 | WDTC1 | ENSG00000142784 | 0.246 | 0 | 1 | 54 |
| 50 | DDX17 | ENSG00000100201 | 0.252 | 0 | 6 | 54 |
| 51 | RIF1 | ENSG00000080345 | 0.254 | 0 | 3 | 54 |
| 52 | EPHA4 | ENSG00000116106 | 0.256 | 0 | 2 | 50 |
| 53 | STX1A | ENSG00000106089 | 0.257 | 0 | 6 | 53 |
| 54 | GIT1 | ENSG00000108262 | 0.257 | 0 | 2 | 54 |
| 55 | PSMC1 | ENSG00000100764 | 0.258 | 0 | 1 | 54 |
| 56 | CECR2 | ENSG00000099954 | 0.259 | 0 | 3 | 37 |
| 57 | REV3L | ENSG00000009413 | 0.261 | 0 | 4 | 53 |
| 58 | ARHGAP21 | ENSG00000107863 | 0.261 | 0 | 1 | 54 |
| 59 | EFNB2 | ENSG00000125266 | 0.261 | 0 | 1 | 52 |
| 60 | TSHZ3 | ENSG00000121297 | 0.263 | 0 | 1 | 51 |
| 61 | HECTD1 | ENSG00000092148 | 0.265 | 0 | 8 | 54 |
| 62 | ELF2 | ENSG00000109381 | 0.271 | 0 | 1 | 54 |
| 63 | SEC24C | ENSG00000176986 | 0.275 | 0 | 1 | 54 |
| 64 | TEAD1 | ENSG00000187079 | 0.276 | 0 | 1 | 52 |
| 65 | MAP1A | ENSG00000166963 | 0.277 | 0 | 1 | 50 |
| 66 | PRKCE | ENSG00000171132 | 0.278 | 0 | 2 | 54 |
| 67 | FASN | ENSG00000169710 | 0.284 | 0 | 1 | 54 |
| 68 | CASZ1 | ENSG00000130940 | 0.285 | 0 | 16 | 40 |
| 69 | DPYSL2 | ENSG00000092964 | 0.290 | 0 | 1 | 54 |
| 70 | RB1CC1 | ENSG00000023287 | 0.293 | 0 | 5 | 54 |
| 71 | GATAD2A | ENSG00000167491 | 0.295 | 0 | 1 | 54 |
| 72 | TRIM24 | ENSG00000122779 | 0.297 | 0 | 1 | 54 |
| 73 | ANKS1B | ENSG00000185046 | 0.298 | 0 | 2 | 28 |
| 74 | ZNF496 | ENSG00000162714 | 0.300 | 0 | 2 | 54 |
| 75 | ADGRL1 | ENSG00000072071 | 0.302 | 0 | 31 | 51 |
| 76 | SYK | ENSG00000165025 | 0.302 | 0 | 3 | 46 |
| 77 | RICTOR | ENSG00000164327 | 0.305 | 0 | 2 | 54 |
| 78 | GTF2IRD1 | ENSG00000006704 | 0.309 | 0 | 1 | 53 |
| 79 | RAB1A | ENSG00000138069 | 0.311 | 0 | 1 | 54 |
| 80 | UBQLN4 | ENSG00000160803 | 0.314 | 0 | 1 | 54 |
| 81 | NRXN2 | ENSG00000110076 | 0.317 | 0 | 1 | 42 |
| 82 | PTPN1 | ENSG00000196396 | 0.321 | 0 | 2 | 54 |
| 83 | MDM2 | ENSG00000135679 | 0.325 | 0 | 1 | 54 |
| 84 | KDR | ENSG00000128052 | 0.328 | 0 | 1 | 52 |
| 85 | PHF2 | ENSG00000197724 | 0.329 | 0 | 1 | 54 |
| 86 | ADAMTS1 | ENSG00000154734 | 0.329 | 0 | 1 | 52 |
| 87 | PTPRM | ENSG00000173482 | 0.334 | 0 | 1 | 52 |
| 88 | ZMYND8 | ENSG00000101040 | 0.336 | 0 | 5 | 54 |
| 89 | RBFOX2 | ENSG00000100320 | 0.336 | 0 | 5 | 52 |
| 90 | SEMA6A | ENSG00000092421 | 0.337 | 0 | 2 | 52 |
| 91 | TNIK | ENSG00000154310 | 0.338 | 0 | 1 | 54 |
| 92 | NPAS2 | ENSG00000170485 | 0.338 | 0 | 1 | 52 |
| 93 | PRICKLE2 | ENSG00000163637 | 0.340 | 0 | 3 | 52 |
| 94 | CCNK | ENSG00000090061 | 0.340 | 0 | 1 | 54 |
| 95 | ADCY1 | ENSG00000164742 | 0.346 | 0 | 1 | 37 |
| 96 | ARHGAP32 | ENSG00000134909 | 0.349 | 0 | 1 | 52 |
| 97 | MAZ | ENSG00000103495 | 0.350 | 0 | 1 | 54 |
| 98 | HIF1A | ENSG00000100644 | 0.353 | 0 | 5 | 54 |
| 99 | PITPNM3 | ENSG00000091622 | 0.354 | 0 | 3 | 47 |
| 100 | UBQLN1 | ENSG00000135018 | 0.357 | 0 | 1 | 54 |
| 101 | ARPC4 | ENSG00000241553 | 0.362 | 0 | 1 | 54 |
| 102 | UBN2 | ENSG00000157741 | 0.362 | 0 | 1 | 54 |
| 103 | SYCP2 | ENSG00000196074 | 0.365 | 0 | 6 | 44 |
| 104 | JPH1 | ENSG00000104369 | 0.365 | 0 | 4 | 33 |
| 105 | RPL18 | ENSG00000063177 | 0.365 | 0 | 1 | 54 |
| 106 | ANKRD44 | ENSG00000065413 | 0.365 | 0 | 1 | 52 |
| 107 | CAMK1D | ENSG00000183049 | 0.365 | 0 | 1 | 51 |
| 108 | RREB1 | ENSG00000124782 | 0.366 | 0 | 8 | 53 |
| 109 | NUP205 | ENSG00000155561 | 0.366 | 0 | 2 | 54 |
| 110 | AGAP1 | ENSG00000157985 | 0.366 | 0 | 1 | 52 |
| 111 | DENND2B | ENSG00000166444 | 0.373 | 0 | 1 | 52 |
| 112 | SH2B1 | ENSG00000178188 | 0.374 | 0 | 1 | 54 |
| 113 | NDC1 | ENSG00000058804 | 0.375 | 0 | 1 | 54 |
| 114 | NCOA3 | ENSG00000124151 | 0.376 | 0 | 1 | 54 |
| 115 | EPB41L1 | ENSG00000088367 | 0.386 | 0 | 1 | 52 |
| 116 | IREB2 | ENSG00000136381 | 0.389 | 0 | 2 | 54 |
| 117 | LCP1 | ENSG00000136167 | 0.390 | 0 | 1 | 53 |
| 118 | HNRNPUL2 | ENSG00000214753 | 0.391 | 0 | 3 | 54 |
| 119 | KHDRBS1 | ENSG00000121774 | 0.394 | 0 | 1 | 54 |
| 120 | PREP | ENSG00000085377 | 0.395 | 0 | 2 | 54 |
| 121 | AKT2 | ENSG00000105221 | 0.398 | 0 | 2 | 54 |
| 122 | SLC6A6 | ENSG00000131389 | 0.399 | 0 | 3 | 54 |
| 123 | NLGN2 | ENSG00000169992 | 0.401 | 1 | 1 | 53 |
| 124 | DDX23 | ENSG00000174243 | 0.412 | 1 | 5 | 54 |
| 125 | C4B | ENSG00000224389 | 0.413 | 1 | 3 | 40 |
| 126 | MYO9A | ENSG00000066933 | 0.414 | 1 | 2 | 52 |
| 127 | PTBP2 | ENSG00000117569 | 0.420 | 1 | 2 | 54 |
| 128 | PPTC7 | ENSG00000196850 | 0.421 | 1 | 2 | 54 |
| 129 | IGF2R | ENSG00000197081 | 0.421 | 1 | 1 | 54 |
| 130 | CALR | ENSG00000179218 | 0.423 | 1 | 1 | 54 |
| 131 | EDNRA | ENSG00000151617 | 0.424 | 1 | 1 | 48 |
| 132 | CREBL2 | ENSG00000111269 | 0.427 | 1 | 1 | 54 |
| 133 | SV2A | ENSG00000159164 | 0.429 | 1 | 3 | 49 |
| 134 | P4HA1 | ENSG00000122884 | 0.429 | 1 | 1 | 54 |
| 135 | MYO9B | ENSG00000099331 | 0.430 | 1 | 1 | 54 |
| 136 | RYR3 | ENSG00000198838 | 0.430 | 1 | 1 | 34 |
| 137 | FRYL | ENSG00000075539 | 0.434 | 1 | 24 | 54 |
| 138 | SGMS1 | ENSG00000198964 | 0.439 | 1 | 2 | 54 |
| 139 | WDR47 | ENSG00000085433 | 0.440 | 1 | 2 | 54 |
| 140 | NCAPG2 | ENSG00000146918 | 0.441 | 1 | 1 | 51 |
| 141 | DENND4A | ENSG00000174485 | 0.443 | 1 | 1 | 54 |
| 142 | DLG2 | ENSG00000150672 | 0.444 | 1 | 4 | 37 |
| 143 | CAPZA2 | ENSG00000198898 | 0.444 | 1 | 2 | 54 |
| 144 | DOCK4 | ENSG00000128512 | 0.445 | 1 | 4 | 53 |
| 145 | PRPF6 | ENSG00000101161 | 0.445 | 1 | 2 | 54 |
| 146 | THBS2 | ENSG00000186340 | 0.446 | 1 | 1 | 52 |
| 147 | BRWD1 | ENSG00000185658 | 0.447 | 1 | 2 | 54 |
| 148 | ARHGAP29 | ENSG00000137962 | 0.448 | 1 | 8 | 52 |
| 149 | DGKZ | ENSG00000149091 | 0.451 | 1 | 1 | 54 |
| 150 | DDX39A | ENSG00000123136 | 0.455 | 1 | 1 | 54 |
| 151 | VCL | ENSG00000035403 | 0.461 | 1 | 4 | 54 |
| 152 | HDAC3 | ENSG00000171720 | 0.463 | 1 | 1 | 54 |
| 153 | DCAF8 | ENSG00000132716 | 0.466 | 1 | 1 | 54 |
| 154 | PLCB3 | ENSG00000149782 | 0.468 | 1 | 1 | 54 |
| 155 | PLK2 | ENSG00000145632 | 0.471 | 1 | 1 | 53 |
| 156 | RPS20 | ENSG00000008988 | 0.481 | 1 | 9 | 54 |
| 157 | ECE1 | ENSG00000117298 | 0.482 | 1 | 1 | 54 |
| 158 | LEO1 | ENSG00000166477 | 0.488 | 1 | 3 | 54 |
| 159 | NISCH | ENSG00000010322 | 0.488 | 1 | 1 | 54 |
| 160 | TFRC | ENSG00000072274 | 0.489 | 1 | 1 | 54 |
| 161 | SMAP1 | ENSG00000112305 | 0.491 | 1 | 1 | 54 |
| 162 | NUDCD3 | ENSG00000015676 | 0.492 | 1 | 1 | 54 |
| 163 | NMNAT2 | ENSG00000157064 | 0.494 | 1 | 4 | 35 |
| 164 | B2M | ENSG00000166710 | 0.501 | 1 | 5 | 54 |
| 165 | ZNRF3 | ENSG00000183579 | 0.503 | 1 | 7 | 52 |
| 166 | SENP7 | ENSG00000138468 | 0.503 | 1 | 3 | 54 |
| 167 | ERC1 | ENSG00000082805 | 0.503 | 1 | 1 | 54 |
| 168 | GAB1 | ENSG00000109458 | 0.503 | 1 | 1 | 54 |
| 169 | WNK2 | ENSG00000165238 | 0.508 | 1 | 1 | 44 |
| 170 | GSK3B | ENSG00000082701 | 0.512 | 1 | 6 | 54 |
| 171 | TASP1 | ENSG00000089123 | 0.518 | 1 | 5 | 53 |
| 172 | NKTR | ENSG00000114857 | 0.520 | 1 | 1 | 54 |
| 173 | BMP6 | ENSG00000153162 | 0.523 | 1 | 1 | 51 |
| 174 | PLXNA2 | ENSG00000076356 | 0.525 | 1 | 1 | 52 |
| 175 | ASH2L | ENSG00000129691 | 0.529 | 1 | 1 | 54 |
| 176 | HDAC1 | ENSG00000116478 | 0.529 | 1 | 1 | 54 |
| 177 | ATP5F1B | ENSG00000110955 | 0.531 | 1 | 1 | 54 |
| 178 | LRRC37A2 | ENSG00000238083 | 0.534 | 1 | 1 | 53 |
| 179 | IRF5 | ENSG00000128604 | 0.536 | 1 | 1 | 51 |
| 180 | BICC1 | ENSG00000122870 | 0.541 | 1 | 4 | 41 |
| 181 | BIRC3 | ENSG00000023445 | 0.543 | 1 | 7 | 43 |
| 182 | HOXA3 | ENSG00000105997 | 0.544 | 1 | 1 | 33 |
| 183 | CALU | ENSG00000128595 | 0.556 | 1 | 1 | 54 |
| 184 | SLC4A3 | ENSG00000114923 | 0.558 | 1 | 2 | 49 |
| 185 | CCDC186 | ENSG00000165813 | 0.560 | 1 | 2 | 54 |
| 186 | DNER | ENSG00000187957 | 0.561 | 1 | 1 | 37 |
| 187 | NRP2 | ENSG00000118257 | 0.562 | 1 | 1 | 47 |
| 188 | ODF2 | ENSG00000136811 | 0.564 | 1 | 1 | 54 |
| 189 | YPEL3 | ENSG00000090238 | 0.565 | 1 | 1 | 54 |
| 190 | UIMC1 | ENSG00000087206 | 0.566 | 1 | 1 | 54 |
| 191 | SARS1 | ENSG00000031698 | 0.570 | 1 | 2 | 54 |
| 192 | ADIPOR1 | ENSG00000159346 | 0.571 | 1 | 1 | 54 |
| 193 | RIC1 | ENSG00000107036 | 0.577 | 1 | 1 | 54 |
| 194 | FOS | ENSG00000170345 | 0.581 | 1 | 1 | 54 |
| 195 | VPS33A | ENSG00000139719 | 0.583 | 1 | 1 | 54 |
| 196 | FEM1C | ENSG00000145780 | 0.585 | 2 | 1 | 54 |
| 197 | CLPX | ENSG00000166855 | 0.586 | 2 | 1 | 54 |
| 198 | KAT2A | ENSG00000108773 | 0.587 | 2 | 1 | 54 |
| 199 | CEMIP2 | ENSG00000135048 | 0.588 | 2 | 2 | 54 |
| 200 | CITED2 | ENSG00000164442 | 0.588 | 2 | 1 | 54 |
| 201 | PDZRN3 | ENSG00000121440 | 0.591 | 2 | 1 | 51 |
| 202 | MEST | ENSG00000106484 | 0.592 | 2 | 1 | 52 |
