# Somatic versus germline

## Question

What are the counts and top genes for somatic origin and for oncogenicity, held apart from the germline rare-disease ranking?

## Filters

Rows are deduplicated GRCh38 ClinVar alleles in `warehouse/clinvar/allele_dedup.parquet`. One row is one `variant_key`, or one AlleleID when the record has no VCF allele. The query is `warehouse/reports/sql/somatic_versus_germline.sql`. It does not join the gene mart and does not rewrite a mart.

Somatic origin is `germline_somatic = 'somatic'`. That class is `OriginSimple = 'somatic'`. Mixed `germline/somatic` is not called somatic, and it is not called germline. `other` is `not applicable`, `unknown`, or `not provided`. On this file the class column matches `OriginSimple` with no mismatches.

An oncogenicity call is an `Oncogenicity` value other than blank, `-`, `not provided`, or a phrase that starts with `no classification`, `no assertion`, or `no interpretation`. That is the same rule as `onco_has_call` in `scripts/load_clinvar.py`. `Uncertain significance`, benign, likely benign, and conflicting oncogenicity are calls. They are not the gene ranking.

The oncogenicity gene table is sorted by oncogenic direction only: `Oncogenic`, `Likely oncogenic`, or `Oncogenic/Likely oncogenic`.

`germline_pathogenic` stays the existing flag: germline origin only, no oncogenicity call, and a pathogenic or likely pathogenic clinical-significance piece with no conflicting piece. This report does not compute the rare-disease ranking. That ranking also uses two or more review stars, LOEUF, and expression breadth. Oncogenicity is not added to the flag.

A gene is a ClinVar `GeneSymbol` token split on `;`. Empty tokens and `-`, `.`, and `na` are dropped. No symbol in this file contains `|`. An allele with two symbols is counted once for each symbol. The columns inside a gene row overlap. They are not a partition of that gene's alleles.

`review_stars` is the germline review status stored on the dedup row. `ReviewStatusOncogenicity` and `ReviewStatusClinicalImpact` are not columns on `allele_dedup`, so this report does not star-filter oncogenicity or somatic clinical impact.

## Denominator

Both questions use every dedup allele.

| question | numerator | denominator |
| --- | ---: | ---: |
| somatic origin (`germline_somatic = 'somatic'`) | 7001 | 4639253 |
| oncogenicity call | 2405 | 4639253 |

The denominators are the same rows because both fields are present on every allele. The numerators are different sets. Gene shares for somatic origin are out of the 7001. Gene shares for oncogenic direction are out of the 1773 oncogenic-direction alleles, which sit inside the 2405 calls and inside the same 4639253.

Mixed origin is 3560 of the same 4639253. It is not added to the 7001.

Fractions below are rounded to six decimals. The four origin-class fractions sum to 0.999999.

## Result

Origin class. `germline_pathogenic_n` is the flag inside that class. It is zero outside germline.

| germline_somatic | alleles | fraction of 4639253 | germline_pathogenic_n |
| --- | ---: | ---: | ---: |
| germline | 4316349 | 0.930397 | 338587 |
| other | 312343 | 0.067326 | 0 |
| somatic | 7001 | 0.001509 | 0 |
| germline/somatic | 3560 | 0.000767 | 0 |

`other` is `not applicable` 243516, `unknown` 65105, and `not provided` 3722.

Somatic-origin clinical significance, out of 7001. All of these rows have `germline_pathogenic` false, including the rows whose text says pathogenic.

| ClinicalSignificance | alleles |
| --- | ---: |
| Uncertain significance | 2031 |
| `-` | 1838 |
| drug response | 1693 |
| Pathogenic | 731 |
| Likely pathogenic | 357 |
| not provided | 126 |
| Likely benign | 43 |
| other | 40 |
| Benign | 38 |
| association | 34 |
| Pathogenic/Likely pathogenic | 20 |
| Conflicting classifications of pathogenicity | 19 |
| risk factor | 10 |
| Likely pathogenic; drug response | 6 |
| Pathogenic; other | 5 |
| Likely risk allele | 3 |
| protective | 3 |
| Affects | 2 |
| Pathogenic; drug response | 1 |
| association not found | 1 |

Text that contains `pathogenic` and does not contain `conflict` covers 1120 of the 7001. Text that mentions drug response covers 1700, of which 1693 are exactly `drug response`.

Somatic-origin germline review stars, out of 7001. Null stars are `ReviewStatus` `-` (1838) or `no classification provided` (126). No somatic-origin allele has 3 or 4 stars.

| review_stars | alleles |
| --- | ---: |
| 0 | 4261 |
| 1 | 739 |
| 2 | 37 |
| null | 1964 |

Top 15 genes by somatic-origin allele count. Sorted by `somatic_n` only. The fifteen counts sum to 3098 of 7001 alleles. `drug_response_n`, `vus_n`, `review_stars_0_n`, and `oncogenic_direction_n` overlap. `germline_pathogenic_n` counts other alleles of the same symbol. It is not added to `somatic_n`.

| symbol | somatic_n | drug_response_n | vus_n | review_stars_0_n | oncogenic_direction_n | germline_pathogenic_n |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HOXB13 | 620 | 0 | 617 | 617 | 0 | 2 |
| ABCB1 | 595 | 595 | 0 | 595 | 0 | 4 |
| OPRM1 | 460 | 459 | 0 | 459 | 0 | 0 |
| CYP2D6 | 225 | 225 | 0 | 225 | 0 | 0 |
| TP53 | 173 | 2 | 13 | 74 | 56 | 712 |
| SPOP | 166 | 0 | 158 | 159 | 7 | 8 |
| COMT | 161 | 161 | 0 | 161 | 0 | 1 |
| EGFR | 158 | 48 | 41 | 83 | 14 | 68 |
| UGT2B7 | 127 | 127 | 0 | 127 | 0 | 0 |
| ARID1A | 98 | 1 | 3 | 28 | 51 | 126 |
| ATRX | 81 | 0 | 1 | 3 | 51 | 139 |
| PIK3CA | 71 | 1 | 4 | 14 | 9 | 31 |
| RB1 | 59 | 0 | 3 | 25 | 19 | 900 |
| DICER1 | 52 | 0 | 5 | 0 | 2 | 777 |
| PTEN | 52 | 1 | 3 | 9 | 21 | 1046 |

Symbols with at least one somatic-origin allele: 1381. Alleles with a usable symbol: 6996. Alleles with no usable symbol: 5. Alleles with more than one symbol: 20, which is why symbol rows are 7047 rather than 6996.

Mixed origin, against the same 4639253, and not inside the 7001. Of 3560 mixed alleles, 1748 have pathogenic text and no conflict word, and all 3560 have `germline_pathogenic` false. Top 15 by mixed-origin count:

| symbol | mixed_n | pathogenic_text_n | oncogenic_direction_n | germline_pathogenic_n |
| --- | ---: | ---: | ---: | ---: |
| TP53 | 444 | 299 | 247 | 712 |
| ABCC8 | 168 | 2 | 0 | 402 |
| KCNJ11 | 161 | 11 | 0 | 59 |
| HOXB13 | 142 | 0 | 0 | 2 |
| PTEN | 133 | 106 | 86 | 1046 |
| NF1 | 125 | 105 | 46 | 5480 |
| BRCA1 | 100 | 60 | 5 | 4077 |
| BRCA2 | 99 | 62 | 11 | 5406 |
| MLH1 | 78 | 58 | 6 | 1719 |
| MSH6 | 76 | 25 | 1 | 2020 |
| APC | 72 | 66 | 59 | 1998 |
| ATM | 72 | 53 | 22 | 3408 |
| MSH2 | 72 | 50 | 3 | 1918 |
| DICER1 | 61 | 32 | 2 | 777 |
| PIK3CA | 61 | 47 | 26 | 31 |

Mixed-origin symbols: 520. One mixed allele has no usable symbol. Ten have more than one symbol (3570 symbol rows).

Oncogenicity, out of 4639253. `germline_pathogenic_n` on a called label is zero.

| Oncogenicity | alleles | somatic | mixed | germline | other | germline_pathogenic_n |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `-` | 4636174 | 5672 | 2485 | 4316349 | 311668 | 338587 |
| Likely oncogenic | 1299 | 792 | 506 | 0 | 1 | 0 |
| no classification for the single variant | 674 | 0 | 0 | 0 | 674 | 0 |
| Uncertain significance | 617 | 423 | 194 | 0 | 0 | 0 |
| Oncogenic | 472 | 111 | 361 | 0 | 0 | 0 |
| Benign | 6 | 1 | 5 | 0 | 0 | 0 |
| Likely benign | 6 | 1 | 5 | 0 | 0 | 0 |
| Conflicting classifications of oncogenicity | 3 | 1 | 2 | 0 | 0 | 0 |
| Oncogenic/Likely oncogenic | 2 | 0 | 2 | 0 | 0 | 0 |

Calls are 2405 of 4639253 (0.000518). Directional labels are 1773 (0.000382): 1299 likely oncogenic, 472 oncogenic, and 2 `Oncogenic/Likely oncogenic`. The other calls are 617 uncertain, 12 benign or likely benign, and 3 conflicting. `no classification for the single variant` is not a call.

Of the 1773, 903 are somatic origin, 869 are mixed origin, 1 is other origin, and 0 are germline origin. Of the 2405 calls, 1329 are somatic origin, 1075 are mixed, 1 is other, and 0 are germline. The one other-origin call is `chr4:105243651:A:AT`, TET2, `OriginSimple` `unknown`, `Likely oncogenic`, clinical significance `-`, `germline_pathogenic` false.

Top 15 genes by oncogenic-direction allele count, out of 1773 alleles and 236 symbols. Sorted by `oncogenic_direction_n` only. `somatic_origin_n` and `mixed_origin_n` count directional alleles, not every allele of that origin. `germline_pathogenic_n` is not added to the directional count. NF2 and NRAS each include one `Oncogenic/Likely oncogenic` allele, so `oncogenic_n` plus `likely_oncogenic_n` is one short of `oncogenic_direction_n` on those two rows. The other rows add up.

| symbol | oncogenic_direction_n | oncogenic_n | likely_oncogenic_n | oncogenicity_call_n | somatic_origin_n | mixed_origin_n | germline_pathogenic_n |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TP53 | 303 | 147 | 156 | 343 | 56 | 247 | 712 |
| PTEN | 107 | 28 | 79 | 124 | 21 | 86 | 1046 |
| APC | 95 | 19 | 76 | 98 | 36 | 59 | 1998 |
| ARID1A | 56 | 3 | 53 | 65 | 51 | 5 | 126 |
| ATRX | 55 | 0 | 55 | 66 | 51 | 4 | 139 |
| NF1 | 54 | 9 | 45 | 64 | 8 | 46 | 5480 |
| RB1 | 44 | 11 | 33 | 54 | 19 | 25 | 900 |
| KMT2D | 40 | 0 | 40 | 42 | 30 | 10 | 1013 |
| PIK3CA | 35 | 17 | 18 | 48 | 9 | 26 | 31 |
| CDKN2A | 30 | 9 | 21 | 38 | 14 | 16 | 143 |
| SMAD4 | 29 | 3 | 26 | 40 | 14 | 15 | 272 |
| ATM | 28 | 6 | 22 | 48 | 6 | 22 | 3408 |
| NF2 | 27 | 17 | 9 | 30 | 18 | 9 | 277 |
| EGFR | 26 | 14 | 12 | 49 | 14 | 12 | 68 |
| NRAS | 24 | 23 | 0 | 24 | 13 | 11 | 13 |

Symbols with any oncogenicity call: 335. Call symbol-rows: 2409, from 4 multi-symbol alleles. Directional symbol-rows: 1774, from 1773 alleles. The extra directional row is the one multi-symbol allele `MFSD11;SRSF2`. Neither symbol is in the top 15.

The sets overlap, and none of the overlap is germline pathogenic:

| set | alleles | germline_pathogenic |
| --- | ---: | ---: |
| somatic origin | 7001 | 0 |
| mixed origin | 3560 | 0 |
| oncogenicity call | 2405 | 0 |
| oncogenic direction | 1773 | 0 |
| somatic origin and an oncogenicity call | 1329 | 0 |
| somatic origin and oncogenic direction | 903 | 0 |
| mixed origin and an oncogenicity call | 1075 | 0 |
| mixed origin and oncogenic direction | 869 | 0 |
| germline origin and an oncogenicity call | 0 | 0 |
| germline_pathogenic flag | 338587 | 338587 |

Truth-set allele `chr17:31334927:C:T`: GeneSymbol NF1, ClinicalSignificance Pathogenic, OriginSimple `germline/somatic`, Oncogenicity `Likely oncogenic`, SomaticClinicalImpact `-`, review_stars 2, `germline_pathogenic` false. It is one of the 3560 mixed-origin alleles and one of the 1773 directional alleles. It is not one of the 7001 somatic-origin alleles and not one of the 338587 germline-pathogenic alleles. NF1's gene row is 54 directional alleles against 5480 germline-pathogenic alleles. The 54 are not part of the 5480.

Leakage of a non-germline origin, or of an oncogenicity call, into `germline_pathogenic`: 0.

## Caveat

Somatic origin is not oncogenicity, and neither count is the germline rare-disease ranking. Of the 7001 somatic-origin alleles, 2031 are `Uncertain significance`, 1838 have clinical significance `-`, and 1693 are exactly `drug response`. The head of the somatic gene list is HOXB13 (620 alleles, 617 of them zero-star VUS, none oncogenic-direction) and pharmacogenes. ABCB1, OPRM1, CYP2D6, COMT, and UGT2B7 account for 1567 of the 1700 somatic-origin alleles whose clinical significance mentions drug response. Only 903 of the 7001 are oncogenic or likely oncogenic. A Pathogenic label does not move a somatic-origin allele into `germline_pathogenic`: 1120 somatic-origin alleles have that text and still have the flag false.

Mixed origin is the class most easily mistaken for the germline ranking. The NF1 truth allele is Pathogenic at two germline review stars and is still excluded, because its origin is mixed and it has an oncogenicity call. Across mixed origin, 1748 alleles have pathogenic text and no conflict word, and none are `germline_pathogenic`. Mixed top genes are not a cancer ranking either. ABCC8 has 168 mixed alleles, 2 with pathogenic text, and no oncogenic-direction allele. KCNJ11 and HOXB13 are mostly not oncogenic-direction.

Oncogenicity is sparse, and the directional gene table is not a second pathogenic ranking. TP53 has 303 directional alleles and 712 germline-pathogenic alleles. NF1 has 54 and 5480. Adding the directional count to the germline-pathogenic count would put oncogenicity into the germline ranking. This report does not do that, and the flag does not do that. No germline-origin allele has an oncogenicity call. TET2 `chr4:105243651:A:AT` is the one directional allele whose origin is `unknown`, so the oncogenicity numerator is not nested inside somatic or mixed origin.

`review_stars` does not rescue either set. It is the germline review status. Oncogenicity review status is not on this table. Somatic-origin alleles are mostly zero stars (4261) or null (1964). Thirty-seven have two stars, and none have more.

Somatic clinical impact is a third field and is not ranked here. A tier is present on 1316 of 4639253 alleles: Tier I - Strong 753, Tier II - Potential 427, Tier III - Unknown 95, Tier IV - Benign/Likely benign 41. `no classification for the single variant` (674) is not counted as a tier. None of the 1316 are `germline_pathogenic`. The overlap with an oncogenicity call is 229, not 1316. The NF1 truth allele has impact `-`.

These are ClinVar labels on deduplicated alleles. They are not a census of mutations in tumors. Symbol totals are not an allele partition, and a ClinVar symbol was not joined to GENCODE, so this report does not split a symbol that more than one GENCODE gene uses.
