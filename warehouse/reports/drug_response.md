# Drug response

## Question

How small is the ClinVar drug-response set, and what constraint and expression do those genes have?

This is ClinVar's drug-response clinical-significance piece only. It is not a pharmacogenetics catalog. CPIC and PharmGKB are not on disk.

Numbers are from `warehouse/reports/sql/drug_response.sql`. Fractions are DuckDB `round(x, 6)`.

## Filters

- Alleles are rows in `warehouse/clinvar/allele_dedup.parquet`. One row per `variant_key`, or per `AlleleID` when the record has no VCF allele.
- A drug-response allele has a clinical-significance piece equal to `drug response`. Pieces are the lowercased `ClinicalSignificance` string split on `|`, `;`, and `/`, trimmed, with empty pieces and `.`, `-`, `na`, `n/a`, `none`, and `not provided` removed. Every matching piece in this file is that exact string (0 alleles matched only by a longer substring).
- The class is not exclusive. An allele with another piece stays in this set. The gene mart's `clinvar_drug_response_n` uses the same piece rule.
- Origin is `germline_somatic` on the dedup row (`germline`, `somatic`, `germline/somatic`, or `other`), from `OriginSimple`. It is not the oncogenicity column.
- Review stars are the dedup `review_stars` value. The raw `ReviewStatus` strings are reported beside them.
- Genes: `GeneSymbol` split on `;` only. Tokens `.` and `-` are dropped. A token must equal `gene_mart.symbol`. It is not a substring. No drug-response `GeneSymbol` contains `|`.
- Gene summaries keep symbols that match exactly one GENCODE row. This run drops 0 symbols for collision.
- Constraint and expression are the gene-mart columns (`loeuf`, `loeuf_decile`, `lof.pLI`, `mis.z_score`, `max_tissue`, `max_median_tpm`, `n_tissues_tpm_ge_1`, `tau`). Those columns were already left-joined on bare ENSG. This report does not fill a miss by symbol.
- LOEUF &lt; 0.35, LOEUF &lt; 0.6, decile 0, deciles 0–2, pLI ≥ 0.9, missense z ≥ 3.09, and tau ≥ 0.8 are comparison cuts. They do not define the set. Rates use genes with a non-null value of that metric.

## Denominator

- Alleles: 4,639,253 dedup rows. Those rows cover 4,639,270 source rows.
- Drug-response claims are fractions of 4,639,253 unless a smaller filter is named.
- The 1,979 drug-response rows were not collapsed: each has `n_rows` 1, and there are 1,979 distinct `VariationID` values and 1,979 distinct `AlleleID` values. Nine have no `variant_key`.
- Genes: 78,691 rows on `gene_mart.parquet`, of which 20,097 are `protein_coding`.
- The mart column `clinvar_drug_response_n` exists. Its sum across genes is 1,978. That sum is gene links, not alleles. Six alleles match two symbols, so they are counted twice in the column. The allele denominator stays 1,979.
- Recomputed per-symbol allele counts equal `clinvar_drug_response_n` for all 110 matched symbols (0 disagreements). Symbols that match more than one GENCODE gene: 0.

## Result

The set is small: 1,979 / 4,639,253 alleles (0.000427). It touches 110 / 78,691 genes (0.001398), and 107 / 20,097 protein-coding genes (0.005324).

Most of the allele count is one unreviewed phenotype. `PhenotypeList` equals `Tramadol response` for 1,567 / 1,979 alleles (0.791814). All 1,567 are somatic origin and 0 stars (`no assertion criteria provided`). They sit on five genes, and no allele carries two of those symbols: ABCB1 595, OPRM1 459, CYP2D6 225, COMT 161, UGT2B7 127. Those five genes account for 1,583 / 1,979 alleles (0.799899). Dropping the exact `Tramadol response` string leaves 412 alleles. That remainder is still only ClinVar's label.

Review is thin. 1,775 / 1,979 alleles (0.896918) are 0 stars. Alleles at 3 or 4 stars: 161 / 1,979 (0.081354), of which 132 are `reviewed by expert panel` and 29 are `practice guideline`. All 29 practice-guideline alleles are germline CYP2C19. No somatic-origin allele has more than one star. Of the 161 alleles at 3 or 4 stars, 157 are germline-only and 4 are germline/somatic (two EGFR, one `MT-ND1;MT-RNR1;RNR1`, one TNF).

| stars | ReviewStatus | alleles |
| --- | --- | ---: |
| 0 | no assertion criteria provided | 1,775 |
| 1 | criteria provided, single submitter | 25 |
| 1 | criteria provided, conflicting classifications | 3 |
| 2 | criteria provided, multiple submitters, no conflicts | 15 |
| 3 | reviewed by expert panel | 132 |
| 4 | practice guideline | 29 |

| origin | alleles | star 0 | star 1 | star 2 | star 3 | star 4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| somatic | 1,700 | 1,683 | 17 | 0 | 0 | 0 |
| germline | 244 | 68 | 9 | 10 | 128 | 29 |
| germline/somatic | 18 | 7 | 2 | 5 | 4 | 0 |
| other | 17 | 17 | 0 | 0 | 0 | 0 |

Somatic origin is 1,700 / 1,979 (0.859020). Germline origin is 244 / 1,979 (0.123295). Among germline alleles, 157 / 244 (0.643443) are already at 3 or 4 stars. Oncogenicity is empty on 1,970 / 1,979 alleles (`Oncogenic` 5, `Likely oncogenic` 4). Somatic clinical impact is empty on 1,976 / 1,979. The origin label is not an oncogenicity call.

1,889 / 1,979 alleles (0.954522) have no other clinical-significance piece. 90 have at least one other piece. The warehouse `germline_pathogenic` flag is true for 51 of the 1,979; those 51 stay in this set, and the mart counts them in the pathogenic column as well. 1,769 / 1,979 are `single nucleotide variant`.

| ClinicalSignificance | alleles |
| --- | ---: |
| drug response | 1,889 |
| Pathogenic; drug response | 40 |
| Likely pathogenic; drug response | 17 |
| Uncertain significance; drug response | 7 |
| Likely benign; drug response; other | 4 |
| Benign/Likely benign; drug response | 3 |
| Benign/Likely benign; drug response; other | 3 |
| Benign; drug response | 3 |
| Conflicting classifications of pathogenicity; drug response | 3 |
| Likely benign; drug response | 3 |
| Pathogenic/Likely pathogenic; drug response | 2 |
| drug response; risk factor | 2 |
| Conflicting classifications of pathogenicity; drug response; other | 1 |
| association; drug response; risk factor | 1 |
| drug response; other | 1 |

Phenotype strings with at least 10 alleles. Four somatic strings, Tramadol, prednisolone, tyrosine kinase inhibitor, and corticosteroids, account for 1,671 of the 1,700 somatic-origin alleles.

| PhenotypeList | alleles | star 0 | somatic | germline |
| --- | ---: | ---: | ---: | ---: |
| Tramadol response | 1,567 | 1,567 | 1,567 | 0 |
| Prednisolone response | 38 | 38 | 38 | 0 |
| Tyrosine kinase inhibitor response | 36 | 23 | 36 | 0 |
| Corticosteroids response | 30 | 30 | 30 | 0 |
| Thalidomide response | 15 | 15 | 0 | 15 |
| Warfarin response | 11 | 8 | 0 | 7 |

Gene link: 1,972 / 1,979 alleles match at least one GENCODE symbol. The other 7 do not. Three have `GeneSymbol` `-` and `GeneID` `-1` (germline, 0 stars, significance `drug response`). Four have only the token `covers 10 genes, none of which curated to show dosage sensitivity`, which is not a symbol and was not split into ten genes. Two further tokens do not match, `RNR1` (2 alleles) and `LOC129996006` (1 allele), but those alleles also carry a symbol that does match. Of the 110 matched genes, 71 have exactly one drug-response allele. Biotypes: protein-coding 107 genes (1,972 links), Mt_rRNA 1 gene (4 links, MT-RNR1), Mt_tRNA 1 gene (1 link, MT-TS1), lncRNA 1 gene (1 link, CASC6).

Constraint and expression, protein-coding genes. The drug-response genes are not a constrained class. Median LOEUF is 0.913 against 0.968 for all protein-coding genes with a LOEUF. The share below 0.6 is 23 / 104 (0.221154) against 3,895 / 18,098 (0.215217). Decile 0 is 12 / 104 (0.115385) against 1,851 / 17,840 (0.103756). Median pLI is 0 in both. pLI ≥ 0.9 is 16 / 104 (0.153846) against 3,158 / 18,098 (0.174494).

Expression is only mildly shifted. Median tau is 0.881 against 0.810. Tau ≥ 0.8 is 66 / 106 (0.622642) against 9,724 / 18,867 (0.515397). Median tissues at or above 1 TPM is 41.5 against 52.0. Median peak median-TPM is 59.31 against 43.21. That is a somewhat narrower profile, not one tissue program. Liver is the max tissue for several of the named genes below (CYP2D6, UGT2B7, CYP2C19, CYP2C9). It is not the max tissue for the set. RYR1 peaks in skeletal muscle, CFTR in pancreas, and OPRM1's maximum is 2.53 TPM in testis with only 2 tissues at or above 1 TPM.

| metric | drug-response protein-coding | all protein-coding |
| --- | ---: | ---: |
| genes | 107 | 20,097 |
| with LOEUF | 104 | 18,098 |
| median LOEUF | 0.913 | 0.968 |
| LOEUF p25–p75 | 0.652–1.182 | 0.652–1.289 |
| LOEUF &lt; 0.35 | 9 / 104 (0.086538) | 1,395 / 18,098 (0.077080) |
| LOEUF &lt; 0.6 | 23 / 104 (0.221154) | 3,895 / 18,098 (0.215217) |
| with a decile | 104 | 17,840 |
| decile 0 | 12 / 104 (0.115385) | 1,851 / 17,840 (0.103756) |
| deciles 0–2 | 32 / 104 (0.307692) | 5,537 / 17,840 (0.310370) |
| with pLI | 104 | 18,098 |
| median pLI | 0 | 0 |
| pLI ≥ 0.9 | 16 / 104 (0.153846) | 3,158 / 18,098 (0.174494) |
| with missense z | 105 | 18,548 |
| median missense z | 0.838 | 0.750 |
| missense z ≥ 3.09 | 11 / 105 (0.104762) | 1,550 / 18,548 (0.083567) |
| with tau | 106 | 18,867 |
| median tau | 0.881 | 0.810 |
| tau p25–p75 | 0.745–0.963 | 0.637–0.942 |
| tau ≥ 0.8 | 66 / 106 (0.622642) | 9,724 / 18,867 (0.515397) |
| median tissues ≥1 TPM | 41.5 | 52.0 |
| median max median TPM | 59.31 | 43.21 |

Background decile rates use the 17,840 protein-coding genes with a decile, not the 18,098 with a LOEUF. Missense-z rates use genes with a missense z, which is not the same set as genes with a LOEUF.

The germline protein-coding slice is smaller and less constrained: 71 genes, 70 with a LOEUF, median LOEUF 0.993 (p25–p75 0.826–1.235). LOEUF &lt; 0.35 is 2 / 70 (0.028571). LOEUF &lt; 0.6 is 8 / 70 (0.114286). Decile 0 is 5 / 70 (0.071429). Deciles 0–2 are 12 / 70 (0.171429). pLI ≥ 0.9 is 7 / 70 (0.100000). Median missense z is 0.562, and missense z ≥ 3.09 is 5 / 71 (0.070423). All 71 have a tau; median tau is 0.899 (p25–p75 0.773–0.982); tau ≥ 0.8 is 51 / 71 (0.718310). Median tissues at or above 1 TPM is 34.0. Median peak median-TPM is 55.35. CYP2C19, the only practice-guideline gene here, has LOEUF 1.668 (decile 8), missense z −2.28, and a liver maximum of 13.25 TPM in 5 tissues (tau 0.980). RYR1's 38 alleles are all germline and all expert-panel; its LOEUF is 0.644 (decile 2), missense z 4.01, max tissue skeletal muscle, tau 0.994.

Genes with at least 10 drug-response alleles. Counts are alleles carrying that symbol. Stars ≥3 is expert panel or practice guideline. Star 4 is practice guideline and, in this set, only CYP2C19.

| symbol | alleles | germline | somatic | mixed | other | stars ≥3 | star 4 | LOEUF | decile | pLI | mis.z | max tissue | max TPM | tissues ≥1 | tau |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| ABCB1 | 598 | 0 | 595 | 3 | 0 | 0 | 0 | 0.575 | 1 | 0 | 2.77 | Adrenal Gland | 155.95 | 48 | 0.963 |
| OPRM1 | 461 | 0 | 459 | 2 | 0 | 0 | 0 | 1.082 | 5 | 0 | −0.03 | Testis | 2.53 | 2 | 0.979 |
| CYP2D6 | 236 | 2 | 225 | 4 | 5 | 0 | 0 | 1.706 | 8 | 0 | −1.32 | Liver | 215.32 | 48 | 0.986 |
| COMT | 161 | 0 | 161 | 0 | 0 | 0 | 0 | 1.128 | 6 | 0 | 0.74 | Nerve - Tibial | 74.70 | 54 | 0.508 |
| UGT2B7 | 127 | 0 | 127 | 0 | 0 | 0 | 0 | 1.203 | 6 | 0 | −1.18 | Liver | 127.37 | 8 | 0.973 |
| EGFR | 56 | 0 | 48 | 4 | 4 | 2 | 0 | 0.475 | 1 | 0.983 | 3.06 | Skin - Sun Exposed (Lower leg) | 78.34 | 52 | 0.770 |
| RYR1 | 38 | 38 | 0 | 0 | 0 | 38 | 0 | 0.644 | 2 | 0 | 4.01 | Muscle - Skeletal | 423.54 | 22 | 0.994 |
| TRPC6 | 35 | 0 | 35 | 0 | 0 | 0 | 0 | 0.708 | 2 | 0 | 1.91 | Lung | 23.21 | 31 | 0.891 |
| CYP2C19 | 31 | 31 | 0 | 0 | 0 | 30 | 29 | 1.668 | 8 | 0 | −2.28 | Liver | 13.25 | 5 | 0.980 |
| CFTR | 29 | 29 | 0 | 0 | 0 | 29 | 0 | 1.152 | 6 | 0 | −1.07 | Pancreas | 65.79 | 9 | 0.986 |
| DPYD | 20 | 20 | 0 | 0 | 0 | 20 | 0 | 1.090 | 5 | 0 | 0.10 | Cells - Cultured fibroblasts | 43.02 | 52 | 0.745 |
| CYP2C9 | 11 | 11 | 0 | 0 | 0 | 0 | 0 | 1.607 | 8 | 0 | −2.00 | Liver | 433.76 | 7 | 0.998 |
| MYO1E | 11 | 0 | 11 | 0 | 0 | 0 | 0 | 0.557 | 1 | 0 | 0.72 | Cells - EBV-transformed lymphocytes | 52.22 | 52 | 0.770 |
| VKORC1 | 10 | 9 | 0 | 0 | 1 | 7 | 0 | 1.172 | 6 | 0.001 | −0.10 | Cells - Cultured fibroblasts | 306.70 | 54 | 0.669 |

Metric gaps, from a left join back to `gene_constraint` and `gene_expression` on bare ENSG. Six of 110 genes have a null LOEUF. Three of those are protein-coding: DDX3X and MT-ND1 have no constraint row; IFNL4 has a constraint row whose LOEUF is null. CASC6, MT-RNR1, and MT-TS1 also have no constraint row. One gene, H3C2, has no expression row (LOEUF on the constraint row is 1.892). The other 109 matched genes have a tau.

## Caveat

CPIC and PharmGKB are not in this warehouse. Nothing here is a guideline strength, an allele-function assignment, or a dosing recommendation from those catalogs. A large ClinVar row count would still not be that catalog. This count is not large.

The 1,979 rows are not 1,979 independent findings. One phenotype string, at 0 stars and somatic origin, is 0.791814 of the alleles, on five genes. The somatic origin field is not an oncogenicity call: 1,970 alleles have an empty oncogenicity string.

`clinvar_drug_response_n` summed across genes is 1,978 because six alleles match two symbols. Use 1,979 as the allele denominator. The symbol join itself is intact for this set: 0 symbols land on two GENCODE genes, and the recomputed count matches the mart column on every matched symbol. Seven alleles stay out of the gene table because they have no gene token or only a dosage-sensitivity phrase.

LOEUF and tau describe the genes, not the quality of each allele. ABCB1's LOEUF does not review its 598 alleles. Those are 0-star. Medians are across genes, unweighted by allele count, so the tramadol genes do not dominate the constraint and expression summaries. Means of maximum TPM are not used.
