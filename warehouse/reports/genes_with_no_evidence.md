# Genes with no evidence

This is the boundary of what the warehouse has joined onto a GENCODE gene. It is not a claim that the gene does nothing.

Numbers come from `warehouse/reports/sql/genes_with_no_evidence.sql`, read-only against `warehouse/marts/gene_mart.parquet`. The mart was not rewritten. Each source join was recounted. Mismatch counts were 0.

## Question

How many genes have no expression measurement, no constraint call, no ClinVar submission, no GWAS locus by either rule, and no MaveDB assay?

## Filters

A gene is counted only when every predicate below is true. The same predicates are applied to the protein-coding subset.

No expression measurement means the gene is absent from GTEx. `max_median_tpm`, `max_tissue`, `tau`, `n_tissues_tpm_ge_1`, and `n_tissues_tpm_ge_10` are all null together. That is a left-join miss to `warehouse/gtex/gene_expression.parquet` on bare ENSG. Tau present is the opposite state: the maximum median TPM is above 0, so tau is defined. Those genes are excluded, including genes whose median TPM stays below 1 in every tissue.

A gene GTEx tested with all-zero expression is not this state. For those, `max_median_tpm` is 0, `n_tissues_tpm_ge_1` and `n_tissues_tpm_ge_10` are 0, and `max_tissue` and `tau` are null because tau divides by the maximum. There are 11,784 such genes, 342 of them protein-coding. They stay out of the headline count. Counting them as "no expression" would raise the intersection from 908 to 2,137, and the protein-coding intersection from 10 to 28.

The three expression states partition the mart: 24,065 absent, 11,784 tested and all-zero, 42,842 with tau present. Among protein-coding genes the split is 888, 342, and 18,867.

No constraint call means the constraint row is missing. `loeuf`, `loeuf_decile`, `lof.pLI`, `mis.z_score`, `syn.z_score`, and `constraint_flags` are all null. That is a left-join miss to `warehouse/constraint/gene_constraint.parquet` on bare ENSG. It was not filled in by symbol or by MANE transcript.

A constraint row with a null LOEUF is a different state. 492 genes, 469 of them protein-coding, have `constraint_flags` set and a null `loeuf`. Most of those still have a missense or synonymous z-score. They are not counted as "no constraint call." Widening the filter from "row missing" to "LOEUF null" does not change this intersection: it stays 908 and 10, because none of those 492 pass the other filters. Constraint states partition the mart: 60,078 rows missing, 492 rows with a null LOEUF, 18,121 with a LOEUF. Protein-coding: 1,530, 469, and 18,098.

No ClinVar submission means `clinvar_n` is 0 or null. On this mart the column is never null. It is 0 for 57,239 genes and for 1,056 protein-coding genes. The count is an exact match of a ClinVar `GeneSymbol` token, split on `;`, to the GENCODE symbol. It is not an ENSG match. When `clinvar_n` is 0, the pathogenic, likely-pathogenic, VUS, benign, likely-benign, conflicting, and drug-response counts are 0 as well.

No GWAS locus by either rule means `gwas_loci_by_overlap` is 0 or null and `gwas_loci_by_mapped_gene` is 0 or null. Neither column is null. Both are 0 for 3,458 genes and for 380 protein-coding genes. Overlap requires the same chromosome and either an intersection of the locus member span with the gene span, or the lead position inside the gene. The span is the minimum and maximum member positions, not the ±1 Mb window used to build the clump. The mapped-gene rule is an exact token match to the symbol and does not require the locus to sit on the gene.

No MaveDB assay means all four counts are 0 or null: `mavedb_experiments_by_symbol`, `mavedb_mapped_variants_by_symbol`, `mavedb_variant_keys_in_span`, and `mavedb_score_sets_in_span`. None of the four is null. All four are 0 for 78,490 genes and for 19,923 protein-coding genes. The symbol counts are exact matches to `target_gene`, summed from `gene_coverage.parquet`. `mavedb_experiments_by_symbol` is the sum of `n_score_sets`, not a count of experiment URNs. The span counts are distinct scored variant keys and score-set URNs whose GRCh38 position falls inside the gene.

`eqtl_tissues` and `sqtl_tissues` are not filters.

## Denominator

`gene_mart.parquet` has 78,691 rows, one per GENCODE gene, and 78,691 distinct ENSG ids. The protein-coding subset is 20,097 genes, `biotype = 'protein_coding'`.

## Result

908 genes pass all five filters. That is 0.011539 of 78,691.

10 of those are protein-coding. That is 0.000498 of 20,097.

All 908 have `eqtl_tissues` null and `sqtl_tissues` null. None are in the tested-and-negative state, which would be a zero rather than a null. QTL status was not required.

| filters required | genes | protein-coding |
| --- | ---: | ---: |
| all five | 908 | 10 |
| the other four, expression absence not required | 2977 | 52 |
| the other four, constraint-row absence not required | 908 | 10 |
| the other four, ClinVar zero not required | 938 | 21 |
| the other four, GWAS zero not required | 23199 | 174 |
| the other four, MaveDB zero not required | 908 | 10 |
| all five, but all-zero TPM counted as no expression | 2137 | 28 |
| all five, but any null LOEUF counted as no constraint call | 908 | 10 |

Twelve biotypes account for the 908. Every other biotype in the mart contributes 0.

| biotype | in the set | biotype denominator | fraction of that biotype |
| --- | ---: | ---: | ---: |
| lncRNA | 806 | 34880 | 0.023108 |
| transcribed_unprocessed_pseudogene | 19 | 1587 | 0.011972 |
| processed_pseudogene | 18 | 9487 | 0.001897 |
| unprocessed_pseudogene | 17 | 1949 | 0.008722 |
| TEC | 12 | 1019 | 0.011776 |
| misc_RNA | 10 | 2207 | 0.004531 |
| protein_coding | 10 | 20097 | 0.000498 |
| transcribed_processed_pseudogene | 8 | 1149 | 0.006963 |
| miRNA | 3 | 1879 | 0.001597 |
| snRNA | 2 | 1901 | 0.001052 |
| transcribed_unitary_pseudogene | 2 | 201 | 0.009950 |
| rRNA_pseudogene | 1 | 497 | 0.002012 |

The protein-coding table is the whole protein-coding result, not a sample. Coordinates are 1-based and inclusive, as stored on the mart. A blank `mane_enst` is null.

| ensg | symbol | chrom | start | end | strand | span bp | mane_enst | symbol_n_gencode |
| --- | --- | --- | ---: | ---: | --- | ---: | --- | ---: |
| ENSG00000284733 | OR4F29 | chr1 | 450740 | 451678 | - | 939 | ENST00000426406 | 1 |
| ENSG00000293691 | ENSG00000293691 | chr5 | 70487719 | 70491228 | + | 3510 | | 1 |
| ENSG00000293689 | ENSG00000293689 | chr5 | 70508196 | 70512509 | + | 4314 | | 1 |
| ENSG00000291314 | ENSG00000291314 | chrX | 10566888 | 10576955 | - | 10068 | | 1 |
| ENSG00000288706 | ENSG00000288706 | chrX | 23783278 | 23783367 | + | 90 | ENST00000683890 | 1 |
| ENSG00000293661 | ENSG00000293661 | chrX | 135248590 | 135252011 | - | 3422 | | 1 |
| ENSG00000293662 | ENSG00000293662 | chrX | 135305788 | 135306402 | - | 615 | | 1 |
| ENSG00000293663 | ENSG00000293663 | chrX | 135391849 | 135392207 | - | 359 | | 1 |
| ENSG00000287585 | ENSG00000287585 | chrX | 149414877 | 149415506 | - | 630 | ENST00000653668 | 1 |
| ENSG00000292372 | WASH6P | chrY | 57207346 | 57212230 | + | 4885 | | 2 |

One gene from each biotype in the set, the first symbol in that biotype and then the ENSG. Each of these rows has the five null and zero states above.

| biotype | ensg | symbol | chrom | start | end | span bp |
| --- | --- | --- | --- | ---: | ---: | ---: |
| TEC | ENSG00000286012 | ENSG00000286012 | chr21 | 8240736 | 8241047 | 312 |
| lncRNA | ENSG00000291031 | BCORP1 | chrY | 19455431 | 19567280 | 111850 |
| miRNA | ENSG00000278791 | ENSG00000278791 | chr1 | 632325 | 632413 | 89 |
| misc_RNA | ENSG00000273974 | ENSG00000273974 | chrX | 74280402 | 74280495 | 94 |
| processed_pseudogene | ENSG00000292364 | AMD1P2 | chrY | 57015105 | 57016096 | 992 |
| protein_coding | ENSG00000287585 | ENSG00000287585 | chrX | 149414877 | 149415506 | 630 |
| rRNA_pseudogene | ENSG00000292356 | RNA5SP498 | chrY | 1300256 | 1300375 | 120 |
| snRNA | ENSG00000207332 | RNU6-146P | chrX | 3607258 | 3607365 | 108 |
| transcribed_processed_pseudogene | ENSG00000278647 | CDK2AP2P1 | chr9 | 41280460 | 41280787 | 328 |
| transcribed_unitary_pseudogene | ENSG00000310578 | ENSG00000310578 | chrX | 135330476 | 135330703 | 228 |
| transcribed_unprocessed_pseudogene | ENSG00000286156 | BMS1P23 | chr2 | 94196057 | 94203287 | 7231 |
| unprocessed_pseudogene | ENSG00000292371 | DDX11L16 | chrY | 57212184 | 57214397 | 2214 |

## Caveat

The headline definition of "no expression" is a missing GTEx row, not a measured zero. The headline definition of "no constraint call" is a missing constraint row, not a null LOEUF. ClinVar, GWAS, and MaveDB use zero because those mart columns are coalesced and contain no nulls.

GWAS is what makes the intersection small. Dropping only that requirement leaves 23,199 genes, 174 of them protein-coding. Dropping the constraint requirement or the MaveDB requirement leaves the headline unchanged at 908 and 10. MaveDB is already zero for 78,490 of 78,691 genes, so an assay filter does not pick out this set. A gene with no GWAS locus is outside every clumped member span and is not equal to a mapped-gene token. That is a statement about catalog placement, not a test of the gene.

816 of the 908 genes use their ENSG id as the symbol, including 747 of 806 lncRNAs and 8 of the 10 protein-coding genes. ClinVar, the mapped-gene rule, and the MaveDB symbol counts compare that string. They would miss a submission stored under another name. The joins that do not use the symbol also miss these genes: no GTEx row, no constraint row, no locus overlap, and no scored MaveDB variant inside the span. 1,002 of 1,137 MaveDB coverage rows have a `target_gene` that is not a GENCODE symbol, so those rows never attach to any gene.

63 genes in the set share a symbol with another GENCODE gene. `WASH6P` is one of them (`symbol_n_gencode` 2). A symbol-based count of zero is the count for that symbol, copied onto each gene that carries it. It does not hide a count on the other gene.

The set is uneven across chromosomes: 236 of 2,958 chrX genes, 139 of 672 chrY genes, and 151 of 3,188 chr9 genes. None of the 37 chrM genes are in it. Gene span in the set runs from 68 bp to 343,736 bp, with median 8,024 bp. Three of the ten protein-coding genes have a MANE transcript on the mart and still have no constraint row.

GTEx v8 and the gnomAD constraint file do not cover every GENCODE v49 gene. Absence from those tables is the age and scope of the source, not a measurement that the gene is silent or unconstrained.
