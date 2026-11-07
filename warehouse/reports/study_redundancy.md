# Study redundancy

Numbers are from `warehouse/reports/sql/study_redundancy.sql`, read against `warehouse/gwas/associations.parquet`, `warehouse/gwas/loci.parquet`, and `warehouse/gwas/traits.parquet`. The SQL does not rewrite those tables.

## Question

How far do 1,205,644 association rows shrink after trait-and-distance clumping, and which traits are dominated by repeated studies of the same locus?

## Filters

Clumps are the rows already stored in `loci.parquet`. Within one `MAPPED_TRAIT_URI` and one chromosome, an eligible association joins the lead with the smallest `p_value_float`. Ties break by smaller position, then rsID, then file line. Every other still-unassigned association within 1,000,000 bp of that lead joins it. This is not single-linkage: a row within 1,000,000 bp of a member but farther than that from the lead stays out. `start` and `end` are the minimum and maximum member positions, not the search window.

An association is eligible when `chrom` is canonical, `pos` is an integer of at least 1, and `CNV` is not `Y`. A null position stays in `associations.parquet` and stays out of `loci.parquet`. A row with more than one trait URI is clumped once per URI, so `sum(n_associations)` can exceed the eligible row count.

Trait names come from `traits.parquet` on `trait_uri`, including a null URI. That join misses nothing: loci with no trait row are 0, and no trait has more clump memberships than catalog rows. 1,777 trait URIs have associations and no locus. Their trait-table counts sum to 13,479. Those associations had no eligible position. Leads on the same trait and chromosome are at least 1,000,015 bp apart. No locus has a span above 2,000,000 bp, a lead outside its span, or `n_associations` below 1.

For the trait list, memberships are `sum(n_associations)` on that URI. `k50` is the smallest number of loci, largest `n_associations` first, that reach half of those memberships. Ties break by chromosome, lead position, then lead rsID, nulls first. No named trait had two loci tied for the largest count. A trait is named when it has at least 500 memberships and `k50` is 5 or less, so half of the rows is at least 250. The floor is a reporting cut, not a property of the clump.

Study, PubMed, position, and author-trait counts are computed only for each named trait's largest locus. The count is the eligible rows for that URI on that chromosome with `pos` inside 1,000,000 bp of the lead. For every named trait that window count equals the locus `n_associations` (`bad_windows` is 0), so the window is exactly that clump.

## Denominator

| Count | Value |
| --- | ---: |
| Association rows | 1,205,644 |
| Null `pos` (and null `chrom`) | 216,880 |
| `pos` below 1 | 0 |
| `CNV` equal to `Y` | 0 |
| Eligible rows | 988,764 |
| Loci | 311,160 |
| Sum of locus `n_associations` | 1,083,056 |
| Rows with one trait URI | 1,105,561 |
| Rows with more than one trait URI | 100,074 |
| Rows with an empty trait URI | 9 |
| Eligible rows with more than one trait URI | 81,387 |
| Extra memberships from those rows | 94,292 |
| Eligible URI slots, including 8 empty-URI rows | 1,083,056 |
| Unplaced URI slots | 235,766 |
| Trait-table URI counts, placed and unplaced | 1,318,822 |

The 216,880 null-position rows are the only associations excluded from clumping. 988,764 + 216,880 = 1,205,644. The membership sum is the eligible rows plus 94,292 extra URI copies: 988,764 + 94,292 = 1,083,056.

Within one trait, the denominator for "most of the rows" is that trait's clump memberships, not its catalog-row count. Catalog rows also include positions that never entered a clump. Those unplaced counts are shown beside the named traits.

## Result

Eligible rows shrink to 311,160 loci. That is 3.177671 eligible rows per locus, and the loci are 0.314696 of the eligible rows. Setting 1,205,644 raw rows against 311,160 loci gives 3.874675 raw rows per locus, and the loci are 0.258086 of the raw rows. That second ratio mixes a merge with a drop: 0.179887 of raw rows have a null position and are not represented by any locus.

Counted the way the clumper counts, 1,083,056 memberships sit in 311,160 loci, or 3.480704 memberships per locus. The largest single locus has 4,303 rows.

| Rows in the locus | Loci | Memberships | Share of loci | Share of memberships |
| --- | ---: | ---: | ---: | ---: |
| 1 | 181,401 | 181,401 | 0.582983 | 0.167490 |
| 2–4 | 86,545 | 222,942 | 0.278137 | 0.205845 |
| 5–9 | 25,281 | 162,872 | 0.081248 | 0.150382 |
| 10–49 | 16,165 | 292,621 | 0.051951 | 0.270181 |
| 50–99 | 1,196 | 80,359 | 0.003844 | 0.074197 |
| 100–499 | 524 | 91,219 | 0.001684 | 0.084224 |
| 500 or more | 48 | 51,642 | 0.000154 | 0.047682 |

Most loci are one catalog row. Most rows are not. The 181,401 singleton loci are 0.582983 of loci and 0.167490 of memberships. Loci with 10 or more rows are 17,933 loci (0.057633) and 515,841 memberships (0.476283).

Concentration inside a trait is a different cut. Of 15,109 trait URIs that have a locus, 1,015 have at least 100 memberships and 288 have at least 500. Among those 288, the median `k50` is 53. Eighteen have `k50` of 5 or less (0.062500 of the 288). Another 102 traits with 100 to 499 memberships also have `k50` of 5 or less. They are not listed. The 18 named traits hold 39,593 memberships, 0.036557 of all memberships. Their largest loci hold 17,118 of those rows (0.432349). Their five largest loci hold 29,394 (0.742404).

Three traits need only one locus to pass half: milk amount, lipoprotein A measurement, and amount of phosphatidylcholine.

| Trait | Memberships | Unplaced | Loci | k50 | Largest locus | Share | Top 5 rows | Top 5 share | Lead |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| milk amount | 704 | 0 | 76 | 1 | 538 | 0.7642 | 614 | 0.8722 | chr19:48702851 rs679574 |
| lipoprotein A measurement | 613 | 63 | 56 | 1 | 461 | 0.7520 | 535 | 0.8728 | chr6:160344023 rs572889 |
| amount of phosphatidylcholine | 5,974 | 237 | 219 | 1 | 4,303 | 0.7203 | 5,139 | 0.8602 | chr11:61803311 rs174547 |
| glycerophospholipid measurement | 515 | 35 | 23 | 2 | 256 | 0.4971 | 481 | 0.9340 | chr4:114483120 rs10002657 |
| eye measurement | 1,088 | 18 | 156 | 2 | 529 | 0.4862 | 759 | 0.6976 | chr15:28285036 rs1667394 |
| amount of phosphatidylinositol | 1,399 | 110 | 44 | 2 | 664 | 0.4746 | 1,329 | 0.9500 | chr11:61820833 rs174564 |
| degree of unsaturation measurement | 2,324 | 110 | 146 | 2 | 1,099 | 0.4729 | 1,554 | 0.6687 | chr11:61867955 rs181479770 |
| docosahexaenoic acid to total fatty acids percentage | 1,606 | 89 | 101 | 2 | 754 | 0.4695 | 1,079 | 0.6719 | chr11:61820833 rs174564 |
| amount of phosphatidylethanolamine | 2,976 | 36 | 93 | 2 | 1,361 | 0.4573 | 2,785 | 0.9358 | chr15:58431476 rs1800588 |
| memory impairment | 3,536 | 38 | 243 | 2 | 1,594 | 0.4508 | 2,256 | 0.6380 | chr19:44903416 rs10119 |
| amount of diglyceride | 4,715 | 0 | 54 | 2 | 1,586 | 0.3364 | 3,851 | 0.8168 | chr8:11114773 rs4840546 |
| docosahexaenoic acid measurement | 2,260 | 133 | 163 | 3 | 846 | 0.3743 | 1,366 | 0.6044 | chr11:61776027 rs174528 |
| hyperuricemia | 717 | 4 | 45 | 3 | 186 | 0.2594 | 535 | 0.7462 | chr11:64593747 rs121907892 |
| amount of ceramide | 1,752 | 6 | 125 | 3 | 438 | 0.2500 | 1,268 | 0.7237 | chr4:47543124 rs35818294 |
| sphingomyelin measurement | 2,785 | 108 | 170 | 3 | 623 | 0.2237 | 2,035 | 0.7307 | chr20:12982070 rs364585 |
| omega-3 polyunsaturated fatty acid measurement | 3,381 | 163 | 218 | 4 | 1,181 | 0.3493 | 1,892 | 0.5596 | chr11:61776027 rs174528 |
| Inguinal hernia | 2,221 | 160 | 73 | 4 | 410 | 0.1846 | 1,357 | 0.6110 | chr8:25848843 rs4410916 |
| amygdala volume | 1,027 | 55 | 116 | 5 | 289 | 0.2814 | 559 | 0.5443 | chr4:102267552 rs13107325 |

Trait URIs, in the same order: `http://purl.obolibrary.org/obo/OBA_VT0015043`, `http://www.ebi.ac.uk/efo/EFO_0006925`, `http://purl.obolibrary.org/obo/OBA_2040179`, `http://www.ebi.ac.uk/efo/EFO_0007630`, `http://www.ebi.ac.uk/efo/EFO_0004731`, `http://purl.obolibrary.org/obo/OBA_2040181`, `http://www.ebi.ac.uk/efo/EFO_0022261`, `http://www.ebi.ac.uk/efo/EFO_0022262`, `http://purl.obolibrary.org/obo/OBA_2040180`, `http://www.ebi.ac.uk/efo/EFO_0001072`, `http://purl.obolibrary.org/obo/OBA_2040178`, `http://www.ebi.ac.uk/efo/EFO_0007761`, `http://www.ebi.ac.uk/efo/EFO_0009104`, `http://purl.obolibrary.org/obo/OBA_2040177`, `http://www.ebi.ac.uk/efo/EFO_0010118`, `http://www.ebi.ac.uk/efo/EFO_0010119`, `http://purl.obolibrary.org/obo/HP_0000023`, `http://purl.obolibrary.org/obo/OBA_2050034`.

These are not the largest traits. Among traits with at least 5,000 memberships, amount of phosphatidylcholine is the only one with `k50` of 5 or less. The next smallest `k50` values in that set are fatty acid amount (7), memory performance (11), and triglyceride measurement (12). Body height is the largest trait and is not concentrated. Gut microbiome measurement is the most diffuse trait in that set.

| Trait | Memberships | Loci | Rows per locus | k50 | Largest locus | Share | Top 5 share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| amount of phosphatidylcholine | 5,974 | 219 | 27.28 | 1 | 4,303 | 0.7203 | 0.8602 |
| fatty acid amount | 18,611 | 419 | 44.42 | 7 | 4,186 | 0.2249 | 0.4573 |
| memory performance | 9,341 | 969 | 9.64 | 11 | 3,109 | 0.3328 | 0.4472 |
| triglyceride measurement | 19,045 | 869 | 21.92 | 12 | 3,260 | 0.1712 | 0.3713 |
| Alzheimer disease | 5,461 | 1,208 | 4.52 | 82 | 1,035 | 0.1895 | 0.2441 |
| body mass index | 19,849 | 1,334 | 14.88 | 160 | 276 | 0.0139 | 0.0560 |
| body height | 38,386 | 1,752 | 21.91 | 218 | 312 | 0.0081 | 0.0295 |
| gut microbiome measurement | 10,921 | 1,691 | 6.46 | 405 | 120 | 0.0110 | 0.0242 |

Height has 21.91 rows per locus, close to the phosphatidylcholine rate of 27.28, but its largest locus is 312 of 38,386 rows and half of its rows take 218 loci. Fatty acid amount has the highest rows-per-locus figure in this contrast, 44.42. Its largest locus has 4,186 rows at chr11:61854139 (rs73487492). The only larger locus in the file is the phosphatidylcholine locus of 4,303. Seven loci still fall short of the rule used above, and the top five hold 0.4573 of that trait, just under half. Alzheimer disease has a locus of 1,035 rows and still needs 82 loci to cover half of its 5,461 memberships. A large locus is not the same thing as a trait dominated by a few loci.

The largest locus of each named trait is not a pile of independent replications. `STUDY ACCESSION` counts catalog accessions. One PubMed id often has many accessions. No counted row had a blank accession or a blank PubMed id. The lead mapped-gene text is `mapped_genes.raw` on that locus, the lead association's `MAPPED_GENE` string.

| Trait | Rows | Accessions | PubMed ids | Positions | Author strings in the locus | Lead `MAPPED_GENE` |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| milk amount | 538 | 25 | 2 | 47 | 20 | `FUT2` |
| lipoprotein A measurement | 461 | 20 | 18 | 412 | 10 | `SLC22A2 - SLC22A3` |
| amount of phosphatidylcholine | 4,303 | 386 | 11 | 224 | 320 | `FADS2, FADS1` |
| glycerophospholipid measurement | 256 | 1 | 1 | 256 | 1 | `ARSJ - UGT8` |
| eye measurement | 529 | 217 | 4 | 18 | 215 | `HERC2` |
| amount of phosphatidylinositol | 664 | 55 | 4 | 104 | 44 | `FADS2, FADS1` |
| degree of unsaturation measurement | 1,099 | 25 | 3 | 326 | 2 | `FADS2 - FADS3` |
| docosahexaenoic acid to total fatty acids percentage | 754 | 22 | 5 | 205 | 3 | `FADS2, FADS1` |
| amount of phosphatidylethanolamine | 1,361 | 153 | 6 | 174 | 108 | `LIPC, ALDH1A2` |
| memory impairment | 1,594 | 23 | 3 | 311 | 9 | `TOMM40` |
| amount of diglyceride | 1,586 | 21 | 1 | 318 | 21 | `XKR6` |
| docosahexaenoic acid measurement | 846 | 35 | 11 | 232 | 12 | `MYRF, TMEM258` |
| hyperuricemia | 186 | 30 | 2 | 24 | 29 | `SLC22A12` |
| amount of ceramide | 438 | 47 | 3 | 88 | 35 | `ATP10D` |
| sphingomyelin measurement | 623 | 76 | 9 | 136 | 63 | `LINC01723` |
| omega-3 polyunsaturated fatty acid measurement | 1,181 | 42 | 11 | 353 | 14 | `MYRF, TMEM258` |
| Inguinal hernia | 410 | 10 | 6 | 400 | 2 | `EBF2` |
| amygdala volume | 289 | 69 | 2 | 87 | 22 | `SLC39A8` |

Glycerophospholipid measurement's largest locus is one study, one paper, and 256 positions: one SNP list, not a repeated study. Amount of diglyceride's largest locus is one paper (21 accessions, 318 positions, 1,586 rows). Inguinal hernia's largest locus has 400 positions and 6 papers. Lipoprotein A measurement has 412 positions and 18 papers, so the rows are mostly distinct sites from a modest set of papers. Eye measurement is the clearest repeated-accession case in the list: 529 rows, 18 positions, 217 accessions, and 4 papers, with 215 author trait strings in the locus. Amount of phosphatidylcholine, the biggest locus in the file, is 4,303 rows from 11 papers, 386 accessions, 224 positions, and 320 author trait strings.

Six of the eighteen largest loci are on chromosome 11 between 61,776,027 and 61,867,955. Fatty acid amount's 4,186-row locus is in that same interval and is not in the named set. These are separate trait URIs, so the same papers can be counted once under each URI. Clumping does not merge across traits. Two of the six are stored as `MYRF, TMEM258` and the other four as `FADS2, FADS1` or `FADS2 - FADS3`. The strings are not the same, and this report does not decide which gene they belong to.

## Caveat

This is catalog redundancy, not biological importance. A locus with thousands of rows is a place where the file repeats itself. It is not a ranking of genes, traits, or mechanisms.

`n_associations` counts catalog rows inside a 1,000,000 bp lead window. It does not count distinct papers, distinct accessions, or distinct positions. Those three counts come apart in the table above. A study accession is not a publication. One GWAS that reports many SNPs in a window, or many accessions for related metabolite labels, produces the same shrinkage as many papers reporting one SNP.

The 1,000,000 bp rule is distance, not linkage disequilibrium. A window can hold more than one signal, and a signal wider than the window can be split. Leads on one trait and chromosome are forced more than 1,000,000 bp apart.

The lead `MAPPED_GENE` text is not a fine-mapped gene. `FADS2 - FADS3` and `SLC22A2 - SLC22A3` contain no comma, so each stays one catalog string. Lipoprotein A measurement's largest locus is labeled `SLC22A2 - SLC22A3` in the file. That label is reported as stored.

Rows with no position, 216,880 of them, never enter the shrinkage. Traits that look small in the locus file may still have a large unplaced count. The named-trait shares use clump memberships only.

The named list is the large end of the concentrated traits. 102 traits with 100 to 499 memberships also put half their rows in five or fewer loci. The traits that hold most of the catalog, including height, body mass index, and Alzheimer disease, do not.
