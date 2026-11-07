# Dark genome

This is a count of bases, not a claim that the uncounted bases have no function.

Contig lengths were read from the `##contig` header of `D:\variant-workbench\raw\gnomad\genomes\gnomad.genomes.v4.1.sites.chr21.vcf.bgz` by `scripts/vcf_stream.py`.
The table is not hard-coded. The header lists chr1–chr22, chrX, and chrY. It has no chrM contig.
Denominator: 3088269832 bases.

Each category is the union of its own intervals, so overlap inside a category is counted once.
A base that falls in two categories is counted in both. The categories are not an exclusive partition.
The remainder is the bases that are in none of the four categories.

Definitions:

- Protein-coding CDS: `warehouse/intervals/cds.parquet` rows with `feature = CDS` on a gene whose biotype is `protein_coding`. `stop_codon` rows are not included.
- Gene span: `start`–`end` of every GENCODE gene on these contigs, any biotype.
- Tier A neighborhood: inclusive `[pos - 1000, pos + 1000]` around each key in `warehouse/variants/tier_a_keys.tsv`, clipped to the contig.
- GWAS clump: `start`–`end` in `warehouse/gwas/loci.parquet`. Those columns are the minimum and maximum member positions. They are not expanded to the ±1,000,000 bp search radius used to build the clump.

Loci whose start equals end (single association): 202104.
CDS intervals outside the contig set, excluded: 13.
Gene spans outside the contig set, excluded: 37.
GWAS intervals outside the contig set, excluded: 0.
Tier A keys outside the contig set, excluded: 3136.

| chrom | length | protein_coding_cds | gene_span | tier_a_1kb | gwas_clump | in_any | in_none |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chr1 | 248956422 | 3692890 | 170564304 | 129486997 | 223855916 | 227020823 | 21935599 |
| chr2 | 242193529 | 2658580 | 177118711 | 110003022 | 235963619 | 238183713 | 4009816 |
| chr3 | 198295559 | 2071450 | 146899747 | 88101788 | 193774361 | 194518256 | 3777303 |
| chr4 | 190214555 | 1437020 | 126986989 | 71033262 | 184536181 | 186221666 | 3992889 |
| chr5 | 181538259 | 1685904 | 126120766 | 76295036 | 175903606 | 177482682 | 4055577 |
| chr6 | 170805979 | 1825511 | 121753166 | 77243539 | 167603120 | 168444852 | 2361127 |
| chr7 | 159345973 | 1700804 | 113390601 | 73451877 | 154718242 | 155449261 | 3896712 |
| chr8 | 145138636 | 1244481 | 104015875 | 61273338 | 141965069 | 142864995 | 2273641 |
| chr9 | 138394717 | 1473413 | 88434243 | 56866481 | 109327176 | 117253403 | 21141314 |
| chr10 | 133797422 | 1422076 | 94130495 | 64220101 | 130866190 | 131234902 | 2562520 |
| chr11 | 135086622 | 2154852 | 95430824 | 68347957 | 131094289 | 131201055 | 3885567 |
| chr12 | 133275309 | 1910471 | 97254264 | 72458828 | 130214442 | 130560481 | 2714828 |
| chr13 | 114364328 | 659045 | 61873150 | 34386347 | 94779039 | 95607863 | 18756465 |
| chr14 | 107043718 | 1143307 | 65867775 | 46312423 | 86829582 | 87917164 | 19126554 |
| chr15 | 101991189 | 1253436 | 67782683 | 48211424 | 79063552 | 81048570 | 20942619 |
| chr16 | 90338345 | 1551473 | 60525149 | 46408353 | 77680408 | 78804233 | 11534112 |
| chr17 | 83257441 | 2093662 | 62646341 | 55815272 | 78916483 | 79230815 | 4026626 |
| chr18 | 80373285 | 567127 | 50575450 | 30396510 | 74381924 | 74686649 | 5686636 |
| chr19 | 58617616 | 2417239 | 43133746 | 46711405 | 55502997 | 55686031 | 2931585 |
| chr20 | 64444167 | 872672 | 42540548 | 35187014 | 60135582 | 60827471 | 3616696 |
| chr21 | 46709983 | 350140 | 25813629 | 16387396 | 33606329 | 36057514 | 10652469 |
| chr22 | 50818468 | 761342 | 27632143 | 24430787 | 34491393 | 35540028 | 15278440 |
| chrX | 156040895 | 1380227 | 83309230 | 45345544 | 79393453 | 122783577 | 33257318 |
| chrY | 57227415 | 94950 | 10430907 | 498809 | 1729558 | 11450713 | 45776702 |
| total | 3088269832 | 36422072 | 2064230736 | 1378873510 | 2736332511 | 2820076717 | 268193115 |

| set | bases | fraction of denominator |
| --- | ---: | ---: |
| protein_coding_cds | 36422072 | 0.011794 |
| gene_span | 2064230736 | 0.668410 |
| tier_a_1kb | 1378873510 | 0.446487 |
| gwas_clump | 2736332511 | 0.886041 |
| in_any | 2820076717 | 0.913157 |
| in_none | 268193115 | 0.086843 |

The four category fractions are not exclusive, so they do not sum to the in_any fraction.
in_none is the exclusive remainder.
