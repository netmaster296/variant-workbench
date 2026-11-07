# Pleiotropy

Numbers are from `warehouse/reports/sql/pleiotropy.sql`, run with DuckDB against the Parquet files. No mart was rewritten.

## Question

Which genes, and which Tier A variants, are linked to many EFO traits?

A shared gene is not a shared mechanism. A gene named inside a large locus is a weaker link than a shared lead. Both counts are below.

## Filters

Clumped loci are `warehouse/gwas/loci.parquet`. One row is one EFO trait and one lead. The locus span is `end - start`, the distance between the outermost member positions. It is not the ±1,000,000 bp search radius used to build the clump.

- Empty `trait_uri`: 8 loci. They are not an EFO trait and are out of every trait count.
- A gene link is an exact match between a `mapped_genes.tokens` value and a GENCODE symbol in `warehouse/intervals/genes.parquet`. A substring match is not used. `HLA-DQB1` is one token. `FBLN1 - LINC01589` is one token and matches neither symbol.
- A large locus has span >= 1,000,000 bp (10,022 loci). A point locus has span 0 (202,096 loci with a trait). "Span < 1 Mb" means the gene has at least one linking locus below that span on that trait.
- The lead's own mapped-gene field is `mapped_genes.raw`, split on commas. Tokens that appear only on a non-lead member are kept in the gene count and flagged.
- A lead gains a variant key only when `warehouse/clinvar/allele_dedup.parquet` has the same rsID, the same chromosome, and the same position, with a non-null `variant_key`. This is the ledger rule. An rsID match at a different position is not attached (269 loci). `tier_a.source` containing `gwas` is a wider set (non-leads and position mismatches) and is not the set used here.
- The key must be in `warehouse/variants/tier_a.parquet`. Keys with `ref = alt` are not alternate alleles and are dropped from the variant counts.
- A shared lead site is the same chromosome, position, and rsID. It does not require ref or alt. A shared Tier A allele is a `ref <> alt` key from the rule above.
- Trait pairs are unordered. Each pair is stored once.

## Denominator

| set | count |
| --- | ---: |
| Clumped loci | 311,160 |
| Loci with an EFO trait | 311,152 |
| Distinct EFO traits on those loci | 15,108 |
| Point loci (span 0) | 202,096 |
| Loci with span > 0 and < 1 Mb | 99,034 |
| Loci with span >= 1 Mb | 10,022 |
| GENCODE genes with a symbol | 78,691 |
| Distinct GENCODE symbols | 77,078 |
| Distinct mapped-gene tokens | 39,283 |
| Tokens that match a symbol | 22,403 |
| Tokens that do not | 16,880 |
| Symbols linked to at least one trait | 22,403 |
| GENCODE genes carrying those symbols | 23,724 |
| Gene–trait links | 358,827 |
| Traits with at least one matched gene | 14,233 |
| Ledger GWAS rows | 312,875 |
| Ledger rows with no ref/alt key | 290,815 |
| Ledger rows with a four-field key | 22,060 |
| Loci that gained a same-position key | 20,345 |
| Distinct Tier A keys from that join | 5,831 |
| Of those, `ref = alt` | 47 |
| Alternate-allele keys used below | 5,784 |
| Traits with at least one such allele | 6,332 |
| Lead sites (`chrom`, position, rsID) | 187,383 |

Most GWAS loci have no ref/alt key (290,815 of 312,875 ledger GWAS rows). The Catalog stores no alleles. A key is attached later, and only for a lead whose rsID matches a ClinVar allele at the same chromosome and position. Those 5,831 keys are all present in `tier_a.parquet` (0 missing). The join holds. 20,345 loci produced the 22,060 keyed ledger rows; 1,635 of the loci match more than one ClinVar allele, which is why the ledger has more keyed rows than keyed loci.

16,879 of the 16,880 unmatched tokens contain ` - `. They are the catalog's intergenic pair and are not split. One remaining token is the unsplit string `No mapped genes; No mapped genes; HLA-DQA1`. 279 matched symbols belong to more than one GENCODE gene. Trait counts below are per symbol, which is the catalog token. The 23,724 gene count is the same figure as `gene_mart.gwas_loci_by_mapped_gene > 0`.

Of the 47 `ref = alt` keys, 19 rsIDs have no other ClinVar allele at that site. Those 19 leads have no alternate allele in Tier A. The other `ref = alt` keys copy the trait list of a real allele at the same rsID (`rs429358` `chr19:44908684:T:T` has the same 274 traits as `chr19:44908684:T:C`). They are not counted as a second variant.

## Result

"Many" is the tail of the distribution, not a separate biological class. Two traits is ordinary for a mapped gene and above the median for a resolved lead allele.

### Genes

Of 22,403 symbols linked to at least one trait, the median is 7 traits.

| traits | symbols |
| ---: | ---: |
| 1 | 3,242 |
| 2–4 | 5,500 |
| 5–9 | 4,685 |
| 10–19 | 4,147 |
| 20–49 | 3,397 |
| 50–99 | 1,011 |
| >= 100 | 421 |

19,161 symbols are linked to two or more traits. 1,432 are linked to 50 or more, and 1,209 of those are protein-coding. The other biotypes in that tail are lncRNA (205), transcribed pseudogenes (16), and miRNA (2). The lncRNA rows are catalog names, not a claim that the lncRNA is the effector.

For symbols that belong to one GENCODE gene, the lead sits inside that gene for 230,448 of 354,963 gene–trait links. On point loci the lead is inside the gene for 147,871 of 147,873 links. The two exceptions are `XG` (the lead is on chrY at a position inside the numeric span of the chrX gene) and `SPRY3` (the lead is on chrY; the only GENCODE gene of that symbol is on chrX). Wide clumps are different: 50,077 of 358,827 gene–trait links exist only through a locus of span >= 1 Mb, and 126,243 links use a token that is not on the lead's own mapped-gene field.

The top symbols stay at the top if the count is restricted to point loci. `ARHGEF3` has 1,287 traits and 1,237 of those have a point locus. Dropping large windows does not create this list. It does cut some of it: `FADS2` falls from 880 traits to 316 on point loci, `APOE` from 723 to 276, `FADS1` from 675 to 200.

| symbol | ENSG | biotype | traits | span < 1 Mb | point locus | token on the lead |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| ARHGEF3 | ENSG00000163947 | protein_coding | 1,287 | 1,279 | 1,237 | 1,275 |
| ABO | ENSG00000175164 | protein_coding | 1,224 | 1,193 | 687 | 1,133 |
| CFH | ENSG00000000971 | protein_coding | 1,111 | 1,100 | 985 | 1,089 |
| GCKR | ENSG00000084734 | protein_coding | 1,013 | 958 | 719 | 973 |
| JMJD1C | ENSG00000171988 | protein_coding | 973 | 961 | 717 | 914 |
| ATXN2 | ENSG00000204842 | protein_coding | 900 | 819 | 500 | 806 |
| FADS2 | ENSG00000134824 | protein_coding | 880 | 855 | 316 | 742 |
| ZFPM2 | ENSG00000169946 | protein_coding | 802 | 798 | 716 | 800 |
| SARM1 | ENSG00000004139 | protein_coding | 734 | 728 | 643 | 727 |
| APOE | ENSG00000130203 | protein_coding | 723 | 671 | 276 | 521 |
| ZFPM2-AS1 | ENSG00000251003 | lncRNA | 703 | 700 | 643 | 686 |
| FADS1 | ENSG00000149485 | protein_coding | 675 | 651 | 200 | 536 |
| VTN | ENSG00000109072 | protein_coding | 673 | 668 | 592 | 659 |
| TRIB1AL | ENSG00000253111 | lncRNA | 588 | 586 | 261 | 555 |
| PPP1R3B-DT | ENSG00000248538 | lncRNA | 575 | 513 | 312 | 526 |

Among the 148,121 trait labels on symbols with 50 or more traits, 56,322 contain the word "measurement". That is a label check, not an ontology. EFO terms are not independent traits.

### Tier A variants

5,784 alternate-allele keys are leads. The median number of traits is 1. These keys are 5,496 rsIDs. 262 rsIDs have two or three keys, and for every one of those rsIDs the keys carry the same trait set. An extra key is not an extra trait list.

| traits | variant keys |
| ---: | ---: |
| 1 | 3,662 |
| 2–4 | 1,495 |
| 5–9 | 343 |
| 10–19 | 151 |
| 20–49 | 87 |
| >= 50 | 46 |

2,122 keys are linked to two or more traits. 284 are linked to 10 or more.

This is not the list of the most trait-linked leads in the catalog. A lead site has no allele unless ClinVar matched it. Of 187,383 lead sites, 42,963 are linked to two or more traits, and 2,010 of those have an alternate-allele key. Of 2,047 sites linked to 10 or more traits, 268 have such a key. Of 203 sites linked to 50 or more traits, 41 have such a key.

The lead linked to the most traits is `rs1354034` at chr3:56,815,721, inside `ARHGEF3` (chr3:56,727,418–57,079,329). It is the lead for 1,238 traits. It has no ref/alt key, so it is absent from the table below. The gene count for `ARHGEF3` (1,287) is that site plus other loci that also name the gene.

Genome `AF` is the Tier A genome value, shown to 6 decimal places. Call set is genome (G), exome (E), or both.

| traits | rsID | variant key | VEP symbol | consequence | genome AF | call set |
| ---: | --- | --- | --- | --- | ---: | --- |
| 894 | rs1260326 | chr2:27508073:T:C | GCKR | missense_variant&splice_region_variant | 0.670353 | G+E |
| 641 | rs704 | chr17:28367840:G:A | VTN | missense_variant | 0.514958 | G+E |
| 346 | rs13107325 | chr4:102267552:C:T | SLC39A8 | missense_variant | 0.044701 | G+E |
| 285 | rs3184504 | chr12:111446804:T:C | SH2B3 | missense_variant | 0.676306 | G+E |
| 274 | rs429358 | chr19:44908684:T:C | APOE | missense_variant | 0.157359 | G+E |
| 257 | rs28929474 | chr14:94378610:C:T | SERPINA1 | missense_variant | 0.012157 | G+E |
| 233 | rs7412 | chr19:44908822:C:T | APOE | missense_variant | 0.077837 | G+E |
| 208 | rs507666 | chr9:133273983:A:G | ABO | intron_variant | 0.831624 | G |
| 205 | rs1047891 | chr2:210675783:C:A | CPS1 | missense_variant | 0.326239 | G+E |
| 191 | rs11591147 | chr1:55039974:G:T | PCSK9 | missense_variant | 0.012244 | G+E |
| 182 | rs2519093 | chr9:133266456:T:C | ABO | intron_variant | 0.832183 | G |
| 163 | rs1800961 | chr20:44413724:C:T | HNF4A | missense_variant | 0.025100 | G+E |
| 160 | rs738409 | chr22:43928847:C:G | PNPLA3 | missense_variant | 0.228136 | G+E |
| 158 | rs10455872 | chr6:160589086:A:G | LPA | intron_variant | 0.042449 | G |
| 142 | rs12740374 | chr1:109274968:G:T | CELSR2 | 3_prime_UTR_variant | 0.221992 | G+E |

`rs28929474` also matches `chr14:94378610:C:G`. That key has the same 257 traits. It is in the exome call set only (exome AF about 4.8e-6), it is not in the genome call set, and it has no VEP block. It is one of the 262 duplicated rsID trait lists, not a second pleiotropic signal.

5,701 of the 5,784 keys are in the genome call set, the exome call set, or both. 83 are in neither. Among the 5,701, the recorded VEP impact is MODIFIER for 2,848, MODERATE for 1,916, LOW for 762, and HIGH for 62. 113 of those observed keys have no VEP impact. The most common consequence strings are `intron_variant` (2,030) and `missense_variant` (1,872). A keyed lead is often not a coding change.

### Trait–trait table

There are 6,907,749 trait pairs that share at least one matched gene. The median pair shares 1 gene. 484,616 pairs share at least 5 genes, 183,933 share at least 10, and 29,363 share at least 50. The maximum is 2,086.

The same pairs counted with a tighter gene rule, and counted again by the lead:

| shared link | trait pairs |
| --- | ---: |
| Any matched gene | 6,907,749 |
| A gene whose loci on both traits have span < 1 Mb | 6,669,751 |
| A gene whose loci on both traits are point loci | 4,087,595 |
| Genes only through a locus of span >= 1 Mb | 237,998 |
| A lead site (chromosome, position, rsID; alleles not required) | 2,417,256 |
| A gene and a lead site | 2,351,085 |
| A Tier A alternate allele | 938,863 |
| A gene and a Tier A alternate allele | 936,424 |
| A Tier A alternate allele and no matched gene | 2,439 |
| A lead site and no matched gene | 66,171 |

4,556,664 pairs share a mapped gene and do not share a lead site. That is the weak overlap: different leads, overlapping gene labels. 936,424 of 6,907,749 gene-sharing pairs (0.136) also share a Tier A allele. Absence of a Tier A allele is not absence of a shared lead. 1,478,393 of the 2,417,256 site-sharing pairs have no resolved alternate allele. Most GWAS loci have no ref/alt key (290,815 of 312,875 ledger GWAS rows).

The pairs with the most shared genes are broad traits that share thousands of gene labels and few leads. Body height and hematological measurement share 2,086 genes, 2 of them on point loci, 4 lead sites, and 1 Tier A allele.

| trait A | trait B | genes | genes, span < 1 Mb | genes, point | lead sites | Tier A alleles |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| body height | hematological measurement | 2,086 | 247 | 2 | 4 | 1 |
| body height | gut microbiome measurement | 1,869 | 251 | 5 | 1 | 0 |
| body height | body mass index | 1,806 | 319 | 15 | 17 | 2 |
| bone tissue density | body height | 1,677 | 289 | 6 | 3 | 0 |
| body height | protein measurement | 1,541 | 274 | 7 | 6 | 3 |
| body height | body shape measurement | 1,509 | 252 | 2 | 64 | 10 |
| body height | systolic blood pressure | 1,419 | 220 | 7 | 7 | 1 |
| diastolic blood pressure | systolic blood pressure | 1,387 | 750 | 53 | 355 | 25 |
| body height | educational attainment | 1,364 | 266 | 7 | 4 | 0 |
| body height | health trait | 1,353 | 135 | 1 | 46 | 7 |

The pairs with the most shared Tier A alleles are closer measurements. They still share more gene labels than alleles. The allele count can exceed the rsID count when several ClinVar alleles inherit one lead's trait list. 47,843 variant-sharing pairs have more keys than rsIDs.

| trait A | trait B | Tier A alleles | rsIDs | lead sites | genes | genes, point |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| glomerular filtration rate | serum creatinine amount | 30 | 24 | 145 | 798 | 36 |
| diastolic blood pressure | systolic blood pressure | 25 | 25 | 355 | 1,387 | 53 |
| reticulocyte amount | reticulocyte count | 23 | 22 | 253 | 540 | 50 |
| pulse pressure measurement | systolic blood pressure | 22 | 22 | 245 | 939 | 27 |
| cholesterol in very large HDL measurement | total lipids in very large HDL measurement | 22 | 20 | 162 | 169 | 106 |
| phospholipids in very large HDL measurement | total lipids in very large HDL measurement | 22 | 20 | 165 | 177 | 104 |
| HDL cholesterol change measurement | cholesteryl esters in HDL measurement | 22 | 19 | 174 | 158 | 96 |
| hematocrit | hemoglobin measurement | 22 | 19 | 215 | 672 | 46 |
| cystatin C measurement | serum creatinine amount | 22 | 17 | 60 | 345 | 18 |
| linoleic acid measurement | omega-6 polyunsaturated fatty acid measurement | 21 | 20 | 116 | 375 | 32 |

The pair with the most shared lead sites is Alzheimer disease and family history of Alzheimer’s disease: 584 sites, 7 Tier A alleles, 887 genes, 94 of them on point loci. Repeated studies and parent-like EFO terms produce that pattern. This report does not collapse the ontology.

## Caveat

A shared gene is a repeated catalog label. It is not a shared mechanism, a colocalization, or a pathway. On a point locus the matched symbol usually contains the lead, and that is still only the gene the catalog named. On a wide locus the token can come from another member up to the clump radius, and the lead often does not sit in the gene. The height pairs show the difference: thousands of shared genes, a handful of shared leads.

A shared Tier A allele is the stronger count, and it is a minority of leads. Most GWAS loci have no ref/alt key (290,815 of 312,875 ledger GWAS rows). The alleles that do resolve are the ones ClinVar already knew. `rs1354034` is the clearest miss: 1,238 traits, inside `ARHGEF3`, and no key. Do not read the variant table as the catalog's pleiotropic leads. Read it as the pleiotropic leads that gained an alternate allele.

EFO traits repeat one another. Blood-pressure pair, HDL fraction pairs, creatinine and filtration, and disease versus family history are not independent phenotypes. There is no linkage-disequilibrium panel here, so two nearby leads are not tested for tagging the same variant. The distance clump is not a causal set. Nothing in this table is a polygenic score.
