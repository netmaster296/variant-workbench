# Coding versus noncoding

## Question

What share of clumped GWAS loci have a lead association whose catalog `CONTEXT` is missense, intronic, UTR, or intergenic, overall and among the largest trait groups? `INTERGENIC` is reported beside that term. It is not used as a substitute for it.

## Filters

Loci are the rows in `warehouse/gwas/loci.parquet`. Each row is one distance clump within one mapped trait URI: the lead is the member with the smallest p-value, and every member is within 1,000,000 bp of that lead. Associations with no position are not in this file. The class is not recomputed from the clump window.

The lead association is the row in `warehouse/gwas/associations.parquet` with the same trait URI, chromosome, position, rsID, and parsed p-value. A mapped trait URI is one comma-separated token of `MAPPED_TRAIT_URI`, trimmed. An empty URI matches an empty `MAPPED_TRAIT_URI`. The join recovered 311,160 of 311,160 loci. For 3,843 loci the key matched more than one association row (maximum 12). Those rows agreed on `CONTEXT` and on `INTERGENIC`, and the match count never exceeded the locus member count, so the class does not depend on which duplicate is kept. Trait names joined for every non-null URI (0 misses).

`CONTEXT` is split on `;` and on ` x `. A locus counts as:

- missense if a token is `missense_variant`
- intronic if a token is `intron_variant`
- UTR if a token is `3_prime_UTR_variant` or `5_prime_UTR_variant`
- intergenic, by context, if a token is `intergenic_variant`

No lead had more than one token, so these four classes do not overlap. `stop_gained`, frameshift, inframe insertion or deletion, `start_lost`, `stop_lost`, `stop_retained_variant`, and `protein_altering_variant` are not missense. Splice terms are not intronic. `regulatory_region_variant`, `TF_binding_site_variant`, and `TFBS_ablation` are not intergenic.

The `INTERGENIC` flag is a separate column. Flag `1` is intergenic by that flag. Flag `0` is the catalog's other value. Two leads have a blank flag. The flag is not merged into the `CONTEXT` class.

Mapped genes, gene spans, and other members of the clump are not used. The checked lead for trait `http://www.ebi.ac.uk/efo/EFO_0004348` on chr22 is `rs9330814` at 45968311, with `CONTEXT` `intron_variant`, `INTERGENIC` `0`, and 19 members. That is the lead, not another SNP in the same clump.

A trait group is one mapped trait URI. There is no parent-ontology table on disk, so traits are not rolled up further. URIs in the large-trait table include EFO, OBA, and MONDO. The largest groups are the 18 traits with at least 1,000 clumped loci. Eight loci have an empty trait URI. They stay in the overall denominator and in none of those 18.

## Denominator

- 311,160 clumped loci. This is every row in `warehouse/gwas/loci.parquet`.
- 311,150 loci whose lead `CONTEXT` is non-empty. Ten leads have an empty `CONTEXT` string. Eight of those ten have a `SNPS` id beginning with `esv` or `nsv`. All ten have `INTERGENIC` `0`.

Shares of the 311,150 are blank for the empty-context row, because that row is not in the second denominator. Counts are exact. A share is the count divided by the stated denominator, rounded half up to 6 decimal places. Rounded shares of 311,160 sum to 0.999999, short of 1 by 0.000001. The six exclusive counts sum to 311,160.

These 311,160 rows are trait-specific clumps, not distinct variants. They cover 187,396 distinct lead positions. Forty-one leads have no rsID. The same position can be counted once for each trait in which it is the lead.

## Result

Every lead with a `CONTEXT` has exactly one term. The exclusive class and the "contains this term" count are the same number.

| class | lead `CONTEXT` | loci | share of 311,160 | share of 311,150 |
| --- | --- | ---: | ---: | ---: |
| intronic | `intron_variant` | 194523 | 0.625154 | 0.625174 |
| intergenic | `intergenic_variant` | 68777 | 0.221034 | 0.221041 |
| other single term | any other non-empty term | 19290 | 0.061994 | 0.061996 |
| missense | `missense_variant` | 17411 | 0.055955 | 0.055957 |
| UTR | `3_prime_UTR_variant` or `5_prime_UTR_variant` | 11149 | 0.035830 | 0.035832 |
| no `CONTEXT` | empty string | 10 | 0.000032 |  |

UTR is 8,205 `3_prime_UTR_variant` (0.026369 of loci; 0.026370 of loci with a `CONTEXT`) and 2,944 `5_prime_UTR_variant` (0.009461; 0.009462). No lead is both.

Intronic, UTR, and intergenic `CONTEXT` together are 274,449 loci (0.882019 of 311,160; 0.882047 of 311,150). The four requested classes together are 291,860 (0.937974; 0.938004).

The other 19,290 leads are also one term each:

| term | loci | share of 311,160 | share of 311,150 |
| --- | ---: | ---: | ---: |
| `regulatory_region_variant` | 6925 | 0.022255 | 0.022256 |
| `non_coding_transcript_exon_variant` | 6345 | 0.020391 | 0.020392 |
| `synonymous_variant` | 2474 | 0.007951 | 0.007951 |
| `stop_gained` | 1078 | 0.003464 | 0.003465 |
| `TF_binding_site_variant` | 598 | 0.001922 | 0.001922 |
| `splice_region_variant` | 590 | 0.001896 | 0.001896 |
| `splice_polypyrimidine_tract_variant` | 580 | 0.001864 | 0.001864 |
| `splice_donor_region_variant` | 142 | 0.000456 | 0.000456 |
| `frameshift_variant` | 114 | 0.000366 | 0.000366 |
| `splice_donor_variant` | 100 | 0.000321 | 0.000321 |
| `splice_donor_5th_base_variant` | 89 | 0.000286 | 0.000286 |
| `mature_miRNA_variant` | 73 | 0.000235 | 0.000235 |
| `inframe_insertion` | 58 | 0.000186 | 0.000186 |
| `splice_acceptor_variant` | 48 | 0.000154 | 0.000154 |
| `start_lost` | 35 | 0.000112 | 0.000112 |
| `stop_lost` | 19 | 0.000061 | 0.000061 |
| `inframe_deletion` | 16 | 0.000051 | 0.000051 |
| `TFBS_ablation` | 3 | 0.000010 | 0.000010 |
| `stop_retained_variant` | 2 | 0.000006 | 0.000006 |
| `protein_altering_variant` | 1 | 0.000003 | 0.000003 |

Of those, 1,323 are a protein-sequence term other than missense (`stop_gained`, `frameshift_variant`, `inframe_insertion`, `inframe_deletion`, `start_lost`, `stop_lost`, `stop_retained_variant`, `protein_altering_variant`). Synonymous is 2,474. Splice terms are 1,549. Regulatory or transcription-factor terms are 7,526. Noncoding-transcript exon is 6,345. Mature miRNA is 73.

`INTERGENIC` on the same leads:

| lead `INTERGENIC` | loci | share of 311,160 | share of 311,150 |
| --- | ---: | ---: | ---: |
| `0` | 210631 | 0.676922 | 0.676944 |
| `1` | 100527 | 0.323072 | 0.323082 |
| blank | 2 | 0.000006 | 0.000006 |

Flag `1` is 100,527 loci. The intergenic `CONTEXT` term is 68,777. Those sets are not nested.

| `CONTEXT` class | flag `0` | flag `1` | flag blank | loci |
| --- | ---: | ---: | ---: | ---: |
| intronic | 161744 | 32778 | 1 | 194523 |
| intergenic | 9566 | 59210 | 1 | 68777 |
| other single term | 10779 | 8511 | 0 | 19290 |
| missense | 17406 | 5 | 0 | 17411 |
| UTR | 11126 | 23 | 0 | 11149 |
| no `CONTEXT` | 10 | 0 | 0 | 10 |

59,210 leads carry both the intergenic term and flag `1`. 9,566 have the intergenic term and flag `0`, and 1 has the intergenic term and a blank flag. 41,317 have flag `1` without the intergenic term: 32,778 intronic, 8,511 other, 23 UTR, and 5 missense.

There are 15,109 trait groups with at least one clump, including the empty-URI group. Eighteen traits have at least 1,000 loci. They hold 22,050 loci (0.070864 of 311,160). Every one of those leads has a `CONTEXT`, so the two denominators are the same inside this table. Shares below are of that trait's loci.

| trait | loci | missense | intronic | UTR | 3' | 5' | intergenic `CONTEXT` | other | `INTERGENIC` = 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| body height | 1752 | 81 (0.046233) | 1122 (0.640411) | 65 (0.037100) | 58 | 7 | 390 (0.222603) | 94 (0.053653) | 589 (0.336187) |
| gut microbiome measurement | 1691 | 9 (0.005322) | 1050 (0.620934) | 23 (0.013601) | 18 | 5 | 546 (0.322886) | 63 (0.037256) | 757 (0.447664) |
| bone tissue density | 1510 | 28 (0.018543) | 948 (0.627815) | 30 (0.019868) | 24 | 6 | 453 (0.300000) | 51 (0.033775) | 644 (0.426490) |
| protein measurement | 1395 | 149 (0.106810) | 781 (0.559857) | 57 (0.040860) | 36 | 21 | 313 (0.224373) | 95 (0.068100) | 463 (0.331900) |
| educational attainment | 1368 | 20 (0.014620) | 896 (0.654971) | 45 (0.032895) | 33 | 12 | 359 (0.262427) | 48 (0.035088) | 501 (0.366228) |
| body mass index | 1334 | 80 (0.059970) | 817 (0.612444) | 44 (0.032984) | 37 | 7 | 338 (0.253373) | 55 (0.041229) | 458 (0.343328) |
| Alzheimer disease | 1208 | 21 (0.017384) | 772 (0.639073) | 15 (0.012417) | 12 | 3 | 368 (0.304636) | 32 (0.026490) | 519 (0.429636) |
| hematological measurement | 1204 | 27 (0.022425) | 830 (0.689369) | 53 (0.044020) | 36 | 17 | 236 (0.196013) | 58 (0.048173) | 337 (0.279900) |
| breastfeeding duration | 1125 | 4 (0.003556) | 724 (0.643556) | 11 (0.009778) | 9 | 2 | 353 (0.313778) | 33 (0.029333) | 496 (0.440889) |
| COVID-19 | 1125 | 13 (0.011556) | 834 (0.741333) | 20 (0.017778) | 18 | 2 | 228 (0.202667) | 30 (0.026667) | 281 (0.249778) |
| systolic blood pressure | 1100 | 43 (0.039091) | 706 (0.641818) | 43 (0.039091) | 34 | 9 | 250 (0.227273) | 58 (0.052727) | 367 (0.333636) |
| smoking initiation | 1098 | 26 (0.023679) | 726 (0.661202) | 20 (0.018215) | 16 | 4 | 279 (0.254098) | 47 (0.042805) | 405 (0.368852) |
| free cholesterol measurement | 1044 | 24 (0.022989) | 641 (0.613985) | 16 (0.015326) | 13 | 3 | 310 (0.296935) | 53 (0.050766) | 431 (0.412835) |
| diastolic blood pressure | 1028 | 65 (0.063230) | 643 (0.625486) | 36 (0.035019) | 27 | 9 | 229 (0.222763) | 55 (0.053502) | 340 (0.330739) |
| body shape measurement | 1027 | 36 (0.035054) | 666 (0.648491) | 40 (0.038948) | 31 | 9 | 227 (0.221032) | 58 (0.056475) | 329 (0.320351) |
| allergen exposure measurement | 1024 | 5 (0.004883) | 632 (0.617188) | 18 (0.017578) | 14 | 4 | 332 (0.324219) | 37 (0.036133) | 461 (0.450195) |
| glioma pathogenesis-related protein 1 measurement | 1013 | 13 (0.012833) | 668 (0.659427) | 21 (0.020731) | 19 | 2 | 280 (0.276407) | 31 (0.030602) | 402 (0.396841) |
| body weight | 1004 | 35 (0.034861) | 650 (0.647410) | 42 (0.041833) | 32 | 10 | 235 (0.234064) | 42 (0.041833) | 324 (0.322709) |

Pooled, those 18 traits are 679 missense (0.030794), 14,106 intronic (0.639728), 599 UTR (0.027166), 5,726 intergenic by `CONTEXT` (0.259683), 940 other (0.042630), and 8,104 with `INTERGENIC` `1` (0.367528). Intronic is the largest class in every one of the 18. Missense ranges from 0.003556 (breastfeeding duration) to 0.106810 (protein measurement). Intergenic `CONTEXT` ranges from 0.196013 (hematological measurement) to 0.324219 (allergen exposure measurement). Flag `1` is higher than the intergenic `CONTEXT` count in every one of the 18.

Mapped trait URIs for that table, in the same order: `http://purl.obolibrary.org/obo/OBA_VT0001253`, `http://www.ebi.ac.uk/efo/EFO_0007874`, `http://purl.obolibrary.org/obo/OBA_1000110`, `http://www.ebi.ac.uk/efo/EFO_0004747`, `http://www.ebi.ac.uk/efo/EFO_0011015`, `http://www.ebi.ac.uk/efo/EFO_0004340`, `http://purl.obolibrary.org/obo/MONDO_0004975`, `http://www.ebi.ac.uk/efo/EFO_0004503`, `http://www.ebi.ac.uk/efo/EFO_0006864`, `http://purl.obolibrary.org/obo/MONDO_0100096`, `http://www.ebi.ac.uk/efo/EFO_0006335`, `http://www.ebi.ac.uk/efo/EFO_0005670`, `http://www.ebi.ac.uk/efo/EFO_0008591`, `http://www.ebi.ac.uk/efo/EFO_0006336`, `http://www.ebi.ac.uk/efo/EFO_0801227`, `http://www.ebi.ac.uk/efo/EFO_0007944`, `http://www.ebi.ac.uk/efo/EFO_0802567`, `http://www.ebi.ac.uk/efo/EFO_0004338`.

Numbers are from `warehouse/reports/sql/coding_versus_noncoding.sql`, run from the repository root.

## Caveat

This is a label on the reported lead association. It is not a statement that the locus acts through that consequence, and it is not a fine-mapped or causal variant. A missense lead can sit on a haplotype whose relevant base is not the missense base. An intergenic or intronic lead can sit near a coding variant that this count never sees. Clumping is by distance within a trait, not by linkage disequilibrium.

`CONTEXT` and `INTERGENIC` disagree often enough to change the intergenic share. Using the term, intergenic is 0.221034 of loci. Using the flag, it is 0.323072. Most of the extra flag-`1` leads are intronic by `CONTEXT`. This report does not treat either field as a correction of the other.

Missense is 0.055955 of loci, not the whole coding slice. Another 0.004252 are a different protein-sequence term, and 0.007951 are synonymous. Intronic, UTR, and intergenic leads, 0.882019 of loci by `CONTEXT`, are the part of this catalog that a coding clinical record and a gene-constraint score do not classify. This track does not join ClinVar or constraint, so it does not measure how often those leads also have a submission or a constrained gene. A gene constraint score is not a functional class for the lead SNP.

The 18 largest traits are 0.070864 of the loci. Their pooled missense share (0.030794) is lower than the catalog-wide missense share, and their pooled intergenic-`CONTEXT` share (0.259683) is higher. That is a description of this cutoff, not a property of every large trait: protein measurement is 0.106810 missense, and hematological measurement is 0.196013 intergenic by `CONTEXT`.
