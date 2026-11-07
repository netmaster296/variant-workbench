# Splicing

## Question

Which ClinVar splice-region alleles have a SpliceAI or Pangolin score, and which of those genes are GTEx sGenes?

This is an overlap of three evidence classes. A ClinVar label, a predictor score, and an sQTL are not one claim, and none of them shows that the allele is the sQTL.

Numbers are the result sets of `warehouse/reports/sql/splicing.sql`. The variant-level columns are the view `variant_overlap`. The gene-level columns are the view `gene_overlap`. Neither view is written to a mart.

## Filters

Source tables:

- `warehouse/clinvar/allele_dedup.parquet` for the ClinVar allele, gene, and clinical significance.
- `warehouse/variants/tier_a.parquet` for gnomAD v4.1 scores and VEP consequence. Genome columns are `spliceai_ds_max` and `pangolin_largest_ds`. Exome columns are `spliceai_ds_max_exome` and `pangolin_largest_ds_exome`. Consequence strings are `vep_consequence` and `vep_consequence_exome`, terms joined by `&`.
- `warehouse/gtex/sgenes.parquet` for sGene status. `warehouse/gtex/sqtl_pairs.parquet` was not read.
- `warehouse/crosswalk/genes.parquet` only to attach an ENSG. The key is `GeneID` = `ncbi_gene_id_primary`.

`allele_dedup` has no splice-consequence column. The splice-region call is the gnomAD VEP consequence on the same `variant_key`.

An allele is splice-region if either consequence string contains `splice_region_variant`, `splice_donor_region_variant`, or `splice_donor_5th_base_variant`. VEP does not co-emit the parent term on the two child terms (overlap counts are 0 on both datasets), so a parent-only filter would drop the donor-side region alleles. The three terms do not co-occur with each other here.

Not in the set:

- `splice_donor_variant` or `splice_acceptor_variant` with no region term: 37622 alleles.
- `splice_polypyrimidine_tract_variant` with no region term: 92683 alleles. An allele that also carries a region term stays in.
- 786 splice-region alleles also carry a canonical donor or acceptor term on the stored consequence string. They stay in. The set is not "region and nothing else."

Not used as filters: clinical significance, review stars, origin, variant type, or a score cutoff. Score cutoffs below describe SpliceAI. They do not define the set. Pangolin is not put on the SpliceAI scale.

A gene is an sGene when any `sgenes.parquet` row for that bare ENSG has `qval_le_0_05`. That flag is qval <= 0.05 at the intron-cluster grain. It is not collapsed before this count. Absence from the table means the ENSG is not in the GTEx v8 sGene file. It is not a negative sQTL result.

The variant-key join is intact: 0 keyed ClinVar alleles are missing from `tier_a`, and 0 Tier A rows tagged `clinvar` are missing from `allele_dedup`.

## Denominator

| set | n |
| --- | ---: |
| `allele_dedup` rows | 4639253 |
| null `variant_key`, excluded | 20608 |
| keyed GRCh38 alleles, the denominator | 4618645 |
| in the genome sites VCF | 2093198 |
| in the exome sites VCF | 3134607 |
| in neither sites VCF | 1316666 |
| in at least one sites VCF | 3301979 |
| splice-region alleles | 195091 |

The 20608 null keys have no VCF allele, so they cannot take a score or a VEP call. The 1316666 alleles in neither sites VCF are inside the denominator and outside the splice-region count: there is no consequence string to call them. Splice-region alleles are 0.042240 of keyed alleles and 0.059083 of alleles that have a VCF.

Origin is not a filter. Of the 195091: germline 174809, other 20131, germline/somatic 98, somatic 53. Type is not a filter: single nucleotide variant 181376, Deletion 7606, Duplication 3306, Microsatellite 2241, Insertion 562.

## Result

Where both VCFs contain the allele, the splice-region call agrees. Zero alleles are present in both VCFs and called splice-region on only one side. The genome-only and exome-only rows are alleles missing from the other VCF.

| call | n |
| --- | ---: |
| region in both VCFs | 104771 |
| genome region, not in the exome VCF | 6202 |
| exome region, not in the genome VCF | 84118 |
| present in both, region in only one | 0 |

| term | genome | exome |
| --- | ---: | ---: |
| `splice_region_variant` | 93506 | 158656 |
| `splice_donor_region_variant` | 11979 | 20313 |
| `splice_donor_5th_base_variant` | 5488 | 9920 |
| parent co-occurring with either child | 0 | 0 |

Genome labels sum to 110973, the alleles in the genome VCF. Exome labels sum to 188889, the alleles in the exome VCF.

A score means a non-null cell. It does not mean a high score. 195085 of 195091 splice-region alleles have at least one of the four scores. Six have none. No score is stored for a dataset the allele is not in.

| score | non-null | missing among alleles in that VCF | non-null with no call in that VCF |
| --- | ---: | ---: | ---: |
| `spliceai_ds_max` | 109939 | 1034 of 110973 | 0 |
| `spliceai_ds_max_exome` | 187292 | 1597 of 188889 | 0 |
| `pangolin_largest_ds` | 110877 | 96 of 110973 | 0 |
| `pangolin_largest_ds_exome` | 188887 | 2 of 188889 | 0 |

When both values exist they are equal on this set. An allele that is in both VCFs does not have either score on only one side. The columns are still not combined.

| comparison | n |
| --- | ---: |
| in both VCFs | 104771 |
| both SpliceAI scores non-null | 103815 |
| those SpliceAI scores equal | 103815 |
| SpliceAI max absolute difference | 0 |
| both Pangolin scores non-null | 104769 |
| those Pangolin scores equal | 104769 |
| Pangolin max absolute difference | 0 |

Nulls in the bin table include alleles that are not in that VCF. Quantiles are among non-null values. The 0.2 and 0.5 cuts are the usual SpliceAI delta-score descriptions, not a filter. The same numeric bins are shown for Pangolin only so the shape is visible. A negative Pangolin value is not a low SpliceAI score.

| score | null | <0 | 0 | (0, 0.2) | [0.2, 0.5) | [0.5, 0.8) | >=0.8 | min | p50 | p90 | max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| spliceai_genome | 85152 | 0 | 54148 | 39537 | 8437 | 4273 | 3544 | 0 | 0.01 | 0.34 | 1 |
| spliceai_exome | 7799 | 0 | 88163 | 66533 | 15776 | 8878 | 7942 | 0 | 0.01 | 0.44 | 1 |
| pangolin_genome | 84214 | 51794 | 27845 | 27654 | 1954 | 1144 | 486 | -0.96 | 0 | 0.04 | 0.92 |
| pangolin_exome | 6204 | 92153 | 45224 | 44901 | 3509 | 2150 | 950 | -0.97 | 0 | 0.04 | 0.92 |

Among non-null SpliceAI scores, 16254 of 109939 genome scores are >= 0.2 (0.147846) and 7817 are >= 0.5 (0.071103). Exome: 32596 of 187292 are >= 0.2 (0.174038) and 16820 are >= 0.5 (0.089806). The exome counts are larger because more alleles are in the exome VCF, not because the two scores disagree. Of non-null Pangolin scores, 51794 of 110877 genome values are negative (0.467130) and 92153 of 188887 exome values are negative (0.487874). Median Pangolin is 0.

Clinical significance stays in its own column. `-` is the source blank, not a missing join. `other` is the 150 alleles outside the nine strings below (drug response, risk factor, not provided, and the long tail). `>=0.2` and `>=0.5` are counts of non-null scores at those cuts. They are not subtracted from the non-null counts.

| clinical significance | n | spliceai genome | genome >=0.2 | genome >=0.5 | spliceai exome | exome >=0.2 | exome >=0.5 | pangolin genome | pangolin exome |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Uncertain significance | 75079 | 39056 | 8057 | 3554 | 71831 | 16743 | 7895 | 39283 | 72351 |
| Likely benign | 67697 | 36516 | 1290 | 332 | 64911 | 1970 | 444 | 36788 | 65364 |
| `-` | 19109 | 9660 | 3577 | 2012 | 18150 | 7817 | 4585 | 9884 | 18498 |
| Conflicting classifications of pathogenicity | 13940 | 10290 | 1451 | 751 | 13617 | 2217 | 1228 | 10374 | 13718 |
| Benign | 7818 | 7008 | 416 | 141 | 7696 | 447 | 149 | 7087 | 7781 |
| Benign/Likely benign | 5594 | 5337 | 280 | 77 | 5529 | 286 | 77 | 5384 | 5578 |
| Pathogenic | 2207 | 677 | 362 | 268 | 2092 | 1046 | 768 | 679 | 2102 |
| Likely pathogenic | 2063 | 593 | 386 | 334 | 1940 | 1227 | 1002 | 594 | 1962 |
| Pathogenic/Likely pathogenic | 1434 | 727 | 414 | 334 | 1389 | 788 | 636 | 728 | 1393 |
| other | 150 | 75 | 21 | 14 | 137 | 55 | 36 | 76 | 140 |

Among alleles with a non-null score, SpliceAI >= 0.5 is common in pathogenic classes and rare in benign classes. Pathogenic: 268 of 677 genome scores (0.395864) and 768 of 2092 exome scores (0.367113). Likely pathogenic: 334 of 593 (0.563238) and 1002 of 1940 (0.516495). Pathogenic/Likely pathogenic: 334 of 727 (0.459422) and 636 of 1389 (0.457883). Likely benign: 332 of 36516 (0.009092) and 444 of 64911 (0.006840). Uncertain significance: 3554 of 39056 (0.090998) and 7895 of 71831 (0.109911).

sQTL evidence is a different column. The gene link does not use the VEP symbol. 194351 alleles match one ENSG, 238 match two, and 502 match none. Of the unlinked alleles, 435 have `GeneID` `-1` (a multi-gene or non-gene ClinVar row). The other 67 carry an NCBI id that is not `ncbi_gene_id_primary` on any crosswalk gene. For the 238 two-ENSG alleles, 7 have both ENSGs as sGenes, 148 have neither, and 83 disagree. Those 83 are not called sGenes in the one-ENSG column.

| clinical significance | n | one ENSG and sGene | one ENSG, tested, not sGene | one ENSG, not in the sGene table | unlinked | two ENSGs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Uncertain significance | 75079 | 61262 | 10501 | 3091 | 171 | 54 |
| Likely benign | 67697 | 55248 | 9120 | 3072 | 208 | 49 |
| `-` | 19109 | 15352 | 2041 | 1563 | 41 | 112 |
| Conflicting classifications of pathogenicity | 13940 | 11398 | 1932 | 584 | 20 | 6 |
| Benign | 7818 | 6248 | 1078 | 471 | 14 | 7 |
| Benign/Likely benign | 5594 | 4397 | 777 | 392 | 20 | 8 |
| Pathogenic | 2207 | 1751 | 337 | 113 | 5 | 1 |
| Likely pathogenic | 2063 | 1653 | 281 | 108 | 21 | 0 |
| Pathogenic/Likely pathogenic | 1434 | 1166 | 238 | 29 | 0 | 1 |
| other | 150 | 106 | 31 | 11 | 2 | 0 |

The gene question uses genes of scored alleles, not genes of the six unscored alleles. One ClinVar `GeneID` can attach to two ENSGs. Both are counted. There are 14253 ENSGs and 14225 symbols.

| biotype | genes | sGenes | tested, not sGene | not in the sGene table |
| --- | ---: | ---: | ---: | ---: |
| protein_coding | 14140 | 10655 | 2447 | 1038 |
| lncRNA | 105 | 60 | 8 | 37 |
| transcribed_unprocessed_pseudogene | 5 | 4 | 0 | 1 |
| transcribed_unitary_pseudogene | 2 | 1 | 0 | 1 |
| IG_C_gene | 1 | 0 | 0 | 1 |
| total | 14253 | 10720 | 2455 | 1078 |

Of those 10720 sGenes, the median number of tissues with qval <= 0.05 is 4 (minimum 1, maximum 49).

A separate, narrower column asks whether the ClinVar allele itself is the lead variant of a significant sGene phenotype. 229 splice-region alleles are that lead for a linked ENSG. Another 37 are a significant lead only for a different ENSG. 266 are a significant lead for some gene in the sGene table. The pair file was not scanned, so this is not the count of alleles that are significant sQTL pairs.

Six splice-region alleles have no SpliceAI or Pangolin score in either column:

| variant_key | gene | type | in genomes | in exomes | clinical significance | linked ENSGs |
| --- | --- | --- | --- | --- | --- | ---: |
| chr13:94504908:T:C | DCT | single nucleotide variant | true | false | Likely benign | 1 |
| chr14:101760722:AGCTGGTGAGGGGAGGGGCTGGTCGAGGGGAGGG:A | PPP2R5C | Deletion | true | true | Likely benign | 1 |
| chr15:84961522:G:C | SLC28A1 | single nucleotide variant | true | true | Benign | 1 |
| chr21:45990416:CGTGAAGGTGACTGGGGGGAGATAGGATGGACGGGGAGGGACGAGGAGGAATGGGGCGAGATGGGGAGGGACGGAGTGGACGGCGTGAAGGTGACCCGGGGAGGGATGGGGTGGACAGTGTGAAGGTGACCAGGGGAAGGACGGGGAGGGACGGGGAGGGATGGGGTGAGGTGATCCCGGCAGGAGGGACAGGGAGGAGTGGGGTGGACGGT:C | COL6A1 | Deletion | true | false | Uncertain significance | 1 |
| chr2:226179800:G:C | LOC646736 | single nucleotide variant | true | false | `-` | 0 |
| chr2:94755952:C:G | FAM95A | single nucleotide variant | true | false | `-` | 1 |

## Caveat

Most tested genes are already sGenes. In `sgenes.parquet`, 16895 of 23879 genes have qval <= 0.05 in at least one tissue (0.707525). Among protein-coding crosswalk genes that appear in that table, 12765 of 16189 do (0.788498). Among splice-region genes that appear in the table, 10720 of 13175 do (0.813662), and 10655 of 13102 protein-coding genes do (0.813235). The overlap is close to the base rate. It is not evidence that the ClinVar allele explains the sQTL.

For one-ENSG alleles, 27520 of 33688 with a SpliceAI score >= 0.2 on either dataset map to an sGene (0.816908). Of one-ENSG alleles whose SpliceAI scores are all below 0.2, 129739 of 159022 do (0.815856). A higher SpliceAI score does not move the sGene flag. On the one-ENSG rows, pathogenic alleles are an sGene in 1751 of 2201 cases and likely-benign alleles in 55248 of 67440. Clinical class moves the SpliceAI tail, not the sGene column.

229 of 195085 scored splice-region alleles are themselves the lead variant for a linked gene (0.001174). That lead is the top variant for an intron cluster, not a colocalization and not a test that this allele is causal. Significant pairs that are not the lead were not counted.

1078 linked genes are missing from the sGene table. GTEx v8 and the GENCODE crosswalk do not use the same gene catalog. Missing is "not in this file under this ENSG," not "no splicing QTL."

Genome and exome scores were equal wherever both were non-null on this set. That is not a reason to add them or to treat one dataset's coverage as the other's. Alleles absent from a VCF have a null score there even when the other dataset scored them. Pangolin `largest_ds` runs below 0, with a median of 0, so it cannot be read as a SpliceAI probability. A non-null SpliceAI score is the usual state of this set (195085 of 195091). The median non-null SpliceAI score is 0.01.
