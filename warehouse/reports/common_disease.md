# Common-disease gene hypotheses

## Question

For an EFO trait with at least 50 clumped GWAS loci, which genes mapped on those loci are also eGenes in a tissue assigned to the trait?

A match is a shared-gene hypothesis. It is not colocalization.

## Filters

Minimum locus count: 50. A trait is one `trait_uri`. Loci are the rows in `warehouse/gwas/loci.parquet`. The count is the number of clumps, not the number of catalog association rows.

Trait names come from `warehouse/gwas/traits.parquet`, joined on `trait_uri` with nulls matching nulls. That join missed no locus among traits with at least 50 clumps. The gene mart and `egenes.parquet` do not carry the trait name, so the name table is required for the tissue rule.

Tissue assignment is a name rule, and only these four classes are used. A trait that matched two classes would have been left unmapped. None did. Traits that match none are listed in the unmapped table and have no eGene claim.

Blood, matched to QTL tissue `Whole_Blood`: the name matches a hematologic term (hematology, hematocrit, hemoglobin, erythrocyte, red blood cell, red cell, leukocyte, white blood cell, white cell, neutrophil, lymphocyte, monocyte, eosinophil, basophil, platelet, thrombocyte, reticulocyte, corpuscular, granulocyte, myeloid leukocyte, blood cell, or hematopoiesis, including British haemat- spellings) and does not match leukemia, lymphoma, neoplasm, cancer, carcinoma, tumor, esophagitis, platelet-derived, monounsaturated, monoclonal, hemoglobin A1, glycated, glycosylated, A1c, or blood pressure.

Brain, matched to every GTEx v8 eQTL tissue whose name starts with `Brain_`, 13 tissues, including spinal cord and not including pituitary: the name matches a listed neurologic, psychiatric, cognitive, imaging, stroke, or sleep term (Alzheimer, Parkinson, Lewy body, dementia, schizophrenia, bipolar, depressive or depression, attention deficit or ADHD, autism, epilepsy, migraine, neuroticism or neurotic disorder, psychiatric disorder, psychotic, insomnia, cognitive, intelligence, memory, mathematical ability, neuroimaging, neuropsychological, neurofibrillary, the word brain, amygdala volume, cerebellar volume, cerebral cortex, cortical thickness, hippocampal volume, ischemic stroke, stroke disorder, large artery stroke, small vessel stroke, stroke outcome, sleep duration, sleep characteristic, white matter, grey or gray matter, multiple sclerosis, amyotrophic lateral sclerosis, obsessive-compulsive, Tourette, anxiety, narcolepsy, cerebral amyloid) and does not match natriuretic, brain-derived, neuroblastoma, neurofilament, neuropathy, neuropathic, sleep apnea, antidepressant, antipsychotic, or educational attainment.

Heart, matched to `Heart_Atrial_Appendage` and `Heart_Left_Ventricle`: the name matches QRS, QT interval, PR interval, JT interval, atrial fibrillation, atrial flutter, angina pectoris, cardiac troponin, coronary artery, coronary atherosclerosis, heart failure, heart rate, heart malformation, cardiomyopathy, left ventricular, myocardial infarction, myocardial ischemia, imaging of the heart, sudden cardiac arrest, ventricular septal defect, cardiac arrest, aortic stenosis, or aortic valve, and does not match blood pressure, pulse pressure, or stroke.

Liver enzyme, matched to `Liver`: the name matches alanine aminotransferase, alanine transaminase, aspartate aminotransferase, aspartate transaminase, gamma-glutamyl transferase, gamma-glutamyltransferase, or alkaline phosphatase, and does not match bone alkaline, intestinal alkaline, or phenylalanine. Alkaline phosphatase is in the usual liver-enzyme panel. It is not specific to liver.

The exact patterns are in `warehouse/reports/sql/common_disease.sql`. Every assigned name is in the trait table below.

eGene means `qval_le_0_05` is true in `warehouse/gtex/egenes.parquet`. That flag disagrees with `qval <= 0.05` on 0 rows. Rows with a larger q-value are not eGenes. A gene with no eGene row in a tissue was not tested there. Untested is not a negative result.

A mapped gene is a token in `mapped_genes.tokens`. The token must equal one `symbol` in `warehouse/marts/gene_mart.parquet`, that symbol must be on exactly one mart row, and `symbol_n_gencode` must be 1. The eGene join is mart `ensg` to egenes `ensg`. It does not use the gene symbol, and it does not use `egenes.gene_name`.

Assigned traits contribute 94704 mapped-gene tokens. 30507 contain ` - `. None of those distance labels matched a symbol, and this report does not split them. Ambiguous symbols stay out of the hypothesis count. Non-distance tokens that matched no symbol: 0.

Locus `start` and `end` are the minimum and maximum member positions. They are not a 1 Mb window. This report does not overlap them with the eQTL variant, and it does not expand them.

Expression tissue names and QTL tissue names are different strings. The map below is the only translation used. `median_tpm.parquet` was read only to check that each expression label exists. It is not an eGene test. All 17 QTL labels are present in `egenes.parquet`. All 17 expression labels are present in `median_tpm.parquet`. `Liver` is the only pair whose two strings are identical. `Whole Blood` is not `Whole_Blood`. The same kind of split is why `Adipose - Subcutaneous` is not `Adipose_Subcutaneous`. Adipose is not a tissue in this map.

## Denominator

Traits with at least one clumped locus: 15109. Clumped loci: 311160. The median trait has 3 loci. The 90th percentile is 31. Traits with fewer than 50 loci: 14025 / 15109 (0.928).

Traits with at least 50 loci: 1084, covering 221818 loci. Assigned to one tissue class: 151 traits and 47827 loci. Unmapped: 933 traits and 173991 loci. Traits matching two classes: 0. Traits with a null name: 1.

The hypothesis denominator is the resolved trait-gene pair: one assigned trait and one gene whose symbol is on at least one of that trait's clumps. Locus-gene links, which count a gene once per clump, are reported beside the pairs. Distinct genes are also reported, because the same gene can be named on related traits.

Resolved trait-gene pairs: 63380. Shared-gene hypothesis pairs: 39116. Resolved locus-gene links: 63629. Hypothesis links: 39277. Pairs mapped on more than one locus of the same trait: 249. Distinct resolved genes: 13728. Distinct hypothesis genes: 9275. Hypothesis genes that fall in two classes: 2911. Class totals of distinct genes cannot be added.

## Result

Of 63380 resolved trait-gene pairs, 39116 / 63380 (0.617) are shared-gene hypotheses. The gene was tested in at least one assigned tissue for 54638 pairs. Among those tested pairs, 39116 / 54638 (0.716) are eGenes in the class. 8742 pairs were not tested in the class. 60689 / 63380 (0.958) of resolved pairs are eGenes in at least one GTEx tissue of any kind, so an eGene call somewhere is the common case.

Blood: 46 traits, 17468 loci, 16116 / 28816 (0.559) hypothesis pairs, 4206 distinct hypothesis genes. Among pairs tested in Whole Blood, 16116 / 23814 (0.677). Among distinct resolved genes tested in Whole Blood, 4206 / 6337 (0.664). Among all genes tested in Whole Blood, 12360 / 20315 (0.608).

Brain: 73 traits, 22787 loci, 19147 / 25685 (0.745) hypothesis pairs, 6407 distinct hypothesis genes. Among pairs tested in at least one brain tissue, 19147 / 23179 (0.826). Of the distinct brain hypothesis genes, 1489 / 6407 (0.232) are eGenes in exactly one of the 13 brain tissues, and 300 are eGenes in all 13. A brain hypothesis needs an eGene call in any one of those tissues. That is a wider net than the one-tissue blood and liver rules.

Heart: 27 traits, 5039 loci, 2916 / 5162 (0.565) hypothesis pairs, 1610 distinct hypothesis genes. Among pairs tested in at least one heart tissue, 2916 / 4464 (0.653). Of the distinct heart hypothesis genes, 887 are eGenes in both heart tissues and 723 in only one.

Liver enzyme: 5 traits, 2533 loci, 937 / 3717 (0.252) hypothesis pairs, 617 distinct hypothesis genes. Among pairs tested in Liver, 937 / 3181 (0.295). Among distinct resolved genes tested in Liver, 617 / 2115 (0.292). Among all genes tested in Liver, 5734 / 22262 (0.258).

The genes that recur across the most traits in a class are in the gene table. Recurrence is mostly one clump per trait: for the top blood and liver genes, `max_loci_one_trait` is 1. `n_traits` counts related EFO terms. It is not a count of independent loci. Protein-coding genes are the majority of hypothesis pairs, and noncoding genes are kept: blood 14884 / 16116 (0.924), brain 16995 / 19147 (0.888), heart 2615 / 2916 (0.897), liver enzyme 849 / 937 (0.906).

The same rule does not promote every familiar GWAS gene. APOE is a hypothesis for 31 brain traits and is an eGene in 2 of 13 brain tissues. On 14 blood traits it is mapped and tested in Whole Blood and is not an eGene. HBB is mapped for 16 blood traits, tested in Whole Blood, and is not an eGene. GPT is mapped for 2 liver-enzyme traits and is not a liver eGene. PNPLA3 is mapped for 3 liver-enzyme traits and is not a liver eGene. ALPL is a liver eGene and a hypothesis for 2 liver-enzyme traits, including alkaline phosphatase measurement. GGT1 is a hypothesis for 1 liver-enzyme trait, serum gamma-glutamyl transferase measurement. NOS1AP is an eGene in both heart tissues and a hypothesis for QT interval and QRS duration.

## Caveat

A shared-gene hypothesis says only that the catalog named the gene on a clump and that GTEx calls the gene an eGene in an assigned tissue. It does not say the GWAS variant and the eQTL variant are the same variant, are in linkage disequilibrium, or sit inside the locus. No summary statistics and no linkage panel were used.

On the 47827 loci of assigned traits, the member span (`end - start`) has minimum 0 bp, median 2628 bp, and maximum 1987351 bp. 22609 of those loci have `start` equal to `end`. The maximum is under 2,000,000 bp. That is the span of the members, not a window grown out to 1 Mb.

The tissue rule is lexical. It is not an ontology. Blood pressure, blood protein amount, educational attainment, bilirubin measurement, and cardioembolic stroke stay unmapped. Coronary-artery and aortic-valve traits are scored on atrial appendage and left ventricle, which are not coronary artery or aortic valve. Parent terms such as hematological measurement and brain attribute are included and are not independent of the more specific traits. Family-history and biomarker names that matched the brain rule are included. They are not the same claim as a diagnosis measured in the participant's brain.

Raising the minimum from 50 loci to 100, inside the same name rule, leaves the pair fractions nearly unchanged: blood 15938 / 28496 (0.559), brain 18222 / 24412 (0.746), heart 2542 / 4510 (0.564), liver enzyme 937 / 3717 (0.252).

GTEx v8 `gene_name` and the gene-mart symbol differ for 809 hypothesis genes. Each of those ENSG ids has one GTEx gene name. The eGene join does not use either name. Examples are below. Ambiguous mart symbols, 483 symbols in the mart, are excluded rather than assigned to one gene.

## Tissue name map

QTL tissue is the `egenes.parquet` label. Expression tissue is the `median_tpm.parquet` label. Only the QTL label is used to call an eGene.

| class_name | qtl_tissue | expression_tissue |
| --- | --- | --- |
| blood | Whole_Blood | Whole Blood |
| brain | Brain_Amygdala | Brain - Amygdala |
| brain | Brain_Anterior_cingulate_cortex_BA24 | Brain - Anterior cingulate cortex (BA24) |
| brain | Brain_Caudate_basal_ganglia | Brain - Caudate (basal ganglia) |
| brain | Brain_Cerebellar_Hemisphere | Brain - Cerebellar Hemisphere |
| brain | Brain_Cerebellum | Brain - Cerebellum |
| brain | Brain_Cortex | Brain - Cortex |
| brain | Brain_Frontal_Cortex_BA9 | Brain - Frontal Cortex (BA9) |
| brain | Brain_Hippocampus | Brain - Hippocampus |
| brain | Brain_Hypothalamus | Brain - Hypothalamus |
| brain | Brain_Nucleus_accumbens_basal_ganglia | Brain - Nucleus accumbens (basal ganglia) |
| brain | Brain_Putamen_basal_ganglia | Brain - Putamen (basal ganglia) |
| brain | Brain_Spinal_cord_cervical_c-1 | Brain - Spinal cord (cervical c-1) |
| brain | Brain_Substantia_nigra | Brain - Substantia nigra |
| heart | Heart_Atrial_Appendage | Heart - Atrial Appendage |
| heart | Heart_Left_Ventricle | Heart - Left Ventricle |
| liver_enzyme | Liver | Liver |

## Class counts

| class_name | traits | loci | tokens | distance_or_unmatched_tokens | ambiguous_tokens | unique_symbol_tokens | locus_gene_links | hypothesis_links | trait_gene_pairs | hypothesis_pairs | tested_not_egene_pairs | not_tested_pairs | any_tissue_egene_pairs | distinct_resolved_genes | distinct_hypothesis_genes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| blood | 46 | 17468 | 41145 | 11980 | 270 | 28895 | 28895 | 16156 | 28816 | 16116 | 7698 | 5002 | 27711 | 8188 | 4206 |
| brain | 73 | 22787 | 40803 | 14755 | 208 | 25840 | 25840 | 19262 | 25685 | 19147 | 4032 | 2506 | 24413 | 8883 | 6407 |
| heart | 27 | 5039 | 7507 | 2296 | 43 | 5168 | 5168 | 2920 | 5162 | 2916 | 1548 | 698 | 4983 | 2957 | 1610 |
| liver_enzyme | 5 | 2533 | 5249 | 1476 | 47 | 3726 | 3726 | 939 | 3717 | 937 | 2244 | 536 | 3582 | 2503 | 617 |

Unmatched tokens in this table are the distance labels. The non-distance unmatched count is 0.

## Hypothesis genes by number of class tissues

One row per distinct hypothesis gene. Blood and liver enzyme each have one tissue, so every hypothesis gene in those classes has `n_sig_tissues` 1.

| class_name | n_sig_tissues | genes |
| --- | ---: | ---: |
| blood | 1 | 4206 |
| brain | 1 | 1489 |
| brain | 2 | 1285 |
| brain | 3 | 913 |
| brain | 4 | 570 |
| brain | 5 | 423 |
| brain | 6 | 324 |
| brain | 7 | 234 |
| brain | 8 | 215 |
| brain | 9 | 172 |
| brain | 10 | 165 |
| brain | 11 | 163 |
| brain | 12 | 154 |
| brain | 13 | 300 |
| heart | 1 | 723 |
| heart | 2 | 887 |
| liver_enzyme | 1 | 617 |

## eGene rates by QTL tissue

`tested_in_file` and `egenes_in_file` count genes in `egenes.parquet`, not the GWAS gene set. `resolved_genes_tested_here` and `resolved_genes_egene_here` are distinct genes that were mapped for at least one trait in the class. `trait_genes_egene_here` counts trait-gene pairs, so a gene named on several traits is counted several times.

| class_name | qtl_tissue | expression_tissue | tested_in_file | egenes_in_file | file_egene_rate | resolved_genes_tested | resolved_genes_egene | trait_genes_egene |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| blood | Whole_Blood | Whole Blood | 20315 | 12360 | 0.608 | 6337 | 4206 | 16116 |
| brain | Brain_Amygdala | Brain - Amygdala | 24069 | 3726 | 0.155 | 7317 | 1034 | 2937 |
| brain | Brain_Anterior_cingulate_cortex_BA24 | Brain - Anterior cingulate cortex (BA24) | 24342 | 5640 | 0.232 | 7313 | 1635 | 4473 |
| brain | Brain_Caudate_basal_ganglia | Brain - Caudate (basal ganglia) | 24718 | 8362 | 0.338 | 7349 | 2592 | 7376 |
| brain | Brain_Cerebellar_Hemisphere | Brain - Cerebellar Hemisphere | 25144 | 10027 | 0.399 | 7202 | 3086 | 8825 |
| brain | Brain_Cerebellum | Brain - Cerebellum | 25461 | 11240 | 0.441 | 7250 | 3533 | 10108 |
| brain | Brain_Cortex | Brain - Cortex | 24849 | 9082 | 0.365 | 7360 | 2829 | 7927 |
| brain | Brain_Frontal_Cortex_BA9 | Brain - Frontal Cortex (BA9) | 24676 | 7335 | 0.297 | 7362 | 2210 | 6157 |
| brain | Brain_Hippocampus | Brain - Hippocampus | 24420 | 5517 | 0.226 | 7343 | 1690 | 4953 |
| brain | Brain_Hypothalamus | Brain - Hypothalamus | 25096 | 5499 | 0.219 | 7433 | 1603 | 4339 |
| brain | Brain_Nucleus_accumbens_basal_ganglia | Brain - Nucleus accumbens (basal ganglia) | 24890 | 8198 | 0.329 | 7349 | 2464 | 7105 |
| brain | Brain_Putamen_basal_ganglia | Brain - Putamen (basal ganglia) | 23804 | 6902 | 0.290 | 7215 | 2185 | 6281 |
| brain | Brain_Spinal_cord_cervical_c-1 | Brain - Spinal cord (cervical c-1) | 24669 | 4483 | 0.182 | 7392 | 1383 | 4082 |
| brain | Brain_Substantia_nigra | Brain - Substantia nigra | 24044 | 3301 | 0.137 | 7326 | 990 | 2842 |
| heart | Heart_Atrial_Appendage | Heart - Atrial Appendage | 23194 | 10991 | 0.474 | 2465 | 1297 | 2350 |
| heart | Heart_Left_Ventricle | Heart - Left Ventricle | 21353 | 9642 | 0.452 | 2365 | 1200 | 2141 |
| liver_enzyme | Liver | Liver | 22262 | 5734 | 0.258 | 2115 | 617 | 937 |

## Genes mapped on the most traits

Twenty genes per class, ordered by number of assigned traits, then by the sum of loci across those traits, then by symbol. The full gene list is the `result_gene` query in the SQL. `n_sig_tissues` is how many of the class QTL tissues call the gene an eGene.

| class_name | rank | symbol | ensg | biotype | n_traits | n_loci_sum | max_loci_one_trait | n_sig_tissues | n_tested_tissues |
| --- | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| blood | 1 | MIR4435-2HG | ENSG00000172965 | lncRNA | 36 | 36 | 1 | 1 | 1 |
| blood | 2 | FADS2 | ENSG00000134824 | protein_coding | 34 | 34 | 1 | 1 | 1 |
| blood | 3 | IKZF1 | ENSG00000185811 | protein_coding | 33 | 33 | 1 | 1 | 1 |
| blood | 4 | ABO | ENSG00000175164 | protein_coding | 32 | 32 | 1 | 1 | 1 |
| blood | 5 | CCDC162P | ENSG00000203799 | transcribed_unitary_pseudogene | 31 | 31 | 1 | 1 | 1 |
| blood | 6 | MARK4 | ENSG00000007047 | protein_coding | 31 | 31 | 1 | 1 | 1 |
| blood | 7 | HBS1L | ENSG00000112339 | protein_coding | 30 | 30 | 1 | 1 | 1 |
| blood | 8 | GFI1B | ENSG00000165702 | protein_coding | 29 | 29 | 1 | 1 | 1 |
| blood | 9 | CCND3 | ENSG00000112576 | protein_coding | 27 | 27 | 1 | 1 | 1 |
| blood | 10 | HLA-B | ENSG00000234745 | protein_coding | 27 | 27 | 1 | 1 | 1 |
| blood | 11 | PVT1 | ENSG00000249859 | lncRNA | 26 | 26 | 1 | 1 | 1 |
| blood | 12 | SH2B3 | ENSG00000111252 | protein_coding | 26 | 26 | 1 | 1 | 1 |
| blood | 13 | THADA | ENSG00000115970 | protein_coding | 26 | 26 | 1 | 1 | 1 |
| blood | 14 | BRD3 | ENSG00000169925 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 15 | CHD7 | ENSG00000171316 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 16 | FADS1 | ENSG00000149485 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 17 | GCKR | ENSG00000084734 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 18 | JAK2 | ENSG00000096968 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 19 | ZFPM1 | ENSG00000179588 | protein_coding | 25 | 25 | 1 | 1 | 1 |
| blood | 20 | CD36 | ENSG00000135218 | protein_coding | 24 | 24 | 1 | 1 | 1 |
| brain | 1 | RBFOX1 | ENSG00000078328 | protein_coding | 33 | 46 | 2 | 2 | 13 |
| brain | 2 | APOE | ENSG00000130203 | protein_coding | 31 | 31 | 1 | 2 | 13 |
| brain | 3 | MAPT | ENSG00000186868 | protein_coding | 27 | 27 | 1 | 3 | 13 |
| brain | 4 | SLC39A8 | ENSG00000138821 | protein_coding | 26 | 26 | 1 | 2 | 13 |
| brain | 5 | LINC01122 | ENSG00000233723 | lncRNA | 25 | 29 | 2 | 1 | 13 |
| brain | 6 | FOXP2 | ENSG00000128573 | protein_coding | 25 | 26 | 2 | 2 | 13 |
| brain | 7 | DCC | ENSG00000187323 | protein_coding | 25 | 25 | 1 | 2 | 13 |
| brain | 8 | MACROD2 | ENSG00000172264 | protein_coding | 24 | 24 | 1 | 2 | 13 |
| brain | 9 | RERE | ENSG00000142599 | protein_coding | 24 | 24 | 1 | 7 | 13 |
| brain | 10 | ASTN2 | ENSG00000148219 | protein_coding | 23 | 24 | 2 | 5 | 13 |
| brain | 11 | MSRA | ENSG00000175806 | protein_coding | 22 | 23 | 2 | 6 | 13 |
| brain | 12 | ZFHX3 | ENSG00000140836 | protein_coding | 21 | 24 | 2 | 1 | 13 |
| brain | 13 | AUTS2 | ENSG00000158321 | protein_coding | 21 | 22 | 2 | 1 | 13 |
| brain | 14 | SOX5 | ENSG00000134532 | protein_coding | 21 | 22 | 2 | 1 | 13 |
| brain | 15 | TENM2 | ENSG00000145934 | protein_coding | 21 | 22 | 2 | 1 | 13 |
| brain | 16 | BANK1 | ENSG00000153064 | protein_coding | 21 | 21 | 1 | 3 | 13 |
| brain | 17 | RGS6 | ENSG00000182732 | protein_coding | 21 | 21 | 1 | 1 | 13 |
| brain | 18 | HDAC9 | ENSG00000048052 | protein_coding | 20 | 21 | 2 | 1 | 13 |
| brain | 19 | FURIN | ENSG00000140564 | protein_coding | 20 | 20 | 1 | 2 | 13 |
| brain | 20 | KANSL1 | ENSG00000120071 | protein_coding | 20 | 20 | 1 | 13 | 13 |
| heart | 1 | CDKN1A | ENSG00000124762 | protein_coding | 14 | 14 | 1 | 2 | 2 |
| heart | 2 | ZFPM2 | ENSG00000169946 | protein_coding | 12 | 12 | 1 | 1 | 2 |
| heart | 3 | FADS2 | ENSG00000134824 | protein_coding | 11 | 11 | 1 | 2 | 2 |
| heart | 4 | PHACTR1 | ENSG00000112137 | protein_coding | 11 | 11 | 1 | 1 | 2 |
| heart | 5 | PRKCA | ENSG00000154229 | protein_coding | 11 | 11 | 1 | 2 | 2 |
| heart | 6 | APOE | ENSG00000130203 | protein_coding | 10 | 10 | 1 | 1 | 2 |
| heart | 7 | COL4A2 | ENSG00000134871 | protein_coding | 10 | 10 | 1 | 1 | 2 |
| heart | 8 | SYNPO2L-AS1 | ENSG00000271848 | lncRNA | 10 | 10 | 1 | 1 | 2 |
| heart | 9 | TARID | ENSG00000227954 | lncRNA | 10 | 10 | 1 | 2 | 2 |
| heart | 10 | CASQ2 | ENSG00000118729 | protein_coding | 9 | 9 | 1 | 2 | 2 |
| heart | 11 | CDH13 | ENSG00000140945 | protein_coding | 9 | 9 | 1 | 2 | 2 |
| heart | 12 | ZPR1 | ENSG00000109917 | protein_coding | 9 | 9 | 1 | 1 | 2 |
| heart | 13 | ABO | ENSG00000175164 | protein_coding | 8 | 8 | 1 | 2 | 2 |
| heart | 14 | ADAMTS7 | ENSG00000136378 | protein_coding | 8 | 8 | 1 | 2 | 2 |
| heart | 15 | CAV1 | ENSG00000105974 | protein_coding | 8 | 8 | 1 | 1 | 2 |
| heart | 16 | CCDC141 | ENSG00000163492 | protein_coding | 8 | 8 | 1 | 1 | 2 |
| heart | 17 | CELSR2 | ENSG00000143126 | protein_coding | 8 | 8 | 1 | 2 | 2 |
| heart | 18 | FADS1 | ENSG00000149485 | protein_coding | 8 | 8 | 1 | 2 | 2 |
| heart | 19 | FHOD3 | ENSG00000134775 | protein_coding | 8 | 8 | 1 | 1 | 2 |
| heart | 20 | IL6R | ENSG00000160712 | protein_coding | 8 | 8 | 1 | 1 | 2 |
| liver_enzyme | 1 | ABO | ENSG00000175164 | protein_coding | 5 | 5 | 1 | 1 | 1 |
| liver_enzyme | 2 | CHEK2 | ENSG00000183765 | protein_coding | 5 | 5 | 1 | 1 | 1 |
| liver_enzyme | 3 | MLIP | ENSG00000146147 | protein_coding | 5 | 5 | 1 | 1 | 1 |
| liver_enzyme | 4 | NYNRIN | ENSG00000205978 | protein_coding | 5 | 5 | 1 | 1 | 1 |
| liver_enzyme | 5 | TRIB1AL | ENSG00000253111 | lncRNA | 5 | 5 | 1 | 1 | 1 |
| liver_enzyme | 6 | ALDH1A2 | ENSG00000128918 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 7 | APOH | ENSG00000091583 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 8 | ATP8B1-AS1 | ENSG00000267040 | lncRNA | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 9 | CD276 | ENSG00000103855 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 10 | CPN1 | ENSG00000120054 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 11 | DLG5 | ENSG00000151208 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 12 | EFHD1 | ENSG00000115468 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 13 | FCGRT | ENSG00000104870 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 14 | HSD17B13 | ENSG00000170509 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 15 | L3MBTL3 | ENSG00000198945 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 16 | LINC02068 | ENSG00000223387 | lncRNA | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 17 | LIPC | ENSG00000166035 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 18 | MLYCD | ENSG00000103150 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 19 | NCK1 | ENSG00000158092 | protein_coding | 4 | 4 | 1 | 1 | 1 |
| liver_enzyme | 20 | NECTIN2 | ENSG00000130202 | protein_coding | 4 | 4 | 1 | 1 | 1 |

## Biotype of hypothesis pairs

| class_name | biotype | hypothesis_trait_genes | distinct_genes |
| --- | --- | ---: | ---: |
| blood | protein_coding | 14884 | 3863 |
| blood | lncRNA | 1113 | 289 |
| blood | transcribed_unprocessed_pseudogene | 56 | 31 |
| blood | transcribed_unitary_pseudogene | 42 | 8 |
| blood | processed_pseudogene | 12 | 7 |
| blood | transcribed_processed_pseudogene | 4 | 4 |
| blood | unprocessed_pseudogene | 3 | 2 |
| blood | miRNA | 1 | 1 |
| blood | IG_C_pseudogene | 1 | 1 |
| brain | protein_coding | 16995 | 5629 |
| brain | lncRNA | 2061 | 724 |
| brain | transcribed_unprocessed_pseudogene | 47 | 26 |
| brain | transcribed_unitary_pseudogene | 29 | 14 |
| brain | processed_pseudogene | 8 | 8 |
| brain | transcribed_processed_pseudogene | 4 | 3 |
| brain | unprocessed_pseudogene | 1 | 1 |
| brain | misc_RNA | 1 | 1 |
| brain | snoRNA | 1 | 1 |
| heart | protein_coding | 2615 | 1443 |
| heart | lncRNA | 286 | 156 |
| heart | transcribed_unprocessed_pseudogene | 10 | 7 |
| heart | processed_pseudogene | 2 | 1 |
| heart | transcribed_processed_pseudogene | 1 | 1 |
| heart | transcribed_unitary_pseudogene | 1 | 1 |
| heart | unprocessed_pseudogene | 1 | 1 |
| liver_enzyme | protein_coding | 849 | 552 |
| liver_enzyme | lncRNA | 80 | 57 |
| liver_enzyme | transcribed_unprocessed_pseudogene | 4 | 4 |
| liver_enzyme | processed_pseudogene | 2 | 2 |
| liver_enzyme | misc_RNA | 1 | 1 |
| liver_enzyme | transcribed_unitary_pseudogene | 1 | 1 |

## Definition checks

These rows apply the same rule to a fixed symbol list. They are checks, not a ranking. A symbol can be a hypothesis in one class and not in another.

| symbol | class_name | mapped_traits | hypothesis_traits | class_tissues_tested | class_tissues_egene |
| --- | --- | ---: | ---: | ---: | ---: |
| APOE | brain | 31 | 31 | 13 | 2 |
| APOE | blood | 14 | 0 | 1 | 0 |
| HBB | blood | 16 | 0 | 1 | 0 |
| GPT | liver_enzyme | 2 | 0 | 1 | 0 |
| PNPLA3 | liver_enzyme | 3 | 0 | 1 | 0 |
| ALPL | liver_enzyme | 2 | 2 | 1 | 1 |
| GGT1 | liver_enzyme | 1 | 1 | 1 | 1 |
| NOS1AP | heart | 2 | 2 | 2 | 2 |

## Symbol and GTEx gene name

809 hypothesis genes have a gene-mart symbol different from the GTEx v8 `gene_name` on that ENSG. Each example below is two names for one ENSG, including later symbols such as ADA2 for CECR1 and clone or antisense names. The eGene row is that ENSG. It was not chosen by matching the name. The twelve examples are alphabetical by symbol:

| symbol | gene_name | ensg | n_classes |
| --- | --- | --- | ---: |
| ABITRAM | FAM206A | ENSG00000119328 | 1 |
| ABRAXAS2 | FAM175B | ENSG00000165660 | 1 |
| ABTB3 | BTBD11 | ENSG00000151136 | 1 |
| ACSL6-AS1 | AC034228.3 | ENSG00000223548 | 1 |
| ACTMAP | C19orf54 | ENSG00000188493 | 1 |
| ACTN1-DT | ACTN1-AS1 | ENSG00000259062 | 1 |
| ACTN3-AS1 | CTD-3074O7.2 | ENSG00000250105 | 1 |
| ADA2 | CECR1 | ENSG00000093072 | 2 |
| ADAMTSL4-AS2 | RP11-54A4.2 | ENSG00000237781 | 1 |
| ADGRB3-DT | RP3-525N10.2 | ENSG00000230910 | 1 |
| ADGRF2P | ADGRF2 | ENSG00000164393 | 1 |
| ADISSP | C20orf27 | ENSG00000101220 | 2 |

## Assigned traits

Each row is one EFO trait with at least 50 clumps and one tissue class. `hypothesis_genes / resolved_genes` is the pair rate. `not_tested` genes had no eGene row in the class. `loci_with_hypothesis` is clumps with at least one hypothesis gene.

| class_name | n_loci | resolved_genes | hypothesis_genes | pair_rate | tested_not_egene | not_tested | any_tissue_egene | loci_with_hypothesis | loci_with_resolved_gene | mapped_trait | trait_uri |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| blood | 1204 | 3730 | 1927 | 0.517 | 999 | 804 | 3556 | 812 | 1059 | hematological measurement | http://www.ebi.ac.uk/efo/EFO_0004503 |
| blood | 853 | 1514 | 813 | 0.537 | 439 | 262 | 1461 | 487 | 715 | platelet count | http://www.ebi.ac.uk/efo/EFO_0004309 |
| blood | 748 | 1283 | 702 | 0.547 | 327 | 254 | 1235 | 423 | 627 | leukocyte quantity | http://purl.obolibrary.org/obo/OBA_VT0000217 |
| blood | 735 | 1265 | 689 | 0.545 | 338 | 238 | 1222 | 418 | 608 | erythrocyte count | http://www.ebi.ac.uk/efo/EFO_0004305 |
| blood | 665 | 1222 | 594 | 0.486 | 327 | 301 | 1158 | 349 | 558 | hemoglobin measurement | http://www.ebi.ac.uk/efo/EFO_0004509 |
| blood | 651 | 1356 | 747 | 0.551 | 404 | 205 | 1300 | 408 | 567 | erythrocyte volume | http://purl.obolibrary.org/obo/OBA_0003460 |
| blood | 641 | 998 | 494 | 0.495 | 274 | 230 | 962 | 319 | 532 | hematocrit | http://www.ebi.ac.uk/efo/EFO_0004348 |
| blood | 603 | 971 | 547 | 0.563 | 234 | 190 | 934 | 347 | 499 | lymphocyte count | http://www.ebi.ac.uk/efo/EFO_0004587 |
| blood | 563 | 896 | 505 | 0.564 | 209 | 182 | 868 | 318 | 471 | neutrophil count | http://www.ebi.ac.uk/efo/EFO_0004833 |
| blood | 561 | 1013 | 588 | 0.580 | 254 | 171 | 966 | 352 | 473 | eosinophil count | http://www.ebi.ac.uk/efo/EFO_0004842 |
| blood | 540 | 1074 | 576 | 0.536 | 321 | 177 | 1029 | 343 | 474 | mean corpuscular hemoglobin | http://www.ebi.ac.uk/efo/EFO_0004527 |
| blood | 540 | 877 | 502 | 0.572 | 235 | 140 | 843 | 320 | 453 | monocyte count | http://www.ebi.ac.uk/efo/EFO_0005091 |
| blood | 504 | 990 | 570 | 0.576 | 267 | 153 | 957 | 337 | 446 | platelet volume | http://purl.obolibrary.org/obo/OBA_0003277 |
| blood | 500 | 730 | 426 | 0.584 | 213 | 91 | 708 | 290 | 420 | mean corpuscular hemoglobin concentration | http://www.ebi.ac.uk/efo/EFO_0004528 |
| blood | 494 | 494 | 263 | 0.532 | 135 | 96 | 479 | 230 | 383 | red blood cell density | http://purl.obolibrary.org/obo/OBA_2045260 |
| blood | 478 | 961 | 557 | 0.580 | 283 | 121 | 935 | 330 | 434 | Red cell distribution width | http://www.ebi.ac.uk/efo/EFO_0009188 |
| blood | 477 | 503 | 294 | 0.584 | 130 | 79 | 490 | 243 | 366 | myeloid leukocyte count | http://purl.obolibrary.org/obo/OBA_2045223 |
| blood | 470 | 672 | 381 | 0.567 | 186 | 105 | 647 | 259 | 389 | platelet crit | http://www.ebi.ac.uk/efo/EFO_0007985 |
| blood | 447 | 556 | 314 | 0.565 | 139 | 103 | 534 | 233 | 345 | eosinophil percentage of leukocytes | http://www.ebi.ac.uk/efo/EFO_0007991 |
| blood | 447 | 713 | 406 | 0.569 | 211 | 96 | 688 | 265 | 371 | reticulocyte amount | http://purl.obolibrary.org/obo/OBA_0002272 |
| blood | 441 | 688 | 416 | 0.605 | 186 | 86 | 669 | 262 | 364 | reticulocyte count | http://www.ebi.ac.uk/efo/EFO_0007986 |
| blood | 415 | 575 | 344 | 0.598 | 162 | 69 | 559 | 239 | 346 | mean reticulocyte volume | http://www.ebi.ac.uk/efo/EFO_0010701 |
| blood | 394 | 584 | 352 | 0.603 | 157 | 75 | 563 | 240 | 335 | platelet component distribution width | http://www.ebi.ac.uk/efo/EFO_0007984 |
| blood | 353 | 440 | 271 | 0.616 | 100 | 69 | 427 | 201 | 288 | neutrophil percentage of leukocytes | http://www.ebi.ac.uk/efo/EFO_0007990 |
| blood | 343 | 440 | 280 | 0.636 | 102 | 58 | 423 | 207 | 281 | lymphocyte percentage of leukocytes | http://www.ebi.ac.uk/efo/EFO_0007993 |
| blood | 338 | 458 | 272 | 0.594 | 106 | 80 | 443 | 203 | 275 | basophil count | http://www.ebi.ac.uk/efo/EFO_0005090 |
| blood | 337 | 467 | 302 | 0.647 | 105 | 60 | 453 | 206 | 272 | monocyte percentage of leukocytes | http://www.ebi.ac.uk/efo/EFO_0007989 |
| blood | 302 | 327 | 214 | 0.654 | 75 | 38 | 319 | 178 | 245 | neutrophil measurement | http://www.ebi.ac.uk/efo/EFO_0803548 |
| blood | 294 | 764 | 455 | 0.596 | 185 | 124 | 732 | 207 | 263 | neutrophil-to-lymphocyte ratio | http://www.ebi.ac.uk/efo/EFO_0008447 |
| blood | 281 | 292 | 186 | 0.637 | 74 | 32 | 281 | 158 | 222 | lymphocyte amount | http://purl.obolibrary.org/obo/OBA_0003352 |
| blood | 167 | 216 | 136 | 0.630 | 47 | 33 | 206 | 103 | 139 | platelet quantity | http://purl.obolibrary.org/obo/OBA_VT0003179 |
| blood | 153 | 208 | 118 | 0.567 | 55 | 35 | 195 | 88 | 134 | platelet-to-lymphocyte ratio | http://www.ebi.ac.uk/efo/EFO_0008446 |
| blood | 149 | 152 | 80 | 0.526 | 45 | 27 | 143 | 67 | 119 | basophil measurement | http://www.ebi.ac.uk/efo/EFO_0803539 |
| blood | 146 | 129 | 80 | 0.620 | 31 | 18 | 125 | 68 | 99 | granulocyte percentage of myeloid white cells | http://www.ebi.ac.uk/efo/EFO_0007997 |
| blood | 141 | 145 | 86 | 0.593 | 32 | 27 | 142 | 78 | 111 | basophil percentage of leukocytes | http://www.ebi.ac.uk/efo/EFO_0007992 |
| blood | 139 | 122 | 72 | 0.590 | 32 | 18 | 121 | 66 | 98 | granulocyte count | http://www.ebi.ac.uk/efo/EFO_0007987 |
| blood | 134 | 117 | 45 | 0.385 | 36 | 36 | 109 | 38 | 88 | fetal hemoglobin measurement | http://www.ebi.ac.uk/efo/EFO_0004576 |
| blood | 129 | 123 | 79 | 0.642 | 31 | 13 | 119 | 71 | 103 | eosinophil percentage of granulocytes | http://www.ebi.ac.uk/efo/EFO_0007996 |
| blood | 126 | 118 | 77 | 0.653 | 28 | 13 | 114 | 66 | 95 | neutrophil percentage of granulocytes | http://www.ebi.ac.uk/efo/EFO_0007994 |
| blood | 125 | 168 | 99 | 0.589 | 43 | 26 | 159 | 73 | 97 | lymphocyte:monocyte ratio | http://www.ebi.ac.uk/efo/EFO_0600088 |
| blood | 117 | 145 | 79 | 0.545 | 45 | 21 | 134 | 61 | 90 | erythrocyte attribute | http://purl.obolibrary.org/obo/OBA_2045276 |
| blood | 80 | 78 | 42 | 0.538 | 25 | 11 | 74 | 39 | 65 | immature platelet measurement | http://www.ebi.ac.uk/efo/EFO_0803541 |
| blood | 72 | 72 | 40 | 0.556 | 22 | 10 | 69 | 37 | 61 | immature platelet count | http://www.ebi.ac.uk/efo/EFO_0803544 |
| blood | 61 | 60 | 41 | 0.683 | 17 | 2 | 58 | 35 | 47 | eosinophil measurement | http://www.ebi.ac.uk/efo/EFO_0803540 |
| blood | 57 | 54 | 26 | 0.481 | 18 | 10 | 51 | 22 | 44 | platelet larger cell ratio | http://www.ebi.ac.uk/efo/EFO_0009186 |
| blood | 53 | 56 | 29 | 0.518 | 14 | 13 | 51 | 23 | 38 | clonal hematopoiesis | http://www.ebi.ac.uk/efo/EFO_0010819 |
| brain | 1208 | 1740 | 1284 | 0.738 | 290 | 166 | 1655 | 731 | 870 | Alzheimer disease | http://purl.obolibrary.org/obo/MONDO_0004975 |
| brain | 976 | 1224 | 902 | 0.737 | 187 | 135 | 1154 | 610 | 729 | cognitive function measurement | http://www.ebi.ac.uk/efo/EFO_0008354 |
| brain | 969 | 1001 | 732 | 0.731 | 164 | 105 | 946 | 528 | 647 | memory performance | http://www.ebi.ac.uk/efo/EFO_0004874 |
| brain | 949 | 1201 | 884 | 0.736 | 194 | 123 | 1135 | 551 | 663 | major depressive disorder | http://purl.obolibrary.org/obo/MONDO_0002009 |
| brain | 944 | 1431 | 1048 | 0.732 | 241 | 142 | 1346 | 601 | 714 | schizophrenia | http://purl.obolibrary.org/obo/MONDO_0005090 |
| brain | 820 | 1052 | 798 | 0.759 | 157 | 97 | 1006 | 558 | 654 | brain attribute | http://purl.obolibrary.org/obo/OBA_2045210 |
| brain | 791 | 901 | 664 | 0.737 | 146 | 91 | 847 | 476 | 561 | attention deficit-hyperactivity disorder | http://purl.obolibrary.org/obo/MONDO_0007743 |
| brain | 759 | 877 | 656 | 0.748 | 151 | 70 | 844 | 432 | 517 | family history of Alzheimer’s disease | http://www.ebi.ac.uk/efo/EFO_0009268 |
| brain | 745 | 964 | 699 | 0.725 | 171 | 94 | 910 | 491 | 593 | mathematical ability | http://www.ebi.ac.uk/efo/EFO_0004875 |
| brain | 727 | 672 | 518 | 0.771 | 100 | 54 | 648 | 401 | 490 | dementia | http://purl.obolibrary.org/obo/MONDO_0001627 |
| brain | 707 | 1124 | 836 | 0.744 | 180 | 108 | 1065 | 451 | 532 | intelligence | http://www.ebi.ac.uk/efo/EFO_0004337 |
| brain | 646 | 814 | 631 | 0.775 | 117 | 66 | 784 | 427 | 498 | neuroimaging measurement | http://www.ebi.ac.uk/efo/EFO_0004346 |
| brain | 644 | 827 | 642 | 0.776 | 119 | 66 | 795 | 431 | 505 | cerebral cortex area attribute | http://purl.obolibrary.org/obo/OBA_2045267 |
| brain | 633 | 938 | 662 | 0.706 | 152 | 124 | 877 | 394 | 483 | insomnia | http://purl.obolibrary.org/obo/MONDO_0013600 |
| brain | 592 | 680 | 543 | 0.799 | 84 | 53 | 655 | 387 | 447 | cortical thickness | http://www.ebi.ac.uk/efo/EFO_0004840 |
| brain | 590 | 675 | 504 | 0.747 | 101 | 70 | 640 | 344 | 420 | neuroticism measurement | http://www.ebi.ac.uk/efo/EFO_0007660 |
| brain | 589 | 589 | 418 | 0.710 | 108 | 63 | 562 | 326 | 412 | depressive symptom measurement | http://www.ebi.ac.uk/efo/EFO_0007006 |
| brain | 567 | 730 | 558 | 0.764 | 101 | 71 | 700 | 374 | 433 | white matter microstructure measurement | http://www.ebi.ac.uk/efo/EFO_0005674 |
| brain | 548 | 620 | 474 | 0.765 | 90 | 56 | 588 | 342 | 402 | bipolar disorder | http://purl.obolibrary.org/obo/MONDO_0004985 |
| brain | 464 | 462 | 332 | 0.719 | 82 | 48 | 433 | 256 | 328 | sleep duration trait | http://purl.obolibrary.org/obo/OBA_2040171 |
| brain | 387 | 508 | 378 | 0.744 | 81 | 49 | 475 | 241 | 285 | autism spectrum disorder | http://purl.obolibrary.org/obo/MONDO_0005258 |
| brain | 377 | 320 | 242 | 0.756 | 50 | 28 | 309 | 205 | 252 | Ischemic stroke | http://purl.obolibrary.org/obo/HP_0002140 |
| brain | 345 | 404 | 303 | 0.750 | 59 | 42 | 387 | 212 | 255 | brain volume | http://purl.obolibrary.org/obo/OBA_2050009 |
| brain | 342 | 412 | 313 | 0.760 | 63 | 36 | 398 | 217 | 258 | multiple sclerosis | http://purl.obolibrary.org/obo/MONDO_0005301 |
| brain | 325 | 196 | 149 | 0.760 | 32 | 15 | 189 | 131 | 167 | cognitive decline measurement | http://www.ebi.ac.uk/efo/EFO_0007710 |
| brain | 315 | 315 | 231 | 0.733 | 52 | 32 | 299 | 182 | 223 | migraine disorder | http://purl.obolibrary.org/obo/MONDO_0005277 |
| brain | 311 | 318 | 250 | 0.786 | 41 | 27 | 306 | 194 | 236 | brain connectivity attribute | http://purl.obolibrary.org/obo/OBA_2045295 |
| brain | 253 | 298 | 229 | 0.768 | 38 | 31 | 283 | 157 | 185 | Parkinson disease | http://purl.obolibrary.org/obo/MONDO_0005180 |
| brain | 250 | 226 | 154 | 0.681 | 50 | 22 | 217 | 134 | 177 | insomnia measurement | http://www.ebi.ac.uk/efo/EFO_0007876 |
| brain | 243 | 206 | 140 | 0.680 | 33 | 33 | 191 | 119 | 151 | memory impairment | http://www.ebi.ac.uk/efo/EFO_0001072 |
| brain | 232 | 185 | 137 | 0.741 | 32 | 16 | 174 | 122 | 148 | obsessive-compulsive disorder | http://purl.obolibrary.org/obo/MONDO_0008114 |
| brain | 207 | 156 | 119 | 0.763 | 26 | 11 | 154 | 103 | 125 | late-onset Alzheimers disease | http://www.ebi.ac.uk/efo/EFO_1001870 |
| brain | 196 | 154 | 118 | 0.766 | 24 | 12 | 145 | 105 | 129 | cerebellar volume measurement | http://www.ebi.ac.uk/efo/EFO_0020864 |
| brain | 194 | 264 | 201 | 0.761 | 40 | 23 | 251 | 130 | 150 | amount of iron in brain | http://purl.obolibrary.org/obo/OBA_2045560 |
| brain | 186 | 140 | 114 | 0.814 | 16 | 10 | 135 | 99 | 118 | hippocampal volume | http://www.ebi.ac.uk/efo/EFO_0005035 |
| brain | 180 | 135 | 92 | 0.681 | 26 | 17 | 127 | 79 | 106 | psychiatric disorder | http://purl.obolibrary.org/obo/MONDO_0002025 |
| brain | 167 | 150 | 113 | 0.753 | 25 | 12 | 141 | 100 | 119 | Tourette syndrome | http://purl.obolibrary.org/obo/MONDO_0007661 |
| brain | 164 | 148 | 112 | 0.757 | 27 | 9 | 143 | 95 | 112 | stroke disorder | http://purl.obolibrary.org/obo/MONDO_0005098 |
| brain | 164 | 155 | 122 | 0.787 | 18 | 15 | 147 | 90 | 107 | white matter hyperintensity measurement | http://www.ebi.ac.uk/efo/EFO_0005665 |
| brain | 151 | 116 | 85 | 0.733 | 22 | 9 | 109 | 72 | 91 | sleep characteristic | http://purl.obolibrary.org/obo/OBA_2040170 |
| brain | 146 | 144 | 110 | 0.764 | 22 | 12 | 134 | 78 | 95 | Anxiety | http://purl.obolibrary.org/obo/HP_0000739 |
| brain | 138 | 118 | 92 | 0.780 | 14 | 12 | 110 | 82 | 99 | anxiety measurement | http://www.ebi.ac.uk/efo/EFO_0009863 |
| brain | 134 | 114 | 80 | 0.702 | 17 | 17 | 103 | 72 | 94 | obsessive-compulsive symptom measurement | http://www.ebi.ac.uk/efo/EFO_0007802 |
| brain | 133 | 106 | 80 | 0.755 | 18 | 8 | 106 | 73 | 91 | schizophrenia symptom severity measurement | http://www.ebi.ac.uk/efo/EFO_0007927 |
| brain | 121 | 80 | 57 | 0.713 | 16 | 7 | 80 | 54 | 68 | neuropsychological test | http://www.ebi.ac.uk/efo/EFO_0003926 |
| brain | 116 | 121 | 98 | 0.810 | 14 | 9 | 120 | 73 | 84 | amygdala volume | http://purl.obolibrary.org/obo/OBA_2050034 |
| brain | 106 | 84 | 67 | 0.798 | 9 | 8 | 82 | 58 | 70 | amyotrophic lateral sclerosis | http://purl.obolibrary.org/obo/MONDO_0004976 |
| brain | 105 | 95 | 76 | 0.800 | 12 | 7 | 91 | 55 | 65 | Alzheimer's disease biomarker measurement | http://www.ebi.ac.uk/efo/EFO_0006514 |
| brain | 105 | 76 | 58 | 0.763 | 11 | 7 | 73 | 55 | 68 | vascular dementia | http://purl.obolibrary.org/obo/MONDO_0004648 |
| brain | 101 | 75 | 57 | 0.760 | 8 | 10 | 73 | 51 | 65 | anxiety disorder measurement | http://www.ebi.ac.uk/efo/EFO_0007795 |
| brain | 100 | 71 | 60 | 0.845 | 4 | 7 | 68 | 56 | 64 | sporadic amyotrophic lateral sclerosis | http://purl.obolibrary.org/obo/MONDO_0005145 |
| brain | 99 | 71 | 48 | 0.676 | 17 | 6 | 66 | 45 | 57 | major depressive episode | http://www.ebi.ac.uk/efo/EFO_0007634 |
| brain | 96 | 86 | 70 | 0.814 | 5 | 11 | 81 | 61 | 73 | brain aneurysm | http://purl.obolibrary.org/obo/MONDO_0005291 |
| brain | 91 | 99 | 71 | 0.717 | 16 | 12 | 94 | 45 | 59 | cognitive domain measurement | http://www.ebi.ac.uk/efo/EFO_0009098 |
| brain | 88 | 72 | 58 | 0.806 | 7 | 7 | 71 | 46 | 54 | brain physiology trait | http://purl.obolibrary.org/obo/OBA_VT0015058 |
| brain | 87 | 69 | 49 | 0.710 | 13 | 7 | 66 | 45 | 62 | anxiety disorder | http://purl.obolibrary.org/obo/MONDO_0005618 |
| brain | 85 | 53 | 39 | 0.736 | 6 | 8 | 50 | 37 | 47 | short-term memory | http://www.ebi.ac.uk/efo/EFO_0004335 |
| brain | 80 | 53 | 35 | 0.660 | 10 | 8 | 51 | 30 | 44 | cerebral amyloid angiopathy | http://purl.obolibrary.org/obo/MONDO_0005620 |
| brain | 80 | 65 | 44 | 0.677 | 12 | 9 | 60 | 41 | 54 | large artery stroke | http://purl.obolibrary.org/obo/MONDO_0005490 |
| brain | 78 | 73 | 54 | 0.740 | 11 | 8 | 70 | 49 | 56 | bipolar I disorder | http://www.ebi.ac.uk/efo/EFO_0009963 |
| brain | 78 | 41 | 28 | 0.683 | 5 | 8 | 38 | 22 | 34 | epilepsy | http://purl.obolibrary.org/obo/MONDO_0005027 |
| brain | 75 | 54 | 37 | 0.685 | 13 | 4 | 51 | 33 | 48 | Lewy body dementia | http://purl.obolibrary.org/obo/MONDO_0007488 |
| brain | 73 | 79 | 60 | 0.759 | 7 | 12 | 73 | 46 | 56 | white matter integrity | http://www.ebi.ac.uk/efo/EFO_0004641 |
| brain | 72 | 44 | 32 | 0.727 | 7 | 5 | 43 | 31 | 42 | neurofibrillary tangles measurement | http://www.ebi.ac.uk/efo/EFO_0006797 |
| brain | 71 | 51 | 35 | 0.686 | 13 | 3 | 47 | 33 | 47 | neurotic disorder | http://purl.obolibrary.org/obo/MONDO_0005379 |
| brain | 69 | 39 | 27 | 0.692 | 5 | 7 | 35 | 27 | 36 | stroke outcome severity measurement | http://www.ebi.ac.uk/efo/EFO_0009603 |
| brain | 67 | 73 | 57 | 0.781 | 7 | 9 | 68 | 43 | 52 | brain age measurement | http://www.ebi.ac.uk/efo/EFO_0010602 |
| brain | 59 | 43 | 29 | 0.674 | 8 | 6 | 42 | 27 | 36 | small vessel stroke | http://www.ebi.ac.uk/efo/EFO_1001504 |
| brain | 58 | 51 | 34 | 0.667 | 7 | 10 | 48 | 27 | 39 | cerebral amyloid deposition measurement | http://www.ebi.ac.uk/efo/EFO_0007707 |
| brain | 57 | 44 | 32 | 0.727 | 7 | 5 | 42 | 31 | 36 | psychotic symptoms | http://www.ebi.ac.uk/efo/EFO_0005940 |
| brain | 55 | 42 | 33 | 0.786 | 6 | 3 | 40 | 27 | 35 | grey matter volume measurement | http://www.ebi.ac.uk/efo/EFO_0005420 |
| brain | 54 | 33 | 26 | 0.788 | 5 | 2 | 31 | 25 | 30 | narcolepsy | http://purl.obolibrary.org/obo/MONDO_0021107 |
| brain | 53 | 38 | 27 | 0.711 | 10 | 1 | 36 | 26 | 34 | ADHD symptom measurement | http://www.ebi.ac.uk/efo/EFO_0007860 |
| heart | 607 | 831 | 459 | 0.552 | 237 | 135 | 801 | 290 | 446 | coronary artery disorder | http://purl.obolibrary.org/obo/MONDO_0005010 |
| heart | 569 | 687 | 387 | 0.563 | 213 | 87 | 663 | 297 | 453 | atrial fibrillation | http://purl.obolibrary.org/obo/MONDO_0004981 |
| heart | 396 | 393 | 235 | 0.598 | 125 | 33 | 384 | 194 | 294 | heart rate | http://purl.obolibrary.org/obo/OBA_1001087 |
| heart | 329 | 346 | 201 | 0.581 | 99 | 46 | 330 | 153 | 239 | heart failure | http://purl.obolibrary.org/obo/MONDO_0005252 |
| heart | 305 | 296 | 179 | 0.605 | 89 | 28 | 289 | 147 | 220 | QT interval | http://www.ebi.ac.uk/efo/EFO_0004682 |
| heart | 297 | 306 | 171 | 0.559 | 99 | 36 | 293 | 146 | 226 | aortic stenosis | http://www.ebi.ac.uk/efo/EFO_0000266 |
| heart | 295 | 259 | 148 | 0.571 | 74 | 37 | 253 | 133 | 214 | PR interval | http://www.ebi.ac.uk/efo/EFO_0004462 |
| heart | 276 | 376 | 204 | 0.543 | 116 | 56 | 367 | 141 | 212 | myocardial infarction | http://purl.obolibrary.org/obo/MONDO_0005068 |
| heart | 260 | 192 | 93 | 0.484 | 53 | 46 | 178 | 85 | 162 | coronary artery calcification | http://www.ebi.ac.uk/efo/EFO_0004723 |
| heart | 194 | 158 | 87 | 0.551 | 54 | 17 | 155 | 79 | 136 | QRS duration | http://purl.obolibrary.org/obo/OBA_1001086 |
| heart | 183 | 146 | 89 | 0.610 | 50 | 7 | 142 | 85 | 135 | JT interval | http://www.ebi.ac.uk/efo/EFO_0007885 |
| heart | 133 | 83 | 37 | 0.446 | 32 | 14 | 81 | 36 | 79 | cardiac troponin T measurement | http://www.ebi.ac.uk/efo/EFO_0005043 |
| heart | 127 | 113 | 72 | 0.637 | 25 | 16 | 108 | 62 | 93 | magnetic resonance imaging of the heart | http://www.ebi.ac.uk/efo/EFO_0022611 |
| heart | 121 | 110 | 69 | 0.627 | 29 | 12 | 105 | 58 | 89 | hypertrophic cardiomyopathy | http://purl.obolibrary.org/obo/MONDO_0005045 |
| heart | 117 | 130 | 73 | 0.562 | 29 | 28 | 120 | 64 | 105 | ventricular septal defect | http://purl.obolibrary.org/obo/MONDO_0002070 |
| heart | 114 | 84 | 38 | 0.452 | 35 | 11 | 84 | 37 | 75 | cardiac troponin I measurement | http://www.ebi.ac.uk/efo/EFO_0010071 |
| heart | 90 | 80 | 36 | 0.450 | 35 | 9 | 76 | 32 | 64 | left ventricular structural measurement | http://www.ebi.ac.uk/efo/EFO_0008205 |
| heart | 73 | 49 | 22 | 0.449 | 10 | 17 | 45 | 22 | 44 | conotruncal heart malformations | http://purl.obolibrary.org/obo/MONDO_0016581 |
| heart | 72 | 98 | 60 | 0.612 | 22 | 16 | 94 | 42 | 60 | angina pectoris | http://www.ebi.ac.uk/efo/EFO_0003913 |
| heart | 72 | 73 | 39 | 0.534 | 25 | 9 | 69 | 34 | 60 | atrial flutter | http://purl.obolibrary.org/obo/MONDO_0005310 |
| heart | 67 | 60 | 37 | 0.617 | 18 | 5 | 59 | 35 | 52 | QRS-T angle | http://www.ebi.ac.uk/efo/EFO_0020097 |
| heart | 65 | 61 | 35 | 0.574 | 18 | 8 | 60 | 30 | 50 | coronary atherosclerosis | http://purl.obolibrary.org/obo/MONDO_0021661 |
| heart | 64 | 63 | 38 | 0.603 | 16 | 9 | 61 | 32 | 51 | myocardial ischemia | http://purl.obolibrary.org/obo/MONDO_0024644 |
| heart | 58 | 57 | 36 | 0.632 | 17 | 4 | 56 | 31 | 45 | QRS complex | http://www.ebi.ac.uk/efo/EFO_0005054 |
| heart | 53 | 38 | 26 | 0.684 | 11 | 1 | 38 | 26 | 37 | left ventricular function | http://www.ebi.ac.uk/efo/EFO_0004295 |
| heart | 51 | 38 | 24 | 0.632 | 9 | 5 | 37 | 22 | 33 | aortic valve calcification | http://purl.obolibrary.org/obo/MONDO_0005463 |
| heart | 51 | 35 | 21 | 0.600 | 8 | 6 | 35 | 19 | 30 | sudden cardiac arrest | http://www.ebi.ac.uk/efo/EFO_0004278 |
| liver_enzyme | 683 | 951 | 238 | 0.250 | 562 | 151 | 916 | 197 | 538 | serum alanine aminotransferase amount | http://purl.obolibrary.org/obo/OBA_2050062 |
| liver_enzyme | 682 | 946 | 231 | 0.244 | 582 | 133 | 915 | 192 | 545 | aspartate aminotransferase measurement | http://www.ebi.ac.uk/efo/EFO_0004736 |
| liver_enzyme | 659 | 935 | 226 | 0.242 | 574 | 135 | 900 | 185 | 534 | serum gamma-glutamyl transferase measurement | http://www.ebi.ac.uk/efo/EFO_0004532 |
| liver_enzyme | 403 | 741 | 204 | 0.275 | 432 | 105 | 709 | 142 | 344 | alkaline phosphatase measurement | http://www.ebi.ac.uk/efo/EFO_0004533 |
| liver_enzyme | 106 | 144 | 38 | 0.264 | 94 | 12 | 142 | 33 | 85 | aspartate aminotransferase to alanine aminotransferase ratio | http://www.ebi.ac.uk/efo/EFO_0010934 |

## Unmapped traits

933 traits have at least 50 clumped loci and no tissue class. No shared-gene hypothesis is claimed for them.

One of these rows has a null `mapped_trait` in `traits.parquet`: `http://www.ebi.ac.uk/efo/EFO_0010726` (95 loci). The trait row is present. The name is null, so the name rule does not assign it.

| n_loci | mapped_trait | trait_uri |
| ---: | --- | --- |
| 1752 | body height | http://purl.obolibrary.org/obo/OBA_VT0001253 |
| 1691 | gut microbiome measurement | http://www.ebi.ac.uk/efo/EFO_0007874 |
| 1510 | bone tissue density | http://purl.obolibrary.org/obo/OBA_1000110 |
| 1395 | protein measurement | http://www.ebi.ac.uk/efo/EFO_0004747 |
| 1368 | educational attainment | http://www.ebi.ac.uk/efo/EFO_0011015 |
| 1334 | body mass index | http://www.ebi.ac.uk/efo/EFO_0004340 |
| 1125 | COVID-19 | http://purl.obolibrary.org/obo/MONDO_0100096 |
| 1125 | breastfeeding duration | http://www.ebi.ac.uk/efo/EFO_0006864 |
| 1100 | systolic blood pressure | http://www.ebi.ac.uk/efo/EFO_0006335 |
| 1098 | smoking initiation | http://www.ebi.ac.uk/efo/EFO_0005670 |
| 1044 | free cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0008591 |
| 1028 | diastolic blood pressure | http://www.ebi.ac.uk/efo/EFO_0006336 |
| 1027 | body shape measurement | http://www.ebi.ac.uk/efo/EFO_0801227 |
| 1024 | allergen exposure measurement | http://www.ebi.ac.uk/efo/EFO_0007944 |
| 1013 | glioma pathogenesis-related protein 1 measurement | http://www.ebi.ac.uk/efo/EFO_0802567 |
| 1004 | body weight | http://www.ebi.ac.uk/efo/EFO_0004338 |
| 987 | type 2 diabetes mellitus | http://purl.obolibrary.org/obo/MONDO_0005148 |
| 964 | high density lipoprotein cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0004612 |
| 950 | blood protein amount | http://purl.obolibrary.org/obo/OBA_VT0005416 |
| 927 | self reported educational attainment | http://www.ebi.ac.uk/efo/EFO_0004784 |
| 908 | vaginal microbiome measurement | http://www.ebi.ac.uk/efo/EFO_0011013 |
| 887 | lipid measurement | http://www.ebi.ac.uk/efo/EFO_0004529 |
| 869 | triglyceride measurement | http://www.ebi.ac.uk/efo/EFO_0004530 |
| 868 | adolescent idiopathic scoliosis | http://purl.obolibrary.org/obo/MONDO_0005488 |
| 853 | diet measurement | http://www.ebi.ac.uk/efo/EFO_0008111 |
| 849 | FEV/FVC ratio | http://www.ebi.ac.uk/efo/EFO_0004713 |
| 814 | urate measurement | http://www.ebi.ac.uk/efo/EFO_0004531 |
| 810 | Stuttering | http://purl.obolibrary.org/obo/HP_0025268 |
| 810 | metabolic syndrome | http://www.ebi.ac.uk/efo/EFO_0000195 |
| 806 | glomerular filtration rate | http://purl.obolibrary.org/obo/OBA_0003747 |
| 804 | color vision disorder | http://purl.obolibrary.org/obo/MONDO_0001703 |
| 800 | serum creatinine amount | http://purl.obolibrary.org/obo/OBA_2050096 |
| 787 | forced expiratory volume | http://www.ebi.ac.uk/efo/EFO_0004314 |
| 786 | pulse pressure measurement | http://www.ebi.ac.uk/efo/EFO_0005763 |
| 780 | alcohol consumption characteristic | http://purl.obolibrary.org/obo/OBA_1000840 |
| 779 | environmental exposure measurement | http://www.ebi.ac.uk/efo/EFO_0008360 |
| 775 | heel bone mineral density | http://www.ebi.ac.uk/efo/EFO_0009270 |
| 711 | smoking status measurement | http://www.ebi.ac.uk/efo/EFO_0006527 |
| 705 | colorectal cancer | http://purl.obolibrary.org/obo/MONDO_0005575 |
| 692 | sex hormone-binding globulin measurement | http://www.ebi.ac.uk/efo/EFO_0004696 |
| 671 | health trait | http://www.ebi.ac.uk/efo/EFO_0007652 |
| 657 | low density lipoprotein cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0004611 |
| 652 | total cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0004574 |
| 640 | vital capacity | http://www.ebi.ac.uk/efo/EFO_0004312 |
| 628 | sex interaction measurement | http://www.ebi.ac.uk/efo/EFO_0008343 |
| 604 | appendicular lean mass | http://www.ebi.ac.uk/efo/EFO_0004980 |
| 601 | refractive error | http://purl.obolibrary.org/obo/MONDO_0004892 |
| 599 | anthropometric measurement | http://www.ebi.ac.uk/efo/EFO_0004302 |
| 598 | waist-hip ratio | http://www.ebi.ac.uk/efo/EFO_0004343 |
| 589 | testosterone measurement | http://www.ebi.ac.uk/efo/EFO_0004908 |
| 586 | obesity disorder | http://purl.obolibrary.org/obo/MONDO_0011122 |
| 570 | asthma | http://purl.obolibrary.org/obo/MONDO_0004979 |
| 558 | BMI-adjusted waist circumference | http://www.ebi.ac.uk/efo/EFO_0007789 |
| 530 | response to bronchodilator | http://purl.obolibrary.org/obo/GO_0097366 |
| 516 | BMI-adjusted waist-hip ratio | http://www.ebi.ac.uk/efo/EFO_0007788 |
| 514 | BMI-adjusted hip circumference | http://www.ebi.ac.uk/efo/EFO_0008039 |
| 514 | breast carcinoma | http://purl.obolibrary.org/obo/MONDO_0004989 |
| 507 | Crohn disease | http://purl.obolibrary.org/obo/MONDO_0005011 |
| 496 | oligodendroglioma | http://purl.obolibrary.org/obo/MONDO_0016695 |
| 489 | glioblastoma | http://purl.obolibrary.org/obo/MONDO_0018177 |
| 471 | age at onset | http://www.ebi.ac.uk/efo/EFO_0004847 |
| 467 | lean body mass | http://www.ebi.ac.uk/efo/EFO_0004995 |
| 465 | C-reactive protein measurement | http://www.ebi.ac.uk/efo/EFO_0004458 |
| 464 | serum albumin amount | http://purl.obolibrary.org/obo/OBA_2050068 |
| 464 | size | http://www.ebi.ac.uk/efo/EFO_0001695 |
| 459 | ulcerative colitis | http://purl.obolibrary.org/obo/MONDO_0005101 |
| 458 | amount of metabolite in blood serum | http://purl.obolibrary.org/obo/OBA_2050092 |
| 450 | substance abuse | http://purl.obolibrary.org/obo/MONDO_0002491 |
| 448 | Abnormality of refraction | http://purl.obolibrary.org/obo/HP_0000539 |
| 448 | fat pad mass | http://purl.obolibrary.org/obo/OBA_0004528 |
| 448 | socioeconomic status | http://www.ebi.ac.uk/efo/EFO_0022866 |
| 447 | social inhibition characteristic | http://purl.obolibrary.org/obo/OBA_1000278 |
| 446 | body fat percentage | http://www.ebi.ac.uk/efo/EFO_0007800 |
| 437 | breast density | http://purl.obolibrary.org/obo/OBA_2050351 |
| 434 | cystatin C measurement | http://www.ebi.ac.uk/efo/EFO_0004617 |
| 433 | chronic obstructive pulmonary disease | http://purl.obolibrary.org/obo/MONDO_0005002 |
| 431 | hypothyroidism | http://purl.obolibrary.org/obo/MONDO_0005420 |
| 429 | hypertensive disorder | http://purl.obolibrary.org/obo/MONDO_0005044 |
| 426 | prostate carcinoma | http://purl.obolibrary.org/obo/MONDO_0005159 |
| 422 | amino acid measurement | http://www.ebi.ac.uk/efo/EFO_0005134 |
| 420 | HbA1c measurement | http://www.ebi.ac.uk/efo/EFO_0004541 |
| 420 | chronotype measurement | http://www.ebi.ac.uk/efo/EFO_0008328 |
| 419 | fatty acid amount | http://purl.obolibrary.org/obo/OBA_VT0010046 |
| 415 | bilirubin measurement | http://www.ebi.ac.uk/efo/EFO_0004570 |
| 413 | blood glucose amount | http://purl.obolibrary.org/obo/OBA_VT0000188 |
| 407 | wellbeing measurement | http://www.ebi.ac.uk/efo/EFO_0007869 |
| 401 | inflammatory bowel disease | http://purl.obolibrary.org/obo/MONDO_0005265 |
| 398 | metabolite measurement | http://www.ebi.ac.uk/efo/EFO_0004725 |
| 397 | cleft lip | http://www.ebi.ac.uk/efo/EFO_0003959 |
| 395 | lifestyle measurement | http://www.ebi.ac.uk/efo/EFO_0010724 |
| 394 | age at assessment | http://www.ebi.ac.uk/efo/EFO_0008007 |
| 394 | phospholipids:total lipids ratio | http://www.ebi.ac.uk/efo/EFO_0020946 |
| 393 | psoriasis | http://purl.obolibrary.org/obo/MONDO_0005083 |
| 388 | total blood protein measurement | http://www.ebi.ac.uk/efo/EFO_0004536 |
| 380 | apolipoprotein A 1 measurement | http://www.ebi.ac.uk/efo/EFO_0004614 |
| 376 | whole body water mass | http://www.ebi.ac.uk/efo/EFO_0009805 |
| 373 | risk-taking behaviour | http://www.ebi.ac.uk/efo/EFO_0008579 |
| 368 | uric acid measurement | http://www.ebi.ac.uk/efo/EFO_0004761 |
| 367 | gout | http://purl.obolibrary.org/obo/MONDO_0005393 |
| 367 | reaction time measurement | http://www.ebi.ac.uk/efo/EFO_0008393 |
| 361 | mean arterial pressure | http://www.ebi.ac.uk/efo/EFO_0006340 |
| 359 | age at menarche | http://www.ebi.ac.uk/efo/EFO_0004703 |
| 359 | chromosome, telomeric region length | http://purl.obolibrary.org/obo/OBA_0003784 |
| 358 | balding measurement | http://www.ebi.ac.uk/efo/EFO_0007825 |
| 357 | Abnormality of the skeletal system | http://purl.obolibrary.org/obo/HP_0000924 |
| 357 | non-high density lipoprotein cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0005689 |
| 356 | DNA methylation | http://www.ebi.ac.uk/efo/EFO_0022599 |
| 355 | hemoglobin A1 measurement | http://www.ebi.ac.uk/efo/EFO_0007629 |
| 353 | amount of Thyroid stimulating hormone | http://purl.obolibrary.org/obo/OBA_1000983 |
| 351 | body surface area | http://www.ebi.ac.uk/efo/EFO_0022196 |
| 350 | tobacco smoke exposure measurement | http://www.ebi.ac.uk/efo/EFO_0009115 |
| 349 | smoking behavior | http://www.ebi.ac.uk/efo/EFO_0004318 |
| 348 | frailty measurement | http://www.ebi.ac.uk/efo/EFO_0009885 |
| 348 | waist circumference | http://purl.obolibrary.org/obo/OBA_1001085 |
| 346 | Myopia | http://purl.obolibrary.org/obo/HP_0000545 |
| 346 | venous thromboembolism | http://purl.obolibrary.org/obo/MONDO_0005399 |
| 345 | IGF-1 measurement | http://www.ebi.ac.uk/efo/EFO_0004627 |
| 345 | severe acute respiratory syndrome | http://purl.obolibrary.org/obo/MONDO_0005091 |
| 343 | birth weight | http://www.ebi.ac.uk/efo/EFO_0004344 |
| 341 | grip strength measurement | http://www.ebi.ac.uk/efo/EFO_0006941 |
| 338 | carotid artery thickness | http://purl.obolibrary.org/obo/OBA_2050108 |
| 338 | systemic lupus erythematosus | http://purl.obolibrary.org/obo/MONDO_0007915 |
| 335 | hip circumference | http://purl.obolibrary.org/obo/OBA_1000032 |
| 334 | autoimmune thyroid disease | http://purl.obolibrary.org/obo/MONDO_0005623 |
| 334 | facial morphology trait | http://purl.obolibrary.org/obo/OBA_VT0003743 |
| 331 | glucose measurement | http://www.ebi.ac.uk/efo/EFO_0004468 |
| 330 | calcium measurement | http://www.ebi.ac.uk/efo/EFO_0004838 |
| 330 | prostate specific antigen amount | http://purl.obolibrary.org/obo/OBA_2050200 |
| 329 | cholesteryl ester measurement | http://www.ebi.ac.uk/efo/EFO_0010351 |
| 327 | electrocardiography | http://www.ebi.ac.uk/efo/EFO_0004327 |
| 326 | prostate cancer | http://purl.obolibrary.org/obo/MONDO_0008315 |
| 314 | drug use measurement | http://www.ebi.ac.uk/efo/EFO_0007010 |
| 311 | cytokine measurement | http://www.ebi.ac.uk/efo/EFO_0004873 |
| 308 | visceral adipose tissue quantity | http://purl.obolibrary.org/obo/OBA_2045224 |
| 299 | serum IgG glycosylation measurement | http://www.ebi.ac.uk/efo/EFO_0005193 |
| 297 | physical activity measurement | http://www.ebi.ac.uk/efo/EFO_0008002 |
| 295 | strand of hair color | http://purl.obolibrary.org/obo/OBA_0002313 |
| 292 | sexual dimorphism measurement | http://www.ebi.ac.uk/efo/EFO_0021796 |
| 289 | cardiovascular disorder | http://purl.obolibrary.org/obo/MONDO_0004995 |
| 289 | retinal vasculature measurement | http://www.ebi.ac.uk/efo/EFO_0010554 |
| 288 | disease progression measurement | http://www.ebi.ac.uk/efo/EFO_0008336 |
| 287 | metabolic dysfunction-associated steatotic liver disease | http://purl.obolibrary.org/obo/MONDO_0013209 |
| 286 | breast cancer | http://purl.obolibrary.org/obo/MONDO_0007254 |
| 282 | mortality | http://www.ebi.ac.uk/efo/EFO_0004352 |
| 282 | total lipids in large HDL | http://www.ebi.ac.uk/efo/EFO_0022189 |
| 281 | Cleft palate | http://purl.obolibrary.org/obo/HP_0000175 |
| 281 | blood urea nitrogen amount | http://purl.obolibrary.org/obo/OBA_VT0005265 |
| 280 | executive function measurement | http://www.ebi.ac.uk/efo/EFO_0009332 |
| 280 | rheumatoid arthritis | http://purl.obolibrary.org/obo/MONDO_0008383 |
| 279 | open-angle glaucoma | http://purl.obolibrary.org/obo/MONDO_0005338 |
| 278 | polyunsaturated fatty acids to total fatty acids percentage | http://www.ebi.ac.uk/efo/EFO_0022303 |
| 275 | alcohol drinking | http://www.ebi.ac.uk/efo/EFO_0004329 |
| 273 | treatment-resistant hypertension | http://www.ebi.ac.uk/efo/EFO_1002006 |
| 270 | astrocytoma (excluding glioblastoma) | http://purl.obolibrary.org/obo/MONDO_0019781 |
| 269 | ovarian carcinoma | http://purl.obolibrary.org/obo/MONDO_0005140 |
| 266 | peak expiratory flow | http://www.ebi.ac.uk/efo/EFO_0009718 |
| 264 | comparative body size at age 10, self-reported | http://www.ebi.ac.uk/efo/EFO_0009819 |
| 264 | delivery measurement | http://www.ebi.ac.uk/efo/EFO_0006922 |
| 263 | osteoarthritis | http://purl.obolibrary.org/obo/MONDO_0005178 |
| 262 | smoking behavior trait | http://purl.obolibrary.org/obo/OBA_2050116 |
| 261 | endometriosis | http://purl.obolibrary.org/obo/MONDO_0005133 |
| 259 | osteoarthritis, knee | http://purl.obolibrary.org/obo/MONDO_0005416 |
| 255 | blood VLDL cholesterol amount | http://purl.obolibrary.org/obo/OBA_VT0005144 |
| 255 | intraocular pressure measurement | http://www.ebi.ac.uk/efo/EFO_0004695 |
| 249 | cholesteryl esters in large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022159 |
| 247 | omega-6 polyunsaturated fatty acid measurement | http://www.ebi.ac.uk/efo/EFO_0005680 |
| 246 | HDL cholesterol change measurement | http://www.ebi.ac.uk/efo/EFO_0007805 |
| 245 | apolipoprotein B measurement | http://www.ebi.ac.uk/efo/EFO_0004615 |
| 245 | parental longevity | http://www.ebi.ac.uk/efo/EFO_0007796 |
| 244 | age at first sexual intercourse measurement | http://www.ebi.ac.uk/efo/EFO_0009749 |
| 244 | polygenic risk score | http://www.ebi.ac.uk/efo/EFO_0030082 |
| 243 | abdominal adipose tissue measurement | http://www.ebi.ac.uk/efo/EFO_0803316 |
| 243 | parental genotype effect measurement | http://www.ebi.ac.uk/efo/EFO_0005939 |
| 243 | polyunsaturated fatty acids to monounsaturated fatty acids ratio | http://www.ebi.ac.uk/efo/EFO_0022302 |
| 239 | cholesterol in large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0021900 |
| 239 | triglyceride:HDL cholesterol ratio | http://www.ebi.ac.uk/efo/EFO_0007929 |
| 237 | mosquito bite reaction size measurement | http://www.ebi.ac.uk/efo/EFO_0008378 |
| 236 | free androgen index | http://www.ebi.ac.uk/efo/EFO_0007005 |
| 235 | aging | http://www.ebi.ac.uk/efo/EFO_0022597 |
| 235 | base metabolic rate measurement | http://www.ebi.ac.uk/efo/EFO_0007777 |
| 235 | cholesteryl esters in HDL measurement | http://www.ebi.ac.uk/efo/EFO_0920022 |
| 234 | type 1 diabetes mellitus | http://purl.obolibrary.org/obo/MONDO_0005147 |
| 233 | cholesterol to total lipids in large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022234 |
| 232 | gastroesophageal reflux disease | http://purl.obolibrary.org/obo/MONDO_0007186 |
| 231 | amount of phospholipid | http://purl.obolibrary.org/obo/OBA_1001011 |
| 231 | linoleic acid measurement | http://www.ebi.ac.uk/efo/EFO_0006807 |
| 231 | low density lipoprotein triglyceride measurement | http://www.ebi.ac.uk/efo/EFO_0009946 |
| 230 | amount of vitamin D | http://purl.obolibrary.org/obo/OBA_1000968 |
| 230 | facial morphology | http://www.ebi.ac.uk/efo/EFO_0004743 |
| 228 | HDL particle size | http://purl.obolibrary.org/obo/OBA_2050115 |
| 227 | retinal layer thickness | http://purl.obolibrary.org/obo/OBA_2050110 |
| 226 | anorexia nervosa | http://purl.obolibrary.org/obo/MONDO_0005351 |
| 225 | lumbar disc herniation | http://www.ebi.ac.uk/efo/EFO_1002005 |
| 224 | concentration of large HDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022188 |
| 223 | squamous cell lung carcinoma | http://purl.obolibrary.org/obo/MONDO_0005097 |
| 222 | osteoarthritis, hip | http://purl.obolibrary.org/obo/MONDO_0006629 |
| 221 | life span determination trait | http://purl.obolibrary.org/obo/OBA_VT0005372 |
| 220 | cholesteryl esters to total lipids in large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022248 |
| 220 | pyruvate measurement | http://www.ebi.ac.uk/efo/EFO_0010117 |
| 219 | amount of phosphatidylcholine | http://purl.obolibrary.org/obo/OBA_2040179 |
| 219 | chronic kidney disease | http://purl.obolibrary.org/obo/MONDO_0005300 |
| 219 | free cholesterol to total lipids in large LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022280 |
| 219 | response to antineoplastic agent | http://purl.obolibrary.org/obo/GO_0097327 |
| 219 | substance-related disorder | http://purl.obolibrary.org/obo/MONDO_0002494 |
| 219 | triglycerides to phosphoglycerides ratio | http://www.ebi.ac.uk/efo/EFO_0022327 |
| 218 | free cholesterol in large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022157 |
| 218 | lung carcinoma | http://purl.obolibrary.org/obo/MONDO_0005138 |
| 218 | omega-3 polyunsaturated fatty acid measurement | http://www.ebi.ac.uk/efo/EFO_0010119 |
| 217 | Tinnitus | http://purl.obolibrary.org/obo/HP_0000360 |
| 216 | amount of cystatin-C (human) in blood | http://purl.obolibrary.org/obo/OBA_2052375 |
| 216 | corneal resistance factor | http://www.ebi.ac.uk/efo/EFO_0010067 |
| 216 | phospholipids in very large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022298 |
| 215 | total lipids in very large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022312 |
| 214 | cholesterol in very large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022229 |
| 213 | saturated fatty acids to total fatty acids percentage | http://www.ebi.ac.uk/efo/EFO_0022305 |
| 212 | ankylosing spondylitis | http://purl.obolibrary.org/obo/MONDO_0005306 |
| 212 | triglycerides to total lipids in very large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022339 |
| 211 | cigarettes per day measurement | http://www.ebi.ac.uk/efo/EFO_0006525 |
| 211 | glioma | http://purl.obolibrary.org/obo/MONDO_0021042 |
| 210 | free cholesterol to total lipids in medium LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022283 |
| 208 | lipoprotein measurement | http://www.ebi.ac.uk/efo/EFO_0004732 |
| 207 | diverticular disease | http://www.ebi.ac.uk/efo/EFO_0009959 |
| 207 | triglycerides to total lipids in small HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022336 |
| 206 | free cholesterol:total lipids ratio | http://www.ebi.ac.uk/efo/EFO_0020945 |
| 206 | response to xenobiotic stimulus | http://purl.obolibrary.org/obo/GO_0009410 |
| 205 | triglycerides in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022158 |
| 203 | cholesterol to total lipids in medium HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022237 |
| 202 | cholesterol:total lipids ratio | http://www.ebi.ac.uk/efo/EFO_0020943 |
| 201 | phospholipids in large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022190 |
| 200 | 3-hydroxypropylmercapturic acid measurement | http://www.ebi.ac.uk/efo/EFO_0007014 |
| 200 | alcohol dependence | http://purl.obolibrary.org/obo/MONDO_0007079 |
| 199 | 3-hydroxy-1-methylpropylmercapturic acid measurement | http://www.ebi.ac.uk/efo/EFO_0007015 |
| 199 | dental caries | http://purl.obolibrary.org/obo/MONDO_0005276 |
| 199 | glycoprotein measurement | http://www.ebi.ac.uk/efo/EFO_0004555 |
| 199 | nose morphology trait | http://purl.obolibrary.org/obo/OBA_VT0002233 |
| 198 | atopic eczema | http://purl.obolibrary.org/obo/MONDO_0004980 |
| 198 | cup-to-disc ratio measurement | http://www.ebi.ac.uk/efo/EFO_0006939 |
| 198 | restless legs syndrome | http://purl.obolibrary.org/obo/MONDO_0005391 |
| 197 | diabetes mellitus | http://purl.obolibrary.org/obo/MONDO_0005015 |
| 196 | triglycerides to total lipids in large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022330 |
| 195 | blood immunoglobulin amount | http://purl.obolibrary.org/obo/OBA_VT0002460 |
| 195 | pain measurement | http://www.ebi.ac.uk/efo/EFO_0010639 |
| 194 | S-6-hydroxywarfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803326 |
| 194 | nicotine dependence symptom count | http://www.ebi.ac.uk/efo/EFO_0009262 |
| 193 | cholesteryl esters:total lipids ratio | http://www.ebi.ac.uk/efo/EFO_0020944 |
| 192 | cardiovascular measurement | http://www.ebi.ac.uk/efo/EFO_0004298 |
| 192 | glutamine measurement | http://www.ebi.ac.uk/efo/EFO_0009768 |
| 191 | age at menopause | http://www.ebi.ac.uk/efo/EFO_0004704 |
| 190 | Agents acting on the renin-angiotensin system use measurement | http://www.ebi.ac.uk/efo/EFO_0009931 |
| 190 | lung adenocarcinoma | http://purl.obolibrary.org/obo/MONDO_0005061 |
| 190 | post-traumatic stress disorder | http://purl.obolibrary.org/obo/MONDO_0005146 |
| 189 | Antihypertensive use measurement | http://www.ebi.ac.uk/efo/EFO_0009927 |
| 189 | triglycerides in IDL measurement | http://www.ebi.ac.uk/efo/EFO_0022149 |
| 188 | free cholesterol in very large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022273 |
| 188 | polyunsaturated fatty acid measurement | http://www.ebi.ac.uk/efo/EFO_0010733 |
| 188 | response to antipsychotic drug | http://purl.obolibrary.org/obo/GO_0097332 |
| 187 | total lipids in large VLDL | http://www.ebi.ac.uk/efo/EFO_0022175 |
| 187 | triglycerides to total lipids in medium HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022333 |
| 186 | total lipids in small VLDL | http://www.ebi.ac.uk/efo/EFO_0022148 |
| 185 | VLDL particle size | http://purl.obolibrary.org/obo/OBA_2050114 |
| 185 | cholesteryl esters to total lipids in very small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022259 |
| 185 | free cholesterol to total lipids in medium HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022282 |
| 185 | hip geometry | http://www.ebi.ac.uk/efo/EFO_0004685 |
| 185 | response to antidepressant | http://purl.obolibrary.org/obo/GO_0036276 |
| 184 | triglycerides to total lipids in small LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022337 |
| 184 | uterine corpus leiomyoma | http://purl.obolibrary.org/obo/MONDO_0007886 |
| 183 | concentration of very large HDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0920031 |
| 183 | height growth attribute | http://purl.obolibrary.org/obo/OBA_2045231 |
| 183 | mosquito bite reaction itch intensity measurement | http://www.ebi.ac.uk/efo/EFO_0008377 |
| 183 | soluble triggering receptor expressed on myeloid cells 2 measurement | http://www.ebi.ac.uk/efo/EFO_0010151 |
| 181 | cerebrospinal fluid composition attribute | http://purl.obolibrary.org/obo/OBA_2045291 |
| 180 | cholesterol in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0021903 |
| 179 | cancer | http://purl.obolibrary.org/obo/MONDO_0004992 |
| 179 | free cholesterol in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022267 |
| 179 | triglycerides in small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022145 |
| 178 | Eczematoid dermatitis | http://purl.obolibrary.org/obo/HP_0000964 |
| 178 | daytime rest measurement | http://www.ebi.ac.uk/efo/EFO_0007828 |
| 178 | smoking cessation | http://www.ebi.ac.uk/efo/EFO_0004319 |
| 177 | free cholesterol to total lipids in very large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022288 |
| 177 | gestational age | http://www.ebi.ac.uk/efo/EFO_0005112 |
| 177 | hospitalisation | http://www.ebi.ac.uk/efo/EFO_0020922 |
| 176 | cholesterol to total lipids in medium VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022239 |
| 176 | cholesterol to total lipids in very small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022245 |
| 176 | cholesteryl esters to total lipids in medium VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022253 |
| 175 | alanine measurement | http://www.ebi.ac.uk/efo/EFO_0009765 |
| 175 | free cholesterol in HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022264 |
| 175 | language measurement | http://www.ebi.ac.uk/efo/EFO_0007797 |
| 175 | triglycerides to total lipids in IDL percentage | http://www.ebi.ac.uk/efo/EFO_0022329 |
| 174 | insulin measurement | http://www.ebi.ac.uk/efo/EFO_0004467 |
| 174 | peripheral arterial disease | http://purl.obolibrary.org/obo/MONDO_0005386 |
| 174 | response to radiation | http://purl.obolibrary.org/obo/GO_0009314 |
| 173 | cholesteryl esters in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0920026 |
| 173 | lung cancer | http://purl.obolibrary.org/obo/MONDO_0008903 |
| 173 | total lipids in HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022307 |
| 172 | respiratory system disorder | http://purl.obolibrary.org/obo/MONDO_0005087 |
| 172 | response to selective serotonin reuptake inhibitor | http://www.ebi.ac.uk/efo/EFO_0005658 |
| 170 | concentration of medium HDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0920032 |
| 170 | sclerosing cholangitis | http://purl.obolibrary.org/obo/MONDO_0018646 |
| 170 | sphingomyelin measurement | http://www.ebi.ac.uk/efo/EFO_0010118 |
| 169 | free cholesterol to total lipids in small LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022286 |
| 169 | intermediate density lipoprotein measurement | http://www.ebi.ac.uk/efo/EFO_0008595 |
| 168 | cholesteryl esters to total lipids in medium HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022251 |
| 168 | glycine measurement | http://www.ebi.ac.uk/efo/EFO_0009767 |
| 167 | triglycerides in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022325 |
| 167 | triglycerides in very small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022144 |
| 166 | Drugs used in diabetes use measurement | http://www.ebi.ac.uk/efo/EFO_0009924 |
| 166 | hair color | http://www.ebi.ac.uk/efo/EFO_0003924 |
| 166 | hyperopia | http://purl.obolibrary.org/obo/MONDO_0004891 |
| 166 | irritable bowel syndrome | http://purl.obolibrary.org/obo/MONDO_0005052 |
| 166 | triglycerides:total lipids ratio | http://www.ebi.ac.uk/efo/EFO_0020947 |
| 165 | central corneal thickness | http://www.ebi.ac.uk/efo/EFO_0005213 |
| 165 | exocrine pancreatic carcinoma | http://purl.obolibrary.org/obo/MONDO_0005192 |
| 165 | household income | http://www.ebi.ac.uk/efo/EFO_0009695 |
| 164 | total lipids in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022310 |
| 163 | docosahexaenoic acid measurement | http://www.ebi.ac.uk/efo/EFO_0007761 |
| 163 | phospholipids in HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022293 |
| 163 | social risk factor | http://www.ebi.ac.uk/efo/EFO_0022923 |
| 163 | triglycerides to total lipids in very small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022341 |
| 162 | mitochondrial DNA measurement | http://www.ebi.ac.uk/efo/EFO_0006312 |
| 162 | orofacial cleft | http://purl.obolibrary.org/obo/MONDO_0000358 |
| 161 | coffee consumption measurement | http://www.ebi.ac.uk/efo/EFO_0006781 |
| 161 | free cholesterol to total lipids in small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022287 |
| 160 | body composition measurement | http://www.ebi.ac.uk/efo/EFO_0005106 |
| 160 | body fat distribution | http://www.ebi.ac.uk/efo/EFO_0004341 |
| 160 | income | http://www.ebi.ac.uk/efo/EFO_0022958 |
| 159 | triglycerides to total lipids in medium VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022335 |
| 158 | esterified cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0008589 |
| 158 | eye colour measurement | http://www.ebi.ac.uk/efo/EFO_0009764 |
| 158 | physical activity | http://www.ebi.ac.uk/efo/EFO_0003940 |
| 158 | triglycerides in small LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022323 |
| 157 | cholesterol to total lipids in IDL percentage | http://www.ebi.ac.uk/efo/EFO_0022233 |
| 156 | R-warfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803324 |
| 156 | allergic rhinitis | http://purl.obolibrary.org/obo/MONDO_0011786 |
| 156 | eye measurement | http://www.ebi.ac.uk/efo/EFO_0004731 |
| 155 | mental disorder | http://purl.obolibrary.org/obo/MONDO_0005084 |
| 155 | saturated fatty acids measurement | http://www.ebi.ac.uk/efo/EFO_0022304 |
| 154 | alcohol use disorder measurement | http://www.ebi.ac.uk/efo/EFO_0009458 |
| 154 | cholesteryl esters to total lipids in very large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022258 |
| 154 | phospholipids in chylomicrons and extremely large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022292 |
| 154 | response to vaccine | http://www.ebi.ac.uk/efo/EFO_0004645 |
| 154 | triglycerides in medium LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022322 |
| 154 | triglycerides to total lipids in large LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022331 |
| 153 | amount of serum globulin type protein | http://purl.obolibrary.org/obo/OBA_2045206 |
| 153 | phospholipids in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022295 |
| 152 | circulating fibrinogen amount | http://purl.obolibrary.org/obo/OBA_0000061 |
| 152 | concentration of chylomicrons and extremely large VLDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022260 |
| 152 | interleukin-6 measurement | http://www.ebi.ac.uk/efo/EFO_0004810 |
| 152 | response to diisocyanate | http://www.ebi.ac.uk/efo/EFO_0006995 |
| 151 | overall survival | http://www.ebi.ac.uk/efo/EFO_0000638 |
| 150 | childhood onset asthma | http://purl.obolibrary.org/obo/MONDO_0005405 |
| 150 | free cholesterol in chylomicrons and extremely large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022263 |
| 150 | mood disorder | http://purl.obolibrary.org/obo/MONDO_0005371 |
| 150 | total lipids in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022311 |
| 149 | aortic measurement | http://www.ebi.ac.uk/efo/EFO_0020865 |
| 149 | cholesterol to total lipids in very large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022244 |
| 148 | periodontitis | http://purl.obolibrary.org/obo/MONDO_0005076 |
| 147 | basal cell carcinoma | http://purl.obolibrary.org/obo/MONDO_0020804 |
| 147 | biological sex | http://purl.obolibrary.org/obo/PATO_0000047 |
| 147 | triglycerides in VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022326 |
| 146 | choline measurement | http://www.ebi.ac.uk/efo/EFO_0010116 |
| 146 | degree of unsaturation measurement | http://www.ebi.ac.uk/efo/EFO_0022261 |
| 146 | growth rate attribute | http://purl.obolibrary.org/obo/OBA_2045230 |
| 146 | opiate dependence | http://purl.obolibrary.org/obo/MONDO_0005530 |
| 146 | pelvic organ prolapse | http://purl.obolibrary.org/obo/MONDO_0000082 |
| 146 | taste liking measurement | http://www.ebi.ac.uk/efo/EFO_0010155 |
| 144 | total lipids in medium VLDL | http://www.ebi.ac.uk/efo/EFO_0022153 |
| 144 | total lipids in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022313 |
| 144 | traffic air pollution measurement | http://www.ebi.ac.uk/efo/EFO_0007908 |
| 143 | triglycerides to total lipids in small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022338 |
| 142 | concentration of very large VLDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0920030 |
| 142 | femoral neck bone mineral density | http://www.ebi.ac.uk/efo/EFO_0007785 |
| 142 | glaucoma | http://purl.obolibrary.org/obo/MONDO_0005041 |
| 142 | inflammatory biomarker measurement | http://www.ebi.ac.uk/efo/EFO_0004872 |
| 142 | phospholipids in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022296 |
| 142 | triglycerides in HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022317 |
| 141 | total lipids in chylomicrons and extremely large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022306 |
| 141 | triglycerides in large LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022319 |
| 141 | triglycerides to total lipids in medium LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022334 |
| 140 | allergic disease | http://purl.obolibrary.org/obo/MONDO_0005271 |
| 140 | blood insulin amount | http://purl.obolibrary.org/obo/OBA_VT0001560 |
| 140 | sodium measurement | http://www.ebi.ac.uk/efo/EFO_0009282 |
| 140 | triglycerides in medium HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022321 |
| 139 | Hashimoto thyroiditis | http://purl.obolibrary.org/obo/MONDO_0007699 |
| 139 | abdominal aortic aneurysm | http://purl.obolibrary.org/obo/MONDO_0005350 |
| 139 | cholesterol in chylomicrons and extremely large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0021898 |
| 139 | pilocytic astrocytoma | http://purl.obolibrary.org/obo/MONDO_0016691 |
| 138 | attempted suicide | http://www.ebi.ac.uk/efo/EFO_0004321 |
| 138 | cholesteryl esters in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0920027 |
| 138 | diabetic kidney disease | http://purl.obolibrary.org/obo/MONDO_0005016 |
| 138 | total lipids in large LDL | http://www.ebi.ac.uk/efo/EFO_0022163 |
| 138 | tyrosine measurement | http://www.ebi.ac.uk/efo/EFO_0005058 |
| 137 | APOE carrier status | http://www.ebi.ac.uk/efo/EFO_0007659 |
| 137 | diabetes mellitus biomarker | http://www.ebi.ac.uk/efo/EFO_0006842 |
| 137 | trait in response to paclitaxel | http://purl.obolibrary.org/obo/OBA_2040073 |
| 137 | triglycerides in LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022320 |
| 136 | cholesteryl esters in VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0920021 |
| 136 | free cholesterol to total lipids in medium VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022284 |
| 136 | macula attribute | http://purl.obolibrary.org/obo/OBA_2045216 |
| 135 | cholesterol to total lipids in small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022242 |
| 135 | citrate measurement | http://www.ebi.ac.uk/efo/EFO_0010114 |
| 135 | phospholipids in large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022169 |
| 135 | primary angle-closure glaucoma | http://purl.obolibrary.org/obo/MONDO_0001868 |
| 135 | radiation-induced disorder | http://purl.obolibrary.org/obo/MONDO_0043459 |
| 135 | urinary metabolite measurement | http://www.ebi.ac.uk/efo/EFO_0005116 |
| 134 | cholesterol to total lipids in small HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022240 |
| 134 | free cholesterol in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022270 |
| 134 | information processing speed | http://www.ebi.ac.uk/efo/EFO_0004363 |
| 134 | presbycusis | http://purl.obolibrary.org/obo/MONDO_0043765 |
| 133 | cholesterol in small HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022226 |
| 133 | factor VII measurement | http://www.ebi.ac.uk/efo/EFO_0004619 |
| 133 | pulmonary function measurement | http://www.ebi.ac.uk/efo/EFO_0003892 |
| 133 | triglycerides in large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022178 |
| 133 | triglycerides in medium VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022155 |
| 132 | total lipids in small LDL | http://www.ebi.ac.uk/efo/EFO_0022168 |
| 132 | triglycerides in chylomicrons and extremely large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022316 |
| 131 | Pain | http://purl.obolibrary.org/obo/HP_0012531 |
| 131 | amyloid-beta measurement | http://www.ebi.ac.uk/efo/EFO_0005194 |
| 131 | cholesterol to total lipids in large LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022235 |
| 131 | trait in response to paliperidone | http://purl.obolibrary.org/obo/OBA_2040132 |
| 130 | IgA glomerulonephritis | http://purl.obolibrary.org/obo/MONDO_0005342 |
| 130 | concentration of large VLDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022173 |
| 130 | free cholesterol to total lipids in small HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022285 |
| 130 | serum urea amount | http://purl.obolibrary.org/obo/OBA_2050065 |
| 129 | Thyroid preparation use measurement | http://www.ebi.ac.uk/efo/EFO_0009933 |
| 129 | age-related macular degeneration | http://purl.obolibrary.org/obo/MONDO_0005150 |
| 129 | free cholesterol in large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022265 |
| 128 | concentration of small HDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0920033 |
| 128 | cutaneous melanoma | http://purl.obolibrary.org/obo/MONDO_0005012 |
| 128 | free cholesterol in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022274 |
| 128 | potassium measurement | http://www.ebi.ac.uk/efo/EFO_0009283 |
| 128 | reading and spelling ability | http://www.ebi.ac.uk/efo/EFO_0005301 |
| 128 | total lipids in IDL | http://www.ebi.ac.uk/efo/EFO_0022161 |
| 127 | cholesteryl esters to total lipids in large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022250 |
| 127 | free cholesterol to total lipids in large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022279 |
| 127 | plasminogen activator inhibitor 1 measurement | http://www.ebi.ac.uk/efo/EFO_0004792 |
| 127 | progression free survival | http://www.ebi.ac.uk/efo/EFO_0004920 |
| 127 | rheumatoid factor seropositivity measurement | http://www.ebi.ac.uk/efo/EFO_0007791 |
| 127 | spine bone mineral density | http://www.ebi.ac.uk/efo/EFO_0007701 |
| 126 | autoimmune disease | http://purl.obolibrary.org/obo/MONDO_0007179 |
| 126 | cholesteryl esters to total lipids in small HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022254 |
| 126 | optic disc size trait | http://purl.obolibrary.org/obo/OBA_VT0006216 |
| 126 | phosphoglycerides measurement | http://www.ebi.ac.uk/efo/EFO_0022291 |
| 126 | phospholipids in VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022301 |
| 125 | Calcium channel blocker use measurement | http://www.ebi.ac.uk/efo/EFO_0009930 |
| 125 | amount of ceramide | http://purl.obolibrary.org/obo/OBA_2040177 |
| 124 | response to stimulus | http://purl.obolibrary.org/obo/GO_0050896 |
| 123 | Graves disease | http://purl.obolibrary.org/obo/MONDO_0005364 |
| 123 | anti-citrullinated protein antibody seropositivity | http://www.ebi.ac.uk/efo/EFO_0007837 |
| 123 | colorectal health | http://www.ebi.ac.uk/efo/EFO_0008460 |
| 123 | nicotine dependence | http://purl.obolibrary.org/obo/MONDO_0008575 |
| 123 | phospholipids in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022299 |
| 123 | response to angiotensin-converting enzyme inhibitor | http://www.ebi.ac.uk/efo/EFO_0005325 |
| 123 | urinary bladder cancer | http://purl.obolibrary.org/obo/MONDO_0001187 |
| 121 | cholesteryl esters to total lipids in IDL percentage | http://www.ebi.ac.uk/efo/EFO_0022247 |
| 121 | cholesteryl esters to total lipids in small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022256 |
| 121 | sexual activity behaviour attribute | http://purl.obolibrary.org/obo/OBA_2045293 |
| 120 | diabetic retinopathy | http://purl.obolibrary.org/obo/MONDO_0005266 |
| 120 | eye morphology trait | http://purl.obolibrary.org/obo/OBA_VT0002092 |
| 119 | family history of lung cancer | http://www.ebi.ac.uk/efo/EFO_0006953 |
| 119 | total hip arthroplasty | http://www.ebi.ac.uk/efo/EFO_0009806 |
| 118 | tuberculosis | http://purl.obolibrary.org/obo/MONDO_0018076 |
| 118 | valine measurement | http://www.ebi.ac.uk/efo/EFO_0009792 |
| 117 | amount of thyroxine | http://purl.obolibrary.org/obo/OBA_1001010 |
| 117 | attention function measurement | http://www.ebi.ac.uk/efo/EFO_0007636 |
| 117 | cataract | http://purl.obolibrary.org/obo/MONDO_0005129 |
| 117 | factor XI measurement | http://www.ebi.ac.uk/efo/EFO_0004694 |
| 117 | tonsillectomy risk measurement | http://www.ebi.ac.uk/efo/EFO_0007924 |
| 116 | celiac disease | http://purl.obolibrary.org/obo/MONDO_0005130 |
| 115 | dyslexia | http://purl.obolibrary.org/obo/MONDO_0005489 |
| 115 | p-tau measurement | http://www.ebi.ac.uk/efo/EFO_0004763 |
| 115 | sarcoidosis | http://purl.obolibrary.org/obo/MONDO_0019338 |
| 115 | survival time | http://www.ebi.ac.uk/efo/EFO_0000714 |
| 114 | aging rate | http://purl.obolibrary.org/obo/OBA_0005494 |
| 114 | free cholesterol to total lipids in IDL percentage | http://www.ebi.ac.uk/efo/EFO_0022278 |
| 114 | gestational diabetes | http://purl.obolibrary.org/obo/MONDO_0005406 |
| 114 | hip bone mineral density | http://www.ebi.ac.uk/efo/EFO_0007702 |
| 114 | puberty onset measurement | http://www.ebi.ac.uk/efo/EFO_0005677 |
| 113 | lip morphology trait | http://purl.obolibrary.org/obo/OBA_VT0003769 |
| 113 | triglycerides in very large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022324 |
| 113 | triglycerides to total lipids in very large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022340 |
| 113 | visual perception quality | http://purl.obolibrary.org/obo/OBA_2045263 |
| 112 | Varicose veins | http://purl.obolibrary.org/obo/HP_0002619 |
| 112 | health study participation | http://www.ebi.ac.uk/efo/EFO_0010130 |
| 112 | upper aerodigestive tract neoplasm | http://purl.obolibrary.org/obo/MONDO_0005398 |
| 111 | cholesteryl esters to total lipids in very large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022257 |
| 111 | concentration of small VLDL particles | http://www.ebi.ac.uk/efo/EFO_0022147 |
| 111 | hormone replacement therapy | http://www.ebi.ac.uk/efo/EFO_0003961 |
| 111 | number of siblings | http://www.ebi.ac.uk/efo/EFO_0009713 |
| 111 | von Willebrand factor characteristic | http://purl.obolibrary.org/obo/OBA_1000199 |
| 110 | aggressive behavior characteristic | http://purl.obolibrary.org/obo/OBA_1000376 |
| 110 | cholesteryl esters in LDL measurement | http://www.ebi.ac.uk/efo/EFO_0920023 |
| 110 | interleukin 10 measurement | http://www.ebi.ac.uk/efo/EFO_0004750 |
| 110 | neuropathic pain | http://www.ebi.ac.uk/efo/EFO_0005762 |
| 109 | 25-hydroxyvitamin D3 measurement | http://www.ebi.ac.uk/efo/EFO_0803357 |
| 109 | leucine measurement | http://www.ebi.ac.uk/efo/EFO_0009770 |
| 109 | trauma exposure measurement | http://www.ebi.ac.uk/efo/EFO_0010703 |
| 108 | X-12063 measurement | http://www.ebi.ac.uk/efo/EFO_0021283 |
| 108 | adiponectin measurement | http://www.ebi.ac.uk/efo/EFO_0004502 |
| 107 | handedness | http://www.ebi.ac.uk/efo/EFO_0009902 |
| 107 | hearing loss disorder | http://purl.obolibrary.org/obo/MONDO_0005365 |
| 107 | renal cell carcinoma | http://purl.obolibrary.org/obo/MONDO_0005086 |
| 106 | Vertigo | http://purl.obolibrary.org/obo/HP_0002321 |
| 106 | lobe attachment | http://www.ebi.ac.uk/efo/EFO_0007667 |
| 105 | Diuretic use measurement | http://www.ebi.ac.uk/efo/EFO_0009928 |
| 105 | alopecia | http://purl.obolibrary.org/obo/MONDO_0004907 |
| 105 | chylomicron amount | http://purl.obolibrary.org/obo/OBA_2050296 |
| 105 | free cholesterol to total lipids in very large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022289 |
| 105 | pneumonia | http://purl.obolibrary.org/obo/MONDO_0005249 |
| 105 | response to fenofibrate | http://purl.obolibrary.org/obo/GO_1901557 |
| 105 | trait in response to apixaban | http://purl.obolibrary.org/obo/OBA_2050328 |
| 105 | triglycerides in large HDL measurement | http://www.ebi.ac.uk/efo/EFO_0022318 |
| 104 | creatine kinase measurement | http://www.ebi.ac.uk/efo/EFO_0004534 |
| 104 | response to COVID-19 vaccine | http://www.ebi.ac.uk/efo/EFO_0803362 |
| 103 | Headache | http://purl.obolibrary.org/obo/HP_0002315 |
| 103 | cholesterol in large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0021902 |
| 103 | chronic rhinosinusitis | http://purl.obolibrary.org/obo/MONDO_0006031 |
| 103 | factor VIII measurement | http://www.ebi.ac.uk/efo/EFO_0004630 |
| 103 | juvenile idiopathic arthritis | http://purl.obolibrary.org/obo/MONDO_0011429 |
| 103 | total phospholipids in lipoprotein particles measurement | http://www.ebi.ac.uk/efo/EFO_0022315 |
| 103 | trait in response to escitalopram | http://purl.obolibrary.org/obo/OBA_2040055 |
| 102 | asthma exacerbation measurement | http://www.ebi.ac.uk/efo/EFO_0007614 |
| 102 | cholesterol to total lipids in large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022236 |
| 102 | neurofilament light polypeptide amount | http://purl.obolibrary.org/obo/OBA_2050255 |
| 102 | non-melanoma skin carcinoma | http://www.ebi.ac.uk/efo/EFO_0009260 |
| 102 | pregnancy disorder | http://purl.obolibrary.org/obo/MONDO_0024575 |
| 102 | urinary albumin to creatinine ratio | http://www.ebi.ac.uk/efo/EFO_0007778 |
| 101 | cardioembolic stroke | http://www.ebi.ac.uk/efo/EFO_1001976 |
| 101 | docosahexaenoic acid to total fatty acids percentage | http://www.ebi.ac.uk/efo/EFO_0022262 |
| 101 | eosinophilic esophagitis | http://purl.obolibrary.org/obo/MONDO_0005361 |
| 101 | hemorrhoid | http://purl.obolibrary.org/obo/MONDO_0004872 |
| 101 | nodular goiter | http://purl.obolibrary.org/obo/MONDO_0006869 |
| 101 | plasma cell myeloma | http://purl.obolibrary.org/obo/MONDO_0009693 |
| 101 | response to carboplatin | http://purl.obolibrary.org/obo/GO_0097328 |
| 101 | t-tau measurement | http://www.ebi.ac.uk/efo/EFO_0004760 |
| 101 | tea consumption measurement | http://www.ebi.ac.uk/efo/EFO_0010091 |
| 100 | HMG CoA reductase inhibitor use measurement | http://www.ebi.ac.uk/efo/EFO_0009932 |
| 100 | colorectal adenoma | http://purl.obolibrary.org/obo/MONDO_0005484 |
| 100 | total lipids in VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022314 |
| 100 | total lipids in medium LDL | http://www.ebi.ac.uk/efo/EFO_0022180 |
| 99 | HIV-1 infection | http://www.ebi.ac.uk/efo/EFO_0000180 |
| 99 | RS-6-hydroxywarfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803328 |
| 99 | age at diagnosis | http://www.ebi.ac.uk/efo/EFO_0004918 |
| 99 | corpus callosum volume | http://purl.obolibrary.org/obo/OBA_2045249 |
| 99 | free cholesterol in large LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022176 |
| 99 | omega-6:omega-3 polyunsaturated fatty acid ratio | http://www.ebi.ac.uk/efo/EFO_0010732 |
| 98 | disease recurrence | http://www.ebi.ac.uk/efo/EFO_0004952 |
| 98 | empathy measurement | http://www.ebi.ac.uk/efo/EFO_0009183 |
| 98 | free cholesterol in IDL measurement | http://www.ebi.ac.uk/efo/EFO_0022181 |
| 98 | response to cisplatin | http://purl.obolibrary.org/obo/GO_0072718 |
| 98 | tissue plasminogen activator amount | http://purl.obolibrary.org/obo/OBA_2050207 |
| 97 | Hodgkins lymphoma | http://purl.obolibrary.org/obo/MONDO_0004952 |
| 97 | abdominal fat cell number | http://www.ebi.ac.uk/efo/EFO_0021534 |
| 97 | cholesterol in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022230 |
| 97 | essential hypertension | http://purl.obolibrary.org/obo/MONDO_0001134 |
| 97 | glucagon measurement | http://www.ebi.ac.uk/efo/EFO_0008463 |
| 97 | hepatocellular carcinoma | http://purl.obolibrary.org/obo/MONDO_0007256 |
| 97 | metabolic disease | http://purl.obolibrary.org/obo/MONDO_0005066 |
| 97 | ovarian serous carcinoma | http://www.ebi.ac.uk/efo/EFO_1001516 |
| 97 | response to methotrexate | http://purl.obolibrary.org/obo/GO_0031427 |
| 97 | response to simvastatin | http://purl.obolibrary.org/obo/GO_1903491 |
| 96 | cervical carcinoma | http://purl.obolibrary.org/obo/MONDO_0005131 |
| 96 | gait characteristic | http://purl.obolibrary.org/obo/OBA_1000869 |
| 96 | spontaneous preterm birth | http://www.ebi.ac.uk/efo/EFO_0006917 |
| 96 | trait in response to warfarin | http://purl.obolibrary.org/obo/OBA_2040000 |
| 95 | cholesteryl esters to total lipids in small LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022255 |
| 95 | (null) | http://www.ebi.ac.uk/efo/EFO_0010726 |
| 94 | amount of C-C motif chemokine 2 | http://purl.obolibrary.org/obo/OBA_2045327 |
| 94 | complement factor H-related protein 3 measurement | http://www.ebi.ac.uk/efo/EFO_0600056 |
| 93 | amount of phosphatidylethanolamine | http://purl.obolibrary.org/obo/OBA_2040180 |
| 93 | clostridium difficile infection | http://www.ebi.ac.uk/efo/EFO_0009130 |
| 92 | Astigmatism | http://purl.obolibrary.org/obo/HP_0000483 |
| 92 | total cholesterol change measurement | http://www.ebi.ac.uk/efo/EFO_0007806 |
| 91 | viral load | http://www.ebi.ac.uk/efo/EFO_0010125 |
| 90 | amount of chloride | http://purl.obolibrary.org/obo/OBA_1000959 |
| 90 | cannabis dependence | http://purl.obolibrary.org/obo/MONDO_0005689 |
| 90 | cholesteryl esters in large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022170 |
| 90 | longitudinal BMI measurement | http://www.ebi.ac.uk/efo/EFO_0005937 |
| 90 | neuroblastoma | http://purl.obolibrary.org/obo/MONDO_0005072 |
| 90 | trait in response to platinum atom | http://purl.obolibrary.org/obo/OBA_2040037 |
| 89 | blood sodium bicarbonate amount | http://purl.obolibrary.org/obo/OBA_2050094 |
| 89 | free cholesterol to total lipids in large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022281 |
| 89 | free cholesterol to total lipids in very small VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022290 |
| 89 | gastric cancer | http://purl.obolibrary.org/obo/MONDO_0001056 |
| 89 | photoreceptor cell layer thickness measurement | http://www.ebi.ac.uk/efo/EFO_0803370 |
| 89 | response to corticosteroid | http://purl.obolibrary.org/obo/GO_0031960 |
| 89 | sleep apnea measurement | http://www.ebi.ac.uk/efo/EFO_0007817 |
| 89 | susceptibility to childhood ear infection measurement | http://www.ebi.ac.uk/efo/EFO_0007904 |
| 89 | total lipids in very small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022156 |
| 88 | bitter alcoholic beverage consumption measurement | http://www.ebi.ac.uk/efo/EFO_0010092 |
| 88 | blood phosphate measurement | http://www.ebi.ac.uk/efo/EFO_0010972 |
| 88 | cholesteryl esters in IDL measurement | http://www.ebi.ac.uk/efo/EFO_0920029 |
| 88 | cholesteryl esters in medium VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022182 |
| 88 | emphysema imaging measurement | http://www.ebi.ac.uk/efo/EFO_0007626 |
| 88 | maximum cigarettes per day measurement | http://www.ebi.ac.uk/efo/EFO_0009264 |
| 88 | remnant cholesterol measurement | http://www.ebi.ac.uk/efo/EFO_0010815 |
| 88 | vitiligo | http://purl.obolibrary.org/obo/MONDO_0008661 |
| 87 | complex trait | http://www.ebi.ac.uk/efo/EFO_0010578 |
| 87 | essential tremor | http://purl.obolibrary.org/obo/MONDO_0003233 |
| 87 | gallstones | http://purl.obolibrary.org/obo/MONDO_0005346 |
| 87 | peripheral neuropathy | http://purl.obolibrary.org/obo/MONDO_0005244 |
| 87 | phospholipids in IDL measurement | http://www.ebi.ac.uk/efo/EFO_0022164 |
| 86 | B-cell chronic lymphocytic leukemia | http://purl.obolibrary.org/obo/MONDO_0004948 |
| 86 | ascending aorta diameter | http://purl.obolibrary.org/obo/OBA_2050001 |
| 86 | cholesterol to total lipids in very large HDL percentage | http://www.ebi.ac.uk/efo/EFO_0022243 |
| 86 | facial height measurement | http://www.ebi.ac.uk/efo/EFO_0007856 |
| 86 | head circumference | http://purl.obolibrary.org/obo/OBA_VT0000047 |
| 86 | preeclampsia | http://purl.obolibrary.org/obo/MONDO_0005081 |
| 86 | retinal nerve fibre layer thickness | http://purl.obolibrary.org/obo/OBA_2050111 |
| 86 | trait in response to tenofovir (anhydrous) | http://purl.obolibrary.org/obo/OBA_2040103 |
| 85 | free cholesterol in very small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022275 |
| 85 | peptic ulcer disease | http://purl.obolibrary.org/obo/MONDO_0004247 |
| 85 | response to beta blocker | http://www.ebi.ac.uk/efo/EFO_0007766 |
| 84 | RS-warfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803325 |
| 84 | adipose amount | http://purl.obolibrary.org/obo/OBA_VT0005452 |
| 84 | free cholesterol in small LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022271 |
| 84 | pathological myopia | http://www.ebi.ac.uk/efo/EFO_0004207 |
| 84 | total knee arthroplasty | http://www.ebi.ac.uk/efo/EFO_0020975 |
| 84 | verbal-numerical reasoning measurement | http://www.ebi.ac.uk/efo/EFO_0008394 |
| 83 | amount of alpha-fetoprotein | http://purl.obolibrary.org/obo/OBA_1000942 |
| 83 | behavioural disinhibition measurement | http://www.ebi.ac.uk/efo/EFO_0006946 |
| 83 | body height at birth | http://www.ebi.ac.uk/efo/EFO_0006784 |
| 83 | cholesterol in very small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022231 |
| 83 | cholesteryl esters in large LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022165 |
| 83 | circadian rhythm | http://purl.obolibrary.org/obo/GO_0007623 |
| 83 | response to cranial radiation therapy | http://www.ebi.ac.uk/efo/EFO_0010950 |
| 83 | skin pigmentation | http://purl.obolibrary.org/obo/OBA_VT0002095 |
| 82 | 1,5 anhydroglucitol measurement | http://www.ebi.ac.uk/efo/EFO_0008009 |
| 82 | able to hear with hearing aids | http://www.ebi.ac.uk/efo/EFO_0009720 |
| 82 | childhood trauma measurement | http://www.ebi.ac.uk/efo/EFO_0007979 |
| 82 | cholesterol in IDL measurement | http://www.ebi.ac.uk/efo/EFO_0021899 |
| 82 | cholesterol to total lipids in small LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022241 |
| 82 | free cholesterol in small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022272 |
| 82 | magnesium measurement | http://www.ebi.ac.uk/efo/EFO_0004845 |
| 82 | phospholipids in medium LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022183 |
| 82 | thyroid gland carcinoma | http://purl.obolibrary.org/obo/MONDO_0015075 |
| 81 | FEV change measurement | http://www.ebi.ac.uk/efo/EFO_0005921 |
| 81 | adverse effect | http://www.ebi.ac.uk/efo/EFO_0009658 |
| 81 | blood nickel amount | http://purl.obolibrary.org/obo/OBA_2050064 |
| 81 | cholesteryl esters in medium LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022186 |
| 81 | fibromyalgia | http://purl.obolibrary.org/obo/MONDO_0005546 |
| 81 | senile cataract | http://purl.obolibrary.org/obo/MONDO_0004847 |
| 80 | age at first birth measurement | http://www.ebi.ac.uk/efo/EFO_0009101 |
| 80 | cholesteryl esters in very large VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0920024 |
| 80 | fetal genotype effect measurement | http://www.ebi.ac.uk/efo/EFO_0007959 |
| 80 | free cholesterol in medium LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022268 |
| 80 | non-small cell lung carcinoma | http://purl.obolibrary.org/obo/MONDO_0005233 |
| 80 | phospholipids in small LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022297 |
| 80 | response to statin | http://purl.obolibrary.org/obo/GO_0036273 |
| 79 | bone mineral content measurement | http://www.ebi.ac.uk/efo/EFO_0007621 |
| 79 | cholesteryl esters in small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0920025 |
| 79 | energy intake | http://www.ebi.ac.uk/efo/EFO_0003939 |
| 79 | estrogen-receptor negative breast cancer | http://purl.obolibrary.org/obo/MONDO_0006513 |
| 79 | histidine measurement | http://www.ebi.ac.uk/efo/EFO_0009769 |
| 79 | nephrolithiasis | http://purl.obolibrary.org/obo/MONDO_0008171 |
| 79 | skin disorder | http://purl.obolibrary.org/obo/MONDO_0005093 |
| 79 | sugar consumption measurement | http://www.ebi.ac.uk/efo/EFO_0010158 |
| 78 | Beta blocking agent use measurement | http://www.ebi.ac.uk/efo/EFO_0009929 |
| 78 | LDL particle size | http://purl.obolibrary.org/obo/OBA_2050113 |
| 78 | antihyperlipidemic drug use measurement | http://www.ebi.ac.uk/efo/EFO_0803367 |
| 78 | bone fracture | http://purl.obolibrary.org/obo/MONDO_0005315 |
| 78 | concentration of very small VLDL particles | http://www.ebi.ac.uk/efo/EFO_0022151 |
| 78 | free cholesterol in VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022276 |
| 78 | free cholesterol in medium VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022269 |
| 78 | isoleucine measurement | http://www.ebi.ac.uk/efo/EFO_0009793 |
| 78 | lymphoid neoplasm | http://purl.obolibrary.org/obo/MONDO_0005157 |
| 77 | lung disease severity measurement | http://www.ebi.ac.uk/efo/EFO_0007744 |
| 77 | phospholipids in large LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022174 |
| 77 | phospholipids in very small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022300 |
| 77 | total lipids in lipoprotein particles measurement | http://www.ebi.ac.uk/efo/EFO_0022309 |
| 76 | apolipoprotein B to apolipoprotein A1 ratio | http://www.ebi.ac.uk/efo/EFO_0021897 |
| 76 | cholesteryl esters to total lipids in medium LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022252 |
| 76 | estradiol measurement | http://www.ebi.ac.uk/efo/EFO_0004697 |
| 76 | free cholesterol in LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022266 |
| 76 | human papilloma virus infection | http://purl.obolibrary.org/obo/MONDO_0005161 |
| 76 | lactate measurement | http://www.ebi.ac.uk/efo/EFO_0007745 |
| 76 | milk amount | http://purl.obolibrary.org/obo/OBA_VT0015043 |
| 76 | obstructive sleep apnea syndrome | http://purl.obolibrary.org/obo/MONDO_0007147 |
| 76 | subcutaneous adipose tissue quality | http://purl.obolibrary.org/obo/OBA_1000023 |
| 76 | velopharyngeal dysfunction | http://www.ebi.ac.uk/efo/EFO_0009336 |
| 75 | alcoholic liver disease | http://purl.obolibrary.org/obo/MONDO_0043693 |
| 75 | cholesterol in large LDL measurement | http://www.ebi.ac.uk/efo/EFO_0021901 |
| 75 | cholesterol in small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022228 |
| 75 | electroencephalogram measurement | http://www.ebi.ac.uk/efo/EFO_0004357 |
| 75 | interleukin 12 measurement | http://www.ebi.ac.uk/efo/EFO_0004753 |
| 75 | keratinocyte carcinoma | http://www.ebi.ac.uk/efo/EFO_0010176 |
| 75 | polycystic ovary syndrome | http://purl.obolibrary.org/obo/MONDO_0008487 |
| 75 | trait in response to hydrochlorothiazide | http://purl.obolibrary.org/obo/OBA_2040093 |
| 75 | triglycerides to total lipids in large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022332 |
| 74 | calvaria morphology trait | http://purl.obolibrary.org/obo/OBA_VT0000076 |
| 74 | creatinine clearance measurement | http://www.ebi.ac.uk/efo/EFO_0007934 |
| 74 | cytotoxicity measurement | http://www.ebi.ac.uk/efo/EFO_0006952 |
| 74 | hormone measurement | http://www.ebi.ac.uk/efo/EFO_0004730 |
| 74 | osteonecrosis | http://purl.obolibrary.org/obo/MONDO_0005380 |
| 74 | phenylalanine measurement | http://www.ebi.ac.uk/efo/EFO_0005001 |
| 74 | primary biliary cholangitis | http://purl.obolibrary.org/obo/MONDO_0005388 |
| 74 | rosacea severity measurement | http://www.ebi.ac.uk/efo/EFO_0009180 |
| 74 | serum IgA amount | http://purl.obolibrary.org/obo/OBA_2050074 |
| 74 | social deprivation | http://www.ebi.ac.uk/efo/EFO_0009696 |
| 74 | total cholesterol in medium LDL | http://www.ebi.ac.uk/efo/EFO_0022185 |
| 73 | Dupuytren Contracture | http://www.ebi.ac.uk/efo/EFO_0004229 |
| 73 | Inguinal hernia | http://purl.obolibrary.org/obo/HP_0000023 |
| 73 | melanoma | http://purl.obolibrary.org/obo/MONDO_0005105 |
| 73 | optic nerve glioma | http://purl.obolibrary.org/obo/MONDO_0003235 |
| 73 | phospholipids in small VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022146 |
| 72 | R-6-hydroxywarfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803327 |
| 72 | acylcarnitine measurement | http://www.ebi.ac.uk/efo/EFO_0005059 |
| 72 | beta-amyloid 1-42 measurement | http://www.ebi.ac.uk/efo/EFO_0004670 |
| 72 | blood pressure trait | http://purl.obolibrary.org/obo/OBA_VT0000183 |
| 72 | cholesterol to total lipids in medium LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022238 |
| 72 | diffusing capacity of the lung for carbon monoxide | http://www.ebi.ac.uk/efo/EFO_0009369 |
| 72 | peptide measurement | http://www.ebi.ac.uk/efo/EFO_0010520 |
| 72 | response to anticonvulsant | http://purl.obolibrary.org/obo/GO_0036277 |
| 71 | Drug-induced agranulocytosis | http://purl.obolibrary.org/obo/HP_0012235 |
| 71 | cholesteryl esters to total lipids in large LDL percentage | http://www.ebi.ac.uk/efo/EFO_0022249 |
| 71 | concentration of medium VLDL particles | http://www.ebi.ac.uk/efo/EFO_0022152 |
| 71 | corneal topography | http://www.ebi.ac.uk/efo/EFO_0004345 |
| 71 | glucose tolerance test | http://www.ebi.ac.uk/efo/EFO_0004307 |
| 71 | malaria | http://purl.obolibrary.org/obo/MONDO_0005136 |
| 71 | mood instability measurement | http://www.ebi.ac.uk/efo/EFO_0008475 |
| 71 | pursuit maintenance gain measurement | http://www.ebi.ac.uk/efo/EFO_0008433 |
| 70 | albuminuria | http://www.ebi.ac.uk/efo/EFO_0004285 |
| 70 | carotid plaque build | http://www.ebi.ac.uk/efo/EFO_0006501 |
| 70 | cholelithiasis | http://purl.obolibrary.org/obo/MONDO_0012672 |
| 70 | keratoconus | http://purl.obolibrary.org/obo/MONDO_0015486 |
| 70 | phospholipids in medium VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022154 |
| 70 | total cholesterol in small LDL | http://www.ebi.ac.uk/efo/EFO_0022179 |
| 70 | triglyceride change measurement | http://www.ebi.ac.uk/efo/EFO_0007681 |
| 69 | acne | http://purl.obolibrary.org/obo/MONDO_0011438 |
| 69 | gluteofemoral adipose tissue measurement | http://www.ebi.ac.uk/efo/EFO_0803317 |
| 69 | insulin sensitivity measurement | http://www.ebi.ac.uk/efo/EFO_0004471 |
| 69 | interferon gamma measurement | http://www.ebi.ac.uk/efo/EFO_0008165 |
| 69 | squamous cell carcinoma | http://purl.obolibrary.org/obo/MONDO_0005096 |
| 68 | S-7-hydroxywarfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803329 |
| 68 | nevus count | http://www.ebi.ac.uk/efo/EFO_0004632 |
| 68 | renal sinus volume | http://purl.obolibrary.org/obo/OBA_2100005 |
| 68 | seasonal gut microbiome measurement | http://www.ebi.ac.uk/efo/EFO_0007753 |
| 67 | RS-10-hydroxywarfarin to RS-warfarin ratio measurement | http://www.ebi.ac.uk/efo/EFO_0803335 |
| 67 | amount of vascular endothelial growth factor A | http://purl.obolibrary.org/obo/OBA_2050131 |
| 67 | atopic asthma | http://www.ebi.ac.uk/efo/EFO_0010638 |
| 67 | blood cobalt amount | http://purl.obolibrary.org/obo/OBA_2050077 |
| 67 | eye disorder | http://purl.obolibrary.org/obo/MONDO_0005328 |
| 67 | family history | http://www.ebi.ac.uk/efo/EFO_0000493 |
| 67 | ferritin measurement | http://www.ebi.ac.uk/efo/EFO_0004459 |
| 67 | glucose homeostasis trait | http://purl.obolibrary.org/obo/OBA_1000528 |
| 67 | ocular sarcoidosis | http://www.ebi.ac.uk/efo/EFO_0010723 |
| 67 | quality of life cycle | http://purl.obolibrary.org/obo/OBA_2045261 |
| 67 | systemic sclerosis | http://purl.obolibrary.org/obo/MONDO_0005100 |
| 67 | trait in response to citalopram | http://purl.obolibrary.org/obo/OBA_2040057 |
| 66 | S-warfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803323 |
| 66 | androgenetic alopecia | http://purl.obolibrary.org/obo/MONDO_0005339 |
| 66 | endometrial cancer | http://purl.obolibrary.org/obo/MONDO_0011962 |
| 66 | facial width measurement | http://www.ebi.ac.uk/efo/EFO_0007855 |
| 66 | leprosy | http://purl.obolibrary.org/obo/MONDO_0005124 |
| 66 | loneliness measurement | http://www.ebi.ac.uk/efo/EFO_0007865 |
| 66 | stressful life event measurement | http://www.ebi.ac.uk/efo/EFO_0007781 |
| 66 | thyroid cancer | http://purl.obolibrary.org/obo/MONDO_0002108 |
| 65 | age of onset of childhood onset asthma | http://purl.obolibrary.org/obo/OBA_2001015 |
| 65 | cholesterol in small LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022227 |
| 65 | response to TNF antagonist | http://www.ebi.ac.uk/efo/EFO_0004653 |
| 65 | response to clozapine | http://purl.obolibrary.org/obo/GO_0097338 |
| 64 | HOMA-B | http://www.ebi.ac.uk/efo/EFO_0004469 |
| 64 | HOMA-IR | http://www.ebi.ac.uk/efo/EFO_0004501 |
| 64 | cirrhosis of liver | http://purl.obolibrary.org/obo/MONDO_0005155 |
| 64 | osteoporosis | http://purl.obolibrary.org/obo/MONDO_0005298 |
| 63 | Moderate albuminuria | http://purl.obolibrary.org/obo/HP_0012594 |
| 63 | Progressive visual loss | http://purl.obolibrary.org/obo/HP_0000529 |
| 63 | concentration of large LDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022160 |
| 63 | estrogen-receptor positive breast cancer | http://purl.obolibrary.org/obo/MONDO_0006512 |
| 63 | interleukin 18 measurement | http://www.ebi.ac.uk/efo/EFO_0004581 |
| 63 | response to ketamine | http://www.ebi.ac.uk/efo/EFO_0009748 |
| 63 | total lipids in LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022308 |
| 63 | worry measurement | http://www.ebi.ac.uk/efo/EFO_0009589 |
| 62 | X-11315 measurement | http://www.ebi.ac.uk/efo/EFO_0021223 |
| 62 | influenza A (H1N1) | http://www.ebi.ac.uk/efo/EFO_1001488 |
| 62 | interleukin-1 beta measurement | http://www.ebi.ac.uk/efo/EFO_0004812 |
| 62 | response to sulfonylurea | http://www.ebi.ac.uk/efo/EFO_0007922 |
| 61 | Hypercholesterolemia | http://purl.obolibrary.org/obo/HP_0003124 |
| 61 | R-6-hydroxywarfarin to R-warfarin ratio measurement | http://www.ebi.ac.uk/efo/EFO_0803333 |
| 61 | cholesterol to total lipids in chylomicrons and extremely large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022232 |
| 61 | conduct disorder | http://purl.obolibrary.org/obo/MONDO_0005352 |
| 61 | endometrial carcinoma | http://purl.obolibrary.org/obo/MONDO_0002447 |
| 61 | head size trait | http://purl.obolibrary.org/obo/OBA_VT0000002 |
| 61 | renal carcinoma | http://purl.obolibrary.org/obo/MONDO_0005206 |
| 60 | S-warfarin to R-warfarin ratio measurement | http://www.ebi.ac.uk/efo/EFO_0803331 |
| 60 | aseptic loosening | http://www.ebi.ac.uk/efo/EFO_0010725 |
| 60 | economic and social preference | http://www.ebi.ac.uk/efo/EFO_0004827 |
| 60 | interleukin 7 measurement | http://www.ebi.ac.uk/efo/EFO_0008189 |
| 60 | phospholipids in LDL measurement | http://www.ebi.ac.uk/efo/EFO_0022294 |
| 60 | serum dimethylarginine amount | http://purl.obolibrary.org/obo/OBA_2050060 |
| 59 | RS-10-hydroxywarfarin measurement | http://www.ebi.ac.uk/efo/EFO_0803330 |
| 59 | axial length measurement | http://www.ebi.ac.uk/efo/EFO_0005318 |
| 59 | cholesteryl esters to total lipids in chylomicrons and extremely large VLDL percentage | http://www.ebi.ac.uk/efo/EFO_0022246 |
| 59 | drug-induced liver injury | http://purl.obolibrary.org/obo/MONDO_0005359 |
| 59 | dry eye syndrome | http://purl.obolibrary.org/obo/MONDO_0006733 |
| 59 | neuritic plaque measurement | http://www.ebi.ac.uk/efo/EFO_0006798 |
| 59 | skin cancer | http://purl.obolibrary.org/obo/MONDO_0002898 |
| 58 | Bell's palsy | http://purl.obolibrary.org/obo/MONDO_0005665 |
| 58 | Self-injurious behavior | http://purl.obolibrary.org/obo/HP_0100716 |
| 58 | Takayasu arteritis | http://purl.obolibrary.org/obo/MONDO_0017991 |
| 58 | fish consumption measurement | http://www.ebi.ac.uk/efo/EFO_0010139 |
| 58 | hearing threshold trait | http://purl.obolibrary.org/obo/OBA_2045234 |
| 58 | neck of femur size | http://purl.obolibrary.org/obo/OBA_0002725 |
| 58 | oropharynx cancer | http://purl.obolibrary.org/obo/MONDO_0004608 |
| 58 | response to opioid | http://www.ebi.ac.uk/efo/EFO_0008541 |
| 58 | total cholesterol in IDL | http://www.ebi.ac.uk/efo/EFO_0022184 |
| 58 | vascular endothelial function measurement | http://www.ebi.ac.uk/efo/EFO_0803369 |
| 57 | acute lymphoblastic leukemia | http://purl.obolibrary.org/obo/MONDO_0004967 |
| 57 | developmental dysplasia of the hip | http://purl.obolibrary.org/obo/MONDO_0000158 |
| 57 | interleukin-8 measurement | http://www.ebi.ac.uk/efo/EFO_0004811 |
| 57 | kit ligand amount | http://purl.obolibrary.org/obo/OBA_2050260 |
| 57 | leptin measurement | http://www.ebi.ac.uk/efo/EFO_0005000 |
| 57 | non-Hodgkin lymphoma | http://purl.obolibrary.org/obo/MONDO_0018908 |
| 57 | response to exercise | http://www.ebi.ac.uk/efo/EFO_0007768 |
| 57 | serum IgM amount | http://purl.obolibrary.org/obo/OBA_2050070 |
| 57 | social interaction measurement | http://www.ebi.ac.uk/efo/EFO_0009592 |
| 57 | sunburn | http://purl.obolibrary.org/obo/MONDO_0005326 |
| 57 | thalamus volume | http://www.ebi.ac.uk/efo/EFO_0006935 |
| 57 | urinary uromodulin measurement | http://www.ebi.ac.uk/efo/EFO_0005663 |
| 56 | Behcet disease | http://purl.obolibrary.org/obo/MONDO_0007191 |
| 56 | cholesteryl esters in small LDL measurement | http://www.ebi.ac.uk/efo/EFO_0920028 |
| 56 | concentration of small LDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022166 |
| 56 | idiopathic pulmonary fibrosis | http://www.ebi.ac.uk/efo/EFO_0000768 |
| 56 | lipoprotein A measurement | http://www.ebi.ac.uk/efo/EFO_0006925 |
| 56 | nephrotic syndrome | http://purl.obolibrary.org/obo/MONDO_0005377 |
| 56 | polychlorinated biphenyls measurement | http://www.ebi.ac.uk/efo/EFO_0007042 |
| 56 | skin neoplasm | http://purl.obolibrary.org/obo/MONDO_0002531 |
| 56 | total cholesterol in small VLDL | http://www.ebi.ac.uk/efo/EFO_0022171 |
| 56 | trauma complication | http://www.ebi.ac.uk/efo/EFO_0009888 |
| 55 | Abnormal skin pigmentation | http://purl.obolibrary.org/obo/HP_0001000 |
| 55 | Back pain | http://purl.obolibrary.org/obo/HP_0003418 |
| 55 | PCSK9 protein measurement | http://www.ebi.ac.uk/efo/EFO_0006899 |
| 55 | amount of galanin peptides (human) in blood | http://purl.obolibrary.org/obo/OBA_2052182 |
| 55 | anti-hepatitis E virus antibody measurement | http://www.ebi.ac.uk/efo/EFO_0803690 |
| 55 | complement factor D measurement | http://www.ebi.ac.uk/efo/EFO_0020283 |
| 55 | digestive system disorder | http://purl.obolibrary.org/obo/MONDO_0004335 |
| 55 | dimethylglycine measurement | http://www.ebi.ac.uk/efo/EFO_0010476 |
| 55 | erectile dysfunction | http://purl.obolibrary.org/obo/MONDO_0005362 |
| 55 | hyperlipidemia | http://purl.obolibrary.org/obo/MONDO_0021187 |
| 55 | initial pursuit acceleration | http://www.ebi.ac.uk/efo/EFO_0008434 |
| 55 | joint damage measurement | http://www.ebi.ac.uk/efo/EFO_0005413 |
| 55 | obsolete amount of tumor necrosis factor alpha | http://purl.obolibrary.org/obo/OBA_1001013 |
| 54 | 3-hydroxybutyrate measurement | http://www.ebi.ac.uk/efo/EFO_0010982 |
| 54 | Inhalant adrenergic use measurement | http://www.ebi.ac.uk/efo/EFO_0009941 |
| 54 | amount of diglyceride | http://purl.obolibrary.org/obo/OBA_2040178 |
| 54 | beverage consumption measurement | http://www.ebi.ac.uk/efo/EFO_0010088 |
| 54 | cholesterol in medium VLDL measurement | http://www.ebi.ac.uk/efo/EFO_0022225 |
| 54 | eye color | http://www.ebi.ac.uk/efo/EFO_0003949 |
| 54 | ghrelin measurement | http://www.ebi.ac.uk/efo/EFO_0600001 |
| 54 | nucleotide measurement | http://www.ebi.ac.uk/efo/EFO_0010513 |
| 54 | opioid use disorder | http://www.ebi.ac.uk/efo/EFO_0010702 |
| 54 | poisoning | http://purl.obolibrary.org/obo/MONDO_0029000 |
| 54 | response to antihypertensive drug | http://www.ebi.ac.uk/efo/EFO_0005405 |
| 54 | response to diuretic | http://purl.obolibrary.org/obo/GO_0036270 |
| 54 | sensory perception of smell | http://purl.obolibrary.org/obo/GO_0007608 |
| 54 | snoring measurement | http://www.ebi.ac.uk/efo/EFO_0008341 |
| 54 | suicide behaviour measurement | http://www.ebi.ac.uk/efo/EFO_0006882 |
| 53 | 5-HIAA measurement | http://www.ebi.ac.uk/efo/EFO_0005132 |
| 53 | adult onset asthma | http://www.ebi.ac.uk/efo/EFO_1002011 |
| 53 | amount of C-C motif chemokine 4 | http://purl.obolibrary.org/obo/OBA_2045334 |
| 53 | amount of creatine | http://purl.obolibrary.org/obo/OBA_1001005 |
| 53 | carnitine measurement | http://www.ebi.ac.uk/efo/EFO_0010469 |
| 53 | carpal tunnel syndrome | http://purl.obolibrary.org/obo/MONDO_0007275 |
| 53 | concentration of medium LDL particles measurement | http://www.ebi.ac.uk/efo/EFO_0022172 |
| 53 | descending aorta diameter | http://purl.obolibrary.org/obo/OBA_2050002 |
| 53 | diabetic neuropathy | http://purl.obolibrary.org/obo/MONDO_0006626 |
| 53 | hearing process quality | http://purl.obolibrary.org/obo/OBA_2045262 |
| 53 | interleukin 17 measurement | http://www.ebi.ac.uk/efo/EFO_0008174 |
| 53 | mean fractional anisotropy measurement | http://www.ebi.ac.uk/efo/EFO_0008399 |
| 52 | DDT metabolite measurement | http://www.ebi.ac.uk/efo/EFO_0007886 |
| 52 | S-7-hydroxywarfarin to S-warfarin ratio measurement | http://www.ebi.ac.uk/efo/EFO_0803334 |
| 52 | amount of C-X-C motif chemokine 6 | http://purl.obolibrary.org/obo/OBA_2045346 |
| 52 | amount of C-type lectin domain family 11 member A | http://purl.obolibrary.org/obo/OBA_2045342 |
| 52 | amount of platelet-derived growth factor complex BB dimer (human) | http://purl.obolibrary.org/obo/OBA_2050167 |
| 52 | clear cell renal carcinoma | http://purl.obolibrary.org/obo/MONDO_0005005 |
| 52 | liver fibrosis measurement | http://www.ebi.ac.uk/efo/EFO_0010576 |
| 52 | mammographic density measurement | http://www.ebi.ac.uk/efo/EFO_0005941 |
| 52 | monoclonal gammopathy | http://purl.obolibrary.org/obo/MONDO_0004960 |
| 52 | otosclerosis | http://purl.obolibrary.org/obo/MONDO_0005349 |
| 52 | radiation injury | http://purl.obolibrary.org/obo/MONDO_0043458 |
| 52 | response to combination chemotherapy | http://www.ebi.ac.uk/efo/EFO_0007965 |
| 52 | response to lithium ion | http://purl.obolibrary.org/obo/GO_0010226 |
| 52 | seropositivity measurement | http://www.ebi.ac.uk/efo/EFO_0007034 |
| 52 | testicular germ cell tumor | http://purl.obolibrary.org/obo/MONDO_0010108 |
| 51 | 2,4,5-trichlorophenol measurement | http://www.ebi.ac.uk/efo/EFO_0022052 |
| 51 | 4,6-dinitro-o-cresol measurement | http://www.ebi.ac.uk/efo/EFO_0022053 |
| 51 | Cannabis use | http://www.ebi.ac.uk/efo/EFO_0007585 |
| 51 | Corneal astigmatism | http://www.ebi.ac.uk/efo/EFO_1002040 |
| 51 | MHPG measurement | http://www.ebi.ac.uk/efo/EFO_0005133 |
| 51 | YKL40 measurement | http://www.ebi.ac.uk/efo/EFO_0004869 |
| 51 | age of onset of asthma | http://purl.obolibrary.org/obo/OBA_2001001 |
| 51 | aldrin measurement | http://www.ebi.ac.uk/efo/EFO_0022054 |
| 51 | amount of advanced glycosylation end product-specific receptor in blood | http://purl.obolibrary.org/obo/OBA_2051078 |
| 51 | amount of desmoglein-4 (human) in blood | http://purl.obolibrary.org/obo/OBA_2052063 |
| 51 | azinphos methyl measurement | http://www.ebi.ac.uk/efo/EFO_0022055 |
| 51 | cadmium chloride measurement | http://www.ebi.ac.uk/efo/EFO_0022056 |
| 51 | chlorpyrifos measurement | http://www.ebi.ac.uk/efo/EFO_0022057 |
| 51 | diazinon measurement | http://www.ebi.ac.uk/efo/EFO_0022059 |
| 51 | dibutyl phthalate measurement | http://www.ebi.ac.uk/efo/EFO_0022058 |
| 51 | dicofol measurement | http://www.ebi.ac.uk/efo/EFO_0022060 |
| 51 | dieldrin measurement | http://www.ebi.ac.uk/efo/EFO_0022061 |
| 51 | disulfoton measurement | http://www.ebi.ac.uk/efo/EFO_0022062 |
| 51 | endosulfan measurement | http://www.ebi.ac.uk/efo/EFO_0022063 |
| 51 | endrin measurement | http://www.ebi.ac.uk/efo/EFO_0022064 |
| 51 | erosive tooth wear attribute | http://purl.obolibrary.org/obo/OBA_2045285 |
| 51 | ethion measurement | http://www.ebi.ac.uk/efo/EFO_0022065 |
| 51 | heptachlor epoxide measurement | http://www.ebi.ac.uk/efo/EFO_0022067 |
| 51 | heptachlor measurement | http://www.ebi.ac.uk/efo/EFO_0022066 |
| 51 | lymphotoxin-alpha amount | http://purl.obolibrary.org/obo/OBA_2050291 |
| 51 | mercuric chloride measurement | http://www.ebi.ac.uk/efo/EFO_0022068 |
| 51 | methoxychlor measurement | http://www.ebi.ac.uk/efo/EFO_0022069 |
| 51 | multisite chronic pain | http://www.ebi.ac.uk/efo/EFO_0010100 |
| 51 | parathion measurement | http://www.ebi.ac.uk/efo/EFO_0022070 |
| 51 | pentachlorophenol measurement | http://www.ebi.ac.uk/efo/EFO_0022071 |
| 51 | post-traumatic stress disorder symptom measurement | http://www.ebi.ac.uk/efo/EFO_0008535 |
| 51 | potassium chromate measurement | http://www.ebi.ac.uk/efo/EFO_0022072 |
| 51 | reading | http://www.ebi.ac.uk/efo/EFO_0005229 |
| 51 | response to trauma exposure | http://www.ebi.ac.uk/efo/EFO_0008483 |
| 51 | rostrum of corpus callosum volume | http://purl.obolibrary.org/obo/OBA_2045228 |
| 51 | sensory perception of sweet taste | http://purl.obolibrary.org/obo/GO_0050916 |
| 51 | systemic-onset juvenile idiopathic arthritis | http://purl.obolibrary.org/obo/MONDO_0019434 |
| 51 | total cholesterol in large LDL | http://www.ebi.ac.uk/efo/EFO_0022167 |
| 51 | trait in response to thiazide | http://purl.obolibrary.org/obo/OBA_2040083 |
| 51 | wet macular degeneration | http://purl.obolibrary.org/obo/MONDO_0005417 |
| 50 | Epstein-Barr virus infection | http://purl.obolibrary.org/obo/MONDO_0005111 |
| 50 | amount of fms-related tyrosine kinase 3 ligand (human) in blood | http://purl.obolibrary.org/obo/OBA_2052944 |
| 50 | autosomal dominant compelling helio-ophthalmic outburst syndrome | http://www.ebi.ac.uk/efo/EFO_0007887 |
| 50 | chronic rhinosinusitis with nasal polyps | http://www.ebi.ac.uk/efo/EFO_1002029 |
| 50 | diabetic eye disease | http://www.ebi.ac.uk/efo/EFO_0009486 |
| 50 | fms-related tyrosine kinase 3 ligand measurement | http://www.ebi.ac.uk/efo/EFO_0010785 |
| 50 | polyp of colon | http://purl.obolibrary.org/obo/MONDO_0021400 |
| 50 | renal transplant outcome measurement | http://www.ebi.ac.uk/efo/EFO_0005199 |
| 50 | serum iron amount | http://purl.obolibrary.org/obo/OBA_2050055 |
| 50 | systolic blood pressure change measurement | http://www.ebi.ac.uk/efo/EFO_0006944 |
