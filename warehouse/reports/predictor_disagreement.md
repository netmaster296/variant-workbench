# Predictor disagreement

## Question

On Tier A coding variants, where do CADD, REVEL, SpliceAI, a MaveDB score, and ClinVar stars disagree?

Agreement is the set to read first. Disagreement is the set worth looking at. A star count is confidence in a ClinVar label, not a deleteriousness score, so stars disagree with a predictor only when a high-star pathogenic or benign label and that predictor's call point opposite ways.

## Filters

Coding is taken from the stored VEP consequence, not from `vep_impact`. `scripts/annotate_gnomad.py` kept one VEP block per site: a non-empty `MANE_SELECT` if any block has one, otherwise `CANONICAL=YES`, otherwise the highest `IMPACT`. The column is `vep_consequence` on the genome side and `vep_consequence_exome` on the exome side. A site is coding on that side when the site is in that callset (`in_genomes` or `in_exomes`) and the consequence, split on `&`, contains any of:

`missense_variant`, `synonymous_variant`, `stop_gained`, `stop_lost`, `start_lost`, `start_retained_variant`, `stop_retained_variant`, `frameshift_variant`, `inframe_insertion`, `inframe_deletion`, `protein_altering_variant`, `incomplete_terminal_codon_variant`, `coding_sequence_variant`.

Splice acceptor and splice donor sites with no coding-sequence term stay out. So do intron, UTR, upstream, downstream, and intergenic terms. `NMD_transcript_variant` does not make a site coding by itself.

Genome columns and exome columns are never added and never averaged. A key that is coding on both sides is counted once in each denominator. Predictor names in the files are `cadd_phred`, `revel_max`, `spliceai_ds_max`, `pangolin_largest_ds`, `sift_max`, and `polyphen_max`, with `_exome` on the exome copy.

Calls, used only when the score is non-null:

- REVEL pathogenic-supporting: `revel_max` >= 0.644. Benign-supporting: <= 0.290. The open interval is not a call. These are the supporting cuts Pejaver et al. (Am J Hum Genet, 2022) published for REVEL and that the UCSC REVEL track attributes to that paper's Table 2. Moderate (>= 0.773) and strong (>= 0.932) pathogenic cuts are reported only as a sensitivity check; Hopkins et al. (Hum Mutat, 2023) restate those two pathogenic cuts from the same paper.
- CADD high: phred >= 20, the top 1% of the phred ranking. CADD low: phred < 10, outside the top 10%. The band from 10 up to but not including 20 is not a call. Phred >= 25 and >= 30 are sensitivity cuts, not the primary call.
- SpliceAI uses the three operating points gnomAD documents for `spliceai_ds_max`: 0.2 (high recall), 0.5 (recommended), 0.8 (high precision). The primary splice call is >= 0.5. A score under 0.2 is not treated as contradicting a high REVEL, because the tools are not answering the same question.
- Pangolin values in this file go below zero (splice loss). Agreement with SpliceAI uses `abs(pangolin_largest_ds)` at 0.5 and 0.2. That is the SpliceAI numeric point applied to the absolute value. It is not a published Pangolin cutoff.
- SIFT and PolyPhen are counted for missingness only. No SIFT or PolyPhen call is made here. gnomAD describes `sift_max` as 0 to 1, with scores below 0.05 predicted to affect protein function. That cut is not applied.

REVEL is compared only on rows whose consequence contains `missense_variant`. CADD and SpliceAI comparisons that do not use REVEL use all coding rows.

ClinVar is `warehouse/clinvar/allele_dedup.parquet` joined on `variant_key`. The join did not multiply rows. `review_stars` is the warehouse star bucket. Star 1 is split on the stored `ReviewStatus` text into single submitter versus conflicting classifications. A directional class requires germline origin (`germline_somatic = germline`), no oncogenicity call, and no conflict piece:

- P/LP: a `;`-separated piece contains `pathogenic` and does not contain `conflict`.
- B/LB: a piece contains `benign` and contains neither `conflict` nor `pathogenic`.
- A conflict piece, or both a pathogenic piece and a benign piece, is `conflict_or_mixed`.

This matched the stored `germline_pathogenic` flag on every coding row (0 rows where the flag and the P/LP class differed). High-star means stars >= 2, the review text does not contain `conflict`, and the class is P/LP or B/LB. Every coding ClinVar key had `n_rows = 1`, so the deduplicated star is not a max across mixed records.

MaveDB uses `mapped_variants_scored.parquet` allele rows only (`post_mapped_type = Allele`) with a non-null numeric `score`. A haplotype score is not copied onto a member allele. A `(score_set_urn, variant_key)` pair is kept only when every numeric score on that pair is the same value. Scores are never pooled, ranked, or subtracted across score sets. Inside one score set, a tail "disagrees" when its score sits on the other tail's side of that other tail's median. A set is treated as separated only when each tail has enough variants and at most 25% of each tail is on the wrong side. The median comparison is inside that urn only.

## Denominator

Tier A has 5,928,738 keys, each once. A key absent from a callset is not in that side's denominator. It is not a missing score.

| side | in that callset | coding | of which missense | on the chosen block's MANE transcript | also coding on the other side |
| --- | ---: | ---: | ---: | ---: | ---: |
| genome | 3,220,597 | 1,624,724 | 1,206,859 | 1,601,015 | 1,542,356 |
| exome | 3,283,322 | 2,552,394 | 1,862,583 | 2,520,283 | 1,542,356 |

The coding denominators are the two populations below. They are not added. Of the keys coding on both sides, CADD differed on 0 of 1,542,356 where both scores exist, REVEL differed on 0 of 981,190, and SpliceAI differed on 0 of 1,530,214. One key, `chr22:50697556:G:GGCCCCGGCCCCGCGCCCGGC` (SHANK3), has the same three consequence terms in a different order. The callsets agree where they overlap. The exome coding set is larger: 1,542,356 / 2,552,394 of exome coding keys are also genome coding, against 1,542,356 / 1,624,724 of genome coding keys.

HIGH or MODERATE impact that failed the coding-term rule: 12,497 genome rows and 35,701 exome rows. Splice acceptor or donor with no coding-sequence term: 12,486 and 35,701. Those sites are not in the denominator.

## Result

### Missingness

gnomAD does not put every tool on every site. On Tier A rows that are simply present in the callset, CADD is the exception: `cadd_phred` and `cadd_phred_exome` are populated on every in-callset row (0 missing of 3,220,597 genome and 0 of 3,283,322 exome). The other tools are not. Missing on the in-callset rows:

| side | REVEL missing | SpliceAI missing | Pangolin missing | SIFT missing | PolyPhen missing |
| --- | ---: | ---: | ---: | ---: | ---: |
| genome | 2,185,785 / 3,220,597 | 657,704 / 3,220,597 | 386,462 / 3,220,597 | 2,049,245 / 3,220,597 | 2,042,919 / 3,220,597 |
| exome | 1,684,796 / 3,283,322 | 83,366 / 3,283,322 | 12,135 / 3,283,322 | 1,429,870 / 3,283,322 | 1,419,993 / 3,283,322 |

Coding rows are a more completely annotated slice. CADD is still complete. REVEL is missing mostly because it is a missense score. SIFT and PolyPhen behave the same way. Missing counts:

| side | denominator | CADD | REVEL | SpliceAI | Pangolin | SIFT | PolyPhen | no ClinVar row | ClinVar row, stars null | any MaveDB allele score |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| genome | 1,624,724 coding | 0 | 593,138 | 12,768 | 1,523 | 456,434 | 450,199 | 13,196 | 16,891 | 4,246 present |
| exome | 2,552,394 coding | 0 | 958,338 | 18,988 | 116 | 703,295 | 693,539 | 19,810 | 33,496 | 13,337 present |
| genome | 1,206,859 missense | 0 | 177,270 | 8,240 | 948 | 39,596 | 33,425 |  |  |  |
| exome | 1,862,583 missense | 0 | 272,273 | 12,202 | 99 | 16,039 | 6,397 |  |  |  |

REVEL is missing on 177,270 / 1,206,859 genome missense rows and 272,273 / 1,862,583 exome missense rows (0.1469 and 0.1462). A later disagreement rate uses only the rows both tools actually scored.

ClinVar stars on coding rows:

| bucket | genome | exome |
| --- | ---: | ---: |
| no ClinVar row | 13,196 | 19,810 |
| stars null | 16,891 | 33,496 |
| 0 | 24,477 | 39,855 |
| 1, single submitter | 1,098,222 | 1,820,698 |
| 1, conflicting classifications | 110,651 | 137,278 |
| 2 | 356,173 | 492,414 |
| 3 | 5,080 | 8,807 |
| 4 | 34 | 36 |

No coding row had both the single-submitter phrase and the conflicting phrase. Star 1 single submitter is 1,098,222 / 1,624,724 of genome coding rows and 1,820,698 / 2,552,394 of exome coding rows. Stars >= 2 are common, but a high star is not a pathogenic or benign label. The high-star directional classes are much smaller: genome 1,397 P/LP and 2,524 B/LB; exome 3,404 P/LP and 3,262 B/LB.

MaveDB allele scores reach 4,246 / 1,624,724 genome coding keys (137 score sets, 8,020 key-set rows) and 13,337 / 2,552,394 exome coding keys (152 score sets, 24,969 key-set rows). Most coding variants have no assay score. Across the whole scored file, 304,093 `(urn, key)` pairs had one numeric allele score, 397 pairs had conflicting allele scores and were dropped, and 22,892 keys had a haplotype score and no unique allele score. Those 22,892 are not given a MaveDB value in this report.

### Agreement, the set to trust first

On missense rows with both CADD and REVEL, the two calls agree here:

| side | both scored | REVEL >= 0.644 and CADD >= 20 | REVEL <= 0.290 and CADD < 10 |
| --- | ---: | ---: | ---: |
| genome | 1,029,589 | 100,823 | 166,732 |
| exome | 1,590,310 | 187,085 | 231,955 |

The same pathogenic-side agreement at stricter REVEL cuts, still with CADD >= 20: moderate REVEL, 54,220 genome and 105,708 exome; strong REVEL, 9,588 genome and 21,821 exome.

High-star labels on top of that agreement are the clinical set to trust first. Denominator: high-star missense with both scores.

| side | class | both scored | agrees with both calls |
| --- | --- | ---: | ---: |
| genome | P/LP | 688 | 607 are REVEL >= 0.644 and CADD >= 20 |
| genome | B/LB | 1,156 | 137 are REVEL <= 0.290 and CADD < 10 |
| exome | P/LP | 1,483 | 1,329 are REVEL >= 0.644 and CADD >= 20 |
| exome | B/LB | 1,244 | 144 are REVEL <= 0.290 and CADD < 10 |

Read the 607 and 1,329 first for pathogenic, and the 137 and 144 first for benign. The benign trust cell is a minority of high-star B/LB missense (137 / 1,156 and 144 / 1,244). Most high-star benign rows are in a middle band on at least one tool, which is not agreement and not the disagreement cell below.

Inside one MaveDB score set, the same rule is: trust the sets where high-star P/LP and high-star B/LB sit at opposite ends, and the variants that sit with their own class. Six exome sets and one genome set had at least 10 high-star P/LP and 10 high-star B/LB. In every one of them the P/LP median is the lower score. Five exome sets and the one genome set also meet the 25% wrong-side rule. Wrong-side counts are the interesting remainder; the rest of those sets is the agreement.

| side | score set | target | P/LP n | B/LB n | median P/LP | median B/LB | P/LP on the B/LB side of the B/LB median | B/LB on the P/LP side of the P/LP median |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| exome | `urn:mavedb:00001225-a-1` | BRCA2, Hap1 saturation editing | 48 | 137 | -0.873 | 0.134 | 0 | 0 |
| exome | `urn:mavedb:00001242-a-1` | BRCA2, humanized mouse ES saturation editing | 45 | 141 | -1.843 | -0.061 | 1 | 2 |
| exome | `urn:mavedb:00000097-0-2` | BRCA1 saturation editing | 37 | 76 | -2.050 | -0.025 | 0 | 1 |
| exome | `urn:mavedb:00000003-a-1` | BRCA1 RING, E3 assay | 12 | 59 | -2.737 | -0.481 | 0 | 1 |
| exome | `urn:mavedb:00001224-a-1` | BRCA2 homology-directed repair | 19 | 23 | 1.087 | 4.479 | 1 | 0 |
| genome | `urn:mavedb:00001225-a-1` | BRCA2, Hap1 saturation editing | 11 | 76 | -0.901 | 0.132 | 0 | 0 |

A lower score is the P/LP side only inside these sets. That direction is not exported to any other urn.

SpliceAI and Pangolin mostly agree at the low end. Among coding rows with both scores, both are under 0.2 on 1,556,775 / 1,610,506 genome rows and 2,439,220 / 2,533,300 exome rows, and both are at least 0.5 on 5,187 and 11,161. Those two corners are agreement. The absolute-value cut is an analysis choice, as noted above.

### Disagreement

**REVEL against CADD, missense, both scores present.** The rare opposing cell is REVEL >= 0.644 and CADD < 10: 80 genome, 135 exome. The large opposing cell is REVEL <= 0.290 and CADD >= 20: 304,861 / 1,029,589 genome and 460,648 / 1,590,310 exome.

That large cell is mostly the CADD cut. The median CADD among REVEL-benign missense is 18.58 on the genome side and 18.93 on the exome side, just under 20. Of REVEL-benign missense, 304,861 / 685,577 genome and 460,648 / 1,008,149 exome are still CADD >= 20. Raising the CADD cut shrinks the cell:

| side | REVEL <= 0.290 and CADD >= 20 | and CADD >= 25 | and CADD >= 30 |
| --- | ---: | ---: | ---: |
| genome | 304,861 | 60,420 | 4,754 |
| exome | 460,648 | 91,382 | 7,783 |

The disagreement that survives a top-0.1% CADD rank is those 4,754 and 7,783, plus the 80 and 135 where REVEL is pathogenic-supporting and CADD is under 10. CADD under 15 with REVEL >= 0.644 is 202 genome and 338 exome. CADD under 20 with REVEL >= 0.644 is 755 and 1,369, which includes the middle CADD band. The full 3-by-3, including bands that are not calls:

| side | REVEL call | CADD < 10 | CADD 10 to < 20 | CADD >= 20 |
| --- | --- | ---: | ---: | ---: |
| genome | <= 0.290 | 166,732 | 213,984 | 304,861 |
| genome | between | 5,620 | 22,142 | 214,672 |
| genome | >= 0.644 | 80 | 675 | 100,823 |
| exome | <= 0.290 | 231,955 | 315,546 | 460,648 |
| exome | between | 9,215 | 36,501 | 347,991 |
| exome | >= 0.644 | 135 | 1,234 | 187,085 |

Three genome examples of REVEL >= 0.644 and CADD < 10, the sharp rare cell: `chr9:128940377:C:T` PHYHD1 (CADD 8.633, REVEL 0.948, single-submitter benign), `chr12:6330660:G:C` TNFRSF1A (CADD 8.849, REVEL 0.820, two stars, not a P/LP or B/LB class), `chr17:43091733:G:T` BRCA1 (CADD 8.217, REVEL 0.806, conflicting ClinVar). The exome side repeats PHYHD1 and adds `chr15:40415470:C:G` IVD (CADD 6.412, REVEL 0.907, single-submitter P/LP) and `chr16:81259710:C:A` BCO1 (CADD 5.062, REVEL 0.850, single submitter).

**High-star labels against the missense scores.** There is no high-star missense P/LP with REVEL <= 0.290 and CADD < 10 (0 genome, 0 exome). High-star P/LP with REVEL <= 0.290 at any CADD is 8 / 688 genome and 18 / 1,483 exome. Of those, 7 genome and 15 exome are CADD >= 20. The other 1 and 3 sit in the band from 10 to 20, which is not a CADD call. None are CADD < 10. REVEL is the score that splits from the label.

Those rows are expert-panel (3-star) missense, and several are splice-aware rather than a clean protein-score failure. Lowest REVEL on the exome side:

| variant | gene | CADD | REVEL | SpliceAI | label |
| --- | --- | ---: | ---: | ---: | --- |
| `chr19:17840276:C:T` | JAK3 | 22.8 | 0.068 | 0.01 | likely pathogenic |
| `chr17:44389286:C:T` | ITGA2B | 25.4 | 0.089 | 0.61 | likely pathogenic, also splice-region |
| `chr17:44385570:A:C` | ITGA2B | 16.71 | 0.134 | 0.89 | likely pathogenic |
| `chr19:17834571:C:T` | JAK3 | 33.0 | 0.153 | 0.37 | pathogenic, also splice-region |
| `chr12:88109135:A:G` | CEP290 | 26.5 | 0.157 | 0.01 | likely pathogenic |

On the genome side the same JAK3 and CEP290 rows recur, and the list adds `chr11:22255485:C:G` ANO5 (CADD 33.0, REVEL 0.176, SpliceAI 1.0, pathogenic), `chr2:151696713:T:C` NEB (CADD 24.9, REVEL 0.182, SpliceAI 0.28), and `chr7:44146531:G:C` GCK (CADD 16.41, REVEL 0.191, SpliceAI 0.34). For ITGA2B and ANO5, SpliceAI agrees with the panel and REVEL does not. That is a disagreement between a missense score and a splice score, not a reason to demote the 3-star label.

The other clinical disagreement is larger. High-star B/LB missense with REVEL >= 0.644 and CADD >= 20: 153 / 1,156 genome and 164 / 1,244 exome. One more on each side is REVEL >= 0.644 and CADD < 10. Highest REVEL among the exome 3-star benign or likely benign rows: `chr19:38543551:A:G` RYR1 (CADD 32.0, REVEL 0.983, SpliceAI 0), `chr7:44153306:C:T` GCK (26.1, 0.983, 0), `chr2:165130164:C:T` SCN3A (27.0, 0.977, 0.13), `chr1:68439059:T:G` RPE65 (26.6, 0.974, 0), `chr12:120996569:C:A` HNF1A (27.2, 0.969, 0.02, benign). The genome side repeats the first four keys with the same scores and, instead of that HNF1A allele, has `chr12:120996568:C:G` HNF1A (CADD 25.0, REVEL 0.964, SpliceAI 0.02, benign). Here CADD and REVEL agree with each other and disagree with the expert panel. SpliceAI is low, so splicing is not the explanation in these five. The 3-star label is the one to trust first. The rows are the disagreement to inspect, not a reclassification.

**SpliceAI against REVEL, missense, both present.** A negative SpliceAI is not counted as disagreement. The opposing pair is SpliceAI >= 0.5 and REVEL <= 0.290: 2,266 genome and 4,139 exome. That is 2,266 / 4,161 of genome missense with SpliceAI >= 0.5, and 4,139 / 8,458 of the exome equivalent. About half of missense variants at the recommended splice point look benign to REVEL. SpliceAI >= 0.8 and REVEL <= 0.290 is 840 genome and 1,592 exome. Examples with SpliceAI 1.0 and REVEL under 0.04, all single-submitter and not high-star: exome `chr2:110666273:T:A` BUB1 (REVEL 0.022, CADD 17.29), `chr7:147639222:A:G` CNTNAP2 (0.028, 22.2), `chr8:89958731:T:C` NBN (0.028, 25.8); genome repeats BUB1 and adds `chr12:49599505:T:C` FAM186B (REVEL 0.030, CADD 14.54) and `chr7:50543899:T:C` DDC (0.039, 19.75).

Among high-star rows with a SpliceAI score, SpliceAI >= 0.5 and REVEL <= 0.290 is 1 of 1,397 genome P/LP, 3 of 2,523 genome B/LB, 5 of 3,402 exome P/LP, and 3 of 3,261 exome B/LB. High-star P/LP with SpliceAI < 0.2, REVEL <= 0.290, and CADD < 10: 0 on both sides. The panel-versus-all-three-tools cell is empty.

**Pangolin against SpliceAI**, coding rows with both scores, using the absolute Pangolin value. Opposing corners, not the whole discordant band: SpliceAI >= 0.5 and absolute Pangolin < 0.2 is 1,024 genome and 1,734 exome. Absolute Pangolin >= 0.5 and SpliceAI < 0.2 is 1,319 genome and 2,259 exome. Pangolin itself is negative on 455,581 genome and 704,943 exome of those both-scored rows, which is why the absolute value is required.

**MaveDB against CADD, inside one score set.** Cross-set comparison of the numbers was not done. A variant in two urns is not averaged. Sets with at least 25 coding variants at CADD >= 20 and 25 at CADD < 10, and unequal medians: 24 genome, 37 exome. Sets that also keep each wrong-side fraction at or under 0.25: 2 genome, 9 exome. The other sets overlap. Overlap, not a variant list, is the usual CADD-versus-assay result.

Where a set does separate, the direction is not shared, which is why the scores cannot be pooled. On exome coding rows in `urn:mavedb:00001213-a-1` (TP53 exons 5–8 in HCT116), the CADD >= 20 median is 0.832 and the CADD < 10 median is -1.018; 14 of 239 high-CADD variants sit at or below the low-CADD median, and 0 of 92 low-CADD variants sit at or above the high-CADD median. On `urn:mavedb:00000097-0-2` (BRCA1 saturation editing) the medians go the other way: -0.439 versus 0.056, with 51 of 294 and 14 of 111 on the wrong side. Two exome pairs of urns have identical counts and medians on this slice (`urn:mavedb:00001264-a-1` and `-a-2`, XRCC2; `urn:mavedb:00001263-a-1` and `-a-2`, BRCA2), so the nine exome sets are not nine independent assays. The per-urn table is in the SQL result `mave_cadd_sets`.

**MaveDB against high-star labels, the set that does not separate.** `urn:mavedb:00000003-b-1` (BRCA1 RING, yeast two-hybrid) has 12 high-star P/LP and 59 high-star B/LB on the exome coding slice. Medians are -0.665 and 0.060, so the P/LP median is still lower, but 5 of 12 P/LP variants sit at or below the B/LB median. That fraction is above the 0.25 rule. Those five are the within-set disagreement. The E3 assay on the same domain (`urn:mavedb:00000003-a-1`) does separate. Same gene, different score set, different relationship to the labels. The scores were not compared across those two urns.

## Caveat

CADD >= 20 is a rank among all substitutions, and missense variants pile up above it. Treating that line as a pathogenic call makes REVEL look wrong on hundreds of thousands of benign-supporting missense variants. The median CADD in that REVEL bin is under 20. The cell that remains at CADD >= 30 is the computational disagreement to look at. The empty cell matters too: no high-star missense P/LP is REVEL-benign and CADD-low together.

REVEL's supporting cuts are a calibration for missense evidence strength, not a definition of truth. They are applied here only to missense rows. About 15% of missense coding rows have no REVEL score, and those rows are outside every REVEL rate above. SpliceAI >= 0.5 with a benign REVEL is a mechanism split (splice versus protein), not proof that one tool failed. Canonical splice acceptor and donor sites were excluded by the coding-term rule, so this is not a splice-site report.

ClinVar star 1 mixes nothing in this table: single submitter and conflicting classifications were split, and they did not co-occur. Most coding rows are single-submitter. A two-star VUS is still not a directional label. High-star counts above are germline only, with oncogenicity calls removed. Somatic and drug-response rows are in `other` and are not used as either side of a disagreement.

MaveDB scores are not on one scale. No cross-set comparison was done. A higher score is not pathogenic. The lower-score P/LP pattern is reported only for the seven sets that had enough high-star labels of both classes, and it is not a license to flip the sign on the TP53 set or on any untested urn. Haplotype rows were not used. An assay covers a few genes; 13,337 exome coding keys with an allele score is not a genome-wide functional map.

Genome and exome copies agreed wherever both existed, and they were still not added. A rate in one table is not a rate for the other callset. SIFT and PolyPhen were not turned into calls. The SQL that produced these counts is `warehouse/reports/sql/predictor_disagreement.sql`.
