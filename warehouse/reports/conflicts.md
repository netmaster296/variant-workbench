# Conflicts

Question: where do ClinVar, gnomAD, the GWAS Catalog, and MaveDB disagree?

These rows audit the warehouse as much as they describe biology. A count is kept only when the join columns agree. Genome allele frequency is `AF`. Exome allele frequency is `AF_exome`. Those two columns are never added. The population frequency columns are not used, and ancestry groups are not ranked.

The statements that produced the numbers are in `warehouse/reports/sql/conflicts.sql`.

## Germline pathogenic alleles with a high genome frequency

Question: how many germline ClinVar pathogenic or likely pathogenic alleles are common in gnomAD genomes?

Filters:

- `germline_pathogenic` is true on `warehouse/clinvar/allele_dedup.parquet`. That flag is already limited to germline origin, with no oncogenicity call, and it requires a pathogenic or likely pathogenic piece and no conflicting piece. Likely pathogenic is included. Low-penetrance text is included when the flag is true.
- The allele has a `variant_key`, and that key is in `warehouse/variants/tier_a.parquet`.
- `in_genomes` is true and `AF` is not null.
- Three thresholds, all on genome `AF`: at or above 0.05, at or above 0.01, and at or above 0.001.
- `AF_exome` is reported beside `AF` and is not part of the threshold.
- The cuts are descriptive frequency cuts. They are not an ACMG benign stand-alone classification.

Denominator: 333,246 keyed germline pathogenic or likely pathogenic alleles. Another 5,341 flagged alleles have no `variant_key` and are not in this join. All 333,246 keyed alleles join to Tier A, and on every join the chromosome, position, ref, and alt agree. 295,503 are not in the genome sites file, so they have no genome AF. Absence is not an AF of zero. One genome row has a null `AF`. The measured denominator is the other 37,742.

Result:

| Genome `AF` | Alleles | Share of 37,742 | `PASS` and stars at least 2 |
| --- | ---: | ---: | ---: |
| at or above 0.05 | 8 | 0.000212 | 0 |
| at or above 0.01 | 27 | 0.000715 | 12 |
| at or above 0.001 | 187 | 0.004954 | 139 |

`AF` matches `AC / AN` within 0.0001 on all 37,742 measured rows. No row was dropped for an allele mismatch.

Of the 27 alleles at or above 0.01, 18 have `filter = PASS`, 16 have `review_stars` of at least 2, and 12 have both. Six of the 27 already carry "low penetrance", "risk factor", or "susceptibility" in `ClinicalSignificance`. Five of those six are inside the 12. The other seven of the 12 do not: `chr7:117548628:GTT:G` (CFTR), `chrX:154536002:C:T` (G6PD), `chr11:5227002:T:A` (HBB), `chr1:152312600:CACTG:C` (FLG), `chr1:152313385:G:A` (FLG), `chr3:165830741:T:C` (BCHE), and `chr6:32040110:G:T` (CYP21A2).

Three alleles at or above 0.01 have a genome `AF` at least ten times the exome `AF`: `chr10:48066:C:A` (TUBB8), `chr19:1037767:G:A` (CNN2), and `chr19:1037757:G:A` (CNN2). All three have `review_stars` 0 and a filter other than `PASS`. Two alleles in the same list are not in the exomes at all (`chrX:2748343:G:C` and `chr11:75561608:C:T`). A missing exome row is not an exome AF of zero.

The eight alleles at or above 0.05 are the first eight rows below. None of them are both `PASS` and two or more stars. The two `PASS` rows in that bin have zero stars: a ClinVar gene symbol of `-` whose phenotype text is an Xg blood-group allele, and an NCR3 row whose significance text already says "risk factor".

| variant_key | GeneSymbol | stars | filter | genome AF | exome AF | ClinicalSignificance | consequence |
| --- | --- | ---: | --- | ---: | ---: | --- | --- |
| chr7:142750561:C:T | PRSS1 | 2 | InbreedingCoeff | 0.3862 | 0.202123 | Pathogenic/Likely pathogenic | missense_variant |
| chrX:2748343:G:C | - | 0 | PASS | 0.354497 |  | Pathogenic | upstream_gene_variant |
| chr4:102901325:AAC:A | SLC9B1 | 2 | AS_VQSR;InbreedingCoeff | 0.341192 | 0.10245 | Pathogenic | frameshift_variant |
| chr8:71362155:C:CT | EYA1 | 1 | AS_VQSR | 0.190688 | 0.15792 | Likely pathogenic | upstream_gene_variant |
| chr6:31593133:C:G | NCR3 | 0 | PASS | 0.129762 | 0.0792867 | Pathogenic; risk factor | upstream_gene_variant |
| chr8:99832368:G:T | VPS13B | 2 | AS_VQSR | 0.0809806 | 0.110904 | Pathogenic | splice_acceptor_variant |
| chr10:48066:C:A | TUBB8 | 0 | AS_VQSR;InbreedingCoeff | 0.067775 | 0.000111699 | Likely pathogenic | missense_variant |
| chr19:1037767:G:A | CNN2 | 0 | InbreedingCoeff | 0.0511315 | 0.000109274 | Pathogenic | missense_variant |
| chr19:1037757:G:A | CNN2 | 0 | AS_VQSR;InbreedingCoeff | 0.0468217 | 0.000190874 | Pathogenic | missense_variant |
| chr12:11267400:A:AG | PRB3 | 0 | AS_VQSR | 0.0413653 | 0.0144371 | Pathogenic | frameshift_variant |
| chr6:26092913:G:A | HFE | 2 | PASS | 0.0372927 | 0.0590333 | Pathogenic/Pathogenic, low penetrance; risk factor | missense_variant |
| chr7:117548628:GTT:G | CFTR | 2 | PASS | 0.0356948 | 0.0286172 | Pathogenic/Likely pathogenic | splice_polypyrimidine_tract_variant&intron_variant |
| chrX:154536002:C:T | G6PD | 2 | PASS | 0.0334892 | 0.00368199 | Pathogenic/Likely pathogenic | missense_variant |
| chr8:93755776:A:T | TMEM67 | 2 | AS_VQSR | 0.0331748 | 0.00782628 | Pathogenic/Likely pathogenic | splice_acceptor_variant |
| chrX:106036488:C:T | SERPINA7 | 0 | PASS | 0.0314199 | 0.00344396 | Pathogenic | missense_variant |
| chr14:94380925:T:A | SERPINA1 | 2 | PASS | 0.028729 | 0.0371547 | Pathogenic/Pathogenic, low penetrance; other | missense_variant |
| chr6:31592893:C:T | NCR3 | 0 | PASS | 0.0232547 | 0.0211588 | Pathogenic | 5_prime_UTR_variant |
| chr11:75561608:C:T | SERPINH1 | 0 | PASS | 0.0220249 |  | Pathogenic | upstream_gene_variant |
| chr1:145927447:C:T | RBM8A | 2 | PASS | 0.0184438 | 0.0260167 | Pathogenic/Likely pathogenic/Pathogenic, low penetrance; other | 5_prime_UTR_variant |
| chr2:218890289:T:A | WNT10A | 2 | PASS | 0.0140464 | 0.0214924 | Pathogenic/Likely pathogenic/Pathogenic, low penetrance | missense_variant |
| chr11:27368192:T:C | LGR4 | 0 | PASS | 0.013546 | 0.0203732 | Pathogenic | missense_variant |
| chr11:5227002:T:A | HBB;LOC106099062;LOC107133510 | 2 | PASS | 0.0127188 | 0.00160112 | Pathogenic | missense_variant |
| chr1:152312600:CACTG:C | FLG | 2 | PASS | 0.0125877 | 0.0196365 | Pathogenic/Likely pathogenic | frameshift_variant |
| chr3:165830741:T:C | BCHE | 2 | PASS | 0.0124371 | 0.0163971 | Pathogenic/Likely pathogenic | missense_variant |
| chr14:94378610:C:T | SERPINA1 | 2 | PASS | 0.0121565 | 0.016249 | Pathogenic/Likely pathogenic; risk factor | missense_variant |
| chr6:32040110:G:T | CYP21A2 | 2 | PASS | 0.0108544 | 0.00438799 | Pathogenic/Likely pathogenic | missense_variant |
| chr1:152313385:G:A | FLG | 2 | PASS | 0.0104279 | 0.0176366 | Pathogenic/Likely pathogenic | stop_gained |

Frequencies in that table are the stored doubles as DuckDB prints them. Empty exome AF means `in_exomes` is false. The HBB row's `GeneSymbol` is the dedup value, three symbols joined by semicolons.

Caveat: most of the extreme tail is a review or filter problem, not a high-confidence pathogenic allele that is also common on both frequency columns. The alleles that remain after `PASS` and two stars, from HFE downward, are real stored disagreements between a pathogenic label and a genome AF of at least 0.01. Several of those labels already say low penetrance or risk factor. The CFTR row is an intron polypyrimidine-tract allele, not a coding missense. HBB's phenotype text includes "HEMOGLOBIN S"; its genome AF is 0.0127 and its exome AF is 0.0016, about eightfold, which is under the tenfold cut used above. Do not read 8 out of 333,246. The 295,503 alleles with no genome row are not known to be rare.

## GWAS risk alleles that are the common allele

Question: after a ref/alt key exists, how often is the reported risk allele the common gnomAD allele?

Filters:

- The association has an rsid, a canonical chromosome, a position, and one `risk_allele` token. Rows with a semicolon in `risk_allele` are out. There are 755,765 such association rows and 320,040 distinct rsids.
- `associations.variant_key` is null on all 1,205,644 catalog rows. It is not the join key.
- The key is a ClinVar `variant_key` with the same rsid, chromosome, and position, whose `ref` or `alt` equals the risk-allele token.
- A claim is one `(rsid, risk allele)` pair. If that pair matches more than one ClinVar key, it is ambiguous and is dropped. That is 207 pairs.
- The common-allele call uses only the case where the risk allele equals `alt` and does not equal `ref`. Genome `AF` is the frequency of that alt. `AN` must be at least 1,000.
- Two thresholds: genome `AF` greater than 0.5, and genome `AF` greater than 0.9.
- The reference-allele frequency is not a stored column. `1 - AF` is not used, because another alternate at the same site would make that quantity larger than the reference frequency. Risk alleles that match `ref` are counted and then not called common or rare.

Denominator: 7,661 unambiguous claims whose risk allele is the alt, with a measured genome `AF` and `AN` of at least 1,000. One further alt claim has `AN` below 1,000 and is outside this denominator. Behind the 1,139 claims above 0.5 sit 6,337 association rows.

Result:

| Risk allele is the alt, genome `AF` | Claims | Share of 7,661 |
| --- | ---: | ---: |
| greater than 0.9 | 48 | 0.00627 |
| greater than 0.5 | 1,139 | 0.1487 |
| 0.05 to 0.5, inclusive | 3,955 | 0.5162 |
| below 0.05 | 2,567 | 0.3351 |

Of the 1,139 claims with genome `AF` above 0.5:

| Catalog median `RISK ALLELE FREQUENCY` | Claims |
| --- | ---: |
| also greater than 0.5 | 701 |
| 0.5 or below | 176 |
| not a plain number | 262 |

The 701 are cases where the catalog's own frequency and the genome `AF` both say the risk allele is common. The 176 are the frequency disagreement. Median RAF on these claims runs from 0.027 to 0.992. No claim has a median RAF above 1, so the split is not a percent-versus-fraction parse. One of the 1,139 has at least one study RAF above 1.

The match set, before the alt restriction, is 11,584 `(rsid, risk)` pairs that hit a ClinVar key on rsid, position, and allele letter. 11,377 hit exactly one key. Of those, 7,684 are alt and not ref, 3,666 are ref and not alt, and 27 have `ref` equal to `alt`, so the letter matches both columns and the claim is unused.

Caveat: a risk allele with frequency above 0.5 is often the coded effect allele for a quantitative trait, not a rare-disease contradiction. Where the catalog RAF is also high, the catalogs agree. The claim that was not made is the reference-matched half: 3,666 claims name the reference base as the risk allele, and this table cannot prove that base is the major allele. Ambiguous rsids were not assigned a key.

## High missense-z genes and MaveDB scores

Question: among genes with a high missense z score, are the MaveDB scores mostly mild?

Filters:

- Protein-coding gene mart, `warehouse/marts/gene_mart_protein_coding.parquet`. High means `mis.z_score` at or above 3.09. That threshold is applied here. It is not a column in the file.
- The gene's symbol maps to one GENCODE gene (`symbol_n_gencode = 1`). The MaveDB join is exact `target_gene` = `symbol`.
- Scores are the `score` column inside one `score_set_urn`. `raw_score` is not used. Score sets are not pooled, and a score is not compared with a score from another set. The sets are not on one scale.
- A variant counts once per score set (`urn`). It must have one numeric `score`.
- `protein_change` is classified only when it is a single three-letter HGVS protein change. `Ter` is nonsense, not missense. `=` or the same amino acid repeated is synonymous. A different amino acid is missense.
- A score set can be called only when it has at least 10 synonymous, 10 nonsense, and 20 missense variants under that parse.
- It is oriented only when the nonsense median is outside the synonymous 10th to 90th percentile and that percentile width is smaller than the gap between the two medians. Otherwise the synonymous band already covers the nonsense scores, and "mild" would be circular.
- A missense score is mild when it falls inside that synonymous band. Mostly mild means more than half of the missense variants in that set.

Denominator: 20,097 protein-coding genes. 18,548 have a `mis.z_score`. 1,550 gene rows are at or above 3.09. Four of those rows share a symbol with a second gene and are excluded: HERC3 (`ENSG00000138641`, `ENSG00000287542`) and PDE8B (`ENSG00000113231`, `ENSG00000284762`). The symbol denominator is 1,546 genes. 37 of them have `mavedb_experiments_by_symbol` greater than 0. Nine of those 37 have no row in `mapped_variants_scored.parquet` (the mart's mapped-variant count is already 0): BAP1, FYN, HECTD1, HNRNPUL1, ITCH, PTPN11, RAF1, SHOC2, SPTAN1. 28 genes have at least one scored row. 13 of those 28 have no `protein_change` on any row, so missense cannot be separated from synonymous. The remaining genes are not pooled into one fraction.

Result: no oriented score set is mostly mild.

Fourteen score sets in five genes meet the count floor. Eight of the fourteen are oriented. In those eight, the share of missense scores inside the synonymous band is between 0.286 and 0.464.

| Gene | z | Score set | Title | Missense | Synonymous | Nonsense | Inside synonymous band | Oriented |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| HMGCR | 4.235 | urn:mavedb:00000035-a-1 | HMGCR rosuvastatin imputed and refined | 16,853 | 748 | 847 | 16,192 (0.961) | no |
| HMGCR | 4.235 | urn:mavedb:00000035-a-2 | HMGCR no statin imputed and refined | 16,853 | 748 | 847 | 16,082 (0.954) | no |
| HMGCR | 4.235 | urn:mavedb:00000035-a-3 | HMGCR atorvastatin imputed and refined | 16,853 | 748 | 847 | 16,216 (0.962) | no |
| KRAS | 3.090 | urn:mavedb:00000115-a-8 | binding fitness from KRAS-DARPin K55 bindingPCA of KRAS block3 | 4,675 | 10 | 37 | 1,889 (0.404) | yes |
| KRAS | 3.090 | urn:mavedb:00000115-a-24 | binding fitness from KRAS-RAF1RDB bindingPCA of KRAS block3 | 4,833 | 18 | 39 | 1,768 (0.366) | yes |
| KRAS | 3.090 | urn:mavedb:00000115-a-28 | abundance fitness from abundancePCA of KRAS block3 | 4,404 | 14 | 36 | 2,895 (0.657) | no |
| PTEN | 4.120 | urn:mavedb:00000013-a-1 | PTEN VAMP-seq | 4,112 | 156 | 140 | 1,847 (0.449) | yes |
| PTEN | 4.120 | urn:mavedb:00000101-a-1 | PTEN VAMP-seq Fill-in | 4,185 | 162 | 167 | 3,279 (0.784) | no |
| PTEN | 4.120 | urn:mavedb:00000102-0-1 | PTEN VAMP-seq Combined | 5,083 | 193 | 201 | 2,356 (0.464) | yes |
| PTEN | 4.120 | urn:mavedb:00001244-a-1 | PTEN Variant Imaging in Neurons | 960 | 48 | 45 | 349 (0.364) | yes |
| PTEN | 4.120 | urn:mavedb:00001244-b-1 | PTEN VIS-seq in hiPSCs | 960 | 48 | 45 | 275 (0.286) | yes |
| SCN5A | 4.845 | urn:mavedb:00000098-a-1 | SCN5A response to veratridine+brevetoxin+ouabain mixture | 224 | 12 | 12 | 93 (0.415) | yes |
| SCN5A | 4.845 | urn:mavedb:00000098-a-2 | SCN5A response to veratridine+brevetoxin+ouabain mixture | 224 | 12 | 12 | 93 (0.415) | yes |
| SPOP | 5.538 | urn:mavedb:00001258-a-1 | Selection of loss-of-function SPOP via yeast survival assay | 6,261 | 354 | 344 | 3,772 (0.602) | no |

z is `mis.z_score`, rounded. The two SCN5A sets have the same title and the same summary counts. They are not added.

HMGCR is the set that would have looked "mostly mild" (about 0.96) if the synonymous band were trusted. It is not. The synonymous 10th percentile is 0, and the nonsense median is inside the band (about 0.06 to 0.09). Synonymous and nonsense scores overlap, including scores of 0 on both. SPOP, the PTEN fill-in set, and KRAS `00000115-a-28` do put the nonsense median outside the synonymous band, but the synonymous band is wider than the gap between the medians, so the separation is too weak to call. SPOP's nonsense median is on the high side of its synonymous median. That direction is left as a property of that score set. It is not applied to any other set.

The oriented PTEN sets are the cleanest comparison: synonymous medians sit near 1, or near 0.2 on the two imaging sets, and the nonsense medians sit on the other side of a narrower band. Fewer than half of the missense scores fall in the synonymous band. That agrees with a high missense z. It is not a conflict. KRAS `00000115-a-8` rests on only 10 synonymous variants, so its band is thin. It is still not mostly mild. KRAS has many more score sets, but the others do not contain 10 synonymous and 10 nonsense in the same set.

Caveat: a genome-wide statement about MaveDB and missense constraint is not available from these files. Most high-z genes have no assay, and most assays that exist cannot be called mild or severe without pretending the scores share a scale. The gene-level answer is that none of the score sets that can be oriented are mostly mild. Variant counts above are assayed variants inside that score set, not a sum across sets.

## Conflicting ClinVar labels with a high REVEL or SpliceAI score

Question: how often does a conflicting ClinVar label sit on a variant with a high REVEL or SpliceAI score?

Filters:

- `ClinicalSignificance` contains "conflict". That is ClinVar's conflicting label, including "Conflicting classifications of pathogenicity".
- The allele has a `variant_key`. Every such key is in Tier A.
- High REVEL is genome `revel_max` at or above 0.75. The same column is also counted at 0.5.
- High SpliceAI is genome `spliceai_ds_max` at or above 0.5. The same column is also counted at 0.2 and at 0.8.
- `revel_max_exome` and `spliceai_ds_max_exome` are not copied into the genome columns and are not used in the high count. They are counted only where the genome score is null.

Denominator: 183,646 conflicting alleles with a `variant_key`. A further 168 conflicting rows have no key. 127,850 of the keyed alleles are in the genome sites file. Genome `revel_max` is present for 69,836. Genome `spliceai_ds_max` is present for 126,271. At least one of those two genome scores is present for 126,674.

Result:

| Score | Threshold | Alleles |
| --- | --- | ---: |
| `revel_max` | at or above 0.50 | 16,574 |
| `revel_max` | at or above 0.75 | 6,430 |
| `spliceai_ds_max` | at or above 0.20 | 5,491 |
| `spliceai_ds_max` | at or above 0.50 | 2,405 |
| `spliceai_ds_max` | at or above 0.80 | 1,371 |
| either high cut | `revel_max` >= 0.75 or `spliceai_ds_max` >= 0.50 | 8,786 |
| both high cuts | both of those | 49 |

6,430 / 69,836 is 0.092 of conflicting alleles that have a genome REVEL score. 2,405 / 126,271 is 0.019 of those that have a genome SpliceAI score. 8,786 / 126,674 is 0.069 of conflicting alleles with either genome score.

Where the genome score is null, `revel_max_exome` is at or above 0.75 for 3,871 alleles, and `spliceai_ds_max_exome` is at or above 0.5 for 1,513. Those counts are not added to the 6,430 or the 2,405. Where both datasets have a REVEL score, none differ by more than 0.05. The same is true for SpliceAI.

Caveat: a conflicting label with a high score is a disagreement between submitters, recorded by ClinVar, next to a predictor that is high. It is not evidence that the predictor columns were joined to the wrong allele. The sampled high-REVEL rows below are rare on genome `AF` and have the same score on the exome columns. REVEL is missing for 58,014 of the 127,850 conflicting alleles that are in the genome file, which is expected for variants the site VCF did not score. gnomAD does not put every predictor on every site.

## What a re-read of the source columns showed

The checks below are a second read of the stored columns, not a new join rule.

Germline frequency join. For all 333,246 keyed germline pathogenic alleles, Tier A ref and alt equal the dedup ref and alt. For all 37,742 with a genome AF, `AF` is within 0.0001 of `AC / AN`.

`chr6:26092913:G:A`, HFE, rs1800562. Dedup: `ClinicalSignificance` is "Pathogenic/Pathogenic, low penetrance; risk factor", `review_stars` 2, `germline_pathogenic` true, ref G, alt A. Tier A: ref G, alt A, `filter` PASS, `AC` 5,679, `AN` 152,282, `AF` 0.0372927, which is 5679/152282. `AF_exome` is 0.0590333. The label already qualifies penetrance. The two frequency columns agree that the allele is not rare. This is not a bad join.

`chr11:5227002:T:A`, HBB, rs334. Dedup: `ClinicalSignificance` "Pathogenic", stars 2, `germline_pathogenic` true, ref T, alt A. The phenotype text includes "HEMOGLOBIN S". Tier A: `filter` PASS, `AC` 1,937, `AN` 152,294, `AF` 0.0127188, equal to 1937/152294. `AF_exome` is 0.00160112. The pathogenic label and a genome AF above 0.01 are both stored. The exome column is lower and was not folded in.

`chr10:48066:C:A`, TUBB8. Dedup: "Likely pathogenic", stars 0, `germline_pathogenic` true, no rsid. Tier A: `filter` "AS_VQSR;InbreedingCoeff", `AC` 5,104, `AN` 75,308, `AF` 0.067775, equal to 5104/75308. `AF_exome` is 0.000111699. The genome and exome frequencies disagree by about 600-fold, and the filter is not PASS. Averaging them would hide that.

`chr17:50600964:A:T`, rs198530. The association row has `variant_key` null, `risk_allele` T, `STRONGEST SNP-RISK ALLELE` "rs198530-T", `RISK ALLELE FREQUENCY` "0.992059662", chromosome chr17, position 50600964, trait text "WFIKKN2 protein levels". ClinVar has the same rsid, chromosome, and position, ref A, alt T, key `chr17:50600964:A:T`, `GeneSymbol` CACNA1G, `ClinicalSignificance` Benign, stars 2. Tier A: `AC` 151,520, `AN` 152,256, `AF` 0.995166, equal to 151520/152256, `filter` PASS, `AF_exome` 0.993302. The risk allele is the alt. The catalog RAF, the genome AF, and the exome AF all say it is common, and the ClinVar label on that key is Benign.

`chr14:54902784:A:G`, rs1753589. The association row has `variant_key` null, `risk_allele` G, raw allele "rs1753589-G", `RISK ALLELE FREQUENCY` "0.027", chromosome chr14, position 54902784. ClinVar has the same rsid and position, ref A, alt G, `GeneSymbol` GCH1, `ClinicalSignificance` Benign, stars 2, `germline_pathogenic` false. Tier A genome `AF` is 0.959873 (`AC` 146,085, `AN` 152,192), and `AF_exome` is 0.990208. The letters and the position match. The frequency columns do not: the catalog cell is 0.027 and both gnomAD frequency columns are above 0.95. That is one of the 176. The join did not assign the wrong base. The catalog frequency is a study field, not a gnomAD field, and this report does not decide which number a reader should trust.

`chr10:13288508:T:C`, PHYH. Dedup: "Conflicting classifications of pathogenicity", stars 1, `germline_pathogenic` false. Tier A: `revel_max` 0.997, `revel_max_exome` 0.997, `spliceai_ds_max` 0.04, `spliceai_ds_max_exome` 0.04, `filter` PASS, genome `AF` about 3.9e-5 (`AC` 6, `AN` 152,134). The conflicting label is the ClinVar text. The predictor columns agree across the two datasets, and the allele is rare. The same pattern holds for the other top `revel_max` rows that were re-read with this one (TYR, COL4A3, CYP27A1, COL5A2): conflict text, stars 1, `revel_max` equal to `revel_max_exome` at 0.994 to 0.997, and a small genome AF.

MaveDB, PTEN score set `urn:mavedb:00000013-a-1`, re-read from `mapped_variants_scored.parquet`:

| urn | protein_change | score |
| --- | --- | ---: |
| urn:mavedb:00000013-a-1#1821 | NP_000305.3:p.Pro231= | 1.064 |
| urn:mavedb:00000013-a-1#1880 | NP_000305.3:p.Asp236= | 0.942 |
| urn:mavedb:00000013-a-1#1071 | NP_000305.3:p.Val166Ter | -0.038 |
| urn:mavedb:00000013-a-1#1185 | NP_000305.3:p.Tyr177Ter | 0.512 |
| urn:mavedb:00000013-a-1#1756 | NP_000305.3:p.Ser227Asn | -0.223 |

The two synonymous scores sit inside this set's synonymous band (0.775 to 1.145). The missense score does not, so it is not called mild. One nonsense score is near zero and one is 0.512. The call uses the nonsense median of this set, 0.064, which is outside the band. It does not assume every nonsense score is low, and it does not compare 0.512 with a score from another set.

HMGCR score set `urn:mavedb:00000035-a-1`: `urn:mavedb:00000035-a-1#1859` is `NP_000850.1:p.Ser799=` with `score` 0, and `urn:mavedb:00000035-a-1#210` is `NP_000850.1:p.Met867Ter` with `score` 0. A synonymous change and a nonsense change share a score. That is why the high fraction of HMGCR missense scores inside the synonymous band is not reported as "mostly mild".

No join in these four comparisons was broken. Ambiguous GWAS keys, reference-only risk alleles, and MaveDB sets without an internal synonymous and nonsense scale were left uncalled instead of being patched.
