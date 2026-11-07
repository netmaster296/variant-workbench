# ClinVar against frequency

Frequencies are gnomAD v4.1 global allele frequencies on `warehouse/variants/tier_a.parquet`. Genome AF is the column `AF`. Exome AF is `AF_exome`. The two are never added, averaged, or otherwise combined. An allele that is in both call sets is counted once in the genome table and once in the exome table.

The SQL is `warehouse/reports/sql/clinvar_versus_frequency.sql`.

## Question

Are pathogenic and likely pathogenic alleles mostly rare?

## Filters

Alleles are rows in `warehouse/clinvar/allele_dedup.parquet`: one row per `variant_key`, or one row per `AlleleID` when the record has no VCF allele. Phenotypes are already packed into that row. This is not a count of submissions.

The class is the first match on the lowered `ClinicalSignificance` text:

- `conflicting` if the text contains `conflict`. This wins even when the word pathogenic is also present.
- `pathogenic` if `pathogenic` is still present after the phrase `likely pathogenic` is removed. `Pathogenic/Likely pathogenic` is this class, not both classes.
- `likely_pathogenic` if the text contains `likely pathogenic` and the row was not already classified.
- `benign` if `benign` remains after `likely benign` is removed. `Benign/Likely benign` is this class.
- `likely_benign`, then `vus` (`uncertain significance` or `vus`), then `drug_response`, then `risk_or_other_effect` (risk factor, risk allele, association, protective, affects, or confers sensitivity), then `no_classification` (blank, `-`, not provided, no classification, no assertion, no interpretation), then `other`.

Low penetrance stays inside pathogenic or likely pathogenic. A pathogenic row that also says drug response or risk factor stays pathogenic. Those words are flags on the row, not a reason to move it.

`germline_pathogenic` is the loaded flag. It is not recomputed here. It is the germline-origin subset of the pathogenic and likely pathogenic classes, with no oncogenicity call and no conflicting piece.

AF bins are the same for both datasets and are half-open: null, `<0.0001`, `0.0001–<0.001`, `0.001–<0.01`, `0.01–<0.05`, `>=0.05`. A stored AF of 0 is inside `<0.0001`. It is not a null. Bins use the stored AF even when `filter` is not `PASS`. The exception section separates `PASS`.

The common-allele check is every pathogenic or likely pathogenic allele with genome AF `>= 0.01` or exome AF `>= 0.01`. That is the whole tail, 34 alleles, not a draw. Names are `Name` from `warehouse/clinvar/alleles.parquet`.

Ancestry-specific AF columns are not used.

## Denominator

Pathogenic or likely pathogenic alleles, excluding conflicting labels: 376,350. Of those, 243,482 are pathogenic and 132,868 are likely pathogenic. The germline flag is 338,587. Text pathogenic or likely pathogenic outside that flag: 37,763.

| set | genome AF non-null | genome AF null | exome AF non-null | exome AF null | alleles |
| --- | ---: | ---: | ---: | ---: | ---: |
| pathogenic | 28365 | 215117 | 70884 | 172598 | 243482 |
| likely pathogenic | 10545 | 122323 | 34630 | 98238 | 132868 |
| pathogenic or likely pathogenic | 38910 | 337440 | 105514 | 270836 | 376350 |
| germline flag | 37742 | 300845 | 100336 | 238251 | 338587 |
| text P/LP, flag false | 1168 | 36595 | 5178 | 32585 | 37763 |

Among the 376,350, 33,204 have both AFs, 5,706 have only a genome AF, 72,310 have only an exome AF, and 265,130 have neither. Those are presence counts. The frequencies are not combined.

Why the genome AF is null for the 337,440: 9,227 have no `variant_key`, 328,212 have a key but are not in the genome call set, and 1 is in the genome call set with a null AF. Why the exome AF is null for the 270,836: the same 9,227 have no key, 261,608 have a key but are not in the exome call set, and 1 is in the exome call set with a null AF. The two in-call-set nulls are AC 0 and AN 0: `CBS` `chr21:43063074:A:G` (genome filter `AC0;AS_VQSR`; exome AF is a measured 0) and `HGF` `chr7:81755200:TTCA:T` (exome filter `AC0;AS_VQSR`; genome AF is 0.000007, AC 1, AN 152,112, `PASS`).

Of the 337,440 with a null genome AF, 171,369 are single nucleotide variants and 166,071 are not (deletions 95,847, duplications 37,974, microsatellites 10,712, indels 9,071, insertions 7,090, copy number 5,219, other types 158). Of the 270,836 with a null exome AF, 124,805 are single nucleotide variants. A null is not an AF of zero. Separately, 787 pathogenic or likely pathogenic alleles have genome AF exactly 0, and 12,052 have exome AF exactly 0. Each of those has AC 0 and a filter containing `AC0`, and each is inside `<0.0001`.

chrM accounts for 191 of the null genome AFs. There is no chrM sites VCF.

## Result

Yes, among alleles that have a frequency. It is not a statement about the nulls.

Of pathogenic or likely pathogenic alleles with a non-null genome AF, 37,179 of 38,910 are below 0.0001 (0.955513), 38,720 of 38,910 are below 0.001 (0.995117), and 38,883 of 38,910 are below 0.01 (0.999306). Twenty-seven are at or above 0.01 (0.000694). Eight of those are at or above 0.05.

Of those with a non-null exome AF, 103,680 of 105,514 are below 0.0001 (0.982618), 105,302 of 105,514 are below 0.001 (0.997991), and 105,490 of 105,514 are below 0.01 (0.999773). Twenty-four are at or above 0.01 (0.000227). Six of those are at or above 0.05.

The germline flag is the same shape. Genome: 37,715 of 37,742 below 0.01 (0.999285), and the same 27 at or above 0.01. Exome: 100,314 of 100,336 below 0.01 (0.999781), and 22 at or above 0.01. The two text P/LP alleles at or above 0.01 that the flag excludes are both somatic, and both cross 0.01 on the exome AF only.

Split by class, still among non-null AFs:

| set | dataset | non-null | below 0.001 | below 0.01 | >= 0.01 | >= 0.05 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| pathogenic | genome | 28365 | 28200 | 28340 | 25 | 6 |
| pathogenic | exome | 70884 | 70702 | 70863 | 21 | 5 |
| likely pathogenic | genome | 10545 | 10520 | 10543 | 2 | 2 |
| likely pathogenic | exome | 34630 | 34600 | 34627 | 3 | 1 |

No likely pathogenic allele has a genome AF in `0.01–<0.05`. Both likely pathogenic genome AFs at or above 0.01 are at or above 0.05.

The same bins for the other classes show that this is not true of every label. Among non-null genome AFs, 134,170 of 249,057 benign alleles are at or above 0.01 (0.538712), 175 of 440 drug-response alleles are (0.397727), and 318 of 563 risk-or-other-effect alleles are (0.564831). Likely benign is 14,419 of 531,644 (0.027122). VUS is 984 of 1,068,153 (0.000921). Conflicting is 655 of 127,848 (0.005123). Drug response and risk labels were not pulled into the pathogenic bins. Conflicting text that also says pathogenic is 183,682 rows, and those rows stay in conflicting.

Every one of the 27 genome AFs and 24 exome AFs at or above 0.01 has AN of at least 41,628 (genomes) or 139,418 (exomes). None of those frequencies is a handful of alleles over a tiny denominator.

### What the 34 alleles at or above 0.01 actually were

Thirty-four alleles have a genome AF or an exome AF at or above 0.01 (27 genome, 24 exome, 17 both). Nine have either AF at or above 0.05. The name join returned 34 rows and 34 keys, so it did not multiply. All 34 have `n_rows` 1. All 34 have the same ref and alt on the dedup row and on `tier_a`. None has drug response in the significance. None has conflicting significance, and none has a review status of conflicting classifications. Thirty-two are germline. Two are somatic. Stars: 13 with 0, one with 1, 18 with 2, one with 3, one with 4.

The files have no founder column. Founder status was not assigned.

**PASS on the call that crosses 0.01, and two or more stars (15).** All 15 are germline pathogenic, not likely-pathogenic-only. Five also say low penetrance. Three also say risk factor or risk allele. None says drug response. AF values below are the stored frequencies rounded to 6 decimal places. AC and AN are unrounded. The stored name and phenotype are:

| gene | name | stars | label extras | genome AF (AC/AN, filter) | exome AF (AC/AN, filter) |
| --- | --- | ---: | --- | --- | --- |
| HFE | p.Cys282Tyr | 2 | low penetrance; risk factor | 0.037293 (5679/152282, PASS) | 0.059033 (86298/1461852, PASS) |
| CFTR | c.1210-12T[5] | 2 | phenotype includes a bronchiectasis modifier | 0.035695 (5090/142598, PASS) | 0.028617 (40213/1405206, PASS) |
| G6PD | p.Val98Met | 2 | chrX | 0.033489 (3754/112096, PASS) | 0.003682 (4043/1098048, PASS) |
| SERPINA1 | p.Glu288Val | 2 | low penetrance; also `other`; phenotype PI S | 0.028729 (4376/152320, PASS) | 0.037155 (54315/1461862, PASS) |
| RBM8A | c.-21G>A | 2 | low penetrance; also `other` | 0.018444 (2808/152246, PASS) | 0.026017 (37795/1452720, PASS) |
| WNT10A | p.Phe228Ile | 2 | low penetrance | 0.014046 (2138/152210, PASS) | 0.021492 (31289/1455814, PASS) |
| HBB | p.Glu7Val | 2 | phenotype hemoglobin S; symbol is `HBB;LOC106099062;LOC107133510` | 0.012719 (1937/152294, PASS) | 0.001601 (2335/1458356, PASS) |
| FLG | p.Ser761fs | 2 |  | 0.012588 (1913/151974, PASS) | 0.019637 (28705/1461820, AS_VQSR) |
| BCHE | p.Asp98Gly | 2 | phenotype includes postanesthetic apnea; the label does not say drug response | 0.012437 (1894/152286, PASS) | 0.016397 (23967/1461658, PASS) |
| SERPINA1 | p.Glu366Lys | 2 | risk factor; phenotype PI Z | 0.012157 (1851/152264, PASS) | 0.016249 (23754/1461878, PASS) |
| CYP21A2 | p.Val282Leu | 2 |  | 0.010854 (1639/150998, PASS) | 0.004388 (6337/1444170, PASS) |
| FLG | p.Arg501Ter | 2 |  | 0.010428 (1580/151516, PASS) | 0.017637 (25781/1461786, PASS) |
| F2 | c.*97G>A | 2 | low penetrance; established risk allele; risk factor | 0.009566 (1457/152304, PASS) | 0.011594 (15393/1327714, PASS) |
| GJB2 | p.Met34Thr | 3 | expert panel | 0.009406 (1433/152352, PASS) | 0.012316 (18002/1461676, PASS) |
| CFTR | p.Phe508del | 4 | practice guideline | 0.007884 (1200/152206, PASS) | 0.012353 (18037/1460114, PASS) |

HFE p.Cys282Tyr is the only row in this group with an AF at or above 0.05. CFTR p.Phe508del, GJB2 p.Met34Thr, and F2 c.*97G>A are in the tail because the exome AF is at or above 0.01. Their genome AFs are below 0.01. HBB and G6PD are the other way around: the genome AF is at or above 0.01 and the exome AF is not. G6PD is on chrX, so the stored AN is the X denominator (112,096 genomes, not the autosomal ~152,000). The FLG frameshift is `PASS` in genomes and `AS_VQSR` in exomes. Both of those AFs are shown because they were not collapsed.

**PASS and stars 0, no assertion criteria (6).** Germline. No drug response. No conflict. One risk-factor label.

| gene | name | phenotype in the column | genome AF (filter) | exome AF (filter) |
| --- | --- | --- | --- | --- |
| `-` | NC_000023.11:g.2748343G>C | Xg(a-) blood group | 0.354497 (PASS) | not in exomes |
| NCR3 | NC_000006.12:g.31593133C>G | malaria susceptibility | 0.129762 (PASS) | 0.079287 (PASS) |
| SERPINA7 | p.Asp191Asn | thyroxine-binding globulin, slow | 0.031420 (PASS) | 0.003444 (PASS) |
| NCR3 | c.-172G>A | malaria, severe, susceptibility | 0.023255 (PASS) | 0.021159 (PASS) |
| SERPINH1 | NC_000011.10:g.75561608C>T | preterm premature rupture of membranes | 0.022025 (PASS) | not in exomes |
| LGR4 | p.Asp844Gly | delayed puberty, self-limited | 0.013546 (PASS) | 0.020373 (PASS) |

The first NCR3 row is `Pathogenic; risk factor`. The second is `Pathogenic` only. The Xg allele and SERPINA7 are on chrX. The gene symbol stored for the Xg allele is `-`.

**PASS exome AF, one star, single submitter (1).** `EYA1` c.34-5645dup, likely pathogenic, germline, duplication, phenotype `EYA1-related disorders`. No drug response, no risk factor, no low penetrance, no conflict. Genome AF 0.190688 (16,395/85,978) with filter `AS_VQSR`. Exome AF 0.157920 (111,975/709,060) with filter `PASS`. The exome frequency is a PASS call and is common. There is no rsID on the row.

**The crossing AF is not PASS (10).** The allele string matched. gnomAD did not mark that call `PASS`. Stars and labels do not repair the filter.

| gene | name | stars | genome AF (filter) | exome AF (filter) |
| --- | --- | ---: | --- | --- |
| PRSS1 | p.Ala16Val | 2 | 0.386200 (InbreedingCoeff) | 0.202123 (InbreedingCoeff) |
| SLC9B1 | p.Leu447fs | 2 | 0.341192 (AS_VQSR;InbreedingCoeff) | 0.102450 (AS_VQSR;InbreedingCoeff) |
| VPS13B | c.9331-1G>T | 2 | 0.080981 (AS_VQSR) | 0.110904 (AS_VQSR) |
| TUBB8 | p.Gly109Val | 0 | 0.067775 (AS_VQSR;InbreedingCoeff) | 0.000112 (AS_VQSR) |
| CNN2 | p.Arg266Gln | 0 | 0.051132 (InbreedingCoeff) | 0.000109 (AS_VQSR) |
| CNN2 | p.Gly263Ser | 0 | 0.046822 (AS_VQSR;InbreedingCoeff) | 0.000191 (AS_VQSR) |
| PRB3 | p.His284fs | 0 | 0.041365 (AS_VQSR) | 0.014437 (AS_VQSR) |
| TMEM67 | c.224-2A>T | 2 | 0.033175 (AS_VQSR) | 0.007826 (AS_VQSR) |
| DNAH5 | c.5272-1G>T | 2 | 0.000274 (AS_VQSR) | 0.011380 (AS_VQSR) |
| BBS4 | p.Cys264Phe | 0 | 0.000053 (AS_VQSR) | 0.010898 (AS_VQSR) |

PRSS1, SLC9B1, and VPS13B are two-star germline labels (hereditary pancreatitis, Usher syndrome, Cohen syndrome) on failed filters, with both AFs high. SLC9B1's ClinVar type is Microsatellite; the stored name is a two-base deletion, and the matched alleles are 3 bases to 1 base. TUBB8 and both CNN2 rows are stars 0, and the genome AF is about 0.05 while the exome AF is about 0.0001. DNAH5 and BBS4 go the other way. TMEM67's genome AN is 41,628, the smallest AN in this tail, and the filter is `AS_VQSR`. PRB3's phenotype is `PRB3M(NULL)`. None of these ten says drug response, low penetrance, or risk factor.

**Not germline (2).** Both are somatic, stars 0, no assertion criteria, significance `Pathogenic`. `ETV5` p.His179Pro, myoepithelial tumor, genome AF 0.009641 (`AS_VQSR`, AN 14,106), exome AF 0.017473 (`AS_VQSR`, AN 204,310). `PPP2R1B` p.Gly90Asp, lung carcinoma, both filters `PASS`, genome AF 0.006914, exome AF 0.010873. These are the only two alleles in the tail with the germline flag false.

The band `0.001–<0.01` was counted and not opened allele by allele. Genome: 163 alleles, 160 germline, 3 low penetrance, 0 drug response, 3 risk, 0 conflicting review status, 133 with two or more stars, 152 with filter `PASS`. Exome: 188 alleles, 183 germline, 2 low penetrance, 0 drug response, 2 risk, 0 conflicting review status, 153 with two or more stars, 162 `PASS`.

### Bins

A zero means the SQL returned no row for that cell. Row totals match the class counts above.

Genome AF:

| class | null | <0.0001 | 0.0001–<0.001 | 0.001–<0.01 | 0.01–<0.05 | >=0.05 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pathogenic | 215117 | 26912 | 1288 | 140 | 19 | 6 |
| likely pathogenic | 122323 | 10267 | 253 | 23 | 0 | 2 |
| conflicting | 55966 | 84462 | 33677 | 9054 | 414 | 241 |
| vus | 1330101 | 957378 | 102766 | 7025 | 429 | 555 |
| benign | 30332 | 29797 | 32178 | 52912 | 45677 | 88493 |
| likely benign | 612805 | 409357 | 65152 | 42716 | 13244 | 1175 |
| drug response | 1453 | 190 | 42 | 33 | 47 | 128 |
| risk or other effect | 442 | 171 | 44 | 30 | 28 | 290 |
| other | 1171 | 137 | 30 | 19 | 17 | 145 |
| no classification | 176423 | 29993 | 13117 | 14360 | 8314 | 10373 |

Exome AF:

| class | null | <0.0001 | 0.0001–<0.001 | 0.001–<0.01 | 0.01–<0.05 | >=0.05 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pathogenic | 172598 | 69346 | 1356 | 161 | 16 | 5 |
| likely pathogenic | 98238 | 34334 | 266 | 27 | 2 | 1 |
| conflicting | 25042 | 123290 | 28640 | 6213 | 398 | 231 |
| vus | 699657 | 1621004 | 72695 | 4193 | 352 | 353 |
| benign | 54675 | 50727 | 51270 | 46933 | 18091 | 57693 |
| likely benign | 309467 | 748259 | 59087 | 20753 | 5464 | 1419 |
| drug response | 1190 | 505 | 35 | 28 | 31 | 104 |
| risk or other effect | 512 | 304 | 43 | 24 | 21 | 101 |
| other | 1335 | 165 | 5 | 3 | 0 | 11 |
| no classification | 145066 | 72361 | 14075 | 10247 | 4153 | 6678 |

Germline-flag bins, pathogenic and likely pathogenic together: genome null 300,845, `<0.0001` 36,040, `0.0001–<0.001` 1,515, `0.001–<0.01` 160, `0.01–<0.05` 19, `>=0.05` 8. Exome null 238,251, `<0.0001` 98,532, `0.0001–<0.001` 1,599, `0.001–<0.01` 183, `0.01–<0.05` 16, `>=0.05` 6.

### Join check

The join is not broken, so the frequency claim is not stopped. Dedup rows 4,639,253, joined rows 4,639,253. Distinct keyed alleles 4,618,645, extra key rows 0, null keys 20,608. Tier A rows 5,928,738, extra key rows 0. Keyed alleles missing from tier A: 0. Ref, alt, chrom, or pos mismatches: 0. An AF is never stored when `in_genomes` or `in_exomes` is false. Across every ClinVar allele, not just pathogenic ones, 78 genome matches and 3,134 exome matches have a null AF. That is a null frequency on a matched site, not a duplicated join. The pathogenic slice of those nulls is the CBS and HGF rows above.

## Caveat

Most pathogenic and likely pathogenic alleles have no AF in one or both call sets. Below 0.01 describes 38,883 of the 38,910 with a genome AF, not 38,883 of all 376,350. The 337,440 null genome AFs were not counted as rare. Not being in the sites VCF is what a rare allele often looks like, and it is also what a structural allele with no matching VCF record looks like. Both are in the null. The SNV share of the null is above, so the null is not only structural.

The bin tables trust the stored AF, including sites whose filter is `AS_VQSR` or `InbreedingCoeff`. The alleles above 0.05 are mostly those failed filters, or stars-0 labels with no assertion criteria. Treating the `>=0.05` bin as common pathogenic disease alleles would be wrong. HFE p.Cys282Tyr is the one two-or-more-star PASS allele in that bin.

Genome AF and exome AF disagree on several of the 34, including HBB, G6PD, TUBB8, CNN2, DNAH5, and BBS4. Both numbers are reported. They were not averaged. Three of the 34 are on chrX (G6PD, SERPINA7, and the Xg allele). Their AF uses the AN stored for that site.

Low penetrance and risk-factor wording are on some of the PASS, high-star exceptions. They are not most of the tail, and they are rare in the `0.001–<0.01` band. Drug-response rows and conflicting labels were expected and were not what this tail showed. The drug-response class is common when it has an AF, and it sits in its own bin.
