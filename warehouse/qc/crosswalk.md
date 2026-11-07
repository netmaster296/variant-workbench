# Crosswalk

Lane G. One row per versioned GENCODE gene id. Counts below were queried from the Parquet files after they were written.

Genes: `D:\variant-workbench\warehouse\crosswalk\genes.parquet` (3,395,821 bytes).
Unmatched: `D:\variant-workbench\warehouse\crosswalk\unmatched.parquet` (120,929 bytes).

## Row count

Status: pass.

| check | count |
| --- | --- |
| genes.parquet rows | 78,691 |
| distinct ensg_version | 78,691 |
| distinct ensg | 78,691 |
| interval genes.parquet rows | 78,691 |
| expected | 78,691 |

The crosswalk does not add a row for an NCBI id or a RefSeq id that failed to join. Those ids stay on the gene they matched, or in `unmatched.parquet`.

## Join rules

Primary key: bare ENSG (`coords.bare_ensg`). Symbol is not a gene key.

ClinVar `GeneID` and `HGNC_ID` attach by exact `GeneSymbol` on single-gene rows (`GeneID` is digits, symbol has no space, comma, or semicolon). `GeneID` `-1` is not an NCBI gene id. Semicolon-separated symbols on those rows are split only to see which symbols failed to match. They do not assign an NCBI id.

If one symbol has several NCBI GeneIDs, `ncbi_gene_id` keeps every id, semicolon-separated, sorted numerically. `ncbi_gene_id_primary` is the chosen id: the most allele rows, then the smaller numeric id on a tie. `ncbi_gene_id_rule` is `sole`, `most_alleles`, or `smallest_id`. The gene is not dropped. Every GENCODE gene with that symbol gets the same ClinVar ids. `symbol_n_gencode` is how many GENCODE genes share the symbol. `ncbi_clinvar_symbols` lists every ClinVar symbol that uses the primary NCBI id.

`hgnc_id` is the GTF `hgnc_id` when the gene has one. Otherwise it is the ClinVar HGNC id when that symbol has exactly one. If ClinVar has several and the GTF has none, `hgnc_id` stays null and `hgnc_id_clinvar` keeps the list. No id is dropped. `hgnc_disagree` is true when both sources are present and the sets differ. Ids are compared as `HGNC:` plus the integer, so a bare number or leading zeros would not count as a disagreement.

`mane_enst` is the bare ENST. The GENCODE tag `MANE_Select` wins when the gene has one. Otherwise the id comes from a constraint transcript with `mane_select` and an ENST id on that bare ENSG. `mane_enst_gencode` and `mane_enst_constraint` both stay populated. `mane_enst_disagree` is true when both are present and the sets differ. If a source has several ENST ids, the stored `mane_enst` is the lexicographically smallest and the full set stays in the source column. `mane_enst_version` is the GENCODE versioned id when the chosen source is GENCODE. Constraint ENST ids in this file have no version, so the version is null when the source is constraint.

`mane_refseq` is the constraint transcript with `mane_select` whose id starts with `NM_` or `NR_`. RefSeq rows have a null bare ENSG. Pairing stays inside the constraint file: if the symbol has exactly one MANE ENST, the RefSeq id attaches to that transcript's bare ENSG; if it has no MANE ENST and exactly one bare ENSG on any constraint transcript, it attaches to that ENSG; otherwise it is unattached. It is not joined to GENCODE by symbol. Null when none attach.

`constraint_ensg_matched` is true when the gene's bare ENSG equals `gene_constraint.bare_ensg`. A constraint row whose `bare_ensg` is null is unmatched on id, not a failed symbol join.

`MANE_Plus_Clinical` is not MANE Select.

Identifier check against coords.bare_ensg and coords.bare_enst: GENCODE gene rows differing 0; MANE transcripts differing 0; gene-constraint rows differing 0; transcript-constraint rows differing 0.

ClinVar single-gene symbols: 20,327. Symbols with more than one NCBI GeneID: 0.
Single-gene ClinVar rows rejected as non-symbols: 0. Multigene rows with an empty or prose token: 0.
MANE RefSeq transcripts left unattached: 11. They have no unique bare ENSG inside the constraint file, so they are not symbol-joined to GENCODE.
Unattached MANE RefSeq transcripts (up to 15): EEF1AKNMT NM_015935.5, FAM207A NM_058190.4, LOC100133315 NM_001375847.2, LOC100421372 NM_001378785.1, LOC100505841 NM_001195535.4, LOC102724488 NM_001370184.1, LOC112577516 NM_001364461.3, LOC112694756 NM_001365304.2, LOC114841035 NM_001370414.1, LOC728392 NM_001162371.3, METTL11B NM_001136107.2.

## Protein-coding constraint match

Match is `biotype = protein_coding` and `constraint_ensg_matched`. An earlier partial join saw 18,567 of 20,097 protein-coding genes and about 1,530 misses. That figure is not a target. These counts are recomputed from `genes.parquet`.

| set | genes | rate of protein-coding |
| --- | --- | --- |
| protein-coding | 20,097 | 100% |
| matched to a constraint bare ENSG | 18,567 | 92.39% |
| unmatched | 1,530 | 7.61% |
| all biotypes matched | 18,613 |  |
| all biotypes unmatched | 60,078 |  |

Protein-coding genes that missed on bare ENSG were not rescued by symbol. The next table only describes what a symbol comparison would have seen. It is not the join.

| protein-coding unmatched on ENSG | genes |
| --- | --- |
| total | 1,530 |
| symbol hits a constraint row with a different bare ENSG | 5 |
| symbol hits only a constraint row whose bare ENSG is null | 6 |
| no constraint row with that symbol | 1,519 |

These symbol hits were not used as a join.

| symbol | GENCODE ENSG | constraint ENSG | constraint transcript |
| --- | --- | --- | --- |
| CAST | ENSG00000310517 | ENSG00000153113 | ENST00000675179 |
| FAM174C | ENSG00000310576 | ENSG00000228300 | ENST00000409293 |
| HOMEZ | ENSG00000290292 | ENSG00000215271 | ENST00000357460 |
| PAXX | ENSG00000310560 | ENSG00000148362 | ENST00000371620 |
| SOD2 | ENSG00000291237 | ENSG00000112096 | ENST00000538183 |

Protein-coding symbols that meet only a constraint row with a null bare ENSG. That row is unmatched on id.

| symbol | GENCODE ENSG | constraint gene_id | constraint transcript |
| --- | --- | --- | --- |
| NBPF8 | ENSG00000270231 | 728841 | NM_001037501.4 |
| NPIPA3 | ENSG00000224712 | 642778 | XM_011522584.2 |
| POLR2J2 | ENSG00000228049 | 246721 | NM_032959.6 |
| SPATA31C1 | ENSG00000230246 | 441452 | XM_011518702.1 |
| TARP | ENSG00000289746 | 445347 | XM_017012226.1 |
| ZNF883 | ENSG00000228623 | 169834 | NM_001101338.1 |

## Duplicate symbols

Symbols on more than one GENCODE gene: 483, covering 2,096 genes. Symbols that also occur on a protein-coding gene: 36. Symbols shared by more than one protein-coding gene: 27 (54 protein-coding genes).

Every gene is kept. ClinVar ids for that symbol are copied onto each gene. `symbol_n_gencode` records the fan-out.

| genes with the symbol | symbols |
| --- | --- |
| 2 | 453 |
| 3 | 9 |
| 4 | 4 |
| 5 | 1 |
| 6 | 1 |
| 7 | 6 |
| 8 | 1 |
| 9 | 1 |
| 19 | 1 |
| 22 | 1 |
| 27 | 1 |
| 33 | 1 |
| 50 | 1 |
| 170 | 1 |
| 756 | 1 |

Largest duplicate symbols:

| symbol | genes | protein_coding | biotypes |
| --- | --- | --- | --- |
| Y_RNA | 756 | 0 | misc_RNA |
| Metazoa_SRP | 170 | 0 | misc_RNA |
| U3 | 50 | 0 | snoRNA |
| U6 | 33 | 0 | snRNA |
| SNORA70 | 27 | 0 | snoRNA |
| U8 | 22 | 0 | snoRNA |
| U2 | 19 | 0 | snRNA |
| 5S_rRNA | 9 | 0 | rRNA, rRNA_pseudogene |
| SNORA72 | 8 | 0 | snoRNA |
| 7SK | 7 | 0 | misc_RNA |
| SNORA62 | 7 | 0 | snoRNA |
| SNORA63 | 7 | 0 | snoRNA |
| SNORA75 | 7 | 0 | snoRNA |
| U4 | 7 | 0 | snRNA |
| U7 | 7 | 0 | snRNA |

Protein-coding symbols that are not unique:

| symbol | genes | ensgs |
| --- | --- | --- |
| AKAP17A | 2 | ENSG00000197976, ENSG00000292343 |
| ASMT | 2 | ENSG00000196433, ENSG00000292336 |
| ASMTL | 2 | ENSG00000169093, ENSG00000292339 |
| C4orf36 | 2 | ENSG00000163633, ENSG00000285458 |
| CD99 | 2 | ENSG00000002586, ENSG00000292348 |
| CRLF2 | 2 | ENSG00000205755, ENSG00000292363 |
| CSF2RA | 2 | ENSG00000198223, ENSG00000292357 |
| DHRSX | 2 | ENSG00000169084, ENSG00000292338 |
| DUSP13B | 2 | ENSG00000079393, ENSG00000293542 |
| FAM174C | 2 | ENSG00000228300, ENSG00000310576 |
| GTPBP6 | 2 | ENSG00000178605, ENSG00000292358 |
| HERC3 | 2 | ENSG00000138641, ENSG00000287542 |
| IL3RA | 2 | ENSG00000185291, ENSG00000292332 |
| IL9R | 2 | ENSG00000124334, ENSG00000292373 |
| KYAT1 | 2 | ENSG00000171097, ENSG00000286112 |
| P2RY8 | 2 | ENSG00000182162, ENSG00000292333 |
| PDE4C | 2 | ENSG00000105650, ENSG00000285188 |
| PDE8B | 2 | ENSG00000113231, ENSG00000284762 |
| PINX1 | 2 | ENSG00000254093, ENSG00000258724 |
| PLCXD1 | 2 | ENSG00000182378, ENSG00000292344 |
| POLR2J3 | 2 | ENSG00000168255, ENSG00000285437 |
| PPP2R3B | 2 | ENSG00000167393, ENSG00000292327 |
| SHOX | 2 | ENSG00000185960, ENSG00000292354 |
| SLC25A6 | 2 | ENSG00000169100, ENSG00000292334 |
| VAMP7 | 2 | ENSG00000124333, ENSG00000292366 |
| WASH6P | 2 | ENSG00000182484, ENSG00000292372 |
| ZBED1 | 2 | ENSG00000214717, ENSG00000292345 |

## ClinVar symbols with no GENCODE symbol

Rows are in `unmatched.parquet` with `source = clinvar`. `n` is allele rows. For `no_gencode_symbol`, `n_single` is single-gene rows and `n_multi` is multigene rows whose semicolon-separated symbol list contains the symbol. A symbol that matches any GENCODE gene is absent here, including when several genes share it.

| reason | keys | allele rows | single-gene rows | multigene rows |
| --- | --- | --- | --- | --- |
| missing_symbol | 1 | 1,655 |  |  |
| no_gencode_symbol | 7,534 | 21,507 | 4,924 | 16,583 |
| not_a_gene_symbol | 3,701 | 9,143 |  |  |

`no_gencode_symbol` is the ClinVar symbol miss. `not_a_gene_symbol` is prose such as a dosage-coverage label. `missing_symbol` is GeneSymbol `-`. Those two are not symbol joins.

## Constraint ids with no GENCODE gene

`gene_constraint` rows with null `bare_ensg`: 88. A null bare ENSG is unmatched on id. It is not counted as a symbol miss.

| reason | gene rows | transcript rows with that gene_id | mane_select gene rows |
| --- | --- | --- | --- |
| bare_ensg_not_in_gencode | 15 | 51 | 4 |
| bare_ensg_null | 88 | 231 | 11 |

`n` on a constraint unmatched row is the number of `transcript_constraint` rows with that `gene_id`. NCBI ids that the constraint loader collapsed onto an Ensembl gene are not listed again here. The gene id that remains is the bare ENSG, and it is unmatched only when that ENSG is absent from GENCODE.

## NCBI and HGNC

| metric | count |
| --- | --- |
| genes with an NCBI GeneID | 19,891 |
| genes with no NCBI GeneID | 58,800 |
| genes whose symbol has more than one NCBI GeneID | 0 |
| primary rule sole | 19,891 |
| primary rule most_alleles | 0 |
| primary rule smallest_id | 0 |
| genes whose primary NCBI id is shared by more than one ClinVar symbol | 82 |
| distinct primary NCBI ids shared by more than one ClinVar symbol | 82 |
| genes with GTF and ClinVar HGNC | 19,870 |
| genes with GTF HGNC only | 21,707 |
| genes with ClinVar HGNC only | 14 |
| genes whose resolved hgnc_id is null | 37,100 |
| genes whose ClinVar HGNC list has more than one id | 0 |
| HGNC disagreements | 0 |

NCBI ids used by more than one ClinVar symbol (up to 20):

| ncbi_gene_id | clinvar symbols | gencode genes carrying that id |
| --- | --- | --- |
| 100129250 | SMIM27;TZMP1 | 1 |
| 100130933 | ERLN;SMIM6 | 1 |
| 100506243 | KRABD1;KRBOX1 | 1 |
| 10227 | MFSD10;SLC75A1 | 1 |
| 10994 | HACL2;ILVBL | 1 |
| 11149 | BVES;POPDC1 | 1 |
| 11318 | ACKR5;GPR182 | 1 |
| 113655 | MFSD3;SLC33A2 | 1 |
| 113828 | FAM83F;SACK1F | 1 |
| 115811 | DRC10;IQCD | 1 |
| 116151 | FAM210B;MIMS2 | 1 |
| 124411 | KRABD5;KRBOX5 | 1 |
| 125228 | FAM210A;MIMS1 | 1 |
| 146712 | B3GNTL1;QTGAL | 1 |
| 147007 | TMEM199;VMA12 | 1 |
| 148808 | MFSD4A;SLC60A1 | 1 |
| 151254 | C2CD6;CATSPERT | 1 |
| 160419 | C12orf50;ZC3H11D | 1 |
| 161291 | CDC50B;TMEM30B | 1 |
| 161823 | ADAL;MAPDA | 1 |

No gene has both a GTF HGNC id and a different ClinVar HGNC id.

## MANE

GENCODE MANE Select is the transcript tag `MANE_Select`, not `MANE_Plus_Clinical`. Constraint MANE is `gene_constraint.mane_select`. Overlap is a constraint MANE row whose `bare_ensg` is a GENCODE MANE gene. The two difference sets are GENCODE-only and constraint-only. Neither is dropped from the note.

| set | genes |
| --- | --- |
| GENCODE genes with mane_select | 19,276 |
| GENCODE genes with a MANE_Select transcript tag | 19,276 |
| GENCODE MANE_Select transcripts | 19,276 |
| constraint gene rows with mane_select | 17,494 |
| overlap on bare ENSG | 17,471 |
| GENCODE only | 1,805 |
| constraint only | 23 |
| constraint only, bare ENSG set, not a GENCODE MANE gene | 12 |
| constraint only, bare ENSG null (unmatched on id) | 11 |
| constraint MANE bare ENSG absent from GENCODE entirely | 4 |

Set identity (overlap + GENCODE only = GENCODE; overlap + constraint only = constraint): pass.

On the crosswalk row, a GENCODE MANE transcript is preferred for `mane_enst` and the constraint id is still stored. A constraint-only MANE gene that joins on bare ENSG keeps the constraint ENST when it has one, and the RefSeq id when the pair attached.

| crosswalk count | genes |
| --- | --- |
| mane_enst_gencode set | 19,276 |
| mane_enst_constraint set | 17,477 |
| both ENST ids present and equal | 17,441 |
| both ENST ids present and different | 28 |
| mane_enst chosen from GENCODE | 19,276 |
| mane_enst chosen from constraint | 8 |
| mane_enst null | 59,407 |
| mane_refseq set | 17,466 |
| more than one mane_refseq | 0 |
| more than one GENCODE MANE ENST | 0 |
| more than one constraint MANE ENST | 0 |
| constraint NM_/NR_ mane transcripts in the transcript file | 17,481 |
| those RefSeq ids stored on unmatched constraint rows | 15 |
| crosswalk plus unmatched RefSeq ids | 17,481 |

Every constraint `NM_` or `NR_` MANE transcript is either `mane_refseq` on the GENCODE gene whose bare ENSG the constraint file paired it to, or `mane_refseq` on the unmatched constraint row when that ENSG is missing or null. Accounting: pass (17,466 + 15 = 17,481; transcript file 17,481).

Transcript id disagreements: 28. `mane_enst` uses the GENCODE id. The constraint id stays in `mane_enst_constraint`.

| symbol | ensg | gencode | constraint | chosen | mane_refseq |
| --- | --- | --- | --- | --- | --- |
| ACYP2 | ENSG00000170634 | ENST00000394666 | ENST00000607452 | ENST00000394666 | NM_001320586.2 |
| ATAT1 | ENSG00000137343 | ENST00000376485 | ENST00000330083 | ENST00000376485 | NM_001031722.4 |
| BRAF | ENSG00000157764 | ENST00000646891 | ENST00000644969 | ENST00000646891 | NM_001374258.1 |
| CCL4L2 | ENSG00000276070 | ENST00000620576 | ENST00000617405 | ENST00000620576 | NM_001291475.2 |
| CCL5 | ENSG00000271503 | ENST00000605140 | ENST00000651122 | ENST00000605140 | NM_001278736.2 |
| CTBP2 | ENSG00000175029 | ENST00000337195 | ENST00000309035 | ENST00000337195 | NM_022802.3 |
| FAM149A | ENSG00000109794 | ENST00000706927 | ENST00000389354 | ENST00000706927 | NM_001367768.3 |
| GLT1D1 | ENSG00000151948 | ENST00000442111 | ENST00000281703 | ENST00000442111 | NM_144669.3 |
| GP6 | ENSG00000088053 | ENST00000417454 | ENST00000310373 | ENST00000417454 | NM_001083899.2 |
| GUK1 | ENSG00000143774 | ENST00000453943 | ENST00000312726 | ENST00000453943 | NM_000858.7 |
| HNF4A | ENSG00000101076 | ENST00000316673 | ENST00000316099 | ENST00000316673 | NM_000457.6 |
| IL12A | ENSG00000168811 | ENST00000699704 | ENST00000305579 | ENST00000699704 | NM_000882.4 |
| IL36B | ENSG00000136696 | ENST00000327407 | ENST00000259213 | ENST00000327407 | NM_014438.5 |
| KLF12 | ENSG00000118922 | ENST00000703967 | ENST00000377669 | ENST00000703967 | NM_007249.5 |
| MGA | ENSG00000174197 | ENST00000703841 | ENST00000219905 | ENST00000703841 | NM_001164273.2 |
| NMNAT3 | ENSG00000163864 | ENST00000704800 | ENST00000643695 | ENST00000704800 | NM_001320510.2 |
| PLCXD2 | ENSG00000240891 | ENST00000636933 | ENST00000477665 | ENST00000636933 | NM_001185106.1 |
| RAF1 | ENSG00000132155 | ENST00000251849 | ENST00000442415 | ENST00000251849 | NM_001354689.3 |
| RGS20 | ENSG00000147509 | ENST00000276500 | ENST00000297313 | ENST00000276500 | NM_170587.4 |
| SLC35B3 | ENSG00000124786 | ENST00000710437 | ENST00000644923 | ENST00000710437 | NM_001370476.2 |
| SSBP2 | ENSG00000145687 | ENST00000615665 | ENST00000320672 | ENST00000615665 | NM_012446.5 |
| TBC1D1 | ENSG00000065882 | ENST00000698857 | ENST00000261439 | ENST00000698857 | NM_015173.4 |
| TCF7L2 | ENSG00000148737 | ENST00000355995 | ENST00000627217 | ENST00000355995 | NM_001146274.2 |
| TGFBR3L | ENSG00000260001 | ENST00000713907 | ENST00000565886 | ENST00000713907 | NM_001195259.2 |
| TMEM230 | ENSG00000089063 | ENST00000202834 | ENST00000342308 | ENST00000202834 | NM_001009923.2 |
| TRIM61 | ENSG00000183439 | ENST00000710271 | ENST00000329314 | ENST00000710271 | NM_001012414.3 |
| TXNDC8 | ENSG00000204193 | ENST00000374511 | ENST00000423740 | ENST00000374511 | NM_001286946.2 |
| UFSP1 | ENSG00000176125 | ENST00000672365 | ENST00000388761 | ENST00000672365 | NM_001015072.4 |

Constraint MANE rows that are not GENCODE MANE genes: 23. Null `bare_ensg` rows are unmatched on id. `in_gencode` false means that bare ENSG is not a GENCODE gene, so the row is not in `genes.parquet`. A same symbol in GENCODE was not given this constraint id.

| symbol | bare_ensg | in_gencode | gene_id | transcript | biotype | gencode ENST | constraint ENST | refseq |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EEF1AKNMT |  | no | 51603 | NM_015935.5 |  |  |  | NM_015935.5 |
| FAM207A |  | no | 85395 | NM_058190.4 |  |  |  | NM_058190.4 |
| LOC100133315 |  | no | 100133315 | NM_001375847.2 |  |  |  | NM_001375847.2 |
| LOC100421372 |  | no | 100421372 | NM_001378785.1 |  |  |  | NM_001378785.1 |
| LOC100505841 |  | no | 100505841 | NM_001195535.4 |  |  |  | NM_001195535.4 |
| LOC102724488 |  | no | 102724488 | NM_001370184.1 |  |  |  | NM_001370184.1 |
| LOC112577516 |  | no | 112577516 | NM_001364461.3 |  |  |  | NM_001364461.3 |
| LOC112694756 |  | no | 112694756 | NM_001365304.2 |  |  |  | NM_001365304.2 |
| LOC114841035 |  | no | 114841035 | NM_001370414.1 |  |  |  | NM_001370414.1 |
| LOC728392 |  | no | 728392 | NM_001162371.3 |  |  |  | NM_001162371.3 |
| METTL11B |  | no | 149281 | NM_001136107.2 |  |  |  | NM_001136107.2 |
| SOD2 | ENSG00000112096 | no | ENSG00000112096 | ENST00000538183 |  |  |  | NM_000636.4 |
| SMIM2 | ENSG00000139656 | no | ENSG00000139656 | ENST00000400419 |  |  |  | NM_024058.3 |
| KAAG1 | ENSG00000146049 | yes | ENSG00000146049 | ENST00000274766 | lncRNA |  | ENST00000274766 | NM_181337.4 |
| PAXX | ENSG00000148362 | no | ENSG00000148362 | ENST00000371620 |  |  |  | NM_183241.3 |
| CAST | ENSG00000153113 | yes | ENSG00000153113 | ENST00000675179 | lncRNA |  | ENST00000675179 | NM_001750.7 |
| GPATCH4 | ENSG00000160818 | yes | ENSG00000160818 | ENST00000368232 | protein_coding |  | ENST00000368232 | NM_182679.3 |
| ZNF781 | ENSG00000196381 | yes | ENSG00000196381 | ENST00000358582 | transcribed_unprocessed_pseudogene |  | ENST00000358582 | NM_152605.4 |
| HOMEZ | ENSG00000215271 | no | ENSG00000215271 | ENST00000357460 |  |  |  | NM_020834.3 |
| FAM174C | ENSG00000228300 | yes | ENSG00000228300 | ENST00000409293 | protein_coding |  | ENST00000409293 | NM_017914.4 |
| ZNF559-ZNF177 | ENSG00000270011 | yes | ENSG00000270011 | ENST00000683217 | protein_coding |  | ENST00000683217 |  |
| TSTD3 | ENSG00000279170 | yes | ENSG00000279170 | ENST00000452647 | transcribed_unitary_pseudogene |  | ENST00000452647 | NM_001195131.2 |
| IQCJ-SCHIP1 | ENSG00000283154 | yes | ENSG00000283154 | ENST00000638749 | protein_coding |  | ENST00000638749 |  |

The same MANE ENST is on different ENSG ids in GENCODE and in constraint: 7. Each source stays on its own gene. The transcript was not copied across.

| ENST | GENCODE symbol | GENCODE ENSG | constraint symbol | constraint ENSG | constraint ENSG in GENCODE |
| --- | --- | --- | --- | --- | --- |
| ENST00000357460 | HOMEZ | ENSG00000290292 | HOMEZ | ENSG00000215271 | no |
| ENST00000371620 | PAXX | ENSG00000310560 | PAXX | ENSG00000148362 | no |
| ENST00000409293 | FAM174C | ENSG00000310576 | FAM174C | ENSG00000228300 | yes |
| ENST00000538183 | SOD2 | ENSG00000291237 | SOD2 | ENSG00000112096 | no |
| ENST00000638749 | SCHIP1 | ENSG00000151967 | IQCJ-SCHIP1 | ENSG00000283154 | yes |
| ENST00000675179 | CAST | ENSG00000310517 | CAST | ENSG00000153113 | yes |
| ENST00000683217 | ZNF177 | ENSG00000188629 | ZNF559-ZNF177 | ENSG00000270011 | yes |

## Column notes

`ncbi_gene_id` is every ClinVar GeneID for the symbol. `ncbi_gene_id_primary` is the chosen one. They are equal when the symbol has one GeneID.

`hgnc_id` is the resolved id. `hgnc_id_gtf` and `hgnc_id_clinvar` are the sources, including when they agree.

`mane_enst` is bare. `mane_enst_version` has the GENCODE version suffix when the source is GENCODE.
