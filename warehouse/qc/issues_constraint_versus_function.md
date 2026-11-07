# Span join is not the assayed gene

`mavedb_variant_keys_in_span` counts a GRCh38 key whose position falls inside the gene. That gene is not the assay target.

`mapped_variants_scored.parquet` has 212102 distinct keys with a chromosome and position, on 141 target names (238232 distinct key, target, and position rows). Those keys fall in 124 gene spans: 98 protein-coding and 26 lncRNA. 201746 of the keys fall in at least one protein-coding span. 7609 of those fall in two or three protein-coding genes.

Of the 98 protein-coding genes, 53 contain at least one key whose `target_gene` equals the symbol. 6 of those 53 also contain keys under a different name. 45 contain no key whose `target_gene` equals the symbol.

Some of the 45 are the same gene under another string. DDX3X has 12332 keys named `DDX3X Exon 1` through `DDX3X Exon 17`, plus `NC_000023.11`, and none named `DDX3X`. NPC1 has 975 keys named `ENST00000269228.10`, which is the versioned MANE transcript `ENST00000269228` on NPC1. SMN1 has 592 keys named `SMN Tudor domain`. YAP1 has 306 keys named `hYAP65 WW domain`.

Others are a different gene's assay, or an enhancer, sitting in this span:

| symbol | keys | target_gene on the key |
| --- | ---: | --- |
| ENSG00000267618 | 4606 | RAD51D |
| DHRS9 | 1860 | ECR11 enhancer |
| CELSR2 | 1800 | SORT1 enhancer |
| POU5F1B | 1792 | MYC enhancer (rs6983267) |
| LMBR1 | 1455 | ZRS enhancer |
| SEPTIN5 | 1154 | GP1BB promoter |
| ENSG00000284874 | 1154 | GP1BB promoter |
| ENSG00000284931 | 822 | HBG1 promoter |
| ENSG00000276490 | 625 | CYP2C19 |
| CCNH | 589 | RASA1 |
| MRNIP | 476 | SQSTM1 |
| EPHA8 | 29 | EPHB2 SAM domain |

UBC's 682 keys mix `UBC (L15S)` with `UBB`. The full list of the 45 is the last query in `D:\variant-workbench\warehouse\reports\sql\constraint_versus_function.sql`.

A score joined through the span would be credited to the overlapped gene. That claim is stopped. `constraint_versus_function.md` does not estimate a MaveDB effect size from this join. An exact `target_name` = symbol count is a separate filter. It does not repair the span join.

The mart column was not changed.
