# Constraint

Lane B. Source `D:\variant-workbench\raw\gnomad\constraint\gnomad.v4.1.constraint_metrics.tsv`.

Truth set `D:\variant-workbench\warehouse\qc\truth_set.tsv` was read and was not rewritten. Constraint cases in that file are `mane_select`, `loeuf_low_1`, `loeuf_low_2`, `loeuf_low_3`, `loeuf_near_1_1`, `loeuf_near_1_2`, and `loeuf_near_1_3`.

Required checks: pass.

## Counts

| Metric | Count |
|---|---:|
| Transcript rows | 211,523 |
| Gene rows | 18,716 |
| Source data rows read | 211,523 |
| Raw distinct `gene_id` values | 36,721 |
| `canonical=true` transcript rows | 36,315 |
| `mane_select=true` transcript rows | 34,962 |
| Genes with a MANE transcript | 17,494 |
| Genes with only a canonical transcript | 1,201 |
| Genes with neither MANE nor canonical | 21 |
| Genes with no MANE transcript | 1,222 |
| Genes dropped by the preference rule | 18,005 |
| Dropped Ensembl gene ids | 0 |
| Dropped NCBI gene ids | 18,005 |
| Transcript rows not kept on the gene table | 192,807 |

`canonical=true` is not one row per gene. The transcript table has 36,315 canonical rows and the gene table has 18,716 rows. A1BG is stored twice as canonical in the transcript file, `NM_130786.4` and `ENST00000263100`.

DuckDB reread the Parquet files. Transcript rows 211,523. Gene rows 18,716. Gene `gene_id` values are unique: pass.

## Preference

One gene row is the winning transcript, not every raw `gene_id`. RefSeq rows carry an NCBI gene id. Ensembl rows carry a bare ENSG. Those are the same gene when they share a symbol and that symbol has one Ensembl id, so they compete for the metric. The stored `gene_id` stays the Ensembl id, even when the winning transcript is RefSeq.

Order: `mane_select` true, then `canonical` true, then longest `cds_length`. A null `cds_length` does not beat a measured length. When a RefSeq transcript and an Ensembl transcript still tie, the Ensembl transcript is kept (`bare_enst` is set).

Symbols with more than one Ensembl gene id stay one row per Ensembl id. The symbol `NA` in the source is a missing symbol, stored null, and is not a group key. Gene rows with a null symbol: 508.

Symbols with two Ensembl gene rows: `MATR3`, `PINX1`, `POLR2J3`, `SIGLEC5`, `TBCE`.

Genes dropped by the preference rule are gene ids present in the TSV that are not the `gene_id` of a gene-table row. 18,005 are NCBI ids that lost to an Ensembl transcript of the same symbol. 0 are Ensembl ids that lost because a RefSeq transcript outranked them:

- None.

## A1BG

Gene-table transcript is `ENST00000263100` (`gene_id` `ENSG00000121410`, LOEUF 1.34, decile 7). Check that this is `ENST00000263100`: pass.

Transcript-file canonical copies include both `NM_130786.4` and `ENST00000263100`: pass.

| Transcript | Canonical | MANE | gene_id | cds_length | LOEUF | Decile |
|---|---|---|---|---:|---:|---:|
| `ENST00000263100` | True | True | `ENSG00000121410` | 1485 | 1.34 | 7 |
| `ENST00000600966` | False | False | `ENSG00000121410` | 917 | 1.256 | null |
| `NM_130786.4` | True | True | `1` | null | 1.34 | null |

## Truth set

Cited from `D:\variant-workbench\warehouse\qc\truth_set.tsv`. LOEUF must be within 0.001 of the truth-set value. Decile must match as an integer. The transcript must match the Ensembl id named in that row. "About 0.04", "about 0.041", and "about 0.048" are those truth-set values.

| Case | Gene | Expected | Observed | Result |
|---|---|---|---|---|
| `mane_select` | A1BG | ENST00000263100 | gene_id ENSG00000121410; transcript ENST00000263100; bare_ensg ENSG00000121410; LOEUF 1.34; decile 7; mane_select True; canonical True | pass |
| `loeuf_low_1` | XPO6 | LOEUF 0.04; decile 0; transcript ENST00000304658 | gene_id ENSG00000169180; transcript ENST00000304658; bare_ensg ENSG00000169180; LOEUF 0.04; decile 0; mane_select True; canonical True | pass |
| `loeuf_low_2` | MED13L | LOEUF 0.041; decile 0; transcript ENST00000281928 | gene_id ENSG00000123066; transcript ENST00000281928; bare_ensg ENSG00000123066; LOEUF 0.041; decile 0; mane_select True; canonical True | pass |
| `loeuf_low_3` | EIF4G1 | LOEUF 0.048; decile 0; transcript ENST00000346169 | gene_id ENSG00000114867; transcript ENST00000346169; bare_ensg ENSG00000114867; LOEUF 0.048; decile 0; mane_select True; canonical True | pass |
| `loeuf_near_1_1` | CABLES2 | LOEUF 1.0; decile 5; transcript ENST00000279101 | gene_id ENSG00000149679; transcript ENST00000279101; bare_ensg ENSG00000149679; LOEUF 1; decile 5; mane_select True; canonical True | pass |
| `loeuf_near_1_2` | CFAP69 | LOEUF 1.0; decile 5; transcript ENST00000389297 | gene_id ENSG00000105792; transcript ENST00000389297; bare_ensg ENSG00000105792; LOEUF 1; decile 5; mane_select True; canonical True | pass |
| `loeuf_near_1_3` | CHI3L1 | LOEUF 1.0; decile 5; transcript ENST00000255409 | gene_id ENSG00000133048; transcript ENST00000255409; bare_ensg ENSG00000133048; LOEUF 1; decile 5; mane_select True; canonical True | pass |

## Null rates

Transcript file:

| Column | Nulls | Null rate |
|---|---:|---:|
| `lof.oe_ci.upper` | 19,712 | 9.32% |
| `lof.oe_ci.upper_bin_decile` | 193,654 | 91.55% |
| `lof.pLI` | 19,712 | 9.32% |
| `mis.z_score` | 238 | 0.11% |
| `syn.z_score` | 305 | 0.14% |
| `constraint_flags` | 0 | 0.00% |
| `bare_ensg` | 110,556 | 52.27% |
| `bare_enst` | 110,556 | 52.27% |
| `chromosome` | 110,556 | 52.27% |
| `cds_length` | 110,556 | 52.27% |
| `mane_select` | 0 | 0.00% |
| `canonical` | 0 | 0.00% |
| `gene` | 802 | 0.38% |
| `gene_id` | 0 | 0.00% |
| `transcript` | 0 | 0.00% |
| `lof_hc_lc.obs` | 18,796 | 8.89% |
| `lof_hc_lc.exp` | 18,796 | 8.89% |
| `lof_hc_lc.possible` | 18,796 | 8.89% |
| `lof_hc_lc.oe` | 18,796 | 8.89% |
| `lof_hc_lc.mu` | 18,796 | 8.89% |
| `lof_hc_lc.pLI` | 18,796 | 8.89% |
| `lof_hc_lc.pNull` | 18,796 | 8.89% |
| `lof_hc_lc.pRec` | 18,796 | 8.89% |
| `lof.obs` | 19,712 | 9.32% |
| `lof.exp` | 19,712 | 9.32% |
| `lof.possible` | 19,712 | 9.32% |
| `lof.oe` | 19,712 | 9.32% |
| `lof.mu` | 19,712 | 9.32% |
| `lof.pNull` | 19,712 | 9.32% |
| `lof.pRec` | 19,712 | 9.32% |
| `lof.oe_ci.lower` | 19,712 | 9.32% |
| `lof.oe_ci.upper_rank` | 193,654 | 91.55% |
| `lof.z_raw` | 19,712 | 9.32% |
| `lof.z_score` | 19,712 | 9.32% |
| `mis.obs` | 238 | 0.11% |
| `mis.exp` | 238 | 0.11% |
| `mis.possible` | 238 | 0.11% |
| `mis.oe` | 238 | 0.11% |
| `mis.mu` | 238 | 0.11% |
| `mis.oe_ci.lower` | 238 | 0.11% |
| `mis.oe_ci.upper` | 238 | 0.11% |
| `mis.z_raw` | 238 | 0.11% |
| `mis_pphen.obs` | 49,254 | 23.29% |
| `mis_pphen.exp` | 49,254 | 23.29% |
| `mis_pphen.possible` | 49,254 | 23.29% |
| `mis_pphen.oe` | 49,254 | 23.29% |
| `syn.obs` | 305 | 0.14% |
| `syn.exp` | 305 | 0.14% |
| `syn.possible` | 305 | 0.14% |
| `syn.oe` | 305 | 0.14% |
| `syn.mu` | 305 | 0.14% |
| `syn.oe_ci.lower` | 305 | 0.14% |
| `syn.oe_ci.upper` | 305 | 0.14% |
| `syn.z_raw` | 305 | 0.14% |
| `level` | 110,556 | 52.27% |
| `transcript_type` | 110,556 | 52.27% |
| `num_coding_exons` | 110,556 | 52.27% |

Columns with zero nulls: 5 of 57.

Gene file:

| Column | Nulls | Null rate |
|---|---:|---:|
| `loeuf` | 505 | 2.70% |
| `loeuf_decile` | 849 | 4.54% |
| `lof.oe_ci.upper` | 505 | 2.70% |
| `lof.oe_ci.upper_bin_decile` | 849 | 4.54% |
| `lof.pLI` | 505 | 2.70% |
| `mis.z_score` | 22 | 0.12% |
| `syn.z_score` | 28 | 0.15% |
| `constraint_flags` | 0 | 0.00% |
| `bare_ensg` | 88 | 0.47% |
| `bare_enst` | 91 | 0.49% |
| `chromosome` | 91 | 0.49% |
| `cds_length` | 91 | 0.49% |
| `mane_select` | 0 | 0.00% |
| `canonical` | 0 | 0.00% |
| `gene` | 508 | 2.71% |
| `gene_id` | 0 | 0.00% |
| `transcript` | 0 | 0.00% |
| `lof_hc_lc.obs` | 319 | 1.70% |
| `lof_hc_lc.exp` | 319 | 1.70% |
| `lof_hc_lc.possible` | 319 | 1.70% |
| `lof_hc_lc.oe` | 319 | 1.70% |
| `lof_hc_lc.mu` | 319 | 1.70% |
| `lof_hc_lc.pLI` | 319 | 1.70% |
| `lof_hc_lc.pNull` | 319 | 1.70% |
| `lof_hc_lc.pRec` | 319 | 1.70% |
| `lof.obs` | 505 | 2.70% |
| `lof.exp` | 505 | 2.70% |
| `lof.possible` | 505 | 2.70% |
| `lof.oe` | 505 | 2.70% |
| `lof.mu` | 505 | 2.70% |
| `lof.pNull` | 505 | 2.70% |
| `lof.pRec` | 505 | 2.70% |
| `lof.oe_ci.lower` | 505 | 2.70% |
| `lof.oe_ci.upper_rank` | 849 | 4.54% |
| `lof.z_raw` | 505 | 2.70% |
| `lof.z_score` | 505 | 2.70% |
| `mis.obs` | 22 | 0.12% |
| `mis.exp` | 22 | 0.12% |
| `mis.possible` | 22 | 0.12% |
| `mis.oe` | 22 | 0.12% |
| `mis.mu` | 22 | 0.12% |
| `mis.oe_ci.lower` | 22 | 0.12% |
| `mis.oe_ci.upper` | 22 | 0.12% |
| `mis.z_raw` | 22 | 0.12% |
| `mis_pphen.obs` | 118 | 0.63% |
| `mis_pphen.exp` | 118 | 0.63% |
| `mis_pphen.possible` | 118 | 0.63% |
| `mis_pphen.oe` | 118 | 0.63% |
| `syn.obs` | 28 | 0.15% |
| `syn.exp` | 28 | 0.15% |
| `syn.possible` | 28 | 0.15% |
| `syn.oe` | 28 | 0.15% |
| `syn.mu` | 28 | 0.15% |
| `syn.oe_ci.lower` | 28 | 0.15% |
| `syn.oe_ci.upper` | 28 | 0.15% |
| `syn.z_raw` | 28 | 0.15% |
| `level` | 91 | 0.49% |
| `transcript_type` | 91 | 0.49% |
| `num_coding_exons` | 91 | 0.49% |

Columns with zero nulls: 5 of 59.

## Join

Joined `gene_constraint.bare_ensg` to `D:\variant-workbench\warehouse\intervals\genes.parquet` column `ensg`. Gene rows: 18,716. Rows with a bare ENSG: 18,628. Those missing from intervals: 15. Gene rows with no bare ENSG (RefSeq gene id only): 88. Interval protein-coding genes: 20,097. Matched to a constraint gene row: 18,567. Protein-coding join misses: 1,530.

## File notes

- Header has 55 columns and matches the gnomAD constraint schema in plan.md. No column was renamed.
- LOEUF is `lof.oe_ci.upper`. On the gene file that float is also stored as `loeuf`. The decile is `lof.oe_ci.upper_bin_decile`, also stored as `loeuf_decile`.
- `constraint_flags` is kept as text. `[]` means no flag and is not null.
- `level` is a number 1, 2, or 3, or null. Counts: 1=6,844, 2=89,443, 3=4,680, null=110,556. It is not a gene-versus-transcript split. Every source row is a transcript row.
- Ensembl-id rows: 100,967. Rows whose `gene_id` is not Ensembl (NCBI Gene id on the RefSeq transcript): 110,556.
- On those non-Ensembl rows, `chromosome` is null in 110,556, `cds_length` is null in 110,556, and the LOEUF decile is null in 110,556. The decile used for the truth set is on the Ensembl transcript.
- `chromosome` is passed through `coords.canon_chrom`. Values present: `chr1`, `chr10`, `chr11`, `chr12`, `chr13`, `chr14`, `chr15`, `chr16`, `chr17`, `chr18`, `chr19`, `chr2`, `chr20`, `chr21`, `chr22`, `chr3`, `chr4`, `chr5`, `chr6`, `chr7`, `chr8`, `chr9`. `chrX` and `chrY` are absent when they are not in that list.
- Numeric source fields are float64. `canonical` and `mane_select` are boolean. `NA` tokens are null.
- `bare_ensg` and `bare_enst` come from `scripts/coords.py` and are null when the id is not Ensembl.
