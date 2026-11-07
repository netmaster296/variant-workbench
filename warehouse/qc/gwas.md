# GWAS Catalog

This table is a literature index of reported associations. It is not a file of genome-wide summary statistics, and it does not contain individual genotypes. Do not compute a polygenic score, or any other score of a person, from these rows. Lead SNPs are the associations an author reported, not an effect-size vector that can be summed.

## Source

- Zip: `raw/gwas/gwas-catalog-associations_ontology-annotated-full.zip`
- Member: `gwas-catalog-download-associations-alt-full.tsv`
- Compressed bytes: 72,739,186
- Uncompressed bytes: 749,786,313
- Catalog columns kept: all 38, in file order, as text. Empty cells stay empty strings.
- `DATE ADDED TO CATALOG` ranges from `2008-06-10` to `2026-10-23` (the values are `YYYY-MM-DD`).
- There is no assembly column. Current `CHR_POS` values are GRCh38. An empty or multi-valued `CHR_POS` is an unplaced association (haplotype, interaction, or similar), not a second assembly. Those cells are not split and are not lifted over.

## Columns added

| Column | Rule |
|---|---|
| `chrom` | Canonical chromosome from `CHR_ID` via `canon_chrom`. Null when `CHR_ID` is empty or is not one chromosome (`6 x 6`, `2;2`, and the same kind of cell). |
| `chrom_raw` | `CHR_ID` with surrounding whitespace removed. Null when that is empty. Multi-chromosome text is kept here and does not become `chrom`. |
| `pos` | `CHR_POS` as an integer >= 1. Null when the cell is empty, zero, or not a single integer. Multi-position cells stay null. |
| `rsid` | `norm_rsid(SNP_ID_CURRENT)`. The file stores a bare integer (`4823308` becomes `rs4823308`). Empty ids and suffixed junk such as `7658266b` stay null. The `SNPS` column is not a fallback. |
| `p_value_float` | `P-VALUE` when the whole cell is a plain finite number >= 0. The raw text stays in `P-VALUE`. |
| `effect_float` | `OR or BETA` when the whole cell is a plain finite number. `Infinity`, blanks, and prose stay null. The raw text stays in `OR or BETA`. |
| `risk_allele` | Allele token from `STRONGEST SNP-RISK ALLELE`. `rs128494-C` stores `C`. `rs145611121-GC` stores `GC`. This is not a REF/ALT pair. `?` is null. Several SNP tokens stay in order, joined by `;`. |
| `variant_key` | Null on every row. This file has no ref or alt. Alleles are not guessed from `risk_allele`. A later rsID match to ClinVar or gnomAD is how a row can gain a `chrom:pos:ref:alt` key. |
| `locus_key` | `chrom:pos` when both `chrom` and `pos` exist. Otherwise null. |

Sample size and ancestry stay as prose in `INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE`. They are not parsed.

P-values that underflow float64: 6,279. Those ranks are stored as the smallest positive float so they remain more significant than any representable p-value. Negative p-values left null: 0. Non-representable effects left null: 0.

## Counts

| Measure | Value |
|---|---|
| Association rows | 1,205,644 (equal to the 1,205,644 baseline) |
| Null `pos` | 216,880 / 1,205,644 (0.179887) |
| Null `chrom` | 216,880 / 1,205,644 (0.179887) |
| `chrom` null but `chrom_raw` set (multi-chromosome or other non-canonical `CHR_ID`) | 5,238 |
| `CONTEXT` contains `intergenic` | 213,162 / 1,205,644 (0.176803) |
| `INTERGENIC` flag equal to `1` | 512,000 |
| Distinct non-empty `MAPPED_TRAIT_URI` | 16,885 |
| Trait-table rows, including one empty-URI bucket | 16,886 |
| Loci (clumps) | 311,160 |
| Raw rows per locus | 3.87 |
| Eligible rows | 988,764 |
| Eligible rows per locus | 3.18 |
| Sum of locus `n_associations` | 1,083,056 |
| Singleton loci | 181,401 / 311,160 (0.582983) |
| `CNV` equal to `Y` (kept in associations, excluded from clumps) | 0 |
| Rows with a `locus_key` | 988,764 |
| Rows with an `rsid` | 996,237 |
| Rows with a parsed `effect_float` | 1,015,663 |
| Rows with a `risk_allele` | 878,626 |
| Rows whose `risk_allele` has more than one token | 1,985 |
| Rows with more than one mapped trait URI | 100,074 |
| Rows with an empty `MAPPED_TRAIT_URI` | 9 |
| Multi-URI rows whose trait-name list could not be aligned | 1,395 |
| URIs with more than one mapped name | 0 |
| Non-empty URIs with no mapped name | 4 |

`CONTEXT` and the `INTERGENIC` flag are not the same claim. The intergenic fraction above is the `CONTEXT` substring. Some rows have `INTERGENIC` of `1` and a `CONTEXT` that does not contain `intergenic`.

The sum of locus `n_associations` is larger than the eligible row count when a row lists more than one `MAPPED_TRAIT_URI`: that row is clumped once inside each URI. It is smaller than the raw row count because null positions never enter a locus. File check of that sum: 1,083,056.

Parquet row counts read back from the files: associations 1,205,644, traits 16,886, loci 311,160. Null `pos` read back from associations: 216,880. `CONTEXT` containing `intergenic`, read back: 213,162. `variant_key` nulls read back: 1,205,644.

## Clumping

Within one `MAPPED_TRAIT_URI`, and within one chromosome, associations with a canonical `chrom`, a `pos` >= 1, and `CNV` not equal to `Y` are clumped as follows.

1. The lead is the still-unassigned association with the smallest `p_value_float`. Ties break by smaller position, then `rsid`, then file line.
2. Every other still-unassigned association on that chromosome within 1,000,000 bp of the lead joins the clump.
3. Repeat until none remain.

A missing p-value sorts after every real p-value. This is not single-linkage: a SNP within 1 Mb of a member but more than 1 Mb from the lead stays out and can seed its own clump. Leads on the same trait and chromosome are therefore more than 1,000,000 bp apart, and every member is within 1,000,000 bp of its own lead.

`start` and `end` are the minimum and maximum member positions (1-based, inclusive). They are the span of the collapsed associations, not the ±1 Mb search window. A one-SNP locus has `start` = `end` = `lead_pos`.

Rows with a null position stay in `associations.parquet` and stay out of `loci.parquet`. `CNV` of `Y` would do the same; this release has 0 such rows. Associations with an empty `MAPPED_TRAIT_URI` are clumped together under a null `trait_uri` when they have a position, and they are the empty-URI row in `traits.parquet`.

## Traits and mapped genes

`traits.parquet` has one row per `MAPPED_TRAIT_URI`:

- `trait_uri`
- `mapped_trait`, the ontology label. A URI that is alone on a row keeps the whole `MAPPED_TRAIT` text, commas included (`osteoarthritis, hand`). Several URIs on one row are paired with comma-separated labels only when the counts match. The stored name is the most common label, then the alphabetical label on a tie.
- `n_associations`, catalog rows that list that URI. A row with two URIs counts in both.
- `n_author_traits`, distinct non-empty `DISEASE/TRAIT` strings on those rows.

`loci.parquet` column `mapped_genes` is a struct, not a plain string:

- `mapped_genes.raw` is the lead association's `MAPPED_GENE` text, unchanged. Null when the lead gene is blank.
- `mapped_genes.tokens` is the comma-split token list from every member of the clump, lead first, then the other members in position order, with duplicates removed. Tokens are trimmed. Empty tokens and a token that is exactly `-` are removed. `HLA-DQB1` stays one token.

The intergenic catalog form `FBLN1 - LINC01589` has no comma, so it stays one token. A symbol matches a token, not a substring: `FBLN1` does not match that token. The raw string is kept so a later step can see the unsplit text. `BRCA1` does match a token of `BRCA1, NBR2`.

## Truth set

`warehouse/qc/truth_set.tsv`.

### gwas_intergenic (PASS)

File line 8 (header is line 1). Read back from `associations.parquet` in file order: chrom `chr22`, pos `45604058`, rsid `rs4823308`, locus_key `chr22:45604058`, variant_key `None`, CONTEXT `intergenic_variant`, trait URI `http://www.ebi.ac.uk/efo/EFO_0004348`, CNV `N`, risk_allele `A`.

The position is not null, so the row is outside the null-position set. `CNV` is not `Y`, and both chrom and pos exist, so the row is eligible for a clump. The clump that contains it is trait `http://www.ebi.ac.uk/efo/EFO_0004348`, chrom `chr22`, lead `rs9330814` at `45968311` (p `2e-69`), span `45589577`–`45981128`, members `19`. That locus row is present in `loci.parquet`.

### gwas_empty_pos (PASS)

File line 5. Read back from `associations.parquet`: study `GCST004604`, pubmed `27863252`, SNPS `chr22:37462936`, CHR_POS `''`, chrom `None`, pos `None`, locus_key `None`, variant_key `None`, risk_allele `G`, trait URI `http://www.ebi.ac.uk/efo/EFO_0004348`.

The row is in the association table. `pos` and `locus_key` are null. `STRONGEST SNP-RISK ALLELE` on this row looks like a coordinate plus an allele (`chr22:37462936-G`); only the allele token is stored, and it is not written into `pos`. The clumper never sees this row. Every locus has a non-null `lead_pos` and `chrom` (null lead positions in the locus file: 0), so this association is absent from `loci.parquet`.
