# MaveDB QC

Dump title: MaveDB public data. asOf: 2026-06-24T18:13:15.647594+00:00.
Archive: `D:\variant-workbench\raw\mavedb\mavedb-dump.2026062418131.tar.gz`.
The tar was streamed. It was not extracted.

## AppleDouble

Skipped 6943 members whose names contain `/._` or whose filename starts with `._`. Their bytes were not read.

## Experiments and score sets

Experiment sets: 1254. Experiments: 1992. Score sets: 2777.
Mapped score sets (a `mapped/*.mapped-variants.json` file is in the tar): 1170.
Unmapped score sets (no mapped file): 1607.
Experiments with at least one mapped score set: 918.
Experiments with no mapped score set: 1074.
Mapped JSON files read: 1170.

MaveDB has no assay-type field. `assay_text` is `methodText`.

Assay text source:

- scoreSets.methodText: 2746
- experiments.methodText: 23
- missing: 8

mappingState versus mapped file:

- None|no_file: 3
- complete|file: 831
- failed|file: 60
- failed|no_file: 1604
- incomplete|file: 279

Score sets whose `datasetColumns.scoreColumns` lists `score`: 2777/2777.
That name is metadata on the score set. It is not a value in the mapped JSON.

## Scores

Mapped JSON has no score field. Scores were loaded from `csv/*.scores.csv`. The join key both sides share is the CSV `accession`, the same string as the mapped variant urn. Scores are not normalized across score sets. Counts are in Score join.

## Mapped records

Raw records in mapped JSON files: 3212230.
`current` is not true, excluded: 0.
Exact duplicate records of the same variant urn in one file, excluded: 11898.
Same variant urn with a different postMapped payload, kept: 22550.
Records kept: 3200332.
Allele rows written: 6350898. A haplotype is one record and one row per member allele.
Records with postMapped null: 184159. These stay in the table with no coordinate.
Records placed on a GRCh38 chromosome: 715755. Of those, 649568 have a variant_key and 66187 have a chromosome but no key.
Mapped records that are non-GRCh38 (postMapped present, and no allele label is a GRCh38 chromosome accession): 2300418.
- protein labels (`NP_`): 2294797
- no sequence label: 5621
- other labels: 0

Allele labels:

- on a GRCh38 ChromosomeAccession: 3590197
- `NP_` protein, variant_key left null: 2570921
- missing label: 5621
- any other label, variant_key left null: 0

A chromosome is built only when the label is a GRCh38 `ChromosomeAccession` in ClinVar `variant_summary.txt.gz` (`Assembly == GRCh38`). The VCF position is `start + 1` when `end - start` equals the `vrs_ref_allele_seq` length. `variant_key` comes from `coords.variant_key`, so an empty reference or alternate stays null. No anchor base was added.

GRCh38 alleles whose span matched but `variant_key` was still null: 75553 (empty ref 53537, empty alt 22016).
GRCh38 alleles whose span did not equal the reference length, position left null: 0.
Accessions whose ClinVar chromosome did not canonicalize: 0.
Ambiguous `vrs_ref_allele_seq` values: 0. Extra copies of that extension: 0.
Variant urns that do not start with the score-set urn from the filename: 19100. These are ids stored on the mapped record. The examples are temporary ids, not a filename parse failure.
Examples:

- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#1`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#2`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#3`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#35`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#36`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#37`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#38`
- `urn:mavedb:00001239-a-6` vs `tmp:992af924-8d46-4bc4-8f9f-f8b089677731#39`

Mapped files with no main.json score set: 0.
Mapped JSON members that were not arrays: 0.
Mapped records that were not objects: 0.

ClinVar GRCh38 accessions with more than one chromosome value. The majority value was used:

- NC_000024.10: {'Y': 3305, 'X': 1}; used Y

Multi-target score sets are not copied onto every gene. Their mapped records stay on the score-set row with the joined target names.
Multi-target score sets that had a mapped file: 0. Records in those files: 0.

## Null rates

mapped_variants rows: 6350898.
- variant_key null: 2836254/6350898 (44.6591%)
- chrom null: 2760701/6350898 (43.4695%)
- pos null: 2760701/6350898 (43.4695%)
- ref null: 184174/6350898 (2.9000%)
- alt null: 184174/6350898 (2.9000%)
- protein_change null: 3774357/6350898 (59.4303%)
- sequence_label null: 189780/6350898 (2.9882%)
- target_gene null: 0/6350898 (0.0000%)

experiments rows: 2777.
- ensembl_id null: 2727/2777 (98.1995%)
- assay_text null: 8/2777 (0.2881%)
- target_name null: 0/2777 (0.0000%)
- mapping_state null: 3/2777 (0.1080%)

## Gene coverage

Distinct raw target names: 1137.
Names MaveDB categories as `protein_coding`: 1101. This category includes phrases and non-human targets. It is not a GENCODE biotype.
Names that are Homo sapiens: 685.
Target names that exactly match a GENCODE protein-coding symbol: 134.
Of those, Homo sapiens: 126.
GENCODE protein-coding genes in `db/variants.sqlite` (`biotype == protein_coding`): 20097.
Protein-coding gene ids touched by an exact symbol match: 134.
Protein-coding gene ids with no MaveDB score set under that rule: 19963.
Target names with no exact protein-coding symbol match: 1003.
- names that contain a space or punctuation and do not match: 701
- names that match a GENCODE symbol of another biotype: 1

The unmatched human protein-coding gene list is the rest of the catalog. It is not written out. MaveDB covers a small slice of proteins.

Sample of target names that are not protein-coding symbols:

- hYAP65 WW domain
- BRCA1 RING domain
- E4B
- ALDOB enhancer
- ECR11 enhancer
- LTV1 enhancer
- SUL1 promoter
- PAB1
- HSP90
- Gal4
- BCL11A enhancer
- F9 promoter
- FOXE1 promoter
- GP1BB promoter
- HBB promoter
- HBG1 promoter
- HNF4A promoter
- IRF4 enhancer
- IRF6 enhancer
- LDLR promoter

Sample of names that match a non-protein-coding GENCODE symbol:

- IGHG1: IG_C_gene

Target genes with at least one mapped variant record: 651.
Target genes with any GRCh38 genomic mapping: 141. This flag is true when a chromosome was built, including alleles that have no variant_key.
Target genes with any variant_key: 141.
Sum of per-gene mapped-variant counts: 3200332. This counts kept variant records for score sets with one target name, not allele rows.

## Truth set

Mapped case `mavedb_mapped` expects variant `urn:mavedb:00000352-0-1#464` to produce `chr5:146469691:A:T`.
Hits written for that urn: 1.
- label NC_000005.10, chrom_raw 5, chrom chr5, pos 146469691, ref A, alt T, key chr5:146469691:A:T, target TCERG1 (Y438R mutant)

Unmapped case `mavedb_unmapped` expects score set `urn:mavedb:00000004-a-1` gene `E4B` to stay without a coordinate.
On-disk checks: mapped key ok True; target name ok True; unmapped score set absent from mapped_variants (0 rows) and experiments say no file and gene E4B coverage [{'target_gene': 'E4B', 'n_mapped_variants': 0, 'has_genomic_mapping': False}]; parquet rows variants 6350898, experiments 2777, genes 1137.

## Result

PASS

## Score join

Mapped JSON still has no score. Scores were read from `csv/*.scores.csv` without extracting the tar.
AppleDouble members skipped in this pass: 6943. Score files read: 2777. Files missing `accession` or `scores.score`: 0.

The column both sides share is `accession` in the CSV and `urn` in `mapped_variants.parquet` (the mapped `variantUrn`). It is not HGVS. For `urn:mavedb:00000352-0-1#464` the CSV `hgvs_pro` is `p.Asn22Ile` and the mapped `protein_change` is null, while `accession` equals the mapped urn.
Score set `urn:mavedb:00001239-a-6` uses `tmp:` ids on both sides.
That set has 19100 mapped urns, 19100 score rows, and 19100 accessions equal to a mapped urn.
README says score values are not normalized across score sets. No normalized score was invented.

`scores.score` is a float only when the cell is a plain number. Numeric: 7999960. Not numeric, float left null, raw text kept: 205780.
Files with a sigma column (`scores.sigma` or `scores.Sigma`): 581.
Files with a raw score column (`scores.raw_score` or `scores.score_raw`): 573.
Files with `scores.raw_sigma`: 521.
Files with other `scores.*` columns, stored as JSON in `extra_scores`: 2118.
Those other columns are not given their own fields. Empty accession cells: 0.

Score rows: 8205740. Distinct `(score_set_urn, accession)` keys: 8205740.
Duplicate keys: 0. Duplicate keys whose `scores.score` text disagrees: 0.
Keys used for the attach (unique, or duplicates with identical score text): 8205740. Conflicting keys were not attached.

Score rows whose accession equals a mapped urn in the same score set: 3177782.
Score rows with no mapped urn: 5027958.
Distinct score keys that match a mapped urn: 3177782.
Distinct score keys with no mapped urn: 5027958.

Score rows in score sets that have a mapped file: 3180836.
Of those, accession matched a mapped urn: 3177782.
Of those, accession missed: 3054.

Mapped rows: 6350898. Rows with a score cell attached: 6350898. Rows with no score attached: 0.
Of those attached cells, 6204960 are plain numbers and 145938 are not, so `score` is null and `score_text` keeps the cell.
Distinct mapped urns: 3177782. Urns with a score attached: 3177782. Urns with no score: 0.
A miss leaves `score` and `score_text` null. It is not filled from another column.

`mapped_variants.parquet` was not rewritten (6350898 rows). `mapped_variants_scored.parquet` has the same rows and the same columns, then the score columns. Scored rows: 6350898.

Truth variant `urn:mavedb:00000352-0-1#464` in the scored file: [('chr5:146469691:A:T', -0.1538729821505993, '-0.1538729821505993', 1)].
CSV score text seen while reading: -0.1538729821505993.
Unmapped score set `urn:mavedb:00000004-a-1` (E4B) has 98298 score rows (98297 numeric) and 0 mapped rows. It still has no coordinate.
