# Release notes

Open-source genetic variant database (GRCh38). Research use only. Built by AI agents (Grok ran the build and wrote the reports; Claude prepared this publication). Source terms apply: see NOTICE.md in the repository.

## Assets

Each asset is a zip that keeps the `warehouse/<folder>/` layout (member paths are `<folder>/<file>`). Every asset is under 1.9 GiB (GitHub's limit is 2 GiB). `SHA256SUMS` has the sha256 of each zip; `MANIFEST.tsv` lists every table with row counts and bytes.

| asset | bytes | contents |
|---|---:|---|
| variant-workbench-clinvar.zip | 562,360,623 | ClinVar variant_summary derived tables (alleles, allele_dedup). |
| variant-workbench-constraint.zip | 32,237,879 | gnomAD v4.1 constraint metrics (gene and transcript). |
| variant-workbench-crosswalk.zip | 3,469,927 | Gene identifier crosswalk (GENCODE, NCBI, HGNC, MANE) and unmatched rows. |
| variant-workbench-gtex.zip | 1,000,184,695 | GTEx v8 open-access tables: median TPM, gene expression, eGenes, sGenes, sQTL pairs (the eQTL pairs table is in its own asset). |
| variant-workbench-gtex-eqtl_pairs.zip | 1,969,181,641 | GTEx v8 eQTL significant pairs (71,478,479 rows), alone because of its size. |
| variant-workbench-gwas.zip | 88,815,755 | GWAS Catalog associations, loci, traits. |
| variant-workbench-intervals.zip | 38,852,112 | GENCODE v49 genes, transcripts, exons, CDS, UTRs. |
| variant-workbench-ledger.zip | 381,884,271 | Evidence ledger (one row per piece of evidence, six classes). |
| variant-workbench-marts.zip | 5,155,988 | Per-gene mart (all genes, protein-coding subset). |
| variant-workbench-mavedb.zip | 907,432,269 | MaveDB experiments, scores, mapped variants, gene coverage (6 score sets held back, see below). |
| variant-workbench-variants.zip | 952,078,910 | Tier A variant table from gnomAD v4.1 (genome, exome, combined), keys, unresolved GWAS rsIDs, field map. Third-party predictor scores removed. |

Total: 5,941,654,070 bytes in 11 assets.

## Exclusions and changes from the working warehouse

- Third-party predictor scores are not included (CADD, REVEL, SpliceAI, Pangolin, phyloP, SIFT, PolyPhen; genome and exome forms, 14 columns). They were removed from `variants/tier_a.parquet` (14 columns), `tier_a_genome.parquet` (7) and `tier_a_exome.parquet` (7), and their 14 rows were removed from `variants/field_map.tsv`. Row order and all other columns are unchanged. Regenerate them locally from gnomAD with `scripts/annotate_gnomad.py`.
- Six MaveDB score sets are held back: urn:mavedb:00000044-a-1, -a-2, -b-1, -b-2 (SARS-CoV-2 RBD, BSD-3 notice) and urn:mavedb:00001255-a-1, -a-2 (AIRE, data must be pre-printed prior to use). They are removed from `experiments` (2,777 -> 2,771 rows), `scores` (8,205,740 -> 7,777,824), `mapped_variants` and `mapped_variants_scored` (6,350,898 -> 6,329,224 each) and `gene_coverage` (1,137 -> 1,135; the two target names used only by those sets are gone). The ledger, the gene marts, the crosswalk and the Tier A tables never contained rows from these sets, so they are unchanged.
- The 28 Fort Lauderdale score sets (PMID 31395865) and the RHO set urn:mavedb:00001275-a-1 are kept; cite the publications named in each set's notice.
- Not published: `raw/`, `db/` (including `variants.sqlite`, rebuildable with `scripts/load_core.py`), `warehouse/tmp`, and the `reports`, `qc` and `manifest` folders (these live in the git repository).
- Tier A `vep_*` columns come from gnomAD's own VEP v105 / GENCODE Release 39 annotation, not from GENCODE v49.

## Source terms

The project's own contributions are CC0 1.0. The source data keep their providers' terms: gnomAD v4.1 (CC0, no re-identification, trademark notice), ClinVar/NCBI (no NCBI restrictions, see the NCBI disclaimer), GENCODE v49 and GWAS Catalog (EMBL-EBI Terms of Use), GTEx v8 open-access (acknowledge GTEx), MaveDB (CC0, per-score-set notices). See NOTICE.md in the repository.

## Using the assets

```
# verify (from the folder with the downloads)
sha256sum -c SHA256SUMS
# unpack everything into one warehouse folder
for z in variant-workbench-*.zip; do unzip -o "$z" -d warehouse; done
```

`gtex/eqtl_pairs.parquet` is in its own asset and is not split; no reassembly is needed. All tables are Parquet except `variants/tier_a_keys.tsv`, `tier_a_unresolved.tsv` and `field_map.tsv`.
