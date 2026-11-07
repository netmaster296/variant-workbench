# variant-workbench

Open-source genetic variant database: a GRCh38 warehouse of genes, studied variants and their evidence from gnomAD, ClinVar, GENCODE, GTEx, the GWAS Catalog and MaveDB. Research use only.

> **Research use only.** This is not medical advice, not a diagnostic tool and
> not for clinical decisions. It is not a personal genome tool: it holds no
> individual genotypes and cannot read anyone's genome. Behavioral and
> cognitive GWAS hits in it are weak observational claims. The reports were
> written by an AI agent and are published as written.

This repository holds the code that builds the warehouse, the build plan, the
QC notes and the written reports. The warehouse tables themselves are
published as GitHub Release assets (see [Release assets](#release-assets)).
The raw downloads (about 768 GB) are not redistributed; the scripts fetch them
from the providers.

## What this is, and is not

- It is a reproducible build: about 20 single-file Python scripts that turn
  six public sources into Parquet tables queried with DuckDB, plus the derived
  tables.
- It models **public evidence** about genes and about variants someone has
  already studied or measured: ClinVar labels, GWAS Catalog associations,
  GTEx expression and QTLs, MaveDB assay scores, gnomAD frequencies and gnomAD
  constraint. Each kind of evidence stays a separate claim.
- It is not a complete map of what the genome does. Most of GRCh38 falls
  outside gene exons, lead QTLs and GWAS loci; the build measures that
  fraction (`warehouse/reports/dark_genome.md`) instead of explaining every
  base.
- It contains no DNA sequence and no individual genotypes, only variant keys,
  annotations and aggregate statistics.
- It is not a polygenic-score tool, a linkage-disequilibrium resource or a
  colocalization analysis (see [Limitations](#limitations)).

## Sources and versions

| Source | Version used | What was used | Terms (details in `NOTICE.md`) |
|---|---|---|---|
| gnomAD | v4.1 | Sites VCFs for genomes and exomes, chr1-22, X and Y (48 VCFs with indexes); v4.1 constraint metrics | Primary data CC0; do not attempt to re-identify participants; "gnomAD" is a trademark of The Broad Institute, Inc. |
| ClinVar (NCBI) | Weekly `variant_summary.txt.gz`, file dated 2026-11-01 | GRCh38 rows | NCBI places no restrictions; NCBI disclaimer applies |
| GENCODE | Release 49 (Ensembl 115) | Comprehensive annotation GTF | Open access; EMBL-EBI Terms of Use |
| NHGRI-EBI GWAS Catalog | `releases/latest` retrieved 2026-11-03; TSV dated 2026-10-29 | Associations, ontology-annotated full download | EMBL-EBI Terms of Use |
| GTEx | v8 (RNA-SeQC v1.1.9) | Median TPM by tissue; single-tissue cis-eQTL and cis-sQTL results | Open-access files; acknowledge GTEx |
| MaveDB | Bulk dump as of 2026-06-24 (`mavedb-dump.2026062418131`, Zenodo DOI [10.5281/zenodo.20840937](https://doi.org/10.5281/zenodo.20840937)) | Score sets, scores, mapped variants | CC0, with per-score-set notes |

All files were downloaded on 2026-11-03. Exact URLs are in
`scripts/fetch_sources.py`. This tool includes data from the gnomAD v4.1
release. This project is not affiliated with or endorsed by gnomAD, the Broad
Institute, NCBI, EMBL-EBI, GENCODE, GTEx or MaveDB.

## Repository layout

```text
plan.md                     the build plan the agent followed (rules, schemas, waves, limits)
scripts/                    the build: fetch, load, build, annotate, QC and report helpers
warehouse/
  reports/                  17 written reports (Markdown), as written by the build agent
    sql/                    the DuckDB SQL behind 16 of them (run with scripts/run_report_sql.py)
    figures/                cognitive_loci_map.png
  qc/                       QC notes per lane, the truth-set report, wave notes
  manifest/                 the raw-file manifest (files.tsv) and disk note, paths relative to ROOT
logs/fetch.status           the download's final status line (DONE)
docs/RELEASE_NOTES.md       the text of the GitHub Release notes
LICENSE                     GPL-3.0-or-later (code)
DATA-LICENSE.md             CC0 1.0 (this project's data contributions)
NOTICE.md                   source terms, required notices, exclusions
CITATIONS.md, CITATION.cff  how to cite this project and its sources
```

Not in git: `raw/` (downloads), `db/` (SQLite catalog), `logs/*.log`, and
every warehouse data table (`.parquet`, `.tsv`). See `.gitignore`.

## Warehouse layout

All tables are Parquet under `warehouse/`, GRCh38, queried with DuckDB.
Column-level detail is in `plan.md` and in each lane's QC note in
`warehouse/qc/`.

| Folder | Tables | Contents |
|---|---|---|
| `intervals/` | `genes`, `transcripts`, `exons`, `cds`, `utrs` | GENCODE v49 gene models, 1-based inclusive coordinates |
| `constraint/` | `transcript_constraint`, `gene_constraint` | gnomAD v4.1 constraint; one row per gene picked by MANE, then canonical, then longest CDS |
| `clinvar/` | `alleles`, `allele_dedup` | ClinVar GRCh38 records with VCF-style variant keys and review-star buckets; one row per variant key |
| `gwas/` | `associations`, `traits`, `loci` | All 38 catalog columns plus parsed fields; traits by EFO URI; associations clumped within 1 Mb per trait |
| `gtex/` | `median_tpm`, `gene_expression`, `egenes`, `sgenes`, `eqtl_pairs`, `sqtl_pairs` | Expression by tissue with tau; lead and significant eQTL and sQTL pairs |
| `mavedb/` | `experiments`, `scores`, `mapped_variants`, `mapped_variants_scored`, `gene_coverage` | Score-set inventory, scores, GRCh38-mapped variants, coverage per gene |
| `crosswalk/` | `genes`, `unmatched` | Gene identifier crosswalk (Ensembl, symbol, NCBI, HGNC, MANE) and join misses |
| `marts/` | `gene_mart`, `gene_mart_protein_coding` | One row per GENCODE gene: identity, constraint, expression, ClinVar, GWAS, MaveDB, QTL counts |
| `variants/` | `tier_a`, `tier_a_genome`, `tier_a_exome`, `tier_a_keys.tsv`, `tier_a_unresolved.tsv`, `field_map.tsv` | Tier A variants (ClinVar, GWAS with a resolved key, MaveDB mapped, lead QTLs) with gnomAD genome and exome frequencies in separate columns and gnomAD's VEP consequences |
| `ledger/` | `evidence` | Evidence ledger: one row per claim, with the evidence class on every row (clinvar, gwas, eqtl, sqtl, mavedb, constraint) |

The `variant_key` is `chrom:pos:ref:alt` with a 1-based VCF position and
uppercase alleles. Ensembl ids are joined without their version. Genome and
exome allele counts are never added together.

## Release assets

The warehouse tables are attached to the GitHub Release as 11 zip files, each
under 2 GiB. Each zip keeps the `warehouse/<folder>/` layout (members are
`<folder>/<file>`), so unpacking them all into one `warehouse` folder under your
ROOT gives the paths the scripts and the report SQL expect. Check the downloads
against `SHA256SUMS`, and `MANIFEST.tsv` lists every table with bytes and rows.
Both `SHA256SUMS` and `MANIFEST.tsv` are attached to the release.

```bash
sha256sum -c SHA256SUMS
for z in variant-workbench-*.zip; do unzip -o "$z" -d warehouse; done
```

PowerShell equivalent (run in the folder that holds the zips and `SHA256SUMS`;
set `$root` to your ROOT, the folder that should contain `warehouse`):

```powershell
# Compare each zip's SHA-256 with SHA256SUMS (prints False for any mismatch)
Get-Content SHA256SUMS | ForEach-Object {
  $hash, $name = $_ -split '\s+', 2
  $name = $name.TrimStart('*').Trim()
  (Get-FileHash -Algorithm SHA256 $name).Hash.ToLower() -eq $hash.ToLower()
}
# Members are <folder>/<file>, so extract into the warehouse folder
$root = 'E:\vw-data'
Get-ChildItem variant-workbench-*.zip | ForEach-Object {
  Expand-Archive -LiteralPath $_.FullName -DestinationPath "$root\warehouse" -Force
}
```

`gtex/eqtl_pairs.parquet` is in its own zip and is not split. All tables are
Parquet except `variants/tier_a_keys.tsv`, `tier_a_unresolved.tsv` and
`field_map.tsv`. The reports, QC notes and manifests are in this repository, not
in the release.

| Asset | Contents (member: rows) | Bytes | SHA-256 |
|---|---|---:|---|
| `variant-workbench-clinvar.zip` | `clinvar/allele_dedup.parquet`: 4,639,253, `clinvar/alleles.parquet`: 4,639,270 | 562,360,623 | `9e4d5846cd3bb81c7b0219845505064f0078bf6fa19ce489e273c9105dc6308f` |
| `variant-workbench-constraint.zip` | `constraint/gene_constraint.parquet`: 18,716, `constraint/transcript_constraint.parquet`: 211,523 | 32,237,879 | `12556c3108b9e638fbd18c900f8ac16b482d779325f42cd0681209cbd846019d` |
| `variant-workbench-crosswalk.zip` | `crosswalk/genes.parquet`: 78,691, `crosswalk/unmatched.parquet`: 11,339 | 3,469,927 | `de765d0015af52009d9af0fcf17f8769b2b0bf9f2acd300dba757b306a899ac9` |
| `variant-workbench-gtex-eqtl_pairs.zip` | `gtex/eqtl_pairs.parquet`: 71,478,479 | 1,969,181,641 | `994bce0b68f9f7d3a6b6931b9ce8f1b34f7b02efe006cf57a5f19ea52d70f3b2` |
| `variant-workbench-gtex.zip` | `gtex/egenes.parquet`: 1,207,976, `gtex/gene_expression.parquet`: 56,200, `gtex/median_tpm.parquet`: 3,034,800, `gtex/sgenes.parquet`: 618,932, `gtex/sqtl_pairs.parquet`: 24,445,206 | 1,000,184,695 | `34cf2c3fa82e7676956adc8303c996cc900ecf94df21b4cf5f176446b59d0d9a` |
| `variant-workbench-gwas.zip` | `gwas/associations.parquet`: 1,205,644, `gwas/loci.parquet`: 311,160, `gwas/traits.parquet`: 16,886 | 88,815,755 | `5ac842c98813b44d5a653e482ed4410ba791b41a661696f16b6a43bab492eca2` |
| `variant-workbench-intervals.zip` | `intervals/cds.parquet`: 2,498,156, `intervals/exons.parquet`: 3,673,949, `intervals/genes.parquet`: 78,691, `intervals/transcripts.parquet`: 507,365, `intervals/utrs.parquet`: 770,688 | 38,852,112 | `6e3979d3adcc84178a47e1e9f63a3f0d39af462785656feba280719fdb981199` |
| `variant-workbench-ledger.zip` | `ledger/evidence.parquet`: 10,866,762 | 381,884,271 | `0d6684151bbe512f81e83bb3ef1e6299cf3515ffa085fb7dd7e5dd4616d3e68a` |
| `variant-workbench-marts.zip` | `marts/gene_mart.parquet`: 78,691, `marts/gene_mart_protein_coding.parquet`: 20,097 | 5,155,988 | `33e6bd8393f8102f39465f150bc16af6f231576b0f6118daec70baee88935161` |
| `variant-workbench-mavedb.zip` | `mavedb/experiments.parquet`: 2,771, `mavedb/gene_coverage.parquet`: 1,135, `mavedb/mapped_variants.parquet`: 6,329,224, `mavedb/mapped_variants_scored.parquet`: 6,329,224, `mavedb/scores.parquet`: 7,777,824 | 907,432,269 | `1d090d66605877192b65479658e150e89793580a5ad9aad7911d334a6b5b85c9` |
| `variant-workbench-variants.zip` | `variants/field_map.tsv`: 43, `variants/tier_a.parquet`: 5,928,738, `variants/tier_a_exome.parquet`: 5,928,738, `variants/tier_a_genome.parquet`: 5,928,738, `variants/tier_a_keys.tsv`: 5,928,738, `variants/tier_a_unresolved.tsv`: 416,792 | 952,078,910 | `e38123251628c451da31af29844f9f601ef83591398ea2c5f5354d5e2c2c4384` |

Total: 5,941,654,070 bytes (5.534 GiB) in 11 assets, each under 2 GiB. Row counts are rows in the table, not lines.

The release tables differ from the original build in two ways (see
[What was left out](#what-was-left-out-and-why)): the 14 predictor-score
columns are removed from the Tier A tables (and their 14 rows from
`field_map.tsv`), and six MaveDB score sets are removed from the MaveDB tables.

## Rebuilding from the sources

Nothing here has to be rebuilt to use the release tables. A full rebuild
downloads about 768 GB and takes a large machine.

### Requirements

- Python 3.13 with the packages in `requirements.txt`
  (`pip install -r requirements.txt`).
- Disk: about 768 GB for the downloads (761.6 GB of it is the gnomAD VCFs),
  about 7 GB for the warehouse and about 1 GB for the SQLite catalog.
  `fetch_sources.py` refuses a download that would leave less than 80 GiB
  free, and `load_gtex.py` stops before writing the QTL pair tables if its
  size estimate would leave less than 120 GiB free. Plan on roughly 1 TB free.
- Memory: the heaviest steps (`build_tier_a.py`, `annotate_gnomad.py`,
  `build_gene_mart.py`, `build_ledger.py`) set DuckDB to `memory_limit='64GB'`
  and 8 threads; the report SQL uses 12-16 GB. The original build ran on a
  Windows 11 desktop with about 127 GB of RAM. On a smaller machine, lower
  those settings; DuckDB spills to `warehouse/tmp/`. Peak memory and run times
  were not measured.
- No `bcftools`, `htslib` or `pysam` needed: `scripts/vcf_query.py` and
  `scripts/vcf_stream.py` read the bgzip VCFs and tabix indexes in pure Python
  (with Biopython's BGZF reader).

### Where the data go: ROOT and VW_ROOT

Every script reads and writes under one folder, ROOT. ROOT is the `VW_ROOT`
environment variable when it is set; otherwise it is the folder above
`scripts/`, so by default the data go inside your clone. Inside ROOT the
scripts use `raw/`, `db/`, `logs/` and `warehouse/`.

```powershell
# Windows PowerShell: keep the code in the clone, put the data on a big drive
$env:VW_ROOT = 'E:\vw-data'
```

```bash
# Linux or macOS
export VW_ROOT=/data/vw
```

The scripts import each other from `scripts/`, so run them as
`python scripts/<name>.py` from anywhere; they also work under `python -I`.

### Windows notes

The project was built and run on Windows 11 only.

- `fetch_sources.py` calls `curl.exe`, which ships with Windows 10 and later.
  On Linux or macOS, change `"curl.exe"` to `"curl"` in that script (untested).
- The pure-Python VCF readers should run on Linux and macOS but were not
  tested there.
- `plot_cognitive_loci.py` uses Arial from `C:\Windows\Fonts` and falls back
  to Pillow's default font elsewhere.
- Disk-space checks use ROOT's drive on Windows and the filesystem that holds
  ROOT elsewhere.

### Build order

This is the order in `plan.md` (Step 0 to Step 4). Each step reads the outputs
of the steps before it.

```bash
# Step 0: download (about 768 GB) into raw/; "core" and "gnomad" run the two halves
python scripts/fetch_sources.py all          # writes logs/fetch.log and logs/fetch.status (must end DONE)

# Wave 0: catalog, manifest, truth set
python scripts/load_core.py                  # db/variants.sqlite (its gene table is used by later steps)
python scripts/build_manifest.py             # warehouse/manifest/files.tsv and disk.md; SQLite source_file
python scripts/build_truth_set.py            # warehouse/qc/truth_set.tsv (18 cases looked up in the sources)

# Wave 1: one loader per source (independent of each other)
python scripts/load_intervals.py             # warehouse/intervals/
python scripts/load_constraint.py            # warehouse/constraint/
python scripts/load_clinvar.py               # warehouse/clinvar/
python scripts/load_gwas.py                  # warehouse/gwas/
python scripts/load_gtex.py                  # warehouse/gtex/
python scripts/load_mavedb.py                # warehouse/mavedb/ (inventory and mapped variants)
python scripts/load_mavedb.py --scores-only  # warehouse/mavedb/ scores joined to mapped variants
python scripts/build_crosswalk.py            # warehouse/crosswalk/ (after intervals, constraint, ClinVar)
python scripts/truth_set.py                  # warehouse/qc/truth_set_report.md

# Wave 2: integration
python scripts/build_gene_mart.py            # warehouse/marts/
python scripts/build_tier_a.py               # warehouse/variants/tier_a_keys.tsv, tier_a_unresolved.tsv
python scripts/annotate_gnomad.py            # warehouse/variants/tier_a*.parquet (one pass over the 48 VCFs)
python scripts/build_ledger.py               # warehouse/ledger/evidence.parquet and reports/dark_genome.md
python scripts/truth_set.py                  # rerun; all 18 cases passed in the original build

# Wave 3: reports (optional)
python scripts/run_report_sql.py <report>    # see "Running the report SQL"
python scripts/plot_cognitive_loci.py        # warehouse/reports/figures/cognitive_loci_map.png
```

Each loader writes a QC note to `warehouse/qc/`. Rebuilding overwrites the QC
notes, `dark_genome.md` and the figure in your ROOT.

A rebuild from the sources reproduces the **original** tables, including the
predictor columns and the six MaveDB score sets this release leaves out. If
you rebuild, the terms in `NOTICE.md` for those data apply to you.

### Regenerating the predictor columns locally

The release removes CADD, REVEL, SpliceAI, Pangolin, phyloP, SIFT and PolyPhen
(genome and exome forms, 14 columns) from the Tier A tables. They come from
the gnomAD v4.1 sites VCFs, and `scripts/annotate_gnomad.py` extracts them
again:

1. Download the gnomAD VCFs and indexes into ROOT:
   `python scripts/fetch_sources.py gnomad` (about 761.6 GB).
2. Put `warehouse/variants/tier_a_keys.tsv` in place: it is in
   `variant-workbench-variants.zip`, or rebuild it with `build_tier_a.py`
   (which needs the ClinVar, GWAS, MaveDB and GTEx tables).
3. Run `python scripts/annotate_gnomad.py`. It writes `tier_a_genome.parquet`,
   `tier_a_exome.parquet` and `tier_a.parquet` with every extracted field,
   predictors included, and `field_map.tsv`. It overwrites the release copies
   of those files in your ROOT.

The regenerated scores are under their providers' terms (CADD, REVEL and
SpliceAI are free for non-commercial use only; see `NOTICE.md`). Do not
redistribute them without checking those terms.

### Running the report SQL

The report SQL files are published exactly as the build wrote them, so they
hard-code the original project folder `D:/variant-workbench`. Do not edit
them; run them through the helper, which substitutes your ROOT at run time:

```bash
python scripts/run_report_sql.py study_redundancy             # prints each result table
python scripts/run_report_sql.py ancestry_audit --max-rows 0  # print every row
python scripts/run_report_sql.py splicing --dry-run           # show the rewritten SQL only
```

Expect differences against the release tables: `conflicts.sql`,
`predictor_disagreement.sql` and `splicing.sql` read the predictor columns
and fail until you regenerate them, and counts that involve MaveDB change
slightly because six score sets are excluded. The report text itself was
written by the agent from these results; no script regenerates it.

## What was left out, and why

- **Third-party predictor scores (14 columns).** CADD, REVEL and SpliceAI
  scores are offered for non-commercial use only, and none of their terms
  grants redistribution in an unrestricted public release. The terms for
  Pangolin, phyloP, SIFT and PolyPhen-2 values as distributed by gnomAD were
  not found, so they are left out as unclear. Regenerate them locally (above).
- **Six MaveDB score sets.** `urn:mavedb:00000044-a-1`, `-a-2`, `-b-1` and
  `-b-2` are SARS-CoV-2 receptor-binding-domain scans (viral, not human) whose
  usage note is a BSD 3-Clause licence text. `urn:mavedb:00001255-a-1` and
  `-a-2` (AIRE) say the data "must be pre-printed prior to usage", and no
  preprint is listed. They are removed from the published MaveDB tables
  (`experiments` 2,777 to 2,771 rows, `scores` 8,205,740 to 7,777,824,
  `mapped_variants` and `mapped_variants_scored` 6,350,898 to 6,329,224 each,
  `gene_coverage` 1,137 to 1,135). The evidence ledger, the gene marts, the
  crosswalk and the Tier A tables never held rows from these sets (checked: 0
  ledger rows reference them), so they are unchanged. The 28 Fort
  Lauderdale score sets (PMID 31395865) and the RHO score set
  `urn:mavedb:00001275-a-1` are kept, and their publications are cited in
  `CITATIONS.md`.
- **`raw/` and `db/`.** The raw downloads are too large to host, and they
  include the GTEx sample-level TPM file (17,382 donor-sample columns) and
  MaveDB metadata with researcher names. The SQLite catalog is a derived
  convenience copy; `load_core.py` and `build_manifest.py` rebuild it.
- **Local details.** Personal paths, hardware model names and process IDs were
  redacted from `plan.md` and `warehouse/qc/wave0.md`, and the manifest paths
  were made relative to ROOT. The substance of the plan is unchanged.

## Limitations

From `plan.md` ("What this disk can and cannot support") and the build's QC:

- **Not every base is explained.** Most of GRCh38 falls outside gene exons,
  lead QTLs and GWAS loci. `dark_genome.md` publishes the fraction.
- **No linkage disequilibrium and no real polygenic score.** gnomAD sites VCFs
  carry no individual genotypes, and the GWAS Catalog stores reported
  significant associations, not genome-wide summary statistics. Do not build a
  score from lead SNPs.
- **No personal genome report.** There are no individual genotypes anywhere in
  the project.
- **No statistical colocalization.** Overlap between a GWAS mapped gene and an
  eGene in a relevant tissue is a hypothesis; colocalization needs full summary
  statistics.
- **Not a complete functional map.** MaveDB covers a small slice of proteins
  (in the original build, assays for 134 of 20,097 protein-coding genes by
  symbol).
- **GWAS rows carry no alleles.** The GWAS Catalog gives no ref/alt, so
  `variant_key` is null in `gwas/associations`; GWAS variants enter Tier A only
  through an rsID match to ClinVar. 416,791 GWAS rsIDs with a position stayed
  unresolved (`tier_a_unresolved.tsv`) because the genome-wide gnomAD rsID scan
  was deferred (416,792 rows in `tier_a_unresolved.tsv`; one rsID appears twice).
  GWAS-versus-frequency results cover only a small share of GWAS
  variants.
- **Ancestry.** Most catalog discoveries come from European-ancestry samples.
  A frequency in the discovery group is not a frequency for everyone. Reports
  give coverage and missingness and do not rank groups.
- **Annotation versions differ.** Gene models are GENCODE v49, but the
  `vep_*` consequence columns come from gnomAD's own VEP v105 run on
  GENCODE 39 transcripts.
- **Sources change.** ClinVar classifications are revised weekly; this is a
  snapshot.
- **One claim was stopped.** `constraint_versus_function.md` does not estimate
  a MaveDB effect size, because the span join credited scores to overlapped
  genes rather than assayed genes (`warehouse/qc/issues_constraint_versus_function.md`).
  `plan.md` asked for a `warehouse/qc/issues.md`; the build did not write one.
- **Other gaps.** The plan's `bcftools` acceptance check was replaced by the
  pure-Python reader because `bcftools` has no win-64 conda build. The build also
  found an extra GENCODE file in `raw/gencode/` (`pc_translations`) that no
  script fetched or used.

## The design rule (rule 11)

`plan.md` rule 11 governs the whole project. In plain words:

- The warehouse describes public biology. It is not a tool for ranking people,
  and it must not be used to score anyone.
- A behavioral or cognitive GWAS association is not a reading of a person's
  nature. Those hits stay in the catalog only as weak observational claims,
  shaped by who was studied, by environment and by social confounding.
- Population-frequency work reports coverage and missingness. It does not rank
  groups or infer disease risk for a population from a lead SNP.

## Reading the reports

The 17 reports in `warehouse/reports/` are AI-generated research notes,
published exactly as written. During publication, 87 of their numbers
(spread across all 17 reports and the figure) were spot-checked against the
warehouse and all matched; the rest were not re-derived. Read them with
these points in mind:

- `common_disease.md` lists cognitive, behavioral and social traits
  (intelligence, cognitive function, income, risk-taking and others) next to
  diseases, under a "Common-disease" heading, and does not repeat rule 11's
  caveat. Those rows are weak observational associations, not disease
  findings and not statements about people.
- `figures/cognitive_loci_map.png` maps catalog lead SNPs for educational
  attainment, intelligence and cognitive function. A highlighted region is
  where a study reported a lead variant, not a segment that causes a trait.
  The discovery cohorts were mostly of European ancestry, the associations are
  confounded by environment and social factors, and they do not predict
  anything about an individual or a group.
- `predictor_disagreement.md` uses phrases such as "the one to trust first".
  That is a statement about the warehouse's evidence, not guidance for
  interpreting anyone's variant.
- `clinvar_versus_frequency.md` and `conflicts.md` discuss common alleles with
  pathogenic labels (for example HBB, HFE, CFTR and G6PD). Carriers of one copy
  of a recessive pathogenic allele are usually unaffected.
- One external citation in `rare_disease.md` (gnomAD LOEUF guidance dates)
  could be only partly verified.
- The numbers describe the original build, which still included the
  predictor columns and the six excluded MaveDB score sets.

## Disclaimers

- **Research use only.** Not medical advice, not a diagnostic tool, and not
  for clinical decisions about you or anyone else. If you have a health
  question, ask a qualified clinician.
- **Not a personal genome tool.** It holds no individual genotypes and is not
  designed to interpret a personal genome.
- **Weak evidence stays weak.** Behavioral and cognitive GWAS hits are weak
  observational claims. Nothing here ranks people or groups.
- **No warranty.** The code is provided under GPL-3.0 without warranty (see
  `LICENSE`), and the data as is. Verify anything you rely on against the
  original sources.

## How this was made: AI assistance

This project was built by AI agents.

- A **Grok** agent ran the build on 2026-11-03 and 2026-11-04, following
  `plan.md`: it wrote the scripts, downloaded the sources, built the warehouse,
  ran the QC and truth set, and wrote the 17 reports and the figure.
- **Claude** (Anthropic) prepared the publication: a licence, privacy and
  ethics audit; the portable paths (`VW_ROOT`) and `run_report_sql.py`; the
  redactions; and this README, the licence files, `NOTICE.md` and the
  citations.

The owner chose the description and the licences; the exclusions follow the
audit's licence and privacy findings. The code and reports have not been independently audited by a
human expert. Treat them accordingly.

## Licences

- **Code** (`scripts/`, `warehouse/reports/sql/`): GNU General Public License
  version 3 or, at your option, any later version (`LICENSE`).
- **This project's data contributions** (schemas, derived keys, aggregates,
  reports, QC notes and documentation): CC0 1.0 Universal (`DATA-LICENSE.md`).
- **Source data** keep their providers' terms (`NOTICE.md`). In particular,
  users of the gnomAD-derived values agree not to attempt to re-identify gnomAD
  participants.

```text
variant-workbench
Copyright (C) 2026 netmaster296

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

## Citing

See `CITATIONS.md` and `CITATION.cff`. Please cite the data sources as well:
gnomAD, ClinVar, GENCODE, the GWAS Catalog, GTEx (with its acknowledgement)
and MaveDB, plus the original publications of the MaveDB score sets you use.
