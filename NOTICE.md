# Notices: sources, terms and exclusions

variant-workbench is built from six public data sources. This project's own
code is GPL-3.0 (`LICENSE`) and its own data contributions are CC0 1.0
(`DATA-LICENSE.md`). **Source data keep their providers' terms.** This file
records those terms, quoted exactly from the providers' pages as they were
fetched during the publication audit in November 2026, and lists what was left
out of the release and why. Terms change; check the linked pages before you
rely on them.

This tool includes data from the gnomAD v4.1 release.

## gnomAD v4.1

Used: the v4.1 sites VCFs (genomes and exomes, chr1-22, X and Y) and the v4.1
constraint metrics TSV, from
`https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/`.
Browser: <https://gnomad.broadinstitute.org/>.

Terms, quoted from the gnomAD browser's terms page
(source: <https://github.com/broadinstitute/gnomad-browser/blob/main/browser/about/policies/terms.md>):

> All users of gnomAD data agree to not attempt to reidentify participants.

> The primary data from the gnomAD exomes and genomes are available free of
> restrictions under the Creative Commons Zero Public Domain Dedication. This
> means that you can use it for any purpose without legally having to give
> attribution. However, we request that you actively acknowledge and give
> attribution to the gnomAD project, and link back to the relevant page,
> wherever possible.

> Some annotations may have restrictions on usage. For instance, SpliceAI
> annotations have been computed by Illumina and are provided with permission
> under a CC BY NC 4.0 license for academic and non-commercial use. It is the
> responsibility of users to abide by all relevant licensing requirements.

> We request that developers integrating gnomAD data in their tools include a
> statement acknowledging the inclusion of gnomAD data (e.g., "This tool
> includes data from the gnomAD v4.1 release."). However, to avoid confusion
> and misattribution, we ask that you refrain from incorporating "gnomAD" or
> "Genome Aggregation Database'' into the name of your tool and from using the
> gnomAD logo without permission. gnomAD™ is a trademark owned by The Broad
> Institute, Inc.

What this means here:

- **Anyone who uses the gnomAD-derived values in this release takes on the
  same term: do not attempt to re-identify gnomAD participants.**
- gnomAD™ is a trademark owned by The Broad Institute, Inc. This project is not
  affiliated with or endorsed by the gnomAD project or the Broad Institute. It
  does not use "gnomAD" in its name and does not use the gnomAD logo.
- The CC0 sentence names "the primary data from the gnomAD exomes and genomes".
  It does not name the constraint metrics explicitly; treating them the same
  way is a reasonable reading, not a stated grant.
- The `vep_*` columns in the Tier A tables come from gnomAD's own annotation.
  The VCF headers record `##vep_version=v105` and `##gencode_version=Release 39`,
  so those consequences follow GENCODE 39 transcripts, not the GENCODE v49 gene
  models used elsewhere in the warehouse.
- Third-party predictor scores that gnomAD distributes in the same VCFs are not
  in this release (see "Excluded from this release").

## ClinVar (NCBI)

Used: the weekly `variant_summary.txt.gz` from
`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/` (gzip timestamp
2026-11-01), GRCh38 rows only.

NCBI's "Molecular Data Usage" policy, quoted from
<https://www.ncbi.nlm.nih.gov/home/about/policies/>:

> NCBI itself places no restrictions on the use or distribution of the data
> contained therein. Nor do we accept data when the submitter has requested
> restrictions on reuse or redistribution. However, some submitters of the
> original data (or the country of origin of such data) may claim patent,
> copyright, or other intellectual property rights in all or a portion of the
> data (that has been submitted). NCBI is not in a position to assess the
> validity of such claims and since there is no transfer of rights from
> submitters to NCBI, NCBI has no rights to transfer to a third party.
> Therefore, NCBI cannot provide comment or unrestricted permission concerning
> the use, copying, or distribution of the information contained in the
> molecular databases.

NCBI disclaimer: see "Website and Data Usage Policies and Disclaimers" at
<https://www.ncbi.nlm.nih.gov/home/about/policies/>. From its medical
disclaimer:

> It is not the intention of NLM to provide specific medical advice but rather
> to provide users with information to better understand their health and their
> diagnosed disorders. Specific medical advice will not be provided, and NLM
> urges you to consult with a qualified physician for diagnosis and for answers
> to your personal questions.

Notes: the ClinVar-specific documentation pages returned a bot check during
the audit, so ClinVar-specific terms beyond NCBI's general policy were not
seen. ClinVar `PhenotypeIDS` carry OMIM identifiers only; no OMIM text is
included.

## GENCODE v49

Used: `gencode.v49.annotation.gtf.gz` (GENCODE release 49, Ensembl 115) from
`https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/`.

From <https://www.gencodegenes.org/pages/data_access.html>:

> All the GENCODE project data is open access

GENCODE is distributed by EMBL-EBI under the EMBL-EBI Terms of Use (below).

## NHGRI-EBI GWAS Catalog

Used: `gwas-catalog-associations_ontology-annotated-full.zip` from
`https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/`, downloaded
2026-11-03. The TSV inside is dated 2026-10-29.

From <https://www.ebi.ac.uk/gwas/docs/about>:

> The NHGRI-EBI GWAS Catalog and all its contents are available under the
> general Terms of Use for EMBL-EBI Services. Summary statistics are made
> available under CC0 unless otherwise stated. We advise consumers of data
> hosted by the GWAS Catalog to note the license terms of individual datasets,
> if applicable to their specific use case.

This release uses the associations download, not summary statistics, so the
general EMBL-EBI Terms of Use apply. The association rows include public
bibliographic fields (first author, PubMed id, study accession).

## EMBL-EBI Terms of Use (GENCODE and the GWAS Catalog)

From <https://www.ebi.ac.uk/about/terms-of-use/>:

> Where EMBL-EBI presents scientific data generated by others, EMBL-EBI imposes
> no additional restriction on the use of the contributed data than those
> provided by the data owner, unless otherwise specified in these Terms of Use.
> EMBL-EBI expects attribution (e.g. in publications, services or products) for
> the use of any of its Data Resources and Tools in accordance with good
> scientific practice.

> EMBL-EBI itself places no additional restrictions on the use or
> redistribution of the data available via its Data Resources and Tools other
> than those provided by the original data owners, unless otherwise specified
> in these Terms of Use. The original data may be subject to rights claimed by
> third parties, including but not limited to, patent, copyright, other
> intellectual property rights and Access and Benefit Sharing

## GTEx v8

Used: open-access v8 files from the public `adult-gtex` Google Cloud bucket,
downloaded 2026-11-03: the median TPM by tissue
(`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`) and the
single-tissue cis-eQTL and cis-sQTL tars. The sample-level TPM file was
downloaded but never loaded, and it is not published. Raw sequence data and
protected phenotypes are only on dbGaP and were not used.

No explicit licence for the open-access files was found (the GTEx portal pages
are rendered by JavaScript and could not be read in full). GTEx asks to be
acknowledged. From <https://sites.google.com/broadinstitute.org/gtex-faqs/home>:

> Please acknowledge the GTEx Project and/or Portal. An example acknowledgement
> statement follows: The Genotype-Tissue Expression (GTEx) Project was
> supported by the Common Fund of the Office of the Director of the National
> Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The
> data used for the analyses described in this manuscript were obtained from:
> [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP
> accession number phs000424.vN.pN on MM/DD/YYYY.

**GTEx acknowledgement for this project:** The Genotype-Tissue Expression
(GTEx) Project was supported by the Common Fund of the Office of the Director
of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and
NINDS. The data used here were obtained from the GTEx v8 open-access files in
the GTEx public Google Cloud bucket (`adult-gtex`) on 11/03/2026.

## MaveDB

Used: the MaveDB bulk dump `mavedb-dump.2026062418131.tar.gz` (as of
2026-06-24T18:13:15Z), Zenodo record 20840937, DOI
[10.5281/zenodo.20840937](https://doi.org/10.5281/zenodo.20840937). The dump's
`LICENSE.txt` is the CC0 1.0 Universal legal code, and every score set in this
release has the licence short name `CC0` in the dump's `main.json`.

MaveDB asks to be cited, together with each dataset's original publication.
From <https://mavedb.org/docs/mavedb/citation.html>:

> When citing data downloaded from MaveDB, please refer to the MaveDB
> publication (above) as well as the original publications where the datasets
> were reported. These citations are usually available on the associated
> experiment and score set pages.

The original publication for each score set is listed on its MaveDB page
(search the score-set URN at <https://mavedb.org/>) and in `main.json` in the
dump. Researcher names and ORCID identifiers from `main.json` are not
republished in this release.

Some CC0 score sets carry an extra data-usage note in `main.json`. The notes
and how this release handles them:

### Fort Lauderdale notice: 28 score sets, kept

Score sets: `urn:mavedb:00000014-a-1`, `00000015-a-1`, `00000016-a-1`,
`00000017-a-1`, `00000018-a-1`, `00000019-a-1`, `00000020-a-1`,
`00000021-a-1`, `00000023-a-1`, `00000023-a-2`, `00000024-a-1`,
`00000025-a-1`, `00000026-a-1`, `00000027-a-1`, `00000027-b-1`,
`00000028-a-1`, `00000029-a-1`, `00000029-a-2`, `00000029-b-1`,
`00000030-a-1`, `00000031-a-1`, `00000031-b-1`, `00000031-c-1`,
`00000031-d-1`, `00000032-a-1`, `00000033-a-1`, `00000034-a-1`,
`00000034-b-1` (saturation mutagenesis MPRAs of regulatory elements).

Their note reads:

> We are making our data available prior to publication in line with Fort
> Lauderdale principle, allowing others to use the data but allowing the data
> producers to make the first presentations and to publish the first paper with
> global analyses of the data. In addition, we also reserve the right to
> publish the first analysis of the differences seen in the TERT knock-down
> experiments and alternative cell-type experiments. Studies that do not
> overlap with these intentions may be submitted for publication at any time,
> but must appropriately cite the data source. After publication of the data,
> the first publication of the data producers should be cited for any use of
> these data.

The data were published in 2019, so the citation condition applies. **Cite:**
Kircher M, et al. Saturation mutagenesis of twenty disease-associated
regulatory elements at single base-pair resolution. *Nat Commun.* 2019.
PMID 31395865. doi:[10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w).
Data DOI listed by MaveDB: 10.17605/OSF.IO/75B2M.

### Citation request: RHO, kept

Score set `urn:mavedb:00001275-a-1` (RHO surface trafficking DMS scores,
Method 1). Its note reads:

> Users are requested to cite the associated publication and dataset DOI when
> using these data.

**Cite:** Manian KV, et al. A comprehensive map of missense trafficking
variants in rhodopsin and their response to pharmacologic correction. bioRxiv
preprint, 2025. doi:[10.1101/2025.02.27.640335](https://doi.org/10.1101/2025.02.27.640335).
Dataset: doi:[10.5281/zenodo.19558926](https://doi.org/10.5281/zenodo.19558926).

### Descriptive notes only, kept

`urn:mavedb:00000009-a-2` (SUL1 promoter; PMID 26936925) and
`urn:mavedb:00001257-0-1` (ADSL; preprint doi 10.1101/2025.08.15.670494) have
notes that describe the study and add no condition.

### Excluded score sets

See below.

## Excluded from this release

### Third-party predictor scores (14 columns)

These columns are removed from every published table (`tier_a.parquet`,
`tier_a_genome.parquet`, `tier_a_exome.parquet`, and their rows in
`field_map.tsv`):

| Score | Genome column | Exome column |
|---|---|---|
| CADD (v1.6) | `cadd_phred` | `cadd_phred_exome` |
| REVEL (v1.3) | `revel_max` | `revel_max_exome` |
| SpliceAI (v1.3) | `spliceai_ds_max` | `spliceai_ds_max_exome` |
| Pangolin | `pangolin_largest_ds` | `pangolin_largest_ds_exome` |
| phyloP | `phylop` | `phylop_exome` |
| SIFT | `sift_max` | `sift_max_exome` |
| PolyPhen-2 | `polyphen_max` | `polyphen_max_exome` |

Why:

- CADD, from <https://cadd.gs.washington.edu/>: "CADD scores are freely
  available for all non-commercial applications. If you are planning on using
  them in a commercial application, you can obtain a license through the UW
  CoMotion Express Licensing System."
- REVEL, from <https://sites.google.com/site/revelgenomics/downloads>: "REVEL
  scores are freely available for non-commercial use. For other uses, please
  contact Weiva Sieh."
- SpliceAI: gnomAD's terms (quoted above) give SpliceAI annotations under CC BY
  NC 4.0 for academic and non-commercial use. Illumina's README
  (<https://github.com/Illumina/SpliceAI>) says of the precomputed scores:
  "These annotations are free for academic and not-for-profit use; other use
  requires a commercial license from Illumina, Inc."
- Pangolin, phyloP, SIFT and PolyPhen-2: the audit found no stated terms for
  these values as distributed in the gnomAD VCFs, so they are left out as
  unclear rather than assumed free.

None of those pages grants redistribution in an unrestricted public release.
The README explains how to regenerate the columns locally from the gnomAD VCFs
with `scripts/annotate_gnomad.py`. If you do, the providers' terms above apply
to you.

### Six MaveDB score sets

Removed from the published MaveDB tables (experiment inventory, scores, mapped
variants, gene coverage). The evidence ledger, the gene marts, the crosswalk
and the Tier A tables never contained rows from these sets (0 ledger rows
reference them), so they needed no change:

- `urn:mavedb:00000044-a-1`, `00000044-a-2`, `00000044-b-1`, `00000044-b-2`:
  SARS-CoV-2 receptor-binding-domain scans (viral, not human; 406,242 score
  rows, no mapped rows). Their note is a BSD 3-Clause licence text beginning
  "Copyright (c) 2020, Tyler N. Starr, Allison J. Greaney, and Jesse D. Bloom",
  including a non-endorsement clause. They are left out rather than shipped
  with a software licence attached to data.
- `urn:mavedb:00001255-a-1` and `00001255-a-2`: AIRE deep mutational scan
  (21,674 score rows, no variant keys). Their notes read "Data must be
  pre-printed prior to usage." and "Score data must be pre-printed prior to
  usage." The dump lists no preprint for them, so they are held back.

### Files never published

- `raw/`: every source download (about 768 GB). Re-fetch it with
  `scripts/fetch_sources.py`. It includes the GTEx sample-level TPM file
  (17,382 donor-sample columns) and the MaveDB `main.json` with researcher
  names and ORCID identifiers.
- `db/variants.sqlite`: the build's SQLite catalog (about 1 GB). Rebuild it
  with `scripts/load_core.py` and `scripts/build_manifest.py`.
- `logs/*.log`: local download and install logs. `logs/fetch.status` is kept.

No individual genotypes exist anywhere in the project. The gnomAD files are
sites-only VCFs with no sample columns.

### Reports and QC notes describe the full build

The reports in `warehouse/reports/` and the QC notes in `warehouse/qc/` were
written before these exclusions and are published unchanged. Their numbers
include the excluded predictor columns and MaveDB score sets. Some reports
quote aggregate counts and thresholds based on CADD, REVEL and SpliceAI, and
two of them quote individual scores as worked examples: 26 variants in all
(25 in `predictor_disagreement.md`, 1 in `conflicts.md`). Those values remain
under their providers' terms (non-commercial for CADD, REVEL and SpliceAI).
No table of predictor scores is published as data.
