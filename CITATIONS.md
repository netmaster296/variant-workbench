# Citations

## Citing this project

netmaster296. *variant-workbench: an open-source GRCh38 warehouse of genes,
studied variants and their evidence from public sources.* 2026.
<https://github.com/netmaster296/variant-workbench>. Code: GPL-3.0. Project
data contributions: CC0 1.0.

Machine-readable metadata is in `CITATION.cff`. If you use the data, also cite
the sources below. Each source's terms are in `NOTICE.md`.

"Confirmed" means the citation was read on the provider's own page during the
publication audit (November 2026). "Not confirmed" means the official page
could not be read, so the citation comes from another source; check it before
you use it.

## Data sources

### gnomAD (confirmed)

gnomAD asks for the flagship paper and a link to the browser
(<https://gnomad.broadinstitute.org/>), and gives this version format:

> gnomAD (v4.1, 807,162 individuals, released on 04/19/2024) (Chen et al., 2023).

- Chen S, et al. The gnomAD flagship paper.
  doi:[10.1038/s41586-023-06045-0](https://doi.org/10.1038/s41586-023-06045-0).
  (The gnomAD citation page links the DOI; it does not print the title.)

Statement for tools: "This tool includes data from the gnomAD v4.1 release."

### ClinVar (not confirmed)

The ClinVar documentation pages returned a bot check during the audit, so
NCBI's preferred citation was not seen. Commonly cited ClinVar papers:

- Landrum MJ, et al. *Nucleic Acids Res.* 2014.
  doi:[10.1093/nar/gkt1113](https://doi.org/10.1093/nar/gkt1113).
- Landrum MJ, et al. *Nucleic Acids Res.* 2016.
  doi:[10.1093/nar/gkv1222](https://doi.org/10.1093/nar/gkv1222).

Check <https://www.ncbi.nlm.nih.gov/clinvar/> for the current preferred
citation. Data: ClinVar weekly `variant_summary.txt.gz`, file dated 2026-11-01.

### GENCODE (confirmed)

- Mudge JM, Carbonell-Sala S, Diekhans M, et al. GENCODE 2025: reference gene
  annotation for human and mouse. *Nucleic Acids Res.* 2025;53(D1):D966-D975.
  PMID 39565199.
  doi:[10.1093/nar/gkae1078](https://doi.org/10.1093/nar/gkae1078).

Data: GENCODE release 49 (Ensembl 115).

### NHGRI-EBI GWAS Catalog (confirmed)

- Cerezo M, Sollis E, Ji Y, et al. *Nucleic Acids Res.* 2025;53(D1):D998-D1005.
  doi:[10.1093/nar/gkae1070](https://doi.org/10.1093/nar/gkae1070).
  (The GWAS Catalog citation page lists the authors, journal and DOI.)

Data: associations download from `releases/latest`, retrieved 2026-11-03; the
TSV inside is dated 2026-10-29.

### GTEx (acknowledgement confirmed; paper not confirmed)

Acknowledgement (wording from the GTEx FAQ; see `NOTICE.md`): The
Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of
the Office of the Director of the National Institutes of Health, and by NCI,
NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used here were obtained from the
GTEx v8 open-access files in the GTEx public Google Cloud bucket on 11/03/2026.

The GTEx portal's how-to-cite page is rendered by JavaScript and could not be
read. The v8 paper usually cited (seen only on a secondary page):

- GTEx Consortium. *Science.* 2020.
  doi:[10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776).

### MaveDB (confirmed)

From <https://mavedb.org/docs/mavedb/citation.html>:

- Rubin AF, Stone J, Bianchi AH, et al. MaveDB 2024: a curated community
  database with over seven million variant effects from multiplexed functional
  assays. *Genome Biol.* 2025;26(1):13. PMID 39838450.
  doi:[10.1186/s13059-025-03476-y](https://doi.org/10.1186/s13059-025-03476-y).
- Previous citation: Esposito D, Weile J, Shendure J, et al. MaveDB: an
  open-source platform to distribute and interpret data from multiplexed
  assays of variant effect. *Genome Biol.* 2019;20(1):223. PMID 31679514.
  doi:[10.1186/s13059-019-1845-6](https://doi.org/10.1186/s13059-019-1845-6).

Dataset: MaveDB bulk dump `mavedb-dump.2026062418131` (as of 2026-06-24),
Zenodo, doi:[10.5281/zenodo.20840937](https://doi.org/10.5281/zenodo.20840937).

MaveDB also asks you to cite the original publication of each dataset you use.
Each score set's publication is on its MaveDB page and in the dump's
`main.json`. Two groups ask for it explicitly (publication details read from
`main.json`):

- The 28 Fort Lauderdale score sets (`urn:mavedb:00000014` to `00000034`, listed
  in `NOTICE.md`): Kircher M, et al. Saturation mutagenesis of twenty
  disease-associated regulatory elements at single base-pair resolution.
  *Nat Commun.* 2019. PMID 31395865.
  doi:[10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w).
- RHO score set `urn:mavedb:00001275-a-1`: Manian KV, et al. A comprehensive map
  of missense trafficking variants in rhodopsin and their response to
  pharmacologic correction. bioRxiv preprint, 2025.
  doi:[10.1101/2025.02.27.640335](https://doi.org/10.1101/2025.02.27.640335).
  Dataset: doi:[10.5281/zenodo.19558926](https://doi.org/10.5281/zenodo.19558926).

## Predictor scores (only if you regenerate them locally)

The published tables do not contain these scores. If you regenerate them from
gnomAD (see `README.md`), cite them and follow their terms:

- CADD v1.6 (confirmed on <https://cadd.gs.washington.edu/>): Rentzsch P,
  Schubach M, Shendure J, Kircher M. CADD-Splice: improving genome-wide variant
  effect prediction using deep learning-derived splice scores. *Genome Med.*
  2021. PMID 33618777.
  doi:[10.1186/s13073-021-00835-9](https://doi.org/10.1186/s13073-021-00835-9).
- SpliceAI (confirmed during the audit): Jaganathan et al. *Cell.* 2019.
  doi:[10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015).
- Pangolin (confirmed during the audit): Zeng and Li. *Genome Biol.* 2022;23:103.
- REVEL (title confirmed during the audit; authors, journal and year not confirmed): Ioannidis NM, et al. REVEL: an ensemble
  method for predicting the pathogenicity of rare missense variants. *Am J Hum
  Genet.* 2016. (Authors, journal and year are not confirmed; the audit recorded
  only the title.)

The audit did not record citations for phyloP, SIFT or PolyPhen-2.

## Software

The build depends on these free and open-source projects. Thank you to their
developers.

| Project | Used for | Licence |
|---|---|---|
| Python 3.13 | everything | PSF License |
| DuckDB | SQL over Parquet | MIT |
| Apache Arrow (pyarrow) | Parquet reading and writing | Apache-2.0 |
| NumPy | `load_gtex.py` | BSD-3-Clause |
| Biopython | BGZF reads in `vcf_query.py` | Biopython License Agreement / BSD-3-Clause |
| Pillow | the cognitive-loci figure | MIT-CMU |
