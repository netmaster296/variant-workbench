# Truth-set report

| status | cases |
| --- | ---: |
| pass | 18 |
| fail | 0 |
| waiting | 0 |

18 cases from `D:\variant-workbench\warehouse\qc\truth_set.tsv`.
A case is fail when a file that exists does not match.
A case is waiting when every check against an existing file matched and a table that case still needs is not built.
Waiting is not a failure. No Parquet table was edited.

Failing cases: none.

Waiting cases: none.

## Not built

- `D:\variant-workbench\warehouse\marts\gene_mart.parquet` — present
- `D:\variant-workbench\warehouse\variants\tier_a.parquet` — present
- `D:\variant-workbench\warehouse\variants\tier_a_genome.parquet` — present
- `D:\variant-workbench\warehouse\variants\tier_a_exome.parquet` — present
- `D:\variant-workbench\warehouse\ledger\evidence.parquet` — present

`evidence.parquet` is not a case result. The somatic case is read from `allele_dedup.parquet`.

## Cases

| case | status |
| --- | --- |
| `coding_snp` | pass |
| `cftr_indel` | pass |
| `gwas_intergenic` | pass |
| `gwas_empty_pos` | pass |
| `mismatched_ref` | pass |
| `af_pop_diff` | pass |
| `mane_select` | pass |
| `loeuf_low_1` | pass |
| `loeuf_low_2` | pass |
| `loeuf_low_3` | pass |
| `loeuf_near_1_1` | pass |
| `loeuf_near_1_2` | pass |
| `loeuf_near_1_3` | pass |
| `gtex_egene` | pass |
| `mavedb_mapped` | pass |
| `mavedb_unmapped` | pass |
| `chrom_crosswalk` | pass |
| `somatic_or_oncogenicity` | pass |

### coding_snp: pass

Requirement: ClinVar, GWAS if present, and gnomAD agree on variant_key

Truth-set expected: chr11:126275389:C:T

| check | status | observed | file |
| --- | --- | --- | --- |
| ClinVar alleles | pass | variant_key chr11:126275389:C:T; GeneSymbol FOXRED1; rsid rs267606829; ref C; alt T; review_stars 2; ReviewStatus criteria provided, multiple submitters, no conflicts; ClinicalSignificance Pathogenic | `D:\variant-workbench\warehouse\clinvar\alleles.parquet` |
| GWAS rsid | pass | rsid rs267606829 rows 0; SNP_ID_CURRENT 267606829 rows 0; SNP_ID_CURRENT rs267606829 rows 0 | `D:\variant-workbench\warehouse\gwas\associations.parquet` |
| gnomAD genomes | pass | ref C alts T id rs267606829 contains T true | `D:\variant-workbench\raw\gnomad\genomes\gnomad.genomes.v4.1.sites.chr11.vcf.bgz` |
| gnomAD exomes | pass | ref C alts T id rs267606829 contains T true | `D:\variant-workbench\raw\gnomad\exomes\gnomad.exomes.v4.1.sites.chr11.vcf.bgz` |
| Tier A | pass | rows 1; ref C; alt T | `D:\variant-workbench\warehouse\variants\tier_a.parquet` |

### cftr_indel: pass

Requirement: Indel alleles match gnomAD when gnomAD has that allele, or the miss is recorded

Truth-set expected: stars 4; significance Pathogenic

| check | status | observed | file |
| --- | --- | --- | --- |
| ClinVar stars | pass | variant_key chr7:117559590:ATCT:A; GeneSymbol CFTR; ref ATCT; alt A; review_stars 4; ReviewStatus practice guideline; ClinicalSignificance Pathogenic | `D:\variant-workbench\warehouse\clinvar\alleles.parquet` |
| gnomAD genomes | pass | ref ATCT alts A id rs113993960 contains A true | `D:\variant-workbench\raw\gnomad\genomes\gnomad.genomes.v4.1.sites.chr7.vcf.bgz` |
| gnomAD exomes | pass | ref ATCT alts A id rs113993960 contains A true | `D:\variant-workbench\raw\gnomad\exomes\gnomad.exomes.v4.1.sites.chr7.vcf.bgz` |

### gwas_intergenic: pass

Requirement: Lands outside coding exons and still has a trait URI

Truth-set expected: http://www.ebi.ac.uk/efo/EFO_0004348

| check | status | observed | file |
| --- | --- | --- | --- |
| associations | pass | source line 8; chrom chr22; pos 45604058; rsid rs4823308; MAPPED_TRAIT_URI http://www.ebi.ac.uk/efo/EFO_0004348; CONTEXT intergenic_variant; INTERGENIC 1; locus_key chr22:45604058 | `D:\variant-workbench\warehouse\gwas\associations.parquet` |
| CDS overlap | pass | CDS intervals covering chr22:45604058: 0 | `D:\variant-workbench\warehouse\intervals\cds.parquet` |

### gwas_empty_pos: pass

Requirement: Stays in the association table and stays out of the clump

Truth-set expected: in associations; absent from loci

| check | status | observed | file |
| --- | --- | --- | --- |
| association row | pass | source line 5; STUDY ACCESSION GCST004604; PUBMEDID 27863252; SNPS chr22:37462936; CHR_POS ''; CHR_ID ''; chrom null; pos null; rsid null; locus_key null; MAPPED_TRAIT_URI http://www.ebi.ac.uk/efo/EFO_0004348 | `D:\variant-workbench\warehouse\gwas\associations.parquet` |
| absent from loci | pass | no study column on loci; this row pos null; locus_key null; loci with null lead_pos or chrom 0; row rsid is null, so it cannot be a locus lead_rsid | `D:\variant-workbench\warehouse\gwas\loci.parquet` |

### mismatched_ref: pass

Requirement: The key uses the VCF ref and alt, not the ClinVar display allele

Truth-set expected: chr15:84799209:G:A

| check | status | observed | file |
| --- | --- | --- | --- |
| ClinVar key | pass | variant_key chr15:84799209:G:A; ref G; alt A; ReferenceAlleleVCF G; AlternateAlleleVCF A; ReferenceAllele na; AlternateAllele na; Name NM_014630.3(ZNF592):c.3136G>A (p.Gly1046Arg) | `D:\variant-workbench\warehouse\clinvar\alleles.parquet` |

### af_pop_diff: pass

Requirement: AF_nfe and AF_afr are both stored and neither overwrites the other

Truth-set expected: AF_nfe 0.00000; AF_afr 0.00229358

| check | status | observed | file |
| --- | --- | --- | --- |
| AF_nfe and AF_afr | pass | ref A alts G AF_nfe 0.00000 AF_afr 0.00229358 raw AF_nfe 0.00000 raw AF_afr 0.00229358 | `D:\variant-workbench\raw\gnomad\genomes\gnomad.genomes.v4.1.sites.chr22.vcf.bgz` |
| Tier A frequencies | pass | AF_nfe 0; AF_nfe_exome null; AF_afr 0.00229358; AF_afr_exome null | `D:\variant-workbench\warehouse\variants\tier_a.parquet` |

### mane_select: pass

Requirement: The gene-constraint row uses the MANE Select transcript

Truth-set expected: ENST00000263100

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene A1BG; gene_id ENSG00000121410; transcript ENST00000263100; bare_ensg ENSG00000121410; LOEUF 1.34; decile 7; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol A1BG; loeuf_decile 7; loeuf 1.34; mane_enst ENST00000263100 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_low_1: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 0.04; decile 0; transcript ENST00000304658

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene XPO6; gene_id ENSG00000169180; transcript ENST00000304658; bare_ensg ENSG00000169180; LOEUF 0.04; decile 0; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol XPO6; loeuf_decile 0; loeuf 0.04; mane_enst ENST00000304658 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_low_2: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 0.041; decile 0; transcript ENST00000281928

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene MED13L; gene_id ENSG00000123066; transcript ENST00000281928; bare_ensg ENSG00000123066; LOEUF 0.041; decile 0; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol MED13L; loeuf_decile 0; loeuf 0.041; mane_enst ENST00000281928 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_low_3: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 0.048; decile 0; transcript ENST00000346169

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene EIF4G1; gene_id ENSG00000114867; transcript ENST00000346169; bare_ensg ENSG00000114867; LOEUF 0.048; decile 0; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol EIF4G1; loeuf_decile 0; loeuf 0.048; mane_enst ENST00000346169 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_near_1_1: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 1.0; decile 5; transcript ENST00000279101

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene CABLES2; gene_id ENSG00000149679; transcript ENST00000279101; bare_ensg ENSG00000149679; LOEUF 1; decile 5; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol CABLES2; loeuf_decile 5; loeuf 1; mane_enst ENST00000279101 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_near_1_2: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 1.0; decile 5; transcript ENST00000389297

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene CFAP69; gene_id ENSG00000105792; transcript ENST00000389297; bare_ensg ENSG00000105792; LOEUF 1; decile 5; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol CFAP69; loeuf_decile 5; loeuf 1; mane_enst ENST00000389297 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### loeuf_near_1_3: pass

Requirement: Gene mart carries this LOEUF decile

Truth-set expected: LOEUF 1.0; decile 5; transcript ENST00000255409

| check | status | observed | file |
| --- | --- | --- | --- |
| gene_constraint | pass | gene CHI3L1; gene_id ENSG00000133048; transcript ENST00000255409; bare_ensg ENSG00000133048; LOEUF 1; decile 5; mane_select true | `D:\variant-workbench\warehouse\constraint\gene_constraint.parquet` |
| gene mart | pass | symbol CHI3L1; loeuf_decile 5; loeuf 1; mane_enst ENST00000255409 | `D:\variant-workbench\warehouse\marts\gene_mart.parquet` |

### gtex_egene: pass

Requirement: chr_pos_ref_alt_b38 becomes the same variant key as the split columns

Truth-set expected: chr1:64764:C:T

| check | status | observed | file |
| --- | --- | --- | --- |
| eGene keys | pass | 1 rows; tissue Adipose_Subcutaneous; gene_id ENSG00000227232.5; chr chr1; variant_pos 64764; ref C; alt T; split_key chr1:64764:C:T; id_key chr1:64764:C:T; variant_key chr1:64764:C:T | `D:\variant-workbench\warehouse\gtex\egenes.parquet` |

### mavedb_mapped: pass

Requirement: A mapped GRCh38 variant gets a variant key

Truth-set expected: chr5:146469691:A:T

| check | status | observed | file |
| --- | --- | --- | --- |
| mapped score | pass | 1 rows; urn urn:mavedb:00000352-0-1#464; variant_key chr5:146469691:A:T; chrom chr5; pos 146469691; ref A; alt T; score -0.153873; score_text -0.1538729821505993; target_gene TCERG1 (Y438R mutant) | `D:\variant-workbench\warehouse\mavedb\mapped_variants_scored.parquet` |

### mavedb_unmapped: pass

Requirement: An unmapped score stays at gene level and does not get a coordinate

Truth-set expected: no variant_key

| check | status | observed | file |
| --- | --- | --- | --- |
| experiment | pass | score_set_urn urn:mavedb:00000004-a-1; rows 1; target_name ['E4B']; mapped_file_exists ['false']; mapping_state ['failed']; score_set_title ['Enrich scores for E4B'] | `D:\variant-workbench\warehouse\mavedb\experiments.parquet` |
| scores | pass | score rows 98298; non-null score 98297; non-null score_text 98298 | `D:\variant-workbench\warehouse\mavedb\scores.parquet` |
| coordinates | pass | mapped rows 0; non-null chrom 0; non-null pos 0; non-null variant_key 0 | `D:\variant-workbench\warehouse\mavedb\mapped_variants.parquet` |

### chrom_crosswalk: pass

Requirement: ClinVar chromosome 1 joins a GENCODE gene on chr1

Truth-set expected: chr1

| check | status | observed | file |
| --- | --- | --- | --- |
| ClinVar chromosome | pass | variant_key chr1:243344276:AG:A; Chromosome 1; chrom_raw 1; chrom chr1; GeneSymbol SDCCAG8; pos 243344276; ref AG; alt A | `D:\variant-workbench\warehouse\clinvar\alleles.parquet` |
| crosswalk gene | pass | ensg ENSG00000054282; ensg_version ENSG00000054282.17; symbol SDCCAG8; chrom chr1; start 243255349; end 243500095; biotype protein_coding; pos inside span true | `D:\variant-workbench\warehouse\crosswalk\genes.parquet` |

### somatic_or_oncogenicity: pass

Requirement: Does not count as germline pathogenic

Truth-set expected: exclude from germline pathogenic counts

| check | status | observed | file |
| --- | --- | --- | --- |
| allele_dedup | pass | variant_key chr17:31334927:C:T; GeneSymbol NF1; germline_pathogenic false; Oncogenicity Likely oncogenic; ClinicalSignificance Pathogenic; OriginSimple germline/somatic; germline_somatic germline/somatic | `D:\variant-workbench\warehouse\clinvar\allele_dedup.parquet` |
