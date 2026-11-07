# Wave 2

Free space on D: at the end of this integration: 669,961,555,968 bytes (`Get-PSDrive`). At the start, before the mart and Tier A files were written, it was 671,699,435,520 bytes. The difference is the Parquet and TSV written here. No gnomAD VCF was copied, and exome allele counts were not added to genome allele counts. bcftools was not installed.

## Exome chromosomes

Every exome sites VCF from chr1–chr22, chrX, and chrY is on disk with a `.tbi`. None of the bgzip files are empty. chrY is 113,561,604 bytes and its index is 5,364 bytes. Genome sites for the same 24 chromosomes are also present; there is no chrM sites VCF.

| chromosome | vcf bytes | tbi bytes |
| --- | ---: | ---: |
| chr1 | 18,789,775,480 | 130,696 |
| chr2 | 14,282,176,877 | 114,639 |
| chr3 | 11,586,472,818 | 92,706 |
| chr4 | 7,712,255,367 | 78,479 |
| chr5 | 8,494,344,938 | 78,966 |
| chr6 | 9,132,083,331 | 77,127 |
| chr7 | 9,675,275,184 | 77,730 |
| chr8 | 7,275,385,426 | 66,417 |
| chr9 | 8,041,521,990 | 59,035 |
| chr10 | 7,858,233,309 | 63,934 |
| chr11 | 11,641,720,279 | 72,428 |
| chr12 | 10,632,124,827 | 73,346 |
| chr13 | 3,540,233,028 | 36,968 |
| chr14 | 6,669,068,987 | 47,202 |
| chr15 | 7,374,290,946 | 51,155 |
| chr16 | 9,774,275,728 | 47,912 |
| chr17 | 12,017,042,594 | 55,983 |
| chr18 | 3,390,375,115 | 33,492 |
| chr19 | 12,620,525,657 | 45,979 |
| chr20 | 4,760,834,034 | 33,684 |
| chr21 | 2,250,319,809 | 14,118 |
| chr22 | 5,060,347,554 | 22,536 |
| chrX | 5,746,768,174 | 57,299 |
| chrY | 113,561,604 | 5,364 |

## Where the files disagree with the plan

The loaders' columns were used. These are the disagreements that changed a name or a count.

- `gene_constraint.parquet` stores LOEUF as `loeuf` (the same value as `lof.oe_ci.upper`) and the decile as `loeuf_decile`. There is no column named `LOEUF`. The mart copies `loeuf`, `loeuf_decile`, `lof.pLI`, `mis.z_score`, `syn.z_score`, and `constraint_flags`.
- `gene_expression.parquet` already has `max_tissue`, `max_median_tpm`, `n_tissues_tpm_ge_1`, and `tau`. The mart also keeps `n_tissues_tpm_ge_10`, which is in that file. Tau is the GTEx lane's formula, not recomputed.
- `gene_coverage.parquet` counts score sets (`n_score_sets`), not experiment URNs. `mavedb_experiments_by_symbol` is that score-set count. Exact symbol match is a second, weaker count beside variants whose GRCh38 position falls inside the gene.
- `loci.parquet` `mapped_genes` is a struct. Tokens are already split. `FBLN1 - LINC01589` is one token. A symbol must equal a token. It does not match a substring, and the overlap count does not use the token.
- ClinVar `GeneSymbol` separates multiple symbols with `;`. Each token is matched exactly. A comma phrase is not split. An allele with two matching tokens is counted on each gene.
- The gnomAD `vep` description is `Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|UNIPROT_ISOFORM|SOURCE|DOMAINS|miRNA|HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info`. There is no `Existing_variation` subfield. SIFT and PolyPhen are the INFO fields `sift_max` and `polyphen_max`, not VEP subfields. The chosen block prefers a non-empty `MANE_SELECT`, then `CANONICAL=YES`, then the highest IMPACT. `vep_on_mane` is true only when that block's `MANE_SELECT` is non-empty. The raw VEP string is not stored. The full map is `warehouse/variants/field_map.tsv`.
- The exome chrY header has no `AF_ami`. `AF_ami_exome` is still a column and is null. The five ancestry fields the audit needs (`AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, `AF_amr`) are in both headers. Sex-split `AF_*_XX` / `AF_*_XY`, `faf*`, and histogram fields were not kept.
- eGene and sGene effect columns are `slope`, not `beta`. The ledger uses `pval_nominal`. Only rows with `qval_le_0_05` are evidence. The pair tables were not read.

## Gene mart

`warehouse/marts/gene_mart.parquet` has 78,691 rows, one per GENCODE gene. `gene_mart_protein_coding.parquet` has 20,097 rows, protein-coding only. Coordinates and biotype come from `intervals/genes.parquet`. They match the crosswalk on chrom, start, end, symbol, and biotype for every gene (0 mismatches). `mane_enst` and `mane_enst_constraint` are both kept. `mane_enst` is set on 19,284 genes and `mane_enst_constraint` on 17,477.

Constraint is a left join on bare ENSG. It was not rescued by symbol. Protein-coding genes with no constraint row: 1,530, the same miss the crosswalk reported. A further 469 protein-coding genes match a constraint row whose LOEUF is null, so 1,999 protein-coding mart rows have a null `loeuf`. Across all biotypes, 60,078 genes have no constraint row and 60,570 have a null `loeuf`.

Expression joins on bare ENSG. GTEx v8 is older than GENCODE v49. 24,065 GENCODE genes have no expression row. 1,530 GTEx genes with a bare ENSG are not in GENCODE. 44 GTEx rows have a null `ensg` (the `_PAR_Y` ids). `eqtl_tissues` is null for 40,148 genes and `sqtl_tissues` is null for 55,564. A null means the bare ENSG is absent from that table. A zero means the gene was tested and no tissue has `qval_le_0_05`.

ClinVar counts come from `allele_dedup.parquet` by symbol token. 4,623,299 alleles matched at least one GENCODE symbol. 15,954 did not. 17,285 alleles have more than one `;` token. 2,096 genes share a symbol with another gene, and 119 of those also have a ClinVar count. Those genes each received the full count for the symbol. The pathogenic buckets follow `germline_pathogenic`, so a somatic or oncogenicity row is not germline pathogenic. `Pathogenic/Likely pathogenic` counts as germline pathogenic, not also as likely pathogenic. `Benign/Likely benign` counts as benign. Drug response is counted in addition to any other class. Benign, VUS, and conflicting are not limited to germline origin. Summed across genes, germline pathogenic is 222,278 and likely pathogenic is 116,316. That sum is a few alleles higher than the dedup flag because a two-symbol allele is counted once per symbol. Conflicting 184,177. VUS 2,399,580. Drug response 1,978. 21,452 genes have at least one matched allele.

GWAS keeps both rules. Overlap (same chromosome, and either the locus interval or the lead position meets the gene span): 75,174 genes. Mapped-gene token equals the symbol, with no chromosome requirement and no substring match: 23,724 genes. The locus `start`/`end` are the member positions, not the ±1 Mb search radius.

MaveDB symbol match hits 135 genes. That is the weak key. Distinct GRCh38 variant keys inside the gene span hit 124 genes. No MaveDB symbol collision landed on two GENCODE genes (`mavedb_symbol_collision_genes` is 0).

Truth-set genes on the mart: A1BG `ENST00000263100` LOEUF 1.34 decile 7; XPO6 0.04 decile 0; MED13L 0.041 decile 0; EIF4G1 0.048 decile 0; CABLES2, CFAP69, and CHI3L1 LOEUF 1.0 decile 5. Each symbol is one row, and `mane_enst_constraint` matches `mane_enst` for all seven.

## Tier A

`warehouse/variants/tier_a_keys.tsv` has 5,928,738 keys, sorted by chromosome and position. Sources, as distinct keys carrying that label: ClinVar 4,618,645, eQTL leads 830,649, sQTL leads 459,332, MaveDB GRCh38 212,102, GWAS 13,698, and one `truth_set` key. The significant-pair tables were not used.

GWAS has no ref or alt. A key was attached only when the association rsID matched a ClinVar rsID that already had a variant key. 12,759 GWAS rsIDs with a position matched. 364 of those positions do not match the ClinVar key's chromosome and position; the ClinVar key was still the one tagged, and it was not moved. 416,792 rows (416,791 rsIDs) with a position and no ClinVar key are in `tier_a_unresolved.tsv`. 7,624 GWAS rows have an rsID and no position; they are not in that file. The genome-wide gnomAD ID scan was deferred: bcftools is not available on Windows, and this pass does not stream every sites VCF a second time just to recover alleles for rsIDs that have no ref or alt. The annotation stream matches positions of keys that already have alleles.

`chr22:10510183:A:G` was not in ClinVar, MaveDB, or the lead QTL tables. It was added with source `truth_set` so the frequency case is in the annotated table. `chr11:126275389:C:T` is ClinVar. `chr7:117559590:ATCT:A` is ClinVar and GWAS.

Annotation streamed each chromosome once per dataset through `scripts/vcf_stream.py` (zlib BGZF, with the TBI linear index used only to skip a gap of at least 256 kb). chrY was read first. Its header lists `AF_nfe` and `AF_afr` as separate INFO fields, and the first records that carry both parse them as separate keys. There is no chrM VCF. The 3,136 chrM keys stay in the table with `in_genomes` false and `in_exomes` false.

`tier_a_genome.parquet` and `tier_a_exome.parquet` each have 5,928,738 rows. `tier_a.parquet` joins them on `variant_key` and also has 5,928,738 distinct keys. `in_genomes` is true for 3,220,597 keys. `in_exomes` is true for 3,283,322. Both are true for 2,061,870. A missing site stays in the table. Exome columns are suffixed `_exome`. `AC` and `AC_exome` are not summed.

Checked against the sites VCFs: `chr22:10510183:A:G` has genome `AF_nfe` 0 and `AF_afr` 0.00229358, genome `AC` 2, and a null exome side. `chr11:126275389:C:T` is C/T with genome `AC` 5 and exome `AC` 53. `chr7:117559590:ATCT:A` is ATCT/A with genome `AC` 1200 and exome `AC` 18037.

## Ledger

`warehouse/ledger/evidence.parquet` has 10,866,762 rows and 10,866,762 distinct `evidence_id` values. All six classes are present.

| class | rows |
| --- | ---: |
| clinvar | 4,639,253 |
| constraint | 18,628 |
| eqtl | 475,829 |
| gwas | 312,875 |
| mavedb | 5,277,857 |
| sqtl | 142,320 |

ClinVar is one row per deduplicated allele. `object_id` is the `PhenotypeIDS` string, not one row per cross-reference. 4,618,645 use a variant key. 5,981 with no key and a unique symbol use the bare ENSG. 14,627 with no key and no unique symbol use the dedup key. Strength is review stars.

GWAS is one row per clumped locus, plus one extra row when several ClinVar alleles share the lead rsID at the lead position (1,635 loci). 20,345 loci got a variant key at the same chromosome and position. 290,815 did not. Their `subject_id` is `chrom:pos`, which is not a ref/alt key, because the Catalog has no alleles. 269 leads matched a ClinVar rsID only at a different position; those were not attached. Sample-size prose is in `n_note` for 311,097 loci. Ancestry was not parsed out of that prose, so `ancestry_note` is null. Strength is the lead p-value.

eQTL and sQTL are lead rows with `qval_le_0_05` only (475,829 and 142,320). `object_id` is `ensg|tissue`. The sQTL natural key also includes `phenotype_id`, so intron clusters do not collapse. Pair tables stayed out.

MaveDB rows are mapped scores, not the unmapped score file. The signed score is stored. It is not an absolute value and it is not comparable across score sets. 145,938 mapped rows had no score. 927,103 had a score but neither a variant key nor a unique GENCODE symbol, so they stay in `mapped_variants_scored.parquet` and are not in the ledger. 3,474,823 ledger rows are variant subjects and 1,803,034 are gene subjects.

Constraint is one row per `gene_constraint` bare ENSG (18,628). The 88 rows with a null bare ENSG were not symbol-joined. `object_id` is `loeuf`. Strength is the LOEUF value and is null when the source LOEUF is null.

## Dark genome

`warehouse/reports/dark_genome.md`. Contig lengths were read from the chr21 genome VCF `##contig` header by `scripts/vcf_stream.py`, not from a hard-coded table. The denominator is 3,088,269,832 bases (chr1–chr22, chrX, chrY). chrM is not in that header. 37 gene spans, 13 CDS intervals, and 3,136 Tier A keys fall outside it and are excluded.

| set | bases | fraction |
| --- | ---: | ---: |
| protein-coding CDS | 36,422,072 | 0.011794 |
| any gene span | 2,064,230,736 | 0.668410 |
| within 1 kb of a Tier A key | 1,378,873,510 | 0.446487 |
| clumped GWAS member span | 2,736,332,511 | 0.886041 |
| in any of those sets | 2,820,076,717 | 0.913157 |
| in none of them | 268,193,115 | 0.086843 |

The four category counts are unions inside themselves and are not exclusive of each other. The remainder is exclusive. CDS is `feature = CDS` on protein-coding genes. `stop_codon` rows are not included. The 1 kb window is inclusive `[pos - 1000, pos + 1000]`. The GWAS interval is `loci.start`–`loci.end`, the member span, not the ±1 Mb search radius. That span covers most of the primary assembly because the catalog is dense: 202,104 loci are a single base, and the rest tile a large fraction of the genome. This is not a claim that those bases are explained.

## Truth set

`scripts/truth_set.py` was rerun after the mart, the Tier A join, and the ledger existed. Result: 18 pass, 0 fail, 0 waiting. The nine cases that had already passed were not edited. The gene-mart and frequency cases that had been waiting now pass.

`af_pop_diff` needed one harness change. The joined table has `AF_nfe` and `AF_afr` from the genome VCF and `AF_nfe_exome` and `AF_afr_exome` from the exome VCF. This site is not in the exome VCF, so the exome columns are null. The old check treated every column whose name contains `AF_nfe` or `AF_afr` as required. A null exome column is no longer a failure when a genome `AF_nfe` and a genome `AF_afr` are both stored and no stored `AF_nfe` equals a stored `AF_afr`. The values on this row are genome 0 and 0.00229358.
