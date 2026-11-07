-- ClinVar splice-region alleles, gnomAD SpliceAI and Pangolin, and GTEx sGenes.
-- Genome and exome scores stay in separate columns. They are not averaged.
-- sQTL pairs are not read. sGene status comes from sgenes.parquet only.
-- DuckDB. Statements end with a semicolon. The connection starts empty.

SET threads = 8;
SET memory_limit = '16GB';
SET preserve_insertion_order = false;

CREATE TEMP VIEW allele AS
SELECT
  d.variant_key,
  d."GeneID" AS gene_id,
  d."GeneSymbol" AS gene_symbol,
  d."Type" AS var_type,
  d."ClinicalSignificance" AS clinical_significance,
  d.germline_somatic,
  t.in_genomes,
  t.in_exomes,
  t.spliceai_ds_max AS spliceai_genome,
  t.spliceai_ds_max_exome AS spliceai_exome,
  t.pangolin_largest_ds AS pangolin_genome,
  t.pangolin_largest_ds_exome AS pangolin_exome,
  list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_region_variant') AS sr_g,
  list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_region_variant') AS sr_e,
  list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_donor_region_variant') AS dreg_g,
  list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_donor_region_variant') AS dreg_e,
  list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_donor_5th_base_variant') AS fifth_g,
  list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_donor_5th_base_variant') AS fifth_e,
  (
    list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_donor_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_acceptor_variant')
  ) AS canon_g,
  (
    list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_donor_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_acceptor_variant')
  ) AS canon_e,
  list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_polypyrimidine_tract_variant') AS poly_g,
  list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_polypyrimidine_tract_variant') AS poly_e,
  (
    list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_region_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_donor_region_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence, ''), '&'), 'splice_donor_5th_base_variant')
  ) AS region_g,
  (
    list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_region_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_donor_region_variant')
    OR list_contains(string_split(coalesce(t.vep_consequence_exome, ''), '&'), 'splice_donor_5th_base_variant')
  ) AS region_e
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') d
JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
  ON d.variant_key = t.variant_key
WHERE d.variant_key IS NOT NULL;

CREATE TEMP VIEW region AS
SELECT * FROM allele WHERE region_g OR region_e;

CREATE TEMP VIEW gene_status AS
SELECT
  ensg,
  count(*) FILTER (WHERE qval_le_0_05) AS n_sig_phenotypes,
  count(DISTINCT tissue) FILTER (WHERE qval_le_0_05) AS n_sig_tissues,
  min(qval) FILTER (WHERE qval_le_0_05) AS min_sig_qval
FROM read_parquet('D:/variant-workbench/warehouse/gtex/sgenes.parquet')
GROUP BY ensg;

CREATE TEMP VIEW sig_leads AS
SELECT DISTINCT variant_key, ensg
FROM read_parquet('D:/variant-workbench/warehouse/gtex/sgenes.parquet')
WHERE qval_le_0_05;

CREATE TEMP VIEW variant_overlap AS
SELECT
  r.variant_key,
  r.gene_id,
  r.gene_symbol,
  r.var_type,
  r.clinical_significance,
  r.germline_somatic,
  r.in_genomes,
  r.in_exomes,
  r.region_g,
  r.region_e,
  r.sr_g,
  r.dreg_g,
  r.fifth_g,
  r.sr_e,
  r.dreg_e,
  r.fifth_e,
  r.canon_g,
  r.canon_e,
  r.spliceai_genome,
  r.spliceai_exome,
  r.pangolin_genome,
  r.pangolin_exome,
  count(g.ensg) AS n_linked_ensg,
  count(*) FILTER (WHERE s.n_sig_phenotypes > 0) AS n_linked_sgene,
  count(*) FILTER (
    WHERE s.ensg IS NOT NULL AND coalesce(s.n_sig_phenotypes, 0) = 0
  ) AS n_linked_tested_not_sgene,
  count(*) FILTER (WHERE g.ensg IS NOT NULL AND s.ensg IS NULL) AS n_linked_not_in_sgene_table,
  count(*) FILTER (WHERE lead.ensg IS NOT NULL) AS n_linked_ensg_allele_is_sig_lead
FROM region r
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/crosswalk/genes.parquet') g
  ON g.ncbi_gene_id_primary = r.gene_id
LEFT JOIN gene_status s ON s.ensg = g.ensg
LEFT JOIN sig_leads lead
  ON lead.variant_key = r.variant_key AND lead.ensg = g.ensg
GROUP BY
  r.variant_key, r.gene_id, r.gene_symbol, r.var_type, r.clinical_significance,
  r.germline_somatic, r.in_genomes, r.in_exomes, r.region_g, r.region_e,
  r.sr_g, r.dreg_g, r.fifth_g, r.sr_e, r.dreg_e, r.fifth_e, r.canon_g, r.canon_e,
  r.spliceai_genome, r.spliceai_exome, r.pangolin_genome, r.pangolin_exome;

CREATE TEMP VIEW gene_overlap AS
SELECT
  g.ensg,
  any_value(g.symbol) AS symbol,
  any_value(g.biotype) AS biotype,
  count(DISTINCT r.variant_key) AS n_region_variants,
  count(DISTINCT r.variant_key) FILTER (
    WHERE r.spliceai_genome IS NOT NULL
      OR r.spliceai_exome IS NOT NULL
      OR r.pangolin_genome IS NOT NULL
      OR r.pangolin_exome IS NOT NULL
  ) AS n_scored_variants,
  s.n_sig_phenotypes,
  s.n_sig_tissues,
  s.min_sig_qval,
  coalesce(s.n_sig_phenotypes, 0) > 0 AS is_sgene,
  s.ensg IS NOT NULL AS in_sgene_table
FROM region r
JOIN read_parquet('D:/variant-workbench/warehouse/crosswalk/genes.parquet') g
  ON g.ncbi_gene_id_primary = r.gene_id
LEFT JOIN gene_status s ON s.ensg = g.ensg
GROUP BY g.ensg, s.n_sig_phenotypes, s.n_sig_tissues, s.min_sig_qval, s.ensg;

-- Join audit. A non-zero variant_key miss is a broken join.
SELECT
  'clinvar_key_missing_from_tier_a' AS check_name,
  count(*) AS n
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') d
ANTI JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
  ON d.variant_key = t.variant_key
WHERE d.variant_key IS NOT NULL
UNION ALL
SELECT
  'tier_a_clinvar_source_missing_from_dedup',
  count(*)
FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
ANTI JOIN read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') d
  ON t.variant_key = d.variant_key
WHERE list_contains(string_split(t.source, ','), 'clinvar');

SELECT
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')) AS dedup_rows,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') WHERE variant_key IS NULL) AS null_variant_key,
  (SELECT count(*) FROM allele) AS keyed_in_tier_a,
  (SELECT count(*) FROM allele WHERE in_genomes) AS in_genomes,
  (SELECT count(*) FROM allele WHERE in_exomes) AS in_exomes,
  (SELECT count(*) FROM allele WHERE NOT in_genomes AND NOT in_exomes) AS in_neither_vcf,
  (SELECT count(*) FROM allele WHERE in_genomes OR in_exomes) AS in_either_vcf,
  (SELECT count(*) FROM region) AS splice_region,
  (SELECT count(*) FROM allele WHERE (canon_g OR canon_e) AND NOT (region_g OR region_e)) AS canonical_only_excluded,
  (SELECT count(*) FROM allele WHERE (poly_g OR poly_e) AND NOT (region_g OR region_e)) AS polypyrimidine_only_excluded,
  (SELECT count(*) FROM region WHERE canon_g OR canon_e) AS region_also_canonical;

SELECT
  count(*) FILTER (WHERE region_g AND region_e) AS both_datasets,
  count(*) FILTER (WHERE region_g AND NOT in_exomes) AS genome_only_not_in_exome,
  count(*) FILTER (WHERE region_e AND NOT in_genomes) AS exome_only_not_in_genome,
  count(*) FILTER (WHERE region_g AND in_exomes AND NOT region_e) AS disagree_genome_only,
  count(*) FILTER (WHERE region_e AND in_genomes AND NOT region_g) AS disagree_exome_only,
  count(*) FILTER (WHERE sr_g AND NOT dreg_g AND NOT fifth_g) AS genome_splice_region_variant,
  count(*) FILTER (WHERE dreg_g AND NOT fifth_g) AS genome_donor_region,
  count(*) FILTER (WHERE fifth_g) AS genome_donor_5th_base,
  count(*) FILTER (WHERE fifth_g AND dreg_g) AS genome_fifth_and_donor_region,
  count(*) FILTER (WHERE sr_g AND fifth_g) AS genome_region_and_fifth,
  count(*) FILTER (WHERE sr_g AND dreg_g) AS genome_region_and_donor_region,
  count(*) FILTER (WHERE sr_e AND NOT dreg_e AND NOT fifth_e) AS exome_splice_region_variant,
  count(*) FILTER (WHERE dreg_e AND NOT fifth_e) AS exome_donor_region,
  count(*) FILTER (WHERE fifth_e) AS exome_donor_5th_base,
  count(*) FILTER (WHERE fifth_e AND dreg_e) AS exome_fifth_and_donor_region,
  count(*) FILTER (WHERE sr_e AND fifth_e) AS exome_region_and_fifth,
  count(*) FILTER (WHERE sr_e AND dreg_e) AS exome_region_and_donor_region
FROM variant_overlap;

SELECT
  count(*) AS n_region,
  count(*) FILTER (
    WHERE spliceai_genome IS NOT NULL
      OR spliceai_exome IS NOT NULL
      OR pangolin_genome IS NOT NULL
      OR pangolin_exome IS NOT NULL
  ) AS any_score,
  count(*) FILTER (WHERE in_genomes) AS in_genomes,
  count(*) FILTER (WHERE in_exomes) AS in_exomes,
  count(spliceai_genome) AS spliceai_genome,
  count(spliceai_exome) AS spliceai_exome,
  count(pangolin_genome) AS pangolin_genome,
  count(pangolin_exome) AS pangolin_exome,
  count(*) FILTER (WHERE in_genomes AND spliceai_genome IS NULL) AS spliceai_genome_missing,
  count(*) FILTER (WHERE in_exomes AND spliceai_exome IS NULL) AS spliceai_exome_missing,
  count(*) FILTER (WHERE in_genomes AND pangolin_genome IS NULL) AS pangolin_genome_missing,
  count(*) FILTER (WHERE in_exomes AND pangolin_exome IS NULL) AS pangolin_exome_missing,
  count(*) FILTER (WHERE NOT in_genomes AND spliceai_genome IS NOT NULL) AS spliceai_genome_without_call,
  count(*) FILTER (WHERE NOT in_exomes AND spliceai_exome IS NOT NULL) AS spliceai_exome_without_call,
  count(*) FILTER (WHERE NOT in_genomes AND pangolin_genome IS NOT NULL) AS pangolin_genome_without_call,
  count(*) FILTER (WHERE NOT in_exomes AND pangolin_exome IS NOT NULL) AS pangolin_exome_without_call
FROM variant_overlap;

SELECT
  score,
  n_null,
  n_negative,
  n_zero,
  n_gt_0_lt_0_2,
  n_ge_0_2_lt_0_5,
  n_ge_0_5_lt_0_8,
  n_ge_0_8,
  round(min_v, 4) AS min_v,
  round(p50, 4) AS p50,
  round(p90, 4) AS p90,
  round(max_v, 4) AS max_v
FROM (
  SELECT
    'spliceai_genome' AS score,
    count(*) FILTER (WHERE spliceai_genome IS NULL) AS n_null,
    count(*) FILTER (WHERE spliceai_genome < 0) AS n_negative,
    count(*) FILTER (WHERE spliceai_genome = 0) AS n_zero,
    count(*) FILTER (WHERE spliceai_genome > 0 AND spliceai_genome < 0.2) AS n_gt_0_lt_0_2,
    count(*) FILTER (WHERE spliceai_genome >= 0.2 AND spliceai_genome < 0.5) AS n_ge_0_2_lt_0_5,
    count(*) FILTER (WHERE spliceai_genome >= 0.5 AND spliceai_genome < 0.8) AS n_ge_0_5_lt_0_8,
    count(*) FILTER (WHERE spliceai_genome >= 0.8) AS n_ge_0_8,
    min(spliceai_genome) AS min_v,
    quantile_cont(spliceai_genome, 0.5) AS p50,
    quantile_cont(spliceai_genome, 0.9) AS p90,
    max(spliceai_genome) AS max_v
  FROM variant_overlap
  UNION ALL
  SELECT
    'spliceai_exome',
    count(*) FILTER (WHERE spliceai_exome IS NULL),
    count(*) FILTER (WHERE spliceai_exome < 0),
    count(*) FILTER (WHERE spliceai_exome = 0),
    count(*) FILTER (WHERE spliceai_exome > 0 AND spliceai_exome < 0.2),
    count(*) FILTER (WHERE spliceai_exome >= 0.2 AND spliceai_exome < 0.5),
    count(*) FILTER (WHERE spliceai_exome >= 0.5 AND spliceai_exome < 0.8),
    count(*) FILTER (WHERE spliceai_exome >= 0.8),
    min(spliceai_exome),
    quantile_cont(spliceai_exome, 0.5),
    quantile_cont(spliceai_exome, 0.9),
    max(spliceai_exome)
  FROM variant_overlap
  UNION ALL
  SELECT
    'pangolin_genome',
    count(*) FILTER (WHERE pangolin_genome IS NULL),
    count(*) FILTER (WHERE pangolin_genome < 0),
    count(*) FILTER (WHERE pangolin_genome = 0),
    count(*) FILTER (WHERE pangolin_genome > 0 AND pangolin_genome < 0.2),
    count(*) FILTER (WHERE pangolin_genome >= 0.2 AND pangolin_genome < 0.5),
    count(*) FILTER (WHERE pangolin_genome >= 0.5 AND pangolin_genome < 0.8),
    count(*) FILTER (WHERE pangolin_genome >= 0.8),
    min(pangolin_genome),
    quantile_cont(pangolin_genome, 0.5),
    quantile_cont(pangolin_genome, 0.9),
    max(pangolin_genome)
  FROM variant_overlap
  UNION ALL
  SELECT
    'pangolin_exome',
    count(*) FILTER (WHERE pangolin_exome IS NULL),
    count(*) FILTER (WHERE pangolin_exome < 0),
    count(*) FILTER (WHERE pangolin_exome = 0),
    count(*) FILTER (WHERE pangolin_exome > 0 AND pangolin_exome < 0.2),
    count(*) FILTER (WHERE pangolin_exome >= 0.2 AND pangolin_exome < 0.5),
    count(*) FILTER (WHERE pangolin_exome >= 0.5 AND pangolin_exome < 0.8),
    count(*) FILTER (WHERE pangolin_exome >= 0.8),
    min(pangolin_exome),
    quantile_cont(pangolin_exome, 0.5),
    quantile_cont(pangolin_exome, 0.9),
    max(pangolin_exome)
  FROM variant_overlap
);

SELECT
  count(*) FILTER (WHERE in_genomes AND in_exomes) AS in_both_vcfs,
  count(*) FILTER (WHERE spliceai_genome IS NOT NULL AND spliceai_exome IS NOT NULL) AS spliceai_both,
  count(*) FILTER (
    WHERE spliceai_genome IS NOT NULL AND spliceai_exome IS NOT NULL
      AND spliceai_genome = spliceai_exome
  ) AS spliceai_equal,
  count(*) FILTER (WHERE pangolin_genome IS NOT NULL AND pangolin_exome IS NOT NULL) AS pangolin_both,
  count(*) FILTER (
    WHERE pangolin_genome IS NOT NULL AND pangolin_exome IS NOT NULL
      AND pangolin_genome = pangolin_exome
  ) AS pangolin_equal,
  count(*) FILTER (
    WHERE in_genomes AND in_exomes AND spliceai_genome IS NOT NULL AND spliceai_exome IS NULL
  ) AS spliceai_genome_only_when_both_present,
  count(*) FILTER (
    WHERE in_genomes AND in_exomes AND spliceai_genome IS NULL AND spliceai_exome IS NOT NULL
  ) AS spliceai_exome_only_when_both_present,
  count(*) FILTER (
    WHERE in_genomes AND in_exomes AND pangolin_genome IS NOT NULL AND pangolin_exome IS NULL
  ) AS pangolin_genome_only_when_both_present,
  count(*) FILTER (
    WHERE in_genomes AND in_exomes AND pangolin_genome IS NULL AND pangolin_exome IS NOT NULL
  ) AS pangolin_exome_only_when_both_present,
  max(abs(spliceai_genome - spliceai_exome)) FILTER (
    WHERE spliceai_genome IS NOT NULL AND spliceai_exome IS NOT NULL
  ) AS spliceai_max_abs_diff,
  max(abs(pangolin_genome - pangolin_exome)) FILTER (
    WHERE pangolin_genome IS NOT NULL AND pangolin_exome IS NOT NULL
  ) AS pangolin_max_abs_diff
FROM variant_overlap;

SELECT germline_somatic, count(*) AS n
FROM variant_overlap
GROUP BY 1
ORDER BY n DESC;

SELECT var_type, count(*) AS n
FROM variant_overlap
GROUP BY 1
ORDER BY n DESC;

SELECT
  CASE
    WHEN clinical_significance IN (
      'Uncertain significance',
      'Likely benign',
      '-',
      'Conflicting classifications of pathogenicity',
      'Benign',
      'Benign/Likely benign',
      'Pathogenic',
      'Likely pathogenic',
      'Pathogenic/Likely pathogenic'
    ) THEN clinical_significance
    ELSE 'other'
  END AS clinical_group,
  count(*) AS n,
  count(*) FILTER (WHERE spliceai_genome IS NOT NULL) AS spliceai_genome,
  count(*) FILTER (WHERE (spliceai_genome >= 0.2) IS TRUE) AS spliceai_genome_ge_0_2,
  count(*) FILTER (WHERE (spliceai_genome >= 0.5) IS TRUE) AS spliceai_genome_ge_0_5,
  count(*) FILTER (WHERE spliceai_exome IS NOT NULL) AS spliceai_exome,
  count(*) FILTER (WHERE (spliceai_exome >= 0.2) IS TRUE) AS spliceai_exome_ge_0_2,
  count(*) FILTER (WHERE (spliceai_exome >= 0.5) IS TRUE) AS spliceai_exome_ge_0_5,
  count(*) FILTER (WHERE pangolin_genome IS NOT NULL) AS pangolin_genome,
  count(*) FILTER (WHERE pangolin_exome IS NOT NULL) AS pangolin_exome,
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_sgene = 1) AS sgene_one_ensg,
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_tested_not_sgene = 1) AS tested_not_sgene_one_ensg,
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_not_in_sgene_table = 1) AS not_in_sgene_table_one_ensg,
  count(*) FILTER (WHERE n_linked_ensg = 0) AS unlinked,
  count(*) FILTER (WHERE n_linked_ensg > 1) AS multi_ensg
FROM variant_overlap
GROUP BY 1
ORDER BY clinical_group = 'other', n DESC;

SELECT
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_sgene = 1) AS one_sgene,
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_tested_not_sgene = 1) AS one_tested_not,
  count(*) FILTER (WHERE n_linked_ensg = 1 AND n_linked_not_in_sgene_table = 1) AS one_absent,
  count(*) FILTER (WHERE n_linked_ensg = 0) AS unlinked,
  count(*) FILTER (WHERE n_linked_ensg = 0 AND gene_id = '-1') AS unlinked_gene_id_minus1,
  count(*) FILTER (WHERE n_linked_ensg > 1 AND n_linked_sgene = n_linked_ensg) AS multi_all_sgene,
  count(*) FILTER (WHERE n_linked_ensg > 1 AND n_linked_sgene = 0) AS multi_no_sgene,
  count(*) FILTER (
    WHERE n_linked_ensg > 1 AND n_linked_sgene > 0 AND n_linked_sgene < n_linked_ensg
  ) AS multi_disagree,
  count(*) FILTER (WHERE n_linked_ensg_allele_is_sig_lead > 0) AS allele_is_sig_lead_for_linked_gene,
  count(*) AS n
FROM variant_overlap;

SELECT
  count(*) FILTER (
    WHERE n_linked_ensg = 1
      AND ((spliceai_genome >= 0.2) IS TRUE OR (spliceai_exome >= 0.2) IS TRUE)
      AND n_linked_sgene = 1
  ) AS hi_sgene,
  count(*) FILTER (
    WHERE n_linked_ensg = 1
      AND ((spliceai_genome >= 0.2) IS TRUE OR (spliceai_exome >= 0.2) IS TRUE)
  ) AS hi_n,
  count(*) FILTER (
    WHERE n_linked_ensg = 1
      AND (spliceai_genome IS NOT NULL OR spliceai_exome IS NOT NULL)
      AND NOT ((spliceai_genome >= 0.2) IS TRUE OR (spliceai_exome >= 0.2) IS TRUE)
      AND n_linked_sgene = 1
  ) AS lo_sgene,
  count(*) FILTER (
    WHERE n_linked_ensg = 1
      AND (spliceai_genome IS NOT NULL OR spliceai_exome IS NOT NULL)
      AND NOT ((spliceai_genome >= 0.2) IS TRUE OR (spliceai_exome >= 0.2) IS TRUE)
  ) AS lo_n
FROM variant_overlap;

SELECT
  count(*) AS genes_with_a_scored_allele,
  count(*) FILTER (WHERE biotype = 'protein_coding') AS protein_coding,
  count(*) FILTER (WHERE is_sgene) AS sgenes,
  count(*) FILTER (WHERE in_sgene_table AND NOT is_sgene) AS tested_not_sgene,
  count(*) FILTER (WHERE NOT in_sgene_table) AS not_in_sgene_table,
  count(DISTINCT symbol) AS symbols
FROM gene_overlap
WHERE n_scored_variants > 0;

SELECT
  biotype,
  count(*) AS genes,
  count(*) FILTER (WHERE is_sgene) AS sgenes,
  count(*) FILTER (WHERE in_sgene_table AND NOT is_sgene) AS tested_not_sgene,
  count(*) FILTER (WHERE NOT in_sgene_table) AS not_in_sgene_table
FROM gene_overlap
WHERE n_scored_variants > 0
GROUP BY 1
ORDER BY genes DESC;

SELECT
  'sgenes_table' AS set_name,
  count(*) AS tested_genes,
  count(*) FILTER (WHERE n_sig_phenotypes > 0) AS sgenes,
  round(count(*) FILTER (WHERE n_sig_phenotypes > 0) / count(*), 6) AS sgene_fraction
FROM gene_status
UNION ALL
SELECT
  'protein_coding_in_crosswalk_and_sgenes_table',
  count(*),
  count(*) FILTER (WHERE s.n_sig_phenotypes > 0),
  round(count(*) FILTER (WHERE s.n_sig_phenotypes > 0) / count(*), 6)
FROM read_parquet('D:/variant-workbench/warehouse/crosswalk/genes.parquet') g
JOIN gene_status s USING (ensg)
WHERE g.biotype = 'protein_coding'
UNION ALL
SELECT
  'scored_splice_region_genes_in_sgenes_table',
  count(*),
  count(*) FILTER (WHERE is_sgene),
  round(count(*) FILTER (WHERE is_sgene) / count(*), 6)
FROM gene_overlap
WHERE n_scored_variants > 0 AND in_sgene_table
UNION ALL
SELECT
  'scored_splice_region_protein_coding_in_sgenes_table',
  count(*),
  count(*) FILTER (WHERE is_sgene),
  round(count(*) FILTER (WHERE is_sgene) / count(*), 6)
FROM gene_overlap
WHERE n_scored_variants > 0 AND in_sgene_table AND biotype = 'protein_coding';

SELECT
  count(*) AS sgenes_linked_to_scored_alleles,
  quantile_cont(n_sig_tissues, 0.5) AS median_sig_tissues,
  min(n_sig_tissues) AS min_sig_tissues,
  max(n_sig_tissues) AS max_sig_tissues
FROM gene_overlap
WHERE n_scored_variants > 0 AND is_sgene;

SELECT
  count(*) AS region_alleles_that_are_sig_leads_for_some_gene,
  count(*) FILTER (WHERE linked.variant_key IS NOT NULL) AS of_which_linked_gene,
  count(*) FILTER (WHERE linked.variant_key IS NULL) AS lead_for_a_different_gene
FROM (
  SELECT DISTINCT r.variant_key
  FROM region r
  JOIN sig_leads lead ON lead.variant_key = r.variant_key
) hits
LEFT JOIN (
  SELECT DISTINCT variant_key
  FROM variant_overlap
  WHERE n_linked_ensg_allele_is_sig_lead > 0
) linked ON linked.variant_key = hits.variant_key;

SELECT
  variant_key,
  gene_symbol,
  var_type,
  in_genomes,
  in_exomes,
  clinical_significance,
  n_linked_ensg
FROM variant_overlap
WHERE spliceai_genome IS NULL
  AND spliceai_exome IS NULL
  AND pangolin_genome IS NULL
  AND pangolin_exome IS NULL
ORDER BY variant_key;
