-- Coverage of gnomAD genome ancestry frequencies on Tier A GWAS variants.
-- Genome columns only. Do not add AF_afr_exome (or any AF_*_exome) to AF_*.
-- AN_afr, AN_nfe, AN_eas, AN_sas, and AN_amr are not columns in tier_a.parquet.
-- The allele-number filter is global genome AN >= 1000, not a per-group AN.
-- AF = 0 is non-null. FILTER = PASS is not required.
-- Run with DuckDB (multiple statements). The SELECTs are the numbers in
-- warehouse/reports/ancestry_audit.md.

SET memory_limit = '16GB';
SET threads = 8;
SET preserve_insertion_order = false;
SET temp_directory = 'D:/variant-workbench/warehouse/tmp';

CREATE OR REPLACE TEMP TABLE gwas AS
SELECT
  variant_key,
  source,
  chrom,
  pos,
  ref,
  alt,
  in_genomes,
  filter,
  AN,
  AF_afr,
  AF_nfe,
  AF_eas,
  AF_sas,
  AF_amr
FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet')
WHERE list_contains(string_split(source, ','), 'gwas');

CREATE OR REPLACE TEMP TABLE clin AS
SELECT variant_key, rsid, chrom, pos, ref, alt
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')
WHERE variant_key IN (SELECT variant_key FROM gwas);

CREATE OR REPLACE TEMP TABLE coord AS
SELECT
  c.variant_key,
  c.rsid,
  count(*) FILTER (WHERE a.chrom = c.chrom AND a.pos = c.pos) AS n_exact,
  min(abs(a.pos - c.pos)) FILTER (WHERE a.chrom = c.chrom AND a.pos IS NOT NULL) AS nearest,
  count(*) FILTER (WHERE a.pos IS NULL) AS n_null_pos,
  count(*) FILTER (WHERE a.chrom IS DISTINCT FROM c.chrom AND a.pos IS NOT NULL) AS n_other_chrom,
  (length(c.ref) > 1 OR length(c.alt) > 1) AS indel
FROM clin c
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet') a
  ON a.rsid = c.rsid
GROUP BY c.variant_key, c.rsid, c.ref, c.alt;

-- Join and denominator. associations.variant_key is null on every catalog row.
SELECT
  (SELECT count(*) FROM gwas) AS tier_a_gwas_keys,
  (SELECT count(DISTINCT variant_key) FROM gwas) AS distinct_keys,
  (SELECT count(*) FROM clin) AS clinvar_rows,
  (SELECT count(DISTINCT rsid) FROM clin) AS distinct_rsids,
  (SELECT count(*) FROM clin WHERE rsid IS NULL) AS null_rsid,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')) AS association_rows,
  (SELECT count(variant_key) FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')) AS association_variant_keys_nonnull,
  (SELECT count(*) FROM coord) AS keys_with_rsid_match,
  (SELECT count(*) FROM coord WHERE n_exact > 0) AS keys_same_chrom_pos,
  (SELECT count(*) FROM coord WHERE n_exact = 0 AND nearest = 1 AND indel) AS keys_indel_off_by_one,
  (SELECT count(*) FROM coord WHERE n_exact = 0 AND nearest = 1 AND NOT indel) AS keys_snv_off_by_one,
  (SELECT count(*) FROM coord WHERE n_exact = 0 AND nearest IS NULL AND n_null_pos > 0 AND n_other_chrom = 0) AS keys_null_catalog_pos_only,
  (SELECT count(*) FROM coord WHERE n_exact = 0 AND nearest IS NULL AND n_other_chrom > 0) AS keys_other_chrom_only;

-- Unfiltered coverage. Denominator is every Tier A GWAS key.
SELECT
  count(*) AS denominator,
  count(*) FILTER (WHERE NOT in_genomes) AS not_in_genomes,
  count(*) FILTER (WHERE AF_afr IS NULL AND AF_nfe IS NULL AND AF_eas IS NULL AND AF_sas IS NULL AND AF_amr IS NULL) AS missing_all_five,
  count(*) FILTER (WHERE
    ((AF_afr IS NULL)::int + (AF_nfe IS NULL)::int + (AF_eas IS NULL)::int + (AF_sas IS NULL)::int + (AF_amr IS NULL)::int)
    BETWEEN 1 AND 4) AS partial_missing,
  count(*) FILTER (WHERE AF_afr IS NOT NULL) AS afr_nonnull,
  count(*) FILTER (WHERE AF_afr IS NULL) AS afr_missing,
  count(*) FILTER (WHERE AF_nfe IS NOT NULL) AS nfe_nonnull,
  count(*) FILTER (WHERE AF_nfe IS NULL) AS nfe_missing,
  count(*) FILTER (WHERE AF_eas IS NOT NULL) AS eas_nonnull,
  count(*) FILTER (WHERE AF_eas IS NULL) AS eas_missing,
  count(*) FILTER (WHERE AF_sas IS NOT NULL) AS sas_nonnull,
  count(*) FILTER (WHERE AF_sas IS NULL) AS sas_missing,
  count(*) FILTER (WHERE AF_amr IS NOT NULL) AS amr_nonnull,
  count(*) FILTER (WHERE AF_amr IS NULL) AS amr_missing,
  count(*) FILTER (WHERE isnan(AF_afr) OR isnan(AF_nfe) OR isnan(AF_eas) OR isnan(AF_sas) OR isnan(AF_amr)) AS nan_rows,
  count(*) FILTER (WHERE AF_afr IS NOT NULL AND AF_nfe IS NOT NULL AND AF_eas IS NOT NULL AND AF_sas IS NOT NULL AND AF_amr IS NOT NULL) AS all_five_nonnull,
  round(count(*) FILTER (WHERE AF_afr IS NOT NULL) / count(*)::double, 6) AS afr_nonnull_fraction,
  round(count(*) FILTER (WHERE AF_afr IS NULL) / count(*)::double, 6) AS afr_missing_fraction,
  round(count(*) FILTER (WHERE AF_nfe IS NOT NULL) / count(*)::double, 6) AS nfe_nonnull_fraction,
  round(count(*) FILTER (WHERE AF_nfe IS NULL) / count(*)::double, 6) AS nfe_missing_fraction,
  round(count(*) FILTER (WHERE AF_eas IS NOT NULL) / count(*)::double, 6) AS eas_nonnull_fraction,
  round(count(*) FILTER (WHERE AF_eas IS NULL) / count(*)::double, 6) AS eas_missing_fraction,
  round(count(*) FILTER (WHERE AF_sas IS NOT NULL) / count(*)::double, 6) AS sas_nonnull_fraction,
  round(count(*) FILTER (WHERE AF_sas IS NULL) / count(*)::double, 6) AS sas_missing_fraction,
  round(count(*) FILTER (WHERE AF_amr IS NOT NULL) / count(*)::double, 6) AS amr_nonnull_fraction,
  round(count(*) FILTER (WHERE AF_amr IS NULL) / count(*)::double, 6) AS amr_missing_fraction,
  round(count(*) FILTER (WHERE AF_afr IS NOT NULL AND AF_nfe IS NOT NULL AND AF_eas IS NOT NULL AND AF_sas IS NOT NULL AND AF_amr IS NOT NULL) / count(*)::double, 6) AS all_five_nonnull_fraction
FROM gwas;

-- Filtered coverage. in_genomes and global AN >= 1000.
SELECT
  count(*) FILTER (WHERE in_genomes AND AN >= 1000) AS denominator,
  count(*) FILTER (WHERE NOT in_genomes OR AN IS NULL OR AN < 1000) AS excluded,
  count(*) FILTER (WHERE in_genomes AND AN < 1000) AS in_genomes_an_below_1000,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND filter IS DISTINCT FROM 'PASS') AS not_pass,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_afr IS NOT NULL) AS afr_nonnull,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_afr IS NULL) AS afr_missing,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_nfe IS NOT NULL) AS nfe_nonnull,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_nfe IS NULL) AS nfe_missing,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_eas IS NOT NULL) AS eas_nonnull,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_eas IS NULL) AS eas_missing,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_sas IS NOT NULL) AS sas_nonnull,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_sas IS NULL) AS sas_missing,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_amr IS NOT NULL) AS amr_nonnull,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_amr IS NULL) AS amr_missing,
  count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_afr IS NOT NULL AND AF_nfe IS NOT NULL AND AF_eas IS NOT NULL AND AF_sas IS NOT NULL AND AF_amr IS NOT NULL) AS all_five_nonnull,
  round(count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_sas IS NOT NULL) / count(*) FILTER (WHERE in_genomes AND AN >= 1000)::double, 6) AS sas_nonnull_fraction,
  round(count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_sas IS NULL) / count(*) FILTER (WHERE in_genomes AND AN >= 1000)::double, 6) AS sas_missing_fraction,
  round(count(*) FILTER (WHERE in_genomes AND AN >= 1000 AND AF_afr IS NOT NULL AND AF_nfe IS NOT NULL AND AF_eas IS NOT NULL AND AF_sas IS NOT NULL AND AF_amr IS NOT NULL) / count(*) FILTER (WHERE in_genomes AND AN >= 1000)::double, 6) AS all_five_nonnull_fraction
FROM gwas;

-- Genome FILTER values inside the allele-number filter. Not an exclusion.
SELECT coalesce(filter, '<null>') AS genome_filter, count(*) AS n
FROM gwas
WHERE in_genomes AND AN >= 1000 AND filter IS DISTINCT FROM 'PASS'
GROUP BY 1
ORDER BY 2 DESC, 1;

-- In-genome keys with a null ancestry AF. AF values are not selected.
SELECT
  variant_key,
  filter AS genome_filter,
  AN,
  (AN >= 1000) AS passes_an_filter,
  (AF_afr IS NULL) AS af_afr_null,
  (AF_nfe IS NULL) AS af_nfe_null,
  (AF_eas IS NULL) AS af_eas_null,
  (AF_sas IS NULL) AS af_sas_null,
  (AF_amr IS NULL) AS af_amr_null
FROM gwas
WHERE in_genomes
  AND (AF_afr IS NULL OR AF_nfe IS NULL OR AF_eas IS NULL OR AF_sas IS NULL OR AF_amr IS NULL)
ORDER BY AN, variant_key;

-- rsID matches that are not on the Tier A chromosome. Counts are in the join SELECT.
SELECT
  c.variant_key,
  c.rsid,
  c.chrom AS tier_a_chrom,
  c.pos AS tier_a_pos,
  string_agg(DISTINCT a.chrom, ',') FILTER (WHERE a.pos IS NOT NULL) AS catalog_chroms
FROM clin c
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet') a
  ON a.rsid = c.rsid
WHERE c.variant_key IN (
  SELECT variant_key
  FROM coord
  WHERE n_exact = 0 AND nearest IS NULL AND n_other_chrom > 0
)
GROUP BY c.variant_key, c.rsid, c.chrom, c.pos
ORDER BY c.variant_key;
