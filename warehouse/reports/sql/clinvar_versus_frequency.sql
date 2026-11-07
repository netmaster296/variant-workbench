-- ClinVar clinical significance against gnomAD v4.1 frequency.
-- Alleles: warehouse/clinvar/allele_dedup.parquet (one row per variant_key,
-- or per AlleleID when there is no VCF allele).
-- Frequencies: warehouse/variants/tier_a.parquet.
--   AF is the genome allele frequency. AF_exome is the exome allele frequency.
--   They are never added, averaged, or otherwise combined.
-- Names for the common-allele list: warehouse/clinvar/alleles.parquet.
--   That join is checked for row multiplication.
--
-- Clinical class, first match on the lowered ClinicalSignificance text:
--   conflicting          text contains 'conflict' (wins over pathogenic)
--   pathogenic           'pathogenic' remains after removing 'likely pathogenic'
--                        so Pathogenic/Likely pathogenic counts here, not in both
--   likely_pathogenic    text contains 'likely pathogenic' and the row is not above
--   benign               'benign' remains after removing 'likely benign'
--   likely_benign        text contains 'likely benign'
--   vus                  uncertain significance, or vus
--   drug_response        drug response, and none of the above
--   risk_or_other_effect risk factor, risk allele, association, protective,
--                        affects, or confers sensitivity
--   no_classification    blank, '-', not provided, no classification,
--                        no assertion, no interpretation
--   other                anything else
-- Low penetrance, drug response, and risk factor are also flags. A pathogenic
-- row that also says drug response stays in pathogenic.
-- germline_pathogenic is the column loaded with the dedup table. It is not
-- recomputed here. It is true only for a germline origin, no oncogenicity
-- call, and a pathogenic or likely pathogenic piece with no conflicting piece.
--
-- AF bins are half-open and identical for the two datasets:
--   null, <0.0001, 0.0001–<0.001, 0.001–<0.01, 0.01–<0.05, >=0.05.
-- AF = 0 is inside <0.0001. It is a measured zero, not a missing AF.
-- The common-allele review is every pathogenic or likely pathogenic allele
-- with AF >= 0.01 or AF_exome >= 0.01. That set is the whole tail, not a draw.

CREATE OR REPLACE MACRO af_bin(x) AS (
  CASE
    WHEN x IS NULL THEN 'null'
    WHEN x < 0.0001 THEN '<0.0001'
    WHEN x < 0.001 THEN '0.0001_to_0.001'
    WHEN x < 0.01 THEN '0.001_to_0.01'
    WHEN x < 0.05 THEN '0.01_to_0.05'
    ELSE '>=0.05'
  END
);

CREATE OR REPLACE TEMP VIEW clinvar_af AS
WITH d AS (
  SELECT
    variant_key,
    n_rows,
    "AlleleID" AS allele_id,
    "GeneSymbol" AS gene,
    "Type" AS var_type,
    chrom,
    pos,
    ref,
    alt,
    rsid,
    "ClinicalSignificance" AS sig_raw,
    "ReviewStatus" AS review_status,
    review_stars,
    "OriginSimple" AS origin_raw,
    germline_somatic AS origin_class,
    "Oncogenicity" AS oncogenicity,
    germline_pathogenic,
    substr(replace(replace(coalesce("PhenotypeList", ''), chr(10), ' '), chr(9), ' '), 1, 180) AS phenotype,
    lower(coalesce("ClinicalSignificance", '')) AS sig
  FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')
),
t AS (
  SELECT
    variant_key,
    chrom AS tier_chrom,
    pos AS tier_pos,
    ref AS tier_ref,
    alt AS tier_alt,
    in_genomes,
    in_exomes,
    AF AS af_g,
    AF_exome AS af_e,
    AC AS ac_g,
    AN AS an_g,
    AC_exome AS ac_e,
    AN_exome AS an_e,
    filter AS filt_g,
    filter_exome AS filt_e
  FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet')
)
SELECT
  d.*,
  t.tier_chrom,
  t.tier_pos,
  t.tier_ref,
  t.tier_alt,
  t.in_genomes,
  t.in_exomes,
  t.af_g,
  t.af_e,
  t.ac_g,
  t.an_g,
  t.ac_e,
  t.an_e,
  t.filt_g,
  t.filt_e
FROM d
LEFT JOIN t USING (variant_key);

CREATE OR REPLACE TEMP VIEW classed AS
SELECT
  *,
  (sig LIKE '%low penetrance%') AS low_pen,
  (sig LIKE '%drug response%') AS drug,
  (sig LIKE '%risk factor%' OR sig LIKE '%risk allele%') AS risk,
  (lower(coalesce(review_status, '')) LIKE '%conflicting classifications%') AS review_conflict,
  CASE
    WHEN sig LIKE '%conflict%' THEN 'conflicting'
    WHEN regexp_replace(sig, 'likely pathogenic', '', 'g') LIKE '%pathogenic%' THEN 'pathogenic'
    WHEN sig LIKE '%likely pathogenic%' THEN 'likely_pathogenic'
    WHEN regexp_replace(sig, 'likely benign', '', 'g') LIKE '%benign%' THEN 'benign'
    WHEN sig LIKE '%likely benign%' THEN 'likely_benign'
    WHEN sig LIKE '%uncertain significance%' OR sig LIKE '%vus%' THEN 'vus'
    WHEN sig LIKE '%drug response%' THEN 'drug_response'
    WHEN sig LIKE '%risk factor%'
      OR sig LIKE '%risk allele%'
      OR sig LIKE '%association%'
      OR sig LIKE '%protective%'
      OR sig LIKE '%affects%'
      OR sig LIKE '%confers sensitivity%' THEN 'risk_or_other_effect'
    WHEN sig IN ('', '-', 'not provided')
      OR sig LIKE 'no classification%'
      OR sig LIKE 'no assertion%'
      OR sig LIKE 'no interpretation%' THEN 'no_classification'
    ELSE 'other'
  END AS sig_class,
  (
    variant_key IS NOT NULL
    AND (
      ref IS DISTINCT FROM tier_ref
      OR alt IS DISTINCT FROM tier_alt
      OR chrom IS DISTINCT FROM tier_chrom
      OR pos IS DISTINCT FROM tier_pos
    )
  ) AS allele_mismatch
FROM clinvar_af;

CREATE OR REPLACE TEMP VIEW plp AS
SELECT *
FROM classed
WHERE sig_class IN ('pathogenic', 'likely_pathogenic');

-- Join audit. joined_rows must equal dedup_rows. A keyed miss, an allele
-- mismatch, or a duplicated variant_key stops the claim.
SELECT
  'join_audit' AS result,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')) AS dedup_rows,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') WHERE variant_key IS NULL) AS dedup_null_keys,
  (SELECT count(*) - count(DISTINCT variant_key) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') WHERE variant_key IS NOT NULL) AS dedup_extra_key_rows,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet')) AS tier_rows,
  (SELECT count(*) - count(DISTINCT variant_key) FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet')) AS tier_extra_key_rows,
  (SELECT count(*) FROM classed) AS joined_rows,
  (SELECT count(*) FROM classed WHERE variant_key IS NOT NULL AND in_genomes IS NULL) AS keyed_missing_from_tier,
  (SELECT count(*) FROM classed WHERE allele_mismatch) AS allele_mismatch,
  (SELECT count(*) FROM classed WHERE in_genomes AND af_g IS NULL) AS in_genomes_af_null,
  (SELECT count(*) FROM classed WHERE in_exomes AND af_e IS NULL) AS in_exomes_af_null,
  (SELECT count(*) FROM classed WHERE NOT coalesce(in_genomes, false) AND af_g IS NOT NULL) AS af_without_genome_hit,
  (SELECT count(*) FROM classed WHERE NOT coalesce(in_exomes, false) AND af_e IS NOT NULL) AS af_without_exome_hit;

SELECT
  'flag_agreement' AS result,
  count(*) FILTER (WHERE germline_pathogenic AND sig_class NOT IN ('pathogenic', 'likely_pathogenic')) AS flag_outside_plp,
  count(*) FILTER (WHERE sig_class IN ('pathogenic', 'likely_pathogenic')) AS text_plp,
  count(*) FILTER (WHERE germline_pathogenic) AS germ_flag,
  count(*) FILTER (WHERE sig_class = 'conflicting' AND sig LIKE '%pathogenic%') AS conflicting_text_also_says_pathogenic
FROM classed;

SELECT
  'class_bin' AS result,
  dataset,
  sig_class,
  bin,
  n
FROM (
  SELECT 'genome' AS dataset, sig_class, af_bin(af_g) AS bin, count(*) AS n
  FROM classed
  GROUP BY 1, 2, 3
  UNION ALL
  SELECT 'exome', sig_class, af_bin(af_e), count(*)
  FROM classed
  GROUP BY 1, 2, 3
  UNION ALL
  SELECT 'genome', 'pathogenic_or_likely', af_bin(af_g), count(*)
  FROM plp
  GROUP BY 1, 2, 3
  UNION ALL
  SELECT 'exome', 'pathogenic_or_likely', af_bin(af_e), count(*)
  FROM plp
  GROUP BY 1, 2, 3
  UNION ALL
  SELECT 'genome', 'germline_pathogenic_flag', af_bin(af_g), count(*)
  FROM classed
  WHERE germline_pathogenic
  GROUP BY 1, 2, 3
  UNION ALL
  SELECT 'exome', 'germline_pathogenic_flag', af_bin(af_e), count(*)
  FROM classed
  WHERE germline_pathogenic
  GROUP BY 1, 2, 3
);

SELECT
  'thresholds' AS result,
  label,
  dataset,
  n,
  n_null,
  n_af,
  n_lt_0_0001,
  n_lt_0_001,
  n_lt_0_01,
  n_ge_0_01,
  n_ge_0_05,
  round(n_lt_0_0001 * 1.0 / nullif(n_af, 0), 6) AS frac_af_lt_0_0001,
  round(n_lt_0_001 * 1.0 / nullif(n_af, 0), 6) AS frac_af_lt_0_001,
  round(n_lt_0_01 * 1.0 / nullif(n_af, 0), 6) AS frac_af_lt_0_01,
  round(n_ge_0_01 * 1.0 / nullif(n_af, 0), 6) AS frac_af_ge_0_01
FROM (
  SELECT
    label,
    'genome' AS dataset,
    count(*) AS n,
    count(*) FILTER (WHERE af_g IS NULL) AS n_null,
    count(*) FILTER (WHERE af_g IS NOT NULL) AS n_af,
    count(*) FILTER (WHERE af_g IS NOT NULL AND af_g < 0.0001) AS n_lt_0_0001,
    count(*) FILTER (WHERE af_g IS NOT NULL AND af_g < 0.001) AS n_lt_0_001,
    count(*) FILTER (WHERE af_g IS NOT NULL AND af_g < 0.01) AS n_lt_0_01,
    count(*) FILTER (WHERE af_g >= 0.01) AS n_ge_0_01,
    count(*) FILTER (WHERE af_g >= 0.05) AS n_ge_0_05
  FROM (
    SELECT 'pathogenic' AS label, af_g, af_e, sig_class, germline_pathogenic FROM classed WHERE sig_class = 'pathogenic'
    UNION ALL
    SELECT 'likely_pathogenic', af_g, af_e, sig_class, germline_pathogenic FROM classed WHERE sig_class = 'likely_pathogenic'
    UNION ALL
    SELECT 'pathogenic_or_likely', af_g, af_e, sig_class, germline_pathogenic FROM plp
    UNION ALL
    SELECT 'germline_flag', af_g, af_e, sig_class, germline_pathogenic FROM classed WHERE germline_pathogenic
    UNION ALL
    SELECT 'text_plp_not_germline_flag', af_g, af_e, sig_class, germline_pathogenic FROM plp WHERE NOT germline_pathogenic
    UNION ALL
    SELECT sig_class, af_g, af_e, sig_class, germline_pathogenic FROM classed WHERE sig_class NOT IN ('pathogenic', 'likely_pathogenic')
  )
  GROUP BY 1
  UNION ALL
  SELECT
    label,
    'exome',
    count(*),
    count(*) FILTER (WHERE af_e IS NULL),
    count(*) FILTER (WHERE af_e IS NOT NULL),
    count(*) FILTER (WHERE af_e IS NOT NULL AND af_e < 0.0001),
    count(*) FILTER (WHERE af_e IS NOT NULL AND af_e < 0.001),
    count(*) FILTER (WHERE af_e IS NOT NULL AND af_e < 0.01),
    count(*) FILTER (WHERE af_e >= 0.01),
    count(*) FILTER (WHERE af_e >= 0.05)
  FROM (
    SELECT 'pathogenic' AS label, af_g, af_e FROM classed WHERE sig_class = 'pathogenic'
    UNION ALL
    SELECT 'likely_pathogenic', af_g, af_e FROM classed WHERE sig_class = 'likely_pathogenic'
    UNION ALL
    SELECT 'pathogenic_or_likely', af_g, af_e FROM plp
    UNION ALL
    SELECT 'germline_flag', af_g, af_e FROM classed WHERE germline_pathogenic
    UNION ALL
    SELECT 'text_plp_not_germline_flag', af_g, af_e FROM plp WHERE NOT germline_pathogenic
    UNION ALL
    SELECT sig_class, af_g, af_e FROM classed WHERE sig_class NOT IN ('pathogenic', 'likely_pathogenic')
  )
  GROUP BY 1
);

SELECT
  'null_reason' AS result,
  sig_class,
  count(*) AS n,
  count(*) FILTER (WHERE variant_key IS NULL) AS no_variant_key,
  count(*) FILTER (WHERE variant_key IS NOT NULL AND NOT coalesce(in_genomes, false)) AS not_in_genomes,
  count(*) FILTER (WHERE in_genomes AND af_g IS NULL) AS in_genomes_af_null,
  count(*) FILTER (WHERE af_g IS NOT NULL) AS genome_af,
  count(*) FILTER (WHERE variant_key IS NOT NULL AND NOT coalesce(in_exomes, false)) AS not_in_exomes,
  count(*) FILTER (WHERE in_exomes AND af_e IS NULL) AS in_exomes_af_null,
  count(*) FILTER (WHERE af_e IS NOT NULL) AS exome_af,
  count(*) FILTER (WHERE af_g IS NULL AND chrom = 'chrM') AS genome_null_chrm,
  count(*) FILTER (WHERE germline_pathogenic) AS germ_flag
FROM plp
GROUP BY sig_class
ORDER BY sig_class;

SELECT
  'null_reason_combined' AS result,
  count(*) AS n,
  count(*) FILTER (WHERE variant_key IS NULL) AS no_variant_key,
  count(*) FILTER (WHERE variant_key IS NOT NULL AND NOT coalesce(in_genomes, false)) AS not_in_genomes,
  count(*) FILTER (WHERE in_genomes AND af_g IS NULL) AS in_genomes_af_null,
  count(*) FILTER (WHERE af_g IS NOT NULL) AS genome_af,
  count(*) FILTER (WHERE af_g IS NULL) AS genome_null,
  count(*) FILTER (WHERE variant_key IS NOT NULL AND NOT coalesce(in_exomes, false)) AS not_in_exomes,
  count(*) FILTER (WHERE in_exomes AND af_e IS NULL) AS in_exomes_af_null,
  count(*) FILTER (WHERE af_e IS NOT NULL) AS exome_af,
  count(*) FILTER (WHERE af_e IS NULL) AS exome_null,
  count(*) FILTER (WHERE af_g IS NOT NULL AND af_e IS NOT NULL) AS both_af,
  count(*) FILTER (WHERE af_g IS NOT NULL AND af_e IS NULL) AS genome_only,
  count(*) FILTER (WHERE af_g IS NULL AND af_e IS NOT NULL) AS exome_only,
  count(*) FILTER (WHERE af_g IS NULL AND af_e IS NULL) AS neither_af,
  count(*) FILTER (WHERE af_g IS NULL AND chrom = 'chrM') AS genome_null_chrm
FROM plp;

SELECT
  'null_type' AS result,
  which,
  coalesce(var_type, '(null)') AS var_type,
  n
FROM (
  SELECT 'genome_null' AS which, var_type, count(*) AS n
  FROM plp
  WHERE af_g IS NULL
  GROUP BY 1, 2
  UNION ALL
  SELECT 'genome_af', var_type, count(*)
  FROM plp
  WHERE af_g IS NOT NULL
  GROUP BY 1, 2
  UNION ALL
  SELECT 'exome_null', var_type, count(*)
  FROM plp
  WHERE af_e IS NULL
  GROUP BY 1, 2
  UNION ALL
  SELECT 'exome_af', var_type, count(*)
  FROM plp
  WHERE af_e IS NOT NULL
  GROUP BY 1, 2
)
ORDER BY which, n DESC;

SELECT
  'measured_zero' AS result,
  count(*) FILTER (WHERE af_g = 0) AS genome_af_zero,
  count(*) FILTER (WHERE af_g = 0 AND ac_g = 0) AS genome_ac_zero,
  count(*) FILTER (WHERE af_g = 0 AND filt_g LIKE '%AC0%') AS genome_filter_ac0,
  count(*) FILTER (WHERE af_e = 0) AS exome_af_zero,
  count(*) FILTER (WHERE af_e = 0 AND ac_e = 0) AS exome_ac_zero,
  count(*) FILTER (WHERE af_e = 0 AND filt_e LIKE '%AC0%') AS exome_filter_ac0,
  count(*) FILTER (WHERE af_g >= 0.01 AND an_g < 10000) AS genome_ge_0_01_an_lt_10000,
  count(*) FILTER (WHERE af_e >= 0.01 AND an_e < 10000) AS exome_ge_0_01_an_lt_10000,
  min(an_g) FILTER (WHERE af_g >= 0.01) AS min_an_genome_ge_0_01,
  min(an_e) FILTER (WHERE af_e >= 0.01) AS min_an_exome_ge_0_01
FROM plp;

SELECT
  'in_dataset_af_null' AS result,
  variant_key,
  gene,
  sig_raw,
  var_type,
  germline_pathogenic,
  review_stars,
  in_genomes,
  af_g,
  ac_g,
  an_g,
  filt_g,
  in_exomes,
  af_e,
  ac_e,
  an_e,
  filt_e
FROM plp
WHERE (in_genomes AND af_g IS NULL) OR (in_exomes AND af_e IS NULL);

SELECT
  'mid_band' AS result,
  dataset,
  n,
  n_germ,
  n_low_pen,
  n_drug,
  n_risk,
  n_review_conflict,
  n_stars_0,
  n_stars_1,
  n_stars_ge2,
  n_not_germline_origin,
  n_pass
FROM (
  SELECT
    'genome' AS dataset,
    count(*) AS n,
    count(*) FILTER (WHERE germline_pathogenic) AS n_germ,
    count(*) FILTER (WHERE low_pen) AS n_low_pen,
    count(*) FILTER (WHERE drug) AS n_drug,
    count(*) FILTER (WHERE risk) AS n_risk,
    count(*) FILTER (WHERE review_conflict) AS n_review_conflict,
    count(*) FILTER (WHERE review_stars = 0) AS n_stars_0,
    count(*) FILTER (WHERE review_stars = 1) AS n_stars_1,
    count(*) FILTER (WHERE review_stars >= 2) AS n_stars_ge2,
    count(*) FILTER (WHERE origin_class <> 'germline') AS n_not_germline_origin,
    count(*) FILTER (WHERE filt_g = 'PASS') AS n_pass
  FROM plp
  WHERE af_g >= 0.001 AND af_g < 0.01
  UNION ALL
  SELECT
    'exome',
    count(*),
    count(*) FILTER (WHERE germline_pathogenic),
    count(*) FILTER (WHERE low_pen),
    count(*) FILTER (WHERE drug),
    count(*) FILTER (WHERE risk),
    count(*) FILTER (WHERE review_conflict),
    count(*) FILTER (WHERE review_stars = 0),
    count(*) FILTER (WHERE review_stars = 1),
    count(*) FILTER (WHERE review_stars >= 2),
    count(*) FILTER (WHERE origin_class <> 'germline'),
    count(*) FILTER (WHERE filt_e = 'PASS')
  FROM plp
  WHERE af_e >= 0.001 AND af_e < 0.01
);

CREATE OR REPLACE TEMP VIEW exceptions AS
SELECT
  *,
  (
    (af_g >= 0.01 AND filt_g = 'PASS')
    OR (af_e >= 0.01 AND filt_e = 'PASS')
  ) AS pass_cross,
  CASE
    WHEN NOT germline_pathogenic THEN 'not_germline_flag'
    WHEN (
      (af_g >= 0.01 AND filt_g = 'PASS')
      OR (af_e >= 0.01 AND filt_e = 'PASS')
    ) AND review_stars >= 2 THEN 'pass_stars_ge2'
    WHEN (
      (af_g >= 0.01 AND filt_g = 'PASS')
      OR (af_e >= 0.01 AND filt_e = 'PASS')
    ) AND review_stars = 1 THEN 'pass_stars_1'
    WHEN (
      (af_g >= 0.01 AND filt_g = 'PASS')
      OR (af_e >= 0.01 AND filt_e = 'PASS')
    ) THEN 'pass_stars_0_or_null'
    ELSE 'nonpass_only'
  END AS exception_group
FROM plp
WHERE af_g >= 0.01 OR af_e >= 0.01;

SELECT
  'exception_summary' AS result,
  count(*) AS n,
  count(DISTINCT variant_key) AS n_keys,
  count(*) FILTER (WHERE af_g >= 0.01) AS genome_ge_0_01,
  count(*) FILTER (WHERE af_e >= 0.01) AS exome_ge_0_01,
  count(*) FILTER (WHERE af_g >= 0.01 AND af_e >= 0.01) AS both_ge_0_01,
  count(*) FILTER (WHERE af_g >= 0.05 OR af_e >= 0.05) AS either_ge_0_05,
  count(*) FILTER (WHERE germline_pathogenic) AS germ_flag,
  count(*) FILTER (WHERE NOT germline_pathogenic) AS not_germ_flag,
  count(*) FILTER (WHERE origin_class = 'germline') AS origin_germline,
  count(*) FILTER (WHERE origin_class = 'somatic') AS origin_somatic,
  count(*) FILTER (WHERE origin_class = 'germline/somatic') AS origin_mixed,
  count(*) FILTER (WHERE origin_class = 'other') AS origin_other,
  count(*) FILTER (WHERE low_pen) AS low_pen,
  count(*) FILTER (WHERE drug) AS drug,
  count(*) FILTER (WHERE risk) AS risk,
  count(*) FILTER (WHERE review_conflict) AS review_conflict,
  count(*) FILTER (WHERE sig_class = 'conflicting') AS conflict_class,
  count(*) FILTER (WHERE allele_mismatch) AS allele_mismatch,
  count(*) FILTER (WHERE n_rows > 1) AS multi_source_rows,
  count(*) FILTER (WHERE review_stars = 0) AS stars_0,
  count(*) FILTER (WHERE review_stars = 1) AS stars_1,
  count(*) FILTER (WHERE review_stars = 2) AS stars_2,
  count(*) FILTER (WHERE review_stars = 3) AS stars_3,
  count(*) FILTER (WHERE review_stars = 4) AS stars_4,
  count(*) FILTER (WHERE review_stars IS NULL) AS stars_null,
  count(*) FILTER (WHERE pass_cross) AS pass_on_a_crossing_call,
  count(*) FILTER (WHERE chrom = 'chrX') AS chrom_x,
  count(*) FILTER (WHERE chrom = 'chrY') AS chrom_y
FROM exceptions;

SELECT
  'exception_group' AS result,
  exception_group,
  count(*) AS n,
  count(*) FILTER (WHERE low_pen) AS low_pen,
  count(*) FILTER (WHERE drug) AS drug,
  count(*) FILTER (WHERE risk) AS risk,
  count(*) FILTER (WHERE review_conflict) AS review_conflict,
  count(*) FILTER (WHERE sig_class = 'likely_pathogenic') AS likely_only
FROM exceptions
GROUP BY exception_group
ORDER BY exception_group;

SELECT
  'exception_names' AS result,
  (SELECT count(*) FROM exceptions) AS exception_rows,
  count(*) AS name_join_rows,
  count(DISTINCT e.variant_key) AS distinct_keys,
  count(*) FILTER (WHERE a.hgvs_name IS NULL) AS missing_name
FROM exceptions e
LEFT JOIN (
  SELECT variant_key, "Name" AS hgvs_name
  FROM read_parquet('D:/variant-workbench/warehouse/clinvar/alleles.parquet')
) a USING (variant_key);

SELECT
  'exception_row' AS result,
  e.exception_group,
  e.sig_class,
  e.germline_pathogenic,
  e.origin_class,
  e.review_stars,
  e.review_status,
  e.gene,
  e.chrom,
  e.rsid,
  a.hgvs_name,
  e.sig_raw,
  e.low_pen,
  e.drug,
  e.risk,
  e.var_type,
  e.af_g,
  e.ac_g,
  e.an_g,
  e.filt_g,
  e.af_e,
  e.ac_e,
  e.an_e,
  e.filt_e,
  (e.ref IS NOT DISTINCT FROM e.tier_ref AND e.alt IS NOT DISTINCT FROM e.tier_alt) AS allele_match,
  e.n_rows,
  e.phenotype
FROM exceptions e
LEFT JOIN (
  SELECT variant_key, "Name" AS hgvs_name
  FROM read_parquet('D:/variant-workbench/warehouse/clinvar/alleles.parquet')
) a USING (variant_key)
ORDER BY e.exception_group, coalesce(e.af_g, -1) DESC, coalesce(e.af_e, -1) DESC;
