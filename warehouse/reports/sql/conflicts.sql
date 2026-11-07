-- Where the catalogs disagree.
-- Read-only. Does not write a mart.
-- Genome AF is tier_a.AF. Exome AF is tier_a.AF_exome. They are never added.
-- Population AF columns are not used.
-- associations.variant_key is null on every GWAS row. This file does not join on it.
-- A GWAS ref/alt key is a ClinVar variant_key whose rsid, chrom, and pos match
-- the association and whose ref or alt equals the single risk-allele token.
-- MaveDB "mild" is computed inside one score set. Scores are not pooled.

-- ---------------------------------------------------------------------------
-- 1. Germline pathogenic or likely pathogenic alleles with a high genome AF
-- ---------------------------------------------------------------------------

-- germline_pathogenic is the dedup flag. It is true for pathogenic and for
-- likely pathogenic, and it is false for somatic origin, mixed origin,
-- an oncogenicity call, and any significance text that contains a conflict.

SELECT
  count(*) FILTER (WHERE germline_pathogenic AND variant_key IS NOT NULL) AS gplp_keyed,
  count(*) FILTER (WHERE germline_pathogenic AND variant_key IS NULL) AS gplp_no_key
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet');

SELECT
  count(*) AS joined,
  count(*) FILTER (WHERE g.ref = t.ref AND g.alt = t.alt) AS ref_alt_equal,
  count(*) FILTER (WHERE g.ref IS DISTINCT FROM t.ref OR g.alt IS DISTINCT FROM t.alt) AS ref_alt_differ,
  count(*) FILTER (WHERE t.in_genomes) AS in_genomes,
  count(*) FILTER (WHERE t.in_genomes AND t.AF IS NOT NULL) AS genome_af_measured,
  count(*) FILTER (WHERE t.in_genomes AND t.AF IS NULL) AS in_genomes_af_null,
  count(*) FILTER (WHERE NOT t.in_genomes) AS not_in_genomes,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.05) AS af_ge_0_05,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.01) AS af_ge_0_01,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.001) AS af_ge_0_001,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.05 AND t.filter = 'PASS' AND g.review_stars >= 2) AS pass_star2_ge_0_05,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.01 AND t.filter = 'PASS' AND g.review_stars >= 2) AS pass_star2_ge_0_01,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.001 AND t.filter = 'PASS' AND g.review_stars >= 2) AS pass_star2_ge_0_001,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.01 AND (
      g.ClinicalSignificance ILIKE '%low penetrance%'
      OR g.ClinicalSignificance ILIKE '%risk factor%'
      OR g.ClinicalSignificance ILIKE '%susceptibility%'
  )) AS qualified_label_ge_0_01,
  count(*) FILTER (WHERE t.in_genomes AND t.AF >= 0.01 AND t.in_exomes AND t.AF_exome IS NOT NULL AND t.AF >= 10 * t.AF_exome) AS exome_af_10x_lower_ge_0_01,
  count(*) FILTER (
    WHERE t.in_genomes AND t.AF IS NOT NULL AND t.AN > 0
      AND abs(t.AF - t.AC * 1.0 / t.AN) > 1e-4
  ) AS af_disagrees_with_ac_over_an
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') g
JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
  ON g.variant_key = t.variant_key
WHERE g.germline_pathogenic
  AND g.variant_key IS NOT NULL;

-- Alleles at or above genome AF 0.01. Exome AF is a separate column.
SELECT
  g.variant_key,
  g.GeneSymbol,
  g.ClinicalSignificance,
  g.review_stars,
  t.filter,
  t.AF AS genome_af,
  t.AC,
  t.AN,
  t.in_exomes,
  t.AF_exome AS exome_af,
  t.vep_consequence
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') g
JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
  ON g.variant_key = t.variant_key
WHERE g.germline_pathogenic
  AND t.in_genomes
  AND t.AF >= 0.01
ORDER BY t.AF DESC;

-- ---------------------------------------------------------------------------
-- 2. GWAS risk alleles that are the common allele, only with a ref/alt key
-- ---------------------------------------------------------------------------

-- Ready rows: one risk-allele token, an rsid, and a position.
-- variant_key on the association stays null and is not a join key.

WITH a AS (
  SELECT
    rsid,
    chrom,
    pos,
    upper(risk_allele) AS risk,
    try_cast("RISK ALLELE FREQUENCY" AS DOUBLE) AS raf
  FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')
  WHERE rsid IS NOT NULL
    AND risk_allele IS NOT NULL
    AND position(';' IN risk_allele) = 0
    AND chrom IS NOT NULL
    AND pos IS NOT NULL
),
c AS (
  SELECT rsid, variant_key, chrom, pos, ref, alt
  FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')
  WHERE rsid IS NOT NULL AND variant_key IS NOT NULL
),
pairs AS (
  SELECT
    a.rsid,
    a.risk,
    a.raf,
    c.variant_key,
    c.ref,
    c.alt,
    (a.risk = c.alt) AS is_alt,
    (a.risk = c.ref) AS is_ref
  FROM a
  JOIN c
    ON a.rsid = c.rsid
   AND a.chrom = c.chrom
   AND a.pos = c.pos
   AND (a.risk = c.ref OR a.risk = c.alt)
),
key_n AS (
  SELECT rsid, risk, count(DISTINCT variant_key) AS n_keys
  FROM pairs
  GROUP BY 1, 2
),
one AS (
  SELECT
    p.rsid,
    p.risk,
    any_value(p.variant_key) AS variant_key,
    any_value(p.ref) AS ref,
    any_value(p.alt) AS alt,
    bool_or(p.is_alt) AS is_alt,
    bool_or(p.is_ref) AS is_ref,
    count(*) AS n_assoc,
    median(p.raf) AS median_raf,
    max(p.raf) AS max_raf
  FROM pairs p
  JOIN key_n k ON p.rsid = k.rsid AND p.risk = k.risk
  WHERE k.n_keys = 1
  GROUP BY p.rsid, p.risk
),
j AS (
  SELECT o.*, t.AF, t.AN, t.in_genomes
  FROM one o
  LEFT JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
    ON o.variant_key = t.variant_key
)
SELECT
  (SELECT count(*) FROM a) AS assoc_ready,
  (SELECT count(DISTINCT rsid) FROM a) AS rsids_ready,
  (SELECT count(*) FROM key_n) AS matched_claims,
  (SELECT count(*) FROM key_n WHERE n_keys = 1) AS unambiguous_claims,
  (SELECT count(*) FROM key_n WHERE n_keys > 1) AS ambiguous_claims,
  count(*) FILTER (WHERE is_alt AND is_ref) AS ref_equals_alt,
  count(*) FILTER (WHERE is_alt AND NOT is_ref) AS alt_claims,
  count(*) FILTER (WHERE is_ref AND NOT is_alt) AS ref_claims,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000) AS alt_measured,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5) AS alt_af_gt_0_5,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.9) AS alt_af_gt_0_9,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF >= 0.05 AND AF <= 0.5) AS alt_af_0_05_to_0_5,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF < 0.05) AS alt_af_lt_0_05,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5 AND median_raf > 0.5) AS both_say_common,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5 AND median_raf <= 0.5) AS raf_not_common,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5 AND median_raf IS NULL) AS raf_missing,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5 AND max_raf > 1) AS claims_with_a_raf_above_1,
  sum(n_assoc) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN >= 1000 AND AF > 0.5) AS assoc_rows_behind_af_gt_0_5,
  count(*) FILTER (WHERE is_alt AND NOT is_ref AND in_genomes AND AF IS NOT NULL AND AN < 1000) AS alt_an_below_1000
FROM j;

-- ---------------------------------------------------------------------------
-- 3. High missense-z genes and MaveDB scores, one score set at a time
-- ---------------------------------------------------------------------------

SELECT
  count(*) AS protein_coding_rows,
  count(*) FILTER (WHERE "mis.z_score" IS NOT NULL) AS z_present,
  count(*) FILTER (WHERE "mis.z_score" >= 3.09) AS z_ge_3_09,
  count(*) FILTER (WHERE "mis.z_score" >= 3.09 AND symbol_n_gencode = 1) AS z_ge_3_09_unique_symbol,
  count(*) FILTER (WHERE "mis.z_score" >= 3.09 AND symbol_n_gencode > 1) AS z_ge_3_09_shared_symbol,
  count(*) FILTER (WHERE "mis.z_score" >= 3.09 AND symbol_n_gencode = 1 AND mavedb_experiments_by_symbol > 0) AS z_with_mavedb_symbol_count
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart_protein_coding.parquet');

-- Mild is inside the synonymous 10th-90th percentile of the same score set.
-- A set is oriented only when the nonsense median falls outside that band
-- and the band is narrower than the gap between the synonymous and nonsense medians.
-- Column used: score. raw_score is not mixed in. Sets are not pooled.

WITH genes AS (
  SELECT symbol, ensg, "mis.z_score" AS z
  FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart_protein_coding.parquet')
  WHERE "mis.z_score" >= 3.09
    AND symbol_n_gencode = 1
),
raw AS (
  SELECT
    g.symbol,
    g.z,
    s.score_set_urn,
    s.urn,
    s.score,
    CASE
      WHEN s.protein_change IS NULL THEN 'missing'
      WHEN position('|' IN s.protein_change) > 0 THEN 'multi'
      WHEN regexp_matches(s.protein_change, 'p\.[A-Z][a-z]{2}[0-9]+Ter$') THEN 'nonsense'
      WHEN regexp_matches(s.protein_change, 'p\.[A-Z][a-z]{2}[0-9]+\*$') THEN 'nonsense'
      WHEN regexp_matches(s.protein_change, 'p\.[A-Z][a-z]{2}[0-9]+=$') THEN 'synonymous'
      WHEN regexp_matches(s.protein_change, 'p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}$')
           AND regexp_extract(s.protein_change, 'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})$', 1)
             = regexp_extract(s.protein_change, 'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})$', 2)
        THEN 'synonymous'
      WHEN regexp_matches(s.protein_change, 'p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}$')
           AND regexp_extract(s.protein_change, 'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})$', 2) <> 'Ter'
        THEN 'missense'
      ELSE 'other'
    END AS hgvs_class
  FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet') s
  JOIN genes g ON s.target_gene = g.symbol
),
urns AS (
  SELECT
    symbol,
    any_value(z) AS z,
    score_set_urn,
    urn,
    any_value(hgvs_class) AS hgvs_class,
    any_value(score) AS score
  FROM raw
  GROUP BY symbol, score_set_urn, urn
  HAVING count(DISTINCT hgvs_class) = 1
     AND count(DISTINCT score) FILTER (WHERE score IS NOT NULL) = 1
     AND count(*) FILTER (WHERE score IS NULL) = 0
     AND any_value(hgvs_class) IN ('missense', 'synonymous', 'nonsense')
),
sets AS (
  SELECT
    symbol,
    any_value(z) AS z,
    score_set_urn,
    count(*) FILTER (WHERE hgvs_class = 'missense') AS n_missense,
    count(*) FILTER (WHERE hgvs_class = 'synonymous') AS n_syn,
    count(*) FILTER (WHERE hgvs_class = 'nonsense') AS n_nonsense,
    quantile_cont(score, 0.1) FILTER (WHERE hgvs_class = 'synonymous') AS syn_p10,
    quantile_cont(score, 0.9) FILTER (WHERE hgvs_class = 'synonymous') AS syn_p90,
    median(score) FILTER (WHERE hgvs_class = 'synonymous') AS syn_med,
    median(score) FILTER (WHERE hgvs_class = 'nonsense') AS non_med,
    median(score) FILTER (WHERE hgvs_class = 'missense') AS mis_med
  FROM urns
  GROUP BY symbol, score_set_urn
),
lim AS (
  SELECT *
  FROM sets
  WHERE n_syn >= 10
    AND n_nonsense >= 10
    AND n_missense >= 20
),
scored AS (
  SELECT
    l.*,
    (l.non_med < l.syn_p10 OR l.non_med > l.syn_p90) AS separated,
    ((l.syn_p90 - l.syn_p10) < abs(l.syn_med - l.non_med)) AS width_lt_gap,
    count(*) FILTER (
      WHERE u.hgvs_class = 'missense'
        AND u.score BETWEEN l.syn_p10 AND l.syn_p90
    ) AS n_in_syn_band
  FROM lim l
  JOIN urns u
    ON l.symbol = u.symbol
   AND l.score_set_urn = u.score_set_urn
  GROUP BY ALL
)
SELECT
  s.symbol,
  s.z,
  s.score_set_urn,
  e.score_set_title,
  s.n_missense,
  s.n_syn,
  s.n_nonsense,
  s.syn_med,
  s.non_med,
  s.mis_med,
  s.syn_p10,
  s.syn_p90,
  s.separated,
  s.width_lt_gap,
  (s.separated AND s.width_lt_gap) AS oriented,
  s.n_in_syn_band,
  s.n_in_syn_band * 1.0 / s.n_missense AS frac_in_syn_band
FROM scored s
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
  ON s.score_set_urn = e.score_set_urn
ORDER BY s.symbol, s.score_set_urn;

-- ---------------------------------------------------------------------------
-- 4. Conflicting ClinVar labels with high REVEL or high SpliceAI
-- ---------------------------------------------------------------------------

-- High REVEL: revel_max >= 0.75. Also counted at 0.5.
-- High SpliceAI: spliceai_ds_max >= 0.5. Also counted at 0.2 and 0.8.
-- Exome score columns are counted only when the genome score is null.
-- They are not copied into revel_max or spliceai_ds_max.

SELECT
  count(*) AS conflict_keyed,
  count(*) FILTER (WHERE t.variant_key IS NULL) AS missing_from_tier_a,
  count(*) FILTER (WHERE t.in_genomes) AS in_genomes,
  count(*) FILTER (WHERE t.revel_max IS NOT NULL) AS revel_max_present,
  count(*) FILTER (WHERE t.spliceai_ds_max IS NOT NULL) AS spliceai_ds_max_present,
  count(*) FILTER (WHERE t.revel_max IS NOT NULL OR t.spliceai_ds_max IS NOT NULL) AS either_genome_score,
  count(*) FILTER (WHERE t.revel_max >= 0.5) AS revel_ge_0_5,
  count(*) FILTER (WHERE t.revel_max >= 0.75) AS revel_ge_0_75,
  count(*) FILTER (WHERE t.spliceai_ds_max >= 0.2) AS spliceai_ge_0_2,
  count(*) FILTER (WHERE t.spliceai_ds_max >= 0.5) AS spliceai_ge_0_5,
  count(*) FILTER (WHERE t.spliceai_ds_max >= 0.8) AS spliceai_ge_0_8,
  count(*) FILTER (WHERE t.revel_max >= 0.75 OR t.spliceai_ds_max >= 0.5) AS either_high,
  count(*) FILTER (WHERE t.revel_max >= 0.75 AND t.spliceai_ds_max >= 0.5) AS both_high,
  count(*) FILTER (
    WHERE t.revel_max IS NOT NULL AND t.revel_max_exome IS NOT NULL
      AND abs(t.revel_max - t.revel_max_exome) > 0.05
  ) AS revel_datasets_differ_gt_0_05,
  count(*) FILTER (
    WHERE t.spliceai_ds_max IS NOT NULL AND t.spliceai_ds_max_exome IS NOT NULL
      AND abs(t.spliceai_ds_max - t.spliceai_ds_max_exome) > 0.05
  ) AS spliceai_datasets_differ_gt_0_05,
  count(*) FILTER (WHERE t.revel_max IS NULL AND t.revel_max_exome >= 0.75) AS exome_only_revel_ge_0_75,
  count(*) FILTER (WHERE t.spliceai_ds_max IS NULL AND t.spliceai_ds_max_exome >= 0.5) AS exome_only_spliceai_ge_0_5
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') g
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet') t
  ON g.variant_key = t.variant_key
WHERE g.variant_key IS NOT NULL
  AND g.ClinicalSignificance ILIKE '%conflict%';

-- ---------------------------------------------------------------------------
-- Spot checks. These re-read stored columns. They do not impute a frequency.
-- ---------------------------------------------------------------------------

SELECT variant_key, ref, alt, AF, AC, AN, AC * 1.0 / AN AS ac_over_an,
       filter, in_exomes, AF_exome, revel_max, revel_max_exome,
       spliceai_ds_max, spliceai_ds_max_exome
FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet')
WHERE variant_key IN (
  'chr6:26092913:G:A',
  'chr11:5227002:T:A',
  'chr10:48066:C:A',
  'chr17:50600964:A:T',
  'chr14:54902784:A:G',
  'chr10:13288508:T:C'
)
ORDER BY variant_key;

SELECT variant_key, GeneSymbol, ClinicalSignificance, review_stars,
       germline_pathogenic, rsid, ref, alt
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet')
WHERE variant_key IN (
  'chr6:26092913:G:A',
  'chr11:5227002:T:A',
  'chr10:48066:C:A',
  'chr17:50600964:A:T',
  'chr14:54902784:A:G',
  'chr10:13288508:T:C'
)
ORDER BY variant_key;

SELECT rsid, chrom, pos, risk_allele, "STRONGEST SNP-RISK ALLELE",
       "RISK ALLELE FREQUENCY", variant_key, "DISEASE/TRAIT"
FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')
WHERE rsid IN ('rs198530', 'rs1753589')
ORDER BY rsid;

SELECT score_set_urn, urn, target_gene, protein_change, score
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
WHERE urn IN (
  'urn:mavedb:00000013-a-1#1821',
  'urn:mavedb:00000013-a-1#1880',
  'urn:mavedb:00000013-a-1#1071',
  'urn:mavedb:00000013-a-1#1185',
  'urn:mavedb:00000013-a-1#1756',
  'urn:mavedb:00000035-a-1#1859',
  'urn:mavedb:00000035-a-1#210'
)
ORDER BY score_set_urn, urn;
