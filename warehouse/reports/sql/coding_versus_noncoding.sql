-- Coding versus noncoding shares of clumped GWAS loci.
-- Run from the repository root with DuckDB 2.x.
-- Class comes from the lead association CONTEXT and INTERGENIC.
-- Mapped genes and interval overlap are not used.

CREATE OR REPLACE TEMP MACRO share6(n, d) AS (
  CASE
    WHEN d IS NULL OR d = 0 OR n IS NULL THEN NULL
    ELSE printf(
      '%d.%06d',
      CASE
        WHEN ((n * 1000000) % d) * 2 >= d THEN (n * 1000000) // d + 1
        ELSE (n * 1000000) // d
      END // 1000000,
      CASE
        WHEN ((n * 1000000) % d) * 2 >= d THEN (n * 1000000) // d + 1
        ELSE (n * 1000000) // d
      END % 1000000
    )
  END
);

CREATE OR REPLACE TEMP TABLE lead_class AS
WITH assoc AS (
  SELECT
    chrom,
    pos,
    rsid,
    p_value_float,
    CONTEXT,
    INTERGENIC,
    SNPS,
    CASE
      WHEN "MAPPED_TRAIT_URI" IS NULL OR trim("MAPPED_TRAIT_URI") = '' THEN [CAST(NULL AS VARCHAR)]
      ELSE
        list_filter(
          list_transform(string_split("MAPPED_TRAIT_URI", ','), lambda x: trim(x)),
          lambda x: x <> ''
        )
    END AS uris
  FROM read_parquet('warehouse/gwas/associations.parquet')
),
assoc_uri AS (
  SELECT
    a.chrom,
    a.pos,
    a.rsid,
    a.p_value_float,
    a.CONTEXT,
    a.INTERGENIC,
    a.SNPS,
    u.uri AS trait_uri
  FROM assoc a,
  LATERAL unnest(
    CASE WHEN len(a.uris) = 0 THEN [CAST(NULL AS VARCHAR)] ELSE a.uris END
  ) AS u(uri)
),
matched AS (
  SELECT
    l.trait_uri,
    l.chrom,
    l.lead_pos,
    l.lead_rsid,
    l.lead_p,
    l.n_associations,
    a.CONTEXT,
    a.INTERGENIC,
    a.SNPS
  FROM read_parquet('warehouse/gwas/loci.parquet') l
  JOIN assoc_uri a
    ON a.chrom IS NOT DISTINCT FROM l.chrom
   AND a.pos IS NOT DISTINCT FROM l.lead_pos
   AND a.rsid IS NOT DISTINCT FROM l.lead_rsid
   AND a.p_value_float IS NOT DISTINCT FROM l.lead_p
   AND a.trait_uri IS NOT DISTINCT FROM l.trait_uri
),
one AS (
  SELECT
    trait_uri,
    chrom,
    lead_pos,
    lead_rsid,
    lead_p,
    any_value(n_associations) AS n_associations,
    count(*) AS n_match,
    count(DISTINCT CONTEXT) AS n_ctx,
    count(DISTINCT INTERGENIC) AS n_flag,
    min(CONTEXT) AS CONTEXT,
    min(INTERGENIC) AS INTERGENIC,
    min(SNPS) AS SNPS
  FROM matched
  GROUP BY 1, 2, 3, 4, 5
),
tok AS (
  SELECT
    *,
    list_filter(
      list_transform(
        string_split(replace(coalesce(CONTEXT, ''), ' x ', ';'), ';'),
        lambda x: trim(x)
      ),
      lambda x: x <> ''
    ) AS toks
  FROM one
)
SELECT
  trait_uri,
  chrom,
  lead_pos,
  lead_rsid,
  lead_p,
  n_associations,
  n_match,
  n_ctx,
  n_flag,
  CONTEXT,
  INTERGENIC,
  SNPS,
  toks,
  len(toks) AS n_tok,
  list_contains(toks, 'missense_variant') AS has_missense,
  list_contains(toks, 'intron_variant') AS has_intron,
  list_contains(toks, '3_prime_UTR_variant') AS has_utr3,
  list_contains(toks, '5_prime_UTR_variant') AS has_utr5,
  (list_contains(toks, '3_prime_UTR_variant') OR list_contains(toks, '5_prime_UTR_variant')) AS has_utr,
  list_contains(toks, 'intergenic_variant') AS has_ig,
  len(list_filter(toks, lambda x: x NOT IN (
    'missense_variant',
    'intron_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'intergenic_variant'
  ))) AS n_other,
  CASE
    WHEN trim(coalesce(INTERGENIC, '')) = '1' THEN '1'
    WHEN trim(coalesce(INTERGENIC, '')) = '0' THEN '0'
    ELSE 'blank'
  END AS flag
FROM tok;

CREATE OR REPLACE TEMP TABLE classed AS
SELECT
  *,
  (has_missense::INTEGER + has_intron::INTEGER + has_utr::INTEGER + has_ig::INTEGER) AS n_requested,
  CASE
    WHEN n_tok = 0 THEN 'no_context'
    WHEN n_other = 0 AND has_missense AND NOT has_intron AND NOT has_utr AND NOT has_ig THEN 'missense'
    WHEN n_other = 0 AND has_intron AND NOT has_missense AND NOT has_utr AND NOT has_ig THEN 'intronic'
    WHEN n_other = 0 AND has_utr AND NOT has_missense AND NOT has_intron AND NOT has_ig THEN 'utr'
    WHEN n_other = 0 AND has_ig AND NOT has_missense AND NOT has_intron AND NOT has_utr THEN 'intergenic_context'
    WHEN n_other = 0 AND n_requested > 1 THEN 'mixed_requested'
    WHEN n_other > 0 AND n_requested = 0 THEN 'other'
    ELSE 'mixed_with_other'
  END AS exclusive_class
FROM lead_class;

-- Join and class checks. A non-zero disagreement or an unmatched locus means the class claim stops.
SELECT
  'qc' AS section,
  metric,
  n
FROM (
  SELECT 'loci' AS metric, count(*) AS n FROM read_parquet('warehouse/gwas/loci.parquet')
  UNION ALL
  SELECT 'leads_recovered', count(*) FROM classed
  UNION ALL
  SELECT 'context_disagreement', count(*) FROM classed WHERE n_ctx <> 1
  UNION ALL
  SELECT 'intergenic_flag_disagreement', count(*) FROM classed WHERE n_flag <> 1
  UNION ALL
  SELECT 'match_exceeds_members', count(*) FROM classed WHERE n_match > n_associations
  UNION ALL
  SELECT 'loci_with_duplicate_lead_rows', count(*) FROM classed WHERE n_match > 1
  UNION ALL
  SELECT 'max_lead_row_matches', max(n_match) FROM classed
  UNION ALL
  SELECT 'multi_token_context', count(*) FROM classed WHERE n_tok > 1
  UNION ALL
  SELECT 'trait_name_miss_nonnull_uri', count(*)
  FROM classed c
  LEFT JOIN read_parquet('warehouse/gwas/traits.parquet') t
    ON t.trait_uri IS NOT DISTINCT FROM c.trait_uri
  WHERE c.trait_uri IS NOT NULL AND t.trait_uri IS NULL
  UNION ALL
  SELECT 'empty_trait_uri_loci', count(*) FROM classed WHERE trait_uri IS NULL
  UNION ALL
  SELECT 'trait_groups_including_empty_uri', count(*) FROM (SELECT trait_uri FROM classed GROUP BY 1)
  UNION ALL
  SELECT 'distinct_lead_positions', count(*) FROM (SELECT DISTINCT chrom, lead_pos FROM classed)
  UNION ALL
  SELECT 'null_lead_rsid', count(*) FROM classed WHERE lead_rsid IS NULL
  UNION ALL
  SELECT 'blank_context', count(*) FROM classed WHERE n_tok = 0
  UNION ALL
  SELECT 'blank_context_esv_or_nsv', count(*)
  FROM classed
  WHERE n_tok = 0 AND (SNPS ILIKE 'esv%' OR SNPS ILIKE 'nsv%')
  UNION ALL
  SELECT 'truth_lead_is_intron_flag0', count(*)
  FROM classed
  WHERE trait_uri = 'http://www.ebi.ac.uk/efo/EFO_0004348'
    AND chrom = 'chr22'
    AND lead_rsid = 'rs9330814'
    AND lead_pos = 45968311
    AND CONTEXT = 'intron_variant'
    AND INTERGENIC = '0'
    AND n_associations = 19
);

SELECT
  'exclusive' AS section,
  exclusive_class,
  count(*) AS n,
  share6(count(*), (SELECT count(*) FROM classed)) AS share_of_loci,
  share6(
    count(*),
    (SELECT count(*) FROM classed WHERE n_tok > 0)
  ) AS share_of_loci_with_context
FROM classed
GROUP BY exclusive_class
ORDER BY n DESC;

SELECT
  'contains' AS section,
  label,
  n,
  share6(n, (SELECT count(*) FROM classed)) AS share_of_loci,
  share6(n, (SELECT count(*) FROM classed WHERE n_tok > 0)) AS share_of_loci_with_context
FROM (
  SELECT 'missense_variant' AS label, count(*) FILTER (WHERE has_missense) AS n FROM classed
  UNION ALL
  SELECT 'intron_variant', count(*) FILTER (WHERE has_intron) FROM classed
  UNION ALL
  SELECT 'utr', count(*) FILTER (WHERE has_utr) FROM classed
  UNION ALL
  SELECT '3_prime_UTR_variant', count(*) FILTER (WHERE has_utr3) FROM classed
  UNION ALL
  SELECT '5_prime_UTR_variant', count(*) FILTER (WHERE has_utr5) FROM classed
  UNION ALL
  SELECT 'utr_both_ends', count(*) FILTER (WHERE has_utr3 AND has_utr5) FROM classed
  UNION ALL
  SELECT 'intergenic_variant', count(*) FILTER (WHERE has_ig) FROM classed
  UNION ALL
  SELECT 'intronic_or_utr_or_intergenic_context', count(*) FILTER (WHERE has_intron OR has_utr OR has_ig) FROM classed
  UNION ALL
  SELECT 'requested_four_any', count(*) FILTER (WHERE has_missense OR has_intron OR has_utr OR has_ig) FROM classed
  UNION ALL
  SELECT 'context_present', count(*) FILTER (WHERE n_tok > 0) FROM classed
  UNION ALL
  SELECT 'context_absent', count(*) FILTER (WHERE n_tok = 0) FROM classed
);

SELECT
  'intergenic_flag' AS section,
  flag,
  count(*) AS n,
  share6(count(*), (SELECT count(*) FROM classed)) AS share_of_loci,
  share6(count(*), (SELECT count(*) FROM classed WHERE n_tok > 0)) AS share_of_loci_with_context
FROM classed
GROUP BY flag
ORDER BY flag;

SELECT
  'flag_by_class' AS section,
  exclusive_class,
  count(*) FILTER (WHERE flag = '0') AS flag_0,
  count(*) FILTER (WHERE flag = '1') AS flag_1,
  count(*) FILTER (WHERE flag = 'blank') AS flag_blank,
  count(*) AS n
FROM classed
GROUP BY exclusive_class
ORDER BY n DESC;

SELECT
  'other_term' AS section,
  u.tok AS term,
  count(*) AS n,
  share6(count(*), (SELECT count(*) FROM classed)) AS share_of_loci,
  share6(count(*), (SELECT count(*) FROM classed WHERE n_tok > 0)) AS share_of_loci_with_context
FROM classed, LATERAL unnest(toks) AS u(tok)
WHERE exclusive_class = 'other'
GROUP BY u.tok
ORDER BY n DESC, term;

SELECT
  'other_group' AS section,
  grp,
  count(*) AS n,
  share6(count(*), (SELECT count(*) FROM classed)) AS share_of_loci,
  share6(count(*), (SELECT count(*) FROM classed WHERE n_tok > 0)) AS share_of_loci_with_context
FROM (
  SELECT
    CASE CONTEXT
      WHEN 'stop_gained' THEN 'other_protein_sequence'
      WHEN 'frameshift_variant' THEN 'other_protein_sequence'
      WHEN 'inframe_insertion' THEN 'other_protein_sequence'
      WHEN 'inframe_deletion' THEN 'other_protein_sequence'
      WHEN 'start_lost' THEN 'other_protein_sequence'
      WHEN 'stop_lost' THEN 'other_protein_sequence'
      WHEN 'stop_retained_variant' THEN 'other_protein_sequence'
      WHEN 'protein_altering_variant' THEN 'other_protein_sequence'
      WHEN 'synonymous_variant' THEN 'synonymous'
      WHEN 'splice_region_variant' THEN 'splice'
      WHEN 'splice_polypyrimidine_tract_variant' THEN 'splice'
      WHEN 'splice_donor_region_variant' THEN 'splice'
      WHEN 'splice_donor_variant' THEN 'splice'
      WHEN 'splice_donor_5th_base_variant' THEN 'splice'
      WHEN 'splice_acceptor_variant' THEN 'splice'
      WHEN 'regulatory_region_variant' THEN 'regulatory_or_tf'
      WHEN 'TF_binding_site_variant' THEN 'regulatory_or_tf'
      WHEN 'TFBS_ablation' THEN 'regulatory_or_tf'
      WHEN 'non_coding_transcript_exon_variant' THEN 'noncoding_transcript_exon'
      WHEN 'mature_miRNA_variant' THEN 'mature_mirna'
      ELSE 'unlisted_other'
    END AS grp
  FROM classed
  WHERE exclusive_class = 'other'
)
GROUP BY grp
ORDER BY n DESC;

SELECT
  'largest_trait' AS section,
  t.mapped_trait,
  c.trait_uri,
  count(*) AS n_loci,
  count(*) FILTER (WHERE c.n_tok > 0) AS n_with_context,
  count(*) FILTER (WHERE c.exclusive_class = 'missense') AS n_missense,
  share6(count(*) FILTER (WHERE c.exclusive_class = 'missense'), count(*)) AS share_missense,
  count(*) FILTER (WHERE c.exclusive_class = 'intronic') AS n_intronic,
  share6(count(*) FILTER (WHERE c.exclusive_class = 'intronic'), count(*)) AS share_intronic,
  count(*) FILTER (WHERE c.exclusive_class = 'utr') AS n_utr,
  share6(count(*) FILTER (WHERE c.exclusive_class = 'utr'), count(*)) AS share_utr,
  count(*) FILTER (WHERE c.has_utr3) AS n_utr3,
  count(*) FILTER (WHERE c.has_utr5) AS n_utr5,
  count(*) FILTER (WHERE c.exclusive_class = 'intergenic_context') AS n_intergenic_context,
  share6(count(*) FILTER (WHERE c.exclusive_class = 'intergenic_context'), count(*)) AS share_intergenic_context,
  count(*) FILTER (WHERE c.exclusive_class = 'other') AS n_other,
  share6(count(*) FILTER (WHERE c.exclusive_class = 'other'), count(*)) AS share_other,
  count(*) FILTER (WHERE c.flag = '1') AS n_flag_1,
  share6(count(*) FILTER (WHERE c.flag = '1'), count(*)) AS share_flag_1
FROM classed c
JOIN read_parquet('warehouse/gwas/traits.parquet') t
  ON t.trait_uri IS NOT DISTINCT FROM c.trait_uri
GROUP BY t.mapped_trait, c.trait_uri
HAVING count(*) >= 1000
ORDER BY n_loci DESC, lower(t.mapped_trait), c.trait_uri;

SELECT
  'largest_trait_pool' AS section,
  count(*) AS n_traits,
  sum(n_loci) AS n_loci,
  sum(n_with_context) AS n_with_context,
  sum(n_missense) AS n_missense,
  share6(sum(n_missense), sum(n_loci)) AS share_missense,
  sum(n_intronic) AS n_intronic,
  share6(sum(n_intronic), sum(n_loci)) AS share_intronic,
  sum(n_utr) AS n_utr,
  share6(sum(n_utr), sum(n_loci)) AS share_utr,
  sum(n_intergenic_context) AS n_intergenic_context,
  share6(sum(n_intergenic_context), sum(n_loci)) AS share_intergenic_context,
  sum(n_other) AS n_other,
  share6(sum(n_other), sum(n_loci)) AS share_other,
  sum(n_flag_1) AS n_flag_1,
  share6(sum(n_flag_1), sum(n_loci)) AS share_flag_1,
  share6(sum(n_loci), (SELECT count(*) FROM classed)) AS share_of_all_loci
FROM (
  SELECT
    count(*) AS n_loci,
    count(*) FILTER (WHERE n_tok > 0) AS n_with_context,
    count(*) FILTER (WHERE exclusive_class = 'missense') AS n_missense,
    count(*) FILTER (WHERE exclusive_class = 'intronic') AS n_intronic,
    count(*) FILTER (WHERE exclusive_class = 'utr') AS n_utr,
    count(*) FILTER (WHERE exclusive_class = 'intergenic_context') AS n_intergenic_context,
    count(*) FILTER (WHERE exclusive_class = 'other') AS n_other,
    count(*) FILTER (WHERE flag = '1') AS n_flag_1
  FROM classed
  GROUP BY trait_uri
  HAVING count(*) >= 1000
);
