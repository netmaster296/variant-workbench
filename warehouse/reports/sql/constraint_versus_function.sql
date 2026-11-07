-- LOEUF decile versus ClinVar pathogenic fraction, and MaveDB assay coverage.
-- Read-only. Does not write a mart.
--
-- Denominator: protein-coding genes in gene_mart_protein_coding.parquet.
-- LOEUF decile is loeuf_decile (gnomAD lof.oe_ci.upper_bin_decile). Decile 0 is
-- the most constrained bin. constraint_flags is null only when the gene has no
-- constraint row; a matched row stores [] when it has no flag.
-- Pathogenic at two or more stars is the mart's germline_pathogenic bucket
-- (not likely pathogenic alone) in star columns 2, 3, and 4.
-- The allele fraction is sum(those alleles) / sum(clinvar_n) inside the decile.
-- clinvar_n is every matched ClinVar allele, not only the pathogenic ones.
--
-- MaveDB effect sizes are not pooled. Score sets do not share a scale.
-- The span join is diagnosed at the bottom and is not used as an assay join.
-- Exact assay: experiments.target_name equals the protein-coding symbol,
-- organism_name is Homo sapiens or null, and mapped_variants_scored.score
-- is numeric. One score per (score_set_urn, urn). Gallus gallus and
-- Escherichia coli score sets are excluded.

CREATE TEMP TABLE coding AS
SELECT
  ensg,
  symbol,
  loeuf,
  loeuf_decile,
  constraint_flags,
  symbol_n_gencode,
  clinvar_n,
  mavedb_experiments_by_symbol,
  mavedb_variant_keys_in_span,
  (
    clinvar_germline_pathogenic_star_2
    + clinvar_germline_pathogenic_star_3
    + clinvar_germline_pathogenic_star_4
  ) AS p_ge2,
  (
    clinvar_likely_pathogenic_star_2
    + clinvar_likely_pathogenic_star_3
    + clinvar_likely_pathogenic_star_4
  ) AS lp_ge2
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart_protein_coding.parquet');

-- Denominator.
SELECT
  count(*) AS protein_coding_genes,
  count(*) FILTER (WHERE constraint_flags IS NULL) AS no_constraint_row,
  count(*) FILTER (WHERE constraint_flags IS NOT NULL AND loeuf IS NULL) AS row_null_loeuf,
  count(*) FILTER (WHERE loeuf IS NOT NULL AND loeuf_decile IS NULL) AS loeuf_without_decile,
  count(*) FILTER (WHERE constraint_flags IS NOT NULL AND loeuf_decile IS NULL) AS row_null_decile,
  count(*) FILTER (WHERE loeuf_decile IS NULL) AS null_decile,
  count(*) FILTER (WHERE loeuf_decile IS NOT NULL) AS with_decile,
  count(*) FILTER (WHERE clinvar_n > 0) AS with_clinvar,
  count(*) FILTER (WHERE p_ge2 > 0) AS with_p_ge2,
  count(*) FILTER (WHERE mavedb_experiments_by_symbol > 0) AS mart_symbol_assay,
  count(*) FILTER (WHERE mavedb_variant_keys_in_span > 0) AS mart_span_key
FROM coding;

-- Same symbol in two deciles. Both have p_ge2 = 0, so they do not move numerators.
SELECT
  symbol,
  count(*) AS gene_rows,
  string_agg(ensg || ' decile ' || decile_text, ', ' ORDER BY ensg) AS copies,
  any_value(clinvar_n) AS clinvar_n,
  any_value(p_ge2) AS p_ge2
FROM (
  SELECT
    symbol,
    ensg,
    clinvar_n,
    p_ge2,
    cast(loeuf_decile AS VARCHAR) AS decile_text
  FROM coding
  WHERE symbol IN (
    SELECT symbol
    FROM coding
    WHERE loeuf_decile IS NOT NULL
    GROUP BY symbol
    HAVING count(DISTINCT loeuf_decile) > 1
  )
)
GROUP BY symbol
ORDER BY symbol;

-- Allele-weighted pathogenic fraction by LOEUF decile.
-- Within a decile each symbol occurs once. Median of the gene fraction is 0
-- in every decile; mean_frac is the unweighted mean among genes with clinvar_n > 0.
SELECT
  loeuf_decile,
  count(*) AS genes,
  count(*) FILTER (WHERE clinvar_n > 0) AS genes_with_clinvar,
  sum(clinvar_n) AS alleles,
  sum(p_ge2) AS p_ge2,
  sum(lp_ge2) AS lp_ge2,
  round(sum(p_ge2)::DOUBLE / nullif(sum(clinvar_n), 0), 6) AS frac_p_ge2,
  round((sum(p_ge2) + sum(lp_ge2))::DOUBLE / nullif(sum(clinvar_n), 0), 6) AS frac_p_or_lp_ge2,
  round(avg(p_ge2::DOUBLE / clinvar_n) FILTER (WHERE clinvar_n > 0), 6) AS mean_gene_frac,
  round(median(p_ge2::DOUBLE / clinvar_n) FILTER (WHERE clinvar_n > 0), 6) AS median_gene_frac,
  count(*) FILTER (WHERE p_ge2 > 0) AS genes_with_p_ge2,
  round(median(loeuf), 3) AS median_loeuf,
  min(loeuf) AS min_loeuf,
  max(loeuf) AS max_loeuf
FROM coding
WHERE loeuf_decile IS NOT NULL
GROUP BY loeuf_decile
ORDER BY loeuf_decile;

-- Spearman, average ranks for ties. y is the gene-level fraction.
WITH base AS (
  SELECT ensg, loeuf_decile::DOUBLE AS x, (p_ge2::DOUBLE / clinvar_n) AS y
  FROM coding
  WHERE loeuf_decile IS NOT NULL AND clinvar_n > 0
),
rx AS (
  SELECT ensg, x, y, avg(rn) OVER (PARTITION BY x) AS rank_x
  FROM (
    SELECT ensg, x, y, row_number() OVER (ORDER BY x, ensg) AS rn
    FROM base
  )
),
ranks AS (
  SELECT rank_x, avg(rn) OVER (PARTITION BY y) AS rank_y
  FROM (
    SELECT ensg, x, y, rank_x, row_number() OVER (ORDER BY y, ensg) AS rn
    FROM rx
  )
)
SELECT 'frac_clinvar_gt0' AS contrast, round(corr(rank_x, rank_y), 4) AS spearman, count(*) AS n
FROM ranks;

WITH base AS (
  SELECT ensg, loeuf_decile::DOUBLE AS x, (p_ge2::DOUBLE / clinvar_n) AS y
  FROM coding
  WHERE loeuf_decile IS NOT NULL AND clinvar_n >= 10
),
rx AS (
  SELECT ensg, x, y, avg(rn) OVER (PARTITION BY x) AS rank_x
  FROM (
    SELECT ensg, x, y, row_number() OVER (ORDER BY x, ensg) AS rn
    FROM base
  )
),
ranks AS (
  SELECT rank_x, avg(rn) OVER (PARTITION BY y) AS rank_y
  FROM (
    SELECT ensg, x, y, rank_x, row_number() OVER (ORDER BY y, ensg) AS rn
    FROM rx
  )
)
SELECT 'frac_clinvar_ge10' AS contrast, round(corr(rank_x, rank_y), 4) AS spearman, count(*) AS n
FROM ranks;

WITH base AS (
  SELECT
    ensg,
    loeuf_decile::DOUBLE AS x,
    CASE WHEN p_ge2 > 0 THEN 1.0 ELSE 0.0 END AS y
  FROM coding
  WHERE loeuf_decile IS NOT NULL
),
rx AS (
  SELECT ensg, x, y, avg(rn) OVER (PARTITION BY x) AS rank_x
  FROM (
    SELECT ensg, x, y, row_number() OVER (ORDER BY x, ensg) AS rn
    FROM base
  )
),
ranks AS (
  SELECT rank_x, avg(rn) OVER (PARTITION BY y) AS rank_y
  FROM (
    SELECT ensg, x, y, rank_x, row_number() OVER (ORDER BY y, ensg) AS rn
    FROM rx
  )
)
SELECT 'any_p_ge2' AS contrast, round(corr(rank_x, rank_y), 4) AS spearman, count(*) AS n
FROM ranks;

-- Largest pathogenic numerators. These sit in the middle deciles as well as decile 0.
SELECT
  symbol,
  loeuf_decile,
  round(loeuf, 3) AS loeuf,
  clinvar_n,
  p_ge2,
  round(p_ge2::DOUBLE / clinvar_n, 4) AS frac
FROM coding
WHERE clinvar_n > 0
ORDER BY p_ge2 DESC
LIMIT 10;

-- Exact target-name match, from gene_coverage and from experiments.
SELECT
  count(*) AS coverage_names,
  count(DISTINCT g.ensg) AS genes
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet') c
JOIN coding g ON c.target_gene = g.symbol;

SELECT
  count(*) AS names_without_own_target_name
FROM (
  SELECT c.target_gene
  FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet') c
  JOIN coding g ON c.target_gene = g.symbol
  ANTI JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
    ON e.target_name = c.target_gene
);

-- Sole-target score sets whose organism is human or unrecorded.
CREATE TEMP TABLE assay_sets AS
SELECT
  e.score_set_urn,
  e.target_name AS symbol,
  e.organism_name,
  g.ensg,
  g.loeuf_decile
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
JOIN coding g ON e.target_name = g.symbol
WHERE e.organism_name IS NULL OR e.organism_name = 'Homo sapiens';

SELECT
  count(*) AS score_sets,
  count(DISTINCT symbol) AS genes,
  count(*) FILTER (WHERE organism_name IS NULL) AS null_organism_sets,
  count(DISTINCT symbol) FILTER (WHERE organism_name IS NULL) AS genes_with_null_organism_set
FROM assay_sets;

-- Symbols whose only exact score sets are not human.
SELECT e.target_name AS symbol, e.organism_name, count(*) AS score_sets
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
JOIN coding g ON e.target_name = g.symbol
GROUP BY e.target_name, e.organism_name
HAVING e.target_name NOT IN (SELECT symbol FROM assay_sets)
ORDER BY e.target_name, e.organism_name;

-- Coverage names with no score set of their own. n_score_sets still counts
-- each of them on a multi-target score set, so the mart symbol flag is not
-- "this gene has its own assay".
SELECT c.target_gene
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet') c
JOIN coding g ON c.target_gene = g.symbol
ANTI JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
  ON e.target_name = c.target_gene
ORDER BY c.target_gene;

-- Exact symbol, organism not recorded, and no Homo sapiens or other-organism set.
SELECT e.target_name AS symbol
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
JOIN coding g ON e.target_name = g.symbol
GROUP BY e.target_name
HAVING count(*) FILTER (WHERE e.organism_name = 'Homo sapiens') = 0
   AND count(*) FILTER (WHERE e.organism_name IS NOT NULL AND e.organism_name <> 'Homo sapiens') = 0
   AND count(*) FILTER (WHERE e.organism_name IS NULL) > 0
ORDER BY e.target_name;

-- Numeric scores. One row per variant urn inside a score set.
CREATE TEMP TABLE assay_scores AS
SELECT
  a.score_set_urn,
  a.symbol,
  a.loeuf_decile,
  a.ensg,
  s.urn,
  any_value(s.score) AS score
FROM assay_sets a
JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet') s
  ON s.score_set_urn = a.score_set_urn
WHERE s.score IS NOT NULL
GROUP BY a.score_set_urn, a.symbol, a.loeuf_decile, a.ensg, s.urn;

SELECT
  count(*) AS scored_urns,
  count(DISTINCT score_set_urn) AS score_sets_with_score,
  count(DISTINCT symbol) AS genes_with_score
FROM assay_scores;

SELECT
  loeuf_decile,
  count(DISTINCT symbol) AS genes,
  count(DISTINCT score_set_urn) AS score_sets,
  count(*) AS scored_urns
FROM assay_scores
GROUP BY loeuf_decile
ORDER BY loeuf_decile NULLS LAST;

-- Per score set, then the spread of those medians. Not an effect size.
WITH per AS (
  SELECT
    score_set_urn,
    symbol,
    count(*) AS n,
    min(score) AS min_score,
    quantile_cont(score, 0.5) AS median_score,
    max(score) AS max_score
  FROM assay_scores
  GROUP BY score_set_urn, symbol
)
SELECT
  count(*) AS score_sets,
  count(DISTINCT symbol) AS genes,
  count(*) FILTER (WHERE n >= 20) AS sets_with_at_least_20,
  min(median_score) AS min_of_medians,
  quantile_cont(median_score, 0.25) AS p25_of_medians,
  quantile_cont(median_score, 0.5) AS p50_of_medians,
  quantile_cont(median_score, 0.75) AS p75_of_medians,
  max(median_score) AS max_of_medians,
  min(min_score) AS min_score,
  max(max_score) AS max_score,
  count(*) FILTER (WHERE median_score < 0) AS median_below_0,
  count(*) FILTER (WHERE median_score > 1) AS median_above_1,
  count(*) FILTER (WHERE min_score >= 0 AND max_score <= 1) AS entirely_inside_0_1,
  count(*) FILTER (WHERE min_score < 0 AND max_score > 0) AS crosses_0,
  count(*) FILTER (WHERE max_score - min_score > 10) AS width_above_10
FROM per;

WITH per AS (
  SELECT symbol, score_set_urn, quantile_cont(score, 0.5) AS median_score
  FROM assay_scores
  GROUP BY symbol, score_set_urn
)
SELECT
  count(*) FILTER (WHERE n_sets >= 2) AS genes_with_two_sets,
  count(*) FILTER (WHERE n_sets >= 2 AND max_median - min_median > 1) AS genes_median_gap_above_1
FROM (
  SELECT symbol, count(*) AS n_sets, min(median_score) AS min_median, max(median_score) AS max_median
  FROM per
  GROUP BY symbol
);

-- Same gene, two scales. KCNQ2 medians are not one effect size.
SELECT
  score_set_urn,
  count(*) AS n,
  min(score) AS min_score,
  quantile_cont(score, 0.5) AS median_score,
  max(score) AS max_score
FROM assay_scores
WHERE symbol = 'KCNQ2'
GROUP BY score_set_urn
ORDER BY median_score;

SELECT
  score_set_urn,
  symbol,
  loeuf_decile,
  count(*) AS n,
  min(score) AS min_score,
  quantile_cont(score, 0.5) AS median_score,
  max(score) AS max_score
FROM assay_scores
WHERE score_set_urn IN (
  'urn:mavedb:00001244-b-1',
  'urn:mavedb:00001231-a-1',
  'urn:mavedb:00000062-a-1'
)
GROUP BY score_set_urn, symbol, loeuf_decile
ORDER BY symbol;

-- Span join. A key inside a gene is not that gene's assay.
-- Stopped as an effect-size join. See warehouse/qc/issues_constraint_versus_function.md.
CREATE TEMP TABLE keyed AS
SELECT DISTINCT variant_key, target_gene, chrom, pos
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
WHERE variant_key IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL;

SELECT
  count(*) AS key_target_positions,
  count(DISTINCT variant_key) AS variant_keys,
  count(DISTINCT target_gene) AS target_names
FROM keyed;

CREATE TEMP TABLE span_hit AS
SELECT
  g.ensg,
  g.symbol,
  g.biotype,
  k.target_gene,
  k.variant_key,
  k.target_gene = g.symbol AS exact_name
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet') g
JOIN keyed k
  ON k.chrom = g.chrom
 AND k.pos BETWEEN g.start AND g."end";

SELECT biotype, count(DISTINCT ensg) AS genes
FROM span_hit
GROUP BY biotype
ORDER BY genes DESC;

SELECT
  count(DISTINCT variant_key) AS keys_in_a_coding_span,
  count(DISTINCT variant_key) FILTER (WHERE n_genes > 1) AS keys_in_multiple_coding_genes,
  max(n_genes) AS max_coding_genes
FROM (
  SELECT variant_key, count(DISTINCT ensg) AS n_genes
  FROM span_hit
  WHERE biotype = 'protein_coding'
  GROUP BY variant_key
);

SELECT
  count(*) AS coding_genes,
  count(*) FILTER (WHERE exact_keys > 0) AS genes_with_exact_name_key,
  count(*) FILTER (WHERE exact_keys = 0) AS genes_with_no_exact_name_key,
  count(*) FILTER (WHERE exact_keys > 0 AND exact_keys < keys) AS genes_mixed
FROM (
  SELECT
    ensg,
    count(DISTINCT variant_key) AS keys,
    count(DISTINCT variant_key) FILTER (WHERE exact_name) AS exact_keys
  FROM span_hit
  WHERE biotype = 'protein_coding'
  GROUP BY ensg
);

-- Coding spans whose keys never use the gene symbol as target_gene.
-- A target that merely contains the symbol, or is that gene's MANE transcript,
-- is still not an exact symbol match. Listed so the failed join can be checked.
SELECT
  s.symbol,
  s.ensg,
  count(DISTINCT s.variant_key) AS keys,
  string_agg(DISTINCT s.target_gene, ' | ') AS target_names
FROM span_hit s
WHERE s.biotype = 'protein_coding'
  AND s.ensg NOT IN (
    SELECT ensg FROM span_hit WHERE biotype = 'protein_coding' AND exact_name
  )
GROUP BY s.symbol, s.ensg
ORDER BY keys DESC;
