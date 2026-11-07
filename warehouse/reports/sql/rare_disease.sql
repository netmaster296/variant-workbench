-- Rare-disease ranking from the protein-coding gene mart.
-- Read-only. This file does not write a mart or fill a null LOEUF by symbol.
--
-- Filters
--   universe:
--     warehouse/marts/gene_mart_protein_coding.parquet
--     every row has biotype protein_coding
--   loss-of-function intolerance:
--     loeuf IS NOT NULL AND loeuf < 0.6
--     loeuf is gnomAD v4.1 lof.oe_ci.upper
--     null LOEUF stays null; there is no symbol join to gene_constraint
--   high-confidence germline pathogenic alleles:
--     plp_stars_ge_2 >= 1
--     plp_stars_ge_2 is germline pathogenic plus likely pathogenic
--     at review stars 2, 3, and 4 only
--     star 0, star 1, and star null are not in that count
--     Pathogenic/Likely pathogenic is already only in the pathogenic class
--   broad expression:
--     n_tissues_tpm_ge_1 >= 27
--     that column counts GTEx v8 tissues with median TPM >= 1
--     the median file has 54 tissues, so 27 is at least half
--     a null expression value fails
--     tau is not a filter
--   order of the high-confidence set:
--     loeuf ascending, plp_stars_ge_2 descending,
--     n_tissues_tpm_ge_1 descending, symbol, ensg
--   split-out, only single-submitter or conflicting review status:
--     same universe, LOEUF, and breadth filters
--     plp_stars_ge_2 = 0
--     plp_star_1 >= 1
--     plp_star_0 = 0 and plp_star_null = 0
--     star 1 is the review status "single submitter" or
--     "conflicting classifications"; the mart does not separate those phrases
--     order: loeuf ascending, plp_star_1 descending,
--     n_tissues_tpm_ge_1 descending, symbol, ensg

CREATE OR REPLACE TEMP VIEW coding AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart_protein_coding.parquet');

CREATE OR REPLACE TEMP VIEW constraint_gene AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/constraint/gene_constraint.parquet');

CREATE OR REPLACE TEMP VIEW expression_gene AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/gtex/gene_expression.parquet');

CREATE OR REPLACE TEMP VIEW scored AS
SELECT
  c.ensg,
  c.symbol,
  c.chrom,
  c.biotype,
  c.loeuf,
  c.loeuf_decile,
  c.constraint_flags,
  c.symbol_n_gencode,
  c.n_tissues_tpm_ge_1,
  c.tau,
  c.clinvar_germline_pathogenic_n,
  c.clinvar_likely_pathogenic_n,
  c.clinvar_conflicting_n,
  c.clinvar_germline_pathogenic_star_0,
  c.clinvar_germline_pathogenic_star_1,
  c.clinvar_germline_pathogenic_star_2,
  c.clinvar_germline_pathogenic_star_3,
  c.clinvar_germline_pathogenic_star_4,
  c.clinvar_germline_pathogenic_star_null,
  c.clinvar_likely_pathogenic_star_0,
  c.clinvar_likely_pathogenic_star_1,
  c.clinvar_likely_pathogenic_star_2,
  c.clinvar_likely_pathogenic_star_3,
  c.clinvar_likely_pathogenic_star_4,
  c.clinvar_likely_pathogenic_star_null,
  (
    c.clinvar_germline_pathogenic_star_2
    + c.clinvar_germline_pathogenic_star_3
    + c.clinvar_germline_pathogenic_star_4
    + c.clinvar_likely_pathogenic_star_2
    + c.clinvar_likely_pathogenic_star_3
    + c.clinvar_likely_pathogenic_star_4
  ) AS plp_stars_ge_2,
  (
    c.clinvar_germline_pathogenic_star_2
    + c.clinvar_germline_pathogenic_star_3
    + c.clinvar_germline_pathogenic_star_4
  ) AS pathogenic_stars_ge_2,
  (
    c.clinvar_likely_pathogenic_star_2
    + c.clinvar_likely_pathogenic_star_3
    + c.clinvar_likely_pathogenic_star_4
  ) AS likely_pathogenic_stars_ge_2,
  (
    c.clinvar_germline_pathogenic_star_1
    + c.clinvar_likely_pathogenic_star_1
  ) AS plp_star_1,
  (
    c.clinvar_germline_pathogenic_star_0
    + c.clinvar_likely_pathogenic_star_0
  ) AS plp_star_0,
  (
    c.clinvar_germline_pathogenic_star_null
    + c.clinvar_likely_pathogenic_star_null
  ) AS plp_star_null,
  (
    c.clinvar_germline_pathogenic_star_3
    + c.clinvar_likely_pathogenic_star_3
  ) AS plp_star_3,
  (
    c.clinvar_germline_pathogenic_star_4
    + c.clinvar_likely_pathogenic_star_4
  ) AS plp_star_4,
  (
    c.clinvar_germline_pathogenic_n
    + c.clinvar_likely_pathogenic_n
  ) AS plp_any_star,
  (
    c.clinvar_germline_pathogenic_star_0
    + c.clinvar_germline_pathogenic_star_1
    + c.clinvar_germline_pathogenic_star_2
    + c.clinvar_germline_pathogenic_star_3
    + c.clinvar_germline_pathogenic_star_4
    + c.clinvar_germline_pathogenic_star_null
  ) AS pathogenic_star_sum,
  (
    c.clinvar_likely_pathogenic_star_0
    + c.clinvar_likely_pathogenic_star_1
    + c.clinvar_likely_pathogenic_star_2
    + c.clinvar_likely_pathogenic_star_3
    + c.clinvar_likely_pathogenic_star_4
    + c.clinvar_likely_pathogenic_star_null
  ) AS likely_pathogenic_star_sum,
  (
    c.clinvar_conflicting_star_0
    + c.clinvar_conflicting_star_1
    + c.clinvar_conflicting_star_2
    + c.clinvar_conflicting_star_3
    + c.clinvar_conflicting_star_4
    + c.clinvar_conflicting_star_null
  ) AS conflicting_star_sum
FROM coding c;

-- query: audit_identity
SELECT
  count(*) AS protein_coding_genes,
  count(DISTINCT ensg) AS distinct_ensg,
  count(*) FILTER (WHERE biotype = 'protein_coding') AS biotype_protein_coding,
  count(*) FILTER (WHERE biotype IS DISTINCT FROM 'protein_coding') AS other_biotype,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS symbol_collision_genes,
  count(*) FILTER (
    WHERE clinvar_germline_pathogenic_n <> pathogenic_star_sum
       OR clinvar_likely_pathogenic_n <> likely_pathogenic_star_sum
       OR clinvar_conflicting_n <> conflicting_star_sum
  ) AS star_bucket_mismatches
FROM scored;

-- query: audit_constraint
SELECT
  count(*) AS protein_coding_genes,
  count(*) FILTER (WHERE k.bare_ensg IS NULL) AS no_constraint_row,
  count(*) FILTER (WHERE k.bare_ensg IS NULL AND s.chrom = 'chrX') AS no_constraint_row_chrx,
  count(*) FILTER (WHERE k.bare_ensg IS NULL AND s.chrom = 'chrY') AS no_constraint_row_chry,
  count(*) FILTER (WHERE k.bare_ensg IS NULL AND s.chrom = 'chrM') AS no_constraint_row_chrm,
  count(*) FILTER (
    WHERE k.bare_ensg IS NULL
      AND s.chrom NOT IN ('chrX', 'chrY', 'chrM')
  ) AS no_constraint_row_other_chrom,
  count(*) FILTER (WHERE k.bare_ensg IS NOT NULL AND s.loeuf IS NULL) AS constraint_row_loeuf_null,
  count(*) FILTER (WHERE s.loeuf IS NOT NULL) AS loeuf_not_null,
  count(*) FILTER (
    WHERE s.loeuf IS NOT NULL AND s.loeuf IS DISTINCT FROM k.loeuf
  ) AS loeuf_value_disagreements,
  count(*) FILTER (WHERE s.loeuf IS NOT NULL AND k.bare_ensg IS NULL) AS loeuf_without_constraint_row,
  (
    SELECT count(*)
    FROM (
      SELECT bare_ensg
      FROM constraint_gene
      WHERE bare_ensg IS NOT NULL
      GROUP BY bare_ensg
      HAVING count(*) > 1
    )
  ) AS duplicate_constraint_bare_ensg,
  (SELECT count(*) FROM constraint_gene WHERE chromosome = 'chrX') AS constraint_rows_chrx,
  (SELECT count(*) FROM constraint_gene WHERE chromosome = 'chrY') AS constraint_rows_chry,
  (SELECT count(*) FROM constraint_gene WHERE chromosome = 'chrM') AS constraint_rows_chrm
FROM scored s
LEFT JOIN constraint_gene k ON s.ensg = k.bare_ensg;

-- query: audit_expression
SELECT
  count(*) FILTER (WHERE e.ensg IS NULL) AS no_expression_row,
  count(*) FILTER (WHERE e.ensg IS NOT NULL) AS expression_row,
  count(*) FILTER (
    WHERE s.n_tissues_tpm_ge_1 IS DISTINCT FROM e.n_tissues_tpm_ge_1
  ) AS breadth_disagreements,
  count(*) FILTER (WHERE s.tau IS DISTINCT FROM e.tau) AS tau_disagreements,
  min(e.n_tissues) AS min_n_tissues,
  max(e.n_tissues) AS max_n_tissues,
  count(DISTINCT e.n_tissues) AS distinct_n_tissues
FROM scored s
LEFT JOIN expression_gene e ON s.ensg = e.ensg;

-- query: denominators
SELECT
  count(*) AS protein_coding_genes,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_not_null,
  count(*) FILTER (WHERE loeuf IS NULL) AS loeuf_null,
  count(*) FILTER (WHERE loeuf < 0.6) AS loeuf_lt_0_6,
  count(*) FILTER (WHERE loeuf < 0.45) AS loeuf_lt_0_45,
  count(*) FILTER (WHERE loeuf = 0.6) AS loeuf_eq_0_6,
  count(*) FILTER (WHERE loeuf IS NOT NULL AND loeuf_decile IS NULL) AS loeuf_present_decile_null,
  count(*) FILTER (WHERE loeuf < 0.6 AND loeuf_decile IS NULL) AS loeuf_lt_0_6_decile_null,
  count(*) FILTER (
    WHERE loeuf < 0.6 AND loeuf_decile IS NULL AND plp_stars_ge_2 >= 1
  ) AS loeuf_lt_0_6_decile_null_and_hc,
  count(*) FILTER (WHERE n_tissues_tpm_ge_1 IS NULL) AS expression_null,
  count(*) FILTER (WHERE n_tissues_tpm_ge_1 IS NOT NULL) AS expression_present,
  count(*) FILTER (WHERE n_tissues_tpm_ge_1 >= 27) AS breadth_ge_27,
  count(*) FILTER (WHERE n_tissues_tpm_ge_1 = 54) AS breadth_eq_54,
  count(*) FILTER (WHERE plp_stars_ge_2 >= 1) AS hc_genes,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS symbol_collision_genes
FROM scored;

-- query: breadth_median
SELECT n_tissues_tpm_ge_1 AS breadth_median
FROM (
  SELECT
    n_tissues_tpm_ge_1,
    row_number() OVER (ORDER BY n_tissues_tpm_ge_1) AS rn,
    count(*) OVER () AS n
  FROM scored
  WHERE n_tissues_tpm_ge_1 IS NOT NULL
)
WHERE rn = (n + 1) / 2;

-- query: funnel
SELECT
  count(*) FILTER (WHERE loeuf < 0.6 AND plp_stars_ge_2 >= 1) AS intolerant_and_hc,
  count(*) FILTER (
    WHERE loeuf < 0.6 AND plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 IS NULL
  ) AS intolerant_hc_expression_null,
  count(*) FILTER (
    WHERE loeuf < 0.6 AND plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 < 27
  ) AS intolerant_hc_breadth_below_27,
  count(*) FILTER (
    WHERE loeuf < 0.6 AND plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 >= 27
  ) AS result_genes,
  count(*) FILTER (
    WHERE loeuf < 0.45 AND plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 >= 27
  ) AS result_if_loeuf_lt_0_45,
  count(*) FILTER (
    WHERE loeuf = 0.6 AND plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 >= 27
  ) AS excluded_loeuf_equal_0_6,
  count(*) FILTER (
    WHERE loeuf < 0.6 AND n_tissues_tpm_ge_1 >= 27
  ) AS intolerant_and_broad,
  count(*) FILTER (
    WHERE plp_stars_ge_2 >= 1 AND n_tissues_tpm_ge_1 >= 27 AND loeuf IS NULL
  ) AS hc_broad_loeuf_null,
  count(*) FILTER (
    WHERE loeuf < 0.6
      AND n_tissues_tpm_ge_1 >= 27
      AND plp_stars_ge_2 = 0
      AND plp_star_1 >= 1
      AND plp_star_0 = 0
      AND plp_star_null = 0
  ) AS split_only_star_1,
  count(*) FILTER (
    WHERE loeuf < 0.6
      AND n_tissues_tpm_ge_1 >= 27
      AND plp_stars_ge_2 = 0
      AND plp_star_1 >= 1
      AND (plp_star_0 > 0 OR plp_star_null > 0)
  ) AS star_1_plus_weaker,
  count(*) FILTER (
    WHERE loeuf < 0.6
      AND n_tissues_tpm_ge_1 >= 27
      AND plp_any_star = 0
      AND clinvar_conflicting_n > 0
  ) AS conflicting_class_only,
  count(*) FILTER (
    WHERE loeuf < 0.6
      AND n_tissues_tpm_ge_1 >= 27
      AND plp_stars_ge_2 = 0
      AND plp_star_1 = 0
      AND (plp_star_0 > 0 OR plp_star_null > 0)
  ) AS weaker_review_only
FROM scored;

-- query: null_loeuf_hc_broad
SELECT
  count(*) FILTER (WHERE k.bare_ensg IS NULL) AS no_constraint_row,
  count(*) FILTER (WHERE k.bare_ensg IS NOT NULL AND s.loeuf IS NULL) AS constraint_row_loeuf_null
FROM scored s
LEFT JOIN constraint_gene k ON s.ensg = k.bare_ensg
WHERE s.plp_stars_ge_2 >= 1
  AND s.n_tissues_tpm_ge_1 >= 27
  AND s.loeuf IS NULL;

-- query: result_summary
SELECT
  count(*) AS result_genes,
  min(loeuf) AS loeuf_min,
  max(loeuf) AS loeuf_max,
  min(plp_stars_ge_2) AS hc_min,
  max(plp_stars_ge_2) AS hc_max,
  min(n_tissues_tpm_ge_1) AS breadth_min,
  max(n_tissues_tpm_ge_1) AS breadth_max,
  sum(plp_stars_ge_2) AS hc_alleles,
  sum(pathogenic_stars_ge_2) AS pathogenic_alleles,
  sum(likely_pathogenic_stars_ge_2) AS likely_pathogenic_alleles,
  count(*) FILTER (WHERE pathogenic_stars_ge_2 = 0) AS likely_pathogenic_only_genes,
  count(*) FILTER (WHERE plp_star_4 > 0) AS genes_with_star_4,
  count(*) FILTER (WHERE plp_star_3 > 0) AS genes_with_star_3,
  count(*) FILTER (WHERE plp_star_3 = 0 AND plp_star_4 = 0) AS genes_star_2_only,
  count(*) FILTER (WHERE plp_star_1 > 0) AS genes_also_star_1,
  count(*) FILTER (WHERE clinvar_conflicting_n > 0) AS genes_also_conflicting,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS symbol_collision_genes,
  count(*) FILTER (
    WHERE constraint_flags IS NOT NULL AND constraint_flags <> '[]'
  ) AS flagged_genes,
  count(*) FILTER (WHERE constraint_flags LIKE '%outlier_lof%') AS outlier_lof_genes,
  count(*) FILTER (WHERE loeuf_decile IS NULL) AS null_decile_genes,
  count(*) FILTER (WHERE loeuf_decile = 0) AS decile_0,
  count(*) FILTER (WHERE loeuf_decile = 1) AS decile_1,
  count(*) FILTER (WHERE loeuf_decile = 2) AS decile_2,
  count(*) FILTER (WHERE chrom = 'chrX') AS chrx,
  count(*) FILTER (WHERE chrom = 'chrY') AS chry
FROM scored
WHERE loeuf < 0.6
  AND plp_stars_ge_2 >= 1
  AND n_tissues_tpm_ge_1 >= 27;

-- query: result_medians
SELECT
  max(loeuf) FILTER (WHERE loeuf_rn = (n + 1) / 2) AS loeuf_median,
  max(plp_stars_ge_2) FILTER (WHERE hc_rn = (n + 1) / 2) AS hc_median,
  max(n_tissues_tpm_ge_1) FILTER (WHERE breadth_rn = (n + 1) / 2) AS breadth_median
FROM (
  SELECT
    loeuf,
    plp_stars_ge_2,
    n_tissues_tpm_ge_1,
    count(*) OVER () AS n,
    row_number() OVER (ORDER BY loeuf) AS loeuf_rn,
    row_number() OVER (ORDER BY plp_stars_ge_2) AS hc_rn,
    row_number() OVER (ORDER BY n_tissues_tpm_ge_1) AS breadth_rn
  FROM scored
  WHERE loeuf < 0.6
    AND plp_stars_ge_2 >= 1
    AND n_tissues_tpm_ge_1 >= 27
);

-- query: result_flags
SELECT symbol, ensg, loeuf, plp_stars_ge_2, n_tissues_tpm_ge_1, constraint_flags
FROM scored
WHERE loeuf < 0.6
  AND plp_stars_ge_2 >= 1
  AND n_tissues_tpm_ge_1 >= 27
  AND constraint_flags IS NOT NULL
  AND constraint_flags <> '[]'
ORDER BY loeuf, plp_stars_ge_2 DESC, n_tissues_tpm_ge_1 DESC, symbol, ensg;

-- query: hc_broad_null_loeuf_genes
SELECT s.symbol, s.ensg, s.chrom, s.plp_stars_ge_2, s.n_tissues_tpm_ge_1,
       k.constraint_flags,
       k."lof.exp" AS lof_exp,
       k.loeuf AS constraint_loeuf,
       CASE WHEN k.bare_ensg IS NULL THEN 'no_constraint_row' ELSE 'constraint_row_loeuf_null' END AS reason
FROM scored s
LEFT JOIN constraint_gene k ON s.ensg = k.bare_ensg
WHERE s.plp_stars_ge_2 >= 1
  AND s.n_tissues_tpm_ge_1 >= 27
  AND s.loeuf IS NULL
  AND k.bare_ensg IS NOT NULL
ORDER BY s.symbol;

-- query: examples_low_loeuf_excluded
SELECT symbol, ensg, loeuf, plp_stars_ge_2, plp_any_star, n_tissues_tpm_ge_1
FROM scored
WHERE symbol IN ('XPO6', 'EIF4G1', 'MED13L')
ORDER BY symbol;

-- query: intolerant_hc_no_expression
SELECT symbol, ensg, loeuf, plp_stars_ge_2, n_tissues_tpm_ge_1
FROM scored
WHERE loeuf < 0.6
  AND plp_stars_ge_2 >= 1
  AND n_tissues_tpm_ge_1 IS NULL
ORDER BY loeuf, symbol;

-- query: loeuf_equal_0_6_otherwise_passing
SELECT symbol, ensg, loeuf, plp_stars_ge_2, n_tissues_tpm_ge_1
FROM scored
WHERE loeuf = 0.6
  AND plp_stars_ge_2 >= 1
  AND n_tissues_tpm_ge_1 >= 27
ORDER BY plp_stars_ge_2 DESC, symbol;

-- query: split_summary
SELECT
  count(*) AS split_genes,
  min(loeuf) AS loeuf_min,
  max(loeuf) AS loeuf_max,
  sum(plp_star_1) AS star_1_alleles,
  min(plp_star_1) AS star_1_min,
  max(plp_star_1) AS star_1_max,
  count(*) FILTER (WHERE clinvar_conflicting_n > 0) AS also_conflicting,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS symbol_collision_genes,
  count(*) FILTER (WHERE loeuf_decile = 0) AS decile_0,
  count(*) FILTER (WHERE loeuf_decile = 1) AS decile_1,
  count(*) FILTER (WHERE loeuf_decile = 2) AS decile_2,
  count(*) FILTER (WHERE loeuf_decile IS NULL) AS decile_null
FROM scored
WHERE loeuf < 0.6
  AND n_tissues_tpm_ge_1 >= 27
  AND plp_stars_ge_2 = 0
  AND plp_star_1 >= 1
  AND plp_star_0 = 0
  AND plp_star_null = 0;

-- query: result_genes
SELECT
  row_number() OVER (
    ORDER BY loeuf, plp_stars_ge_2 DESC, n_tissues_tpm_ge_1 DESC, symbol, ensg
  ) AS rank,
  symbol,
  ensg,
  loeuf,
  loeuf_decile,
  plp_stars_ge_2,
  n_tissues_tpm_ge_1
FROM scored
WHERE loeuf < 0.6
  AND plp_stars_ge_2 >= 1
  AND n_tissues_tpm_ge_1 >= 27
ORDER BY loeuf, plp_stars_ge_2 DESC, n_tissues_tpm_ge_1 DESC, symbol, ensg;

-- query: split_genes
SELECT
  row_number() OVER (
    ORDER BY loeuf, plp_star_1 DESC, n_tissues_tpm_ge_1 DESC, symbol, ensg
  ) AS rank,
  symbol,
  ensg,
  loeuf,
  loeuf_decile,
  plp_star_1,
  n_tissues_tpm_ge_1
FROM scored
WHERE loeuf < 0.6
  AND n_tissues_tpm_ge_1 >= 27
  AND plp_stars_ge_2 = 0
  AND plp_star_1 >= 1
  AND plp_star_0 = 0
  AND plp_star_null = 0
ORDER BY loeuf, plp_star_1 DESC, n_tissues_tpm_ge_1 DESC, symbol, ensg;
