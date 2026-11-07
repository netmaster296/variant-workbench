-- Genes with no expression measurement, no constraint row, no ClinVar
-- submission, no GWAS locus by either rule, and no MaveDB assay.
-- Read-only against warehouse/marts/gene_mart.parquet. Does not write a mart.
--
-- Headline "no expression" is a missing GTEx row (all five expression columns
-- null). A gene GTEx tested with maximum median TPM 0 has tau null as well,
-- but max_median_tpm is 0, not null, and is not in the headline set.
-- Headline "no constraint call" is a missing constraint row (all six constraint
-- columns null). A row with a null LOEUF is not that state.
-- ClinVar, both GWAS counts, and all four MaveDB counts are 0. Those columns
-- are never null on this mart; the predicates still accept null.

SET threads = 8;
SET preserve_insertion_order = false;

CREATE OR REPLACE TEMP VIEW gene AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet');

CREATE OR REPLACE TEMP VIEW called AS
SELECT
  ensg,
  symbol,
  chrom,
  start,
  "end",
  strand,
  biotype,
  mane_enst,
  symbol_n_gencode,
  ("end" - start + 1) AS span_bp,
  (symbol = ensg) AS symbol_is_ensg,
  max_median_tpm,
  max_tissue,
  tau,
  n_tissues_tpm_ge_1,
  n_tissues_tpm_ge_10,
  loeuf,
  loeuf_decile,
  "lof.pLI",
  "mis.z_score",
  "syn.z_score",
  constraint_flags,
  clinvar_n,
  clinvar_germline_pathogenic_n,
  clinvar_likely_pathogenic_n,
  clinvar_vus_n,
  clinvar_benign_n,
  clinvar_likely_benign_n,
  clinvar_conflicting_n,
  clinvar_drug_response_n,
  gwas_loci_by_overlap,
  gwas_traits_by_overlap,
  gwas_loci_by_mapped_gene,
  gwas_traits_by_mapped_gene,
  mavedb_experiments_by_symbol,
  mavedb_mapped_variants_by_symbol,
  mavedb_variant_keys_in_span,
  mavedb_score_sets_in_span,
  eqtl_tissues,
  sqtl_tissues,
  (
    max_median_tpm IS NULL
    AND max_tissue IS NULL
    AND tau IS NULL
    AND n_tissues_tpm_ge_1 IS NULL
    AND n_tissues_tpm_ge_10 IS NULL
  ) AS no_expression_row,
  (
    max_median_tpm = 0
    AND max_tissue IS NULL
    AND tau IS NULL
    AND n_tissues_tpm_ge_1 = 0
    AND n_tissues_tpm_ge_10 = 0
  ) AS expression_all_zero,
  (tau IS NOT NULL) AS tau_present,
  (
    loeuf IS NULL
    AND loeuf_decile IS NULL
    AND "lof.pLI" IS NULL
    AND "mis.z_score" IS NULL
    AND "syn.z_score" IS NULL
    AND constraint_flags IS NULL
  ) AS no_constraint_row,
  (loeuf IS NULL AND constraint_flags IS NOT NULL) AS constraint_row_loeuf_null,
  (clinvar_n IS NULL OR clinvar_n = 0) AS no_clinvar,
  (
    (gwas_loci_by_overlap IS NULL OR gwas_loci_by_overlap = 0)
    AND (gwas_loci_by_mapped_gene IS NULL OR gwas_loci_by_mapped_gene = 0)
  ) AS no_gwas,
  (
    (mavedb_experiments_by_symbol IS NULL OR mavedb_experiments_by_symbol = 0)
    AND (mavedb_mapped_variants_by_symbol IS NULL OR mavedb_mapped_variants_by_symbol = 0)
    AND (mavedb_variant_keys_in_span IS NULL OR mavedb_variant_keys_in_span = 0)
    AND (mavedb_score_sets_in_span IS NULL OR mavedb_score_sets_in_span = 0)
  ) AS no_mavedb
FROM gene;

CREATE OR REPLACE TEMP VIEW in_set AS
SELECT *
FROM called
WHERE no_expression_row
  AND no_constraint_row
  AND no_clinvar
  AND no_gwas
  AND no_mavedb;

-- Denominators.
SELECT
  'denominator' AS result,
  count(*) AS genes,
  count(*) FILTER (WHERE biotype = 'protein_coding') AS protein_coding,
  count(DISTINCT ensg) AS distinct_ensg
FROM gene;

-- Each state, both denominators. Expression states partition the mart.
-- Constraint states partition the mart. ClinVar, GWAS, and MaveDB nulls are
-- counted beside the zeros.
SELECT
  'margins' AS result,
  count(*) FILTER (WHERE no_expression_row) AS expression_row_absent,
  count(*) FILTER (WHERE expression_all_zero) AS expression_tested_all_zero,
  count(*) FILTER (WHERE tau_present) AS tau_present,
  count(*) FILTER (WHERE no_constraint_row) AS constraint_row_absent,
  count(*) FILTER (WHERE constraint_row_loeuf_null) AS constraint_row_loeuf_null,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_present,
  count(*) FILTER (WHERE clinvar_n IS NULL) AS clinvar_null,
  count(*) FILTER (WHERE clinvar_n = 0) AS clinvar_zero,
  count(*) FILTER (WHERE gwas_loci_by_overlap IS NULL) AS gwas_overlap_null,
  count(*) FILTER (WHERE gwas_loci_by_mapped_gene IS NULL) AS gwas_mapped_null,
  count(*) FILTER (WHERE no_gwas) AS gwas_neither,
  count(*) FILTER (WHERE mavedb_experiments_by_symbol IS NULL) AS mavedb_experiments_null,
  count(*) FILTER (WHERE mavedb_variant_keys_in_span IS NULL) AS mavedb_span_null,
  count(*) FILTER (WHERE no_mavedb) AS mavedb_all_zero,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND no_expression_row) AS coding_expression_row_absent,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND expression_all_zero) AS coding_expression_tested_all_zero,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND tau_present) AS coding_tau_present,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND no_constraint_row) AS coding_constraint_row_absent,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND constraint_row_loeuf_null) AS coding_constraint_row_loeuf_null,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND loeuf IS NOT NULL) AS coding_loeuf_present,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND clinvar_n = 0) AS coding_clinvar_zero,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND no_gwas) AS coding_gwas_neither,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND no_mavedb) AS coding_mavedb_all_zero
FROM called;

-- Headline intersection, the all-zero-TPM alternative, and the null-LOEUF alternative.
SELECT
  'intersection' AS result,
  count(*) AS no_evidence,
  count(*) FILTER (WHERE biotype = 'protein_coding') AS protein_coding_no_evidence,
  round(count(*) * 1.0 / (SELECT count(*) FROM gene), 6) AS fraction_of_genes,
  round(
    count(*) FILTER (WHERE biotype = 'protein_coding') * 1.0
    / (SELECT count(*) FROM gene WHERE biotype = 'protein_coding'),
    6
  ) AS fraction_of_protein_coding,
  (
    SELECT count(*)
    FROM called
    WHERE (no_expression_row OR expression_all_zero)
      AND no_constraint_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS if_all_zero_counts_as_no_expression,
  (
    SELECT count(*)
    FROM called
    WHERE biotype = 'protein_coding'
      AND (no_expression_row OR expression_all_zero)
      AND no_constraint_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS coding_if_all_zero_counts_as_no_expression,
  (
    SELECT count(*)
    FROM called
    WHERE no_expression_row
      AND loeuf IS NULL
      AND no_clinvar AND no_gwas AND no_mavedb
  ) AS if_any_null_loeuf_counts_as_no_constraint,
  (
    SELECT count(*)
    FROM called
    WHERE biotype = 'protein_coding'
      AND no_expression_row
      AND loeuf IS NULL
      AND no_clinvar AND no_gwas AND no_mavedb
  ) AS coding_if_any_null_loeuf_counts_as_no_constraint
FROM in_set;

-- QTL columns are not filters. On this set they are null, not zero.
SELECT
  'qtl_on_set' AS result,
  count(*) FILTER (WHERE eqtl_tissues IS NULL) AS eqtl_null,
  count(*) FILTER (WHERE eqtl_tissues = 0) AS eqtl_zero,
  count(*) FILTER (WHERE eqtl_tissues > 0) AS eqtl_positive,
  count(*) FILTER (WHERE sqtl_tissues IS NULL) AS sqtl_null,
  count(*) FILTER (WHERE sqtl_tissues = 0) AS sqtl_zero,
  count(*) FILTER (WHERE sqtl_tissues > 0) AS sqtl_positive
FROM in_set;

-- How many genes each filter removes from the other four.
SELECT
  'drop_one_filter' AS result,
  count(*) FILTER (
    WHERE no_constraint_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS without_requiring_expression_absent,
  count(*) FILTER (
    WHERE no_expression_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS without_requiring_constraint_row_absent,
  count(*) FILTER (
    WHERE no_expression_row AND no_constraint_row AND no_gwas AND no_mavedb
  ) AS without_requiring_clinvar_zero,
  count(*) FILTER (
    WHERE no_expression_row AND no_constraint_row AND no_clinvar AND no_mavedb
  ) AS without_requiring_gwas_zero,
  count(*) FILTER (
    WHERE no_expression_row AND no_constraint_row AND no_clinvar AND no_gwas
  ) AS without_requiring_mavedb_zero,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND no_constraint_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS coding_without_requiring_expression_absent,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND no_expression_row AND no_clinvar AND no_gwas AND no_mavedb
  ) AS coding_without_requiring_constraint_row_absent,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND no_expression_row AND no_constraint_row AND no_gwas AND no_mavedb
  ) AS coding_without_requiring_clinvar_zero,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND no_expression_row AND no_constraint_row AND no_clinvar AND no_mavedb
  ) AS coding_without_requiring_gwas_zero,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND no_expression_row AND no_constraint_row AND no_clinvar AND no_gwas
  ) AS coding_without_requiring_mavedb_zero
FROM called;

SELECT
  'set_shape' AS result,
  count(*) AS n,
  count(*) FILTER (WHERE symbol_is_ensg) AS symbol_equals_ensg,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND symbol_is_ensg) AS coding_symbol_equals_ensg,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS shared_symbol,
  count(*) FILTER (WHERE biotype = 'lncRNA' AND symbol_is_ensg) AS lncrna_symbol_equals_ensg,
  min(span_bp) AS min_span_bp,
  median(span_bp) AS median_span_bp,
  max(span_bp) AS max_span_bp,
  count(*) FILTER (WHERE chrom = 'chrX') AS chrX,
  count(*) FILTER (WHERE chrom = 'chrY') AS chrY,
  count(*) FILTER (WHERE chrom = 'chr9') AS chr9,
  count(*) FILTER (WHERE chrom = 'chrM') AS chrM,
  (SELECT count(*) FROM gene WHERE chrom = 'chrX') AS chrX_genes,
  (SELECT count(*) FROM gene WHERE chrom = 'chrY') AS chrY_genes,
  (SELECT count(*) FROM gene WHERE chrom = 'chr9') AS chr9_genes,
  (SELECT count(*) FROM gene WHERE chrom = 'chrM') AS chrM_genes
FROM in_set;

-- Biotypes in the set. Biotypes with zero rows are omitted here.
SELECT
  d.biotype,
  h.n_in_set,
  d.n_biotype,
  round(h.n_in_set * 1.0 / d.n_biotype, 6) AS fraction_of_biotype
FROM (
  SELECT biotype, count(*) AS n_biotype
  FROM gene
  GROUP BY biotype
) d
JOIN (
  SELECT biotype, count(*) AS n_in_set
  FROM in_set
  GROUP BY biotype
) h ON d.biotype = h.biotype
ORDER BY h.n_in_set DESC, d.biotype;

-- Every protein-coding gene in the set.
SELECT
  'protein_coding' AS result,
  ensg,
  symbol,
  chrom,
  start,
  "end",
  strand,
  span_bp,
  mane_enst,
  symbol_n_gencode
FROM in_set
WHERE biotype = 'protein_coding'
ORDER BY chrom, start, ensg;

-- One example from each biotype: first symbol, then ENSG.
SELECT
  biotype,
  ensg,
  symbol,
  chrom,
  start,
  "end",
  span_bp
FROM (
  SELECT
    *,
    row_number() OVER (PARTITION BY biotype ORDER BY symbol, ensg) AS rn
  FROM in_set
)
WHERE rn = 1
ORDER BY biotype;

-- Column consistency the headline depends on. Expect zeros.
SELECT
  'column_consistency' AS result,
  count(*) FILTER (
    WHERE clinvar_n = 0 AND (
      clinvar_germline_pathogenic_n <> 0
      OR clinvar_likely_pathogenic_n <> 0
      OR clinvar_vus_n <> 0
      OR clinvar_benign_n <> 0
      OR clinvar_likely_benign_n <> 0
      OR clinvar_conflicting_n <> 0
      OR clinvar_drug_response_n <> 0
    )
  ) AS clinvar_subclass_mismatch,
  count(*) FILTER (
    WHERE gwas_loci_by_overlap = 0 AND gwas_traits_by_overlap <> 0
  ) AS gwas_overlap_trait_mismatch,
  count(*) FILTER (
    WHERE gwas_loci_by_mapped_gene = 0 AND gwas_traits_by_mapped_gene <> 0
  ) AS gwas_mapped_trait_mismatch,
  count(*) FILTER (
    WHERE no_expression_row AND NOT (
      max_tissue IS NULL AND tau IS NULL
      AND n_tissues_tpm_ge_1 IS NULL AND n_tissues_tpm_ge_10 IS NULL
      AND max_median_tpm IS NULL
    )
  ) AS expression_flag_mismatch
FROM called;

-- Join checks against the source parquet. Mismatch columns must be 0.
-- A source row that has no GENCODE gene is an orphan, not a mart miss.
SELECT
  'join_expression' AS result,
  count(*) FILTER (WHERE m.max_median_tpm IS NULL AND e.ensg IS NOT NULL) AS mart_null_source_hit,
  count(*) FILTER (WHERE m.max_median_tpm IS NOT NULL AND e.ensg IS NULL) AS mart_hit_source_miss,
  count(*) FILTER (WHERE e.ensg IS NOT NULL AND m.max_median_tpm IS DISTINCT FROM e.max_median_tpm) AS tpm_diff,
  count(*) FILTER (WHERE e.ensg IS NOT NULL AND m.tau IS DISTINCT FROM e.tau) AS tau_diff,
  count(*) FILTER (WHERE e.ensg IS NOT NULL AND m.max_tissue IS DISTINCT FROM e.max_tissue) AS tissue_diff,
  count(*) FILTER (WHERE e.ensg IS NOT NULL AND m.n_tissues_tpm_ge_1 IS DISTINCT FROM e.n_tissues_tpm_ge_1) AS ge1_diff,
  count(*) FILTER (WHERE e.ensg IS NOT NULL AND m.n_tissues_tpm_ge_10 IS DISTINCT FROM e.n_tissues_tpm_ge_10) AS ge10_diff
FROM gene m
LEFT JOIN (
  SELECT *
  FROM read_parquet('D:/variant-workbench/warehouse/gtex/gene_expression.parquet')
  WHERE ensg IS NOT NULL
) e ON m.ensg = e.ensg;

SELECT
  'join_constraint' AS result,
  count(*) FILTER (WHERE m.constraint_flags IS NULL AND c.bare_ensg IS NOT NULL) AS mart_miss_source_hit,
  count(*) FILTER (WHERE m.constraint_flags IS NOT NULL AND c.bare_ensg IS NULL) AS mart_hit_source_miss,
  count(*) FILTER (WHERE c.bare_ensg IS NOT NULL AND m.loeuf IS DISTINCT FROM c.loeuf) AS loeuf_diff,
  count(*) FILTER (
    WHERE c.bare_ensg IS NOT NULL
      AND m.loeuf_decile IS DISTINCT FROM CAST(c.loeuf_decile AS INTEGER)
  ) AS decile_diff,
  count(*) FILTER (WHERE c.bare_ensg IS NOT NULL AND m."lof.pLI" IS DISTINCT FROM c."lof.pLI") AS pli_diff,
  count(*) FILTER (WHERE c.bare_ensg IS NOT NULL AND m."mis.z_score" IS DISTINCT FROM c."mis.z_score") AS mis_diff,
  count(*) FILTER (WHERE c.bare_ensg IS NOT NULL AND m."syn.z_score" IS DISTINCT FROM c."syn.z_score") AS syn_diff,
  count(*) FILTER (WHERE c.bare_ensg IS NOT NULL AND m.constraint_flags IS DISTINCT FROM c.constraint_flags) AS flags_diff
FROM gene m
LEFT JOIN (
  SELECT *
  FROM read_parquet('D:/variant-workbench/warehouse/constraint/gene_constraint.parquet')
  WHERE bare_ensg IS NOT NULL
) c ON m.ensg = c.bare_ensg;

WITH tokens AS (
  SELECT DISTINCT dedup_key, trim(token) AS token
  FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet'),
    unnest(string_split(coalesce(GeneSymbol, ''), ';')) AS u(token)
),
counted AS (
  SELECT token AS symbol, count(*)::BIGINT AS clinvar_n
  FROM tokens
  WHERE token <> ''
    AND token NOT IN ('.', '-')
    AND token IN (
      SELECT symbol
      FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet')
      WHERE symbol IS NOT NULL AND symbol <> ''
    )
  GROUP BY token
)
SELECT
  'join_clinvar' AS result,
  (
    SELECT count(*)
    FROM gene m
    LEFT JOIN counted c ON m.symbol = c.symbol
    WHERE m.clinvar_n IS DISTINCT FROM coalesce(c.clinvar_n, 0)
  ) AS genes_count_mismatch,
  (SELECT count(*) FROM gene WHERE clinvar_n IS NULL) AS mart_clinvar_null;

WITH recomputed AS (
  SELECT
    g.ensg,
    count(DISTINCT concat_ws(
      '|',
      coalesce(l.trait_uri, ''),
      l.chrom,
      cast(l.lead_pos AS VARCHAR),
      coalesce(l.lead_rsid, '')
    ))::BIGINT AS gwas_loci_by_overlap
  FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet') g
  JOIN read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') l
    ON g.chrom = l.chrom
   AND (
        (l.start <= g."end" AND l."end" >= g.start)
        OR (l.lead_pos BETWEEN g.start AND g."end")
   )
  GROUP BY g.ensg
)
SELECT
  'join_gwas_overlap' AS result,
  (
    SELECT count(*)
    FROM gene m
    LEFT JOIN recomputed r ON m.ensg = r.ensg
    WHERE m.gwas_loci_by_overlap IS DISTINCT FROM coalesce(r.gwas_loci_by_overlap, 0)
  ) AS genes_count_mismatch,
  (SELECT count(*) FROM gene WHERE gwas_loci_by_overlap IS NULL) AS mart_null;

WITH recomputed AS (
  SELECT
    g.ensg,
    count(DISTINCT concat_ws(
      '|',
      coalesce(t.trait_uri, ''),
      coalesce(t.chrom, ''),
      cast(t.lead_pos AS VARCHAR),
      coalesce(t.lead_rsid, '')
    ))::BIGINT AS gwas_loci_by_mapped_gene
  FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet') g
  JOIN (
    SELECT trait_uri, chrom, lead_pos, lead_rsid, unnest(mapped_genes.tokens) AS token
    FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')
  ) t ON t.token = g.symbol
  GROUP BY g.ensg
)
SELECT
  'join_gwas_mapped' AS result,
  (
    SELECT count(*)
    FROM gene m
    LEFT JOIN recomputed r ON m.ensg = r.ensg
    WHERE m.gwas_loci_by_mapped_gene IS DISTINCT FROM coalesce(r.gwas_loci_by_mapped_gene, 0)
  ) AS genes_count_mismatch,
  (SELECT count(*) FROM gene WHERE gwas_loci_by_mapped_gene IS NULL) AS mart_null;

WITH recomputed AS (
  SELECT
    target_gene AS symbol,
    sum(n_score_sets)::BIGINT AS mavedb_experiments_by_symbol,
    sum(n_mapped_variants)::BIGINT AS mavedb_mapped_variants_by_symbol
  FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet')
  GROUP BY target_gene
)
SELECT
  'join_mavedb_symbol' AS result,
  (
    SELECT count(*)
    FROM gene m
    LEFT JOIN recomputed r ON m.symbol = r.symbol
    WHERE m.mavedb_experiments_by_symbol IS DISTINCT FROM coalesce(r.mavedb_experiments_by_symbol, 0)
       OR m.mavedb_mapped_variants_by_symbol IS DISTINCT FROM coalesce(r.mavedb_mapped_variants_by_symbol, 0)
  ) AS genes_count_mismatch,
  (
    SELECT count(*)
    FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet') c
    WHERE c.target_gene IS NOT NULL
      AND c.target_gene NOT IN (
        SELECT symbol
        FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet')
        WHERE symbol IS NOT NULL
      )
  ) AS coverage_rows_symbol_not_in_gencode,
  (
    SELECT count(*)
    FROM read_parquet('D:/variant-workbench/warehouse/mavedb/gene_coverage.parquet')
  ) AS coverage_rows;

WITH variants AS (
  SELECT DISTINCT variant_key, score_set_urn, chrom, pos
  FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
  WHERE variant_key IS NOT NULL AND chrom IS NOT NULL AND pos IS NOT NULL
),
recomputed AS (
  SELECT
    g.ensg,
    count(DISTINCT v.variant_key)::BIGINT AS mavedb_variant_keys_in_span,
    count(DISTINCT v.score_set_urn)::BIGINT AS mavedb_score_sets_in_span
  FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet') g
  JOIN variants v
    ON v.chrom = g.chrom
   AND v.pos BETWEEN g.start AND g."end"
  GROUP BY g.ensg
)
SELECT
  'join_mavedb_span' AS result,
  (
    SELECT count(*)
    FROM gene m
    LEFT JOIN recomputed r ON m.ensg = r.ensg
    WHERE m.mavedb_variant_keys_in_span IS DISTINCT FROM coalesce(r.mavedb_variant_keys_in_span, 0)
       OR m.mavedb_score_sets_in_span IS DISTINCT FROM coalesce(r.mavedb_score_sets_in_span, 0)
  ) AS genes_count_mismatch;
