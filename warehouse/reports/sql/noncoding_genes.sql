-- Noncoding genes. Read-only. Does not write or edit a mart.
-- Question, filters, and caveats: warehouse/reports/noncoding_genes.md
--
-- Source is the full gene mart, not gene_mart_protein_coding.parquet.
-- Included biotypes are every biotype except protein_coding and the
-- immunoglobulin / T-cell receptor gene segments, which encode polypeptides.
-- Nearby is the same chrom and a span distance of at most 10,000 bp.
-- Distance is 0 when the inclusive spans share a base. Otherwise it is
-- the number of bases strictly between them (right_start - left_end - 1).
-- The neighbor used for "nearby" has biotype protein_coding only.
-- "Has expression" means max_median_tpm > 0. A null TPM is not in the
-- GTEx v8 median table. A zero was in the table and was not detected.
-- GWAS overlap and mapped-gene counts are the mart columns. The checks
-- below recompute them from loci.parquet.

PRAGMA threads=4;
SET memory_limit='12GB';

CREATE OR REPLACE TEMP TABLE mart AS
SELECT
  ensg,
  symbol,
  chrom,
  start,
  "end",
  strand,
  biotype,
  symbol_n_gencode,
  loeuf,
  max_tissue,
  max_median_tpm,
  n_tissues_tpm_ge_1,
  n_tissues_tpm_ge_10,
  tau,
  gwas_loci_by_overlap,
  gwas_traits_by_overlap,
  gwas_loci_by_mapped_gene,
  gwas_traits_by_mapped_gene,
  CASE
    WHEN biotype = 'protein_coding' THEN 'excluded_protein_coding'
    WHEN biotype IN (
      'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene',
      'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene'
    ) THEN 'excluded_ig_tr_segment'
    WHEN biotype = 'lncRNA' THEN 'lncRNA'
    WHEN biotype LIKE '%pseudogene' THEN 'pseudogene'
    WHEN biotype = 'TEC' THEN 'TEC'
    WHEN biotype = 'artifact' THEN 'artifact'
    ELSE 'other_ncrna'
  END AS cls
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet');

CREATE OR REPLACE TEMP TABLE nc AS
SELECT *
FROM mart
WHERE cls IN ('lncRNA', 'other_ncrna', 'pseudogene', 'TEC', 'artifact');

-- Identity: mart coordinates against interval genes. Expression: mart
-- columns against gene_expression.parquet. A non-zero mismatch stops the claim.
SELECT
  'qc_identity' AS section,
  (SELECT count(*) FROM mart) AS mart_n,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet')) AS interval_n,
  (
    SELECT count(*)
    FROM mart m
    FULL OUTER JOIN read_parquet('D:/variant-workbench/warehouse/intervals/genes.parquet') g
      ON m.ensg = g.ensg
    WHERE m.ensg IS NULL
       OR g.ensg IS NULL
       OR m.chrom IS DISTINCT FROM g.chrom
       OR m.start IS DISTINCT FROM g.start
       OR m."end" IS DISTINCT FROM g."end"
       OR m.symbol IS DISTINCT FROM g.symbol
       OR m.biotype IS DISTINCT FROM g.biotype
       OR m.strand IS DISTINCT FROM g.strand
  ) AS identity_mismatch;

SELECT
  'qc_expression' AS section,
  count(*) AS mart_n,
  count(*) FILTER (WHERE m.max_median_tpm IS NOT NULL) AS mart_in_gtex,
  count(*) FILTER (WHERE e.ensg IS NOT NULL) AS file_in_gtex,
  count(*) FILTER (WHERE m.max_median_tpm IS NOT NULL AND e.ensg IS NULL) AS mart_only,
  count(*) FILTER (WHERE m.max_median_tpm IS NULL AND e.ensg IS NOT NULL) AS file_only,
  count(*) FILTER (
    WHERE e.ensg IS NOT NULL AND (
      m.max_median_tpm IS DISTINCT FROM e.max_median_tpm
      OR m.max_tissue IS DISTINCT FROM e.max_tissue
      OR m.n_tissues_tpm_ge_1 IS DISTINCT FROM e.n_tissues_tpm_ge_1
      OR m.n_tissues_tpm_ge_10 IS DISTINCT FROM e.n_tissues_tpm_ge_10
      OR m.tau IS DISTINCT FROM e.tau
    )
  ) AS value_mismatch
FROM mart m
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/gtex/gene_expression.parquet') e
  ON m.ensg = e.ensg;

CREATE OR REPLACE TEMP TABLE nearest AS
SELECT
  ensg,
  min(distance_bp) AS nearest_bp,
  bool_or(span_overlap) AS span_overlap,
  count(*) FILTER (WHERE distance_bp <= 10000) AS n_coding_10kb,
  count(*) FILTER (WHERE distance_bp <= 100000) AS n_coding_100kb,
  count(*) FILTER (WHERE distance_bp = 0 AND NOT span_overlap) AS n_abut
FROM (
  SELECT
    nc.ensg,
    (c."end" >= nc.start AND c.start <= nc."end") AS span_overlap,
    CASE
      WHEN c."end" >= nc.start AND c.start <= nc."end" THEN 0
      WHEN c."end" < nc.start THEN nc.start - c."end" - 1
      ELSE c.start - nc."end" - 1
    END AS distance_bp
  FROM nc
  JOIN mart c
    ON nc.chrom = c.chrom
   AND c.biotype = 'protein_coding'
)
GROUP BY ensg;

CREATE OR REPLACE TEMP TABLE nearest_igtr AS
SELECT
  nc.ensg,
  min(
    CASE
      WHEN g."end" >= nc.start AND g.start <= nc."end" THEN 0
      WHEN g."end" < nc.start THEN nc.start - g."end" - 1
      ELSE g.start - nc."end" - 1
    END
  ) AS nearest_igtr_bp
FROM nc
JOIN mart g
  ON nc.chrom = g.chrom
 AND g.cls = 'excluded_ig_tr_segment'
GROUP BY nc.ensg;

CREATE OR REPLACE TEMP TABLE nc_near AS
SELECT
  nc.*,
  nr.nearest_bp,
  coalesce(nr.span_overlap, false) AS span_overlap,
  coalesce(nr.n_coding_10kb, 0) AS n_coding_10kb,
  coalesce(nr.n_coding_100kb, 0) AS n_coding_100kb,
  coalesce(nr.n_abut, 0) AS n_abut,
  ig.nearest_igtr_bp,
  (nc.max_median_tpm IS NOT NULL) AS in_gtex,
  (nc.max_median_tpm = 0) AS gtex_zero,
  (nc.max_median_tpm > 0) AS has_expr,
  (coalesce(nc.n_tissues_tpm_ge_1, 0) > 0) AS tpm_ge_1,
  (coalesce(nc.n_tissues_tpm_ge_10, 0) > 0) AS tpm_ge_10,
  (nc.gwas_loci_by_overlap > 0) AS has_overlap,
  (nc.gwas_loci_by_mapped_gene > 0) AS has_mapped,
  (coalesce(nr.nearest_bp, 1000000000) <= 10000) AS within_10kb,
  (coalesce(nr.span_overlap, false) = false AND coalesce(nr.nearest_bp, 1000000000) <= 10000) AS flank_10kb,
  (coalesce(nr.nearest_bp, 1000000000) <= 100000) AS within_100kb
FROM nc
LEFT JOIN nearest nr USING (ensg)
LEFT JOIN nearest_igtr ig USING (ensg);

SELECT
  'qc_nearest' AS section,
  (SELECT count(*) FROM nc) AS included_n,
  (SELECT count(*) FROM nearest) AS nearest_n,
  (SELECT count(*) FROM nc ANTI JOIN nearest USING (ensg)) AS no_coding_on_chrom,
  (SELECT min(nearest_bp) FROM nearest) AS min_distance_bp,
  (SELECT max(nearest_bp) FROM nearest) AS max_distance_bp,
  (SELECT count(*) FROM nc_near WHERE n_abut > 0) AS genes_that_abut,
  (
    SELECT count(*)
    FROM nc_near
    WHERE nearest_igtr_bp <= 10000
      AND (nearest_bp IS NULL OR nearest_bp > 10000)
  ) AS igtr_10kb_but_not_coding_10kb;

-- Same locus key the mart used. Overlap predicate matches build_gene_mart.py.
CREATE OR REPLACE TEMP TABLE gwas_overlap_re AS
SELECT
  ensg,
  count(DISTINCT locus_id) AS n_overlap,
  count(DISTINCT locus_id) FILTER (WHERE lead_inside) AS n_lead_inside
FROM (
  SELECT
    g.ensg,
    concat_ws(
      '|',
      coalesce(l.trait_uri, ''),
      l.chrom,
      cast(l.lead_pos AS VARCHAR),
      coalesce(l.lead_rsid, '')
    ) AS locus_id,
    (l.lead_pos BETWEEN g.start AND g."end") AS lead_inside
  FROM mart g
  JOIN read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') l
    ON g.chrom = l.chrom
   AND (
        (l.start <= g."end" AND l."end" >= g.start)
        OR (l.lead_pos BETWEEN g.start AND g."end")
   )
)
GROUP BY ensg;

CREATE OR REPLACE TEMP TABLE gwas_mapped_re AS
SELECT
  g.ensg,
  count(DISTINCT concat_ws(
    '|',
    coalesce(t.trait_uri, ''),
    coalesce(t.chrom, ''),
    cast(t.lead_pos AS VARCHAR),
    coalesce(t.lead_rsid, '')
  )) AS n_mapped
FROM mart g
JOIN (
  SELECT
    trait_uri,
    chrom,
    lead_pos,
    lead_rsid,
    unnest(mapped_genes.tokens) AS token
  FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')
) t
  ON t.token = g.symbol
GROUP BY g.ensg;

SELECT
  'qc_gwas_overlap' AS section,
  count(*) FILTER (WHERE coalesce(r.n_overlap, 0) <> g.gwas_loci_by_overlap) AS gene_mismatch,
  count(*) FILTER (
    WHERE g.cls IN ('lncRNA', 'other_ncrna', 'pseudogene', 'TEC', 'artifact')
      AND coalesce(r.n_overlap, 0) <> g.gwas_loci_by_overlap
  ) AS included_mismatch,
  sum(g.gwas_loci_by_overlap) AS mart_sum,
  sum(coalesce(r.n_overlap, 0)) AS recompute_sum
FROM mart g
LEFT JOIN gwas_overlap_re r USING (ensg);

SELECT
  'qc_gwas_mapped' AS section,
  count(*) FILTER (WHERE coalesce(r.n_mapped, 0) <> g.gwas_loci_by_mapped_gene) AS gene_mismatch,
  count(*) FILTER (
    WHERE g.cls IN ('lncRNA', 'other_ncrna', 'pseudogene', 'TEC', 'artifact')
      AND coalesce(r.n_mapped, 0) <> g.gwas_loci_by_mapped_gene
  ) AS included_mismatch,
  sum(g.gwas_loci_by_mapped_gene) AS mart_sum,
  sum(coalesce(r.n_mapped, 0)) AS recompute_sum
FROM mart g
LEFT JOIN gwas_mapped_re r USING (ensg);

CREATE OR REPLACE TEMP TABLE nc_full AS
SELECT
  n.*,
  coalesce(o.n_lead_inside, 0) AS n_lead_inside
FROM nc_near n
LEFT JOIN gwas_overlap_re o USING (ensg);

SELECT
  'class_summary' AS section,
  cls,
  count(*) AS n,
  count(*) FILTER (WHERE in_gtex) AS in_gtex,
  count(*) FILTER (WHERE gtex_zero) AS gtex_zero,
  count(*) FILTER (WHERE NOT in_gtex) AS not_in_gtex,
  count(*) FILTER (WHERE has_expr) AS has_expr,
  count(*) FILTER (WHERE tpm_ge_1) AS tpm_ge_1,
  count(*) FILTER (WHERE tpm_ge_10) AS tpm_ge_10,
  count(*) FILTER (WHERE has_overlap) AS gwas_overlap,
  count(*) FILTER (WHERE n_lead_inside > 0) AS gwas_lead_inside,
  count(*) FILTER (WHERE has_overlap AND n_lead_inside = 0) AS gwas_overlap_no_lead,
  count(*) FILTER (WHERE has_mapped) AS gwas_mapped,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS gwas_mapped_unique_symbol,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode > 1) AS gwas_mapped_shared_symbol,
  count(*) FILTER (WHERE has_overlap AND has_mapped) AS overlap_and_mapped,
  count(*) FILTER (WHERE has_mapped AND NOT has_overlap) AS mapped_not_overlap,
  count(*) FILTER (WHERE span_overlap) AS coding_span_overlap,
  count(*) FILTER (WHERE within_10kb) AS coding_within_10kb,
  count(*) FILTER (WHERE flank_10kb) AS coding_flank_10kb,
  count(*) FILTER (WHERE within_100kb) AS coding_within_100kb,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_not_null,
  count(*) FILTER (WHERE symbol_n_gencode > 1) AS shared_symbol,
  median(max_median_tpm) FILTER (WHERE has_expr) AS median_max_tpm,
  quantile_cont(tau, 0.25) FILTER (WHERE has_expr) AS tau_p25,
  median(tau) FILTER (WHERE has_expr) AS tau_p50,
  quantile_cont(tau, 0.75) FILTER (WHERE has_expr) AS tau_p75,
  quantile_cont(n_tissues_tpm_ge_1, 0.25) FILTER (WHERE has_expr) AS tissues_ge_1_p25,
  median(n_tissues_tpm_ge_1) FILTER (WHERE has_expr) AS tissues_ge_1_p50,
  quantile_cont(n_tissues_tpm_ge_1, 0.75) FILTER (WHERE has_expr) AS tissues_ge_1_p75,
  quantile_cont(nearest_bp, 0.25) AS nearest_bp_p25,
  median(nearest_bp) AS nearest_bp_p50,
  quantile_cont(nearest_bp, 0.75) AS nearest_bp_p75,
  median(nearest_bp) FILTER (WHERE NOT span_overlap) AS nearest_bp_p50_if_no_span_overlap,
  median(gwas_loci_by_overlap) FILTER (WHERE has_overlap) AS median_overlap_loci,
  median(gwas_traits_by_overlap) FILTER (WHERE has_overlap) AS median_overlap_traits,
  median(gwas_loci_by_mapped_gene) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS median_mapped_loci_unique_symbol,
  median(gwas_traits_by_mapped_gene) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS median_mapped_traits_unique_symbol
FROM nc_full
GROUP BY cls
ORDER BY CASE cls
  WHEN 'lncRNA' THEN 1
  WHEN 'other_ncrna' THEN 2
  WHEN 'pseudogene' THEN 3
  WHEN 'TEC' THEN 4
  ELSE 5
END;

SELECT
  'included_total' AS section,
  count(*) AS n,
  count(*) FILTER (WHERE in_gtex) AS in_gtex,
  count(*) FILTER (WHERE gtex_zero) AS gtex_zero,
  count(*) FILTER (WHERE NOT in_gtex) AS not_in_gtex,
  count(*) FILTER (WHERE has_expr) AS has_expr,
  count(*) FILTER (WHERE tpm_ge_1) AS tpm_ge_1,
  count(*) FILTER (WHERE tpm_ge_10) AS tpm_ge_10,
  count(*) FILTER (WHERE has_overlap) AS gwas_overlap,
  count(*) FILTER (WHERE n_lead_inside > 0) AS gwas_lead_inside,
  count(*) FILTER (WHERE has_overlap AND n_lead_inside = 0) AS gwas_overlap_no_lead,
  count(*) FILTER (WHERE has_mapped) AS gwas_mapped,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS gwas_mapped_unique_symbol,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode > 1) AS gwas_mapped_shared_symbol,
  count(DISTINCT symbol) FILTER (WHERE has_mapped) AS mapped_symbols,
  count(DISTINCT symbol) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS mapped_unique_symbols,
  count(*) FILTER (WHERE has_overlap AND has_mapped) AS overlap_and_mapped,
  count(*) FILTER (WHERE has_mapped AND NOT has_overlap) AS mapped_not_overlap,
  count(*) FILTER (WHERE has_expr AND has_overlap) AS expr_and_overlap,
  count(*) FILTER (WHERE has_expr AND has_mapped) AS expr_and_mapped,
  count(*) FILTER (WHERE has_expr AND has_overlap AND has_mapped) AS expr_overlap_mapped,
  count(*) FILTER (WHERE span_overlap) AS coding_span_overlap,
  count(*) FILTER (WHERE within_10kb) AS coding_within_10kb,
  count(*) FILTER (WHERE flank_10kb) AS coding_flank_10kb,
  count(*) FILTER (WHERE within_100kb) AS coding_within_100kb,
  count(*) FILTER (WHERE has_overlap AND within_10kb) AS overlap_and_10kb,
  count(*) FILTER (WHERE has_overlap AND span_overlap) AS overlap_and_coding_span,
  count(*) FILTER (WHERE has_overlap AND NOT within_10kb) AS overlap_and_not_10kb,
  count(*) FILTER (WHERE has_mapped AND within_10kb) AS mapped_and_10kb,
  count(*) FILTER (WHERE has_mapped AND span_overlap) AS mapped_and_coding_span,
  count(*) FILTER (WHERE has_mapped AND NOT within_10kb) AS mapped_and_not_10kb,
  count(*) FILTER (WHERE has_expr AND within_10kb) AS expr_and_10kb,
  count(*) FILTER (WHERE has_expr AND tau < 0.5) AS expr_tau_lt_0_5,
  count(*) FILTER (WHERE has_expr AND tau >= 0.5 AND tau < 0.85) AS expr_tau_mid,
  count(*) FILTER (WHERE has_expr AND tau >= 0.85) AS expr_tau_ge_0_85,
  count(*) FILTER (WHERE has_expr AND n_tissues_tpm_ge_1 = 0) AS expr_never_tpm_ge_1,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_not_null,
  count(*) FILTER (WHERE loeuf IS NULL) AS loeuf_null,
  median(max_median_tpm) FILTER (WHERE has_expr) AS median_max_tpm,
  quantile_cont(tau, 0.25) FILTER (WHERE has_expr) AS tau_p25,
  median(tau) FILTER (WHERE has_expr) AS tau_p50,
  quantile_cont(tau, 0.75) FILTER (WHERE has_expr) AS tau_p75,
  quantile_cont(n_tissues_tpm_ge_1, 0.25) FILTER (WHERE has_expr) AS tissues_ge_1_p25,
  median(n_tissues_tpm_ge_1) FILTER (WHERE has_expr) AS tissues_ge_1_p50,
  quantile_cont(n_tissues_tpm_ge_1, 0.75) FILTER (WHERE has_expr) AS tissues_ge_1_p75,
  quantile_cont(nearest_bp, 0.25) AS nearest_bp_p25,
  median(nearest_bp) AS nearest_bp_p50,
  quantile_cont(nearest_bp, 0.75) AS nearest_bp_p75,
  median(nearest_bp) FILTER (WHERE NOT span_overlap) AS nearest_bp_p50_if_no_span_overlap,
  median(gwas_loci_by_overlap) FILTER (WHERE has_overlap) AS median_overlap_loci,
  median(gwas_traits_by_overlap) FILTER (WHERE has_overlap) AS median_overlap_traits,
  median(gwas_loci_by_mapped_gene) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS median_mapped_loci_unique_symbol,
  median(gwas_traits_by_mapped_gene) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS median_mapped_traits_unique_symbol
FROM nc_full;

SELECT
  'biotype' AS section,
  cls,
  biotype,
  count(*) AS n,
  count(*) FILTER (WHERE in_gtex) AS in_gtex,
  count(*) FILTER (WHERE gtex_zero) AS gtex_zero,
  count(*) FILTER (WHERE NOT in_gtex) AS not_in_gtex,
  count(*) FILTER (WHERE has_expr) AS has_expr,
  count(*) FILTER (WHERE tpm_ge_1) AS tpm_ge_1,
  count(*) FILTER (WHERE tpm_ge_10) AS tpm_ge_10,
  count(*) FILTER (WHERE has_overlap) AS gwas_overlap,
  count(*) FILTER (WHERE has_mapped) AS gwas_mapped,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS gwas_mapped_unique_symbol,
  count(*) FILTER (WHERE span_overlap) AS coding_span_overlap,
  count(*) FILTER (WHERE within_10kb) AS coding_within_10kb,
  count(*) FILTER (WHERE flank_10kb) AS coding_flank_10kb,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_not_null
FROM nc_full
GROUP BY cls, biotype
ORDER BY CASE cls
  WHEN 'lncRNA' THEN 1
  WHEN 'other_ncrna' THEN 2
  WHEN 'pseudogene' THEN 3
  WHEN 'TEC' THEN 4
  ELSE 5
END, n DESC, biotype;

SELECT
  'lnc_conditional' AS section,
  count(*) AS n,
  count(*) FILTER (WHERE has_expr) AS has_expr,
  count(*) FILTER (WHERE has_overlap) AS gwas_overlap,
  count(*) FILTER (WHERE n_lead_inside > 0) AS gwas_lead_inside,
  count(*) FILTER (WHERE has_mapped) AS gwas_mapped,
  count(*) FILTER (WHERE has_mapped AND symbol_n_gencode = 1) AS gwas_mapped_unique_symbol,
  count(*) FILTER (WHERE has_overlap AND within_10kb) AS overlap_and_10kb,
  count(*) FILTER (WHERE has_overlap AND span_overlap) AS overlap_and_coding_span,
  count(*) FILTER (WHERE has_overlap AND has_mapped) AS overlap_and_mapped,
  count(*) FILTER (WHERE has_overlap AND NOT within_10kb) AS overlap_and_not_10kb,
  count(*) FILTER (WHERE has_mapped AND within_10kb) AS mapped_and_10kb,
  count(*) FILTER (WHERE has_mapped AND span_overlap) AS mapped_and_coding_span,
  count(*) FILTER (WHERE has_mapped AND NOT within_10kb) AS mapped_and_not_10kb,
  count(*) FILTER (WHERE has_expr AND has_overlap AND has_mapped) AS expr_overlap_mapped,
  count(*) FILTER (WHERE has_expr AND tau >= 0.85) AS expr_tau_ge_0_85,
  count(*) FILTER (WHERE has_expr AND tau < 0.5) AS expr_tau_lt_0_5,
  count(*) FILTER (WHERE has_expr AND n_tissues_tpm_ge_1 = 0) AS expr_never_tpm_ge_1,
  count(*) FILTER (WHERE in_gtex AND has_expr) AS in_gtex_and_expr,
  count(*) FILTER (WHERE in_gtex AND tpm_ge_1) AS in_gtex_and_tpm_ge_1
FROM nc_full
WHERE cls = 'lncRNA';

SELECT
  'top_tissue' AS section,
  cls,
  max_tissue,
  n
FROM (
  SELECT
    cls,
    max_tissue,
    count(*) AS n,
    row_number() OVER (PARTITION BY cls ORDER BY count(*) DESC, max_tissue) AS rk
  FROM nc_full
  WHERE has_expr AND max_tissue IS NOT NULL
  GROUP BY cls, max_tissue
)
WHERE rk <= 5
ORDER BY CASE cls
  WHEN 'lncRNA' THEN 1
  WHEN 'other_ncrna' THEN 2
  WHEN 'pseudogene' THEN 3
  WHEN 'TEC' THEN 4
  ELSE 5
END, rk;

SELECT
  'repeated_mapped_symbols' AS section,
  symbol,
  any_value(cls) AS cls,
  any_value(biotype) AS biotype,
  count(*) AS genes,
  count(*) FILTER (WHERE has_mapped) AS mapped_genes,
  max(symbol_n_gencode) AS symbol_n_gencode,
  min(gwas_loci_by_mapped_gene) FILTER (WHERE has_mapped) AS min_mapped_loci,
  max(gwas_loci_by_mapped_gene) FILTER (WHERE has_mapped) AS max_mapped_loci
FROM nc_full
WHERE symbol_n_gencode > 1
GROUP BY symbol
HAVING count(*) FILTER (WHERE has_mapped) >= 10
ORDER BY mapped_genes DESC, symbol;

SELECT
  'symbol_shared_with_protein_coding' AS section,
  nc.symbol,
  nc.ensg,
  nc.gwas_loci_by_mapped_gene,
  nc.gwas_loci_by_overlap,
  nc.max_median_tpm
FROM nc_full nc
SEMI JOIN mart pc
  ON nc.symbol = pc.symbol
 AND pc.biotype = 'protein_coding'
ORDER BY nc.symbol;

-- Scale only. Not part of the noncoding denominator and not a ranking.
SELECT
  'protein_coding_scale' AS section,
  count(*) AS n,
  count(*) FILTER (WHERE max_median_tpm IS NOT NULL) AS in_gtex,
  count(*) FILTER (WHERE max_median_tpm = 0) AS gtex_zero,
  count(*) FILTER (WHERE max_median_tpm IS NULL) AS not_in_gtex,
  count(*) FILTER (WHERE max_median_tpm > 0) AS has_expr,
  count(*) FILTER (WHERE n_tissues_tpm_ge_1 > 0) AS tpm_ge_1,
  count(*) FILTER (WHERE n_tissues_tpm_ge_10 > 0) AS tpm_ge_10,
  count(*) FILTER (WHERE gwas_loci_by_overlap > 0) AS gwas_overlap,
  count(*) FILTER (WHERE gwas_loci_by_mapped_gene > 0) AS gwas_mapped,
  count(*) FILTER (WHERE gwas_loci_by_mapped_gene > 0 AND symbol_n_gencode = 1) AS gwas_mapped_unique_symbol,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS loeuf_not_null,
  median(max_median_tpm) FILTER (WHERE max_median_tpm > 0) AS median_max_tpm,
  median(tau) FILTER (WHERE max_median_tpm > 0) AS median_tau,
  median(n_tissues_tpm_ge_1) FILTER (WHERE max_median_tpm > 0) AS median_tissues_ge_1
FROM mart
WHERE biotype = 'protein_coding';

SELECT
  'excluded_counts' AS section,
  cls,
  biotype,
  count(*) AS n
FROM mart
WHERE cls IN ('excluded_protein_coding', 'excluded_ig_tr_segment')
GROUP BY cls, biotype
ORDER BY cls, n DESC;
