-- ClinVar drug-response size, constraint, and expression.
-- Read-only. Does not write a mart.
-- DuckDB 2.x (lambda x: ...). Multiple statements, in order.
-- Paths are the warehouse parquet files.

CREATE OR REPLACE TEMP VIEW alleles AS
SELECT
  dedup_key,
  variant_key,
  n_rows,
  "AlleleID",
  "VariationID",
  "GeneID",
  "GeneSymbol",
  "Type",
  "ClinicalSignificance",
  "ReviewStatus",
  review_stars,
  "PhenotypeList",
  "OriginSimple",
  germline_somatic,
  germline_pathogenic,
  "Oncogenicity",
  "SomaticClinicalImpact",
  list_filter(
    list_transform(
      regexp_split_to_array(lower(coalesce("ClinicalSignificance", '')), '[|;/]'),
      lambda x: trim(x)
    ),
    lambda x: x <> ''
      AND x NOT IN ('.', '-', 'na', 'n/a', 'none', 'not provided')
  ) AS pieces,
  list_filter(
    list_distinct(
      list_transform(string_split(coalesce("GeneSymbol", ''), ';'), lambda x: trim(x))
    ),
    lambda x: x <> '' AND x NOT IN ('.', '-')
  ) AS tokens
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet');

-- A drug-response allele has the exact piece "drug response".
-- The mart's clinvar_drug_response_n uses the same piece rule and is not
-- exclusive of pathogenic, benign, VUS, or conflicting pieces.
CREATE OR REPLACE TEMP VIEW drug AS
SELECT *
FROM alleles
WHERE len(list_filter(pieces, lambda x: contains(x, 'drug response'))) > 0;

CREATE OR REPLACE TEMP VIEW drug_tokens AS
SELECT
  d.dedup_key,
  d.review_stars,
  d.germline_somatic,
  d.germline_pathogenic,
  d."PhenotypeList",
  unnest(d.tokens) AS symbol
FROM drug d;

-- Constraint and expression come from the gene mart (bare-ENSG joins already
-- done there). This symbol join is the only gene key. Symbols that match
-- more than one GENCODE row are excluded from the gene summaries below.
CREATE OR REPLACE TEMP VIEW symbol_mart AS
SELECT
  symbol,
  count(*) AS n_genes,
  count(DISTINCT ensg) AS n_ensg,
  any_value(ensg) AS ensg,
  any_value(biotype) AS biotype,
  any_value(symbol_n_gencode) AS symbol_n_gencode,
  any_value(loeuf) AS loeuf,
  any_value(loeuf_decile) AS loeuf_decile,
  any_value("lof.pLI") AS pli,
  any_value("mis.z_score") AS mis_z,
  any_value(max_tissue) AS max_tissue,
  any_value(max_median_tpm) AS max_median_tpm,
  any_value(n_tissues_tpm_ge_1) AS n_tissues_tpm_ge_1,
  any_value(tau) AS tau,
  any_value(clinvar_drug_response_n) AS clinvar_drug_response_n
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')
GROUP BY symbol;

CREATE OR REPLACE TEMP VIEW drug_genes AS
SELECT
  t.symbol,
  m.ensg,
  m.biotype,
  m.n_genes,
  m.loeuf,
  m.loeuf_decile,
  m.pli,
  m.mis_z,
  m.max_tissue,
  m.max_median_tpm,
  m.n_tissues_tpm_ge_1,
  m.tau,
  m.clinvar_drug_response_n,
  count(DISTINCT t.dedup_key) AS alleles,
  count(DISTINCT t.dedup_key) FILTER (WHERE t.germline_somatic = 'germline') AS alleles_germline,
  count(DISTINCT t.dedup_key) FILTER (WHERE t.germline_somatic = 'somatic') AS alleles_somatic,
  count(DISTINCT t.dedup_key) FILTER (WHERE t.germline_somatic = 'germline/somatic') AS alleles_mixed,
  count(DISTINCT t.dedup_key) FILTER (
    WHERE t.germline_somatic NOT IN ('germline', 'somatic', 'germline/somatic')
  ) AS alleles_other_origin,
  count(DISTINCT t.dedup_key) FILTER (WHERE t.review_stars >= 3) AS alleles_star_ge3,
  count(DISTINCT t.dedup_key) FILTER (WHERE t.review_stars = 4) AS alleles_star4
FROM drug_tokens t
JOIN symbol_mart m ON t.symbol = m.symbol
WHERE m.n_genes = 1
GROUP BY ALL;

-- Join audit. A non-zero disagreement or a symbol collision stops the gene claim.
SELECT
  'join_audit' AS section,
  (SELECT count(*) FROM drug) AS drug_alleles,
  (SELECT count(*) FROM drug d WHERE len(d.tokens) = 0) AS alleles_no_token,
  (SELECT count(*) FROM drug d
    WHERE len(d.tokens) > 0
      AND NOT EXISTS (
        SELECT 1 FROM drug_tokens t
        JOIN symbol_mart m ON t.symbol = m.symbol
        WHERE t.dedup_key = d.dedup_key
      )) AS alleles_no_matching_symbol,
  (SELECT count(DISTINCT symbol) FROM drug_tokens) AS symbols,
  (SELECT count(*) FROM (
      SELECT t.symbol
      FROM drug_tokens t
      ANTI JOIN symbol_mart m ON t.symbol = m.symbol
      GROUP BY t.symbol
    )) AS unmatched_symbols,
  (SELECT count(*) FROM (
      SELECT t.symbol
      FROM drug_tokens t
      JOIN symbol_mart m ON t.symbol = m.symbol
      WHERE m.n_genes > 1
      GROUP BY t.symbol
    )) AS colliding_symbols,
  (SELECT count(*) FROM drug_genes) AS genes_unique_symbol,
  (SELECT count(*) FROM drug_genes WHERE alleles <> clinvar_drug_response_n) AS count_disagreements,
  (SELECT sum(clinvar_drug_response_n) FROM drug_genes) AS sum_matched_gene_links,
  (SELECT sum(clinvar_drug_response_n)
    FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')) AS mart_column_sum,
  (SELECT count(*) FROM drug
    WHERE len(list_filter(pieces, lambda x: x = 'drug response'))
      <> len(list_filter(pieces, lambda x: contains(x, 'drug response')))
  ) AS piece_not_exact;

SELECT
  'allele_size' AS section,
  (SELECT count(*) FROM alleles) AS allele_denominator,
  (SELECT sum(n_rows) FROM alleles) AS source_rows_in_denominator,
  (SELECT count(*) FROM drug) AS drug_alleles,
  round((SELECT count(*) FROM drug)::DOUBLE / (SELECT count(*) FROM alleles), 6) AS drug_fraction,
  (SELECT count(DISTINCT "VariationID") FROM drug) AS distinct_variation_id,
  (SELECT count(DISTINCT "AlleleID") FROM drug) AS distinct_allele_id,
  (SELECT count(*) FROM drug WHERE n_rows > 1) AS collapsed_drug_alleles,
  (SELECT sum(n_rows) FROM drug) AS drug_source_rows,
  (SELECT count(*) FROM drug WHERE variant_key IS NULL) AS drug_alleles_without_variant_key,
  (SELECT count(*) FROM drug WHERE len(tokens) = 0) AS no_symbol_token,
  (SELECT count(*) FROM drug WHERE contains(coalesce("GeneSymbol", ''), '|')) AS gene_symbol_contains_pipe,
  (SELECT count(*) FROM drug WHERE len(tokens) > 1) AS multi_token_alleles,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')) AS mart_genes,
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')
    WHERE biotype = 'protein_coding') AS mart_protein_coding;

SELECT
  'class_split' AS section,
  count(*) AS drug_alleles,
  count(*) FILTER (WHERE len(list_filter(pieces, lambda x: x <> 'drug response')) = 0) AS only_drug_piece,
  count(*) FILTER (WHERE len(list_filter(pieces, lambda x: x <> 'drug response')) > 0) AS drug_plus_other_piece,
  count(*) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_flag,
  count(*) FILTER (WHERE "PhenotypeList" = 'Tramadol response') AS tramadol_exact,
  count(*) FILTER (WHERE "PhenotypeList" <> 'Tramadol response' OR "PhenotypeList" IS NULL) AS not_tramadol_exact,
  count(*) FILTER (WHERE review_stars = 0) AS star_0,
  count(*) FILTER (WHERE review_stars >= 2) AS star_ge2,
  count(*) FILTER (WHERE review_stars >= 3) AS star_ge3,
  count(*) FILTER (WHERE review_stars = 4) AS star_4,
  count(*) FILTER (WHERE "Type" = 'single nucleotide variant') AS snv,
  count(*) FILTER (
    WHERE "Oncogenicity" IS NULL
      OR lower(trim("Oncogenicity")) IN ('', '-', '.', 'na', 'n/a', 'none', 'not provided')
  ) AS oncogenicity_empty,
  count(*) FILTER (WHERE "Oncogenicity" = 'Oncogenic') AS oncogenic,
  count(*) FILTER (WHERE "Oncogenicity" = 'Likely oncogenic') AS likely_oncogenic,
  count(*) FILTER (
    WHERE "SomaticClinicalImpact" IS NULL
      OR lower(trim("SomaticClinicalImpact")) IN ('', '-', '.', 'na', 'n/a', 'none', 'not provided')
  ) AS somatic_impact_empty
FROM drug;

SELECT 'stars' AS section, review_stars, "ReviewStatus", count(*) AS alleles
FROM drug
GROUP BY review_stars, "ReviewStatus"
ORDER BY review_stars, alleles DESC;

SELECT 'origin' AS section, germline_somatic, review_stars, count(*) AS alleles
FROM drug
GROUP BY germline_somatic, review_stars
ORDER BY germline_somatic, review_stars;

SELECT 'significance' AS section, "ClinicalSignificance", count(*) AS alleles
FROM drug
GROUP BY "ClinicalSignificance"
ORDER BY alleles DESC, "ClinicalSignificance";

SELECT
  'phenotype_ge10' AS section,
  "PhenotypeList",
  count(*) AS alleles,
  count(*) FILTER (WHERE review_stars = 0) AS star_0,
  count(*) FILTER (WHERE germline_somatic = 'somatic') AS somatic,
  count(*) FILTER (WHERE germline_somatic = 'germline') AS germline
FROM drug
GROUP BY "PhenotypeList"
HAVING count(*) >= 10
ORDER BY alleles DESC;

SELECT 'tramadol_genes' AS section, symbol, count(DISTINCT dedup_key) AS alleles
FROM drug_tokens
WHERE "PhenotypeList" = 'Tramadol response'
GROUP BY symbol
ORDER BY alleles DESC;

SELECT 'multi_token' AS section, "GeneSymbol", count(*) AS alleles
FROM drug
WHERE len(tokens) > 1
GROUP BY "GeneSymbol"
ORDER BY alleles DESC, "GeneSymbol";

SELECT 'star4_genes' AS section, "GeneSymbol", germline_somatic, count(*) AS alleles
FROM drug
WHERE review_stars = 4
GROUP BY "GeneSymbol", germline_somatic
ORDER BY alleles DESC;

SELECT
  'mixed_star_ge3' AS section,
  "GeneSymbol",
  "ClinicalSignificance",
  "PhenotypeList",
  review_stars,
  count(*) AS alleles
FROM drug
WHERE germline_somatic = 'germline/somatic' AND review_stars >= 3
GROUP BY "GeneSymbol", "ClinicalSignificance", "PhenotypeList", review_stars
ORDER BY alleles DESC, "GeneSymbol";

SELECT 'unmatched_tokens' AS section, symbol, count(DISTINCT dedup_key) AS alleles
FROM drug_tokens t
ANTI JOIN symbol_mart m ON t.symbol = m.symbol
GROUP BY symbol
ORDER BY alleles DESC, symbol;

SELECT
  'no_token_alleles' AS section,
  "GeneSymbol",
  "GeneID",
  "ClinicalSignificance",
  germline_somatic,
  review_stars,
  count(*) AS alleles
FROM drug
WHERE len(tokens) = 0
GROUP BY "GeneSymbol", "GeneID", "ClinicalSignificance", germline_somatic, review_stars;

-- Protein-coding comparison. Cuts are descriptive. They do not define the set.
-- Denominators are genes with a non-null value of that metric.
CREATE OR REPLACE TEMP VIEW coding_background AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')
WHERE biotype = 'protein_coding';

SELECT
  'constraint_expression' AS section,
  'drug_protein_coding' AS set_name,
  count(*) AS genes,
  count(*) FILTER (WHERE loeuf IS NOT NULL) AS with_loeuf,
  round(median(loeuf), 3) AS median_loeuf,
  round(quantile_cont(loeuf, 0.25), 3) AS loeuf_p25,
  round(quantile_cont(loeuf, 0.75), 3) AS loeuf_p75,
  count(*) FILTER (WHERE loeuf < 0.35) AS loeuf_lt_0_35,
  count(*) FILTER (WHERE loeuf < 0.6) AS loeuf_lt_0_6,
  count(*) FILTER (WHERE loeuf_decile IS NOT NULL) AS with_decile,
  count(*) FILTER (WHERE loeuf_decile = 0) AS decile_0,
  count(*) FILTER (WHERE loeuf_decile <= 2) AS decile_0_to_2,
  count(*) FILTER (WHERE pli IS NOT NULL) AS with_pli,
  round(median(pli), 3) AS median_pli,
  count(*) FILTER (WHERE pli >= 0.9) AS pli_ge_0_9,
  count(*) FILTER (WHERE mis_z IS NOT NULL) AS with_mis_z,
  round(median(mis_z), 3) AS median_mis_z,
  count(*) FILTER (WHERE mis_z >= 3.09) AS mis_z_ge_3_09,
  count(*) FILTER (WHERE tau IS NOT NULL) AS with_tau,
  round(median(tau), 3) AS median_tau,
  round(quantile_cont(tau, 0.25), 3) AS tau_p25,
  round(quantile_cont(tau, 0.75), 3) AS tau_p75,
  count(*) FILTER (WHERE tau >= 0.8) AS tau_ge_0_8,
  round(median(n_tissues_tpm_ge_1), 1) AS median_tissues_ge1,
  round(median(max_median_tpm), 2) AS median_max_tpm
FROM drug_genes
WHERE biotype = 'protein_coding'
UNION ALL
SELECT
  'constraint_expression',
  'all_protein_coding',
  count(*),
  count(*) FILTER (WHERE loeuf IS NOT NULL),
  round(median(loeuf), 3),
  round(quantile_cont(loeuf, 0.25), 3),
  round(quantile_cont(loeuf, 0.75), 3),
  count(*) FILTER (WHERE loeuf < 0.35),
  count(*) FILTER (WHERE loeuf < 0.6),
  count(*) FILTER (WHERE loeuf_decile IS NOT NULL),
  count(*) FILTER (WHERE loeuf_decile = 0),
  count(*) FILTER (WHERE loeuf_decile <= 2),
  count(*) FILTER (WHERE "lof.pLI" IS NOT NULL),
  round(median("lof.pLI"), 3),
  count(*) FILTER (WHERE "lof.pLI" >= 0.9),
  count(*) FILTER (WHERE "mis.z_score" IS NOT NULL),
  round(median("mis.z_score"), 3),
  count(*) FILTER (WHERE "mis.z_score" >= 3.09),
  count(*) FILTER (WHERE tau IS NOT NULL),
  round(median(tau), 3),
  round(quantile_cont(tau, 0.25), 3),
  round(quantile_cont(tau, 0.75), 3),
  count(*) FILTER (WHERE tau >= 0.8),
  round(median(n_tissues_tpm_ge_1), 1),
  round(median(max_median_tpm), 2)
FROM coding_background
UNION ALL
SELECT
  'constraint_expression',
  'drug_germline_protein_coding',
  count(*),
  count(*) FILTER (WHERE loeuf IS NOT NULL),
  round(median(loeuf), 3),
  round(quantile_cont(loeuf, 0.25), 3),
  round(quantile_cont(loeuf, 0.75), 3),
  count(*) FILTER (WHERE loeuf < 0.35),
  count(*) FILTER (WHERE loeuf < 0.6),
  count(*) FILTER (WHERE loeuf_decile IS NOT NULL),
  count(*) FILTER (WHERE loeuf_decile = 0),
  count(*) FILTER (WHERE loeuf_decile <= 2),
  count(*) FILTER (WHERE pli IS NOT NULL),
  round(median(pli), 3),
  count(*) FILTER (WHERE pli >= 0.9),
  count(*) FILTER (WHERE mis_z IS NOT NULL),
  round(median(mis_z), 3),
  count(*) FILTER (WHERE mis_z >= 3.09),
  count(*) FILTER (WHERE tau IS NOT NULL),
  round(median(tau), 3),
  round(quantile_cont(tau, 0.25), 3),
  round(quantile_cont(tau, 0.75), 3),
  count(*) FILTER (WHERE tau >= 0.8),
  round(median(n_tissues_tpm_ge_1), 1),
  round(median(max_median_tpm), 2)
FROM drug_genes
WHERE biotype = 'protein_coding' AND alleles_germline > 0;

SELECT
  'genes_ge10' AS section,
  symbol,
  biotype,
  alleles,
  alleles_germline,
  alleles_somatic,
  alleles_mixed,
  alleles_other_origin,
  alleles_star_ge3,
  alleles_star4,
  round(loeuf, 3) AS loeuf,
  loeuf_decile,
  round(pli, 3) AS pli,
  round(mis_z, 2) AS mis_z,
  max_tissue,
  round(max_median_tpm, 2) AS max_median_tpm,
  n_tissues_tpm_ge_1,
  round(tau, 3) AS tau
FROM drug_genes
WHERE alleles >= 10
ORDER BY alleles DESC, symbol;

SELECT
  'gene_span' AS section,
  count(*) AS matched_genes,
  count(*) FILTER (WHERE biotype = 'protein_coding') AS protein_coding_genes,
  count(*) FILTER (WHERE alleles = 1) AS genes_with_one_allele,
  sum(alleles) AS gene_links,
  count(*) FILTER (WHERE loeuf IS NULL) AS loeuf_null,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND loeuf IS NULL) AS coding_loeuf_null,
  count(*) FILTER (WHERE tau IS NULL) AS tau_null
FROM drug_genes;

SELECT
  'null_metric_genes' AS section,
  g.symbol,
  g.ensg,
  g.biotype,
  g.alleles,
  g.loeuf IS NULL AS loeuf_null,
  g.tau IS NULL AS tau_null,
  c.bare_ensg IS NOT NULL AS constraint_row,
  c.loeuf AS constraint_loeuf,
  e.ensg IS NOT NULL AS expression_row
FROM drug_genes g
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/constraint/gene_constraint.parquet') c
  ON g.ensg = c.bare_ensg
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/gtex/gene_expression.parquet') e
  ON g.ensg = e.ensg
WHERE g.loeuf IS NULL OR g.tau IS NULL
ORDER BY g.symbol;

SELECT
  'biotype' AS section,
  biotype,
  count(*) AS genes,
  sum(alleles) AS gene_links
FROM drug_genes
GROUP BY biotype
ORDER BY genes DESC, biotype;
