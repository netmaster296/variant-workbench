-- Tissue programs. Read-only. Does not write a mart.
-- Numbers for warehouse/reports/tissue_programs.md.
-- Tau is gene_mart.tau. This file does not recompute it.
-- Expression tissues are GTEx v8 median-TPM column names ("Whole Blood").
-- QTL tissue names ("Whole_Blood") are not read and are not treated as the same strings.
-- Run in a fresh DuckDB connection. Tables are temporary.

CREATE OR REPLACE TEMP TABLE group_rule AS
SELECT * FROM (VALUES
  ('blood_cell',
    '\b(?:platelets?|erythrocytes?|leukocytes?|lymphocytes?|monocytes?|neutrophils?|eosinophils?|basophils?|granulocytes?|reticulocytes?|red blood cells?|white blood cells?|corpuscular\w*|h[ae]matocrit\w*|h[ae]matolog\w*|hemoglobins?|haemoglobins?)\b',
    'hemoglobin a1|haemoglobin a1|hba1c|glycated hemoglobin|glycated haemoglobin|amount of |protein|antigen|glycoprotein|molecule| in blood| in serum| in cerebrospinal|folate'),
  ('brain',
    '\b(?:brain(?:stem)?|alzheimer\w*|parkinson\w*|schizophreni\w*|hippocamp\w*|cerebell(?:um|ar)\w*|amygdala\w*|hypothalam\w*|substantia nigra|intracerebral|cerebral)\b',
    'neurotrophic|amount of | in blood|protein measurement| in serum| in cerebrospinal|protease'),
  ('liver',
    '\bliver\b|hepatic|\b(?:alanine aminotransferase|aspartate aminotransferase|bilirubin|gamma-glutamyl transferase|gamma-glutamyltransferase)\b',
    'amount of | in blood|alkaline phosphatase|protein'),
  ('heart',
    '\b(?:heart|cardiac|atrial|ventricular|myocard\w*|cardiomyopath\w*)\b',
    'amount of |lateral ventricle|cerebral ventricle|protein|use measurement'),
  ('lung',
    '\b(?:lung|pulmonary|asthma|chronic obstructive)\b|\bfev\b|\bfvc\b|\bforced expiratory\b|\bforced vital\b',
    'amount of | in blood|surfactant-associated protein|scarlet fever'),
  ('kidney',
    '\b(?:kidney|renal|glomerul\w*)\b|\bserum creatinine\b|\bcreatinine clearance\b|\balbumin to creatinine\b|\bglomerular filtration\b',
    'amount of | in blood|hyperparathyroid|molecule|amyloid|adipose'),
  ('adipose',
    '\badipose\b|\bbody fat\b',
    'amount of |cardiac enriched|renal sinus'),
  ('thyroid',
    'thyroid|\b(?:thyroxin\w*|triiodothyronine|tri-iodothyronine|thyrotropin|hashimoto\w*)\b',
    'parathyroid| in blood|binding globulin|preparation use|anti-thyroid drug'),
  ('skin',
    '\b(?:skin|psoriasis|dermatitis|eczema|cutaneous)\b',
    'dyskinesia|uveal|choroidal|amount of | in blood|microbiome|tuberculin'),
  ('pancreas',
    '\bpancrea',
    'amount of | in blood| in serum| in cerebrospinal|lipase-related|prohormone|ribonuclease|trypsin inhibitor'),
  ('breast',
    '\bbreast\b',
    'breastfeeding|protein'),
  ('prostate',
    '\bprostate\b',
    'amount of | in blood| in cerebrospinal|expressed protein|microseminoprotein'),
  ('ovary',
    '\bovar(?:y|ian)\b',
    ''),
  ('gut',
    '\b(?:colon|colonic|colorectal|crohn|ulcerative colitis|inflammatory bowel)\b',
    'colony|amount of | in blood'),
  ('artery',
    '\b(?:artery|arterial|atherosclerosis)\b',
    ''),
  ('skeletal_muscle',
    '\b(?:skeletal muscle|thigh muscle)\b',
    'troponin|amount of |fat infiltration')
) AS v(grp, include_re, exclude_re);

-- Median-TPM names only. One row per index tissue.
CREATE OR REPLACE TEMP TABLE tissue_map AS
SELECT * FROM (VALUES
  ('blood_cell', 'Whole Blood'),
  ('brain', 'Brain - Amygdala'),
  ('brain', 'Brain - Anterior cingulate cortex (BA24)'),
  ('brain', 'Brain - Caudate (basal ganglia)'),
  ('brain', 'Brain - Cerebellar Hemisphere'),
  ('brain', 'Brain - Cerebellum'),
  ('brain', 'Brain - Cortex'),
  ('brain', 'Brain - Frontal Cortex (BA9)'),
  ('brain', 'Brain - Hippocampus'),
  ('brain', 'Brain - Hypothalamus'),
  ('brain', 'Brain - Nucleus accumbens (basal ganglia)'),
  ('brain', 'Brain - Putamen (basal ganglia)'),
  ('brain', 'Brain - Spinal cord (cervical c-1)'),
  ('brain', 'Brain - Substantia nigra'),
  ('liver', 'Liver'),
  ('heart', 'Heart - Atrial Appendage'),
  ('heart', 'Heart - Left Ventricle'),
  ('lung', 'Lung'),
  ('kidney', 'Kidney - Cortex'),
  ('kidney', 'Kidney - Medulla'),
  ('adipose', 'Adipose - Subcutaneous'),
  ('adipose', 'Adipose - Visceral (Omentum)'),
  ('thyroid', 'Thyroid'),
  ('skin', 'Skin - Not Sun Exposed (Suprapubic)'),
  ('skin', 'Skin - Sun Exposed (Lower leg)'),
  ('pancreas', 'Pancreas'),
  ('breast', 'Breast - Mammary Tissue'),
  ('prostate', 'Prostate'),
  ('ovary', 'Ovary'),
  ('gut', 'Colon - Sigmoid'),
  ('gut', 'Colon - Transverse'),
  ('gut', 'Small Intestine - Terminal Ileum'),
  ('artery', 'Artery - Aorta'),
  ('artery', 'Artery - Coronary'),
  ('artery', 'Artery - Tibial'),
  ('skeletal_muscle', 'Muscle - Skeletal')
) AS v(grp, tissue);

CREATE OR REPLACE TEMP TABLE trait_dim AS
SELECT
  tr.trait_uri,
  tr.mapped_trait,
  lower(tr.mapped_trait) AS name,
  count(*)::BIGINT AS n_loci
FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') l
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
  ON l.trait_uri IS NOT DISTINCT FROM tr.trait_uri
GROUP BY 1, 2, 3;

CREATE OR REPLACE TEMP TABLE hits AS
SELECT d.trait_uri, d.mapped_trait, d.n_loci, r.grp
FROM trait_dim d
JOIN group_rule r
  ON d.name IS NOT NULL
 AND d.name <> ''
 AND regexp_matches(d.name, r.include_re)
 AND (r.exclude_re = '' OR NOT regexp_matches(d.name, r.exclude_re));

CREATE OR REPLACE TEMP TABLE single_trait AS
SELECT h.grp, h.trait_uri, h.mapped_trait, h.n_loci
FROM hits h
WHERE h.trait_uri IN (
  SELECT trait_uri FROM hits GROUP BY 1 HAVING count(*) = 1
);

CREATE OR REPLACE TEMP TABLE ambiguous_trait AS
SELECT
  h.mapped_trait,
  any_value(h.n_loci) AS n_loci,
  string_agg(h.grp, ', ' ORDER BY h.grp) AS groups
FROM hits h
WHERE h.trait_uri IN (
  SELECT trait_uri FROM hits GROUP BY 1 HAVING count(*) > 1
)
GROUP BY h.mapped_trait, h.trait_uri;

CREATE OR REPLACE TEMP TABLE assigned AS
SELECT s.grp, s.trait_uri, s.mapped_trait, l.mapped_genes
FROM single_trait s
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') l
  ON l.trait_uri IS NOT DISTINCT FROM s.trait_uri;

CREATE OR REPLACE TEMP TABLE symbol_n AS
SELECT symbol, count(*)::BIGINT AS n_genes
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')
GROUP BY symbol;

CREATE OR REPLACE TEMP TABLE token_rows AS
SELECT a.grp, unnest(a.mapped_genes.tokens) AS token
FROM assigned a
WHERE a.mapped_genes IS NOT NULL
  AND a.mapped_genes.tokens IS NOT NULL;

CREATE OR REPLACE TEMP TABLE token_class AS
SELECT
  tk.grp,
  tk.token,
  CASE
    WHEN s.n_genes = 1 THEN 'unique'
    WHEN s.n_genes > 1 THEN 'ambiguous'
    WHEN tk.token LIKE '% - %' THEN 'pair'
    ELSE 'other_unmatched'
  END AS kind
FROM (
  SELECT DISTINCT grp, token FROM token_rows
) tk
LEFT JOIN symbol_n s ON tk.token = s.symbol;

CREATE OR REPLACE TEMP TABLE genes AS
SELECT
  tk.grp,
  g.ensg,
  g.symbol,
  g.biotype,
  g.tau,
  g.max_tissue,
  g.max_median_tpm,
  g.n_tissues_tpm_ge_1,
  CASE
    WHEN g.max_median_tpm IS NULL THEN 'no_gtex_row'
    WHEN g.tau IS NULL THEN 'tau_undefined'
    WHEN g.tau < 0.5 THEN 'broad'
    WHEN g.tau >= 0.95 THEN 'private'
    ELSE 'intermediate'
  END AS tau_class
FROM token_class tk
JOIN read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet') g
  ON g.symbol = tk.token
WHERE tk.kind = 'unique';

CREATE OR REPLACE TEMP TABLE qc_check AS
WITH mart AS (
  SELECT * FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet')
),
expr AS (
  SELECT * FROM read_parquet('D:/variant-workbench/warehouse/gtex/gene_expression.parquet')
),
med AS (
  SELECT * FROM read_parquet('D:/variant-workbench/warehouse/gtex/median_tpm.parquet')
),
med_ensg AS (
  SELECT ensg, count(*)::BIGINT AS n_rows, count(DISTINCT tissue)::BIGINT AS n_tissues
  FROM med
  WHERE ensg IS NOT NULL
  GROUP BY ensg
),
checks AS (
  SELECT 'mart_rows' AS name, 78691::DOUBLE AS expected,
    (SELECT count(*) FROM mart)::DOUBLE AS actual
  UNION ALL
  SELECT 'mart_ensg_unique', 78691,
    (SELECT count(DISTINCT ensg) FROM mart WHERE ensg IS NOT NULL)
  UNION ALL
  SELECT 'mart_ensg_null', 0,
    (SELECT count(*) FROM mart WHERE ensg IS NULL)
  UNION ALL
  SELECT 'no_gtex_row_on_mart', 24065,
    (SELECT count(*) FROM mart WHERE max_median_tpm IS NULL)
  UNION ALL
  SELECT 'no_median_ensg_match', 24065,
    (SELECT count(*) FROM mart g ANTI JOIN med_ensg e ON g.ensg = e.ensg)
  UNION ALL
  SELECT 'expressed_genes_not_54_tissues', 0,
    (SELECT count(*) FROM mart g JOIN med_ensg e ON g.ensg = e.ensg WHERE e.n_rows <> 54 OR e.n_tissues <> 54)
  UNION ALL
  SELECT 'median_ensg_tissue_dups', 0,
    (SELECT count(*) FROM (
      SELECT ensg, tissue FROM med WHERE ensg IS NOT NULL
      GROUP BY 1, 2 HAVING count(*) > 1
    ))
  UNION ALL
  SELECT 'tau_mart_vs_expression_mismatch', 0,
    (SELECT count(*) FROM mart g JOIN expr e ON g.ensg = e.ensg
     WHERE g.tau IS DISTINCT FROM e.tau)
  UNION ALL
  SELECT 'tau_set_when_max_tpm_null', 0,
    (SELECT count(*) FROM mart WHERE max_median_tpm IS NULL AND tau IS NOT NULL)
  UNION ALL
  SELECT 'tau_set_when_max_tpm_zero', 0,
    (SELECT count(*) FROM mart WHERE max_median_tpm = 0 AND tau IS NOT NULL)
  UNION ALL
  SELECT 'map_tissue_missing_from_median', 0,
    (SELECT count(*) FROM tissue_map m
     ANTI JOIN (SELECT DISTINCT tissue FROM med) t ON m.tissue = t.tissue)
  UNION ALL
  SELECT 'map_tissue_with_underscore', 0,
    (SELECT count(*) FROM tissue_map WHERE position('_' IN tissue) > 0)
  UNION ALL
  SELECT 'map_group_missing', 0,
    (SELECT count(*) FROM group_rule r
     ANTI JOIN (SELECT DISTINCT grp FROM tissue_map) m ON r.grp = m.grp)
  UNION ALL
  SELECT 'loci_without_trait_row', 0,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') l
     ANTI JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
       ON l.trait_uri IS NOT DISTINCT FROM tr.trait_uri)
  UNION ALL
  SELECT 'genes_grp_ensg_dups', 0,
    (SELECT count(*) FROM (
      SELECT grp, ensg FROM genes GROUP BY 1, 2 HAVING count(*) > 1
    ))
  UNION ALL
  SELECT 'unique_token_not_one_gene', 0,
    (SELECT count(*) FROM token_class tk
     JOIN symbol_n s ON tk.token = s.symbol
     WHERE tk.kind = 'unique' AND s.n_genes <> 1)
  UNION ALL
  SELECT 'median_tissue_count', 54,
    (SELECT count(DISTINCT tissue) FROM med)
)
SELECT name, expected, actual, actual = expected AS ok FROM checks;

CREATE OR REPLACE TEMP TABLE funnel AS
SELECT
  s.grp,
  count(DISTINCT s.trait_uri)::BIGINT AS n_traits,
  sum(s.n_loci)::BIGINT AS n_loci,
  (SELECT count(*) FROM token_class t WHERE t.grp = s.grp) AS n_tokens,
  (SELECT count(*) FROM token_class t WHERE t.grp = s.grp AND t.kind = 'pair') AS n_pair_tokens,
  (SELECT count(*) FROM token_class t WHERE t.grp = s.grp AND t.kind = 'other_unmatched') AS n_other_unmatched,
  (SELECT count(*) FROM token_class t WHERE t.grp = s.grp AND t.kind = 'ambiguous') AS n_ambiguous_symbols,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp) AS n_genes,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding') AS n_protein_coding,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.tau_class = 'no_gtex_row') AS n_no_gtex,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.tau_class = 'no_gtex_row') AS n_no_gtex_pc,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.tau_class = 'tau_undefined') AS n_tau_undefined,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.tau_class = 'tau_undefined') AS n_tau_undefined_pc,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.tau_class = 'broad') AS n_broad,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.tau_class = 'broad') AS n_broad_pc,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.tau_class = 'intermediate') AS n_intermediate,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.tau_class = 'intermediate') AS n_intermediate_pc,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.tau_class = 'private') AS n_private,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.tau_class = 'private') AS n_private_pc,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'lncRNA') AS n_lncrna,
  (SELECT count(*) FROM genes g WHERE g.grp = s.grp AND g.biotype = 'protein_coding' AND g.max_median_tpm IS NOT NULL) AS n_pc_with_gtex
FROM single_trait s
GROUP BY s.grp;

CREATE OR REPLACE TEMP TABLE class_sum_check AS
SELECT
  grp,
  n_genes = n_no_gtex + n_tau_undefined + n_broad + n_intermediate + n_private AS all_ok,
  n_protein_coding = n_no_gtex_pc + n_tau_undefined_pc + n_broad_pc + n_intermediate_pc + n_private_pc AS pc_ok
FROM funnel;

CREATE OR REPLACE TEMP TABLE genome_tau AS
SELECT
  count(*)::BIGINT AS n_genes,
  count(*) FILTER (WHERE biotype = 'protein_coding')::BIGINT AS n_pc,
  count(*) FILTER (WHERE max_median_tpm IS NULL)::BIGINT AS n_no_gtex,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND max_median_tpm IS NULL)::BIGINT AS n_no_gtex_pc,
  count(*) FILTER (WHERE max_median_tpm IS NOT NULL AND tau IS NULL)::BIGINT AS n_tau_undefined,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND max_median_tpm IS NOT NULL AND tau IS NULL)::BIGINT AS n_tau_undefined_pc,
  count(*) FILTER (WHERE tau < 0.5)::BIGINT AS n_broad,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND tau < 0.5)::BIGINT AS n_broad_pc,
  count(*) FILTER (WHERE tau >= 0.5 AND tau < 0.95)::BIGINT AS n_intermediate,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND tau >= 0.5 AND tau < 0.95)::BIGINT AS n_intermediate_pc,
  count(*) FILTER (WHERE tau >= 0.95)::BIGINT AS n_private,
  count(*) FILTER (WHERE biotype = 'protein_coding' AND tau >= 0.95)::BIGINT AS n_private_pc,
  count(*) FILTER (WHERE tau < 0.3)::BIGINT AS n_tau_lt_0_3,
  min(tau) AS tau_min,
  max(tau) AS tau_max
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet');

CREATE OR REPLACE TEMP TABLE tissue_stat AS
SELECT
  g.grp,
  t.tissue,
  max(CASE WHEN m.tissue IS NOT NULL THEN 1 ELSE 0 END) AS is_index,
  count(*)::BIGINT AS n_all,
  count(*) FILTER (WHERE t.median_tpm >= 1)::BIGINT AS n_all_ge1,
  count(*) FILTER (WHERE g.biotype = 'protein_coding')::BIGINT AS n_pc,
  count(*) FILTER (WHERE g.biotype = 'protein_coding' AND t.median_tpm >= 1)::BIGINT AS n_pc_ge1,
  median(t.median_tpm) FILTER (WHERE g.biotype = 'protein_coding') AS median_tpm_pc
FROM genes g
JOIN read_parquet('D:/variant-workbench/warehouse/gtex/median_tpm.parquet') t
  ON g.ensg = t.ensg
LEFT JOIN tissue_map m
  ON m.grp = g.grp AND m.tissue = t.tissue
GROUP BY g.grp, t.tissue;

CREATE OR REPLACE TEMP TABLE tissue_balance AS
SELECT
  grp,
  min(n_all) AS min_all,
  max(n_all) AS max_all,
  min(n_pc) AS min_pc,
  max(n_pc) AS max_pc,
  count(*)::BIGINT AS n_tissues
FROM tissue_stat
GROUP BY grp;

CREATE OR REPLACE TEMP TABLE bg_pc AS
SELECT
  t.tissue,
  count(*)::BIGINT AS n_pc,
  count(*) FILTER (WHERE t.median_tpm >= 1)::BIGINT AS n_ge1
FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet') g
JOIN read_parquet('D:/variant-workbench/warehouse/gtex/median_tpm.parquet') t
  ON g.ensg = t.ensg
WHERE g.biotype = 'protein_coding'
  AND g.max_median_tpm IS NOT NULL
GROUP BY t.tissue;

CREATE OR REPLACE TEMP TABLE tissue_ranked AS
SELECT
  s.grp,
  s.tissue,
  s.is_index,
  s.n_all,
  s.n_all_ge1,
  s.n_pc,
  s.n_pc_ge1,
  s.n_pc_ge1::DOUBLE / nullif(s.n_pc, 0) AS frac_pc,
  s.n_all_ge1::DOUBLE / nullif(s.n_all, 0) AS frac_all,
  s.median_tpm_pc,
  b.n_pc AS bg_n_pc,
  b.n_ge1 AS bg_n_ge1,
  b.n_ge1::DOUBLE / nullif(b.n_pc, 0) AS bg_frac_pc,
  rank() OVER (PARTITION BY s.grp ORDER BY s.n_pc_ge1::DOUBLE / nullif(s.n_pc, 0) DESC, s.tissue) AS rank_pc,
  row_number() OVER (PARTITION BY s.grp ORDER BY s.n_pc_ge1::DOUBLE / nullif(s.n_pc, 0) DESC, s.tissue) AS rn_pc
FROM tissue_stat s
JOIN bg_pc b ON s.tissue = b.tissue;

CREATE OR REPLACE TEMP TABLE any_index AS
SELECT
  g.grp,
  count(*)::BIGINT AS n_pc,
  count(*) FILTER (WHERE x.hit)::BIGINT AS n_pc_any_index_ge1,
  count(*) FILTER (WHERE g.tau_class = 'broad')::BIGINT AS n_broad_pc,
  count(*) FILTER (WHERE g.tau_class = 'broad' AND x.hit)::BIGINT AS n_broad_pc_any_index_ge1,
  count(*) FILTER (WHERE g.tau_class = 'private')::BIGINT AS n_private_pc,
  count(*) FILTER (WHERE g.tau_class = 'private' AND g.max_tissue IN (SELECT tissue FROM tissue_map m WHERE m.grp = g.grp))::BIGINT AS n_private_pc_max_in_index,
  median(g.n_tissues_tpm_ge_1) FILTER (WHERE g.tau_class = 'broad') AS broad_median_n_tissues_ge1,
  median(g.max_median_tpm) FILTER (WHERE g.tau_class = 'broad') AS broad_median_max_tpm,
  median(g.max_median_tpm) FILTER (WHERE g.tau_class = 'private') AS private_median_max_tpm,
  count(*) FILTER (WHERE g.tau_class = 'private' AND g.max_median_tpm >= 1)::BIGINT AS n_private_pc_max_ge1,
  count(*) FILTER (WHERE g.tau_class = 'private' AND g.max_median_tpm >= 10)::BIGINT AS n_private_pc_max_ge10,
  count(*) FILTER (WHERE g.tau_class = 'private' AND g.max_tissue = 'Testis')::BIGINT AS n_private_pc_max_testis
FROM (
  SELECT *
  FROM genes
  WHERE biotype = 'protein_coding' AND max_median_tpm IS NOT NULL
) g
LEFT JOIN (
  SELECT g2.grp, g2.ensg, bool_or(t.median_tpm >= 1) AS hit
  FROM genes g2
  JOIN tissue_map m ON m.grp = g2.grp
  JOIN read_parquet('D:/variant-workbench/warehouse/gtex/median_tpm.parquet') t
    ON t.ensg = g2.ensg AND t.tissue = m.tissue
  WHERE g2.biotype = 'protein_coding' AND g2.max_median_tpm IS NOT NULL
  GROUP BY g2.grp, g2.ensg
) x ON x.grp = g.grp AND x.ensg = g.ensg
GROUP BY g.grp;

CREATE OR REPLACE TEMP TABLE private_max_top AS
SELECT grp, max_tissue, n, rn
FROM (
  SELECT
    grp,
    max_tissue,
    count(*)::BIGINT AS n,
    row_number() OVER (PARTITION BY grp ORDER BY count(*) DESC, max_tissue) AS rn
  FROM genes
  WHERE biotype = 'protein_coding' AND tau_class = 'private'
  GROUP BY grp, max_tissue
)
WHERE rn <= 5;

CREATE OR REPLACE TEMP TABLE multi_group_genes AS
SELECT
  count(*)::BIGINT AS n_ensg_in_2plus,
  count(*) FILTER (WHERE biotype = 'protein_coding')::BIGINT AS n_pc_in_2plus
FROM (
  SELECT ensg, any_value(biotype) AS biotype, count(DISTINCT grp) AS n_grp
  FROM genes
  GROUP BY ensg
  HAVING count(DISTINCT grp) > 1
);

CREATE OR REPLACE TEMP TABLE case_mismatch AS
SELECT count(*)::BIGINT AS n_tokens
FROM (
  SELECT DISTINCT token FROM token_rows
) tk
LEFT JOIN symbol_n s ON tk.token = s.symbol
JOIN symbol_n s2 ON lower(tk.token) = lower(s2.symbol)
WHERE s.symbol IS NULL
  AND tk.token NOT LIKE '% - %';

CREATE OR REPLACE TEMP TABLE other_unmatched_top AS
SELECT token, count(*)::BIGINT AS n_locus_rows
FROM token_rows tk
LEFT JOIN symbol_n s ON tk.token = s.symbol
WHERE s.symbol IS NULL
  AND tk.token NOT LIKE '% - %'
GROUP BY token
ORDER BY n_locus_rows DESC, token
LIMIT 15;

CREATE OR REPLACE TEMP TABLE coverage AS
SELECT
  (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet'))::BIGINT AS n_locus_rows,
  (SELECT count(*) FROM trait_dim)::BIGINT AS n_traits_with_loci,
  (SELECT coalesce(sum(n_loci), 0) FROM single_trait)::BIGINT AS n_assigned_loci,
  (SELECT count(*) FROM single_trait)::BIGINT AS n_assigned_traits,
  (SELECT coalesce(sum(n_loci), 0) FROM (
     SELECT trait_uri, any_value(n_loci) AS n_loci FROM hits
     WHERE trait_uri IN (SELECT trait_uri FROM hits GROUP BY 1 HAVING count(*) > 1)
     GROUP BY trait_uri
  ))::BIGINT AS n_ambiguous_loci,
  (SELECT count(DISTINCT trait_uri) FROM hits
   WHERE trait_uri IN (SELECT trait_uri FROM hits GROUP BY 1 HAVING count(*) > 1))::BIGINT AS n_ambiguous_traits;

CREATE OR REPLACE TEMP TABLE left_out AS
SELECT * FROM (VALUES
  ('bone_density_or_fracture', '\bbone mineral\b|\bbone tissue density\b|\bbone fracture\b|\bheel bone\b'),
  ('cholesterol_or_triglyceride', '\bcholesterol\b|\btriglyceride\b'),
  ('blood_pressure', '\bblood pressure\b'),
  ('depression', '\bdepress')
) AS v(label, re);

CREATE OR REPLACE TEMP TABLE left_out_counts AS
SELECT
  x.label,
  count(*)::BIGINT AS n_traits,
  coalesce(sum(d.n_loci), 0)::BIGINT AS n_loci
FROM left_out x
JOIN trait_dim d ON regexp_matches(d.name, x.re)
WHERE d.trait_uri NOT IN (SELECT trait_uri FROM single_trait)
  AND d.trait_uri NOT IN (
    SELECT trait_uri FROM hits GROUP BY 1 HAVING count(*) > 1
  )
GROUP BY x.label;

CREATE OR REPLACE TEMP TABLE lab_share AS
SELECT * FROM (VALUES
  ('heart', 'name contains troponin', 'troponin'),
  ('prostate', 'prostate specific antigen', 'prostate specific antigen'),
  ('thyroid', 'TSH, thyroxine, triiodothyronine, or thyrotropin', 'thyroxine|thyroid stimulating|triiodothyronine|tri-iodothyronine|thyrotropin'),
  ('liver', 'aminotransferase, bilirubin, or gamma-glutamyl transferase', 'aminotransferase|bilirubin|gamma-glutamyl transferase'),
  ('kidney', 'filtration or creatinine', 'glomerular filtration|serum creatinine|creatinine clearance|albumin to creatinine'),
  ('artery', 'arterial pressure', 'arterial pressure')
) AS v(grp, label, re);

CREATE OR REPLACE TEMP TABLE lab_share_counts AS
SELECT
  a.grp,
  a.label,
  count(DISTINCT s.trait_uri)::BIGINT AS n_traits,
  coalesce(sum(s.n_loci), 0)::BIGINT AS n_loci
FROM lab_share a
JOIN single_trait s
  ON s.grp = a.grp
 AND regexp_matches(lower(s.mapped_trait), a.re)
GROUP BY a.grp, a.label;

-- Excess is the protein-coding detection fraction minus the same fraction
-- for every protein-coding mart gene that has a GTEx row. It is not a test.
-- bg_rank is the tissue's rank in that background, 1 = highest fraction.
CREATE OR REPLACE TEMP TABLE bg_rank AS
SELECT
  tissue,
  rank() OVER (ORDER BY n_ge1::DOUBLE / nullif(n_pc, 0) DESC, tissue) AS bg_rank
FROM bg_pc;

CREATE OR REPLACE TEMP TABLE tissue_excess AS
SELECT
  s.grp,
  s.tissue,
  s.is_index,
  s.n_pc,
  s.n_pc_ge1,
  s.frac_pc,
  s.median_tpm_pc,
  s.bg_n_pc,
  s.bg_frac_pc,
  s.frac_pc - s.bg_frac_pc AS excess,
  s.rank_pc,
  rank() OVER (
    PARTITION BY s.grp
    ORDER BY (s.frac_pc - s.bg_frac_pc) DESC, s.tissue
  ) AS rank_excess,
  b.bg_rank
FROM tissue_ranked s
JOIN bg_rank b ON s.tissue = b.tissue;

CREATE OR REPLACE TEMP TABLE broad_span AS
SELECT
  grp,
  count(*) FILTER (
    WHERE biotype = 'protein_coding' AND tau_class = 'broad'
  )::BIGINT AS n_broad_pc,
  count(*) FILTER (
    WHERE biotype = 'protein_coding'
      AND tau_class = 'broad'
      AND n_tissues_tpm_ge_1 = 54
  )::BIGINT AS n_broad_pc_all_54
FROM genes
GROUP BY grp;

CREATE OR REPLACE TEMP TABLE gene_union AS
SELECT
  count(DISTINCT ensg)::BIGINT AS n_genes,
  count(DISTINCT ensg) FILTER (WHERE tau_class = 'no_gtex_row')::BIGINT AS n_no_gtex,
  count(DISTINCT ensg) FILTER (WHERE biotype = 'protein_coding')::BIGINT AS n_pc,
  count(DISTINCT ensg) FILTER (
    WHERE biotype = 'protein_coding' AND tau_class = 'no_gtex_row'
  )::BIGINT AS n_no_gtex_pc,
  count(DISTINCT ensg) FILTER (WHERE biotype = 'protein_coding' AND tau < 0.5)::BIGINT AS n_broad_pc,
  count(DISTINCT ensg) FILTER (
    WHERE biotype = 'protein_coding' AND tau >= 0.95
  )::BIGINT AS n_private_pc
FROM genes;

CREATE OR REPLACE TEMP TABLE private_where AS
SELECT
  g.grp,
  count(*)::BIGINT AS n_private_pc,
  count(*) FILTER (WHERE g.max_tissue = 'Testis')::BIGINT AS n_testis,
  count(*) FILTER (WHERE g.max_tissue = 'Liver')::BIGINT AS n_liver,
  count(*) FILTER (
    WHERE g.max_tissue IN (
      'Cells - Cultured fibroblasts',
      'Cells - EBV-transformed lymphocytes'
    )
  )::BIGINT AS n_cell_culture,
  count(*) FILTER (
    WHERE g.max_tissue IN (SELECT tissue FROM tissue_map m WHERE m.grp = g.grp)
  )::BIGINT AS n_index
FROM genes g
WHERE g.biotype = 'protein_coding'
  AND g.tau_class = 'private'
GROUP BY g.grp;

CREATE OR REPLACE TEMP TABLE top_tissue_all AS
SELECT grp, tissue, n_all, n_all_ge1, frac_all
FROM (
  SELECT
    grp,
    tissue,
    n_all,
    n_all_ge1,
    n_all_ge1::DOUBLE / nullif(n_all, 0) AS frac_all,
    row_number() OVER (
      PARTITION BY grp
      ORDER BY n_all_ge1::DOUBLE / nullif(n_all, 0) DESC, tissue
    ) AS rn
  FROM tissue_stat
)
WHERE rn = 1;
