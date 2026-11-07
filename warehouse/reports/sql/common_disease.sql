-- Shared-gene hypotheses. Not colocalization.
-- A mapped gene is a hypothesis when it is an eGene (qval_le_0_05) in a QTL
-- tissue assigned to the trait. Locus start/end are the member span. They are
-- not expanded to 1 Mb and they are not an overlap window.
--
-- QTL tissue strings are the egenes.parquet labels (underscores). Expression
-- tissue strings are the median_tpm.parquet labels (spaces). tissue_map is the
-- only translation. The two strings are not treated as equal.
--
-- Trait names come from traits.parquet. The eGene test uses egenes.parquet.
-- Gene identity uses the gene mart (symbol -> one ENSG).

CREATE OR REPLACE TEMP VIEW loci AS
SELECT * FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet');

CREATE OR REPLACE TEMP VIEW traits AS
SELECT * FROM read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet');

CREATE OR REPLACE TEMP VIEW egenes AS
SELECT * FROM read_parquet('D:/variant-workbench/warehouse/gtex/egenes.parquet');

CREATE OR REPLACE TEMP VIEW mart AS
SELECT * FROM read_parquet('D:/variant-workbench/warehouse/marts/gene_mart.parquet');

CREATE OR REPLACE TEMP VIEW expression_tissues AS
SELECT DISTINCT tissue
FROM read_parquet('D:/variant-workbench/warehouse/gtex/median_tpm.parquet');

CREATE OR REPLACE TEMP TABLE tissue_map (
  class_name VARCHAR,
  qtl_tissue VARCHAR,
  expression_tissue VARCHAR
);

INSERT INTO tissue_map VALUES
  ('blood', 'Whole_Blood', 'Whole Blood'),
  ('liver_enzyme', 'Liver', 'Liver'),
  ('heart', 'Heart_Atrial_Appendage', 'Heart - Atrial Appendage'),
  ('heart', 'Heart_Left_Ventricle', 'Heart - Left Ventricle'),
  ('brain', 'Brain_Amygdala', 'Brain - Amygdala'),
  ('brain', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain - Anterior cingulate cortex (BA24)'),
  ('brain', 'Brain_Caudate_basal_ganglia', 'Brain - Caudate (basal ganglia)'),
  ('brain', 'Brain_Cerebellar_Hemisphere', 'Brain - Cerebellar Hemisphere'),
  ('brain', 'Brain_Cerebellum', 'Brain - Cerebellum'),
  ('brain', 'Brain_Cortex', 'Brain - Cortex'),
  ('brain', 'Brain_Frontal_Cortex_BA9', 'Brain - Frontal Cortex (BA9)'),
  ('brain', 'Brain_Hippocampus', 'Brain - Hippocampus'),
  ('brain', 'Brain_Hypothalamus', 'Brain - Hypothalamus'),
  ('brain', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain - Nucleus accumbens (basal ganglia)'),
  ('brain', 'Brain_Putamen_basal_ganglia', 'Brain - Putamen (basal ganglia)'),
  ('brain', 'Brain_Spinal_cord_cervical_c-1', 'Brain - Spinal cord (cervical c-1)'),
  ('brain', 'Brain_Substantia_nigra', 'Brain - Substantia nigra');

-- Name rule for traits with at least 50 clumped loci. A trait that matches
-- two classes is left unmapped. The minimum locus count is 50.
CREATE OR REPLACE TEMP TABLE trait_loci AS
SELECT
  l.trait_uri,
  count(*) AS n_loci,
  count(DISTINCT struct_pack(chrom := l.chrom, lead_pos := l.lead_pos)) AS n_locus_keys,
  any_value(t.mapped_trait) AS mapped_trait,
  count(DISTINCT t.mapped_trait) AS n_names,
  count(*) FILTER (WHERE t.trait_uri IS NULL AND l.trait_uri IS NOT NULL) AS loci_missing_trait_row
FROM loci l
LEFT JOIN traits t
  ON l.trait_uri IS NOT DISTINCT FROM t.trait_uri
GROUP BY l.trait_uri;

CREATE OR REPLACE TEMP TABLE classified AS
SELECT
  trait_uri,
  n_loci,
  n_locus_keys,
  mapped_trait,
  n_names,
  loci_missing_trait_row,
  (is_blood::INTEGER + is_brain::INTEGER + is_heart::INTEGER + is_liver::INTEGER) AS n_classes,
  CASE
    WHEN (is_blood::INTEGER + is_brain::INTEGER + is_heart::INTEGER + is_liver::INTEGER) <> 1 THEN NULL
    WHEN is_blood THEN 'blood'
    WHEN is_brain THEN 'brain'
    WHEN is_heart THEN 'heart'
    WHEN is_liver THEN 'liver_enzyme'
  END AS class_name
FROM (
  SELECT
    *,
    (
      regexp_matches(lower(coalesce(mapped_trait, '')),
        '(hematolog|haematolog|hematocrit|haematocrit|hemoglobin|haemoglobin|erythrocyte|red blood cell|red cell|leukocyte|leucocyte|white blood cell|white cell|neutrophil|lymphocyte|monocyte|eosinophil|basophil|platelet|thrombocyte|reticulocyte|corpuscular|granulocyte|myeloid leukocyte|blood cell|hematopoiesis|haematopoiesis)')
      AND NOT regexp_matches(lower(mapped_trait),
        '(leukemia|leukaemia|lymphoma|neoplasm|cancer|carcinoma|tumor|tumour|esophagitis|platelet-derived|monounsaturated|monoclonal|hemoglobin a1|haemoglobin a1|glycated|glycosylated|a1c|blood pressure)')
    ) AS is_blood,
    (
      regexp_matches(lower(coalesce(mapped_trait, '')),
        '(alzheimer|parkinson|lewy body|dementia|schizophrenia|bipolar|depressive|depression|attention deficit|(^|[^a-z])adhd([^a-z]|$)|autism|epilepsy|migraine|neuroticism|neurotic disorder|psychiatric disorder|psychotic|insomnia|cognitive|intelligence|(^|[^a-z])memory([^a-z]|$)|mathematical ability|neuroimaging|neuropsycholog|neurofibrillary|(^|[^a-z])brain([^a-z]|$)|amygdala volume|cerebellar volume|cerebral cortex|cortical thickness|hippocampal volume|ischemic stroke|stroke disorder|large artery stroke|small vessel stroke|stroke outcome|sleep duration|sleep characteristic|white matter|grey matter|gray matter|multiple sclerosis|amyotrophic lateral sclerosis|obsessive-compulsive|tourette|(^|[^a-z])anxiety([^a-z]|$)|narcolepsy|cerebral amyloid)')
      AND NOT regexp_matches(lower(mapped_trait),
        '(natriuretic|brain-derived|brain derived|neuroblastoma|neurofilament|neuropathy|neuropathic|sleep apnea|antidepressant|antipsychotic|educational attainment)')
    ) AS is_brain,
    (
      regexp_matches(lower(coalesce(mapped_trait, '')),
        '(qrs|qt interval|pr interval|jt interval|atrial fibrillation|atrial flutter|angina pectoris|cardiac troponin|coronary artery|coronary atherosclerosis|heart failure|heart rate|heart malformation|cardiomyopathy|left ventricular|myocardial infarction|myocardial ischemia|imaging of the heart|sudden cardiac arrest|ventricular septal defect|cardiac arrest|aortic stenosis|aortic valve)')
      AND NOT regexp_matches(lower(mapped_trait),
        '(blood pressure|pulse pressure|stroke)')
    ) AS is_heart,
    (
      regexp_matches(lower(coalesce(mapped_trait, '')),
        '(alanine aminotransferase|alanine transaminase|aspartate aminotransferase|aspartate transaminase|gamma-glutamyl transferase|gamma-glutamyltransferase|alkaline phosphatase)')
      AND NOT regexp_matches(lower(mapped_trait),
        '(bone alkaline|intestinal alkaline|phenylalanine)')
    ) AS is_liver
  FROM trait_loci
  WHERE n_loci >= 50
) flagged;

CREATE OR REPLACE TEMP TABLE symbol_unique AS
SELECT symbol, any_value(ensg) AS ensg, any_value(biotype) AS biotype,
       any_value(symbol_n_gencode) AS symbol_n_gencode
FROM mart
GROUP BY symbol
HAVING count(*) = 1;

CREATE OR REPLACE TEMP TABLE symbol_ambiguous AS
SELECT symbol, count(*) AS n_genes
FROM mart
GROUP BY symbol
HAVING count(*) > 1;

-- Tested and significant rows in the assigned QTL tissues only.
CREATE OR REPLACE TEMP TABLE class_tested AS
SELECT
  m.class_name,
  e.ensg,
  count(DISTINCT e.tissue) AS n_tested_tissues,
  count(DISTINCT e.tissue) FILTER (WHERE e.qval_le_0_05) AS n_sig_tissues
FROM egenes e
JOIN tissue_map m ON e.tissue = m.qtl_tissue
GROUP BY m.class_name, e.ensg;

CREATE OR REPLACE TEMP TABLE any_egene AS
SELECT ensg, count(DISTINCT tissue) AS n_sig_tissues_any
FROM egenes
WHERE qval_le_0_05
GROUP BY ensg;

CREATE OR REPLACE TEMP TABLE tokens AS
SELECT
  c.class_name,
  c.trait_uri,
  c.mapped_trait,
  c.n_loci,
  l.chrom,
  l.lead_pos,
  l.lead_rsid,
  l.start,
  l."end" AS end_pos,
  unnest(l.mapped_genes.tokens) AS token
FROM classified c
JOIN loci l ON l.trait_uri IS NOT DISTINCT FROM c.trait_uri
WHERE c.class_name IS NOT NULL;

CREATE OR REPLACE TEMP TABLE token_status AS
SELECT
  t.*,
  u.ensg,
  u.biotype,
  u.symbol_n_gencode,
  (t.token LIKE '% - %') AS distance_label,
  (u.ensg IS NOT NULL) AS unique_symbol,
  (a.symbol IS NOT NULL) AS ambiguous_symbol
FROM tokens t
LEFT JOIN symbol_unique u ON t.token = u.symbol
LEFT JOIN symbol_ambiguous a ON t.token = a.symbol;

-- One row per assigned locus and uniquely resolved gene.
-- A gene with no class_tested row was not tested in that class. That is not an eGene call.
CREATE OR REPLACE TEMP TABLE links AS
SELECT
  s.class_name,
  s.trait_uri,
  s.mapped_trait,
  s.chrom,
  s.lead_pos,
  s.ensg,
  any_value(s.biotype) AS biotype,
  any_value(s.token) AS symbol,
  any_value(s.symbol_n_gencode) AS symbol_n_gencode,
  any_value(s.start) AS start_pos,
  any_value(s.end_pos) AS end_pos,
  coalesce(max(ct.n_tested_tissues), 0) AS n_tested_tissues,
  coalesce(max(ct.n_sig_tissues), 0) AS n_sig_tissues,
  coalesce(max(ct.n_sig_tissues), 0) > 0 AS is_hypothesis,
  count(*) FILTER (WHERE ae.ensg IS NOT NULL) > 0 AS egene_any
FROM token_status s
LEFT JOIN class_tested ct
  ON s.ensg = ct.ensg AND s.class_name = ct.class_name
LEFT JOIN any_egene ae ON s.ensg = ae.ensg
WHERE s.unique_symbol
  AND s.symbol_n_gencode = 1
GROUP BY s.class_name, s.trait_uri, s.mapped_trait, s.chrom, s.lead_pos, s.ensg;

CREATE OR REPLACE TEMP TABLE trait_gene AS
SELECT
  class_name,
  trait_uri,
  any_value(mapped_trait) AS mapped_trait,
  ensg,
  any_value(symbol) AS symbol,
  any_value(biotype) AS biotype,
  count(*) AS n_loci_mapped,
  any_value(n_tested_tissues) AS n_tested_tissues,
  any_value(n_sig_tissues) AS n_sig_tissues,
  bool_or(is_hypothesis) AS is_hypothesis,
  bool_or(egene_any) AS egene_any
FROM links
GROUP BY class_name, trait_uri, ensg;

CREATE OR REPLACE TEMP TABLE result_denominator AS
SELECT
  (SELECT count(*) FROM trait_loci) AS traits_with_loci,
  (SELECT count(*) FROM loci) AS loci_all,
  50 AS min_loci,
  (SELECT count(*) FROM classified) AS traits_ge_min,
  (SELECT coalesce(sum(n_loci), 0) FROM classified) AS loci_ge_min,
  (SELECT count(*) FROM classified WHERE class_name IS NOT NULL) AS traits_assigned,
  (SELECT coalesce(sum(n_loci), 0) FROM classified WHERE class_name IS NOT NULL) AS loci_assigned,
  (SELECT count(*) FROM classified WHERE class_name IS NULL) AS traits_unmapped,
  (SELECT coalesce(sum(n_loci), 0) FROM classified WHERE class_name IS NULL) AS loci_unmapped,
  (SELECT count(*) FROM classified WHERE n_classes > 1) AS traits_multi_class,
  (SELECT count(*) FROM classified WHERE mapped_trait IS NULL) AS traits_null_name,
  (SELECT coalesce(sum(loci_missing_trait_row), 0) FROM classified) AS loci_missing_trait_row,
  (SELECT count(*) FROM traits) AS trait_rows,
  (SELECT count(DISTINCT trait_uri) FROM traits) AS trait_uris,
  (SELECT quantile_cont(n_loci, 0.5) FROM trait_loci) AS loci_p50,
  (SELECT quantile_cont(n_loci, 0.9) FROM trait_loci) AS loci_p90,
  (SELECT count(*) FROM trait_loci WHERE n_loci < 50) AS traits_below_min;

CREATE OR REPLACE TEMP TABLE result_class AS
SELECT
  c.class_name,
  count(*) AS n_traits,
  sum(c.n_loci) AS n_loci,
  (SELECT count(*) FROM tokens t WHERE t.class_name = c.class_name) AS token_rows,
  (SELECT count(*) FROM token_status t WHERE t.class_name = c.class_name AND NOT t.unique_symbol AND NOT t.ambiguous_symbol) AS tokens_unmatched,
  (SELECT count(*) FROM token_status t WHERE t.class_name = c.class_name AND t.ambiguous_symbol) AS tokens_ambiguous,
  (SELECT count(*) FROM token_status t WHERE t.class_name = c.class_name AND t.unique_symbol AND t.symbol_n_gencode = 1) AS tokens_unique,
  (SELECT count(*) FROM links k WHERE k.class_name = c.class_name) AS locus_gene_links,
  (SELECT count(*) FROM links k WHERE k.class_name = c.class_name AND k.is_hypothesis) AS hypothesis_links,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name) AS trait_genes,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name AND g.is_hypothesis) AS hypothesis_trait_genes,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name AND g.n_tested_tissues > 0) AS tested_trait_genes,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name AND g.n_tested_tissues > 0 AND NOT g.is_hypothesis) AS tested_not_egene_trait_genes,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name AND g.n_tested_tissues = 0) AS not_tested_trait_genes,
  (SELECT count(*) FROM trait_gene g WHERE g.class_name = c.class_name AND g.egene_any) AS any_tissue_egene_trait_genes,
  (SELECT count(DISTINCT g.ensg) FROM trait_gene g WHERE g.class_name = c.class_name) AS distinct_resolved_genes,
  (SELECT count(DISTINCT g.ensg) FROM trait_gene g WHERE g.class_name = c.class_name AND g.is_hypothesis) AS distinct_hypothesis_genes
FROM classified c
WHERE c.class_name IS NOT NULL
GROUP BY c.class_name
ORDER BY c.class_name;

CREATE OR REPLACE TEMP TABLE trait_gene_counts AS
SELECT
  class_name,
  trait_uri,
  count(*) AS resolved_genes,
  count(*) FILTER (WHERE is_hypothesis) AS hypothesis_genes,
  count(*) FILTER (WHERE n_tested_tissues > 0 AND NOT is_hypothesis) AS tested_not_egene,
  count(*) FILTER (WHERE n_tested_tissues = 0) AS not_tested,
  count(*) FILTER (WHERE egene_any) AS any_tissue_egene
FROM trait_gene
GROUP BY class_name, trait_uri;

CREATE OR REPLACE TEMP TABLE trait_locus_counts AS
SELECT
  class_name,
  trait_uri,
  count(DISTINCT struct_pack(chrom := chrom, lead_pos := lead_pos)) FILTER (WHERE is_hypothesis) AS loci_with_hypothesis,
  count(DISTINCT struct_pack(chrom := chrom, lead_pos := lead_pos)) AS loci_with_resolved_gene
FROM links
GROUP BY class_name, trait_uri;

CREATE OR REPLACE TEMP TABLE result_trait AS
SELECT
  c.class_name,
  c.n_loci,
  c.mapped_trait,
  c.trait_uri,
  coalesce(g.resolved_genes, 0) AS resolved_genes,
  coalesce(g.hypothesis_genes, 0) AS hypothesis_genes,
  coalesce(g.tested_not_egene, 0) AS tested_not_egene,
  coalesce(g.not_tested, 0) AS not_tested,
  coalesce(g.any_tissue_egene, 0) AS any_tissue_egene,
  coalesce(k.loci_with_hypothesis, 0) AS loci_with_hypothesis,
  coalesce(k.loci_with_resolved_gene, 0) AS loci_with_resolved_gene
FROM classified c
LEFT JOIN trait_gene_counts g
  ON c.trait_uri IS NOT DISTINCT FROM g.trait_uri
 AND c.class_name = g.class_name
LEFT JOIN trait_locus_counts k
  ON c.trait_uri IS NOT DISTINCT FROM k.trait_uri
 AND c.class_name = k.class_name
WHERE c.class_name IS NOT NULL
ORDER BY c.class_name, c.n_loci DESC, c.mapped_trait;

CREATE OR REPLACE TEMP TABLE result_tissue AS
SELECT
  m.class_name,
  m.qtl_tissue,
  m.expression_tissue,
  (SELECT count(*) FROM egenes e WHERE e.tissue = m.qtl_tissue) AS tested_in_file,
  (SELECT count(*) FROM egenes e WHERE e.tissue = m.qtl_tissue AND e.qval_le_0_05) AS egenes_in_file,
  (
    SELECT count(DISTINCT g.ensg)
    FROM trait_gene g
    JOIN egenes e ON e.ensg = g.ensg AND e.tissue = m.qtl_tissue
    WHERE g.class_name = m.class_name
  ) AS resolved_genes_tested_here,
  (
    SELECT count(DISTINCT g.ensg)
    FROM trait_gene g
    JOIN egenes e ON e.ensg = g.ensg AND e.tissue = m.qtl_tissue AND e.qval_le_0_05
    WHERE g.class_name = m.class_name
  ) AS resolved_genes_egene_here,
  (
    SELECT count(*)
    FROM trait_gene g
    JOIN egenes e ON e.ensg = g.ensg AND e.tissue = m.qtl_tissue AND e.qval_le_0_05
    WHERE g.class_name = m.class_name
  ) AS trait_genes_egene_here
FROM tissue_map m
ORDER BY m.class_name, m.qtl_tissue;

CREATE OR REPLACE TEMP TABLE result_gene AS
SELECT
  class_name,
  symbol,
  ensg,
  any_value(biotype) AS biotype,
  count(*) AS n_traits,
  sum(n_loci_mapped) AS n_loci_sum,
  max(n_loci_mapped) AS max_loci_one_trait,
  any_value(n_sig_tissues) AS n_sig_tissues,
  any_value(n_tested_tissues) AS n_tested_tissues
FROM trait_gene
WHERE is_hypothesis
GROUP BY class_name, symbol, ensg
ORDER BY class_name, n_traits DESC, n_loci_sum DESC, symbol;

CREATE OR REPLACE TEMP TABLE result_gene_top AS
SELECT * FROM (
  SELECT
    *,
    row_number() OVER (
      PARTITION BY class_name
      ORDER BY n_traits DESC, n_loci_sum DESC, symbol
    ) AS rk
  FROM result_gene
) ranked
WHERE rk <= 20
ORDER BY class_name, rk;

CREATE OR REPLACE TEMP TABLE result_biotype AS
SELECT class_name, biotype,
       count(*) AS hypothesis_trait_genes,
       count(DISTINCT ensg) AS distinct_genes
FROM trait_gene
WHERE is_hypothesis
GROUP BY class_name, biotype
ORDER BY class_name, hypothesis_trait_genes DESC;

CREATE OR REPLACE TEMP TABLE result_unmapped AS
SELECT n_loci, mapped_trait, trait_uri
FROM classified
WHERE class_name IS NULL
ORDER BY n_loci DESC, mapped_trait;

CREATE OR REPLACE TEMP TABLE result_diag AS
SELECT
  (SELECT count(*) FROM mart) AS mart_genes,
  (SELECT count(*) - count(DISTINCT ensg) FROM mart) AS mart_dup_ensg,
  (SELECT count(*) FROM symbol_ambiguous) AS ambiguous_symbols,
  (SELECT count(*) FROM tokens) AS token_rows,
  (SELECT count(*) FROM token_status WHERE distance_label) AS distance_token_rows,
  (SELECT count(*) FROM token_status WHERE distance_label AND (unique_symbol OR ambiguous_symbol)) AS distance_tokens_matched,
  (SELECT count(*) FROM token_status WHERE NOT distance_label AND NOT unique_symbol AND NOT ambiguous_symbol) AS other_unmatched_rows,
  (SELECT count(*) FROM symbol_unique WHERE symbol_n_gencode <> 1) AS unique_symbol_n_not_1,
  (SELECT count(*) FROM links WHERE symbol_n_gencode <> 1) AS links_symbol_n_not_1,
  (SELECT count(*) FROM classified WHERE n_loci <> n_locus_keys) AS traits_locus_key_collision,
  (SELECT count(*) FROM egenes WHERE qval_le_0_05 IS DISTINCT FROM (qval <= 0.05)) AS qval_flag_disagreements,
  (SELECT min(l."end" - l.start) FROM loci l JOIN classified c ON l.trait_uri IS NOT DISTINCT FROM c.trait_uri WHERE c.class_name IS NOT NULL) AS min_member_span_bp,
  (SELECT median(l."end" - l.start) FROM loci l JOIN classified c ON l.trait_uri IS NOT DISTINCT FROM c.trait_uri WHERE c.class_name IS NOT NULL) AS median_member_span_bp,
  (SELECT max(l."end" - l.start) FROM loci l JOIN classified c ON l.trait_uri IS NOT DISTINCT FROM c.trait_uri WHERE c.class_name IS NOT NULL) AS max_member_span_bp,
  (SELECT count(*) FROM loci l JOIN classified c ON l.trait_uri IS NOT DISTINCT FROM c.trait_uri WHERE c.class_name IS NOT NULL AND l.start = l."end") AS point_loci,
  (SELECT count(*) FROM tissue_map m WHERE NOT EXISTS (
      SELECT 1 FROM egenes e WHERE e.tissue = m.qtl_tissue
  )) AS qtl_tissues_missing,
  (SELECT count(*) FROM tissue_map m WHERE m.expression_tissue NOT IN (
      SELECT tissue FROM expression_tissues
  )) AS expression_labels_missing,
  (SELECT count(*) FROM tissue_map m WHERE m.qtl_tissue = m.expression_tissue) AS tissues_same_string,
  (SELECT count(*) FROM tissue_map) AS tissues_in_map,
  (
    SELECT count(*) FROM (
      SELECT l.ensg
      FROM links l
      JOIN egenes e ON e.ensg = l.ensg AND e.qval_le_0_05
      GROUP BY l.ensg
      HAVING count(DISTINCT e.gene_name) > 1
    )
  ) AS hypothesis_ensg_many_gene_names,
  (
    SELECT count(*) FROM (
      SELECT l.ensg, any_value(l.symbol) AS symbol, any_value(e.gene_name) AS gene_name
      FROM links l
      JOIN egenes e ON e.ensg = l.ensg AND e.qval_le_0_05 AND l.is_hypothesis
      GROUP BY l.ensg
      HAVING any_value(l.symbol) IS DISTINCT FROM any_value(e.gene_name)
         AND count(DISTINCT e.gene_name) = 1
    )
  ) AS hypothesis_symbol_ne_gene_name;

CREATE OR REPLACE TEMP TABLE result_unmatched_token AS
SELECT token, count(*) AS n
FROM token_status
WHERE NOT unique_symbol AND NOT ambiguous_symbol AND NOT distance_label
GROUP BY token
ORDER BY n DESC, token;

CREATE OR REPLACE TEMP TABLE result_ambiguous_token AS
SELECT token, count(*) AS n
FROM token_status
WHERE ambiguous_symbol
GROUP BY token
ORDER BY n DESC, token
LIMIT 15;

CREATE OR REPLACE TEMP TABLE result_spot AS
SELECT
  v.symbol,
  g.class_name,
  g.mapped_trait,
  g.n_loci_mapped,
  g.n_tested_tissues,
  g.n_sig_tissues,
  g.is_hypothesis,
  g.egene_any,
  g.biotype
FROM (VALUES
  ('APOE'), ('BIN1'), ('CR1'), ('GPT'), ('GGT1'), ('PNPLA3'),
  ('HBB'), ('TET2'), ('NOS1AP'), ('SORT1'), ('LDLR'), ('ALPL')
) v(symbol)
LEFT JOIN trait_gene g ON g.symbol = v.symbol
ORDER BY v.symbol, g.class_name, g.n_loci_mapped DESC;

CREATE OR REPLACE TEMP TABLE result_subset_100 AS
SELECT
  class_name,
  count(*) AS n_traits,
  sum(n_loci) AS n_loci,
  sum(resolved_genes) AS resolved_genes,
  sum(hypothesis_genes) AS hypothesis_genes
FROM result_trait
WHERE n_loci >= 100
GROUP BY class_name
ORDER BY class_name;

SELECT 'denominator' AS section, * FROM result_denominator;
SELECT 'class' AS section, * FROM result_class;
SELECT 'trait' AS section, * FROM result_trait;
SELECT 'tissue' AS section, * FROM result_tissue;
SELECT 'gene_top' AS section, * FROM result_gene_top;
SELECT 'biotype' AS section, * FROM result_biotype;
SELECT 'unmapped' AS section, * FROM result_unmapped;
SELECT 'diag' AS section, * FROM result_diag;
SELECT 'unmatched_token' AS section, * FROM result_unmatched_token;
SELECT 'ambiguous_token' AS section, * FROM result_ambiguous_token;
SELECT 'spot' AS section, * FROM result_spot;
SELECT 'subset_100' AS section, * FROM result_subset_100;
SELECT 'gene_n' AS section, class_name, count(*) AS genes, max(n_traits) AS max_traits
FROM result_gene GROUP BY class_name;

CREATE OR REPLACE TEMP TABLE result_sig_hist AS
SELECT class_name, n_sig_tissues, count(*) AS genes
FROM result_gene
GROUP BY class_name, n_sig_tissues
ORDER BY class_name, n_sig_tissues;

CREATE OR REPLACE TEMP TABLE result_name_mismatch AS
SELECT
  l.ensg,
  any_value(l.symbol) AS symbol,
  any_value(e.gene_name) AS gene_name,
  count(DISTINCT l.class_name) AS n_classes
FROM links l
JOIN egenes e ON e.ensg = l.ensg AND e.qval_le_0_05 AND l.is_hypothesis
GROUP BY l.ensg
HAVING count(DISTINCT e.gene_name) = 1
   AND any_value(l.symbol) IS DISTINCT FROM any_value(e.gene_name)
ORDER BY symbol
LIMIT 12;

CREATE OR REPLACE TEMP TABLE result_gene_union AS
SELECT
  (SELECT count(DISTINCT ensg) FROM trait_gene) AS resolved_genes,
  (SELECT count(DISTINCT ensg) FROM trait_gene WHERE is_hypothesis) AS hypothesis_genes,
  (SELECT count(*) FROM trait_gene) AS trait_gene_pairs,
  (SELECT count(*) FROM trait_gene WHERE is_hypothesis) AS hypothesis_pairs,
  (SELECT count(*) FROM links) AS locus_gene_links,
  (SELECT count(*) FROM links WHERE is_hypothesis) AS hypothesis_links,
  (SELECT count(*) FROM trait_gene WHERE n_loci_mapped > 1) AS pairs_on_multiple_loci,
  (
    SELECT count(*) FROM (
      SELECT ensg FROM trait_gene WHERE is_hypothesis GROUP BY ensg HAVING count(DISTINCT class_name) > 1
    )
  ) AS hypothesis_genes_in_two_classes;

SELECT 'gene' AS section, * FROM result_gene;
SELECT 'sig_hist' AS section, * FROM result_sig_hist;
SELECT 'name_mismatch' AS section, * FROM result_name_mismatch;
SELECT 'gene_union' AS section, * FROM result_gene_union;
