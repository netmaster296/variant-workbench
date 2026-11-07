-- Counts and top genes for somatic origin and for oncogenicity.
-- Read-only. Does not write a mart and does not rank genes by LOEUF,
-- expression, or germline review stars.
--
-- Denominator for both questions: every row of allele_dedup.
-- Somatic origin is germline_somatic = 'somatic' only.
-- Mixed germline/somatic stays its own class.
-- Oncogenicity is not copied into germline_pathogenic.

CREATE OR REPLACE TEMP VIEW allele AS
SELECT
  dedup_key,
  variant_key,
  GeneSymbol,
  ClinicalSignificance,
  ReviewStatus,
  OriginSimple,
  germline_somatic,
  Oncogenicity,
  SomaticClinicalImpact,
  review_stars,
  germline_pathogenic,
  contains(lower(coalesce(ClinicalSignificance, '')), 'drug response') AS drug_response,
  contains(lower(coalesce(ClinicalSignificance, '')), 'pathogenic')
    AND NOT contains(lower(coalesce(ClinicalSignificance, '')), 'conflict') AS pathogenic_text,
  (
    contains(lower(coalesce(ClinicalSignificance, '')), 'uncertain significance')
    OR lower(coalesce(ClinicalSignificance, '')) LIKE 'vus%'
  ) AS vus_text,
  NOT (
    Oncogenicity IS NULL
    OR lower(trim(Oncogenicity)) IN ('', '-', '.', 'na', 'n/a', 'none', 'not provided')
    OR lower(Oncogenicity) LIKE 'no classification%'
    OR lower(Oncogenicity) LIKE 'no assertion%'
    OR lower(Oncogenicity) LIKE 'no interpretation%'
  ) AS onco_call,
  Oncogenicity IN (
    'Oncogenic',
    'Likely oncogenic',
    'Oncogenic/Likely oncogenic'
  ) AS onco_positive,
  SomaticClinicalImpact IN (
    'Tier I - Strong',
    'Tier II - Potential',
    'Tier III - Unknown',
    'Tier IV - Benign/Likely benign'
  ) AS impact_tier
FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet');

CREATE OR REPLACE TEMP VIEW token AS
SELECT
  a.*,
  trim(u.token) AS symbol
FROM allele a,
unnest(string_split(coalesce(a.GeneSymbol, ''), chr(59))) AS u(token)
WHERE trim(u.token) NOT IN ('', '-', '.', 'na', 'n/a', 'none');

-- Origin class must match OriginSimple. Mixed stays mixed.
SELECT
  count(*) AS rows,
  count(*) FILTER (
    WHERE germline_somatic IS DISTINCT FROM CASE OriginSimple
      WHEN 'somatic' THEN 'somatic'
      WHEN 'germline/somatic' THEN 'germline/somatic'
      WHEN 'germline' THEN 'germline'
      ELSE 'other'
    END
  ) AS origin_class_mismatches,
  count(*) FILTER (
    WHERE germline_pathogenic
      AND (
        germline_somatic <> 'germline'
        OR onco_call
      )
  ) AS germline_pathogenic_leakage
FROM allele;

SELECT
  germline_somatic,
  count(*) AS alleles,
  round(count(*) / (SELECT count(*) FROM allele), 6) AS fraction_of_alleles,
  count(*) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM allele
GROUP BY germline_somatic
ORDER BY alleles DESC;

SELECT OriginSimple, germline_somatic, count(*) AS alleles
FROM allele
WHERE germline_somatic = 'other'
GROUP BY 1, 2
ORDER BY alleles DESC;

SELECT
  'somatic_origin' AS set_name,
  count(*) FILTER (WHERE germline_somatic = 'somatic') AS alleles
FROM allele
UNION ALL
SELECT 'mixed_origin', count(*) FILTER (WHERE germline_somatic = 'germline/somatic')
FROM allele
UNION ALL
SELECT 'oncogenicity_call', count(*) FILTER (WHERE onco_call)
FROM allele
UNION ALL
SELECT 'oncogenic_direction', count(*) FILTER (WHERE onco_positive)
FROM allele
UNION ALL
SELECT 'impact_tier', count(*) FILTER (WHERE impact_tier)
FROM allele
UNION ALL
SELECT 'germline_pathogenic', count(*) FILTER (WHERE germline_pathogenic)
FROM allele
UNION ALL
SELECT 'somatic_and_onco_call', count(*) FILTER (WHERE germline_somatic = 'somatic' AND onco_call)
FROM allele
UNION ALL
SELECT 'somatic_and_onco_positive', count(*) FILTER (WHERE germline_somatic = 'somatic' AND onco_positive)
FROM allele
UNION ALL
SELECT 'mixed_and_onco_call', count(*) FILTER (WHERE germline_somatic = 'germline/somatic' AND onco_call)
FROM allele
UNION ALL
SELECT 'mixed_and_onco_positive', count(*) FILTER (WHERE germline_somatic = 'germline/somatic' AND onco_positive)
FROM allele
UNION ALL
SELECT 'germline_and_onco_call', count(*) FILTER (WHERE germline_somatic = 'germline' AND onco_call)
FROM allele
UNION ALL
SELECT 'other_and_onco_call', count(*) FILTER (WHERE germline_somatic = 'other' AND onco_call)
FROM allele
UNION ALL
SELECT 'impact_and_onco_call', count(*) FILTER (WHERE impact_tier AND onco_call)
FROM allele
UNION ALL
SELECT 'somatic_drug_response_exact', count(*) FILTER (
  WHERE germline_somatic = 'somatic' AND ClinicalSignificance = 'drug response'
)
FROM allele
UNION ALL
SELECT 'somatic_drug_response_any', count(*) FILTER (
  WHERE germline_somatic = 'somatic' AND drug_response
)
FROM allele
UNION ALL
SELECT 'somatic_vus_text', count(*) FILTER (WHERE germline_somatic = 'somatic' AND vus_text)
FROM allele
UNION ALL
SELECT 'somatic_significance_blank', count(*) FILTER (
  WHERE germline_somatic = 'somatic' AND ClinicalSignificance = '-'
)
FROM allele
UNION ALL
SELECT 'somatic_pathogenic_text', count(*) FILTER (
  WHERE germline_somatic = 'somatic' AND pathogenic_text
)
FROM allele
UNION ALL
SELECT 'mixed_pathogenic_text', count(*) FILTER (
  WHERE germline_somatic = 'germline/somatic' AND pathogenic_text
)
FROM allele
UNION ALL
SELECT 'somatic_no_symbol', count(*) FILTER (
  WHERE germline_somatic = 'somatic'
    AND (
      GeneSymbol IS NULL
      OR trim(GeneSymbol) IN ('', '-', '.', 'na', 'n/a', 'none')
    )
)
FROM allele
UNION ALL
SELECT 'somatic_multi_symbol', count(*) FILTER (
  WHERE germline_somatic = 'somatic' AND contains(GeneSymbol, chr(59))
)
FROM allele
UNION ALL
SELECT 'mixed_no_symbol', count(*) FILTER (
  WHERE germline_somatic = 'germline/somatic'
    AND (
      GeneSymbol IS NULL
      OR trim(GeneSymbol) IN ('', '-', '.', 'na', 'n/a', 'none')
    )
)
FROM allele
UNION ALL
SELECT 'mixed_multi_symbol', count(*) FILTER (
  WHERE germline_somatic = 'germline/somatic' AND contains(GeneSymbol, chr(59))
)
FROM allele
UNION ALL
SELECT 'onco_call_multi_symbol', count(*) FILTER (WHERE onco_call AND contains(GeneSymbol, chr(59)))
FROM allele
UNION ALL
SELECT 'onco_positive_multi_symbol', count(*) FILTER (
  WHERE onco_positive AND contains(GeneSymbol, chr(59))
)
FROM allele
UNION ALL
SELECT 'gene_symbol_contains_pipe', count(*) FILTER (WHERE contains(GeneSymbol, '|'))
FROM allele;

SELECT ClinicalSignificance, count(*) AS alleles
FROM allele
WHERE germline_somatic = 'somatic'
GROUP BY 1
ORDER BY alleles DESC, ClinicalSignificance;

SELECT review_stars, count(*) AS alleles
FROM allele
WHERE germline_somatic = 'somatic'
GROUP BY 1
ORDER BY review_stars NULLS LAST;

SELECT ReviewStatus, count(*) AS alleles
FROM allele
WHERE germline_somatic = 'somatic' AND review_stars IS NULL
GROUP BY 1
ORDER BY alleles DESC, ReviewStatus;

SELECT Oncogenicity, count(*) AS alleles,
  count(*) FILTER (WHERE germline_somatic = 'somatic') AS somatic_n,
  count(*) FILTER (WHERE germline_somatic = 'germline/somatic') AS mixed_n,
  count(*) FILTER (WHERE germline_somatic = 'germline') AS germline_n,
  count(*) FILTER (WHERE germline_somatic = 'other') AS other_n,
  count(*) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM allele
GROUP BY Oncogenicity
ORDER BY alleles DESC, Oncogenicity;

SELECT SomaticClinicalImpact, count(*) AS alleles,
  count(*) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM allele
GROUP BY 1
ORDER BY alleles DESC, SomaticClinicalImpact;

SELECT
  variant_key,
  GeneSymbol,
  ClinicalSignificance,
  OriginSimple,
  germline_somatic,
  Oncogenicity,
  SomaticClinicalImpact,
  review_stars,
  germline_pathogenic
FROM allele
WHERE variant_key = 'chr17:31334927:C:T';

SELECT
  variant_key,
  GeneSymbol,
  OriginSimple,
  germline_somatic,
  Oncogenicity,
  ClinicalSignificance,
  germline_pathogenic
FROM allele
WHERE onco_positive AND germline_somatic = 'other';

SELECT
  count(DISTINCT symbol) FILTER (WHERE germline_somatic = 'somatic') AS somatic_genes,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic') AS somatic_alleles_with_symbol,
  count(*) FILTER (WHERE germline_somatic = 'somatic') AS somatic_token_rows,
  count(DISTINCT symbol) FILTER (WHERE germline_somatic = 'germline/somatic') AS mixed_genes,
  count(*) FILTER (WHERE germline_somatic = 'germline/somatic') AS mixed_token_rows,
  count(DISTINCT symbol) FILTER (WHERE onco_call) AS onco_call_genes,
  count(*) FILTER (WHERE onco_call) AS onco_call_token_rows,
  count(DISTINCT symbol) FILTER (WHERE onco_positive) AS onco_positive_genes,
  count(*) FILTER (WHERE onco_positive) AS onco_positive_token_rows
FROM token;

-- Top genes by somatic-origin allele count. Not sorted by germline_pathogenic.
-- drug_response_n, vus_n, review_stars_0_n, and oncogenic_direction_n overlap.
SELECT
  symbol,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic') AS somatic_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic' AND drug_response) AS drug_response_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic' AND vus_text) AS vus_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic' AND review_stars = 0) AS review_stars_0_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic' AND onco_positive) AS oncogenic_direction_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM token
GROUP BY symbol
HAVING count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic') > 0
ORDER BY somatic_n DESC, symbol
LIMIT 15;

-- Mixed origin, held out of both the somatic-origin count and the germline flag.
-- pathogenic_text_n is clinical-significance text, not germline_pathogenic.
SELECT
  symbol,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'germline/somatic') AS mixed_n,
  count(DISTINCT dedup_key) FILTER (
    WHERE germline_somatic = 'germline/somatic' AND pathogenic_text
  ) AS pathogenic_text_n,
  count(DISTINCT dedup_key) FILTER (
    WHERE germline_somatic = 'germline/somatic' AND onco_positive
  ) AS oncogenic_direction_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM token
GROUP BY symbol
HAVING count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'germline/somatic') > 0
ORDER BY mixed_n DESC, symbol
LIMIT 15;

-- Top genes by oncogenic or likely oncogenic alleles. Uncertain, benign,
-- and conflicting oncogenicity stay in the call count and out of this sort.
SELECT
  symbol,
  count(DISTINCT dedup_key) FILTER (WHERE onco_positive) AS oncogenic_direction_n,
  count(DISTINCT dedup_key) FILTER (WHERE Oncogenicity = 'Oncogenic') AS oncogenic_n,
  count(DISTINCT dedup_key) FILTER (WHERE Oncogenicity = 'Likely oncogenic') AS likely_oncogenic_n,
  count(DISTINCT dedup_key) FILTER (WHERE onco_call) AS oncogenicity_call_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_somatic = 'somatic' AND onco_positive) AS somatic_origin_n,
  count(DISTINCT dedup_key) FILTER (
    WHERE germline_somatic = 'germline/somatic' AND onco_positive
  ) AS mixed_origin_n,
  count(DISTINCT dedup_key) FILTER (WHERE germline_pathogenic) AS germline_pathogenic_n
FROM token
GROUP BY symbol
HAVING count(DISTINCT dedup_key) FILTER (WHERE onco_positive) > 0
ORDER BY oncogenic_direction_n DESC, symbol
LIMIT 15;
