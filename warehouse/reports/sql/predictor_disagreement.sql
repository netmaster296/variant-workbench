-- Predictor disagreement on Tier A coding variants.
-- Reads warehouse parquet only. Does not write a mart.
-- Genome columns and exome columns stay in separate rows (side = genome or exome).
-- Allele counts are not added. Predictor scores are not averaged across callsets.
-- MaveDB scores are never pooled across score sets. Haplotype rows are not used.

SET threads = 4;

CREATE TEMP TABLE coding_term (term VARCHAR);
INSERT INTO coding_term VALUES
    ('missense_variant'),
    ('synonymous_variant'),
    ('stop_gained'),
    ('stop_lost'),
    ('start_lost'),
    ('start_retained_variant'),
    ('stop_retained_variant'),
    ('frameshift_variant'),
    ('inframe_insertion'),
    ('inframe_deletion'),
    ('protein_altering_variant'),
    ('incomplete_terminal_codon_variant'),
    ('coding_sequence_variant');

CREATE TEMP TABLE tier AS
SELECT *
FROM read_parquet('D:/variant-workbench/warehouse/variants/tier_a.parquet');

CREATE TEMP TABLE coding AS
SELECT
    'genome' AS side,
    variant_key,
    chrom,
    pos,
    ref,
    alt,
    vep_symbol AS symbol,
    vep_consequence AS consequence,
    vep_impact AS impact,
    vep_on_mane AS on_mane,
    cadd_phred AS cadd,
    revel_max AS revel,
    spliceai_ds_max AS spliceai,
    pangolin_largest_ds AS pangolin,
    sift_max AS sift,
    polyphen_max AS polyphen,
    list_contains(string_split(vep_consequence, '&'), 'missense_variant') AS missense
FROM tier
WHERE in_genomes
  AND len(list_intersect(string_split(vep_consequence, '&'), (SELECT list(term) FROM coding_term))) > 0
UNION ALL
SELECT
    'exome',
    variant_key,
    chrom,
    pos,
    ref,
    alt,
    vep_symbol_exome,
    vep_consequence_exome,
    vep_impact_exome,
    vep_on_mane_exome,
    cadd_phred_exome,
    revel_max_exome,
    spliceai_ds_max_exome,
    pangolin_largest_ds_exome,
    sift_max_exome,
    polyphen_max_exome,
    list_contains(string_split(vep_consequence_exome, '&'), 'missense_variant')
FROM tier
WHERE in_exomes
  AND len(list_intersect(string_split(vep_consequence_exome, '&'), (SELECT list(term) FROM coding_term))) > 0;

-- ClinVar is one row per variant_key. Flags follow the piece split used in load_clinvar.py:
-- phrases are separated by ' | ', and a phrase is split on ';'.
-- A piece that contains 'conflict' is not also counted as pathogenic.
-- Pathogenic and benign pieces on the same key are mixed. That is stricter than
-- germline_pathogenic, which ignores a benign piece. The QC query counts the difference.
CREATE TEMP TABLE cv_flags AS
WITH subset AS (
    SELECT d.*
    FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') d
    SEMI JOIN (SELECT DISTINCT variant_key FROM coding) k
        ON d.variant_key = k.variant_key
),
pieces AS (
    SELECT
        s.variant_key,
        lower(trim(piece)) AS txt
    FROM subset s
    CROSS JOIN unnest(string_split(coalesce(s."ClinicalSignificance", ''), ' | ')) AS u(phrase)
    CROSS JOIN unnest(string_split(phrase, ';')) AS v(piece)
    WHERE length(trim(piece)) > 0
),
sig AS (
    SELECT
        variant_key,
        bool_or(contains(txt, 'conflict')) AS has_conflict,
        bool_or(contains(txt, 'pathogenic') AND NOT contains(txt, 'conflict')) AS has_path,
        bool_or(
            contains(txt, 'benign')
            AND NOT contains(txt, 'conflict')
            AND NOT contains(txt, 'pathogenic')
        ) AS has_benign
    FROM pieces
    GROUP BY variant_key
),
onco_pieces AS (
    SELECT
        s.variant_key,
        lower(trim(piece)) AS txt
    FROM subset s
    CROSS JOIN unnest(string_split(coalesce(s."Oncogenicity", ''), ' | ')) AS u(phrase)
    CROSS JOIN unnest(string_split(phrase, ';')) AS v(piece)
    WHERE length(trim(piece)) > 0
      AND lower(trim(piece)) NOT IN ('', '.', '-', 'na', 'n/a', 'none', 'not provided')
      AND NOT starts_with(lower(trim(piece)), 'no classification')
      AND NOT starts_with(lower(trim(piece)), 'no assertion')
      AND NOT starts_with(lower(trim(piece)), 'no interpretation')
),
onco AS (
    SELECT variant_key, true AS has_onco
    FROM onco_pieces
    GROUP BY variant_key
)
SELECT
    s.variant_key,
    s.n_rows,
    s."ClinicalSignificance" AS clinical_significance,
    s."ReviewStatus" AS review_status,
    s.review_stars,
    s.germline_somatic,
    s."Oncogenicity" AS oncogenicity,
    s.germline_pathogenic,
    s."GeneSymbol" AS clinvar_symbol,
    coalesce(sig.has_conflict, false) AS has_conflict,
    coalesce(sig.has_path, false) AS has_path,
    coalesce(sig.has_benign, false) AS has_benign,
    coalesce(onco.has_onco, false) AS has_onco
FROM subset s
LEFT JOIN sig ON s.variant_key = sig.variant_key
LEFT JOIN onco ON s.variant_key = onco.variant_key;

CREATE TEMP TABLE coding_cv AS
SELECT
    c.*,
    f.n_rows AS clinvar_n_rows,
    f.clinical_significance,
    f.review_status,
    f.review_stars,
    f.germline_somatic,
    f.oncogenicity,
    f.germline_pathogenic,
    f.clinvar_symbol,
    f.has_conflict,
    f.has_path,
    f.has_benign,
    f.has_onco,
    CASE
        WHEN f.variant_key IS NULL THEN 'no_clinvar'
        WHEN f.has_conflict OR (f.has_path AND f.has_benign) THEN 'conflict_or_mixed'
        WHEN f.has_path AND f.germline_somatic = 'germline' AND NOT f.has_onco THEN 'p_lp'
        WHEN f.has_benign AND f.germline_somatic = 'germline' AND NOT f.has_onco THEN 'b_lb'
        ELSE 'other'
    END AS cv_class,
    (
        f.variant_key IS NOT NULL
        AND f.review_stars >= 2
        AND NOT contains(lower(coalesce(f.review_status, '')), 'conflict')
        AND NOT (f.has_conflict OR (f.has_path AND f.has_benign))
        AND (
            (f.has_path AND f.germline_somatic = 'germline' AND NOT f.has_onco)
            OR (f.has_benign AND f.germline_somatic = 'germline' AND NOT f.has_onco)
        )
    ) AS high_star
FROM coding c
LEFT JOIN cv_flags f ON c.variant_key = f.variant_key;

-- One numeric score per (score set, allele key). Disagreeing duplicate scores are dropped.
-- post_mapped_type Haplotype is not an allele score and is not loaded here.
CREATE TEMP TABLE mave_unique AS
SELECT
    variant_key,
    score_set_urn,
    any_value(score) AS score,
    count(*) AS n_rows
FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
WHERE variant_key IS NOT NULL
  AND post_mapped_type = 'Allele'
  AND score IS NOT NULL
GROUP BY variant_key, score_set_urn
HAVING count(DISTINCT score) = 1;

CREATE TEMP TABLE mave_coding AS
SELECT
    c.side,
    c.variant_key,
    c.symbol,
    c.missense,
    c.cadd,
    c.revel,
    c.spliceai,
    c.cv_class,
    c.high_star,
    c.review_stars,
    m.score_set_urn,
    m.score,
    e.target_name,
    e.score_set_title
FROM coding_cv c
JOIN mave_unique m ON c.variant_key = m.variant_key
LEFT JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e
    ON m.score_set_urn = e.score_set_urn;

-- Key and join checks. A multiplied ClinVar join would make coding_cv larger than coding.
SELECT 'qc_keys' AS result_set,
    (SELECT count(*) FROM tier) AS tier_rows,
    (SELECT count(DISTINCT variant_key) FROM tier) AS tier_keys,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') WHERE variant_key IS NOT NULL) AS clinvar_keyed_rows,
    (SELECT count(DISTINCT variant_key) FROM read_parquet('D:/variant-workbench/warehouse/clinvar/allele_dedup.parquet') WHERE variant_key IS NOT NULL) AS clinvar_keys,
    (SELECT count(*) FROM coding) AS coding_rows,
    (SELECT count(*) FROM coding_cv) AS coding_cv_rows,
    (SELECT count(*) FROM mave_unique m ANTI JOIN tier t ON m.variant_key = t.variant_key) AS mave_keys_missing_from_tier,
    (SELECT count(*) FROM mave_unique m ANTI JOIN read_parquet('D:/variant-workbench/warehouse/mavedb/experiments.parquet') e ON m.score_set_urn = e.score_set_urn) AS mave_sets_missing_experiment;

SELECT 'qc_clinvar_flag' AS result_set,
    side,
    count(*) FILTER (WHERE germline_pathogenic AND cv_class = 'p_lp') AS flag_and_class,
    count(*) FILTER (WHERE germline_pathogenic AND cv_class <> 'p_lp') AS flag_not_class,
    count(*) FILTER (WHERE cv_class = 'p_lp' AND NOT germline_pathogenic) AS class_not_flag
FROM coding_cv
GROUP BY side
ORDER BY side;

SELECT 'denominator' AS result_set,
    side,
    count(*) AS coding_n,
    count(*) FILTER (WHERE missense) AS missense_n,
    count(*) FILTER (WHERE on_mane) AS on_mane_n,
    count(*) FILTER (WHERE length(ref) = 1 AND length(alt) = 1) AS snv_n,
    count(DISTINCT variant_key) AS keys_n
FROM coding_cv
GROUP BY side
ORDER BY side;

SELECT 'coding_both_sides' AS result_set,
    count(*) AS keys_coding_both,
    count(*) FILTER (WHERE g.cadd IS NOT NULL AND e.cadd IS NOT NULL AND g.cadd <> e.cadd) AS cadd_both_differ,
    count(*) FILTER (WHERE g.cadd IS NOT NULL AND e.cadd IS NOT NULL) AS cadd_both_present,
    count(*) FILTER (WHERE g.revel IS NOT NULL AND e.revel IS NOT NULL AND g.revel <> e.revel) AS revel_both_differ,
    count(*) FILTER (WHERE g.revel IS NOT NULL AND e.revel IS NOT NULL) AS revel_both_present,
    count(*) FILTER (WHERE g.spliceai IS NOT NULL AND e.spliceai IS NOT NULL AND g.spliceai <> e.spliceai) AS spliceai_both_differ,
    count(*) FILTER (WHERE g.spliceai IS NOT NULL AND e.spliceai IS NOT NULL) AS spliceai_both_present,
    count(*) FILTER (WHERE g.consequence <> e.consequence) AS consequence_differ
FROM coding_cv g
JOIN coding_cv e
    ON g.variant_key = e.variant_key
   AND g.side = 'genome'
   AND e.side = 'exome';

SELECT 'excluded_splice_without_coding_term' AS result_set,
    'genome' AS side,
    count(*) AS in_callset_n,
    count(*) FILTER (
        WHERE (
            list_contains(string_split(vep_consequence, '&'), 'splice_acceptor_variant')
            OR list_contains(string_split(vep_consequence, '&'), 'splice_donor_variant')
        )
        AND len(list_intersect(string_split(vep_consequence, '&'), (SELECT list(term) FROM coding_term))) = 0
    ) AS splice_acceptor_or_donor_not_coding_n,
    count(*) FILTER (
        WHERE vep_impact IN ('HIGH', 'MODERATE')
          AND len(list_intersect(string_split(vep_consequence, '&'), (SELECT list(term) FROM coding_term))) = 0
    ) AS high_or_moderate_not_coding_n
FROM tier
WHERE in_genomes
UNION ALL
SELECT
    'excluded_splice_without_coding_term',
    'exome',
    count(*),
    count(*) FILTER (
        WHERE (
            list_contains(string_split(vep_consequence_exome, '&'), 'splice_acceptor_variant')
            OR list_contains(string_split(vep_consequence_exome, '&'), 'splice_donor_variant')
        )
        AND len(list_intersect(string_split(vep_consequence_exome, '&'), (SELECT list(term) FROM coding_term))) = 0
    ),
    count(*) FILTER (
        WHERE vep_impact_exome IN ('HIGH', 'MODERATE')
          AND len(list_intersect(string_split(vep_consequence_exome, '&'), (SELECT list(term) FROM coding_term))) = 0
    )
FROM tier
WHERE in_exomes
ORDER BY side;

SELECT 'missingness_coding' AS result_set,
    side,
    count(*) AS coding_n,
    count(cadd) AS cadd_n,
    count(*) - count(cadd) AS cadd_missing,
    count(revel) AS revel_n,
    count(*) - count(revel) AS revel_missing,
    count(spliceai) AS spliceai_n,
    count(*) - count(spliceai) AS spliceai_missing,
    count(pangolin) AS pangolin_n,
    count(*) - count(pangolin) AS pangolin_missing,
    count(sift) AS sift_n,
    count(*) - count(sift) AS sift_missing,
    count(polyphen) AS polyphen_n,
    count(*) - count(polyphen) AS polyphen_missing,
    count(*) FILTER (WHERE cv_class <> 'no_clinvar') AS clinvar_n,
    count(*) FILTER (WHERE cv_class = 'no_clinvar') AS clinvar_missing,
    count(*) FILTER (WHERE cv_class <> 'no_clinvar' AND review_stars IS NULL) AS stars_null_n,
    count(DISTINCT variant_key) FILTER (WHERE variant_key IN (SELECT variant_key FROM mave_unique)) AS mave_allele_scored_n
FROM coding_cv
GROUP BY side
ORDER BY side;

SELECT 'missingness_missense' AS result_set,
    side,
    count(*) AS missense_n,
    count(cadd) AS cadd_n,
    count(*) - count(cadd) AS cadd_missing,
    count(revel) AS revel_n,
    count(*) - count(revel) AS revel_missing,
    count(spliceai) AS spliceai_n,
    count(*) - count(spliceai) AS spliceai_missing,
    count(pangolin) AS pangolin_n,
    count(*) - count(pangolin) AS pangolin_missing,
    count(sift) AS sift_n,
    count(*) - count(sift) AS sift_missing,
    count(polyphen) AS polyphen_n,
    count(*) - count(polyphen) AS polyphen_missing
FROM coding_cv
WHERE missense
GROUP BY side
ORDER BY side;

-- Coding sites are fully CADD-scored. The callset, not the coding subset, is where CADD is missing.
SELECT 'missingness_in_callset' AS result_set,
    'genome' AS side,
    count(*) AS in_callset_n,
    count(*) - count(cadd_phred) AS cadd_missing,
    count(*) - count(revel_max) AS revel_missing,
    count(*) - count(spliceai_ds_max) AS spliceai_missing,
    count(*) - count(pangolin_largest_ds) AS pangolin_missing,
    count(*) - count(sift_max) AS sift_missing,
    count(*) - count(polyphen_max) AS polyphen_missing
FROM tier
WHERE in_genomes
UNION ALL
SELECT
    'missingness_in_callset',
    'exome',
    count(*),
    count(*) - count(cadd_phred_exome),
    count(*) - count(revel_max_exome),
    count(*) - count(spliceai_ds_max_exome),
    count(*) - count(pangolin_largest_ds_exome),
    count(*) - count(sift_max_exome),
    count(*) - count(polyphen_max_exome)
FROM tier
WHERE in_exomes
ORDER BY side;

SELECT 'consequence_differ_key' AS result_set,
    g.variant_key,
    g.symbol AS genome_symbol,
    e.symbol AS exome_symbol,
    g.consequence AS genome_consequence,
    e.consequence AS exome_consequence
FROM coding_cv g
JOIN coding_cv e
    ON g.variant_key = e.variant_key
   AND g.side = 'genome'
   AND e.side = 'exome'
WHERE g.consequence <> e.consequence;

SELECT 'mave_coverage' AS result_set,
    side,
    count(DISTINCT variant_key) AS coding_keys,
    count(DISTINCT variant_key) FILTER (WHERE score IS NOT NULL) AS keys_with_allele_score,
    count(*) AS key_set_rows,
    count(DISTINCT score_set_urn) AS score_sets
FROM mave_coding
GROUP BY side
ORDER BY side;

SELECT 'mave_allele_ambiguity' AS result_set,
    (SELECT count(*) FROM mave_unique) AS unique_scored_pairs,
    (SELECT count(DISTINCT variant_key) FROM mave_unique) AS unique_scored_keys,
    (SELECT count(DISTINCT score_set_urn) FROM mave_unique) AS unique_scored_sets,
    (
        SELECT count(*)
        FROM (
            SELECT variant_key, score_set_urn
            FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
            WHERE variant_key IS NOT NULL
              AND post_mapped_type = 'Allele'
              AND score IS NOT NULL
            GROUP BY 1, 2
            HAVING count(DISTINCT score) > 1
        )
    ) AS allele_pairs_score_conflict,
    (
        SELECT count(DISTINCT variant_key)
        FROM read_parquet('D:/variant-workbench/warehouse/mavedb/mapped_variants_scored.parquet')
        WHERE variant_key IS NOT NULL
          AND post_mapped_type = 'Haplotype'
          AND score IS NOT NULL
          AND variant_key NOT IN (SELECT variant_key FROM mave_unique)
    ) AS keys_haplotype_score_only;

SELECT 'star_bucket' AS result_set,
    side,
    CASE
        WHEN cv_class = 'no_clinvar' THEN 'no_clinvar'
        WHEN review_stars IS NULL THEN 'stars_null'
        WHEN review_stars = 1 AND contains(lower(review_status), 'conflicting')
             AND contains(lower(review_status), 'single submitter') THEN 'star1_both_phrases'
        WHEN review_stars = 1 AND contains(lower(review_status), 'conflicting') THEN 'star1_conflicting'
        WHEN review_stars = 1 AND contains(lower(review_status), 'single submitter') THEN 'star1_single'
        WHEN review_stars = 1 THEN 'star1_other'
        ELSE 'star_' || cast(review_stars AS VARCHAR)
    END AS bucket,
    count(*) AS n
FROM coding_cv
GROUP BY 2, 3
ORDER BY 1, 2;

SELECT 'cv_class' AS result_set,
    side,
    cv_class,
    count(*) AS n,
    count(*) FILTER (WHERE high_star) AS high_star_n,
    count(*) FILTER (WHERE clinvar_n_rows > 1) AS multi_row_n
FROM coding_cv
GROUP BY 2, 3
ORDER BY 1, 2;

-- Missense pathogenicity calls. REVEL supporting cuts are Pejaver et al. 2022
-- as restated by the UCSC REVEL track from that paper's Table 2:
--   >= 0.644 pathogenic supporting, <= 0.290 benign supporting.
-- The open interval is not a call.
-- CADD phred >= 20 is the top 1% rank and >= 10 is the top 10% rank.
-- A middle CADD score is not a call. Opposition requires both tools to call.
SELECT 'revel_cadd_3x3' AS result_set,
    side,
    CASE
        WHEN revel >= 0.644 THEN 'revel_p'
        WHEN revel <= 0.290 THEN 'revel_b'
        ELSE 'revel_mid'
    END AS revel_call,
    CASE
        WHEN cadd >= 20 THEN 'cadd_high'
        WHEN cadd < 10 THEN 'cadd_low'
        ELSE 'cadd_mid'
    END AS cadd_call,
    count(*) AS n
FROM coding_cv
WHERE missense AND revel IS NOT NULL AND cadd IS NOT NULL
GROUP BY 2, 3, 4
ORDER BY 1, 2, 3;

SELECT 'revel_cadd_3x3_high_star' AS result_set,
    side,
    cv_class,
    CASE
        WHEN revel >= 0.644 THEN 'revel_p'
        WHEN revel <= 0.290 THEN 'revel_b'
        ELSE 'revel_mid'
    END AS revel_call,
    CASE
        WHEN cadd >= 20 THEN 'cadd_high'
        WHEN cadd < 10 THEN 'cadd_low'
        ELSE 'cadd_mid'
    END AS cadd_call,
    count(*) AS n
FROM coding_cv
WHERE missense
  AND high_star
  AND cv_class IN ('p_lp', 'b_lb')
  AND revel IS NOT NULL
  AND cadd IS NOT NULL
GROUP BY 2, 3, 4, 5
ORDER BY 1, 2, 3, 4;

-- Stricter REVEL tails (Pejaver moderate >= 0.773, strong >= 0.932), same CADD cuts.
SELECT 'revel_stricter' AS result_set,
    side,
    count(*) FILTER (WHERE revel >= 0.644 AND cadd >= 20) AS agree_p_supporting,
    count(*) FILTER (WHERE revel >= 0.773 AND cadd >= 20) AS agree_p_moderate,
    count(*) FILTER (WHERE revel >= 0.932 AND cadd >= 20) AS agree_p_strong,
    count(*) FILTER (WHERE revel <= 0.290 AND cadd < 10) AS agree_b,
    count(*) FILTER (WHERE revel >= 0.644 AND cadd < 10) AS disagree_revel_p_cadd_low,
    count(*) FILTER (WHERE revel >= 0.773 AND cadd < 10) AS disagree_revel_mod_cadd_low,
    count(*) FILTER (WHERE revel <= 0.290 AND cadd >= 20) AS disagree_revel_b_cadd_high,
    count(*) AS both_scored_missense
FROM coding_cv
WHERE missense AND revel IS NOT NULL AND cadd IS NOT NULL
GROUP BY side
ORDER BY side;

-- SpliceAI operating points are the three gnomAD documents for spliceai_ds_max:
-- 0.2 high recall, 0.5 recommended, 0.8 high precision.
-- A negative SpliceAI does not contradict a high REVEL. The opposing pair is a
-- splice call against a benign protein call, or a splice call against a low CADD rank.
SELECT 'spliceai_vs_revel' AS result_set,
    side,
    CASE
        WHEN spliceai >= 0.8 THEN 'splice_0.8'
        WHEN spliceai >= 0.5 THEN 'splice_0.5'
        WHEN spliceai >= 0.2 THEN 'splice_0.2'
        ELSE 'splice_below_0.2'
    END AS splice_bin,
    CASE
        WHEN revel >= 0.644 THEN 'revel_p'
        WHEN revel <= 0.290 THEN 'revel_b'
        ELSE 'revel_mid'
    END AS revel_call,
    count(*) AS n
FROM coding_cv
WHERE missense AND spliceai IS NOT NULL AND revel IS NOT NULL
GROUP BY 2, 3, 4
ORDER BY 1, 2, 3;

SELECT 'spliceai_vs_cadd' AS result_set,
    side,
    CASE
        WHEN spliceai >= 0.8 THEN 'splice_0.8'
        WHEN spliceai >= 0.5 THEN 'splice_0.5'
        WHEN spliceai >= 0.2 THEN 'splice_0.2'
        ELSE 'splice_below_0.2'
    END AS splice_bin,
    CASE
        WHEN cadd >= 20 THEN 'cadd_high'
        WHEN cadd < 10 THEN 'cadd_low'
        ELSE 'cadd_mid'
    END AS cadd_call,
    count(*) AS n
FROM coding_cv
WHERE spliceai IS NOT NULL AND cadd IS NOT NULL
GROUP BY 2, 3, 4
ORDER BY 1, 2, 3;

SELECT 'spliceai_high_star' AS result_set,
    side,
    cv_class,
    count(*) AS scored_n,
    count(*) FILTER (WHERE spliceai >= 0.2) AS ge_0_2,
    count(*) FILTER (WHERE spliceai >= 0.5) AS ge_0_5,
    count(*) FILTER (WHERE spliceai >= 0.8) AS ge_0_8,
    count(*) FILTER (WHERE spliceai >= 0.5 AND revel <= 0.290) AS splice_and_revel_b,
    count(*) FILTER (WHERE spliceai >= 0.5 AND revel >= 0.644) AS splice_and_revel_p,
    count(*) FILTER (WHERE spliceai < 0.2 AND revel <= 0.290 AND cadd < 10) AS all_lean_benign,
    count(*) FILTER (WHERE spliceai >= 0.5 AND revel <= 0.290 AND cadd < 10) AS splice_only
FROM coding_cv
WHERE high_star
  AND cv_class IN ('p_lp', 'b_lb')
  AND spliceai IS NOT NULL
GROUP BY 2, 3
ORDER BY 1, 2;

-- Pangolin loss scores are negative. |pangolin_largest_ds| uses the SpliceAI 0.5
-- point as a shared numeric cut. It is not a published Pangolin threshold.
SELECT 'pangolin_vs_spliceai' AS result_set,
    side,
    count(*) AS both_n,
    count(*) FILTER (WHERE spliceai >= 0.5 AND abs(pangolin) >= 0.5) AS both_ge_0_5,
    count(*) FILTER (WHERE spliceai >= 0.5 AND abs(pangolin) < 0.2) AS spliceai_only,
    count(*) FILTER (WHERE spliceai < 0.2 AND abs(pangolin) >= 0.5) AS pangolin_only,
    count(*) FILTER (WHERE spliceai < 0.2 AND abs(pangolin) < 0.2) AS both_low,
    count(*) FILTER (WHERE pangolin < 0) AS pangolin_negative_n
FROM coding_cv
WHERE spliceai IS NOT NULL AND pangolin IS NOT NULL
GROUP BY side
ORDER BY side;

-- Within one score set, CADD tails. Medians are not compared across sets.
-- A set separates only when each tail has at least 25 variants and at most 25%
-- of each tail sits on the other tail's side of that other tail's median.
CREATE TEMP TABLE mave_cadd_set AS
SELECT
    m.side,
    m.score_set_urn,
    any_value(m.target_name) AS target_name,
    any_value(m.score_set_title) AS score_set_title,
    count(*) FILTER (WHERE m.cadd >= 20) AS n_h,
    count(*) FILTER (WHERE m.cadd < 10) AS n_l,
    median(m.score) FILTER (WHERE m.cadd >= 20) AS med_h,
    median(m.score) FILTER (WHERE m.cadd < 10) AS med_l,
    count(*) FILTER (
        WHERE m.cadd >= 20
          AND s.med_h > s.med_l
          AND m.score <= s.med_l
    ) AS n_h_below_low_median,
    count(*) FILTER (
        WHERE m.cadd < 10
          AND s.med_h > s.med_l
          AND m.score >= s.med_h
    ) AS n_l_above_high_median,
    count(*) FILTER (
        WHERE m.cadd >= 20
          AND s.med_h < s.med_l
          AND m.score >= s.med_l
    ) AS n_h_above_low_median,
    count(*) FILTER (
        WHERE m.cadd < 10
          AND s.med_h < s.med_l
          AND m.score <= s.med_h
    ) AS n_l_below_high_median
FROM mave_coding m
JOIN (
    SELECT
        side,
        score_set_urn,
        count(*) FILTER (WHERE cadd >= 20) AS n_h,
        count(*) FILTER (WHERE cadd < 10) AS n_l,
        median(score) FILTER (WHERE cadd >= 20) AS med_h,
        median(score) FILTER (WHERE cadd < 10) AS med_l
    FROM mave_coding
    WHERE cadd IS NOT NULL
    GROUP BY 1, 2
) s
    ON m.side = s.side
   AND m.score_set_urn = s.score_set_urn
WHERE m.cadd IS NOT NULL
  AND s.n_h >= 25
  AND s.n_l >= 25
  AND s.med_h IS DISTINCT FROM s.med_l
GROUP BY m.side, m.score_set_urn;

SELECT 'mave_cadd_sets' AS result_set,
    side,
    score_set_urn,
    target_name,
    score_set_title,
    n_h,
    n_l,
    med_h,
    med_l,
    CASE WHEN med_h > med_l THEN n_h_below_low_median ELSE n_h_above_low_median END AS n_h_wrong,
    CASE WHEN med_h > med_l THEN n_l_above_high_median ELSE n_l_below_high_median END AS n_l_wrong,
    (med_h > med_l) AS high_cadd_has_higher_score
FROM mave_cadd_set
ORDER BY side, n_h + n_l DESC;

SELECT 'mave_cadd_summary' AS result_set,
    side,
    count(*) AS sets_tested,
    count(*) FILTER (
        WHERE (CASE WHEN med_h > med_l THEN n_h_below_low_median ELSE n_h_above_low_median END) <= 0.25 * n_h
          AND (CASE WHEN med_h > med_l THEN n_l_above_high_median ELSE n_l_below_high_median END) <= 0.25 * n_l
    ) AS sets_separated,
    sum(n_h) AS sum_n_h,
    sum(n_l) AS sum_n_l
FROM mave_cadd_set
GROUP BY side
ORDER BY side;

-- Same shape for high-star germline classes. Minimum 10 on each side.
-- Cross-set numeric comparison is not done.
CREATE TEMP TABLE mave_cv_set AS
SELECT
    m.side,
    m.score_set_urn,
    any_value(m.target_name) AS target_name,
    any_value(m.score_set_title) AS score_set_title,
    count(*) FILTER (WHERE m.high_star AND m.cv_class = 'p_lp') AS n_p,
    count(*) FILTER (WHERE m.high_star AND m.cv_class = 'b_lb') AS n_b,
    median(m.score) FILTER (WHERE m.high_star AND m.cv_class = 'p_lp') AS med_p,
    median(m.score) FILTER (WHERE m.high_star AND m.cv_class = 'b_lb') AS med_b,
    count(*) FILTER (
        WHERE m.high_star AND m.cv_class = 'p_lp'
          AND s.med_p > s.med_b
          AND m.score <= s.med_b
    ) AS n_p_below_b,
    count(*) FILTER (
        WHERE m.high_star AND m.cv_class = 'b_lb'
          AND s.med_p > s.med_b
          AND m.score >= s.med_p
    ) AS n_b_above_p,
    count(*) FILTER (
        WHERE m.high_star AND m.cv_class = 'p_lp'
          AND s.med_p < s.med_b
          AND m.score >= s.med_b
    ) AS n_p_above_b,
    count(*) FILTER (
        WHERE m.high_star AND m.cv_class = 'b_lb'
          AND s.med_p < s.med_b
          AND m.score <= s.med_p
    ) AS n_b_below_p
FROM mave_coding m
JOIN (
    SELECT
        side,
        score_set_urn,
        count(*) FILTER (WHERE high_star AND cv_class = 'p_lp') AS n_p,
        count(*) FILTER (WHERE high_star AND cv_class = 'b_lb') AS n_b,
        median(score) FILTER (WHERE high_star AND cv_class = 'p_lp') AS med_p,
        median(score) FILTER (WHERE high_star AND cv_class = 'b_lb') AS med_b
    FROM mave_coding
    GROUP BY 1, 2
) s
    ON m.side = s.side
   AND m.score_set_urn = s.score_set_urn
WHERE s.n_p >= 10
  AND s.n_b >= 10
  AND s.med_p IS DISTINCT FROM s.med_b
GROUP BY m.side, m.score_set_urn;

SELECT 'mave_cv_sets' AS result_set,
    side,
    score_set_urn,
    target_name,
    score_set_title,
    n_p,
    n_b,
    med_p,
    med_b,
    CASE WHEN med_p > med_b THEN n_p_below_b ELSE n_p_above_b END AS n_p_wrong,
    CASE WHEN med_p > med_b THEN n_b_above_p ELSE n_b_below_p END AS n_b_wrong,
    (med_p > med_b) AS plp_has_higher_score
FROM mave_cv_set
ORDER BY side, n_p + n_b DESC;

SELECT 'mave_cv_largest' AS result_set,
    side,
    max(n_p) AS max_n_p,
    max(n_b) AS max_n_b,
    count(*) FILTER (WHERE n_p >= 10 AND n_b >= 10) AS sets_with_10_each
FROM (
    SELECT
        side,
        score_set_urn,
        count(*) FILTER (WHERE high_star AND cv_class = 'p_lp') AS n_p,
        count(*) FILTER (WHERE high_star AND cv_class = 'b_lb') AS n_b
    FROM mave_coding
    GROUP BY 1, 2
)
GROUP BY side
ORDER BY side;

-- CADD phred 20 is a low bar for missense. Show whether the REVEL-benign /
-- CADD-high cell survives stricter CADD cuts. Still missense, both scores present.
SELECT 'revel_cadd_sensitivity' AS result_set,
    side,
    count(*) AS both_scored_missense,
    count(*) FILTER (WHERE revel <= 0.290) AS revel_b,
    median(cadd) FILTER (WHERE revel <= 0.290) AS median_cadd_when_revel_b,
    count(*) FILTER (WHERE revel <= 0.290 AND cadd >= 20) AS revel_b_cadd_ge_20,
    count(*) FILTER (WHERE revel <= 0.290 AND cadd >= 25) AS revel_b_cadd_ge_25,
    count(*) FILTER (WHERE revel <= 0.290 AND cadd >= 30) AS revel_b_cadd_ge_30,
    count(*) FILTER (WHERE revel >= 0.644) AS revel_p,
    median(cadd) FILTER (WHERE revel >= 0.644) AS median_cadd_when_revel_p,
    count(*) FILTER (WHERE revel >= 0.644 AND cadd < 20) AS revel_p_cadd_lt_20,
    count(*) FILTER (WHERE revel >= 0.644 AND cadd < 15) AS revel_p_cadd_lt_15,
    count(*) FILTER (WHERE revel >= 0.644 AND cadd < 10) AS revel_p_cadd_lt_10
FROM coding_cv
WHERE missense AND revel IS NOT NULL AND cadd IS NOT NULL
GROUP BY side
ORDER BY side;

SELECT 'example_plp_revel_b' AS result_set,
    side, variant_key, symbol, consequence, cadd, revel, spliceai, review_stars, review_status, clinical_significance
FROM coding_cv
WHERE missense
  AND high_star
  AND cv_class = 'p_lp'
  AND revel <= 0.290
  AND clinvar_n_rows = 1
QUALIFY row_number() OVER (
    PARTITION BY side
    ORDER BY revel, cadd DESC, variant_key
) <= 5
ORDER BY side, revel, cadd DESC;

SELECT 'example_plp_predictors_benign' AS result_set,
    side, variant_key, symbol, clinvar_symbol, consequence,
    cadd, revel, spliceai, review_stars, review_status, clinical_significance
FROM coding_cv
WHERE missense
  AND high_star
  AND cv_class = 'p_lp'
  AND revel <= 0.290
  AND cadd < 10
  AND clinvar_n_rows = 1
QUALIFY row_number() OVER (
    PARTITION BY side
    ORDER BY revel, cadd, variant_key
) <= 5
ORDER BY side, revel, cadd, variant_key;

SELECT 'example_blb_predictors_path' AS result_set,
    side, variant_key, symbol, clinvar_symbol, consequence,
    cadd, revel, spliceai, review_stars, review_status, clinical_significance
FROM coding_cv
WHERE missense
  AND high_star
  AND cv_class = 'b_lb'
  AND revel >= 0.644
  AND cadd >= 20
  AND clinvar_n_rows = 1
QUALIFY row_number() OVER (
    PARTITION BY side
    ORDER BY revel DESC, cadd DESC, variant_key
) <= 5
ORDER BY side, revel DESC, cadd DESC, variant_key;

SELECT 'example_revel_p_cadd_low' AS result_set,
    side, variant_key, symbol, consequence, cadd, revel, spliceai, cv_class, review_stars
FROM coding_cv
WHERE missense
  AND revel >= 0.644
  AND cadd < 10
QUALIFY row_number() OVER (
    PARTITION BY side
    ORDER BY revel DESC, cadd, variant_key
) <= 3
ORDER BY side, revel DESC, cadd;

SELECT 'example_splice_revel_b' AS result_set,
    side, variant_key, symbol, consequence, cadd, revel, spliceai, cv_class, review_stars, high_star
FROM coding_cv
WHERE missense
  AND spliceai >= 0.5
  AND revel <= 0.290
QUALIFY row_number() OVER (
    PARTITION BY side
    ORDER BY spliceai DESC, revel, variant_key
) <= 3
ORDER BY side, spliceai DESC, revel;
