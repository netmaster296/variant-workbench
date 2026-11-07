-- Study redundancy: how far GWAS Catalog association rows shrink after
-- trait-and-distance clumping, and which traits pile those rows into a few loci.
--
-- Run with DuckDB from any working directory. The file is several statements.
-- The DuckDB CLI runs the whole file. The Python client runs one statement per execute() call.
-- Numbers in warehouse/reports/study_redundancy.md come from this file.
--
-- Paths:
--   D:/variant-workbench/warehouse/gwas/associations.parquet
--   D:/variant-workbench/warehouse/gwas/loci.parquet
--   D:/variant-workbench/warehouse/gwas/traits.parquet
--
-- Clumps are the rows already stored in loci.parquet. This file does not rebuild them.
-- Eligible association: canonical chrom, pos >= 1, and CNV other than Y.
-- A row with several MAPPED_TRAIT_URI values is one membership per URI.
-- k50 is the smallest number of a trait's loci, largest n_associations first,
-- whose memberships reach half of that trait's clump memberships.
-- Ties break by chrom, lead_pos, then lead_rsid (nulls first).

CREATE OR REPLACE TEMP TABLE assoc AS
SELECT
    count(*) AS n_rows,
    count(*) FILTER (WHERE pos IS NULL) AS null_pos,
    count(*) FILTER (WHERE chrom IS NULL) AS null_chrom,
    count(*) FILTER (WHERE pos IS NOT NULL AND pos < 1) AS pos_below_1,
    count(*) FILTER (WHERE upper(trim(CNV)) = 'Y') AS n_cnv,
    count(*) FILTER (
        WHERE chrom IS NOT NULL
          AND pos IS NOT NULL
          AND pos >= 1
          AND coalesce(upper(trim(CNV)), '') <> 'Y'
    ) AS n_eligible
FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet');

CREATE OR REPLACE TEMP TABLE eligible AS
SELECT
    chrom,
    pos,
    "STUDY ACCESSION" AS study,
    PUBMEDID AS pubmed,
    "DISEASE/TRAIT" AS author_trait,
    list_filter(
        list_transform(
            string_split(coalesce("MAPPED_TRAIT_URI", ''), ','),
            lambda x: trim(x)
        ),
        lambda x: x <> ''
    ) AS uris
FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')
WHERE chrom IS NOT NULL
  AND pos IS NOT NULL
  AND pos >= 1
  AND coalesce(upper(trim(CNV)), '') <> 'Y';

CREATE OR REPLACE TEMP TABLE uri_slots AS
SELECT
    count(*) AS n_rows,
    count(*) FILTER (WHERE n_uri = 0) AS empty_uri_rows,
    count(*) FILTER (WHERE n_uri = 1) AS one_uri_rows,
    count(*) FILTER (WHERE n_uri > 1) AS multi_uri_rows,
    count(*) FILTER (WHERE eligible AND n_uri = 0) AS eligible_empty,
    count(*) FILTER (WHERE eligible AND n_uri > 1) AS eligible_multi,
    sum(CASE WHEN eligible AND n_uri > 1 THEN n_uri - 1 ELSE 0 END) AS extra_eligible_slots,
    sum(CASE WHEN eligible THEN slots ELSE 0 END) AS eligible_slots,
    sum(CASE WHEN NOT eligible THEN slots ELSE 0 END) AS unplaced_slots,
    sum(slots) AS all_slots
FROM (
    SELECT
        len(uris) AS n_uri,
        true AS eligible,
        CASE WHEN len(uris) = 0 THEN 1 ELSE len(uris) END AS slots
    FROM eligible
    UNION ALL
    SELECT
        len(uris) AS n_uri,
        false AS eligible,
        CASE WHEN len(uris) = 0 THEN 1 ELSE len(uris) END AS slots
    FROM (
        SELECT list_filter(
            list_transform(
                string_split(coalesce("MAPPED_TRAIT_URI", ''), ','),
                lambda x: trim(x)
            ),
            lambda x: x <> ''
        ) AS uris
        FROM read_parquet('D:/variant-workbench/warehouse/gwas/associations.parquet')
        WHERE NOT (
            chrom IS NOT NULL
            AND pos IS NOT NULL
            AND pos >= 1
            AND coalesce(upper(trim(CNV)), '') <> 'Y'
        )
    )
);

CREATE OR REPLACE TEMP TABLE loci AS
SELECT
    trait_uri,
    chrom,
    lead_pos,
    lead_rsid,
    lead_p,
    n_associations,
    start,
    "end",
    mapped_genes.raw AS lead_mapped_gene,
    row_number() OVER (
        PARTITION BY trait_uri
        ORDER BY n_associations DESC, chrom, lead_pos, lead_rsid NULLS FIRST
    ) AS rk,
    sum(n_associations) OVER (PARTITION BY trait_uri) AS memberships,
    count(*) OVER (PARTITION BY trait_uri) AS n_loci
FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet');

CREATE OR REPLACE TEMP TABLE trait_half AS
SELECT trait_uri, min(rk) AS k50
FROM (
    SELECT
        trait_uri,
        rk,
        memberships,
        sum(n_associations) OVER (
            PARTITION BY trait_uri
            ORDER BY rk
            ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
        ) AS cum_n
    FROM loci
)
WHERE cum_n * 1.0 / memberships >= 0.5
GROUP BY trait_uri;

CREATE OR REPLACE TEMP TABLE trait_roll AS
SELECT
    l.trait_uri,
    any_value(l.n_loci) AS n_loci,
    any_value(l.memberships) AS memberships,
    max(CASE WHEN l.rk = 1 THEN l.n_associations END) AS top1,
    max(CASE WHEN l.rk = 1 THEN l.chrom END) AS top1_chrom,
    max(CASE WHEN l.rk = 1 THEN l.lead_pos END) AS top1_pos,
    max(CASE WHEN l.rk = 1 THEN l.lead_rsid END) AS top1_rsid,
    max(CASE WHEN l.rk = 1 THEN l.lead_mapped_gene END) AS top1_gene,
    max(CASE WHEN l.rk = 1 THEN l.start END) AS top1_start,
    max(CASE WHEN l.rk = 1 THEN l."end" END) AS top1_end,
    sum(CASE WHEN l.rk <= 5 THEN l.n_associations ELSE 0 END) AS top5
FROM loci l
GROUP BY l.trait_uri;

-- A tie for the largest n_associations would make one top1 row the wrong summary.
-- n_tied_for_top counts loci that share the winning n_associations.
CREATE OR REPLACE TEMP TABLE named AS
SELECT
    tr.mapped_trait,
    tr.trait_uri,
    tr.n_associations AS catalog_rows,
    tr.n_author_traits,
    r.n_loci,
    r.memberships,
    h.k50,
    r.top1,
    r.top1 * 1.0 / r.memberships AS top1_share,
    r.top5,
    r.top5 * 1.0 / r.memberships AS top5_share,
    r.top1_chrom,
    r.top1_pos,
    r.top1_rsid,
    r.top1_gene,
    r.top1_start,
    r.top1_end,
    tr.n_associations - r.memberships AS unplaced,
    (
        SELECT count(*)
        FROM trait_roll r2
        JOIN loci l2
          ON l2.trait_uri IS NOT DISTINCT FROM r2.trait_uri
         AND l2.n_associations = r2.top1
        WHERE r2.trait_uri IS NOT DISTINCT FROM r.trait_uri
    ) AS loci_tied_for_top
FROM trait_roll r
JOIN trait_half h
  ON r.trait_uri IS NOT DISTINCT FROM h.trait_uri
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
  ON r.trait_uri IS NOT DISTINCT FROM tr.trait_uri
WHERE r.memberships >= 500
  AND h.k50 <= 5;

-- 1. Shrinkage.
SELECT
    a.n_rows,
    a.null_pos,
    a.null_chrom,
    a.pos_below_1,
    a.n_cnv,
    a.n_eligible,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS n_loci,
    (SELECT sum(n_associations) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS memberships,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') WHERE n_associations = 1) AS singleton_loci,
    (SELECT max(n_associations) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS max_locus_rows,
    u.empty_uri_rows,
    u.one_uri_rows,
    u.multi_uri_rows,
    u.eligible_empty,
    u.eligible_multi,
    u.extra_eligible_slots,
    u.eligible_slots,
    u.unplaced_slots,
    u.all_slots,
    a.null_pos * 1.0 / a.n_rows AS null_pos_fraction,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) * 1.0 / a.n_rows AS loci_per_raw_row,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) * 1.0 / a.n_eligible AS loci_per_eligible_row,
    a.n_eligible * 1.0 / (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS eligible_rows_per_locus,
    (SELECT sum(n_associations) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) * 1.0
        / (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS memberships_per_locus
FROM assoc a
CROSS JOIN uri_slots u;

-- 2. Locus size. Bins partition the 311,160 loci.
SELECT
    CASE
        WHEN n_associations = 1 THEN '1'
        WHEN n_associations BETWEEN 2 AND 4 THEN '2-4'
        WHEN n_associations BETWEEN 5 AND 9 THEN '5-9'
        WHEN n_associations BETWEEN 10 AND 49 THEN '10-49'
        WHEN n_associations BETWEEN 50 AND 99 THEN '50-99'
        WHEN n_associations BETWEEN 100 AND 499 THEN '100-499'
        ELSE '500+'
    END AS n_associations_bin,
    count(*) AS n_loci,
    sum(n_associations) AS memberships
FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')
GROUP BY 1
ORDER BY min(n_associations);

-- 3. Traits with at least 500 clump memberships whose half-coverage needs 5 loci or fewer.
SELECT
    mapped_trait,
    trait_uri,
    n_loci,
    memberships,
    catalog_rows,
    unplaced,
    n_author_traits,
    k50,
    top1,
    round(top1_share, 4) AS top1_share,
    top5,
    round(top5_share, 4) AS top5_share,
    top1_chrom,
    top1_pos,
    top1_rsid,
    top1_gene,
    top1_start,
    top1_end,
    loci_tied_for_top
FROM named
ORDER BY k50, top1_share DESC, memberships DESC;

-- 4. Largest locus of each named trait, counted inside the lead +/- 1,000,000 bp window.
-- window_rows must equal top1. A mismatch means another clump of the same trait
-- contributed rows, and the study counts must not be used.
SELECT
    n.mapped_trait,
    n.k50,
    n.memberships,
    n.top1,
    count(*) AS window_rows,
    count(*) - n.top1 AS window_minus_top1,
    count(DISTINCT e.study) AS n_study_accessions,
    count(DISTINCT e.pubmed) AS n_pubmed,
    count(DISTINCT e.pos) AS n_positions,
    count(DISTINCT e.author_trait) AS n_author_traits_in_locus,
    count(*) FILTER (WHERE e.study IS NULL OR e.study = '') AS blank_study_rows,
    count(*) FILTER (WHERE e.pubmed IS NULL OR e.pubmed = '') AS blank_pubmed_rows
FROM named n
JOIN eligible e
  ON e.chrom = n.top1_chrom
 AND e.pos BETWEEN n.top1_pos - 1000000 AND n.top1_pos + 1000000
 AND list_contains(e.uris, n.trait_uri)
GROUP BY n.mapped_trait, n.k50, n.memberships, n.top1, n.top1_share
ORDER BY n.k50, n.top1_share DESC, n.memberships DESC;

-- 5. How common this concentration is, and the large traits that miss the rule.
SELECT
    (SELECT count(*) FROM trait_half) AS traits_with_loci,
    (SELECT count(*) FROM trait_roll WHERE memberships >= 100) AS traits_ge_100,
    (SELECT count(*) FROM trait_roll r JOIN trait_half h ON r.trait_uri IS NOT DISTINCT FROM h.trait_uri WHERE r.memberships >= 100 AND h.k50 <= 5) AS ge_100_k50_le_5,
    (SELECT count(*) FROM trait_roll r JOIN trait_half h ON r.trait_uri IS NOT DISTINCT FROM h.trait_uri WHERE r.memberships >= 100 AND r.memberships < 500 AND h.k50 <= 5) AS from_100_to_499_k50_le_5,
    (SELECT count(*) FROM trait_roll WHERE memberships >= 500) AS traits_ge_500,
    (SELECT median(h.k50) FROM trait_roll r JOIN trait_half h ON r.trait_uri IS NOT DISTINCT FROM h.trait_uri WHERE r.memberships >= 500) AS median_k50_ge_500,
    (SELECT count(*) FROM named) AS named_traits,
    (SELECT sum(memberships) FROM named) AS named_memberships,
    (SELECT sum(top1) FROM named) AS named_top1_rows,
    (SELECT sum(top5) FROM named) AS named_top5_rows,
    (SELECT sum(memberships) FROM named) * 1.0 / (SELECT sum(n_associations) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')) AS named_membership_fraction;

SELECT
    tr.mapped_trait,
    r.n_loci,
    r.memberships,
    h.k50,
    r.top1,
    round(r.top1 * 1.0 / r.memberships, 4) AS top1_share,
    r.top5,
    round(r.top5 * 1.0 / r.memberships, 4) AS top5_share,
    round(r.memberships * 1.0 / r.n_loci, 2) AS rows_per_locus
FROM trait_roll r
JOIN trait_half h
  ON r.trait_uri IS NOT DISTINCT FROM h.trait_uri
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
  ON r.trait_uri IS NOT DISTINCT FROM tr.trait_uri
WHERE r.memberships >= 5000
ORDER BY h.k50, r.memberships DESC;

-- Loci with at least 4,000 rows, so the fatty-acid lead can be placed against the file maximum.
SELECT
    tr.mapped_trait,
    lo.chrom,
    lo.lead_pos,
    lo.lead_rsid,
    lo.n_associations
FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') lo
JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
  ON lo.trait_uri IS NOT DISTINCT FROM tr.trait_uri
WHERE lo.n_associations >= 4000
ORDER BY lo.n_associations DESC, tr.mapped_trait;

-- 6. Join and clump checks. The report uses the study counts only if bad_windows is 0.
WITH lead_gaps AS (
    SELECT
        lead_pos - lag(lead_pos) OVER (
            PARTITION BY trait_uri, chrom
            ORDER BY lead_pos
        ) AS gap
    FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')
),
trait_gap AS (
    SELECT
        tr.trait_uri,
        tr.n_associations AS catalog_rows,
        coalesce(sum(lo.n_associations), 0) AS memberships
    FROM read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
    LEFT JOIN read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') lo
      ON lo.trait_uri IS NOT DISTINCT FROM tr.trait_uri
    GROUP BY tr.trait_uri, tr.n_associations
)
SELECT
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet') lo
      ANTI JOIN read_parquet('D:/variant-workbench/warehouse/gwas/traits.parquet') tr
      ON lo.trait_uri IS NOT DISTINCT FROM tr.trait_uri) AS loci_missing_trait,
    (SELECT count(*) FROM trait_gap WHERE memberships > catalog_rows) AS traits_memberships_exceed_catalog,
    (SELECT count(*) FROM trait_gap WHERE memberships = 0) AS traits_with_no_locus,
    (SELECT sum(catalog_rows) FROM trait_gap WHERE memberships = 0) AS catalog_rows_on_traits_with_no_locus,
    (SELECT count(*) FROM read_parquet('D:/variant-workbench/warehouse/gwas/loci.parquet')
      WHERE lead_pos < start OR lead_pos > "end" OR "end" - start > 2000000 OR n_associations < 1) AS loci_breaking_span_rule,
    (SELECT count(*) FROM lead_gaps WHERE gap IS NOT NULL AND gap <= 1000000) AS lead_pairs_within_1mb,
    (SELECT min(gap) FROM lead_gaps WHERE gap IS NOT NULL) AS min_lead_gap,
    (SELECT count(*) FROM (
        SELECT n.trait_uri
        FROM named n
        JOIN eligible e
          ON e.chrom = n.top1_chrom
         AND e.pos BETWEEN n.top1_pos - 1000000 AND n.top1_pos + 1000000
         AND list_contains(e.uris, n.trait_uri)
        GROUP BY n.trait_uri, n.top1
        HAVING count(*) <> n.top1
    )) AS bad_windows,
    (SELECT count(*) FROM named WHERE loci_tied_for_top <> 1) AS named_traits_with_tied_top_locus;
