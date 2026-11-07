-- Pleiotropy: genes and Tier A variants linked to many EFO traits.
-- Run from the repository root. The DuckDB CLI prints every statement.
-- Python duckdb returns only the last result of a multi-statement string;
-- execute statements separately if you need each table.
--
-- A shared mapped gene is not a shared mechanism.
-- A gene named inside a wide clump is a weaker link than a shared lead.
-- Most ledger GWAS rows have no ref/alt key.

SET preserve_insertion_order = false;

CREATE TEMP TABLE genes AS
SELECT ensg, symbol, biotype, chrom, start, "end"
FROM read_parquet('warehouse/intervals/genes.parquet')
WHERE symbol IS NOT NULL AND symbol <> '';

CREATE TEMP TABLE sym AS
SELECT
  symbol,
  count(*) AS n_ensg,
  arg_min(ensg, ensg) AS ensg,
  arg_min(biotype, ensg) AS biotype
FROM genes
GROUP BY symbol;

CREATE TEMP TABLE loci AS
SELECT
  trait_uri,
  chrom,
  lead_pos,
  lead_rsid,
  lead_p,
  n_associations,
  start,
  "end",
  ("end" - start) AS span,
  mapped_genes
FROM read_parquet('warehouse/gwas/loci.parquet');

CREATE TEMP TABLE traits AS
SELECT trait_uri, mapped_trait, n_associations, n_author_traits
FROM read_parquet('warehouse/gwas/traits.parquet');

-- Catalog tokens. Intergenic "GENE - GENE" phrases stay one token.
CREATE TEMP TABLE member_tok AS
SELECT trait_uri, chrom, lead_pos, lead_rsid, span, token
FROM (
  SELECT
    trait_uri,
    chrom,
    lead_pos,
    lead_rsid,
    span,
    unnest(mapped_genes.tokens) AS token
  FROM loci
)
WHERE trait_uri IS NOT NULL
  AND trait_uri <> ''
  AND token IS NOT NULL
  AND token <> ''
  AND token <> '-';

-- Lead association only. mapped_genes.raw is that row, split on commas.
CREATE TEMP TABLE lead_tok AS
SELECT trait_uri, chrom, lead_pos, lead_rsid, span, trim(tok) AS token
FROM (
  SELECT
    trait_uri,
    chrom,
    lead_pos,
    lead_rsid,
    span,
    unnest(string_split(mapped_genes.raw, ',')) AS tok
  FROM loci
  WHERE mapped_genes.raw IS NOT NULL
)
WHERE trait_uri IS NOT NULL
  AND trait_uri <> ''
  AND trim(tok) <> ''
  AND trim(tok) <> '-';

CREATE TEMP TABLE gene_locus AS
SELECT
  m.trait_uri,
  m.chrom,
  m.lead_pos,
  m.lead_rsid,
  m.span,
  m.token AS symbol,
  EXISTS (
    SELECT 1
    FROM lead_tok d
    WHERE d.trait_uri IS NOT DISTINCT FROM m.trait_uri
      AND d.chrom IS NOT DISTINCT FROM m.chrom
      AND d.lead_pos IS NOT DISTINCT FROM m.lead_pos
      AND d.lead_rsid IS NOT DISTINCT FROM m.lead_rsid
      AND d.token = m.token
  ) AS on_lead
FROM member_tok m
SEMI JOIN sym s ON s.symbol = m.token;

CREATE TEMP TABLE gene_trait AS
SELECT
  g.symbol,
  g.trait_uri,
  any_value(s.n_ensg) AS n_ensg,
  any_value(s.ensg) AS ensg,
  any_value(s.biotype) AS biotype,
  count(*) AS n_loci,
  bool_or(g.span < 1000000) AS any_small,
  bool_or(g.span = 0) AS any_point,
  bool_or(g.span >= 1000000) AS any_large,
  bool_or(g.on_lead) AS any_lead
FROM gene_locus g
JOIN sym s ON s.symbol = g.symbol
GROUP BY g.symbol, g.trait_uri;

-- Same rule as the ledger: exact rsID, same chromosome, same position.
-- A match at another position is not attached. ref = alt is not an allele.
CREATE TEMP TABLE clin AS
SELECT rsid, variant_key, chrom, pos
FROM read_parquet('warehouse/clinvar/allele_dedup.parquet')
WHERE variant_key IS NOT NULL
  AND rsid IS NOT NULL;

CREATE TEMP TABLE lead_var AS
SELECT
  l.trait_uri,
  l.chrom,
  l.lead_pos,
  l.lead_rsid,
  l.span,
  c.variant_key
FROM loci l
JOIN clin c
  ON c.rsid = l.lead_rsid
 AND c.chrom = l.chrom
 AND c.pos = l.lead_pos
WHERE l.trait_uri IS NOT NULL
  AND l.trait_uri <> '';

CREATE TEMP TABLE tier AS
SELECT
  t.variant_key,
  t.ref,
  t.alt,
  t.in_genomes,
  t.in_exomes,
  t.AF,
  t.AF_exome,
  t.vep_symbol,
  t.vep_consequence,
  t.source
FROM read_parquet('warehouse/variants/tier_a.parquet') t
SEMI JOIN (SELECT DISTINCT variant_key FROM lead_var) k
  ON k.variant_key = t.variant_key;

CREATE TEMP TABLE lead_alt AS
SELECT v.trait_uri, v.chrom, v.lead_pos, v.lead_rsid, v.variant_key
FROM lead_var v
JOIN tier t ON t.variant_key = v.variant_key
WHERE t.ref <> t.alt;

-- Checks. missing_tier_keys must be 0 or the Tier A claim stops.
SELECT 'check_ledger' AS label, gwas_rows, ref_alt_keys, no_ref_alt
FROM (
  SELECT
    count(*) AS gwas_rows,
    count(*) FILTER (
      WHERE split_part(subject_id, ':', 4) <> ''
        AND split_part(subject_id, ':', 5) = ''
    ) AS ref_alt_keys,
    count(*) FILTER (WHERE split_part(subject_id, ':', 4) = '') AS no_ref_alt
  FROM read_parquet('warehouse/ledger/evidence.parquet')
  WHERE evidence_class = 'gwas'
);

SELECT 'check_tier_join' AS label,
       (SELECT count(*) FROM (SELECT DISTINCT variant_key FROM lead_var)) AS lead_keys,
       (
         SELECT count(*)
         FROM (SELECT DISTINCT variant_key FROM lead_var) k
         ANTI JOIN tier t ON t.variant_key = k.variant_key
       ) AS missing_tier_keys;

SELECT 'check_other_position' AS label, count(*) AS loci_rsid_other_position_only
FROM loci l
WHERE l.lead_rsid IS NOT NULL
  AND EXISTS (SELECT 1 FROM clin c WHERE c.rsid = l.lead_rsid)
  AND NOT EXISTS (
    SELECT 1 FROM clin c
    WHERE c.rsid = l.lead_rsid
      AND c.chrom = l.chrom
      AND c.pos = l.lead_pos
  );

SELECT 'check_mart_ensg' AS label,
       (
         SELECT count(DISTINCT g.ensg)
         FROM member_tok m
         JOIN genes g ON g.symbol = m.token
       ) AS ensg_from_tokens,
       (
         SELECT count(*)
         FROM read_parquet('warehouse/marts/gene_mart.parquet')
         WHERE gwas_loci_by_mapped_gene > 0
       ) AS mart_mapped_genes;

SELECT 'denominator_loci' AS label,
       count(*) AS loci,
       count(*) FILTER (WHERE trait_uri IS NULL OR trait_uri = '') AS empty_trait_loci,
       count(DISTINCT trait_uri) FILTER (
         WHERE trait_uri IS NOT NULL AND trait_uri <> ''
       ) AS traits,
       count(*) FILTER (
         WHERE trait_uri IS NOT NULL AND trait_uri <> '' AND span = 0
       ) AS point_loci,
       count(*) FILTER (
         WHERE trait_uri IS NOT NULL AND trait_uri <> ''
           AND span > 0 AND span < 1000000
       ) AS mid_loci,
       count(*) FILTER (
         WHERE trait_uri IS NOT NULL AND trait_uri <> '' AND span >= 1000000
       ) AS large_loci
FROM loci;

SELECT 'denominator_tokens' AS label,
       (SELECT count(DISTINCT token) FROM member_tok) AS tokens,
       (SELECT count(*) FROM sym) AS gencode_symbols,
       (SELECT count(*) FROM genes) AS gencode_genes,
       (
         SELECT count(*)
         FROM (SELECT DISTINCT token FROM member_tok) t
         ANTI JOIN sym s ON s.symbol = t.token
       ) AS unmatched_tokens,
       (
         SELECT count(*)
         FROM (SELECT DISTINCT token FROM member_tok) t
         ANTI JOIN sym s ON s.symbol = t.token
         WHERE t.token LIKE '% - %'
       ) AS unmatched_dash_tokens;

SELECT 'denominator_genes' AS label,
       (
         SELECT count(*)
         FROM (SELECT symbol FROM gene_trait GROUP BY symbol)
       ) AS symbols,
       (
         SELECT count(*)
         FROM (
           SELECT symbol
           FROM gene_trait
           GROUP BY symbol
           HAVING max(n_ensg) > 1
         )
       ) AS colliding_symbols,
       (
         SELECT count(DISTINCT g.ensg)
         FROM gene_trait t
         JOIN genes g ON g.symbol = t.symbol
       ) AS ensg_rows,
       (SELECT count(DISTINCT trait_uri) FROM gene_trait) AS traits_with_gene,
       (SELECT count(*) FROM gene_trait) AS gene_trait_edges;

SELECT 'denominator_variants' AS label,
       (SELECT count(*) FROM lead_var) AS keyed_rows_including_ref_eq,
       (SELECT count(DISTINCT (trait_uri, chrom, lead_pos, lead_rsid)) FROM lead_var) AS keyed_loci,
       (SELECT count(*) FROM tier) AS tier_keys,
       (SELECT count(*) FROM tier WHERE ref = alt) AS ref_eq_alt_keys,
       (SELECT count(*) FROM (SELECT DISTINCT variant_key FROM lead_alt)) AS alt_keys,
       (SELECT count(*) FROM lead_alt) AS alt_rows,
       (SELECT count(DISTINCT trait_uri) FROM lead_alt) AS traits_with_alt,
       (SELECT count(DISTINCT lead_rsid) FROM lead_alt) AS alt_rsids;

SELECT 'ref_eq_only_rsids' AS label, count(*) AS n
FROM (
  SELECT v.lead_rsid
  FROM lead_var v
  JOIN tier t ON t.variant_key = v.variant_key
  GROUP BY v.lead_rsid
  HAVING bool_and(t.ref = t.alt)
);

SELECT 'multi_allele_loci' AS label, count(*) AS loci
FROM (
  SELECT trait_uri, chrom, lead_pos, lead_rsid
  FROM lead_var
  GROUP BY trait_uri, chrom, lead_pos, lead_rsid
  HAVING count(DISTINCT variant_key) > 1
);

SELECT 'vep_impact_observed' AS label, impact, n
FROM (
  SELECT coalesce(vep_impact, '(none)') AS impact, count(*) AS n
  FROM read_parquet('warehouse/variants/tier_a.parquet') t
  SEMI JOIN (SELECT DISTINCT variant_key FROM lead_alt) k
    ON k.variant_key = t.variant_key
  WHERE t.ref <> t.alt
    AND (coalesce(t.in_genomes, false) OR coalesce(t.in_exomes, false))
  GROUP BY 1
)
ORDER BY n DESC;

SELECT 'vep_consequence_observed' AS label, consequence, n
FROM (
  SELECT coalesce(t.vep_consequence, '(no vep)') AS consequence, count(*) AS n
  FROM read_parquet('warehouse/variants/tier_a.parquet') t
  SEMI JOIN (SELECT DISTINCT variant_key FROM lead_alt) k
    ON k.variant_key = t.variant_key
  WHERE t.ref <> t.alt
    AND (coalesce(t.in_genomes, false) OR coalesce(t.in_exomes, false))
  GROUP BY 1
)
ORDER BY n DESC
LIMIT 5;

SELECT 'alt_not_in_gnomad' AS label, count(*) AS n
FROM tier
WHERE ref <> alt
  AND NOT coalesce(in_genomes, false)
  AND NOT coalesce(in_exomes, false);

SELECT 'measurement_label' AS label,
       count(*) FILTER (WHERE n_traits >= 50) AS symbols_ge50,
       sum(n_traits) FILTER (WHERE n_traits >= 50) AS trait_labels,
       sum(n_meas) FILTER (WHERE n_traits >= 50) AS labels_with_measurement
FROM (
  SELECT
    g.symbol,
    count(*) AS n_traits,
    count(*) FILTER (WHERE lower(t.mapped_trait) LIKE '%measurement%') AS n_meas
  FROM gene_trait g
  JOIN traits t ON t.trait_uri = g.trait_uri
  GROUP BY g.symbol
);

SELECT 'biotype_ge50' AS label, biotype, count(*) AS symbols
FROM (
  SELECT symbol, any_value(biotype) AS biotype, count(*) AS n_traits
  FROM gene_trait
  GROUP BY symbol
)
WHERE n_traits >= 50
GROUP BY biotype
ORDER BY symbols DESC, biotype;

SELECT 'gene_degree' AS label,
       count(*) AS symbols,
       count(*) FILTER (WHERE n_traits = 1) AS t1,
       count(*) FILTER (WHERE n_traits BETWEEN 2 AND 4) AS t2_4,
       count(*) FILTER (WHERE n_traits BETWEEN 5 AND 9) AS t5_9,
       count(*) FILTER (WHERE n_traits BETWEEN 10 AND 19) AS t10_19,
       count(*) FILTER (WHERE n_traits BETWEEN 20 AND 49) AS t20_49,
       count(*) FILTER (WHERE n_traits BETWEEN 50 AND 99) AS t50_99,
       count(*) FILTER (WHERE n_traits >= 100) AS t100,
       quantile_cont(n_traits, 0.5) AS median_traits,
       max(n_traits) AS max_traits,
       count(*) FILTER (WHERE n_point >= 1) AS symbols_point,
       count(*) FILTER (WHERE n_point >= 10) AS point_ge10,
       count(*) FILTER (WHERE n_point >= 50) AS point_ge50,
       max(n_point) AS max_point
FROM (
  SELECT
    symbol,
    count(*) AS n_traits,
    count(*) FILTER (WHERE any_point) AS n_point
  FROM gene_trait
  GROUP BY symbol
);

SELECT 'gene_edges' AS label,
       count(*) AS edges,
       count(*) FILTER (WHERE any_large AND NOT any_small) AS only_large,
       count(*) FILTER (WHERE any_point) AS point_edges,
       count(*) FILTER (WHERE any_lead) AS lead_field_edges,
       count(*) FILTER (WHERE NOT any_lead) AS member_only_edges
FROM gene_trait;

SELECT 'lead_inside_unique_symbol' AS label,
       count(*) AS edges,
       count(*) FILTER (WHERE inside) AS lead_inside,
       count(*) FILTER (WHERE any_point) AS point_edges,
       count(*) FILTER (WHERE point_inside) AS point_inside
FROM (
  SELECT
    gl.symbol,
    gl.trait_uri,
    bool_or(gl.chrom = g.chrom AND gl.lead_pos BETWEEN g.start AND g."end") AS inside,
    bool_or(gl.span = 0) AS any_point,
    bool_or(
      gl.span = 0
      AND gl.chrom = g.chrom
      AND gl.lead_pos BETWEEN g.start AND g."end"
    ) AS point_inside
  FROM gene_locus gl
  JOIN genes g ON g.symbol = gl.symbol
  WHERE gl.symbol IN (SELECT symbol FROM sym WHERE n_ensg = 1)
  GROUP BY gl.symbol, gl.trait_uri
);

SELECT 'variant_degree' AS label,
       count(*) AS variants,
       count(*) FILTER (WHERE n_traits = 1) AS t1,
       count(*) FILTER (WHERE n_traits BETWEEN 2 AND 4) AS t2_4,
       count(*) FILTER (WHERE n_traits BETWEEN 5 AND 9) AS t5_9,
       count(*) FILTER (WHERE n_traits BETWEEN 10 AND 19) AS t10_19,
       count(*) FILTER (WHERE n_traits BETWEEN 20 AND 49) AS t20_49,
       count(*) FILTER (WHERE n_traits >= 50) AS t50,
       quantile_cont(n_traits, 0.5) AS median_traits,
       max(n_traits) AS max_traits
FROM (
  SELECT variant_key, count(DISTINCT trait_uri) AS n_traits
  FROM lead_alt
  GROUP BY variant_key
);

SELECT 'multi_key_rsids' AS label,
       count(*) AS rsids,
       count(*) FILTER (WHERE n_sets = 1) AS identical_trait_sets,
       max(n_keys) AS max_keys
FROM (
  SELECT lead_rsid, count(*) AS n_keys, count(DISTINCT traits) AS n_sets
  FROM (
    SELECT
      lead_rsid,
      variant_key,
      string_agg(trait_uri, '|' ORDER BY trait_uri) AS traits
    FROM (
      SELECT DISTINCT lead_rsid, variant_key, trait_uri
      FROM lead_alt
    )
    GROUP BY lead_rsid, variant_key
  )
  GROUP BY lead_rsid
  HAVING count(*) > 1
);

SELECT 'site_degree' AS label,
       count(*) AS sites,
       count(*) FILTER (WHERE n_traits >= 2) AS ge2,
       count(*) FILTER (WHERE n_traits >= 10) AS ge10,
       count(*) FILTER (WHERE n_traits >= 50) AS ge50,
       max(n_traits) AS max_traits
FROM (
  SELECT chrom, lead_pos, lead_rsid, count(DISTINCT trait_uri) AS n_traits
  FROM loci
  WHERE trait_uri IS NOT NULL
    AND trait_uri <> ''
    AND lead_rsid IS NOT NULL
  GROUP BY chrom, lead_pos, lead_rsid
);

SELECT 'sites_with_alt_key' AS label,
       count(*) FILTER (WHERE n_traits >= 2) AS sites_ge2,
       count(*) FILTER (WHERE n_traits >= 2 AND keyed) AS keyed_ge2,
       count(*) FILTER (WHERE n_traits >= 10) AS sites_ge10,
       count(*) FILTER (WHERE n_traits >= 10 AND keyed) AS keyed_ge10,
       count(*) FILTER (WHERE n_traits >= 50) AS sites_ge50,
       count(*) FILTER (WHERE n_traits >= 50 AND keyed) AS keyed_ge50
FROM (
  SELECT
    s.n_traits,
    EXISTS (
      SELECT 1
      FROM lead_alt a
      WHERE a.chrom = s.chrom
        AND a.lead_pos = s.lead_pos
        AND a.lead_rsid = s.lead_rsid
    ) AS keyed
  FROM (
    SELECT chrom, lead_pos, lead_rsid, count(DISTINCT trait_uri) AS n_traits
    FROM loci
    WHERE trait_uri IS NOT NULL
      AND trait_uri <> ''
      AND lead_rsid IS NOT NULL
    GROUP BY chrom, lead_pos, lead_rsid
  ) s
);

SELECT 'top_site' AS label, n_traits, lead_rsid, chrom, lead_pos, keyed
FROM (
  SELECT
    s.chrom,
    s.lead_pos,
    s.lead_rsid,
    s.n_traits,
    EXISTS (
      SELECT 1
      FROM lead_alt a
      WHERE a.chrom = s.chrom
        AND a.lead_pos = s.lead_pos
        AND a.lead_rsid = s.lead_rsid
    ) AS keyed
  FROM (
    SELECT chrom, lead_pos, lead_rsid, count(DISTINCT trait_uri) AS n_traits
    FROM loci
    WHERE trait_uri IS NOT NULL
      AND trait_uri <> ''
      AND lead_rsid IS NOT NULL
    GROUP BY chrom, lead_pos, lead_rsid
  ) s
  ORDER BY s.n_traits DESC, s.chrom, s.lead_pos
  LIMIT 1
);

SELECT 'top_genes' AS label,
       symbol, ensg, biotype, n_traits, n_small, n_point, n_lead
FROM (
  SELECT
    symbol,
    any_value(ensg) AS ensg,
    any_value(biotype) AS biotype,
    count(*) AS n_traits,
    count(*) FILTER (WHERE any_small) AS n_small,
    count(*) FILTER (WHERE any_point) AS n_point,
    count(*) FILTER (WHERE any_lead) AS n_lead
  FROM gene_trait
  GROUP BY symbol
)
ORDER BY n_traits DESC, symbol
LIMIT 15;

SELECT 'top_variants' AS label,
       d.n_traits, d.rsid, d.variant_key, t.vep_symbol, t.vep_consequence,
       t.in_genomes, t.in_exomes, t.AF
FROM (
  SELECT
    v.variant_key,
    count(DISTINCT v.trait_uri) AS n_traits,
    any_value(v.lead_rsid) AS rsid
  FROM lead_alt v
  GROUP BY v.variant_key
) d
JOIN tier t ON t.variant_key = d.variant_key
ORDER BY d.n_traits DESC, d.variant_key
LIMIT 16;

CREATE TEMP TABLE gene_pairs AS
SELECT
  a.trait_uri AS trait_a,
  b.trait_uri AS trait_b,
  count(*) AS n_shared_genes,
  count(*) FILTER (WHERE a.any_small AND b.any_small) AS n_shared_genes_small,
  count(*) FILTER (WHERE a.any_point AND b.any_point) AS n_shared_genes_point
FROM gene_trait a
JOIN gene_trait b
  ON a.symbol = b.symbol
 AND a.trait_uri < b.trait_uri
GROUP BY a.trait_uri, b.trait_uri;

CREATE TEMP TABLE var_pairs AS
SELECT
  a.trait_uri AS trait_a,
  b.trait_uri AS trait_b,
  count(DISTINCT a.variant_key) AS n_shared_variants,
  count(DISTINCT a.lead_rsid) AS n_shared_rsids
FROM lead_alt a
JOIN lead_alt b
  ON a.variant_key = b.variant_key
 AND a.trait_uri < b.trait_uri
GROUP BY a.trait_uri, b.trait_uri;

CREATE TEMP TABLE site_pairs AS
SELECT
  a.trait_uri AS trait_a,
  b.trait_uri AS trait_b,
  count(*) AS n_shared_sites
FROM (
  SELECT DISTINCT trait_uri, chrom, lead_pos, lead_rsid
  FROM loci
  WHERE trait_uri IS NOT NULL
    AND trait_uri <> ''
    AND lead_rsid IS NOT NULL
) a
JOIN (
  SELECT DISTINCT trait_uri, chrom, lead_pos, lead_rsid
  FROM loci
  WHERE trait_uri IS NOT NULL
    AND trait_uri <> ''
    AND lead_rsid IS NOT NULL
) b
  ON a.chrom = b.chrom
 AND a.lead_pos = b.lead_pos
 AND a.lead_rsid = b.lead_rsid
 AND a.trait_uri < b.trait_uri
GROUP BY a.trait_uri, b.trait_uri;

SELECT 'pair_summary' AS label,
       (SELECT count(*) FROM gene_pairs) AS gene_pairs,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes_small > 0) AS gene_pairs_small,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes_point > 0) AS gene_pairs_point,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes_small = 0) AS gene_pairs_only_large,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes >= 5) AS genes_ge5,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes >= 10) AS genes_ge10,
       (SELECT count(*) FROM gene_pairs WHERE n_shared_genes >= 50) AS genes_ge50,
       (SELECT max(n_shared_genes) FROM gene_pairs) AS max_genes,
       (SELECT max(n_shared_genes_point) FROM gene_pairs) AS max_point_genes,
       (SELECT count(*) FROM site_pairs) AS site_pairs,
       (SELECT count(*) FROM site_pairs WHERE n_shared_sites >= 10) AS sites_ge10,
       (SELECT max(n_shared_sites) FROM site_pairs) AS max_sites,
       (SELECT count(*) FROM var_pairs) AS variant_pairs,
       (SELECT count(*) FROM var_pairs WHERE n_shared_variants >= 2) AS variants_ge2,
       (SELECT count(*) FROM var_pairs WHERE n_shared_variants >= 10) AS variants_ge10,
       (SELECT max(n_shared_variants) FROM var_pairs) AS max_variants,
       (
         SELECT count(*)
         FROM var_pairs
         WHERE n_shared_variants > n_shared_rsids
       ) AS pairs_more_keys_than_rsids,
       (
         SELECT count(*) FROM site_pairs
       ) - (
         SELECT count(*) FROM var_pairs
       ) AS site_pairs_without_alt,
       (
         SELECT count(*)
         FROM gene_pairs g
         SEMI JOIN site_pairs s
           ON s.trait_a = g.trait_a AND s.trait_b = g.trait_b
       ) AS gene_pairs_with_site,
       (
         SELECT count(*)
         FROM gene_pairs g
         SEMI JOIN var_pairs v
           ON v.trait_a = g.trait_a AND v.trait_b = g.trait_b
       ) AS gene_pairs_with_variant,
       (
         SELECT count(*)
         FROM var_pairs v
         ANTI JOIN gene_pairs g
           ON g.trait_a = v.trait_a AND g.trait_b = v.trait_b
       ) AS variant_pairs_without_gene,
       (
         SELECT count(*)
         FROM site_pairs s
         ANTI JOIN gene_pairs g
           ON g.trait_a = s.trait_a AND g.trait_b = s.trait_b
       ) AS site_pairs_without_gene;

SELECT 'top_gene_pairs' AS label,
       g.n_shared_genes,
       g.n_shared_genes_small,
       g.n_shared_genes_point,
       coalesce(s.n_shared_sites, 0) AS n_shared_sites,
       coalesce(v.n_shared_variants, 0) AS n_shared_variants,
       ta.mapped_trait AS trait_a,
       tb.mapped_trait AS trait_b
FROM gene_pairs g
LEFT JOIN site_pairs s
  ON s.trait_a = g.trait_a AND s.trait_b = g.trait_b
LEFT JOIN var_pairs v
  ON v.trait_a = g.trait_a AND v.trait_b = g.trait_b
JOIN traits ta ON ta.trait_uri = g.trait_a
JOIN traits tb ON tb.trait_uri = g.trait_b
ORDER BY g.n_shared_genes DESC, ta.mapped_trait, tb.mapped_trait
LIMIT 10;

SELECT 'top_variant_pairs' AS label,
       v.n_shared_variants,
       v.n_shared_rsids,
       coalesce(s.n_shared_sites, 0) AS n_shared_sites,
       coalesce(g.n_shared_genes, 0) AS n_shared_genes,
       coalesce(g.n_shared_genes_small, 0) AS n_shared_genes_small,
       coalesce(g.n_shared_genes_point, 0) AS n_shared_genes_point,
       ta.mapped_trait AS trait_a,
       tb.mapped_trait AS trait_b
FROM var_pairs v
LEFT JOIN site_pairs s
  ON s.trait_a = v.trait_a AND s.trait_b = v.trait_b
LEFT JOIN gene_pairs g
  ON g.trait_a = v.trait_a AND g.trait_b = v.trait_b
JOIN traits ta ON ta.trait_uri = v.trait_a
JOIN traits tb ON tb.trait_uri = v.trait_b
ORDER BY v.n_shared_variants DESC, v.n_shared_rsids DESC, ta.mapped_trait, tb.mapped_trait
LIMIT 10;
