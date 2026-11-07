# Ancestry audit

## Question

For Tier A variants that also came from GWAS, what fraction have a non-null gnomAD genome `AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, and `AF_amr`, and what fraction are missing each?

## Filters

The set is rows in `warehouse/variants/tier_a.parquet` whose `source` list contains `gwas`. That label is an rsID match. It is not a `variant_key` match. `warehouse/gwas/associations.parquet` has 1,205,644 rows, and `variant_key` is null on every row, because the catalog has no ref or alt. Each of the 13,698 keys has one row in `warehouse/clinvar/allele_dedup.parquet`, the rsID on that row is non-null, and the rsID occurs in the association file. Those keys use 12,902 distinct rsIDs.

13,107 keys have an association on the same chromosome and position. 430 keys are indels whose nearest catalog position on that chromosome is 1 bp away. No SNV in this set is off by one base. 158 keys match an rsID whose catalog position is null and that has no positioned row on another chromosome. 3 keys match an rsID only on the other sex chromosome: `chrX:1294784:C:G` (`rs28602228`) and `chrX:2748343:G:C` (`rs311103`) are on chrY in the catalog, and `chrY:1304146:T:C` (`rs78468042`) is on chrX. The frequencies below are the gnomAD genome frequencies of the Tier A allele.

The columns are the genome fields `AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, and `AF_amr`. In the gnomAD v4.1 genome sites header those names are the African/African-American, Non-Finnish European, East Asian, South Asian, and Admixed American genetic ancestry groups. Exome columns are separate. `AF_afr_exome` and the other `AF_*_exome` fields are not added to the genome fields, and `AN_exome` is not added to `AN`. Exome values are not in these fractions.

A stored 0 is non-null. No ancestry AF in this set is NaN (0 rows). Genome `FILTER = PASS` is not required. Of the keys that pass the allele-number filter, 66 have another genome FILTER: `AS_VQSR` 46, `InbreedingCoeff` 13, `AC0` 4, `AS_VQSR;InbreedingCoeff` 2, and `AC0;AS_VQSR` 1.

The allele-number filter is `in_genomes` and global genome `AN >= 1000`. One thousand is a floor on alleles in the whole genome callset. It is not a floor inside each ancestry group. The genome VCF header defines `AN_afr`, `AN_nfe`, `AN_eas`, `AN_sas`, and `AN_amr`, and `tier_a.parquet` does not contain those columns, so a per-group allele number was not applied.

Unfiltered counts drop that allele-number rule. They keep every Tier A GWAS key, including keys absent from the genome sites VCF.

The queries are `warehouse/reports/sql/ancestry_audit.sql`. A fraction is the count divided by that table's denominator, rounded to 6 decimal places. `1 / 13111` is `0.000076`.

## Denominator

Unfiltered: 13,698 distinct `variant_key` values.

Filtered: 13,111 of those keys. The filter leaves out 587 keys. 584 have `in_genomes` false. 3 are in the genome callset with global `AN` of 2, 146, and 306.

## Result

With `in_genomes` and global `AN >= 1000`, the denominator is 13,111. For one column, missing means that column is null. For the all-five row, missing means at least one of the five is null.

| column | non-null | missing | non-null fraction | missing fraction |
| --- | ---: | ---: | ---: | ---: |
| AF_afr | 13111 | 0 | 1.000000 | 0.000000 |
| AF_nfe | 13111 | 0 | 1.000000 | 0.000000 |
| AF_eas | 13111 | 0 | 1.000000 | 0.000000 |
| AF_sas | 13110 | 1 | 0.999924 | 0.000076 |
| AF_amr | 13111 | 0 | 1.000000 | 0.000000 |
| all five | 13110 | 1 | 0.999924 | 0.000076 |

The one gap inside the filter is `chr3:48565086:C:T` (global `AN` 2304, genome FILTER `PASS`). `AF_sas` is null. `AF_afr`, `AF_nfe`, `AF_eas`, and `AF_amr` are non-null.

Without the allele-number filter, the denominator is 13,698.

| column | non-null | missing | non-null fraction | missing fraction |
| --- | ---: | ---: | ---: | ---: |
| AF_afr | 13113 | 585 | 0.957293 | 0.042707 |
| AF_nfe | 13114 | 584 | 0.957366 | 0.042634 |
| AF_eas | 13113 | 585 | 0.957293 | 0.042707 |
| AF_sas | 13111 | 587 | 0.957147 | 0.042853 |
| AF_amr | 13113 | 585 | 0.957293 | 0.042707 |
| all five | 13111 | 587 | 0.957147 | 0.042853 |

584 keys are missing every one of the five columns. Those are the 584 keys with `in_genomes` false. Three keys are in the genome callset and are missing only some columns. Allele-frequency values are not shown.

| variant_key | genome FILTER | AN | passes AN >= 1000 | null columns |
| --- | --- | ---: | --- | --- |
| chr15:88855936:C:T | AC0;AS_VQSR | 2 | no | AF_afr, AF_eas, AF_sas, AF_amr |
| chr16:88731880:C:G | PASS | 146 | no | AF_sas |
| chr3:48565086:C:T | PASS | 2304 | yes | AF_sas |

The third key with global `AN < 1000` is `chr2:151593174:A:G` (`AN` 306, FILTER `AS_VQSR`). All five columns are non-null, so it is not in that table. The unfiltered all-five missing count, 587, is keys missing at least one column. It is not the 584 keys missing every column.

## Caveat

Most catalog discoveries were made in European-ancestry samples. A frequency in the discovery group is not a frequency for everyone. These fractions are coverage of stored gnomAD genome frequencies on the Tier A alleles above. They are not disease rates. The columns stay in the order `AF_afr`, `AF_nfe`, `AF_eas`, `AF_sas`, `AF_amr`. A lead SNP's frequency in one group is not a risk for that group, and it is not a frequency for any other group. A non-null cell does not mean that group had enough alleles for a stable frequency. This table has no per-group allele number. Global `AN` mixes the whole genome callset.
