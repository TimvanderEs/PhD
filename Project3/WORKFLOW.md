# HELIOS 10k cognitive GWAS workflow

This repository documents the final Chinese 10k GWAS freeze used in the
downstream Project 2 analyses.

## Workflow

| Step | Record |
|---|---|
| Cognitive factor | Reproduced from the historical script and saved score for 7,403 participants |
| Genotype QC | Portable version of the recovered PLINK commands and retained output naming |
| REGENIE Step 1 | Completed REGENIE v4.1 log and recovered command |
| REGENIE Step 2 | Completed REGENIE v4.1 log and recovered command |
| Summary-statistic preparation | Recovered build conversion and rsID-mapping scripts |
| Final freeze | EC2 checksums, row counts and cross-format comparisons |
| SG100K LD reference | Recovered scripts and completed chromosome logs |

## Cognitive phenotype

Six cognitive tasks were directionally aligned, residualised and combined by
principal-component analysis. The reproduced participant set exactly matches
the 7,403 IDs in the historical score file. PC1 explained 34.489% of the total
variance.

The stored score has negative coefficients for the direction-aligned tasks, so
higher stored `g` values indicate poorer performance. This is an arbitrary PCA
sign. The aggregate tables also provide the sign-equivalent
higher-performance orientation for presentation.

## Chinese association test

The completed 28 November 2025 REGENIE v4.1 Step 1 log records 5,700
participants and 279,546 LD-pruned variants. It produced the
`CHI_10k_step1_pred.list` prediction list used by Step 2. Step 1 used the
following phenotype and covariates:

- phenotype column `g_raw`;
- age, sex, age squared, age-by-sex, age-squared-by-sex and 20 genetic
  principal components; and
- categorical treatment of sex.

The completed Step 2 log from the same date records:

- phenotype column `g_raw`;
- age, sex, age squared, age-by-sex, age-squared-by-sex and 20 genetic
  principal components;
- categorical treatment of sex;
- `--apply-rint`, `--minMAC 0.5` and block size 100;
- 5,700 participants entering Step 2; and
- 7,043,723 variants in the genotype input.

The Step 1 and Step 2 logs use the same 25 covariate columns and report
covariate data for 6,353 participants.

The modal variant-level sample size was 5,592. The corresponding post-QC LDSC
input had mean chi-square 1.063.

## Summary-statistic freeze

The retained build-37 Chinese summary-statistic file,
`HELIOS_10k_lifted_clean_with_stats_FINAL.txt.gz`, contains 5,615,100 rows. All
rows passed checks for field count, chromosome, position, allele coding, sample
size, frequency, P value and consistency of Z with beta/SE. The modal
variant-level sample size is 5,592.

The complete downstream position-based table contains 4,875,695 rows. Every
coordinate, allele and statistic combination was matched to the final freeze.
The LAVA and MiXeR exports have identical normalized chromosome, position,
allele, Z and N fields. The 4,867,848-row MTAG and PLEIO exports likewise have
identical normalized rsID, allele, Z and N fields. Exact file fingerprints are
listed in [`results/summary_statistics/freeze_manifest.tsv`](results/summary_statistics/freeze_manifest.tsv).

Two earlier Chinese–Indian–Malay METAL outputs were also checked. They contain
6,686,830 and 6,688,950 variants, with maximum sample sizes of 7,789 and 7,664.
These are historical analyses of earlier cognitive-phenotype versions and are
not inputs to the retained Chinese-only freeze.

## Downstream analyses

The April 2026 `10k_clean` LDSC analysis used an East Asian LD-score reference.
It reported observed-scale heritability of 0.3048 (SE 0.0949) and genetic
correlation with East Asian educational attainment of -0.0502 (SE 0.0835,
P=0.5476).

The successful MTAG v1.0.8 run completed on 2 December 2025 and combined the
10k HELIOS result with East Asian educational attainment. It retained
3,403,649 shared variants and is currently classified as exploratory.

The repository contains code, configuration, aggregate results and
non-identifying provenance. Participant-level HELIOS and SG100K data remain
controlled.
