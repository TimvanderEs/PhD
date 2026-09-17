# HELIOS cognitive GWAS workflow

This repository documents the HELIOS cognitive GWAS. Its main analysis is the
ancestry-stratified GWAS in Chinese, Indian and Malay participants followed by
trans-ancestry meta-analysis. The original 10k analysis, expanded 22k
validation and a later Chinese-only freeze used in Project 2 are kept as
separate analysis records.

## Workflow

| Step | Record |
|---|---|
| Cognitive factor | Reproduced from the historical script and saved score for 7,403 participants |
| Cognitive validation | Measurement summaries, 1,000 bootstrap PCAs, leave-one-task-out and ancestry-specific checks |
| Expanded-cohort validation | Exact 12k and combined 22k PCA checks, residual reconstruction, measurement and structural sensitivity analyses |
| Genotype QC | Portable ancestry-specific PLINK workflow for Chinese, Indian and Malay participants |
| Ancestry-specific GWAS | REGENIE v4.1 commands and recovered Step 2 run records |
| Trans-ancestry meta-analysis | METAL inverse-variance model and the exact post-analysis filters |
| Summary-statistic preparation | Recovered build conversion and rsID-mapping scripts |
| Chinese-only freeze | EC2 checksums, row counts and cross-format comparisons |
| SG100K LD reference | Recovered scripts and completed chromosome logs |

## Cognitive phenotype

Six cognitive tasks were directionally aligned, residualised and combined by
principal-component analysis. The reproduced participant set exactly matches
the 7,403 IDs in the historical score file. PC1 explained 34.489% of the total
variance. This is one of the 10k phenotype iterations; the retained
meta-analysis outputs are identified separately by checksum because the
phenotype was updated during the project.

The stored score has negative coefficients for the direction-aligned tasks, so
higher stored `g` values indicate poorer performance. This is an arbitrary PCA
sign. The aggregate tables also provide the sign-equivalent
higher-performance orientation for presentation.

## Expanded 22k validation

The archived 12k extension PCA was reproduced in 13,774 participants. The
combined phenotype residualisation and PCA were reproduced in 22,422
participants; PC1 explained 37.04% of variance. All validation outputs use the
sign-equivalent orientation in which higher values indicate better cognitive
performance.

The retained task residualisation model included age, sex, age squared,
age-by-sex terms, ancestry indicators and education as a categorical
covariate. Omitting education preserved the PC1 loading pattern (Tucker
congruence 0.9999), although the resulting participant scores correlated
r=0.981 with the retained score. The genetic effect of this choice requires a
GWAS sensitivity analysis using the no-education phenotype.

The three latency indicators remained the strongest PC1 contributors in 1,000
bootstrap replicates. Loadings were consistent across Chinese, Indian and
Malay participants and across the HELIOS10K and HELIOS20K freezes. The
four-indicator PCA with one collapsed speed domain correlated r=0.931 with the
six-task PC1. Full results are in [`validation/22k`](validation/22k).

The ancestry-specific expanded GWAS used 17,966 Chinese, 2,310 Indian and
1,242 Malay participants. Its LDSC record reported h2=0.0392 (SE 0.0202),
intercept=1.0008 (SE 0.0058), mean chi-square=1.0175 and lambda GC=1.0245.

## Trans-ancestry association analysis

GWAS was run separately in the Chinese, Indian and Malay ancestry groups using
REGENIE v4.1. The recovered 10k Step 2 records use the ancestry-specific
genotype data, a precomputed `g` phenotype, 20 ancestry-specific principal
components, block size 200 and minimum MAC 0.5. The genotype inputs contained
7,043,723 Chinese, 8,387,428 Indian and 7,660,514 Malay variants.

The three association results were combined with METAL using the standard-error
scheme. The retained variants were present in at least two ancestry groups, had
heterogeneity P > 0.05 and had a maximum between-ancestry allele-frequency
difference < 0.5. The executable workflow is in
[`scripts/meta_analysis`](scripts/meta_analysis).

Two 10k trans-ancestry iterations were recovered and validated. The later file
contains 6,686,830 variants and reaches a maximum combined N of 7,789; the
earlier file contains 6,688,950 variants and reaches N = 7,664. Their
fingerprints and ancestry-specific maximum sample sizes are listed in
[`results/trans_ancestry/meta_analysis_manifest.tsv`](results/trans_ancestry/meta_analysis_manifest.tsv).

## Refined Chinese-only association test

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

The modal variant-level sample size was 5,592.

## Chinese-only summary-statistic freeze

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

## Downstream analyses

The separately retained April 2026 `10k_clean` LDSC analysis used an East Asian LD-score reference.
It reported observed-scale heritability of 0.3048 (SE 0.0949) and genetic
correlation with East Asian educational attainment of -0.0502 (SE 0.0835,
P=0.5476).

The successful MTAG v1.0.8 run completed on 2 December 2025 and combined the
10k HELIOS result with East Asian educational attainment. It retained
3,403,649 shared variants and is currently classified as exploratory.

The repository contains code, configuration, aggregate results and
non-identifying provenance. Participant-level HELIOS and SG100K data remain
controlled.
