# HELIOS 10k workflow inventory

This inventory records which parts of the 10k cognitive GWAS are supported by
an original script, log or result. Draft methods are used to identify the
intended workflow but are not treated as evidence that an analysis ran.

## Analysis status

| Stage | Evidence | Status |
|---|---|---|
| Cognitive factor | Historical code, phenotype extracts and saved score | Reproduced for the same 7,403 participants |
| Genotype QC | Historical QC scripts and matching genotype prefix | Clean reconstruction available; completed logs not recovered |
| Chinese REGENIE Step 1 | Completed REGENIE v4.1 log | Curated script available |
| Chinese REGENIE Step 2 | Completed REGENIE v4.1 log | Curated script available |
| Indian and Malay association tests | Multiple completed logs with differing phenotype and sample sets | Final runs not established |
| Trans-ancestry meta-analysis | No verified final 10k METAL output | Not recovered |
| Summary-statistic conversion | Executed shell and Python scripts | Curated scripts available |
| SG100K covariate-adjusted LD scores | Source scripts and completed chromosome logs | Curated scripts available |
| LDSC | Completed April 2026 10k analysis | Recovered |
| MTAG with EAS educational attainment | Completed MTAG v1.0.8 log | Exploratory |

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

The Step 1 and Step 2 logs name their covariate files
`CHI_gwas_10k_SESfix.cov` and `CHI_gwas_10k_nodup.cov`, respectively. Both
contain the same 25 requested covariate columns and report covariate data for
6,353 participants. The retained files are not available for a byte-level
comparison.

The modal variant-level sample size was 5,592. The corresponding post-QC LDSC
input had mean chi-square 1.063.

## Downstream analyses

The April 2026 `10k_clean` LDSC analysis used an East Asian LD-score reference.
It reported observed-scale heritability of 0.3048 (SE 0.0949) and genetic
correlation with East Asian educational attainment of -0.0502 (SE 0.0835,
P=0.5476).

The successful MTAG v1.0.8 run completed on 2 December 2025 and combined the
10k HELIOS result with East Asian educational attainment. It retained
3,403,649 shared variants and is currently classified as exploratory.

## Input integrity

The EC2 file named
`HELIOS_10k_COG_Chinese_QCrem_miss02_maf005_hwe1e06_g.regenie` contains the
expected 7,043,723 association rows followed by a malformed row and appended
chromosome 22 records from a higher-N run. It must not be released or analysed
as a single result. Clean downstream exports and the original execution log
are the valid records for the 10k analysis.

## Work needed before release

1. Verify the reconstructed Chinese genotype QC script against retained PLINK
   files or logs.
2. Confirm whether Indian and Malay 10k association tests were completed; rerun
   them if necessary.
3. Produce the retained trans-ancestry meta-analysis from the verified
   ancestry-specific results.
4. Confirm whether the 10k LDSC and exploratory MTAG analyses are retained in
   the manuscript.

Only code, configuration, aggregate results and non-identifying provenance
will be committed. Participant-level HELIOS and SG100K data remain controlled.
