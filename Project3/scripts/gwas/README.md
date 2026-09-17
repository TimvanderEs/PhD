# Ancestry-specific GWAS

The Project 3 discovery analysis ran REGENIE separately in Chinese, Indian and
Malay participants before meta-analysis. The recovered 10k Step 2 records use a
precomputed `g` phenotype, ancestry-specific PC covariate files, block size 200
and minimum MAC 0.5. `run_regenie_step2_10k_ancestry.sh` provides that command
for any of the three groups.

The genotype inputs recorded in the logs contain 7,043,723 Chinese, 8,387,428
Indian and 7,660,514 Malay variants. Full genotype, phenotype, covariate and
LOCO prediction files remain controlled.

## Refined Chinese-only freeze

`run_regenie_step1_10k_chinese.sh` and
`run_regenie_step2_10k_chinese.sh` reproduce the completed REGENIE v4.1
commands for the Chinese 10k cognitive-factor analysis. Both runs completed on
28 November 2025. Step 1 used 279,546 LD-pruned variants and 5,700 participants;
Step 2 tested 7,043,723 variants in the genotype input for the same 5,700
participants.

The input-preparation script provides the required REGENIE schema:

- `prepare_regenie_10k_inputs.R` creates the `g_raw` phenotype and covariate
  files with age, sex, age interactions and 20 ancestry-specific PCs.

The completed logs use `CHI_gwas_10k_SESfix.cov` for Step 1 and
`CHI_gwas_10k_nodup.cov` for Step 2. Both runs request the same covariate columns
and report 6,353 participants with covariate data.

All scripts take controlled-data paths as arguments. Participant-level inputs
and REGENIE predictions are not distributed publicly.
