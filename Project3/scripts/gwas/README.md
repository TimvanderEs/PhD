# Ancestry-specific GWAS

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
and REGENIE predictions are not distributed publicly. The retained workflow is
the final Chinese-only 10k freeze.
