# Ancestry-specific GWAS

## Expanded 22k analysis

The analysis presented at WCPG 2026 used ancestry-specific REGENIE 4.1 Step 2
runs in 17,966 Chinese, 2,310 Indian and 1,242 Malay participants. The exact
Step 2 options recovered from the three archived logs are represented in
[`scripts/gwas_22k`](../gwas_22k). The input `g` phenotype had already been
Blom inverse-normal transformed within ancestry, so Step 2 did not use
`--apply-rint`. The original prediction lists were not retained with immutable
Step 1 lineage.

## Historical 10k analysis

The Project 3 discovery analysis ran REGENIE separately in Chinese, Indian and
Malay participants before meta-analysis. The recovered 10k Step 2 records use a
precomputed `g` phenotype, ancestry-specific PC covariate files, block size 200
and minimum MAC 0.5. The score was prepared by a within-ancestry 5-SD mask and
external Blom RINT. `run_regenie_step2_10k_ancestry.sh` provides the association
command for any of the three groups.

The genotype inputs recorded in the logs contain 7,043,723 Chinese, 8,387,428
Indian and 7,660,514 Malay variants. Full genotype, phenotype, covariate and
LOCO prediction files remain controlled.

## Later refined Chinese-only freeze

`run_regenie_step1_10k_chinese.sh` and
`run_regenie_step2_10k_chinese.sh` reproduce the completed REGENIE v4.1
commands for the Chinese 10k cognitive-factor analysis. Both runs completed on
28 November 2025. Step 1 used 279,546 LD-pruned variants and 5,700 participants;
Step 2 tested 7,043,723 variants in the genotype input for the same 5,700
participants.

The input-preparation script provides the required REGENIE schema:

- `prepare_regenie_10k_inputs.R` accepts only an explicitly named `g_raw`
  column and creates the phenotype and covariate schemas with age, sex, age
  interactions and 20 ancestry-specific PCs. It does not infer the upstream
  score definition.

The completed logs use `CHI_gwas_10k_SESfix.cov` for Step 1 and
`CHI_gwas_10k_nodup.cov` for Step 2. Both runs request the same covariate columns
and report 6,353 participants with covariate data.

The historical Step 1 command did not use `--apply-rint`, while Step 2 did.
These wrappers deliberately preserve that record. They should not be described
as a matched RINT implementation; a new analysis would apply the option in
both steps.

All scripts take controlled-data paths as arguments. Participant-level inputs
and REGENIE predictions are not distributed publicly.
