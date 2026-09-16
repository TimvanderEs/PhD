# Project 3: HELIOS cognitive GWAS

Analysis code and non-sensitive reproducibility material for the HELIOS
cognitive GWAS project.

## Contents

- [`scripts/phenotype`](scripts/phenotype): reproduction of the six-task 10k
  cognitive factor.
- [`scripts/gwas`](scripts/gwas): REGENIE Step 1 and Step 2 commands recovered
  from the completed Chinese association run, plus a reconstruction of the
  input-preparation step.
- [`scripts/summary_statistics`](scripts/summary_statistics): conversion and
  build/identifier mapping of association results.
- [`scripts/ld_reference`](scripts/ld_reference): preparation of the SG100K
  covariates and covariate-adjusted LD scores.
- [`results/phenotype`](results/phenotype): aggregate PCA loadings,
  correlations and validation statistics.
- [`provenance`](provenance): software versions, source-file checksums and run
  status, including the source-archive audit.

The current audit of retained and missing 10k analysis records is in
[`WORKFLOW_INVENTORY.md`](WORKFLOW_INVENTORY.md).

## Data availability

Individual-level HELIOS and SG100K data cannot be deposited in this
repository. Access is subject to the relevant cohort consent, data-access and
institutional governance procedures. The public repository contains code,
small aggregate results and provenance only.

## Current scope

This archive is restricted to the 10k analysis. The cognitive-factor
calculation has been reproduced against the saved historical score for all
7,403 participants, and the completed Chinese REGENIE Step 1 and Step 2
commands have been recovered. The final Indian and Malay runs and the 10k
trans-ancestry meta-analysis will be added after their records are identified
or the analyses are rerun.
