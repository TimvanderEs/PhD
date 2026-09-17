# Project 3: HELIOS cognitive GWAS

Code and aggregate reproducibility material for the HELIOS cognitive GWAS.

## Contents

- [`scripts/phenotype`](scripts/phenotype): construction of the six-task
  cognitive factor.
- [`scripts/qc`](scripts/qc): ancestry-specific genotype quality control.
- [`scripts/gwas`](scripts/gwas): REGENIE input preparation and the completed
  Chinese Step 1 and Step 2 commands.
- [`scripts/summary_statistics`](scripts/summary_statistics): conversion and
  build/identifier mapping of association results.
- [`scripts/ld_reference`](scripts/ld_reference): preparation of the SG100K
  covariates and covariate-adjusted LD scores.
- [`results`](results): aggregate phenotype results and the summary-statistic
  freeze manifest.
- [`provenance`](provenance): run records, software versions and source-file
  checksums.

The retained analysis is described in [`WORKFLOW.md`](WORKFLOW.md).

## Analysis scope

This release documents the final Chinese 10k GWAS freeze used in the downstream
Project 2 analyses. Earlier Chinese–Indian–Malay meta-analyses used previous
versions of the cognitive phenotype and are not part of the retained workflow.

## Data availability

Individual-level HELIOS and SG100K data cannot be deposited in this
repository. Access is subject to the relevant cohort consent, data-access and
institutional governance procedures. The public repository contains code,
small aggregate results and provenance only. Summary-statistic access and
release are handled separately under the relevant study approvals.
