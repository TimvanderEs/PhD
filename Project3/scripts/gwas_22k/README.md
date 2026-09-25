# Expanded 22k ancestry-specific GWAS

The WCPG 2026 analysis ran REGENIE 4.1 Step 2 separately in Chinese, Indian and
Malay participants. The three archived logs record the same options represented
by `run_regenie_step2_22k_ancestry.sh`:

- one ancestry-specific PLINK bed prefix;
- one externally transformed `g` phenotype column;
- 20 ancestry-specific principal components;
- block size 200; and
- minimum MAC 0.5.

The completed analyses included 17,966 Chinese, 2,310 Indian and 1,242 Malay
participants. Input variant counts and log checksums are in
[`baseline_gwas_qc.tsv`](../../validation/22k/00_baseline/baseline_gwas_qc.tsv).

The input phenotype was prepared after PCA by splitting the pooled score into
the three ancestry analysis groups, masking values beyond five within-group
standard deviations and applying a Blom inverse-normal transformation within
group. The Step 2 command therefore does not use `--apply-rint`. The portable
preparation script is
[`prepare_ancestry_g_for_regenie.R`](../phenotype/prepare_ancestry_g_for_regenie.R).

Example:

```bash
bash run_regenie_step2_22k_ancestry.sh \
  /path/to/regenie \
  /path/to/ancestry_genotypes \
  /path/to/ancestry_phenotype.txt \
  /path/to/ancestry_covariates.txt \
  /path/to/ancestry_predictions.list \
  /path/to/output_prefix
```

The original prediction-list files were referred to by generic names and were
not retained with immutable Step 1 lineage. This wrapper therefore reproduces
the verified Step 2 settings only. `--phenoCol g` and the 20-column
`--covarColList` are stated explicitly here; the archived commands omitted
them because their input files contained only those analysis columns.
Controlled genotype, phenotype, covariate and prediction files are not
distributed.
