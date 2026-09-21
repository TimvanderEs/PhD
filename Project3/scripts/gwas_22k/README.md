# Expanded 22k ancestry-specific GWAS

The WCPG 2026 analysis ran REGENIE 4.1 Step 2 separately in Chinese, Indian and
Malay participants. The three archived logs record the same options represented
by `run_regenie_step2_22k_ancestry.sh`:

- one ancestry-specific PLINK bed prefix;
- one `g` phenotype column;
- 20 ancestry-specific principal components;
- block size 200; and
- minimum MAC 0.5.

The completed analyses included 17,966 Chinese, 2,310 Indian and 1,242 Malay
participants. Input variant counts and log checksums are in
[`baseline_gwas_qc.tsv`](../../validation/22k/00_baseline/baseline_gwas_qc.tsv).

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
the verified Step 2 command only. Controlled genotype, phenotype, covariate and
prediction files are not distributed.
