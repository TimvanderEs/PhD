# HELIOS 12k extension and combined 22k validation

This directory contains aggregate checks of the expanded HELIOS cognitive phenotype. Participant-level inputs are not included.

The validation covers:

- exact reproduction of the archived 12k extension score;
- exact reconstruction of all six combined-sample residuals;
- exact reproduction and direction check of the combined 22k PC1;
- a residualisation sensitivity analysis omitting education;
- raw task distributions, including ancestry and freeze summaries;
- 1,000 participant bootstrap replicates;
- leave-one-task-out and collapsed-speed PCAs;
- ancestry- and freeze-specific correlation and loading checks;
- sample sizes and settings recovered from the three REGENIE logs;
- baseline LDSC quality-control metrics.

The combined PCA includes 22,422 participants. Of these, 21,565 have recorded ancestry and freeze labels, and 21,518 entered the ancestry-specific GWASs after genotype/covariate intersection. See [`sample_flow.tsv`](00_baseline/sample_flow.tsv) and [`22k_summary.md`](22k_summary.md).

## Reproduction

Python 3 with NumPy, pandas, and Pillow is required. Run with authorised local copies of the participant-level files:

```bash
Project3/validation/run_22k_validation.sh \
  --cognitive /path/to/merged.txt \
  --covariates /path/to/merged_covariates_fixed.txt \
  --metadata /path/to/Full_data_HELIOS_22k.txt \
  --pca-input /path/to/22k_pre_norm_post_resid_Edu.txt \
  --extension-residuals /path/to/residualized_cognitive_phenotypes_12k.txt \
  --extension-g /path/to/g_scores_12k.txt \
  --validation-plan /path/to/HELIOS_cognition_validation_plan_for_codex.md \
  --gwas-log /path/to/HELIOS_22k_COG_Chinese.log \
  --gwas-log /path/to/HELIOS_22k_COG_Indian.log \
  --gwas-log /path/to/HELIOS_22k_COG_Malay.log \
  --ldsc-log /path/to/current_22k_ldsc.log \
  --output-root /path/to/output/22k
```

Only aggregate tables, reports, and figures are written.
