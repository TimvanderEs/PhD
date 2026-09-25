# Expanded 22k phenotype validation and GWAS

This directory documents the expanded HELIOS analysis presented at WCPG 2026.
It includes aggregate phenotype checks, ancestry-specific GWAS QC, the METAL
record and cross-freeze sensitivity analyses. Participant-level inputs are not
included.

See the [WCPG 2026 conference page](../../WCPG_2026.md) for a short guide.

The validation covers:

- exact reproduction of the archived 12k extension score;
- exact reconstruction of all six combined-sample residuals;
- exact reproduction and direction check of the combined 22k PC1;
- reconstruction of the ancestry split, 5-SD mask and external Blom RINT used
  to prepare the GWAS phenotype;
- a residualisation sensitivity analysis omitting education;
- raw task distributions, including ancestry and freeze summaries;
- 1,000 participant bootstrap replicates;
- leave-one-task-out and collapsed-speed PCAs;
- ancestry- and freeze-specific correlation and loading checks;
- sample sizes and settings recovered from the three REGENIE logs;
- baseline LDSC quality-control metrics; and
- a checksum-linked comparison of the 10k, 12k and 22k genetic signal.

The combined PCA includes 22,422 participants. Of these, 21,565 have recorded ancestry and freeze labels, and 21,518 entered the ancestry-specific GWASs after genotype/covariate intersection. See [`sample_flow.tsv`](00_baseline/sample_flow.tsv) and [`22k_summary.md`](22k_summary.md).

## Phenotype orientation

PCA sign is arbitrary. The archived 22k PC1 is stored in the orientation in
which higher values indicate poorer cognitive performance. Figures and public
interpretation reverse the sign so that higher values indicate better
performance. This does not alter P values, heritability or model fit;
association-effect signs reverse accordingly.

## GWAS phenotype

The saved pooled PC1 was not passed directly to REGENIE. It was split into the
Chinese, Indian and Malay genetic-analysis groups, values beyond five
within-group standard deviations were set to missing, and the remaining values
were Blom inverse-normal transformed within group. The completed Step 2 runs
therefore did not use `--apply-rint`. Aggregate pooled-versus-within-ancestry
checks are in [`../phenotype_pipeline`](../phenotype_pipeline).

## Education covariate sensitivity

Omitting education from residualisation preserved the PCA loading pattern
(Tucker congruence 0.9999) and produced a participant score correlated r =
0.981 with the retained phenotype. These results show phenotypic robustness,
but do not establish genetic equivalence; a matched GWAS sensitivity analysis
is still required.

## LDSC record

The archived poster estimate was h² = 0.0392 (SE 0.0202). The retained log does
not identify its LD-score path, so the exact reference for that run requires
confirmation. A later checksum-linked rerun used the frozen SG100K 20 cM
covariate-adjusted LD-score reference and gave h² = 0.0394 (SE 0.0188).

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
  --validation-plan /path/to/HELIOS_cognition_validation_plan.md \
  --gwas-log /path/to/HELIOS_22k_COG_Chinese.log \
  --gwas-log /path/to/HELIOS_22k_COG_Indian.log \
  --gwas-log /path/to/HELIOS_22k_COG_Malay.log \
  --ldsc-log /path/to/current_22k_ldsc.log \
  --output-root /path/to/output/22k
```

Only aggregate tables, reports, and figures are written.

The cross-freeze LDSC and Z-score checks are in
[`03_genetic_signal`](03_genetic_signal). They use the same SG100K cov-LDSC
reference for each freeze and keep the trans-ancestry and refined Chinese-only
10k analyses separate.

The verified expanded-cohort REGENIE Step 2 wrapper is in
[`scripts/gwas_22k`](../../scripts/gwas_22k). Exact Step 1 prediction-file
lineage was not retained; the repository therefore does not claim that the
Step 1 transformation state was independently verified.
