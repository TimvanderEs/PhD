# HELIOS cognition validation

This directory contains sensitivity checks for the six-task cognitive phenotype. The historical phenotypes and GWAS outputs remain unchanged.

The original 10k validation covers:

- exact reproduction of the historical PCA;
- task-level measurement summaries;
- Pearson and Spearman correlations;
- 1,000 participant bootstrap replicates;
- leave-one-task-out PCA;
- a four-indicator PCA in which the three speed tasks are collapsed;
- ancestry-specific PCA.

Participant-level phenotype and genotype files are not included. Run the analysis with authorised local copies:

```bash
Project3/validation/run_phase1.sh \
  --cognitive /path/to/Cognition_phenotypes_10K_Pritesh.txt \
  --covariates /path/to/merged_pheno_cov_filtered_10k.txt \
  --historical-g /path/to/g_scores_10k.txt \
  --validation-plan /path/to/HELIOS_cognition_validation_plan.md \
  --historical-r-script /path/to/Final_Pheno_script_fixed.R \
  --step1-log /path/to/CHI_10k_step1.log \
  --step2-log /path/to/HELIOS_10k_Chinese_g.log
```

The aggregate results and figures are version-controlled. Participant scores are neither written nor committed. See [phase1_summary.md](phase1_summary.md) for the findings.

The expanded 22k phenotype validation and GWAS presented at WCPG 2026 are in
[`22k`](22k). This record reproduces the archived 12k extension score and the
combined 22k residualisation/PCA, repeats the measurement and PCA-structure
checks, and compares the genetic signal across the retained 10k, 12k and 22k
GWAS files. See [`22k/22k_summary.md`](22k/22k_summary.md) for the findings.

The WCPG-focused audit of the near-zero HELIOS 10k cognition–KoGES education
genetic correlation is in
[`edu_local_cancellation`](edu_local_cancellation). The historical control and
corrected quantitative-education LAVA runs are complete. The corrected result
found 47 positive and 45 negative covariance estimates across 92 eligible
regions, with 90.8% descriptive cancellation; no region survived FDR
correction. The report explains why this is consistent with, but does not
prove, genome-wide cancellation.
