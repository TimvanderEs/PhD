# HELIOS 10k cognition validation

This directory contains sensitivity checks for the six-task cognitive phenotype. The historical phenotype and GWAS remain unchanged.

Phase 1 covers:

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
