# Cognitive phenotype

`reproduce_10k_cognitive_factor.R` is the cleaned implementation of the
historical six-task PCA used for the HELIOS 10k analysis. It replaces a working
script that contained repeated and superseded blocks.

The workflow:

1. applies reciprocal transformations to reaction time, Stroop box time,
   Stroop ink time and pairing guesses;
2. masks values more than four standard deviations from the task mean;
3. residualises each task for age, sex, age-by-sex, age squared,
   age-squared-by-sex, ancestry and household-income terms; and
4. derives the first principal component from the standardised residuals.

The historical saved score is oriented so that higher `g` values indicate
poorer cognitive performance. The sign-equivalent higher-performance loadings
are also reported for figures and interpretation.

Run with:

```bash
Rscript reproduce_10k_cognitive_factor.R \
  Cognition_phenotypes_10K_Pritesh.txt \
  merged_pheno_cov_filtered_10k.txt \
  private_output \
  g_scores_10k.txt
```

The fourth argument is optional and is used only to compare the reproduced
score with the historical score. The three aggregate output tables may be
shared; `cognitive_factor_scores.tsv` is participant-level data and must remain
outside the repository.
