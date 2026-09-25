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

## Final GWAS phenotype preparation

`prepare_ancestry_g_for_regenie.R` represents the post-PCA preparation used by
the historical 10k and expanded 22k ancestry-specific analyses. For one
ancestry group it:

1. intersects the pooled cognitive-factor score with the ancestry analysis IDs;
2. masks values beyond five within-group standard deviations; and
3. applies the Blom inverse-normal transformation within that group.

```bash
Rscript prepare_ancestry_g_for_regenie.R \
  /controlled/path/pooled_scores.tsv \
  /controlled/path/Chinese_ids.tsv \
  /controlled/path/Chinese_g.pheno \
  g 5
```

Run the script separately for Chinese, Indian and Malay participants. Its
output is participant-level controlled data and must not be committed.
