# Phenotype-to-GWAS lineage

This record separates the cognitive-factor construction from the final
phenotype preparation used for ancestry-specific GWAS.

## Retained expanded 22k analysis

The WCPG 2026 analysis used the following sequence:

1. residualise the six cognitive tasks in the pooled sample for age, age
   squared, sex, age-by-sex and age-squared-by-sex terms, Indian and Malay
   indicators, and education as a categorical covariate;
2. centre and scale the residuals and extract PC1 in 22,422 participants;
3. split PC1 into the Chinese, Indian and Malay genetic-analysis groups;
4. set values beyond five within-group standard deviations to missing;
5. apply a Blom inverse-normal transformation within group; and
6. run REGENIE Step 2 with 20 ancestry-specific PCs, block size 200 and minimum
   MAC 0.5, without applying a second RINT.

The archived PC1 is oriented so that higher values indicate poorer cognitive
performance. Figures reverse this arbitrary sign so that higher values indicate
better performance.

The transformation lineage is supported by the recovered phenotype-preparation
script, the phenotype filenames in the completed REGENIE logs, and exact sample
count reconstruction for the historical 10k and nested 12k runs. The exact
participant-level September 2025 phenotype files are controlled data and were
not retained in the public repository.

## Analysis versions

[`pipeline_freezes.tsv`](pipeline_freezes.tsv) distinguishes the primary 22k
analysis from the historical 10k, the nested 12k diagnostic analysis and the
later refined Chinese-only freeze. These are different analysis versions and
should not be combined into a sample-size progression.

## RINT sensitivity

[`rint_scope_comparison.tsv`](rint_scope_comparison.tsv) compares pooled and
within-ancestry Blom transformations of the same saved factor. Within each
ancestry the transformations retained identical rank order (Spearman rho = 1).
Pearson correlations were 0.9989–0.9999 in the exact historical 10k sample and
0.9990–1.0000 in the reconstructed expanded sample.

This establishes robustness of phenotype ordering only. No completed GWAS pair
was retained in which RINT scope was the sole difference, so the repository
does not claim GWAS-level equivalence.

## Code

The controlled-input implementation of the post-PCA step is
[`prepare_ancestry_g_for_regenie.R`](../../scripts/phenotype/prepare_ancestry_g_for_regenie.R).
Its participant-level output must remain outside the repository.

