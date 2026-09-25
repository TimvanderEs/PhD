# 12k and combined 22k baseline reproduction

The archived 12k extension score was reproduced in **13,774 participants** (r = **1.000000000000**; PC1 variance explained = **38.93%**).

The current combined 22k score was reproduced in **22,422 participants** (r = **1.000000000000**; PC1 variance explained = **37.04%**). The maximum score difference after sign and scale alignment was **1.07e-13**.

All six task residuals were independently reconstructed from the retained analysis table. The model used age, sex as a factor, age squared, age-by-sex and age-squared-by-sex terms, Indian and Malay ancestry indicators, and education as a factor. The minimum task-wise residual correlation was **1.000000000000** and the largest absolute numerical difference was **1.39e-10**. The **57** records without a complete residualisation covariate set account for the difference between the 22,479-row analysis table and the 22,422-person PCA.

Removing education from the residualisation model left the loading vector nearly unchanged (Tucker congruence **0.9999**) but changed participant scores enough that the two PC1s correlated **r = 0.981**. This is a phenotype-construction sensitivity only; its genetic effect requires a GWAS rerun with the no-education score.

The stored combined score is oriented so that higher values indicate poorer performance. It is multiplied by -1 in all higher-performance sensitivity outputs. PCA sign was fixed from the task directions before any genetic result was considered.

For the ancestry-specific GWAS, the pooled score was split by genetic-analysis group, values beyond five within-group standard deviations were set to missing and the remaining values were Blom inverse-normal transformed within group. The archived Step 2 commands therefore used a pre-transformed `g` column and did not request a second RINT.

The matching ancestry-specific REGENIE analyses included 21,518 participants. Their METAL output contained 10,286,285 variants before filtering and 6,518,064 variants after retaining results represented in at least two ancestry groups, with heterogeneity P > 0.05 and maximum allele-frequency difference < 0.5. No variant reached P < 5 × 10⁻⁸ before or after filtering. The corresponding minima were 8.17 × 10⁻⁸ and 6.91 × 10⁻⁷, respectively. Aggregate checks are recorded in `baseline_meta_analysis_qc.tsv`.

Baseline status: **PASS**.
