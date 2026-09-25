# Expanded 22k phenotype validation and GWAS

## Baseline

- Archived 12k extension PCA: **N = 13,774**, PC1 variance **38.93%**, exact saved-score reproduction.
- Combined 22k PCA: **N = 22,422**, PC1 variance **37.04%**, exact saved-score reproduction.
- All six combined-sample residuals were reproduced from the retained covariates (minimum task-wise r = **1.000000000000**).
- Omitting education from residualisation preserved the loading pattern (congruence **0.9999**) but produced a PC1 correlated **r = 0.981** with the current score.
- Loading congruence between the archived 12k and combined 22k solutions: **1.000**.
- The stored 22k score has higher values for poorer performance. Poster-facing
  loadings reverse the arbitrary PCA sign so higher values indicate better
  performance.
- Before GWAS, PC1 was split by genetic-ancestry analysis group, values beyond
  five within-group standard deviations were set to missing and the remaining
  scores were Blom inverse-normal transformed within group. REGENIE Step 2 did
  not apply a second RINT.

### Education covariate sensitivity

Omitting education from residualisation preserved the loading pattern (Tucker
congruence **0.9999**) and produced a score correlated **r = 0.981** with the
retained phenotype. This is evidence of phenotypic stability only; genetic
equivalence requires a matched GWAS sensitivity analysis.

## Genome-wide analysis

- Ancestry-specific REGENIE analyses included **17,966 Chinese**, **2,310 Indian** and **1,242 Malay** participants (total **N = 21,518**).
- The matching METAL output contained **10,286,285** variants before post-analysis filters and no variant reached P < 5 × 10⁻⁸ (minimum P = **8.17 × 10⁻⁸**).
- After retaining variants with `HetDf >= 1`, `HetPVal > 0.05` and
  `MaxFreq - MinFreq < 0.5`, **6,518,064** variants remained. None reached
  P < 5 × 10⁻⁸ (minimum P = **6.906 × 10⁻⁷**).
- The raw top marker, `chr2:90389726:T:G`, had `HetDf = 0` and direction
  `+??`, so it was represented in only one ancestry and failed the first
  retained filter.
- The archived poster record reported h² = **0.0392** (SE **0.0202**),
  intercept = **1.0008** (SE **0.0058**), mean χ² = **1.0175** and λGC =
  **1.0245**. Its log does not retain the LD-score path, so the exact reference
  for that run requires confirmation.
- A later checksum-linked rerun using the frozen SG100K 20 cM cov-LDSC
  reference reported h² = **0.0394** (SE **0.0188**), intercept = **1.0006**
  (SE **0.0061**), mean χ² = **1.0158** and λGC = **1.0177**.

## Measurement

- All six raw tasks are complete in the combined 22k PCA sample.
- Effective category counts were Quiz **6.1**, Working Memory **6.3**, and Pairing **22.3**.
- The combined construction retains Pairing zero-error participants; it does not repeat the reciprocal-Pairing exclusion found in the historical 10k pipeline.

## PCA stability

- Reaction Time loading 95% interval: **[0.428, 0.450]**.
- Stroop Box loading 95% interval: **[0.545, 0.557]**.
- Stroop Ink loading 95% interval: **[0.479, 0.494]**.
- Leave-one-task-out score correlations: **r = 0.947–0.986**.
- Recorded ancestry loading congruence: **0.999–1.000**.
- HELIOS10K/HELIOS20K loading congruence: **0.998–1.000**.
- Collapsed-speed versus six-task PC1: **r = 0.931**.

The speed-heavy loading pattern persists in the larger sample, across the two recorded freezes, and across ancestry groups. It is therefore unlikely to be solely a small-pilot-sample artefact. Collapsing the three latency measures still changes the score materially, supporting the interpretation that giving speed three indicators affects the composition of PC1.

## Genetic signal across freezes

The 12k extension had no detectable cov-LDSC signal (h² = **-0.0036**, SE
**0.0276**). The combined 22k Z scores reconstructed from the retained 10k and
12k components with r = **0.937**; reversing the 12k effects gave r =
**-0.321**. The low 22k estimate is therefore not explained by a global sign
error or the retained post-meta-analysis filters.

The approximately 30% result belongs to the later refined Chinese-only 10k
analysis, not to the 10k trans-ancestry iteration. See
[`03_genetic_signal/genetic_signal_audit.md`](03_genetic_signal/genetic_signal_audit.md).

The 10k and 22k estimates arise from distinct phenotype and analysis
iterations. They should not be interpreted as a clean sample-size scaling
experiment without matching the phenotype definition, SNP set, ancestry
composition and LD-score reference.

At phenotype level, pooled and within-ancestry Blom scores preserved identical
within-ancestry rank order (Spearman rho = 1) and correlated at r =
0.9990–1.0000 in the reconstructed expanded sample. This is not a matched GWAS
comparison; see [`../phenotype_pipeline`](../phenotype_pipeline).
