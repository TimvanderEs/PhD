# Genetic signal across the HELIOS freezes

This check compares the retained 10k trans-ancestry meta-analysis, the 12k
extension, the combined 22k meta-analysis and the later refined Chinese-only
10k GWAS. The files were aligned to the same frozen SG100K 20 cM cov-LDSC
reference. Input fingerprints and retained variant counts are in
[`input_manifest.tsv`](input_manifest.tsv); reference fingerprints are in
[`reference_manifest.tsv`](reference_manifest.tsv).

The preparation, cov-LDSC and concordance steps are implemented in
[`prepare_ldsc_iteration_audit.py`](../../scripts/prepare_ldsc_iteration_audit.py),
[`run_ldsc_iteration_audit.sh`](../../scripts/run_ldsc_iteration_audit.sh) and
[`summarize_ldsc_iteration_concordance.py`](../../scripts/summarize_ldsc_iteration_concordance.py).

## Results

The combined 22k estimate was h² = 0.0394 (SE 0.0188) with a free intercept of
1.0006. The 12k extension had no detectable signal (h² = -0.0036, SE 0.0276).
These results were essentially unchanged by the retained meta-analysis filters
or by using the stricter allele-frequency-difference threshold of 0.05.
The retained workflow uses a threshold of 0.50; the historical 10k filename
ending in `AFdiff05` is therefore not used as a definition of the filter.

The 10k trans-ancestry free-intercept estimate was h² = 0.2023 (SE 0.0530), but
its intercept was 0.9751 while its mean chi-square was only 1.0078. Constraining
the intercept to one reduced the estimate to 0.0541 (SE 0.0407), close to the
22k constrained estimate of 0.0432 (SE 0.0150). The apparent 10k-to-22k drop is
therefore not robust to the intercept specification.

The approximately 30% estimate belongs to a different analysis: the later
refined Chinese-only 10k GWAS. Its matched free-intercept result was h² =
0.2672 (SE 0.0911). It should not be presented as the baseline of a
trans-ancestry 10k-to-22k series.

The 22k Z scores closely matched an N-weighted reconstruction from the 10k and
12k components (r = 0.937). Reversing the 12k effects gave r = -0.321. This
rules out a global sign reversal in the retained combined result. The same
pattern held under both allele-frequency-difference thresholds.

Phenotypic checks also argue against collapse of the cognitive factor. PC1
explained 38.93% of variance in the 12k extension and 37.04% in the combined
sample, compared with 34.49% in the historical 10k construction. The 12k and
22k loading vectors had Tucker congruence 0.9996, and the saved 12k score and
combined-PCA score correlated r = 0.979 in the same 13,774 participants.

## Interpretation

The expanded phenotype is phenotypically coherent, but the 12k extension adds
little detectable common-SNP signal. Stable PCA loadings do not imply that the
between-person variance is genetic: stable environmental, cohort and testing
effects can preserve the factor structure while lowering SNP heritability.
With mean chi-square values near one, LDSC estimates are also sensitive to the
intercept model.

The 12k extension GWAS was predominantly Chinese (12,363 of 13,751
participants), so the result is unlikely to be explained solely by adding the
smaller Indian and Malay groups. Sample counts are in
[`gwas_sample_sizes.tsv`](gwas_sample_sizes.tsv).

The refined Chinese-only 10k and combined 22k results remain genetically
compatible: on their shared variant set, cov-LDSC gave r_g = 0.908 (SE 0.281,
P = 0.0012). This is not an independent replication because the participant
sets overlap. It supports compatible effect directions, but the 22k
heritability is low and the uncertainty in r_g is substantial.

The exact Step 1 prediction files used by the September 2025 12k and 22k Step
2 runs were not retained with immutable names. Later logs using the same generic
prediction-list names have different dates and participant counts. This does
not explain the results above, but it prevents an exact Step 1 audit from the
recovered files alone.

## Reporting decision

Use h² = 0.039 (SE 0.019) as the combined 22k cov-LDSC estimate. Describe the
12k extension as having no detectable LDSC signal, not as having zero true
heritability. Keep the refined Chinese-only 10k estimate separate and avoid a
claim that heritability fell from approximately 30% to 4% after sample
expansion.

Full estimates are in [`ldsc_iteration_comparison.tsv`](ldsc_iteration_comparison.tsv)
and the sign checks are in [`zscore_concordance.tsv`](zscore_concordance.tsv).
The audit used cov-LDSC v1.0.0 at software commit
`926bbbcdde8bac30a8b53c1e9ebf283b5934bda8`.
