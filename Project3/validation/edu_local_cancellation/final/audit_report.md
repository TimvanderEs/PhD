# Final audit report

## Question

Could the small genome-wide genetic correlation between HELIOS cognition and
KoGES educational attainment reflect cancellation between positive and
negative local genetic covariance?

## Validation

The historical four-trait LAVA analysis was reproduced before the corrected
result was used. The control matched the archived row sets, all univariate
quantities and the bivariate rho and r2 estimates exactly. The archived and
control confidence intervals and P values differed only by the expected
unseeded Monte Carlo variation in LAVA 0.1.5; P-value correlation was 0.999994
and no nominal or FDR decision changed.

Changing only KoGES education from binary to quantitative retained all 92
HELIOS–education regions. Control and corrected rho estimates correlated
0.9997, with no sign changes. The correction changed one nominal decision and
no FDR decision.

## Findings

The genome-wide LDSC estimate was rg = 0.0502 (SE = 0.0835, P = 0.5476) after
orienting HELIOS so that positive values indicate better cognitive
performance. Of 3,064 LAVA blocks, 92 were eligible for bivariate analysis.
Their local covariance estimates were nearly evenly divided: 47 positive and
45 negative. The positive contribution summed to 0.03560, the absolute
negative contribution to 0.04282, and the signed total to -0.00722. Thus only
9.2% of the absolute local covariance remained after summation, corresponding
to a descriptive cancellation fraction of 90.8%.

The pattern was stable to leave-one-chromosome-out analysis (cancellation
86.9–97.6%), MHC exclusion (91.7%) and removal of the largest absolute local
contribution (93.8%). Fifteen regions were nominally significant, but none
survived BH correction across the 92 HELIOS–education tests.

## Assessment

The data support the presence of opposing local covariance within the eligible
subset. They provide suggestive, not definitive, evidence that local
cancellation contributes to the small genome-wide correlation. The main
constraints are that only 3.0% of all LAVA blocks were eligible, the subset
contained 11.1% of summed estimable HELIOS local h2 and 6.48% of the education
total, and no local correlation survived FDR correction. The aggregate rg of
the eligible subset (-0.0171) is therefore descriptive and cannot replace the
genome-wide LDSC estimate.

## Suggested manuscript wording

> The genome-wide genetic correlation between HELIOS cognitive performance
> and KoGES educational attainment was small and non-significant (rg = 0.050,
> SE = 0.084, P = 0.548). Among 92 regions eligible for local analysis, genetic
> covariance was bidirectional (47 positive and 45 negative estimates), and
> positive and negative contributions largely offset (descriptive cancellation
> fraction = 90.8%). This pattern was stable in chromosome, MHC and outlier
> sensitivity analyses, although no local association survived FDR correction
> and most genomic regions were not eligible. The findings are therefore
> consistent with local cancellation but do not establish it as the cause of
> the near-zero genome-wide estimate.
