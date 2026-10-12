# Final controlled-analysis report

## Question

Does a large Singapore Chinese LD reference, or a pooled Singapore covariate-adjusted reference, improve LDSC results relative to 1000 Genomes EAS when all other analytical choices are held constant?

## Result

No measurable advantage was demonstrated for the six Chinese SG100K GWAS examined here.

- Standard Chinese Singapore versus exact 1000 Genomes EAS: mean h2 0.1618 versus 0.1601; mean SE 0.02035 versus 0.02022.
- Trait-level paired h2 differences: 0.0002 to 0.0037; no nominal P < 0.05 and no BH-adjusted P < 0.05.
- Genetic correlation: r = 0.99945 across 15 pairs, mean absolute difference 0.00736, 15/15 directions concordant.
- Reference size: h2 and SE were nearly unchanged from n = 500 to n = 27,178, even though LD-score agreement with the full panel improved as expected.
- Held-out simulation: mean h2 0.2489 for 1000 Genomes EAS and 0.2464 for both local panels. Five replicates are informative as a pipeline check but not powered to rank panels.

The larger local sample improves the numerical stability of estimated LD scores, but sampling variation in the LD reference was not the limiting source of uncertainty for these GWAS. The 1000 Genomes EAS panel appears sufficiently well matched for the variants and Chinese GWAS considered.

The pooled analysis produced a separate methodological finding: adjusting for 20 PCs changed its h2 scale materially, whereas the same adjustment had little effect within the Chinese-only sample. Pooled multi-ancestry LD scores therefore require an explicitly justified covariate model and should not be substituted for an ancestry-specific panel without sensitivity analyses.

## Statistical note

An initial summary incorrectly used bias-corrected pseudo-values built from four-decimal LDSC log estimates. Multiplication by 200 blocks amplified rounding error and produced implausible point differences. The corrected analysis retains the full-sample h2 contrast as the point estimate and calculates the paired jackknife SE from the full-precision delete-block differences. The corrected script and regenerated tables are included in this repository.

## Supported chapter conclusion

Large local LD references are scientifically valuable, especially for populations lacking an appropriate public panel, but reference size alone does not guarantee better LDSC performance. Under a controlled 1 cM specification, the Singapore Chinese and 1000 Genomes EAS panels gave effectively equivalent h2 and rg results. The chapter should focus on reference construction, population representation, and robustness rather than claim a universal increase in heritability, precision, or power.
