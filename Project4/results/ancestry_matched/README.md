# SG100K ancestry-matched LDSC results

This directory contains aggregate results from the preliminary six-trait ancestry-matched validation. All 18 univariate and 45 bivariate matched analyses completed. Six additional univariate and 15 bivariate runs used a translated-row 1KG EAS comparator.

The Chinese comparator here is retained for provenance and is superseded for performance claims by the exact raw-genotype analysis in [`../../validation/controlled_reference_experiment/`](../../validation/controlled_reference_experiment/).

## Main findings

In the final controlled experiment, mean h2 was 0.1618 with the standard Singapore Chinese panel and 0.1601 with exact 1KG EAS, and no paired trait difference was nominally significant. Those values should be used in the chapter and public summary.

Chinese genetic-correlation estimates were effectively unchanged by panel choice: Pearson r = 0.9990, mean absolute difference = 0.0111, and direction agreement = 15/15. Six of 15 trait pairs passed Benjamini--Hochberg correction with both panels.

The Indian and Malay GWAS were much less informative. Mean chi-square was close to one, only height and BMI had h2/SE >= 1.96 in each ancestry, and no rg survived correction across 15 pairs. These results are consistent with limited GWAS power at N approximately 7,400 and 5,600; they should not be used to rank the local panels.

Type 2 diabetes h2 is on the observed scale. No liability-scale conversion was performed because population prevalence was not prespecified.

## Tables

- [`ancestry_matched_h2.tsv`](ancestry_matched_h2.tsv): all 24 univariate results, diagnostics, confidence intervals, and log paths.
- [`ancestry_matched_rg.tsv`](ancestry_matched_rg.tsv): all 60 bivariate results with BH correction within each ancestry/reference analysis.
- [`h2_summary_by_analysis.tsv`](h2_summary_by_analysis.tsv): compact h2 summaries.
- [`rg_summary_by_analysis.tsv`](rg_summary_by_analysis.tsv): compact rg summaries.
- [`chinese_h2_panel_comparison.tsv`](chinese_h2_panel_comparison.tsv): trait-level Chinese panel comparison.
- [`chinese_rg_panel_comparison.tsv`](chinese_rg_panel_comparison.tsv): Chinese rg concordance across panels.
- [`CHECKSUMS.sha256`](CHECKSUMS.sha256): checksums for this result directory.

Raw logs and run provenance are stored under [`../../validation/ancestry_matched_sg100k/`](../../validation/ancestry_matched_sg100k/).
Trait labels, transformations, and sample sizes are defined in [`../../metadata/historical_benchmark_traits.tsv`](../../metadata/historical_benchmark_traits.tsv).
