# SG100K ancestry-matched validation

This validation uses the Chinese, Indian, and Malay association columns retained in the six SG100K release 2.3.1 provider packages. It is separate from the historical analysis in which the same trans-ancestry GWAS was paired with three different Singapore LD panels.

## Design

- Six traits: systolic blood pressure, waist-to-hip ratio, height, BMI, LDL cholesterol, and type 2 diabetes.
- Three matched analyses: Chinese GWAS with the Singapore Chinese panel, Indian GWAS with the Singapore Indian panel, and Malay GWAS with the Singapore Malay panel.
- Per matched analysis: six univariate h2 runs and all 15 pairwise rg runs.
- Chinese comparator: the same Chinese GWAS was also analysed with the 1KG EAS 1 cM panel.
- Benjamini--Hochberg correction is applied across the 15 rg tests within each GWAS-ancestry by LD-panel analysis.

The 1KG EAS LD-score files originally used rsIDs, whereas the Singapore inputs use `chr:position:alleles` identifiers. A complete 2,628,602-row crosswalk was used to translate the identifiers in a copy of the 1KG EAS panel. LD-score values were not changed. All 867,103 panel variants mapped uniquely; no variant was missing or duplicated.

## Status and interpretation limits

The Chinese/1KG comparison in this directory is an earlier validation that translated and row-subset existing 1KG LD-score files. It is retained for provenance. The primary Chinese reference comparison is now the exact raw-genotype experiment in [`../controlled_reference_experiment/`](../controlled_reference_experiment/), which uses 874,025 variants in identical order across all panels.

The final controlled experiment also retained full-precision delete-block values, allowing covariance-aware paired jackknife SEs for h2 differences.

No suitable public South Asian comparator was present on the EC2 instance. Indian and Malay outputs therefore document ancestry-matched local-panel estimates but do not test whether those panels outperform a public reference.

Type 2 diabetes h2 is reported on the observed scale. Population prevalence was not specified, so no liability-scale conversion was performed.

## Files

- [`input_sample_sizes.tsv`](input_sample_sizes.tsv): ancestry- and trait-specific GWAS sample sizes.
- `qc/`: preparation counts, run status, session information, and checksums.
- `h2/` and `rg/`: complete LDSC logs for the retained analyses. Paths in the aggregate result tables are relative to this validation directory.
- [`CHECKSUMS.sha256`](CHECKSUMS.sha256): checksums for the retained validation bundle.
- [`../../results/ancestry_matched/`](../../results/ancestry_matched/): reader-facing aggregate results and comparisons.

Provider-controlled GWAS files and individual-level genotype data are not included.

In `qc/input_counts.tsv`, `munged_rows` is the fixed 867,103-row merge-alleles output grid; variants without usable statistics are represented as missing. The `n_snps` field in the h2 result table gives the number retained in each regression.
