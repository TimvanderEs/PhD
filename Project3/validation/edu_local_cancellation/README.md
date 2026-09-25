# HELIOS cognition–KoGES education local-cancellation audit

This audit tests whether opposing local genetic effects could contribute to the
small genome-wide genetic correlation between HELIOS cognition and KoGES
educational attainment. Results are presented with higher HELIOS values
meaning better cognitive performance.

## Result

The audit was completed on 19 September 2026.

- Global LDSC: rg = 0.0502, SE = 0.0835, P = 0.5476 after orienting HELIOS to
  higher cognitive performance.
- LAVA tested 3,064 East Asian LD blocks. The corrected run identified 150
  locally heritable HELIOS blocks, 1,120 locally heritable education blocks
  and 92 blocks eligible for the HELIOS–education analysis.
- Among the 92 eligible blocks, 47 local covariance estimates were positive
  and 45 were negative. Fifteen local correlations had nominal P < 0.05; none
  survived BH correction across the 92 tests.
- Positive covariance summed to 0.03560 and the absolute negative contribution
  to 0.04282. The signed sum was -0.00722, giving a descriptive cancellation
  fraction of 90.8%.
- Leave-one-chromosome-out cancellation fractions ranged from 86.9% to 97.6%.
  Excluding the MHC gave 91.7%; excluding the largest absolute contribution
  gave 93.8%.

These results are consistent with substantial cancellation among the analysed
regions. They do not establish cancellation as the cause of the near-zero
genome-wide estimate: only 92 of 3,064 blocks were eligible, they represented
11.1% of the summed estimable HELIOS local h2 and 6.48% of the corresponding
education total, and no local association survived FDR correction. The
eligible-region aggregate rg was -0.0171 and should not be treated as a
genome-wide estimate.

The concise report and suggested manuscript wording are in
[`final/audit_report.md`](final/audit_report.md).

## Corrected LAVA run

The recovered historical configuration had declared quantitative KoGES
education as binary. Two isolated four-trait LAVA 0.1.5 runs were therefore
performed with identical summary statistics, overlap matrix, LD reference,
3,064 blocks and univariate threshold (P < 2e-5):

1. a historical control retaining the binary declaration; and
2. a corrected run declaring education as quantitative.

The control reproduced all 11,057 univariate rows, 726 bivariate rows,
univariate estimates and bivariate point estimates exactly. LAVA 0.1.5 derives
bivariate confidence intervals and P values from unseeded Monte Carlo draws;
those fields were not byte-identical, but their correlations exceeded 0.9986
and all nominal and FDR decisions were unchanged.

The corrected run produced 11,051 univariate and 717 all-pair bivariate rows.
All 92 HELIOS–education blocks remained eligible; their rho estimates
correlated 0.9997 with the control, no sign changed, one nominal P-value
decision changed, and no FDR decision changed. All reported results use the
corrected run.

## Directory guide

- [`00_discovery`](00_discovery): source inventory and input concordance.
- [`01_global_ldsc`](01_global_ldsc): global LDSC logs and summary.
- [`02_lava_inputs`](02_lava_inputs): historical configuration, controlled
  rerun inputs, exact outputs, logs and R session information.
- [`03_local_covariance`](03_local_covariance): corrected region-level results,
  summaries and ranked regions.
- [`04_cancellation`](04_cancellation): covariance decomposition.
- [`05_sensitivity`](05_sensitivity): chromosome, MHC and outlier exclusions.
- [`06_plots`](06_plots): figures generated from the corrected output.
- [`final`](final): key values, post-processing log, software session and final
  report.
- [`provenance`](provenance): comparisons, run record and SHA-256 manifests.
- [`scripts`](scripts): validation, comparison and post-processing code.

The summary-statistic inputs and LD reference are controlled or too large for
this repository. Their paths and hashes are recorded in
[`provenance/exact_run_input_checksums.sha256`](provenance/exact_run_input_checksums.sha256).

## Reproduction

On the analysis server, run:

```bash
bash 02_lava_inputs/rerun/run_control_and_corrected.sh /path/to/new_run
bash scripts/run_postprocessing.sh /path/to/new_run /path/to/new_analysis
```

`run_postprocessing.sh` stops unless the historical control passes the
reproduction gate. The scripts contain the server-side input and reference
paths used for this audit.
