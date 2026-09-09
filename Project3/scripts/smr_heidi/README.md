# SMR/HEIDI reanalysis

This workflow reanalyses merged webSMR/HEIDI outputs using a trait-specific
Bonferroni threshold and retains the original `P_SMR < 5e-8` definition as an
ultra-stringent sensitivity analysis. It uses only the Python standard library.

## Primary definition

For each ancestry-trait analysis, a primary hit satisfies:

- `P_SMR < 0.05 / N`, where `N` is the number of valid merged
  probe-by-QTL-resource tests for that ancestry and trait; and
- `P_HEIDI > 0.01`.

The output also reports a study-wide Bonferroni tier, `P_HEIDI >= 0.05`, at
least 10 HEIDI SNPs, SMR-multi support, and a unique-probe Bonferroni
sensitivity tier. These classifications aid sensitivity assessment; they do
not by themselves establish expression mediation or a shared causal variant.

## Input layout

`<smr_root>` must contain these seven merged results at the following relative
paths:

```text
EU_webSMR/EU_EDU_SMR_merged_raw.zip
EU_webSMR/EU_COG_SMR_merged_raw.zip
EU_webSMR/EU_MDD_SMR_merged_raw.zip
EU_webSMR/EU_SCZ_SMR_merged_raw.zip
EAS_SMR_web/EDU_EAS_eSMR.merged.tsv
EAS_SMR_web/EAS_MDD_eSMR.merged.tsv
EAS_SMR_web/SCZ_EAS_eSMR.merged.tsv
```

The optional supplementary TSV directory may contain:

```text
12_PLEIO_locus_master.tsv
14_PLEIO_gene_direction_map.tsv
11_LAVA_PLEIO_overlap.tsv
```

When the first two optional files are absent, the SMR/HEIDI reanalysis still
completes and records the skipped integration in the QC table.

## Run

```bash
bash scripts/smr_heidi/run_smr_reanalysis.sh \
  /path/to/SMR \
  /path/to/output/SMR_REANALYSIS \
  /path/to/canonical_supplementary_tsv
```

The third argument can be omitted when the PLEIO/LAVA integration tables are
not required. The Python program can also be called directly; use `--help` for
its command-line options.

## Main outputs

- `00_SMR_reanalysis_QC.tsv`: input, duplication, and integration checks.
- `20_SMR_HEIDI_all_results_reanalysed.tsv.gz`: harmonised results and all
  threshold/sensitivity flags.
- `21_SMR_threshold_comparison.tsv`: counts under the primary and sensitivity
  definitions.
- `22A`–`22D`: hit-level tables for each significance tier.
- `23A`–`23C`: unique-gene and cross-ancestry summaries.
- `24A`–`24C`: optional PLEIO/SMR/LAVA convergence summaries.
- `99_SMR_reanalysis_checks.tsv`: copy of the complete run checks.

The EC2 analysis run completed all recorded QC checks, including confirmation
that the EAS MDD and SCZ source files were not duplicates.
