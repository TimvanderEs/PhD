# Project 3: multivariate pleiotropy and directional locus classification

[![Project 3 synthetic test](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml/badge.svg)](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml)

This directory contains analysis code supporting an upcoming molecular psychiatry manuscript. The current release includes the post-hoc directional classification and single-trait overlap audit of PLEIO/FUMA loci, together with the SMR/HEIDI threshold and sensitivity reanalysis. Additional analysis scripts will be added as the manuscript workflow is finalised.

## Contents

| Path | Purpose |
| --- | --- |
| `scripts/classify_PLEIO_FUMA_loci_v2.R` | Classifies independent significant SNPs and collapses those labels to PLEIO loci; compares corrected and legacy classifications; audits overlap with model-relevant single-trait FUMA loci. |
| `scripts/run_PLEIO_reclassification_audit.sh` | Command-line wrapper that runs the R workflow and prints the main summaries and QC results. |
| `scripts/smr_heidi/` | Reanalyses merged webSMR/HEIDI outputs using trait-specific and study-wide multiple-testing thresholds, sensitivity definitions, QC checks, and optional PLEIO/LAVA integration. |
| `examples/synthetic/` | Small, artificial input files that document the required directory layout and column names. They contain no study data. |
| `examples/run_synthetic_example.sh` | Runs the complete workflow on the synthetic example. |
| `tests/validate_synthetic_output.R` | Checks the synthetic run against its expected locus labels, overlaps, and QC values. |

## Directional classification

For each independent significant SNP, the script applies an absolute Z-score threshold (default: `1.96`) and determines the sign represented in each domain:

- cognitive domain: educational attainment (`EDU_Z`) and cognitive performance (`G_Z` or an accepted alias);
- psychiatric domain: MDD, SCZ, or both, according to the PLEIO model;
- opposite cognitive and psychiatric signs are labelled `concordant_raw` under the phenotype/sign convention used for these analyses;
- matching signs are labelled `discordant_raw`;
- contradictory signs within a domain are labelled `mixed_within_domain`;
- a domain without a qualifying Z score is labelled `trait_dominant_or_unclear`.

SNP-level labels are collapsed to each locus as `concordant`, `discordant`, `mixed`, `mixed_dual`, or `trait_dominant_or_unclear`. Full definitions and the collapse precedence are implemented in the R script.

The corrected prioritisation audit compares a PLEIO locus only with the single-trait FUMA results included in that model. The output retains the earlier four-trait/same-sign classification so that changes can be inspected explicitly.

## Requirements

- R with the `data.table` package (synthetic workflow validated with R 4.2.3)
- Python 3.8 or later (SMR/HEIDI workflow; no third-party Python packages)
- Bash for the convenience wrapper

Install the R dependency with:

```r
install.packages("data.table")
```

## Quick validation with synthetic data

From this `Project3` directory, run:

```bash
bash examples/run_synthetic_example.sh
```

The example writes ignored output files to `examples/synthetic/output/`. The input identifiers and values are entirely synthetic and are provided only to demonstrate the workflow.

The same example and its expected outputs are checked automatically on GitHub whenever Project 3 files change.

## Running the workflow on study outputs

```bash
bash scripts/run_PLEIO_reclassification_audit.sh \
  EUR \
  FOURTRAIT \
  /path/to/FUMA/EUR/PLEIO \
  /path/to/FUMA/EUR \
  /path/to/raw_merged_PLEIO.tsv.gz \
  /path/to/output/EUR_FOURTRAIT \
  G
```

Arguments, in order, are:

1. ancestry label;
2. model: `FOURTRAIT`, `MDD_3TRAIT`, or `SCZ_3TRAIT`;
3. PLEIO FUMA directory containing one `GenomicRiskLoci*.txt` and one `IndSigSNPs*.txt` file;
4. single-trait FUMA root containing `EDU`, `G` (or the name supplied as argument 7), `MDD`, and `SCZ` subdirectories;
5. raw merged PLEIO summary-statistics file;
6. output directory;
7. optional cognitive-factor folder name (default: `G`).

The R script can also be called directly. Run the following for its argument summary:

```bash
Rscript scripts/classify_PLEIO_FUMA_loci_v2.R --help
```

## SMR/HEIDI reanalysis

The SMR workflow replaces a single fixed primary threshold with a
trait-specific Bonferroni threshold while retaining `P_SMR < 5e-8` and other
HEIDI/instrument definitions as sensitivity analyses. It can optionally rebuild
PLEIO–SMR and PLEIO–LAVA convergence tables from the canonical supplementary
TSVs. See [`scripts/smr_heidi/README.md`](scripts/smr_heidi/README.md) for the
input layout, command, outputs, and interpretation caveat.

### Required input fields

| Input | Required fields |
| --- | --- |
| PLEIO and single-trait `GenomicRiskLoci` files | chromosome, start, and end; a genomic-locus identifier is recommended |
| PLEIO `IndSigSNPs` file | SNP identifier; a genomic-locus identifier or genomic coordinates |
| Raw merged PLEIO file | SNP identifier, `EDU_Z`, a cognitive Z-score column, and the psychiatric Z-score column(s) used by the selected model |

Several common FUMA and analysis-pipeline column aliases are accepted. If a required field is absent or ambiguous, the script stops with a diagnostic message.

## Main outputs

Each filename begins with `<ANCESTRY>_<MODEL>`:

- `_PLEIO_locus_master.tsv`: corrected locus labels, single-trait overlaps, and retained legacy/audit fields;
- `_PLEIO_IndSigSNPs_with_direction.tsv`: SNP-level Z scores and direction labels;
- `_PLEIO_locus_direction_summary.tsv`: corrected counts by direction and prioritisation status;
- `_locus_direction_transition.tsv`: legacy-to-corrected locus transitions;
- `_prioritisation_transition.tsv`: changes caused by using model-relevant comparators;
- `_classification_QC.tsv`: run-level counts and missingness checks;
- `_run.log`: console output captured by the Bash wrapper.

## Data availability and privacy

No individual-level participant data, credentials, or study summary statistics are included here. The synthetic example is not derived from the study. Access to the source GWAS data and downstream FUMA files remains subject to the terms of their original providers.

## Citation and licence

The manuscript citation and DOI will be added when available. Until then, please cite this repository and the specific Git commit used for an analysis. The repository is distributed under the [MIT License](../LICENSE).
