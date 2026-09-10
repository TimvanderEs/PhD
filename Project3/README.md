# Project 3: multivariate pleiotropy and directional locus classification

[![Project 3 checks](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml/badge.svg)](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml)

This directory contains analysis code supporting an upcoming Molecular Psychiatry manuscript. It includes the ancestry-specific LDSC, LAVA, and PLEIO preparation/run workflows; the post-hoc directional classification and single-trait overlap audit of PLEIO/FUMA loci; FUMA/MAGMA, GTEx v8, and g:Profiler enrichment analyses; and the SMR/HEIDI threshold and sensitivity reanalysis.

## Reproducibility scope

The repository separates analyses that can be rebuilt entirely from the public
archive from steps that require provider-controlled GWAS or webSMR inputs.

| Analysis | Public record included here | What a reader can reproduce |
| --- | --- | --- |
| Directional PLEIO/FUMA classification | Complete R workflow, command-line wrapper, synthetic inputs, and expected outputs | Run the full synthetic example immediately; rerun the study analysis after supplying the authorised PLEIO/FUMA and summary-statistics inputs. |
| LDSC, LAVA, and PLEIO | Executed scripts, trait/model configurations, selected non-individual-level matrices and calibrations, and final-output fingerprints | Recreate the workflows after obtaining the source GWAS and required reference resources. |
| SMR/HEIDI | Complete threshold/sensitivity reanalysis code, input contract, and run provenance | Rerun after supplying the seven merged webSMR files listed in its guide. |
| FUMA/MAGMA and GTEx v8 | Exact FUMA parameters, raw term-level results, reviewed manifests, and processed tables | Rebuild and compare all 13 archived runs directly from this repository. |
| g:Profiler | Exact query lists, available custom backgrounds, original exports, checksums, and processed tables | Validate all 12 archived queries and rebuild the 61 significant-term table without calling a changing live database. |

Two gaps are documented rather than silently reconstructed: no EAS
cognitive-performance single-trait FUMA run was recovered, and the original
g:Profiler database release and EAS background option were not recorded. The
relevant manifests and workflow guides mark these limitations explicitly.

The badge above runs the synthetic workflow and rebuilds the archived
FUMA/MAGMA and g:Profiler results on every Project 3 change.

## Analysis guides

- [Directional PLEIO/FUMA classification](scripts/classify_PLEIO_FUMA_loci_v2.R)
- [LDSC](scripts/ldsc/README.md)
- [LAVA](scripts/lava/README.md)
- [PLEIO preparation and model runs](scripts/pleio/README.md)
- [SMR/HEIDI](scripts/smr_heidi/README.md)
- [FUMA/MAGMA and GTEx v8 enrichment](scripts/fuma_magma/README.md)
- [g:Profiler over-representation analysis](scripts/gprofiler/README.md)

## Repository layout

| Path | Contents |
| --- | --- |
| `scripts/` | Analysis code, per-workflow instructions, upstream configurations, and selected workflow inputs. |
| `inputs/` | FUMA/MAGMA run manifests and g:Profiler query, background, and settings records. |
| `results/` | Archived FUMA/MAGMA and g:Profiler raw exports and compact processed results. |
| `provenance/` | Software versions, checksums, final-output fingerprints, and expected run summaries. |
| `examples/` and `tests/` | Privacy-safe synthetic inputs plus automated integrity and workflow checks. |

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

## Software requirements

The privacy-safe synthetic check requires Bash, R, and the R package
`data.table`. The archived enrichment checks additionally require Python 3.9
or later, NumPy, and pandas. Install these public-validation dependencies with:

```bash
Rscript -e 'install.packages("data.table")'
python3 -m pip install -r scripts/fuma_magma/requirements.txt
```

Rebuilding Supplementary Figure S8 also requires the R packages `ggplot2` and
`patchwork`. Upstream reruns require GenomicSEM, LAVA 0.1.5, and PLEIO at the
revision recorded in its guide. Each analysis guide lists its additional
reference files and workflow-specific requirements.

## Quick validation with synthetic data

From this `Project3` directory, run:

```bash
bash examples/run_synthetic_example.sh
```

The example writes ignored output files to `examples/synthetic/output/`. The input identifiers and values are entirely synthetic and are provided only to demonstrate the workflow.

The same example and its expected outputs are checked automatically on GitHub whenever Project 3 files change.

## Reproducing the upstream analyses

Each upstream workflow has a self-contained guide:

- [LDSC](scripts/ldsc/README.md)
- [LAVA](scripts/lava/README.md)
- [PLEIO preparation and model runs](scripts/pleio/README.md)
- [SMR/HEIDI](scripts/smr_heidi/README.md)
- [FUMA/MAGMA and GTEx v8 enrichment](scripts/fuma_magma/README.md)
- [g:Profiler over-representation analysis](scripts/gprofiler/README.md)

The ancestry/model configurations preserve the executed trait order, disease
prevalences, thresholds, overlap matrices, and output names. The
[provenance record](REPRODUCIBILITY.md) gives software revisions, successful
run dates, row counts, and SHA-256 fingerprints for the retained final outputs.

## Running the directional classifier on study outputs

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

### Required input fields

| Input | Required fields |
| --- | --- |
| PLEIO and single-trait `GenomicRiskLoci` files | chromosome, start, and end; a genomic-locus identifier is recommended |
| PLEIO `IndSigSNPs` file | SNP identifier; a genomic-locus identifier or genomic coordinates |
| Raw merged PLEIO file | SNP identifier, `EDU_Z`, a cognitive Z-score column, and the psychiatric Z-score column(s) used by the selected model |

Several common FUMA and analysis-pipeline column aliases are accepted. If a required field is absent or ambiguous, the script stops with a diagnostic message.

## Directional-classification outputs

Each filename begins with `<ANCESTRY>_<MODEL>`:

- `_PLEIO_locus_master.tsv`: corrected locus labels, single-trait overlaps, and retained legacy/audit fields;
- `_PLEIO_IndSigSNPs_with_direction.tsv`: SNP-level Z scores and direction labels;
- `_PLEIO_locus_direction_summary.tsv`: corrected counts by direction and prioritisation status;
- `_locus_direction_transition.tsv`: legacy-to-corrected locus transitions;
- `_prioritisation_transition.tsv`: changes caused by using model-relevant comparators;
- `_classification_QC.tsv`: run-level counts and missingness checks;
- `_run.log`: console output captured by the Bash wrapper.

## SMR/HEIDI reanalysis

The SMR workflow replaces a single fixed primary threshold with a
trait-specific Bonferroni threshold while retaining `P_SMR < 5e-8` and other
HEIDI/instrument definitions as sensitivity analyses. It can optionally rebuild
PLEIO–SMR and PLEIO–LAVA convergence tables from the canonical supplementary
TSVs. See [`scripts/smr_heidi/README.md`](scripts/smr_heidi/README.md) for the
input layout, command, outputs, and interpretation caveat.

## Data availability and privacy

No individual-level participant data, credentials, or source GWAS summary statistics are included here. The synthetic example is not derived from the study. The included LAVA overlap matrices, locus definitions, PLEIO inverse-CDF calibrations, FUMA parameter files, MAGMA/GTEx enrichment outputs, g:Profiler gene lists/exports, and compact summaries are non-individual-level analysis metadata or derived results. Access to source GWAS data and any unarchived FUMA output remains subject to the terms of the original providers.

## Citation and licence

The manuscript citation and DOI will be added when available. Until then, please cite this repository and the specific Git commit used for an analysis. The repository is distributed under the [MIT License](../LICENSE).
