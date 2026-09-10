# Project 3: multivariate pleiotropy and directional locus classification

[![Project 3 checks](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml/badge.svg)](https://github.com/TimvanderEs/PhD/actions/workflows/project3-synthetic-test.yml)

This directory contains analysis code supporting an upcoming Molecular Psychiatry manuscript. It includes the ancestry-specific LDSC, LAVA, and PLEIO preparation/run workflows; the post-hoc directional classification and single-trait overlap audit of PLEIO/FUMA loci; FUMA/MAGMA, GTEx v8, and g:Profiler enrichment analyses; and the SMR/HEIDI threshold and sensitivity reanalysis.

## Reproducibility scope

The repository separates analyses that can be rebuilt entirely from the public
archive from steps that require provider-controlled GWAS or webSMR inputs.

| Analysis | Public record included here | What a reader can reproduce |
| --- | --- | --- |
| Directional PLEIO/FUMA classification | Complete R workflow, command-line wrapper, synthetic inputs, and expected outputs | Run the full synthetic example immediately; rerun the study analysis after supplying the authorised PLEIO/FUMA and summary-statistics inputs. |
| Summary-statistics QC | Recovered settings and checksummed run-log provenance for five of eight ancestry-by-trait runs | Audit recovered thresholds and counts; three missing logs and two EAS sample-definition discrepancies are explicitly identified. |
| LDSC and PLEIO | Executed scripts, trait/model configurations, selected non-individual-level matrices and calibrations, and final-output fingerprints | Recreate the workflows after obtaining the source GWAS and required reference resources. |
| LAVA | Portable runner, manuscript-target and recovered input-info files, exact-threshold/FDR post-processing, and archived-output fingerprints | Audit the historical result now; verify the exact EAS SCZ input identity and rerun. The archived 726, exact-threshold 698, and workbook 697 EAS row sets are not conflated. |
| SMR/HEIDI | Complete threshold/sensitivity reanalysis code, input contract, and run provenance | Rerun after supplying the seven merged webSMR files listed in its guide. |
| FUMA/MAGMA and GTEx v8 | Exact FUMA parameters, raw term-level results, reviewed manifests, and processed tables | Rebuild and compare all 13 archived runs directly from this repository. |
| g:Profiler | Exact query lists, available custom backgrounds, original exports, checksums, and processed tables | Validate all 12 archived queries and rebuild the 61 significant-term table without calling a changing live database. |

Unrecovered or conflicting provenance is documented rather than silently
reconstructed. This includes the EAS cognitive-function single-trait FUMA run,
three SumstatsQC logs, two EAS sample-definition discrepancies, the original
g:Profiler database release and EAS background option, and the webSMR portal
job/version identifiers. The LAVA guide also records the clean-rerun blocker.

The badge above runs the directional and LAVA synthetic workflows, repository
consistency checks, and archived FUMA/MAGMA and g:Profiler rebuilds on every
Project 3 change.

## Analysis guides

- [Directional PLEIO/FUMA classification](scripts/classify_pleio_locus_direction.R)
- [Summary-statistics QC](scripts/sumstats_qc/README.md)
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
| `inputs/` | FUMA/MAGMA, g:Profiler, and SMR/HEIDI run manifests plus retained query/background records. |
| `results/` | Archived FUMA/MAGMA and g:Profiler raw exports and compact processed results. |
| `provenance/` | Software versions, checksums, final-output fingerprints, and expected run summaries. |
| `examples/` and `tests/` | Privacy-safe synthetic inputs plus automated integrity and workflow checks. |

## Directional classification

For each independent significant SNP, the script applies an absolute Z-score threshold (default: `1.96`) and determines the sign represented in each domain:

- cognitive/educational domain: educational attainment (`EA_Z`) and cognitive function (`CF_Z`);
- psychiatric domain: MDD, SCZ, or both, according to the PLEIO model;
- opposite cognitive and psychiatric signs are labelled `concordant_raw` under the phenotype/sign convention used for these analyses;
- matching signs are labelled `discordant_raw`;
- contradictory signs within a domain are labelled `mixed_within_domain`;
- a domain without a qualifying Z score is labelled `Unassigned`.

SNP-level labels are collapsed to the manuscript locus classes `Concordant`,
`Discordant`, `Dual`, `Mixed`, and `Unassigned`. A Mixed SNP takes precedence;
otherwise, a locus containing both Concordant and Discordant SNPs is Dual, a
locus containing only one informative class inherits that class, and all
remaining loci are Unassigned.

The prioritisation step compares a PLEIO locus only with the single-trait FUMA
results included in that model. A locus with no overlap is marked
`pleio_prioritised`; this does not imply complete biological novelty.

## Terminology and supplementary material

The repository uses the manuscript abbreviations throughout: `EA` means
educational attainment and `CF` means cognitive function. Some original source
filenames and upstream phenotype identifiers contain `Edu`, `G`, `COGENT`, or
`HEL10k`; those identifiers are retained only where needed to locate or verify
the executed source files. Reader-facing outputs use `EA` and `CF`.

The directional locus outputs correspond to Supplementary Tables S1-S3, the
sensitivity and audit outputs are reported in Supplementary Tables S7-S9, and
Supplementary Table S8 is specifically the locus-level directional audit.
FUMA/MAGMA and GENE2FUNC results are reported in Supplementary Tables S10-S12,
and g:Profiler results in Supplementary Tables S13-S14 and Supplementary Figure
S7.

## Software requirements

The privacy-safe synthetic check requires Bash, R, and the R package
`data.table`. The archived enrichment checks additionally require Python 3.9
or later, NumPy, and pandas. Install these public-validation dependencies with:

```bash
Rscript -e 'install.packages("data.table")'
python3 -m pip install -r scripts/fuma_magma/requirements.txt
```

Rebuilding Supplementary Figure S7 also requires the R packages `ggplot2` and
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
- [Summary-statistics QC](scripts/sumstats_qc/README.md)
- [LAVA](scripts/lava/README.md)
- [PLEIO preparation and model runs](scripts/pleio/README.md)
- [SMR/HEIDI](scripts/smr_heidi/README.md)
- [FUMA/MAGMA and GTEx v8 enrichment](scripts/fuma_magma/README.md)
- [g:Profiler over-representation analysis](scripts/gprofiler/README.md)

The ancestry/model configurations preserve the executed trait order, overlap
matrices, and output names. Where an executed configuration conflicts with the
manuscript, the historical and corrected-target files are kept separately. The
[provenance record](REPRODUCIBILITY.md) gives software revisions, run dates,
row counts, SHA-256 fingerprints, and verification status.

## Running the directional classifier on study outputs

```bash
bash scripts/run_pleio_locus_classification.sh \
  EUR \
  FOURTRAIT \
  /path/to/FUMA/EUR/PLEIO \
  /path/to/FUMA/EUR \
  /path/to/raw_merged_PLEIO.tsv.gz \
  /path/to/output/EUR_FOURTRAIT
```

Arguments, in order, are:

1. ancestry label;
2. model: `FOURTRAIT`, `MDD_3TRAIT`, or `SCZ_3TRAIT`;
3. PLEIO FUMA directory containing one `GenomicRiskLoci*.txt` and one `IndSigSNPs*.txt` file;
4. single-trait FUMA root containing the model-relevant `EA`, `CF`, `MDD`, and/or `SCZ` subdirectories;
5. raw merged PLEIO summary-statistics file;
6. output directory;
7. optional source folder name for EA (default: `EA`);
8. optional source folder name for CF (default: `CF`).

The optional folder arguments support historical archives whose directory
names differ from the publication labels; they do not change output variable
names.

The R script can also be called directly. Run the following for its argument summary:

```bash
Rscript scripts/classify_pleio_locus_direction.R --help
```

### Required input fields

| Input | Required fields |
| --- | --- |
| PLEIO and single-trait `GenomicRiskLoci` files | chromosome, start, and end; a genomic-locus identifier is recommended |
| PLEIO `IndSigSNPs` file | SNP identifier; a genomic-locus identifier or genomic coordinates |
| Raw merged PLEIO file | SNP identifier, `EA_Z`, `CF_Z`, and the psychiatric Z-score column(s) used by the selected model |

Several source-pipeline aliases are accepted and normalised to `EA_Z` and
`CF_Z`. If a required field is absent or ambiguous, the script stops with a
diagnostic message.

## Directional-classification outputs

Each filename begins with `<ANCESTRY>_<MODEL>`:

- `_PLEIO_locus_master.tsv`: locus classes, class counts, single-trait overlaps, and PLEIO-prioritisation fields;
- `_PLEIO_IndSigSNPs_with_direction.tsv`: SNP-level Z scores and direction labels;
- `_PLEIO_locus_direction_summary.tsv`: counts and proportions for each manuscript locus class;
- `_PLEIO_locus_prioritisation_summary.tsv`: counts of overlaps and PLEIO-prioritised loci;
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

See the manuscript-ready [Code Availability statement](CODE_AVAILABILITY.md).
The manuscript citation and archival DOI will be added when available. Until
then, cite this repository and the specific Git commit used for an analysis.
The repository is distributed under the [MIT License](../LICENSE).
