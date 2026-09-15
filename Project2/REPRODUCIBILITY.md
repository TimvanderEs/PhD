# Reproducibility record

This document describes which analyses can be run from the public repository
and which require controlled or third-party inputs. Software versions, run
manifests and file fingerprints are retained under [`provenance/`](provenance/).

## Analysis status

| Workflow | Material retained | Reproduction status |
| --- | --- | --- |
| SumstatsQC | Settings and five recovered run logs | Requires source summary statistics; three logs were not recovered |
| LDSC/GenomicSEM | Scripts, trait manifests and selected covariance/correlation matrices | Requires source summary statistics and LD-score reference files |
| LAVA | Runner, overlap matrices, locus definitions, input manifests and post-processing code | Corrected rerun required; see below |
| PLEIO | Four- and three-trait preparation scripts, model manifests and importance-sampling files | Requires source summary statistics and the cited PLEIO software revision |
| Directional locus classification | Complete workflow, synthetic data and expected outputs | Fully testable from this repository |
| SMR/HEIDI | Reanalysis code and checksummed input manifest | Requires the seven merged webSMR result files listed in the workflow guide |
| FUMA/MAGMA/GTEx | Parameters and raw term-level results for 13 archived runs | Rebuildable from the archived files |
| g:Profiler | Query lists, retained backgrounds and 12 raw exports | Rebuildable from the archived files; 61 significant terms are expected |

## LAVA status

The archived East Asian bivariate file contains 726 rows. It was generated with
a rounded `2e-5` univariate screen and a recovered input configuration that does
not match the manuscript sample definitions. Applying the exact threshold
`P <= 0.05/3064` to the archived univariate and bivariate files retains 698
rows. The current manuscript workbook contains 697 rows because one eligible
EA-SCZ row with a zero univariate P value is absent.

These are distinct records rather than alternative labels for the same output.
The exact East Asian SCZ input must be confirmed and the LAVA analysis rerun
with the manuscript-target configuration before the final denominator is
reported. The corresponding files are documented in
[`scripts/lava/README.md`](scripts/lava/README.md).

## Other provenance limitations

- Five of eight SumstatsQC logs were recovered. EUR EA, EUR MDD and EAS CF logs
  are not available.
- The manuscript reports 22,778 East Asian SCZ cases, whereas a recovered QC
  log reports 27,888. The exact input identity remains to be confirmed.
- A single-trait East Asian cognitive-function FUMA run was not recovered;
  comparisons that require it are marked incomplete.
- The original g:Profiler database release and East Asian background setting
  were not retained. Raw exports are archived so the reported tables can be
  reconstructed without querying the current live database.
- Original webSMR job identifiers, portal version and access dates were not
  recovered. Available file timestamps are not presented as portal run dates.

## Provenance files

| File | Contents |
| --- | --- |
| [`software_versions.tsv`](provenance/software_versions.tsv) | Recorded software and runtime versions |
| [`final_output_fingerprints.tsv`](provenance/final_output_fingerprints.tsv) | Output names, row counts, completion times, checksums and verification status |
| [`lava_pair_test_counts.tsv`](provenance/lava_pair_test_counts.tsv) | LAVA row counts by trait pair and analysis set |
| [`sumstatsqc_runs.tsv`](provenance/sumstatsqc_runs.tsv) | Recovered QC settings, counts and log hashes |
| [`eur_pleio_alignment_fingerprints.tsv`](provenance/eur_pleio_alignment_fingerprints.tsv) | Variant counts and hashes for the aligned EUR PLEIO inputs |
| [`pleio_isf_fingerprints.tsv`](provenance/pleio_isf_fingerprints.tsv) | Fingerprints for PLEIO importance-sampling files |
| [`pleio_gws_counts.tsv`](provenance/pleio_gws_counts.tsv) | Genome-wide significant PLEIO SNP counts |
| [`fuma_magma_raw_files.sha256`](provenance/fuma_magma_raw_files.sha256) | Checksums for archived FUMA/MAGMA inputs and outputs |

Historical identifiers such as `Edu`, `G`, `COGENT` and `HEL10k` remain in a
small number of source filenames. Reader-facing labels use `EA` and `CF`.
