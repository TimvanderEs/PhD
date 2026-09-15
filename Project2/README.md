# Project 2: cross-ancestry multivariate psychiatric genetics

[![Reproducibility checks](https://github.com/TimvanderEs/PhD/actions/workflows/project2-reproducibility.yml/badge.svg)](https://github.com/TimvanderEs/PhD/actions/workflows/project2-reproducibility.yml)

This directory contains the analysis code and non-sensitive derived material
for a study of genetic overlap between major depression, schizophrenia,
educational attainment and cognitive function in European- and East
Asian-ancestry data. The work includes genome-wide and local genetic
correlation analyses, multivariate association testing, post-hoc directional
classification and functional follow-up.

## Analyses

- summary-statistics quality control;
- LDSC and GenomicSEM;
- LAVA local genetic correlation;
- four- and three-trait PLEIO models;
- directional classification and single-trait overlap of PLEIO/FUMA loci;
- SMR/HEIDI sensitivity analyses;
- FUMA, MAGMA and GTEx enrichment analyses; and
- g:Profiler over-representation analyses.

## Repository contents

| Path | Description |
| --- | --- |
| [`scripts/`](scripts/) | Analysis scripts, configuration files and workflow-specific instructions |
| [`inputs/`](inputs/) | Run manifests and non-sensitive query/background files |
| [`results/`](results/) | Archived FUMA/MAGMA and g:Profiler exports and processed tables |
| [`provenance/`](provenance/) | Software versions, checksums and output fingerprints |
| [`examples/`](examples/) | Privacy-safe synthetic input data |
| [`tests/`](tests/) | Automated workflow and archive-integrity checks |

The principal custom analysis is
[`scripts/classify_pleio_locus_direction.R`](scripts/classify_pleio_locus_direction.R).
It assigns SNP-level directional labels and collapses them to the locus classes
reported in the manuscript. Detailed instructions are provided alongside each
workflow:

- [Summary-statistics QC](scripts/sumstats_qc/README.md)
- [LDSC](scripts/ldsc/README.md)
- [LAVA](scripts/lava/README.md)
- [PLEIO](scripts/pleio/README.md)
- [SMR/HEIDI](scripts/smr_heidi/README.md)
- [FUMA/MAGMA and GTEx](scripts/fuma_magma/README.md)
- [g:Profiler](scripts/gprofiler/README.md)

Directional locus outputs support Supplementary Tables S1-S3 and S7-S9;
Supplementary Table S8 contains the locus-level directional classification. FUMA/MAGMA and
GENE2FUNC results support Supplementary Tables S10-S12, while the g:Profiler
results support Supplementary Tables S13-S14 and Supplementary Figure S7.

## Reproducibility

The repository can be tested without access to study data:

```bash
bash Project2/tests/run_synthetic_test.sh
python3 Project2/tests/test_lava_postprocessing.py
python3 Project2/tests/test_repository_consistency.py
```

The continuous-integration workflow additionally rebuilds the archived
FUMA/MAGMA and g:Profiler summaries. Re-running the study analyses requires the
source GWAS summary statistics and reference panels described in the individual
workflow guides.

The archived East Asian LAVA run requires a corrected rerun before its test
denominator is treated as final. The repository therefore records the original
726-row output, the 698 rows retained by the exact threshold, and the 697 rows
in the current manuscript workbook as separate objects. See
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the full status record.

Publication-facing files use `EA` for educational attainment and `CF` for
cognitive function. Historical labels are retained only in source filenames
where they are needed to identify an executed analysis.

## Data and code availability

Individual-level data and provider-controlled GWAS summary statistics are not
redistributed. The repository contains custom code, synthetic examples,
analysis configurations, run manifests, checksums and eligible derived results.
The full availability statement is in
[`DATA_AND_CODE_AVAILABILITY.md`](DATA_AND_CODE_AVAILABILITY.md).

## Citation and licence

Until a versioned archive and DOI are available, cite this repository together
with the Git commit used for the analysis. Code in this repository is released
under the [MIT License](../LICENSE); third-party data and software remain subject
to their original terms.
