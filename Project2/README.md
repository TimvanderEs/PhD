# Cross-ancestry multivariate psychiatric genetics

[![Reproducibility checks](https://github.com/TimvanderEs/PhD/actions/workflows/project2-reproducibility.yml/badge.svg)](https://github.com/TimvanderEs/PhD/actions/workflows/project2-reproducibility.yml)

> **WCPG 2026:** [poster results and code guide](WCPG_2026.md)

Code and non-sensitive results for analyses of genetic overlap between major
depression, schizophrenia, educational attainment and cognitive function in
European- and East Asian-ancestry data.

## Analysis code

- [Summary-statistics quality control](scripts/sumstats_qc/README.md)
- [LDSC and GenomicSEM](scripts/ldsc/README.md)
- [LAVA](scripts/lava/README.md)
- [PLEIO](scripts/pleio/README.md)
- [Directional locus classification](scripts/classify_pleio_locus_direction.R)
- [SMR and HEIDI](scripts/smr_heidi/README.md)
- [FUMA, MAGMA and GTEx](scripts/fuma_magma/README.md)
- [g:Profiler](scripts/gprofiler/README.md)

The directional classification script assigns SNP-level labels and summarises
them as Concordant, Discordant, Dual, Mixed or Unassigned loci. Its outputs
support Supplementary Tables S1-S3 and S7-S9; Supplementary Table S8 contains
the locus-level classifications. FUMA, MAGMA and GENE2FUNC outputs
support Supplementary Tables S10-S12; g:Profiler outputs support Supplementary
Tables S13-S14 and Supplementary Figure S7.

## Repository structure

- `scripts/`: analysis code and configuration files
- `inputs/`: run manifests and non-sensitive inputs
- `results/`: retained functional-analysis outputs
- `provenance/`: software versions and file checksums
- `examples/`: synthetic data for the directional classifier
- `tests/`: workflow and file-integrity checks

## Tests

From the repository root:

```bash
bash Project2/tests/run_synthetic_test.sh
python3 Project2/tests/test_lava_postprocessing.py
python3 Project2/tests/test_repository_consistency.py
```

Re-running the study analyses requires the source GWAS summary statistics and
reference data listed in the workflow guides.

## Data and code availability

No individual-level data or provider-controlled GWAS summary statistics are
stored here. See [Data and code availability](DATA_AND_CODE_AVAILABILITY.md).

Reader-facing files use `EA` for educational attainment and `CF` for cognitive
function. Older labels remain only where they form part of an original input
or output filename.

## Citation

The publication release will be archived with a DOI. Until then, cite the
repository commit used for the analysis. Code is released under the
[MIT License](../LICENSE).
