# Reproducibility

The directional locus classifier and the LAVA post-processing code can be
tested with the synthetic data in this repository. The retained FUMA, MAGMA,
GTEx and g:Profiler outputs can also be checked and rebuilt locally. These
tests run automatically through GitHub Actions.

LDSC, LAVA, PLEIO and the other study-level analyses require the original GWAS
summary statistics and ancestry-matched reference data. These inputs are not
redistributed. Required file formats and commands are documented with each
workflow under [`scripts/`](scripts/).

Recorded software versions and checksums are provided in
[`provenance/`](provenance/).
