# LDSC and GenomicSEM

The primary genetic covariance matrices were estimated with
`GenomicSEM::munge()` and `GenomicSEM::ldsc()`. The analysis used GenomicSEM
0.0.5 in R 4.4.3, with LDSC as the underlying method.

## Run the primary analysis

East Asian ancestry:

```bash
Rscript run_genomicsem_ldsc.R \
  --trait-config config/EAS_primary_traits.tsv \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --hm3 /path/to/eas_ldscores/w_hm3.snplist \
  --ld /path/to/eas_ldscores \
  --out-dir /path/to/output/EAS \
  --out-prefix EAS_MDD_EA_SCZ_CF
```

For EUR, use `config/EUR_primary_traits.tsv` with the EUR summary statistics
and LD scores. File paths in the trait configuration are resolved relative to
`--sumstats-dir`. `N=NA` uses the per-variant sample-size field.

`--validate-only` checks paths and input fields. `--skip-munge` reuses existing
munged files.

## Outputs

The script writes the covariance and correlation matrices returned by
GenomicSEM and the corresponding nearest-positive-definite matrices used for
plotting. The filenames end in `_raw.csv` and `_nearPD.csv`, respectively.
`*_matrix_adjustment_QC.tsv` records the size of the matrix adjustment.

Benjamini-Hochberg correction was applied across the 12 ancestry-specific
cross-trait genetic correlations: six EUR and six EAS tests.

## EAS sensitivity analysis

```bash
bash run_eas_clean_ldsc.sh \
  /path/to/ldsc \
  /path/to/eas_ldscores \
  /path/to/SCZ_EAS.sumstats.gz \
  /path/to/MDD_EAS.sumstats.gz \
  /path/to/EA_EAS.sumstats.gz \
  /path/to/CF_EAS.sumstats.gz \
  /path/to/output/EAS_clean
```

This pairwise sensitivity analysis used the standalone Python-compatible LDSC
implementation. It is separate from the primary GenomicSEM matrix analysis.

## Software

- R 4.4.3
- GenomicSEM 0.0.5
- Matrix 1.7-3
- ggplot2 4.0.3

GWAS summary statistics, HapMap3 lists and LD-score reference files are not
redistributed.
