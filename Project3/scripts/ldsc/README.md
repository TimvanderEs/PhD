# LDSC workflows

The primary ancestry-specific covariance matrices were estimated with
`GenomicSEM::munge()` and `GenomicSEM::ldsc()`. The portable R script retains
the trait order, disease prevalences, standardisation, and near-positive-
definite correlation conversion used on EC2.

The manuscript-facing EAS estimates were later checked with a clean pairwise
LDSC rerun. That rerun remunged SCZ, used the public-only munged MDD file, and
used the actual case fractions for the disease-trait liability conversions.
It is therefore provided separately in `run_eas_clean_ldsc.sh`.

## Primary matrix runs

Example EAS command:

```bash
Rscript run_genomicsem_ldsc.R \
  --trait-config config/EAS_primary_traits.tsv \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --hm3 /path/to/eas_ldscores/w_hm3.snplist \
  --ld /path/to/eas_ldscores \
  --out-dir /path/to/output/EAS \
  --out-prefix EAS_MDD_EDU_SCZ_HEL10k
```

Use `config/EUR_primary_traits.tsv` with the EUR summary-statistics and LD-score
paths for the EUR run. `config/EAS_HEL12k_sensitivity_traits.tsv` records the
five-trait HELIOS-12k sensitivity/QC iteration; it was not the primary
manuscript matrix.

The `file` values in each config are resolved relative to `--sumstats-dir`.
`N=NA` tells GenomicSEM to use the per-variant `N` column in the file. Use
`--validate-only` for a path/schema check or `--skip-munge` to reuse existing
munged files.

## Clean EAS pairwise rerun

```bash
bash run_eas_clean_ldsc.sh \
  /path/to/ldsc \
  /path/to/eas_ldscores \
  /path/to/SCZ_EAS_simple_for_munge.sumstats.gz \
  /path/to/MDD_EAS_public_noWHI_no23andMe.sumstats.gz \
  /path/to/Edu_EAS.sumstats.gz \
  /path/to/HEL10k_EAS.sumstats.gz \
  /path/to/output/EAS_clean
```

The exact rerun settings were sample/population prevalence `0.4503/0.01` for
SCZ and `0.1278/0.15` for MDD. All six pairwise correlations among MDD, SCZ,
educational attainment, and HELIOS-10k cognitive performance are run.

## EC2 software environment

- R 4.4.3
- GenomicSEM 0.0.5
- Matrix 1.7-3
- ggplot2 4.0.3

The clean shell workflow used the Python-2-compatible Bulik-Sullivan LDSC
implementation installed on EC2. Reference LD scores and HapMap3 lists are
external resources and are not redistributed.
