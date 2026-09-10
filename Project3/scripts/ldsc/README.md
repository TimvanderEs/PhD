# LDSC workflows

The primary ancestry-specific matrices were estimated by
`GenomicSEM::munge()` and `GenomicSEM::ldsc()`, not by a direct invocation of
the standalone LDSC command-line program. The portable R workflow preserves
the executed trait order, disease prevalences, and standardisation settings.
It records GenomicSEM 0.0.5 and R 4.4.3 in the EC2 environment.

Calling this workflow only “LDSC v1.0.1” is therefore incomplete. Manuscript
Methods and software reporting should identify GenomicSEM and its version,
while citing LDSC as the underlying method.

## Raw and adjusted matrices

`run_genomicsem_ldsc.R` writes the inferential covariance estimate returned by
GenomicSEM and the `Matrix::nearPD()` transformation separately:

- `*_genetic_covariance_raw.csv` and `*_genetic_correlations_raw.csv` contain
  values derived directly from `ldsc_result$S`;
- `*_genetic_covariance_nearPD.csv` and
  `*_genetic_correlations_nearPD.csv` contain the transformed matrix used for
  the reconstructed heatmap/correlation export; and
- `*_matrix_adjustment_QC.tsv` quantifies the transformation.

The older archived correlation files in
[`final_output_fingerprints.tsv`](../../provenance/final_output_fingerprints.tsv)
were generated from the nearPD matrix but had ambiguous filenames. Their
fingerprints remain unchanged as historical provenance. A rerun uses the
explicit filenames above.

## Primary matrix runs

Example EAS command:

```bash
Rscript run_genomicsem_ldsc.R \
  --trait-config config/EAS_primary_traits.tsv \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --hm3 /path/to/eas_ldscores/w_hm3.snplist \
  --ld /path/to/eas_ldscores \
  --out-dir /path/to/output/EAS \
  --out-prefix EAS_MDD_EA_SCZ_CF
```

Use `config/EUR_primary_traits.tsv` with the EUR summary statistics and LD
scores for EUR. `config/EAS_HEL12k_sensitivity_traits.tsv` records the
five-trait HELIOS-12k QC/sensitivity iteration; it is not the primary matrix.

The `file` values are resolved relative to `--sumstats-dir`. `N=NA` tells
GenomicSEM to use the per-variant `N` field. Use `--validate-only` for a
path/schema check or `--skip-munge` to reuse existing munged files.

## Multiple-testing family

The publication workbook applies Benjamini-Hochberg correction once across all
12 cross-trait genetic-correlation tests (six EUR plus six EAS). It does not
apply two separate six-test corrections. Text, captions, and future scripts
should use the wording “BH correction across the 12 ancestry-specific genetic
correlations.”

## Clean EAS pairwise sensitivity rerun

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

This sensitivity run remunged SCZ, used the public-only MDD file, and used
sample/population prevalences `0.4503/0.01` for SCZ and `0.1278/0.15` for MDD.
It ran the six pairwise correlations among MDD, SCZ, educational attainment
(EA), and HELIOS-10k cognitive function (CF). The current publication workbook
traces to the GenomicSEM matrix, so this pairwise result must remain labelled a
sensitivity analysis unless the manuscript results are deliberately replaced.

## Software and external resources

- R 4.4.3
- GenomicSEM 0.0.5
- Matrix 1.7-3
- ggplot2 4.0.3

The EAS sensitivity shell workflow used the Python-2-compatible
Bulik-Sullivan LDSC implementation installed on EC2; its exact repository
revision was not retained. Reference LD scores and HapMap3 lists are external
resources and are not redistributed.
