# HELIOS genetic ancestry and recorded ethnicity

This analysis distinguishes three variables that were previously described interchangeably:

- `grids_meta.ancestry`: GRID genetic ancestry assignment;
- `FREG5_Race`: race recorded from the Singapore NRIC; and
- `FIAQ10_Demo1`: participant-described ethnic heritage.

The main findings are in [`HELIOS_ancestry_selfreport_analysis.md`](HELIOS_ancestry_selfreport_analysis.md), with draft manuscript wording in [`THESIS_TEXT.md`](THESIS_TEXT.md). The central result is that the ancestry-specific 10k and 22k GWAS inputs use the core GRID genetic groups, not the NRIC or questionnaire variables. Self-described ethnicity and GRID ancestry are strongly related but are not interchangeable.

The follow-up [`phenotype_sensitivity`](phenotype_sensitivity) analysis tests the effect of correcting the original pooled phenotype ancestry coding, restricting construction to the three core GRID groups, adding participant-described ethnic heritage and excluding discordant records. All versions produced essentially unchanged PC1 loadings and scores.

## Contents

- `figures/`: PCA and concordance figures in PNG and PDF format;
- `tables/`: aggregate counts, concordance estimates, statistical tests and lineage checks;
- `scripts/analyse_ancestry_selfreport.py`: complete analysis and figure-generation script;
- `scripts/run_phenotype_ethnicity_sensitivity.R`: EC2 phenotype sensitivity script;
- `phenotype_sensitivity/`: aggregate results from the phenotype sensitivity analysis;
- `validation_checks.tsv`: automated integrity checks; and
- `software_and_settings.tsv`: software versions and disclosure-control setting.

No participant identifiers, individual ancestry records or individual PCA coordinates are written to this directory. Detailed cross-tabulations apply primary and row-complementary suppression with a minimum cell size of five.

## Reproduction

Run with authorised local copies of the controlled SG100K/HELIOS files:

```bash
python3 scripts/analyse_ancestry_selfreport.py \
  --annotation /path/to/Sample_Annotation_r231.csv \
  --core /path/to/SG100K_Core_r1_v1_1.csv \
  --interview /path/to/SG100K_IAQ_anon_r1_v1_1.csv \
  --data-dictionary "/path/to/SG100K Data Dictionary - v46_r1_v1_0_3.xlsx" \
  --full22 /path/to/Full_data_HELIOS_22k.txt \
  --covariates10 /path/to/merged_pheno_cov_filtered_10k.txt \
  --sg100k-covariates /path/to/sg100k_r231_covariate_file.txt \
  --gwas22-chinese /path/to/CHI_gwas.cov \
  --gwas22-indian /path/to/IND_gwas.cov \
  --gwas22-malay /path/to/MAL_gwas.cov \
  --output-root .
```

The script requires Python, NumPy, pandas, Pillow and openpyxl. It fails if participant identifiers are duplicated in the genetic annotation, GRID PC1–PC5 are incomplete, or the independent NRIC variables disagree after duplicate resolution.
