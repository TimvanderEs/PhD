# Project 3: HELIOS cognitive GWAS

> **WCPG 2026:** [poster results, methods, QC and references](WCPG_2026.md)

This repository documents cognitive GWAS in Chinese, Indian and Malay HELIOS
participants and their trans-ancestry meta-analysis. It preserves the
historical 10k analysis, the expanded 22k analysis presented at WCPG 2026 and
a later refined Chinese-only freeze used in downstream work.

## September 2026 update

This release clarifies the phenotype-to-GWAS lineage and separates the retained
22k analysis from the historical 10k, nested 12k and refined Chinese-only
versions. It adds the post-PCA ancestry-specific phenotype preparation,
ancestry-versus-recorded-ethnicity sensitivity analysis, HELIOS–KoGES local
covariance audit and six-task Genomic SEM validation. The accompanying code,
aggregate results and provenance records contain no participant-level data.

## WCPG 2026 analysis

| Metric | Result |
|---|---:|
| Cognitive PCA sample | 22,422 |
| PC1 variance explained | 37.04% |
| Chinese GWAS | 17,966 |
| Indian GWAS | 2,310 |
| Malay GWAS | 1,242 |
| Total GWAS sample | 21,518 |
| Retained METAL variants | 6,518,064 |
| Genome-wide significant variants | 0 |
| Minimum retained P | 6.906 × 10⁻⁷ |

The PCA loadings shown on the poster use the orientation in which higher scores
indicate better performance. The archived score has the opposite, sign-equivalent
orientation.

Before GWAS, the retained PC1 was split by genetic-ancestry analysis group,
values beyond five within-group standard deviations were set to missing and
the remaining scores were Blom inverse-normal transformed within group. The
REGENIE Step 2 commands therefore did not apply a second RINT.

The trans-ancestry analysis is a fixed-effect inverse-variance METAL analysis
of the three ancestry-specific GWASs. It targets shared effects and does not
model ancestry-heterogeneous effects explicitly. See the
[expanded 22k summary](validation/22k/22k_summary.md).

## Analysis versions

| Analysis | PCA sample | Final phenotype preparation | GWAS sample | Role |
|---|---:|---|---:|---|
| Expanded 22k | 22,422 | Within-ancestry 5-SD mask and Blom RINT | 21,518 | Primary WCPG 2026 analysis |
| Historical 10k PCA reconstruction | 7,403 | Within-ancestry 5-SD mask and Blom RINT | 7,403 | Historical factor/Step 2 lineage |
| 12k extension | 13,774 | Within-ancestry 5-SD mask and Blom RINT | 13,751 | Nested diagnostic analysis |
| Refined Chinese-only 10k | — | REGENIE RINT in Step 2 only | 5,700 | Separate downstream freeze |

## Analysis records

- [Expanded 22k phenotype validation and GWAS](validation/22k/README.md), the
  analysis presented at WCPG 2026
- [Historical 10k phenotype and trans-ancestry workflow](WORKFLOW.md)
- [Refined Chinese-only summary-statistic freeze](results/summary_statistics/freeze_manifest.tsv)
- [Cross-freeze genetic-signal audit](validation/22k/03_genetic_signal/genetic_signal_audit.md)
- [Phenotype-to-GWAS lineage and RINT sensitivity](validation/phenotype_pipeline/README.md)
- [Six-task Genomic SEM validation](validation/genomic_sem_components/FINAL_REPORT.md)

## Code and documentation

- [Phenotype construction](scripts/phenotype)
- [Ancestry-specific genotype QC](scripts/qc)
- [Historical 10k GWAS](scripts/gwas)
- [Expanded 22k GWAS](scripts/gwas_22k)
- [METAL meta-analysis](scripts/meta_analysis)
- [Summary-statistic preparation](scripts/summary_statistics)
- [SG100K covariate-adjusted LD reference](scripts/ld_reference)
- [Software and cohort references](REFERENCES.md)
- [Full workflow](WORKFLOW.md)

## Data availability

Individual-level HELIOS and SG100K data cannot be deposited in this
repository. Access is subject to the relevant cohort consent, data-access and
institutional governance procedures. The public repository contains code,
small aggregate results and provenance only. Summary-statistic access and
release are handled separately under the relevant study approvals.

## Citation

The publication release will be archived with a DOI. Until then, cite the
repository commit used for the analysis. Code is released under the
[MIT License](../LICENSE).
