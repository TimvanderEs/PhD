# Project 4: Singapore LD reference panels

This project contains the LD-reference work from **Chapter 3, Study 2** of my
PhD: construction and evaluation of linkage disequilibrium (LD) reference panels
from a PRECISE-SG100K-derived Singapore genomic resource. It asks whether ancestry
matching, reference size, or adjustment for population structure changes
LD-score regression (LDSC) estimates of SNP heritability and genetic correlation.

The work includes Chinese, Indian and Malay ancestry-specific panels, pooled
panels, and a controlled comparison with 1000 Genomes East Asian data. It
complements the [HELIOS cognitive GWAS in Project 3](../Project3/).

## Where to start

- **Research overview:** the main findings and study design below.
- **Chapter 3 material:** the [chapter-to-code guide](CHAPTER3_GUIDE.md) connects
  the study to construction scripts, aggregate results and validation records.
- **Try the code:** the [portable workflow](portable/README.md) includes a
  synthetic example that requires no controlled participant data or LDSC installation.
- **Reproduce the research analysis:** the [workflow](WORKFLOW.md) documents
  required controlled inputs, software and execution order.

The research outputs and portable prototype have different scopes: the archived
results come from the controlled analysis, while the synthetic example exercises
input validation, command planning and paired-comparison arithmetic.

## Main finding

Under a controlled 1 cM specification, the standard Singapore Chinese reference
and an exact 1000 Genomes East Asian reference produced very similar results for
six Chinese SG100K GWAS: systolic blood pressure, waist-to-hip ratio, height,
body mass index, LDL cholesterol and type 2 diabetes.

| Reference | N | Mean h2 | Mean SE | Mean h2/SE |
| --- | ---: | ---: | ---: | ---: |
| 1000 Genomes EAS | 504 | 0.1601 | 0.02022 | 7.23 |
| Singapore Chinese | 27,178 | 0.1618 | 0.02035 | 7.29 |

None of the six paired h2 differences was nominally significant, and none survived Benjamini–Hochberg correction. Across 15 trait pairs, genetic-correlation estimates correlated at r = 0.99945, had a mean absolute difference of 0.00736, and agreed in direction for all pairs. Non-significance does not establish statistical equivalence.

The values above are documented in the [panel summary](results/controlled_reference/controlled_h2_panel_summary.tsv),
[paired comparisons](results/controlled_reference/paired_h2_differences.tsv) and
[genetic-correlation concordance table](results/controlled_reference/controlled_rg_concordance.tsv).

The much larger local reference therefore improved LD-score stability but did not measurably improve h2 precision or calibration for these GWAS. A five-replicate held-out simulation likewise did not distinguish the panels and is retained as a pipeline calibration check rather than definitive evidence of equivalence.

## Controlled experiment

The final experiment held the following constant across five panels: GRCh38 map, 1 cM LD window, 874,025-variant universe and order, Chinese GWAS inputs, 1000 Genomes EAS regression weights, and 200 jackknife blocks.

The five panels were:

1. exact 1000 Genomes EAS calculated from 504 raw-genotype samples;
2. Singapore Chinese standard LD scores, N = 27,178 unrelated participants;
3. Singapore Chinese LD scores adjusted for 20 PCs;
4. pooled Singapore standard LD scores, N = 37,351 unrelated participants; and
5. pooled Singapore LD scores adjusted for 20 PCs.

PC adjustment had little effect within the Chinese reference but materially changed the pooled reference. This makes covariate specification an important design choice for pooled multi-ancestry panels; the larger estimate is not, by itself, evidence of greater accuracy.

Complete results and run records are in [`validation/controlled_reference_experiment/`](validation/controlled_reference_experiment/). The concise interpretation is in its [`FINAL_REPORT.md`](validation/controlled_reference_experiment/FINAL_REPORT.md).

## Repository contents

- [`portable/`](portable/): Python 3 configuration-driven validation, command
  planning, panel construction from prepared genotypes, LDSC execution and
  paired-comparison prototype, with a synthetic example and tests. No controlled
  inputs are required for the example.

- [`scripts/`](scripts/): panel construction, exact 1KG recomputation, LDSC, reference-size, simulation, and summary scripts.
- [`WORKFLOW.md`](WORKFLOW.md): required inputs, execution order, and paired-jackknife definition.
- [`results/controlled_reference/`](results/controlled_reference/): reader-facing final aggregate tables.
- [`validation/controlled_reference_experiment/`](validation/controlled_reference_experiment/): final controlled outputs, logs, QC records, and software information.
- [`results/ancestry_matched/`](results/ancestry_matched/): earlier Chinese, Indian, and Malay ancestry-matched validation.
- [`results/recovered/`](results/recovered/): recovered historical Singapore analyses.
- [`PROJECT4_AUDIT.md`](PROJECT4_AUDIT.md): provenance audit and interpretation boundaries.
- [`CHAPTER3_GUIDE.md`](CHAPTER3_GUIDE.md): guide to the Chapter 3 study and its repository materials.
- [`CHAPTER_DRAFT.md`](CHAPTER_DRAFT.md): extended working text for the LD-reference study.

## Earlier analyses

The repository retains two earlier analysis stages for provenance:

- a pooled 20 cM cov-LDSC analysis; and
- Chinese, Indian, and Malay 1 cM panels applied to trans-ancestry SG100K GWAS, followed by ancestry-matched GWAS analyses.

These analyses used different samples, variant sets, or specifications and should not be combined as a single performance benchmark. The final controlled experiment supersedes the earlier translated-row 1KG comparison for claims about Chinese panel performance.

## Scope

This repository addresses post-GWAS LDSC. REGENIE uses individual-level study genotypes to model polygenic effects and does not use an external LDSC reference panel. The Singapore panels are therefore not replacement inputs for the original REGENIE association analyses.

Individual-level genotypes, sample identifiers, covariates, and provider-controlled GWAS files are not distributed here.
