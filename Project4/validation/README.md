# External-cohort input audit

This directory records an audit of the summary statistics used in the early external-cohort LDSC benchmark. That benchmark must not be used as evidence about LD-panel performance.

The archived download script assigned trait names to GWAS Catalog accessions that describe different phenotypes. It also applied fixed sample sizes when a file did not contain a recognised sample-size column. The resulting cohort/trait labels and LDSC estimates are invalid.

The audit has three components:

- [`cohort_accession_audit.tsv`](cohort_accession_audit.tsv): intended labels versus official GWAS Catalog traits;
- [`bbj_height_validation.md`](bbj_height_validation.md): a controlled validation using the genuine Biobank Japan height GWAS; and
- [`scripts/run_bbj_height_validation.sh`](scripts/run_bbj_height_validation.sh): the exact reusable munging and LDSC commands.

Small historical and corrected logs are retained under [`../results/validation/`](../results/validation/). The scripts that generated the invalid labels are preserved under [`provenance/source_archive/`](provenance/source_archive/) for audit only and must not be run as a data-download workflow.

The later six-trait SG100K comparison in `Project4/results/historical/` is a separate analysis. Its semantic phenotype identifiers, sample sizes and source release were verified independently; the accession error described here does not apply to those six files.

The design of that later analysis is reconstructed in [`analysis_design_reconstruction.md`](analysis_design_reconstruction.md), with a compact machine-readable summary in [`analysis_matrix.tsv`](analysis_matrix.tsv). The Chinese, Indian, and Malay result directories identify the LD reference ancestry; the GWAS inputs were the same fixed-effect trans-ancestry files in every directory.

The historical five-dataset benchmark remains unresolved apart from the corrected BioBank Japan height analysis. Dataset-level status is recorded in [`five_cohort_rebuild_status.tsv`](five_cohort_rebuild_status.tsv).

The new six-trait ancestry-matched SG100K validation is documented separately in [`ancestry_matched_sg100k/`](ancestry_matched_sg100k/). It does not rehabilitate the historical external-cohort benchmark.

The final reference-panel comparison is documented in [`controlled_reference_experiment/`](controlled_reference_experiment/). It recalculates 1KG EAS and Singapore LD scores on one 874,025-variant universe and contains the primary h2, rg, reference-size, and held-out simulation results.
