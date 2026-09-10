# FUMA SNP2GENE, MAGMA, and GTEx v8 enrichment

This directory contains the Python scripts used to audit and consolidate the
model-wide PLEIO and constituent single-trait FUMA/MAGMA results.

## What is archived

`results/fuma_magma/raw/<ANCESTRY>/<MODEL>/` contains, for each selected run:

- `params.config`: the complete FUMA web-run settings and software versions;
- `magma.gsa.out`: raw MAGMA competitive gene-set results;
- `magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out`: raw GTEx v8 tissue-specific
  gene-property results;
- `magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out`: raw GTEx v8 general-tissue
  gene-property results.

These are the exact EC2 files used to build the publication tables. The source
GWAS summary statistics are not redistributed. The large
`magma.gsa.sets.genes.out` database-membership files are also omitted because
they were not used by the reported term-level comparisons and may reproduce
third-party gene-set database content.

The seven available single-trait runs used FUMA v1.5.2. The six final PLEIO
runs used FUMA v1.8.2. All used MAGMA v1.08, 1000 Genomes Phase 3 ancestry-
matched LD, GRCh37 coordinates, protein-coding Ensembl v102 genes, MHC
exclusion, a 10 kb positional-mapping window, no eQTL or chromatin-interaction
mapping, and GTEx v8 expression panels. Exact per-run differences, including
input field mappings and sample-size fields, remain in each `params.config`.

## Rebuild the extracted tables

Create an environment with Python 3, NumPy, and pandas, then run these commands
from the `Project3` directory:

```bash
python3 scripts/fuma_magma/extract_fuma_magma_gtex.py \
  --run-manifest inputs/fuma_magma/single_trait_run_manifest.tsv \
  --scope single-trait \
  --outdir work/fuma_single_trait

python3 scripts/fuma_magma/extract_fuma_magma_gtex.py \
  --run-manifest inputs/fuma_magma/pleio_run_manifest.tsv \
  --scope multivariate \
  --outdir work/fuma_pleio

python3 scripts/fuma_magma/compare_magma_single_vs_pleio.py \
  --single-dir work/fuma_single_trait \
  --pleio-dir work/fuma_pleio \
  --outdir work/fuma_single_vs_pleio
```

The extractor applies Bonferroni and Benjamini-Hochberg corrections within each
source file. Reported competitive gene-set and GTEx enrichments require both a
corrected-significant result and a positive MAGMA beta.

Compact retained outputs are in `results/fuma_magma/processed/`. Absolute EC2
paths in their provenance columns identify the original files; a local rerun
will naturally record local paths instead.

## Interpretation limit

No EAS cognitive-performance single-trait SNP2GENE run was available. The EAS
four-trait and nested three-trait comparisons are therefore explicitly marked
incomplete in
`results/fuma_magma/processed/comparison/11_missing_or_incomplete_comparisons.tsv`.
`PLEIO_only_among_available_constituents` is a descriptive significance pattern,
not evidence that an enrichment is biologically unique or causal.
