# FUMA MAGMA and GTEx enrichment

The archived FUMA results are under
`results/fuma_magma/raw/<ANCESTRY>/<MODEL>/`. Each run includes the FUMA
settings, MAGMA competitive gene-set results and GTEx v8 tissue results used
for the supplementary tables.

The seven single-trait runs used FUMA 1.5.2. The six PLEIO runs used FUMA
1.8.2. All runs used MAGMA 1.08, GRCh37 coordinates, ancestry-matched 1000
Genomes Phase 3 LD, Ensembl v102 protein-coding genes, MHC exclusion and a 10 kb
positional-mapping window. The settings for each run are stored in
`params.config`.

GWAS summary statistics and third-party gene-set membership files are not
included.

## Rebuild the tables

Install Python, NumPy and pandas, then run from `Project2`:

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

Multiple-testing correction is applied within each source file. Significant
enrichments require a positive MAGMA coefficient and a corrected P value below
0.05.

The retained tables are under `results/fuma_magma/processed/`. They support
Supplementary Tables S10-S12. The EAS cognitive-function single-trait result
was unavailable, so comparisons that require it are marked incomplete.
