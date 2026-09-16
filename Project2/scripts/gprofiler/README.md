# gProfiler

The 12 g:Profiler exports used in the paper are stored as compressed CSV files
under `results/gprofiler/raw/`. Their query files, settings and checksums are
listed in `inputs/gprofiler/run_manifest.tsv`.

Queries used the human organism, g:SCS correction and an adjusted-P threshold
of 0.05. EUR directional gene lists used all genes mapped in the corresponding
PLEIO model as the custom background. EAS analyses used the complete mapped
gene list for each model. The original g:Profiler database release was not
recorded, so the archived exports are used to reproduce the reported results.

## Rebuild the tables

Run from `Project2`:

```bash
python3 scripts/gprofiler/consolidate_gprofiler_results.py \
  --project-root . \
  --outdir results/gprofiler/processed
```

The script checks the input and export hashes, then writes:

- `significant_terms.tsv`: all terms with adjusted P < 0.05
- `run_summary.tsv`: query-level counts
- `source_summary.tsv`: counts by query and annotation source
- `manifest_validation.tsv`: file and row-count checks

The outputs support Supplementary Tables S13-S14.

## Rebuild Supplementary Figure S7

```bash
Rscript scripts/gprofiler/plot_gprofiler_enrichment.R \
  results/gprofiler/processed/significant_terms.tsv \
  work/gprofiler_enrichment \
  8 false
```

This command requires `data.table`, `ggplot2` and `patchwork`.
