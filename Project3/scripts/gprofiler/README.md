# g:Profiler over-representation analysis

This directory makes the publication g:Profiler analysis auditable without
silently replacing it with a query against a newer database release.

## Archived analysis

The 12 original CSV exports are stored as deterministic gzip files in
`results/gprofiler/raw/`. Their uncompressed and compressed SHA-256 hashes,
input files, settings, row counts, and local download timestamps are recorded
in `inputs/gprofiler/run_manifest.tsv`.

The analysis used `hsapiens`, g:SCS multiple-testing correction, a significance
threshold of adjusted P < 0.05, and the following sources: GO biological
process, molecular function and cellular component; Reactome; KEGG;
WikiPathways; Human Phenotype Ontology; Human Protein Atlas; CORUM; TRANSFAC;
and miRTarBase.

For EUR, each directional gene list was tested using all mapped genes from the
corresponding PLEIO model as a custom background. The publication label `Dual`
corresponds to the archived input lists whose source filenames contain
`mixed_dual`; this is an archived source filename only, and its publication
class is `Dual`. For EAS, the complete model-wide mapped-gene list was queried.
No EAS custom-background file or original web request payload was recovered,
so the EAS background option is explicitly marked `not_recorded`.

The original web exports also do not identify the g:Profiler database release,
ordered-query option, or electronic-GO-annotation option. Because g:Profiler's
underlying databases change, a new live query may not return byte-identical
results. The archived exports are the primary record of the reported analysis.

## Validate and rebuild the processed tables

Run from the `Project3` directory:

```bash
python3 scripts/gprofiler/consolidate_gprofiler_results.py \
  --project-root . \
  --outdir results/gprofiler/processed
```

The script validates every input and raw-export checksum, verifies the expected
row and significant-term counts, and writes:

- `significant_terms.tsv`: all 61 terms with g:SCS-adjusted P < 0.05;
- `run_summary.tsv`: one row for each of the 12 queries;
- `source_summary.tsv`: counts and minimum adjusted P by query and source;
- `manifest_validation.tsv`: checksum, row-count, uniqueness, and background
  containment checks.

`run_summary.tsv` and `source_summary.tsv` supply the g:Profiler overview in
Supplementary Table S13. The 61 significant terms reproduce Supplementary
Table S14 by ancestry, model, direction class, source label, term identifier,
adjusted P value, and intersection size.

## Rebuild Supplementary Figure S7

With R packages `data.table`, `ggplot2`, and `patchwork` installed:

```bash
Rscript scripts/gprofiler/plot_gprofiler_enrichment.R \
  results/gprofiler/processed/significant_terms.tsv \
  work/gprofiler_enrichment \
  8 false
```

The publication default selects the eight lowest adjusted-P terms per query
without ties and retains the complete 61-term table alongside the figure.
