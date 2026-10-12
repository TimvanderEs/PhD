# LD-reference analysis scripts

These scripts construct and compare Singapore LD-score reference panels for
Chapter 3, Study 2. Start with the [research workflow](../WORKFLOW.md) for input
requirements and execution order, or the [portable example](../portable/README.md)
to try validation and command planning without controlled data.

## Primary controlled comparison

| Step | Entry points | Purpose |
| --- | --- | --- |
| Prepare reference inputs | `prepare_controlled_reference_inputs.sh`, `prepare_covldsc_eigenvectors.py`, `interpolate_bim_genetic_map.py` | Prepare unrelated samples, PC covariates and genetic-map positions |
| Construct Singapore panels | `build_controlled_ldscore_panels.sh` | Build Chinese and pooled panels with and without PC adjustment |
| Recalculate the comparator | `build_exact_1kgeas_panel.sh`, `prepare_exact_1kgeas_mappings.py` | Build the exact 1000 Genomes EAS comparator |
| Match reference panels | `match_ldscore_panels.py`, `align_ldscore_panel_order.py` | Match variant membership and order |
| Estimate heritability and genetic correlation | `run_controlled_panel_ldsc.sh` | Apply five reference panels to the same Chinese GWAS inputs |
| Summarise comparisons | `summarise_paired_h2_jackknife.py`, `summarise_controlled_rg.py` | Produce paired heritability contrasts and cross-panel genetic-correlation summaries |
| Check exported results | `validate_controlled_experiment.py` | Check table sizes, run completion, the common variant count and identical reader-facing copies |

## Sensitivity analyses

- `build_reference_size_panels.sh` and `run_reference_size_h2.sh` evaluate
  smaller Chinese reference samples. The `summarise_reference_size*.py`
  scripts summarise heritability and LD-score stability.
- `run_heldout_simulation.sh` and `summarise_simulation.py` implement and
  summarise the five-replicate calibration experiment.
- `prepare_sg100k_ancestry_sumstats.sh`, `run_sg100k_ancestry_matched_ldsc.sh`
  and `summarise_sg100k_ancestry_matched.py` support the earlier
  [ancestry-matched analyses](../validation/ancestry_matched_sg100k/).

Earlier recovery and translated-1KG scripts are retained to explain archived
outputs. The primary comparison uses the exact 1KG reconstruction, as described
in the [controlled-analysis report](../validation/controlled_reference_experiment/FINAL_REPORT.md).

## Execution and data access

The shell scripts retain default paths from the research environment. Configure
the environment variables documented at the top of each script before running
it. The [workflow](../WORKFLOW.md) lists the required PLINK, Python and
LDSC/cov-LDSC environments and controlled input files.

You can validate the exported aggregate results locally from `Project4/`:

```sh
python3 scripts/validate_controlled_experiment.py \
  validation/controlled_reference_experiment \
  --reader-results results/controlled_reference
```

This checks the saved result bundle; it does not rerun LDSC or verify the
underlying controlled genotypes. Participant-level inputs remain subject to
their original access conditions.
