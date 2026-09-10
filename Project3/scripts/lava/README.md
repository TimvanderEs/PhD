# LAVA local genetic-correlation workflow

`run_lava.R` is a path-portable reconstruction of the four-trait EAS and EUR
workflows used on EC2. By default, it screens univariate results at exactly
`0.05 / number of regions` and then runs `run.univ.bivar()` at every processed
locus. `postprocess_lava.py` independently reapplies that eligibility rule and
performs Benjamini-Hochberg correction within each ancestry-by-trait-pair
family.

## Important audit status

The archived outputs are evidence of the executed historical analysis, not a
clean final rerun:

- the historical runner used the rounded threshold `2e-5`;
- the archived EAS bivariate file has 726 rows;
- applying the stated `P <= 0.05/3064` rule to the archived EAS univariate and
  bivariate files retains 698 rows;
- the current publication workbook has 697 EAS rows because one EA--SCZ row at
  locus 605, whose univariate P value is numerically zero, is absent; and
- the recovered input-info files contain shifted or incorrect disease
  case/control fields.

For those reasons, the repository does **not** relabel 726 or 697 as the
verified final denominator. The clean resolution is to verify that the exact
EAS SCZ file is the Lam et al. dataset, rerun both ancestries with the corrected
input-info files, and regenerate the manuscript tables with
`postprocess_lava.py`. The three count sets are preserved separately in
[`provenance/lava_pair_test_counts.tsv`](../../provenance/lava_pair_test_counts.tsv).

## Input configurations

The candidate `EAS_manuscript_target_input.info.tsv` and
`EUR_manuscript_target_input.info.tsv` files use the sample counts reported in
the current manuscript. They are clearly named as targets, not executed-run
records. Continuous traits use `cases=1` and `controls=0`, as required by LAVA
input processing.

| Ancestry | Trait | Active cases / controls | Audit note |
| --- | --- | --- | --- |
| EUR | MDD | 525,197 / 3,362,335 | Matches manuscript Table 1. |
| EUR | SCZ | 53,386 / 77,258 | Matches manuscript Table 1 and recovered QC log. |
| EAS | MDD | 15,771 / 178,777 | Matches manuscript Table 1. |
| EAS | SCZ | 22,778 / 35,362 | Matches manuscript Table 1 and the primary GWAS publication; the recovered SumstatsQC log instead reports 27,888 / 35,362, so verify the exact file identity/header before rerun. |

The EAS SCZ target is supported by the primary GWAS report,
[Lam et al., Nature Genetics 2019](https://doi.org/10.1038/s41588-019-0512-x),
which states 22,778 cases and 35,362 controls. The conflicting QC-log value is
retained because file-specific provenance should not be overwritten by the
publication-level count.

The exact configurations inferred from the archived run are retained as
`EAS_recovered_run_input.info.tsv` and `EUR_recovered_run_input.info.tsv`.
They must not be used for a clean rerun. They explain why the archived EAS
univariate output has a liability-scale estimate for EA but not MDD.

## Run EAS

```bash
Rscript run_lava.R \
  --input-info config/EAS_manuscript_target_input.info.tsv \
  --sample-overlap config/EAS_sample_overlap.tsv \
  --locus-file config/EAS_LAVA.locfile \
  --ref-prefix /path/to/g1000_eas \
  --phenos MDD,EA,SCZ,CF \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --out-dir /path/to/results_EAS \
  --out-prefix EAS_MDD_EA_SCZ_CF
```

## Run EUR

```bash
Rscript run_lava.R \
  --input-info config/EUR_manuscript_target_input.info.tsv \
  --sample-overlap config/EUR_sample_overlap.tsv \
  --locus-file config/blocks_s2500_m25_f1_w200.GRCh37_hg19.LOCFILE \
  --ref-prefix /path/to/g1000_eur \
  --phenos MDD,EA,SCZ,CF \
  --sumstats-dir /path/to/EUR/summary_statistics \
  --out-dir /path/to/results_EUR \
  --out-prefix EUR_MDD_EA_SCZ_CF
```

The input-info filenames are resolved relative to `--sumstats-dir`. Summary
statistics must follow the format required by `LAVA::process.input()`. The
PLINK reference prefix must resolve to `.bed`, `.bim`, and `.fam` files. Use
`--univ-threshold` only to document an intentional sensitivity analysis.

## Consolidate and apply FDR

For a new EAS run:

```bash
python3 postprocess_lava.py \
  --univ /path/to/results_EAS/EAS_MDD_EA_SCZ_CF.univ.lava \
  --bivar /path/to/results_EAS/EAS_MDD_EA_SCZ_CF.bivar.lava \
  --ancestry EAS \
  --n-regions 3064 \
  --trait-aliases config/EAS_trait_aliases.tsv \
  --out-dir /path/to/results_EAS/consolidated
```

Use `--n-regions 2495` and `config/EUR_trait_aliases.tsv` for EUR. The program
stops if an eligible bivariate row is missing, retains valid P values equal to
zero, records all thresholds and row counts, and writes:

- `lava_bivariate_eligible.tsv`, the eligible row-level result;
- `lava_pair_summary.tsv`, test and FDR-significant counts by pair; and
- `lava_postprocess_qc.tsv`, the machine-readable audit record.

The synthetic check in `tests/test_lava_postprocessing.py` verifies the exact
threshold, zero-P retention, reversed trait order, and within-pair BH logic.

## Retained inputs and software

The overlap/intercept matrices and locus definitions are the exact
non-individual-level inputs retained with the successful historical runs.
Older EAS outputs without HELIOS-10k pairings were QC iterations and are not
represented as final.

- R 4.2.2
- LAVA 0.1.5
- LAVA commit `4738b097bf929ec8af40225c196c57d46d3d8a22`
- data.table 1.17.8

The upstream software is available from
[josefin-werme/LAVA](https://github.com/josefin-werme/LAVA). Summary statistics
and PLINK reference panels are not redistributed.
