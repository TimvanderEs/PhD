# LAVA local genetic-correlation workflow

`run_lava.R` is a path-portable version of the four-trait EAS and EUR scripts
used on EC2. It applies a univariate P-value screen of `2e-5` and then runs
`run.univ.bivar()` for every locus that passes LAVA processing.

## EAS

```bash
Rscript run_lava.R \
  --input-info config/EAS_input.info.tsv \
  --sample-overlap config/EAS_sample_overlap.tsv \
  --locus-file config/EAS_LAVA.locfile \
  --ref-prefix /path/to/g1000_eas \
  --phenos MDD_EAS,Edu_EAS,SCZ_EAS,HEL10k_EAS \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --out-dir /path/to/results_EAS \
  --out-prefix EAS_MDD_Edu_SCZ_HELIOS
```

## EUR

```bash
Rscript run_lava.R \
  --input-info config/EUR_input.info.tsv \
  --sample-overlap config/EUR_sample_overlap.tsv \
  --locus-file config/blocks_s2500_m25_f1_w200.GRCh37_hg19.LOCFILE \
  --ref-prefix /path/to/g1000_eur \
  --phenos MDD_EUR,Edu_EUR,SCZ_EUR,G_EUR \
  --sumstats-dir /path/to/EUR/summary_statistics \
  --out-dir /path/to/results_EUR \
  --out-prefix EUR_MDD_Edu_SCZ_G
```

The input-info filenames are resolved relative to `--sumstats-dir`. LAVA input
files must follow the format required by `process.input()`. The PLINK reference
prefix must resolve to `.bed`, `.bim`, and `.fam` files.

The overlap/intercept matrices and locus definitions are the exact,
non-individual-level inputs retained with the successful runs. The older EAS
files without HELIOS-10k pairings are excluded because they are QC iterations,
not the manuscript result.

## EC2 software environment

- R 4.2.2
- LAVA 0.1.5
- LAVA Git commit `4738b097bf929ec8af40225c196c57d46d3d8a22`
- data.table 1.17.8

The upstream software is available from
[josefin-werme/LAVA](https://github.com/josefin-werme/LAVA). Summary statistics
and PLINK reference panels are not redistributed.
