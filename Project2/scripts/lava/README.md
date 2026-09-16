# LAVA

`run_lava.R` runs ancestry-specific local heritability and bivariate local
genetic-correlation analyses. Univariate significance is defined as 0.05
divided by the number of regions. `postprocess_lava.py` selects eligible
bivariate tests and applies Benjamini-Hochberg correction separately to each
ancestry and trait pair.

## Requirements

- R 4.2.2
- LAVA 0.1.5 at commit `4738b097bf929ec8af40225c196c57d46d3d8a22`
- data.table 1.17.8
- ancestry-matched PLINK reference panels
- GWAS summary statistics formatted for `LAVA::process.input()`

LAVA is available from [josefin-werme/LAVA](https://github.com/josefin-werme/LAVA).
Summary statistics and reference panels are not included in this repository.

## Run the analysis

East Asian ancestry:

```bash
Rscript run_lava.R \
  --input-info config/EAS_input.info.tsv \
  --sample-overlap config/EAS_sample_overlap.tsv \
  --locus-file config/EAS_LAVA.locfile \
  --ref-prefix /path/to/g1000_eas \
  --phenos MDD,EA,SCZ,CF \
  --sumstats-dir /path/to/EAS/summary_statistics \
  --out-dir /path/to/results_EAS \
  --out-prefix EAS_MDD_EA_SCZ_CF
```

European ancestry:

```bash
Rscript run_lava.R \
  --input-info config/EUR_input.info.tsv \
  --sample-overlap config/EUR_sample_overlap.tsv \
  --locus-file config/blocks_s2500_m25_f1_w200.GRCh37_hg19.LOCFILE \
  --ref-prefix /path/to/g1000_eur \
  --phenos MDD,EA,SCZ,CF \
  --sumstats-dir /path/to/EUR/summary_statistics \
  --out-dir /path/to/results_EUR \
  --out-prefix EUR_MDD_EA_SCZ_CF
```

The PLINK reference prefix must resolve to `.bed`, `.bim` and `.fam` files.
Input-info filenames are resolved relative to `--sumstats-dir`.

## Apply FDR correction

```bash
python3 postprocess_lava.py \
  --univ /path/to/results_EAS/EAS_MDD_EA_SCZ_CF.univ.lava \
  --bivar /path/to/results_EAS/EAS_MDD_EA_SCZ_CF.bivar.lava \
  --ancestry EAS \
  --n-regions 3064 \
  --trait-aliases config/EAS_trait_aliases.tsv \
  --out-dir /path/to/results_EAS/consolidated
```

For EUR, use `--n-regions 2495` and `config/EUR_trait_aliases.tsv`.

The final eligible set contained 698 EAS tests and 3,142 EUR tests. Counts by
trait pair are provided in
[`provenance/lava_pair_test_counts.tsv`](../../provenance/lava_pair_test_counts.tsv).

The post-processing script writes the eligible bivariate results, a summary by
trait pair and a QC file containing thresholds and row counts. The synthetic
test checks the eligibility threshold, zero P values, trait order and
within-pair FDR correction:

```bash
python3 ../../tests/test_lava_postprocessing.py
```
