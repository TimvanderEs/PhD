# Trans-ancestry meta-analysis

The Chinese, Malay and Indian REGENIE results were combined with METAL using a
fixed-effect inverse-variance model (`SCHEME STDERR`). The study order is kept
as Malay, Chinese, Indian to reproduce the `Direction` field in the retained
outputs.

Run the meta-analysis with:

```bash
bash run_metal_10k_trans_ancestry.sh \
  /path/to/metal \
  Malay_g.regenie \
  Chinese_g.regenie \
  Indian_g.regenie \
  OFF \
  results/HELIOS_10k_trans_ancestry_g
```

The fifth argument records whether METAL genomic control is `ON` or `OFF`.
Both settings were evaluated during the project; the revised analysis used
`OFF` after the ancestry-specific results had been checked for calibration.

Apply the retained variant filters and append chromosome and position parsed
from the marker identifier with:

```bash
python3 filter_metal_10k_trans_ancestry.py \
  results/HELIOS_10k_trans_ancestry_g1.sumstats.txt \
  results/HELIOS_10k_trans_ancestry_g.filtered.txt.gz
```

The filter retains variants represented in at least two ancestry groups
(`HetDf >= 1`), with heterogeneity P > 0.05 and maximum allele-frequency
difference < 0.5.
