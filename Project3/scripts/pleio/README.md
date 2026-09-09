# PLEIO input preparation and model runs

These wrappers reproduce the four-trait and nested three-trait PLEIO models.
They do not redistribute PLEIO itself. The EC2 checkout reports that upstream
software as unlicensed, so obtain it from
[cuelee/pleio](https://github.com/cuelee/pleio) and check out the revision used:

```bash
git clone https://github.com/cuelee/pleio.git
cd pleio
git checkout b18d507533ae39dcb6284573d3b36b91cdc13603
```

This upstream revision has a split runtime: `ldsc_preprocess.py` uses Python 2
syntax, while `pleio.py` uses Python 3 syntax. Create both recorded environments:

```bash
conda env create -f environment_preprocess.yml
conda env create -f environment.yml
```

## 0. Align the EUR analysis-ready files

The recovered EC2 alignment workflow standardises the source columns, finds
the shared SNP set, uses educational attainment as the allele anchor, removes
A/T and C/G variants, and reverses Z scores when effect alleles are swapped:

```bash
bash align_eur_pleio_inputs.sh \
  /path/to/MDD_EUR_LDSC.sumstats.TAB.gz \
  /path/to/EA4_Edu.qc.txt_noINFO.Z.TAB.gz \
  /path/to/SCZ_EUR_for_LDSC.neff.sumstats.TAB.gz \
  /path/to/Cognition_MTAG.qc.txt_noINFO.Z.TAB.gz \
  /path/to/EUR/pleio_full_aligned_hdr
```

The script validates the exact source headers before applying the original
column mappings. In the retained EC2 run, the subsequent PLEIO preprocessing
reported 4,930,902 shared variants and no allele mismatches.

## 1. Prepare each four-trait model

The input GWAS files must contain `SNP`, `A1`, `A2`, `Z`, and `N` columns.
The EUR manifest names the outputs from step 0.

```bash
bash prepare_four_trait_inputs.sh \
  --pleio-dir /path/to/pleio \
  --trait-config config/EAS_four_trait_manifest.tsv \
  --sumstats-dir /path/to/EAS/analysis_ready_sumstats \
  --ref-ld-prefix /path/to/eas_ldscores/ \
  --out-dir /path/to/analysis/EAS_four_trait \
  --python /path/to/project3-pleio-preprocess/bin/python

bash prepare_four_trait_inputs.sh \
  --pleio-dir /path/to/pleio \
  --trait-config config/EUR_four_trait_manifest.tsv \
  --sumstats-dir /path/to/EUR/pleio_full_aligned_hdr \
  --ref-ld-prefix /path/to/eur_ldscores/@ \
  --out-dir /path/to/analysis/EUR_four_trait \
  --python /path/to/project3-pleio-preprocess/bin/python
```

The wrapper saves the fully resolved manifest beside the output directory and
then invokes PLEIO's `ldsc_preprocess.py`. The expected outputs are
`metain.txt.gz`, `sg.txt.gz`, and `ce.txt.gz`.

## 2. Derive the nested three-trait inputs

```bash
Rscript derive_three_trait_inputs.R \
  --input-dir /path/to/analysis/EAS_four_trait \
  --model-config config/EAS_three_trait_models.tsv \
  --out-root /path/to/analysis/EAS_three_trait

Rscript derive_three_trait_inputs.R \
  --input-dir /path/to/analysis/EUR_four_trait \
  --model-config config/EUR_three_trait_models.tsv \
  --out-root /path/to/analysis/EUR_three_trait
```

The nested inputs are exact column/submatrix selections from the corresponding
four-trait `metain`, `sg`, and `ce` files. No model parameters are re-estimated
at this stage.

## 3. Run PLEIO

Run the wrapper once for each prepared directory:

```bash
bash run_pleio.sh \
  --pleio-py /path/to/pleio/pleio.py \
  --input-dir /path/to/analysis/EAS_four_trait \
  --out-prefix /path/to/results/HELIOS_PLEIO_results \
  --python /path/to/project3-pleio-model/bin/python

bash run_pleio.sh \
  --pleio-py /path/to/pleio/pleio.py \
  --input-dir /path/to/analysis/EAS_three_trait/MDD_EA_G \
  --out-prefix /path/to/results/HELIOS_PLEIO_EAS_MDD_EA_G_results \
  --python /path/to/project3-pleio-model/bin/python
```

Repeat the second command for `SCZ_EA_G` and for the EUR `MDD_EA_CF` and
`SCZ_EA_CF` directories. The EUR four-trait output prefix used on EC2 was
`HELIOS_PLEIO_EUR_results`. All six final models used `--parallel --create` and
100,000 importance samples (the PLEIO default, made explicit by the wrapper).

For an exact recalculation, pass the matching retained calibration with
`--isf isf/<MODEL>.isf`. This replaces `--create` and avoids a new Monte Carlo
importance-sampling approximation:

| Analysis | Calibration file |
| --- | --- |
| EAS four-trait | `isf/EAS_FOURTRAIT.isf` |
| EAS MDD-EA-G | `isf/EAS_MDD_EA_G.isf` |
| EAS SCZ-EA-G | `isf/EAS_SCZ_EA_G.isf` |
| EUR four-trait | `isf/EUR_FOURTRAIT.isf` |
| EUR MDD-EA-CF | `isf/EUR_MDD_EA_CF.isf` |
| EUR SCZ-EA-CF | `isf/EUR_SCZ_EA_CF.isf` |

Use `--ncores` to fix the worker count.

## EC2 software provenance

- PLEIO Git commit `b18d507533ae39dcb6284573d3b36b91cdc13603`
- PLEIO log version 2.0; `ldsc_preprocess` log version 1.0
- preprocessing environment: Python 2.7.18, NumPy 1.16.5, pandas 0.24.2,
  SciPy 1.2.1
- model-run environment: Python 3.9.25, NumPy 2.0.2, pandas 2.3.3,
  SciPy 1.13.1

The raw/analysis-ready GWAS files and LD-score panels are governed by their
source studies and are not included here.
