# Genotype quality control

The recovered Chinese genotype prefix records the filters
`--geno 0.02`, `--maf 0.005` and `--hwe 1e-6`. Equivalent Chinese, Indian and
Malay prefixes were named in the historical scripts.

`run_ancestry_genotype_qc.sh` consolidates the recovered commands into a
path-independent workflow. It applies:

- genotype and sample missingness below 2%;
- MAF at least 0.005 and HWE P greater than 1e-6;
- F-statistic between -0.2 and 0.2;
- removal of participants in pairs with PI-HAT greater than 0.75; and
- 20 ancestry-specific PCs from common, LD-pruned variants outside the
  supplied high-LD regions.

These settings reproduce the retained genotype prefix used by the completed
REGENIE runs.
