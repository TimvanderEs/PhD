# Genotype quality control

The recovered Chinese genotype prefix records the filters
`--geno 0.02`, `--maf 0.005` and `--hwe 1e-6`. Equivalent Chinese, Indian and
Malay prefixes were named in the historical scripts.

`reconstruct_ancestry_genotype_qc.sh` is a cleaned, path-independent
reconstruction of those scripts. It applies:

- genotype and sample missingness below 2%;
- MAF at least 0.005 and HWE P greater than 1e-6;
- F-statistic between -0.2 and 0.2;
- removal of participants in pairs with PI-HAT greater than 0.75; and
- 20 ancestry-specific PCs from common, LD-pruned variants outside the
  supplied high-LD regions.

The manuscript drafts also describe a MAF threshold of 0.01 and sample filters
for missingness, heterozygosity, sex concordance and relatedness. These are
draft method descriptions for a later dataset and do not replace the recovered
10k thresholds. The original completed PLINK logs were not present on the
restored EC2 volume, so the reconstruction must be checked against the retained
10k genotype files before it is labelled as the executed QC workflow. No sex
concordance command was found in the supplied archive.
