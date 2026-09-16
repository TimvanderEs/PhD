# SG100K covariate-adjusted LD reference

These scripts capture the recovered workflow used to construct the SG100K
covariates and calculate covariate-adjusted LD scores.

`prepare_sg100k_covariates.sh` performs the recorded preparation steps:

- MAF at least 0.01;
- exclusion of predefined high-LD regions;
- LD pruning with `--indep-pairphase 500 50 0.2`;
- removal of related samples with `--king-cutoff 0.0884`; and
- calculation of 50 principal components in unrelated participants.

`compute_covariate_adjusted_ld_scores.sh` runs cov-LDSC v1.0.0 for chromosomes
1–22 with a 20 cM window and six concurrent chromosome jobs. The EC2 logs
record 37,364 individuals in the chromosome-specific inputs and successful
completion of all chromosome jobs in June 2025.

The scripts require controlled SG100K genotype inputs and therefore cannot be
run from the public repository alone.
