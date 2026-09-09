# Reproducibility and EC2 provenance

The public workflows were reconstructed from the scripts, internal program
logs, final-output timestamps, manuscript consolidation manifests, and output
checksums retained on the analysis EC2 instance. Failed preprocessing attempts
and explicitly labelled QC/older outputs were excluded.

## What can be reproduced

Readers with access to the governed GWAS summary statistics and the cited
reference panels can rerun:

1. EAS and EUR genome-wide LDSC/GenomicSEM covariance analyses;
2. the clean EAS pairwise LDSC rerun used for final figure-facing estimates;
3. four-trait EAS and EUR LAVA analyses;
4. EAS and EUR four-trait PLEIO input preparation and association tests;
5. the nested MDD and SCZ three-trait PLEIO models;
6. downstream PLEIO directional classification and SMR/HEIDI reanalysis.

The repository cannot make controlled or third-party GWAS data public. It
therefore supports computational reproduction conditional on lawful data
access, plus immediate execution of the fully synthetic directional-label
example.

## Decisions made during the audit

- The primary EAS LAVA files are the four-trait `HELIOS` outputs. They contain
  all six pairings and 726 eligible bivariate tests. Older EAS outputs with
  three pairings were marked as QC and are not represented as final.
- The primary EAS LDSC matrix contains HELIOS-10k. A later five-trait file that
  also added HELIOS-12k is retained only as a sensitivity/QC configuration.
- The first EUR PLEIO preprocessing attempts failed because the input schema
  was not valid. The successful model used the `full_aligned_hdr` inputs with
  4,930,902 shared variants and no allele mismatches.
- The nested three-trait PLEIO files were derived by exact column and covariance
  submatrix selection from each successful four-trait preparation.
- PLEIO itself is referenced by commit rather than copied because the checked
  out software reports no redistribution licence.

## Verification files

- `provenance/software_versions.tsv` records the EC2 runtime versions.
- `provenance/eur_pleio_alignment_fingerprints.tsv` records the exact row
  counts and decompressed hashes reproduced by the recovered alignment script.
- `provenance/final_output_fingerprints.tsv` records final filenames, completion
  times, data-row counts, and SHA-256 hashes.
- `provenance/pleio_gws_counts.tsv` records the post-run genome-wide significant
  SNP counts for the six PLEIO models.
- `provenance/pleio_isf_fingerprints.tsv` verifies the six retained
  importance-sampling calibration files.
- `provenance/lava_pair_test_counts.tsv` records the bivariate-test counts by
  trait pair.

These fingerprints establish which retained EC2 files were treated as final.
They are not checksums of raw participant or GWAS input data.
