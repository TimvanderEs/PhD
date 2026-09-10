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
6. downstream PLEIO directional classification and SMR/HEIDI reanalysis;
7. consolidation and comparison of the archived FUMA/MAGMA and GTEx v8
   enrichment outputs;
8. checksum-validated reconstruction of the significant g:Profiler tables and
   Supplementary Figure S8 source data.

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
- Six model-wide PLEIO and seven available single-trait FUMA SNP2GENE runs were
  selected with explicit reviewed manifests. EAS cognitive performance was not
  available as a single-trait FUMA run, so every affected comparison remains
  labelled incomplete.
- The 12 g:Profiler raw exports produce exactly 61 g:SCS-adjusted P < 0.05
  terms and reconcile with the final publication workbook. EUR query and custom
  background lists were recovered. The EAS background option and exact
  g:Profiler database release were not retained and are marked `not_recorded`.

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
- `inputs/fuma_magma/*.tsv` identify the exact archived FUMA runs used in the
  enrichment comparisons; their complete web settings are stored in the
  corresponding `results/fuma_magma/raw/*/*/params.config` files.
- `provenance/fuma_magma_raw_files.sha256` fingerprints all 52 archived FUMA
  parameter and raw MAGMA/GTEx files.
- `inputs/gprofiler/run_manifest.tsv` records the query/background files, raw
  row counts, analysis settings, timestamps, and both compressed and original
  CSV SHA-256 hashes.

These fingerprints establish which retained EC2 files were treated as final.
They are not checksums of raw participant or GWAS input data.
