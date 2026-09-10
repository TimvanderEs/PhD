# Reproducibility and EC2 provenance

The public workflows were reconstructed from the scripts, internal program
logs, final-output timestamps, manuscript consolidation manifests, and output
checksums retained on the analysis EC2 instance. Failed preprocessing attempts
and explicitly labelled QC/older outputs were excluded.

## What can be reproduced

Readers with access to the governed GWAS summary statistics and the cited
reference panels can rerun:

1. the five recovered SumstatsQC runs and the documented QC settings for the
   remaining ancestry-trait combinations, conditional on source data;
2. EAS and EUR genome-wide LDSC/GenomicSEM covariance analyses, with raw and
   nearPD-adjusted matrices kept distinct;
3. the clean EAS pairwise LDSC sensitivity analysis;
4. four-trait EAS and EUR LAVA analyses after the exact EAS SCZ input identity
   is verified;
5. EAS and EUR four-trait PLEIO input preparation and association tests;
6. the nested MDD and SCZ three-trait PLEIO models;
7. downstream PLEIO directional classification and SMR/HEIDI reanalysis;
8. consolidation and comparison of the archived FUMA/MAGMA and GTEx v8
   enrichment outputs;
9. checksum-validated reconstruction of the significant g:Profiler tables and
   Supplementary Figure S7 source data.

The repository cannot make controlled or third-party GWAS data public. It
therefore supports computational reproduction conditional on lawful data
access, plus immediate execution of the fully synthetic directional-label
example.

## Decisions made during the audit

- The four-trait EAS `HELIOS` LAVA files are the relevant archived outputs, but
  their 726 bivariate rows were generated with a rounded `2e-5` screen and an
  incorrect input-info configuration. They are historical provenance, not a
  verified publication denominator. Applying `P <= 0.05/3064` to the archived
  files yields 698 rows. The publication workbook has 697 because one valid
  EA--SCZ row with a univariate P value of zero is absent. A corrected rerun is
  required; all three count sets remain separately labelled.
- Manuscript-target LAVA configurations use the disease counts in manuscript
  Table 1 and `1/0` for continuous traits. Exact recovered configurations are
  retained alongside them. The manuscript and Lam et al. primary GWAS report
  22,778 EAS SCZ cases, whereas the recovered SumstatsQC log reports 27,888;
  the exact summary-statistics file identity/header must be verified before
  rerunning.
- The primary EAS LDSC matrix contains HELIOS-10k. A later five-trait file that
  also added HELIOS-12k is retained only as a sensitivity/QC configuration.
- The primary LDSC workflow is `GenomicSEM::ldsc()` (GenomicSEM 0.0.5), and the
  publication workbook applies BH correction across all 12 ancestry-specific
  genetic correlations. The later standalone pairwise EAS run remains a
  sensitivity analysis. New outputs explicitly distinguish the raw covariance
  matrix from its nearPD transformation.
- The first EUR PLEIO preprocessing attempts failed because the input schema
  was not valid. The successful model used the `full_aligned_hdr` inputs with
  4,930,902 shared variants and no allele mismatches.
- The nested three-trait PLEIO files were derived by exact column and covariance
  submatrix selection from each successful four-trait preparation.
- PLEIO itself is referenced by commit rather than copied because the checked
  out software reports no redistribution licence.
- Six model-wide PLEIO and seven available single-trait FUMA SNP2GENE runs were
  selected with explicit reviewed manifests. EAS cognitive function was not
  available as a single-trait FUMA run, so every affected comparison remains
  labelled incomplete.
- The 12 g:Profiler raw exports produce exactly 61 g:SCS-adjusted P < 0.05
  terms and reconcile with the final publication workbook. EUR query and custom
  background lists were recovered. The EAS background option and exact
  g:Profiler database release were not retained and are marked `not_recorded`.
- Five SumstatsQC logs were recovered. EUR EA, EUR MDD, and EAS CF logs were
  not. The EAS MDD and SCZ sample-definition differences are preserved in the
  manifest rather than reconciled by assumption.
- All seven SMR/HEIDI merged inputs were recovered and checksummed. Original
  webSMR portal job IDs, software version, and access dates were not recovered;
  available file timestamps are not represented as portal run times.

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
  trait pair for the archived runner, exact-threshold audit, and publication
  workbook as separate analysis sets.
- `provenance/sumstatsqc_runs.tsv` records recovered QC settings, row counts,
  log hashes, and missing/conflicting sample metadata.
- `inputs/fuma_magma/*.tsv` identify the exact archived FUMA runs used in the
  enrichment comparisons; their complete web settings are stored in the
  corresponding `results/fuma_magma/raw/*/*/params.config` files.
- `provenance/fuma_magma_raw_files.sha256` fingerprints all 52 archived FUMA
  parameter and raw MAGMA/GTEx files.
- `inputs/gprofiler/run_manifest.tsv` records the query/background files, raw
  row counts, analysis settings, timestamps, and both compressed and original
  CSV SHA-256 hashes.
- `inputs/smr_heidi/run_manifest.tsv` records the seven input hashes and the
  limits of retained portal provenance.

These fingerprints establish which EC2 files were retained and how each is now
interpreted. The `provenance_status` field prevents an archived output from
being mistaken for a verified final analysis. They are not checksums of raw
participant or GWAS input data.

Publication-facing traits are labelled `EA` and `CF`. Historical identifiers
such as `Edu`, `G`, `COGENT`, and `HEL10k` remain only in source filenames or
upstream input identifiers when changing them would obscure exact provenance.
