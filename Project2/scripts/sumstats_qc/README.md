# Summary-statistics quality control

The summary statistics were processed with SumstatsQC v0.1 before
ancestry-specific harmonisation. Per-trait settings are recorded in
[`provenance/sumstats_qc_settings.tsv`](../../provenance/sumstats_qc_settings.tsv).

Original GWAS files and full logs are not redistributed because their access
and contents are governed by the source studies.

## Command pattern

The workflow used the following SumstatsQC stages for each input:

```bash
03_SumstatsQC_sumstats_munge.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX --multicpu=Y
04_SumstatsQC_refsumstats_merge.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX --multicpu=Y
05_SumstatsQC_QCparameters.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX \
  --AF=0.005 --AFB=0.15 --AMB=AMBIGUOUS_AF --INFO_score=0.3 \
  --qt=TRAIT_TYPE --multicpu=Y
06_SumstatsQC_Postprocessing.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX \
  --qt=TRAIT_TYPE --AF=0.005 --AFB=0.15 --AMB=AMBIGUOUS_AF --INFO_score=0.3
07_SumstatsQC_extractqcinfo.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX \
  --qt=TRAIT_TYPE --AF=0.005 --AFB=0.15 --AMB=AMBIGUOUS_AF --INFO_score=0.3
08_SumstatsQC_visualization.sh --sumstats=INPUT --pop=REFERENCE --prefix=PREFIX \
  --qt=TRAIT_TYPE --AF=0.005 --AFB=0.15 --AMB=AMBIGUOUS_AF \
  --INFO_score=0.3 --multicpu=Y
```

`REFERENCE` was `1000g_eur` or `1000g_eas`. `AMBIGUOUS_AF` was 0.35 for
EUR runs and 0.40 for EAS runs. Quantitative and binary
settings are recorded per trait in the manifest.
