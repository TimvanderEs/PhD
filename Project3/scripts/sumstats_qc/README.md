# Summary-statistics quality control

The manuscript reports SumstatsQC v0.1 before ancestry-specific harmonisation.
Five of the eight ancestry-by-trait run logs were recovered locally. This
directory records the exact settings visible in those logs and marks the three
missing logs explicitly; it does not invent commands or counts for them.

The machine-readable record is
[`provenance/sumstatsqc_runs.tsv`](../../provenance/sumstatsqc_runs.tsv).
Original GWAS files and full logs are not redistributed because their access
and contents are governed by the source studies.

## Recovered command pattern

The retained logs record the following SumstatsQC stages for each input:

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
recovered EUR runs and 0.40 for recovered EAS runs. Quantitative and binary
settings are recorded per trait in the manifest.

## Audit findings that require manuscript clarification

- The recovered EAS MDD log treats the file as quantitative with `N=194,542`,
  whereas the manuscript describes 15,771 cases and 178,777 controls
  (`N=194,548`).
- The recovered EAS SCZ log reports 27,888 cases and 35,362 controls, whereas
  the manuscript reports 22,778 cases and 35,362 controls.
- Logs for EUR EA, EUR MDD, and EAS CF were not recovered from the local
  circulation package or the accessible archive.
- The manuscript's post-QC shared variant sets (5,601,717 EUR and 4,717,676
  EAS) are not the same quantity as the variants ultimately tested by PLEIO
  (4,930,902 EUR and 4,556,091 EAS). The Methods should name the latter as the
  post-alignment PLEIO analysis universe rather than presenting both numbers
  as interchangeable.

Before submission, resolve the two EAS phenotype/count discrepancies against
the source GWAS documentation and recover the three missing logs if possible.
The repository remains transparent and executable without those controlled
inputs, but the missing provenance prevents a claim of complete raw-data
reproduction.
