# Chinese-only summary-statistic freeze

`freeze_manifest.tsv` identifies the final Chinese 10k build-37 file and the
exact derivatives used by downstream analyses. Full summary statistics are not
stored in this code repository.

The final file passed structural and numeric validation across all 5,615,100
rows. Cross-format checks were performed on normalized common columns so that
the LAVA/MiXeR and MTAG/PLEIO derivatives could be compared independently of
their analysis-specific headers.

The 10k trans-ancestry results are documented separately in
[`../trans_ancestry`](../trans_ancestry).
