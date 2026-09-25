# Controlled LAVA validation design

The historical output cannot be repaired by relabelling or by reversing the arbitrary HELIOS PCA sign. KoGES EDU was passed to LAVA as a binary phenotype, and LAVA used the supplied case fraction when reconstructing marginal SNP correlations.

The minimum controlled rerun has two stages, both in a new directory.

## 1. Historical-control run

Run all four historical phenotypes (`MDD_EAS`, `Edu_EAS`, `SCZ_EAS`, and
`HEL10k_EAS`) with:

- the recovered LAVA v0.1.5 environment;
- the same four summary-statistic files;
- the same 3,064-locus file;
- the recovered 1000 Genomes EAS PLINK reference;
- the complete historical 4 x 4 sample-overlap matrix;
- `univ.thresh = 2e-5`; and
- the historical, incorrect EDU case/control declaration.

All four phenotypes are required because `process.input()` constructs a common
SNP universe across every supplied GWAS. A two-trait control could therefore
change the locus SNP sets even though only the HELIOS–EDU pair is of interest.
The resulting target-pair rows must reproduce the 92 archived HELIOS–EDU rows.
This control verifies that the moved reference path and reconstructed runner
are equivalent to the historical run.

## 2. Corrected quantitative-EDU run

Keep all four phenotypes and change only the EDU declaration:

```text
phenotype  cases  controls  filename
MDD_EAS    NA     NA        /home/ec2-user/project2/inputfiles/EAS_raw/MDD_EAS_for_LAVA.sumstats.gz
Edu_EAS    NA     NA        /home/ec2-user/project2/inputfiles/EAS_raw/Edu_EAS_for_LAVA.sumstats.gz
SCZ_EAS    13305  16244     /home/ec2-user/project2/inputfiles/EAS_raw/SCZ_EAS_for_LAVA.sumstats.gz
HEL10k_EAS NA     NA        /home/ec2-user/project2/inputfiles/EAS_raw/HELIOS_10k_for_LAVA.sumstats.gz
```

Retain every other file, overlap value, locus, reference, LAVA version and
threshold from the historical-control run. Write all univariate and bivariate
rows to the new audit directory. Do not overwrite the December 2025 outputs.

Only if the historical-control run reproduces the archived target rows should the corrected output be used for local covariance decomposition.

The control and corrected configurations may run concurrently, but each uses
the historical single-worker locus loop. This preserves the original
per-locus execution path while avoiding duplicate wall-clock time.

## Completion

Both runs completed on 19 September 2026. The historical control reproduced
the complete row sets and all deterministic quantities exactly. LAVA's
unseeded Monte Carlo intervals and P values were highly concordant and produced
identical nominal and FDR decisions. The corrected run retained all 92
HELIOS–education regions and is the source of the final analysis.

Exact outputs and logs are in [`rerun/results`](rerun/results). The validation
statistics are in
[`../provenance/historical_control_comparison.tsv`](../provenance/historical_control_comparison.tsv)
and
[`../provenance/control_corrected_impact.tsv`](../provenance/control_corrected_impact.tsv).
