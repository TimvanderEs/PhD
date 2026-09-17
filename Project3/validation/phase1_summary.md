# Phase 1 summary

## Scope

Phase 1 reproduced the historical HELIOS 10k cognitive factor and examined measurement quality and PCA structure. No phenotype was redefined for association testing, and no GWAS was rerun.

## Baseline

- PCA sample: **7,403 participants**.
- PC1 variance explained: **34.49%**.
- Reconstructed versus saved historical score: **r > 0.999999999**; maximum aligned difference below **1e-8**.
- The normalized-phenotype table contains the same 7,403 saved `g` values with no numerical differences.
- All six loadings and task-PC1 correlations matched the archived reference within the pre-specified tolerance.
- The two validated trans-ancestry outputs remain separately identified by checksum (maximum variant-level N = 7,789 and 7,664); the refined Chinese freeze and historical `10k_clean` LDSC record are not conflated with them.

## Measurement

- Recorded source missingness was 0% across all six tasks, but the discrete tasks had much lower resolution than the latency measures.
- Effective categories: Quiz **6.0**, Working Memory **6.3**, Pairing **22.3**.
- The historical reciprocal Pairing transform converts **346 zero-error scores** to non-finite values. These participants are consequently excluded from the historical complete-case PCA.
- Reliability coefficients could not be estimated because repeated or item-level data were not available.
- The available covariate file contains one cohort/freeze label, so between-freeze comparisons are not estimable.

## PCA structure

- The bootstrap confirmed that the largest loadings are the three latency-derived indicators. Their 95% loading intervals were: Reaction Time **[0.482, 0.501]**, Stroop Box **[0.571, 0.586]**, and Stroop Ink **[0.489, 0.510]**.
- Leave-one-task-out solutions correlated **0.932–0.996** with the full PC1.
- Ancestry-specific loading congruence with the full solution ranged from **0.993 to 0.999**.
- The four-indicator PCA using one collapsed speed domain correlated **r = 0.881** with the original PC1.

## Interpretation

The historical phenotype is exactly reproducible and its speed-heavy structure is not a bootstrap artefact. The weighting partly reflects the covariance among three separately entered latency measures. The collapsed-speed analysis provides the planned domain-balance sensitivity test for a later phase, but it does not by itself establish that the historical score is invalid.

The Pairing zero-score transformation is the clearest construction issue identified in Phase 1. Any alternative phenotype work should pre-specify a transformation that retains zero-error participants, while preserving the historical phenotype unchanged as the primary archived analysis.

Phase 1 is complete. The next planned step is Phase 2 phenotype sensitivity work; it should begin only after review of these results.
