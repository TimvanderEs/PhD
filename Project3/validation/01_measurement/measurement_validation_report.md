# Task-level measurement validation

The source cognitive table contains **8,299 unique participants**. Recorded missingness is 0% for all six source measures, but their information content differs substantially.

Quiz, Working Memory, and Pairing are bounded/discrete measures with **12**, **18**, and **36** observed values, respectively. Their effective category counts are reported in `task_measurement_summary.tsv`.

The historical reciprocal transformation of Pairing guesses converts zero errors to infinity and then missing. This affects **346 participants (4.17%)**, including the best-performing end of that source scale. This is a property of the archived construction and was retained only to reproduce the historical phenotype.

Repeated measurements and item-level responses were not available, so test-retest reliability, alpha, omega, and split-half reliability were not estimated. Only one cohort/freeze label (`HELIOS10K`) was present; between-freeze distribution comparisons are therefore not estimable from these files.
