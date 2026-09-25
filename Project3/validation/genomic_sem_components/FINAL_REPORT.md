# Six-task Genomic SEM validation

## Question

This analysis asked whether the six task-level HELIOS GWAS from the 10k and
22k freezes support either a single genetic cognitive factor or correlated
latency and non-latency factors. The two freezes contain overlapping
participants and are not independent replications.

## Inputs and methods

All six 10k trans-ancestry meta-analysis files and all six 22k release 2.3.1
archives passed compression, schema and row-width checks. Timing measures and
pairing guesses were sign-reversed so that every effect represented better
cognitive performance. Summary statistics were munged to EAS HapMap3 variants
with LDSC v1.0.1. Genetic covariance matrices were estimated with
`GenomicSEM::ldsc()` and fitted by DWLS in GenomicSEM 0.0.5.

The one-factor model loaded all six tasks on `g`. The two-factor model assigned
reaction time, Stroop box and Stroop ink to a latency factor, and quiz, working
memory and pairing guesses to a non-latency factor. The factors were allowed
to correlate.

## Heritability estimates

| Task | 10k h² (SE) | 10k P | 22k h² (SE) | 22k P |
|---|---:|---:|---:|---:|
| Reaction time | 0.070 (0.050) | 0.167 | 0.027 (0.020) | 0.177 |
| Stroop box | 0.146 (0.051) | 0.0045 | 0.045 (0.022) | 0.040 |
| Stroop ink | 0.048 (0.055) | 0.389 | 0.020 (0.020) | 0.329 |
| Quiz | 0.060 (0.051) | 0.243 | 0.035 (0.021) | 0.094 |
| Working memory | 0.123 (0.062) | 0.0477 | 0.063 (0.021) | 0.0026 |
| Pairing guesses | 0.054 (0.054) | 0.318 | 0.057 (0.019) | 0.0025 |

These are univariate LDSC estimates within each multivariable run. The 10k
standard errors are large, and four of its six estimates are not different
from zero at P < 0.05. The larger 22k freeze improves precision, although only
working memory and pairing guesses are clearly different from zero and Stroop
box is nominally significant.

## Structural model comparison

| Freeze | Model | χ² | df | AIC | CFI | SRMR |
|---|---|---:|---:|---:|---:|---:|
| 10k | One factor | 4.983 | 9 | 28.983 | 1.000 | 0.160 |
| 10k | Two factors | 4.290 | 8 | 30.290 | 1.000 | 0.148 |
| 22k | One factor | 8.899 | 9 | 32.899 | 1.000 | 0.142 |
| 22k | Two factors | 7.069 | 8 | 33.069 | 1.000 | 0.141 |

For 10k, the extra factor gave a descriptive Δχ²(1) = 0.693, approximate
P = 0.405, and increased AIC by 1.307. For 22k, the corresponding values were
Δχ²(1) = 1.830, approximate P = 0.176, and an AIC increase of 0.170. These are
not treated as formal likelihood-ratio tests because the models were fitted
by DWLS and the one-factor case lies on the two-factor correlation boundary.
The small SRMR changes also show that the two-factor specification did not
materially improve residual fit. CFI should not be interpreted alone here,
particularly because SRMR remained high and the 10k matrix was unstable.

The 22k standardized one-factor loadings were 0.589 for reaction time, 0.799
for Stroop box, 0.706 for Stroop ink, 0.781 for quiz, 0.471 for working memory
and 0.081 for pairing guesses. Only Stroop box reached P < 0.05. In the 22k
two-factor model, the factor correlation was 0.726 (SE 0.413, P = 0.079).
Pairing guesses had a near-zero loading on the non-latency factor. This is task
heterogeneity, but it is not a clean latency-versus-non-latency separation.

## Quality-control findings

The 10k genetic covariance matrix had a minimum eigenvalue of -0.0586.
GenomicSEM substantially smoothed the matrix, and both fitted models contained
negative residual variances. Some genetic correlations exceeded the valid
range before smoothing. The 10k factor estimates are therefore exploratory
and should not be interpreted substantively.

The 22k matrix had a minimum eigenvalue of 0.00223 and required no smoothing.
Neither fitted model produced a negative residual variance. This makes the
22k comparison technically more reliable, but it remains weakly powered at
the task level and neither model achieved good residual fit.

## Conclusion

The rerun explains why the early Genomic SEM work did not yield a defensible
common-factor GWAS. The 10k covariance structure is too unstable for latent
interpretation, and the better-behaved 22k data do not favour two factors over
one. The later 60k result may support a latency/non-latency structure because
of its greater power, but it is independent evidence and is not validated by
this analysis.

Suggested thesis wording:

> Exploratory Genomic SEM analyses of six task-level GWAS did not identify a
> stable latent structure in the early HELIOS freezes. The 10k genetic
> covariance matrix was non-positive-semidefinite and required substantial
> smoothing. In the 22k freeze, the covariance matrix was estimable without
> smoothing, but a correlated latency/non-latency model did not improve on a
> single-factor model (descriptive Δχ²(1) = 1.83, approximate P = 0.176), and
> residual fit remained poor. We therefore did not use these models to
> construct a multivariate cognitive GWAS.

## Reproducibility record

- [Source audit](source_audit.tsv)
- [Direction and factor coding](input_orientation.tsv)
- [Pipeline scripts](scripts)
- [10k machine-readable results](results/10k)
- [22k machine-readable results](results/22k)
- [Run logs](logs)
- [Compact archive checksums](archive_checksums.sha256)
- [EC2 intermediate-file checksums](ec2_intermediate_checksums.sha256)

Large source GWAS and intermediate munged files remain in the controlled EC2
workspace and are not deposited in GitHub.
