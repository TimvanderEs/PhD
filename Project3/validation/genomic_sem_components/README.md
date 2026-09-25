# Cognitive-component Genomic SEM validation

This analysis tests whether the six HELIOS cognitive-task GWAS are adequately
described by one genetic factor or are better represented by correlated latency
and non-latency factors. The 10k and 22k freezes are analysed separately. They
are nested releases and are not treated as independent replications.

## Inputs

- The 10k analysis uses the six Chinese–Indian–Malay fixed-effect meta-analysis
  files ending in `filt2studies_hetp05_AFdiff05`.
- The 22k analysis uses the six release 2.3.1 provider archives, including the
  trans-ancestry effect, standard-error, P-value, allele-frequency and sample-
  size fields.
- [source_audit.tsv](source_audit.tsv) records checksums, row counts, sample-size
  ranges and schema checks for all 12 source files.

All effects are oriented toward higher cognitive performance before munging.
The timing traits and pairing-guesses count are therefore multiplied by -1;
quiz and working-memory scores retain their original direction. The complete
coding is in [input_orientation.tsv](input_orientation.tsv).

## Processing

The recovered allele-aware SG100K marker map is linked to rsIDs and hg19
coordinates before HapMap3 munging. Reproduction against the archived 22k
working-memory conversion yielded the same 2,627,227 mapped variants and an
identical SHA-256 over every field consumed by LDSC (`SNP`, `A1`, `A2`, `N`,
`Z`, and `P`). Of these markers, 1,580 have an unresolved or conflicting hg19
position in the broader legacy build bridge. They retain their validated rsID
and allele mapping; the position is not used by `munge_sumstats.py` or LDSC.

Summary statistics are munged against the EAS HapMap3 list with LDSC v1.0.1.
Multivariable LDSC and structural models use GenomicSEM 0.0.5 under R 4.4.3.
Two models are fitted by DWLS:

1. one factor loading on all six tasks;
2. correlated latency (reaction time, Stroop box and Stroop ink) and
   non-latency (quiz, working memory and pairing guesses) factors.

The scripts in [scripts](scripts) reconstruct the mapping, prepare signed
summary statistics, run munging and LDSC, and fit both models. Source GWAS and
large intermediate files are retained on the controlled EC2 workspace and are
not added to GitHub.

## Results

Both analyses completed on 24 September 2026. The task-level SNP-heritability
estimates were:

| Task | 10k h² (SE) | 22k h² (SE) |
|---|---:|---:|
| Reaction time | 0.070 (0.050) | 0.027 (0.020) |
| Stroop box | 0.146 (0.051) | 0.045 (0.022) |
| Stroop ink | 0.048 (0.055) | 0.020 (0.020) |
| Quiz | 0.060 (0.051) | 0.035 (0.021) |
| Working memory | 0.123 (0.062) | 0.063 (0.021) |
| Pairing guesses | 0.054 (0.054) | 0.057 (0.019) |

The two-factor model did not improve on the one-factor model in either freeze:

| Freeze | Δχ² (1 df) | Approx. P | ΔAIC, two minus one | One-factor SRMR | Two-factor SRMR |
|---|---:|---:|---:|---:|---:|
| 10k | 0.693 | 0.405 | +1.307 | 0.160 | 0.148 |
| 22k | 1.830 | 0.176 | +0.170 | 0.142 | 0.141 |

The P values are descriptive chi-square-difference approximations. They are
not treated as formal likelihood-ratio tests because the models were fitted
by DWLS and the one-factor case lies on the two-factor correlation boundary.

The 10k genetic covariance matrix was not positive semidefinite and required
substantial smoothing. Both 10k models also produced negative residual
variances. Its parameter estimates should therefore not be interpreted as a
valid latent structure. The 22k covariance matrix was positive semidefinite
and the models completed without these warnings, but the high SRMR, imprecise
loadings and non-significant nested-model test provide no evidence that the
two-factor model is preferable. In the 22k two-factor model, the estimated
factor correlation was 0.726 (SE 0.413, P = 0.079).

The result supports a limited conclusion: the early HELIOS component GWAS did
not establish a stable common genetic factor, and splitting the tasks into
latency and non-latency groups did not solve that problem. It does not test or
validate the later 60k two-factor model, which uses a larger freeze and must be
evaluated from its own inputs and results.

See [FINAL_REPORT.md](FINAL_REPORT.md) for the interpretation and thesis-ready
wording. Machine-readable results are in [`results/10k`](results/10k) and
[`results/22k`](results/22k); the run logs are retained in [`logs`](logs).
