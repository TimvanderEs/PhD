# HELIOS phenotype sensitivity to ancestry and self-described ethnicity

## Question

This analysis tested whether the 22k cognitive factor changed when:

1. the supplied seven-category GRID genetic-ancestry assignment replaced the original Chinese-reference ancestry dummies;
2. phenotype construction was restricted to the three core GRID groups used for the ancestry-specific GWASs;
3. participant-described ethnic heritage was added to the phenotype model; and
4. records whose core GRID assignment differed from their self-description were excluded.

All models retained the original age, age-squared, sex, age-by-sex, age-squared-by-sex and education terms. Participant-level inputs and results remained on the controlled EC2 instance. Only aggregate outputs and non-identifiable figures were exported.

## Reproduction of the retained phenotype

The retained six task residuals and cognitive-factor score were reproduced before any sensitivity analysis. Every task-level residual correlation was **r = 1.000000000000**, with a largest absolute numerical difference of **1.19 × 10^-10**. The reproduced PC1 correlated **r = 1.000000000000** with the retained score after accounting for the arbitrary PCA sign.

See [`baseline_reproduction.tsv`](tables/baseline_reproduction.tsv) and [`validation_checks.tsv`](validation_checks.tsv).

## Ancestry-coding audit

The retained PCA included **22,422** participants:

- **21,596 (96.32%)** belonged to the core GRID Chinese, Indian or Malay groups;
- **826 (3.68%)** belonged to a pairwise intermediate or Other GRID category; and
- all 826 intermediate/Other records had `indian = 0` and `malay = 0` in the retained phenotype table and therefore entered the original regression under its Chinese reference coding.

An additional **7 (0.03%)** core GRID Indian or Malay participants also had both original indicators set to zero. These coding details affected pooled phenotype construction, not the ancestry-specific GWAS membership: the intermediate/Other categories did not enter the Chinese, Indian or Malay GWAS strata.

See [`ancestry_coding_audit.tsv`](tables/ancestry_coding_audit.tsv).

## Sensitivity results

| Sensitivity | N | Loading congruence | PC1 score correlation | Interpretation |
|---|---:|---:|---:|---|
| Seven-category GRID adjustment versus retained phenotype | 22,422 | 0.99999994 | 0.999851 | No material change |
| Core GRID groups only versus retained phenotype | 21,596 | 0.99998224 | 0.999882 | No material change |
| Add self-described ethnicity within the same questionnaire subset | 13,102 | 0.99999997 | 0.999729 | No material change |
| Exclude discordant records within the questionnaire subset | 12,989 | 0.99999884 | 0.999982 | No material change |

The maximum loading difference was **0.00384** in the core-only refit and **0.00021** when self-described ethnicity was added. PC1 variance explained was 37.04% in the retained model, 37.04% with seven GRID categories and 36.77% in the core-only refit. The questionnaire-subset models explained 38.81% of task variance before and after self-described ethnicity was added.

Full estimates are in [`model_comparisons.tsv`](tables/model_comparisons.tsv) and [`loading_comparisons.tsv`](tables/loading_comparisons.tsv).

## Incremental contribution of self-described ethnicity

Within the same 13,102-person core-GRID questionnaire subset, adding participant-described ethnic heritage produced very small increments in task-level explained variance:

- reaction time: **0.060 percentage points** (FDR q = 0.074);
- Stroop box: **0.025 percentage points** (q = 0.346);
- Stroop ink: **0.026 percentage points** (q = 0.346);
- quiz: **0.006 percentage points** (q = 0.842);
- working memory: **0.075 percentage points** (q = 0.049); and
- pairing guesses: **0.108 percentage points** (q = 0.014).

Working memory and pairing guesses met the six-test FDR threshold, but the associated increments in explained variance were less than 0.11 percentage points and did not materially alter the multivariate cognitive factor. Statistical detectability should therefore be distinguished from practical importance.

See [`self_description_incremental_fit.tsv`](tables/self_description_incremental_fit.tsv).

## Interpretation

The original pooled phenotype model did not fully distinguish core, pairwise intermediate and Other GRID categories because the latter entered the Chinese reference coding. This should be described transparently. However, directly correcting the coding or excluding the intermediate/Other categories produced essentially the same loading pattern and participant scores. Likewise, adding participant-described ethnic heritage or excluding genetic/self-description-discordant records had negligible effects on PC1.

These results support retaining the published 22k cognitive factor for the primary genetic analyses. A GWAS rerun is not indicated solely by the distinction between genetic ancestry and self-described ethnicity. The thesis should nevertheless define the constructs separately and should not interpret the genetic ancestry strata as homogeneous social or cultural groups.

## Reproducibility

The analysis was run on 23 September 2026 under:

```text
/home/ec2-user/HELIOS_validation/ancestry_ethnicity_sensitivity_20260923
```

The executable script is [`../scripts/run_phenotype_ethnicity_sensitivity.R`](../scripts/run_phenotype_ethnicity_sensitivity.R). Input checksums, R session information and automated checks are retained in this directory. No participant identifiers, individual scores, task values, ethnic-heritage responses or PCA coordinates are included in the exported outputs.
