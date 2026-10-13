# Project 4 evidence audit

## Scope

This audit covers:

- the current five-page manuscript draft;
- the `precise-covldsc-main` code archive;
- the documented pooled 20 cM cov-LDSC outputs;
- the invalid 40- to 45-row external-cohort benchmark and its input provenance; and
- the later Chinese, Indian, and Malay 1 cM LD-score panels and LDSC results on EC2; and
- the final controlled five-panel comparison, reference-size sensitivity, and held-out simulation completed on 22 September 2026.

The purpose of this document is to distinguish results that are ready to report from preliminary analyses and unresolved provenance.

## Main conclusion

The historical manuscript headline must be replaced because its main quantitative claims came from an invalid external-cohort benchmark. The archived download script assigned trait names to unrelated GWAS Catalog accessions, and its formatters substituted fixed sample sizes. Those rows are excluded. A controlled rerun using genuine Biobank Japan height summary statistics produced h2 = 0.3875 (SE = 0.0188), compared with 0.0067 (SE = 0.0026) for the historical file labelled as BBJ height. The input audit is in [`validation/`](validation/).

The latest documented SG100K analysis is separate from the invalid accession-based benchmark. It constructs Chinese, Indian, and Malay Singapore panels with a 1 cM window and uses six semantically identified SG100K release 2.3.1 GWAS. The GWAS inputs in the documented pipeline are fixed-effect trans-ancestry estimates.

A preliminary ancestry-matched validation completed 63 local-panel runs and 21 Chinese-GWAS/1KG-EAS comparisons. Its translated-row 1KG comparator was subsequently superseded by an exact comparator recalculated from raw 1KG EAS genotypes.

The final controlled experiment used one 874,025-variant universe and identical analytical settings across exact 1KG EAS, Singapore Chinese, PC-adjusted Chinese, pooled Singapore, and PC-adjusted pooled panels. Standard Chinese and 1KG EAS estimates were effectively equivalent: mean h2 was 0.1618 versus 0.1601, mean SE was 0.02035 versus 0.02022, none of six paired h2 differences was nominally significant, and rg estimates correlated at 0.99945. Reference-size and pilot simulation results also did not demonstrate a local-panel precision or accuracy advantage.

The strongest defensible chapter narrative is therefore:

> A population-specific LD-reference development study showing that pooled multi-ancestry LD scores are sensitive to covariate specification, whereas a standard Singapore Chinese panel and 1000 Genomes EAS give nearly equivalent LDSC results under a controlled 1 cM design.

## Historical analyses

### Pooled Singapore 20 cM cov-LDSC panel

- Initial merged sample: 46,508 participants.
- Recorded unrelated sample in the chromosome-level cov-LDSC logs: 37,364 participants.
- LD-score window: 20 cM.
- Covariates supplied to cov-LDSC: 50 principal components.
- The runner does not include age or sex.
- The archived code and README disagree about GRCh37 versus GRCh38.
- The external benchmark compared this 20 cM panel with 1 cM 1000 Genomes panels.

This panel is suitable as an exploratory or sensitivity analysis until window-size and covariate effects are separated.

### Ancestry-specific Singapore 1 cM panels

The latest documented workflow created complete autosomal panels for three groups:

| Panel | Unrelated participants | LD-score SNPs | Sumstat–LD overlap | Files complete |
| --- | ---: | ---: | ---: | --- |
| Chinese | 27,181 | 878,093 | 827,024 | 22/22 chromosomes |
| Indian | 5,848 | 912,313 | 859,224 | 22/22 chromosomes |
| Malay | 4,327 | 909,907 | 856,937 | 22/22 chromosomes |

Construction used ancestry-specific genotype files whose filenames contain `HELIOS`, the KING unrelated-sample keep list, SNP-only variants, MAF >= 0.01, genotype missingness <= 0.02, HWE P >= 1e-6, genetic-map cM positions, and a 1 cM LDSC window. No covariate file was supplied, so these are standard rather than covariate-adjusted LD scores, even though the cov-LDSC fork of `ldsc.py` was used to calculate them.

The chromosome logs and outputs are complete. The original panel-construction script and exact extraction list for this historical workflow were not available; a similarly named surviving file contains 925,225 unique variants despite the `867k` label. This limitation applies to the documented panels only. The final controlled experiment uses reconstructed scripts and an exact 874,025-variant manifest.

### GWAS-to-panel pairing

The provider files contain separate Chinese, Indian, Malay, and trans-ancestry result columns. However, `/home/ec2-user/sg100K_sumstats/prep_for_munge.sh` explicitly selects only the `FREQ_trans`, `BETA_trans`, `SE_trans`, `P_trans`, and `N_trans` fields. All documented h2 and rg logs then reuse those six trans-ancestry files with the Chinese, Indian, and Malay LD panels. No ancestry-specific munged file or matching LDSC run was available.

The documented results are therefore a panel-sensitivity analysis for trans-ancestry GWAS, not an ancestry-matched benchmark. The new validation is reported separately under [`results/ancestry_matched/`](results/ancestry_matched/) and [`validation/ancestry_matched_sg100k/`](validation/ancestry_matched_sg100k/). The full design reconstruction is in [`validation/analysis_design_reconstruction.md`](validation/analysis_design_reconstruction.md).

### Preliminary ancestry-matched validation

The ancestry-specific columns were harmonised to the 867,103-variant Singapore consensus set. Quantitative-trait sample sizes were approximately 33,200 Chinese, 7,400 Indian, and 5,600 Malay participants. T2D comprised 3,037/30,214 Chinese, 1,546/5,892 Indian, and 903/4,750 Malay cases/controls.

For Chinese GWAS, mean h2 was 0.1452 with the Singapore Chinese panel and 0.1231 with 1KG EAS. Trait-level Singapore estimates were 10.1%--22.6% larger, but mean SE was 0.0183 versus 0.0149 and mean h2/SE was 7.28 versus 7.55. Mean intercepts were 1.0050 and 1.0052. Larger h2 therefore reflects a panel-dependent scale shift and is not evidence of better accuracy or precision.

Chinese rg estimates were highly concordant across panels (r = 0.9990, mean absolute difference = 0.0111, 15/15 directions concordant). Six trait pairs passed BH correction within the 15-pair family under both panels: SBP--BMI, SBP--T2D, WHR--BMI, WHR--T2D, height--BMI, and BMI--T2D.

Only height and BMI had h2/SE >= 1.96 in the Indian and Malay analyses. No Indian or Malay rg survived BH correction. These analyses document the matched estimates but are too weak to determine whether the local panels improve inference. No suitable public South Asian comparator was available on EC2.

The Chinese/1KG comparison in this preliminary stage used a translated row subset of existing 1KG LD-score files. It remains useful as provenance but is not the final performance comparison.

### Final controlled reference experiment

The final experiment recalculated the 1KG EAS comparator from the raw genotypes and regenerated four Singapore panels under the same GRCh38 map, 1 cM window, SNP set, variant order, regression weights, and 200 jackknife blocks. The raw 1KG release lacked one requested chromosome 5 marker, so all panels were intersected to 874,025 variants.

| Panel | Reference N | Mean h2 | Mean SE | Mean h2/SE | Mean intercept |
| --- | ---: | ---: | ---: | ---: | ---: |
| Exact 1KG EAS | 504 | 0.16008 | 0.02022 | 7.2326 | 1.00720 |
| Singapore Chinese, standard | 27,178 | 0.16180 | 0.02035 | 7.2935 | 1.00525 |
| Singapore Chinese, 20 PCs | 27,178 | 0.16202 | 0.02040 | 7.2682 | 1.00535 |
| Singapore pooled, standard | 37,351 | 0.14752 | 0.02005 | 6.7603 | 1.00693 |
| Singapore pooled, 20 PCs | 37,351 | 0.17053 | 0.02162 | 7.2221 | 1.00605 |

The standard Chinese-versus-1KG h2 differences ranged from 0.0002 to 0.0037. None was nominally significant; BH-adjusted P values ranged from 0.322 to 0.924. Across 15 rg pairs, r = 0.99945, the mean absolute difference was 0.00736, and all directions agreed.

Downsampling the Chinese panel from 27,178 participants to 10,000, 2,500, or 500 changed mean h2 and SE very little. LD-score correlation with the full panel improved monotonically with N, confirming that reference size affects LD-score sampling noise, but that noise was not the limiting uncertainty for these six GWAS.

Five held-out simulated traits produced almost identical results across 1KG EAS, the full overlapping Chinese panel, and a disjoint 10,000-person Chinese panel. Because five replicates do not provide adequate power for small differences in bias or coverage, this is a calibration check rather than proof of equal performance.

PC adjustment had little effect within the Chinese panel: mean h2 changed by 0.00022 and no trait survived BH correction. It had a larger effect in the pooled panel: mean h2 increased by 0.0230 and five of six differences survived correction. This supports specification sensitivity in pooled multi-ancestry LD scores, not greater accuracy of either pooled model.

### Cohort and naming provenance

The evidence supports identifying the 46,508-person reference as PRECISE-SG100K-derived rather than HELIOS-only:

- the published HELIOS cohort contains 10,004 participants;
- the merged genotype log contains 46,508 participants and records the prefix `precise_sg100k`;
- the ancestry-specific genotype sources contain 33,392 Chinese, 7,449 Indian, and 5,666 Malay participants, totalling 46,507 before unrelatedness filtering;
- the six benchmark GWAS directories are stored under `sg100K_sumstats`, identify genomic data release 2.3.1, and report approximately 46,200--46,300 participants per quantitative trait; and
- those GWAS were run within ancestry and combined by fixed-effect meta-analysis.

The internal genotype prefix `HELIOS_*` is therefore not evidence that the panels use HELIOS participants alone. The manuscript should use **PRECISE-SG100K-derived reference panel** and identify the data release. HELIOS should be reserved for analyses that use the 10,004-participant HELIOS cohort.

## Historical trans-ancestry 1 cM results

### Heritability

Six SG100K release 2.3.1 traits were evaluated with the 1000 Genomes EAS 1 cM panel and each of the three Singapore 1 cM panels.

| LD panel | Mean h2 | Mean SE(h2) | Mean h2/SE | Mean intercept |
| --- | ---: | ---: | ---: | ---: |
| 1000 Genomes EAS, 1 cM | 0.12737 | 0.01290 | 8.9717 | 0.98028 |
| Singapore Chinese, 1 cM | 0.12525 | 0.01343 | 8.5026 | 0.98158 |
| Singapore Indian, 1 cM | 0.12920 | 0.01560 | 7.5371 | 0.99277 |
| Singapore Malay, 1 cM | 0.12372 | 0.01355 | 8.2985 | 0.98188 |
| Pooled Singapore cov-LDSC, 20 cM | 0.05257 | 0.00605 | 7.8012 | 0.99558 |

Relative to 1000 Genomes EAS, mean h2 was 1.7% lower with the Chinese panel, 1.4% higher with the Indian panel, and 2.9% lower with the Malay panel. These differences are small. Mean h2/SE was highest with 1000 Genomes EAS, so the latest results do not demonstrate a precision or power advantage for the Singapore panels.

The pooled 20 cM panel produced much smaller h2 estimates for all six traits. This pattern is consistent with its larger LD scores changing the LDSC slope scale. It must not be interpreted as a direct population-panel effect because window size and panel construction changed simultaneously.

T2D requires separate caution: mean chi-square was below one and intercepts were approximately 0.897 across the 1 cM panels. This suggests under-dispersion or an input/statistical-scale issue and weakens interpretation of its h2 and rg results.

### Genetic correlation

All 45 planned runs were available: 15 trait pairs evaluated with each of the three ancestry-specific Singapore panels.

| Panel comparison | Pearson correlation of rg | Same direction |
| --- | ---: | ---: |
| Chinese versus Indian | 0.9837 | 14/15 |
| Chinese versus Malay | 0.9965 | 15/15 |
| Indian versus Malay | 0.9850 | 14/15 |

Six trait pairs were nominally significant with all three panels: SBP–BMI, SBP–T2D, T2D–BMI, WHR–BMI, WHR–SBP, and WHR–T2D. The significance labels are currently based on unadjusted P < 0.05 and should be replaced with a stated multiple-testing procedure or described as nominal.

The defensible conclusion is that the rg structure is robust to the choice among the three ancestry-specific Singapore LD panels. The pooled trans-ancestry GWAS does not have one perfectly matched ancestry-specific reference, so agreement across panels is more informative than selecting one panel as the best.

## Why the manuscript headline must change

The draft reports mean h2 of 0.166 for Singapore and 0.085 for 1000 Genomes EAS. These values can be reproduced from two early 40-row univariate result tables, but the underlying trait/cohort identities are invalid. For example, the file called BBJ height was downloaded from `GCST90018907`, which is psoriasis vulgaris. The other JBB labels are similarly incorrect; most CBK and KOGES/TWB labels are also unrelated to their Catalog phenotypes, and several UKB labels point to vitamin D studies. In addition, one formatter used N = 160,000 for JBB and another used N = 400,000 when N was not recognised.

Even before the accession failure was identified, the aggregate diagnostics did not support greater accuracy:

- one Singapore estimate is greater than 1 and strongly affects the mean;
- median h2 is 0.0881 for Singapore and 0.0576 for 1000 Genomes EAS;
- mean SE(h2) is 0.02092 for Singapore and 0.00789 for 1000 Genomes EAS;
- mean h2/SE is 8.92 for Singapore and 9.81 for 1000 Genomes EAS; and
- mean intercept is 1.299 for Singapore and 1.291 for 1000 Genomes EAS.

The external-cohort estimates must now be excluded completely. They cannot be rehabilitated by relabelling because phenotype identity, cohort identity and N are analysis inputs. The six-trait SG100K comparison remains usable because it was traced independently to semantic phenotype identifiers and provider metadata. It should be described as a trans-ancestry GWAS panel-sensitivity analysis, not as an ancestry-matched comparison.

The terms `95.3% increase in heritability`, `70.5% improvement in discovery power`, `superior`, and `attenuation-bias correction` should be removed unless a controlled comparison or simulation supports them. h2/SE may be reported as a precision diagnostic, but it should not be labelled discovery power.

## Manuscript inconsistencies

| Topic | Current draft/archive | Required resolution |
| --- | --- | --- |
| Resource name | Draft calls a 46,508-person resource HELIOS-derived | Merge provenance identifies `precise_sg100k`; use PRECISE-SG100K-derived reference panel and reserve HELIOS for the 10,004-person cohort |
| GWAS-panel matching | Result directories are named Chinese, Indian, and Malay | Those labels identify the LD panel only; every run used the same trans-ancestry SG100K GWAS files |
| Primary method | Draft describes a pooled cov-LDSC resource; an author comment says the analysis switched back to standard LDSC | Use the final controlled 1 cM five-panel experiment as the primary comparison; retain the historical 20 cM result as development provenance |
| Genome build | GRCh38 in the title/methods; GRCh37 in part of the README; hg19 high-LD file named in logs | Verify genotype, genetic-map, exclusion-region, and GWAS builds separately and record conversions |
| Final unrelated N | Different workflows produced different counts | Use 37,364 for the June 2025 historical pooled panel; 27,181/5,848/4,327 for the documented ancestry panels; and 27,178 Chinese plus 37,351 pooled for the final controlled experiment |
| Covariates | Draft says age, sex, and 50 PCs | Documented 20 cM runner supplies 50 PCs only; 1 cM panels supply no covariates |
| Window | Draft presents 20 cM as the panel setting and compares it with standard panels | Treat 20 cM versus 1 cM as a confounded historical comparison; use the matched 1 cM experiment for panel claims |
| External benchmark | Draft says 45 trait–cohort combinations | Accession audit shows incorrect phenotype labels and fixed-N substitution; exclude all accession-derived UKB/JBB/CBK/KOGES_TWB rows and rebuild only from independently verified sources |
| Headline h2 | 0.166 versus 0.085 | Derived from mislabeled inputs and invalid; remove rather than reinterpret |
| Intercepts | Draft reports means below one | Documented external benchmark means are approximately 1.30; the latest six-trait analysis is near one except T2D |
| Software citation | README cites Luo and Li, Bioinformatics 2018 | Cite Luo et al., Human Molecular Genetics 2021, DOI 10.1093/hmg/ddab130; 2018 is the code copyright date |
| Code completeness | The source archive was incomplete | Use the reconstructed, parameterised scripts and retained controlled run records in this repository |

## Evaluation framework

| Metric | What supports a useful panel | What raises concern | Direction of improvement |
| --- | --- | --- | --- |
| h2 | Agreement with an external benchmark, simulation truth, or genotype-based estimate; stability across reasonable specifications | Large shifts caused by window size, outliers, negative or greater-than-one estimates | Neither higher nor lower is inherently better |
| SE(h2) | Smaller SE without worse calibration or reduced SNP coverage | Smaller SE produced by an invalid scale, constrained model, or fewer difficult traits | Lower is preferable only with comparable estimands and calibration |
| h2/SE | Stable or larger values under an otherwise matched analysis | Relabelling the statistic as discovery power; increases driven only by h2 inflation | Larger can indicate precision, not correctness or discovery power |
| rg | Similar direction and magnitude across plausible panels; agreement with known or independently estimated relationships | Sign changes, estimates outside the valid range, or dependence on one panel | Agreement and stability matter more than magnitude |
| SE(rg) | Smaller paired SE with stable rg and adequate univariate h2 | Large SE, failed rg caused by weak h2, or selective exclusion of failed traits | Lower is preferable for the same valid estimand |
| LDSC intercept | Values compatible with the GWAS processing and mean chi-square | Strong deviation from one, unexplained panel shifts, or treating values below one as proof of no confounding | Closer to one is often reassuring but not a stand-alone ranking |
| Ratio | Interpretable alongside mean chi-square and intercept | Undefined, negative, or unstable values; comparison when mean chi-square is close to or below one | No universal better direction |
| Mean chi-square | Sufficient signal for LDSC and consistency across panels using the same summary statistics | Values below one, minimal signal, or panel-dependent values that should be invariant | Context dependent |
| SNP count/overlap | High, documented, comparable overlap with correct alleles and build | Panel differences caused by unequal SNP universes or unresolved identifiers | More is useful only after QC and harmonisation |
| Convergence/failures | All prespecified runs complete, with failures retained in the denominator | Silent omission of failed or implausible results | Higher completion is preferable |
| Cross-cohort replication | Similar panel effects in cohorts of matching ancestry | Apparent benefit restricted to one cohort or ancestry-mismatched analyses | Greater consistency is preferable |

## Required analyses before a definitive manuscript

### Essential

- [x] Use standard ancestry-specific 1 cM LD scores as the primary estimand; treat pooled covariate-adjusted scores as a sensitivity analysis.
- [x] Reconstruct and archive the exact controlled panel build scripts and 874,025-variant intersection manifest.
- [x] Archive parameterised h2 and rg runners and generate the new validation tables directly from LDSC logs.
- [ ] Verify the build of the genotype data, genetic map, high-LD exclusion regions, and every benchmark summary-statistic set.
- [ ] Run a controlled window analysis at 1, 2, 5, 10, and 20 cM while holding samples, SNPs, and covariates fixed.
- [x] Run a controlled 1 cM covariate analysis comparing standard and 20-PC scores in Chinese and pooled reference samples.
- [x] Audit the historical external inputs. The accession-based UKB/JBB/CBK/KOGES_TWB benchmark is invalid and excluded.
- [ ] If an external benchmark remains a chapter aim, define a new prespecified cohort–trait manifest using independently verified source files, phenotype definitions and sample sizes, then rerun it from the beginning.
- [x] Apply BH correction across the 15 rg tests within each ancestry/reference analysis in the new matched validation.
- [ ] Investigate T2D under-dispersion before treating T2D results as confirmatory.
- [x] Prepare the retained Chinese, Indian, and Malay GWAS columns separately and run all matched local-panel h2 and rg analyses.
- [ ] Add prespecified public comparators for Indian and Malay; the current EC2 archive contains only the valid 1KG EAS comparator used for Chinese GWAS.

### Strongly recommended

- [x] Run a held-out known-h2 simulation as a pipeline calibration check; expand beyond five replicates before using it to rank panels.
- [x] Report trait-level Chinese panel differences using paired delete-block jackknife SEs under an exactly matched design.
- [x] Report medians alongside means in the new validation summaries.
- [ ] Predefine rules for impossible, negative, non-convergent, or weak-signal estimates.
- [x] Stratify the new validation by Chinese, Indian, and Malay GWAS and exclude the invalid UK/European benchmark from interpretation.
- [ ] Verify whether pooled trans-ancestry SG100K summary statistics require a mixture-weighted or cohort-specific reference rather than choosing one ancestry panel.

## Proposed thesis chapter structure

### 1. Chapter overview

State the LD-reference problem, the Singapore reference resource, the two-stage evaluation, and the link to the HELIOS cognitive GWAS chapter.

### 2. Introduction

1. LDSC and dependence on reference LD.
2. Bias from ancestry and population mismatch.
3. Genetic heterogeneity in Singapore.
4. Rationale and assumptions of cov-LDSC.
5. Difference between a pooled admixed reference and ancestry-specific references.
6. Study aims and falsifiable evaluation criteria.

### 3. Methods

1. Confirmed Singapore genotype resource, data release, and ancestry groups.
2. Sample QC, KING filtering, and ancestry definitions.
3. Variant QC and genome-build verification.
4. Genetic-map preparation.
5. Pooled 20 cM cov-LDSC panel.
6. Ancestry-specific 1 cM standard panels.
7. GWAS benchmark datasets and trait manifest.
8. Summary-statistic harmonisation and common SNP universes.
9. Univariate and bivariate LDSC.
10. Paired comparison metrics, failure rules, and multiple testing.
11. Window and covariate sensitivity analyses.
12. Code, metadata, and LD-score release.

### 4. Results

1. Sample and variant flow.
2. Properties of the pooled 20 cM panel.
3. Identification of the window-size sensitivity.
4. Construction and QC of the ancestry-specific 1 cM panels.
5. Heritability agreement with 1000 Genomes EAS.
6. Genetic-correlation stability across the ancestry-specific Singapore panels.
7. Trait-specific diagnostics, with T2D separated.
8. External-cohort benchmark after a controlled rerun.

### 5. Discussion

1. Main finding: reference construction changes LDSC scale, while matched 1 cM panels yield similar h2.
2. Stable rg estimates across the Singapore ancestry panels.
3. Why higher h2 is not evidence of greater accuracy.
4. Appropriate uses of pooled cov-LDSC and ancestry-specific standard panels.
5. Strengths, including large local reference samples and complete chromosome-level logs.
6. Limitations, including pooled GWAS ancestry, smaller Indian and Malay panels, and unavailable individual-level data for public release.
7. Resource release and implications for non-European analyses.

### 6. Conclusion

Population-matched LD resources are valuable, but they must be evaluated under matched windows and SNP universes. The current evidence supports robustness of the ancestry-specific Singapore panels, not a universal performance advantage.

## Figure and table plan

### Main figures

1. Sample and variant QC flowchart for the pooled and ancestry-specific panels.
2. Mean LD score versus window size, stratified by covariate specification.
3. Paired h2 and SE(h2) comparisons across panels for each trait.
4. rg concordance across Chinese, Indian, and Malay panels, with confidence intervals.
5. Diagnostic plot of intercept, ratio, mean chi-square, and SNP overlap by trait and panel.

### Main tables

1. Reference-panel construction and QC summary.
2. Benchmark cohort and trait manifest.
3. Trait-level h2 results with SE, intercept, ratio, mean chi-square, and SNP count.
4. Pair-level rg results with SE, adjusted P value, and cross-panel stability.

### Supplementary material

1. Per-chromosome LD-score file manifest and checksums.
2. Window and covariate sensitivity results.
3. Complete external-cohort benchmark including failures.
4. Software versions and exact commands.
5. Data-access and permitted-use table.

## Proposed public repository layout

```text
Project4/
├── README.md
├── WORKFLOW.md
├── metadata/
│   ├── panel_manifest.tsv
│   ├── benchmark_manifest.tsv
│   └── software_versions.tsv
├── scripts/
│   ├── reference_panel/
│   ├── sumstats/
│   ├── ldsc/
│   └── summarise/
├── results/
│   ├── primary/
│   └── sensitivity/
└── provenance/
    ├── checksums/
    └── example_logs/
```

The repository should contain parameterised scripts, aggregate result tables, compact representative logs, checksums, and data-access instructions. Controlled genotype data, participant-level covariates, and provider-restricted GWAS files should remain outside GitHub.
