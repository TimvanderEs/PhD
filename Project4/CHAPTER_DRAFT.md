# Population-specific LD reference panels for multi-ancestry Singapore genomic analyses

## Draft status

This evidence-constrained working draft replaces the unsupported performance claims in the earlier manuscript with conclusions traced to scripts, logs, and result tables. The merge log identifies the 46,508-person genotype resource as PRECISE-SG100K-derived rather than HELIOS-only.

## Abstract

Linkage disequilibrium score regression (LDSC) depends on LD estimates that adequately represent the ancestry and variant structure of the genome-wide association study under analysis. Public panels are convenient, but broad continental groupings may not fully represent genetically heterogeneous populations such as Singapore, where people of Chinese, Malay, and Indian ancestry live within the same national setting. We developed LD-score resources from a large multi-ancestry Singapore genomic dataset and evaluated how reference construction affected SNP-heritability and genetic-correlation estimates.

The earlier work comprised a pooled 20 cM cov-LDSC panel and standard 1 cM panels calculated separately in Chinese, Indian, and Malay participants. Six SG100K release 2.3.1 traits were evaluated: systolic blood pressure, waist-to-hip ratio, height, body mass index, LDL cholesterol, and type 2 diabetes. A final controlled experiment analysed the Chinese GWAS with five 1 cM panels: exact 1000 Genomes East Asian LD scores, standard and 20-PC-adjusted Singapore Chinese scores, and standard and 20-PC-adjusted pooled Singapore scores. All panels used the same GRCh38 map, 874,025 variants in the same order, 1000 Genomes regression weights, and 200 jackknife blocks.

In the final controlled experiment, mean SNP-heritability was 0.1618 with the standard Singapore Chinese panel and 0.1601 with exact 1000 Genomes EAS; mean SE was 0.02035 and 0.02022. None of the six paired h2 differences was nominally significant. Genetic-correlation estimates were nearly identical across the two panels (r = 0.99945; mean absolute difference 0.00736; 15/15 directions concordant). Increasing the Chinese reference from 500 to 27,178 participants improved LD-score stability but had little effect on h2 or SE. Five held-out simulated traits also did not distinguish the panels. PC adjustment had little effect within Chinese participants but materially shifted estimates in the pooled multi-ancestry panel.

These results show that reference construction can materially affect LDSC estimates in pooled multi-ancestry data, but they do not show that reference size or ancestry matching alone improves Chinese h2 accuracy or precision over 1000 Genomes EAS. The value of the local resource lies in population representation, robustness checks, and coverage of groups without suitable public panels rather than a universal increase in heritability or power.

## 1. Chapter overview

LDSC is widely used to estimate SNP heritability, quantify genetic correlations, and distinguish polygenic signal from confounding in GWAS summary statistics. Its regression model relies on LD scores obtained from a reference sample. When the reference does not approximate the ancestry composition, allele-frequency spectrum, genome build, and variant universe of the target GWAS, the resulting estimates may be shifted or less precise.

This chapter evaluates locally derived LD-score panels for genetic studies in Singapore. It has two linked aims: to document the construction of the panels and to test whether their use changes LDSC results relative to a commonly used 1000 Genomes East Asian reference. The chapter also documents an important methodological lesson from the development process: a panel comparison is interpretable only when window size, SNP universe, genome build, and analytical model are held constant.

The work complements the preceding cognitive GWAS chapter. That analysis illustrates the practical need for ancestry-aware genomic resources; this chapter examines one such resource directly. It does not assume that a larger heritability estimate is more accurate. Instead, panel performance is evaluated using estimate stability, standard errors, intercepts, SNP coverage, model completion, and cross-panel concordance.

## 2. Introduction

### 2.1 LD-score regression and reference LD

Under the LDSC framework, variants with higher LD scores tag more neighbouring variation and are therefore expected to have larger association statistics for a polygenic trait. Regressing GWAS association statistics on LD scores provides estimates of SNP heritability and, in the bivariate setting, genetic correlation. The reference LD scores are not merely a technical input: they define the predictor used in the regression and influence the scale of the estimated slope.

Reference mismatch can arise through differences in ancestry, allele frequencies, LD structure, variant filters, genome build, or imputation quality. It can also arise when panels are calculated with different LD windows. A fair comparison must therefore distinguish the effect of population matching from the effect of panel construction.

### 2.2 Singapore as a multi-ancestry setting

Singapore includes substantial Chinese, Malay, and Indian ancestry groups. These groups broadly connect to East, South-East, and South Asian genetic diversity, but continental labels do not capture all local structure. A single public East Asian panel is not an obvious match for Indian-ancestry participants and may also incompletely represent Malay or locally admixed haplotypes.

The genotype resource initially contained 46,508 participants. This exceeds the published 10,004-participant HELIOS cohort, and the original merge log records the prefix `precise_sg100k`. This chapter therefore uses *PRECISE-SG100K-derived reference resource* rather than HELIOS-only.

### 2.3 Covariate-adjusted and ancestry-specific LD scores

Covariate-adjusted LDSC was developed to estimate LD after accounting for ancestry-related covariates in admixed or structured samples. This provides one possible strategy for constructing a pooled reference. A second strategy is to calculate standard LD scores within more homogeneous ancestry strata. These approaches target related but not identical reference representations and should be evaluated separately.

The initial analysis used a pooled Singapore panel, principal-component adjustment, and a 20 cM window. A later analysis constructed Chinese, Indian, and Malay panels without covariate adjustment using a 1 cM window. The latter provides a more direct comparison with standard LDSC reference panels, although it does not test the added value of covariate adjustment.

### 2.4 Study aims

The chapter addresses five questions:

1. Can complete ancestry-specific Singapore LD-score panels be constructed after relatedness and variant filtering?
2. How similar are SNP-heritability estimates from the Singapore 1 cM panels and 1000 Genomes East Asian LD scores?
3. How stable are genetic-correlation estimates across the Chinese, Indian, and Malay Singapore panels?
4. What is observed when ancestry-specific GWAS are paired with their corresponding Singapore panels?
5. Which additional controlled analyses are required to determine whether covariate-adjusted pooled LD scores improve inference?

## 3. Methods

### 3.1 Study design

The evaluation was reconstructed from code, execution logs, aggregate result tables, and data-provider README files. Historical analyses were retained for provenance, followed by an exact controlled experiment in which window, map, SNP universe, order, regression weights, and GWAS inputs were fixed across panels.

### 3.2 Reference samples and relatedness filtering

The pooled genotype merge log contained 46,508 participants. A June 2025 pooled cov-LDSC workflow retained 37,364 unrelated participants. A later documented ancestry-specific workflow retained 27,181 Chinese, 5,848 Indian, and 4,327 Malay participants. The final controlled rebuild started from 33,392 Chinese, 7,450 Indian, and 5,666 Malay participants and retained 27,178 Chinese and 37,351 pooled participants after KING filtering. Counts are reported by workflow rather than silently reconciled.

The eight-person difference between the pooled and summed ancestry-specific unrelated samples reflects different documented workflows and should not be silently reconciled. The pooled analysis is reported with N = 37,364; the ancestry-specific panels are reported with their individual counts.

### 3.3 Variant quality control

For the ancestry-specific panels, the recorded construction commands retained biallelic SNPs with minor-allele frequency at least 0.01, genotype missingness at most 0.02, and Hardy-Weinberg equilibrium P at least 1 x 10^-6. The final controlled rebuild applied those SNP, frequency, missingness, and Hardy-Weinberg filters within Chinese participants; Hardy-Weinberg filtering was not applied to the pooled multi-ancestry sample to avoid a Wahlund-effect filter. Variants were aligned to a GRCh38 genetic map. One marker requested through the public-panel crosswalk was absent from the raw 1KG genotype release, leaving an exact final intersection of 874,025 variants across all five panels.

### 3.4 Pooled covariate-adjusted panel

The exploratory pooled panel used the cov-LDSC implementation and a 20 cM window. The runner supplied 50 principal components as covariates; it did not supply age or sex. Statements in the earlier manuscript that age and sex were included are therefore unsupported by the command. Genome-build labels in the archived documentation are inconsistent, so this historical panel is retained as development provenance and is not used for the primary controlled comparison.

### 3.5 Ancestry-specific standard panels

The later workflow calculated LD scores separately in the Chinese, Indian, and Malay samples using a 1 cM window and no covariate file. Although the command used the `ldsc.py` file from a cov-LDSC installation, this invocation corresponds to standard, unadjusted LD-score calculation. Complete `.l2.ldscore.gz`, `.l2.M`, and `.l2.M_5_50` outputs were available for all 22 autosomes in each group.

The Chinese, Indian, and Malay panels contained 878,093, 912,313, and 909,907 LD-score variants, respectively. Their overlaps with the munged benchmark summary statistics were 827,024, 859,224, and 856,937 variants.

### 3.6 Benchmark GWAS

The controlled-window comparison used six SG100K genomic data release 2.3.1 GWAS. Quantitative traits were analysed within ancestry using REGENIE v3.4 after the transformations stated by the provider, with age, sex, cohort, subcohort, and the top 20 ancestry-specific principal components as covariates. Ancestry-specific results were combined by fixed-effect meta-analysis using METAL. The traits were average systolic blood pressure, waist-to-hip ratio, height, body mass index, log10-transformed LDL cholesterol, and type 2 diabetes.

Although the provider files retain separate Chinese, Indian, and Malay association columns, the LDSC preparation script selected only the fixed-effect trans-ancestry columns. The same six trans-ancestry files were used with each ancestry-specific LD panel. These analyses therefore assess sensitivity to LD-panel choice; they are not ancestry-matched GWAS analyses.

The five quantitative-trait sample sizes ranged from 46,216 to 46,314. The type 2 diabetes GWAS included 5,486 cases and 40,856 controls. These are SG100K trans-ancestry results and must not be described as HELIOS-only GWAS.

### 3.7 LDSC analyses and comparison metrics

Univariate LDSC was run with each of the three Singapore 1 cM panels and with the 1000 Genomes East Asian comparator. For each trait and panel, the evaluation should retain SNP heritability, its standard error, the LDSC intercept and standard error, mean chi-square, ratio where defined, SNP count, and completion status.

Bivariate LDSC was run for all 15 pairwise combinations of the six traits with each Singapore panel, producing 45 completed runs. The summary compared point estimates across panels using Pearson correlation and agreement in direction. These documented trans-ancestry results are described as nominal where multiple-testing correction was not retained.

Panel usefulness was not defined by the absolute magnitude of h2. The primary criteria were agreement under matched analytical conditions, precision conditional on calibration, stability across panels, and transparent retention of failed or implausible results.

### 3.8 Ancestry-matched and final controlled validation

The Chinese, Indian, and Malay association columns were extracted separately from the six provider packages, aligned to a common allele key, restricted to the 867,103-variant consensus set, and munged using the provider-supplied variant-level N. Quantitative-trait N was approximately 33,200 Chinese, 7,400 Indian, and 5,600 Malay participants. T2D included 3,037 cases and 30,214 controls in Chinese, 1,546 and 5,892 in Indian, and 903 and 4,750 in Malay participants.

Each ancestry-specific GWAS was first paired with its corresponding 1 cM Singapore panel. Six h2 and 15 rg analyses were run per ancestry. This stage documented the available ancestry-specific signal, but its preliminary Chinese comparator translated and row-subset existing 1KG LD-score files.

The final experiment recalculated 1KG EAS LD scores from 504 raw-genotype samples and regenerated Singapore Chinese and pooled LD scores, with and without 20-PC adjustment. The five panels shared the same 1 cM window, genetic map, 874,025 variants and order, and 1KG EAS regression weights. Thirty h2 and 75 rg analyses were completed. Full-precision delete-one-block values were retained to estimate paired SEs for h2 differences; Benjamini--Hochberg correction was applied across the six traits within each contrast. T2D was retained on the observed scale because population prevalence was not prespecified.

Reference-size sensitivity panels were generated from reproducible Chinese subsamples of 500, 2,500, and 10,000 participants and compared with the full 27,178-person panel. A held-out calibration simulation split the unrelated Chinese sample into 10,000 reference participants and 17,178 non-overlapping GWAS participants. Five quantitative traits with realized h2 = 0.20 were analysed with exact 1KG EAS, an overlapping full Chinese panel, and the disjoint 10,000-person Chinese panel.

### 3.9 Excluded external-cohort benchmark

An earlier benchmark evaluated 45 intended trait-dataset combinations using the pooled Singapore 20 cM panel and 1000 Genomes panels calculated with a 1 cM window. A subsequent provenance audit found that the download script had assigned intended phenotype labels to unrelated GWAS Catalog accessions. The formatters also supplied fixed sample sizes when a recognised sample-size field was absent. These failures invalidate the phenotype identities and LDSC estimates. The accession-derived UKB, Pan-UKB, JBB, CKB and KoGES/TWB rows were therefore excluded from all quantitative comparisons. They are retained only as an audit trail and were not used to motivate or support the six-trait Singapore analysis.

### 3.10 Reproducibility

The code archive includes parameterised panel-construction scripts, exact 1KG recomputation, summary-statistic harmonisation, univariate and bivariate LDSC runners, reference-size and simulation workflows, aggregate result tables, software versions, compact logs, and checksums. Individual-level genotypes, covariates, sample identifiers, and provider-controlled GWAS files are not deposited publicly.

## 4. Results

### 4.1 Reference-panel construction

The ancestry-specific workflow produced complete autosomal LD-score panels for 27,181 Chinese, 5,848 Indian, and 4,327 Malay participants. SNP coverage was broadly similar across the three panels, with 878,093--912,313 LD-score variants and 827,024--859,224 variants overlapping the benchmark summary statistics. The smaller Indian and Malay sample sizes remain relevant when interpreting precision and the stability of local LD estimates.

### 4.2 Trans-ancestry heritability under aligned 1 cM specifications

Across the six traits, mean h2 was 0.127 with 1000 Genomes East Asian LD scores, compared with 0.125 for the Singapore Chinese panel, 0.129 for the Singapore Indian panel, and 0.124 for the Singapore Malay panel. Relative differences from the public comparator were -1.7%, +1.4%, and -2.9%, respectively. These differences were small at the aggregate level.

Mean standard errors were 0.0129, 0.0134, 0.0156, and 0.0136 for the 1000 Genomes East Asian, Singapore Chinese, Indian, and Malay panels. Corresponding mean h2/SE ratios were 8.97, 8.50, 7.54, and 8.30. The controlled-window results therefore do not show a general precision advantage for the Singapore panels.

Mean intercepts were close to one for all four panel choices. Type 2 diabetes was an exception: its mean chi-square was below one and its intercept was approximately 0.897 across the 1 cM analyses. This panel-invariant pattern points to the GWAS input or statistical scale rather than a particular LD reference and should be investigated before the type 2 diabetes results are used as confirmatory evidence.

### 4.3 Genetic-correlation stability

Genetic-correlation estimates were highly concordant across the ancestry-specific Singapore panels. Pearson correlations were 0.984 for Chinese versus Indian, 0.997 for Chinese versus Malay, and 0.985 for Indian versus Malay. Fourteen of 15 trait-pair estimates had the same direction for Chinese versus Indian and Indian versus Malay; all 15 agreed in direction for Chinese versus Malay.

Six pairs were nominally significant in all three panels: systolic blood pressure with body mass index, systolic blood pressure with type 2 diabetes, type 2 diabetes with body mass index, waist-to-hip ratio with body mass index, waist-to-hip ratio with systolic blood pressure, and waist-to-hip ratio with type 2 diabetes. These labels remain descriptive until the multiple-testing procedure is specified.

### 4.4 Final controlled Chinese comparison

For the Chinese GWAS, mean observed-scale h2 was 0.1618 with the standard Singapore Chinese panel and 0.1601 with exact 1KG EAS. Mean SE was 0.02035 and 0.02022, and mean h2/SE was 7.29 and 7.23. Mean intercepts were 1.00525 and 1.00720. Trait-level h2 differences ranged from 0.0002 to 0.0037. Covariance-aware paired jackknife tests found no nominally significant difference across the six traits, and all BH-adjusted P values exceeded 0.32.

The Chinese rg estimates were highly stable: Pearson r was 0.99945 across the 15 pairs, the mean absolute difference was 0.00736, the maximum absolute difference was 0.0197, and all directions agreed.

The PC-adjusted Chinese panel yielded a similar mean h2 of 0.1620. Relative to the standard Chinese panel, mean h2 changed by 0.00022 and no trait survived BH correction. PC adjustment therefore had little practical effect within the Chinese-only reference.

### 4.5 Pooled-panel and reference-size sensitivity

Under the final 1 cM design, mean h2 was 0.1475 with standard pooled Singapore LD scores and 0.1705 after adjustment for 20 PCs. Five of six within-pooled differences survived BH correction. The same adjustment had little effect in the Chinese-only panel, indicating that ancestry-related covariance materially changes pooled multi-ancestry LD scores. The result identifies specification sensitivity; it does not determine which pooled model is more accurate.

Increasing the Chinese reference size improved LD-score stability but did not materially improve the real-GWAS estimates. Mean h2 was 0.1601, 0.1614, 0.1618, and 0.1618 at N = 500, 2,500, 10,000, and 27,178, respectively; corresponding mean SEs were 0.02022, 0.02033, 0.02035, and 0.02035. Across chromosomes 1, 6, and 22, mean correlation with full-panel LD scores increased from 0.99546 at N = 500 to 0.99986 at N = 10,000, while RMSE fell from 1.252 to 0.217.

In five held-out simulations with realized h2 = 0.20, mean estimated h2 was 0.2489 with 1KG EAS and 0.2464 with both the full overlapping and disjoint 10,000-person Chinese panels. Mean SE was 0.0319 for all three panels and 95% coverage was 0.60. The common upward bias and low coverage indicate that this small simulation series is useful for checking the workflow but is not sufficient to rank reference panels.

### 4.6 External-input audit and BBJ validation

The external-input audit showed that the low estimates labelled as Biobank Japan results did not arise from the intended phenotypes. The file labelled as BBJ height was `GCST90018907`, which the GWAS Catalog identifies as psoriasis vulgaris, and it had been formatted with a substituted sample size. Equivalent label failures affected the remaining JBB rows and most of the other accession-derived cohorts.

To test whether genuine BBJ height data lacked signal, the Japanese height GWAS of Akiyama et al. was munged with N = 159,095 and analysed using LDSC v1.0.1 and a complete 1000 Genomes East Asian reference. The corrected run retained 941,939 regression variants and estimated observed-scale h2 = 0.3875 (SE = 0.0188), with mean chi-square = 2.6432 and intercept = 1.2478 (SE = 0.0250). The archived file labelled as BBJ height had produced h2 = 0.0067 (SE = 0.0026). The corrected diagnostics reproduce the strong polygenic signal reported for the source GWAS and demonstrate that the historical near-zero value was an input-provenance failure.

No quantitative conclusion is drawn from the invalid external benchmark. Any future external comparison must begin with a new, independently verified cohort–trait manifest.

## 5. Discussion

### 5.1 Principal findings

This reconstruction produced four main findings. First, complete Singapore LD-score panels can be generated at substantially larger sample sizes than public ancestry-specific references commonly used in LDSC. Second, under an exactly matched 1 cM design, standard Singapore Chinese and 1KG EAS panels produced nearly identical h2 and rg estimates. Third, reference-size increases beyond 500 participants improved LD-score stability but did not reduce GWAS-level SEs in these traits. Fourth, PC adjustment had little effect within Chinese participants but materially altered the pooled multi-ancestry panel.

The evidence therefore supports comparability and robustness rather than universal superiority. The local panels reproduce the public-reference results for Chinese GWAS while providing a framework for populations and ancestry mixtures that are poorly represented by existing resources.

### 5.2 Interpreting population matching

Population matching is expected to matter when LD and allele-frequency patterns differ between reference and target data. In this Chinese comparison, however, 1KG EAS was sufficiently well matched that a 54-fold larger local panel did not materially change estimates or precision. Reference size was therefore not the limiting component of uncertainty for these GWAS and variants.

This conclusion should not be generalized to Indian, Malay, admixed, or trans-ancestry GWAS. Only the Chinese analysis had an appropriate public comparator available locally. Prespecified South Asian and South-East Asian comparators, or methods designed explicitly for cross-ancestry LD, are still required for those settings.

### 5.3 Pooled covariate-adjusted LD scores

The pooled cov-LDSC approach remains scientifically relevant because it attempts to represent structure within a heterogeneous sample. The final experiment isolated the effect of 20-PC adjustment at 1 cM and showed that it materially shifts the pooled h2 scale. This sensitivity was not present within the Chinese-only panel.

The direction of that shift does not identify the more accurate pooled model. Additional simulations representing the Singapore ancestry mixture, and ideally independent individual-level heritability estimates, are needed for calibration. Neither 1 cM nor 20 cM should be treated as universally correct without reference to the target estimand.

### 5.4 Strengths and limitations

Strengths include the large local reference sample, representation of three major Singapore ancestry groups, an exact raw-genotype 1KG comparator, identical variant order across five panels, paired 200-block uncertainty estimates, reference-size sensitivity, and held-out simulation. The reconstruction also prevents an analytical development artifact from becoming an overstated biological conclusion.

Several limitations remain. The retained merge log identifies the resource as PRECISE-SG100K-derived, but its exact data release and contributing-cohort composition still require confirmation. Only six provider packages were available; ALT, smoking, and triglycerides were not available. T2D estimates remain on the observed scale. The Indian and Malay GWAS were too weak for stable rg inference, and no suitable public comparators were available for them. The simulation used only five replicates and showed panel-invariant upward bias, so it cannot establish small performance differences. Finally, the early external benchmark is invalid and the historical 20 cM comparison remains confounded by multiple design changes.

### 5.5 Implications and next steps

The immediate priorities are to expand the simulation replicates and architectures, add public South Asian and South-East Asian comparators for Indian and Malay GWAS, and evaluate a prespecified cross-ancestry strategy for pooled summary statistics. If locally derived panels offer an advantage, it should appear consistently across calibration, precision, SNP overlap, and completion, not only as a larger mean h2.

Once validated, the Singapore LD-score panels could provide a useful resource for studies of Asian populations that are poorly represented by current public references. Their release must define the appropriate ancestry contexts, genome build, variant universe, construction parameters, and limitations.

## 6. Conclusion

Population-specific LD references are important for representation and for analyses without an appropriate public panel, but greater reference size does not automatically improve LDSC performance. Under the final controlled 1 cM design, standard Singapore Chinese and 1KG EAS panels produced effectively equivalent h2 and rg estimates. Covariate adjustment mattered mainly in the pooled multi-ancestry reference. The chapter should therefore emphasize reference construction, robustness, and population coverage rather than claim higher heritability, precision, or discovery power.

## References to verify and format

- Bulik-Sullivan BK et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics*. 2015;47:291--295. https://doi.org/10.1038/ng.3211
- Bulik-Sullivan B et al. An atlas of genetic correlations across human diseases and traits. *Nature Genetics*. 2015;47:1236--1241. https://doi.org/10.1038/ng.3406
- Luo Y et al. Estimating heritability and its enrichment in tissue-specific gene sets in admixed populations. *Human Molecular Genetics*. 2021;30:1521--1534. https://doi.org/10.1093/hmg/ddab130
- The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015;526:68--74. https://doi.org/10.1038/nature15393
- Khoo CM et al. The Health for Life in Singapore (HELIOS) Study: delivering precision medicine research for Asian populations. *Nature Communications*. 2025. https://doi.org/10.1038/s41467-025-65774-0
