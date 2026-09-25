# HELIOS ancestry and self-reported ethnicity analysis

## Data sources

This analysis separates three variables previously described under the same ancestry label:

- grids_meta.ancestry: GRID genetic ancestry assignment and GRID PC1–PC5;
- FREG5_Race: race recorded from the NRIC; and
- FIAQ10_Demo1: participant-described ethnic heritage from the interviewer-administered questionnaire.

Parental ethnic heritage comes from FIAQ14_1_Demo6 and FIAQ14_2_Demo7. The analysis includes **23,615** genetically annotated HELIOS participants. Self-described ethnic heritage is available for **13,906 (58.9%)**.

## Existing ancestry pipeline

The SG100K annotation provides seven genetic categories. The HELIOS counts are:

- Chinese: **18,905**
- Indian: **2,481**
- Malay: **1,367**
- Chinese–Malay: **386**
- Chinese–Indian: **179**
- Indian–Malay: **181**
- Other: **116**

These supplied genetic categories, rather than a new PCA threshold, are the primary ancestry classification.

## Self-reported ethnicity variables

FREG5_Race should be described as **NRIC-recorded race**, not self-reported ancestry. FIAQ10_Demo1 is the appropriate self-described ethnic-heritage variable. Among participants present in both sources, the annotation race and FREG5_Race agree exactly after removal of one conflicting duplicate core record.

Questionnaire availability differs substantially by NRIC-recorded race:

- Chinese: 12,463/19,046 (65.4%).
- Indian: 921/2,680 (34.4%).
- Malay: 461/1,773 (26.0%).
- Other: 61/116 (52.6%).

Comparisons using FIAQ10_Demo1 therefore describe the questionnaire subset and may be affected by differential availability.

## Genetic ancestry classification

The older script that produced PCA_inferred.png copied the recorded label for every participant with a non-missing label and inferred only missing values. It is not an independent genetic-classification result and should not be used as evidence of concordance.

## Admixed/intermediate ancestry

The GRID annotation directly identifies Chinese–Malay, Chinese–Indian and Indian–Malay genetic categories. In participants with both parental ethnic-heritage variables recorded:

- Chinese–Indian: **27/65 (41.5%)** reported one Chinese and one Indian parent.
- Chinese–Malay: **21/152 (13.8%)** reported one Chinese and one Malay parent.
- Indian–Malay: **14/61 (23.0%)** reported one Indian and one Malay parent.

The Chinese–Indian category has the strongest direct support from parental heritage. The lower matching-parent proportions for Chinese–Malay and Indian–Malay support more cautious language: the genetic categories may reflect recent mixture, older admixture, continuous population structure, or limitations of categorical parental responses.

## Self-report/genetic concordance

Within the questionnaire subset:

- Chinese: exact core match **98.87%**; pairwise intermediate involving Chinese **1.00%**; genetic label not involving Chinese **0.13%**.
- Indian: exact core match **92.99%**; pairwise intermediate involving Indian **5.37%**; genetic label not involving Indian **1.64%**.
- Malay: exact core match **75.84%**; pairwise intermediate involving Malay **17.46%**; genetic label not involving Malay **6.70%**.

## Direction of discordance

The apparent excess of Chinese discordance in the earlier PCA is not supported after using rates. Self-described Chinese participants have the **lowest**, not the highest, rate of a genetic label that does not involve their recorded group. The same ordering is present when NRIC-recorded race is used on the full annotated sample.

Rates of a strict non-involving genetic label differ across the three self-described groups (Pearson chi-square = **447.61**, df = 2, p = **6.35e-98**). The risk among self-described Chinese participants is **0.078 times** that among self-described Indian participants (95% CI 0.038–0.160; Fisher p = 1.01e-09) and **0.019 times** that among self-described Malay participants (95% CI 0.010–0.035; Fisher p = 2.13e-31). These comparisons are descriptive of the questionnaire subset because ethnic-heritage availability is differential.

## Chinese–Indian intermediate ancestry

The GRID Chinese–Indian category contains **179** HELIOS participants. Its parental-heritage match is higher than for the other two pairwise genetic categories. This is compatible with a recent mixed-ancestry contribution in a subset, but the genetic label should not be interpreted as proof of parental background for every participant.

## Chinese–Malay intermediate ancestry

The GRID Chinese–Malay category contains **386** participants. Because Chinese and Malay populations are genetically closer and only a minority of participants with parental data report one Chinese and one Malay parent, describe this category as **Chinese–Malay intermediate genetic ancestry** unless additional ancestry-proportion evidence is available.

## Indian–Malay intermediate ancestry

The GRID Indian–Malay category contains **181** participants. Within-Indian heterogeneity and incomplete parental data remain relevant limitations.

## GWAS inclusion consequences

The combined 22k retained analysis table assigned ancestry-specific labels only to GRID core Chinese, Indian and Malay participants. Intermediate/other participants assigned to an ancestry-specific label: **0**.

The retained 10k merged covariate table also maps its Chinese, Indian and Malay Ancestry values exclusively to the corresponding GRID core category. Thus the retained analysis Ancestry variable is a genetic-analysis grouping, not a self-report variable.

The pooled 22k phenotype PCA nevertheless contained 826 pairwise intermediate or Other GRID participants (3.68% of N=22,422). Because the retained phenotype regression used only Indian and Malay indicators, all 826 entered its Chinese reference coding; seven core GRID Indian or Malay participants were coded the same way. This affected phenotype construction but not membership of the ancestry-specific GWASs. A direct seven-category GRID sensitivity and a core-only refit showed that the coding had a negligible effect on PC1.

## Sensitivity analyses

- Nearest-centroid descriptions were repeated using PC1–PC2, PC1–PC3, PC1–PC5 and standardised PC1–PC5. With standardised PC1–PC5, the original core group remained the nearest centroid for **100.0%** of core Chinese, **99.6%** of core Indian and **99.8%** of core Malay participants. Allocation of intermediate groups changed with PC scaling, confirming that this descriptive procedure cannot reconstruct the original GRID classification rule.
- Relatedness/duplicate flags were summarised by genetic category. These flags do not establish whether a visible bridge is family-driven; a kinship matrix would be required for that test.
- The available SG100K covariate file contains subcohort values for core groups but not the pairwise intermediate categories, so a valid intermediate-group batch/subcohort comparison could not be completed from this file.
- Replacing the original phenotype ancestry dummies with all seven GRID categories produced loading congruence 0.99999994 and score correlation r=0.999851. Restricting phenotype construction to the core GRID groups produced congruence 0.99998224 and r=0.999882.
- Within the same 13,102-person questionnaire subset, adding participant-described ethnic heritage produced loading congruence 0.99999997 and r=0.999729. Excluding discordant records produced congruence 0.99999884 and r=0.999982. Full results are in [`phenotype_sensitivity`](phenotype_sensitivity/PHENOTYPE_ETHNICITY_SENSITIVITY.md).

## Interpretation

Self-described ethnic heritage, NRIC-recorded race and genetic ancestry are related but non-equivalent variables. The major genetic clusters align closely with both recorded variables, while GRID also identifies pairwise intermediate structure. The data do not support the hypothesis that self-described Chinese participants are disproportionately discordant after accounting for group size. Correcting the pooled phenotype ancestry coding and modelling self-described ethnicity separately did not materially change the cognitive factor.

## Limitations

- Self-described ethnic heritage is missing non-randomly across NRIC race groups.
- GRID ancestry proportions and the exact classification thresholds were not present in the retained local archive.
- Parental ethnic heritage uses broad categories and cannot distinguish all forms or generations of admixture.
- PCA and categorical ancestry assignments cannot identify why a participant uses a particular ethnic label.
- A formal batch analysis and family-cluster analysis require source variables not present for the intermediate groups in the retained files.
