# Analysis-design reconstruction

Verified against the EC2 archive on 21 September 2026; ancestry-matched and final controlled analyses completed on 22 September 2026.

## Conclusion

The recovered archive did not contain an ancestry-matched GWAS analysis for the Chinese, Indian, and Malay Singapore LD panels. Each panel had been applied to the same six fixed-effect trans-ancestry SG100K summary-statistic files.

A new validation used the Chinese, Indian, and Malay association columns from the provider files. All 63 matched analyses completed, together with 21 preliminary Chinese-GWAS/1KG-EAS comparator analyses. A subsequent controlled experiment recalculated the 1KG comparator from raw genotypes and matched five reference panels on window, genetic map, SNP universe, order, and regression weights. These are explicitly labelled as new validation work rather than recovered analyses.

## What was actually run

### Historical external benchmark

- Five dataset labels: UKB, Pan-UKB, JBB, CKB, and KoGES/TWB.
- Nine intended traits per dataset: ALT, BMI, height, LDL cholesterol, systolic blood pressure, smoking, type 2 diabetes, triglycerides, and waist-to-hip ratio.
- Three LD references: the pooled Singapore panel, 1000 Genomes EAS, and 1000 Genomes EUR.
- The 45 prepared files and five logs per reference are present.
- HELIOS/SG100K was not one of these five dataset labels.
- This benchmark is excluded because the archived accession-to-phenotype mapping and substituted sample sizes were invalid.

The five-dataset issue has not been resolved. Only genuine BioBank Japan height has been rebuilt and validated. Pan-UKB phenotype field IDs are plausible, but the recovered formatter selected the trans-ancestry `beta_meta_hq` and `se_meta_hq` fields rather than the EAS fields and substituted a fixed N of 420,531. The current status of every dataset is recorded in [`five_cohort_rebuild_status.tsv`](five_cohort_rebuild_status.tsv).

### Six-trait Singapore benchmark

- Six SG100K release 2.3.1 traits: systolic blood pressure, waist-to-hip ratio, height, BMI, LDL cholesterol, and type 2 diabetes.
- The provider files contain `FREQ`, `BETA`, `SE`, `P`, and `N` fields for Chinese, Indian, Malay, and trans-ancestry results.
- The recovered preparation script selects only `FREQ_trans`, `BETA_trans`, `SE_trans`, `P_trans`, and `N_trans`.
- The resulting six munged files have total sample sizes of approximately 46,200--46,300, rather than ancestry-specific sample sizes.
- The same six munged files were analysed with the Chinese, Indian, and Malay 1 cM LD panels. The same trans-ancestry files were also used in the earlier pooled-Singapore and 1000 Genomes EAS comparisons.

Only these six SG100K provider packages were recovered on EC2 and in the local project archives. SG100K packages for ALT, smoking, and triglycerides were not found. Files with those trait labels exist only in the excluded external-dataset archive. A nine-trait SG100K analysis therefore requires the three additional provider packages, or a separately specified and documented new GWAS for those phenotypes.

## Evidence on EC2

| Evidence | EC2 path | Interpretation |
| --- | --- | --- |
| Provider header | `/home/ec2-user/sg100K_sumstats/c8_bmi/c8_bmi_chr1.tsv` | Contains Chinese, Indian, Malay, and trans-ancestry result columns |
| Provider metadata | `/home/ec2-user/sg100K_sumstats/*/*.readme` | GWAS were run within ancestry and combined using fixed-effect METAL |
| Preparation script | `/home/ec2-user/sg100K_sumstats/prep_for_munge.sh` | Explicitly selects the five `*_trans` fields |
| Prepared inputs | `/home/ec2-user/sg100K_sumstats/munged/*.sumstats.gz` | Six trans-ancestry files only |
| Heritability log | `/home/ec2-user/project3/helios_ldsc_new_1cm/logs/nohup_h2_all.log` | Reuses each trans-ancestry file with all three ancestry LD panels |
| Genetic-correlation runner | `/home/ec2-user/project3/helios_ldsc_new_1cm/scripts/run_rg_all_helios_1cm.sh` | Reuses the same trans-ancestry inputs with all three ancestry LD panels |
| Result directories | `/home/ec2-user/project3/helios_ldsc_new_1cm/{h2,rg}/{Chinese,Indian,Malay}` | Directory names identify LD-panel ancestry, not GWAS ancestry |
| New matched validation | `/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched` | Separately prepared Chinese, Indian, and Malay GWAS columns and complete retained outputs |
| 1KG EAS identifier crosswalk | `/home/ec2-user/COV_LDSC_paper/rsid_to_marker.map` | Used to identify the raw 1KG genotype variants before recomputing the final comparator |
| Final controlled experiment | `/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment` | Exact raw-genotype 1KG EAS and four Singapore panels on 874,025 shared variants |

A targeted search of the pre-validation EC2 archive found no munged SG100K files, preparation scripts, or LDSC logs selecting `BETA_chinese`, `BETA_indian`, or `BETA_malay`.

## Correct interpretation

The recovered comparison tests the sensitivity of trans-ancestry SG100K LDSC estimates to LD-reference choice. The ancestry-matched validation adds:

1. Chinese GWAS with Singapore Chinese and 1000 Genomes EAS panels;
2. Indian GWAS with the Singapore Indian panel; and
3. Malay GWAS with the Singapore Malay panel.

The preliminary Chinese comparator translated and row-subset existing LD-score files. It is superseded for panel-performance claims by the final controlled experiment in [`controlled_reference_experiment/`](controlled_reference_experiment/), which recomputed 1KG EAS LD scores from raw genotypes. In that experiment, the standard Singapore Chinese and 1KG EAS panels gave mean h2 of 0.1618 and 0.1601, respectively; none of six paired differences was nominally significant, and rg estimates correlated at 0.99945. No suitable public South Asian comparator was available locally, so the Indian and Malay results cannot rank panel performance.
