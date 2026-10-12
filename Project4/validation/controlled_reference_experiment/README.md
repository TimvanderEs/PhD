# Controlled Singapore LD-reference experiment

This directory contains the compact outputs and quality-control records for the final controlled comparison completed on 22 September 2026. Individual-level genotype data, sample identifiers, covariates, and provider-controlled GWAS files are not included.

## Design

Six Chinese-ancestry SG100K GWAS were analysed: systolic blood pressure, waist-to-hip ratio, height, body mass index, LDL cholesterol, and type 2 diabetes. Five LD-score panels were calculated with the same 1 cM window, GRCh38 genetic map, final variant set, variant order, regression-weight panel, and 200-block jackknife:

| Panel | Reference participants | LD-score model |
| --- | ---: | --- |
| 1000 Genomes EAS | 504 | Standard |
| Singapore Chinese | 27,178 unrelated | Standard |
| Singapore Chinese | 27,178 unrelated | Adjusted for 20 PCs |
| Singapore pooled | 37,351 unrelated | Standard |
| Singapore pooled | 37,351 unrelated | Adjusted for 20 PCs |

The exact 1000 Genomes comparator was recalculated from raw genotypes. One requested chromosome 5 marker was absent from that genotype release, leaving 874,025 variants shared in identical order across all five panels.

The workflow completed 30 univariate and 75 bivariate real-GWAS LDSC runs. It also evaluated Chinese reference sizes of 500, 2,500, 10,000, and 27,178 participants, and ran five held-out simulated traits with target h2 = 0.20.

## Main results

The standard Singapore Chinese and exact 1000 Genomes EAS panels gave very similar mean h2 estimates (0.1618 and 0.1601) and mean SEs (0.02035 and 0.02022). Trait-level differences ranged from 0.0002 to 0.0037. Covariance-aware paired jackknife tests found no nominally significant difference across the six traits; all BH-adjusted P values were greater than 0.32.

Genetic-correlation estimates were also stable. Across the 15 trait pairs, the Singapore Chinese and 1000 Genomes estimates correlated at r = 0.99945, had a mean absolute difference of 0.00736, and agreed in direction for all 15 pairs.

Increasing the Chinese reference from 500 to 27,178 participants improved LD-score agreement with the full panel, but had little effect on the six real-GWAS estimates. Mean h2 was 0.1601 at n = 500 and 0.1618 in the full panel; mean SE was 0.02022 and 0.02035, respectively. On chromosomes 1, 6, and 22, mean LD-score correlation with the full panel increased from 0.99546 at n = 500 to 0.99986 at n = 10,000, while mean RMSE decreased from 1.252 to 0.217.

The five simulated traits did not distinguish the panels. Mean estimated h2 was 0.2489 with 1000 Genomes EAS and 0.2464 with both Singapore panels. The simulation is a calibration check, not a definitive performance comparison: five replicates are insufficient to establish small differences in bias or coverage, and the common upward bias occurred with every panel.

The pooled panel was sensitive to PC adjustment. Mean h2 increased from 0.1475 without adjustment to 0.1705 with 20 PCs, with five of six paired differences surviving BH correction. In contrast, PC adjustment changed mean h2 by only 0.00022 in the Chinese-only reference and no trait survived BH correction. This indicates that covariate adjustment matters mainly when ancestries are pooled; it does not demonstrate that either pooled specification is more accurate.

## Interpretation

The controlled experiment does not show that a larger, ancestry-matched Singapore reference improves LDSC precision or accuracy over 1000 Genomes EAS for these six Chinese GWAS. It shows that the two standard 1 cM references produce nearly equivalent h2 and rg estimates once the analysis is genuinely matched. The local reference remains useful for population coverage and for analyses without a suitable public comparator, but larger h2 estimates alone should not be described as improved performance.

This comparison concerns post-GWAS LDSC. REGENIE estimates marker effects from the study genotypes and does not use an external 1KG or Singapore LD-score panel in its association model. A Singapore LDSC panel is therefore not an alternative LD reference for the original REGENIE GWAS.

## Files

- `results/controlled_h2.tsv`: trait-level univariate estimates.
- `results/controlled_h2_panel_summary.tsv`: panel-level h2 summary.
- `results/paired_h2_differences.tsv`: paired jackknife contrasts.
- `results/controlled_rg.tsv`: all 75 genetic-correlation estimates.
- `results/controlled_rg_concordance.tsv`: concordance with exact 1000 Genomes EAS.
- `results/reference_size_*.tsv`: reference-size sensitivity.
- `results/simulation_*.tsv`: held-out simulation results.
- `qc/`: completion status, panel matching, software information, and checksums from the EC2 run.
- `logs/`: compact execution and LDSC logs.
- [`../../scripts/`](../../scripts/): reproducible preparation, analysis, and summarisation scripts.

The paired h2 summary uses the full LDSC estimate as the point estimate and paired delete-one-block values to estimate its SE. This avoids multiplying the four-decimal rounding of the printed LDSC estimate by the 200 jackknife blocks.

`qc/ec2_final_result_checksums.sha256` and `qc/ec2_final_script_checksums.sha256` preserve the absolute-path checksum manifests written on EC2 before export. `CHECKSUMS.sha256` validates the compact local bundle after the corrected summaries were regenerated.
