# Analysis lineage

## HELIOS cognition

```text
Historical six-task residualised PCA score
(stored orientation: higher score = poorer cognitive performance)
    -> refined Chinese-only REGENIE GWAS, phenotype g_raw
       (5700 Step 2 participants; modal per-variant N = 5592)
    -> HELIOS_10k_lifted_clean_with_stats_FINAL.txt.gz
       SHA-256 80d5047fc4b7f3e62723e913ff5e6d3f57780460c4e1b4fcf8dc1a4c73a452b7
    -> HELIOS_10k_for_LAVA.sumstats.gz
       SHA-256 d5f68ced2f301bd8cdf7cc7a402677d1877a0fe3e519982f62060cbb1b18dc6c
    -> EAS_MDD_Edu_SCZ_HELIOS.univ.lava
    -> EAS_MDD_Edu_SCZ_HELIOS.bivar.lava
```

The April 2026 global LDSC record used a separate HapMap3/build-37 export:

```text
10k_clean_g.rsid37.tsv.gz
    -> 10k_clean.sumstats.gz
       SHA-256 dd0a9a4f0aba9a19d7ccce8da623100dcefa1ff37db9977e239900b698f3cd22
    -> 10k_clean_vs_EDU_KOGES_1kgeas.log
```

The global and local inputs are not byte-identical. Across 809,717 allele-aligned variants their Z statistics correlate 0.9929, with 99.51% sign agreement and no allele mismatches. They represent the same refined Chinese-only phenotype signal but separate formatting stages. This is not the trans-ancestry 10k meta-analysis.

## KoGES educational attainment

```text
KoGES educational-attainment GWAS
    -> EDU_KOGES.sumstats.gz
       SHA-256 53fba8c2467c9af5b5f66bd68277bec5151563a4f378cb6ee0ade3107e007548
    -> Edu_EAS_for_LAVA.sumstats.gz
       SHA-256 0765f10ba451f68bb6b91e5d758137afeb6ecaca590fa676f6002444bee85196
    -> EAS_MDD_Edu_SCZ_HELIOS.univ.lava
    -> EAS_MDD_Edu_SCZ_HELIOS.bivar.lava
```

All 1,031,845 LDSC variants aligned to the LAVA export have effectively identical Z statistics (r > 0.999999999999; sign agreement 100%). The GWAS file lineage is therefore clear.

## Historical LAVA specification and correction

The input files are identifiable, but the analysis specification is not valid for the target question. `input.info_EAS.txt` declared `Edu_EAS` using `15,771` cases and `178,777` controls. The LAVA run log confirms that LAVA v0.1.5 treated EDU as binary. Inspection of the installed `process.locus()` implementation confirms that the case fraction entered `process.binary()` and changed the reconstructed marginal SNP correlations. Consequently, it changed local h2, eligibility, rho and P values.

The output contains every row expected under that specification (726/726
overall; 92 HELIOS–EDU). Completeness was not the problem; the phenotype-type
setting was.

The exact historical control and corrected quantitative-education run were
completed on 19 September 2026. The control reproduced all deterministic
outputs exactly. Because LAVA 0.1.5 uses unseeded Monte Carlo draws for
bivariate intervals and P values, those simulated fields were highly
concordant rather than byte-identical; inferential decisions were unchanged.

Correcting education changed the all-pair output from 726 to 717 bivariate
rows, but all 92 HELIOS–education regions remained eligible. Target-pair rho
correlation was 0.9997, no direction changed, one nominal decision changed and
no FDR decision changed. The corrected output is used for all final summaries
and figures.
