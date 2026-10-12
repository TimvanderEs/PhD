# Project 1: Local genetic overlap of depression and cognition

This project examines the regional genetic relationships of **major depressive
disorder (MDD)** with **intelligence (INT)** and **educational attainment (EA)**.
It asks whether their negative genome-wide genetic correlations conceal
different patterns of sharing across the genome, and whether locally correlated
regions support shared or distinct association signals.

Using European-ancestry GWAS summary statistics, the study combines LD-score
regression, LAVA local genetic correlation, SuSiE-RSS fine-mapping and
colocalisation. It forms the basis of Chapter 2 of my PhD. This overview reflects
the preprint draft dated 7 October 2026.

## Main findings

- The primary analysis identified **19 FDR-significant region–trait-pair
  comparisons**: 16 MDD–EA results (12 negative and 4 positive) and 3 negative
  MDD–INT results, after correction across 71 eligible comparisons.
- SuSiE-based colocalisation supported a shared MDD–INT signal in the major
  histocompatibility complex (MHC; PP.H4 = 0.987), but distinct MDD and EA signals
  in the same interval (PP.H3 = 0.997).
- A secondary single-signal analysis supported another MDD–INT signal on
  chromosome 9. The MDD–EA chromosome 5 result remained inconclusive.

The findings show regional heterogeneity and reinforce that EA and intelligence
are related but distinct phenotypes. Local genetic correlation does not by itself
establish colocalisation or causal direction. The complex LD structure of the
MHC and the European-only analysis warrant cautious interpretation and replication.

## Analysis software

The selected-locus fine-mapping and colocalisation analyses used
[COLOC-reporter, maintained by Thomas Spargo](https://github.com/ThomasPSpargo/COLOC-reporter/tree/master).
That repository provides the analysis scripts, environment specification,
input preparation guidance, fine-mapping and colocalisation options, and
summary-statistic plotting tools.

This directory provides the study overview; COLOC-reporter hosts the
selected-locus implementation and setup instructions. Study-specific scripts,
complete outputs and the manuscript are not currently distributed here. Source
GWAS data remain subject to the original studies' access conditions.
