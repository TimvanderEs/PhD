# Baseline reproduction

The historical six-task PCA was reproduced in **7,403 participants**. PC1 explained **34.49%** of variance. The reconstructed stored-orientation score correlated **r = 1.000000000000** with the saved score; the largest absolute difference after scale alignment was **2.02e-12**.

All six loadings and task-PC1 correlations were within the pre-specified tolerance of the historical reference. The `g` column in the retained normalized-phenotype table also matched the separate saved-score file exactly. Baseline status: **PASS**.

The principal Project 3 genetic record is the Chinese–Indian–Malay meta-analysis. The later validated output contains 6,686,830 variants and reaches a maximum variant-level N of 7,789; the earlier retained iteration contains 6,688,950 variants and reaches N = 7,664. Both are identified by checksum in the trans-ancestry manifest.

Two related records are retained separately. The refined Chinese-only REGENIE freeze entered 5,700 participants in Steps 1 and 2 (modal per-variant N = 5,592). The LDSC record is labelled `historical_10k_clean`; it gave h2 = 0.3048 (SE 0.0949), intercept = 1.0257 (SE 0.0081), and rg with KoGES EDU = -0.0502 (SE 0.0835, P = 0.5476). These values are not presented as if they came from the same output file.

The participant-level REGENIE phenotype is not distributed in this repository. Reproduction of its construction is supported by the exact PCA match, while the final Chinese GWAS settings are supported by its Step 1 and Step 2 logs.
