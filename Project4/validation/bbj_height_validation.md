# Biobank Japan height validation

## Question

The archived external benchmark showed very low heritability for nearly every file labelled as a Biobank Japan trait. This audit tested whether that reflected genuinely weak BBJ signal or an input/provenance error.

## Historical input failure

The archived download script labelled `GCST90018907` as BBJ height. The GWAS Catalog identifies that accession as **psoriasis vulgaris**, not height. The same problem affects every JBB label in the script; see [`cohort_accession_audit.tsv`](cohort_accession_audit.tsv).

Two further issues affected the historical benchmark:

1. the Singapore-panel formatter supplied `N = 160000` to JBB files irrespective of the phenotype-specific sample size;
2. the 1000 Genomes EAS formatter used `N = 400000` whenever a recognised N column was absent.

The archived JBB height-labelled run therefore analysed the wrong phenotype with the wrong sample size. Its low estimate is not evidence that BBJ height lacks polygenic signal.

## Correct source

The validation used:

```text
GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz
SHA-256: 345ad177a9f25884879b70fa7af6f43a7dd19f103a851db281507abf1c2535f6
N: 159095
effect: BETA
P value: P_BOLT_LMM_INF
effect allele: ALLELE1
other allele: ALLELE0
allele frequency: A1FREQ
INFO: INFO
```

This is the Japanese height GWAS described by Akiyama et al. (2019), DOI: [10.1038/s41467-019-12276-5](https://doi.org/10.1038/s41467-019-12276-5).

## Reproduction

The exact parameterised commands are in [`scripts/run_bbj_height_validation.sh`](scripts/run_bbj_height_validation.sh). LDSC v1.0.1 was used. Summary statistics were restricted to HapMap3 variants during munging and analysed with the complete 22-chromosome 1000 Genomes East Asian LD-score panel at:

```text
/home/ec2-user/input_files/reference/eas_ldscores/
```

The alternative `1kgeas_consensus_867k` and `1kgeas_consensus_930k` directories were not used because chromosome 10 contained only a header (34-byte compressed file). The selected reference contained 64,243 chromosome 10 variants plus its header.

## Results

| Diagnostic | Historical file labelled BBJ height | Correct BBJ height rerun |
| --- | ---: | ---: |
| Munged SNPs | 1,217,311 | 943,742 |
| SNPs in LDSC regression | 1,088,621 | 941,939 |
| Observed-scale h2 | 0.0067 | 0.3875 |
| SE(h2) | 0.0026 | 0.0188 |
| Mean chi-square | 1.0792 | 2.6432 |
| Intercept | 1.0230 | 1.2478 |
| SE(intercept) | 0.0104 | 0.0250 |
| Ratio | 0.2908 | 0.1508 |

The corrected munging log reports mean chi-square 2.644, lambda GC 1.767 and 9,572 genome-wide significant variants. These diagnostics closely reproduce the published mean chi-square of approximately 2.63 and intercept of approximately 1.21 for the source height GWAS.

## Conclusion

The suspicion was justified: the near-zero historical result arose from a mislabeled accession and incorrect sample-size handling. Genuine BBJ height summary statistics show strong polygenic signal and a precise LDSC heritability estimate. All rows derived from the archived UKB/JBB/CBK/KOGES_TWB accession map are excluded from substantive interpretation. They may be replaced only after each phenotype, cohort, genome build, effect scale and sample size has been independently verified.

