# Chapter 3, Study 2: from research question to code

Chapter 3 brings together the HELIOS cognitive GWAS (Study 1, in
[Project 3](../Project3/)) and the evaluation of Singapore LD reference panels
(Study 2, in this directory). This guide maps the LD-reference study to the
available code and aggregate evidence.

## Research question

Does a large local LD reference change SNP-heritability or genetic-correlation
estimates relative to a public reference when the analysis settings are matched?
The study separates ancestry composition, reference size and principal-component
(PC) adjustment, rather than treating a higher heritability estimate as evidence
of better performance.

## Chapter material and repository contents

| Chapter 3 component | Available material |
| --- | --- |
| Reference-panel construction | [Input preparation](scripts/prepare_controlled_reference_inputs.sh), [Singapore panel construction](scripts/build_controlled_ldscore_panels.sh), [exact 1000 Genomes EAS reconstruction](scripts/build_exact_1kgeas_panel.sh) |
| Primary controlled LDSC comparison | [Analysis runner](scripts/run_controlled_panel_ldsc.sh), [aggregate results](results/controlled_reference/), [design and validation records](validation/controlled_reference_experiment/) |
| Paired heritability and genetic-correlation comparisons | [Paired jackknife summary](scripts/summarise_paired_h2_jackknife.py), [genetic-correlation summary](scripts/summarise_controlled_rg.py), [final report](validation/controlled_reference_experiment/FINAL_REPORT.md) |
| Reference-size sensitivity | [Panel subsampling and construction](scripts/build_reference_size_panels.sh), [heritability analyses](scripts/run_reference_size_h2.sh), [LD-score stability results](results/controlled_reference/reference_size_ldscore_summary.tsv) |
| Held-out simulation | [Simulation runner](scripts/run_heldout_simulation.sh), [replicate results](results/controlled_reference/simulation_replicates.tsv) |
| Earlier ancestry-specific and HELIOS analyses | [Ancestry-matched results](results/ancestry_matched/), [recovered HELIOS comparisons](results/recovered/), [provenance audit](PROJECT4_AUDIT.md) |

## Current interpretation

The primary comparison uses five panels with the same GRCh38 genetic map,
1 cM window, 874,025 variants in identical order, regression weights and
200 jackknife blocks. It applies each panel to the same six Chinese SG100K
GWAS. The local Chinese reference contains 27,178 unrelated participants;
the exact 1000 Genomes East Asian comparator contains 504.

The two standard references gave similar heritability and genetic-correlation
estimates. Increasing reference size improved LD-score stability but did not
demonstrate better heritability precision for these GWAS. Adjustment for 20 PCs
had little effect within the Chinese reference and a larger effect in the pooled
reference. Five held-out simulation replicates provide a calibration check,
not a definitive comparison of panel accuracy.

Earlier pooled 20 cM analyses and ancestry-specific analyses used different
specifications. Their outputs remain available for provenance and should be
interpreted separately from the primary controlled comparison. Likewise,
comparisons among HELIOS phenotype or sample freezes should not be read as a
pure sample-size experiment.

## Running the work

Start with the [portable synthetic example](portable/README.md) to explore
validation, command planning and paired-comparison calculations using Python
3.9+ and its standard library. The example does not perform scientific LDSC
inference. The portable wrapper is a local prototype; its README records the
remaining validation needed for a publication release.

The [research workflow](WORKFLOW.md) describes execution with real inputs.
Rebuilding the Singapore panels requires approved access to individual-level
genotypes and covariates, and rerunning the benchmark requires the relevant
provider-controlled GWAS files. These inputs are not distributed here.

The repository supplies scripts, aggregate result tables and compact validation
records. The [extended LD-reference draft](CHAPTER_DRAFT.md) supplies additional
methods and interpretation; it is working text rather than a final thesis or
published manuscript.
