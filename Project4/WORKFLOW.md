# Controlled-analysis workflow

The final analysis was run on EC2 because panel construction requires controlled individual-level genotypes. The public repository contains the scripts, aggregate outputs, and compact logs, but not those restricted inputs.

## Required inputs

- PLINK BED/BIM/FAM genotype data for the Singapore reference resource.
- Chinese, Indian, and Malay two-column PLINK keep files.
- Raw 1000 Genomes EAS BED/BIM/FAM genotypes.
- GRCh38 genetic maps for chromosomes 1–22.
- The rsID-to-Singapore marker crosswalk.
- Six munged Chinese SG100K GWAS files.
- PLINK 1.9, PLINK 2, Python 3, Python 2.7 with the LDSC environment, and the cov-LDSC fork.

Paths in the shell scripts are defaults from the archived EC2 run. They can be overridden with the environment variables declared at the top of each script. Set `ROOT` to a writable analysis directory containing `scripts/`, `input/`, and `maps/`.

## Run order

1. Prepare unrelated samples and 20-PC covariates.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   GENO=/controlled/sg100k_consensus \
   bash scripts/prepare_controlled_reference_inputs.sh
   ```

2. Build standard and 20-PC-adjusted Chinese and pooled 1 cM panels.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   bash scripts/build_controlled_ldscore_panels.sh
   ```

3. Recalculate 1KG EAS from raw genotypes and form the exact five-panel intersection.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   RAW_1KG=/controlled/g1000_eas \
   bash scripts/build_exact_1kgeas_panel.sh
   ```

4. Run the real-GWAS analyses and regenerate the aggregate summaries.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   SOURCE_ROOT=/controlled/SG100K_ancestry_matched \
   bash scripts/run_controlled_panel_ldsc.sh

   python3 scripts/summarise_paired_h2_jackknife.py /analysis/controlled_reference_experiment
   python3 scripts/summarise_controlled_rg.py /analysis/controlled_reference_experiment
   ```

5. Build the Chinese reference-size panels and analyse the six GWAS.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   bash scripts/build_reference_size_panels.sh

   ROOT=/analysis/controlled_reference_experiment \
   bash scripts/run_reference_size_h2.sh
   ```

6. Run the held-out calibration simulation.

   ```bash
   ROOT=/analysis/controlled_reference_experiment \
   bash scripts/run_heldout_simulation.sh
   ```

7. Validate the compact exported bundle.

   ```bash
   python3 scripts/validate_controlled_experiment.py \
     validation/controlled_reference_experiment \
     --reader-results results/controlled_reference \
     --report validation/controlled_reference_experiment/qc/local_validation.json
   ```

## Statistical summaries

The h2 comparison uses the full-sample LDSC contrast as its point estimate. The paired SE is calculated from differences between corresponding delete-one-block estimates:

`SE = sqrt((B - 1) / B * sum((d_i - mean(d))^2))`, with `B = 200`.

This retains the covariance between estimates generated from the same GWAS and avoids treating the two h2 estimates as independent. BH correction is applied across the six traits within each prespecified contrast.

The five held-out simulations are a calibration check. A performance study should use substantially more replicates, more than one causal architecture, and prespecified error and coverage criteria.
