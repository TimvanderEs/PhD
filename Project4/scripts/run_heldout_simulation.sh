#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
GENO="${GENO:-${ROOT}/input/sg100k_consensus}"
PLINK="${PLINK:-/home/ec2-user/software/plink/plink}"
PLINK2="${PLINK2:-/home/ec2-user/software/plink/plink2}"
LDSC_PYTHON="${LDSC_PYTHON:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
LDSC="${LDSC:-/home/ec2-user/software/cov-ldsc/ldsc.py}"
REPLICATES="${REPLICATES:-5}"
TARGET_H2="${TARGET_H2:-0.2}"
CAUSAL_SNPS="${CAUSAL_SNPS:-10000}"

SIM_ROOT="${ROOT}/simulation"
KEEP_ROOT="${ROOT}/input/downsamples"
PANEL_ROOT="${ROOT}/panels/simulation_matched"
mkdir -p "${SIM_ROOT}" "${PANEL_ROOT}" "${ROOT}/logs/simulation" "${ROOT}/qc"
test -s "${KEEP_ROOT}/simulation_reference.keep"
test -s "${KEEP_ROOT}/simulation_gwas.keep"

python3 "${ROOT}/scripts/match_ldscore_panels.py" \
  --panel "1kgeas=${ROOT}/panels/1kgeas_matched" \
  --panel "chinese_full_overlap=${ROOT}/panels/chinese_standard" \
  --panel "chinese_disjoint_n10000=${ROOT}/panels/reference_size_raw/chinese_n10000" \
  --output-root "${PANEL_ROOT}" \
  --report "${ROOT}/qc/simulation_panel_matching.tsv"

common_snps="${SIM_ROOT}/simulation_common.snplist"
: > "${common_snps}"
for chromosome in $(seq 1 22); do
  zcat "${PANEL_ROOT}/1kgeas/${chromosome}.l2.ldscore.gz" | awk 'NR > 1 {print $2}' >> "${common_snps}"
done
sort -u "${common_snps}" -o "${common_snps}"

"${PLINK2}" \
  --bfile "${GENO}" \
  --keep "${KEEP_ROOT}/simulation_gwas.keep" \
  --extract "${ROOT}/qc/king_prune.prune.in" \
  --pca approx 10 \
  --threads 8 \
  --out "${SIM_ROOT}/heldout_gwas_pca"

panels=(1kgeas chinese_full_overlap chinese_disjoint_n10000)
printf 'replicate\tpanel\tstatus\tlog\n' > "${ROOT}/qc/simulation_ldsc_status.tsv"
for replicate in $(seq 1 "${REPLICATES}"); do
  rep="${SIM_ROOT}/replicate_${replicate}"
  mkdir -p "${rep}"
  python3 "${ROOT}/scripts/generate_simulation_effects.py" \
    --bim "${GENO}.bim" \
    --extract "${common_snps}" \
    --causal "${CAUSAL_SNPS}" \
    --seed "$((2026092200 + replicate))" \
    --output "${rep}/effects.tsv"
  "${PLINK}" \
    --bfile "${GENO}" \
    --keep "${KEEP_ROOT}/simulation_gwas.keep" \
    --extract "${common_snps}" \
    --score "${rep}/effects.tsv" 1 2 3 header sum \
    --out "${rep}/polygenic_score" \
    > "${ROOT}/logs/simulation/replicate_${replicate}.score.stdout.log" 2>&1
  python3 "${ROOT}/scripts/prepare_simulated_phenotype.py" \
    --score "${rep}/polygenic_score.profile" \
    --eigenvec "${SIM_ROOT}/heldout_gwas_pca.eigenvec" \
    --h2 "${TARGET_H2}" \
    --seed "$((2026092300 + replicate))" \
    --output "${rep}/phenotype.tsv" \
    --report "${rep}/phenotype_qc.json"
  "${PLINK2}" \
    --bfile "${GENO}" \
    --keep "${KEEP_ROOT}/simulation_gwas.keep" \
    --extract "${common_snps}" \
    --pheno "${rep}/phenotype.tsv" \
    --pheno-name PHENO \
    --covar "${SIM_ROOT}/heldout_gwas_pca.eigenvec" \
    --covar-name PC1-PC10 \
    --glm hide-covar omit-ref \
    --threads 8 \
    --out "${rep}/gwas" \
    > "${ROOT}/logs/simulation/replicate_${replicate}.gwas.stdout.log" 2>&1
  python3 "${ROOT}/scripts/format_plink_glm_for_ldsc.py" \
    --glm "${rep}/gwas.PHENO.glm.linear" \
    --output "${rep}/simulation.sumstats.gz"

  for panel in "${panels[@]}"; do
    out="${rep}/h2_${panel}"
    if "${LDSC_PYTHON}" "${LDSC}" \
        --h2 "${rep}/simulation.sumstats.gz" \
        --ref-ld-chr "${PANEL_ROOT}/${panel}/" \
        --w-ld-chr "${PANEL_ROOT}/1kgeas/" \
        --not-M-5-50 \
        --n-blocks 200 \
        --print-delete-vals \
        --out "${out}" \
        > "${ROOT}/logs/simulation/replicate_${replicate}.${panel}.ldsc.stdout.log" 2>&1; then
      status=complete
    else
      status=failed
    fi
    printf '%s\t%s\t%s\t%s\n' "${replicate}" "${panel}" "${status}" "${out}.log" \
      >> "${ROOT}/qc/simulation_ldsc_status.tsv"
  done
done

if grep -q $'\tfailed\t' "${ROOT}/qc/simulation_ldsc_status.tsv"; then
  echo "One or more simulated-trait LDSC analyses failed; see simulation_ldsc_status.tsv" >&2
  exit 1
fi

python3 "${ROOT}/scripts/summarise_simulation.py" \
  --root "${ROOT}" \
  --replicates "${REPLICATES}" \
  --target-h2 "${TARGET_H2}"
date --iso-8601=seconds > "${ROOT}/qc/simulation_complete.timestamp"
