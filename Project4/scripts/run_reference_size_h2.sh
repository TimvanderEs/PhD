#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
SOURCE_ROOT="${SOURCE_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
MUNGED_DIR="${MUNGED_DIR:-${SOURCE_ROOT}/munged/chinese}"
PANEL_ROOT="${ROOT}/panels/reference_size_matched"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
LDSC="${LDSC:-/home/ec2-user/software/cov-ldsc/ldsc.py}"

traits=(c1_sbp_avg c11_whr dbi13_height c8_bmi c37_dlab_ldl_mmol_l_log10 c20_t2d)
panels=(1kgeas chinese_n500 chinese_n2500 chinese_n10000 chinese_full)
mkdir -p "${ROOT}/results/reference_size_h2" "${ROOT}/logs/reference_size_h2" "${ROOT}/qc"
printf 'panel\ttrait\tstatus\tlog\n' > "${ROOT}/qc/reference_size_h2_status.tsv"

for panel in "${panels[@]}"; do
  mkdir -p "${ROOT}/results/reference_size_h2/${panel}"
  for trait in "${traits[@]}"; do
    sumstats="${MUNGED_DIR}/${trait}.chinese.sumstats.gz"
    output="${ROOT}/results/reference_size_h2/${panel}/${trait}"
    if "${PYTHON_BIN}" "${LDSC}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${PANEL_ROOT}/${panel}/" \
        --w-ld-chr "${PANEL_ROOT}/1kgeas/" \
        --not-M-5-50 \
        --n-blocks 200 \
        --print-delete-vals \
        --out "${output}" \
        > "${ROOT}/logs/reference_size_h2/${panel}.${trait}.stdout.log" 2>&1; then
      status=complete
    else
      status=failed
    fi
    printf '%s\t%s\t%s\t%s\n' "${panel}" "${trait}" "${status}" "${output}.log" \
      >> "${ROOT}/qc/reference_size_h2_status.tsv"
  done
done
if grep -q $'\tfailed\t' "${ROOT}/qc/reference_size_h2_status.tsv"; then
  echo "One or more reference-size LDSC analyses failed; see reference_size_h2_status.tsv" >&2
  exit 1
fi
python3 "${ROOT}/scripts/summarise_reference_size.py" "${ROOT}"
date --iso-8601=seconds > "${ROOT}/qc/reference_size_h2_complete.timestamp"
