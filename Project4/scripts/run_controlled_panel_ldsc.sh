#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
SOURCE_ROOT="${SOURCE_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
MUNGED_DIR="${MUNGED_DIR:-${SOURCE_ROOT}/munged/chinese}"
WEIGHT_PANEL="${WEIGHT_PANEL:-${ROOT}/panels/1kgeas_matched}"
LDSC_SCRIPT="${LDSC_SCRIPT:-/home/ec2-user/software/cov-ldsc/ldsc.py}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
N_BLOCKS="${N_BLOCKS:-200}"

traits=(
  c1_sbp_avg
  c11_whr
  dbi13_height
  c8_bmi
  c37_dlab_ldl_mmol_l_log10
  c20_t2d
)
panels=(1kgeas_matched chinese_standard chinese_cov20pc pooled_standard pooled_cov20pc)

mkdir -p "${ROOT}/results/h2" "${ROOT}/results/rg" "${ROOT}/logs/ldsc" "${ROOT}/qc"
printf 'analysis\tpanel\ttrait_1\ttrait_2\tstatus\tlog\n' > "${ROOT}/qc/controlled_ldsc_status.tsv"

for panel in "${panels[@]}"; do
  reference="${ROOT}/panels/${panel}/"
  if [[ $(find "${reference}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l) -ne 22 ]]; then
    echo "Incomplete panel: ${reference}" >&2
    exit 1
  fi
  mkdir -p "${ROOT}/results/h2/${panel}" "${ROOT}/results/rg/${panel}"
  sumstat_files=()
  for trait in "${traits[@]}"; do
    sumstats="${MUNGED_DIR}/${trait}.chinese.sumstats.gz"
    test -s "${sumstats}"
    sumstat_files+=("${sumstats}")
    output="${ROOT}/results/h2/${panel}/${trait}"
    if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${reference}" \
        --w-ld-chr "${WEIGHT_PANEL}/" \
        --not-M-5-50 \
        --n-blocks "${N_BLOCKS}" \
        --print-delete-vals \
        --out "${output}" \
        > "${ROOT}/logs/ldsc/${panel}.${trait}.h2.stdout.log" 2>&1; then
      status=complete
    else
      status=failed
    fi
    printf 'h2\t%s\t%s\tNA\t%s\t%s\n' \
      "${panel}" "${trait}" "${status}" "${output}.log" \
      >> "${ROOT}/qc/controlled_ldsc_status.tsv"
  done

  for ((i = 0; i < ${#traits[@]}; i++)); do
    for ((j = i + 1; j < ${#traits[@]}; j++)); do
      first_trait="${traits[$i]}"
      second_trait="${traits[$j]}"
      output="${ROOT}/results/rg/${panel}/${first_trait}__${second_trait}"
      if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
          --rg "${sumstat_files[$i]},${sumstat_files[$j]}" \
          --ref-ld-chr "${reference}" \
          --w-ld-chr "${WEIGHT_PANEL}/" \
          --not-M-5-50 \
          --n-blocks "${N_BLOCKS}" \
          --print-delete-vals \
          --out "${output}" \
          > "${ROOT}/logs/ldsc/${panel}.${first_trait}__${second_trait}.rg.stdout.log" 2>&1; then
        status=complete
      else
        status=failed
      fi
      printf 'rg\t%s\t%s\t%s\t%s\t%s\n' \
        "${panel}" "${first_trait}" "${second_trait}" "${status}" "${output}.log" \
        >> "${ROOT}/qc/controlled_ldsc_status.tsv"
    done
  done
done

if grep -q $'\tfailed\t' "${ROOT}/qc/controlled_ldsc_status.tsv"; then
  echo "One or more controlled LDSC analyses failed; see controlled_ldsc_status.tsv" >&2
  exit 1
fi

date --iso-8601=seconds > "${ROOT}/qc/controlled_ldsc_complete.timestamp"
