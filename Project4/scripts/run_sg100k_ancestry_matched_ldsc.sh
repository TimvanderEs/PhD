#!/usr/bin/env bash
set -euo pipefail

OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
LD_ROOT="${LD_ROOT:-/home/ec2-user/project3/sg100k_ancestry_ldscores_1cm}"
LDSC_SCRIPT="${LDSC_SCRIPT:-/home/ec2-user/software/cov-ldsc/ldsc.py}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"

MUNGED_DIR="${OUTPUT_ROOT}/munged"
H2_DIR="${OUTPUT_ROOT}/h2"
RG_DIR="${OUTPUT_ROOT}/rg"
LOG_DIR="${OUTPUT_ROOT}/logs"
QC_DIR="${OUTPUT_ROOT}/qc"

traits=(
  c1_sbp_avg
  c11_whr
  dbi13_height
  c8_bmi
  c37_dlab_ldl_mmol_l_log10
  c20_t2d
)
ancestries=(Chinese Indian Malay)

mkdir -p "${H2_DIR}" "${RG_DIR}" "${LOG_DIR}" "${QC_DIR}"
printf 'analysis\tancestry\ttrait_1\ttrait_2\tstatus\tlog\n' > "${QC_DIR}/run_status.tsv"

for ancestry in "${ancestries[@]}"; do
  ancestry_lower=$(printf '%s' "${ancestry}" | tr '[:upper:]' '[:lower:]')
  reference="${LD_ROOT}/${ancestry}/"
  if [[ $(find "${reference}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l) -ne 22 ]]; then
    echo "Incomplete LD panel: ${reference}" >&2
    exit 1
  fi

  mkdir -p "${H2_DIR}/${ancestry}" "${RG_DIR}/${ancestry}"
  files=()
  for trait in "${traits[@]}"; do
    sumstats="${MUNGED_DIR}/${ancestry_lower}/${trait}.${ancestry_lower}.sumstats.gz"
    if [[ ! -s "${sumstats}" ]]; then
      echo "Missing munged input: ${sumstats}" >&2
      exit 1
    fi
    files+=("${sumstats}")

    out="${H2_DIR}/${ancestry}/${trait}.${ancestry_lower}_matched"
    stdout_log="${LOG_DIR}/${trait}.${ancestry_lower}.h2.stdout.log"
    echo "[$(date --iso-8601=seconds)] h2 ${ancestry}: ${trait}"
    if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --h2 "${sumstats}" \
        --ref-ld-chr "${reference}" \
        --w-ld-chr "${reference}" \
        --out "${out}" \
        > "${stdout_log}" 2>&1; then
      status=complete
    else
      status=failed
    fi
    printf 'h2\t%s\t%s\tNA\t%s\t%s\n' "${ancestry}" "${trait}" "${status}" "${out}.log" >> "${QC_DIR}/run_status.tsv"
  done

  for ((i = 0; i < ${#files[@]}; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do
      first="${files[$i]}"
      second="${files[$j]}"
      trait_1="${traits[$i]}"
      trait_2="${traits[$j]}"
      out="${RG_DIR}/${ancestry}/${trait_1}__${trait_2}.${ancestry_lower}_matched"
      stdout_log="${LOG_DIR}/${trait_1}__${trait_2}.${ancestry_lower}.rg.stdout.log"
      echo "[$(date --iso-8601=seconds)] rg ${ancestry}: ${trait_1} vs ${trait_2}"
      if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
          --rg "${first},${second}" \
          --ref-ld-chr "${reference}" \
          --w-ld-chr "${reference}" \
          --out "${out}" \
          > "${stdout_log}" 2>&1; then
        status=complete
      else
        status=failed
      fi
      printf 'rg\t%s\t%s\t%s\t%s\t%s\n' "${ancestry}" "${trait_1}" "${trait_2}" "${status}" "${out}.log" >> "${QC_DIR}/run_status.tsv"
    done
  done
done

echo "[$(date --iso-8601=seconds)] Matched LDSC analyses complete"
