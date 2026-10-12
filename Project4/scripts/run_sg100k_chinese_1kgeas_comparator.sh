#!/usr/bin/env bash
set -euo pipefail

# Public-panel comparator for the Chinese-specific SG100K GWAS.

OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
EAS_LD_ROOT="${EAS_LD_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched/1kgeas_singapore_ids}"
LDSC_SCRIPT="${LDSC_SCRIPT:-/home/ec2-user/software/cov-ldsc/ldsc.py}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"

MUNGED_DIR="${OUTPUT_ROOT}/munged/chinese"
H2_DIR="${OUTPUT_ROOT}/h2/Chinese_1KGEAS"
RG_DIR="${OUTPUT_ROOT}/rg/Chinese_1KGEAS"
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

mkdir -p "${H2_DIR}" "${RG_DIR}" "${LOG_DIR}" "${QC_DIR}"
status_file="${QC_DIR}/chinese_1kgeas_run_status.tsv"
printf 'analysis\ttrait_1\ttrait_2\tstatus\tlog\n' > "${status_file}"

if [[ $(find "${EAS_LD_ROOT}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l) -ne 22 ]]; then
  echo "Incomplete 1000 Genomes EAS LD panel: ${EAS_LD_ROOT}" >&2
  exit 1
fi

first_reference_snp=$(python3 -c 'import gzip, sys; handle = gzip.open(sys.argv[1], "rt"); next(handle); print(next(handle).split()[1])' "${EAS_LD_ROOT}/1.l2.ldscore.gz")
if [[ "${first_reference_snp}" != chr*:*:*:* ]]; then
  echo "1KG EAS panel has not been translated to Singapore variant identifiers: ${first_reference_snp}" >&2
  exit 1
fi

files=()
for trait in "${traits[@]}"; do
  sumstats="${MUNGED_DIR}/${trait}.chinese.sumstats.gz"
  [[ -s "${sumstats}" ]] || { echo "Missing input: ${sumstats}" >&2; exit 1; }
  files+=("${sumstats}")
  out="${H2_DIR}/${trait}.chinese_1kgeas"
  echo "[$(date --iso-8601=seconds)] h2 Chinese GWAS with 1KG EAS: ${trait}"
  if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
      --h2 "${sumstats}" \
      --ref-ld-chr "${EAS_LD_ROOT}/" \
      --w-ld-chr "${EAS_LD_ROOT}/" \
      --out "${out}" \
      > "${LOG_DIR}/${trait}.chinese_1kgeas.h2.stdout.log" 2>&1; then
    status=complete
  else
    status=failed
  fi
  printf 'h2\t%s\tNA\t%s\t%s\n' "${trait}" "${status}" "${out}.log" >> "${status_file}"
done

for ((i = 0; i < ${#files[@]}; i++)); do
  for ((j = i + 1; j < ${#files[@]}; j++)); do
    first="${files[$i]}"
    second="${files[$j]}"
    trait_1="${traits[$i]}"
    trait_2="${traits[$j]}"
    out="${RG_DIR}/${trait_1}__${trait_2}.chinese_1kgeas"
    echo "[$(date --iso-8601=seconds)] rg Chinese GWAS with 1KG EAS: ${trait_1} vs ${trait_2}"
    if "${PYTHON_BIN}" "${LDSC_SCRIPT}" \
        --rg "${first},${second}" \
        --ref-ld-chr "${EAS_LD_ROOT}/" \
        --w-ld-chr "${EAS_LD_ROOT}/" \
        --out "${out}" \
        > "${LOG_DIR}/${trait_1}__${trait_2}.chinese_1kgeas.rg.stdout.log" 2>&1; then
      status=complete
    else
      status=failed
    fi
    printf 'rg\t%s\t%s\t%s\t%s\n' "${trait_1}" "${trait_2}" "${status}" "${out}.log" >> "${status_file}"
  done
done

echo "[$(date --iso-8601=seconds)] Chinese 1KG EAS comparison complete"
