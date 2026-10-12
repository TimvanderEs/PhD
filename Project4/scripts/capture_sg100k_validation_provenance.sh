#!/usr/bin/env bash
set -euo pipefail

OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
LD_ROOT="${LD_ROOT:-/home/ec2-user/project3/sg100k_ancestry_ldscores_1cm}"
EAS_LD_ROOT="${EAS_LD_ROOT:-/home/ec2-user/project3/1kgeas_consensus_867k}"
EAS_TRANSLATED_ROOT="${EAS_TRANSLATED_ROOT:-${OUTPUT_ROOT}/1kgeas_singapore_ids}"
LDSC_ROOT="${LDSC_ROOT:-/home/ec2-user/software/cov-ldsc}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"

QC_DIR="${OUTPUT_ROOT}/qc"
mkdir -p "${QC_DIR}"

{
  echo "captured_utc=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "kernel=$(uname -srmo)"
  echo "python=$(${PYTHON_BIN} --version 2>&1)"
  echo "ldsc_commit=$(git -C "${LDSC_ROOT}" rev-parse HEAD 2>/dev/null || echo unavailable)"
  echo "ldsc_status=$(git -C "${LDSC_ROOT}" status --porcelain 2>/dev/null | wc -l | tr -d ' ') modified_or_untracked_files"
  for ancestry in Chinese Indian Malay; do
    count=$(find "${LD_ROOT}/${ancestry}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l | tr -d ' ')
    echo "Singapore_${ancestry}_ldscore_chromosomes=${count}"
  done
  eas_count=$(find "${EAS_LD_ROOT}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l | tr -d ' ')
  echo "1KG_EAS_ldscore_chromosomes=${eas_count}"
  eas_translated_count=$(find "${EAS_TRANSLATED_ROOT}" -maxdepth 1 -name '*.l2.ldscore.gz' | wc -l | tr -d ' ')
  echo "1KG_EAS_translated_ldscore_chromosomes=${eas_translated_count}"
} > "${QC_DIR}/session_info.txt"

sha256sum \
  "${LDSC_ROOT}/ldsc.py" \
  "${LDSC_ROOT}/munge_sumstats.py" \
  > "${QC_DIR}/software_and_script_checksums.sha256"

find "${OUTPUT_ROOT}/scripts" -maxdepth 1 -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  >> "${QC_DIR}/software_and_script_checksums.sha256"

find "${OUTPUT_ROOT}/h2" "${OUTPUT_ROOT}/rg" -type f -name '*.log' -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "${QC_DIR}/result_log_checksums.sha256"
