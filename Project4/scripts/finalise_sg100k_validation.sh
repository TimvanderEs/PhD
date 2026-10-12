#!/usr/bin/env bash
set -euo pipefail

OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
SCRIPT_DIR="${OUTPUT_ROOT}/scripts"

python3 "${SCRIPT_DIR}/summarise_sg100k_ancestry_matched.py" "${OUTPUT_ROOT}"
python3 "${SCRIPT_DIR}/summarise_sg100k_validation_results.py" "${OUTPUT_ROOT}/tables"
bash "${SCRIPT_DIR}/capture_sg100k_validation_provenance.sh"

matched_complete=$(awk -F '\t' 'NR > 1 && $5 == "complete" {count++} END {print count + 0}' "${OUTPUT_ROOT}/qc/run_status.tsv")
comparator_complete=$(awk -F '\t' 'NR > 1 && $4 == "complete" {count++} END {print count + 0}' "${OUTPUT_ROOT}/qc/chinese_1kgeas_run_status.tsv")

if [[ "${matched_complete}" -ne 63 ]]; then
  echo "Expected 63 completed matched analyses; found ${matched_complete}" >&2
  exit 1
fi
if [[ "${comparator_complete}" -ne 21 ]]; then
  echo "Expected 21 completed comparator analyses; found ${comparator_complete}" >&2
  exit 1
fi

echo "Validated 63 matched and 21 comparator analyses."
