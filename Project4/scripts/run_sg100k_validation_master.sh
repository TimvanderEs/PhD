#!/usr/bin/env bash
set -euo pipefail

OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
SCRIPT_DIR="${OUTPUT_ROOT}/scripts"
LOG_DIR="${OUTPUT_ROOT}/logs"
EAS_SOURCE_ROOT="${EAS_SOURCE_ROOT:-/home/ec2-user/project3/1kgeas_consensus_867k}"
RSID_MAPPING="${RSID_MAPPING:-/home/ec2-user/COV_LDSC_paper/rsid_to_marker.map}"
EAS_TRANSLATED_ROOT="${EAS_TRANSLATED_ROOT:-${OUTPUT_ROOT}/1kgeas_singapore_ids}"

mkdir -p "${LOG_DIR}" "${OUTPUT_ROOT}/tables"

expected_inputs=18
prepared_inputs=$(find "${OUTPUT_ROOT}/munged" -type f -name '*.sumstats.gz' -size +0c | wc -l)
if [[ "${prepared_inputs}" -ne "${expected_inputs}" ]]; then
  echo "Expected ${expected_inputs} prepared inputs; found ${prepared_inputs}" >&2
  exit 1
fi

echo "[$(date --iso-8601=seconds)] Starting ancestry-matched SG100K LDSC validation"
bash "${SCRIPT_DIR}/run_sg100k_ancestry_matched_ldsc.sh"

echo "[$(date --iso-8601=seconds)] Translating 1KG EAS rsIDs to Singapore identifiers"
python3 "${SCRIPT_DIR}/translate_1kgeas_to_singapore_ids.py" \
  --source "${EAS_SOURCE_ROOT}" \
  --mapping "${RSID_MAPPING}" \
  --output "${EAS_TRANSLATED_ROOT}"

echo "[$(date --iso-8601=seconds)] Starting Chinese GWAS / 1KG EAS comparator"
bash "${SCRIPT_DIR}/run_sg100k_chinese_1kgeas_comparator.sh"

echo "[$(date --iso-8601=seconds)] Finalising validation tables and provenance"
bash "${SCRIPT_DIR}/finalise_sg100k_validation.sh"

echo "[$(date --iso-8601=seconds)] Validation workflow complete"
