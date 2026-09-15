#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/synthetic" && pwd)
PROJECT_DIR=$(cd -- "${EXAMPLE_DIR}/../.." && pwd)

bash "${PROJECT_DIR}/scripts/run_pleio_locus_classification.sh" \
  SYNTHETIC \
  FOURTRAIT \
  "${EXAMPLE_DIR}/fuma/PLEIO" \
  "${EXAMPLE_DIR}/fuma" \
  "${EXAMPLE_DIR}/raw_merged_PLEIO.tsv" \
  "${EXAMPLE_DIR}/output"
