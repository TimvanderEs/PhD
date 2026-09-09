#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

bash "${PROJECT_DIR}/examples/run_synthetic_example.sh"
Rscript "${PROJECT_DIR}/tests/validate_synthetic_output.R" \
  "${PROJECT_DIR}/examples/synthetic/output"
