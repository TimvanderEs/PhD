#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON:-python3}"

exec "${python_bin}" "${script_dir}/scripts/run_22k_validation.py" "$@"
