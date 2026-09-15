#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_smr_reanalysis.sh <smr_root> <output_dir> [supplementary_tsv_dir]

smr_root must contain the EU_webSMR and EAS_SMR_web directories documented
in README.md. supplementary_tsv_dir is optional; when supplied, it enables
the PLEIO/LAVA convergence tables.
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SMR_ROOT="$1"
OUTPUT_DIR="$2"
SUPP_TSV_DIR="${3:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if [[ ! -d "$SMR_ROOT" ]]; then
  echo "SMR input directory not found: $SMR_ROOT" >&2
  exit 1
fi

args=(
  --smr-root "$SMR_ROOT"
  --output "$OUTPUT_DIR"
)

if [[ -n "$SUPP_TSV_DIR" ]]; then
  if [[ ! -d "$SUPP_TSV_DIR" ]]; then
    echo "Supplementary TSV directory not found: $SUPP_TSV_DIR" >&2
    exit 1
  fi
  args+=(--supp-tsv "$SUPP_TSV_DIR")
fi

"$PYTHON_BIN" "$SCRIPT_DIR/reanalyse_smr_heidi.py" "${args[@]}"

echo
echo "SMR/HEIDI reanalysis complete: $OUTPUT_DIR"
echo "Primary summary: $OUTPUT_DIR/21_SMR_threshold_comparison.tsv"
