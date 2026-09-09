#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash prepare_four_trait_inputs.sh \
    --pleio-dir DIR --trait-config FILE --sumstats-dir DIR \
    --ref-ld-prefix PREFIX --out-dir DIR [--python CMD] [--validate-only]

The PLEIO directory must contain ldsc_preprocess.py. The trait config must have
these tab-separated columns: file, type, sample_prev, population_prev, name.
EOF
}

PLEIO_DIR=""
TRAIT_CONFIG=""
SUMSTATS_DIR=""
REF_LD_PREFIX=""
OUT_DIR=""
PYTHON_BIN="python"
VALIDATE_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pleio-dir) PLEIO_DIR="$2"; shift 2 ;;
    --trait-config) TRAIT_CONFIG="$2"; shift 2 ;;
    --sumstats-dir) SUMSTATS_DIR="$2"; shift 2 ;;
    --ref-ld-prefix) REF_LD_PREFIX="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --validate-only) VALIDATE_ONLY=1; shift ;;
    --help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value_name in PLEIO_DIR TRAIT_CONFIG SUMSTATS_DIR REF_LD_PREFIX OUT_DIR; do
  if [[ -z "${!value_name}" ]]; then
    echo "ERROR: missing required option for $value_name" >&2
    usage >&2
    exit 2
  fi
done

for path in "$PLEIO_DIR/ldsc_preprocess.py" "$TRAIT_CONFIG" "$SUMSTATS_DIR"; do
  if [[ ! -e "$path" ]]; then
    echo "ERROR: missing required path: $path" >&2
    exit 1
  fi
done
if [[ ! -d "$(dirname "$REF_LD_PREFIX")" ]]; then
  echo "ERROR: reference LD-score directory does not exist: $(dirname "$REF_LD_PREFIX")" >&2
  exit 1
fi
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: Python command not found: $PYTHON_BIN" >&2
  exit 1
fi

header="$(head -n 1 "$TRAIT_CONFIG")"
expected=$'file\ttype\tsample_prev\tpopulation_prev\tname'
if [[ "$header" != "$expected" ]]; then
  echo "ERROR: unexpected trait-config header: $header" >&2
  echo "Expected: $expected" >&2
  exit 1
fi

SUMSTATS_DIR="$(cd "$SUMSTATS_DIR" && pwd)"
PLEIO_DIR="$(cd "$PLEIO_DIR" && pwd)"
TRAIT_CONFIG="$(cd "$(dirname "$TRAIT_CONFIG")" && pwd)/$(basename "$TRAIT_CONFIG")"
mkdir -p "$(dirname "$OUT_DIR")"
OUT_DIR="$(cd "$(dirname "$OUT_DIR")" && pwd)/$(basename "$OUT_DIR")"
RESOLVED_MANIFEST="${OUT_DIR}.input_manifest.tsv"

printf 'FILE\tTYPE\tSPREV\tPPREV\tNAME\n' > "$RESOLVED_MANIFEST"
while IFS=$'\t' read -r file trait_type sample_prev population_prev trait_name; do
  [[ -n "$file" ]] || continue
  if [[ "$file" = /* ]]; then
    resolved="$file"
  else
    resolved="$SUMSTATS_DIR/$file"
  fi
  if [[ ! -f "$resolved" ]]; then
    echo "ERROR: missing summary-statistics file: $resolved" >&2
    exit 1
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$resolved" "$trait_type" "$sample_prev" "$population_prev" "$trait_name" \
    >> "$RESOLVED_MANIFEST"
done < <(tail -n +2 "$TRAIT_CONFIG")

trait_count=$(( $(wc -l < "$RESOLVED_MANIFEST") - 1 ))
if [[ "$trait_count" -lt 2 ]]; then
  echo "ERROR: at least two traits are required." >&2
  exit 1
fi

echo "Validated $trait_count traits"
echo "Resolved manifest: $RESOLVED_MANIFEST"
echo "Reference LD prefix: $REF_LD_PREFIX"
echo "Output directory: $OUT_DIR"
if [[ "$VALIDATE_ONLY" -eq 1 ]]; then
  exit 0
fi

(
  cd "$PLEIO_DIR"
  "$PYTHON_BIN" ./ldsc_preprocess.py \
    --ref-ld-chr "$REF_LD_PREFIX" \
    --out "$OUT_DIR" \
    --input "$RESOLVED_MANIFEST" \
    --w-ld-chr "$REF_LD_PREFIX"
)

for output in metain.txt.gz sg.txt.gz ce.txt.gz; do
  if [[ ! -s "$OUT_DIR/$output" ]]; then
    echo "ERROR: expected PLEIO preparation output is missing: $OUT_DIR/$output" >&2
    exit 1
  fi
done
echo "Completed four-trait PLEIO preparation: $OUT_DIR"
