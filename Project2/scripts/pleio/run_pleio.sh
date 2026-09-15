#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash run_pleio.sh --pleio-py FILE --input-dir DIR --out-prefix PATH
    [--python CMD] [--nis N] [--ncores N] [--isf FILE]

By default the script uses --create, --parallel, and 100000 importance samples.
Supplying --isf reuses an existing inverse-CDF file instead of --create.
EOF
}

PLEIO_PY=""
INPUT_DIR=""
OUT_PREFIX=""
PYTHON_BIN="python"
NIS="100000"
NCORES=""
ISF=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pleio-py) PLEIO_PY="$2"; shift 2 ;;
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --nis) NIS="$2"; shift 2 ;;
    --ncores) NCORES="$2"; shift 2 ;;
    --isf) ISF="$2"; shift 2 ;;
    --help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$PLEIO_PY" || -z "$INPUT_DIR" || -z "$OUT_PREFIX" ]]; then
  usage >&2
  exit 2
fi
for path in "$PLEIO_PY" "$INPUT_DIR/metain.txt.gz" "$INPUT_DIR/sg.txt.gz" "$INPUT_DIR/ce.txt.gz"; do
  if [[ ! -s "$path" ]]; then
    echo "ERROR: missing required file: $path" >&2
    exit 1
  fi
done
if ! [[ "$NIS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --nis must be a positive integer." >&2
  exit 2
fi
if [[ -n "$NCORES" ]] && ! [[ "$NCORES" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --ncores must be a positive integer." >&2
  exit 2
fi
if [[ -n "$ISF" && ! -s "$ISF" ]]; then
  echo "ERROR: --isf file does not exist: $ISF" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"
PLEIO_PY="$(cd "$(dirname "$PLEIO_PY")" && pwd)/$(basename "$PLEIO_PY")"
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUT_PREFIX="$(cd "$(dirname "$OUT_PREFIX")" && pwd)/$(basename "$OUT_PREFIX")"

command=(
  "$PYTHON_BIN" "$PLEIO_PY"
  --metain "$INPUT_DIR/metain.txt.gz"
  --sg "$INPUT_DIR/sg.txt.gz"
  --ce "$INPUT_DIR/ce.txt.gz"
  --nis "$NIS"
  --parallel
  --out "$OUT_PREFIX"
)
if [[ -n "$NCORES" ]]; then
  command+=(--ncores "$NCORES")
fi
if [[ -n "$ISF" ]]; then
  command+=(--isf "$ISF")
else
  command+=(--create)
fi

printf 'Running:'
printf ' %q' "${command[@]}"
printf '\n'
(
  cd "$(dirname "$PLEIO_PY")"
  "${command[@]}"
) > "${OUT_PREFIX}.stdout.log" 2>&1

if [[ ! -s "${OUT_PREFIX}.txt.gz" ]]; then
  echo "ERROR: PLEIO did not create ${OUT_PREFIX}.txt.gz" >&2
  exit 1
fi
echo "Completed PLEIO run: $OUT_PREFIX"
