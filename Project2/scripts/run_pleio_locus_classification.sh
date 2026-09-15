#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 8 ]]; then
  cat >&2 <<'USAGE'
Usage:
  bash run_pleio_locus_classification.sh \
    <ANCESTRY> <MODEL> <PLEIO_FUMA_DIR> <SINGLE_TRAIT_FUMA_ROOT> \
    <RAW_MERGED_PLEIO_FILE> <OUTPUT_DIR> [EA_FOLDER] [CF_FOLDER]

Example with canonical source-directory names:
  bash run_pleio_locus_classification.sh \
    EUR FOURTRAIT \
    /path/to/FUMA/EUR/PLEIO \
    /path/to/FUMA/EUR \
    /path/to/raw_merged_PLEIO.tsv.gz \
    /path/to/output/EUR_FOURTRAIT

MODEL must be FOURTRAIT, MDD_3TRAIT, or SCZ_3TRAIT.
EA_FOLDER and CF_FOLDER default to EA and CF. Supply different values only when
the source archive uses historical directory names; output fields remain EA/CF.
USAGE
  exit 2
fi

ANCESTRY=$1
MODEL=$2
PLEIO_DIR=$3
SINGLE_ROOT=$4
RAW_FILE=$5
OUT_DIR=$6
EA_FOLDER=${7:-EA}
CF_FOLDER=${8:-CF}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
R_SCRIPT="${SCRIPT_DIR}/classify_pleio_locus_direction.R"
PREFIX="${ANCESTRY}_${MODEL}"

if [[ ! -f "$R_SCRIPT" ]]; then
  echo "ERROR: Cannot find $R_SCRIPT" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

Rscript "$R_SCRIPT" \
  --ancestry "$ANCESTRY" \
  --model "$MODEL" \
  --pleio-dir "$PLEIO_DIR" \
  --single-trait-root "$SINGLE_ROOT" \
  --raw-file "$RAW_FILE" \
  --out-dir "$OUT_DIR" \
  --ea-folder "$EA_FOLDER" \
  --cf-folder "$CF_FOLDER" \
  --zthr 1.96 \
  --prefix "$PREFIX" \
  2>&1 | tee "${OUT_DIR}/${PREFIX}_run.log"

echo
echo "===== Locus direction summary ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_PLEIO_locus_direction_summary.tsv" || \
  sed -n '1,40p' "${OUT_DIR}/${PREFIX}_PLEIO_locus_direction_summary.tsv"

echo
echo "===== Locus prioritisation summary ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_PLEIO_locus_prioritisation_summary.tsv" || \
  sed -n '1,20p' "${OUT_DIR}/${PREFIX}_PLEIO_locus_prioritisation_summary.tsv"

echo
echo "===== QC ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_classification_QC.tsv" || \
  sed -n '1,20p' "${OUT_DIR}/${PREFIX}_classification_QC.tsv"
