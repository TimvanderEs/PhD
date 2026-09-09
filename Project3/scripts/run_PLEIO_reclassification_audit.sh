#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
  cat >&2 <<'USAGE'
Usage:
  bash run_PLEIO_reclassification_audit.sh \
    <ANCESTRY> <MODEL> <PLEIO_FUMA_DIR> <SINGLE_TRAIT_FUMA_ROOT> \
    <RAW_MERGED_PLEIO_FILE> <OUTPUT_DIR> [CF_FOLDER]

Example:
  bash run_PLEIO_reclassification_audit.sh \
    EUR FOURTRAIT \
    /home/ec2-user/project2/RESULTS/FUMA/FUMA/EUR/PLEIO \
    /home/ec2-user/project2/RESULTS/FUMA/FUMA/EUR \
    /home/ec2-user/project2/RESULTS/PLEIO/EU/HELIOS_PLEIO_EUR_summary_all_merged.tsv \
    /home/ec2-user/project2/RESULTS/PLEIO/reclassification_audit/EUR_FOURTRAIT \
    G

MODEL must be FOURTRAIT, MDD_3TRAIT, or SCZ_3TRAIT.
USAGE
  exit 2
fi

ANCESTRY=$1
MODEL=$2
PLEIO_DIR=$3
SINGLE_ROOT=$4
RAW_FILE=$5
OUT_DIR=$6
CF_FOLDER=${7:-G}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
R_SCRIPT="${SCRIPT_DIR}/classify_PLEIO_FUMA_loci_v2.R"
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
  --cf-folder "$CF_FOLDER" \
  --zthr 1.96 \
  --prefix "$PREFIX" \
  2>&1 | tee "${OUT_DIR}/${PREFIX}_run.log"

echo
echo "===== Corrected direction summary ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_PLEIO_locus_direction_summary.tsv" || \
  cat "${OUT_DIR}/${PREFIX}_PLEIO_locus_direction_summary.tsv"

echo
echo "===== Legacy -> corrected locus transitions ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_locus_direction_transition.tsv" || \
  cat "${OUT_DIR}/${PREFIX}_locus_direction_transition.tsv"

echo
echo "===== Prioritisation transitions ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_prioritisation_transition.tsv" || \
  cat "${OUT_DIR}/${PREFIX}_prioritisation_transition.tsv"

echo
echo "===== QC ====="
column -t -s $'\t' "${OUT_DIR}/${PREFIX}_classification_QC.tsv" || \
  cat "${OUT_DIR}/${PREFIX}_classification_QC.tsv"
