#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --ldsc /path/to/ldsc.py --reference-prefix /path/to/prefix@ --prepared-dir DIR --output-dir DIR" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ldsc) LDSC=$2; shift 2 ;;
    --reference-prefix) REFERENCE=$2; shift 2 ;;
    --prepared-dir) PREPARED=$2; shift 2 ;;
    --output-dir) OUTPUT=$2; shift 2 ;;
    *) usage ;;
  esac
done

[[ -n "${LDSC:-}" && -n "${REFERENCE:-}" && -n "${PREPARED:-}" && -n "${OUTPUT:-}" ]] || usage
mkdir -p "$OUTPUT"

for rule in reference_maf01 archived_qc_afdiff050 strict_qc_afdiff005; do
  "$LDSC" \
    --rg "$PREPARED/trans10k_n7789.${rule}.sumstats.gz,$PREPARED/extension12k_n13751.${rule}.sumstats.gz,$PREPARED/combined22k_n21518.${rule}.sumstats.gz" \
    --ref-ld-chr "$REFERENCE" \
    --w-ld-chr "$REFERENCE" \
    --out "$OUTPUT/$rule"
done

for analysis in trans10k_n7789 extension12k_n13751 combined22k_n21518; do
  "$LDSC" \
    --h2 "$PREPARED/${analysis}.reference_maf01.sumstats.gz" \
    --ref-ld-chr "$REFERENCE" \
    --w-ld-chr "$REFERENCE" \
    --no-intercept \
    --out "$OUTPUT/${analysis}.reference_maf01.intercept1"
done
