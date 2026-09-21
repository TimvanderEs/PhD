#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <regenie> <bed-prefix> <phenotype-file> <covariate-file> <prediction-list> <output-prefix>" >&2
  exit 2
fi

regenie=$1
bed_prefix=$2
phenotype_file=$3
covariate_file=$4
prediction_list=$5
output_prefix=$6

"${regenie}" \
  --step 2 \
  --bed "${bed_prefix}" \
  --covarFile "${covariate_file}" \
  --phenoFile "${phenotype_file}" \
  --bsize 200 \
  --minMAC 0.5 \
  --pred "${prediction_list}" \
  --out "${output_prefix}"
