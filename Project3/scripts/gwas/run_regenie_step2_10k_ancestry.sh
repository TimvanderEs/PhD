#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
  echo "Usage: $0 <regenie> <bed-prefix> <phenotype-file> <covariate-file> <prediction-list> <output-prefix> [threads]" >&2
  exit 2
fi

regenie=$1
bed_prefix=$2
phenotype_file=$3
covariate_file=$4
prediction_list=$5
output_prefix=$6
threads=${7:-8}

"${regenie}" \
  --step 2 \
  --bed "${bed_prefix}" \
  --phenoFile "${phenotype_file}" \
  --phenoCol g \
  --covarFile "${covariate_file}" \
  --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
  --bsize 200 \
  --minMAC 0.5 \
  --threads "${threads}" \
  --pred "${prediction_list}" \
  --out "${output_prefix}"
