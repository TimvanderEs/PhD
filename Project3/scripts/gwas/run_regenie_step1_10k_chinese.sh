#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <regenie> <bed-prefix> <pruned-snps> <phenotype-file> <covariate-file> <output-prefix>" >&2
  exit 2
fi

regenie=$1
bed_prefix=$2
pruned_snps=$3
phenotype_file=$4
covariate_file=$5
output_prefix=$6

"${regenie}" \
  --step 1 \
  --extract "${pruned_snps}" \
  --bed "${bed_prefix}" \
  --phenoFile "${phenotype_file}" \
  --phenoCol g_raw \
  --covarFile "${covariate_file}" \
  --covarColList Age,Sex,Age2,AgeSex,Age2Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
  --catCovarList Sex \
  --bsize 1000 \
  --out "${output_prefix}"
