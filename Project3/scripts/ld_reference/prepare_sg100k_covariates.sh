#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <plink1.9> <plink2> <bfile-prefix> <high-LD-exclusion-file> <output-directory>" >&2
  exit 2
fi

plink1=$1
plink2=$2
bfile_prefix=$3
high_ld_exclusion=$4
output_directory=$5

mkdir -p "${output_directory}"

pruned_prefix="${output_directory}/independent_snps"
unrelated_prefix="${output_directory}/sg100k_unrelated"
pca_prefix="${output_directory}/sg100k_unrelated_pca"

"${plink1}" \
  --bfile "${bfile_prefix}" \
  --maf 0.01 \
  --exclude "${high_ld_exclusion}" \
  --indep-pairphase 500 50 0.2 \
  --out "${pruned_prefix}"

"${plink2}" \
  --bfile "${bfile_prefix}" \
  --extract "${pruned_prefix}.prune.in" \
  --king-cutoff 0.0884 \
  --out "${unrelated_prefix}"

"${plink2}" \
  --bfile "${bfile_prefix}" \
  --keep "${unrelated_prefix}.king.cutoff.in.id" \
  --extract "${pruned_prefix}.prune.in" \
  --freq counts \
  --threads 35 \
  --pca 50 \
  --out "${pca_prefix}"
