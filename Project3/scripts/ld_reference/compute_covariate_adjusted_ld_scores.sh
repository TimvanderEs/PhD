#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <python> <cov-ldsc.py> <chromosome-bfile-template> <covariate-file> <output-directory>" >&2
  echo "Template example: /data/sg100k_chr{CHR}_matched" >&2
  exit 2
fi

python_executable=$1
cov_ldsc_script=$2
bfile_template=$3
covariate_file=$4
output_directory=$5

mkdir -p "${output_directory}"

batch_size=6
running=0

for chromosome in $(seq 1 22); do
  bfile_prefix=${bfile_template//\{CHR\}/${chromosome}}
  output_prefix="${output_directory}/sg100k_chr${chromosome}_20cm_covldsc"

  "${python_executable}" "${cov_ldsc_script}" \
    --bfile "${bfile_prefix}" \
    --l2 \
    --ld-wind-cm 20 \
    --cov "${covariate_file}" \
    --out "${output_prefix}" &

  running=$((running + 1))
  if (( running % batch_size == 0 )); then
    wait
  fi
done

wait
