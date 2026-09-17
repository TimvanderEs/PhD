#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <plink1.9> <input-bfile-prefix> <ancestry-label> <high-LD-range-file> <output-directory>" >&2
  exit 2
fi

plink=$1
input_prefix=$2
ancestry=$3
high_ld_ranges=$4
output_directory=$5

mkdir -p "${output_directory}"

qc1_prefix="${output_directory}/HELIOS_${ancestry}_QC1"
prune_prefix="${output_directory}/Indep_SNP_${ancestry}"
report_prefix="${output_directory}/HELIOS_${ancestry}_reports"
heterozygosity_exclusions="${output_directory}/${ancestry}_heterozygosity_exclusions.txt"
relatedness_exclusions="${output_directory}/${ancestry}_relatedness_exclusions.txt"
sample_exclusions="${output_directory}/${ancestry}_sample_exclusions.txt"
final_prefix="${output_directory}/HELIOS_${ancestry}_QCrem_miss02_maf005_hwe1e06"

# Variant and sample missingness, allele frequency and HWE filters recovered
# from the historical QC scripts.
"${plink}" \
  --bfile "${input_prefix}" \
  --geno 0.02 \
  --mind 0.02 \
  --maf 0.005 \
  --hwe 1e-6 \
  --make-bed \
  --out "${qc1_prefix}"

# Common, LD-pruned variants outside the supplied high-LD regions were used
# for heterozygosity, relatedness and ancestry-specific PCA.
"${plink}" \
  --bfile "${qc1_prefix}" \
  --maf 0.05 \
  --exclude range "${high_ld_ranges}" \
  --indep-pairwise 200 100 0.1 \
  --out "${prune_prefix}"

"${plink}" \
  --bfile "${qc1_prefix}" \
  --extract "${prune_prefix}.prune.in" \
  --het \
  --genome \
  --out "${report_prefix}"

awk 'NR > 1 && ($6 < -0.2 || $6 > 0.2) {print $1, $2}' \
  "${report_prefix}.het" > "${heterozygosity_exclusions}"

awk 'NR > 1 && $10 > 0.75 {print $1, $2; print $3, $4}' \
  "${report_prefix}.genome" > "${relatedness_exclusions}"

sort -u "${heterozygosity_exclusions}" "${relatedness_exclusions}" \
  > "${sample_exclusions}"

"${plink}" \
  --bfile "${qc1_prefix}" \
  --remove "${sample_exclusions}" \
  --geno 0.02 \
  --maf 0.005 \
  --hwe 1e-6 \
  --make-bed \
  --out "${final_prefix}"

"${plink}" \
  --bfile "${final_prefix}" \
  --extract "${prune_prefix}.prune.in" \
  --pca 20 header \
  --out "${final_prefix}"

echo "QC complete: ${final_prefix}"
