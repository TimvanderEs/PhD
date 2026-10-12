#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
GENO="${GENO:-${ROOT}/input/sg100k_consensus}"
PLINK="${PLINK:-/home/ec2-user/software/plink/plink}"
PLINK2="${PLINK2:-/home/ec2-user/software/plink/plink2}"
THREADS="${THREADS:-8}"
PCS="${PCS:-20}"

mkdir -p "${ROOT}/qc" "${ROOT}/logs" "${ROOT}/work"
for extension in bed bim fam; do
  test -s "${GENO}.${extension}"
done
for ancestry in chinese indian malay; do
  test -s "${ROOT}/input/${ancestry}.keep"
done

sort -u \
  "${ROOT}/input/chinese.keep" \
  "${ROOT}/input/indian.keep" \
  "${ROOT}/input/malay.keep" \
  > "${ROOT}/input/pooled.keep"

"${PLINK}" \
  --bfile "${GENO}" \
  --keep "${ROOT}/input/pooled.keep" \
  --maf 0.01 \
  --geno 0.02 \
  --indep-pairwise 500 50 0.2 \
  --threads "${THREADS}" \
  --out "${ROOT}/qc/king_prune"

"${PLINK2}" \
  --bfile "${GENO}" \
  --extract "${ROOT}/qc/king_prune.prune.in" \
  --king-cutoff 0.0884 \
  --threads "${THREADS}" \
  --out "${ROOT}/qc/king"

KING_KEEP="${ROOT}/qc/king.king.cutoff.in.id"
test -s "${KING_KEEP}"
for ancestry in chinese indian malay; do
  awk 'NR==FNR { keep[$1 SUBSEP $2]=1; next } (($1 SUBSEP $2) in keep)' \
    "${KING_KEEP}" "${ROOT}/input/${ancestry}.keep" \
    > "${ROOT}/input/${ancestry}_unrelated.keep"
done
sort -u \
  "${ROOT}/input/chinese_unrelated.keep" \
  "${ROOT}/input/indian_unrelated.keep" \
  "${ROOT}/input/malay_unrelated.keep" \
  > "${ROOT}/input/pooled_unrelated.keep"

printf 'sample_set\tcount\n' > "${ROOT}/qc/reference_sample_counts.tsv"
for sample_set in chinese indian malay pooled chinese_unrelated indian_unrelated malay_unrelated pooled_unrelated; do
  printf '%s\t%s\n' "${sample_set}" "$(wc -l < "${ROOT}/input/${sample_set}.keep")" \
    >> "${ROOT}/qc/reference_sample_counts.tsv"
done

for sample_set in pooled_unrelated chinese_unrelated; do
  "${PLINK2}" \
    --bfile "${GENO}" \
    --keep "${ROOT}/input/${sample_set}.keep" \
    --extract "${ROOT}/qc/king_prune.prune.in" \
    --pca approx "${PCS}" \
    --threads "${THREADS}" \
    --out "${ROOT}/qc/pca_${sample_set}"
  python3 "${ROOT}/scripts/prepare_covldsc_eigenvectors.py" \
    --eigenvec "${ROOT}/qc/pca_${sample_set}.eigenvec" \
    --components "${PCS}" \
    --output "${ROOT}/input/pca_${sample_set}_${PCS}pc.eigenvec" \
    --report "${ROOT}/qc/pca_${sample_set}_${PCS}pc.json"
done

sha256sum \
  "${GENO}.bim" \
  "${GENO}.fam" \
  "${ROOT}/qc/king_prune.prune.in" \
  "${KING_KEEP}" \
  "${ROOT}/input/pooled_unrelated.keep" \
  "${ROOT}/input/chinese_unrelated.keep" \
  "${ROOT}/input/pca_pooled_unrelated_${PCS}pc.eigenvec" \
  "${ROOT}/input/pca_chinese_unrelated_${PCS}pc.eigenvec" \
  > "${ROOT}/qc/input_checksums.sha256"

date --iso-8601=seconds > "${ROOT}/qc/inputs_complete.timestamp"
