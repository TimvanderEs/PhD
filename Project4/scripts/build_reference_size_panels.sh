#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
GENO="${GENO:-${ROOT}/input/sg100k_consensus}"
PLINK="${PLINK:-/home/ec2-user/software/plink/plink}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
LDSC="${LDSC:-/home/ec2-user/software/cov-ldsc/ldsc.py}"
MAP_ROOT="${MAP_ROOT:-${ROOT}/maps/GRCh38/no_chr_in_chrom_field}"

KEEP_ROOT="${ROOT}/input/downsamples"
RAW_ROOT="${ROOT}/panels/reference_size_raw"
MATCHED_ROOT="${ROOT}/panels/reference_size_matched"
SENTINEL_ROOT="${ROOT}/panels/reference_size_sentinel"
mkdir -p "${KEEP_ROOT}" "${RAW_ROOT}" "${MATCHED_ROOT}" "${SENTINEL_ROOT}" \
  "${ROOT}/work/reference_size" "${ROOT}/logs/reference_size" "${ROOT}/qc"

python3 "${ROOT}/scripts/generate_reference_downsamples.py" \
  --keep "${ROOT}/input/chinese_unrelated.keep" \
  --output "${KEEP_ROOT}" \
  --sizes 500 2500 10000 \
  --replicates 5 \
  --simulation-reference-size 10000

run_ldscore() {
  local bfile="$1"
  local keep="$2"
  local output="$3"
  "${PYTHON_BIN}" "${LDSC}" \
    --bfile "${bfile}" \
    --keep "${keep}" \
    --l2 \
    --ld-wind-cm 1 \
    --out "${output}"
}

for size in 500 2500 10000; do
  mkdir -p "${RAW_ROOT}/chinese_n${size}"
done
for replicate in 2 3 4 5; do
  for size in 500 2500 10000; do
    mkdir -p "${SENTINEL_ROOT}/chinese_n${size}_r${replicate}"
  done
done

for chromosome in $(seq 1 22); do
  work="${ROOT}/work/reference_size/chr${chromosome}"
  mkdir -p "${work}"
  shared="${ROOT}/work/panel_build/chr${chromosome}/shared.snplist"
  test -s "${shared}"
  "${PLINK}" \
    --bfile "${GENO}" \
    --chr "${chromosome}" \
    --keep "${ROOT}/input/chinese_unrelated.keep" \
    --extract "${shared}" \
    --make-bed \
    --out "${work}/chinese"
  python3 "${ROOT}/scripts/interpolate_bim_genetic_map.py" \
    --bim "${work}/chinese.bim" \
    --map "${MAP_ROOT}/plink.chr${chromosome}.GRCh38.map" \
    --output "${work}/chinese_cm.bim" \
    --report "${ROOT}/qc/reference_size_chr${chromosome}_map.json"
  ln -sf "${work}/chinese.bed" "${work}/chinese_cm.bed"
  ln -sf "${work}/chinese.fam" "${work}/chinese_cm.fam"

  pids=()
  for size in 500 2500 10000; do
    run_ldscore \
      "${work}/chinese_cm" \
      "${KEEP_ROOT}/chinese_n${size}_r1.keep" \
      "${RAW_ROOT}/chinese_n${size}/${chromosome}" \
      > "${ROOT}/logs/reference_size/chr${chromosome}.n${size}.r1.log" 2>&1 &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do wait "${pid}"; done

  if [[ "${chromosome}" == 1 || "${chromosome}" == 6 || "${chromosome}" == 22 ]]; then
    for replicate in 2 3 4 5; do
      pids=()
      for size in 500 2500 10000; do
        run_ldscore \
          "${work}/chinese_cm" \
          "${KEEP_ROOT}/chinese_n${size}_r${replicate}.keep" \
          "${SENTINEL_ROOT}/chinese_n${size}_r${replicate}/${chromosome}" \
          > "${ROOT}/logs/reference_size/chr${chromosome}.n${size}.r${replicate}.log" 2>&1 &
        pids+=("$!")
      done
      for pid in "${pids[@]}"; do wait "${pid}"; done
    done
  fi
  rm -f "${work}"/*.bed "${work}"/*.bim "${work}"/*.fam
done

python3 "${ROOT}/scripts/match_ldscore_panels.py" \
  --panel "1kgeas=${ROOT}/panels/1kgeas_matched" \
  --panel "chinese_full=${ROOT}/panels/chinese_standard" \
  --panel "chinese_n10000=${RAW_ROOT}/chinese_n10000" \
  --panel "chinese_n2500=${RAW_ROOT}/chinese_n2500" \
  --panel "chinese_n500=${RAW_ROOT}/chinese_n500" \
  --output-root "${MATCHED_ROOT}" \
  --report "${ROOT}/qc/reference_size_panel_matching.tsv"

python3 "${ROOT}/scripts/summarise_reference_size_ldscores.py" "${ROOT}"

find "${RAW_ROOT}" "${MATCHED_ROOT}" "${SENTINEL_ROOT}" -type f -print0 \
  | sort -z | xargs -0 sha256sum > "${ROOT}/qc/reference_size_checksums.sha256"
date --iso-8601=seconds > "${ROOT}/qc/reference_size_panels_complete.timestamp"
