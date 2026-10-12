#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
GENO="${GENO:-${ROOT}/input/sg100k_consensus}"
PLINK="${PLINK:-/home/ec2-user/software/plink/plink}"
COVLDSC_PYTHON="${COVLDSC_PYTHON:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
COVLDSC="${COVLDSC:-${ROOT}/software/cov-ldsc-vectorized/ldsc.py}"
MAP_ROOT="${MAP_ROOT:-${ROOT}/maps/GRCh38/no_chr_in_chrom_field}"
PUBLIC_EAS="${PUBLIC_EAS:-${ROOT}/input/1kgeas_singapore_ids_complete}"
THREADS="${THREADS:-4}"
PCS="${PCS:-20}"
MAX_CHROMOSOMES="${MAX_CHROMOSOMES:-4}"

# Each chromosome launches two independent LDSC jobs. Prevent BLAS from
# oversubscribing the eight-core host inside each process.
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CHINESE_KEEP="${ROOT}/input/chinese_unrelated.keep"
POOLED_KEEP="${ROOT}/input/pooled_unrelated.keep"
CHINESE_COV="${ROOT}/input/pca_chinese_unrelated_${PCS}pc.eigenvec"
POOLED_COV="${ROOT}/input/pca_pooled_unrelated_${PCS}pc.eigenvec"

for required in "${GENO}.bed" "${GENO}.bim" "${GENO}.fam" "${CHINESE_KEEP}" "${POOLED_KEEP}" "${CHINESE_COV}" "${POOLED_COV}"; do
  test -s "${required}"
done
mkdir -p "${ROOT}/work/panel_build" "${ROOT}/logs/panel_build" "${ROOT}/qc/maps"
for panel in chinese_standard chinese_cov20pc pooled_standard pooled_cov20pc 1kgeas_matched; do
  mkdir -p "${ROOT}/panels/${panel}"
done
printf 'chromosome\tshared_variants\tchinese_samples\tpooled_samples\tstatus\n' > "${ROOT}/qc/panel_build_status.tsv"

run_ldscore() {
  local bfile="$1"
  local output="$2"
  local covariates="${3:-}"
  if [[ -n "${covariates}" ]]; then
    "${COVLDSC_PYTHON}" "${COVLDSC}" \
      --bfile "${bfile}" \
      --l2 \
      --ld-wind-cm 1 \
      --cov "${covariates}" \
      --out "${output}"
  else
    "${COVLDSC_PYTHON}" "${COVLDSC}" \
      --bfile "${bfile}" \
      --l2 \
      --ld-wind-cm 1 \
      --out "${output}"
  fi
}

process_chromosome() {
  local chromosome="$1"
  chromosome_work="${ROOT}/work/panel_build/chr${chromosome}"
  mkdir -p "${chromosome_work}"
  map_file="${MAP_ROOT}/plink.chr${chromosome}.GRCh38.map"
  test -s "${map_file}"
  test -s "${PUBLIC_EAS}/${chromosome}.l2.ldscore.gz"

  # HWE is evaluated only within the ancestry-matched Chinese sample. Applying
  # an HWE filter to the pooled multi-ancestry sample would remove valid loci
  # through the Wahlund effect and make the comparison ancestry-dependent.
  "${PLINK}" \
    --bfile "${GENO}" \
    --chr "${chromosome}" \
    --keep "${CHINESE_KEEP}" \
    --snps-only just-acgt \
    --maf 0.01 \
    --geno 0.02 \
    --hwe 1e-6 \
    --write-snplist \
    --threads "${THREADS}" \
    --out "${chromosome_work}/chinese_qc"
  "${PLINK}" \
    --bfile "${GENO}" \
    --chr "${chromosome}" \
    --keep "${POOLED_KEEP}" \
    --snps-only just-acgt \
    --maf 0.01 \
    --geno 0.02 \
    --write-snplist \
    --threads "${THREADS}" \
    --out "${chromosome_work}/pooled_qc"

  awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' \
    "${chromosome_work}/pooled_qc.snplist" \
    "${chromosome_work}/chinese_qc.snplist" \
    > "${chromosome_work}/local_shared.snplist"
  zcat "${PUBLIC_EAS}/${chromosome}.l2.ldscore.gz" \
    | awk 'NR > 1 {print $2}' \
    > "${chromosome_work}/public_eas.snplist"
  awk 'NR==FNR { keep[$1]=1; next } ($1 in keep)' \
    "${chromosome_work}/public_eas.snplist" \
    "${chromosome_work}/local_shared.snplist" \
    > "${chromosome_work}/shared.snplist"
  test -s "${chromosome_work}/shared.snplist"

  # Preserve the public EAS LD scores but restrict them to the same regression
  # variants used to construct every local panel.
  awk 'NR==FNR { keep[$1]=1; next } FNR == 1 || ($2 in keep)' \
    "${chromosome_work}/shared.snplist" \
    <(zcat "${PUBLIC_EAS}/${chromosome}.l2.ldscore.gz") \
    | gzip -c > "${ROOT}/panels/1kgeas_matched/${chromosome}.l2.ldscore.gz"
  matched_rows=$(wc -l < "${chromosome_work}/shared.snplist")
  printf '%s\n' "${matched_rows}" > "${ROOT}/panels/1kgeas_matched/${chromosome}.l2.M"
  printf '%s\n' "${matched_rows}" > "${ROOT}/panels/1kgeas_matched/${chromosome}.l2.M_5_50"

  for sample_set in chinese pooled; do
    keep_variable="${sample_set^^}_KEEP"
    keep_file="${!keep_variable}"
    "${PLINK}" \
      --bfile "${GENO}" \
      --chr "${chromosome}" \
      --keep "${keep_file}" \
      --extract "${chromosome_work}/shared.snplist" \
      --make-bed \
      --threads "${THREADS}" \
      --out "${chromosome_work}/${sample_set}"
    python3 "${ROOT}/scripts/interpolate_bim_genetic_map.py" \
      --bim "${chromosome_work}/${sample_set}.bim" \
      --map "${map_file}" \
      --output "${chromosome_work}/${sample_set}_cm.bim" \
      --report "${ROOT}/qc/maps/chr${chromosome}_${sample_set}.json"
    ln -sf "${chromosome_work}/${sample_set}.bed" "${chromosome_work}/${sample_set}_cm.bed"
    ln -sf "${chromosome_work}/${sample_set}.fam" "${chromosome_work}/${sample_set}_cm.fam"
  done

  run_ldscore "${chromosome_work}/chinese_cm" "${ROOT}/panels/chinese_standard/${chromosome}" \
    > "${ROOT}/logs/panel_build/chr${chromosome}.chinese_standard.stdout.log" 2>&1 &
  pid_chinese_standard=$!
  run_ldscore "${chromosome_work}/pooled_cm" "${ROOT}/panels/pooled_standard/${chromosome}" \
    > "${ROOT}/logs/panel_build/chr${chromosome}.pooled_standard.stdout.log" 2>&1 &
  pid_pooled_standard=$!
  wait_for_batch "${pid_chinese_standard}" "${pid_pooled_standard}"

  run_ldscore "${chromosome_work}/chinese_cm" "${ROOT}/panels/chinese_cov20pc/${chromosome}" "${CHINESE_COV}" \
    > "${ROOT}/logs/panel_build/chr${chromosome}.chinese_cov20pc.stdout.log" 2>&1 &
  pid_chinese_cov=$!
  run_ldscore "${chromosome_work}/pooled_cm" "${ROOT}/panels/pooled_cov20pc/${chromosome}" "${POOLED_COV}" \
    > "${ROOT}/logs/panel_build/chr${chromosome}.pooled_cov20pc.stdout.log" 2>&1 &
  pid_pooled_cov=$!
  wait_for_batch "${pid_chinese_cov}" "${pid_pooled_cov}"

  expected_rows=$(wc -l < "${chromosome_work}/shared.snplist")
  for panel in chinese_standard chinese_cov20pc pooled_standard pooled_cov20pc 1kgeas_matched; do
    gzip -t "${ROOT}/panels/${panel}/${chromosome}.l2.ldscore.gz"
    observed_rows=$(($(zcat "${ROOT}/panels/${panel}/${chromosome}.l2.ldscore.gz" | wc -l) - 1))
    if [[ "${observed_rows}" -ne "${expected_rows}" ]]; then
      echo "Row-count mismatch for chr${chromosome} ${panel}: ${observed_rows} != ${expected_rows}" >&2
      exit 1
    fi
  done

  # All five references represent the identical per-chromosome SNP universe.
  # Set M explicitly to that shared universe for the controlled comparison.
  for panel in chinese_standard chinese_cov20pc pooled_standard pooled_cov20pc; do
    printf '%s\n' "${expected_rows}" > "${ROOT}/panels/${panel}/${chromosome}.l2.M"
    printf '%s\n' "${expected_rows}" > "${ROOT}/panels/${panel}/${chromosome}.l2.M_5_50"
  done

  chinese_samples=$(wc -l < "${chromosome_work}/chinese.fam")
  pooled_samples=$(wc -l < "${chromosome_work}/pooled.fam")
  printf '%s\t%s\t%s\t%s\tcomplete\n' \
    "${chromosome}" "${expected_rows}" "${chinese_samples}" "${pooled_samples}" \
    >> "${ROOT}/qc/panel_build_status.tsv"

  rm -f "${chromosome_work}"/*.bed "${chromosome_work}"/*.bim "${chromosome_work}"/*.fam
}

wait_for_batch() {
  local status=0
  local pid
  for pid in "$@"; do
    wait "${pid}" || status=1
  done
  return "${status}"
}

chromosome_pids=()
for chromosome in $(seq 1 22); do
  process_chromosome "${chromosome}" &
  chromosome_pids+=("$!")
  if [[ "${#chromosome_pids[@]}" -ge "${MAX_CHROMOSOMES}" ]]; then
    wait_for_batch "${chromosome_pids[@]}"
    chromosome_pids=()
  fi
done
if [[ "${#chromosome_pids[@]}" -gt 0 ]]; then
  wait_for_batch "${chromosome_pids[@]}"
fi

python3 "${ROOT}/scripts/align_ldscore_panel_order.py" "${ROOT}"
find "${ROOT}/panels" -type f -print0 | sort -z | xargs -0 sha256sum > "${ROOT}/qc/panel_checksums.sha256"
date --iso-8601=seconds > "${ROOT}/qc/panels_complete.timestamp"
