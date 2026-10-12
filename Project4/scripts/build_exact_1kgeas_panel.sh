#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
RAW_1KG="${RAW_1KG:-/home/ec2-user/input_files/reference/1KGEAS/g1000_eas}"
SOURCE_LDSCORES="${SOURCE_LDSCORES:-/home/ec2-user/input_files/reference/eas_ldscores}"
CROSSWALK="${CROSSWALK:-/home/ec2-user/COV_LDSC_paper/rsid_to_marker.map}"
PLINK="${PLINK:-/home/ec2-user/software/plink/plink}"
PLINK2="${PLINK2:-/home/ec2-user/software/plink/plink2}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"
LDSC="${LDSC:-${ROOT}/software/cov-ldsc-vectorized/ldsc.py}"
MAP_ROOT="${MAP_ROOT:-${ROOT}/maps/GRCh38/no_chr_in_chrom_field}"
MAX_CHROMOSOMES="${MAX_CHROMOSOMES:-4}"

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for suffix in bed bim fam; do test -s "${RAW_1KG}.${suffix}"; done
test -s "${CROSSWALK}"
mapping_root="${ROOT}/work/exact_1kgeas_mappings"
work_root="${ROOT}/work/exact_1kgeas"
panel_root="${ROOT}/panels/1kgeas_matched"
archive_root="${ROOT}/archive/row_subset_1kgeas_panel"
mkdir -p "${mapping_root}" "${work_root}" "${panel_root}" "${ROOT}/logs/exact_1kgeas" "${ROOT}/qc"
printf 'chromosome\trequested_markers\textracted_raw_markers\tmissing_from_raw_genotypes\n' \
  > "${ROOT}/qc/exact_1kgeas_build_status.tsv"

python3 "${ROOT}/scripts/prepare_exact_1kgeas_mappings.py" \
  --crosswalk "${CROSSWALK}" \
  --source-ldscore-root "${SOURCE_LDSCORES}" \
  --panel-work-root "${ROOT}/work/panel_build" \
  --output-root "${mapping_root}" \
  --report "${ROOT}/qc/exact_1kgeas_mapping.tsv"

# Retain the row-subset comparator for auditability, then replace the active
# comparator with LD scores recomputed from raw 1KG EAS genotypes on exactly
# the same SNP set used by the Singapore panels.
if [[ ! -d "${archive_root}" ]]; then
  mkdir -p "${archive_root}"
  cp -a "${panel_root}/." "${archive_root}/"
fi

process_chromosome() {
  local chromosome="$1"
  local work="${work_root}/chr${chromosome}"
  mkdir -p "${work}"
  "${PLINK}" \
    --bfile "${RAW_1KG}" \
    --chr "${chromosome}" \
    --extract "${mapping_root}/chr${chromosome}.extract" \
    --make-bed \
    --out "${work}/raw"
  python3 "${ROOT}/scripts/rewrite_1kgeas_bim.py" \
    --bim "${work}/raw.bim" \
    --mapping "${mapping_root}/chr${chromosome}.rsid_to_marker.tsv" \
    --output "${work}/renamed.bim" \
    > "${ROOT}/logs/exact_1kgeas/chr${chromosome}.rewrite.log"
  ln -sf "${work}/raw.bed" "${work}/renamed.bed"
  ln -sf "${work}/raw.fam" "${work}/renamed.fam"
  # GRCh37-to-GRCh38 rsID translation can change local marker order. Reorder
  # the BED together with the rewritten BIM before applying the GRCh38 map.
  "${PLINK2}" \
    --bfile "${work}/renamed" \
    --sort-vars \
    --make-pgen \
    --threads 1 \
    --out "${work}/sorted_pgen" \
    > "${ROOT}/logs/exact_1kgeas/chr${chromosome}.sort.log" 2>&1
  "${PLINK2}" \
    --pfile "${work}/sorted_pgen" \
    --make-bed \
    --threads 1 \
    --out "${work}/sorted" \
    > "${ROOT}/logs/exact_1kgeas/chr${chromosome}.bed_conversion.log" 2>&1
  python3 "${ROOT}/scripts/interpolate_bim_genetic_map.py" \
    --bim "${work}/sorted.bim" \
    --map "${MAP_ROOT}/plink.chr${chromosome}.GRCh38.map" \
    --output "${work}/exact_cm.bim" \
    --report "${ROOT}/qc/exact_1kgeas_chr${chromosome}_map.json"
  ln -sf "${work}/sorted.bed" "${work}/exact_cm.bed"
  ln -sf "${work}/sorted.fam" "${work}/exact_cm.fam"
  "${PYTHON_BIN}" "${LDSC}" \
    --bfile "${work}/exact_cm" \
    --l2 \
    --ld-wind-cm 1 \
    --out "${panel_root}/${chromosome}" \
    > "${ROOT}/logs/exact_1kgeas/chr${chromosome}.stdout.log" 2>&1
  requested=$(wc -l < "${mapping_root}/chr${chromosome}.extract")
  expected=$(wc -l < "${work}/sorted.bim")
  observed=$(($(zcat "${panel_root}/${chromosome}.l2.ldscore.gz" | wc -l) - 1))
  if [[ "${expected}" -ne "${observed}" ]]; then
    echo "Exact 1KG EAS row mismatch on chromosome ${chromosome}: ${observed} != ${expected}" >&2
    return 1
  fi
  printf '%s\n' "${expected}" > "${panel_root}/${chromosome}.l2.M"
  printf '%s\n' "${expected}" > "${panel_root}/${chromosome}.l2.M_5_50"
  printf '%s\t%s\t%s\t%s\n' \
    "${chromosome}" "${requested}" "${expected}" "$((requested - expected))" \
    >> "${ROOT}/qc/exact_1kgeas_build_status.tsv"
  rm -f "${work}"/*.bed "${work}"/*.bim "${work}"/*.fam \
    "${work}"/*.pgen "${work}"/*.pvar "${work}"/*.psam
}

wait_batch() {
  local status=0
  local pid
  for pid in "$@"; do wait "${pid}" || status=1; done
  return "${status}"
}

pids=()
for chromosome in $(seq 1 22); do
  process_chromosome "${chromosome}" &
  pids+=("$!")
  if [[ "${#pids[@]}" -ge "${MAX_CHROMOSOMES}" ]]; then
    wait_batch "${pids[@]}"
    pids=()
  fi
done
if [[ "${#pids[@]}" -gt 0 ]]; then wait_batch "${pids[@]}"; fi

# The public LD-score files and raw 1KG genotype archive differ by a small
# number of release-specific markers. Re-intersect all active panels after the
# raw-genotype build so their regression rows, LD-calculation variants, order,
# and M values are genuinely identical.
matched_root="${ROOT}/work/final_five_panel_intersection_$(date -u +%Y%m%dT%H%M%SZ)"
python3 "${ROOT}/scripts/match_ldscore_panels.py" \
  --panel "1kgeas_matched=${ROOT}/panels/1kgeas_matched" \
  --panel "chinese_standard=${ROOT}/panels/chinese_standard" \
  --panel "chinese_cov20pc=${ROOT}/panels/chinese_cov20pc" \
  --panel "pooled_standard=${ROOT}/panels/pooled_standard" \
  --panel "pooled_cov20pc=${ROOT}/panels/pooled_cov20pc" \
  --output-root "${matched_root}" \
  --report "${ROOT}/qc/final_five_panel_intersection.tsv"
for panel in 1kgeas_matched chinese_standard chinese_cov20pc pooled_standard pooled_cov20pc; do
  cp -a "${matched_root}/${panel}/." "${ROOT}/panels/${panel}/"
done

date --iso-8601=seconds > "${ROOT}/qc/exact_1kgeas_complete.timestamp"
