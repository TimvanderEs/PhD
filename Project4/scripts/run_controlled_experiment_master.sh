#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260922/controlled_reference_experiment}"
PANEL_BUILD_PID="${PANEL_BUILD_PID:?Set PANEL_BUILD_PID to the active build process}"
MASTER_LOG="${ROOT}/logs/controlled_experiment_master.log"

log() {
  printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$*" | tee -a "${MASTER_LOG}"
}

log "Waiting for controlled panel build PID ${PANEL_BUILD_PID}"
while ps -p "${PANEL_BUILD_PID}" -o cmd= 2>/dev/null | grep -q 'build_controlled_ldscore_panels.sh'; do
  sleep 30
done
test -s "${ROOT}/qc/panels_complete.timestamp"
log "Recomputing the 1000 Genomes EAS comparator from raw genotypes on the exact shared SNP set"
bash "${ROOT}/scripts/build_exact_1kgeas_panel.sh" \
  > "${ROOT}/logs/build_exact_1kgeas_panel.log" 2>&1
test -s "${ROOT}/qc/exact_1kgeas_complete.timestamp"
python3 "${ROOT}/scripts/align_ldscore_panel_order.py" "${ROOT}"
find "${ROOT}/panels" -maxdepth 2 -type f -print0 | sort -z | xargs -0 sha256sum \
  > "${ROOT}/qc/panel_checksums.sha256"
completed_chromosomes=$(awk -F '\t' 'NR > 1 && $5 == "complete" {count++} END {print count + 0}' "${ROOT}/qc/panel_build_status.tsv")
if [[ "${completed_chromosomes}" -ne 22 ]]; then
  log "Panel validation failed: ${completed_chromosomes}/22 chromosomes complete"
  exit 1
fi
log "Controlled five-panel build validated"

log "Starting controlled real-GWAS LDSC and reference-size panel build"
bash "${ROOT}/scripts/build_reference_size_panels.sh" \
  > "${ROOT}/logs/build_reference_size_panels.log" 2>&1 &
reference_pid=$!
controlled_status=0
bash "${ROOT}/scripts/run_controlled_panel_ldsc.sh" \
  > "${ROOT}/logs/run_controlled_panel_ldsc.log" 2>&1 || controlled_status=1
reference_status=0
wait "${reference_pid}" || reference_status=1
if [[ "${controlled_status}" -ne 0 || "${reference_status}" -ne 0 ]]; then
  log "Concurrent stage failed: controlled_status=${controlled_status}; reference_status=${reference_status}"
  exit 1
fi
python3 "${ROOT}/scripts/summarise_paired_h2_jackknife.py" "${ROOT}"
python3 "${ROOT}/scripts/summarise_controlled_rg.py" "${ROOT}"
log "Controlled real-GWAS estimates, paired jackknife, and reference-size panels complete"

log "Starting real-GWAS reference-size h2 sensitivity"
bash "${ROOT}/scripts/run_reference_size_h2.sh" \
  > "${ROOT}/logs/run_reference_size_h2.log" 2>&1
log "Real-GWAS reference-size h2 sensitivity complete"

log "Starting held-out genotype simulation"
bash "${ROOT}/scripts/run_heldout_simulation.sh" \
  > "${ROOT}/logs/run_heldout_simulation.log" 2>&1
log "Held-out genotype simulation complete"

{
  echo "completed_utc=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "kernel=$(uname -srmo)"
  echo "ldsc_python=$(/home/ec2-user/miniconda3/envs/ldsc/bin/python --version 2>&1)"
  echo "plink=$(/home/ec2-user/software/plink/plink --version 2>&1 | head -1)"
  echo "plink2=$(/home/ec2-user/software/plink/plink2 --version 2>&1 | head -1)"
  df -h "${ROOT}"
} > "${ROOT}/qc/final_session_info.txt"
find "${ROOT}/scripts" -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum \
  > "${ROOT}/qc/final_script_checksums.sha256"
find "${ROOT}/results" -type f -print0 | sort -z | xargs -0 sha256sum \
  > "${ROOT}/qc/final_result_checksums.sha256"
date --iso-8601=seconds > "${ROOT}/qc/controlled_experiment_complete.timestamp"
log "All controlled reference analyses complete"
