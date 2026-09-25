#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-/home/ec2-user/HELIOS_validation/genomic_sem_components_20260924}
RAW=${ROOT}/raw
PREPARED=${ROOT}/prepared
MUNGED=${ROOT}/munged
LOGS=${ROOT}/logs
RESULTS=${ROOT}/results
SCRIPTS=${ROOT}/scripts

mkdir -p "${PREPARED}" "${MUNGED}" "${LOGS}" "${RESULTS}"

ROSETTA=${PREPARED}/helios_component_rosetta.tsv.gz
if [[ ! -s "${ROSETTA}" ]]; then
  python3 "${SCRIPTS}/build_component_rosetta.py" \
    --rsid-to-marker /home/ec2-user/COV_LDSC_paper/rsid_to_marker.map \
    --retained-markers /home/ec2-user/COV_LDSC_paper/sg_master_bp38_unordered_allele.map \
    --build-map /home/ec2-user/input_files/misc/new_header_build38_merge.txt \
    --output "${ROSETTA}" \
    >"${LOGS}/build_rosetta.stdout.log" \
    2>"${LOGS}/build_rosetta.stderr.log"
fi

declare -A TENK_INPUTS=(
  [reaction_time]="${RAW}/10k/HELIOS_COG_All_DC7R6_React_Avg_filt2studies_hetp05_AFdiff05.summary.txt.gz"
  [stroop_box]="${RAW}/10k/HELIOS_COG_All_DC7R8_Stroopbox_Avg_filt2studies_hetp05_AFdiff05.summary.txt.gz"
  [stroop_ink]="${RAW}/10k/HELIOS_COG_All_DC7R10_Stroopink_Avg_filt2studies_hetp05_AFdiff05.summary.txt.gz"
  [quiz_score]="${RAW}/10k/HELIOS_COG_All_DC7R5_Quiz_Score_filt2studies_hetp05_AFdiff05.summary.txt.gz"
  [working_memory]="${RAW}/10k/HELIOS_COG_All_DC7R12_Wm_Score_filt2studies_hetp05_AFdiff05.summary.txt.gz"
  [pairing_guesses]="${RAW}/10k/HELIOS_COG_All_DC7R4_Pairing7_Guesses_filt2studies_hetp05_AFdiff05.summary.txt.gz"
)

declare -A TWENTYTWO_INPUTS=(
  [reaction_time]="${RAW}/22k/dc7r6_react_avg.tar"
  [stroop_box]="${RAW}/22k/dc7r8_stroopbox_avg.tar"
  [stroop_ink]="/home/ec2-user/input_files/sumstats/dc7r10_stroopink_avg.tar"
  [quiz_score]="${RAW}/22k/dc7r5_quiz_score.tar"
  [working_memory]="/home/ec2-user/input_files/sumstats/dc7r12_wm_score.tar"
  [pairing_guesses]="${RAW}/22k/dc7r4_pairing7_guesses.tar"
)

sign_for_trait() {
  case "$1" in
    reaction_time|stroop_box|stroop_ink|pairing_guesses) echo -1 ;;
    quiz_score|working_memory) echo 1 ;;
    *) return 2 ;;
  esac
}

prepare_freeze() {
  local freeze=$1
  local -n inputs=$2
  for trait in reaction_time stroop_box stroop_ink quiz_score working_memory pairing_guesses; do
    local output=${PREPARED}/${freeze}_${trait}.ldsc_input.tsv.gz
    if [[ ! -s "${output}" ]]; then
      python3 "${SCRIPTS}/prepare_component_sumstats.py" \
        --input "${inputs[$trait]}" \
        --rosetta "${ROSETTA}" \
        --output "${output}" \
        --sign "$(sign_for_trait "${trait}")" \
        >"${LOGS}/${freeze}_${trait}.prepare.stdout.log" \
        2>"${LOGS}/${freeze}_${trait}.prepare.stderr.log"
    fi
  done
}

prepare_freeze 10k TENK_INPUTS
prepare_freeze 22k TWENTYTWO_INPUTS

LDSC_PYTHON=/home/ec2-user/miniconda3/envs/ldsc/bin/python
MUNGE=/home/ec2-user/software/ldsc/munge_sumstats.py
HM3=/home/ec2-user/input_files/reference/eas_ldscores/w_hm3.snplist

for freeze in 10k 22k; do
  for trait in reaction_time stroop_box stroop_ink quiz_score working_memory pairing_guesses; do
    output=${MUNGED}/${freeze}_${trait}.sumstats.gz
    if [[ ! -s "${output}" ]]; then
      "${LDSC_PYTHON}" "${MUNGE}" \
        --signed-sumstats Z,0 \
        --out "${MUNGED}/${freeze}_${trait}" \
        --merge-alleles "${HM3}" \
        --N-col N \
        --chunksize 500000 \
        --a1 A1 \
        --a2 A2 \
        --snp SNP \
        --sumstats "${PREPARED}/${freeze}_${trait}.ldsc_input.tsv.gz" \
        --p P \
        >"${LOGS}/${freeze}_${trait}.munge.stdout.log" \
        2>"${LOGS}/${freeze}_${trait}.munge.stderr.log"
    fi
  done
  /home/ec2-user/miniconda3/envs/r_env/bin/Rscript \
    "${SCRIPTS}/run_genomic_sem_models.R" \
    "${freeze}" "${MUNGED}" "${RESULTS}/${freeze}" \
    >"${LOGS}/${freeze}.genomic_sem.stdout.log" \
    2>"${LOGS}/${freeze}.genomic_sem.stderr.log"
  /home/ec2-user/miniconda3/envs/r_env/bin/Rscript \
    "${SCRIPTS}/summarize_genomic_sem_outputs.R" \
    "${freeze}" "${RESULTS}/${freeze}"
done

sha256sum "${SCRIPTS}"/* "${PREPARED}"/*.gz "${MUNGED}"/*.sumstats.gz >"${ROOT}/checksums.sha256"
df -h /home/ec2-user >"${ROOT}/filesystem_status.txt"
