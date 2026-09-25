#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <new_output_directory>" >&2
    exit 2
fi

output_directory=$1
script_directory=$(cd "$(dirname "$0")" && pwd)
reference_prefix=/home/ec2-user/input_files/reference/1KGEAS/g1000_eas

if [[ -e "$output_directory/control/control.univ.lava" || -e "$output_directory/corrected/corrected.univ.lava" ]]; then
    echo "Refusing to overwrite an existing validation run: $output_directory" >&2
    exit 1
fi

mkdir -p "$output_directory/control" "$output_directory/corrected" "$output_directory/provenance"

sha256sum \
    "$script_directory/input.info_control_binary.tsv" \
    "$script_directory/input.info_corrected_quantitative.tsv" \
    "$script_directory/sample_overlap.tsv" \
    "$script_directory/run_four_trait_lava.R" \
    "$script_directory/run_control_and_corrected.sh" \
    /home/ec2-user/software/COLOC-reporter/EAS_LAVA.locfile \
    /home/ec2-user/project2/inputfiles/EAS_raw/MDD_EAS_for_LAVA.sumstats.gz \
    /home/ec2-user/project2/inputfiles/EAS_raw/Edu_EAS_for_LAVA.sumstats.gz \
    /home/ec2-user/project2/inputfiles/EAS_raw/SCZ_EAS_for_LAVA.sumstats.gz \
    /home/ec2-user/project2/inputfiles/EAS_raw/HELIOS_10k_for_LAVA.sumstats.gz \
    "${reference_prefix}.bed" \
    "${reference_prefix}.bim" \
    "${reference_prefix}.fam" \
    > "$output_directory/provenance/input_checksums.sha256"

Rscript "$script_directory/run_four_trait_lava.R" \
    "$script_directory/input.info_control_binary.tsv" \
    "$script_directory/sample_overlap.tsv" \
    "$reference_prefix" \
    "$output_directory/control/control" \
    > "$output_directory/control/control.log" 2>&1 &
control_pid=$!

Rscript "$script_directory/run_four_trait_lava.R" \
    "$script_directory/input.info_corrected_quantitative.tsv" \
    "$script_directory/sample_overlap.tsv" \
    "$reference_prefix" \
    "$output_directory/corrected/corrected" \
    > "$output_directory/corrected/corrected.log" 2>&1 &
corrected_pid=$!

echo "Control PID: $control_pid"
echo "Corrected PID: $corrected_pid"
echo "Locus workers per run: 1 (historical execution path)"

control_status=0
corrected_status=0
wait "$control_pid" || control_status=$?
wait "$corrected_pid" || corrected_status=$?

if [[ $control_status -ne 0 || $corrected_status -ne 0 ]]; then
    echo "LAVA run failed: control=$control_status corrected=$corrected_status" >&2
    exit 1
fi

echo "Both LAVA runs completed successfully"
