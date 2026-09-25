#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <completed_run_root> <analysis_output_root>" >&2
    exit 2
fi

run_root=$1
analysis_root=$2
config_root=/home/ec2-user/HELIOS_validation/edu_local_cancellation_20260918/config
archived_root=/home/ec2-user/software/COLOC-reporter/results_EAS

source /home/ec2-user/miniconda3/etc/profile.d/conda.sh
conda activate coloc_reporter

mkdir -p "$analysis_root/provenance" "$analysis_root/final"

Rscript "$config_root/compare_lava_control.R" \
    "$archived_root/EAS_MDD_Edu_SCZ_HELIOS.univ.lava" \
    "$archived_root/EAS_MDD_Edu_SCZ_HELIOS.bivar.lava" \
    "$run_root/control/control.univ.lava" \
    "$run_root/control/control.bivar.lava" \
    "$analysis_root/provenance/historical_control_comparison.tsv" \
    > "$analysis_root/provenance/historical_control_comparison.log" 2>&1

Rscript "$config_root/compare_control_corrected.R" \
    "$run_root/control/control.univ.lava" \
    "$run_root/control/control.bivar.lava" \
    "$run_root/corrected/corrected.univ.lava" \
    "$run_root/corrected/corrected.bivar.lava" \
    "$analysis_root/provenance/control_corrected_impact.tsv" \
    "$analysis_root/provenance/target_pair_eligibility_comparison.tsv" \
    > "$analysis_root/provenance/control_corrected_impact.log" 2>&1

Rscript "$config_root/analyze_local_cancellation.R" \
    "$run_root/corrected/corrected.univ.lava" \
    "$run_root/corrected/corrected.bivar.lava" \
    "$config_root/EAS_LAVA.locfile" \
    "$config_root/global_ldsc_summary.tsv" \
    "$analysis_root" \
    > "$analysis_root/final/postprocessing.log" 2>&1

(
    cd "$analysis_root"
    find . -type f ! -path "./provenance/analysis_outputs.sha256" -print0 \
        | sort -z \
        | xargs -0 sha256sum \
        > "./provenance/analysis_outputs.sha256"
)

echo "Historical-control validation and corrected-output analysis completed"
