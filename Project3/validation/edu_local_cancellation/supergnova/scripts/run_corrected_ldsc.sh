#!/usr/bin/env bash
set -euo pipefail

work=${1:-/home/ec2-user/HELIOS_validation/edu_supergnova_20260921}
python=${LDSC_PYTHON:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}
ldsc=${LDSC_SCRIPT:-/home/ec2-user/software/ldsc/ldsc.py}
reference=${EAS_LDSCORES:-/home/ec2-user/input_files/reference/eas_ldscores/}

mkdir -p "$work/results/ldsc_corrected"
"$python" "$ldsc" \
    --rg "$work/input/helios_g_n5592.supergnova.sumstats.gz,$work/input/koges_edu_n71678.supergnova.sumstats.gz" \
    --ref-ld-chr "$reference" \
    --w-ld-chr "$reference" \
    --out "$work/results/ldsc_corrected/helios_g_koges_edu"
