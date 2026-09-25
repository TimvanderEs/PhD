#!/usr/bin/env bash
set -euo pipefail

work=${1:-/home/ec2-user/HELIOS_validation/edu_supergnova_20260921}
software=${SUPERGNOVA_DIR:-/home/ec2-user/software/SUPERGNOVA}
threads=${SUPERGNOVA_THREADS:-8}

cd "$software"
/usr/bin/time -v "$software/.venv/bin/python" supergnova.py \
    "$work/input/helios_g_n5592.supergnova.sumstats.gz" \
    "$work/input/koges_edu_n71678.supergnova.sumstats.gz" \
    --N1 5592 \
    --N2 71678 \
    --bfile "$work/reference/by_chr/eas_chr@_maf05" \
    --partition "$work/input/eas_lava_3064_regions.bed" \
    --thread "$threads" \
    --out "$work/results/helios_g_koges_edu.supergnova.txt"
