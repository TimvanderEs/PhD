#!/usr/bin/env bash
set -euo pipefail

work=${1:-/home/ec2-user/HELIOS_validation/edu_supergnova_20260921}
plink=${PLINK2:-/home/ec2-user/software/plink/plink2}
python=${SUPERGNOVA_PYTHON:-/home/ec2-user/software/SUPERGNOVA/.venv/bin/python}
combined="$work/reference/eas_maf05_combined"
map_dir="$work/resources/genetic_map/extracted"

mkdir -p "$work/reference/by_chr" "$work/provenance/genetic_map_reports" "$work/logs/reference_by_chr"

for chromosome in $(seq 1 22); do
    output="$work/reference/by_chr/eas_chr${chromosome}_maf05"
    "$plink" \
        --bfile "$combined" \
        --chr "$chromosome" \
        --make-bed \
        --out "$output" \
        > "$work/logs/reference_by_chr/chr${chromosome}.log" 2>&1

    "$python" "$work/scripts/interpolate_bim_genetic_map.py" \
        --bim "${output}.bim" \
        --map "$map_dir/plink.chr${chromosome}.GRCh37.map" \
        --output "${output}.bim.with_map" \
        --report "$work/provenance/genetic_map_reports/chr${chromosome}.json"
    mv "${output}.bim.with_map" "${output}.bim"
done

echo "Prepared 22 chromosome-specific MAF >= 0.05 EAS reference files with GRCh37 cM coordinates."
