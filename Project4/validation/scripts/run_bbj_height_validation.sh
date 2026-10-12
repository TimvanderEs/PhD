#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <bbj_height_raw.gz> <w_hm3.snplist> <eas_ldscore_prefix> <new_output_dir>" >&2
    exit 2
fi

raw_sumstats=$1
hm3_snplist=$2
eas_ld_prefix=${3%/}/
output_dir=$4

ldsc_dir=${LDSC_DIR:-/home/ec2-user/software/ldsc}
ldsc_python=${LDSC_PYTHON:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}

for required_file in "$raw_sumstats" "$hm3_snplist" "${eas_ld_prefix}1.l2.ldscore.gz" "${eas_ld_prefix}22.l2.ldscore.gz"; do
    if [[ ! -s "$required_file" ]]; then
        echo "Missing or empty input: $required_file" >&2
        exit 1
    fi
done

if [[ -e "$output_dir/bbj_height_correct_hm3.sumstats.gz" || -e "$output_dir/bbj_height_correct_1kgeas.log" ]]; then
    echo "Refusing to overwrite an existing validation run: $output_dir" >&2
    exit 1
fi

mkdir -p "$output_dir"
sha256sum "$raw_sumstats" "$hm3_snplist" "${eas_ld_prefix}"*.l2.ldscore.gz

"$ldsc_python" "$ldsc_dir/munge_sumstats.py" \
    --sumstats "$raw_sumstats" \
    --out "$output_dir/bbj_height_correct_hm3" \
    --merge-alleles "$hm3_snplist" \
    --snp SNP \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --N 159095 \
    --p P_BOLT_LMM_INF \
    --signed-sumstats BETA,0 \
    --frq A1FREQ \
    --info INFO \
    --chunksize 500000

"$ldsc_python" "$ldsc_dir/ldsc.py" \
    --h2 "$output_dir/bbj_height_correct_hm3.sumstats.gz" \
    --ref-ld-chr "$eas_ld_prefix" \
    --w-ld-chr "$eas_ld_prefix" \
    --out "$output_dir/bbj_height_correct_1kgeas"

