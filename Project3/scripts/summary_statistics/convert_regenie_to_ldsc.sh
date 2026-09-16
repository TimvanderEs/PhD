#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input.regenie> <output.sumstats.tsv.gz>" >&2
    exit 2
fi

input_file=$1
output_file=$2
temporary_output="${output_file%.gz}"

if [[ "${temporary_output}" == "${output_file}" ]]; then
    echo "Output filename must end in .gz" >&2
    exit 2
fi

# Required REGENIE columns:
# CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
awk 'BEGIN {
         OFS="\t"
     }
     NR == 1 {
         print "SNP", "CHR", "BP", "A1", "A2", "N", "BETA", "SE", "P"
         next
     }
     {
         if ($12 == "NA" || $12 == "nan" || $12 == "NaN") next

         # REGENIE reports -log10(P) in column 12. ALLELE1 is the
         # effect allele for the BETA in column 9.
         p = exp(-$12 * log(10))
         print $3, $1, $2, $5, $4, $7, $9, $10, p
     }' "${input_file}" > "${temporary_output}"

gzip -f "${temporary_output}"

if [[ "${temporary_output}.gz" != "${output_file}" ]]; then
    mv "${temporary_output}.gz" "${output_file}"
fi

echo "Wrote ${output_file}"
