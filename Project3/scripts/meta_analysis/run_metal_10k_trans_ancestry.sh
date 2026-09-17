#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <metal> <Malay-regenie> <Chinese-regenie> <Indian-regenie> <genomic-control:ON|OFF> <output-prefix>" >&2
  exit 2
fi

metal=$1
malay_result=$2
chinese_result=$3
indian_result=$4
genomic_control=$(printf '%s' "$5" | tr '[:lower:]' '[:upper:]')
output_prefix=$6

if [[ "${genomic_control}" != "ON" && "${genomic_control}" != "OFF" ]]; then
  echo "Genomic control must be ON or OFF" >&2
  exit 2
fi

for input_file in "${malay_result}" "${chinese_result}" "${indian_result}"; do
  if [[ ! -f "${input_file}" ]]; then
    echo "Missing input: ${input_file}" >&2
    exit 2
  fi
done

mkdir -p "$(dirname "${output_prefix}")"
config=$(mktemp "${TMPDIR:-/tmp}/helios_metal.XXXXXX")
trap 'rm -f "${config}"' EXIT

{
  printf 'GENOMICCONTROL %s\n' "${genomic_control}"
  printf '%s\n' \
    'SCHEME STDERR' \
    'STDERR_PRINT_PRECISION 8' \
    'AVERAGEFREQ ON' \
    'MINMAXFREQ ON' \
    'MARKER ID' \
    'ALLELE ALLELE1 ALLELE0' \
    'EFFECT BETA' \
    'STDERR SE' \
    'FREQLABEL A1FREQ' \
    'WEIGHTLABEL N' \
    'SEPARATOR WHITESPACE'
  printf 'PROCESS %s\n' "${malay_result}" "${chinese_result}" "${indian_result}"
  printf 'OUTFILE %s .sumstats.txt\n' "${output_prefix}"
  printf '%s\n' 'ANALYZE HETEROGENEITY' 'QUIT'
} > "${config}"

"${metal}" "${config}"
