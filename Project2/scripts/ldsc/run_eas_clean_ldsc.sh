#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:"
  echo "  bash run_eas_clean_ldsc.sh LDSC_DIR REF_DIR SCZ_INPUT MDD_MUNGED EA_MUNGED CF_MUNGED OUT_DIR"
  echo
  echo "LDSC_DIR must contain ldsc.py and munge_sumstats.py."
  echo "REF_DIR must contain EAS LD-score files and w_hm3.snplist."
}

if [[ ${1:-} == "--help" ]]; then
  usage
  exit 0
fi
if [[ $# -ne 7 ]]; then
  usage >&2
  exit 2
fi

LDSC_DIR="$1"
REF_DIR="$2"
SCZ_INPUT="$3"
MDD_MUNGED="$4"
EA_MUNGED="$5"
CF_MUNGED="$6"
OUT_DIR="$7"

LDSC="$LDSC_DIR/ldsc.py"
MUNGE="$LDSC_DIR/munge_sumstats.py"
for path in "$LDSC" "$MUNGE" "$REF_DIR/w_hm3.snplist" "$SCZ_INPUT" "$MDD_MUNGED" "$EA_MUNGED" "$CF_MUNGED"; do
  if [[ ! -e "$path" ]]; then
    echo "ERROR: missing required path: $path" >&2
    exit 1
  fi
done

MUNGE_DIR="$OUT_DIR/munged"
H2_DIR="$OUT_DIR/h2"
RG_DIR="$OUT_DIR/rg"
mkdir -p "$MUNGE_DIR" "$H2_DIR" "$RG_DIR"

python "$MUNGE" \
  --sumstats "$SCZ_INPUT" \
  --out "$MUNGE_DIR/SCZ" \
  --merge-alleles "$REF_DIR/w_hm3.snplist" \
  --snp SNP --a1 A1 --a2 A2 --N-col N \
  --chunksize 500000 --signed-sumstats Z,0

SCZ_MUNGED="$MUNGE_DIR/SCZ.sumstats.gz"

python "$LDSC" --h2 "$SCZ_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --samp-prev 0.4503 --pop-prev 0.01 --out "$H2_DIR/SCZ_h2_liab"
python "$LDSC" --h2 "$MDD_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --samp-prev 0.1278 --pop-prev 0.15 --out "$H2_DIR/MDD_h2_liab"
python "$LDSC" --h2 "$EA_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --out "$H2_DIR/EA_h2_obs"
python "$LDSC" --h2 "$CF_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --out "$H2_DIR/CF_h2_obs"

declare -a PAIRS=(
  "SCZ_MDD|$SCZ_MUNGED,$MDD_MUNGED"
  "SCZ_EA|$SCZ_MUNGED,$EA_MUNGED"
  "SCZ_CF|$SCZ_MUNGED,$CF_MUNGED"
  "MDD_EA|$MDD_MUNGED,$EA_MUNGED"
  "MDD_CF|$MDD_MUNGED,$CF_MUNGED"
  "EA_CF|$EA_MUNGED,$CF_MUNGED"
)

for item in "${PAIRS[@]}"; do
  name="${item%%|*}"
  pair="${item#*|}"
  python "$LDSC" --rg "$pair" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
    --out "$RG_DIR/$name"
done

echo "Completed clean EAS LDSC workflow: $OUT_DIR"
