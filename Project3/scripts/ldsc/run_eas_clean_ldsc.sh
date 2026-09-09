#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:"
  echo "  bash run_eas_clean_ldsc.sh LDSC_DIR REF_DIR SCZ_INPUT MDD_MUNGED EDU_MUNGED G_MUNGED OUT_DIR"
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
EDU_MUNGED="$5"
G_MUNGED="$6"
OUT_DIR="$7"

LDSC="$LDSC_DIR/ldsc.py"
MUNGE="$LDSC_DIR/munge_sumstats.py"
for path in "$LDSC" "$MUNGE" "$REF_DIR/w_hm3.snplist" "$SCZ_INPUT" "$MDD_MUNGED" "$EDU_MUNGED" "$G_MUNGED"; do
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
  --out "$MUNGE_DIR/SCZ_EAS" \
  --merge-alleles "$REF_DIR/w_hm3.snplist" \
  --snp SNP --a1 A1 --a2 A2 --N-col N \
  --chunksize 500000 --signed-sumstats Z,0

SCZ_MUNGED="$MUNGE_DIR/SCZ_EAS.sumstats.gz"

python "$LDSC" --h2 "$SCZ_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --samp-prev 0.4503 --pop-prev 0.01 --out "$H2_DIR/SCZ_EAS_h2_liab"
python "$LDSC" --h2 "$MDD_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --samp-prev 0.1278 --pop-prev 0.15 --out "$H2_DIR/MDD_EAS_h2_liab"
python "$LDSC" --h2 "$EDU_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --out "$H2_DIR/EDU_EAS_h2_obs"
python "$LDSC" --h2 "$G_MUNGED" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
  --out "$H2_DIR/HEL10k_EAS_h2_obs"

declare -a PAIRS=(
  "SCZ_MDD|$SCZ_MUNGED,$MDD_MUNGED"
  "SCZ_EDU|$SCZ_MUNGED,$EDU_MUNGED"
  "SCZ_HEL10k|$SCZ_MUNGED,$G_MUNGED"
  "MDD_EDU|$MDD_MUNGED,$EDU_MUNGED"
  "MDD_HEL10k|$MDD_MUNGED,$G_MUNGED"
  "EDU_HEL10k|$EDU_MUNGED,$G_MUNGED"
)

for item in "${PAIRS[@]}"; do
  name="${item%%|*}"
  pair="${item#*|}"
  python "$LDSC" --rg "$pair" --ref-ld-chr "$REF_DIR/" --w-ld-chr "$REF_DIR/" \
    --out "$RG_DIR/$name"
done

echo "Completed clean EAS LDSC workflow: $OUT_DIR"
