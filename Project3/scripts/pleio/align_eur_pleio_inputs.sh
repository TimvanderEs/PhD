#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash align_eur_pleio_inputs.sh [--validate-only] \
    MDD_TAB.gz EA_TAB.gz SCZ_TAB.gz CF_TAB.gz OUT_DIR

Expected source columns:
  MDD/EA/CF: SNP=1, A1=4, A2=5, Z=6, N=7
  SCZ:         SNP=1, A1=4, A2=5, N=11, Z=12
EOF
}

VALIDATE_ONLY=0
if [[ ${1:-} == "--validate-only" ]]; then
  VALIDATE_ONLY=1
  shift
fi
if [[ ${1:-} == "--help" ]]; then
  usage
  exit 0
fi
if [[ $# -ne 5 ]]; then
  usage >&2
  exit 2
fi

MDD="$1"
EA="$2"
SCZ="$3"
CF="$4"
OUT_DIR="$5"
for path in "$MDD" "$EA" "$SCZ" "$CF"; do
  if [[ ! -s "$path" ]]; then
    echo "ERROR: missing input: $path" >&2
    exit 1
  fi
  gzip -t "$path"
done

first_line() {
  local path="$1"
  local line
  set +o pipefail
  line="$(gzip -cd "$path" | head -n 1)"
  set -o pipefail
  printf '%s\n' "$line"
}

check_header() {
  local label="$1"
  local path="$2"
  local expected="$3"
  local observed
  observed="$(first_line "$path")"
  if [[ "$observed" != "$expected" ]]; then
    echo "ERROR: unexpected $label header" >&2
    echo "Expected: $expected" >&2
    echo "Observed: $observed" >&2
    exit 1
  fi
}

standard_header=$'SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP'
scz_header=$'SNP\tCHR\tBP\tA1\tA2\tMAF\tINFO\tBETA\tSE\tP\tN\tZ'
check_header MDD "$MDD" "$standard_header"
check_header EA "$EA" "$standard_header"
check_header CF "$CF" "$standard_header"
check_header SCZ "$SCZ" "$scz_header"

echo "Validated EUR PLEIO source headers."
if [[ "$VALIDATE_ONLY" -eq 1 ]]; then
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARMONIZE_AWK="$SCRIPT_DIR/harmonize_to_anchor.awk"
HEADER_AWK="$SCRIPT_DIR/add_header_and_5cols.awk"
for path in "$HARMONIZE_AWK" "$HEADER_AWK"; do
  [[ -s "$path" ]] || { echo "ERROR: missing helper: $path" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/project3_eur_pleio_align.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT
STD_DIR="$WORK_DIR/std"
mkdir -p "$STD_DIR"

gzip -cd "$MDD" | awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","A1","A2","Z","N";next} {print $1,$4,$5,$6,$7}' | gzip -c > "$STD_DIR/MDD.full.std.tsv.gz"
gzip -cd "$EA" | awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","A1","A2","Z","N";next} {print $1,$4,$5,$6,$7}' | gzip -c > "$STD_DIR/EA.full.std.tsv.gz"
gzip -cd "$SCZ" | awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","A1","A2","Z","N";next} {print $1,$4,$5,$12,$11}' | gzip -c > "$STD_DIR/SCZ.full.std.tsv.gz"
gzip -cd "$CF" | awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","A1","A2","Z","N";next} {print $1,$4,$5,$6,$7}' | gzip -c > "$STD_DIR/CF.full.std.tsv.gz"

for trait in MDD EA SCZ CF; do
  gzip -cd "$STD_DIR/${trait}.full.std.tsv.gz" |
    awk 'NR>1 {print $1}' |
    LC_ALL=C sort -u > "$WORK_DIR/${trait}.snps"
done

comm -12 "$WORK_DIR/MDD.snps" "$WORK_DIR/EA.snps" > "$WORK_DIR/common.1"
comm -12 "$WORK_DIR/common.1" "$WORK_DIR/SCZ.snps" > "$WORK_DIR/common.2"
comm -12 "$WORK_DIR/common.2" "$WORK_DIR/CF.snps" > "$WORK_DIR/common.snps"

gzip -cd "$STD_DIR/EA.full.std.tsv.gz" |
  awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,toupper($2),toupper($3)}' |
  LC_ALL=C sort -k1,1 > "$WORK_DIR/ea.alleles.tsv"
join -t $'\t' -1 1 -2 1 "$WORK_DIR/common.snps" "$WORK_DIR/ea.alleles.tsv" \
  > "$WORK_DIR/anchor_alleles.tsv"

harmonize() {
  local anchor="$1"
  local input="$2"
  local no_header="$3"
  gzip -cd "$input" | awk -f "$HARMONIZE_AWK" "$anchor" - | gzip -c > "$no_header"
}

harmonize "$WORK_DIR/anchor_alleles.tsv" "$STD_DIR/MDD.full.std.tsv.gz" "$WORK_DIR/MDD.no_header.gz"
harmonize "$WORK_DIR/ea.alleles.tsv" "$STD_DIR/EA.full.std.tsv.gz" "$WORK_DIR/EA.no_header.gz"
harmonize "$WORK_DIR/anchor_alleles.tsv" "$STD_DIR/SCZ.full.std.tsv.gz" "$WORK_DIR/SCZ.no_header.gz"
harmonize "$WORK_DIR/anchor_alleles.tsv" "$STD_DIR/CF.full.std.tsv.gz" "$WORK_DIR/CF.no_header.gz"

for trait in MDD EA SCZ CF; do
  {
    printf 'SNP\tA1\tA2\tZ\tN\n'
    gzip -cd "$WORK_DIR/${trait}.no_header.gz" | awk -f "$HEADER_AWK"
  } | gzip -c > "$OUT_DIR/${trait}.tmp.gz"
done

mv "$OUT_DIR/MDD.tmp.gz" "$OUT_DIR/MDD_EUR.aligned.tsv.gz"
mv "$OUT_DIR/EA.tmp.gz" "$OUT_DIR/EA_EUR.aligned.tsv.gz"
mv "$OUT_DIR/SCZ.tmp.gz" "$OUT_DIR/SCZ_EUR.aligned.tsv.gz"
mv "$OUT_DIR/CF.tmp.gz" "$OUT_DIR/CF_EUR.aligned.tsv.gz"

echo "Completed EUR PLEIO alignment: $OUT_DIR"
