#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash test_eur_pleio_alignment.sh /path/to/Project3" >&2
  exit 2
fi
PROJECT_DIR="$(cd "$1" && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/project3_align_test.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

printf 'SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\nrs1\t1\t1\tA\tG\t1\t100\t0.3\nrs2\t1\t2\tC\tT\t2\t100\t0.2\nrs3\t1\t3\tA\tT\t3\t100\t0.1\nrs4\t1\t4\tA\tC\t4\t100\t0.05\n' | gzip -c > "$TMP/mdd.gz"
printf 'SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\nrs1\t1\t1\tA\tG\t0.1\t200\t0.3\nrs2\t1\t2\tT\tC\t0.2\t200\t0.2\nrs3\t1\t3\tA\tT\t0.3\t200\t0.1\nrs4\t1\t4\tA\tC\t0.4\t200\t0.05\n' | gzip -c > "$TMP/edu.gz"
printf 'SNP\tCHR\tBP\tA1\tA2\tMAF\tINFO\tBETA\tSE\tP\tN\tZ\nrs1\t1\t1\tA\tG\t0.1\t1\t0\t1\t0.3\t300\t-1\nrs2\t1\t2\tT\tC\t0.1\t1\t0\t1\t0.2\t300\t-2\nrs3\t1\t3\tA\tT\t0.1\t1\t0\t1\t0.1\t300\t-3\nrs4\t1\t4\tC\tA\t0.1\t1\t0\t1\t0.05\t300\t-4\n' | gzip -c > "$TMP/scz.gz"
printf 'SNP\tCHR\tBP\tA1\tA2\tZ\tN\tP\nrs1\t1\t1\tG\tA\t5\t400\t0.3\nrs2\t1\t2\tT\tC\t6\t400\t0.2\nrs3\t1\t3\tA\tT\t7\t400\t0.1\nrs4\t1\t4\tA\tC\t8\t400\t0.05\n' | gzip -c > "$TMP/cog.gz"

bash "$PROJECT_DIR/scripts/pleio/align_eur_pleio_inputs.sh" \
  "$TMP/mdd.gz" "$TMP/edu.gz" "$TMP/scz.gz" "$TMP/cog.gz" "$TMP/out"

for file in "$TMP/out"/*.aligned.tsv.gz; do
  [[ "$(gzip -cd "$file" | wc -l | tr -d ' ')" == "4" ]]
  [[ "$(gzip -cd "$file" | head -n 1)" == $'SNP\tA1\tA2\tZ\tN' ]]
done

gzip -cd "$TMP/out/MDD_EUR.aligned.tsv.gz" | grep -F $'rs2\tT\tC\t-2\t100' >/dev/null
gzip -cd "$TMP/out/SCZ_EUR.aligned.tsv.gz" | grep -F $'rs4\tA\tC\t4\t300' >/dev/null
gzip -cd "$TMP/out/COGENT_EUR.aligned.tsv.gz" | grep -F $'rs1\tA\tG\t-5\t400' >/dev/null
if gzip -cd "$TMP/out/MDD_EUR.aligned.tsv.gz" | grep -F 'rs3' >/dev/null; then
  echo "ERROR: strand-ambiguous rs3 was retained." >&2
  exit 1
fi

echo "EUR PLEIO alignment test passed."
