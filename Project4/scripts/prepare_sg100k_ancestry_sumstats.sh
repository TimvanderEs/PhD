#!/usr/bin/env bash
set -euo pipefail

# Prepare Chinese-, Indian-, and Malay-specific SG100K summary statistics for LDSC.
# The provider files are controlled inputs and are not distributed with this repository.

SOURCE_ROOT="${SOURCE_ROOT:-/home/ec2-user/sg100K_sumstats}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/home/ec2-user/COV_LDSC_paper/VALIDATION_20260921/SG100K_ancestry_matched}"
MASTER_MAP="${MASTER_MAP:-/home/ec2-user/COV_LDSC_paper/sg_master_bp38_unordered_allele.map}"
MERGE_ALLELES="${MERGE_ALLELES:-/home/ec2-user/COV_LDSC_paper/sg_consensus_867k.merge_alleles.txt}"
MUNGE_SCRIPT="${MUNGE_SCRIPT:-/home/ec2-user/software/cov-ldsc/munge_sumstats.py}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/miniconda3/envs/ldsc/bin/python}"

PREPARED_DIR="${OUTPUT_ROOT}/prepared"
MUNGED_DIR="${OUTPUT_ROOT}/munged"
QC_DIR="${OUTPUT_ROOT}/qc"
LOG_DIR="${OUTPUT_ROOT}/logs"

traits=(
  c1_sbp_avg
  c11_whr
  dbi13_height
  c8_bmi
  c37_dlab_ldl_mmol_l_log10
  c20_t2d
)
ancestries=(chinese indian malay)

mkdir -p "${PREPARED_DIR}" "${MUNGED_DIR}" "${QC_DIR}" "${LOG_DIR}"

for required in "${MASTER_MAP}" "${MERGE_ALLELES}" "${MUNGE_SCRIPT}"; do
  if [[ ! -s "${required}" ]]; then
    echo "Required input is missing or empty: ${required}" >&2
    exit 1
  fi
done

printf 'trait\tancestry\tprepared_rows\tmunged_rows\n' > "${QC_DIR}/input_counts.tsv"

for trait in "${traits[@]}"; do
  echo "[$(date --iso-8601=seconds)] Preparing ${trait}"

  files=()
  for chr in {1..22}; do
    file="${SOURCE_ROOT}/${trait}/${trait}_chr${chr}.tsv"
    if [[ ! -s "${file}" ]]; then
      echo "Missing chromosome input: ${file}" >&2
      exit 1
    fi
    files+=("${file}")
  done

  chinese_out="${PREPARED_DIR}/${trait}.chinese.ldsc_input.tsv"
  indian_out="${PREPARED_DIR}/${trait}.indian.ldsc_input.tsv"
  malay_out="${PREPARED_DIR}/${trait}.malay.ldsc_input.tsv"
  extraction_marker="${QC_DIR}/${trait}.extraction_complete"

  if [[ ! -s "${extraction_marker}" ]]; then
    rm -f "${chinese_out}" "${chinese_out}.gz" "${indian_out}" "${indian_out}.gz" "${malay_out}" "${malay_out}.gz"
    for out in "${chinese_out}" "${indian_out}" "${malay_out}"; do
      printf 'SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tZ\tP\tN\n' > "${out}"
    done

    LC_ALL=C awk -v OFS='\t' \
    -v MAP="${MASTER_MAP}" \
    -v OUT_CHINESE="${chinese_out}" \
    -v OUT_INDIAN="${indian_out}" \
    -v OUT_MALAY="${malay_out}" '
    function missing(value) {
      return value == "" || value == "NA" || value == "NaN" || value == "nan" || value == "."
    }
    BEGIN {
      while ((getline line < MAP) > 0) {
        split(line, map_fields, "\t")
        variant_map[map_fields[1]] = map_fields[2]
      }
      close(MAP)
      ancestry[1] = "chinese"
      ancestry[2] = "indian"
      ancestry[3] = "malay"
      output[1] = OUT_CHINESE
      output[2] = OUT_INDIAN
      output[3] = OUT_MALAY
    }
    $1 == "CHR" && $2 == "POS" {
      delete column
      for (i = 1; i <= NF; i++) column[$i] = i
      next
    }
    {
      chr = $(column["CHR"])
      bp = $(column["POS"])
      a1 = $(column["A1"])
      a2 = $(column["A2"])

      if (missing(chr) || missing(bp) || missing(a1) || missing(a2)) next
      if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/) next

      if (a1 < a2) pair = a1 ":" a2
      else pair = a2 ":" a1
      lookup = chr ":" bp ":" pair
      if (!(lookup in variant_map)) next

      snp = variant_map[lookup]
      split(snp, snp_fields, ":")
      ref_a1 = snp_fields[3]
      ref_a2 = snp_fields[4]

      for (group = 1; group <= 3; group++) {
        suffix = ancestry[group]
        frq_index = column["FREQ_" suffix]
        beta_index = column["BETA_" suffix]
        se_index = column["SE_" suffix]
        p_index = column["P_" suffix]
        n_index = column["N_" suffix]

        if (!frq_index || !beta_index || !se_index || !p_index || !n_index) {
          print "Required ancestry column missing for " suffix > "/dev/stderr"
          exit 2
        }

        frq = $frq_index
        beta = $beta_index
        se = $se_index
        p = $p_index
        n = $n_index
        if (missing(frq) || missing(beta) || missing(se) || missing(p) || missing(n)) continue
        if ((se + 0) <= 0 || (p + 0) <= 0 || (p + 0) > 1 || (n + 0) <= 0) continue

        beta += 0
        se += 0
        frq += 0
        if (a1 == ref_a1 && a2 == ref_a2) {
          beta_out = beta
          frq_out = frq
        } else if (a1 == ref_a2 && a2 == ref_a1) {
          beta_out = -beta
          frq_out = 1 - frq
        } else {
          continue
        }
        z = beta_out / se
        print snp, chr, bp, ref_a1, ref_a2, frq_out, beta_out, se, z, p, n >> output[group]
      }
    }
    ' "${files[@]}"
    printf 'completed\t%s\n' "$(date --iso-8601=seconds)" > "${extraction_marker}"
  else
    echo "[$(date --iso-8601=seconds)] Reusing completed extraction for ${trait}"
  fi

  for ancestry in "${ancestries[@]}"; do
    prepared="${PREPARED_DIR}/${trait}.${ancestry}.ldsc_input.tsv"
    if [[ -s "${prepared}" ]]; then
      gzip -f "${prepared}"
    fi
    prepared="${prepared}.gz"
    if [[ ! -s "${prepared}" ]]; then
      echo "Prepared file is missing after extraction: ${prepared}" >&2
      exit 1
    fi
    prepared_rows=$(gzip -cd "${prepared}" | awk 'END {print NR - 1}')

    mkdir -p "${MUNGED_DIR}/${ancestry}"
    prefix="${MUNGED_DIR}/${ancestry}/${trait}.${ancestry}"
    if [[ ! -s "${prefix}.sumstats.gz" ]]; then
      "${PYTHON_BIN}" "${MUNGE_SCRIPT}" \
        --sumstats "${prepared}" \
        --out "${prefix}" \
        --merge-alleles "${MERGE_ALLELES}" \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --N-col N \
        --p P \
        --signed-sumstats Z,0 \
        --chunksize 500000 \
        > "${LOG_DIR}/${trait}.${ancestry}.munge.stdout.log" 2>&1
    fi

    if [[ ! -s "${prefix}.sumstats.gz" ]]; then
      echo "Munging failed for ${trait} ${ancestry}" >&2
      exit 1
    fi
    munged_rows=$(gzip -cd "${prefix}.sumstats.gz" | awk 'END {print NR - 1}')
    printf '%s\t%s\t%s\t%s\n' "${trait}" "${ancestry}" "${prepared_rows}" "${munged_rows}" >> "${QC_DIR}/input_counts.tsv"
  done
done

sha256sum "${PREPARED_DIR}"/*.gz "${MUNGED_DIR}"/*/*.sumstats.gz > "${QC_DIR}/input_checksums.sha256"
echo "[$(date --iso-8601=seconds)] Ancestry-specific preparation complete"
