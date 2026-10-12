#!/bin/bash
# ------------------------------------------------------------------------------
# Script: master_anchor_download_QC_2026.sh
# ------------------------------------------------------------------------------

BASE_DIR="/home/ec2-user/COV_LDSC_paper/GWAS"
LOG_FILE="$BASE_DIR/download_qc_report.log"
echo "--- Download & QC Report $(date) ---" > "$LOG_FILE"

# Function to perform QC on a downloaded file
check_file() {
    local folder=$1
    local filename=$2
    local path="$BASE_DIR/$folder/$filename"

    echo -n "[QC] Checking $folder/$filename... " | tee -a "$LOG_FILE"

    # 1. Existence check
    if [ ! -f "$path" ]; then
        echo "FAIL (Missing)" | tee -a "$LOG_FILE"
        return 1
    fi

    # 2. Size check (HTML error pages are usually < 50KB, GWAS files are > 10MB)
    local size=$(stat -c%s "$path")
    if [ "$size" -lt 100000 ]; then
        echo "FAIL (File too small - likely HTML error)" | tee -a "$LOG_FILE"
        rm "$path" # Remove the garbage file
        return 1
    fi

    # 3. Gzip Integrity check
    if ! gzip -t "$path" 2>/dev/null; then
        echo "FAIL (Corrupt Gzip)" | tee -a "$LOG_FILE"
        rm "$path"
        return 1
    fi

    # 4. Header peek (Look for 'hm_' prefix or standard GWAS cols)
    local header=$(zcat "$path" | head -n 1)
    if [[ "$header" == *"html"* || "$header" == *"HTML"* ]]; then
        echo "FAIL (Contains HTML tags)" | tee -a "$LOG_FILE"
        rm "$path"
        return 1
    fi

    echo "PASS" | tee -a "$LOG_FILE"
    return 0
}

# --- 1. UKB (EUR) ---
mkdir -p "$BASE_DIR/UKB" && cd "$BASE_DIR/UKB"
declare -A ukb_map=( ["height"]="GCST90000612" ["bmi"]="GCST90000615" ["whr"]="GCST90000620" ["ldl"]="GCST90000617" ["tg"]="GCST90000618" ["sbp"]="GCST90000613" ["alt"]="GCST90000624" ["t2d"]="GCST90000614" ["smoking"]="GCST90000621" )
for t in "${!ukb_map[@]}"; do
    fname="ukb_${t}.h.tsv.gz"
    if [ ! -f "$fname" ]; then
        wget -q -c -O "$fname" "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/${ukb_map[$t]}/harmonised/${ukb_map[$t]}.h.tsv.gz"
    fi
    check_file "UKB" "$fname"
done

# --- 2. PAN-UKB (EAS) ---
mkdir -p "$BASE_DIR/PAN_UKB" && cd "$BASE_DIR/PAN_UKB"
declare -A pan_files=( ["height"]="continuous-50-both_sexes-irnt.tsv.bgz" ["bmi"]="continuous-21001-both_sexes-irnt.tsv.bgz" ["whr"]="continuous-48-both_sexes-irnt.tsv.bgz" ["ldl"]="biomarkers-30780-both_sexes-irnt.tsv.bgz" ["tg"]="biomarkers-30870-both_sexes-irnt.tsv.bgz" ["sbp"]="continuous-4080-both_sexes-irnt.tsv.bgz" ["alt"]="biomarkers-30620-both_sexes-irnt.tsv.bgz" ["t2d"]="categorical-2443-both_sexes-2443.tsv.bgz" ["smoking"]="categorical-20116-both_sexes-2.tsv.bgz" )
for t in "${!pan_files[@]}"; do
    fname="${pan_files[$t]}"
    if [ ! -f "$fname" ]; then
        aws s3 cp --no-sign-request "s3://pan-ukb-us-east-1/sumstats_flat_files/$fname" . > /dev/null
    fi
    check_file "PAN_UKB" "$fname"
done

# --- 3. JBB / 4. CBK / 5. KOGES_TWB (All follow the same GWAS Catalog range) ---
# Logic consolidated for brevity:
for bank in JBB CBK KOGES_TWB; do
    mkdir -p "$BASE_DIR/$bank" && cd "$BASE_DIR/$bank"
    if [ "$bank" == "JBB" ]; then declare -A map=( ["height"]="GCST90018907" ["bmi"]="GCST90018908" ["whr"]="GCST90018909" ["ldl"]="GCST90018917" ["tg"]="GCST90018919" ["sbp"]="GCST90018914" ["alt"]="GCST90018921" ["t2d"]="GCST90018925" ["smoking"]="GCST90018928" ); fi
    if [ "$bank" == "CBK" ]; then declare -A map=( ["height"]="GCST90018961" ["bmi"]="GCST90018963" ["whr"]="GCST90018967" ["ldl"]="GCST90018973" ["tg"]="GCST90018975" ["sbp"]="GCST90018970" ["alt"]="GCST90018977" ["t2d"]="GCST90018981" ["smoking"]="GCST90018984" ); fi
    if [ "$bank" == "KOGES_TWB" ]; then declare -A map=( ["height"]="GCST90018949" ["bmi"]="GCST90018933" ["whr"]="GCST90018956" ["ldl"]="GCST90018938" ["tg"]="GCST90018941" ["sbp"]="GCST90018944" ["t2d"]="GCST90018946" ["smoking"]="GCST90018950" ["alt"]="GCST90018953" ); fi

    for t in "${!map[@]}"; do
        id=${map[$t]}
        fname="${bank,,}_${t}.h.tsv.gz"
        if [ ! -f "$fname" ]; then
            wget -q -c -O "$fname" "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/$id/harmonised/$id.h.tsv.gz"
        fi
        check_file "$bank" "$fname"
    done
done

echo "QC complete. See $LOG_FILE for details."
