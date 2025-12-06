#!/bin/bash
# filepath: /workspaces/BF527_Project/resume_from_step7.sh
#
# resume_from_step7.sh — Resume pipeline from ADMIXTURE analysis
#
# This script continues the pipeline from Step 7 onwards, skipping completed preprocessing steps.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

run_and_log() {
    local log_file="$1"
    shift
    mkdir -p "$(dirname "$log_file")"
    : > "$log_file"
    echo "[$(timestamp)] Running: $*" | tee -a "$log_file"
    "$@" 2>&1 | tee -a "$log_file"
}

# ==============================================================================
# Tool discovery
# ==============================================================================
echo "=== Checking required tools ==="

if [ -x "$SCRIPT_DIR/bin/plink" ]; then
    PLINK_BIN="$SCRIPT_DIR/bin/plink"
else
    PLINK_BIN="${PLINK_BIN:-$(command -v plink || true)}"
fi

if [ -z "$PLINK_BIN" ] || [ ! -x "$PLINK_BIN" ]; then
    echo "ERROR: PLINK binary not found."
    exit 1
fi
echo "Using PLINK at: $PLINK_BIN"

ADMIXTURE_BIN="${ADMIXTURE_BIN:-$(command -v admixture || true)}"
if [ -z "$ADMIXTURE_BIN" ] || [ ! -x "$ADMIXTURE_BIN" ]; then
    echo "ERROR: ADMIXTURE binary not found."
    exit 1
fi
echo "Using ADMIXTURE at: $ADMIXTURE_BIN"
echo ""

# ==============================================================================
# Paths and constants
# ==============================================================================
RESULTS_DIR="$SCRIPT_DIR/results"
LOG_DIR="$RESULTS_DIR/logs"
PRUNED_PREFIX="$RESULTS_DIR/chr22_pruned"
CLEAN_PANEL="$RESULTS_DIR/panel_cleaned.tsv"
SUPER_CLUSTERS="$RESULTS_DIR/clusters_superpop.txt"
SUBPOP_CLUSTERS="$RESULTS_DIR/clusters_subpop.txt"
CV_SUMMARY="$RESULTS_DIR/admixture_cv_errors.tsv"
FST_SUPER_PREFIX="$RESULTS_DIR/chr22_fst_superpop"
FST_SUBPOP_PREFIX="$RESULTS_DIR/chr22_fst_subpop"
PAIRWISE_SUPER_DIR="$RESULTS_DIR/fst_pairs_superpop"
PAIRWISE_SUBPOP_DIR="$RESULTS_DIR/fst_pairs_subpop"
PAIRWISE_SUPER_TSV="$RESULTS_DIR/chr22_fst_pairwise_superpop.tsv"
PAIRWISE_SUBPOP_TSV="$RESULTS_DIR/chr22_fst_pairwise_subpop.tsv"

ADMIXTURE_K_RANGE=(2 3 4 5 6)
ADMIXTURE_SEED=123
ADMIXTURE_CV_FOLDS=10
SUPER_POPS=(AFR EUR EAS AMR SAS)

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

echo "Resuming pipeline from Step 7..."
echo "Working directory: $RESULTS_DIR"
echo ""

# ==============================================================================
# Verify prerequisites from Steps 1-6
# ==============================================================================
echo "=== Verifying prerequisite files from Steps 1-6 ==="
MISSING_FILES=()

[ ! -f "$PRUNED_PREFIX.bed" ] && MISSING_FILES+=("$PRUNED_PREFIX.bed")
[ ! -f "$PRUNED_PREFIX.bim" ] && MISSING_FILES+=("$PRUNED_PREFIX.bim")
[ ! -f "$PRUNED_PREFIX.fam" ] && MISSING_FILES+=("$PRUNED_PREFIX.fam")
[ ! -f "$CLEAN_PANEL" ] && MISSING_FILES+=("$CLEAN_PANEL")
[ ! -f "$SUPER_CLUSTERS" ] && MISSING_FILES+=("$SUPER_CLUSTERS")
[ ! -f "$SUBPOP_CLUSTERS" ] && MISSING_FILES+=("$SUBPOP_CLUSTERS")

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo "ERROR: Missing prerequisite files from Steps 1-6:"
    printf '  • %s\n' "${MISSING_FILES[@]}"
    echo ""
    echo "Please run the full pipeline first: ./run_analysis.sh"
    exit 1
fi

echo "✓ All prerequisite files found"
echo ""

# 提取亚群体列表
mapfile -t SUBPOPS < <(awk 'NR>1 {print $2}' "$CLEAN_PANEL" | sort -u)
echo "Found ${#SUBPOPS[@]} sub-populations"

# ==============================================================================
# STEP 7: ADMIXTURE (K = 2–6) with CV errors + seed
# ==============================================================================
echo "=== Step 7: ADMIXTURE K=2-6 with CV tracking ==="

# 初始化 CV summary 文件
echo -e "K\tCV_Error" > "$CV_SUMMARY"

# 需要在 results 目录中运行 ADMIXTURE
pushd "$RESULTS_DIR" > /dev/null

for K in "${ADMIXTURE_K_RANGE[@]}"; do
    LOG_FILE="$LOG_DIR/admixture_K${K}.log"
    
    # 检查是否已完成
    if [ -f "chr22_pruned.${K}.Q" ] && [ -f "chr22_pruned.${K}.P" ]; then
        echo "K=$K already completed, extracting CV error..."
        CV_LINE=$(grep -i "CV error" "$LOG_FILE" 2>/dev/null | tail -n 1 || echo "")
        if [ -n "$CV_LINE" ]; then
            CV_VALUE=$(echo "$CV_LINE" | grep -oP '(?<=:\s)\d+\.\d+' || echo "NA")
            [ -z "$CV_VALUE" ] && CV_VALUE=$(echo "$CV_LINE" | sed -n 's/.*: *\([0-9]*\.[0-9]*\).*/\1/p')
        else
            CV_VALUE="NA"
        fi
        [ -z "$CV_VALUE" ] && CV_VALUE="NA"
        echo -e "${K}\t${CV_VALUE}" >> "$CV_SUMMARY"
        echo "K=$K CV error: $CV_VALUE (from existing log)"
        continue
    fi
    
    echo "Running ADMIXTURE for K=$K (seed=$ADMIXTURE_SEED, cv=$ADMIXTURE_CV_FOLDS)"
    
    "$ADMIXTURE_BIN" --cv="$ADMIXTURE_CV_FOLDS" --seed="$ADMIXTURE_SEED" \
        "chr22_pruned.bed" "$K" 2>&1 | tee "$LOG_FILE"

    # 提取 CV error
    CV_LINE=$(grep -i "CV error" "$LOG_FILE" | tail -n 1 || true)
    if [ -n "$CV_LINE" ]; then
        CV_VALUE=$(echo "$CV_LINE" | grep -oP '(?<=:\s)\d+\.\d+' || echo "")
        [ -z "$CV_VALUE" ] && CV_VALUE=$(echo "$CV_LINE" | sed -n 's/.*: *\([0-9]*\.[0-9]*\).*/\1/p')
    else
        CV_VALUE="NA"
    fi
    [ -z "$CV_VALUE" ] && CV_VALUE="NA"
    
    echo -e "${K}\t${CV_VALUE}" >> "$CV_SUMMARY"
    echo "K=$K CV error: $CV_VALUE"
done

popd > /dev/null

echo "ADMIXTURE runs complete. CV summary: $CV_SUMMARY"
echo ""

# ==============================================================================
# Helper function: Parse FST file and calculate statistics
# ==============================================================================
calculate_fst_stats() {
    local fst_file="$1"
    
    if [ ! -f "$fst_file" ]; then
        echo "NA\tNA\t0"
        return
    fi
    
    local header
    header=$(head -1 "$fst_file")
    
    local nmiss_col=4
    local fst_col=5
    
    if echo "$header" | grep -qi "NMISS"; then
        nmiss_col=$(echo "$header" | tr '\t' '\n' | grep -n -i "NMISS" | cut -d: -f1)
    fi
    if echo "$header" | grep -qi "FST"; then
        fst_col=$(echo "$header" | tr '\t' '\n' | grep -n -i "FST" | cut -d: -f1)
    fi
    
    awk -v nmiss_col="$nmiss_col" -v fst_col="$fst_col" '
        NR>1 {
            fst = $fst_col
            nmiss = $nmiss_col
            if (fst == "nan" || fst == "NA" || fst == "." || fst == "") next
            if (nmiss == "nan" || nmiss == "NA" || nmiss == "." || nmiss == "") nmiss = 1
            
            sum_fst += fst
            sum_weight += nmiss
            weighted_sum += fst * nmiss
            n++
        }
        END {
            if (n > 0) {
                mean = sum_fst / n
                weighted = (sum_weight > 0) ? weighted_sum / sum_weight : 0
                printf "%.6f\t%.6f\t%d", mean, weighted, n
            } else {
                print "NA\tNA\t0"
            }
        }
    ' "$fst_file"
}

# ==============================================================================
# STEP 8: Global FST calculations
# ==============================================================================
echo "=== Step 8: Global FST (super-populations and sub-populations) ==="

run_and_log "$LOG_DIR/step8_fst_superpop.log" \
    "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
        --fst \
        --within "$SUPER_CLUSTERS" \
        --out "$FST_SUPER_PREFIX" \
        --allow-extra-chr

run_and_log "$LOG_DIR/step8_fst_subpop.log" \
    "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
        --fst \
        --within "$SUBPOP_CLUSTERS" \
        --out "$FST_SUBPOP_PREFIX" \
        --allow-extra-chr

echo "Global FST results saved"
echo ""

# ==============================================================================
# STEP 9: Pairwise FST for super-populations
# ==============================================================================
echo "=== Step 9: Pairwise FST for super-populations ==="
mkdir -p "$PAIRWISE_SUPER_DIR"
echo -e "pop1\tpop2\tmean_fst\tweighted_fst\tn_snps" > "$PAIRWISE_SUPER_TSV"

for ((i = 0; i < ${#SUPER_POPS[@]}; i++)); do
    for ((j = i + 1; j < ${#SUPER_POPS[@]}; j++)); do
        SP1="${SUPER_POPS[$i]}"
        SP2="${SUPER_POPS[$j]}"
        PAIR_PREFIX="$PAIRWISE_SUPER_DIR/fst_${SP1}_${SP2}"
        LOG_FILE="$LOG_DIR/fst_super_${SP1}_${SP2}.log"
        
        echo "  Calculating FST: $SP1 vs $SP2"
        run_and_log "$LOG_FILE" \
            "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
                --fst \
                --within "$SUPER_CLUSTERS" \
                --keep-cluster-names "$SP1" "$SP2" \
                --out "$PAIR_PREFIX" \
                --allow-extra-chr

        STATS=$(calculate_fst_stats "$PAIR_PREFIX.fst")
        echo -e "${SP1}\t${SP2}\t${STATS}" >> "$PAIRWISE_SUPER_TSV"
    done
done

echo "Super-population pairwise FST: $PAIRWISE_SUPER_TSV"
echo ""

# ==============================================================================
# STEP 10: Pairwise FST for sub-populations
# ==============================================================================
echo "=== Step 10: Pairwise FST for sub-populations (${#SUBPOPS[@]} groups) ==="
mkdir -p "$PAIRWISE_SUBPOP_DIR"
echo -e "pop1\tpop2\tsuper_pop1\tsuper_pop2\tmean_fst\tweighted_fst\tn_snps" > "$PAIRWISE_SUBPOP_TSV"

# 创建 pop -> super_pop 映射
declare -A POP_TO_SUPER
while IFS=$'\t' read -r sample pop super_pop gender; do
    [ "$sample" = "sample" ] && continue
    POP_TO_SUPER["$pop"]="$super_pop"
done < "$CLEAN_PANEL"

TOTAL_PAIRS=$(( ${#SUBPOPS[@]} * (${#SUBPOPS[@]} - 1) / 2 ))
CURRENT_PAIR=0

for ((i = 0; i < ${#SUBPOPS[@]}; i++)); do
    for ((j = i + 1; j < ${#SUBPOPS[@]}; j++)); do
        POP1="${SUBPOPS[$i]}"
        POP2="${SUBPOPS[$j]}"
        SUPER1="${POP_TO_SUPER[$POP1]:-UNKNOWN}"
        SUPER2="${POP_TO_SUPER[$POP2]:-UNKNOWN}"
        PAIR_PREFIX="$PAIRWISE_SUBPOP_DIR/fst_${POP1}_${POP2}"
        LOG_FILE="$LOG_DIR/fst_subpop_${POP1}_${POP2}.log"
        
        CURRENT_PAIR=$((CURRENT_PAIR + 1))
        echo -ne "\r  Processing pair $CURRENT_PAIR / $TOTAL_PAIRS: $POP1 vs $POP2          "
        
        "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
            --fst \
            --within "$SUBPOP_CLUSTERS" \
            --keep-cluster-names "$POP1" "$POP2" \
            --out "$PAIR_PREFIX" \
            --allow-extra-chr > "$LOG_FILE" 2>&1

        STATS=$(calculate_fst_stats "$PAIR_PREFIX.fst")
        echo -e "${POP1}\t${POP2}\t${SUPER1}\t${SUPER2}\t${STATS}" >> "$PAIRWISE_SUBPOP_TSV"
    done
done

echo ""
echo "Sub-population pairwise FST: $PAIRWISE_SUBPOP_TSV"
echo ""

# ==============================================================================
# STEP 11: Generate summary statistics
# ==============================================================================
echo "=== Step 11: Generating summary statistics ==="

# 找到最佳 K
BEST_K=$(awk -F'\t' 'NR>1 && $2!="NA" && $2!="" {
    if (best=="" || $2+0 < best_cv+0) {best=$1; best_cv=$2}
} END {print best}' "$CV_SUMMARY")
BEST_CV=$(awk -F'\t' 'NR>1 && $2!="NA" && $2!="" {
    if (best_cv=="" || $2+0 < best_cv+0) {best_cv=$2}
} END {print best_cv}' "$CV_SUMMARY")

# 找到最大 FST
MAX_FST_SUPER=$(awk -F'\t' 'NR>1 && $4!="NA" && $4!="" {
    if ($4+0 > max+0) {max=$4; p1=$1; p2=$2}
} END {
    if (max != "") printf "%s vs %s: %.6f", p1, p2, max
    else print "N/A"
}' "$PAIRWISE_SUPER_TSV")

MAX_FST_SUBPOP=$(awk -F'\t' 'NR>1 && $6!="NA" && $6!="" {
    if ($6+0 > max+0) {max=$6; p1=$1; p2=$2}
} END {
    if (max != "") printf "%s vs %s: %.6f", p1, p2, max
    else print "N/A"
}' "$PAIRWISE_SUBPOP_TSV")

cat > "$RESULTS_DIR/pipeline_stats_resumed.txt" << EOF
================================================================================
           Population Genetics Pipeline - Resumed from Step 7
================================================================================
Date: $(date)
Resumed at: $(timestamp)

--------------------------------------------------------------------------------
STEP 7: ADMIXTURE ANALYSIS
--------------------------------------------------------------------------------
  K range tested: ${ADMIXTURE_K_RANGE[*]}
  Cross-validation folds: $ADMIXTURE_CV_FOLDS
  Random seed: $ADMIXTURE_SEED
  Best K: $BEST_K (CV error: $BEST_CV)
  
  CV Error Summary:
$(awk -F'\t' 'NR>1 {printf "    K=%s: %s\n", $1, $2}' "$CV_SUMMARY")

--------------------------------------------------------------------------------
STEP 8-10: FST ANALYSIS
--------------------------------------------------------------------------------
  Super-populations: ${SUPER_POPS[*]}
  Sub-populations: ${#SUBPOPS[@]} groups
  
  Maximum pairwise FST (super-pop): $MAX_FST_SUPER
  Maximum pairwise FST (sub-pop): $MAX_FST_SUBPOP

================================================================================
OUTPUT FILES
================================================================================
ADMIXTURE:
  • $RESULTS_DIR/chr22_pruned.{2-6}.{Q,P}
  • $CV_SUMMARY

FST Analysis:
  • $FST_SUPER_PREFIX.fst
  • $FST_SUBPOP_PREFIX.fst
  • $PAIRWISE_SUPER_TSV
  • $PAIRWISE_SUBPOP_TSV

Logs:
  • $LOG_DIR/

================================================================================
EOF

echo "Statistics saved to: $RESULTS_DIR/pipeline_stats_resumed.txt"

# ==============================================================================
# Complete
# ==============================================================================
echo ""
echo "=========================================="
echo "Pipeline resumed and completed!"
echo "=========================================="
echo ""
echo "Key outputs:"
echo "  • ADMIXTURE results: $RESULTS_DIR/chr22_pruned.{2-6}.{Q,P}"
echo "  • CV errors: $CV_SUMMARY"
echo "  • FST (super-pop): $PAIRWISE_SUPER_TSV"
echo "  • FST (sub-pop): $PAIRWISE_SUBPOP_TSV"
echo "  • Statistics: $RESULTS_DIR/pipeline_stats_resumed.txt"
echo ""
echo "=========================================="
echo "Completed at $(timestamp)"
echo "=========================================="