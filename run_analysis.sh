#!/bin/bash
#
# run_analysis.sh — Population Genetics Analysis Pipeline (Chr22)
#
# This refactored pipeline performs end-to-end preprocessing, QC, population structure
# inference, and pairwise FST calculations for the 1000 Genomes Project Chromosome 22 data.
# All downstream Python utilities now parse the logs and outputs produced here.

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
# Tool discovery (prefer repo-local binaries when available)
# ==============================================================================
echo "=== Checking required tools ==="

if [ -x "$SCRIPT_DIR/bin/plink" ]; then
    PLINK_BIN="$SCRIPT_DIR/bin/plink"
else
    PLINK_BIN="${PLINK_BIN:-$(command -v plink || true)}"
fi

if [ -z "$PLINK_BIN" ] || [ ! -x "$PLINK_BIN" ]; then
    echo "ERROR: PLINK binary not found. Install PLINK 1.9+ or set PLINK_BIN."
    exit 1
fi
echo "Using PLINK at: $PLINK_BIN"

ADMIXTURE_BIN="${ADMIXTURE_BIN:-$(command -v admixture || true)}"
if [ -z "$ADMIXTURE_BIN" ] || [ ! -x "$ADMIXTURE_BIN" ]; then
    echo "ERROR: ADMIXTURE binary not found. Install ADMIXTURE or set ADMIXTURE_BIN."
    exit 1
fi
echo "Using ADMIXTURE at: $ADMIXTURE_BIN"
echo ""

# ==============================================================================
# Inputs, outputs, and constants
# ==============================================================================
VCF_FILE="$SCRIPT_DIR/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_FILE="$SCRIPT_DIR/integrated_call_samples_v3.20130502.ALL.panel"
RESULTS_DIR="$SCRIPT_DIR/results"
LOG_DIR="$RESULTS_DIR/logs"

# ADMIXTURE parameters
ADMIXTURE_K_RANGE=(2 3 4 5 6)
ADMIXTURE_SEED=123
ADMIXTURE_CV_FOLDS=10

# Population definitions
SUPER_POPS=(AFR EUR EAS AMR SAS)

if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: Missing VCF file at $VCF_FILE"
    exit 1
fi

if [ ! -f "$PANEL_FILE" ]; then
    echo "ERROR: Missing panel metadata at $PANEL_FILE"
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

# Output file prefixes (all absolute paths)
CHR22_PREFIX="$RESULTS_DIR/chr22"
FILTER_PREFIX="$RESULTS_DIR/chr22_filtered"
PRUNED_PREFIX="$RESULTS_DIR/chr22_pruned"
PCA_PREFIX="$RESULTS_DIR/chr22_pca"
FST_SUPER_PREFIX="$RESULTS_DIR/chr22_fst_superpop"
FST_SUBPOP_PREFIX="$RESULTS_DIR/chr22_fst_subpop"
CLEAN_PANEL="$RESULTS_DIR/panel_cleaned.tsv"
SUPER_CLUSTERS="$RESULTS_DIR/clusters_superpop.txt"
SUBPOP_CLUSTERS="$RESULTS_DIR/clusters_subpop.txt"
CV_SUMMARY="$RESULTS_DIR/admixture_cv_errors.tsv"
PAIRWISE_SUPER_DIR="$RESULTS_DIR/fst_pairs_superpop"
PAIRWISE_SUBPOP_DIR="$RESULTS_DIR/fst_pairs_subpop"
PAIRWISE_SUPER_TSV="$RESULTS_DIR/chr22_fst_pairwise_superpop.tsv"
PAIRWISE_SUBPOP_TSV="$RESULTS_DIR/chr22_fst_pairwise_subpop.tsv"

echo "Outputs will be written to: $RESULTS_DIR"
echo "Logs will be written to: $LOG_DIR"
echo ""

# ==============================================================================
# STEP 1: Convert VCF to PLINK binary format
# ==============================================================================
echo "=== Step 1: VCF → PLINK binary ==="
run_and_log "$LOG_DIR/step1_convert.log" \
    "$PLINK_BIN" --vcf "$VCF_FILE" \
        --make-bed \
        --out "$CHR22_PREFIX" \
        --allow-extra-chr

# ==============================================================================
# STEP 2: Variant + sample QC (MAF, missingness)
# ==============================================================================
echo "=== Step 2: Variant & sample QC (MAF=0.01, geno=0.1, mind=0.1) ==="
run_and_log "$LOG_DIR/step2_filter.log" \
    "$PLINK_BIN" --bfile "$CHR22_PREFIX" \
        --maf 0.01 \
        --geno 0.1 \
        --mind 0.1 \
        --make-bed \
        --out "$FILTER_PREFIX" \
        --allow-extra-chr

# ==============================================================================
# STEP 3: LD pruning with expanded window
# ==============================================================================
echo "=== Step 3: LD pruning (window=200, step=10, r²=0.1) ==="
run_and_log "$LOG_DIR/step3_prune.log" \
    "$PLINK_BIN" --bfile "$FILTER_PREFIX" \
        --indep-pairwise 200 10 0.1 \
        --out "$FILTER_PREFIX" \
        --allow-extra-chr

run_and_log "$LOG_DIR/step3_extract.log" \
    "$PLINK_BIN" --bfile "$FILTER_PREFIX" \
        --extract "$FILTER_PREFIX.prune.in" \
        --make-bed \
        --out "$PRUNED_PREFIX" \
        --allow-extra-chr

# ==============================================================================
# STEP 4: Metadata cleaning + alignment
# ==============================================================================
echo "=== Step 4: Cleaning panel metadata ==="
if [ ! -f "$CHR22_PREFIX.fam" ]; then
    echo "ERROR: Expected $CHR22_PREFIX.fam but it was not created."
    exit 1
fi

# 创建干净的 panel 文件（仅包含 VCF 中存在的样本）
awk 'BEGIN{OFS="\t"} 
    NR==FNR {samples[$2]; next} 
    FNR==1 {print "sample","pop","super_pop","gender"; next} 
    ($1 in samples) {print $1,$2,$3,$4}' \
    "$CHR22_PREFIX.fam" "$PANEL_FILE" > "$CLEAN_PANEL"

if [ ! -s "$CLEAN_PANEL" ]; then
    echo "ERROR: Cleaned panel file is empty."
    exit 1
fi

# 提取所有亚群体列表（26个）
mapfile -t SUBPOPS < <(awk 'NR>1 {print $2}' "$CLEAN_PANEL" | sort -u)
echo "Found ${#SUBPOPS[@]} sub-populations: ${SUBPOPS[*]}"

echo "Clean panel saved to $CLEAN_PANEL"

# ==============================================================================
# STEP 5: Build cluster files for super-populations and sub-populations
# ==============================================================================
echo "=== Step 5: Building cluster files ==="

# 超群体聚类文件 (5 groups)
awk 'BEGIN{OFS=" "} 
    NR==FNR {sp[$1]=$3; next} 
    {iid=$2; if (iid in sp) {print $1, iid, sp[iid]}}' \
    "$CLEAN_PANEL" "$PRUNED_PREFIX.fam" > "$SUPER_CLUSTERS"

# 亚群体聚类文件 (26 groups)
awk 'BEGIN{OFS=" "} 
    NR==FNR {pop[$1]=$2; next} 
    {iid=$2; if (iid in pop) {print $1, iid, pop[iid]}}' \
    "$CLEAN_PANEL" "$PRUNED_PREFIX.fam" > "$SUBPOP_CLUSTERS"

if [ ! -s "$SUPER_CLUSTERS" ]; then
    echo "ERROR: Super-population cluster file is empty."
    exit 1
fi

if [ ! -s "$SUBPOP_CLUSTERS" ]; then
    echo "ERROR: Sub-population cluster file is empty."
    exit 1
fi

echo "Super-population clusters: $SUPER_CLUSTERS"
echo "Sub-population clusters: $SUBPOP_CLUSTERS"

# ==============================================================================
# STEP 6: Principal component analysis
# ==============================================================================
echo "=== Step 6: PCA on pruned genotype set ==="
run_and_log "$LOG_DIR/step6_pca.log" \
    "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
        --pca 20 \
        --out "$PCA_PREFIX" \
        --allow-extra-chr

# ==============================================================================
# STEP 7: ADMIXTURE (K = 2–6) with CV errors + seed
# ==============================================================================
echo "=== Step 7: ADMIXTURE K=2-6 with CV tracking ==="

# 初始化 CV summary 文件
echo -e "K\tCV_Error" > "$CV_SUMMARY"

# 需要在 results 目录中运行 ADMIXTURE（因为它会在当前目录生成输出）
pushd "$RESULTS_DIR" > /dev/null

for K in "${ADMIXTURE_K_RANGE[@]}"; do
    # 使用绝对路径存储日志
    LOG_FILE="$LOG_DIR/admixture_K${K}.log"
    echo "Running ADMIXTURE for K=$K (seed=$ADMIXTURE_SEED, cv=$ADMIXTURE_CV_FOLDS)"
    
    # 运行 ADMIXTURE
    "$ADMIXTURE_BIN" --cv="$ADMIXTURE_CV_FOLDS" --seed="$ADMIXTURE_SEED" \
        "$(basename "$PRUNED_PREFIX").bed" "$K" 2>&1 | tee "$LOG_FILE"

    # 提取 CV error（修复：使用更精确的正则表达式）
    # ADMIXTURE 输出格式: "CV error (K=X): Y.YYYYY"
    CV_LINE=$(grep -i "CV error" "$LOG_FILE" | tail -n 1 || true)
    if [ -n "$CV_LINE" ]; then
        # 提取冒号后的数值
        CV_VALUE=$(echo "$CV_LINE" | grep -oP '(?<=:\s)\d+\.\d+' || echo "")
        # 如果上述失败，尝试备用方法
        if [ -z "$CV_VALUE" ]; then
            CV_VALUE=$(echo "$CV_LINE" | sed -n 's/.*: *\([0-9]*\.[0-9]*\).*/\1/p')
        fi
    else
        CV_VALUE="NA"
    fi
    
    # 验证提取的值
    if [ -z "$CV_VALUE" ]; then
        CV_VALUE="NA"
    fi
    
    echo -e "${K}\t${CV_VALUE}" >> "$CV_SUMMARY"
    echo "K=$K CV error: $CV_VALUE"
done

popd > /dev/null

echo "ADMIXTURE runs complete. CV summary stored at $CV_SUMMARY"

# ==============================================================================
# STEP 8: Global FST calculations
# ==============================================================================
echo "=== Step 8: Global FST (super-populations and sub-populations) ==="

# 8a: Super-population FST
run_and_log "$LOG_DIR/step8_fst_superpop.log" \
    "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
        --fst \
        --within "$SUPER_CLUSTERS" \
        --out "$FST_SUPER_PREFIX" \
        --allow-extra-chr

# 8b: Sub-population FST (26 groups)
run_and_log "$LOG_DIR/step8_fst_subpop.log" \
    "$PLINK_BIN" --bfile "$PRUNED_PREFIX" \
        --fst \
        --within "$SUBPOP_CLUSTERS" \
        --out "$FST_SUBPOP_PREFIX" \
        --allow-extra-chr

echo "Global FST results saved to:"
echo "  Super-pop: $FST_SUPER_PREFIX.fst"
echo "  Sub-pop: $FST_SUBPOP_PREFIX.fst"

# ==============================================================================
# Helper function: Parse FST file and calculate statistics
# ==============================================================================
calculate_fst_stats() {
    local fst_file="$1"
    
    if [ ! -f "$fst_file" ]; then
        echo "NA\tNA\t0"
        return
    fi
    
    # 验证 FST 文件格式（PLINK 1.9: CHR SNP POS NMISS FST）
    local header
    header=$(head -1 "$fst_file")
    
    # 确定 FST 和 NMISS 列的位置
    local nmiss_col=4
    local fst_col=5
    
    # 检查是否为预期格式
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
            # 跳过无效值
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
# STEP 9: Pairwise FST for super-populations (5 choose 2 = 10 pairs)
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

        # 使用辅助函数计算统计数据
        STATS=$(calculate_fst_stats "$PAIR_PREFIX.fst")
        echo -e "${SP1}\t${SP2}\t${STATS}" >> "$PAIRWISE_SUPER_TSV"
    done
done

echo "Super-population pairwise FST saved to: $PAIRWISE_SUPER_TSV"

# ==============================================================================
# STEP 10: Pairwise FST for sub-populations (26 choose 2 = 325 pairs)
# ==============================================================================
echo "=== Step 10: Pairwise FST for sub-populations (26 groups) ==="
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

        # 使用辅助函数计算统计数据
        STATS=$(calculate_fst_stats "$PAIR_PREFIX.fst")
        echo -e "${POP1}\t${POP2}\t${SUPER1}\t${SUPER2}\t${STATS}" >> "$PAIRWISE_SUBPOP_TSV"
    done
done

echo ""
echo "Sub-population pairwise FST saved to: $PAIRWISE_SUBPOP_TSV"

# ==============================================================================
# STEP 11: Generate summary statistics
# ==============================================================================
echo "=== Step 11: Generating summary statistics ==="

# 从日志中提取关键统计数据（使用更健壮的解析方法）
extract_stat() {
    local file="$1"
    local pattern="$2"
    local default="${3:-N/A}"
    grep -oP "$pattern" "$file" 2>/dev/null | head -1 || echo "$default"
}

INITIAL_VARIANTS=$(extract_stat "$LOG_DIR/step1_convert.log" '\d+(?= variants loaded)')
INITIAL_SAMPLES=$(extract_stat "$LOG_DIR/step1_convert.log" '\d+(?= people)')

# Step 2 统计（包括因缺失率移除的样本）
REMOVED_VARIANTS_MAF=$(extract_stat "$LOG_DIR/step2_filter.log" '\d+(?= variants removed due to minor allele)' "0")
REMOVED_VARIANTS_GENO=$(extract_stat "$LOG_DIR/step2_filter.log" '\d+(?= variants removed due to missing genotype)' "0")
REMOVED_SAMPLES_MIND=$(extract_stat "$LOG_DIR/step2_filter.log" '\d+(?= people removed due to missing genotype)' "0")
FILTERED_VARIANTS=$(extract_stat "$LOG_DIR/step2_filter.log" '\d+(?= variants and)' | tail -1)
FILTERED_SAMPLES=$(grep -oP '(?<=variants and )\d+' "$LOG_DIR/step2_filter.log" 2>/dev/null | tail -1 || echo "N/A")

# Step 3 统计
PRUNED_IN=$(wc -l < "$FILTER_PREFIX.prune.in" 2>/dev/null || echo "N/A")
PRUNED_OUT=$(wc -l < "$FILTER_PREFIX.prune.out" 2>/dev/null || echo "N/A")
PRUNED_VARIANTS=$(extract_stat "$LOG_DIR/step3_extract.log" '\d+(?= variants and)' | tail -1)

# 找到最佳 K（最低 CV error）
BEST_K=$(awk -F'\t' 'NR>1 && $2!="NA" && $2!="" {
    if (best=="" || $2+0 < best_cv+0) {best=$1; best_cv=$2}
} END {print best}' "$CV_SUMMARY")
BEST_CV=$(awk -F'\t' 'NR>1 && $2!="NA" && $2!="" {
    if (best_cv=="" || $2+0 < best_cv+0) {best_cv=$2}
} END {print best_cv}' "$CV_SUMMARY")

# 找到最大 FST（超群体）
MAX_FST_SUPER=$(awk -F'\t' 'NR>1 && $4!="NA" && $4!="" {
    if ($4+0 > max+0) {max=$4; p1=$1; p2=$2}
} END {
    if (max != "") printf "%s vs %s: %.6f", p1, p2, max
    else print "N/A"
}' "$PAIRWISE_SUPER_TSV")

# 找到最大 FST（亚群体）
MAX_FST_SUBPOP=$(awk -F'\t' 'NR>1 && $6!="NA" && $6!="" {
    if ($6+0 > max+0) {max=$6; p1=$1; p2=$2}
} END {
    if (max != "") printf "%s vs %s: %.6f", p1, p2, max
    else print "N/A"
}' "$PAIRWISE_SUBPOP_TSV")

# 生成详细统计文件
cat > "$RESULTS_DIR/pipeline_stats.txt" << EOF
================================================================================
                    Population Genetics Pipeline Statistics
================================================================================
Date: $(date)
Pipeline Version: 2.0

--------------------------------------------------------------------------------
INPUT DATA
--------------------------------------------------------------------------------
  VCF File: $(basename "$VCF_FILE")
  Panel File: $(basename "$PANEL_FILE")
  Initial variants: $INITIAL_VARIANTS
  Initial samples: $INITIAL_SAMPLES

--------------------------------------------------------------------------------
STEP 2: QUALITY CONTROL (MAF=0.01, geno=0.1, mind=0.1)
--------------------------------------------------------------------------------
  Variants removed (MAF < 0.01): $REMOVED_VARIANTS_MAF
  Variants removed (missing rate > 0.1): $REMOVED_VARIANTS_GENO
  Samples removed (missing rate > 0.1): $REMOVED_SAMPLES_MIND
  Remaining variants: $FILTERED_VARIANTS
  Remaining samples: $FILTERED_SAMPLES

--------------------------------------------------------------------------------
STEP 3: LD PRUNING (window=200, step=10, r²<0.1)
--------------------------------------------------------------------------------
  Variants kept (prune.in): $PRUNED_IN
  Variants removed (prune.out): $PRUNED_OUT
  Final pruned variants: $PRUNED_VARIANTS

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
  Super-populations analyzed: ${SUPER_POPS[*]}
  Sub-populations analyzed: ${#SUBPOPS[@]} groups
  
  Maximum pairwise FST (super-pop): $MAX_FST_SUPER
  Maximum pairwise FST (sub-pop): $MAX_FST_SUBPOP

================================================================================
OUTPUT FILES
================================================================================
QC & Preprocessing:
  • $FILTER_PREFIX.{bed,bim,fam}
  • $PRUNED_PREFIX.{bed,bim,fam}

Population Structure:
  • $PCA_PREFIX.{eigenvec,eigenval}
  • $RESULTS_DIR/chr22_pruned.{2-6}.{Q,P}
  • $CV_SUMMARY

FST Analysis:
  • $PAIRWISE_SUPER_TSV
  • $PAIRWISE_SUBPOP_TSV

Metadata:
  • $CLEAN_PANEL
  • $SUPER_CLUSTERS
  • $SUBPOP_CLUSTERS

Logs:
  • $LOG_DIR/

================================================================================
EOF

echo "Pipeline statistics saved to: $RESULTS_DIR/pipeline_stats.txt"

# ==============================================================================
# Pipeline Complete
# ==============================================================================
echo ""
echo "=========================================="
echo "Pipeline complete! Key outputs:"
echo "=========================================="
echo ""
echo "QC & Preprocessing:"
echo "  • Logs: $LOG_DIR/"
echo "  • Filtered data: $FILTER_PREFIX.{bed,bim,fam}"
echo "  • Pruned data: $PRUNED_PREFIX.{bed,bim,fam}"
echo ""
echo "Population Structure:"
echo "  • PCA: $PCA_PREFIX.{eigenvec,eigenval}"
echo "  • ADMIXTURE (K=2-6): $RESULTS_DIR/chr22_pruned.{K}.{Q,P}"
echo "  • CV errors: $CV_SUMMARY"
echo ""
echo "FST Analysis:"
echo "  • Super-pop global: $FST_SUPER_PREFIX.fst"
echo "  • Sub-pop global: $FST_SUBPOP_PREFIX.fst"
echo "  • Super-pop pairwise: $PAIRWISE_SUPER_TSV"
echo "  • Sub-pop pairwise: $PAIRWISE_SUBPOP_TSV"
echo ""
echo "Metadata:"
echo "  • Cleaned panel: $CLEAN_PANEL"
echo "  • Super-pop clusters: $SUPER_CLUSTERS"
echo "  • Sub-pop clusters: $SUBPOP_CLUSTERS"
echo "  • Statistics: $RESULTS_DIR/pipeline_stats.txt"
echo ""
echo "=========================================="
echo "Script completed successfully at $(timestamp)"
echo "=========================================="