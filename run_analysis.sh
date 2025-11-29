#!/bin/bash
#
# run_analysis.sh
# Population Genetics Analysis Script
# 
# This script performs population genetics analysis using 1000 Genomes Project data (Chromosome 22).
# It includes VCF preprocessing, PCA analysis, ADMIXTURE analysis, and Fst calculations.
#
# Input Files:
#   - VCF: ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#   - Metadata: integrated_call_samples_v3.20130502.ALL.panel
#
# Tools Required:
#   - PLINK 1.9
#   - ADMIXTURE
#

set -e  # Exit on error

# ==============================================================================
# STEP 0: Check if required tools are installed
# ==============================================================================
echo "=== Checking for required tools ==="

# Check if PLINK is installed
if ! command -v plink &> /dev/null; then
    echo "ERROR: plink is not installed or not in PATH."
    echo "Please install PLINK 1.9 and ensure it is accessible."
    exit 1
fi
echo "plink found: $(command -v plink)"

# Check if ADMIXTURE is installed
if ! command -v admixture &> /dev/null; then
    echo "ERROR: admixture is not installed or not in PATH."
    echo "Please install ADMIXTURE and ensure it is accessible."
    exit 1
fi
echo "admixture found: $(command -v admixture)"

echo "All required tools are available."
echo ""

# ==============================================================================
# Define input files and output directory
# ==============================================================================
VCF_FILE="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_FILE="integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR="results"

# Check if input files exist
if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: VCF file not found: $VCF_FILE"
    exit 1
fi

if [ ! -f "$PANEL_FILE" ]; then
    echo "ERROR: Panel file not found: $PANEL_FILE"
    exit 1
fi

# ==============================================================================
# Create output directory
# ==============================================================================
echo "=== Creating output directory ==="
mkdir -p "$OUTPUT_DIR"
echo "Output directory created: $OUTPUT_DIR"
echo ""

# ==============================================================================
# STEP 1: Convert VCF to PLINK binary format (.bed, .bim, .fam)
# ==============================================================================
echo "=== Step 1: Converting VCF to PLINK binary format ==="
# This step converts the VCF file to PLINK binary format for downstream analysis.
# Output files: chr22.bed, chr22.bim, chr22.fam

plink --vcf "$VCF_FILE" \
      --make-bed \
      --out "$OUTPUT_DIR/chr22" \
      --allow-extra-chr

echo "VCF converted to PLINK binary format."
echo ""

# ==============================================================================
# STEP 2: Filter variants (MAF and genotype missingness)
# ==============================================================================
echo "=== Step 2: Filtering variants ==="
# Apply quality control filters:
#   --maf 0.05: Remove variants with minor allele frequency < 5%
#   --geno 0.1: Remove variants with >10% missing genotypes
# Output files: chr22_filtered.bed, chr22_filtered.bim, chr22_filtered.fam

plink --bfile "$OUTPUT_DIR/chr22" \
      --maf 0.05 \
      --geno 0.1 \
      --make-bed \
      --out "$OUTPUT_DIR/chr22_filtered" \
      --allow-extra-chr

echo "Variant filtering complete."
echo ""

# ==============================================================================
# STEP 3: LD pruning
# ==============================================================================
echo "=== Step 3: Performing LD pruning ==="
# LD pruning removes SNPs that are in high linkage disequilibrium.
# --indep-pairwise 50 10 0.1:
#   - Window size: 50 SNPs
#   - Step size: 10 SNPs
#   - r^2 threshold: 0.1
# This generates two files:
#   - chr22_filtered.prune.in: SNPs to keep
#   - chr22_filtered.prune.out: SNPs to remove

plink --bfile "$OUTPUT_DIR/chr22_filtered" \
      --indep-pairwise 50 10 0.1 \
      --out "$OUTPUT_DIR/chr22_filtered" \
      --allow-extra-chr

echo "LD pruning SNP list generated."

# Extract the pruned SNPs to create a new dataset
# Output files: chr22_pruned.bed, chr22_pruned.bim, chr22_pruned.fam
plink --bfile "$OUTPUT_DIR/chr22_filtered" \
      --extract "$OUTPUT_DIR/chr22_filtered.prune.in" \
      --make-bed \
      --out "$OUTPUT_DIR/chr22_pruned" \
      --allow-extra-chr

echo "LD-pruned dataset created."
echo ""

# ==============================================================================
# STEP 4: PCA Analysis
# ==============================================================================
echo "=== Step 4: Running PCA analysis ==="
# Run PCA on the pruned dataset to extract the top 20 principal components.
# This helps identify population structure in the data.
# Output files: chr22_pca.eigenval, chr22_pca.eigenvec

plink --bfile "$OUTPUT_DIR/chr22_pruned" \
      --pca 20 \
      --out "$OUTPUT_DIR/chr22_pca" \
      --allow-extra-chr

echo "PCA analysis complete. Results saved to:"
echo "  - $OUTPUT_DIR/chr22_pca.eigenval (eigenvalues)"
echo "  - $OUTPUT_DIR/chr22_pca.eigenvec (eigenvectors)"
echo ""

# ==============================================================================
# STEP 5: ADMIXTURE Analysis
# ==============================================================================
echo "=== Step 5: Running ADMIXTURE analysis ==="
# Run ADMIXTURE on the pruned dataset for different values of K (number of ancestral populations).
# ADMIXTURE requires the .bed file to be in the same directory for output.
# We run K values from 3 to 5.

# Change to output directory for ADMIXTURE (it outputs to current directory)
cd "$OUTPUT_DIR"

for K in 3 4 5; do
    echo "Running ADMIXTURE with K=$K..."
    admixture chr22_pruned.bed "$K"
    echo "ADMIXTURE K=$K complete."
done

# Return to original directory
cd ..

echo "ADMIXTURE analysis complete for K=3, 4, 5."
echo "Results saved to:"
echo "  - chr22_pruned.3.Q, chr22_pruned.3.P (K=3)"
echo "  - chr22_pruned.4.Q, chr22_pruned.4.P (K=4)"
echo "  - chr22_pruned.5.Q, chr22_pruned.5.P (K=5)"
echo ""

# ==============================================================================
# STEP 6: Convert metadata panel to PLINK cluster format
# ==============================================================================
echo "=== Step 6: Converting metadata panel to PLINK cluster format ==="
# Convert the panel file to a format PLINK accepts for cluster assignment.
# PLINK expects a cluster file with three columns: Family ID, Individual ID, Cluster
# The panel file has columns: sample, pop, super_pop, gender
# We use the 'pop' column as the cluster assignment.

# Create cluster file from panel file (skip header line)
# Format: FID IID Cluster
# Since .fam files from VCF typically use 0 as FID and sample ID as IID,
# we need to match that format.

awk 'NR>1 {print "0", $1, $2}' "$PANEL_FILE" > "$OUTPUT_DIR/clusters.txt"

echo "Cluster file created: $OUTPUT_DIR/clusters.txt"
echo ""

# ==============================================================================
# STEP 7: Fst Calculation
# ==============================================================================
echo "=== Step 7: Calculating pairwise Fst between populations ==="
# Calculate pairwise Fst (fixation index) between populations using PLINK.
# Fst measures genetic differentiation between populations.
# The --fst flag calculates Fst using the cluster assignments.

plink --bfile "$OUTPUT_DIR/chr22_pruned" \
      --fst \
      --within "$OUTPUT_DIR/clusters.txt" \
      --out "$OUTPUT_DIR/chr22_fst" \
      --allow-extra-chr

echo "Fst calculation complete."
echo "Results saved to: $OUTPUT_DIR/chr22_fst.fst"
echo ""

# ==============================================================================
# Analysis Complete
# ==============================================================================
echo "=========================================="
echo "Population genetics analysis complete!"
echo "=========================================="
echo ""
echo "Output files in $OUTPUT_DIR/:"
echo "  - PLINK binary files: chr22.bed/bim/fam, chr22_filtered.bed/bim/fam, chr22_pruned.bed/bim/fam"
echo "  - LD pruning: chr22_filtered.prune.in, chr22_filtered.prune.out"
echo "  - PCA: chr22_pca.eigenval, chr22_pca.eigenvec"
echo "  - ADMIXTURE: chr22_pruned.{3,4,5}.Q, chr22_pruned.{3,4,5}.P"
echo "  - Fst: chr22_fst.fst"
echo "  - Cluster file: clusters.txt"
echo ""
echo "Script completed successfully."
