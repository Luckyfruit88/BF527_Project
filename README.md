# Population Structure and Genetic Differentiation in the 1000 Genomes Project

A Chromosome 22 Analysis using PLINK and ADMIXTURE

## Description

Human populations exhibit substantial genetic diversity shaped by demographic history, migration, genetic drift, and natural selection. The 1000 Genomes Project (Phase 3) provides one of the most comprehensive open-source datasets for studying global population structure. Chromosome 22, although small, contains dense SNP variation and is commonly used in teaching-level population genetics analyses due to its manageable size and representative diversity patterns.

### Objectives

This project investigates global human population structure and genetic differentiation using SNP data from chromosome 22 of the 1000 Genomes Project. We focus on:

- **Principal Component Analysis (PCA)** of individuals across worldwide populations
- **Ancestry estimation** through ADMIXTURE
- **Genetic differentiation** between populations using pairwise FST
- **Visualization and interpretation** of major clustering patterns

## Data Sources

We use publicly available data from the 1000 Genomes Project (Phase 3):

| File | Description |
|------|-------------|
| `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | Chromosome 22 genotype VCF |
| `integrated_call_samples_v3.20130502.ALL.panel` | Sample population metadata |

These data contain genotype calls for **2,504 individuals** across **26 global populations**.

## Prerequisites

### Tools

- **PLINK 1.9** - For genotype data processing and PCA
- **ADMIXTURE** - For ancestry estimation

### Python Libraries

- `pandas` - Data manipulation
- `matplotlib` - Plotting
- `seaborn` - Statistical visualization

Install Python dependencies:

```bash
pip install pandas matplotlib seaborn
```

## Usage

### Step 1: Run the Analysis Pipeline

```bash
./run_analysis.sh
```

This script performs:
1. VCF to PLINK binary format conversion
2. Variant filtering (MAF ≥ 0.05, genotype missingness ≤ 10%)
3. LD pruning for ADMIXTURE analysis
4. PCA analysis (top 20 components)
5. ADMIXTURE analysis for K=3, 4, 5
6. Pairwise FST calculation between populations

### Step 2: Generate Visualizations

```bash
python plot_results.py
```

This script creates:
- PCA scatter plot (PC1 vs PC2)
- ADMIXTURE bar plot for K=5

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                        Input Data                               │
│  VCF File + Sample Metadata Panel                               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    1. Data Processing                           │
│  • Convert VCF → PLINK binary format (.bed/.bim/.fam)           │
│  • Filter variants (MAF ≥ 0.05, --geno 0.1)                     │
│  • LD-prune SNPs (--indep-pairwise 50 10 0.1)                   │
└─────────────────────────────────────────────────────────────────┘
                              │
          ┌───────────────────┼───────────────────┐
          ▼                   ▼                   ▼
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│    2. PCA       │ │  3. ADMIXTURE   │ │    4. FST       │
│                 │ │                 │ │                 │
│ PLINK --pca 20  │ │ K = 3, 4, 5     │ │ Pairwise FST    │
│                 │ │ ancestry est.   │ │ between pops    │
└─────────────────┘ └─────────────────┘ └─────────────────┘
          │                   │                   │
          ▼                   ▼                   ▼
┌─────────────────────────────────────────────────────────────────┐
│                    5. Visualization                             │
│  • PCA scatter plot (colored by super-population)               │
│  • ADMIXTURE bar plot (ancestry proportions)                    │
└─────────────────────────────────────────────────────────────────┘
```

## Expected Outputs

All output files are saved to the `results/` directory.

### Analysis Files

| File | Description |
|------|-------------|
| `chr22.bed/bim/fam` | PLINK binary files (raw) |
| `chr22_filtered.bed/bim/fam` | Filtered PLINK binary files |
| `chr22_pruned.bed/bim/fam` | LD-pruned PLINK binary files |
| `chr22_pca.eigenval` | PCA eigenvalues |
| `chr22_pca.eigenvec` | PCA eigenvectors (principal components) |
| `chr22_pruned.3.Q` | ADMIXTURE ancestry proportions (K=3) |
| `chr22_pruned.4.Q` | ADMIXTURE ancestry proportions (K=4) |
| `chr22_pruned.5.Q` | ADMIXTURE ancestry proportions (K=5) |
| `chr22_fst.fst` | Pairwise FST values |
| `clusters.txt` | Population cluster assignments |

### Visualization Files

| File | Description |
|------|-------------|
| `pca_plot.png` | PCA scatter plot showing PC1 vs PC2, colored by super-population (AFR, EUR, EAS, AMR, SAS) |
| `admixture_k5_plot.png` | Stacked bar chart showing ancestry proportions for K=5, sorted by population |

### Expected Results

- **PCA plots**: Clear clustering of continental groups (AFR, EUR, EAS, AMR, SAS)
- **ADMIXTURE bar plots**: Patterns of shared ancestry and population mixing
- **FST values**: Levels of genetic differentiation between major populations

## Project Structure

```
BF527_Project/
├── README.md                 # This file
├── run_analysis.sh           # Main analysis pipeline script
├── plot_results.py           # Visualization script
├── .gitignore               # Git ignore file
└── results/                  # Output directory (created by scripts)
    ├── chr22*.bed/bim/fam   # PLINK binary files
    ├── chr22_pca.*          # PCA results
    ├── chr22_pruned.*.Q     # ADMIXTURE results
    ├── chr22_fst.fst        # FST results
    ├── pca_plot.png         # PCA visualization
    └── admixture_k5_plot.png # ADMIXTURE visualization
```

## References

- 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.
- Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655-1664.
- Purcell, S., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *The American Journal of Human Genetics*, 81(3), 559-575.
