#!/usr/bin/env python3
"""
create_flowchart.py
Generate a flowchart for the Population Genetics Analysis Pipeline

This script creates a visual flowchart showing all steps in the analysis pipeline.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np


def create_flowchart():
    """Create a flowchart of the population genetics analysis pipeline."""
    
    fig, ax = plt.subplots(figsize=(16, 20))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 20)
    ax.axis('off')
    
    # Define colors
    colors = {
        'input': '#E8F5E9',      # Light green - input files
        'tool': '#E3F2FD',       # Light blue - tools
        'process': '#FFF3E0',    # Light orange - processing steps
        'output': '#FCE4EC',     # Light pink - output files
        'analysis': '#F3E5F5',   # Light purple - analysis steps
        'viz': '#FFFDE7',        # Light yellow - visualization
        'border': '#37474F'      # Dark gray - borders
    }
    
    def draw_box(x, y, width, height, text, color, fontsize=9, bold=False):
        """Draw a rounded rectangle box with text."""
        box = FancyBboxPatch(
            (x - width/2, y - height/2), width, height,
            boxstyle="round,pad=0.03,rounding_size=0.2",
            facecolor=color,
            edgecolor=colors['border'],
            linewidth=1.5
        )
        ax.add_patch(box)
        weight = 'bold' if bold else 'normal'
        ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
                fontweight=weight, wrap=True)
    
    def draw_arrow(start, end, color='#37474F'):
        """Draw an arrow between two points."""
        ax.annotate(
            '', xy=end, xytext=start,
            arrowprops=dict(
                arrowstyle='-|>',
                color=color,
                lw=1.5,
                connectionstyle='arc3,rad=0'
            )
        )
    
    # Title
    ax.text(8, 19.5, 'Population Genetics Analysis Pipeline', 
            ha='center', va='center', fontsize=18, fontweight='bold')
    ax.text(8, 19.0, '1000 Genomes Project - Chromosome 22', 
            ha='center', va='center', fontsize=12, style='italic', color='#666666')
    
    # =========================================================================
    # INPUT FILES (Top)
    # =========================================================================
    draw_box(4, 17.8, 5.5, 0.7, 'VCF File\n(Chr22 genotypes, ~196MB)', colors['input'], bold=True)
    draw_box(12, 17.8, 5.5, 0.7, 'Panel File\n(Sample metadata, pop/super_pop)', colors['input'], bold=True)
    
    # =========================================================================
    # STEP 0: Tool Check
    # =========================================================================
    draw_box(8, 16.5, 4, 0.6, 'Step 0: Check Tools\n(PLINK, ADMIXTURE)', colors['tool'])
    draw_arrow((4, 17.4), (7, 16.8))
    draw_arrow((12, 17.4), (9, 16.8))
    
    # =========================================================================
    # STEP 1: VCF to PLINK
    # =========================================================================
    draw_box(8, 15.3, 5, 0.8, 'Step 1: Convert VCF to PLINK Format\nplink --vcf → .bed/.bim/.fam', colors['process'])
    draw_arrow((8, 16.2), (8, 15.7))
    
    # Output
    draw_box(13, 15.3, 3, 0.6, 'chr22.bed\nchr22.bim\nchr22.fam', colors['output'], fontsize=8)
    draw_arrow((10.5, 15.3), (11.5, 15.3))
    
    # =========================================================================
    # STEP 2: Filter Variants
    # =========================================================================
    draw_box(8, 14.0, 5.5, 0.8, 'Step 2: Filter Variants\n--maf 0.05 --geno 0.1', colors['process'])
    draw_arrow((8, 14.9), (8, 14.4))
    
    # Output
    draw_box(13, 14.0, 3.2, 0.6, 'chr22_filtered.*\n(112,015 variants)', colors['output'], fontsize=8)
    draw_arrow((10.75, 14.0), (11.4, 14.0))
    
    # =========================================================================
    # STEP 3: LD Pruning
    # =========================================================================
    draw_box(8, 12.7, 5.5, 0.8, 'Step 3: LD Pruning\n--indep-pairwise 50 5 0.1', colors['process'])
    draw_arrow((8, 13.6), (8, 13.1))
    
    # Output
    draw_box(13, 12.7, 3.2, 0.6, 'chr22_pruned.*\n(LD-independent SNPs)', colors['output'], fontsize=8)
    draw_arrow((10.75, 12.7), (11.4, 12.7))
    
    # =========================================================================
    # BRANCHING POINT - Three parallel analyses
    # =========================================================================
    # Draw branching lines
    ax.plot([8, 8], [12.3, 11.6], color=colors['border'], lw=1.5)
    ax.plot([4, 12], [11.6, 11.6], color=colors['border'], lw=1.5)
    ax.plot([4, 4], [11.6, 11.2], color=colors['border'], lw=1.5)
    ax.plot([8, 8], [11.6, 11.2], color=colors['border'], lw=1.5)
    ax.plot([12, 12], [11.6, 11.2], color=colors['border'], lw=1.5)
    
    # Arrow heads
    draw_arrow((4, 11.3), (4, 10.9))
    draw_arrow((8, 11.3), (8, 10.9))
    draw_arrow((12, 11.3), (12, 10.9))
    
    # =========================================================================
    # STEP 4: PCA Analysis (Left Branch)
    # =========================================================================
    draw_box(4, 10.5, 4, 0.8, 'Step 4: PCA Analysis\n--pca 20', colors['analysis'])
    
    # Output
    draw_box(4, 9.3, 3.5, 0.7, 'chr22_pca.eigenvec\nchr22_pca.eigenval', colors['output'], fontsize=8)
    draw_arrow((4, 10.1), (4, 9.7))
    
    # =========================================================================
    # STEP 5: ADMIXTURE Analysis (Middle Branch)
    # =========================================================================
    draw_box(8, 10.5, 4, 0.8, 'Step 5: ADMIXTURE\nK=3, K=4, K=5', colors['analysis'])
    
    # Output
    draw_box(8, 9.3, 3.5, 0.7, 'chr22_pruned.3.Q\nchr22_pruned.4.Q\nchr22_pruned.5.Q', colors['output'], fontsize=8)
    draw_arrow((8, 10.1), (8, 9.7))
    
    # =========================================================================
    # STEP 6 & 7: Cluster & FST (Right Branch)
    # =========================================================================
    draw_box(12, 10.5, 4, 0.8, 'Step 6: Create Clusters\n(from Panel file)', colors['process'])
    
    draw_box(12, 9.3, 3.5, 0.6, 'clusters.txt\n(26 populations)', colors['output'], fontsize=8)
    draw_arrow((12, 10.1), (12, 9.6))
    
    draw_box(12, 8.2, 4, 0.8, 'Step 7: FST Calculation\n--fst --within clusters.txt', colors['analysis'])
    draw_arrow((12, 8.95), (12, 8.6))
    
    draw_box(12, 7.1, 3.5, 0.7, 'chr22_fst.fst\nMean Fst: 0.092', colors['output'], fontsize=8)
    draw_arrow((12, 7.8), (12, 7.5))
    
    # =========================================================================
    # VISUALIZATION SECTION
    # =========================================================================
    # Merge arrows to visualization
    ax.plot([4, 4], [8.95, 5.8], color=colors['border'], lw=1.5, linestyle='--')
    ax.plot([8, 8], [8.95, 5.8], color=colors['border'], lw=1.5, linestyle='--')
    ax.plot([4, 12], [5.8, 5.8], color=colors['border'], lw=1.5, linestyle='--')
    ax.plot([12, 12], [6.75, 5.8], color=colors['border'], lw=1.5, linestyle='--')
    
    draw_box(8, 5.3, 5, 0.7, 'Visualization (plot_results.py)\nPython: pandas, matplotlib, seaborn', colors['viz'], bold=True)
    draw_arrow((8, 5.8), (8, 5.65))
    
    # =========================================================================
    # OUTPUT PLOTS
    # =========================================================================
    draw_box(4, 4.0, 4.2, 1.0, 'Plot 1:\nPCA by Super Pop\n(5 continental groups)', colors['viz'], fontsize=8)
    draw_box(8, 4.0, 4.2, 1.0, 'Plot 2:\nPCA by Detailed Pop\n(26 populations)', colors['viz'], fontsize=8)
    draw_box(12, 4.0, 4.2, 1.0, 'Plot 3:\nADMIXTURE K=5\n(Sorted bar chart)', colors['viz'], fontsize=8)
    
    draw_arrow((6, 4.95), (4.5, 4.5))
    draw_arrow((8, 4.95), (8, 4.5))
    draw_arrow((10, 4.95), (11.5, 4.5))
    
    # Output file names
    draw_box(4, 2.8, 3.8, 0.5, 'pca_super_pop.png', colors['output'], fontsize=8)
    draw_box(8, 2.8, 3.8, 0.5, 'pca_detailed_pop.png', colors['output'], fontsize=8)
    draw_box(12, 2.8, 4, 0.5, 'admixture_k5_sorted.png', colors['output'], fontsize=8)
    
    draw_arrow((4, 3.5), (4, 3.05))
    draw_arrow((8, 3.5), (8, 3.05))
    draw_arrow((12, 3.5), (12, 3.05))
    
    # =========================================================================
    # LEGEND
    # =========================================================================
    legend_y = 1.8
    legend_items = [
        (1.5, 'Input Files', colors['input']),
        (4.5, 'Tools', colors['tool']),
        (7.5, 'Processing', colors['process']),
        (10.5, 'Analysis', colors['analysis']),
        (13.5, 'Output', colors['output']),
    ]
    
    ax.text(8, 2.2, 'Legend:', ha='center', fontsize=10, fontweight='bold')
    
    for x, label, color in legend_items:
        box = FancyBboxPatch(
            (x - 0.6, legend_y - 0.2), 1.2, 0.4,
            boxstyle="round,pad=0.02,rounding_size=0.1",
            facecolor=color,
            edgecolor=colors['border'],
            linewidth=1
        )
        ax.add_patch(box)
        ax.text(x + 1, legend_y, label, ha='left', va='center', fontsize=8)
    
    # =========================================================================
    # STATISTICS BOX
    # =========================================================================
    stats_text = (
        "Pipeline Statistics:\n"
        "• Samples: 2,504 individuals\n"
        "• Populations: 26 (5 super-populations)\n"
        "• Variants after filtering: 112,015\n"
        "• Mean Fst: 0.092 (Weighted: 0.094)"
    )
    
    props = dict(boxstyle='round,pad=0.5', facecolor='#ECEFF1', edgecolor=colors['border'], alpha=0.9)
    ax.text(14.5, 17.8, stats_text, ha='left', va='top', fontsize=8, 
            bbox=props, family='monospace')
    
    plt.tight_layout()
    plt.savefig('results/pipeline_flowchart.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print("Flowchart saved to: results/pipeline_flowchart.png")


if __name__ == '__main__':
    create_flowchart()
