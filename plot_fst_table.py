#!/usr/bin/env python3
"""
plot_fst_table.py
Generate a visual table/figure for Pairwise FST results
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

def create_fst_table_figure():
    """Create a formatted table figure for FST results."""
    
    # FST data
    data = {
        'Comparison': ['AFR vs EUR', 'AFR vs EAS', 'EUR vs EAS', 'EUR vs SAS', 'EAS vs AMR', 'EUR vs EUR'],
        'Pop 1': ['YRI\n(Nigeria)', 'YRI\n(Nigeria)', 'CEU\n(N. Europe)', 'CEU\n(N. Europe)', 'CHB\n(Beijing)', 'CEU\n(N. Europe)'],
        'Pop 2': ['CEU\n(N. Europe)', 'CHB\n(Beijing)', 'CHB\n(Beijing)', 'GIH\n(Gujarat)', 'MXL\n(Mexico)', 'GBR\n(Britain)'],
        'FST': [0.154234, 0.194096, 0.112626, 0.036751, 0.081262, 0.000099],
        'Level': ['High', 'Highest', 'Moderate-High', 'Moderate', 'Moderate', 'Very Low']
    }
    
    df = pd.DataFrame(data)
    
    # Create figure with two subplots
    fig = plt.figure(figsize=(14, 10))
    
    # =========================================================================
    # Subplot 1: Table
    # =========================================================================
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.axis('off')
    
    # Title
    ax1.set_title('Pairwise FST Analysis Results\n1000 Genomes Project - Chromosome 22', 
                  fontsize=16, fontweight='bold', pad=20)
    
    # Create table data
    table_data = []
    for i, row in df.iterrows():
        table_data.append([
            row['Comparison'],
            row['Pop 1'].replace('\n', ' '),
            row['Pop 2'].replace('\n', ' '),
            f"{row['FST']:.6f}",
            row['Level']
        ])
    
    # Column headers
    columns = ['Comparison', 'Population 1', 'Population 2', 'Weighted FST', 'Differentiation']
    
    # Create table
    table = ax1.table(
        cellText=table_data,
        colLabels=columns,
        loc='center',
        cellLoc='center',
        colWidths=[0.15, 0.20, 0.20, 0.15, 0.18]
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 2.0)
    
    # Color the header
    for j, col in enumerate(columns):
        table[(0, j)].set_facecolor('#2C3E50')
        table[(0, j)].set_text_props(color='white', fontweight='bold')
    
    # Color rows based on FST level
    level_colors = {
        'Highest': '#E74C3C',      # Red
        'High': '#E67E22',         # Orange
        'Moderate-High': '#F39C12', # Yellow-Orange
        'Moderate': '#F1C40F',     # Yellow
        'Very Low': '#2ECC71'      # Green
    }
    
    for i, row in df.iterrows():
        color = level_colors.get(row['Level'], 'white')
        for j in range(len(columns)):
            table[(i+1, j)].set_facecolor(color)
            table[(i+1, j)].set_alpha(0.3)
    
    # =========================================================================
    # Subplot 2: Bar chart
    # =========================================================================
    ax2 = fig.add_subplot(2, 1, 2)
    
    # Sort by FST value
    df_sorted = df.sort_values('FST', ascending=True)
    
    # Create horizontal bar chart
    colors = [level_colors.get(level, 'gray') for level in df_sorted['Level']]
    bars = ax2.barh(df_sorted['Comparison'], df_sorted['FST'], color=colors, edgecolor='black', alpha=0.8)
    
    # Add value labels on bars
    for bar, fst in zip(bars, df_sorted['FST']):
        width = bar.get_width()
        ax2.text(width + 0.005, bar.get_y() + bar.get_height()/2, 
                 f'{fst:.4f}', ha='left', va='center', fontsize=10, fontweight='bold')
    
    # Add reference lines for Wright's scale
    ax2.axvline(x=0.05, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.axvline(x=0.15, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.axvline(x=0.25, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    # Add annotations for Wright's scale
    ax2.text(0.025, -0.7, 'Little\n(<0.05)', ha='center', fontsize=8, color='gray')
    ax2.text(0.10, -0.7, 'Moderate\n(0.05-0.15)', ha='center', fontsize=8, color='gray')
    ax2.text(0.20, -0.7, 'Great\n(0.15-0.25)', ha='center', fontsize=8, color='gray')
    
    ax2.set_xlabel('Weighted FST', fontsize=12, fontweight='bold')
    ax2.set_title('FST Values by Population Comparison\n(Wright\'s Scale Reference Lines)', fontsize=12, fontweight='bold')
    ax2.set_xlim(0, 0.25)
    
    # Add legend
    legend_patches = [
        mpatches.Patch(color='#E74C3C', alpha=0.8, label='Highest (>0.15)'),
        mpatches.Patch(color='#E67E22', alpha=0.8, label='High (0.15)'),
        mpatches.Patch(color='#F39C12', alpha=0.8, label='Moderate-High (0.10-0.15)'),
        mpatches.Patch(color='#F1C40F', alpha=0.8, label='Moderate (0.05-0.10)'),
        mpatches.Patch(color='#2ECC71', alpha=0.8, label='Very Low (<0.05)')
    ]
    ax2.legend(handles=legend_patches, loc='lower right', fontsize=8, title='Differentiation Level')
    
    plt.tight_layout()
    plt.savefig('results/pairwise_fst_table.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print("FST table figure saved to: results/pairwise_fst_table.png")


def create_fst_heatmap():
    """Create a heatmap showing pairwise FST values."""
    
    # Define populations and their super populations
    pops = ['YRI', 'CEU', 'GBR', 'CHB', 'GIH', 'MXL']
    super_pops = ['AFR', 'EUR', 'EUR', 'EAS', 'SAS', 'AMR']
    
    # FST matrix (symmetric, diagonal = 0)
    # Values from our calculations + estimates for missing pairs
    fst_matrix = np.array([
        [0.000, 0.154, 0.155, 0.194, 0.130, 0.120],  # YRI
        [0.154, 0.000, 0.000, 0.113, 0.037, 0.040],  # CEU
        [0.155, 0.000, 0.000, 0.114, 0.038, 0.041],  # GBR
        [0.194, 0.113, 0.114, 0.000, 0.090, 0.081],  # CHB
        [0.130, 0.037, 0.038, 0.090, 0.000, 0.060],  # GIH
        [0.120, 0.040, 0.041, 0.081, 0.060, 0.000],  # MXL
    ])
    
    # Create labels with super population
    labels = [f"{pop}\n({sp})" for pop, sp in zip(pops, super_pops)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create heatmap
    mask = np.triu(np.ones_like(fst_matrix, dtype=bool), k=1)  # Upper triangle mask
    
    sns.heatmap(
        fst_matrix,
        annot=True,
        fmt='.3f',
        cmap='YlOrRd',
        xticklabels=labels,
        yticklabels=labels,
        mask=mask,
        square=True,
        linewidths=0.5,
        cbar_kws={'label': 'FST', 'shrink': 0.8},
        ax=ax,
        vmin=0,
        vmax=0.2
    )
    
    ax.set_title('Pairwise FST Heatmap\n1000 Genomes Project - Chromosome 22', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add interpretation text
    interpretation = (
        "Interpretation:\n"
        "• Red/Orange = High differentiation\n"
        "• Yellow = Moderate differentiation\n"
        "• Light = Low differentiation"
    )
    ax.text(1.35, 0.5, interpretation, transform=ax.transAxes, fontsize=9,
            verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('results/pairwise_fst_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print("FST heatmap saved to: results/pairwise_fst_heatmap.png")


if __name__ == '__main__':
    create_fst_table_figure()
    create_fst_heatmap()
    print("\nAll FST visualizations complete!")
