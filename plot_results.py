#!/usr/bin/env python3
"""Generate PCA, ADMIXTURE, and CV diagnostic plots for the chr22 workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_hex

sns.set(style="whitegrid")

RESULTS_DIR = Path("results")
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

PANEL_FILE = Path("integrated_call_samples_v3.20130502.ALL.panel")
EIGENVEC_FILE = RESULTS_DIR / "chr22_pca.eigenvec"
EIGENVAL_FILE = RESULTS_DIR / "chr22_pca.eigenval"
FAM_FILE = RESULTS_DIR / "chr22_pruned.fam"
CV_FILE = RESULTS_DIR / "admixture_cv_errors.tsv"
K_RANGE = range(2, 7)

# Super-population colors (colorblind-friendly)
SUPER_COLORS = {
    "AFR": "#d73027",  # Red
    "EUR": "#4575b4",  # Blue
    "EAS": "#1a9850",  # Green
    "SAS": "#ff7f00",  # Orange
    "AMR": "#8c6bb1",  # Purple
}

SUPER_PALETTES = {
    "AFR": "Reds",
    "EUR": "Blues",
    "EAS": "Greens",
    "SAS": "Oranges",
    "AMR": "Purples",
}

SUPER_ORDER = ["AFR", "EUR", "EAS", "SAS", "AMR"]

# Improved ancestry color palettes for ADMIXTURE
# Using distinct, colorblind-friendly palettes for different K values
ANCESTRY_PALETTES = {
    2: ["#e41a1c", "#377eb8"],
    3: ["#e41a1c", "#377eb8", "#4daf4a"],
    4: ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"],
    5: ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"],
    6: ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"],
    7: ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf"],
    8: ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"],
}


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file missing: {path}")


def load_metadata() -> pd.DataFrame:
    ensure_exists(PANEL_FILE)
    panel = pd.read_csv(PANEL_FILE, sep="\t")
    required = {"sample", "pop", "super_pop"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"Panel file missing columns: {', '.join(sorted(missing))}")
    return panel


def load_pca_data() -> pd.DataFrame:
    ensure_exists(EIGENVEC_FILE)
    columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 21)]
    return pd.read_csv(EIGENVEC_FILE, sep=r"\s+", names=columns, header=None)


def load_eigenvalues() -> np.ndarray:
    ensure_exists(EIGENVAL_FILE)
    return pd.read_csv(EIGENVAL_FILE, sep=r"\s+", header=None)[0].values


def compute_pc_variance(eigenvals: np.ndarray) -> np.ndarray:
    total = eigenvals.sum()
    if total == 0:
        return np.zeros_like(eigenvals)
    return eigenvals / total


def build_population_palette(metadata: pd.DataFrame) -> Dict[str, str]:
    palette: Dict[str, str] = {}
    for super_pop, subset in metadata.groupby("super_pop"):
        pops = sorted(subset["pop"].unique())
        base_palette = SUPER_PALETTES.get(super_pop, "husl")
        colors = sns.color_palette(base_palette, len(pops))
        for pop, color in zip(pops, colors):
            palette[pop] = to_hex(color)
    return palette


def plot_pca_super_pop(
    pca_df: pd.DataFrame, 
    metadata: pd.DataFrame, 
    var_pct: np.ndarray, 
    output: Path
) -> None:
    """Plot PCA colored by super-population."""
    merged = pca_df.merge(metadata, left_on="IID", right_on="sample", how="left")
    plt.figure(figsize=(8, 6))
    for super_pop in SUPER_ORDER:
        color = SUPER_COLORS.get(super_pop, "#7f7f7f")
        subset = merged[merged["super_pop"] == super_pop]
        if subset.empty:
            continue
        plt.scatter(subset["PC1"], subset["PC2"], c=color, label=super_pop, s=18, alpha=0.7)
    plt.xlabel(f"PC1 ({var_pct[0] * 100:.1f}% variance explained)")
    plt.ylabel(f"PC2 ({var_pct[1] * 100:.1f}% variance explained)")
    plt.title("PCA colored by super population")
    plt.legend(title="Super population", loc="best")
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def plot_pca_detailed(
    pca_df: pd.DataFrame, 
    metadata: pd.DataFrame, 
    var_pct: np.ndarray, 
    palette: Dict[str, str], 
    output: Path
) -> None:
    """Plot PCA colored by sub-population."""
    merged = pca_df.merge(metadata, left_on="IID", right_on="sample", how="left")
    merged = merged.sort_values(["super_pop", "pop"])
    plt.figure(figsize=(10, 8))
    for pop, subset in merged.groupby("pop"):
        color = palette.get(pop, "#7f7f7f")
        plt.scatter(subset["PC1"], subset["PC2"], c=[color], s=15, alpha=0.7, label=pop)
    plt.xlabel(f"PC1 ({var_pct[0] * 100:.1f}% variance explained)")
    plt.ylabel(f"PC2 ({var_pct[1] * 100:.1f}% variance explained)")
    plt.title("PCA colored by population")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def plot_pca_multi_pc(
    pca_df: pd.DataFrame,
    metadata: pd.DataFrame,
    var_pct: np.ndarray,
    output: Path,
    pc_pairs: List[Tuple[int, int]] = [(1, 2), (1, 3), (2, 3)]
) -> None:
    """
    Plot multiple PC combinations in a single figure.
    
    Args:
        pc_pairs: List of (PCx, PCy) tuples to plot
    """
    merged = pca_df.merge(metadata, left_on="IID", right_on="sample", how="left")
    
    n_plots = len(pc_pairs)
    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5))
    if n_plots == 1:
        axes = [axes]
    
    for ax, (pc_x, pc_y) in zip(axes, pc_pairs):
        for super_pop in SUPER_ORDER:
            color = SUPER_COLORS.get(super_pop, "#7f7f7f")
            subset = merged[merged["super_pop"] == super_pop]
            if subset.empty:
                continue
            ax.scatter(
                subset[f"PC{pc_x}"], 
                subset[f"PC{pc_y}"], 
                c=color, 
                label=super_pop, 
                s=15, 
                alpha=0.7
            )
        
        ax.set_xlabel(f"PC{pc_x} ({var_pct[pc_x - 1] * 100:.1f}%)")
        ax.set_ylabel(f"PC{pc_y} ({var_pct[pc_y - 1] * 100:.1f}%)")
        ax.set_title(f"PC{pc_x} vs PC{pc_y}")
        ax.grid(alpha=0.3)
    
    # Single legend for all subplots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.99, 0.99), 
               title="Super Pop", frameon=True, fontsize=9)
    
    plt.suptitle("PCA - Multiple Principal Component Views", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def plot_pca_3d_projection(
    pca_df: pd.DataFrame,
    metadata: pd.DataFrame,
    var_pct: np.ndarray,
    output: Path
) -> None:
    """Create a 2x2 grid showing PC1-2, PC1-3, PC2-3, and scree plot."""
    merged = pca_df.merge(metadata, left_on="IID", right_on="sample", how="left")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    pc_combos = [(1, 2), (1, 3), (2, 3)]
    
    for ax, (pc_x, pc_y) in zip(axes.flat[:3], pc_combos):
        for super_pop in SUPER_ORDER:
            color = SUPER_COLORS.get(super_pop, "#7f7f7f")
            subset = merged[merged["super_pop"] == super_pop]
            if subset.empty:
                continue
            ax.scatter(
                subset[f"PC{pc_x}"], 
                subset[f"PC{pc_y}"], 
                c=color, 
                label=super_pop, 
                s=12, 
                alpha=0.6
            )
        ax.set_xlabel(f"PC{pc_x} ({var_pct[pc_x - 1] * 100:.1f}%)", fontsize=10)
        ax.set_ylabel(f"PC{pc_y} ({var_pct[pc_y - 1] * 100:.1f}%)", fontsize=10)
        ax.set_title(f"PC{pc_x} vs PC{pc_y}", fontsize=11)
        ax.grid(alpha=0.3)
    
    # Scree plot in the 4th panel
    ax_scree = axes[1, 1]
    n_pcs = min(10, len(var_pct))
    cumulative_var = np.cumsum(var_pct[:n_pcs]) * 100
    individual_var = var_pct[:n_pcs] * 100
    
    x = np.arange(1, n_pcs + 1)
    ax_scree.bar(x, individual_var, color='steelblue', alpha=0.7, label='Individual')
    ax_scree.plot(x, cumulative_var, 'ro-', markersize=6, label='Cumulative')
    ax_scree.set_xlabel("Principal Component", fontsize=10)
    ax_scree.set_ylabel("Variance Explained (%)", fontsize=10)
    ax_scree.set_title("Scree Plot", fontsize=11)
    ax_scree.legend(loc='center right', fontsize=9)
    ax_scree.set_xticks(x)
    ax_scree.grid(alpha=0.3, axis='y')
    
    # Common legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.02),
               ncol=5, title="Super Population", frameon=True, fontsize=9)
    
    plt.suptitle("Principal Component Analysis Summary", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def get_ancestry_colors(k: int) -> List[str]:
    """Get colorblind-friendly colors for K ancestral components."""
    if k in ANCESTRY_PALETTES:
        return ANCESTRY_PALETTES[k]
    # Fallback: generate colors using a perceptually uniform colormap
    return [to_hex(c) for c in sns.color_palette("husl", k)]


def load_admixture(k: int) -> pd.DataFrame:
    q_path = RESULTS_DIR / f"chr22_pruned.{k}.Q"
    ensure_exists(q_path)
    ensure_exists(FAM_FILE)
    q_data = pd.read_csv(q_path, sep=r"\s+", header=None)
    fam = pd.read_csv(FAM_FILE, sep=r"\s+", header=None, 
                      names=["FID", "IID", "PID", "MID", "Sex", "Phenotype"])
    if len(q_data) != len(fam):
        raise ValueError(f"ADMIXTURE Q file for K={k} does not align with FAM file")
    q_data.columns = [f"Ancestry_{i+1}" for i in range(q_data.shape[1])]
    q_data["IID"] = fam["IID"].values
    return q_data


def plot_admixture_sorted(
    q_df: pd.DataFrame, 
    metadata: pd.DataFrame, 
    k: int, 
    output: Path
) -> None:
    """Plot ADMIXTURE results with improved color scheme."""
    merged = q_df.merge(metadata, left_on="IID", right_on="sample", how="left")
    merged["super_order"] = merged["super_pop"].map(
        {sp: idx for idx, sp in enumerate(SUPER_ORDER)}
    ).fillna(99)
    merged = merged.sort_values(["super_order", "pop", "IID"]).reset_index(drop=True)
    
    ancestry_cols = [col for col in merged.columns if col.startswith("Ancestry_")]
    ancestry_colors = get_ancestry_colors(k)

    plt.figure(figsize=(max(14, len(merged) / 60), 4))
    bottom = np.zeros(len(merged))
    
    for i, col in enumerate(ancestry_cols):
        color = ancestry_colors[i] if i < len(ancestry_colors) else "#7f7f7f"
        plt.bar(
            range(len(merged)), 
            merged[col], 
            bottom=bottom, 
            width=1.0, 
            color=color, 
            edgecolor="none",
            label=f"Ancestry {i+1}"
        )
        bottom += merged[col].values

    # Add population separators and labels
    ticks = []
    labels = []
    prev_end = -0.5
    
    for super_pop in SUPER_ORDER:
        group = merged[merged["super_pop"] == super_pop]
        if group.empty:
            continue
        indices = group.index.to_list()
        ticks.append((indices[0] + indices[-1]) / 2)
        labels.append(super_pop)
        # Draw separator line
        plt.axvline(indices[-1] + 0.5, color="white", linewidth=1.5)

    plt.xticks(ticks, labels, fontsize=10)
    plt.ylabel("Ancestry proportion", fontsize=11)
    plt.xlabel("Super population", fontsize=11)
    plt.title(f"ADMIXTURE K={k}", fontsize=12, fontweight='bold')
    plt.xlim(-0.5, len(merged) - 0.5)
    plt.ylim(0, 1)
    
    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(1.12, 1), fontsize=8, frameon=True)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def plot_admixture_combined(
    metadata: pd.DataFrame,
    k_range: range,
    output: Path
) -> None:
    """Plot all K values in a single stacked figure for comparison."""
    valid_k = []
    q_data_list = []
    
    for k in k_range:
        try:
            q_df = load_admixture(k)
            valid_k.append(k)
            q_data_list.append(q_df)
        except (FileNotFoundError, ValueError):
            continue
    
    if not valid_k:
        print("No valid ADMIXTURE results found for combined plot")
        return
    
    n_k = len(valid_k)
    fig, axes = plt.subplots(n_k, 1, figsize=(14, 2.5 * n_k), sharex=True)
    if n_k == 1:
        axes = [axes]
    
    for ax, k, q_df in zip(axes, valid_k, q_data_list):
        merged = q_df.merge(metadata, left_on="IID", right_on="sample", how="left")
        merged["super_order"] = merged["super_pop"].map(
            {sp: idx for idx, sp in enumerate(SUPER_ORDER)}
        ).fillna(99)
        merged = merged.sort_values(["super_order", "pop", "IID"]).reset_index(drop=True)
        
        ancestry_cols = [col for col in merged.columns if col.startswith("Ancestry_")]
        ancestry_colors = get_ancestry_colors(k)
        
        bottom = np.zeros(len(merged))
        for i, col in enumerate(ancestry_cols):
            color = ancestry_colors[i] if i < len(ancestry_colors) else "#7f7f7f"
            ax.bar(
                range(len(merged)), 
                merged[col], 
                bottom=bottom, 
                width=1.0, 
                color=color, 
                edgecolor="none"
            )
            bottom += merged[col].values
        
        # Add separators
        for super_pop in SUPER_ORDER:
            group = merged[merged["super_pop"] == super_pop]
            if group.empty:
                continue
            indices = group.index.to_list()
            ax.axvline(indices[-1] + 0.5, color="white", linewidth=1)
        
        ax.set_ylabel(f"K={k}", fontsize=11, fontweight='bold')
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.5, len(merged) - 0.5)
        ax.set_yticks([0, 0.5, 1])
    
    # Add x-axis labels only on bottom plot
    ticks = []
    labels = []
    merged_last = q_data_list[-1].merge(metadata, left_on="IID", right_on="sample", how="left")
    merged_last["super_order"] = merged_last["super_pop"].map(
        {sp: idx for idx, sp in enumerate(SUPER_ORDER)}
    ).fillna(99)
    merged_last = merged_last.sort_values(["super_order", "pop", "IID"]).reset_index(drop=True)
    
    for super_pop in SUPER_ORDER:
        group = merged_last[merged_last["super_pop"] == super_pop]
        if group.empty:
            continue
        indices = group.index.to_list()
        ticks.append((indices[0] + indices[-1]) / 2)
        labels.append(super_pop)
    
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels, fontsize=10)
    axes[-1].set_xlabel("Super population", fontsize=11)
    
    plt.suptitle("ADMIXTURE Analysis (K=2 to K=6)", fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def plot_cv_errors(cv_file: Path, output: Path) -> None:
    """Plot ADMIXTURE cross-validation errors with best K highlighted."""
    ensure_exists(cv_file)
    cv_data = pd.read_csv(cv_file, sep="\t")
    cv_data = cv_data.replace({"CV_Error": {"NA": np.nan}})
    cv_data["CV_Error"] = cv_data["CV_Error"].astype(float)
    cv_data = cv_data.dropna(subset=["CV_Error"])
    
    if cv_data.empty:
        print("No valid CV error data found")
        return
    
    plt.figure(figsize=(7, 5))
    plt.plot(cv_data["K"], cv_data["CV_Error"], marker="o", color="#2c3e50", 
             linewidth=2, markersize=8, label="CV Error")
    
    # Highlight best K
    best_idx = cv_data["CV_Error"].idxmin()
    best_k = cv_data.loc[best_idx, "K"]
    best_cv = cv_data.loc[best_idx, "CV_Error"]
    plt.scatter([best_k], [best_cv], color="#e74c3c", s=150, zorder=5, 
                edgecolor='white', linewidth=2)
    plt.annotate(
        f"Best K={int(best_k)}\nCV={best_cv:.5f}",
        xy=(best_k, best_cv),
        xytext=(best_k + 0.3, best_cv + 0.002),
        fontsize=10,
        ha='left',
        arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=1.5)
    )
    
    plt.xlabel("K (Number of ancestral populations)", fontsize=11)
    plt.ylabel("Cross-validation error", fontsize=11)
    plt.title("ADMIXTURE Cross-Validation Error", fontsize=13, fontweight='bold')
    plt.xticks(cv_data["K"].astype(int))
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")


def main() -> None:
    print("=" * 60)
    print("Population Genetics Visualization Pipeline")
    print("=" * 60)
    print()
    
    try:
        metadata = load_metadata()
        pca_df = load_pca_data()
        eigenvals = load_eigenvalues()
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return

    variance = compute_pc_variance(eigenvals)
    pop_palette = build_population_palette(metadata)

    print("\n--- Generating PCA plots ---")
    
    # Basic PCA plots
    plot_pca_super_pop(
        pca_df, metadata, variance, 
        FIGURES_DIR / "pca_super_pop.png"
    )
    plot_pca_detailed(
        pca_df, metadata, variance, pop_palette, 
        FIGURES_DIR / "pca_detailed_pop.png"
    )
    
    # Multi-PC view (PC1-2, PC1-3, PC2-3)
    plot_pca_multi_pc(
        pca_df, metadata, variance,
        FIGURES_DIR / "pca_multi_pc.png",
        pc_pairs=[(1, 2), (1, 3), (2, 3)]
    )
    
    # Comprehensive 2x2 grid with scree plot
    plot_pca_3d_projection(
        pca_df, metadata, variance,
        FIGURES_DIR / "pca_summary_grid.png"
    )

    print("\n--- Generating ADMIXTURE plots ---")
    
    # Individual K plots
    for k in K_RANGE:
        try:
            q_df = load_admixture(k)
        except FileNotFoundError:
            print(f"  Skipping ADMIXTURE K={k}: Q file not found")
            continue
        except ValueError as exc:
            print(f"  Skipping ADMIXTURE K={k}: {exc}")
            continue
        plot_admixture_sorted(
            q_df, metadata, k, 
            FIGURES_DIR / f"admixture_k{k}.png"
        )
    
    # Combined ADMIXTURE plot (all K values stacked)
    plot_admixture_combined(
        metadata, K_RANGE,
        FIGURES_DIR / "admixture_combined.png"
    )

    print("\n--- Generating CV error plot ---")
    try:
        plot_cv_errors(CV_FILE, FIGURES_DIR / "admixture_cv_error.png")
    except FileNotFoundError:
        print("  CV summary not found; skipping CV error plot")

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print(f"Outputs stored in: {FIGURES_DIR}")
    print("=" * 60)
    
    # List generated files
    print("\nGenerated files:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  • {f.name}")


if __name__ == "__main__":
    main()
