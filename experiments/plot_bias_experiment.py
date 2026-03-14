#!/usr/bin/env python3
"""
Scatter plots for Experiment 1: L1 vs L2 coordinate bias.

Reads the CSV produced by `cargo run --release --example bias_experiment`
and generates scatter plots comparing first-kind (exponential) and
second-kind (additive) coordinate errors.

Outputs PGF files for direct inclusion in LaTeX documents (text is
rendered by LaTeX using your document fonts) as well as PNG previews.

Usage (from repo root):
    python experiments/plot_bias_experiment.py [path/to/bias_experiment.csv]

Default CSV path: experiments/data/bias_experiment.csv
Output:           paper/figures/bias_scatter_translation.{pgf,png}
                  paper/figures/bias_scatter_rotation.{pgf,png}
"""

import sys
import os
import numpy as np

import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt

# ── Paths (repo-relative, work from any working directory) ──────
REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(REPO_ROOT, "experiments", "data")
FIG_DIR     = os.path.join(REPO_ROOT, "paper", "figures")
DEFAULT_CSV = os.path.join(DATA_DIR, "bias_experiment.csv")

# ── Shared style ────────────────────────────────────────────────
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})


def load_csv(path):
    """Load bias_experiment.csv, return dict of arrays."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return {
        "l1_omega": data[:, 0:3],
        "l1_t":     data[:, 3:6],
        "l2_omega": data[:, 6:9],
        "l2_t":     data[:, 9:12],
    }


def scatter_pair(ax_l1, ax_l2, l1, l2, comp_x, comp_y, labels, title_prefix):
    """Draw matching scatter plots for L1 and L2 on two axes."""
    alpha = max(0.05, min(0.4, 200.0 / len(l1)))

    ax_l1.scatter(l1[:, comp_x], l1[:, comp_y], s=4, alpha=alpha, color="steelblue", rasterized=True)
    ax_l2.scatter(l2[:, comp_x], l2[:, comp_y], s=4, alpha=alpha, color="firebrick", rasterized=True)

    # Mean markers
    l1_mean = l1.mean(axis=0)
    l2_mean = l2.mean(axis=0)
    ax_l1.scatter(l1_mean[comp_x], l1_mean[comp_y], s=80, marker="+", color="black", linewidths=2, zorder=5)
    ax_l2.scatter(l2_mean[comp_x], l2_mean[comp_y], s=80, marker="+", color="black", linewidths=2, zorder=5)

    # Origin crosshair
    for ax in (ax_l1, ax_l2):
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")

    ax_l1.set_xlabel(labels[comp_x])
    ax_l1.set_ylabel(labels[comp_y])
    ax_l2.set_xlabel(labels[comp_x])
    ax_l2.set_ylabel(labels[comp_y])

    l1_bias = np.linalg.norm(l1_mean)
    l2_bias = np.linalg.norm(l2_mean)
    ax_l1.set_title(rf"{title_prefix} --- L1 (exp)" "\n" rf"$\|\mathrm{{bias}}\| = {l1_bias:.4f}$")
    ax_l2.set_title(rf"{title_prefix} --- L2 (additive)" "\n" rf"$\|\mathrm{{bias}}\| = {l2_bias:.4f}$")

    # Use same axis limits for fair comparison
    all_vals = np.concatenate([l1[:, [comp_x, comp_y]], l2[:, [comp_x, comp_y]]])
    margin = 0.05
    for c in range(2):
        lo, hi = np.percentile(all_vals[:, c], [1, 99])
        pad = (hi - lo) * margin
        lo, hi = lo - pad, hi + pad
        if c == 0:
            ax_l1.set_xlim(lo, hi)
            ax_l2.set_xlim(lo, hi)
        else:
            ax_l1.set_ylim(lo, hi)
            ax_l2.set_ylim(lo, hi)


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV
    if not os.path.isfile(csv_path):
        print(f"Error: CSV not found at {csv_path}")
        print("Run the Rust experiment first (from repo root):")
        print("  cd rust && cargo run --release --example bias_experiment && cd ..")
        print(f"Then ensure CSV is at: {DEFAULT_CSV}")
        sys.exit(1)

    os.makedirs(FIG_DIR, exist_ok=True)

    d = load_csv(csv_path)
    n = len(d["l1_t"])
    print(f"Loaded {n} samples from {csv_path}")

    # ── Figure 1: Translation scatter (the main result, shows 21× bias) ──
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    trans_labels = [r"$\delta t_1$", r"$\delta t_2$", r"$\delta t_3$"]
    pairs = [(0, 1), (0, 2), (1, 2)]
    for col, (cx, cy) in enumerate(pairs):
        scatter_pair(axes[0, col], axes[1, col],
                     d["l1_t"], d["l2_t"], cx, cy,
                     trans_labels, "Translation")

    fig.suptitle(
        rf"Experiment 1: Translation Error Scatter \quad ($n={n}$)" "\n"
        r"L1 = first-kind (exponential) $\mid$ L2 = second-kind (additive)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    for ext in ("pgf", "png"):
        path1 = os.path.join(FIG_DIR, f"bias_scatter_translation.{ext}")
        fig.savefig(path1, dpi=150)
        print(f"Saved {path1}")
    plt.close(fig)

    # ── Figure 2: Rotation scatter ──
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    rot_labels = [r"$\delta\omega_1$", r"$\delta\omega_2$", r"$\delta\omega_3$"]
    for col, (cx, cy) in enumerate(pairs):
        scatter_pair(axes[0, col], axes[1, col],
                     d["l1_omega"], d["l2_omega"], cx, cy,
                     rot_labels, "Rotation")

    fig.suptitle(
        rf"Experiment 1: Rotation Error Scatter \quad ($n={n}$)" "\n"
        r"L1 = first-kind (exponential) $\mid$ L2 = second-kind (additive)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    for ext in ("pgf", "png"):
        path2 = os.path.join(FIG_DIR, f"bias_scatter_rotation.{ext}")
        fig.savefig(path2, dpi=150)
        print(f"Saved {path2}")
    plt.close(fig)

    # ── Print summary statistics ──
    l1_t_bias = np.linalg.norm(d["l1_t"].mean(axis=0))
    l2_t_bias = np.linalg.norm(d["l2_t"].mean(axis=0))
    l1_r_bias = np.linalg.norm(d["l1_omega"].mean(axis=0))
    l2_r_bias = np.linalg.norm(d["l2_omega"].mean(axis=0))
    print(f"\nTranslation bias — L1: {l1_t_bias:.6f}  L2: {l2_t_bias:.6f}  ratio: {l2_t_bias/max(l1_t_bias,1e-15):.1f}x")
    print(f"Rotation bias    — L1: {l1_r_bias:.6f}  L2: {l2_r_bias:.6f}  ratio: {l2_r_bias/max(l1_r_bias,1e-15):.1f}x")


if __name__ == "__main__":
    main()
