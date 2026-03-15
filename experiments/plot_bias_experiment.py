#!/usr/bin/env python3
"""
Scatter plots for Experiment 1: L1 vs L2 coordinate bias.

Reads the CSV produced by `cargo run --release --example bias_experiment`
and generates two compact figures (translation, rotation), each showing
the 2D projection that maximises visible L2 bias. L1 (exponential) on
the left, L2 (additive) on the right — matching the paper's Fig 1–2
format.

Outputs PGF files for direct inclusion in LaTeX documents (text is
rendered by LaTeX using your document fonts) as well as PNG previews.

Usage (from repo root):
    python experiments/plot_bias_experiment.py [path/to/data.csv]
    python experiments/plot_bias_experiment.py --all [path/to/data.csv]

Default CSV path: experiments/data/bias_experiment.csv

Without --all (default):
    paper/figures/bias_scatter_translation.{pgf,png}   (1×2, best projection)
    paper/figures/bias_scatter_rotation.{pgf,png}       (1×2, best projection)

With --all:
    Also produces the full 2×3 grids (all axis pairs):
    paper/figures/bias_scatter_translation_all.{pgf,png}
    paper/figures/bias_scatter_rotation_all.{pgf,png}
"""

import sys
import os
import argparse
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

IEEE_COL = 3.5   # inches, single-column
IEEE_DBL = 7.16  # inches, double-column


def load_csv(path):
    """Load bias_experiment.csv, return dict of arrays."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return {
        "l1_omega": data[:, 0:3],
        "l1_t":     data[:, 3:6],
        "l2_omega": data[:, 6:9],
        "l2_t":     data[:, 9:12],
    }


def pick_max_bias_axes(l2_data):
    """Return (comp_x, comp_y) — the two component indices whose L2
    mean has the largest absolute values, so the 2D projection shows
    the most visible bias offset."""
    mean = np.abs(l2_data.mean(axis=0))
    # Pick the two components with largest absolute mean
    ranked = np.argsort(mean)[::-1]
    cx, cy = sorted(ranked[:2])  # keep consistent axis ordering
    return cx, cy


def scatter_panel(ax, data, comp_x, comp_y, labels, color, title, bias_norm):
    """Draw a single scatter panel."""
    alpha = max(0.05, min(0.4, 200.0 / len(data)))

    ax.scatter(data[:, comp_x], data[:, comp_y],
               s=4, alpha=alpha, color=color, rasterized=True)

    # Mean marker
    mean = data.mean(axis=0)
    ax.scatter(mean[comp_x], mean[comp_y],
               s=80, marker="+", color="black", linewidths=2, zorder=5)

    # Origin crosshair
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")

    ax.set_xlabel(labels[comp_x])
    ax.set_ylabel(labels[comp_y])
    ax.set_title(rf"{title}" "\n" rf"$\|\mathrm{{bias}}\| = {bias_norm:.4f}$")


def make_figure(l1, l2, labels, name, fig_name):
    """Create a 1×2 figure (L1 left, L2 right) for the best projection."""
    os.makedirs(FIG_DIR, exist_ok=True)

    cx, cy = pick_max_bias_axes(l2)
    n = len(l1)

    l1_bias = np.linalg.norm(l1.mean(axis=0))
    l2_bias = np.linalg.norm(l2.mean(axis=0))

    fig, (ax_l1, ax_l2) = plt.subplots(1, 2, figsize=(IEEE_DBL, 3.0))

    scatter_panel(ax_l1, l1, cx, cy, labels,
                  color="steelblue",
                  title=rf"{name} --- L1 (exp)",
                  bias_norm=l1_bias)
    scatter_panel(ax_l2, l2, cx, cy, labels,
                  color="firebrick",
                  title=rf"{name} --- L2 (additive)",
                  bias_norm=l2_bias)

    # Matched axis limits for fair comparison
    all_vals = np.concatenate([l1[:, [cx, cy]], l2[:, [cx, cy]]])
    for i, setter_l1, setter_l2 in [
        (0, ax_l1.set_xlim, ax_l2.set_xlim),
        (1, ax_l1.set_ylim, ax_l2.set_ylim),
    ]:
        lo, hi = np.percentile(all_vals[:, i], [1, 99])
        pad = (hi - lo) * 0.05
        setter_l1(lo - pad, hi + pad)
        setter_l2(lo - pad, hi + pad)

    fig.tight_layout()
    for ext in ("pgf", "png"):
        path = os.path.join(FIG_DIR, f"{fig_name}.{ext}")
        fig.savefig(path, dpi=150)
        print(f"Saved {path}")
    plt.close(fig)


def make_figure_all(l1, l2, labels, name, fig_name):
    """Create a 2×3 figure showing all axis pairs. Top row = L1, bottom = L2."""
    os.makedirs(FIG_DIR, exist_ok=True)

    pairs = [(0, 1), (0, 2), (1, 2)]
    l1_bias = np.linalg.norm(l1.mean(axis=0))
    l2_bias = np.linalg.norm(l2.mean(axis=0))

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    for col, (cx, cy) in enumerate(pairs):
        scatter_panel(axes[0, col], l1, cx, cy, labels,
                      color="steelblue",
                      title=rf"{name} --- L1 (exp)",
                      bias_norm=l1_bias)
        scatter_panel(axes[1, col], l2, cx, cy, labels,
                      color="firebrick",
                      title=rf"{name} --- L2 (additive)",
                      bias_norm=l2_bias)

        # Matched axis limits per column
        all_vals = np.concatenate([l1[:, [cx, cy]], l2[:, [cx, cy]]])
        for i, row in enumerate([0, 1]):
            for j, setter in enumerate([
                lambda lo, hi, r=row: axes[r, col].set_xlim(lo, hi),
                lambda lo, hi, r=row: axes[r, col].set_ylim(lo, hi),
            ]):
                lo, hi = np.percentile(all_vals[:, j], [1, 99])
                pad = (hi - lo) * 0.05
                setter(lo - pad, hi + pad)

    n = len(l1)
    fig.suptitle(
        rf"Experiment 1: {name} Error Scatter \quad ($n={n}$)" "\n"
        r"L1 = first-kind (exponential) $\mid$ L2 = second-kind (additive)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    for ext in ("pgf", "png"):
        path = os.path.join(FIG_DIR, f"{fig_name}_all.{ext}")
        fig.savefig(path, dpi=150)
        print(f"Saved {path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Bias scatter plots for Experiment 1.")
    parser.add_argument("csv", nargs="?", default=DEFAULT_CSV,
                        help="Path to bias_experiment.csv")
    parser.add_argument("--all", action="store_true",
                        help="Also produce 2×3 grids with all axis pairs")
    args = parser.parse_args()

    csv_path = args.csv
    if not os.path.isfile(csv_path):
        print(f"Error: CSV not found at {csv_path}")
        print("Run the Rust experiment first (from repo root):")
        print("  cd rust && cargo run --release --example bias_experiment && cd ..")
        print(f"Then ensure CSV is at: {DEFAULT_CSV}")
        sys.exit(1)

    d = load_csv(csv_path)
    n = len(d["l1_t"])
    print(f"Loaded {n} samples from {csv_path}")

    trans_labels = [r"$\delta t_1$", r"$\delta t_2$", r"$\delta t_3$"]
    rot_labels   = [r"$\delta\omega_1$", r"$\delta\omega_2$", r"$\delta\omega_3$"]

    # Always produce the compact 1×2 figures
    make_figure(d["l1_t"], d["l2_t"], trans_labels,
                "Translation", "bias_scatter_translation")
    make_figure(d["l1_omega"], d["l2_omega"], rot_labels,
                "Rotation", "bias_scatter_rotation")

    # Optionally produce the full 2×3 grids
    if args.all:
        make_figure_all(d["l1_t"], d["l2_t"], trans_labels,
                        "Translation", "bias_scatter_translation")
        make_figure_all(d["l1_omega"], d["l2_omega"], rot_labels,
                        "Rotation", "bias_scatter_rotation")

    # ── Summary statistics ──
    l1_t_bias = np.linalg.norm(d["l1_t"].mean(axis=0))
    l2_t_bias = np.linalg.norm(d["l2_t"].mean(axis=0))
    l1_r_bias = np.linalg.norm(d["l1_omega"].mean(axis=0))
    l2_r_bias = np.linalg.norm(d["l2_omega"].mean(axis=0))
    print(f"\nTranslation bias — L1: {l1_t_bias:.6f}  L2: {l2_t_bias:.6f}  ratio: {l2_t_bias/max(l1_t_bias,1e-15):.1f}x")
    print(f"Rotation bias    — L1: {l1_r_bias:.6f}  L2: {l2_r_bias:.6f}  ratio: {l2_r_bias/max(l1_r_bias,1e-15):.1f}x")

    # Report which axes were selected
    cx_t, cy_t = pick_max_bias_axes(d["l2_t"])
    cx_r, cy_r = pick_max_bias_axes(d["l2_omega"])
    print(f"\nSelected translation axes: components {cx_t},{cy_t}")
    print(f"Selected rotation axes:    components {cx_r},{cy_r}")


if __name__ == "__main__":
    main()