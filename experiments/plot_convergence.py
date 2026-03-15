#!/usr/bin/env python3
"""
Convergence basin plots for the pose inference sweep experiment.

Reads the CSV produced by `cargo run --release --example pose_inference --sweep`
and generates figures showing convergence traces across camera configurations,
arc spans, perturbation scales, and Laplace vs saddlepoint methods.

Outputs PGF files for LaTeX inclusion and PNG previews.

Usage (from repo root):
    python experiments/plot_convergence.py [path/to/data.csv]

Default CSV path: experiments/data/convergence_sweep.csv
"""

import sys
import os
import argparse
import numpy as np

import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Paths ────────────────────────────────────────────────────────
REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(REPO_ROOT, "experiments", "data")
FIG_DIR     = os.path.join(REPO_ROOT, "paper", "figures")
DEFAULT_CSV = os.path.join(DATA_DIR, "convergence_sweep.csv")

# ── Style ────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 7,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

IEEE_COL = 3.5
IEEE_DBL = 7.16

METHOD_COLORS = {"Laplace": "steelblue", "SP": "firebrick"}
METHOD_STYLES = {"Laplace": "-", "SP": "--"}
CAM_MARKERS   = {1: "o", 2: "s", 4: "D", 8: "^"}


def load_csv(path):
    """Load convergence_sweep.csv into a structured dict."""
    import csv
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "n_cam": int(r["n_cam"]),
                "arc_deg": float(r["arc_deg"]),
                "pert_scale": float(r["pert_scale"]),
                "method": r["method"],
                "seed": int(r["seed"]),
                "iter": int(r["iter"]),
                "obj": float(r["obj"]),
                "gnorm": float(r["gnorm"]),
                "final_rot_err": float(r["final_rot_err"]),
                "final_trans_err": float(r["final_trans_err"]),
            })
    return rows


def group_traces(rows):
    """Group rows into traces keyed by (n_cam, arc_deg, pert_scale, method, seed)."""
    traces = {}
    for r in rows:
        key = (r["n_cam"], r["arc_deg"], r["pert_scale"], r["method"], r["seed"])
        traces.setdefault(key, []).append(r)
    # Sort each trace by iteration
    for k in traces:
        traces[k].sort(key=lambda r: r["iter"])
    return traces


def fig_convergence_traces(traces, out_name="convergence_traces"):
    """
    Grid: rows = arc_deg, columns = n_cam.
    Each panel shows objective vs iteration for all perturbation scales,
    with Laplace (solid) and SP (dashed).
    """
    os.makedirs(FIG_DIR, exist_ok=True)

    cam_counts = sorted({k[0] for k in traces})
    arc_degs   = sorted({k[1] for k in traces})
    pert_scales = sorted({k[2] for k in traces})

    n_rows = len(arc_degs)
    n_cols = len(cam_counts)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(IEEE_DBL, 2.2 * n_rows),
                             squeeze=False, sharex=True)

    # Color map for perturbation scales
    cmap = plt.cm.viridis
    pert_colors = {p: cmap(i / max(len(pert_scales) - 1, 1))
                   for i, p in enumerate(pert_scales)}

    for ri, arc_deg in enumerate(arc_degs):
        for ci, n_cam in enumerate(cam_counts):
            ax = axes[ri, ci]

            for pert_scale in pert_scales:
                for method in ["Laplace", "SP"]:
                    # Collect all seeds
                    objs_by_seed = []
                    for seed in range(10):  # generous range
                        key = (n_cam, arc_deg, pert_scale, method, seed)
                        if key not in traces:
                            continue
                        tr = traces[key]
                        iters = [r["iter"] for r in tr]
                        objs = [r["obj"] for r in tr]
                        objs_by_seed.append((iters, objs))

                    if not objs_by_seed:
                        continue

                    # Plot mean trace
                    max_len = max(len(o[0]) for o in objs_by_seed)
                    for iters, objs in objs_by_seed:
                        ax.plot(iters, objs,
                                color=pert_colors[pert_scale],
                                linestyle=METHOD_STYLES[method],
                                alpha=0.4, linewidth=0.8)

            ax.set_yscale("symlog", linthresh=1e-2)
            if ri == n_rows - 1:
                ax.set_xlabel("Iteration")
            if ci == 0:
                ax.set_ylabel(rf"Objective ($\mathrm{{arc}}={arc_deg:.0f}^\circ$)")
            if ri == 0:
                ax.set_title(rf"{n_cam} cam{'s' if n_cam > 1 else ''}")

    # Legend
    legend_elements = []
    for p in pert_scales:
        legend_elements.append(Line2D([0], [0], color=pert_colors[p], linewidth=1.5,
                                      label=rf"$\sigma_\mathrm{{pert}}={p}$"))
    legend_elements.append(Line2D([0], [0], color="gray", linestyle="-",
                                  linewidth=1.5, label="Laplace"))
    legend_elements.append(Line2D([0], [0], color="gray", linestyle="--",
                                  linewidth=1.5, label="Saddlepoint"))

    fig.legend(handles=legend_elements, loc="upper center",
               ncol=len(pert_scales) + 2, bbox_to_anchor=(0.5, 1.02),
               frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("pgf", "png"):
        path = os.path.join(FIG_DIR, f"{out_name}.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


def fig_final_error(traces, out_name="convergence_final_error"):
    """
    Bar/scatter plot: final pose error vs perturbation scale,
    grouped by camera count and method.
    """
    os.makedirs(FIG_DIR, exist_ok=True)

    cam_counts  = sorted({k[0] for k in traces})
    arc_degs    = sorted({k[1] for k in traces})
    pert_scales = sorted({k[2] for k in traces})

    fig, axes = plt.subplots(1, len(arc_degs), figsize=(IEEE_DBL, 2.8),
                             squeeze=False, sharey=True)

    for ai, arc_deg in enumerate(arc_degs):
        ax = axes[0, ai]

        for n_cam in cam_counts:
            for method in ["Laplace", "SP"]:
                xs = []
                ys_rot = []
                ys_trans = []
                for pert_scale in pert_scales:
                    errs_rot = []
                    errs_trans = []
                    for seed in range(10):
                        key = (n_cam, arc_deg, pert_scale, method, seed)
                        if key not in traces:
                            continue
                        tr = traces[key]
                        errs_rot.append(tr[-1]["final_rot_err"])
                        errs_trans.append(tr[-1]["final_trans_err"])
                    if errs_rot:
                        xs.append(pert_scale)
                        ys_rot.append(np.mean(errs_rot))
                        ys_trans.append(np.mean(errs_trans))

                if xs:
                    color = METHOD_COLORS[method]
                    ls = METHOD_STYLES[method]
                    marker = CAM_MARKERS.get(n_cam, "x")
                    ax.plot(xs, ys_trans, color=color, linestyle=ls,
                            marker=marker, markersize=4, linewidth=1.2,
                            alpha=0.8, label=f"{n_cam}c {method}")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"Perturbation scale $\sigma_\mathrm{pert}$ (rad)")
        ax.set_title(rf"Arc $= {arc_deg:.0f}^\circ$")
        if ai == 0:
            ax.set_ylabel("Mean translation error (m)")
        ax.grid(True, alpha=0.3)

    # Deduplicated legend
    handles, labels = axes[0, -1].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(),
               loc="upper center", ncol=min(len(by_label), 4),
               bbox_to_anchor=(0.5, 1.05), frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    for ext in ("pgf", "png"):
        path = os.path.join(FIG_DIR, f"{out_name}.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


def fig_gradient_norm(traces, out_name="convergence_gnorm"):
    """
    Grid: rows = arc_deg, columns = n_cam.
    Each panel shows gradient norm vs iteration.
    """
    os.makedirs(FIG_DIR, exist_ok=True)

    cam_counts  = sorted({k[0] for k in traces})
    arc_degs    = sorted({k[1] for k in traces})
    pert_scales = sorted({k[2] for k in traces})

    n_rows = len(arc_degs)
    n_cols = len(cam_counts)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(IEEE_DBL, 2.2 * n_rows),
                             squeeze=False, sharex=True)

    cmap = plt.cm.viridis
    pert_colors = {p: cmap(i / max(len(pert_scales) - 1, 1))
                   for i, p in enumerate(pert_scales)}

    for ri, arc_deg in enumerate(arc_degs):
        for ci, n_cam in enumerate(cam_counts):
            ax = axes[ri, ci]

            for pert_scale in pert_scales:
                for method in ["Laplace", "SP"]:
                    for seed in range(10):
                        key = (n_cam, arc_deg, pert_scale, method, seed)
                        if key not in traces:
                            continue
                        tr = traces[key]
                        iters = [r["iter"] for r in tr]
                        gnorms = [r["gnorm"] for r in tr]
                        ax.plot(iters, gnorms,
                                color=pert_colors[pert_scale],
                                linestyle=METHOD_STYLES[method],
                                alpha=0.4, linewidth=0.8)

            ax.set_yscale("log")
            if ri == n_rows - 1:
                ax.set_xlabel("Iteration")
            if ci == 0:
                ax.set_ylabel(rf"$\|\nabla\|$ ($\mathrm{{arc}}={arc_deg:.0f}^\circ$)")
            if ri == 0:
                ax.set_title(rf"{n_cam} cam{'s' if n_cam > 1 else ''}")

    fig.tight_layout()
    for ext in ("pgf", "png"):
        path = os.path.join(FIG_DIR, f"{out_name}.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Convergence basin plots for pose inference sweep.")
    parser.add_argument("csv", nargs="?", default=DEFAULT_CSV,
                        help="Path to convergence_sweep.csv")
    args = parser.parse_args()

    csv_path = args.csv
    if not os.path.isfile(csv_path):
        print(f"Error: CSV not found at {csv_path}")
        print("Run the Rust sweep first (from repo root):")
        print("  cd rust && cargo run --release --example pose_inference -- --sweep && cd ..")
        print(f"Then ensure CSV is at: {DEFAULT_CSV}")
        sys.exit(1)

    rows = load_csv(csv_path)
    traces = group_traces(rows)
    n_traces = len(traces)
    n_rows = len(rows)
    print(f"Loaded {n_rows} rows, {n_traces} traces from {csv_path}")

    # Summary
    cam_counts = sorted({k[0] for k in traces})
    arc_degs = sorted({k[1] for k in traces})
    pert_scales = sorted({k[2] for k in traces})
    print(f"  Cameras: {cam_counts}")
    print(f"  Arc degs: {arc_degs}")
    print(f"  Perturbation scales: {pert_scales}")

    fig_convergence_traces(traces)
    fig_final_error(traces)
    fig_gradient_norm(traces)


if __name__ == "__main__":
    main()
