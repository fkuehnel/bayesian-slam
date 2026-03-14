# Reproducing and Creating Figures

This document describes the plotting workflow for the paper figures.
All plots use the **Matplotlib → PGF** pipeline so that text is
rendered by LaTeX with the document's own fonts.

---

## 0. What's in This Folder

```
experiments/
├── README.md                      # this file
├── plot_bias_experiment.py        # Experiment 1: bias scatter
├── plot_propagation.py            # Experiment 2: 2nd-order accuracy (planned)
├── plot_convergence.py            # Experiment 3: margMAP vs BA (planned)
├── plot_saddlepoint.py            # Experiments 3–4: SP correction (planned)
└── data/                          # Rust-generated CSVs land here
```

**Data flow convention:**

- Rust examples write CSVs into `experiments/data/`.
- Plot scripts here read from `data/` and write PGF + PNG into `paper/figures/`.
- The manuscript `\input`'s PGF files from `figures/` (relative to `paper/`).

> Create the data directory if it doesn't exist yet:
> ```bash
> mkdir -p experiments/data
> ```

---

## 1. Why PGF

IEEE journals require vector figures with matching typography.
The PGF backend emits LaTeX drawing commands instead of rasterized
text, so `\input{figures/fig.pgf}` inherits whatever font setup
the enclosing document has. No font‑embedding mismatches, no
post‑hoc label replacement.

PNG previews are saved alongside for quick inspection without
a full LaTeX build.

## 2. Shared Style Conventions

Every plotting script must set the Matplotlib PGF preamble
**before creating any figures**. Copy this block verbatim:

```python
import matplotlib
matplotlib.use("pgf")          # must precede pyplot import
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",   # or "lualatex" if the paper moves
    "font.family":   "serif",
    "text.usetex":   True,
    "pgf.rcfonts":   False,        # let LaTeX handle fonts entirely
    "axes.labelsize":  10,
    "font.size":       10,
    "legend.fontsize":  8,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
})
```

> **Key point:** `pgf.rcfonts: False` tells Matplotlib not to inject
> its own font commands into the `.pgf` file. The enclosing LaTeX
> document controls fonts.

### Figure sizes

IEEE single‑column width is ~3.5 in, double‑column ~7.16 in.
Use these as your `figsize` width; height follows from the
aspect ratio of the content.

| Layout | figsize width | Typical use |
|--------|--------------|-------------|
| Single‑column | 3.5 in | One scatter pair |
| 1.5‑column | 5.3 in | Three subplots in a row |
| Double‑column | 7.16 in | 2×3 grid of subplots |

The existing bias experiment uses `figsize=(15, 9)` (cm‑ish units)
because it targets a full‑page supplementary figure. For in‑text
figures, prefer IEEE widths.

### Color palette

| Role | Color | Hex |
|------|-------|-----|
| L1 / first‑kind data | `steelblue` | #4682B4 |
| L2 / second‑kind data | `firebrick` | #B22222 |
| Mean marker | `black` | — |
| Reference lines | `gray` | — |
| Saddlepoint curve | `tab:green` | — |
| Laplace curve | `tab:blue` | — |
| Monte Carlo truth | `black` dashed | — |

Stick to this mapping across all experiments for visual consistency.

### Rasterization

For scatter plots with thousands of points, pass `rasterized=True`
to `ax.scatter()`. This embeds the point cloud as a bitmap inside
the otherwise‑vector PGF file, keeping file sizes manageable while
axes, labels, and annotations remain vector.

## 3. Reproducing Existing Figures

All commands are run from the repo root.

### Experiment 1 — Bias scatter (Figs. 1–2)

```bash
# Step 1: generate the CSV from Rust
cd rust
cargo run --release --example bias_experiment
cd ..
# CSV is now at experiments/data/bias_experiment.csv

# Step 2: plot
python experiments/plot_bias_experiment.py
# PGF + PNG are now in paper/figures/
```

Outputs in `paper/figures/`:
- `bias_scatter_translation.pgf` / `.png`
- `bias_scatter_rotation.pgf` / `.png`

Include in the paper (`paper/SE3_inference_paper.tex`):
```latex
\begin{figure}[t]
    \centering
    \input{figures/bias_scatter_translation.pgf}
    \caption{Translation error scatter…}
    \label{fig:bias_translation}
\end{figure}
```

> **Note:** PGF files are LaTeX source (TikZ commands), so use
> `\input{}`, not `\includegraphics{}`. Make sure `\usepackage{pgf}`
> is in the preamble.

## 4. Adding a New Experiment Plot

Use this template. It follows the same structure as
`plot_bias_experiment.py` and plugs into the shared style.

```python
#!/usr/bin/env python3
"""
<One-line description>.

Reads CSV from:  experiments/data/<n>.csv
Outputs to:      paper/figures/<n>.{pgf,png}

Usage (from repo root):
    python experiments/plot_<n>.py [path/to/data.csv]

Default CSV: experiments/data/<n>.csv
"""

import sys, os
import numpy as np

import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt

# ── Shared style (copy from §2 above) ──────────────────────────
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.family":   "serif",
    "text.usetex":   True,
    "pgf.rcfonts":   False,
    "axes.labelsize":  10,
    "font.size":       10,
    "legend.fontsize":  8,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
})

IEEE_COL = 3.5   # inches, single-column
IEEE_DBL = 7.16  # inches, double-column

# ── Paths (repo-relative, work from any working directory) ──────
REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(REPO_ROOT, "experiments", "data")
FIG_DIR     = os.path.join(REPO_ROOT, "paper", "figures")
DEFAULT_CSV = os.path.join(DATA_DIR, "<n>.csv")


def load_data(path):
    """Load experiment CSV. Adapt column indices to your format."""
    raw = np.loadtxt(path, delimiter=",", skiprows=1)
    return {
        "x": raw[:, 0],
        "y": raw[:, 1],
        # ...
    }


def plot_figure(d):
    """Create and save the figure to paper/figures/."""
    os.makedirs(FIG_DIR, exist_ok=True)

    fig, ax = plt.subplots(figsize=(IEEE_COL, 2.5))

    ax.scatter(d["x"], d["y"], s=4, alpha=0.3,
               color="steelblue", rasterized=True)

    ax.set_xlabel(r"$x$-axis label with \LaTeX")
    ax.set_ylabel(r"$y$-axis label")
    ax.set_title(r"Experiment $N$: descriptive title")

    fig.tight_layout()
    for ext in ("pgf", "png"):
        p = os.path.join(FIG_DIR, f"plot_name.{ext}")
        fig.savefig(p, dpi=150)
        print(f"Saved {p}")
    plt.close(fig)


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV
    if not os.path.isfile(csv_path):
        print(f"CSV not found: {csv_path}")
        print("Run from repo root:")
        print("  cd rust && cargo run --release --example <n> && cd ..")
        sys.exit(1)

    d = load_data(csv_path)
    print(f"Loaded {len(d['x'])} rows from {csv_path}")
    plot_figure(d)


if __name__ == "__main__":
    main()
```

### Checklist for a new plot script

1. **Rust example emits CSV to `data/`.** Comma‑separated,
   header row, one sample per line. Keep columns documented in the
   Rust source.

2. **Python script reads from `data/`, writes to
   `paper/figures/`.** Follow the template above — use `REPO_ROOT`,
   `DATA_DIR`, `FIG_DIR` so the script works from any working directory.

3. **Use shared rcParams verbatim.** Don't tweak font sizes per
   script; consistency across figures matters for the journal.

4. **Use the color palette from §2.** If you need a new role,
   add it to the table and use it consistently.

5. **Rasterize dense point clouds** (`rasterized=True`).
   Keep markers, lines, and text as vector.

6. **Match axis limits across paired subplots** (L1 vs L2,
   Laplace vs saddlepoint) so visual comparison is fair.
   See `scatter_pair()` in `plot_bias_experiment.py` for how
   this is done with `np.percentile`.

7. **LaTeX in labels.** Use raw strings: `r"$\delta t_1$"`.
   These will be compiled by the document's own LaTeX engine.

8. **Save both formats to `paper/figures/`.** PGF for the paper,
   PNG for quick review and the README.

## 5. Planned Figures

| Figure | Plot script | Output in `paper/figures/` | Status |
|--------|-------------|---------------------------|--------|
| Bias scatter (translation) | `plot_bias_experiment.py` | `bias_scatter_translation.{pgf,png}` | ✅ Done |
| Bias scatter (rotation) | `plot_bias_experiment.py` | `bias_scatter_rotation.{pgf,png}` | ✅ Done |
| 2nd-order covariance accuracy | `plot_propagation.py` | `propagation_accuracy.{pgf,png}` | Planned |
| Convergence: margMAP vs BA | `plot_convergence.py` | `convergence_comparison.{pgf,png}` | Planned |
| Saddlepoint c₁ vs depth | `plot_saddlepoint.py` | `sp_c1_vs_depth.{pgf,png}` | Planned |
| c₁ vs number of cameras | `plot_saddlepoint.py` | `sp_c1_vs_cameras.{pgf,png}` | Planned |
| SP validity regime sweep | `plot_sp_regime.py` | `sp_regime_sweep.{pgf,png}` | Planned |

## 6. LaTeX Integration

### Preamble additions (in `paper/SE3_inference_paper.tex`)

```latex
\usepackage{pgf}          % required for .pgf inclusion
% graphicx is already loaded by IEEEtran
```

### Including a PGF figure

Paths are relative to `paper/` (where `SE3_inference_paper.tex` lives):

```latex
\begin{figure}[t]
    \centering
    \input{figures/bias_scatter_translation.pgf}   % \input, not \includegraphics
    \caption{…}
    \label{fig:…}
\end{figure}
```

> **`\input` vs `\includegraphics`:** PGF files are LaTeX source
> code (TikZ commands), so they must be `\input`'d. If the file
> is large, wrap it with `\resizebox{\columnwidth}{!}{\input{…}}`
> to scale. Alternatively, you can use `\includegraphics` if the
> PGF file contains a `\beginpgfgraphicnamed` wrapper — Matplotlib
> does not emit this by default, so prefer `\input`.

### Coexistence with legacy EPS figures

The current paper uses `\getname{…}` to include EPS files from
`paper/eps/`. New PGF figures go in `paper/figures/` and use
`\input{figures/…}`. Both can coexist during the transition —
replace EPS includes one at a time as the new plots are generated.

### Scaling

If the figure was created at the correct `figsize`, no scaling
is needed. If you must rescale:

```latex
\resizebox{\columnwidth}{!}{\input{figures/plot.pgf}}
```

This preserves vector quality but will scale text — which is
why getting `figsize` right at creation time is preferable.

## 7. Troubleshooting

**"Dimension too large" errors during LaTeX compilation.**
PGF scatter plots with >5k points can exceed TeX's memory.
Solutions: (a) `rasterized=True` on the scatter call,
(b) increase TeX memory via `--extra-mem-top=10000000`,
(c) reduce point count by subsampling before plotting.

**Fonts don't match the document.**
Verify `pgf.rcfonts: False` is set and you're using `\input`
not `\includegraphics`. The PGF file should contain no
`\selectfont` commands.

**Slow compilation with many PGF figures.**
Use the TikZ `external` library to cache compiled figures:
```latex
\usetikzlibrary{external}
\tikzexternalize[prefix=tikz-cache/]
```
Each figure compiles once, then reuses the cached PDF.

**Missing LaTeX packages in PGF file.**
Matplotlib's PGF output may require packages not loaded by
IEEEtran. If compilation fails with undefined control sequences,
add `pgf.preamble` to rcParams:
```python
matplotlib.rcParams["pgf.preamble"] = r"\usepackage{amsmath}"
```
