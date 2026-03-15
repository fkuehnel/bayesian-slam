# Experiments, Examples, and Figure Reproduction

This document describes how to run all experiments, reproduce paper figures,
and summarizes key numerical results.

**Data flow:**
`rust/examples/*.rs` → `experiments/data/*.csv` → `experiments/plot_*.py` → `paper/figures/*.{pgf,png}` → `\input{figures/…}` in `paper/SE3_inference_paper.tex`

---

## Getting Started

```bash
cd rust
cargo build --release
cargo test
cargo test -- --nocapture  # see diagnostic output
```

> **Cache trouble?** If a build fails unexpectedly:
> ```bash
> cargo clean && cargo build --release
> ```

---

## Experiment 1: Coordinate Bias (`bias_experiment`)

Demonstrates that Lie-Cartan exponential coordinates (first kind) produce unbiased pose estimates while second-kind (additive) coordinates [Ω, T] show systematic translational bias. This validates Proposition 1 in the paper (§II.3).

**Setup:** A camera at a known pose observes N landmarks with Gaussian noise. The pose is estimated by point cloud alignment (Gauss-Newton on SE(3)). Repeated over 2000 noise realizations, the empirical mean is compared in both coordinate systems.

**Paper config** (N=3, σ=1.0): produces a dramatic 21× bias ratio, demonstrating the effect with few, noisy landmarks.

**Mild config** (N=8, σ=0.5): ~3× ratio, showing the effect persists with more landmarks and lower noise.

```bash
cargo run --release --example bias_experiment          # paper config: N=3, σ=1.0
cargo run --release --example bias_experiment mild     # mild config: N=8, σ=0.5
```

**Output:**
- Console: mean and std in L1 (exponential) and L2 (additive) coordinates, bias ratio
- CSV: `experiments/data/bias_experiment.csv` with per-sample data for scatter plots (paper Figs. 1–2)

**Plotting** (from repo root):

```bash
python experiments/plot_bias_experiment.py
# → paper/figures/bias_scatter_translation.{pgf,png}
# → paper/figures/bias_scatter_rotation.{pgf,png}
```

---

## Experiment 2: Pose Inference (`pose_inference`)

The core use case: estimate camera pose(s) by maximizing the marginalized posterior over landmark positions, comparing Laplace vs saddlepoint-corrected objectives. Supports single-camera and multi-camera configurations.

```bash
cargo run --release --example pose_inference [N_cameras] [N_landmarks] [overlap]
```

**Examples:**
```bash
cargo run --release --example pose_inference                  # 1 camera, 12 landmarks
cargo run --release --example pose_inference 4                # 4 cameras, 12 landmarks
cargo run --release --example pose_inference 4 30 0.5         # 4 cameras, 30 landmarks, 50% overlap
cargo run --release --example pose_inference 2 20 0.7         # stereo, 20 landmarks, 70% overlap
```

**What it does:**
- **Part 1 — Optimization:** Damped Newton on SE(3) with compositive updates f ← f·exp(δξ), run separately for Laplace and saddlepoint objectives. Multi-camera uses Jacobi-style simultaneous updates.
- **Part 2 — Results:** Per-camera pose errors vs ground truth, Laplace–SP differences.
- **Part 3 — Per-landmark detail:** c₁, A/12, −Q₄/8 at the Laplace optimum, confirming closer landmarks drive the correction.

**Parameters:**
| Argument | Default | Description |
|----------|---------|-------------|
| `N_cameras` | 1 | Number of cameras (1 = single, 2+ = ring) |
| `N_landmarks` | 12 | Number of 3D landmarks |
| `overlap` | 1.0 | Fraction of landmarks visible to each camera (0.1–1.0) |

---

## Experiment 3: Multi-Camera Saddlepoint (`multicam_saddlepoint`)

Saddlepoint-corrected landmark marginalization with 2–12 cameras observing a shared 3D point cloud in a ring configuration. Demonstrates how the correction magnitude decreases as more cameras constrain each landmark.

**Setup:** Cameras arranged in a ring at radius 8m, all looking at the origin. 20 landmarks in [−2, 2]³. Each landmark is optimized and marginalized using the multi-camera API.

```bash
cargo run --release --example multicam_saddlepoint             # default: 4 cameras, 20 landmarks
cargo run --release --example multicam_saddlepoint 2 30        # stereo, 30 landmarks
cargo run --release --example multicam_saddlepoint 12 20       # 12 cameras, 20 landmarks
```

**Output:**
- Per-landmark table: views, depth, c₁, A/12, B/8, −Q₄/8
- Camera-count sweep (2–12 cameras) for a test point at the origin
- CSV: `experiments/data/multicam_saddlepoint.csv`

**Key finding:** The stereo case (2 cameras) produces corrections ~200× larger than 4-camera, confirming the correction matters most in the typical multi-view operating regime.

---

## Experiment 4: Extended Multi-Camera Analysis (`multicam_experiment`)

More detailed multi-camera experiment with Monte Carlo integration as ground truth for the marginal integral. Validates the saddlepoint correction against numerical integration.

**Setup:** Similar to Experiment 3 but adds importance-sampling MC integration (200k samples) for 1–2 camera configurations, plus a depth-dependence sweep.

```bash
cargo run --release --example multicam_experiment
```

**Output:**
- Part 1: c₁ vs number of cameras (averaged over point cloud)
- Part 2: Laplace vs saddlepoint vs MC comparison (1 and 2 cameras)
- Part 3: c₁ vs depth at fixed 2-camera baseline
- CSV: `experiments/data/multicam_experiment.csv`

---

## Key Numerical Results

### Multi-camera saddlepoint: camera-count sweep

| N_cam | views | σ_depth | σ/depth | c₁ |
|------:|------:|--------:|--------:|---:|
| 2 | 2 | 0.028 | 0.0037 | 7.5×10⁻³ |
| 3 | 3 | 0.031 | 0.0041 | −1.2×10⁻⁶ |
| 4 | 4 | 0.028 | 0.0037 | −3.5×10⁻⁵ |
| 6 | 6 | 0.023 | 0.0031 | −2.4×10⁻⁵ |
| 8 | 8 | 0.020 | 0.0026 | −1.8×10⁻⁵ |
| 12 | 12 | 0.016 | 0.0022 | −1.2×10⁻⁵ |

### Stereo full point cloud (2 cameras, 30 landmarks)

- 100% saddlepoint validity across all landmarks
- Individual corrections up to c₁ ≈ 9×10⁻³
- Cumulative correction: Σ|Δ| = 4.4×10⁻²
- Correction dominated by negative Q₄ term (depth non-Gaussianity)

### Analytical Q₄: convergence impact

The original FD-based Q₄ created nested finite-difference interference, preventing the saddlepoint optimizer from converging. The closed-form Q₄ eliminates this:

| Method | Iterations | Final |grad| | Step sizes |
|--------|-----------|----------------|------------|
| Laplace (before) | 5 | 1e-10 | h_grad=1e-5, h_hess=1e-4 |
| Saddlepoint (before, FD Q₄) | 60+ | ~2.5 | h_grad=1e-3, h_hess=1e-3 |
| **Saddlepoint (analytical Q₄)** | **5** | **5e-10** | **h_grad=1e-5, h_hess=1e-4** |

The analytical formula at the mode (residual e ≈ 0):

```
f''''_{abcd} = Σ_{mn} Σ⁻¹_{mn} × [
  P_{ma}·D3ⁿ_{bcd} + P_{mb}·D3ⁿ_{acd} + P_{mc}·D3ⁿ_{abd} + P_{md}·D3ⁿ_{abc}
  + H^m_{ab}·Hⁿ_{cd} + H^m_{ac}·Hⁿ_{bd} + H^m_{ad}·Hⁿ_{bc}
]
```

No 4th derivatives of π needed (D4·e vanishes at the mode). The prior contributes zero.

### The depth soft mode

With convergence fixed, single-camera pose inference reveals a fundamental limitation:

```
                           True      Laplace  Saddlepoint
  t₃ (depth)          0.000000     6.451333     6.366716

  Translation error:
    Laplace:      6.454793 m
    Saddlepoint:  6.370215 m
    Difference:   0.084597 m
```

The SP correction shifts the estimate by only ~8.5cm on a ~6.4m error. **This is an information-theoretic limitation, not an approximation artifact.** The projection π(x) = [x/z, y/z] is nearly invariant under depth rescaling, making the posterior flat along the forward translation direction. No per-landmark correction can create absent depth information.

Resolving the depth soft mode requires:
- **Multi-camera / stereo baselines** — parallax provides direct depth triangulation
- **Strong pose priors** — constrain the translation subspace directly
- **Known-scale constraints** — e.g., known inter-landmark distances or IMU-derived scale

Multi-camera results confirm this:
- **1 camera**: trans err ~6.4m, c₁ ~ 2×10⁻²
- **4 cameras**: trans err ~0.5m, c₁ ~ 2×10⁻⁴ (100× smaller)

---

## Reproducing Paper Figures

All paper figures use the **Matplotlib → PGF** pipeline: Rust examples
emit CSV data to `experiments/data/`, Python scripts produce `.pgf`
files in `paper/figures/`, which are `\input`'d into the manuscript.

### Quick reference

| Figure | Rust example | Plot script | Output in `paper/figures/` |
|--------|-------------|-------------|---------------------------|
| Bias scatter (translation) | `bias_experiment.rs` | `plot_bias_experiment.py` | `bias_scatter_translation.{pgf,png}` |
| Bias scatter (rotation) | `bias_experiment.rs` | `plot_bias_experiment.py` | `bias_scatter_rotation.{pgf,png}` |
| 2nd-order covariance accuracy | `propagation.rs` (MC) | `plot_propagation.py` | Planned |
| Convergence: margMAP vs BA | `pose_inference.rs` | `plot_convergence.py` | Planned |
| Saddlepoint c₁ vs depth | `multicam_saddlepoint.rs` | `plot_saddlepoint.py` | Planned |
| c₁ vs number of cameras | `multicam_saddlepoint.rs` | `plot_saddlepoint.py` | Planned |
| SP validity regime sweep | `multicam_experiment.rs` | `plot_sp_regime.py` | Planned |

### End-to-end example (from repo root)

```bash
# 1. Create output directories (once)
mkdir -p experiments/data paper/figures

# 2. Generate data
cd rust
cargo run --release --example bias_experiment
cd ..

# 3. Plot
python experiments/plot_bias_experiment.py

# 4. Include in paper (note: \input, not \includegraphics — PGF is LaTeX source)
#    \begin{figure}[t]
#        \centering
#        \input{figures/bias_scatter_translation.pgf}
#        \caption{…}
#    \end{figure}
```

---

## Plotting Guide

### Why PGF

IEEE journals require vector figures with matching typography.
The PGF backend emits LaTeX drawing commands instead of rasterized
text, so `\input{figures/fig.pgf}` inherits whatever font setup
the enclosing document has. No font-embedding mismatches, no
post-hoc label replacement.

PNG previews are saved alongside for quick inspection without
a full LaTeX build.

### Shared Style Conventions

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

### Figure sizes

| Layout | figsize width | Typical use |
|--------|--------------|-------------|
| Single-column | 3.5 in | One scatter pair |
| 1.5-column | 5.3 in | Three subplots in a row |
| Double-column | 7.16 in | 2×3 grid of subplots |

### Color palette

| Role | Color | Hex |
|------|-------|-----|
| L1 / first-kind data | `steelblue` | #4682B4 |
| L2 / second-kind data | `firebrick` | #B22222 |
| Mean marker | `black` | --- |
| Reference lines | `gray` | --- |
| Saddlepoint curve | `tab:green` | --- |
| Laplace curve | `tab:blue` | --- |
| Monte Carlo truth | `black` dashed | --- |

### Adding a new plot script

Use the template in the codebase (`plot_bias_experiment.py`) as a starting point. Checklist:

1. Rust example emits CSV to `data/` (comma-separated, header row)
2. Python script reads from `data/`, writes to `paper/figures/`
3. Use shared rcParams verbatim
4. Use the color palette above consistently
5. Rasterize dense point clouds (`rasterized=True`)
6. Match axis limits across paired subplots
7. Use raw strings for LaTeX in labels: `r"$\delta t_1$"`
8. Save both PGF and PNG to `paper/figures/`

### LaTeX Integration

```latex
% In preamble:
\usepackage{pgf}

% In document:
\begin{figure}[t]
    \centering
    \input{figures/bias_scatter_translation.pgf}   % \input, not \includegraphics
    \caption{…}
\end{figure}
```

### Troubleshooting

- **"Dimension too large" errors:** Use `rasterized=True` on scatter plots, or increase TeX memory
- **Fonts don't match:** Verify `pgf.rcfonts: False` and use `\input` not `\includegraphics`
- **Slow compilation:** Use TikZ `external` library to cache compiled figures
