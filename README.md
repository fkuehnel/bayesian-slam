<p align="center">
  <img src="assets/hero.png" alt="SU(2), SO(3), and SE(3) geometry with uncertainty propagation" width="1000">
</p>

# Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on SE(3)

[![Rust Tests](https://github.com/fkuehnel/bayesian-slam/actions/workflows/rust-tests.yml/badge.svg)](https://github.com/fkuehnel/bayesian-slam/actions/workflows/rust-tests.yml)

**Authors:** Frank O. Kuehnel, Andre Jalobeanu

---

## Overview

This repository accompanies the paper *Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on the SE(3) Lie Group*, which develops a self-contained framework for Bayesian inference over rigid body poses.

**Paper:** [SE3_inference_paper.pdf](paper/SE3_inference_paper.pdf)

The core ideas originate from work at NASA Ames Research Center (2006) on [robust pose estimation using the SE(3) Lie group](paper/robustEst.pdf) structure. This project refocuses that foundational material into contributions aimed at the broader estimation and inference community. All formulas have been verified symbolically (Mathematica) and numerically (Rust finite differences, Python, Monte Carlo).

## What This Paper Contributes

### 1. Second-Order Uncertainty Propagation on SE(3)

Current filters and optimizers (EKF, IEKF, factor graphs) propagate pose uncertainty using first-order Jacobians. On SE(3), the semi-direct product coupling between rotation and translation introduces systematic bias that first-order methods miss. We derive:

- Composition Jacobians including the novel compact closed-form rotation-translation coupling Jacobian **J_t(Ω, t)**
- Second-order corrections to both the mean and covariance of composed poses
- Explicit formulas using the finite BCH composition via SU(2), avoiding truncated series

**Verified:** Monte Carlo validation (10⁵ samples) confirms a 1.83× improvement in covariance accuracy at σ_ω ≈ 0.2 rad. BCH truncation errors verified across 4 decades (see Table 1 in paper).

### 2. Saddlepoint Marginalization for Projective Observations

When camera observations are projective (x/z, y/z), the joint posterior over pose and 3D landmarks is non-Gaussian in the pose parameters. Standard Laplace approximation (Schur complement in bundle adjustment) misses this asymmetry. We develop:

- A saddlepoint correction with the formula c₁ = (1/12)A + (1/8)B − (1/8)Q₄
- Two distinct cubic contraction types (cross: 6/15 Isserlis pairings, trace: 9/15)
- Quartic term with correct negative sign (from exp(−g₄) ≈ 1 − g₄)
- **Closed-form Q₄**: analytical quartic contraction using only existing projective derivatives (P, H, D3) — no 4th derivatives of π needed (the D4·e term vanishes at the mode)
- Validity guard: correction applied only when |c₁| < 0.5 (σ_depth/depth ≲ 0.3)

**Verified:** Against numerical quadrature, the saddlepoint achieves 6 significant figures (3.2×10⁻⁶ relative error) compared to Laplace's 0.94% error — a ~3000× improvement.

### 3. Complete SE(3) Algebra Toolkit

The paper provides a self-contained reference for the SE(3) machinery, all symbolically verified. This machinery is equivalent to the standard toolkit in **Solà, et. al**, and extends it: 

- Exponential/logarithmic maps with Lie-Cartan coordinates of the first kind (already known)
- Closed-form finite Rodrigues vector composition via the SU(2)/Z₂ ≅ SO(3) double cover
- Phase reflection handling at Θ = π with cutline analysis
- Wei-Norman formula and both left/right Rodrigues Jacobians
- Key identities: S⁻¹ = J_ωl, S⁻¹R = RS⁻¹ = J_ωr (verified to machine precision)
- det S = 2(1−cosΘ)/Θ² (strictly positive, guarantees invertibility)

## Repository Structure

```
RobustPoseEst/
├── README.md                          # This file
├── paper/
│   ├── SE3_inference_paper.tex        # Main manuscript (LaTeX, ~2200 lines)
│   ├── eps/                           # Legacy EPS figures (original plots)
│   ├── figures/                       # PGF + PNG outputs from plotting scripts
│   └── robustEst.tex                  # Original 2006 technical report (reference)
├── verification/
│   ├── README.md                      # Detailed verification results (per-part tables)
│   ├── SE3AlgebraVerification.m       # Mathematica: 24 verified identities (SO3/SE3/BCH/Jacobians)
│   ├── SaddlepointVerification.m      # Mathematica: projective derivs, saddlepoint formula, quadrature
│   └── CouplingJacobianDerivation.m   # Mathematica: symbolic proof of d[S⁻¹T]/dΩ (9/9 Omega subs)
├── rust/
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs                     # Array-based linalg: Vec3/Mat3/Vec6/Mat6
│   │   ├── so3.rs                     # SO(3): Rodrigues exp/log, S matrix, S⁻¹
│   │   ├── se3.rs                     # SE(3): Pose struct, compose/inverse/act/exp/log
│   │   ├── bch.rs                     # Finite BCH via SU(2) quaternions, phase reflection, Jacobians
│   │   ├── jacobians.rs               # J_ωr, J_ωl, J_t coupling (analytic T-form), 6×6 SE(3) Jacobian
│   │   ├── projective.rs              # Pinhole camera, derivatives through 3rd order, third cumulants, analytical Q₄
│   │   ├── saddlepoint.rs             # Landmark optimization, saddlepoint correction (analytical Q₄), validity guard
│   │   └── propagation.rs             # First/second-order covariance transport, MC validation
│   └── examples/
│       ├── bias_experiment.rs         # Experiment 1: L1 vs L2 coordinate bias (21× ratio)
│       ├── pose_inference.rs          # Experiment 2: pose estimation with SP-corrected marginal
│       ├── multicam_saddlepoint.rs    # Experiment 3: multi-camera saddlepoint with camera sweep
│       └── multicam_experiment.rs     # Experiment 4: extended multi-camera analysis with MC truth
└── experiments/                       # ing scripts for paper figures (see ING.md)
    ├── README.md                      # Plotting guide: style, templates, LaTeX integration
    ├── plot_bias_experiment.py        # Experiment 1: L1 vs L2 bias scatter (PGF + PNG)
    └── data/                          # Rust-generated CSVs (created by examples)
```

**Data flow:**
`rust/examples/*.rs` → `experiments/data/*.csv` → `experiments/plot_*.py` → `paper/figures/*.{pgf,png}` → `\input{figures/…}` in `paper/SE3_inference_paper.tex`

See [experiments/README.md](experiments/README.md) for the full plotting guide.

## Rust Implementation

**95 tests passing, 0 failures, 0 warnings.**

| Module | Tests | Status | Description |
|--------|-------|--------|-------------|
| `so3` | 8 | ✅ | Rodrigues exp/log, hat/vee, S matrix, S⁻¹ |
| `se3` | 8 | ✅ | Pose compose/inverse/act, exp/log roundtrip |
| `bch` | 30 | ✅ | Finite BCH via SU(2), phase reflection, composition Jacobians |
| `jacobians` | 13 | ✅ | J_ωr/J_ωl, analytic J_t (T-form), 6×6 SE(3) Jacobian (FD-verified) |
| `projective` | 14 | ✅ | Project, Jacobian, Hessian, 3rd/4th derivs, third cumulants, analytical Q₄ (all FD-verified) |
| `saddlepoint` | 10 | ✅ | Landmark GN optimizer, corrected c₁ formula, validity guard, analytical Q₄, multi-cam |
| `propagation` | 10 | ✅ | First/second-order covariance transport, Levi-Civita, Isserlis correction, MC validation |

### Key verified identities (Rust + Mathematica + Python)

- BCH composition matches matrix log(R_a R_b) to 10⁻¹⁶ across all angle regimes
- Phase reflection at Θ > π produces correct (−r, 2π−Θ) identification
- Composition Jacobians ∂Ω_c/∂Ω_a, ∂Ω_c/∂Ω_b match FD to 10⁻⁸
- S⁻¹R = J_ωr verified symbolically and numerically
- Full 6×6 SE(3) Jacobian block structure [[J_ωr, 0], [J_t, J_ωr]] verified to 10⁻⁹
- **Coupling Jacobian J_t**: analytic T-form verified to 7.3×10⁻¹⁰ against FD; symbolic equivalence proven in Mathematica (9/9 substitutions)
- Projective derivatives through 4th order verified against FD
- Saddlepoint correction matches numerical quadrature to 6 significant figures
- Analytical Q₄ (quartic contraction) verified against FD of exact Hessian

### Design Principles

- **No external dependencies for core algebra**: 3×3 and 6×6 matrices as fixed-size arrays, exploiting skew-symmetric and block-triangular structure
- **Zero heap allocation in inner loops**: SO(3) Jacobians are `[[f64; 3]; 3]`, SE(3) ones are `[[f64; 6]; 6]`
- **Defensive numerics**: Validity guard on saddlepoint (|c₁| > 0.5 → fall back to Laplace), Taylor expansions near singularities (Θ → 0, ε → 0)
- **Exhaustive FD validation**: Every Jacobian and derivative tested against central finite differences


### Getting Started

```bash
git clone https://github.com/fkuehnel/bayesian-slam.git
cd bayesian-slam/rust

cargo build --release
cargo test
cargo test -- --nocapture  # see diagnostic output
```

> **Cache trouble?** Rust's incremental compilation cache can occasionally
> cause stale builds or cryptic errors after switching branches, updating the
> toolchain, or interrupted compilations. If a build fails unexpectedly,
> `cargo clean` is the nuclear option — it wipes the entire `target/`
> directory and forces a full rebuild from scratch:
> ```bash
> cargo clean && cargo build --release
> ```

## Examples

### Experiment 1: Coordinate Bias (`bias_experiment`)

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

See [experiments/README.md](experiments/README.md) for style conventions and the template for new plot scripts.

### Experiment 2: Pose Inference (`pose_inference`)

The core use case: estimate camera pose by maximizing the marginalized posterior over landmark positions, comparing Laplace vs saddlepoint-corrected objectives. This is the full pipeline connecting all three paper contributions.

**Setup:** A camera observes 12 landmarks at varying depths with uncertain 3D priors. The negative log marginalized posterior is optimized over SE(3) using damped Newton with compositive updates f ← f·exp(δξ).

**Part 1 — 1D sweep:** Evaluates both objectives at 21 points along the depth (t_z) direction, revealing where each peaks. The saddlepoint minimum shifts relative to Laplace because close landmarks have depth-dependent non-Gaussian corrections.

**Part 2 — Full 6D optimization:** Runs Newton from a perturbed initial guess, once for Laplace and once for saddlepoint. Both converge, but to slightly different poses.

**Part 3 — Comparison:** Prints converged poses alongside ground truth, pose error metrics, and per-landmark SP correction magnitudes.

**Part 4 — Per-landmark detail:** Shows c₁ vs depth at the Laplace optimum, confirming that closer landmarks drive the correction.

```bash
cargo run --release --example pose_inference           # moderate range (depth≈8, σ_prior=2)
cargo run --release --example pose_inference close     # close range (depth≈4, σ_prior=3)
```

The `close` regime is where the Laplace–saddlepoint divergence is largest — the projective non-Gaussianity is strongest when σ_depth/depth is non-negligible.

### Experiment 3: Multi-Camera Saddlepoint (`multicam_saddlepoint`)

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

### Experiment 4: Extended Multi-Camera Analysis (`multicam_experiment`)

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

## Reproducing Paper Figures

All paper figures use the **Matplotlib → PGF** pipeline: Rust examples
emit CSV data to `experiments/data/`, Python scripts produce `.pgf`
files (LaTeX source for TikZ) in `paper/figures/`, which are `\input`'d
into the manuscript so that figure text is rendered with the document's
own fonts. PNG previews are saved alongside for quick inspection.

The full guide — shared style conventions, color palette, figure sizing
for IEEE column widths, a copy-paste template for new plot scripts, and
LaTeX integration instructions — lives in
**[experiments/README.md](experiments/README.md)**.

### Quick reference

| Figure | Rust example | Plot script | Output in `paper/figures/` |
|--------|-------------|-------------|---------------------------|
| Bias scatter (translation) | `rust/examples/bias_experiment.rs` | `experiments/plot_bias_experiment.py` | `bias_scatter_translation.{pgf,png}` |
| Bias scatter (rotation) | `rust/examples/bias_experiment.rs` | `experiments/plot_bias_experiment.py` | `bias_scatter_rotation.{pgf,png}` |
| 2nd-order covariance accuracy | `rust/src/propagation.rs` (MC) | `experiments/plot_propagation.py` | Planned |
| Convergence: margMAP vs BA | `rust/examples/pose_inference.rs` | `experiments/plot_convergence.py` | Planned |
| Saddlepoint c₁ vs depth | `rust/examples/multicam_saddlepoint.rs` | `experiments/plot_saddlepoint.py` | Planned |
| c₁ vs number of cameras | `rust/examples/multicam_saddlepoint.rs` | `experiments/plot_saddlepoint.py` | Planned |
| SP validity regime sweep | `rust/examples/multicam_experiment.rs` | `experiments/plot_sp_regime.py` | Planned |

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

To add a new figure for another experiment, follow the template and
checklist in [experiments/README.md §4](experiments/README.md).

### Pending Work

| Item | Priority | Description |
|------|----------|-------------|
| Remaining plot scripts | High | `experiments/plot_propagation.py`, `experiments/plot_convergence.py`, `experiments/plot_saddlepoint.py`, `experiments/plot_sp_regime.py` (see [experiments/README.md §5](experiments/README.md)) |
| Good Notations | High | Good consistent notations that the target audience knows |
| Experimental findings | High | Our paper must be based on repeatable experiments |
| More examples | Medium | An example script for pose and 3D landmark inference |
| Verify Solà notation | Medium | Mathematica scripts to confirm notation mappings |
| Depth soft mode analysis | Medium | Investigate multi-camera or prior-based approaches to resolve the translation degeneracy |

### Resolved Items

| Item | Resolution |
|------|------------|
| Plotting workflow | ✅ Documented in [experiments/README.md](experiments/README.md): shared style, PGF pipeline, templates, LaTeX integration |
| Pose inference example | ✅ `rust/examples/pose_inference.rs`: full 6D Newton with SP-corrected marginal, 1D sweep + comparison |
| Closed-form J_t | ✅ Derived T-form with α = (sinΘ−Θ)/(2(1−cosΘ)), verified to 7.3×10⁻¹⁰ |
| T-form symbolic proof | ✅ Proven in Mathematica: 9/9 Ω substitutions with symbolic T (`verification/CouplingJacobianDerivation.m`) |
| Erratum eq 78 | ✅ Second equality was wrong: RHS is β, not α. Corrected in paper |
| Saddlepoint correction formula | ✅ Corrected to c₁ = (1/12)A + (1/8)B − (1/8)Q₄, verified against quadrature |
| SE(3) Jacobian FD test | ✅ Was `#[ignore]`, now passing (bug was in original j_coupling, not the formula) |
| Closed-form Q₄ | ✅ Analytical quartic contraction using P, H, D3 — no 4th derivatives of π needed. Eliminates nested FD; SP optimizer converges in 5 iterations (was 60+). Verified against FD of exact Hessian |
| Pose inference panics | ✅ Fixed NaN/degenerate configurations, added pose prior for depth regularization |

## Notation Correspondence with Solà et al. / Barfoot

The paper includes Appendix F (§F) that maps our conventions to those of Solà, Deray & Atchuthan ("A micro Lie theory," 2018) and Barfoot (*State Estimation for Robotics*, 2024). This assists readers familiar with either reference.

### Block ordering

The fundamental difference is tangent vector ordering:

| | Solà/Barfoot (SDA) | This paper |
|---|---|---|
| Tangent vector | τ = [**ρ**; **θ**] (translation first) | τ = [**Ω**; **t**] (rotation first) |
| Primary Jacobian | Left: J_l | Right: J_r = J_l(−τ) |

All 6×6 matrices are related by the block permutation P = [[0, I₃]; [I₃, 0]], i.e. **M** = P **M**^SDA P.

### Key object correspondence

| Object | SDA / Barfoot | This paper | Relationship |
|---|---|---|---|
| S matrix | **V**(θ), Eq. 174 | **S**(Ω) | Identical formula |
| S⁻¹ | (not identified) | **S**⁻¹ = J_ωl | **Our identity** |
| Adjoint | [[R, [t]×R]; [0, R]] | [[R, 0]; [[T]×R, R]] | Block permutation |
| 6×6 left Jac | [[J_l, **Q**]; [0, J_l]] | [[J_ωl, 0]; [Q̃, J_ωl]] | Block permutation |
| 6×6 right Jac | J_l(−ρ,−θ) | [[J_ωr, 0]; [J_t, J_ωr]] | Block permutation |
| Coupling block | **Q**(ρ,θ) (Barfoot, 4-line expansion) | J_t (compact T-form) | J_t = Q(−t,−Ω)·J_ωr |
| Point action Jac | [R, −R[p]×] | [−R[x]×, I] | Block swap; body frame |

### Three additional identities (not in SDA/Barfoot)

1. **S⁻¹ = J_ωl = J_ωr(−Ω)** — Halves cost: evaluate J_ωl instead of inverting S numerically. Also implies S⁻¹R = RS⁻¹ = J_ωr.

2. **det S = 2(1 − cos Θ)/Θ²** — Exact volume element for Lie-Cartan coordinate chart. Needed for density transformations between [Ω, t] and [Ω, T] coordinates. Continuous at Θ = 0 (det S → 1), vanishes at Θ = 2π.

3. **Compact T-form coupling Jacobian** — Single-line formula for ∂[S⁻¹T]/∂Ω with projector decomposition into axial and transverse components, replacing Barfoot's four-line iterated cross-product expansion. Proven algebraically equivalent in Mathematica (verification/CouplingJacobianDerivation.m, 9/9 Ω substitutions).

## Multi-Camera Saddlepoint Results

### Camera-count sweep (test point at origin)

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

## Pose Inference: Analytical Q₄ and the Depth Soft Mode

### The nested finite-difference problem

The original saddlepoint implementation computed Q₄ via finite differences (FD) of the Hessian. This created a nested FD problem: the outer pose optimizer also uses FD for its gradient and Hessian. The inner FD at h=1e-5 forced the outer FD to use h=1e-3 to avoid interference, producing noisy gradients that prevented quadratic convergence:

| Method | Iterations | Final |grad| | Step sizes |
|--------|-----------|----------------|------------|
| Laplace (before) | 5 | 1e-10 | h_grad=1e-5, h_hess=1e-4 |
| Saddlepoint (before, FD Q₄) | 60+ | ~2.5 | h_grad=1e-3, h_hess=1e-3 |
| **Saddlepoint (analytical Q₄)** | **5** | **5e-10** | **h_grad=1e-5, h_hess=1e-4** |

### The analytical formula

At the mode (residual e ≈ 0), the 4th derivative of the NLL has a clean closed form:

```
f''''_{abcd} = Σ_{mn} Σ⁻¹_{mn} × [
  P_{ma}·D3ⁿ_{bcd} + P_{mb}·D3ⁿ_{acd} + P_{mc}·D3ⁿ_{abd} + P_{md}·D3ⁿ_{abc}
  + H^m_{ab}·Hⁿ_{cd} + H^m_{ac}·Hⁿ_{bd} + H^m_{ad}·Hⁿ_{bc}
]
```

where P = project_jacobian, H = project_hessian, D3 = project_third_deriv. The D4·e term vanishes at the mode (same reason the D3·e term vanishes for third cumulants). The prior contributes zero (quadratic → zero 4th derivative). The contraction Q₄ = Σ f4_{abcd} · S_{ab} · S_{cd} is computed in camera frame where S = R·H⁻¹·Rᵀ.

### The depth soft mode conclusion

With analytical Q₄, the saddlepoint optimizer converges identically to Laplace (5 iterations, |grad|→1e-10). The results reveal a fundamental limitation:

```
                           True      Laplace  Saddlepoint
  t₃ (depth)          0.000000     6.451333     6.366716

  Translation error:
    Laplace:      6.454793 m
    Saddlepoint:  6.370215 m
    Difference:   0.084597 m
```

The saddlepoint correction shifts the estimate by only ~8.5cm on a ~6.4m error. The per-landmark corrections c₁ ≈ 2% are mathematically correct — they capture the non-Gaussian skewness of each landmark marginal — but they cannot resolve the **depth soft mode**, which is a geometric degeneracy: with a single camera, the projection π(x) = [x/z, y/z] is nearly invariant under joint camera-landmark depth rescaling.

**This is an information-theoretic limitation, not an approximation artifact.** The saddlepoint correctly refines the *shape* of each landmark's marginal, but no per-landmark correction can create depth information that isn't present in the observations. Resolving the depth soft mode requires:

- **Multi-camera / stereo baselines** — parallax provides direct depth triangulation
- **Strong pose priors** — constrain the translation subspace directly
- **Known-scale constraints** — e.g., known inter-landmark distances or IMU-derived scale

The multi-camera experiments (Experiments 3–4) confirm this: with 4+ cameras on a ring, depth is well-constrained and both Laplace and saddlepoint give accurate results, with the SP correction properly accounting for residual non-Gaussianity.

## Mathematica Verification

Three Mathematica scripts provide symbolic and numerical verification of every identity and formula. All results are independently confirmed by the Rust test suite (95 tests).

| Script | Scope | Key result |
|--------|-------|------------|
| `SE3AlgebraVerification.m` | 24 SO(3)/SE(3)/BCH identities | Symbolic ✓ and numerical to 10⁻¹⁶ |
| `SaddlepointVerification.m` | Saddlepoint formula vs quadrature | 3000× improvement over Laplace (3.2×10⁻⁶ rel. error) |
| `CouplingJacobianDerivation.m` | d[S⁻¹T]/dΩ symbolic proof | 9/9 Ω substitutions → exact zero |

See **[verification/README.md](verification/README.md)** for detailed per-part results and Mathematica pitfalls.

## Publication Target

Robotics or estimation journal (IEEE Transactions on Robotics, IJRR, or similar). The paper is a methods contribution applicable across:

- Visual-inertial odometry and IMU preintegration
- Spacecraft attitude determination
- Surgical robot registration
- Multi-body dynamics and motion capture

## Relationship to Prior Work

- **Kuehnel (2004)**: *AIP Conf. Proc.* vol. 735, pp. 176–186. [DOI: 10.1063/1.1835212](https://doi.org/10.1063/1.1835212) — First introduction of Lie-Cartan exponential coordinates for unbiased Bayesian pose estimation
- **Kuehnel (2005)**: *AIP Conf. Proc.* vol. 803, pp. 318–329. [DOI: 10.1063/1.2149810](https://doi.org/10.1063/1.2149810) — Local frame junction trees in SLAM
- **Kuehnel (2006)**: Tech. Rep., NASA Ames / USRA-RIACS — [Full SE(3) algebra toolkit](paper/robustEst.pdf) (BCH, Jacobians, phase reflection)
- **Solà, Deray & Atchuthan (2018)**: ["A micro Lie theory for state estimation in robotics"](https://arxiv.org/pdf/1812.01537) — First-order SE(3) machinery; our work extends to second order
- **Ye & Chirikjian (2024)**: ["Uncertainty propagation on unimodular Lie groups"](https://openreview.net/forum?id=duNh060j1J) — Propagation via SDEs; our approach uses closed-form BCH
- **Barfoot (2024)**: [*State Estimation for Robotics*](https://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser24.pdf) — Comprehensive textbook; our coupling Jacobian and saddlepoint go beyond what is covered

## License

Paper content: All rights reserved, Excel Solutions LLC.

Code: MIT License (see `rust/LICENSE` when available).

## Contacts

Frank O. Kuehnel – Excel Solutions LLC<br>
Andre Jalobeanu – Bayesmap Inc.<br>
Email: kuehnelf@gmail.com
