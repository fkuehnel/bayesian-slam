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
└── experiments/                       # Experiments, plotting, and figure reproduction
    ├── README.md                      # Full experiment docs, results, plotting guide
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
```

## Examples and Experiments

Four experiment scripts validate the theory and produce data for paper figures. Full descriptions, usage, and numerical results are in **[experiments/README.md](experiments/README.md)**.

| Example | What it does |
|---------|-------------|
| `bias_experiment` | L1 vs L2 coordinate bias (21× ratio, 2000 MC samples) |
| `pose_inference` | Single/multi-camera pose estimation, Laplace vs SP (1–12 cameras, configurable landmarks and overlap) |
| `multicam_saddlepoint` | Camera-count sweep showing c₁ scaling (2–12 cameras) |
| `multicam_experiment` | SP correction validated against MC integration (200k samples) |

**Quick start:**
```bash
cd rust
cargo run --release --example pose_inference              # single camera
cargo run --release --example pose_inference 4 20 0.8     # 4 cameras, 20 landmarks, 80% overlap
```

**Key findings** (details in [experiments/README.md](experiments/README.md)):
- Analytical Q₄ fixes SP convergence: 60+ iterations → 5 iterations
- Single camera: depth soft mode limits translation accuracy (~6.4m error), SP correction is only ~8.5cm — information-theoretic, not approximation, limitation
- 4 cameras: depth triangulated (~0.5m error), c₁ drops from 2×10⁻² to 2×10⁻⁴

Paper figures use the **Matplotlib → PGF** pipeline. See [experiments/README.md](experiments/README.md) for the full plotting guide, style conventions, and reproduction instructions.

### Pending Work

| Item | Priority | Description |
|------|----------|-------------|
| Remaining plot scripts | High | `plot_propagation.py`, `plot_convergence.py`, `plot_saddlepoint.py`, `plot_sp_regime.py` |
| Good Notations | High | Good consistent notations that the target audience knows |
| Experimental findings | High | Our paper must be based on repeatable experiments |
| Verify Solà notation | Medium | Mathematica scripts to confirm notation mappings |
| Depth soft mode analysis | Medium | Investigate multi-camera or prior-based approaches to resolve the translation degeneracy |

### Resolved Items

| Item | Resolution |
|------|------------|
| Plotting workflow | ✅ Documented in [experiments/README.md](experiments/README.md) |
| Pose inference example | ✅ Single + multi-camera with configurable CLI args |
| Closed-form J_t | ✅ T-form verified to 7.3×10⁻¹⁰; symbolic proof in Mathematica (9/9 subs) |
| Saddlepoint formula | ✅ c₁ = (1/12)A + (1/8)B − (1/8)Q₄, verified against quadrature |
| Closed-form Q₄ | ✅ Analytical quartic contraction; SP converges in 5 iterations (was 60+) |
| Pose inference panics | ✅ Fixed NaN/degenerate configurations, added pose prior |

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
