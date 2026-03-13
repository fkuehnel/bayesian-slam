# Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on SE(3)

**Author:** Frank O. Kuehnel — Excel Solutions LLC

---

## Overview

This repository accompanies the paper *Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on the SE(3) Lie Group*, which develops a self-contained framework for Bayesian inference over rigid body poses.

**Paper:** [SE3_inference_paper.pdf](paper/SE3_inference.pdf)

The core ideas originate from work at NASA Ames Research Center (2008) on robust pose estimation using the SE(3) Lie group structure. This project refocuses that foundational material into three standalone contributions aimed at the broader estimation and inference community. All formulas have been verified symbolically (Mathematica) and numerically (Rust finite differences, Python, Monte Carlo).

## What This Paper Contributes

### 1. Second-Order Uncertainty Propagation on SE(3)

Current filters and optimizers (EKF, IEKF, factor graphs) propagate pose uncertainty using first-order Jacobians. On SE(3), the semi-direct product coupling between rotation and translation introduces systematic bias that first-order methods miss. We derive:

- Closed-form composition Jacobians including the rotation-translation coupling Jacobian **J_t(Ω, t)**
- Second-order corrections to both the mean and covariance of composed poses
- Explicit formulas using the finite BCH composition via SU(2), avoiding truncated series

**Verified:** Monte Carlo validation (10⁵ samples) confirms a 1.83× improvement in covariance accuracy at σ_ω ≈ 0.2 rad. BCH truncation errors verified across 4 decades (see Table 1 in paper).

### 2. Saddlepoint Marginalization for Projective Observations

When camera observations are projective (x/z, y/z), the joint posterior over pose and 3D landmarks is non-Gaussian in the pose parameters. Standard Laplace approximation (Schur complement in bundle adjustment) misses this asymmetry. We develop:

- A saddlepoint correction with the **corrected formula** c₁ = (1/12)A + (1/8)B − (1/8)Q₄
- Two distinct cubic contraction types (cross: 6/15 Isserlis pairings, trace: 9/15)
- Quartic term with correct negative sign (from exp(−g₄) ≈ 1 − g₄)
- Validity guard: correction applied only when |c₁| < 0.5 (σ_depth/depth ≲ 0.3)

**Verified:** Against numerical quadrature, the saddlepoint achieves 6 significant figures (3.2×10⁻⁶ relative error) compared to Laplace's 0.94% error — a ~3000× improvement.

### 3. Complete SE(3) Algebra Toolkit

The paper provides a self-contained reference for the SE(3) machinery, all symbolically verified:

- Exponential/logarithmic maps with Lie-Cartan coordinates of the first kind
- Closed-form finite Rodrigues vector composition via the SU(2)/Z₂ ≅ SO(3) double cover
- Phase reflection handling at Θ = π with cutline analysis
- Wei-Norman formula and both left/right Rodrigues Jacobians
- Key identities: S⁻¹ = J_ωl, S⁻¹R = RS⁻¹ = J_ωr (verified to machine precision)
- det S = 2(1−cosΘ)/Θ² (strictly positive, guarantees invertibility)

## Errata from Original Tech Report

The verification process uncovered errors in the original 2008 formulation:

1. **Third derivative sign** (Appendix E, Eq. thirdderivs): ∂³u/∂x₃'³ = **−6u/x₃'³** (was +6u/x₃'³)
2. **Saddlepoint formula** (Appendix E, Eq. saddlepointfull): The original (5/24)A + (1/8)Q₄ was incorrect. The correct formula is **(1/12)A + (1/8)B − (1/8)Q₄** with two distinct contraction types and a sign flip on the quartic term.
3. **Coupling Jacobian** (Appendix C, Eq. dSinvTdOmega): The original formula used the exponential coordinate **t** where the physical translation **T** was required in the leading skew-symmetric term. A corrected T-form with a new scalar α = (sinΘ − Θ)/(2(1−cosΘ)) has been derived and verified. **Symbolic proof confirms both the T-form and the corrected t-form (with S−S⁻¹(−Ω) structure) are algebraically equivalent** when T = S(Ω)t.

## Repository Structure

```
.
├── README.md                          # This file
├── paper/
│   ├── SE3_inference_paper.tex        # Main manuscript (LaTeX, ~2200 lines)
│   └── robustEst.tex                  # Original 2008 technical report (reference)
├── verification/
│   ├── SE3AlgebraVerification.m       # Mathematica: 24 verified identities (SO3/SE3/BCH/Jacobians)
│   ├── SaddlepointVerification.m      # Mathematica: projective derivs, saddlepoint formula, quadrature
│   └── CouplingJacobianDerivation.m   # Mathematica: symbolic derivation of d[S⁻¹T]/dΩ, equivalence proof
├── rust/
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs                     # Array-based linalg: Vec3/Mat3/Vec6/Mat6
│   │   ├── so3.rs                     # SO(3): Rodrigues exp/log, S matrix, S⁻¹
│   │   ├── se3.rs                     # SE(3): Pose struct, compose/inverse/act/exp/log
│   │   ├── bch.rs                     # Finite BCH via SU(2) quaternions, phase reflection, Jacobians
│   │   ├── jacobians.rs               # J_ωr, J_ωl, J_t coupling (analytic T-form), 6×6 SE(3) Jacobian
│   │   ├── projective.rs              # Pinhole camera, derivatives through 3rd order, third cumulants
│   │   ├── saddlepoint.rs             # Landmark optimization, saddlepoint correction, validity guard
│   │   └── propagation.rs             # First/second-order covariance transport, MC validation
│   └── examples/
│       ├── bias_experiment.rs         # Experiment 1: L1 vs L2 coordinate bias (21× ratio)
│       ├── multicam_saddlepoint.rs    # Multi-camera saddlepoint with point cloud & camera sweep
│       └── multicam_experiment.rs     # Extended multi-camera analysis
└── experiments/                       # Planned: scripts reproducing paper figures
```

## Rust Implementation

**88 tests passing, 0 failures, 0 warnings.**

| Module | Tests | Status | Description |
|--------|-------|--------|-------------|
| `so3` | 8 | ✅ | Rodrigues exp/log, hat/vee, S matrix, S⁻¹ |
| `se3` | 8 | ✅ | Pose compose/inverse/act, exp/log roundtrip |
| `bch` | 30 | ✅ | Finite BCH via SU(2), phase reflection, composition Jacobians |
| `jacobians` | 13 | ✅ | J_ωr/J_ωl, analytic J_t (T-form), 6×6 SE(3) Jacobian (FD-verified) |
| `projective` | 11 | ✅ | Project, Jacobian, Hessian, 3rd derivs, third cumulants (all FD-verified) |
| `saddlepoint` | 7 | ✅ | Landmark GN optimizer, corrected c₁ formula, validity guard, Q₄ by FD |
| `propagation` | 10 | ✅ | First/second-order covariance transport, Levi-Civita, Isserlis correction, MC validation |

### Key verified identities (Rust + Mathematica + Python)

- BCH composition matches matrix log(R_a R_b) to 10⁻¹⁶ across all angle regimes
- Phase reflection at Θ > π produces correct (−r, 2π−Θ) identification
- Composition Jacobians ∂Ω_c/∂Ω_a, ∂Ω_c/∂Ω_b match FD to 10⁻⁸
- S⁻¹R = J_ωr verified symbolically and numerically
- Full 6×6 SE(3) Jacobian block structure [[J_ωr, 0], [J_t, J_ωr]] verified to 10⁻⁹
- **Coupling Jacobian J_t**: analytic T-form verified to 7.3×10⁻¹⁰ against FD; symbolic equivalence with t-form proven in Mathematica
- Projective derivatives through 3rd order verified against FD
- Saddlepoint correction matches numerical quadrature to 6 significant figures

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

### Pending Work

| Item | Priority | Description |
|------|----------|-------------|
| `propagation.rs` | ✅ Done | First/second-order covariance propagation with MC validation (88 tests total) |
| Multi-camera saddlepoint | ✅ Done | `multicam_saddlepoint.rs`: 2–12 cameras, point cloud, camera-count sweep |
| Experiment figures | High | Generate paper figures from Rust (bias scatter, propagation accuracy, convergence) |
| Closed-form Q₄ | Low | Replace FD-based quartic contraction with analytic fourth derivatives of projective model |
| t-form J_t option | Low | Add alternative `j_coupling_tform` using the S−S⁻¹(−Ω) structure (proven equivalent) |

### Resolved Items

| Item | Resolution |
|------|------------|
| Closed-form J_t | ✅ Derived T-form with α = (sinΘ−Θ)/(2(1−cosΘ)), verified to 7.3×10⁻¹⁰ |
| T-form / t-form equivalence | ✅ Proven symbolically in Mathematica (`CouplingJacobianDerivation.m`, Step 5) |
| Saddlepoint correction formula | ✅ Corrected to c₁ = (1/12)A + (1/8)B − (1/8)Q₄, verified against quadrature |
| Third derivative sign | ✅ Fixed: ∂³u/∂x₃'³ = −6u/x₃'³ |
| SE(3) Jacobian FD test | ✅ Was `#[ignore]`, now passing (bug was in original j_coupling, not the formula) |
| det S formula | ✅ Added det S = 2(1−cosΘ)/Θ² to paper Appendix A |

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

2. **det S = 2(1 − cos Θ)/Θ²** — Exact volume element for Lie-Cartan coordinate chart. Needed for density transformations between [Ω, t] and [Ω, T] coordinates. Continuous at Θ = 0 (det S → 1), vanishes at Θ = π.

3. **Compact T-form coupling Jacobian** — Single-line formula for ∂[S⁻¹T]/∂Ω using scalar coefficients x, α and outer products, replacing Barfoot's four-line iterated cross-product expansion. Proven algebraically equivalent in Mathematica (CouplingJacobianDerivation.m, Step 5).

## Multi-Camera Saddlepoint Experiment

The `multicam_saddlepoint` example demonstrates saddlepoint-corrected landmark marginalization with 2–12 cameras observing a shared 3D point cloud in a ring configuration.

### Results: Camera-count sweep (test point at origin)

| N_cam | views | σ_depth | σ/depth | c₁ |
|------:|------:|--------:|--------:|---:|
| 2 | 2 | 0.028 | 0.0037 | 7.5×10⁻³ |
| 3 | 3 | 0.031 | 0.0041 | −1.2×10⁻⁶ |
| 4 | 4 | 0.028 | 0.0037 | −3.5×10⁻⁵ |
| 6 | 6 | 0.023 | 0.0031 | −2.4×10⁻⁵ |
| 8 | 8 | 0.020 | 0.0026 | −1.8×10⁻⁵ |
| 12 | 12 | 0.016 | 0.0022 | −1.2×10⁻⁵ |

**Key finding:** The stereo case (2 cameras) produces the largest saddlepoint correction (~200× larger than 4-camera), confirming that the correction matters most in the typical multi-view operating regime where depth is constrained by triangulation parallax.

### Stereo full point cloud (2 cameras, 30 landmarks)

- 100% saddlepoint validity across all landmarks
- Individual corrections up to c₁ ≈ 9×10⁻³
- Cumulative correction: Σ|Δ| = 4.4×10⁻²
- Correction dominated by negative Q₄ term (depth non-Gaussianity)

### Usage

```bash
cargo run --release --example multicam_saddlepoint [N_cameras] [N_landmarks]
cargo run --release --example multicam_saddlepoint 4 20   # default
cargo run --release --example multicam_saddlepoint 2 30   # stereo
```

## Mathematica Verification Scripts

### SE3AlgebraVerification.m — 24 identities verified

| Part | What | Result |
|------|------|--------|
| 1–2 | H algebra, Rodrigues formula, R^T R = I, det = 1 | Symbolic ✓ |
| 3–4 | S·S⁻¹ = I, S⁻¹ = J_ωl, S⁻¹R = J_ωr, det S, det J_ωr | Symbolic ✓ |
| 5 | BCH finite composition vs matrix log(R_a R_b) | 10⁻¹⁶ (5 cases) |
| 6 | SU(2) quaternions, double cover R(q) = R(−q) | 10⁻¹⁶ |
| 7 | J_ωr from BCH differentiation vs analytic | 6×10⁻¹⁰ |
| 8 | SE(3) exp/log roundtrip, compose, inverse | 10⁻¹⁶ |
| 9–10 | 6×6 SE(3) Jacobian, J_t = d[S⁻¹T]/dΩ · J_ωr | 1.7×10⁻⁹ |
| 11 | Composition Jacobians ∂Ω_c/∂Ω_a, ∂Ω_c/∂Ω_b | 4.7×10⁻⁸ |
| 12 | Ad(f·g) = Ad(f)·Ad(g) | 3×10⁻¹⁶ |
| 13 | BCH truncation order: err₂ ~ s³, err₃ ~ s⁴ | Confirmed across 4 decades |
| 14 | Second-order covariance correction vs MC (10⁵ samples) | 1.83× improvement |
| 15 | Phase reflection at Θ > π | 10⁻¹⁶ |

### SaddlepointVerification.m — Corrected formula derivation

| Part | What | Result |
|------|------|--------|
| 1–2 | Symbolic π derivatives through 3rd order, sign check | −6u/x₃'³ confirmed |
| 3–5 | Numerical quadrature ground truth | log I = −1.8506 |
| 6 | Corrected c₁ = (1/12)A + (1/8)B − (1/8)Q₄ | c₁ = 0.0095 |
| 7 | Laplace 0.94% error, saddlepoint 3.2×10⁻⁶ error | 3000× improvement |
| 7b | Prior strength sweep, regime classification | (σ/depth)² scaling confirmed |
| 9 | Symbolic P×Hess decomposition of f''' | Exact match |
| 10 | Depth scaling analysis | (σ_z/depth)² confirmed |

### CouplingJacobianDerivation.m — d[S⁻¹T]/dΩ from first principles

| Step | What | Result |
|------|------|--------|
| 1–2 | Symbolic S⁻¹, chain rule d/dΩ through (θ, r) | Exact 3×3 expression |
| 3 | Symbolic derivative vs FD | 1.74×10⁻⁸ |
| 4 | Three candidate formulas vs FD | (numerical evaluation has scoping issue) |
| 5 | **T-form vs exact symbolic derivative** | Verified (see symbolic output) |
| 5 | **T-form(T=St) − t-form = 0** | **ZERO (proven)** — algebraic equivalence |
| 6 | Diagnosis of original Rust bug | Formula correct; implementation error |

## Publication Target

Robotics or estimation journal (IEEE Transactions on Robotics, IJRR, or similar). The paper is a methods contribution applicable across:

- Visual-inertial odometry and IMU preintegration
- Spacecraft attitude determination
- Surgical robot registration
- Multi-body dynamics and motion capture

## Relationship to Prior Work

- **Kuehnel (2004)**: *AIP Conf. Proc.* vol. 735, pp. 176–186. [DOI: 10.1063/1.1835212](https://doi.org/10.1063/1.1835212) — First introduction of Lie-Cartan exponential coordinates for unbiased Bayesian pose estimation
- **Kuehnel (2005)**: *AIP Conf. Proc.* vol. 803, pp. 318–329. [DOI: 10.1063/1.2149810](https://doi.org/10.1063/1.2149810) — Local frame junction trees in SLAM
- **Kuehnel (2008)**: Tech. Rep., NASA Ames / USRA-RIACS — Full SE(3) algebra toolkit (BCH, Jacobians, phase reflection)
- **Solà, Deray & Atchuthan (2018)**: ["A micro Lie theory for state estimation in robotics"](https://arxiv.org/pdf/1812.01537) — First-order SE(3) machinery; our work extends to second order
- **Ye & Chirikjian (2024)**: "Uncertainty propagation on unimodular Lie groups" — Propagation via SDEs; our approach uses closed-form BCH
- **Barfoot (2024)**: [*State Estimation for Robotics*](https://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser24.pdf) — Comprehensive textbook; our coupling Jacobian and saddlepoint go beyond what is covered

## License

Paper content: All rights reserved, Frank O. Kuehnel / Excel Solutions LLC.

Code: MIT License (see `rust/LICENSE` when available).

## Contact

Frank O. Kuehnel — Excel Solutions LLC
Email: kuehnelf@gmail.com
