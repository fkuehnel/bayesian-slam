# Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on SE(3)

**Author:** Frank O. Kuehnel

---

## Overview

This repository accompanies the paper *Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on the SE(3) Lie Group*, which develops a self-contained framework for Bayesian inference over rigid body poses.

**Paper:** [SE3_inference_paper.pdf](paper/highorder.pdf)

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
│   ├── SE3_inference_paper.tex        # Main manuscript (LaTeX, ~1980 lines)
│   └── robustEst.tex                  # Original 2008 technical report (reference)
├── verification/
│   ├── SE3AlgebraVerification.m       # Mathematica: 24 verified identities (SO3/SE3/BCH/Jacobians)
│   ├── SaddlepointVerification.m      # Mathematica: projective derivs, saddlepoint formula, quadrature
│   └── CouplingJacobianDerivation.m   # Mathematica: symbolic derivation of d[S⁻¹T]/dΩ, equivalence proof
├── rust/
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs                     # Array-based linalg: Vec3/Mat3/Vec6/Mat6
│       ├── so3.rs                     # SO(3): Rodrigues exp/log, S matrix, S⁻¹
│       ├── se3.rs                     # SE(3): Pose struct, compose/inverse/act/exp/log
│       ├── bch.rs                     # Finite BCH via SU(2) quaternions, phase reflection, Jacobians
│       ├── jacobians.rs              # J_ωr, J_ωl, J_t coupling (analytic T-form), 6×6 SE(3) Jacobian
│       ├── projective.rs             # Pinhole camera, derivatives through 3rd order, third cumulants
│       └── saddlepoint.rs            # Landmark optimization, saddlepoint correction, validity guard
└── experiments/                       # Planned: scripts reproducing paper figures
```

## Rust Implementation

**75 tests passing, 0 failures, 0 warnings.**

| Module | Tests | Status | Description |
|--------|-------|--------|-------------|
| `so3` | 8 | ✅ | Rodrigues exp/log, hat/vee, S matrix, S⁻¹ |
| `se3` | 8 | ✅ | Pose compose/inverse/act, exp/log roundtrip |
| `bch` | 30 | ✅ | Finite BCH via SU(2), phase reflection, composition Jacobians |
| `jacobians` | 13 | ✅ | J_ωr/J_ωl, analytic J_t (T-form), 6×6 SE(3) Jacobian (FD-verified) |
| `projective` | 11 | ✅ | Project, Jacobian, Hessian, 3rd derivs, third cumulants (all FD-verified) |
| `saddlepoint` | 7 | ✅ | Landmark GN optimizer, corrected c₁ formula, validity guard, Q₄ by FD |

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
cd se3-inference/rust

cargo build --release
cargo test
cargo test -- --nocapture  # see diagnostic output
```

### Pending Work

| Item | Priority | Description |
|------|----------|-------------|
| `propagation.rs` | High | First/second-order covariance propagation module (entry point: `bch::compose_jacobians`) |
| Experiment figures | High | Generate paper figures from Rust (bias scatter, propagation accuracy, convergence) |
| Multi-camera saddlepoint | Medium | Extend saddlepoint to landmarks seen from multiple cameras (sum of information matrices) |
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
- **Solà, Deray & Atchuthan (2018)**: "A micro Lie theory for state estimation in robotics" — First-order SE(3) machinery; our work extends to second order
- **Ye & Chirikjian (2024)**: "Uncertainty propagation on unimodular Lie groups" — Propagation via SDEs; our approach uses closed-form BCH
- **Barfoot (2024)**: *State Estimation for Robotics* — Comprehensive textbook; our coupling Jacobian and saddlepoint go beyond what is covered

## License

Paper content: All rights reserved, Frank O. Kuehnel / Excel Solutions LLC.

Code: MIT License (see `rust/LICENSE` when available).

## Contact

Frank O. Kuehnel — Excel Solutions LLC
Email: kuehnelf@gmail.com