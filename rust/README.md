# Rust Implementation

**95 tests passing, 0 failures, 0 warnings.**

| Module | Tests | Status | Description |
|--------|-------|--------|-------------|
| `so3` | 8 | ✅ | Rodrigues exp/log, hat/vee, S matrix, S⁻¹ |
| `se3` | 8 | ✅ | Pose compose/inverse/act, exp/log roundtrip |
| `bch` | 30 | ✅ | Finite BCH via SU(2), phase reflection, composition Jacobians |
| `jacobians` | 13 | ✅ | J_ωr/J_ωl, analytic J_t (T-form), 6×6 SE(3) Jacobian (FD-verified) |
| `projective` | 14 | ✅ | Project, Jacobian, Hessian, 3rd/4th derivs, third cumulants, analytical Q₄ (all FD-verified) |
| `saddlepoint` | 10 | ✅ | Landmark GN optimizer, corrected c₁ formula, validity guard, analytical Q₄, multi-cam |
| `propagation` | 10 | ✅ | First/second-order covariance transport, ε-identity Isserlis correction, MC validation |

## Key Verified Identities (Rust + Mathematica + Python)

- BCH composition matches matrix log(R_a R_b) to 10⁻¹⁶ across all angle regimes
- Phase reflection at Θ > π produces correct (−r, 2π−Θ) identification
- Composition Jacobians ∂Ω_c/∂Ω_a, ∂Ω_c/∂Ω_b match FD to 10⁻⁸
- S⁻¹R = J_ωr verified symbolically and numerically
- Full 6×6 SE(3) Jacobian block structure [[J_ωr, 0], [J_t, J_ωr]] verified to 10⁻⁹
- **Coupling Jacobian J_t**: analytic T-form verified to 7.3×10⁻¹⁰ against FD; symbolic equivalence proven in Mathematica (9/9 substitutions)
- Projective derivatives through 4th order verified against FD
- Saddlepoint correction matches numerical quadrature to 6 significant figures
- Analytical Q₄ (quartic contraction) verified against FD of exact Hessian
- Second-order covariance correction via ε-ε identity matches brute-force Levi-Civita contraction to 10⁻¹⁴

## Design Principles

- **No external dependencies for core algebra**: 3×3 and 6×6 matrices as fixed-size arrays, exploiting skew-symmetric and block-triangular structure
- **Zero heap allocation in inner loops**: SO(3) Jacobians are `[[f64; 3]; 3]`, SE(3) ones are `[[f64; 6]; 6]`
- **Defensive numerics**: Validity guard on saddlepoint (|c₁| > 0.5 → fall back to Laplace), Taylor expansions near singularities (Θ → 0, ε → 0)
- **Algebraic index elimination**: Second-order Isserlis corrections use the ε-ε determinant identity instead of Levi-Civita loops, reducing the covariance correction from O(3⁶) to two 3×3 matrix multiplies (~20× fewer flops)
- **Exhaustive FD validation**: Every Jacobian and derivative tested against central finite differences

## Computational Cost: Saddlepoint vs Laplace

The saddlepoint correction adds ~12,600 FLOPs per landmark on top of Laplace's ~860 FLOPs — a **15.6× overhead** that buys ~3000× better marginal accuracy.

### Per-landmark flop budget

| Step | FLOPs | Path |
|------|------:|------|
| GN optimization (3 iterations) | ~666 | Shared |
| Final Hessian + NLL + det | ~193 | Shared |
| **Laplace total** | **~860** | |
| `third_cumulants()` | ~3,780 | SP only |
| — camera-frame contraction | 864 | |
| — rotation to world (3³×3³ tensor) | 2,916 | |
| `quartic_contraction_analytical()` | ~4,985 | SP only |
| — R·H⁻¹·Rᵀ sandwich | 162 | |
| — 4-nested contraction (3⁴ × 2² × 14) | 4,779 | |
| `saddlepoint_correction_full()` | ~3,778 | SP only |
| — λ tensor (cross-contraction A) | 3,645 | |
| — trace-contraction B + c₁ | 133 | |
| H⁻¹ + ln(1+c₁) | ~33 | SP only |
| **Saddlepoint total** | **~13,440** | |

### Cost breakdown by component

The three saddlepoint-only stages contribute roughly equally:

- **Quartic contraction Q₄** — 40% of overhead (4,985 FLOPs). The 4-nested loop over 3D indices with 2×2 observation dimensions is unavoidable for the full f⁴·H⁻¹·H⁻¹ tensor contraction.
- **Third cumulants κ₃** — 30% (3,780 FLOPs). Dominated by the rotation of the 3×3×3 cumulant tensor from camera to world frame.
- **Cross-contraction A** — 29% (3,645 FLOPs). The λ_{abc} = f³·H⁻¹·H⁻¹·H⁻¹ intermediate tensor requires O(3⁶) work for the full 6-index contraction.

### Multi-camera scaling

With N cameras, cumulant computations scale linearly while the correction formula stays fixed:

| Component | Cost |
|-----------|------|
| Per camera | `third_cumulants` + `quartic_contraction_analytical` = ~8,765 FLOPs |
| Fixed | `saddlepoint_correction_full` = ~3,778 FLOPs |
| **Total SP overhead** | **8,765N + 3,778 FLOPs** |

### Practical impact

In a SLAM system with L landmarks, the saddlepoint adds ~12,600L FLOPs per evaluation — negligible compared to the O(n³) sparse Cholesky factorization that dominates bundle adjustment. The correction is computed per landmark independently, making it trivially parallelizable.

## Getting Started

```bash
git clone https://github.com/fkuehnel/bayesian-slam.git
cd bayesian-slam/rust
cargo build --release
cargo test
```

## Examples and Experiments

Four experiment scripts validate the theory and produce data for paper figures. Full descriptions, usage, and numerical results are in **[../experiments/README.md](../experiments/README.md)**.

| Example | What it does |
|---------|-------------|
| `bias_experiment` | L1 vs L2 coordinate bias (21× ratio, 2000 MC samples) |
| `pose_inference` | Single/multi-camera pose estimation, Laplace vs SP (1–12 cameras, configurable landmarks and overlap) |
| `multicam_saddlepoint` | Camera-count sweep showing c₁ scaling (2–12 cameras) |
| `multicam_experiment` | SP correction validated against MC integration (200k samples) |

**Quick start:**
```bash
cargo run --release --example pose_inference              # single camera
cargo run --release --example pose_inference 4 20 0.8     # 4 cameras, 20 landmarks, 80% overlap
```

**Key findings** (details in [../experiments/README.md](../experiments/README.md)):
- Analytical Q₄ fixes SP convergence: 60+ iterations → 5 iterations
- Single camera: depth soft mode limits translation accuracy (~6.4m error), SP correction is only ~8.5cm — information-theoretic, not approximation, limitation
- 4 cameras: depth triangulated (~0.5m error), c₁ drops from 2×10⁻² to 2×10⁻⁴

Paper figures use the **Matplotlib → PGF** pipeline. See [../experiments/README.md](../experiments/README.md) for the full plotting guide, style conventions, and reproduction instructions.
