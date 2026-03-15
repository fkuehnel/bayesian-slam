# Mathematica Verification Scripts

This directory contains Mathematica notebooks that symbolically and numerically verify every identity and formula in the paper. All results have been independently confirmed by the Rust test suite (95 tests).

## SE3AlgebraVerification.m — 24 identities verified

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

## SaddlepointVerification.m — Corrected formula derivation

| Part | What | Result |
|------|------|--------|
| 1–2 | Symbolic π derivatives through 3rd order, sign check | −6u/x₃'³ confirmed |
| 3–5 | Numerical quadrature ground truth | log I = −1.8506 |
| 6 | Corrected c₁ = (1/12)A + (1/8)B − (1/8)Q₄ | c₁ = 0.0095 |
| 7 | Laplace 0.94% error, saddlepoint 3.2×10⁻⁶ error | 3000× improvement |
| 7b | Prior strength sweep, regime classification | (σ/depth)² scaling confirmed |
| 9 | Symbolic P×Hess decomposition of f''' | Exact match |
| 10 | Depth scaling analysis | (σ_z/depth)² confirmed |

## CouplingJacobianDerivation.m — d[S⁻¹T]/dΩ symbolic proof

| Step | What | Result |
|------|------|--------|
| Setup | Build S⁻¹ directly from Ω = (w₁,w₂,w₃), Θ = √(w·w) implicit | No (θ,r) split |
| Step 1 | Exact derivative via Mathematica D[] through Sqrt[] | 3×3 symbolic |
| Step 2 | T-form as 4 separate terms (cannot add symbolically) | Unsimplified |
| Step 3 | Substitute 9 Ω vectors, simplify each term, add, compare | **9/9 ZERO** |
| Coeffs | α half-angle form, β identity, α ≠ β, erratum disproof | All PROVEN |

### Mathematica pitfalls documented in script header

1. Must parametrize in Ω directly (not θ, r with |r|=1)
2. Must NOT Simplify expressions with symbolic Sqrt[w₁²+w₂²+w₃²]
3. Must NOT add T-form terms symbolically (auto-simplification drops terms)
4. FullSimplify fails on trig(n) vs trig(n/2) — use numerical fallback
