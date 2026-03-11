//! # Baker-Campbell-Hausdorff: composing finite rotations and motions
//!
//! Paper Appendix B. The BCH series for SO(3) can be evaluated in
//! closed form using the SU(2)/Z₂ ≅ SO(3) quaternion representation.
//!
//! ## Key formula (Eq. addSO3finite)
//!
//! Given two rotations R_a = exp(Ĥ·Ω_a) and R_b = exp(Ĥ·Ω_b),
//! the composed Rodrigues vector Ω_c = ln(R_a R_b) is:
//!
//! ```text
//!   Ω_c = (Θ_c / C) · (a r_a + b r_b + ab r_a × r_b)
//!
//!   a = tan(Θ_a/2),  b = tan(Θ_b/2)
//!   A = 1 + a² + b² + a²b²
//!   ε = 1 - ab cos ϑ,  cos ϑ = r_a · r_b
//!   C = √(A - ε²)
//!   tan(Θ_c/2) = C / ε
//! ```
//!
//! ## Phase reflection at Θ = π
//!
//! The logarithmic map loses uniqueness at Θ = π (the topological
//! origin is SU(2)/Z₂). A cutline at cos ϑ = 0 handles the
//! sign flip: vectors crossing the π-sphere boundary reflect
//! to the antipodal direction.

use crate::*;
use crate::so3;

const EPS: f64 = 1e-10;
const PI: f64 = std::f64::consts::PI;

// =========================================================================
// SO(3) finite composition: Ω_c = ln(R_a · R_b) via SU(2) quaternions
// =========================================================================

/// Result of a finite rotation composition.
#[derive(Debug, Clone, Copy)]
pub struct BchResult {
    /// Composed Rodrigues vector Ω_c.
    pub omega_c: Vec3,
    /// Composed angle Θ_c = |Ω_c|.
    pub theta_c: f64,
    /// Normalizing factor C.
    pub big_c: f64,
    /// Epsilon factor ε = 1 - ab cos ϑ.
    pub epsilon: f64,
}

/// Compose two finite Rodrigues vectors: Ω_c = ln(R(Ω_a) · R(Ω_b)).
///
/// Uses the closed-form formula via SU(2) quaternion representation
/// (paper Eq. addSO3finite), with careful handling of limit cases
/// and phase reflection at Θ = π.
///
/// This avoids intermediate rotation matrices entirely — the composition
/// is performed directly in the Rodrigues vector space.
pub fn compose_rodrigues(omega_a: &Vec3, omega_b: &Vec3) -> Vec3 {
    compose_rodrigues_full(omega_a, omega_b).omega_c
}

/// Full composition result including diagnostic quantities.
pub fn compose_rodrigues_full(omega_a: &Vec3, omega_b: &Vec3) -> BchResult {
    let theta_a = norm3(omega_a);
    let theta_b = norm3(omega_b);

    // ─── Trivial cases ───
    if theta_a < EPS {
        return BchResult {
            omega_c: *omega_b,
            theta_c: theta_b,
            big_c: 0.0,
            epsilon: 1.0,
        };
    }
    if theta_b < EPS {
        return BchResult {
            omega_c: *omega_a,
            theta_c: theta_a,
            big_c: 0.0,
            epsilon: 1.0,
        };
    }

    // ─── Unit axis vectors ───
    let r_a = scale3(1.0 / theta_a, omega_a);
    let r_b = scale3(1.0 / theta_b, omega_b);
    let cos_vartheta = dot3(&r_a, &r_b);

    // ─── Half-angle tangents ───
    let a = (theta_a * 0.5).tan();
    let b = (theta_b * 0.5).tan();

    // ─── Normalizing quantities ───
    let a2 = a * a;
    let b2 = b * b;
    let big_a = 1.0 + a2 + b2 + a2 * b2;
    let epsilon = 1.0 - a * b * cos_vartheta;
    let c_sq = big_a - epsilon * epsilon;

    // Ensure non-negative (numerical safety)
    let c_sq = c_sq.max(0.0);
    let big_c = c_sq.sqrt();

    // ─── Composed angle Θ_c via tan(Θ_c/2) = C/ε ───
    // atan2(C, ε) gives half-angle in (-π, π].
    // When ε < 0 (composition crosses π-boundary), θ_c > π.
    let raw_theta_c = 2.0 * big_c.atan2(epsilon);

    // ─── Direction vector: a·r_a + b·r_b + ab·r_a×r_b ───
    let cross_ab = cross3(&r_a, &r_b);
    let dir = [
        a * r_a[0] + b * r_b[0] + a * b * cross_ab[0],
        a * r_a[1] + b * r_b[1] + a * b * cross_ab[1],
        a * r_a[2] + b * r_b[2] + a * b * cross_ab[2],
    ];

    // ─── Prefactor Θ_c / C with limit handling ───
    let prefactor = theta_c_over_c(raw_theta_c.abs(), big_c, epsilon, big_a);

    // Raw Ω_c (may have |Ω_c| > π)
    let sign = if raw_theta_c >= 0.0 { 1.0 } else { -1.0 };
    let mut omega_c = scale3(sign * prefactor, &dir);
    let mut theta_c = norm3(&omega_c);

    // ─── Phase reflection: if |Ω_c| > π, apply (r,Θ) → (-r, 2π-Θ) ───
    // This is the Z₂ identification from SU(2)/Z₂ ≅ SO(3).
    // Opposite points on the π-sphere boundary represent the same rotation.
    if theta_c > PI + EPS {
        let reflected_theta = 2.0 * PI - theta_c;
        omega_c = scale3(-reflected_theta / theta_c, &omega_c);
        theta_c = reflected_theta;
    }

    BchResult { omega_c, theta_c, big_c, epsilon }
}

/// Compute Θ_c/C with proper limit handling.
///
/// Two singular limits (paper §B.1):
///   1. ε → 0 (A ≥ 1): Θ_c/C → π/√A · [1 - (2/π)(ε/√A) + ...]
///   2. C → 0 (ε² ≥ 1): Θ_c/C → (2/ε) · [1 - (1/3)(C/ε)² + ...]
///
/// The general case uses Θ_c/C = 2·atan2(C,ε) / C.
fn theta_c_over_c(theta_c: f64, big_c: f64, epsilon: f64, big_a: f64) -> f64 {
    let eps_abs = epsilon.abs();

    if big_c < EPS && eps_abs < EPS {
        // Both zero: degenerate, return safe value
        return 2.0;
    }

    if eps_abs < 1e-6 * big_a.sqrt() {
        // ─── Limit 1: ε → 0 ───
        // Θ_c/C = π/√A · [1 - (2/π)(ε/√A) + ½(ε/√A)² + ...]
        let sqrt_a = big_a.sqrt();
        let u = epsilon / sqrt_a;
        PI / sqrt_a * (1.0 - (2.0 / PI) * u + 0.5 * u * u)
    } else if big_c < 1e-6 * eps_abs {
        // ─── Limit 2: C → 0 ───
        // Θ_c/C = (2/ε) · [1 - (1/3)(C/ε)² + ...]
        let v = big_c / epsilon;
        (2.0 / epsilon) * (1.0 - v * v / 3.0)
    } else {
        // ─── General case ───
        // Θ_c / C = 2 atan2(C, ε) / C
        theta_c / big_c
    }
}

// =========================================================================
// Unit quaternion helpers: SU(2)/Z₂ ≅ SO(3)
// =========================================================================

/// Unit quaternion q = (w, x, y, z) stored as [w, x, y, z].
pub type Quat = [f64; 4];

/// Rodrigues vector → unit quaternion.
///
/// q = (cos(Θ/2), sin(Θ/2)·r)  where Θ = |Ω|, r = Ω/Θ.
/// Always returns the q with w ≥ 0 (canonical hemisphere).
pub fn rodrigues_to_quat(omega: &Vec3) -> Quat {
    let theta = norm3(omega);
    if theta < EPS {
        // q ≈ (1, Ω/2)  (first-order Taylor)
        return [1.0, omega[0] * 0.5, omega[1] * 0.5, omega[2] * 0.5];
    }
    let half = theta * 0.5;
    let s = half.sin() / theta;
    let mut q = [half.cos(), s * omega[0], s * omega[1], s * omega[2]];
    // Canonical hemisphere: w ≥ 0
    if q[0] < 0.0 {
        for qi in q.iter_mut() { *qi = -*qi; }
    }
    q
}

/// Unit quaternion → Rodrigues vector.
///
/// Ω = 2·arctan(|v|/w) · v/|v|  where q = (w, v).
/// Handles the π-boundary via atan2.
pub fn quat_to_rodrigues(q: &Quat) -> Vec3 {
    let v = [q[1], q[2], q[3]];
    let vn = norm3(&v);
    let w = q[0];
    if vn < EPS {
        // Near identity: Ω ≈ 2v
        return scale3(2.0, &v);
    }
    let theta = 2.0 * vn.atan2(w.abs());
    // If w < 0, the angle is > π, phase-reflect: (r, Θ) → (r, 2π-Θ)
    let theta = if w >= 0.0 { theta } else { 2.0 * PI - theta };
    scale3(theta / vn, &v)
}

/// Quaternion multiplication (Hamilton product).
///
/// The SU(2) group operation. Corresponds to rotation composition:
/// R(q_a q_b) = R(q_a) R(q_b).
pub fn quat_mul(qa: &Quat, qb: &Quat) -> Quat {
    [
        qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3],
        qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2],
        qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1],
        qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0],
    ]
}

/// Quaternion conjugate (inverse for unit quaternions).
pub fn quat_conj(q: &Quat) -> Quat {
    [q[0], -q[1], -q[2], -q[3]]
}

/// Quaternion norm (should be 1 for unit quaternions).
pub fn quat_norm(q: &Quat) -> f64 {
    (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]).sqrt()
}

// =========================================================================
// Composition Jacobians: ∂Ω_c/∂Ω_a and ∂Ω_c/∂Ω_b
// =========================================================================
//
// Given Ω_c = log(exp(Ω_a) · exp(Ω_b)), the chain rule on SO(3) gives:
//
//   ∂Ω_c/∂Ω_a = J_r(Ω_c) · R_b^T · J_r⁻¹(Ω_a)
//   ∂Ω_c/∂Ω_b = J_r(Ω_c) · J_r⁻¹(Ω_b)
//
// Derivation (for ∂Ω_c/∂Ω_b):
//   exp(Ω_b + δ_b) = exp(Ω_b) · exp(J_r⁻¹(Ω_b) δ_b)
//   ⟹ R_c' = R_a R_b exp(J_r⁻¹(Ω_b) δ_b) = R_c exp(η)
//   ⟹ Ω_c + dΩ_c = Ω_c + J_r(Ω_c) η
//
// Derivation (for ∂Ω_c/∂Ω_a):
//   exp(Ω_a + δ_a) = exp(Ω_a) exp(J_r⁻¹(Ω_a) δ_a)
//   ⟹ R_c' = R_a exp(η_a) R_b = R_c · R_b^{-1} exp(η_a) R_b
//           = R_c · exp(R_b^T η_a)
//   ⟹ dΩ_c = J_r(Ω_c) · R_b^T · J_r⁻¹(Ω_a) · δ_a

use crate::jacobians;

/// Composition Jacobian ∂Ω_c/∂Ω_a  (3×3).
///
/// J_r(Ω_c) · R_b^T · J_r⁻¹(Ω_a)
pub fn d_omega_c_d_omega_a(omega_a: &Vec3, omega_b: &Vec3) -> Mat3 {
    let omega_c = compose_rodrigues(omega_a, omega_b);
    let r_b = so3::exp(omega_b);
    let r_bt = transpose3(&r_b);
    let jr_c = jacobians::j_omega_right(&omega_c);
    let jr_inv_a = jacobians::j_omega_right_inv(omega_a);
    mm3(&mm3(&jr_c, &r_bt), &jr_inv_a)
}

/// Composition Jacobian ∂Ω_c/∂Ω_b  (3×3).
///
/// J_r(Ω_c) · J_r⁻¹(Ω_b)
pub fn d_omega_c_d_omega_b(omega_a: &Vec3, omega_b: &Vec3) -> Mat3 {
    let omega_c = compose_rodrigues(omega_a, omega_b);
    let jr_c = jacobians::j_omega_right(&omega_c);
    let jr_inv_b = jacobians::j_omega_right_inv(omega_b);
    mm3(&jr_c, &jr_inv_b)
}

/// Both composition Jacobians at once (avoids redundant BCH evaluation).
///
/// Returns (∂Ω_c/∂Ω_a, ∂Ω_c/∂Ω_b, Ω_c).
pub fn compose_jacobians(omega_a: &Vec3, omega_b: &Vec3)
    -> (Mat3, Mat3, Vec3)
{
    let omega_c = compose_rodrigues(omega_a, omega_b);
    let jr_c = jacobians::j_omega_right(&omega_c);
    let r_bt = transpose3(&so3::exp(omega_b));
    let jr_inv_a = jacobians::j_omega_right_inv(omega_a);
    let jr_inv_b = jacobians::j_omega_right_inv(omega_b);

    let dc_da = mm3(&mm3(&jr_c, &r_bt), &jr_inv_a);
    let dc_db = mm3(&jr_c, &jr_inv_b);
    (dc_da, dc_db, omega_c)
}

// =========================================================================
// SE(3) finite composition in exponential coordinates
// =========================================================================

/// Compose two SE(3) elements in exponential coordinates.
///
/// Given ξ_a = [Ω_a, t_a] and ξ_b = [Ω_b, t_b], compute
/// ξ_c = ln((R_a,T_a) · (R_b,T_b)) = [Ω_c, t_c].
///
/// Paper Eq. (addfiniteSE3):
///   Ω_c = ln(R_a R_b)                          [via compose_rodrigues]
///   t_c = S⁻¹(Ω_c) · (R_a T_b + T_a)
///
/// where T_a = S(Ω_a)·t_a, T_b = S(Ω_b)·t_b.
pub fn compose_se3(xi_a: &Vec6, xi_b: &Vec6) -> Vec6 {
    let omega_a = [xi_a[0], xi_a[1], xi_a[2]];
    let t_a = [xi_a[3], xi_a[4], xi_a[5]];
    let omega_b = [xi_b[0], xi_b[1], xi_b[2]];
    let t_b = [xi_b[3], xi_b[4], xi_b[5]];

    // Physical translations: T = S(Ω)·t
    let big_t_a = mv3(&so3::s_matrix(&omega_a), &t_a);
    let big_t_b = mv3(&so3::s_matrix(&omega_b), &t_b);

    // Rotation composition via BCH
    let omega_c = compose_rodrigues(&omega_a, &omega_b);

    // Composed translation: T_c = R_a T_b + T_a
    let r_a = so3::exp(&omega_a);
    let big_t_c = add3(&mv3(&r_a, &big_t_b), &big_t_a);

    // Back to exponential coordinates: t_c = S⁻¹(Ω_c) T_c
    let t_c = mv3(&so3::s_inv(&omega_c), &big_t_c);

    [omega_c[0], omega_c[1], omega_c[2], t_c[0], t_c[1], t_c[2]]
}

// =========================================================================
// BCH series expansion (for second-order propagation)
// =========================================================================

/// Second-order BCH expansion of the rotational part:
///   Ω_c ≈ Ω_a + Ω_b + ½ Ω_a × Ω_b + O(|Ω|³)
///
/// Used for deriving the second-order mean correction
/// (paper Eq. BCHexpanded, Appendix D).
pub fn bch_second_order(omega_a: &Vec3, omega_b: &Vec3) -> Vec3 {
    let cross = cross3(omega_a, omega_b);
    add3(
        &add3(omega_a, omega_b),
        &scale3(0.5, &cross),
    )
}

/// Third-order BCH expansion:
///   Ω_c ≈ Ω_a + Ω_b + ½ Ω_a×Ω_b
///         + (1/12)(Ω_a×(Ω_a×Ω_b) + Ω_b×(Ω_b×Ω_a)) + O(|Ω|⁴)
pub fn bch_third_order(omega_a: &Vec3, omega_b: &Vec3) -> Vec3 {
    let ab = cross3(omega_a, omega_b);
    let ba = scale3(-1.0, &ab);
    let aab = cross3(omega_a, &ab);
    let bba = cross3(omega_b, &ba);
    let third = scale3(1.0 / 12.0, &add3(&aab, &bba));
    add3(
        &add3(&add3(omega_a, omega_b), &scale3(0.5, &ab)),
        &third,
    )
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::se3::Pose;

    fn approx_eq_vec3(a: &Vec3, b: &Vec3, tol: f64) -> bool {
        (0..3).all(|i| (a[i] - b[i]).abs() < tol)
    }
    fn approx_eq_mat3(a: &Mat3, b: &Mat3, tol: f64) -> bool {
        for i in 0..3 { for j in 0..3 {
            if (a[i][j] - b[i][j]).abs() > tol { return false; }
        }}
        true
    }

    // ─── Basic composition tests ───

    #[test]
    fn compose_identity_left() {
        let w = [0.5, -0.3, 0.8];
        let wc = compose_rodrigues(&[0.0; 3], &w);
        assert!(approx_eq_vec3(&wc, &w, 1e-12));
    }

    #[test]
    fn compose_identity_right() {
        let w = [0.5, -0.3, 0.8];
        let wc = compose_rodrigues(&w, &[0.0; 3]);
        assert!(approx_eq_vec3(&wc, &w, 1e-12));
    }

    #[test]
    fn compose_inverse_gives_zero() {
        let w = [0.5, -0.3, 0.8];
        let w_neg = scale3(-1.0, &w);
        // R(Ω)·R(-Ω) = I, so ln = 0
        let wc = compose_rodrigues(&w, &w_neg);
        assert!(norm3(&wc) < 1e-10,
                "Ω + (-Ω) should be ~0, got {:?}", wc);
    }

    // ─── Agreement with matrix multiplication ───

    #[test]
    fn compose_agrees_with_matrix_mult() {
        let cases: Vec<(Vec3, Vec3)> = vec![
            ([0.3, -0.5, 0.7], [0.1, 0.2, -0.4]),
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),    // orthogonal axes
            ([0.01, 0.02, -0.01], [0.03, -0.01, 0.02]), // small angles
            ([2.5, 0.0, 0.0], [0.0, 0.5, 0.0]),     // near π
            ([0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),     // parallel
        ];
        for (wa, wb) in &cases {
            // Matrix route: log(exp(Ω_a) · exp(Ω_b))
            let r_a = so3::exp(wa);
            let r_b = so3::exp(wb);
            let r_c = mm3(&r_a, &r_b);
            let omega_matrix = so3::log(&r_c);

            // BCH route
            let omega_bch = compose_rodrigues(wa, wb);

            // Compare: the rotation matrices should match
            let r_bch = so3::exp(&omega_bch);
            assert!(approx_eq_mat3(&r_c, &r_bch, 1e-9),
                    "Mismatch at Ω_a={:?}, Ω_b={:?}\n  matrix: {:?}\n  BCH: {:?}",
                    wa, wb, omega_matrix, omega_bch);
        }
    }

    // ─── Collinear case: cos ϑ = ±1 ───

    #[test]
    fn compose_collinear_same_axis() {
        // Same axis: should simply add angles
        let w_a = [0.3, 0.0, 0.0];
        let w_b = [0.5, 0.0, 0.0];
        let wc = compose_rodrigues(&w_a, &w_b);
        assert!(approx_eq_vec3(&wc, &[0.8, 0.0, 0.0], 1e-10),
                "Collinear same: {:?}", wc);
    }

    #[test]
    fn compose_collinear_opposite_axis() {
        // Opposite axis
        let w_a = [0.5, 0.0, 0.0];
        let w_b = [-0.3, 0.0, 0.0];
        let wc = compose_rodrigues(&w_a, &w_b);
        assert!(approx_eq_vec3(&wc, &[0.2, 0.0, 0.0], 1e-10),
                "Collinear opposite: {:?}", wc);
    }

    // ─── Near π boundary ───

    #[test]
    fn compose_near_pi() {
        // Both rotations large but sum < π
        let w_a = [2.0, 0.0, 0.0];
        let w_b = [0.0, 1.0, 0.0];
        let omega_bch = compose_rodrigues(&w_a, &w_b);
        let r_bch = so3::exp(&omega_bch);
        let r_mat = mm3(&so3::exp(&w_a), &so3::exp(&w_b));
        assert!(approx_eq_mat3(&r_mat, &r_bch, 1e-8),
                "Near-π mismatch");
    }

    #[test]
    fn compose_past_pi_reflects() {
        // Two rotations whose sum exceeds π
        let w_a = [2.8, 0.0, 0.0];
        let w_b = [0.8, 0.0, 0.0];
        // The result should be equivalent via matrix mult
        let omega_bch = compose_rodrigues(&w_a, &w_b);
        let r_bch = so3::exp(&omega_bch);
        let r_mat = mm3(&so3::exp(&w_a), &so3::exp(&w_b));
        assert!(approx_eq_mat3(&r_mat, &r_bch, 1e-8),
                "Past-π: BCH={:?}, theta_c={:.4}", omega_bch, norm3(&omega_bch));
    }

    // ─── Small angle BCH series ───

    #[test]
    fn bch_second_order_small_angle() {
        let wa = [0.01, 0.02, -0.01];
        let wb = [0.03, -0.01, 0.02];
        let exact = compose_rodrigues(&wa, &wb);
        let approx = bch_second_order(&wa, &wb);
        // At small angles, second-order should be very accurate
        let err = norm3(&sub3(&exact, &approx));
        assert!(err < 1e-5,
                "2nd-order BCH err={:.2e} at small angle", err);
    }

    #[test]
    fn bch_third_order_better_than_second() {
        let wa = [0.1, 0.2, -0.1];
        let wb = [0.15, -0.1, 0.2];
        let exact = compose_rodrigues(&wa, &wb);
        let e2 = norm3(&sub3(&exact, &bch_second_order(&wa, &wb)));
        let e3 = norm3(&sub3(&exact, &bch_third_order(&wa, &wb)));
        assert!(e3 < e2,
                "3rd order ({:.2e}) should beat 2nd ({:.2e})", e3, e2);
    }

    // ─── SE(3) composition ───

    #[test]
    fn se3_compose_agrees_with_pose() {
        let xi_a: Vec6 = [0.3, -0.2, 0.5, 1.0, -0.5, 0.3];
        let xi_b: Vec6 = [0.1, 0.4, -0.3, 0.5, 1.0, -0.7];

        // Via Pose (matrix route)
        let f_a = Pose::exp(&xi_a);
        let f_b = Pose::exp(&xi_b);
        let f_c = f_a.compose(&f_b);
        let xi_pose = f_c.log();

        // Via BCH (direct exponential coordinate route)
        let xi_bch = compose_se3(&xi_a, &xi_b);

        // Should agree
        for i in 0..6 {
            assert!((xi_pose[i] - xi_bch[i]).abs() < 1e-9,
                    "SE3 compose mismatch [{}]: pose={:.6} bch={:.6}",
                    i, xi_pose[i], xi_bch[i]);
        }
    }

    #[test]
    fn se3_compose_identity() {
        let xi = [0.5, -0.3, 0.7, 2.0, 1.0, -1.0];
        let zero = [0.0; 6];
        let left = compose_se3(&zero, &xi);
        let right = compose_se3(&xi, &zero);
        for i in 0..6 {
            assert!((left[i] - xi[i]).abs() < 1e-12);
            assert!((right[i] - xi[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn se3_compose_inverse() {
        let xi = [0.3, -0.5, 0.7, 1.0, -2.0, 0.5];
        let f = Pose::exp(&xi);
        let f_inv = f.inverse();
        let xi_inv = f_inv.log();
        let composed = compose_se3(&xi, &xi_inv);
        let norm: f64 = composed.iter().map(|x| x*x).sum::<f64>().sqrt();
        assert!(norm < 1e-9, "f·f⁻¹ should be ~0, got norm={:.2e}", norm);
    }

    // ─── Associativity ───

    #[test]
    fn compose_is_associative() {
        let wa = [0.3, -0.2, 0.5];
        let wb = [0.1, 0.4, -0.3];
        let wc_in = [0.2, -0.1, 0.6];

        // (a·b)·c
        let ab = compose_rodrigues(&wa, &wb);
        let abc_left = compose_rodrigues(&ab, &wc_in);

        // a·(b·c)
        let bc = compose_rodrigues(&wb, &wc_in);
        let abc_right = compose_rodrigues(&wa, &bc);

        // Should match (up to phase)
        let r_left = so3::exp(&abc_left);
        let r_right = so3::exp(&abc_right);
        assert!(approx_eq_mat3(&r_left, &r_right, 1e-10),
                "Associativity failed");
    }

    // ─── Diagnostic: BchResult fields ───

    #[test]
    fn bch_result_theta_c_matches_norm() {
        let wa = [0.5, -0.3, 0.8];
        let wb = [0.2, 0.4, -0.1];
        let res = compose_rodrigues_full(&wa, &wb);
        let norm = norm3(&res.omega_c);
        assert!((res.theta_c - norm).abs() < 1e-9,
                "theta_c={:.6} vs |Ω_c|={:.6}", res.theta_c, norm);
    }

    // ─── Quaternion helpers ───

    #[test]
    fn quat_roundtrip_small() {
        let w = [0.1, -0.2, 0.3];
        let q = rodrigues_to_quat(&w);
        assert!((quat_norm(&q) - 1.0).abs() < 1e-14, "not unit");
        let w2 = quat_to_rodrigues(&q);
        assert!(approx_eq_vec3(&w, &w2, 1e-12),
                "roundtrip: {:?} vs {:?}", w, w2);
    }

    #[test]
    fn quat_roundtrip_near_pi() {
        let w = [2.9, 0.0, 0.0];
        let q = rodrigues_to_quat(&w);
        assert!((quat_norm(&q) - 1.0).abs() < 1e-14);
        let w2 = quat_to_rodrigues(&q);
        // Compare rotation matrices (Rodrigues vector may differ by sign near π)
        let r1 = so3::exp(&w);
        let r2 = so3::exp(&w2);
        assert!(approx_eq_mat3(&r1, &r2, 1e-10),
                "near-π roundtrip failed");
    }

    #[test]
    fn quat_roundtrip_identity() {
        let w = [0.0, 0.0, 0.0];
        let q = rodrigues_to_quat(&w);
        assert!((q[0] - 1.0).abs() < 1e-14);
        let w2 = quat_to_rodrigues(&q);
        assert!(norm3(&w2) < 1e-12);
    }

    #[test]
    fn quat_mul_matches_matrix() {
        let wa = [0.5, -0.3, 0.8];
        let wb = [0.2, 0.6, -0.4];
        let qa = rodrigues_to_quat(&wa);
        let qb = rodrigues_to_quat(&wb);
        let qc = quat_mul(&qa, &qb);
        // Normalize (should already be unit)
        let n = quat_norm(&qc);
        assert!((n - 1.0).abs() < 1e-13, "product not unit: {:.15}", n);
        // Compare with matrix product
        let r_mat = mm3(&so3::exp(&wa), &so3::exp(&wb));
        let _wc_mat = so3::log(&r_mat);
        let wc_quat = quat_to_rodrigues(&qc);
        let r_quat = so3::exp(&wc_quat);
        assert!(approx_eq_mat3(&r_mat, &r_quat, 1e-10),
                "quat_mul disagrees with matrix mult");
    }

    #[test]
    fn quat_conj_is_inverse() {
        let w = [0.5, -0.3, 0.8];
        let q = rodrigues_to_quat(&w);
        let qi = quat_conj(&q);
        let prod = quat_mul(&q, &qi);
        // Should be (1, 0, 0, 0)
        assert!((prod[0] - 1.0).abs() < 1e-14);
        assert!(prod[1].abs() < 1e-14);
        assert!(prod[2].abs() < 1e-14);
        assert!(prod[3].abs() < 1e-14);
    }

    #[test]
    fn quat_compose_matches_bch() {
        // Verify that quaternion composition gives the same result
        // as compose_rodrigues (the whole point of the SU(2) approach)
        let wa = [0.7, -0.4, 1.2];
        let wb = [-0.3, 0.9, 0.1];
        let qa = rodrigues_to_quat(&wa);
        let qb = rodrigues_to_quat(&wb);
        let qc = quat_mul(&qa, &qb);
        let wc_quat = quat_to_rodrigues(&qc);
        let wc_bch = compose_rodrigues(&wa, &wb);
        // Compare at the rotation matrix level
        let r_q = so3::exp(&wc_quat);
        let r_b = so3::exp(&wc_bch);
        assert!(approx_eq_mat3(&r_q, &r_b, 1e-10),
                "quat route vs BCH route disagree");
    }

    // ─── Composition Jacobians (finite difference validated) ───

    fn fd_jacobian_a(omega_a: &Vec3, omega_b: &Vec3, h: f64) -> Mat3 {
        let mut jac = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut wp = *omega_a;
            let mut wm = *omega_a;
            wp[j] += h;
            wm[j] -= h;
            let cp = compose_rodrigues(&wp, omega_b);
            let cm = compose_rodrigues(&wm, omega_b);
            for i in 0..3 {
                jac[i][j] = (cp[i] - cm[i]) / (2.0 * h);
            }
        }
        jac
    }

    fn fd_jacobian_b(omega_a: &Vec3, omega_b: &Vec3, h: f64) -> Mat3 {
        let mut jac = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut wp = *omega_b;
            let mut wm = *omega_b;
            wp[j] += h;
            wm[j] -= h;
            let cp = compose_rodrigues(omega_a, &wp);
            let cm = compose_rodrigues(omega_a, &wm);
            for i in 0..3 {
                jac[i][j] = (cp[i] - cm[i]) / (2.0 * h);
            }
        }
        jac
    }

    #[test]
    fn jacobian_a_fd_moderate_angles() {
        let wa = [0.3, -0.5, 0.7];
        let wb = [0.1, 0.2, -0.4];
        let analytic = d_omega_c_d_omega_a(&wa, &wb);
        let numeric = fd_jacobian_a(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-5,
                    "dΩ_c/dΩ_a[{},{}]: analytic={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_b_fd_moderate_angles() {
        let wa = [0.3, -0.5, 0.7];
        let wb = [0.1, 0.2, -0.4];
        let analytic = d_omega_c_d_omega_b(&wa, &wb);
        let numeric = fd_jacobian_b(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-5,
                    "dΩ_c/dΩ_b[{},{}]: analytic={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_a_fd_small_angles() {
        let wa = [0.01, 0.02, -0.01];
        let wb = [0.03, -0.01, 0.02];
        let analytic = d_omega_c_d_omega_a(&wa, &wb);
        let numeric = fd_jacobian_a(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-5,
                    "small dΩ_c/dΩ_a[{},{}]: a={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_b_fd_small_angles() {
        let wa = [0.01, 0.02, -0.01];
        let wb = [0.03, -0.01, 0.02];
        let analytic = d_omega_c_d_omega_b(&wa, &wb);
        let numeric = fd_jacobian_b(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-5,
                    "small dΩ_c/dΩ_b[{},{}]: a={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_a_fd_large_angles() {
        let wa = [1.5, -0.8, 0.3];
        let wb = [0.2, 1.2, -0.7];
        let analytic = d_omega_c_d_omega_a(&wa, &wb);
        let numeric = fd_jacobian_a(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-4,
                    "large dΩ_c/dΩ_a[{},{}]: a={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_b_fd_large_angles() {
        let wa = [1.5, -0.8, 0.3];
        let wb = [0.2, 1.2, -0.7];
        let analytic = d_omega_c_d_omega_b(&wa, &wb);
        let numeric = fd_jacobian_b(&wa, &wb, 1e-7);
        for i in 0..3 { for j in 0..3 {
            assert!((analytic[i][j] - numeric[i][j]).abs() < 1e-4,
                    "large dΩ_c/dΩ_b[{},{}]: a={:.6} fd={:.6}",
                    i, j, analytic[i][j], numeric[i][j]);
        }}
    }

    #[test]
    fn jacobian_a_identity_at_zero_b() {
        // When Ω_b = 0: R_c = R_a, so ∂Ω_c/∂Ω_a = I
        let wa = [0.5, -0.3, 0.7];
        let wb = [0.0, 0.0, 0.0];
        let jac = d_omega_c_d_omega_a(&wa, &wb);
        for i in 0..3 { for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((jac[i][j] - expected).abs() < 1e-9,
                    "dΩ_c/dΩ_a should be I when Ω_b=0: [{},{}]={:.6}",
                    i, j, jac[i][j]);
        }}
    }

    #[test]
    fn jacobian_b_identity_at_zero_a() {
        // When Ω_a = 0: R_c = R_b, so ∂Ω_c/∂Ω_b = I
        let wa = [0.0, 0.0, 0.0];
        let wb = [0.5, -0.3, 0.7];
        let jac = d_omega_c_d_omega_b(&wa, &wb);
        for i in 0..3 { for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((jac[i][j] - expected).abs() < 1e-9,
                    "dΩ_c/dΩ_b should be I when Ω_a=0: [{},{}]={:.6}",
                    i, j, jac[i][j]);
        }}
    }

    #[test]
    fn compose_jacobians_consistent() {
        let wa = [0.4, -0.6, 0.2];
        let wb = [0.3, 0.1, -0.5];
        let (ja, jb, wc) = compose_jacobians(&wa, &wb);
        let ja2 = d_omega_c_d_omega_a(&wa, &wb);
        let jb2 = d_omega_c_d_omega_b(&wa, &wb);
        let wc2 = compose_rodrigues(&wa, &wb);
        assert!(approx_eq_mat3(&ja, &ja2, 1e-12));
        assert!(approx_eq_mat3(&jb, &jb2, 1e-12));
        assert!(approx_eq_vec3(&wc, &wc2, 1e-12));
    }
}
