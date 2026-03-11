//! # SO(3) — Special Orthogonal Group in 3D
//!
//! Rodrigues formula, exponential/logarithmic maps, skew-symmetric matrices.
//!
//! Follows Appendix A of the paper:
//! - Eq. (rodriguesapp): R(Ω) = I + sinΘ H + (1-cosΘ) H²
//! - Eq. (expnotation): Θ = |Ω|, r = Ω/Θ, H = Ĥ·r

use crate::*;

/// Small angle threshold for Taylor expansion.
const EPS: f64 = 1e-10;

// ─── Skew-symmetric (Hat/Vee) ───

/// Hat map: ω ∈ ℝ³ → [ω]× ∈ so(3).
///
/// ```text
/// [ω]× = |  0  -ω₃  ω₂ |
///         |  ω₃  0  -ω₁ |
///         | -ω₂  ω₁  0  |
/// ```
#[inline]
pub fn hat(w: &Vec3) -> Mat3 {
    [[ 0.0, -w[2],  w[1]],
     [ w[2],  0.0, -w[0]],
     [-w[1],  w[0],  0.0]]
}

/// Vee map: [ω]× ∈ so(3) → ω ∈ ℝ³.
#[inline]
pub fn vee(m: &Mat3) -> Vec3 {
    [m[2][1], m[0][2], m[1][0]]
}

// ─── Rodrigues formula ───

/// Rotation matrix from Rodrigues vector Ω ∈ ℝ³.
///
/// R(Ω) = I + sinΘ·H + (1 - cosΘ)·H²
///
/// where Θ = |Ω|, H = [r]× with r = Ω/Θ.
///
/// For small Θ, uses Taylor expansion to avoid division by zero.
pub fn exp(omega: &Vec3) -> Mat3 {
    let theta_sq = dot3(omega, omega);
    let theta = theta_sq.sqrt();

    if theta < EPS {
        // R ≈ I + [ω]× + ½[ω]×²
        let h = hat(omega);
        let h2 = mm3(&h, &h);
        add_mat3(&add_mat3(&I3, &h), &scale_mat3(0.5, &h2))
    } else {
        let h = hat(&scale3(1.0/theta, omega)); // H = [r]×
        let h2 = mm3(&h, &h);
        // R = I + sinΘ·H + (1 - cosΘ)·H²
        add_mat3(
            &add_mat3(&I3, &scale_mat3(theta.sin(), &h)),
            &scale_mat3(1.0 - theta.cos(), &h2),
        )
    }
}

/// Logarithmic map: R ∈ SO(3) → Ω ∈ ℝ³.
///
/// Uses Tr(R) = 1 + 2cosΘ to recover angle,
/// and the antisymmetric part of R to recover axis.
pub fn log(r: &Mat3) -> Vec3 {
    let cos_theta = 0.5 * (trace3(r) - 1.0);
    // Clamp for numerical safety
    let cos_theta = cos_theta.clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    if theta < EPS {
        // Small angle: Ω ≈ vee(R - Rᵀ)/2
        let rt = transpose3(r);
        let skew = sub_mat3(r, &rt);
        [0.5*skew[2][1], 0.5*skew[0][2], 0.5*skew[1][0]]
    } else if (std::f64::consts::PI - theta).abs() < EPS {
        // Near π: special extraction from diagonal
        // r_i = sqrt((R_ii + 1)/2), sign from off-diagonal
        let mut r_axis = [0.0; 3];
        // Find largest diagonal element for numerical stability
        let diag = [r[0][0], r[1][1], r[2][2]];
        let k = if diag[0] >= diag[1] && diag[0] >= diag[2] { 0 }
                else if diag[1] >= diag[2] { 1 }
                else { 2 };
        r_axis[k] = ((diag[k] + 1.0) * 0.5).sqrt();
        let inv = 0.5 / r_axis[k];
        for i in 0..3 {
            if i != k {
                r_axis[i] = r[i][k] * inv;
            }
        }
        scale3(theta, &r_axis)
    } else {
        // General case: Ω = θ/(2sinθ) · vee(R - Rᵀ)
        let rt = transpose3(r);
        let skew = sub_mat3(r, &rt);
        let factor = theta / (2.0 * theta.sin());
        [factor*skew[2][1], factor*skew[0][2], factor*skew[1][0]]
    }
}

/// The S(Ω) matrix connecting exponential coordinates to translation:
/// T = S(Ω)·t.
///
/// S(Ω) = I + H² + (1/Θ)(I - R)H
///
/// where R = R(Ω), H = [r]×, Θ = |Ω|.
pub fn s_matrix(omega: &Vec3) -> Mat3 {
    let theta_sq = dot3(omega, omega);
    let theta = theta_sq.sqrt();

    if theta < EPS {
        // S ≈ I + ½[ω]×
        let h = hat(omega);
        add_mat3(&I3, &scale_mat3(0.5, &h))
    } else {
        let r_unit = scale3(1.0/theta, omega);
        let h = hat(&r_unit);
        let h2 = mm3(&h, &h);
        let r = exp(omega);
        let i_minus_r = sub_mat3(&I3, &r);
        // S = I + H² + (1/Θ)(I - R)H
        add_mat3(
            &add_mat3(&I3, &h2),
            &scale_mat3(1.0/theta, &mm3(&i_minus_r, &h)),
        )
    }
}

/// Inverse of S(Ω), using the closed-form identity:
///
/// S⁻¹(Ω) = (1-x)I + x·r·rᵀ - ½[Ω]×
///
/// where 1-x = (Θ/2)/tan(Θ/2).
///
/// This is identical to J_ωl(Ω) = J_ωr(-Ω).
/// See Eq. (Sinvidentity) in the paper.
pub fn s_inv(omega: &Vec3) -> Mat3 {
    let theta_sq = dot3(omega, omega);
    let theta = theta_sq.sqrt();

    if theta < EPS {
        // S⁻¹ ≈ I - ½[ω]×
        let h = hat(omega);
        sub_mat3(&I3, &scale_mat3(0.5, &h))
    } else {
        let r = scale3(1.0/theta, omega);
        let half_theta = 0.5 * theta;
        let one_minus_x = half_theta / half_theta.tan();
        let x = 1.0 - one_minus_x;
        // S⁻¹ = (1-x)I + x·rrᵀ - ½[Ω]×
        let rrt = outer3(&r, &r);
        let omega_hat = hat(omega);
        add_mat3(
            &add_mat3(
                &scale_mat3(one_minus_x, &I3),
                &scale_mat3(x, &rrt),
            ),
            &scale_mat3(-0.5, &omega_hat),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq_mat3(a: &Mat3, b: &Mat3, tol: f64) -> bool {
        for i in 0..3 { for j in 0..3 {
            if (a[i][j] - b[i][j]).abs() > tol { return false; }
        }}
        true
    }

    fn approx_eq_vec3(a: &Vec3, b: &Vec3, tol: f64) -> bool {
        (0..3).all(|i| (a[i] - b[i]).abs() < tol)
    }

    #[test]
    fn test_exp_identity() {
        let r = exp(&[0.0, 0.0, 0.0]);
        assert!(approx_eq_mat3(&r, &I3, 1e-12));
    }

    #[test]
    fn test_exp_90deg_z() {
        // 90° rotation about z-axis
        let omega = [0.0, 0.0, std::f64::consts::FRAC_PI_2];
        let r = exp(&omega);
        // Should map x→y, y→-x
        let x = mv3(&r, &[1.0, 0.0, 0.0]);
        assert!(approx_eq_vec3(&x, &[0.0, 1.0, 0.0], 1e-10));
    }

    #[test]
    fn test_exp_log_roundtrip() {
        let omega = [0.3, -0.5, 0.7];
        let r = exp(&omega);
        let omega_back = log(&r);
        assert!(approx_eq_vec3(&omega, &omega_back, 1e-10));
    }

    #[test]
    fn test_exp_log_near_pi() {
        // Angle close to π
        let omega = [3.0, 0.1, 0.0];
        let theta = norm3(&omega);
        assert!(theta > 3.0); // close to π ≈ 3.14159
        let r = exp(&omega);
        let omega_back = log(&r);
        let r_back = exp(&omega_back);
        assert!(approx_eq_mat3(&r, &r_back, 1e-8));
    }

    #[test]
    fn test_rotation_is_orthogonal() {
        let omega = [0.5, -1.2, 0.8];
        let r = exp(&omega);
        let rrt = mm3(&r, &transpose3(&r));
        assert!(approx_eq_mat3(&rrt, &I3, 1e-10));
    }

    #[test]
    fn test_s_s_inv_roundtrip() {
        let omega = [0.4, -0.6, 0.3];
        let s = s_matrix(&omega);
        let si = s_inv(&omega);
        let product = mm3(&s, &si);
        assert!(approx_eq_mat3(&product, &I3, 1e-10));
    }

    #[test]
    fn test_hat_vee_roundtrip() {
        let w = [1.0, -2.0, 3.0];
        let m = hat(&w);
        let w_back = vee(&m);
        assert!(approx_eq_vec3(&w, &w_back, 1e-15));
    }

    #[test]
    fn test_hat_is_skew_symmetric() {
        let w = [0.5, -0.3, 0.7];
        let m = hat(&w);
        let mt = transpose3(&m);
        let sum = add_mat3(&m, &mt);
        assert!(approx_eq_mat3(&sum, &Z3, 1e-15));
    }
}
