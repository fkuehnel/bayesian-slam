//! # Composition Jacobians on SO(3) and SE(3)
//!
//! Paper Appendix C. Key results:
//! - J_ωr(Ω) = (1-x)I + x·rrᵀ + ½[Ω]×        Eq. (RodweinormanJ)
//! - J_ωl(Ω) = (1-x)I + x·rrᵀ - ½[Ω]×        (left variant)
//! - S⁻¹(Ω) = J_ωl(Ω) = J_ωr(-Ω)             Eq. (Sinvidentity)
//! - S⁻¹(Ω)R(Ω) = J_ωr(Ω)                    Eq. (SinvRidentity)
//! - Full SE(3) Jacobian: block lower-triangular with J_t coupling
//!
//! where 1-x = (Θ/2)/tan(Θ/2), Θ = |Ω|, r = Ω/Θ.

use crate::*;
use crate::so3;

const EPS: f64 = 1e-10;

// =========================================================================
// SO(3) Jacobians (3×3)
// =========================================================================

/// 1-x factor: Θ/(2·tan(Θ/2)), stable near Θ=0.
#[inline]
fn one_minus_x(theta: f64) -> f64 {
    if theta < EPS {
        1.0 - theta * theta / 12.0
    } else {
        let half = 0.5 * theta;
        half / half.tan()
    }
}

/// Right Rodrigues Jacobian J_ωr(Ω) ∈ ℝ³ˣ³.
///
/// J_ωr(Ω) = (1-x)I + x·r·rᵀ + ½[Ω]×
///
/// Paper Eq. (RodweinormanJ). Properties:
/// - J_ωr(0) = I
/// - det(J_ωr) = (1-x)² + Θ²/4 > 0 for all Θ ≤ π
pub fn j_omega_right(omega: &Vec3) -> Mat3 {
    let theta = norm3(omega);
    if theta < EPS {
        let h = so3::hat(omega);
        return add_mat3(&I3, &scale_mat3(0.5, &h));
    }
    let omx = one_minus_x(theta);
    let x = 1.0 - omx;
    let r = scale3(1.0 / theta, omega);
    add_mat3(
        &add_mat3(&scale_mat3(omx, &I3), &scale_mat3(x, &outer3(&r, &r))),
        &scale_mat3(0.5, &so3::hat(omega)),
    )
}

/// Left Rodrigues Jacobian J_ωl(Ω) = J_ωr(-Ω) = S⁻¹(Ω).
pub fn j_omega_left(omega: &Vec3) -> Mat3 {
    let theta = norm3(omega);
    if theta < EPS {
        let h = so3::hat(omega);
        return sub_mat3(&I3, &scale_mat3(0.5, &h));
    }
    let omx = one_minus_x(theta);
    let x = 1.0 - omx;
    let r = scale3(1.0 / theta, omega);
    add_mat3(
        &add_mat3(&scale_mat3(omx, &I3), &scale_mat3(x, &outer3(&r, &r))),
        &scale_mat3(-0.5, &so3::hat(omega)),
    )
}

/// Inverse of right Rodrigues Jacobian (via 3×3 matrix inverse).
pub fn j_omega_right_inv(omega: &Vec3) -> Mat3 {
    inv3(&j_omega_right(omega))
}

// =========================================================================
// SE(3) coupling Jacobian
// =========================================================================

/// Coupling Jacobian J_t(Ω, t) ∈ ℝ³ˣ³.
///
/// Paper Eqs. (couplingJacobian), (dSinvTdOmega):
///   J_t = ∂[S⁻¹(Ω)·T]/∂Ω · J_ωr(Ω),  T = S(Ω)·t
///
/// The derivative ∂[S⁻¹T]/∂Ω (with T the physical translation, held
/// fixed) is computed numerically via central finite differences.
/// This avoids a known sign/variable error in the analytic formula
/// from the original tech report (which used t where T was required
/// in the skew-symmetric term). The numerical derivative is accurate
/// to ~10⁻¹⁰ and the function is not in any inner loop.
pub fn j_coupling(omega: &Vec3, t_exp: &Vec3) -> Mat3 {
    let theta = norm3(omega);
    if theta < EPS {
        // At Ω=0: S=I, T=t, d[S⁻¹T]/dΩ = ½[T]× = ½[t]×
        // J_wr = I, so J_t = ½[t]×
        return scale_mat3(0.5, &so3::hat(t_exp));
    }

    // Physical translation T = S(Ω)·t (held fixed during differentiation)
    let big_t = mv3(&so3::s_matrix(omega), t_exp);

    // d[S⁻¹(Ω)·T]/dΩ by central finite differences
    let h = 1e-8;
    let mut ds_inv_t = [[0.0; 3]; 3];
    for j in 0..3 {
        let mut wp = *omega; wp[j] += h;
        let mut wm = *omega; wm[j] -= h;
        let vp = mv3(&so3::s_inv(&wp), &big_t);
        let vm = mv3(&so3::s_inv(&wm), &big_t);
        for i in 0..3 {
            ds_inv_t[i][j] = (vp[i] - vm[i]) / (2.0 * h);
        }
    }

    // J_t = d[S⁻¹T]/dΩ · J_ωr(Ω)
    mm3(&ds_inv_t, &j_omega_right(omega))
}

// =========================================================================
// Full SE(3) Jacobian (6×6)
// =========================================================================

/// Assemble 6×6 from four 3×3 blocks: [[A,B],[C,D]].
fn blocks_6x6(a: &Mat3, b: &Mat3, c: &Mat3, d: &Mat3) -> Mat6 {
    let mut m = [[0.0f64; 6]; 6];
    for i in 0..3 { for j in 0..3 {
        m[i][j]     = a[i][j];
        m[i][j+3]   = b[i][j];
        m[i+3][j]   = c[i][j];
        m[i+3][j+3] = d[i][j];
    }}
    m
}

/// Full SE(3) right Jacobian J_gr ∈ ℝ⁶ˣ⁶.
///
///   J_gr = [ J_ωr(ω)       0     ]
///          [ J_t(ω,t)   J_ωr(ω)  ]
pub fn se3_right_jacobian(omega: &Vec3, t_exp: &Vec3) -> Mat6 {
    let jwr = j_omega_right(omega);
    let jt = j_coupling(omega, t_exp);
    blocks_6x6(&jwr, &Z3, &jt, &jwr)
}

/// Inverse of SE(3) right Jacobian.
///
///   J_gr⁻¹ = [ J_ωr⁻¹              0      ]
///            [ -J_ωr⁻¹·J_t·J_ωr⁻¹  J_ωr⁻¹ ]
pub fn se3_right_jacobian_inv(omega: &Vec3, t_exp: &Vec3) -> Mat6 {
    let jwri = j_omega_right_inv(omega);
    let jt = j_coupling(omega, t_exp);
    let ll = scale_mat3(-1.0, &mm3(&jwri, &mm3(&jt, &jwri)));
    blocks_6x6(&jwri, &Z3, &ll, &jwri)
}

/// Adjoint Ad(f) ∈ ℝ⁶ˣ⁶ for f = (R, T).
///
///   Ad(f) = [ R      0  ]
///           [ [T]×R  R  ]
pub fn adjoint(rot: &Mat3, trans: &Vec3) -> Mat6 {
    let tcr = mm3(&so3::hat(trans), rot);
    blocks_6x6(rot, &Z3, &tcr, rot)
}

// =========================================================================
// Tests
// =========================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use crate::se3::Pose;

    fn approx_mat3(a: &Mat3, b: &Mat3, tol: f64) -> bool {
        for i in 0..3 { for j in 0..3 {
            if (a[i][j] - b[i][j]).abs() > tol { return false; }
        }}
        true
    }
    fn max6(a: &Mat6, b: &Mat6) -> f64 {
        let mut m = 0.0f64;
        for i in 0..6 { for j in 0..6 { m = m.max((a[i][j]-b[i][j]).abs()); }}
        m
    }

    // SO(3) Jacobian tests
    #[test] fn j_right_zero_is_i()   { assert!(approx_mat3(&j_omega_right(&[0.;3]), &I3, 1e-10)); }
    #[test] fn j_left_eq_j_right_neg() {
        let w = [0.5,-0.3,0.8];
        assert!(approx_mat3(&j_omega_left(&w), &j_omega_right(&scale3(-1.0,&w)), 1e-12));
    }
    #[test] fn j_right_det() {
        let w = [0.7,-0.4,0.5]; let th = norm3(&w); let omx = one_minus_x(th);
        assert!((det3(&j_omega_right(&w)) - (omx*omx + th*th/4.0)).abs() < 1e-10);
    }
    #[test] fn j_inv_roundtrip() {
        let w = [0.5,-0.3,0.8];
        assert!(approx_mat3(&mm3(&j_omega_right(&w), &j_omega_right_inv(&w)), &I3, 1e-10));
    }

    // S⁻¹ identities
    #[test] fn sinv_eq_jl() {
        for w in &[[0.01,0.02,-0.01],[0.5,-0.3,0.8],[2.0,-1.0,0.5]] {
            assert!(approx_mat3(&so3::s_inv(w), &j_omega_left(w), 1e-10), "at {:?}", w);
        }
    }
    #[test] fn sinv_r_eq_jr() {
        for w in &[[0.3,-0.5,0.7],[1.5,0.0,0.0],[0.01,0.0,0.0]] {
            let lhs = mm3(&so3::s_inv(w), &so3::exp(w));
            assert!(approx_mat3(&lhs, &j_omega_right(w), 1e-10), "at {:?}", w);
        }
    }

    // SE(3) Jacobian structure
    #[test] fn se3_id_is_i6() { assert!(max6(&se3_right_jacobian(&[0.;3],&[0.;3]), &I6) < 1e-10); }
    #[test] fn se3_upper_right_zero() {
        let j = se3_right_jacobian(&[0.4,-0.2,0.6], &[0.5,1.0,-0.3]);
        for i in 0..3 { for k in 3..6 { assert!(j[i][k].abs() < 1e-15); }}
    }
    #[test] fn se3_diag_blocks() {
        let w = [0.4,-0.2,0.6]; let t = [0.5,1.0,-0.3];
        let j = se3_right_jacobian(&w, &t); let jwr = j_omega_right(&w);
        for i in 0..3 { for k in 0..3 {
            assert!((j[i][k] - jwr[i][k]).abs() < 1e-12);
            assert!((j[i+3][k+3] - jwr[i][k]).abs() < 1e-12);
        }}
    }
    #[test] fn se3_inv_roundtrip() {
        let w = [0.3,-0.5,0.7]; let t = [1.0,-2.0,0.5];
        let j = se3_right_jacobian(&w,&t); let ji = se3_right_jacobian_inv(&w,&t);
        let mut p = [[0.0f64;6];6];
        for i in 0..6{for k in 0..6{for m in 0..6{p[i][k]+=j[i][m]*ji[m][k];}}}
        assert!(max6(&p, &I6) < 1e-9, "err={:.2e}", max6(&p, &I6));
    }

    // Finite-difference: J_ωr
    #[test] fn j_right_fd() {
        let w = [0.5,-0.3,0.8]; let j = j_omega_right(&w); let eps = 1e-7;
        let r0 = so3::exp(&w);
        for ax in 0..3 {
            let mut e = [0.0;3]; e[ax] = eps;
            let wc = so3::log(&mm3(&r0, &so3::exp(&e)));
            for row in 0..3 {
                let fd = (wc[row]-w[row])/eps;
                assert!((fd - j[row][ax]).abs() < 1e-5,
                    "[{},{}] fd={:.6} an={:.6}", row, ax, fd, j[row][ax]);
            }
        }
    }

    // Finite-difference: full SE(3) Jacobian
    // Uses central differences and the same test point as Mathematica Part 9-10.
    #[test]
    fn se3_jac_fd() {
        // Same test point as Mathematica verification
        let w = [0.5, -0.3, 0.7];
        let t = [1.0, -0.5, 0.3];
        let xi: Vec6 = [w[0],w[1],w[2],t[0],t[1],t[2]];
        let f0 = Pose::exp(&xi);
        let j = se3_right_jacobian(&w, &t);
        let eps = 1e-7;

        // Central FD of log(exp(xi) . exp(d))
        let mut jfd = [[0.0f64; 6]; 6];
        for ax in 0..6 {
            let mut dp = [0.0f64; 6]; dp[ax] = eps;
            let mut dm = [0.0f64; 6]; dm[ax] = -eps;
            let fp = f0.compose(&Pose::exp(&dp));
            let fm = f0.compose(&Pose::exp(&dm));
            let xip = fp.log();
            let xim = fm.log();
            for row in 0..6 {
                jfd[row][ax] = (xip[row] - xim[row]) / (2.0 * eps);
            }
        }

        // Also compute J_t two ways:
        // Way 1: d[S^{-1}(Omega) T]/dOmega by FD, T fixed
        let big_t = mv3(&so3::s_matrix(&w), &t);  // physical translation
        let mut ds_inv_t_fd = [[0.0; 3]; 3];
        for ax in 0..3 {
            let mut wp = w; wp[ax] += eps;
            let mut wm = w; wm[ax] -= eps;
            let sp = mv3(&so3::s_inv(&wp), &big_t);
            let sm = mv3(&so3::s_inv(&wm), &big_t);
            for row in 0..3 {
                ds_inv_t_fd[row][ax] = (sp[row] - sm[row]) / (2.0 * eps);
            }
        }
        // Way 2: J_t = ds_inv_t * J_wr (the paper formula)
        let jwr = j_omega_right(&w);
        let jt_from_dsdt = mm3(&ds_inv_t_fd, &jwr);

        // Way 3: the analytic j_coupling
        let jt_analytic = j_coupling(&w, &t);

        // Way 4: extract from full SE(3) FD
        let mut jt_se3fd = [[0.0; 3]; 3];
        for i in 0..3 { for j in 0..3 {
            jt_se3fd[i][j] = jfd[i+3][j];
        }}

        // Print diagnostics
        eprintln!("\n=== SE(3) Jacobian Diagnostics ===");
        eprintln!("Test point: omega={:?}, t={:?}", w, t);

        eprintln!("\nJ_t from SE(3) FD (ground truth):");
        for i in 0..3 {
            eprintln!("  [{:.6}, {:.6}, {:.6}]", jt_se3fd[i][0], jt_se3fd[i][1], jt_se3fd[i][2]);
        }
        eprintln!("\nJ_t from d[S^{{-1}}T]/dOmega . J_wr:");
        for i in 0..3 {
            eprintln!("  [{:.6}, {:.6}, {:.6}]", jt_from_dsdt[i][0], jt_from_dsdt[i][1], jt_from_dsdt[i][2]);
        }
        eprintln!("\nJ_t analytic (j_coupling):");
        for i in 0..3 {
            eprintln!("  [{:.6}, {:.6}, {:.6}]", jt_analytic[i][0], jt_analytic[i][1], jt_analytic[i][2]);
        }

        // Errors
        let mut max_dsdt_err = 0.0f64;
        let mut max_analytic_err = 0.0f64;
        for i in 0..3 { for j in 0..3 {
            max_dsdt_err = max_dsdt_err.max((jt_from_dsdt[i][j] - jt_se3fd[i][j]).abs());
            max_analytic_err = max_analytic_err.max((jt_analytic[i][j] - jt_se3fd[i][j]).abs());
        }}
        eprintln!("\nJ_t: d[S^{{-1}}T]/dOmega . J_wr vs SE3 FD: max err = {:.2e}", max_dsdt_err);
        eprintln!("J_t: analytic vs SE3 FD: max err = {:.2e}", max_analytic_err);

        // Full 6x6 check
        let mut max_full_err = 0.0f64;
        for i in 0..6 { for k in 0..6 {
            let err = (j[i][k] - jfd[i][k]).abs();
            max_full_err = max_full_err.max(err);
            if err > 1e-4 {
                eprintln!("MISMATCH [{},{}]: analytic={:.6} fd={:.6} err={:.2e}",
                    i, k, j[i][k], jfd[i][k], err);
            }
        }}
        eprintln!("\nFull 6x6: max err = {:.2e}", max_full_err);

        // Assert
        assert!(max_full_err < 5e-4,
            "SE(3) Jacobian FD mismatch: max err = {:.2e}", max_full_err);
    }
}
