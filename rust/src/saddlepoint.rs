//! # Saddlepoint marginalization for projective observations
//!
//! Paper §IV and Appendix E. Given the joint posterior p(f, {x} | z),
//! we marginalize over 3D landmarks x to obtain the marginal p(f | z)
//! using a saddlepoint approximation that captures non-Gaussian effects
//! from the projective observation model.
//!
//! ## Two-step scheme
//!
//! Step 1: For fixed camera pose f, find optimal landmarks x_opt(f)
//! Step 2: Saddlepoint-approximate ∫ p(f, x | z) dx for each landmark
//!
//! ## Key formula (Eq. saddlepointfull)
//!
//! ∫ e^{ℓ(x)} dx ≈ e^{ℓ(x_opt)} (2π)^{3/2} |H|^{-1/2} (1 + δ^SP)
//!
//! where δ^SP is the saddlepoint correction involving third cumulants.

use crate::*;
use crate::projective;

// =========================================================================
// Landmark optimization (Step 1)
// =========================================================================

/// Result of optimizing a single landmark given the camera pose.
#[derive(Debug, Clone)]
pub struct LandmarkOptResult {
    /// Optimal landmark position in world coordinates.
    pub x_opt: Vec3,
    /// Point in camera frame: x' = f ⋆ x_opt.
    pub xp_opt: Vec3,
    /// Hessian of neg-log-posterior w.r.t. x at x_opt (3×3, positive definite).
    pub hessian: Mat3,
    /// Neg-log-posterior value at x_opt.
    pub nll_opt: f64,
    /// Log-determinant of Hessian: ln|H|.
    pub log_det_h: f64,
}

/// Optimize a single landmark's position given camera pose and observation.
///
/// Solves: x_opt = argmin_x [ ½(z - π(f⋆x))^T Σ_zz^{-1} (z - π(f⋆x))
///                           + ½(x - x_prior)^T Σ_xx^{-1} (x - x_prior) ]
///
/// Uses Gauss-Newton iteration since the projective model is nonlinear.
///
/// Arguments:
/// - `rot`, `trans`: camera pose (R, T) where x' = R·x + T
/// - `z`: 2D observation [u, v]
/// - `sigma_zz_inv`: 2×2 measurement precision
/// - `x_prior`: 3D prior mean
/// - `sigma_xx_inv`: 3×3 prior precision
/// - `max_iter`: maximum GN iterations
pub fn optimize_landmark(
    rot: &Mat3,
    trans: &Vec3,
    z: &[f64; 2],
    sigma_zz_inv: &[[f64; 2]; 2],
    x_prior: &Vec3,
    sigma_xx_inv: &Mat3,
    max_iter: usize,
) -> LandmarkOptResult {
    let mut x = *x_prior;

    for _iter in 0..max_iter {
        let xp = projective::transform_point(rot, trans, &x);

        // Check if point is in front of camera
        if xp[2] < 1e-6 { break; }

        let pi = projective::project(&xp);
        let e = [z[0] - pi[0], z[1] - pi[1]];
        let p = projective::project_jacobian(&xp);

        // Jacobian of π(R·x + T) w.r.t. x is P·R (2×3)
        let mut pr = [[0.0; 3]; 2];
        for i in 0..2 { for j in 0..3 {
            pr[i][j] = p[i][0]*rot[0][j] + p[i][1]*rot[1][j] + p[i][2]*rot[2][j];
        }}

        // Hessian: H = (PR)^T Σ^{-1} (PR) + Σ_xx^{-1}  (3×3)
        let mut h = *sigma_xx_inv;
        for i in 0..3 { for j in 0..3 {
            for m in 0..2 { for n in 0..2 {
                h[i][j] += pr[m][i] * sigma_zz_inv[m][n] * pr[n][j];
            }}
        }}

        // Gradient: g = -(PR)^T Σ^{-1} e + Σ_xx^{-1} (x - x_prior)
        let mut se = [0.0; 2];
        for i in 0..2 {
            se[i] = sigma_zz_inv[i][0] * e[0] + sigma_zz_inv[i][1] * e[1];
        }
        let dx_prior = sub3(&x, x_prior);
        let mut grad = mv3(sigma_xx_inv, &dx_prior);
        for i in 0..3 {
            grad[i] -= pr[0][i] * se[0] + pr[1][i] * se[1];
        }

        // Solve H · delta = -grad
        let h_inv = match inv3_safe(&h) {
            Some(hi) => hi,
            None => break,
        };
        let delta = mv3(&h_inv, &scale3(-1.0, &grad));

        // Update
        x = add3(&x, &delta);

        // Convergence check
        if norm3(&delta) < 1e-10 { break; }
    }

    // Final evaluation at x_opt
    let xp = projective::transform_point(rot, trans, &x);
    let nll = if xp[2] > 1e-6 {
        projective::neg_log_likelihood(z, &xp, sigma_zz_inv)
            + 0.5 * {
                let dx = sub3(&x, x_prior);
                let sdx = mv3(sigma_xx_inv, &dx);
                dot3(&dx, &sdx)
            }
    } else {
        f64::INFINITY
    };

    // Hessian at optimum
    let p = projective::project_jacobian(&xp);
    let mut pr = [[0.0; 3]; 2];
    for i in 0..2 { for j in 0..3 {
        pr[i][j] = p[i][0]*rot[0][j] + p[i][1]*rot[1][j] + p[i][2]*rot[2][j];
    }}
    let mut hessian = *sigma_xx_inv;
    for i in 0..3 { for j in 0..3 {
        for m in 0..2 { for n in 0..2 {
            hessian[i][j] += pr[m][i] * sigma_zz_inv[m][n] * pr[n][j];
        }}
    }}

    let d = det3(&hessian);
    let log_det = if d > 0.0 { d.ln() } else { f64::NEG_INFINITY };

    LandmarkOptResult {
        x_opt: x,
        xp_opt: xp,
        hessian,
        nll_opt: nll,
        log_det_h: log_det,
    }
}

// =========================================================================
// Saddlepoint correction (Step 2)
// =========================================================================

/// Saddlepoint correction for a single landmark marginalization.
///
/// Given the Hessian H and third cumulants κ_{abc} at the optimum,
/// computes the correction δ^SP to the Laplace approximation:
///
///   ∫ e^{ℓ(x)} dx ≈ e^{ℓ(x_opt)} (2π)^{3/2} |H|^{-1/2} (1 + δ^SP)
///
/// The correction involves the contraction:
///   δ^SP = (5/24) Σ_{abc} κ²_{abc} [H⁻³]_{abc}
///
/// where [H⁻³]_{abc} = Σ_{def} H⁻¹_{ad} H⁻¹_{be} H⁻¹_{cf}.
///
/// Paper Eq. (saddlepointfull).
pub fn saddlepoint_correction(
    h_inv: &Mat3,
    kappa: &[[[f64; 3]; 3]; 3],
) -> f64 {
    // Compute [H⁻³]_{abc} = Σ_{def} H⁻¹_{ad} H⁻¹_{be} H⁻¹_{cf}
    // and contract with κ²_{abc}

    // Term 1: (5/24) Σ_{abc} (Σ_{a'b'c'} κ_{a'b'c'} H⁻¹_{aa'} H⁻¹_{bb'} H⁻¹_{cc'})²
    // This is a simplified form. The full correction also has a fourth-cumulant
    // piece (1/8) Σ κ_{abcd} H⁻¹_{ab} H⁻¹_{cd}, but for the projective model
    // the third-cumulant term dominates.

    let mut correction = 0.0;

    // Contract κ with H⁻¹ to get λ_{abc} = Σ κ_{a'b'c'} H⁻¹_{aa'} H⁻¹_{bb'} H⁻¹_{cc'}
    let mut lambda = [[[0.0f64; 3]; 3]; 3];
    for a in 0..3 { for b in 0..3 { for c in 0..3 {
        let mut val = 0.0;
        for ap in 0..3 { for bp in 0..3 { for cp in 0..3 {
            val += kappa[ap][bp][cp] * h_inv[a][ap] * h_inv[b][bp] * h_inv[c][cp];
        }}}
        lambda[a][b][c] = val;
    }}}

    // (5/24) Σ_{abc} λ_{abc}²
    // Actually the correction is:
    // δ^SP = (5/24) Σ_{abc} κ_{abc} λ_{abc}
    for a in 0..3 { for b in 0..3 { for c in 0..3 {
        correction += kappa[a][b][c] * lambda[a][b][c];
    }}}
    correction *= 5.0 / 24.0;

    correction
}

/// Compute the saddlepoint-corrected marginal log-posterior contribution
/// from a single landmark.
///
/// Returns: ℓ(x_opt) - ½ ln|H| + ln(1 + δ^SP) + (3/2) ln(2π)
///
/// Paper Eq. (correctedmarginal).
pub fn landmark_marginal_log_posterior(
    opt: &LandmarkOptResult,
    rot: &Mat3,
    sigma_zz_inv: &[[f64; 2]; 2],
) -> f64 {
    let laplace = -opt.nll_opt - 0.5 * opt.log_det_h
        + 1.5 * (2.0 * std::f64::consts::PI).ln();

    // Saddlepoint correction
    let h_inv = match inv3_safe(&opt.hessian) {
        Some(hi) => hi,
        None => return laplace,
    };

    let kappa = projective::third_cumulants(rot, &opt.xp_opt, sigma_zz_inv);
    let delta_sp = saddlepoint_correction(&h_inv, &kappa);

    laplace + (1.0 + delta_sp).ln()
}

/// Compute the Laplace (standard) marginal log-posterior contribution
/// from a single landmark — no saddlepoint correction.
///
/// This is the baseline used in standard Schur complement approaches.
pub fn landmark_marginal_log_posterior_laplace(opt: &LandmarkOptResult) -> f64 {
    -opt.nll_opt - 0.5 * opt.log_det_h
        + 1.5 * (2.0 * std::f64::consts::PI).ln()
}

// =========================================================================
// Full marginalized MAP (combining Steps 1-3)
// =========================================================================

/// Result of the marginalized MAP evaluation for a set of landmarks.
#[derive(Debug, Clone)]
pub struct MarginalizedMapResult {
    /// Saddlepoint-corrected marginal log-posterior.
    pub log_posterior_sp: f64,
    /// Laplace (uncorrected) marginal log-posterior.
    pub log_posterior_laplace: f64,
    /// Individual landmark optimization results.
    pub landmark_opts: Vec<LandmarkOptResult>,
    /// Saddlepoint correction per landmark (δ^SP_j).
    pub corrections: Vec<f64>,
}

/// Evaluate the marginalized MAP objective for a given camera pose
/// over all landmarks.
///
/// This implements the two-step scheme:
///   Step 1: optimize each landmark x^j given f
///   Step 2: compute saddlepoint-corrected marginal
///
/// The result is the marginal log-posterior p(f | z) up to constants.
pub fn evaluate_marginalized_map(
    rot: &Mat3,
    trans: &Vec3,
    observations: &[([f64; 2], Vec3, Mat3)],  // (z, x_prior, Σ_xx^{-1}) per landmark
    sigma_zz_inv: &[[f64; 2]; 2],
    max_gn_iter: usize,
) -> MarginalizedMapResult {
    let mut log_sp = 0.0;
    let mut log_laplace = 0.0;
    let mut opts = Vec::with_capacity(observations.len());
    let mut corrections = Vec::with_capacity(observations.len());

    for (z, x_prior, sigma_xx_inv) in observations {
        let opt = optimize_landmark(
            rot, trans, z, sigma_zz_inv,
            x_prior, sigma_xx_inv, max_gn_iter,
        );

        let laplace_contrib = landmark_marginal_log_posterior_laplace(&opt);
        let sp_contrib = landmark_marginal_log_posterior(&opt, rot, sigma_zz_inv);

        let delta = sp_contrib - laplace_contrib;

        log_laplace += laplace_contrib;
        log_sp += sp_contrib;
        corrections.push(delta);
        opts.push(opt);
    }

    MarginalizedMapResult {
        log_posterior_sp: log_sp,
        log_posterior_laplace: log_laplace,
        landmark_opts: opts,
        corrections,
    }
}

// =========================================================================
// Helper: safe 3×3 inverse
// =========================================================================

fn inv3_safe(m: &Mat3) -> Option<Mat3> {
    let d = det3(m);
    if d.abs() < 1e-30 { return None; }
    Some(inv3(m))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn isotropic_sigma_inv(sigma: f64) -> [[f64; 2]; 2] {
        let s2 = 1.0 / (sigma * sigma);
        [[s2, 0.0], [0.0, s2]]
    }

    fn isotropic_sigma_inv_3d(sigma: f64) -> Mat3 {
        let s2 = 1.0 / (sigma * sigma);
        [[s2, 0.0, 0.0], [0.0, s2, 0.0], [0.0, 0.0, s2]]
    }

    #[test]
    fn landmark_opt_recovers_true_position() {
        // Camera looking down z-axis
        let rot = I3;
        let trans = [0.0, 0.0, 0.0];
        let x_true = [0.5, -0.3, 5.0];
        let xp = projective::transform_point(&rot, &trans, &x_true);
        let z = projective::project(&xp);

        let sigma_zz_inv = isotropic_sigma_inv(0.01);
        let sigma_xx_inv = isotropic_sigma_inv_3d(10.0);  // weak prior
        let x_prior = [0.6, -0.2, 5.5];  // close initial guess

        let opt = optimize_landmark(
            &rot, &trans, &z, &sigma_zz_inv,
            &x_prior, &sigma_xx_inv, 20,
        );

        // u,v are well-constrained but depth (z) is ill-conditioned
        // from a single camera view. Check lateral accuracy.
        let xp_opt = projective::transform_point(&rot, &trans, &opt.x_opt);
        let pi_opt = projective::project(&xp_opt);
        let pi_true = projective::project(&xp);
        assert!((pi_opt[0] - pi_true[0]).abs() < 0.001,
            "u: opt={:.6} true={:.6}", pi_opt[0], pi_true[0]);
        assert!((pi_opt[1] - pi_true[1]).abs() < 0.001,
            "v: opt={:.6} true={:.6}", pi_opt[1], pi_true[1]);
    }

    #[test]
    fn saddlepoint_correction_is_small_for_distant_point() {
        // Distant point → nearly linear projection → small correction
        let rot = I3;
        let trans = [0.0, 0.0, 0.0];
        let x_true = [0.1, 0.1, 20.0];  // Far away
        let xp = projective::transform_point(&rot, &trans, &x_true);
        let z = projective::project(&xp);

        let sigma_zz_inv = isotropic_sigma_inv(0.01);
        let sigma_xx_inv = isotropic_sigma_inv_3d(1.0);

        let opt = optimize_landmark(
            &rot, &trans, &z, &sigma_zz_inv,
            &x_true, &sigma_xx_inv, 10,
        );

        let h_inv = inv3(&opt.hessian);
        let kappa = projective::third_cumulants(&rot, &opt.xp_opt, &sigma_zz_inv);
        let delta = saddlepoint_correction(&h_inv, &kappa);

        assert!(delta.abs() < 0.1,
            "Distant point should have small SP correction: {:.4}", delta);
    }

    #[test]
    fn saddlepoint_correction_larger_for_close_point() {
        // Close point → more nonlinear projection → larger correction
        let rot = I3;
        let trans = [0.0, 0.0, 0.0];

        let x_far = [0.1, 0.1, 20.0];
        let x_close = [0.1, 0.1, 2.0];

        let sigma_zz_inv = isotropic_sigma_inv(0.01);
        let sigma_xx_inv = isotropic_sigma_inv_3d(1.0);

        let opt_far = optimize_landmark(
            &rot, &trans,
            &projective::project(&x_far),
            &sigma_zz_inv, &x_far, &sigma_xx_inv, 10,
        );
        let opt_close = optimize_landmark(
            &rot, &trans,
            &projective::project(&x_close),
            &sigma_zz_inv, &x_close, &sigma_xx_inv, 10,
        );

        let kf = projective::third_cumulants(&rot, &opt_far.xp_opt, &sigma_zz_inv);
        let kc = projective::third_cumulants(&rot, &opt_close.xp_opt, &sigma_zz_inv);
        let df = saddlepoint_correction(&inv3(&opt_far.hessian), &kf);
        let dc = saddlepoint_correction(&inv3(&opt_close.hessian), &kc);

        assert!(dc.abs() > df.abs(),
            "Close point correction |{:.6}| should exceed far |{:.6}|", dc, df);
    }

    #[test]
    fn marginalized_map_runs_multiple_landmarks() {
        let rot = so3::exp(&[0.1, -0.05, 0.2]);
        let trans = [0.0, 0.0, -5.0];
        let sigma_zz_inv = isotropic_sigma_inv(0.01);

        let landmarks = vec![
            [1.0, 0.5, 8.0],
            [-0.5, 1.0, 7.0],
            [0.3, -0.8, 10.0],
        ];

        let observations: Vec<_> = landmarks.iter().map(|x| {
            let xp = projective::transform_point(&rot, &trans, x);
            let z = projective::project(&xp);
            let sigma_xx_inv = isotropic_sigma_inv_3d(0.5);
            (z, *x, sigma_xx_inv)
        }).collect();

        let result = evaluate_marginalized_map(
            &rot, &trans, &observations, &sigma_zz_inv, 10,
        );

        assert_eq!(result.landmark_opts.len(), 3);
        assert_eq!(result.corrections.len(), 3);
        assert!(result.log_posterior_sp.is_finite());
        assert!(result.log_posterior_laplace.is_finite());
    }

    #[test]
    fn laplace_and_saddlepoint_agree_distant() {
        // For very distant points, SP correction vanishes
        let rot = I3;
        let trans = [0.0, 0.0, 0.0];
        let x = [0.01, 0.01, 100.0];
        let xp = projective::transform_point(&rot, &trans, &x);
        let z = projective::project(&xp);
        let sigma_zz_inv = isotropic_sigma_inv(0.01);
        let sigma_xx_inv = isotropic_sigma_inv_3d(0.01);

        let opt = optimize_landmark(
            &rot, &trans, &z, &sigma_zz_inv,
            &x, &sigma_xx_inv, 10,
        );

        let laplace = landmark_marginal_log_posterior_laplace(&opt);
        let sp = landmark_marginal_log_posterior(&opt, &rot, &sigma_zz_inv);

        assert!((laplace - sp).abs() < 0.01,
            "Distant: Laplace={:.6} SP={:.6}", laplace, sp);
    }
}
