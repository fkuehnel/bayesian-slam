//! # Saddlepoint marginalization for projective observations
//!
//! Paper §IV and Appendix E. Corrected formula (Eq. saddlepointcorrection):
//!
//! I/I_Laplace = 1 + c₁ + O(f'''⁴)
//!
//! c₁ = (1/12)·A + (1/8)·B − (1/8)·Q₄
//!
//! where:
//!   A = cross-contraction: f'''_{ijk} f'''_{lmn} H⁻¹_{il} H⁻¹_{jm} H⁻¹_{kn}
//!   B = trace-contraction: v^T H⁻¹ v,  v_k = f'''_{ijk} H⁻¹_{ij}
//!   Q₄ = quartic: f''''_{ijkl} H⁻¹_{ij} H⁻¹_{kl}
//!
//! ## Validity guard
//!
//! The expansion is valid when |c₁| ≪ 1. When |c₁| > MAX_CORRECTION
//! (default 0.5), we fall back to Laplace and signal `SaddlepointStatus::Invalid`.

use crate::*;
use crate::projective;

/// Maximum |c₁| before falling back to Laplace.
/// 0.5 corresponds to σ_depth/depth ≈ 0.3.
pub const MAX_CORRECTION: f64 = 0.5;

/// Status of the saddlepoint correction for a single landmark.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SaddlepointStatus {
    /// Correction applied, |c₁| < MAX_CORRECTION.
    Valid,
    /// |c₁| exceeded threshold; Laplace returned instead.
    Invalid(f64),
    /// Hessian was singular; Laplace returned.
    SingularHessian,
}

// =========================================================================
// Landmark optimization (Step 1)
// =========================================================================

/// Result of optimizing a single landmark given the camera pose.
#[derive(Debug, Clone)]
pub struct LandmarkOptResult {
    pub x_opt: Vec3,
    pub xp_opt: Vec3,
    pub hessian: Mat3,
    pub nll_opt: f64,
    pub log_det_h: f64,
}

/// Optimize a single landmark via Gauss-Newton.
pub fn optimize_landmark(
    rot: &Mat3, trans: &Vec3,
    z: &[f64; 2], sigma_zz_inv: &[[f64; 2]; 2],
    x_prior: &Vec3, sigma_xx_inv: &Mat3,
    max_iter: usize,
) -> LandmarkOptResult {
    let mut x = *x_prior;

    for _ in 0..max_iter {
        let xp = projective::transform_point(rot, trans, &x);
        if xp[2] < 1e-6 { break; }

        let pi = projective::project(&xp);
        let e = [z[0] - pi[0], z[1] - pi[1]];
        let p = projective::project_jacobian(&xp);

        let mut pr = [[0.0; 3]; 2];
        for i in 0..2 { for j in 0..3 {
            pr[i][j] = p[i][0]*rot[0][j] + p[i][1]*rot[1][j] + p[i][2]*rot[2][j];
        }}

        let mut h = *sigma_xx_inv;
        for i in 0..3 { for j in 0..3 {
            for m in 0..2 { for n in 0..2 {
                h[i][j] += pr[m][i] * sigma_zz_inv[m][n] * pr[n][j];
            }}
        }}

        let mut se = [0.0; 2];
        for i in 0..2 { se[i] = sigma_zz_inv[i][0]*e[0] + sigma_zz_inv[i][1]*e[1]; }
        let dx = sub3(&x, x_prior);
        let mut grad = mv3(sigma_xx_inv, &dx);
        for i in 0..3 { grad[i] -= pr[0][i]*se[0] + pr[1][i]*se[1]; }

        let h_inv = match inv3_safe(&h) { Some(hi) => hi, None => break };
        let delta = mv3(&h_inv, &scale3(-1.0, &grad));
        x = add3(&x, &delta);
        if norm3(&delta) < 1e-10 { break; }
    }

    let xp = projective::transform_point(rot, trans, &x);
    let nll = if xp[2] > 1e-6 {
        projective::neg_log_likelihood(z, &xp, sigma_zz_inv)
            + 0.5 * { let dx = sub3(&x, x_prior); dot3(&dx, &mv3(sigma_xx_inv, &dx)) }
    } else { f64::INFINITY };

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
    LandmarkOptResult {
        x_opt: x, xp_opt: xp, hessian, nll_opt: nll,
        log_det_h: if d > 0.0 { d.ln() } else { f64::NEG_INFINITY },
    }
}

// =========================================================================
// Quartic contraction Q₄ by FD
// =========================================================================

/// Compute Q₄ = Σ f''''_{abcd} H⁻¹_{ab} H⁻¹_{cd} via FD of the Hessian trace.
///
/// Defines T(x) = tr(H_NLL(x) · H⁻¹_opt), then Q₄ = tr(H⁻¹_opt · Hess(T)).
fn quartic_contraction(
    rot: &Mat3, trans: &Vec3, x_opt: &Vec3,
    sigma_zz_inv: &[[f64; 2]; 2], sigma_xx_inv: &Mat3,
    h_inv: &Mat3,
) -> f64 {
    let h = 1e-5;

    // Hessian of NLL at arbitrary world-frame point x
    let hess_at = |x: &Vec3| -> Mat3 {
        let xp = projective::transform_point(rot, trans, x);
        if xp[2] < 1e-6 { return [[1e30; 3]; 3]; }
        let p = projective::project_jacobian(&xp);
        let mut pr = [[0.0; 3]; 2];
        for i in 0..2 { for j in 0..3 {
            pr[i][j] = p[i][0]*rot[0][j] + p[i][1]*rot[1][j] + p[i][2]*rot[2][j];
        }}
        let mut hess = *sigma_xx_inv;
        for i in 0..3 { for j in 0..3 {
            for m in 0..2 { for n in 0..2 {
                hess[i][j] += pr[m][i] * sigma_zz_inv[m][n] * pr[n][j];
            }}
        }}
        hess
    };

    // T(x) = tr(H_NLL(x) · H⁻¹_opt)
    let trace_func = |x: &Vec3| -> f64 {
        let hx = hess_at(x);
        let mut t = 0.0;
        for a in 0..3 { for b in 0..3 { t += hx[a][b] * h_inv[a][b]; }}
        t
    };

    // Hess(T) by central FD, contract with H⁻¹
    let t0 = trace_func(x_opt);
    let mut q4 = 0.0;
    for c in 0..3 { for d in 0..3 {
        let d2t = if c == d {
            let mut xp = *x_opt; xp[c] += h;
            let mut xm = *x_opt; xm[c] -= h;
            (trace_func(&xp) - 2.0*t0 + trace_func(&xm)) / (h*h)
        } else {
            let mut xpp = *x_opt; xpp[c] += h; xpp[d] += h;
            let mut xpm = *x_opt; xpm[c] += h; xpm[d] -= h;
            let mut xmp = *x_opt; xmp[c] -= h; xmp[d] += h;
            let mut xmm = *x_opt; xmm[c] -= h; xmm[d] -= h;
            (trace_func(&xpp) - trace_func(&xpm) - trace_func(&xmp) + trace_func(&xmm))
                / (4.0*h*h)
        };
        q4 += h_inv[c][d] * d2t;
    }}
    q4
}

// =========================================================================
// Saddlepoint correction (Step 2) — corrected formula
// =========================================================================

/// Detailed result of the saddlepoint correction.
#[derive(Debug, Clone)]
pub struct SaddlepointResult {
    /// c₁ = (1/12)A + (1/8)B − (1/8)Q₄.
    pub c1: f64,
    /// Cross-contraction A.
    pub term_a: f64,
    /// Trace-contraction B.
    pub term_b: f64,
    /// Quartic contraction Q₄.
    pub term_q4: f64,
    /// Validity status.
    pub status: SaddlepointStatus,
}

/// Compute c₁ from f''', Q₄, and H⁻¹.
///
/// Paper Eq. (saddlepointcorrection): c₁ = (1/12)A + (1/8)B − (1/8)Q₄
pub fn saddlepoint_correction_full(
    h_inv: &Mat3,
    f3: &[[[f64; 3]; 3]; 3],
    q4: f64,
) -> SaddlepointResult {
    // A: cross-contraction
    let mut lambda = [[[0.0f64; 3]; 3]; 3];
    for a in 0..3 { for b in 0..3 { for c in 0..3 {
        let mut val = 0.0;
        for ap in 0..3 { for bp in 0..3 { for cp in 0..3 {
            val += f3[ap][bp][cp] * h_inv[a][ap] * h_inv[b][bp] * h_inv[c][cp];
        }}}
        lambda[a][b][c] = val;
    }}}
    let mut term_a = 0.0;
    for a in 0..3 { for b in 0..3 { for c in 0..3 {
        term_a += f3[a][b][c] * lambda[a][b][c];
    }}}

    // B: trace-contraction  v_k = Σ_{ij} f'''_{ijk} H⁻¹_{ij}
    let mut v = [0.0; 3];
    for k in 0..3 { for i in 0..3 { for j in 0..3 {
        v[k] += f3[i][j][k] * h_inv[i][j];
    }}}
    let term_b = dot3(&v, &mv3(h_inv, &v));

    let c1 = term_a / 12.0 + term_b / 8.0 - q4 / 8.0;

    let status = if c1.abs() <= MAX_CORRECTION {
        SaddlepointStatus::Valid
    } else {
        SaddlepointStatus::Invalid(c1)
    };

    SaddlepointResult { c1, term_a, term_b, term_q4: q4, status }
}

/// Simplified: cubic terms only (Q₄ = 0).
pub fn saddlepoint_correction(
    h_inv: &Mat3,
    f3: &[[[f64; 3]; 3]; 3],
) -> SaddlepointResult {
    saddlepoint_correction_full(h_inv, f3, 0.0)
}

// =========================================================================
// Marginal log-posterior
// =========================================================================

/// Saddlepoint-corrected marginal log-posterior for one landmark.
///
/// Falls back to Laplace if |c₁| > MAX_CORRECTION.
pub fn landmark_marginal_log_posterior(
    opt: &LandmarkOptResult,
    rot: &Mat3, trans: &Vec3,
    sigma_zz_inv: &[[f64; 2]; 2],
    sigma_xx_inv: &Mat3,
) -> (f64, SaddlepointResult) {
    let laplace = -opt.nll_opt - 0.5 * opt.log_det_h
        + 1.5 * (2.0 * std::f64::consts::PI).ln();

    let h_inv = match inv3_safe(&opt.hessian) {
        Some(hi) => hi,
        None => return (laplace, SaddlepointResult {
            c1: 0.0, term_a: 0.0, term_b: 0.0, term_q4: 0.0,
            status: SaddlepointStatus::SingularHessian,
        }),
    };

    let f3 = projective::third_cumulants(rot, &opt.xp_opt, sigma_zz_inv);
    let q4 = quartic_contraction(rot, trans, &opt.x_opt, sigma_zz_inv, sigma_xx_inv, &h_inv);
    let result = saddlepoint_correction_full(&h_inv, &f3, q4);

    let log_post = match result.status {
        SaddlepointStatus::Valid => laplace + (1.0 + result.c1).ln(),
        _ => laplace,
    };
    (log_post, result)
}

/// Laplace (uncorrected) marginal log-posterior.
pub fn landmark_marginal_log_posterior_laplace(opt: &LandmarkOptResult) -> f64 {
    -opt.nll_opt - 0.5 * opt.log_det_h
        + 1.5 * (2.0 * std::f64::consts::PI).ln()
}

// =========================================================================
// Full marginalized MAP
// =========================================================================

#[derive(Debug, Clone)]
pub struct MarginalizedMapResult {
    pub log_posterior_sp: f64,
    pub log_posterior_laplace: f64,
    pub landmark_opts: Vec<LandmarkOptResult>,
    pub sp_results: Vec<SaddlepointResult>,
    pub n_valid: usize,
    pub n_invalid: usize,
}

pub fn evaluate_marginalized_map(
    rot: &Mat3, trans: &Vec3,
    observations: &[([f64; 2], Vec3, Mat3)],
    sigma_zz_inv: &[[f64; 2]; 2],
    max_gn_iter: usize,
) -> MarginalizedMapResult {
    let mut log_sp = 0.0;
    let mut log_laplace = 0.0;
    let mut opts = Vec::with_capacity(observations.len());
    let mut sp_results = Vec::with_capacity(observations.len());
    let (mut n_valid, mut n_invalid) = (0, 0);

    for (z, x_prior, sigma_xx_inv) in observations {
        let opt = optimize_landmark(rot, trans, z, sigma_zz_inv, x_prior, sigma_xx_inv, max_gn_iter);
        let laplace_c = landmark_marginal_log_posterior_laplace(&opt);
        let (sp_c, sp_r) = landmark_marginal_log_posterior(&opt, rot, trans, sigma_zz_inv, sigma_xx_inv);

        match sp_r.status { SaddlepointStatus::Valid => n_valid += 1, _ => n_invalid += 1 }
        log_laplace += laplace_c;
        log_sp += sp_c;
        sp_results.push(sp_r);
        opts.push(opt);
    }

    MarginalizedMapResult {
        log_posterior_sp: log_sp, log_posterior_laplace: log_laplace,
        landmark_opts: opts, sp_results, n_valid, n_invalid,
    }
}

fn inv3_safe(m: &Mat3) -> Option<Mat3> {
    let d = det3(m);
    if d.abs() < 1e-30 { None } else { Some(inv3(m)) }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::so3;

    fn iso2(s: f64) -> [[f64; 2]; 2] { let v = 1.0/(s*s); [[v,0.0],[0.0,v]] }
    fn iso3(s: f64) -> Mat3 { let v = 1.0/(s*s); [[v,0.0,0.0],[0.0,v,0.0],[0.0,0.0,v]] }

    #[test]
    fn landmark_opt_recovers_projection() {
        let rot = I3; let trans = [0.0; 3];
        let x_true = [0.5, -0.3, 5.0];
        let xp = projective::transform_point(&rot, &trans, &x_true);
        let z = projective::project(&xp);
        let opt = optimize_landmark(&rot, &trans, &z, &iso2(0.01), &[0.6,-0.2,5.5], &iso3(10.0), 20);
        let pi = projective::project(&projective::transform_point(&rot, &trans, &opt.x_opt));
        assert!((pi[0] - z[0]).abs() < 0.001);
        assert!((pi[1] - z[1]).abs() < 0.001);
    }

    #[test]
    fn correction_small_for_distant_well_constrained() {
        let rot = I3; let trans = [0.0; 3];
        let x = [0.5, -0.3, 10.0];
        let z = projective::project(&x);
        let opt = optimize_landmark(&rot, &trans, &z, &iso2(0.01), &x, &iso3(1.0), 10);
        let (_, r) = landmark_marginal_log_posterior(&opt, &rot, &trans, &iso2(0.01), &iso3(1.0));
        assert_eq!(r.status, SaddlepointStatus::Valid);
        assert!(r.c1.abs() < 0.1, "c₁ should be small: {:.4}", r.c1);
        eprintln!("Well-constrained: c₁={:.6}, A={:.4}, B={:.4}, Q₄={:.4}", r.c1, r.term_a, r.term_b, r.term_q4);
    }

    #[test]
    fn correction_invalid_for_weak_prior() {
        let rot = I3; let trans = [0.0; 3];
        let x = [0.5, -0.3, 3.0];
        let z = projective::project(&x);
        let opt = optimize_landmark(&rot, &trans, &z, &iso2(0.01), &x, &iso3(10.0), 10);
        let (_, r) = landmark_marginal_log_posterior(&opt, &rot, &trans, &iso2(0.01), &iso3(10.0));
        assert!(matches!(r.status, SaddlepointStatus::Invalid(_)),
            "Weak prior at close depth should be invalid: c₁={:.4}", r.c1);
    }

    #[test]
    fn correction_larger_for_close_point() {
        let rot = I3; let trans = [0.0; 3];
        let far = optimize_landmark(&rot, &trans, &projective::project(&[0.1,0.1,20.0]),
            &iso2(0.01), &[0.1,0.1,20.0], &iso3(1.0), 10);
        let close = optimize_landmark(&rot, &trans, &projective::project(&[0.1,0.1,5.0]),
            &iso2(0.01), &[0.1,0.1,5.0], &iso3(1.0), 10);
        let (_, rf) = landmark_marginal_log_posterior(&far, &rot, &trans, &iso2(0.01), &iso3(1.0));
        let (_, rc) = landmark_marginal_log_posterior(&close, &rot, &trans, &iso2(0.01), &iso3(1.0));
        assert!(rc.c1.abs() > rf.c1.abs(), "Close |c₁|={:.6} > far |c₁|={:.6}", rc.c1.abs(), rf.c1.abs());
    }

    #[test]
    fn marginalized_map_reports_validity() {
        let rot = so3::exp(&[0.1, -0.05, 0.2]);
        let trans = [0.0, 0.0, -5.0];
        let lms: Vec<([f64;3], f64)> = vec![([1.0,0.5,8.0], 1.0), ([-0.5,1.0,7.0], 1.0), ([0.3,-0.8,8.0], 10.0)];
        let obs: Vec<_> = lms.iter().map(|(x, sx)| {
            let xp = projective::transform_point(&rot, &trans, x);
            (projective::project(&xp), *x, iso3(*sx))
        }).collect();
        let r = evaluate_marginalized_map(&rot, &trans, &obs, &iso2(0.01), 10);
        assert_eq!(r.sp_results.len(), 3);
        assert!(r.log_posterior_sp.is_finite());
        eprintln!("Valid: {}, Invalid: {}", r.n_valid, r.n_invalid);
    }

    #[test]
    fn laplace_and_sp_agree_distant() {
        let rot = I3; let trans = [0.0; 3];
        let x = [0.01, 0.01, 100.0];
        let z = projective::project(&x);
        let opt = optimize_landmark(&rot, &trans, &z, &iso2(0.01), &x, &iso3(0.01), 10);
        let lap = landmark_marginal_log_posterior_laplace(&opt);
        let (sp, r) = landmark_marginal_log_posterior(&opt, &rot, &trans, &iso2(0.01), &iso3(0.01));
        assert_eq!(r.status, SaddlepointStatus::Valid);
        assert!((lap - sp).abs() < 0.01, "Lap={:.6} SP={:.6}", lap, sp);
    }

    #[test]
    fn quartic_contraction_nonzero() {
        let rot = I3; let trans = [0.0; 3];
        let x = [0.5, -0.3, 10.0];
        let z = projective::project(&x);
        let opt = optimize_landmark(&rot, &trans, &z, &iso2(0.01), &x, &iso3(1.0), 10);
        let h_inv = inv3(&opt.hessian);
        let q4 = quartic_contraction(&rot, &trans, &opt.x_opt, &iso2(0.01), &iso3(1.0), &h_inv);
        assert!(q4.is_finite() && q4.abs() > 1e-10, "Q₄={:.2e}", q4);
    }
}
