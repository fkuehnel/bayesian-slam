//! # Pose Inference from Uncertain 3D Point Cloud
//!
//! Estimates camera pose by maximizing the marginalized posterior over
//! landmark positions. Compares Laplace vs saddlepoint-corrected objectives.
//!
//! This is the core use case: given 2D observations and uncertain 3D priors,
//! find the camera pose that best explains the data after integrating out
//! the landmarks.
//!
//! ## Key result
//!
//! At close range with weak landmark priors, the saddlepoint correction
//! shifts the optimal pose relative to Laplace. The shift is along the
//! depth direction — exactly where the projective non-Gaussianity lives.
//!
//! ## Usage
//!
//!   cargo run --release --example pose_inference [close|far]

use se3_inference::*;
use se3_inference::se3::Pose;
use se3_inference::saddlepoint::*;

// ── Tiny PRNG ──────────────────────────────────────────────────────
struct Rng(u64);
impl Rng {
    fn new(s: u64) -> Self { Self(s) }
    fn u(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn n(&mut self) -> f64 {
        let u1 = self.u().max(1e-30);
        let u2 = self.u();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

fn iso2(s: f64) -> [[f64; 2]; 2] { let v = 1.0/(s*s); [[v,0.0],[0.0,v]] }
fn iso3(s: f64) -> Mat3 { let v = 1.0/(s*s); [[v,0.0,0.0],[0.0,v,0.0],[0.0,0.0,v]] }

// ── Marginalized objective ─────────────────────────────────────────

/// Evaluate negative marginalized log-posterior at a given pose.
/// Returns (objective_value, MarginalizedMapResult).
fn eval(
    pose: &Pose,
    obs: &[([f64; 2], Vec3, Mat3)],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
) -> (f64, MarginalizedMapResult) {
    let r = evaluate_marginalized_map(&pose.rot, &pose.trans, obs, sig_zz_inv, 15);
    let v = if use_sp { -r.log_posterior_sp } else { -r.log_posterior_laplace };
    (v, r)
}

/// FD gradient of the objective w.r.t. right perturbation δξ.
fn fd_gradient(
    pose: &Pose,
    obs: &[([f64; 2], Vec3, Mat3)],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
) -> Vec6 {
    // SP objective uses inner FD (h=1e-5) for Q₄, so outer FD must be larger
    // to avoid nested-FD noise amplification. Laplace has no inner FD.
    let h = if use_sp { 1e-3 } else { 1e-5 };
    let (f0, _) = eval(pose, obs, sig_zz_inv, use_sp);
    std::array::from_fn(|j| {
        let mut ep = [0.0; 6]; ep[j] = h;
        let mut em = [0.0; 6]; em[j] = -h;
        let (vp, _) = eval(&pose.compose(&Pose::exp(&ep)), obs, sig_zz_inv, use_sp);
        let (vm, _) = eval(&pose.compose(&Pose::exp(&em)), obs, sig_zz_inv, use_sp);
        // Fall back to one-sided differences when central difference hits a singularity
        match (vp.is_finite(), vm.is_finite()) {
            (true, true) => (vp - vm) / (2.0 * h),
            (true, false) => (vp - f0) / h,
            (false, true) => (f0 - vm) / h,
            (false, false) => 0.0,
        }
    })
}

/// FD Hessian of the objective (symmetric 6×6).
fn fd_hessian(
    pose: &Pose,
    obs: &[([f64; 2], Vec3, Mat3)],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
) -> Mat6 {
    // Use larger step for SP to avoid interference with inner Q₄ FD (h=1e-5)
    let h = if use_sp { 1e-3 } else { 1e-4 };
    let (f0, _) = eval(pose, obs, sig_zz_inv, use_sp);
    let mut hess = [[0.0f64; 6]; 6];

    // Cache f(x ± h*ei), using one-sided FD when a direction hits a singularity
    let mut fp = [0.0f64; 6];
    let mut fm = [0.0f64; 6];
    for i in 0..6 {
        let mut ei = [0.0; 6]; ei[i] = h;
        (fp[i], _) = eval(&pose.compose(&Pose::exp(&ei)), obs, sig_zz_inv, use_sp);
        ei[i] = -h;
        (fm[i], _) = eval(&pose.compose(&Pose::exp(&ei)), obs, sig_zz_inv, use_sp);
        let fpi = if fp[i].is_finite() { fp[i] } else { f0 };
        let fmi = if fm[i].is_finite() { fm[i] } else { f0 };
        fp[i] = fpi;
        fm[i] = fmi;
        hess[i][i] = (fpi - 2.0 * f0 + fmi) / (h * h);
    }

    // Off-diagonal — skip if any of the 4 corners are non-finite
    for i in 0..6 { for j in (i+1)..6 {
        let mut epp = [0.0;6]; epp[i] = h; epp[j] = h;
        let mut epm = [0.0;6]; epm[i] = h; epm[j] = -h;
        let mut emp = [0.0;6]; emp[i] = -h; emp[j] = h;
        let mut emm = [0.0;6]; emm[i] = -h; emm[j] = -h;
        let (fpp, _) = eval(&pose.compose(&Pose::exp(&epp)), obs, sig_zz_inv, use_sp);
        let (fpm, _) = eval(&pose.compose(&Pose::exp(&epm)), obs, sig_zz_inv, use_sp);
        let (fmp, _) = eval(&pose.compose(&Pose::exp(&emp)), obs, sig_zz_inv, use_sp);
        let (fmm, _) = eval(&pose.compose(&Pose::exp(&emm)), obs, sig_zz_inv, use_sp);
        let v = if fpp.is_finite() && fpm.is_finite() && fmp.is_finite() && fmm.is_finite() {
            (fpp - fpm - fmp + fmm) / (4.0 * h * h)
        } else {
            0.0  // treat as uncoupled when we can't compute the cross-derivative
        };
        hess[i][j] = v;
        hess[j][i] = v;
    }}
    hess
}

/// Solve 6×6 linear system via Cholesky.
fn solve6(a: &Mat6, b: &Vec6) -> Vec6 {
    let mut l = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut s = 0.0;
            for k in 0..j { s += l[i][k] * l[j][k]; }
            l[i][j] = if i == j {
                let d = a[i][i] - s;
                if d > 1e-30 { d.sqrt() } else { 1e-15 }
            } else {
                (a[i][j] - s) / l[j][j]
            };
        }
    }
    let mut y = [0.0f64; 6];
    for i in 0..6 {
        let mut s = 0.0;
        for j in 0..i { s += l[i][j] * y[j]; }
        y[i] = (b[i] - s) / l[i][i];
    }
    let mut x = [0.0f64; 6];
    for i in (0..6).rev() {
        let mut s = 0.0;
        for j in (i+1)..6 { s += l[j][i] * x[j]; }
        x[i] = (y[i] - s) / l[i][i];
    }
    x
}

// ── Newton optimizer on SE(3) ──────────────────────────────────────

/// Optimize pose via damped Newton with compositive SE(3) updates.
fn optimize_pose(
    initial: &Pose,
    obs: &[([f64; 2], Vec3, Mat3)],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
    label: &str,
) -> Pose {
    let mut pose = *initial;
    let mut lambda = 1e-3; // LM damping

    let max_iter = if use_sp { 60 } else { 30 };
    for iter in 0..max_iter {
        let (obj, _) = eval(&pose, obs, sig_zz_inv, use_sp);
        let grad = fd_gradient(&pose, obs, sig_zz_inv, use_sp);
        let gnorm: f64 = grad.iter().map(|x| x * x).sum::<f64>().sqrt();

        if iter % 5 == 0 || gnorm < 1e-7 {
            eprintln!("  [{:>12}] iter {:>2}: obj={:>12.6}  |grad|={:.2e}  λ={:.1e}",
                label, iter, obj, gnorm, lambda);
        }
        if gnorm < 1e-7 { break; }

        // Damped Hessian
        let mut hess = fd_hessian(&pose, obs, sig_zz_inv, use_sp);
        for i in 0..6 { hess[i][i] += lambda; }

        let neg_grad: Vec6 = std::array::from_fn(|i| -grad[i]);
        let delta = solve6(&hess, &neg_grad);

        // Try the step — skip if delta is non-finite
        let step_ok = delta.iter().all(|d| d.is_finite());
        let candidate = pose.compose(&Pose::exp(&delta));
        let (new_obj, _) = if step_ok {
            eval(&candidate, obs, sig_zz_inv, use_sp)
        } else {
            (f64::INFINITY, eval(&pose, obs, sig_zz_inv, use_sp).1)
        };

        if new_obj.is_finite() && new_obj < obj {
            pose = candidate;
            lambda = (lambda * 0.3).max(1e-8);
        } else {
            lambda = (lambda * 5.0).min(1e4);
            // When LM is stuck at max damping, try a small gradient descent step
            // with backtracking line search
            if lambda >= 1e4 {
                let mut alpha = 1e-3 / gnorm;
                for _ in 0..10 {
                    let gd_step: Vec6 = std::array::from_fn(|i| -alpha * grad[i]);
                    if !gd_step.iter().all(|d| d.is_finite()) { break; }
                    let gd_candidate = pose.compose(&Pose::exp(&gd_step));
                    let (gd_obj, _) = eval(&gd_candidate, obs, sig_zz_inv, use_sp);
                    if gd_obj.is_finite() && gd_obj < obj {
                        pose = gd_candidate;
                        lambda = 1e-1; // reset damping
                        break;
                    }
                    alpha *= 0.5;
                }
            }
        }
    }
    eprintln!();
    pose
}

/// Pose error: (rotation_rad, translation_m)
fn pose_error(a: &Pose, b: &Pose) -> (f64, f64) {
    let d = a.relative(b).log();
    let rot = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]).sqrt();
    let tr = (d[3]*d[3] + d[4]*d[4] + d[5]*d[5]).sqrt();
    (rot, tr)
}

// ── Main ───────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let close = args.get(1).map(|s| s == "close").unwrap_or(false);

    let (base_depth, sigma_prior, label) = if close {
        (4.0, 3.0, "close range")
    } else {
        (8.0, 2.0, "moderate range")
    };

    println!("═══════════════════════════════════════════════════════════");
    println!("  Pose Inference from Uncertain 3D Point Cloud");
    println!("  Regime: {} (depth≈{}, σ_prior={})", label, base_depth, sigma_prior);
    println!("═══════════════════════════════════════════════════════════\n");

    let sigma_z = 0.01;
    let sig_zz_inv = iso2(sigma_z);
    let sig_xx_inv = iso3(sigma_prior);

    // ── Scene: camera looking along +z, landmarks in front ──
    // True pose: small rotation + translation
    let xi_true: Vec6 = [0.05, -0.03, 0.02, 0.1, -0.05, 0.0];
    let f_true = Pose::exp(&xi_true);

    // Landmarks at varying depths (the key: depth variation makes SP matter)
    let mut rng = Rng::new(12345);
    let n_lm = 12;
    let mut landmarks = Vec::with_capacity(n_lm);
    for i in 0..n_lm {
        let depth = base_depth + 3.0 * (i as f64 / n_lm as f64 - 0.5); // spread ±1.5
        landmarks.push([
            1.5 * (rng.u() - 0.5),
            1.5 * (rng.u() - 0.5),
            depth,
        ]);
    }

    // Generate observations: z_j = π(R·x_j + T) + noise
    let mut observations: Vec<([f64; 2], Vec3, Mat3)> = Vec::new();
    for x in &landmarks {
        let xp = f_true.act(x);
        if xp[2] < 0.5 { continue; }
        let pi = projective::project(&xp);
        let z = [pi[0] + sigma_z * rng.n(), pi[1] + sigma_z * rng.n()];
        // Prior: true position + small offset (simulate uncertainty)
        let x_prior = [
            x[0] + 0.3 * rng.n(),
            x[1] + 0.3 * rng.n(),
            x[2] + 0.3 * rng.n(),
        ];
        observations.push((z, x_prior, sig_xx_inv));
    }

    println!("  {} landmarks, σ_z={}, σ_prior={}", observations.len(), sigma_z, sigma_prior);
    println!("  True pose ξ = [{:.4}, {:.4}, {:.4}, {:.4}, {:.4}, {:.4}]\n",
        xi_true[0], xi_true[1], xi_true[2], xi_true[3], xi_true[4], xi_true[5]);

    // ── Part 1: 1D sweep along translation-z ──
    println!("──────────────────────────────────────────────────────────");
    println!("  Part 1: 1D sweep along z-translation");
    println!("──────────────────────────────────────────────────────────\n");
    println!("{:>8} {:>14} {:>14} {:>12}  {:>6} {:>6}",
        "δt_z", "neg_log_Lap", "neg_log_SP", "Δ(SP-Lap)", "n_val", "n_inv");

    let n_sweep = 21;
    let sweep_range = 0.3;
    let mut best_lap = (f64::MAX, 0.0f64);
    let mut best_sp = (f64::MAX, 0.0f64);

    for k in 0..n_sweep {
        let dt = -sweep_range + 2.0 * sweep_range * (k as f64) / (n_sweep - 1) as f64;
        let mut xi_pert = xi_true;
        xi_pert[5] += dt;
        let pose = Pose::exp(&xi_pert);
        let (obj_lap, r_lap) = eval(&pose, &observations, &sig_zz_inv, false);
        let (obj_sp, _) = eval(&pose, &observations, &sig_zz_inv, true);

        if obj_lap < best_lap.0 { best_lap = (obj_lap, dt); }
        if obj_sp < best_sp.0 { best_sp = (obj_sp, dt); }

        println!("{:>8.4} {:>14.6} {:>14.6} {:>12.6}  {:>6} {:>6}",
            dt, obj_lap, obj_sp, obj_sp - obj_lap,
            r_lap.n_valid, r_lap.n_invalid);
    }

    println!("\n  Laplace minimum at δt_z = {:.4}", best_lap.1);
    println!("  SP minimum at      δt_z = {:.4}", best_sp.1);
    println!("  Shift:              Δt_z = {:.4}", best_sp.1 - best_lap.1);

    // ── Part 2: Full 6D optimization ──
    println!("\n──────────────────────────────────────────────────────────");
    println!("  Part 2: Full 6D pose optimization (Newton on SE(3))");
    println!("──────────────────────────────────────────────────────────\n");

    // Initial guess: perturbed from truth
    let xi_init: Vec6 = [
        xi_true[0] + 0.03,
        xi_true[1] - 0.02,
        xi_true[2] + 0.01,
        xi_true[3] + 0.15,
        xi_true[4] - 0.1,
        xi_true[5] + 0.2,
    ];
    let f_init = Pose::exp(&xi_init);

    let (init_rot_err, init_trans_err) = pose_error(&f_true, &f_init);
    println!("  Initial pose error: rot={:.4} rad, trans={:.4} m\n", init_rot_err, init_trans_err);

    // Optimize with Laplace
    eprintln!("Optimizing with Laplace objective...");
    let f_lap = optimize_pose(&f_init, &observations, &sig_zz_inv, false, "Laplace");

    // Optimize with saddlepoint
    eprintln!("Optimizing with saddlepoint objective...");
    let f_sp = optimize_pose(&f_init, &observations, &sig_zz_inv, true, "Saddlepoint");

    // ── Part 3: Comparison ──
    println!("──────────────────────────────────────────────────────────");
    println!("  Part 3: Results");
    println!("──────────────────────────────────────────────────────────\n");

    let (lap_rot, lap_trans) = pose_error(&f_true, &f_lap);
    let (sp_rot, sp_trans) = pose_error(&f_true, &f_sp);
    let (lap_sp_rot, lap_sp_trans) = pose_error(&f_lap, &f_sp);

    let xi_lap = f_lap.log();
    let xi_sp = f_sp.log();

    println!("  {:>20} {:>12} {:>12} {:>12}", "", "True", "Laplace", "Saddlepoint");
    for (i, name) in ["ω₁","ω₂","ω₃","t₁","t₂","t₃"].iter().enumerate() {
        println!("  {:>20} {:>12.6} {:>12.6} {:>12.6}", name, xi_true[i], xi_lap[i], xi_sp[i]);
    }

    println!("\n  Pose error vs ground truth:");
    println!("    Laplace:      rot={:.6} rad, trans={:.6} m", lap_rot, lap_trans);
    println!("    Saddlepoint:  rot={:.6} rad, trans={:.6} m", sp_rot, sp_trans);
    println!("    Difference:   rot={:.6} rad, trans={:.6} m", lap_sp_rot, lap_sp_trans);

    // Objective values at the optima
    let (obj_lap_at_lap, r_lap) = eval(&f_lap, &observations, &sig_zz_inv, false);
    let (obj_sp_at_sp, r_sp) = eval(&f_sp, &observations, &sig_zz_inv, true);

    println!("\n  Objective at optimum:");
    println!("    Laplace:      {:.6}  (valid: {}, invalid: {})", obj_lap_at_lap, r_lap.n_valid, r_lap.n_invalid);
    println!("    Saddlepoint:  {:.6}  (valid: {}, invalid: {})", obj_sp_at_sp, r_sp.n_valid, r_sp.n_invalid);

    // SP correction magnitudes at the Laplace optimum
    let (_, r_at_lap) = eval(&f_lap, &observations, &sig_zz_inv, true);
    let c1_vals: Vec<f64> = r_at_lap.sp_results.iter().map(|r| r.c1).collect();
    let mean_c1 = c1_vals.iter().sum::<f64>() / c1_vals.len() as f64;
    let max_c1 = c1_vals.iter().fold(0.0f64, |a, b| a.max(b.abs()));

    println!("\n  SP corrections at Laplace optimum:");
    println!("    mean c₁ = {:+.4e},  max |c₁| = {:.4e}", mean_c1, max_c1);

    // ── Part 4: Per-landmark detail ──
    println!("\n──────────────────────────────────────────────────────────");
    println!("  Part 4: Per-landmark saddlepoint corrections");
    println!("──────────────────────────────────────────────────────────\n");
    println!("{:>4} {:>8} {:>12} {:>12} {:>12} {:>4}",
        "LM", "depth", "c₁", "A/12", "-Q₄/8", "ok");

    for (i, r) in r_at_lap.sp_results.iter().enumerate() {
        let xp = projective::transform_point(&f_lap.rot, &f_lap.trans, &observations[i].1);
        let ok = matches!(r.status, SaddlepointStatus::Valid);
        println!("{:>4} {:>8.2} {:>12.4e} {:>12.4e} {:>12.4e} {:>4}",
            i, xp[2], r.c1, r.term_a / 12.0, -r.term_q4 / 8.0,
            if ok { "✓" } else { "✗" });
    }

    println!("\n  Closer landmarks → larger |c₁| → SP correction matters more");
    if lap_sp_trans > 1e-6 {
        println!("  Laplace–SP pose difference: {:.2e} m (translation)", lap_sp_trans);
        println!("  This is the systematic bias that saddlepoint corrects.");
    }
}
