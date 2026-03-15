//! # Pose Inference from Uncertain 3D Point Cloud
//!
//! Estimates camera pose(s) by maximizing the marginalized posterior over
//! landmark positions. Compares Laplace vs saddlepoint-corrected objectives.
//!
//! Supports single-camera and multi-camera configurations. Each landmark
//! can be observed by one or more cameras, controlled by the overlap parameter.
//!
//! ## Usage
//!
//!   cargo run --release --example pose_inference [OPTIONS]
//!
//! Options (positional):
//!   N_cameras     Number of cameras (default: 1)
//!   N_landmarks   Number of landmarks (default: 12)
//!   overlap       Fraction of landmarks visible to each camera, 0.0-1.0 (default: 1.0)
//!
//! Examples:
//!   cargo run --release --example pose_inference                  # 1 camera, 12 landmarks
//!   cargo run --release --example pose_inference 4                # 4 cameras, 12 landmarks
//!   cargo run --release --example pose_inference 4 30 0.5         # 4 cameras, 30 landmarks, 50% overlap
//!   cargo run --release --example pose_inference 2 20 0.7         # stereo, 20 landmarks, 70% overlap

use se3_inference::*;
use se3_inference::se3::Pose;
use se3_inference::saddlepoint::*;
use se3_inference::projective;

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

// ── Camera ring ────────────────────────────────────────────────────
/// Place n cameras on a ring of given radius, all looking at center.
fn ring_cameras(n: usize, radius: f64, center: &Vec3) -> Vec<Pose> {
    (0..n).map(|i| {
        let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let h = 0.5 * (if i % 2 == 0 { 1.0 } else { -1.0 });
        let pos = [center[0] + radius * a.cos(), center[1] + h, center[2] + radius * a.sin()];
        let z_dir = { let d = sub3(center, &pos); scale3(1.0 / norm3(&d), &d) };
        let xr = cross3(&[0.0, 1.0, 0.0], &z_dir);
        let x_dir = if norm3(&xr) > 1e-6 { scale3(1.0 / norm3(&xr), &xr) } else { [1.0, 0.0, 0.0] };
        let y_dir = cross3(&z_dir, &x_dir);
        let r = [[x_dir[0],x_dir[1],x_dir[2]], [y_dir[0],y_dir[1],y_dir[2]], [z_dir[0],z_dir[1],z_dir[2]]];
        Pose::new(r, scale3(-1.0, &mv3(&r, &pos)))
    }).collect()
}

/// Check if a point is visible from a camera (in front and within FOV).
fn visible(cam: &Pose, x: &Vec3, fov_tan: f64) -> bool {
    let p = cam.act(x);
    p[2] > 0.5 && (p[0]/p[2]).abs() < fov_tan && (p[1]/p[2]).abs() < fov_tan
}

// ── Pose prior ───────────────────────────────────────────────────────

/// Gaussian prior on the pose in the Lie algebra.
struct PosePrior {
    pose_ref: Pose,
    sigma_rot: f64,
    sigma_trans: f64,
}

impl PosePrior {
    fn nll(&self, pose: &Pose) -> f64 {
        let delta = self.pose_ref.relative(pose).log();
        let sr2 = 1.0 / (self.sigma_rot * self.sigma_rot);
        let st2 = 1.0 / (self.sigma_trans * self.sigma_trans);
        0.5 * ((delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]) * sr2
             + (delta[3]*delta[3] + delta[4]*delta[4] + delta[5]*delta[5]) * st2)
    }
}

// ── Landmark data ─────────────────────────────────────────────────

/// A landmark with its world-space prior and observations from multiple cameras.
struct Landmark {
    x_prior: Vec3,
    sigma_xx_inv: Mat3,
    /// (camera_index, observation_z)
    observations: Vec<(usize, [f64; 2])>,
}

// ── Marginalized objective (multi-camera) ─────────────────────────

#[allow(dead_code)]
struct EvalResult {
    log_posterior_sp: f64,
    log_posterior_laplace: f64,
    n_valid: usize,
    n_invalid: usize,
    sp_results: Vec<SaddlepointResult>,
    /// Camera-frame depths per landmark (for reporting)
    depths: Vec<f64>,
}

/// Evaluate negative marginalized log-posterior at given camera poses.
fn eval_multicam(
    poses: &[Pose],
    landmarks: &[Landmark],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
    priors: &[PosePrior],
) -> (f64, EvalResult) {
    let mut log_sp = 0.0;
    let mut log_laplace = 0.0;
    let (mut n_valid, mut n_invalid) = (0, 0);
    let mut sp_results = Vec::new();
    let mut depths = Vec::new();

    for lm in landmarks {
        if lm.observations.is_empty() { continue; }

        if lm.observations.len() == 1 {
            // Single-camera path: use original API
            let (ci, z) = &lm.observations[0];
            let pose = &poses[*ci];
            let opt = optimize_landmark(
                &pose.rot, &pose.trans, z, sig_zz_inv,
                &lm.x_prior, &lm.sigma_xx_inv, 15,
            );
            let lap_c = landmark_marginal_log_posterior_laplace(&opt);
            let (sp_c, sp_r) = landmark_marginal_log_posterior(&opt, &pose.rot, sig_zz_inv);
            let xp = projective::transform_point(&pose.rot, &pose.trans, &opt.x_opt);
            depths.push(xp[2]);
            match sp_r.status { SaddlepointStatus::Valid => n_valid += 1, _ => n_invalid += 1 }
            log_laplace += lap_c;
            log_sp += sp_c;
            sp_results.push(sp_r);
        } else {
            // Multi-camera path
            let cam_obs: Vec<CameraObs> = lm.observations.iter().map(|(ci, z)| {
                CameraObs {
                    rot: poses[*ci].rot,
                    trans: poses[*ci].trans,
                    z: *z,
                    sigma_zz_inv: *sig_zz_inv,
                }
            }).collect();
            let opt = optimize_landmark_multicam(&cam_obs, &lm.x_prior, &lm.sigma_xx_inv, 15);
            let lap_c = -opt.nll_opt - 0.5 * opt.log_det_h
                + 1.5 * (2.0 * std::f64::consts::PI).ln();
            let (sp_c, sp_r) = landmark_marginal_multicam(&opt, &cam_obs);
            // Use the min depth across cameras
            let min_depth = opt.xp_opts.iter().map(|xp| xp[2]).fold(f64::MAX, f64::min);
            depths.push(min_depth);
            match sp_r.status { SaddlepointStatus::Valid => n_valid += 1, _ => n_invalid += 1 }
            log_laplace += lap_c;
            log_sp += sp_c;
            sp_results.push(sp_r);
        }
    }

    let obj = if use_sp { -log_sp } else { -log_laplace };
    let prior_nll: f64 = poses.iter().zip(priors.iter()).map(|(p, pr)| pr.nll(p)).sum();

    (obj + prior_nll, EvalResult {
        log_posterior_sp: log_sp, log_posterior_laplace: log_laplace,
        n_valid, n_invalid, sp_results, depths,
    })
}

/// FD gradient of the objective w.r.t. right perturbation of camera `cam_idx`.
fn fd_gradient_cam(
    poses: &[Pose],
    cam_idx: usize,
    landmarks: &[Landmark],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
    priors: &[PosePrior],
) -> Vec6 {
    let h = 1e-5;
    let (f0, _) = eval_multicam(poses, landmarks, sig_zz_inv, use_sp, priors);
    std::array::from_fn(|j| {
        let mut ep = [0.0; 6]; ep[j] = h;
        let mut em = [0.0; 6]; em[j] = -h;
        let mut poses_p = poses.to_vec();
        let mut poses_m = poses.to_vec();
        poses_p[cam_idx] = poses[cam_idx].compose(&Pose::exp(&ep));
        poses_m[cam_idx] = poses[cam_idx].compose(&Pose::exp(&em));
        let (vp, _) = eval_multicam(&poses_p, landmarks, sig_zz_inv, use_sp, priors);
        let (vm, _) = eval_multicam(&poses_m, landmarks, sig_zz_inv, use_sp, priors);
        match (vp.is_finite(), vm.is_finite()) {
            (true, true) => (vp - vm) / (2.0 * h),
            (true, false) => (vp - f0) / h,
            (false, true) => (f0 - vm) / h,
            (false, false) => 0.0,
        }
    })
}

/// FD Hessian of the objective w.r.t. camera `cam_idx` (6×6).
fn fd_hessian_cam(
    poses: &[Pose],
    cam_idx: usize,
    landmarks: &[Landmark],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
    priors: &[PosePrior],
) -> Mat6 {
    let h = 1e-4;
    let perturb = |delta: &Vec6| -> Vec<Pose> {
        let mut ps = poses.to_vec();
        ps[cam_idx] = poses[cam_idx].compose(&Pose::exp(delta));
        ps
    };
    let (f0, _) = eval_multicam(poses, landmarks, sig_zz_inv, use_sp, priors);
    let mut hess = [[0.0f64; 6]; 6];

    let mut fp = [0.0f64; 6];
    let mut fm = [0.0f64; 6];
    for i in 0..6 {
        let mut ei = [0.0; 6]; ei[i] = h;
        (fp[i], _) = eval_multicam(&perturb(&ei), landmarks, sig_zz_inv, use_sp, priors);
        ei[i] = -h;
        (fm[i], _) = eval_multicam(&perturb(&ei), landmarks, sig_zz_inv, use_sp, priors);
        fp[i] = if fp[i].is_finite() { fp[i] } else { f0 };
        fm[i] = if fm[i].is_finite() { fm[i] } else { f0 };
        hess[i][i] = (fp[i] - 2.0 * f0 + fm[i]) / (h * h);
    }

    for i in 0..6 { for j in (i+1)..6 {
        let mut epp = [0.0;6]; epp[i] = h; epp[j] = h;
        let mut epm = [0.0;6]; epm[i] = h; epm[j] = -h;
        let mut emp = [0.0;6]; emp[i] = -h; emp[j] = h;
        let mut emm = [0.0;6]; emm[i] = -h; emm[j] = -h;
        let (fpp, _) = eval_multicam(&perturb(&epp), landmarks, sig_zz_inv, use_sp, priors);
        let (fpm, _) = eval_multicam(&perturb(&epm), landmarks, sig_zz_inv, use_sp, priors);
        let (fmp, _) = eval_multicam(&perturb(&emp), landmarks, sig_zz_inv, use_sp, priors);
        let (fmm, _) = eval_multicam(&perturb(&emm), landmarks, sig_zz_inv, use_sp, priors);
        let v = if fpp.is_finite() && fpm.is_finite() && fmp.is_finite() && fmm.is_finite() {
            (fpp - fpm - fmp + fmm) / (4.0 * h * h)
        } else {
            0.0
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

// ── Block-coordinate Newton optimizer ──────────────────────────────

/// Optimize all camera poses via damped Newton with simultaneous updates.
/// All per-camera Newton steps are computed at the current poses, then
/// applied simultaneously (Jacobi-style), avoiding the zig-zag convergence
/// issue of sequential block-coordinate updates on coupled landmarks.
fn optimize_poses(
    initial: &[Pose],
    landmarks: &[Landmark],
    sig_zz_inv: &[[f64; 2]; 2],
    use_sp: bool,
    priors: &[PosePrior],
    label: &str,
) -> Vec<Pose> {
    let mut poses = initial.to_vec();
    let n_cam = poses.len();
    let max_iter = 30;
    let mut lambda = 1e-3;

    for iter in 0..max_iter {
        let (obj, _) = eval_multicam(&poses, landmarks, sig_zz_inv, use_sp, priors);

        // Compute all per-camera gradients and Newton steps at current poses
        let mut deltas: Vec<Vec6> = Vec::with_capacity(n_cam);
        let mut total_gnorm_sq = 0.0;

        for ci in 0..n_cam {
            let grad = fd_gradient_cam(&poses, ci, landmarks, sig_zz_inv, use_sp, priors);
            total_gnorm_sq += grad.iter().map(|x| x * x).sum::<f64>();

            let mut hess = fd_hessian_cam(&poses, ci, landmarks, sig_zz_inv, use_sp, priors);
            for i in 0..6 { hess[i][i] += lambda; }

            let neg_grad: Vec6 = std::array::from_fn(|i| -grad[i]);
            let delta = solve6(&hess, &neg_grad);
            deltas.push(delta);
        }
        let gnorm = total_gnorm_sq.sqrt();

        if iter % 5 == 0 || gnorm < 1e-7 {
            eprintln!("  [{:>12}] iter {:>2}: obj={:>12.6}  |grad|={:.2e}  λ={:.1e}",
                label, iter, obj, gnorm, lambda);
        }
        if gnorm < 1e-7 { break; }

        // Apply all updates simultaneously
        let mut poses_new = poses.clone();
        let mut all_finite = true;
        for ci in 0..n_cam {
            if !deltas[ci].iter().all(|d| d.is_finite()) { all_finite = false; break; }
            poses_new[ci] = poses[ci].compose(&Pose::exp(&deltas[ci]));
        }

        if all_finite {
            let (new_obj, _) = eval_multicam(&poses_new, landmarks, sig_zz_inv, use_sp, priors);
            if new_obj.is_finite() && new_obj < obj {
                poses = poses_new;
                lambda = (lambda * 0.3).max(1e-8);
            } else {
                lambda = (lambda * 5.0).min(1e4);
                // Gradient descent fallback when Newton fails
                if lambda >= 1e4 {
                    let alpha = 1e-3 / gnorm;
                    let mut gd_poses = poses.clone();
                    for ci in 0..n_cam {
                        let grad = fd_gradient_cam(&poses, ci, landmarks, sig_zz_inv, use_sp, priors);
                        let step: Vec6 = std::array::from_fn(|i| -alpha * grad[i]);
                        if step.iter().all(|d| d.is_finite()) {
                            gd_poses[ci] = poses[ci].compose(&Pose::exp(&step));
                        }
                    }
                    let (gd_obj, _) = eval_multicam(&gd_poses, landmarks, sig_zz_inv, use_sp, priors);
                    if gd_obj.is_finite() && gd_obj < obj {
                        poses = gd_poses;
                        lambda = 1e-1;
                    }
                }
            }
        } else {
            lambda = (lambda * 5.0).min(1e4);
        }
    }
    eprintln!();
    poses
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
    let n_cam: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1);
    let n_lm: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(12);
    let overlap: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let overlap = overlap.clamp(0.1, 1.0);

    let sigma_z = 0.01;
    let sigma_prior = 2.0;
    let sig_zz_inv = iso2(sigma_z);
    let sig_xx_inv = iso3(sigma_prior);
    let fov: f64 = 60.0;
    let fov_tan = (fov.to_radians() / 2.0).tan();

    println!("═══════════════════════════════════════════════════════════");
    println!("  Pose Inference from Uncertain 3D Point Cloud");
    println!("═══════════════════════════════════════════════════════════\n");
    println!("  {} camera(s), {} landmarks, {:.0}% overlap", n_cam, n_lm, overlap * 100.0);
    println!("  σ_z={}, σ_prior={}, FOV={:.0}°\n", sigma_z, sigma_prior, fov);

    let mut rng = Rng::new(12345);

    // ── Camera setup ──
    let radius = 8.0;
    let center = [0.0; 3];

    // True poses: ring cameras with small perturbations
    let nominal_poses = if n_cam == 1 {
        // Single camera: looking along +z from -8m
        vec![Pose::exp(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])]
    } else {
        ring_cameras(n_cam, radius, &center)
    };

    // Apply small true perturbations to each camera
    let true_perturbations: Vec<Vec6> = (0..n_cam).map(|_| {
        [0.05 * rng.n(), 0.05 * rng.n(), 0.05 * rng.n(),
         0.1 * rng.n(), 0.1 * rng.n(), 0.1 * rng.n()]
    }).collect();

    let true_poses: Vec<Pose> = nominal_poses.iter().zip(true_perturbations.iter())
        .map(|(nom, pert)| nom.compose(&Pose::exp(pert)))
        .collect();

    // Print camera positions
    for (i, pose) in true_poses.iter().enumerate() {
        let p = pose.inverse().trans;
        println!("  Cam {}: [{:+.2}, {:+.2}, {:+.2}]", i, p[0], p[1], p[2]);
    }

    // ── Generate landmarks ──
    // For single camera: landmarks in front at varying depths
    // For multi-camera: landmarks near the center of the ring
    let landmarks_world: Vec<Vec3> = (0..n_lm).map(|i| {
        if n_cam == 1 {
            let depth = 8.0 + 3.0 * (i as f64 / n_lm as f64 - 0.5);
            [1.5 * (rng.u() - 0.5), 1.5 * (rng.u() - 0.5), depth]
        } else {
            [3.0 * (rng.u() - 0.5), 3.0 * (rng.u() - 0.5), 3.0 * (rng.u() - 0.5)]
        }
    }).collect();

    // ── Build landmark observations ──
    // Each landmark is observed by cameras that can see it AND pass the overlap filter.
    let mut landmarks: Vec<Landmark> = Vec::new();
    let mut total_obs = 0usize;

    for (li, xw) in landmarks_world.iter().enumerate() {
        let mut obs = Vec::new();
        for (ci, cam) in true_poses.iter().enumerate() {
            if !visible(cam, xw, fov_tan) { continue; }
            // Overlap filter: use a deterministic hash to decide if this
            // camera observes this landmark (ensures reproducibility)
            let hash_val = ((li * 97 + ci * 31) % 1000) as f64 / 1000.0;
            if hash_val >= overlap { continue; }

            let xp = cam.act(xw);
            let pi = projective::project(&xp);
            let z = [pi[0] + sigma_z * rng.n(), pi[1] + sigma_z * rng.n()];
            obs.push((ci, z));
        }
        if obs.is_empty() { continue; }
        total_obs += obs.len();

        let x_prior = [
            xw[0] + 0.3 * rng.n(),
            xw[1] + 0.3 * rng.n(),
            xw[2] + 0.3 * rng.n(),
        ];
        landmarks.push(Landmark {
            x_prior,
            sigma_xx_inv: sig_xx_inv,
            observations: obs,
        });
    }

    println!("\n  {} landmarks visible, {} total observations ({:.1} obs/landmark)\n",
        landmarks.len(), total_obs, total_obs as f64 / landmarks.len().max(1) as f64);

    // ── Initial guesses: perturbed from truth ──
    let init_poses: Vec<Pose> = true_poses.iter().map(|tp| {
        let pert = [0.03 * rng.n(), 0.03 * rng.n(), 0.03 * rng.n(),
                    0.15 * rng.n(), 0.15 * rng.n(), 0.15 * rng.n()];
        tp.compose(&Pose::exp(&pert))
    }).collect();

    let priors: Vec<PosePrior> = init_poses.iter().map(|p| PosePrior {
        pose_ref: *p,
        sigma_rot: 0.3,
        sigma_trans: 2.0,
    }).collect();

    // Print initial errors
    println!("  Initial pose errors:");
    for (i, (tp, ip)) in true_poses.iter().zip(init_poses.iter()).enumerate() {
        let (re, te) = pose_error(tp, ip);
        println!("    Cam {}: rot={:.4} rad, trans={:.4} m", i, re, te);
    }

    // ── Part 1: Optimization ──
    println!("\n──────────────────────────────────────────────────────────");
    println!("  Part 1: Pose optimization (block-coordinate Newton)");
    println!("──────────────────────────────────────────────────────────\n");

    eprintln!("Optimizing with Laplace objective...");
    let lap_poses = optimize_poses(&init_poses, &landmarks, &sig_zz_inv, false, &priors, "Laplace");

    eprintln!("Optimizing with saddlepoint objective...");
    let sp_poses = optimize_poses(&init_poses, &landmarks, &sig_zz_inv, true, &priors, "Saddlepoint");

    // ── Part 2: Results ──
    println!("──────────────────────────────────────────────────────────");
    println!("  Part 2: Results");
    println!("──────────────────────────────────────────────────────────\n");

    println!("  {:>6} {:>12} {:>12} {:>12} {:>12}",
        "Cam", "Lap rot", "Lap trans", "SP rot", "SP trans");
    println!("  {}", "─".repeat(58));

    let mut total_lap_rot = 0.0;
    let mut total_lap_trans = 0.0;
    let mut total_sp_rot = 0.0;
    let mut total_sp_trans = 0.0;

    for i in 0..n_cam {
        let (lr, lt) = pose_error(&true_poses[i], &lap_poses[i]);
        let (sr, st) = pose_error(&true_poses[i], &sp_poses[i]);
        total_lap_rot += lr; total_lap_trans += lt;
        total_sp_rot += sr; total_sp_trans += st;
        println!("  {:>6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}",
            i, lr, lt, sr, st);
    }

    if n_cam > 1 {
        println!("  {}", "─".repeat(58));
        println!("  {:>6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}",
            "mean",
            total_lap_rot / n_cam as f64, total_lap_trans / n_cam as f64,
            total_sp_rot / n_cam as f64, total_sp_trans / n_cam as f64);
    }

    // Pose-space difference between Laplace and SP solutions
    println!("\n  Laplace–SP pose difference:");
    for i in 0..n_cam {
        let (dr, dt) = pose_error(&lap_poses[i], &sp_poses[i]);
        println!("    Cam {}: Δrot={:.6} rad, Δtrans={:.6} m", i, dr, dt);
    }

    // Objective values
    let (obj_lap, r_lap) = eval_multicam(&lap_poses, &landmarks, &sig_zz_inv, false, &priors);
    let (obj_sp, r_sp) = eval_multicam(&sp_poses, &landmarks, &sig_zz_inv, true, &priors);

    println!("\n  Objective at optimum:");
    println!("    Laplace:      {:.6}  (valid: {}, invalid: {})", obj_lap, r_lap.n_valid, r_lap.n_invalid);
    println!("    Saddlepoint:  {:.6}  (valid: {}, invalid: {})", obj_sp, r_sp.n_valid, r_sp.n_invalid);

    // ── Part 3: Per-landmark detail ──
    println!("\n──────────────────────────────────────────────────────────");
    println!("  Part 3: Per-landmark saddlepoint corrections");
    println!("──────────────────────────────────────────────────────────\n");

    let (_, r_at_lap) = eval_multicam(&lap_poses, &landmarks, &sig_zz_inv, true, &priors);

    println!("{:>4} {:>5} {:>8} {:>12} {:>12} {:>12} {:>4}",
        "LM", "views", "depth", "c₁", "A/12", "-Q₄/8", "ok");
    println!("{}", "─".repeat(62));

    for (i, (sp_r, depth)) in r_at_lap.sp_results.iter().zip(r_at_lap.depths.iter()).enumerate() {
        let ok = matches!(sp_r.status, SaddlepointStatus::Valid);
        let n_views = landmarks[i].observations.len();
        println!("{:>4} {:>5} {:>8.2} {:>12.4e} {:>12.4e} {:>12.4e} {:>4}",
            i, n_views, depth, sp_r.c1, sp_r.term_a / 12.0, -sp_r.term_q4 / 8.0,
            if ok { "✓" } else { "✗" });
    }

    // Summary statistics
    let valid_c1: Vec<f64> = r_at_lap.sp_results.iter()
        .filter(|r| matches!(r.status, SaddlepointStatus::Valid))
        .map(|r| r.c1).collect();
    if !valid_c1.is_empty() {
        let mean_c1 = valid_c1.iter().sum::<f64>() / valid_c1.len() as f64;
        let max_c1 = valid_c1.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        println!("\n  SP corrections at Laplace optimum:");
        println!("    mean c₁ = {:+.4e},  max |c₁| = {:.4e}", mean_c1, max_c1);
    }

    // Key insight for single vs multi-camera
    if n_cam == 1 {
        let lap_sp_trans = pose_error(&lap_poses[0], &sp_poses[0]).1;
        if lap_sp_trans > 1e-6 && total_lap_trans > 1.0 {
            println!("\n  Single camera: depth soft mode dominates (trans err ≈ {:.1} m)", total_lap_trans);
            println!("  SP correction is only {:.2e} m — geometric, not approximation, limitation.", lap_sp_trans);
        }
    } else {
        println!("\n  {} cameras: depth triangulated via parallax", n_cam);
        println!("  Mean translation error: Laplace={:.4} m, SP={:.4} m",
            total_lap_trans / n_cam as f64, total_sp_trans / n_cam as f64);
    }
}
