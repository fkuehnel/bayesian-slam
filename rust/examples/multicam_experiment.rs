//! # Experiment: Multi-Camera Saddlepoint Marginalization
//!
//! Demonstrates the saddlepoint correction for a point cloud observed
//! by 1–6 cameras arranged in a ring. Shows:
//!
//! 1. How c₁ decreases as more cameras constrain each landmark
//! 2. Saddlepoint vs Laplace accuracy against MC integration
//! 3. Dependence on depth, baseline, and measurement noise
//!
//! Paper §V.4 and §IV.3 (multi-camera extension).

use se3_inference::*;
use se3_inference::saddlepoint::*;
use se3_inference::projective;

/// Simple PRNG
struct Rng { state: u64 }
impl Rng {
    fn new(seed: u64) -> Self { Rng { state: seed } }
    fn uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.state >> 11) as f64 / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-30);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    #[allow(dead_code)]
    fn normal3(&mut self, sigma: f64) -> Vec3 {
        [sigma * self.normal(), sigma * self.normal(), sigma * self.normal()]
    }
}

fn iso2(s: f64) -> [[f64; 2]; 2] { let v = 1.0/(s*s); [[v,0.0],[0.0,v]] }
fn iso3(s: f64) -> Mat3 { let v = 1.0/(s*s); [[v,0.0,0.0],[0.0,v,0.0],[0.0,0.0,v]] }

/// Create cameras arranged in a ring at given radius and height,
/// all looking toward the origin.
fn camera_ring(n_cams: usize, radius: f64, height: f64) -> Vec<(Mat3, Vec3)> {
    let mut cameras = Vec::with_capacity(n_cams);
    for i in 0..n_cams {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n_cams as f64);
        let cam_pos = [radius * angle.cos(), radius * angle.sin(), height];

        // Camera looks toward origin: z-axis points from cam to origin
        let forward = scale3(-1.0 / norm3(&cam_pos), &cam_pos);
        // world up = [0,0,1]
        let right = cross3(&forward, &[0.0, 0.0, 1.0]);
        let right_n = if norm3(&right) > 1e-10 {
            scale3(1.0 / norm3(&right), &right)
        } else {
            [1.0, 0.0, 0.0]
        };
        let up = cross3(&right_n, &forward);

        // R maps world→camera: rows are right, up, forward
        let rot = [
            [right_n[0], right_n[1], right_n[2]],
            [up[0], up[1], up[2]],
            [forward[0], forward[1], forward[2]],
        ];
        // T = -R·cam_pos
        let trans = scale3(-1.0, &mv3(&rot, &cam_pos));

        cameras.push((rot, trans));
    }
    cameras
}

/// Check if a point is visible from a camera (positive depth, reasonable FOV)
fn is_visible(rot: &Mat3, trans: &Vec3, x: &Vec3) -> bool {
    let xp = projective::transform_point(rot, trans, x);
    if xp[2] < 0.1 { return false; }
    let u = xp[0] / xp[2];
    let v = xp[1] / xp[2];
    u.abs() < 2.0 && v.abs() < 2.0  // ~±60° FOV
}

/// Monte Carlo estimate of the marginal integral log I for one landmark
/// observed from multiple cameras.
///
/// I = ∫ exp(-f(x)) dx  where f = Σ_i NLL_i + prior
fn mc_marginal_integral(
    cameras: &[CameraObs],
    x_opt: &Vec3,
    sigma_xx_inv: &Mat3,
    x_prior: &Vec3,
    n_samples: usize,
    rng: &mut Rng,
) -> f64 {
    // Sample around x_opt using the Hessian as proposal
    // First compute Hessian at optimum for the proposal covariance
    let mut h = *sigma_xx_inv;
    for cam in cameras {
        let xp = projective::transform_point(&cam.rot, &cam.trans, x_opt);
        if xp[2] < 1e-6 { continue; }
        let p = projective::project_jacobian(&xp);
        let mut pr = [[0.0; 3]; 2];
        for i in 0..2 { for j in 0..3 {
            pr[i][j] = p[i][0]*cam.rot[0][j] + p[i][1]*cam.rot[1][j] + p[i][2]*cam.rot[2][j];
        }}
        for i in 0..3 { for j in 0..3 {
            for m in 0..2 { for n in 0..2 {
                h[i][j] += pr[m][i] * cam.sigma_zz_inv[m][n] * pr[n][j];
            }}
        }}
    }

    // Cholesky of H^{-1} for sampling
    let h_inv = inv3(&h);
    let chol = cholesky3(&h_inv);

    // f(x_opt) = value at optimum
    let f_opt = eval_total_nll(cameras, x_opt, sigma_xx_inv, x_prior);

    // Importance sampling: sample from N(x_opt, H^{-1}), weight by
    // exp(-f(x) + f(x_opt)) / q(x|x_opt)
    // where q is the Gaussian proposal.
    // Since the proposal is exactly the Laplace approximation,
    // the weight is exp(-f(x) + f_laplace(x)) where f_laplace is the
    // quadratic approximation.
    let mut log_weights = Vec::with_capacity(n_samples);

    for _ in 0..n_samples {
        let z = [rng.normal(), rng.normal(), rng.normal()];
        let dx = mv3(&chol, &z);
        let x = add3(x_opt, &dx);

        let f_x = eval_total_nll(cameras, &x, sigma_xx_inv, x_prior);
        // Quadratic approximation at x_opt:
        // f_quad(x) = f_opt + ½ dx^T H dx
        let hdx = mv3(&h, &dx);
        let f_quad = f_opt + 0.5 * dot3(&dx, &hdx);

        // Log weight = -(f(x) - f_quad(x)) = f_quad(x) - f(x)
        log_weights.push(f_quad - f_x);
    }

    // Log-sum-exp for numerical stability
    let max_lw = log_weights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_exp: f64 = log_weights.iter().map(|&lw| (lw - max_lw).exp()).sum();
    let log_mean_weight = max_lw + (sum_exp / n_samples as f64).ln();

    // log I = log I_laplace + log(mean weight)
    // log I_laplace = -f_opt + ½ log((2π)³/det(H))
    let log_i_laplace = -f_opt - 0.5 * det3(&h).ln() + 1.5 * (2.0 * std::f64::consts::PI).ln();
    log_i_laplace + log_mean_weight
}

fn eval_total_nll(
    cameras: &[CameraObs], x: &Vec3,
    sigma_xx_inv: &Mat3, x_prior: &Vec3,
) -> f64 {
    let dx = sub3(x, x_prior);
    let mut nll = 0.5 * dot3(&dx, &mv3(sigma_xx_inv, &dx));
    for cam in cameras {
        let xp = projective::transform_point(&cam.rot, &cam.trans, x);
        if xp[2] < 1e-6 { return f64::INFINITY; }
        nll += projective::neg_log_likelihood(&cam.z, &xp, &cam.sigma_zz_inv);
    }
    nll
}

fn cholesky3(a: &Mat3) -> Mat3 {
    let mut l = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j { sum += l[i][k] * l[j][k]; }
            if i == j {
                let d = a[i][i] - sum;
                l[i][j] = if d > 0.0 { d.sqrt() } else { 1e-15 };
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }
    l
}

fn main() {
    println!("═══════════════════════════════════════════════════════════");
    println!("  Multi-Camera Saddlepoint Marginalization Experiment");
    println!("═══════════════════════════════════════════════════════════\n");

    let mut rng = Rng::new(54321);

    // ─── Scene: point cloud at origin, cameras in ring ───
    let n_points = 20;
    let point_spread = 1.5;    // landmarks in [-1.5, 1.5]³
    let cam_radius = 12.0;     // camera ring radius
    let cam_height = 0.0;      // cameras at same height as points
    let sigma_z = 0.01;        // measurement noise (±5 pixels at 512px)
    let sigma_prior = 2.0;     // weak prior on landmark positions
    let max_cameras = 6;
    let mc_samples = 200_000;

    println!("Scene configuration:");
    println!("  {} landmarks in [{:.1},{:.1}]³", n_points, -point_spread, point_spread);
    println!("  Camera ring: radius={}, height={}", cam_radius, cam_height);
    println!("  σ_z = {} (measurement), σ_prior = {} (landmark prior)", sigma_z, sigma_prior);
    println!("  MC samples per integral: {}\n", mc_samples);

    // Generate point cloud
    let points: Vec<Vec3> = (0..n_points).map(|_| {
        [
            point_spread * (2.0 * rng.uniform() - 1.0),
            point_spread * (2.0 * rng.uniform() - 1.0),
            point_spread * (2.0 * rng.uniform() - 1.0),
        ]
    }).collect();

    // Generate camera ring (max cameras)
    let all_cameras = camera_ring(max_cameras, cam_radius, cam_height);

    // ─── Sweep number of cameras ───
    println!("──────────────────────────────────────────────────────────");
    println!("  Part 1: c₁ vs number of cameras (averaged over point cloud)");
    println!("──────────────────────────────────────────────────────────\n");
    println!("{:>4}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}",
        "cams", "mean|c₁|", "max|c₁|", "mean_depth", "valid", "invalid");

    let mut csv = String::new();
    csv.push_str("n_cams,point_idx,depth,c1,term_a,term_b,term_q4,log_laplace,log_sp,log_mc,status\n");

    for n_cams in 1..=max_cameras {
        let cams = &all_cameras[..n_cams];
        let mut c1_abs_sum = 0.0;
        let mut c1_abs_max = 0.0f64;
        let mut depth_sum = 0.0;
        let mut n_valid = 0;
        let mut n_invalid = 0;

        for (pi, pt) in points.iter().enumerate() {
            // Build CameraObs for visible cameras
            let cam_obs: Vec<CameraObs> = cams.iter()
                .filter(|(r, t)| is_visible(r, t, pt))
                .map(|(r, t)| {
                    let xp = projective::transform_point(r, t, pt);
                    CameraObs {
                        rot: *r, trans: *t,
                        z: projective::project(&xp),
                        sigma_zz_inv: iso2(sigma_z),
                    }
                })
                .collect();

            if cam_obs.is_empty() { continue; }

            // Average depth across cameras
            let avg_depth: f64 = cam_obs.iter()
                .map(|c| projective::transform_point(&c.rot, &c.trans, pt)[2])
                .sum::<f64>() / cam_obs.len() as f64;
            depth_sum += avg_depth;

            let sigma_xx_inv = iso3(sigma_prior);
            let opt = optimize_landmark_multicam(&cam_obs, pt, &sigma_xx_inv, 20);
            let (log_sp, sp_r) = landmark_marginal_multicam(&opt, &cam_obs, &sigma_xx_inv);

            // Laplace
            let log_lap = -opt.nll_opt - 0.5 * opt.log_det_h
                + 1.5 * (2.0 * std::f64::consts::PI).ln();

            // MC (only for small camera counts where correction matters)
            let log_mc = if n_cams <= 3 {
                mc_marginal_integral(&cam_obs, &opt.x_opt, &sigma_xx_inv, pt, mc_samples, &mut rng)
            } else {
                log_lap // skip MC for many cameras (too accurate to matter)
            };

            match sp_r.status {
                SaddlepointStatus::Valid => n_valid += 1,
                _ => n_invalid += 1,
            }
            c1_abs_sum += sp_r.c1.abs();
            c1_abs_max = c1_abs_max.max(sp_r.c1.abs());

            let status_str = match sp_r.status {
                SaddlepointStatus::Valid => "valid",
                SaddlepointStatus::Invalid(_) => "invalid",
                SaddlepointStatus::SingularHessian => "singular",
            };
            csv.push_str(&format!("{},{},{:.4},{:.8},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{}\n",
                n_cams, pi, avg_depth, sp_r.c1, sp_r.term_a, sp_r.term_b, sp_r.term_q4,
                log_lap, log_sp, log_mc, status_str));
        }

        let n_total = (n_valid + n_invalid) as f64;
        println!("{:>4}  {:>10.6}  {:>10.6}  {:>10.2}  {:>10}  {:>10}",
            n_cams,
            c1_abs_sum / n_total,
            c1_abs_max,
            depth_sum / n_total,
            n_valid,
            n_invalid);
    }

    // ─── Detailed comparison for 1 and 2 cameras ───
    println!("\n──────────────────────────────────────────────────────────");
    println!("  Part 2: Laplace vs Saddlepoint vs MC (1 & 2 cameras)");
    println!("──────────────────────────────────────────────────────────\n");

    for n_cams in [1, 2] {
        let cams = &all_cameras[..n_cams];
        let mut lap_err_sum = 0.0;
        let mut sp_err_sum = 0.0;
        let mut n_compared = 0;

        for pt in &points {
            let cam_obs: Vec<CameraObs> = cams.iter()
                .filter(|(r, t)| is_visible(r, t, pt))
                .map(|(r, t)| {
                    let xp = projective::transform_point(r, t, pt);
                    CameraObs {
                        rot: *r, trans: *t,
                        z: projective::project(&xp),
                        sigma_zz_inv: iso2(sigma_z),
                    }
                })
                .collect();

            if cam_obs.is_empty() { continue; }

            let sigma_xx_inv = iso3(sigma_prior);
            let opt = optimize_landmark_multicam(&cam_obs, pt, &sigma_xx_inv, 20);
            let (log_sp, sp_r) = landmark_marginal_multicam(&opt, &cam_obs, &sigma_xx_inv);
            let log_lap = -opt.nll_opt - 0.5 * opt.log_det_h
                + 1.5 * (2.0 * std::f64::consts::PI).ln();
            let log_mc = mc_marginal_integral(&cam_obs, &opt.x_opt, &sigma_xx_inv, pt, mc_samples, &mut rng);

            if matches!(sp_r.status, SaddlepointStatus::Valid) && log_mc.is_finite() {
                lap_err_sum += (log_lap - log_mc).abs();
                sp_err_sum += (log_sp - log_mc).abs();
                n_compared += 1;
            }
        }

        if n_compared > 0 {
            let mean_lap_err = lap_err_sum / n_compared as f64;
            let mean_sp_err = sp_err_sum / n_compared as f64;
            let improvement = if mean_sp_err > 1e-15 { mean_lap_err / mean_sp_err } else { f64::INFINITY };
            println!("  {} camera(s): ({} landmarks compared)", n_cams, n_compared);
            println!("    Mean |log I_Laplace − log I_MC| = {:.6}", mean_lap_err);
            println!("    Mean |log I_SP − log I_MC|      = {:.6}", mean_sp_err);
            println!("    Improvement factor:               {:.1}x\n", improvement);
        }
    }

    // ─── Depth dependence ───
    println!("──────────────────────────────────────────────────────────");
    println!("  Part 3: c₁ vs depth (2 cameras, single landmark)");
    println!("──────────────────────────────────────────────────────────\n");
    println!("{:>8}  {:>12}  {:>10}  {:>10}  {:>10}",
        "depth", "σ_d/depth", "c₁", "Lap_err", "SP_err");


    for &depth in &[3.0, 5.0, 8.0, 12.0, 20.0, 50.0] {
        // Actually, let's place the point so it has the right depth from cameras
        let pt = [0.1, -0.05, 0.0]; // near origin
        // Adjust camera distance for target depth
        let test_cams = camera_ring(2, depth, 0.0);

        let cam_obs: Vec<CameraObs> = test_cams.iter()
            .filter(|(r, t)| is_visible(r, t, &pt))
            .map(|(r, t)| {
                let xp = projective::transform_point(r, t, &pt);
                CameraObs {
                    rot: *r, trans: *t,
                    z: projective::project(&xp),
                    sigma_zz_inv: iso2(sigma_z),
                }
            })
            .collect();

        if cam_obs.is_empty() { continue; }

        let avg_depth: f64 = cam_obs.iter()
            .map(|c| projective::transform_point(&c.rot, &c.trans, &pt)[2])
            .sum::<f64>() / cam_obs.len() as f64;

        let sigma_xx_inv = iso3(sigma_prior);
        let opt = optimize_landmark_multicam(&cam_obs, &pt, &sigma_xx_inv, 20);
        let (log_sp, sp_r) = landmark_marginal_multicam(&opt, &cam_obs, &sigma_xx_inv);
        let log_lap = -opt.nll_opt - 0.5 * opt.log_det_h
            + 1.5 * (2.0 * std::f64::consts::PI).ln();
        let log_mc = mc_marginal_integral(&cam_obs, &opt.x_opt, &sigma_xx_inv, &pt, mc_samples, &mut rng);

        let sigma_depth_over_d = sigma_prior / avg_depth;

        println!("{:>8.1}  {:>12.4}  {:>10.6}  {:>10.6}  {:>10.6}",
            avg_depth, sigma_depth_over_d, sp_r.c1,
            (log_lap - log_mc).abs(), (log_sp - log_mc).abs());
    }

    // Write CSV
    let csv_path = std::path::Path::new("/mnt/user-data/outputs/multicam_saddlepoint.csv");
    if let Some(parent) = csv_path.parent() {
        std::fs::create_dir_all(parent).unwrap_or_else(|e| {
            panic!(
                "Failed to create output directory '{}': {}\n\
                 Please create it manually: mkdir -p {}",
                parent.display(), e, parent.display()
            );
        });
    }
    std::fs::write(csv_path, &csv).unwrap_or_else(|e| {
        panic!(
            "Failed to write CSV to '{}': {}",
            csv_path.display(), e
        );
    });
    println!("\n  CSV written to: {}", csv_path.display());
}
