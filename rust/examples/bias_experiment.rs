//! # Experiment 1: Bias in First-Kind vs Second-Kind Coordinates
//!
//! Paper §V.1, Figs. 1-2. Demonstrates that Lie-Cartan exponential
//! coordinates (first kind) produce unbiased pose estimates while
//! second-kind coordinates [Ω, T] show translational bias.
//!
//! ## Setup
//! - Camera at known pose f₀
//! - N landmarks randomly distributed around origin
//! - Gaussian noise added to landmark positions
//! - Pose estimated by minimizing Σ ||R·x + T - p||²
//! - Repeated for many noise realizations
//!
//! ## Output
//! - CSV files for scatter plots in both coordinate systems
//! - Statistics: mean and std in L1 and L2 coordinates
//! - Quantitative bias measurement

use se3_inference::*;
use se3_inference::so3;
use se3_inference::se3::Pose;

/// Simple LCG + Box-Muller PRNG
struct Rng {
    state: u64,
}

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

    fn normal3(&mut self, sigma: f64) -> Vec3 {
        [sigma * self.normal(), sigma * self.normal(), sigma * self.normal()]
    }
}

/// Solve point cloud alignment: find f = argmin Σ ||R·x_j + T - p_j||²
///
/// Uses Gauss-Newton on SE(3) with compositive updates.
/// x_j: source points, p_j: target points (noisy landmarks)
fn align_point_cloud(
    x_world: &[Vec3],  // landmarks in world frame
    p_obs: &[Vec3],    // observed (noisy) landmarks
    f_init: &Pose,     // initial guess
    max_iter: usize,
) -> Pose {
    let mut f = *f_init;
    let n = x_world.len();

    for _iter in 0..max_iter {
        // Build 6×6 normal equations: J^T J δ = -J^T r
        let mut jtj = [[0.0f64; 6]; 6];
        let mut jtr = [0.0f64; 6];

        for j in 0..n {
            // Residual: r_j = f ⋆ x_j - p_j = R·x_j + T - p_j
            let fx = f.act(&x_world[j]);
            let res = sub3(&fx, &p_obs[j]);

            // Jacobian of (f ⋆ x) w.r.t. right perturbation δf:
            //   d(f·exp(δ) ⋆ x)/dδ = R · J_×(x) = R · [-[x]×, I]
            // This is 3×6
            let jx = projective::j_cross(&x_world[j]); // 3×6
            let mut rjx = [[0.0; 6]; 3]; // R · J_×
            for i in 0..3 {
                for k in 0..6 {
                    rjx[i][k] = f.rot[i][0]*jx[0][k] + f.rot[i][1]*jx[1][k] + f.rot[i][2]*jx[2][k];
                }
            }

            // Accumulate J^T J and J^T r
            for a in 0..6 {
                for b in 0..6 {
                    for i in 0..3 {
                        jtj[a][b] += rjx[i][a] * rjx[i][b];
                    }
                }
                for i in 0..3 {
                    jtr[a] += rjx[i][a] * res[i];
                }
            }
        }

        // Solve 6×6 system: jtj · delta = -jtr
        // Using simple Cholesky
        let delta = solve_6x6(&jtj, &[-jtr[0], -jtr[1], -jtr[2], -jtr[3], -jtr[4], -jtr[5]]);

        // Compositive update
        f = f.compose(&Pose::exp(&delta));

        // Check convergence
        let delta_norm: f64 = delta.iter().map(|x| x*x).sum::<f64>().sqrt();
        if delta_norm < 1e-12 { break; }
    }

    f
}

/// Solve 6×6 linear system via Cholesky decomposition
fn solve_6x6(a: &Mat6, b: &Vec6) -> Vec6 {
    // Cholesky: A = L L^T
    let mut l = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j { sum += l[i][k] * l[j][k]; }
            if i == j {
                let diag = a[i][i] - sum;
                l[i][j] = if diag > 1e-30 { diag.sqrt() } else { 1e-15 };
            } else {
                l[i][j] = (a[i][j] - sum) / l[j][j];
            }
        }
    }

    // Forward solve: L y = b
    let mut y = [0.0f64; 6];
    for i in 0..6 {
        let mut sum = 0.0;
        for j in 0..i { sum += l[i][j] * y[j]; }
        y[i] = (b[i] - sum) / l[i][i];
    }

    // Back solve: L^T x = y
    let mut x = [0.0f64; 6];
    for i in (0..6).rev() {
        let mut sum = 0.0;
        for j in (i+1)..6 { sum += l[j][i] * x[j]; }
        x[i] = (y[i] - sum) / l[i][i];
    }
    x
}

fn main() {
    println!("==========================================================");
    println!("Experiment 1: Bias in First-Kind vs Second-Kind Coordinates");
    println!("==========================================================\n");

    // ─── Setup ───
    // Paper §V.1 config: N=3, σ=1.0 gives ~21× bias ratio (dramatic)
    // Mild config: N=8, σ=0.5 gives ~3× ratio
    // Usage: cargo run --example bias_experiment [mild]
    let args: Vec<String> = std::env::args().collect();
    let (n_landmarks, noise_sigma) = if args.len() > 1 && args[1] == "mild" {
        (8, 0.5)
    } else {
        (3, 1.0)
    };
    let n_samples = 2000;

    // True camera pose: looking from z=-10 toward origin
    let omega_true = [0.05, -0.03, 0.02]; // small rotation
    let t_true = [0.1, -0.05, -10.0];     // camera at z=-10
    let xi_true: Vec6 = [omega_true[0], omega_true[1], omega_true[2],
                          t_true[0], t_true[1], t_true[2]];
    let f_true = Pose::exp(&xi_true);

    println!("True pose (exp coords): ω=[{:.4}, {:.4}, {:.4}], t=[{:.4}, {:.4}, {:.4}]",
        xi_true[0], xi_true[1], xi_true[2], xi_true[3], xi_true[4], xi_true[5]);
    println!("True pose T = [{:.4}, {:.4}, {:.4}]",
        f_true.trans[0], f_true.trans[1], f_true.trans[2]);
    println!("Landmarks: {}, Noise σ: {}, Samples: {}\n", n_landmarks, noise_sigma, n_samples);

    // Generate landmarks in world frame (cube around origin)
    let mut rng = Rng::new(12345);
    let mut landmarks = Vec::with_capacity(n_landmarks);
    for _ in 0..n_landmarks {
        landmarks.push([
            2.0 * rng.uniform() - 1.0,
            2.0 * rng.uniform() - 1.0,
            2.0 * rng.uniform() - 1.0,
        ]);
    }

    // Compute true observations: p_j = f_true ⋆ x_j
    let true_obs: Vec<Vec3> = landmarks.iter().map(|x| f_true.act(x)).collect();

    // ─── Monte Carlo sampling ───
    let mut l1_samples: Vec<Vec6> = Vec::with_capacity(n_samples);  // exp coords [Ω, t]
    let mut l2_omega_samples: Vec<Vec3> = Vec::with_capacity(n_samples);  // Ω - Ω₀
    let mut l2_trans_samples: Vec<Vec3> = Vec::with_capacity(n_samples);  // T - T₀

    let mut n_converged = 0;

    for s in 0..n_samples {
        // Add noise to observations
        let noisy_obs: Vec<Vec3> = true_obs.iter().map(|p| {
            add3(p, &rng.normal3(noise_sigma))
        }).collect();

        // Solve alignment problem (initialize at true pose)
        let f_est = align_point_cloud(&landmarks, &noisy_obs, &f_true, 50);

        // Check convergence (residual should be reasonable)
        let residual: f64 = landmarks.iter().zip(noisy_obs.iter())
            .map(|(x, p)| {
                let fx = f_est.act(x);
                let r = sub3(&fx, p);
                dot3(&r, &r)
            })
            .sum::<f64>() / n_landmarks as f64;

        if residual > 100.0 * noise_sigma * noise_sigma {
            continue; // skip non-converged
        }
        n_converged += 1;

        // L1 coordinates: ln(f_true⁻¹ · f_est) = [δΩ, δt]
        let delta_f = f_true.relative(&f_est);
        let xi_delta = delta_f.log();
        l1_samples.push(xi_delta);

        // L2 coordinates: [Ω_est - Ω_true, T_est - T_true]
        let omega_est = so3::log(&f_est.rot);
        let omega_diff = sub3(&omega_est, &omega_true);
        let trans_diff = sub3(&f_est.trans, &f_true.trans);
        l2_omega_samples.push(omega_diff);
        l2_trans_samples.push(trans_diff);

        if (s + 1) % 500 == 0 {
            eprint!("  {} / {} samples...\r", s + 1, n_samples);
        }
    }
    eprintln!();

    println!("Converged: {} / {} ({:.1}%)\n",
        n_converged, n_samples, 100.0 * n_converged as f64 / n_samples as f64);

    // ─── Statistics ───
    let n = l1_samples.len() as f64;

    // L1 mean (exponential coordinates)
    let mut l1_mean = [0.0f64; 6];
    for s in &l1_samples {
        for i in 0..6 { l1_mean[i] += s[i]; }
    }
    for i in 0..6 { l1_mean[i] /= n; }

    // L2 means
    let mut l2_omega_mean = [0.0; 3];
    let mut l2_trans_mean = [0.0; 3];
    for s in &l2_omega_samples { for i in 0..3 { l2_omega_mean[i] += s[i]; } }
    for s in &l2_trans_samples { for i in 0..3 { l2_trans_mean[i] += s[i]; } }
    for i in 0..3 { l2_omega_mean[i] /= n; l2_trans_mean[i] /= n; }

    // Standard deviations
    let mut l1_std = [0.0f64; 6];
    for s in &l1_samples {
        for i in 0..6 { l1_std[i] += (s[i] - l1_mean[i]).powi(2); }
    }
    for i in 0..6 { l1_std[i] = (l1_std[i] / (n - 1.0)).sqrt(); }

    let mut l2_omega_std = [0.0; 3];
    let mut l2_trans_std = [0.0; 3];
    for s in &l2_omega_samples {
        for i in 0..3 { l2_omega_std[i] += (s[i] - l2_omega_mean[i]).powi(2); }
    }
    for s in &l2_trans_samples {
        for i in 0..3 { l2_trans_std[i] += (s[i] - l2_trans_mean[i]).powi(2); }
    }
    for i in 0..3 {
        l2_omega_std[i] = (l2_omega_std[i] / (n - 1.0)).sqrt();
        l2_trans_std[i] = (l2_trans_std[i] / (n - 1.0)).sqrt();
    }

    // ─── Report ───
    println!("══════════════════════════════════════════════════════");
    println!("  FIRST-KIND (Exponential) Coordinates: ln(f₀⁻¹·f)");
    println!("══════════════════════════════════════════════════════");
    println!("  Rotation mean:    [{:+.6}, {:+.6}, {:+.6}]",
        l1_mean[0], l1_mean[1], l1_mean[2]);
    println!("  Rotation std:     [{:.6}, {:.6}, {:.6}]",
        l1_std[0], l1_std[1], l1_std[2]);
    println!("  Translation mean: [{:+.6}, {:+.6}, {:+.6}]",
        l1_mean[3], l1_mean[4], l1_mean[5]);
    println!("  Translation std:  [{:.6}, {:.6}, {:.6}]",
        l1_std[3], l1_std[4], l1_std[5]);
    let l1_rot_bias = (l1_mean[0]*l1_mean[0] + l1_mean[1]*l1_mean[1] + l1_mean[2]*l1_mean[2]).sqrt();
    let l1_trans_bias = (l1_mean[3]*l1_mean[3] + l1_mean[4]*l1_mean[4] + l1_mean[5]*l1_mean[5]).sqrt();
    println!("  |rotation bias|:    {:.6}", l1_rot_bias);
    println!("  |translation bias|: {:.6}", l1_trans_bias);

    println!();
    println!("══════════════════════════════════════════════════════");
    println!("  SECOND-KIND (Additive) Coordinates: [Ω-Ω₀, T-T₀]");
    println!("══════════════════════════════════════════════════════");
    println!("  Rotation mean:    [{:+.6}, {:+.6}, {:+.6}]",
        l2_omega_mean[0], l2_omega_mean[1], l2_omega_mean[2]);
    println!("  Rotation std:     [{:.6}, {:.6}, {:.6}]",
        l2_omega_std[0], l2_omega_std[1], l2_omega_std[2]);
    println!("  Translation mean: [{:+.6}, {:+.6}, {:+.6}]",
        l2_trans_mean[0], l2_trans_mean[1], l2_trans_mean[2]);
    println!("  Translation std:  [{:.6}, {:.6}, {:.6}]",
        l2_trans_std[0], l2_trans_std[1], l2_trans_std[2]);
    let l2_rot_bias = norm3(&l2_omega_mean);
    let l2_trans_bias = norm3(&l2_trans_mean);
    println!("  |rotation bias|:    {:.6}", l2_rot_bias);
    println!("  |translation bias|: {:.6}", l2_trans_bias);

    println!();
    println!("══════════════════════════════════════════════════════");
    println!("  BIAS COMPARISON");
    println!("══════════════════════════════════════════════════════");
    let bias_ratio_rot = if l1_rot_bias > 1e-15 { l2_rot_bias / l1_rot_bias } else { f64::INFINITY };
    let bias_ratio_trans = if l1_trans_bias > 1e-15 { l2_trans_bias / l1_trans_bias } else { f64::INFINITY };
    println!("  Rotation bias ratio (L2/L1):    {:.2}x", bias_ratio_rot);
    println!("  Translation bias ratio (L2/L1): {:.2}x", bias_ratio_trans);
    println!();
    if l2_trans_bias > 2.0 * l1_trans_bias {
        println!("  ✓ Second-kind coordinates show significantly MORE translational bias");
        println!("    than first-kind (exponential) coordinates.");
        println!("    This confirms Proposition 1 (paper §II.3).");
    } else {
        println!("  Note: bias difference may be small at this noise level.");
        println!("  Try increasing noise_sigma or reducing n_landmarks.");
    }

    // ─── Write CSV for plotting ───
    let csv_path = std::path::Path::new("/mnt/user-data/outputs/bias_experiment.csv");
    if let Some(parent) = csv_path.parent() {
        std::fs::create_dir_all(parent).unwrap_or_else(|e| {
            panic!(
                "Failed to create output directory '{}': {}\n\
                 Please create it manually: mkdir -p {}",
                parent.display(), e, parent.display()
            );
        });
    }
    let mut csv = String::new();
    csv.push_str("l1_omega1,l1_omega2,l1_omega3,l1_t1,l1_t2,l1_t3,l2_domega1,l2_domega2,l2_domega3,l2_dt1,l2_dt2,l2_dt3\n");
    for i in 0..l1_samples.len() {
        csv.push_str(&format!("{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8}\n",
            l1_samples[i][0], l1_samples[i][1], l1_samples[i][2],
            l1_samples[i][3], l1_samples[i][4], l1_samples[i][5],
            l2_omega_samples[i][0], l2_omega_samples[i][1], l2_omega_samples[i][2],
            l2_trans_samples[i][0], l2_trans_samples[i][1], l2_trans_samples[i][2],
        ));
    }
    std::fs::write(csv_path, &csv).unwrap_or_else(|e| {
        panic!(
            "Failed to write CSV to '{}': {}",
            csv_path.display(), e
        );
    });
    println!("\n  CSV written to: {}", csv_path.display());
}
