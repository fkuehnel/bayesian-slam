//! # Experiment: Multi-Camera Saddlepoint Marginalization
//!
//! 2+ cameras observing a shared 3D point cloud. Uses the native
//! `optimize_landmark_multicam` / `landmark_marginal_multicam` API.
//!
//! Usage:
//!   cargo run --release --example multicam_saddlepoint [N_cameras] [N_landmarks]

use se3_inference::*;
use se3_inference::se3::Pose;
use se3_inference::saddlepoint::*;
use se3_inference::projective;

// ── tiny PRNG ──────────────────────────────────────────────────────
struct Rng(u64);
impl Rng {
    fn new(s: u64) -> Self { Rng(s) }
    fn u(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn n(&mut self) -> f64 {
        let u1 = self.u().max(1e-30);
        let u2 = self.u();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 { lo + (hi - lo) * self.u() }
}

// ── camera ring ────────────────────────────────────────────────────
fn ring_cameras(n: usize, radius: f64, center: &Vec3) -> Vec<Pose> {
    (0..n).map(|i| {
        let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let h = 0.5 * (if i % 2 == 0 { 1.0 } else { -1.0 });
        let pos = [center[0] + radius * a.cos(), center[1] + h, center[2] + radius * a.sin()];
        let z = { let d = sub3(center, &pos); scale3(1.0 / norm3(&d), &d) };
        let xr = cross3(&[0.0, 1.0, 0.0], &z);
        let x = if norm3(&xr) > 1e-6 { scale3(1.0 / norm3(&xr), &xr) } else { [1.0, 0.0, 0.0] };
        let y = cross3(&z, &x);
        let r = [[x[0],x[1],x[2]], [y[0],y[1],y[2]], [z[0],z[1],z[2]]];
        Pose::new(r, scale3(-1.0, &mv3(&r, &pos)))
    }).collect()
}

fn visible(cam: &Pose, x: &Vec3, ht: f64) -> bool {
    let p = cam.act(x);
    p[2] > 0.5 && (p[0]/p[2]).abs() < ht && (p[1]/p[2]).abs() < ht
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n_cam: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(4);
    let n_pts: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20);

    let radius = 8.0;
    let sigma_z = 0.005;          // ~2.5 px at 512
    let sigma_prior = 2.0;
    let fov: f64 = 60.0;
    let ht = (fov.to_radians() / 2.0).tan();

    let sz_inv = { let v = 1.0/(sigma_z*sigma_z); [[v,0.0],[0.0,v]] };
    let sx_inv = { let v = 1.0/(sigma_prior*sigma_prior);
        [[v,0.0,0.0],[0.0,v,0.0],[0.0,0.0,v]] };

    println!("═══════════════════════════════════════════════════════════════");
    println!("  Multi-Camera Saddlepoint Marginalization");
    println!("═══════════════════════════════════════════════════════════════\n");
    println!("  {} cameras (ring, r={:.1}m) | {} landmarks in [-2,2]³", n_cam, radius, n_pts);
    println!("  σ_z={:.4}  σ_prior={:.1}  FOV={:.0}°\n", sigma_z, sigma_prior, fov);

    let cameras = ring_cameras(n_cam, radius, &[0.0; 3]);
    for (i, c) in cameras.iter().enumerate() {
        let p = c.inverse().trans;
        println!("  Cam {}: [{:+.2}, {:+.2}, {:+.2}]", i, p[0], p[1], p[2]);
    }

    let mut rng = Rng::new(42);
    let landmarks: Vec<Vec3> = (0..n_pts)
        .map(|_| [rng.range(-2.0,2.0), rng.range(-2.0,2.0), rng.range(-2.0,2.0)])
        .collect();

    // ═══════════════════════════════════════════════════════════
    //  Per-landmark analysis
    // ═══════════════════════════════════════════════════════════
    println!("\n{:>4} {:>5} {:>7} {:>11} {:>11} {:>11} {:>11} {:>4}",
        "LM", "views", "depth", "c₁", "A/12", "B/8", "-Q₄/8", "ok");
    println!("{}", "─".repeat(68));

    let mut sum_lap = 0.0f64;
    let mut sum_sp  = 0.0f64;
    let mut nv = 0usize;
    let mut ni = 0usize;
    let mut records: Vec<(usize, usize, f64, f64, bool)> = Vec::new(); // lm, nviews, depth, c1, valid

    for (li, xw) in landmarks.iter().enumerate() {
        // Build CameraObs for every camera that sees this landmark
        let mut obs: Vec<CameraObs> = Vec::new();
        let mut min_depth = f64::MAX;
        for cam in &cameras {
            if visible(cam, xw, ht) {
                let xp = cam.act(xw);
                min_depth = min_depth.min(xp[2]);
                let pi = projective::project(&xp);
                obs.push(CameraObs {
                    rot: cam.rot,
                    trans: cam.trans,
                    z: [pi[0] + sigma_z*rng.n(), pi[1] + sigma_z*rng.n()],
                    sigma_zz_inv: sz_inv,
                });
            }
        }
        if obs.is_empty() { continue; }

        let opt = optimize_landmark_multicam(&obs, xw, &sx_inv, 15);
        let lap = -opt.nll_opt - 0.5*opt.log_det_h
            + 1.5*(2.0*std::f64::consts::PI).ln();
        let (sp, res) = landmark_marginal_multicam(&opt, &obs);

        sum_lap += lap;
        sum_sp  += sp;

        let ok = matches!(res.status, SaddlepointStatus::Valid);
        if ok { nv += 1; } else { ni += 1; }
        records.push((li, obs.len(), min_depth, res.c1, ok));

        println!("{:>4} {:>5} {:>7.2} {:>11.3e} {:>11.3e} {:>11.3e} {:>11.3e} {:>4}",
            li, obs.len(), min_depth, res.c1,
            res.term_a/12.0, res.term_b/8.0, -res.term_q4/8.0,
            if ok { "✓" } else { "✗" });
    }

    // ═══════════════════════════════════════════════════════════
    //  Summary
    // ═══════════════════════════════════════════════════════════
    println!("\n═══════════════════════════════════════════════════════════════");
    println!("  Summary");
    println!("═══════════════════════════════════════════════════════════════");
    let total = nv + ni;
    println!("  Valid: {} / {} ({:.0}%)", nv, total, 100.0*nv as f64/total as f64);

    let vc: Vec<f64> = records.iter().filter(|r| r.4).map(|r| r.3).collect();
    if !vc.is_empty() {
        let mean: f64 = vc.iter().sum::<f64>() / vc.len() as f64;
        let mx = vc.iter().fold(0.0f64, |a,b| a.max(b.abs()));
        println!("  c₁  mean: {:+.4e}   |max|: {:.4e}", mean, mx);
    }
    println!("  Σ log-Laplace:      {:.6}", sum_lap);
    println!("  Σ log-Saddlepoint:  {:.6}", sum_sp);
    println!("  Correction Σ|Δ|:    {:.4e}", (sum_sp - sum_lap).abs());

    // ═══════════════════════════════════════════════════════════
    //  Camera-count sweep  (same test point, varying cameras)
    // ═══════════════════════════════════════════════════════════
    println!("\n═══════════════════════════════════════════════════════════════");
    println!("  Camera-Count Sweep  (point at [0.5, -0.3, 0.0])");
    println!("═══════════════════════════════════════════════════════════════\n");
    println!("{:>6} {:>6} {:>9} {:>9} {:>12}", "N_cam", "views", "σ_depth", "σ/depth", "c₁");
    println!("{}", "─".repeat(48));

    let tp = [0.5, -0.3, 0.0];
    for nc in &[2usize, 3, 4, 6, 8, 12] {
        let cams = ring_cameras(*nc, radius, &[0.0;3]);
        let mut obs: Vec<CameraObs> = Vec::new();
        let mut md = f64::MAX;
        for cam in &cams {
            if visible(cam, &tp, ht) {
                let xp = cam.act(&tp);
                md = md.min(xp[2]);
                let pi = projective::project(&xp);
                obs.push(CameraObs { rot: cam.rot, trans: cam.trans,
                    z: pi, sigma_zz_inv: sz_inv });
            }
        }
        if obs.is_empty() { println!("{:>6} {:>6}", nc, 0); continue; }

        let opt = optimize_landmark_multicam(&obs, &tp, &sx_inv, 15);
        let (_, res) = landmark_marginal_multicam(&opt, &obs);
        let hi = inv3(&opt.hessian);
        let sd = hi[2][2].abs().sqrt();

        let c1s = match res.status {
            SaddlepointStatus::Valid => format!("{:>12.4e}", res.c1),
            _ => format!("{:>12}", "invalid"),
        };
        println!("{:>6} {:>6} {:>9.5} {:>9.5} {}", nc, obs.len(), sd, sd/md, c1s);
    }
    println!("\n  More cameras → smaller σ/depth → smaller c₁");
    println!("  Sweet spot: 2-4 cameras (typical stereo/multi-view regime)\n");

    // CSV
    let csv_path = "../experiments/data/multicam_saddlepoint.csv";
    let mut csv = String::from("landmark,n_views,min_depth,c1,valid\n");
    for (li, nv, d, c, ok) in &records {
        csv.push_str(&format!("{},{},{:.6},{:.8e},{}\n", li, nv, d, c, ok));
    }
    std::fs::write(csv_path, &csv).unwrap_or_else(|e| eprintln!("CSV write: {}", e));
    println!("  CSV written to {}", csv_path);
}
