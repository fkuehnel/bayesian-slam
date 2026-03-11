Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on SE(3)
Author: Frank O. Kuehnel

Overview
This repository accompanies the paper Higher-Order Uncertainty Propagation and Saddlepoint Marginalization on the SE(3) Lie Group, which develops a self-contained framework for Bayesian inference over rigid body poses.
Paper: higherorder.pdf
The core ideas originate from unpublished work at NASA Ames Research Center (2008) on robust pose estimation using the SE(3) Lie group structure. This project refocuses that foundational material—stripping away the SLAM-specific framing—into three standalone contributions aimed at the broader estimation and inference community.
What This Paper Contributes
1. Second-Order Uncertainty Propagation on SE(3)
Current filters and optimizers (EKF, IEKF, factor graphs) propagate pose uncertainty using first-order Jacobians. On SE(3), the semi-direct product coupling between rotation and translation introduces systematic bias that first-order methods miss. We derive:

Closed-form composition Jacobians including the rotation-translation coupling Jacobian J_t(Ω, t)
Second-order corrections to both the mean and covariance of composed poses
Explicit formulas using the finite BCH composition via SU(2), avoiding truncated series

2. Saddlepoint Marginalization for Projective Observations
When camera observations are projective (x/z, y/z), the joint posterior over pose and 3D landmarks is non-Gaussian in the pose parameters. Standard Laplace approximation (Schur complement in bundle adjustment) misses this asymmetry. We develop:

A saddlepoint approximation that captures third-order cumulants of the projective model
A two-step marginalization scheme: optimize landmarks, then integrate with saddlepoint correction
Reduction from (6n + 3m)² joint system to 6n pose-only optimization with improved accuracy

3. Complete SE(3) Algebra Toolkit
The paper provides a self-contained reference for the SE(3) machinery, including:

Exponential/logarithmic maps with Lie-Cartan coordinates of the first kind
Closed-form finite Rodrigues vector composition via the SU(2)/Z₂ ≅ SO(3) double cover
Phase reflection handling at Θ = π with cutline analysis
Wei-Norman formula and both left/right Rodrigues Jacobians
The S⁻¹ identity connecting the translation coupling matrix to the rotation Jacobians

Publication Target
We are preparing this manuscript for submission to a robotics or estimation journal (IEEE Transactions on Robotics, International Journal of Robotics Research, or similar). The paper is structured as a methods contribution applicable across:

Visual-inertial odometry and IMU preintegration
Spacecraft attitude determination
Surgical robot registration
Multi-body dynamics and motion capture

Repository Structure
.
├── README.md                    # This file
├── paper/
│   ├── higherorder.tex          # Main manuscript (LaTeX)
│   └── robustEst.tex            # Original 2008 technical report (reference)
├── rust/
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs               # Library root
│       ├── so3.rs               # SO(3) operations: Rodrigues, exp/log, Jacobians
│       ├── se3.rs               # SE(3) operations: composition, action, S matrix
│       ├── bch.rs               # Finite BCH via SU(2) quaternions, phase reflection
│       ├── jacobians.rs         # J_ωr, J_ωl, J_t coupling Jacobian, S⁻¹ identities
│       ├── propagation.rs       # First and second-order uncertainty propagation
│       ├── saddlepoint.rs       # Saddlepoint marginalization for projective model
│       ├── projective.rs        # Projective observation model, derivatives to 3rd order
│       └── tests/
│           ├── bias_test.rs     # Reproduces bias experiment (L1 vs L2 coordinates)
│           ├── propagation_test.rs  # MC validation of 2nd-order propagation
│           └── saddlepoint_test.rs  # Saddlepoint vs Laplace vs joint MAP
└── experiments/
    ├── bias_parametrization.rs  # Experiment 1: scatter plots in L1 vs L2
    ├── propagation_accuracy.rs  # Experiment 2: 1st vs 2nd order vs MC truth
    └── marginalization.rs       # Experiment 3: convergence basin comparison

Rust Implementation
We provide a Rust implementation of the complete SE(3) inference framework. Rust is chosen for:

Correctness: The type system enforces dimensional consistency (SO(3) vs SE(3) vs ℝ³) at compile time
Performance: Zero-cost abstractions over the matrix algebra, suitable for real-time filtering
Safety: No silent numerical errors from uninitialized memory or index overflows in the Jacobian computations

Planned Modules
ModuleStatusDescriptionso3PlannedRodrigues formula, exp/log on SO(3), axis-angle representationse3PlannedSE(3) group operations, action on points, matrix representationbchPlannedClosed-form finite rotation composition via SU(2), limit cases, phase reflectionjacobiansPlannedJ_ωr, J_ωl, J_t, S⁻¹ identities, additive vs compositive derivativespropagationPlannedFirst-order covariance transport, second-order mean/covariance correctionssaddlepointPlannedProjective model cumulants, saddlepoint correction, marginalized MAPprojectivePlannedCamera model, calibration, derivatives through 3rd order

Design Principles

No dependency on nalgebra for core algebra: The 3×3 and 6×6 matrices arising in SE(3) have enough special structure (skew-symmetric, rotation, block-triangular) that specialized implementations outperform generic linear algebra
Const generics for compile-time dimensions: SO(3) Jacobians are [[f64; 3]; 3], SE(3) ones are [[f64; 6]; 6]—no heap allocation in the inner loop
Property-based testing: Jacobians are validated against finite differences; composition identities (S⁻¹R = J_ωr, etc.) are checked numerically

Getting Started
bash# Clone the repository
git clone https://github.com/[username]/se3-inference.git
cd se3-inference/rust

# Build
cargo build --release

# Run tests (once implemented)
cargo test

# Run bias experiment
cargo run --release --example bias_parametrization
Relationship to Prior Work
This project builds on:

Kuehnel (2004): "Bayesian estimation of nonlinear parameters on SE(3) Lie group," AIP Conf. Proc. vol. 735, pp. 176–186. DOI: 10.1063/1.1835212 — First introduction of Lie-Cartan exponential coordinates for unbiased Bayesian pose estimation and the marginalized MAP estimator
Kuehnel (2005): "Local frame junction trees in SLAM," AIP Conf. Proc. vol. 803, pp. 318–329. DOI: 10.1063/1.2149810 — Graphical model structure of SLAM via junction trees exploiting local frame dependence
Kuehnel (2008): Unpublished technical report at NASA Ames developing the full SE(3) algebra toolkit with BCH composition, Wei-Norman Jacobians, and phase reflection handling, cited as [kuehnel2008] in the manuscript
Solà, Deray & Atchuthan (2018): "A micro Lie theory for state estimation in robotics" — covers similar SE(3) machinery at first order; our work extends to second order
Ye & Chirikjian (2024): "Uncertainty propagation on unimodular Lie groups" — addresses propagation via SDEs; our approach uses the closed-form BCH directly
Barfoot (2024): State Estimation for Robotics — comprehensive textbook; our coupling Jacobian and saddlepoint marginalization go beyond what is covered there

License
Paper content: All rights reserved, Frank O. Kuehnel / Excel Solutions LLC.
Code: MIT License (see rust/LICENSE when available).
Contact
Frank O. Kuehnel — Excel Solutions LLC
Email: kuehnelg@gmail.com