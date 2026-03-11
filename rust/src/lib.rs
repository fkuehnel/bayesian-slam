//! # se3-inference
//!
//! Higher-order uncertainty propagation and saddlepoint marginalization
//! on the SE(3) Lie group.
//!
//! All matrices are stack-allocated fixed-size arrays — no heap, no deps.

pub mod so3;
pub mod se3;
pub mod jacobians;

/// 3×3 matrix, row-major.
pub type Mat3 = [[f64; 3]; 3];
/// 6×6 matrix, row-major.
pub type Mat6 = [[f64; 6]; 6];
/// 3-vector.
pub type Vec3 = [f64; 3];
/// 6-vector: [ω₁, ω₂, ω₃, v₁, v₂, v₃].
pub type Vec6 = [f64; 6];

pub const I3: Mat3 = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
pub const Z3: Mat3 = [[0.0;3];3];

pub const I6: Mat6 = {
    let mut m = [[0.0f64; 6]; 6];
    let mut i = 0;
    while i < 6 { m[i][i] = 1.0; i += 1; }
    m
};

// ─── Vec3 ops ───

#[inline] pub fn dot3(a: &Vec3, b: &Vec3) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}
#[inline] pub fn cross3(a: &Vec3, b: &Vec3) -> Vec3 {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
#[inline] pub fn norm3(v: &Vec3) -> f64 { dot3(v,v).sqrt() }
#[inline] pub fn scale3(s: f64, v: &Vec3) -> Vec3 { [s*v[0], s*v[1], s*v[2]] }
#[inline] pub fn add3(a: &Vec3, b: &Vec3) -> Vec3 { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
#[inline] pub fn sub3(a: &Vec3, b: &Vec3) -> Vec3 { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }

/// Outer product a bᵀ.
#[inline] pub fn outer3(a: &Vec3, b: &Vec3) -> Mat3 {
    [[a[0]*b[0],a[0]*b[1],a[0]*b[2]],
     [a[1]*b[0],a[1]*b[1],a[1]*b[2]],
     [a[2]*b[0],a[2]*b[1],a[2]*b[2]]]
}

// ─── Mat3 ops ───

#[inline] pub fn mv3(m: &Mat3, v: &Vec3) -> Vec3 {
    [m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2],
     m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2],
     m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]]
}
#[inline] pub fn mm3(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = Z3;
    for i in 0..3 { for j in 0..3 {
        c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
    }}
    c
}
#[inline] pub fn transpose3(m: &Mat3) -> Mat3 {
    [[m[0][0],m[1][0],m[2][0]],
     [m[0][1],m[1][1],m[2][1]],
     [m[0][2],m[1][2],m[2][2]]]
}
#[inline] pub fn add_mat3(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = Z3;
    for i in 0..3 { for j in 0..3 { c[i][j] = a[i][j]+b[i][j]; }}
    c
}
#[inline] pub fn sub_mat3(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = Z3;
    for i in 0..3 { for j in 0..3 { c[i][j] = a[i][j]-b[i][j]; }}
    c
}
#[inline] pub fn scale_mat3(s: f64, m: &Mat3) -> Mat3 {
    let mut c = Z3;
    for i in 0..3 { for j in 0..3 { c[i][j] = s*m[i][j]; }}
    c
}
#[inline] pub fn trace3(m: &Mat3) -> f64 { m[0][0]+m[1][1]+m[2][2] }

/// Determinant of 3×3 matrix.
pub fn det3(m: &Mat3) -> f64 {
    m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
  - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
  + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0])
}

/// Inverse of 3×3 matrix (panics if singular).
pub fn inv3(m: &Mat3) -> Mat3 {
    let d = det3(m);
    assert!(d.abs() > 1e-15, "inv3: singular matrix, det={:.2e}", d);
    let id = 1.0 / d;
    [
        [(m[1][1]*m[2][2]-m[1][2]*m[2][1])*id, (m[0][2]*m[2][1]-m[0][1]*m[2][2])*id, (m[0][1]*m[1][2]-m[0][2]*m[1][1])*id],
        [(m[1][2]*m[2][0]-m[1][0]*m[2][2])*id, (m[0][0]*m[2][2]-m[0][2]*m[2][0])*id, (m[0][2]*m[1][0]-m[0][0]*m[1][2])*id],
        [(m[1][0]*m[2][1]-m[1][1]*m[2][0])*id, (m[0][1]*m[2][0]-m[0][0]*m[2][1])*id, (m[0][0]*m[1][1]-m[0][1]*m[1][0])*id],
    ]
}

/// Frobenius norm of 3×3 matrix.
pub fn frob3(m: &Mat3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 { for j in 0..3 { s += m[i][j]*m[i][j]; } }
    s.sqrt()
}
