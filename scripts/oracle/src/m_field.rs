//! `field.tsv` — rust-old/src/field.rs. M1 part 3a: `FieldVector`, the point-charge fields and potential, the
//! superposition, the energy densities and the dipoles. Part 3b: the Gauss's-law fields, the multipoles, the images,
//! the Green's functions and the stress tensor. Part 3c: `trace_field_line`.
//!
//! Rows are class A where a function is + − × ÷, sqrt, copysign and clamp in a fixed order, which the port repeats;
//! the multipole moments (atan2, cos, sin) and the scalar Green's function (cos, sin) are class C, against the
//! corrected oracle, and the scalar Green's function is mapped by ADR 0002 row 5 (the port's value is the conjugate:
//! e^{+jωt}). A 3×3 tensor is nine `f:` tokens, row-major. A vector is
//! three `f:` tokens; a position `[f64; 3]` likewise. A `Result` is its value's tokens on `Ok`, or `E:<Variant>` and
//! the error's Display text on `Err` (ADR 0002 row 4: the port reads both back from the last-error record). The
//! superposition's charges are `v:<4n>` then q, x, y, z per charge.
//!
//! The grids: the specials crossed into vectors, components over every 37th binade with both signs (the rescaled
//! path of `normalized` among them), seeded vectors and SI-scale charges; positions on a lattice, at random, and at
//! offsets from the source around the singular cutoff r² = 1e-30 (√1e-30 = 1e-15 m) on both sides.

use crate::grids::{knife_edge, Rng, pow2, specials};
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::BijliError;
use bijli::field::{
    FieldVector, electric_dipole_field, electric_dipole_potential, electric_energy_density,
    electric_field_charged_sphere, electric_field_disk_axis, electric_field_infinite_line,
    electric_field_infinite_plane, electric_field_inside_cylinder, electric_field_inside_sphere,
    electric_field_outside_cylinder, electric_field_point_charge, electric_field_ring_axis,
    electric_field_superposition, electric_potential_point_charge, electric_quadrupole_potential,
    greens_function_scalar, greens_function_static, image_charge_dielectric, image_charge_plane,
    image_charge_sphere, magnetic_dipole_field, magnetic_energy_density,
    magnetic_field_moving_charge, maxwell_stress_tensor, multipole_moment, poynting_flux_surface,
    spherical_multipole_moment, stress_tensor_force, trace_field_line,
};
use std::hint::black_box;

pub(crate) fn v3(v: FieldVector) -> Vec<Tok> {
    vec![f(v.x), f(v.y), f(v.z)]
}

pub(crate) fn p3(p: [f64; 3]) -> Vec<Tok> {
    vec![f(p[0]), f(p[1]), f(p[2])]
}

pub(crate) fn res_v(r: bijli::error::Result<FieldVector>) -> Vec<Tok> {
    match r {
        Ok(v) => v3(v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

pub(crate) fn res_f(r: bijli::error::Result<f64>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn res_pair(r: bijli::error::Result<(f64, f64)>) -> Vec<Tok> {
    match r {
        Ok((a, b)) => vec![f(a), f(b)],
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn res_image(r: bijli::error::Result<(f64, [f64; 3])>) -> Vec<Tok> {
    match r {
        Ok((q, p)) => cat(&[vec![f(q)], p3(p)]),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// Distances and radii: the specials, negatives, the cutoff 1e-30 from both sides, binades, seeded.
pub(crate) fn distances(rng: &mut Rng) -> Vec<f64> {
    let mut ds = specials();
    ds.extend([
        -2.0, -1e-300, 0.99e-30, 1e-30, 1.01e-30, 1e-15, 0.25, 0.5, 1.0, 2.0, 1e8, 1e200,
    ]);
    for k in (-1074..=1023).step_by(41) {
        ds.push(pow2(k));
    }
    for _ in 0..60 {
        ds.push(rng.log_uniform(-12.0, 6.0));
    }
    ds
}

/// Charge sets for the multipole moments: none, one at the origin, one near the pole, and 1–4 seeded charges (one at
/// the origin or on the z axis in some).
fn charge_sets(rng: &mut Rng) -> Vec<Vec<(f64, [f64; 3])>> {
    let mut sets = vec![
        vec![],
        vec![(4.0, [0.0; 3])],
        vec![(2e-9, [1e-9, 0.0, 1.0])],
    ];
    for n in 0..30 {
        let count = 1 + n % 4;
        let mut cs: Vec<(f64, [f64; 3])> = (0..count)
            .map(|_| {
                (
                    rng.signed_log_uniform(-10.0, -8.0),
                    [
                        rng.uniform(-1.0, 1.0),
                        rng.uniform(-1.0, 1.0),
                        rng.uniform(-1.0, 1.0),
                    ],
                )
            })
            .collect();
        if n % 7 == 3 {
            cs[0].1 = [0.0; 3];
        }
        if n % 5 == 1 {
            cs[0].1 = [0.0, 0.0, rng.uniform(-1.0, 1.0)];
        }
        sets.push(cs);
    }
    sets
}

pub(crate) fn charges_tok(cs: &[(f64, [f64; 3])]) -> Vec<Tok> {
    let mut v = vec![Tok::V(4 * cs.len())];
    for (q, p) in cs {
        v.push(f(*q));
        v.extend(p3(*p));
    }
    v
}

pub(crate) fn cat(parts: &[Vec<Tok>]) -> Vec<Tok> {
    parts.iter().flatten().cloned().collect()
}

/// Vectors: the specials crossed (x, y, 0) and (0, y, z), every 37th binade, and seeded.
pub(crate) fn vectors() -> Vec<FieldVector> {
    let sp = specials();
    let mut vs = Vec::new();
    for &a in &sp {
        for &b in &sp {
            vs.push(FieldVector::new(a, b, 0.0));
            vs.push(FieldVector::new(0.0, a, b));
        }
    }
    let mut rng = Rng::new(0xf1e1d);
    for k in (-1074..=1023).step_by(37) {
        let x = pow2(k);
        vs.push(FieldVector::new(3.0 * x, -4.0 * x, 12.0 * x));
        vs.push(FieldVector::new(x, x, x));
        vs.push(FieldVector::new(
            -x,
            0.0,
            rng.signed_log_uniform(-3.0, 3.0) * x,
        ));
    }
    for _ in 0..300 {
        vs.push(FieldVector::new(
            rng.uniform(-10.0, 10.0),
            rng.uniform(-10.0, 10.0),
            rng.uniform(-10.0, 10.0),
        ));
    }
    for _ in 0..100 {
        vs.push(FieldVector::new(
            rng.signed_log_uniform(-35.0, 35.0),
            rng.signed_log_uniform(-35.0, 35.0),
            rng.signed_log_uniform(-35.0, 35.0),
        ));
    }
    vs
}

/// (source, field point) pairs: a lattice, random points, the cutoff from both sides, specials.
pub(crate) fn pairs(rng: &mut Rng) -> Vec<([f64; 3], [f64; 3])> {
    let mut ps = Vec::new();
    for i in -2..=2 {
        for j in -2..=2 {
            for k in -1..=1 {
                ps.push(([0.0, 0.0, 0.0], [f64::from(i), f64::from(j), f64::from(k)]));
            }
        }
    }
    for _ in 0..200 {
        let a = [
            rng.uniform(-5.0, 5.0),
            rng.uniform(-5.0, 5.0),
            rng.uniform(-5.0, 5.0),
        ];
        let b = [
            rng.uniform(-5.0, 5.0),
            rng.uniform(-5.0, 5.0),
            rng.uniform(-5.0, 5.0),
        ];
        ps.push((a, b));
    }
    // Offsets d·(1, 0, 0) and d·(1, 1, 1)/√3-ish around |d| = 1e-15 (r² = 1e-30), and far.
    for d in [
        1e-16, 9.999e-16, 1e-15, 1.0001e-15, 2e-15, 1e-10, 1e10, 1e150, 1e160,
    ] {
        let base = [1.5, -0.25, 3.0];
        ps.push((base, [base[0] + d, base[1], base[2]]));
        ps.push(([0.0; 3], [d, d, d]));
        ps.push(([0.0; 3], [0.0, -d, 0.0]));
    }
    for x in specials() {
        ps.push(([0.0; 3], [x, 1.0, 0.0]));
        ps.push(([x, 0.0, 0.0], [1.0, 2.0, 3.0]));
    }
    ps
}

/// Charges: SI scales, zero, the specials.
pub(crate) fn charges() -> Vec<f64> {
    let mut qs = vec![1e-6, -1e-6, 1.602_176_634e-19, -3e-9, 1.0, 0.0, -0.0, 1e300];
    qs.extend([f64::NAN, f64::INFINITY]);
    qs
}

pub fn emit() -> Module {
    let mut m = Module::new("field");
    let vs = vectors();
    let mut rng = Rng::new(0xf1e1d2);
    for (n, &v) in vs.iter().enumerate() {
        let ins = v3(v);
        let w = vs[(n * 7 + 3) % vs.len()];
        let both = cat(&[v3(v), v3(w)]);
        let x = [2.0, -0.5, 1e-300, 1e300, 0.0, f64::NAN, 3.7][n % 7];
        m.row(
            "FieldVector::new",
            Class::A,
            &ins,
            &v3(FieldVector::new(black_box(v.x), v.y, v.z)),
        );
        m.row(
            "FieldVector::magnitude",
            Class::A,
            &ins,
            &[f(black_box(v).magnitude())],
        );
        m.row(
            "FieldVector::magnitude_sq",
            Class::A,
            &ins,
            &[f(black_box(v).magnitude_sq())],
        );
        m.row(
            "FieldVector::normalized",
            Class::A,
            &ins,
            &res_v(black_box(v).normalized()),
        );
        m.row(
            "FieldVector::dot",
            Class::A,
            &both,
            &[f(black_box(v).dot(&w))],
        );
        m.row(
            "FieldVector::cross",
            Class::A,
            &both,
            &v3(black_box(v).cross(&w)),
        );
        m.row("FieldVector::add", Class::A, &both, &v3(black_box(v) + w));
        m.row("FieldVector::sub", Class::A, &both, &v3(black_box(v) - w));
        m.row("FieldVector::neg", Class::A, &ins, &v3(-black_box(v)));
        let vx = cat(&[v3(v), vec![f(x)]]);
        m.row(
            "FieldVector::scale",
            Class::A,
            &vx,
            &v3(black_box(v).scale(x)),
        );
        m.row("FieldVector::mul", Class::A, &vx, &v3(black_box(v) * x));
        m.row(
            "f64::mul_field_vector",
            Class::A,
            &cat(&[vec![f(x)], v3(v)]),
            &v3(black_box(x) * v),
        );
        m.row(
            "FieldVector::fmt",
            Class::A,
            &ins,
            &[s(format!("{}", black_box(v)))],
        );
        if n % 5 == 0 {
            m.row(
                "FieldVector::from_array",
                Class::A,
                &ins,
                &v3(FieldVector::from(black_box([v.x, v.y, v.z]))),
            );
            let arr: [f64; 3] = black_box(v).into();
            m.row("FieldVector::into_array", Class::A, &ins, &p3(arr));
        }
    }
    m.row("FieldVector::zero", Class::A, &[], &v3(FieldVector::zero()));
    m.row(
        "FieldVector::default",
        Class::A,
        &[],
        &v3(FieldVector::default()),
    );

    let qs = charges();
    let ps = pairs(&mut rng);
    for (k, &(a, b)) in ps.iter().enumerate() {
        let q = qs[k % qs.len()];
        let ins = cat(&[vec![f(q)], p3(a), p3(b)]);
        m.row(
            "field::electric_field_point_charge",
            Class::A,
            &ins,
            &res_v(electric_field_point_charge(black_box(q), a, b)),
        );
        m.row(
            "field::electric_potential_point_charge",
            Class::A,
            &ins,
            &res_f(electric_potential_point_charge(black_box(q), a, b)),
        );
        let vel = [
            rng.uniform(-1e7, 1e7),
            rng.uniform(-1e7, 1e7),
            rng.uniform(-1e7, 1e7),
        ];
        m.row(
            "field::magnetic_field_moving_charge",
            Class::A,
            &cat(&[vec![f(q)], p3(a), p3(vel), p3(b)]),
            &res_v(magnetic_field_moving_charge(black_box(q), a, vel, b)),
        );
        let dm = vs[(k * 11) % vs.len()];
        let din = cat(&[v3(dm), p3(b)]);
        m.row(
            "field::electric_dipole_field",
            Class::A,
            &din,
            &res_v(electric_dipole_field(&dm, black_box(b))),
        );
        m.row(
            "field::electric_dipole_potential",
            Class::A,
            &din,
            &res_f(electric_dipole_potential(&dm, black_box(b))),
        );
        m.row(
            "field::magnetic_dipole_field",
            Class::A,
            &din,
            &res_v(magnetic_dipole_field(&dm, black_box(b))),
        );
    }
    for n in 0..120 {
        let count = n % 6;
        let mut cs: Vec<(f64, [f64; 3])> = (0..count)
            .map(|_| {
                (
                    rng.signed_log_uniform(-9.0, -5.0),
                    [
                        rng.uniform(-3.0, 3.0),
                        rng.uniform(-3.0, 3.0),
                        rng.uniform(-3.0, 3.0),
                    ],
                )
            })
            .collect();
        let at = [
            rng.uniform(-3.0, 3.0),
            rng.uniform(-3.0, 3.0),
            rng.uniform(-3.0, 3.0),
        ];
        if n % 17 == 5 && !cs.is_empty() {
            cs[0].1 = at; // a charge at the field point: Singularity
        }
        let mut ins = vec![Tok::V(4 * cs.len())];
        for (q, p) in &cs {
            ins.push(f(*q));
            ins.extend(p3(*p));
        }
        ins.extend(p3(at));
        m.row(
            "field::electric_field_superposition",
            Class::A,
            &ins,
            &res_v(electric_field_superposition(black_box(&cs), at)),
        );
    }
    let mut es = specials();
    for k in (-1074..=1023).step_by(23) {
        es.push(pow2(k));
        es.push(-pow2(k) * 1.5);
    }
    for _ in 0..100 {
        es.push(rng.signed_log_uniform(-10.0, 10.0));
    }
    for &e in &es {
        m.row(
            "field::electric_energy_density",
            Class::A,
            &[f(e)],
            &[f(electric_energy_density(black_box(e)))],
        );
        m.row(
            "field::magnetic_energy_density",
            Class::A,
            &[f(e)],
            &[f(magnetic_energy_density(black_box(e)))],
        );
        m.row(
            "field::electric_field_infinite_plane",
            Class::A,
            &[f(e)],
            &[f(electric_field_infinite_plane(black_box(e)))],
        );
    }
    emit_3b(&mut m, &mut rng, &vs);
    emit_3c(&mut m);
    emit_5d(&mut m);
    emit_5d2(&mut m);
    emit_5d3(&mut m);
    emit_5d4(&mut m);
    m
}

/// The dipoles', the quadrupole's, the axes' and the sphere image's range edges (bijli ADR 0002 row 18, M1 part 5d-4),
/// appended: moments and tensors past 2^±600 against near, far and subnormal-offset field points; axes with z and R
/// past 2^±500; and sphere images whose |r|² over- or underflows (one where Rust's |r| of 0 says "inside").
fn emit_5d4(m: &mut Module) {
    let max = f64::MAX;
    let moments = [
        FieldVector::new(0.0, 0.0, 1e-30),
        FieldVector::new(1e300, 0.0, 1e300),
        FieldVector::new(1e-300, 2e-300, 0.0),
        FieldVector::new(max, 0.0, 0.0),
    ];
    let points = [
        [0.0, 0.0, 1.0],
        [1e200, 0.0, 0.0],
        [1e-14, 0.0, 2e-14],
        [3.0, 4.0, 1e-320],
        [1e150, 1e150, 1e150],
    ];
    for dm in moments {
        for b in points {
            let din = cat(&[v3(dm), p3(b)]);
            m.row(
                "field::electric_dipole_field",
                Class::A,
                &din,
                &res_v(electric_dipole_field(&dm, black_box(b))),
            );
            m.row(
                "field::electric_dipole_potential",
                Class::A,
                &din,
                &res_f(electric_dipole_potential(&dm, black_box(b))),
            );
            m.row(
                "field::magnetic_dipole_field",
                Class::A,
                &din,
                &res_v(magnetic_dipole_field(&dm, black_box(b))),
            );
        }
    }
    for scale in [1e300, 1e-300, 1e-20] {
        let qt = [[scale, 0.0, 0.0], [0.0, -0.5 * scale, 0.0], [0.0, 0.0, -0.5 * scale]];
        for b in [[1e100, 0.0, 0.0], [1e-14, 0.0, 1e-14], [1e70, 1e70, 0.0], [2.0, 3.0, 1e-320]] {
            let mut qin: Vec<Tok> = qt
                .iter()
                .flat_map(|row| row.iter().map(|&x| f(x)))
                .collect();
            qin.extend(p3(b));
            m.row(
                "field::electric_quadrupole_potential",
                Class::A,
                &qin,
                &res_f(electric_quadrupole_potential(&qt, black_box(b))),
            );
        }
    }
    for (q, rr, z) in [
        (1e300, 1e-200, 1e-200),
        (1e-300, 1e200, 1e200),
        (1.0, 1e154, 1e154),
        (1e-6, 1e-170, 3e-170),
        (max, 1.0, 1e-100),
    ] {
        m.row(
            "field::electric_field_ring_axis",
            Class::A,
            &[f(q), f(rr), f(z)],
            &res_f(electric_field_ring_axis(black_box(q), rr, z)),
        );
        m.row(
            "field::electric_field_disk_axis",
            Class::A,
            &[f(q), f(rr), f(z)],
            &res_f(electric_field_disk_axis(black_box(q), rr, z)),
        );
    }
    for (q, p, rr) in [
        (1.0, [1e200, 1e200, 0.0], 1.0),
        (1e300, [max, max, 0.0], 1.0),
        (1.0, [1e-200, 0.0, 0.0], 1e-210),
        (1.0, [3e-170, 4e-170, 0.0], 4.9e-170),
        (f64::from_bits(1), [1e300, 0.0, 0.0], 1e-300),
        (1.0, [1e150, 0.0, 0.0], 1e-200),
    ] {
        m.row(
            "field::image_charge_sphere",
            Class::A,
            &cat(&[vec![f(q)], p3(p), vec![f(rr)]]),
            &res_image(image_charge_sphere(black_box(q), p, rr)),
        );
    }
}

/// The point charges' and Green's functions' range edges (bijli ADR 0002 row 18, M1 part 5d-3), appended: results
/// that are representable where an intermediate is not. |E| = kq/r² while the factor kq/r³ under- or overflows (q =
/// 1.1e90 at r = 1e200, q = 1e290 at r = 1e-3); offsets past DBL_MAX; a unit-vector part on the subnormal grid for B;
/// superpositions whose fields overflow but cancel; and for the Green's functions 4πR past DBL_MAX, 1/(4πR) subnormal,
/// kR past DBL_MAX and kR subnormal.
fn emit_5d3(m: &mut Module) {
    let max = f64::MAX;
    let tiny = f64::from_bits(1);
    let qs = [1e-6, 1.1e90, 1e290, 1e300, 1e-300, tiny, 0.0, max];
    let pts: [([f64; 3], [f64; 3]); 9] = [
        ([0.0; 3], [1e200, 0.0, 0.0]),
        ([0.0; 3], [6e199, 8e199, 0.0]),
        ([0.0; 3], [1e-3, 0.0, 0.0]),
        ([0.0; 3], [1e300, 1e300, 0.0]),
        ([-max, 0.0, 0.0], [max, 0.0, 0.0]),
        ([1e160, 0.0, 0.0], [-1e160, 1e160, 0.0]),
        ([0.0; 3], [3.0, 1e-320, 0.0]),
        ([0.0; 3], [1e-14, 0.0, 2e-14]),
        ([1.0, 2.0, 3.0], [1.0 + 1e-15, 2.0, 3.0 + 1e-15]),
    ];
    let vel = [0.0, 2e6, 1e6];
    for &q in &qs {
        for &(a, b) in &pts {
            let ins = cat(&[vec![f(q)], p3(a), p3(b)]);
            m.row(
                "field::electric_field_point_charge",
                Class::A,
                &ins,
                &res_v(electric_field_point_charge(black_box(q), a, b)),
            );
            m.row(
                "field::electric_potential_point_charge",
                Class::A,
                &ins,
                &res_f(electric_potential_point_charge(black_box(q), a, b)),
            );
            m.row(
                "field::magnetic_field_moving_charge",
                Class::A,
                &cat(&[vec![f(q)], p3(a), p3(vel), p3(b)]),
                &res_v(magnetic_field_moving_charge(black_box(q), a, vel, b)),
            );
        }
    }
    let sets: Vec<Vec<(f64, [f64; 3])>> = vec![
        vec![(1e300, [1e-3, 0.0, 0.0]), (1e300, [-1e-3, 0.0, 0.0])],
        vec![(1.1e90, [1e200, 0.0, 0.0]), (-1.1e90, [0.0, 1e200, 0.0]), (2.2e90, [-1e200, 0.0, 0.0])],
        vec![(1e300, [1e-3, 0.0, 0.0]), (-1e300, [1e-3, 1e-18, 0.0]), (1e-6, [0.0, 0.0, 1.0])],
        vec![(1.1e90, [1e200, 0.0, 0.0])],
        vec![(1e-300, [0.0, 0.0, 1e200]), (1e-6, [1.0, 0.0, 0.0])],
    ];
    for cs in &sets {
        let at = [0.0, 0.0, 0.0];
        m.row(
            "field::electric_field_superposition",
            Class::A,
            &cat(&[charges_tok(cs), p3(at)]),
            &res_v(electric_field_superposition(black_box(cs), at)),
        );
    }
    let gp: [([f64; 3], [f64; 3]); 6] = [
        ([1e308, 0.0, 0.0], [0.0; 3]),
        ([1e307, 0.0, 0.0], [0.0; 3]),
        ([1e300, 0.0, 0.0], [0.0; 3]),
        ([1.0, 0.0, 0.0], [0.0; 3]),
        ([3e-29, 4e-29, 0.0], [0.0; 3]),
        ([max, 0.0, 0.0], [-max, 0.0, 0.0]),
    ];
    for &(b, a) in &gp {
        for kk in [2.0, 1e10, 1e-320, 1e300, 0.0] {
            m.row(
                "field::greens_function_scalar",
                Class::C,
                &cat(&[p3(b), p3(a), vec![f(kk)]]),
                &res_pair(greens_function_scalar(black_box(b), a, kk)),
            );
        }
        m.row(
            "field::greens_function_static",
            Class::A,
            &cat(&[p3(b), p3(a)]),
            &res_f(greens_function_static(black_box(b), a)),
        );
    }
}

/// Range-edge vectors (bijli ADR 0002 row 18, M1 part 5d-2): parts past 2^±511 where a square, product or sum leaves
/// the normal range though the result is representable, cancellations of products past DBL_MAX, and the subnormal
/// grid. Exact by construction (powers of two, literals and their products).
pub(crate) fn edge_vectors() -> Vec<FieldVector> {
    let e = f64::EPSILON;
    vec![
        FieldVector::new(pow2(600) * (1.0 + e), pow2(600), 0.0),
        FieldVector::new(pow2(500) * (1.0 - e / 2.0), -pow2(500), 0.0),
        FieldVector::new(pow2(-600), -pow2(-600), pow2(-600)),
        FieldVector::new(f64::MAX, f64::MAX, f64::MAX),
        FieldVector::new(1e300, 1e-300, 0.0),
        FieldVector::new(f64::from_bits(1), 0.0, 1e-310),
        FieldVector::new(pow2(1000), -pow2(999), pow2(998)),
        FieldVector::new(pow2(-1000), pow2(-1001), 0.0),
        FieldVector::new(3.0 * pow2(600), 4.0 * pow2(600), 12.0 * pow2(600)),
        FieldVector::new(pow2(700), pow2(700), 0.0),
        FieldVector::new(-2.5e-200, 1e-200, 7.5e-201),
        FieldVector::new(1.0, pow2(-600), -2.0),
        // Squares on the subnormal grid: Rust rounds each (2.25·2^-1074 to 2), the chain rounds their sum once.
        FieldVector::new(1.5 * pow2(-537), 1.5 * pow2(-537), 1.5 * pow2(-537)),
    ]
}

/// The row-18 edge grid, appended so the rows above keep their places: magnitude, magnitude_sq and normalized of
/// each edge vector, dot and cross of every pair, the plane image past 2^1023 and the dielectric image where its
/// product or sum leaves the range.
fn emit_5d2(m: &mut Module) {
    let vs = edge_vectors();
    for &v in &vs {
        let ins = v3(v);
        m.row("FieldVector::magnitude", Class::A, &ins, &[f(black_box(v).magnitude())]);
        m.row("FieldVector::magnitude_sq", Class::A, &ins, &[f(black_box(v).magnitude_sq())]);
        m.row("FieldVector::normalized", Class::A, &ins, &res_v(black_box(v).normalized()));
        for &w in &vs {
            let both = cat(&[v3(v), v3(w)]);
            m.row("FieldVector::dot", Class::A, &both, &[f(black_box(v).dot(&w))]);
            m.row("FieldVector::cross", Class::A, &both, &v3(black_box(v).cross(&w)));
        }
    }
    for pz in [f64::MAX, -f64::MAX, 1e308, 9e307, 8.9e307, -pow2(1023)] {
        for z in [f64::MAX, 1e308, -1e308, 1.0, 0.0, f64::from_bits(1)] {
            let a = [1.0, -2.0, z];
            let (qi, pi) = image_charge_plane(black_box(-3e-9), a, pz);
            m.row(
                "field::image_charge_plane",
                Class::A,
                &cat(&[vec![f(-3e-9)], p3(a), vec![f(pz)]]),
                &cat(&[vec![f(qi)], p3(pi)]),
            );
        }
    }
    let eps = [
        (1.0, 4.0),
        (1.0, 1.0 + f64::EPSILON),
        (f64::MAX, f64::MAX),
        (1e308, 1.5e308),
        (1e-300, 2e-300),
        (f64::from_bits(1), 1.0),
        (3.0, 3.0),
    ];
    for q in [1e300, f64::MAX, 1e-300, f64::from_bits(1), 1.0, 0.0] {
        for (e1, e2) in eps {
            let b = [0.5, 0.0, 1.0];
            m.row(
                "field::image_charge_dielectric",
                Class::A,
                &cat(&[vec![f(q)], p3(b), vec![f(e1), f(e2)]]),
                &res_image(image_charge_dielectric(black_box(q), b, e1, e2)),
            );
        }
    }
}

/// The singular cutoff's knife edge (bijli ADR 0002 row 17, M1 part 5d-1), appended so the rows above keep their
/// places: grids.rs's `knife_edge` points through every function that tests r² < 1e-30 (the dipoles and the
/// quadrupole on field points about the origin only) and, at 10^-30, the Green's functions' R < 1e-30.
fn emit_5d(m: &mut Module) {
    let q = 1e-6;
    let vel = [3e6, -2e6, 1e6];
    let dm = FieldVector::new(1e-30, -2e-30, 5e-31);
    let qt = [[1e-20, 0.0, 2e-21], [0.0, -5e-21, 0.0], [2e-21, 0.0, -5e-21]];
    for (a, b) in knife_edge(false) {
        let ins = cat(&[vec![f(q)], p3(a), p3(b)]);
        m.row(
            "field::electric_field_point_charge",
            Class::A,
            &ins,
            &res_v(electric_field_point_charge(black_box(q), a, b)),
        );
        m.row(
            "field::electric_potential_point_charge",
            Class::A,
            &ins,
            &res_f(electric_potential_point_charge(black_box(q), a, b)),
        );
        m.row(
            "field::magnetic_field_moving_charge",
            Class::A,
            &cat(&[vec![f(q)], p3(a), p3(vel), p3(b)]),
            &res_v(magnetic_field_moving_charge(black_box(q), a, vel, b)),
        );
        if a != [0.0; 3] {
            continue;
        }
        let din = cat(&[v3(dm), p3(b)]);
        m.row(
            "field::electric_dipole_field",
            Class::A,
            &din,
            &res_v(electric_dipole_field(&dm, black_box(b))),
        );
        m.row(
            "field::electric_dipole_potential",
            Class::A,
            &din,
            &res_f(electric_dipole_potential(&dm, black_box(b))),
        );
        m.row(
            "field::magnetic_dipole_field",
            Class::A,
            &din,
            &res_v(magnetic_dipole_field(&dm, black_box(b))),
        );
        let mut qin: Vec<Tok> = qt
            .iter()
            .flat_map(|row| row.iter().map(|&x| f(x)))
            .collect();
        qin.extend(p3(b));
        m.row(
            "field::electric_quadrupole_potential",
            Class::A,
            &qin,
            &res_f(electric_quadrupole_potential(&qt, black_box(b))),
        );
    }
    for (a, b) in knife_edge(true) {
        m.row(
            "field::greens_function_scalar",
            Class::C,
            &cat(&[p3(b), p3(a), vec![f(2.0)]]),
            &res_pair(greens_function_scalar(black_box(b), a, 2.0)),
        );
        m.row(
            "field::greens_function_static",
            Class::A,
            &cat(&[p3(b), p3(a)]),
            &res_f(greens_function_static(black_box(b), a)),
        );
    }
}

/// Part 3b: the Gauss's-law fields, the multipoles, the images, the Green's functions and the stress tensor.
fn emit_3b(m: &mut Module, rng: &mut Rng, vs: &[FieldVector]) {
    let ds = distances(rng);
    let qs = charges();
    for (k, &d) in ds.iter().enumerate() {
        let q = qs[k % qs.len()];
        let r = ds[(k * 13 + 5) % ds.len()];
        let ins = [f(q), f(d)];
        m.row(
            "field::electric_field_charged_sphere",
            Class::A,
            &ins,
            &res_f(electric_field_charged_sphere(black_box(q), d)),
        );
        m.row(
            "field::electric_field_infinite_line",
            Class::A,
            &ins,
            &res_f(electric_field_infinite_line(black_box(q), d)),
        );
        m.row(
            "field::electric_field_outside_cylinder",
            Class::A,
            &ins,
            &res_f(electric_field_outside_cylinder(black_box(q), d)),
        );
        m.row(
            "field::electric_field_inside_cylinder",
            Class::A,
            &ins,
            &res_f(electric_field_inside_cylinder(black_box(q), d)),
        );
        for (rr, dd) in [(r, d), (d, r), (1.0, d), (d, 0.5 * d), (d, d)] {
            m.row(
                "field::electric_field_inside_sphere",
                Class::A,
                &[f(q), f(rr), f(dd)],
                &res_f(electric_field_inside_sphere(black_box(q), rr, dd)),
            );
        }
        for z in [r, -r, d, 0.0, -0.0, 1e-11] {
            m.row(
                "field::electric_field_ring_axis",
                Class::A,
                &[f(q), f(d), f(z)],
                &res_f(electric_field_ring_axis(black_box(q), d, z)),
            );
            m.row(
                "field::electric_field_disk_axis",
                Class::A,
                &[f(q), f(d), f(z)],
                &res_f(electric_field_disk_axis(black_box(q), d, z)),
            );
        }
    }
    m.row(
        "field::electric_field_ring_axis",
        Class::A,
        &[f(1e-9), f(0.0), f(0.0)],
        &res_f(electric_field_ring_axis(black_box(1e-9), 0.0, 0.0)),
    );

    // Quadrupole tensors from the vector grid, a linear quadrupole among them; positions random and near the origin.
    for (k, &v) in vs.iter().enumerate().step_by(3) {
        let a = vs[(k * 5 + 1) % vs.len()];
        let b = vs[(k * 7 + 2) % vs.len()];
        let mut qt = [[v.x, v.y, v.z], [a.x, a.y, a.z], [b.x, b.y, b.z]];
        if k % 9 == 0 {
            qt = [[-2e-15, 0.0, 0.0], [0.0, -2e-15, 0.0], [0.0, 0.0, 4e-15]];
        }
        let p = if k % 11 == 4 {
            [0.0, 1e-16, 0.0]
        } else {
            [
                rng.uniform(-3.0, 3.0),
                rng.uniform(-3.0, 3.0),
                rng.uniform(-3.0, 3.0),
            ]
        };
        let mut ins: Vec<Tok> = qt
            .iter()
            .flat_map(|row| row.iter().map(|&x| f(x)))
            .collect();
        ins.extend(p3(p));
        m.row(
            "field::electric_quadrupole_potential",
            Class::A,
            &ins,
            &res_f(electric_quadrupole_potential(&qt, black_box(p))),
        );
    }

    // Multipole moments: l ≤ 10 for every m (and |m| = l + 1, refused), then high orders.
    let sets = charge_sets(rng);
    let mut lms: Vec<(u32, i32)> = Vec::new();
    for l in 0..=10_u32 {
        let li = i32::try_from(l).unwrap_or(0);
        for mm in -li - 1..=li + 1 {
            lms.push((l, mm));
        }
    }
    lms.extend([
        (20, 7),
        (20, -20),
        (50, 25),
        (200, 200),
        (200, 0),
        (200, -199),
    ]);
    for (n, cs) in sets.iter().enumerate() {
        for (j, &(l, mm)) in lms.iter().enumerate() {
            if (j + n) % 3 != 0 && l <= 10 {
                continue;
            }
            let ins = cat(&[vec![crate::tsv::i(l), crate::tsv::i(mm)], charges_tok(cs)]);
            m.row(
                "field::multipole_moment",
                Class::C,
                &ins,
                &res_f(multipole_moment(black_box(l), mm, cs)),
            );
            m.row(
                "field::spherical_multipole_moment",
                Class::C,
                &ins,
                &res_pair(spherical_multipole_moment(black_box(l), mm, cs)),
            );
        }
    }

    // Images and Green's functions over the position pairs.
    let ps = pairs(rng);
    for (k, &(a, b)) in ps.iter().enumerate() {
        let q = qs[k % qs.len()];
        let rr = ds[k % ds.len()];
        let e1 = [1.0, 4.0, 2.0, 2.0, 0.0, -1.0, f64::NAN, 80.0][k % 8];
        let e2 = [4.0, 1.0, 2.0, -1.0, 4.0, 1.0, 2.0, f64::INFINITY][(k / 8) % 8];
        let pz = b[2];
        let (qi, pi) = image_charge_plane(black_box(q), a, pz);
        m.row(
            "field::image_charge_plane",
            Class::A,
            &cat(&[vec![f(q)], p3(a), vec![f(pz)]]),
            &cat(&[vec![f(qi)], p3(pi)]),
        );
        m.row(
            "field::image_charge_sphere",
            Class::A,
            &cat(&[vec![f(q)], p3(b), vec![f(rr)]]),
            &res_image(image_charge_sphere(black_box(q), b, rr)),
        );
        m.row(
            "field::image_charge_dielectric",
            Class::A,
            &cat(&[vec![f(q)], p3(b), vec![f(e1), f(e2)]]),
            &res_image(image_charge_dielectric(black_box(q), b, e1, e2)),
        );
        let kk = [
            0.0,
            1.0,
            100.0,
            2.0 * std::f64::consts::PI,
            -3.5,
            1e6,
            f64::NAN,
        ][k % 7];
        m.row(
            "field::greens_function_scalar",
            Class::C,
            &cat(&[p3(b), p3(a), vec![f(kk)]]),
            &res_pair(greens_function_scalar(black_box(b), a, kk)),
        );
        m.row(
            "field::greens_function_static",
            Class::A,
            &cat(&[p3(b), p3(a)]),
            &res_f(greens_function_static(black_box(b), a)),
        );
    }

    // The stress tensor, its force, and the Poynting flux.
    for (k, &e) in vs.iter().enumerate().step_by(2) {
        let b = vs[(k * 3 + 7) % vs.len()].scale(1e-4);
        let nrm = vs[(k * 5 + 1) % vs.len()];
        let t = maxwell_stress_tensor(&black_box(e), &b);
        let tt: Vec<Tok> = t.iter().flat_map(|row| row.iter().map(|&x| f(x))).collect();
        m.row(
            "field::maxwell_stress_tensor",
            Class::A,
            &cat(&[v3(e), v3(b)]),
            &tt,
        );
        m.row(
            "field::stress_tensor_force",
            Class::A,
            &cat(&[tt.clone(), v3(nrm)]),
            &v3(stress_tensor_force(&black_box(t), &nrm)),
        );
    }
    for n in 0..60 {
        let count = n % 5;
        let samples: Vec<(FieldVector, FieldVector, FieldVector, f64)> = (0..count)
            .map(|j| {
                (
                    vs[(n * 7 + j * 3) % vs.len()],
                    vs[(n * 11 + j * 5 + 1) % vs.len()].scale(1e-6),
                    vs[(n * 13 + j) % vs.len()],
                    rng.log_uniform(-4.0, 0.0),
                )
            })
            .collect();
        let mut ins = vec![Tok::V(10 * samples.len())];
        for (e, b, nrm, da) in &samples {
            ins.extend(v3(*e));
            ins.extend(v3(*b));
            ins.extend(v3(*nrm));
            ins.push(f(*da));
        }
        m.row(
            "field::poynting_flux_surface",
            Class::A,
            &ins,
            &[f(poynting_flux_surface(black_box(&samples)))],
        );
    }
}

/// trace_field_line's field functions, by number, from their parameters; the harness builds the same ones.
/// 0 uniform (fx, fy, fz); 1 a point charge q at (x, y, z); 2 an electric dipole p at the origin; 3 a linear
/// quadrupole (+q at ±d ẑ, −2q at the origin); 4 a sink: Singularity once x ≥ a, else x̂; 5 always ZeroCharge;
/// 6 a NaN field once x ≥ a, else x̂; 7 an infinite field once x ≥ a, else x̂; 8 a magnetic dipole m at the origin.
/// A boxed field function, as `trace_field_line` takes one.
type FieldFn = Box<dyn Fn([f64; 3]) -> bijli::error::Result<FieldVector>>;

/// One trace row's inputs: the scenario, the start, the step, max_steps and the scenario's parameters.
type TraceCase = (i64, [f64; 3], f64, usize, Vec<f64>);

fn field_fn(id: i64, p: Vec<f64>) -> FieldFn {
    let ex = FieldVector::new(1.0, 0.0, 0.0);
    match id {
        0 => Box::new(move |_| Ok(FieldVector::new(p[0], p[1], p[2]))),
        1 => Box::new(move |pos| electric_field_point_charge(p[0], [p[1], p[2], p[3]], pos)),
        2 => Box::new(move |pos| electric_dipole_field(&FieldVector::new(p[0], p[1], p[2]), pos)),
        3 => Box::new(move |pos| {
            let cs = [
                (p[0], [0.0, 0.0, p[1]]),
                (p[0], [0.0, 0.0, -p[1]]),
                (-2.0 * p[0], [0.0; 3]),
            ];
            electric_field_superposition(&cs, pos)
        }),
        4 => Box::new(move |pos| {
            if pos[0] >= p[0] {
                Err(BijliError::Singularity)
            } else {
                Ok(ex)
            }
        }),
        5 => Box::new(|_| Err(BijliError::ZeroCharge)),
        6 => Box::new(move |pos| {
            if pos[0] >= p[0] {
                Ok(FieldVector::new(f64::NAN, 0.0, 0.0))
            } else {
                Ok(ex)
            }
        }),
        7 => Box::new(move |pos| {
            if pos[0] >= p[0] {
                Ok(FieldVector::new(f64::INFINITY, 0.0, 0.0))
            } else {
                Ok(ex)
            }
        }),
        _ => Box::new(move |pos| magnetic_dipole_field(&FieldVector::new(p[0], p[1], p[2]), pos)),
    }
}

/// Part 3c: `trace_field_line` over the scenarios, their singular and escaping ends, and the rejected inputs. A row's
/// inputs are the scenario number, the start, the step, max_steps as its 64-bit word (usize::MAX is `i:-1`) and the
/// scenario's parameters (`v:<k>`); its outputs the points as `v:<3n>` and their coordinates, or the error.
fn emit_3c(m: &mut Module) {
    let cases: Vec<TraceCase> = vec![
        (0, [0.0; 3], 0.1, 10, vec![1.0, 0.0, 0.0]),
        (0, [0.5, -1.0, 2.0], 0.5, 7, vec![0.0, 0.0, -2.0]),
        (0, [0.0; 3], 0.1, 50, vec![3.0, 4.0, 12.0]),
        (0, [0.0; 3], 0.5, 2, vec![0.0, 1e-35, 0.0]),
        (0, [0.0; 3], 0.1, 10, vec![0.0, 0.0, 0.0]),
        (0, [0.0; 3], 0.1, 10, vec![-0.0, 0.0, -0.0]),
        (0, [1.0, 1.0, 1.0], 0.25, 9, vec![1e300, 1e300, 0.0]),
        (0, [0.0; 3], 0.1, 0, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], 0.1, 1, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], 0.1, usize::MAX, vec![0.0, 0.0, 0.0]),
        (1, [0.1, 0.0, 0.0], 0.01, 5, vec![1e-6, 0.0, 0.0, 0.0]),
        (1, [0.1, 0.0, 0.0], 0.01, 10, vec![1e-6, 0.0, 0.0, 0.0]),
        (1, [0.3, 0.2, -0.1], 0.05, 60, vec![1e-6, 0.0, 0.0, 0.0]),
        (1, [0.3, 0.2, -0.1], 0.05, 60, vec![-1e-6, 0.0, 0.0, 0.0]),
        (1, [0.1, 0.0, 0.0], 0.03, 40, vec![-1e-6, 0.0, 0.0, 0.0]),
        (1, [2.0, 2.0, 2.0], 0.2, 30, vec![-3e-9, 1.0, 1.0, 1.0]),
        (2, [0.5, 0.0, 0.5], 0.02, 100, vec![0.0, 0.0, 1e-30]),
        (2, [0.0, 0.3, -0.4], 0.05, 80, vec![1e-30, -2e-30, 0.5e-30]),
        (8, [0.5, 0.0, 0.5], 0.02, 100, vec![0.0, 0.0, 1.0]),
        (3, [0.01, 0.0, 0.01], 0.001, 100, vec![1e-9, 1e-3]),
        (3, [0.0, 0.02, 0.0], 0.002, 60, vec![-1e-9, 1e-3]),
        (4, [0.0; 3], 0.1, 100, vec![0.25]),
        (4, [0.0; 3], 0.1, 100, vec![-1.0]),
        (5, [0.0; 3], 0.1, 10, vec![]),
        (5, [0.0; 3], 0.1, 0, vec![]),
        (6, [0.0; 3], 0.1, 100, vec![0.25]),
        (6, [1e-5, -0.0, 1e16], 0.1, 3, vec![-1.0]),
        (7, [0.0; 3], 0.1, 100, vec![0.25]),
        (0, [0.0; 3], -1.0, 10, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], 0.0, 10, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], -0.0, 10, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], f64::NAN, 10, vec![1.0, 0.0, 0.0]),
        (0, [0.0; 3], f64::INFINITY, 10, vec![1.0, 0.0, 0.0]),
        (0, [f64::NAN, 0.0, 0.0], 0.1, 10, vec![1.0, 0.0, 0.0]),
        (
            0,
            [f64::INFINITY, f64::NEG_INFINITY, 1e-5],
            0.1,
            10,
            vec![1.0, 0.0, 0.0],
        ),
        (0, [0.25, 1e20, -3.5e-7], 0.1, 3, vec![f64::NAN, 0.0, 0.0]),
    ];
    for (id, start, step, max, params) in cases {
        let word = i64::from_ne_bytes(max.to_ne_bytes());
        let mut ins = vec![
            crate::tsv::i(id),
            f(start[0]),
            f(start[1]),
            f(start[2]),
            f(step),
            crate::tsv::i(word),
        ];
        ins.push(Tok::V(params.len()));
        ins.extend(params.iter().map(|&x| f(x)));
        let outs = match trace_field_line(black_box(start), step, max, field_fn(id, params)) {
            Ok(points) => {
                let mut v = vec![Tok::V(3 * points.len())];
                for pt in &points {
                    v.extend(p3(*pt));
                }
                v
            }
            Err(e) => vec![err(&e), s(e.to_string())],
        };
        m.row("field::trace_field_line", Class::A, &ins, &outs);
    }
}
