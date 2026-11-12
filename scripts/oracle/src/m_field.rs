//! `field.tsv` — rust-old/src/field.rs, M1 part 3a: `FieldVector`, the point-charge fields and potential, the
//! superposition, the energy densities and the dipoles (3b and 3c add the rest of the module).
//!
//! Every row is class A: the functions are + − × ÷ and sqrt in a fixed order, which the port repeats. A vector is
//! three `f:` tokens; a position `[f64; 3]` likewise. A `Result` is its value's tokens on `Ok`, or `E:<Variant>` and
//! the error's Display text on `Err` (ADR 0002 row 4: the port reads both back from the last-error record). The
//! superposition's charges are `v:<4n>` then q, x, y, z per charge.
//!
//! The grids: the specials crossed into vectors, components over every 37th binade with both signs (the rescaled
//! path of `normalized` among them), seeded vectors and SI-scale charges; positions on a lattice, at random, and at
//! offsets from the source around the singular cutoff r² = 1e-30 (√1e-30 = 1e-15 m) on both sides.

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::field::{
    FieldVector, electric_dipole_field, electric_dipole_potential, electric_energy_density,
    electric_field_point_charge, electric_field_superposition, electric_potential_point_charge,
    magnetic_dipole_field, magnetic_energy_density, magnetic_field_moving_charge,
};
use std::hint::black_box;

fn v3(v: FieldVector) -> Vec<Tok> {
    vec![f(v.x), f(v.y), f(v.z)]
}

fn p3(p: [f64; 3]) -> Vec<Tok> {
    vec![f(p[0]), f(p[1]), f(p[2])]
}

fn res_v(r: bijli::error::Result<FieldVector>) -> Vec<Tok> {
    match r {
        Ok(v) => v3(v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn res_f(r: bijli::error::Result<f64>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn cat(parts: &[Vec<Tok>]) -> Vec<Tok> {
    parts.iter().flatten().cloned().collect()
}

/// Vectors: the specials crossed (x, y, 0) and (0, y, z), every 37th binade, and seeded.
fn vectors() -> Vec<FieldVector> {
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
fn pairs(rng: &mut Rng) -> Vec<([f64; 3], [f64; 3])> {
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
fn charges() -> Vec<f64> {
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
    }
    m
}
