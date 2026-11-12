//! `charge.tsv` — rust-old/src/charge.rs (M1 part 3d): `PointCharge`, Coulomb's law, the Lorentz force, the dipole
//! moment and energies, cyclotron frequency and Larmor radius.
//!
//! Every row is class A: + − × ÷, sqrt and abs in Rust's order, which the port repeats (the `PointCharge` methods
//! call field.rs's functions, whose rows are class A too). A `PointCharge` is eight `f:` tokens: charge, mass, the
//! position's three, the velocity's three. A vector or a position is three `f:` tokens. A `Result` is its value's
//! tokens on `Ok`, or `E:<Variant>` and the error's Display text on `Err` (ADR 0002 row 4).
//!
//! The grids: masses (the specials, both zeros, negatives, the electron and proton masses) crossed with field.rs's
//! charges; field.rs's (source, field point) pairs, the cutoff r² = 1e-30 among them, for the methods and Coulomb's
//! force; field.rs's vectors for the Lorentz force; field.rs's distances for the Coulomb energy; and the cyclotron
//! and Larmor arguments crossed over the specials and SI scales (the weak-field Larmor case among them).

use crate::grids::{Rng, specials};
use crate::m_field::{cat, charges, distances, p3, pairs, res_f, res_v, v3, vectors};
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::charge::{
    ELECTRON_MASS, ELEMENTARY_CHARGE, PROTON_MASS, PointCharge, coulomb_force,
    coulomb_potential_energy, cyclotron_frequency, dipole_moment, larmor_radius, lorentz_force,
    potential_energy,
};
use std::hint::black_box;

fn pc(p: &PointCharge) -> Vec<Tok> {
    cat(&[vec![f(p.charge), f(p.mass)], p3(p.position), p3(p.velocity)])
}

fn res_pc(r: bijli::error::Result<PointCharge>) -> Vec<Tok> {
    match r {
        Ok(p) => pc(&p),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// Masses: the specials, negatives, the particle masses, SI scales.
fn masses() -> Vec<f64> {
    let mut ms = specials();
    ms.extend([-1.0, -1e-300, ELECTRON_MASS, PROTON_MASS, 1e-35, 2.5, 1e30]);
    ms
}

pub fn emit() -> Module {
    let mut m = Module::new("charge");
    let mut rng = Rng::new(0xc4a_29e);
    let qs = charges();
    let ms = masses();
    let ps = pairs(&mut rng);
    let vs = vectors();
    for (k, &mass) in ms.iter().enumerate() {
        for (j, &q) in qs.iter().enumerate() {
            let pos = ps[(k * 13 + j) % ps.len()].1;
            let vel = [
                rng.uniform(-1e7, 1e7),
                rng.uniform(-1e7, 1e7),
                rng.uniform(-1e7, 1e7),
            ];
            m.row(
                "PointCharge::new",
                Class::A,
                &cat(&[vec![f(q), f(mass)], p3(pos), p3(vel)]),
                &res_pc(PointCharge::new(black_box(q), mass, pos, vel)),
            );
        }
    }
    for (k, &(a, b)) in ps.iter().enumerate() {
        let pos = if k % 2 == 0 { a } else { b };
        m.row(
            "PointCharge::electron",
            Class::A,
            &p3(pos),
            &pc(&PointCharge::electron(black_box(pos))),
        );
        m.row(
            "PointCharge::proton",
            Class::A,
            &p3(pos),
            &pc(&PointCharge::proton(black_box(pos))),
        );
    }
    for (k, &(a, b)) in ps.iter().enumerate() {
        let vel = [
            rng.uniform(-1e7, 1e7),
            rng.uniform(-1e7, 1e7),
            rng.uniform(-1e7, 1e7),
        ];
        let p = PointCharge {
            charge: qs[k % qs.len()],
            mass: [1.0, ELECTRON_MASS, PROTON_MASS][k % 3],
            position: a,
            velocity: vel,
        };
        let ins = cat(&[pc(&p), p3(b)]);
        m.row(
            "PointCharge::electric_field_at",
            Class::A,
            &ins,
            &res_v(black_box(&p).electric_field_at(b)),
        );
        m.row(
            "PointCharge::electric_potential_at",
            Class::A,
            &ins,
            &res_f(black_box(&p).electric_potential_at(b)),
        );
        m.row(
            "PointCharge::magnetic_field_at",
            Class::A,
            &ins,
            &res_v(black_box(&p).magnetic_field_at(b)),
        );
        let e = vs[(k * 7 + 1) % vs.len()];
        let bf = vs[(k * 11 + 5) % vs.len()];
        m.row(
            "PointCharge::lorentz_force_in",
            Class::A,
            &cat(&[pc(&p), v3(e), v3(bf)]),
            &v3(black_box(&p).lorentz_force_in(&e, &bf)),
        );
        let q2 = PointCharge {
            charge: qs[(k * 3 + 1) % qs.len()],
            mass: 1.0,
            position: b,
            velocity: [0.0; 3],
        };
        m.row(
            "charge::coulomb_force",
            Class::A,
            &cat(&[pc(&p), pc(&q2)]),
            &res_v(coulomb_force(black_box(&p), &q2)),
        );
    }
    for (k, &v) in vs.iter().enumerate() {
        let q = qs[k % qs.len()];
        let e = vs[(k * 7 + 3) % vs.len()];
        let bf = vs[(k * 13 + 2) % vs.len()];
        m.row(
            "charge::lorentz_force",
            Class::A,
            &cat(&[vec![f(q)], v3(v), v3(e), v3(bf)]),
            &v3(lorentz_force(black_box(q), &v, &e, &bf)),
        );
    }
    let ds = distances(&mut rng);
    for &q in &qs {
        for &x in &ds {
            m.row(
                "charge::dipole_moment",
                Class::A,
                &[f(q), f(x)],
                &[f(dipole_moment(black_box(q), x))],
            );
            m.row(
                "charge::potential_energy",
                Class::A,
                &[f(q), f(x)],
                &[f(potential_energy(black_box(q), x))],
            );
        }
    }
    for (k, &d) in ds.iter().enumerate() {
        for j in 0..3 {
            let q1 = qs[(k + j) % qs.len()];
            let q2 = qs[(k * 3 + j + 1) % qs.len()];
            m.row(
                "charge::coulomb_potential_energy",
                Class::A,
                &[f(q1), f(q2), f(d)],
                &res_f(coulomb_potential_energy(black_box(q1), q2, d)),
            );
        }
    }
    let cq = [
        0.0,
        -0.0,
        ELEMENTARY_CHARGE,
        -ELEMENTARY_CHARGE,
        f64::NAN,
        f64::INFINITY,
        1e-300,
        1.0,
    ];
    let cm = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MIN_POSITIVE,
        ELECTRON_MASS,
        PROTON_MASS,
        1e300,
    ];
    let cb = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        1e-12,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        1e300,
    ];
    let cv = [0.0, -1e6, 1e3, f64::NAN, f64::INFINITY, 1e-300];
    for &q in &cq {
        for &mass in &cm {
            for &bm in &cb {
                m.row(
                    "charge::cyclotron_frequency",
                    Class::A,
                    &[f(q), f(mass), f(bm)],
                    &res_f(cyclotron_frequency(black_box(q), mass, bm)),
                );
                for &v in &cv {
                    m.row(
                        "charge::larmor_radius",
                        Class::A,
                        &[f(mass), f(v), f(q), f(bm)],
                        &res_f(larmor_radius(black_box(mass), v, q, bm)),
                    );
                }
            }
        }
    }
    for _ in 0..200 {
        let q = rng.signed_log_uniform(-19.0, -6.0);
        let mass = rng.log_uniform(-31.0, 0.0);
        let bm = rng.log_uniform(-12.0, 2.0);
        let v = rng.signed_log_uniform(0.0, 8.0);
        m.row(
            "charge::cyclotron_frequency",
            Class::A,
            &[f(q), f(mass), f(bm)],
            &res_f(cyclotron_frequency(black_box(q), mass, bm)),
        );
        m.row(
            "charge::larmor_radius",
            Class::A,
            &[f(mass), f(v), f(q), f(bm)],
            &res_f(larmor_radius(black_box(mass), v, q, bm)),
        );
    }
    m
}
