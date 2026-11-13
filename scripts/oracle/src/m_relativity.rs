//! `relativity.tsv` — rust-old/src/relativity.rs (M2 part 1d-1): the Lorentz factor, the EM tensor, the field
//! transforms, four-vectors, the retarded and Liénard–Wiechert potentials and fields, and the Larmor powers. The
//! retarded-time solver and `lienard_wiechert_fields` (callbacks) are part 1d-2's.
//!
//! Every row is class A: + − × ÷ and sqrt in Rust's order. A `Result` is its value, or `E:<Variant>` and the error's
//! Display text (ADR 0002 row 4). A vector is three `f:` tokens, a four-vector four (ct, x, y, z), an `EmTensor` its
//! sixteen components row-major, a pair of vectors six.
//!
//! The grids: speeds through c (0.1c, 0.9c, 0.99c, (1 − 1e-9)c and c ulp by ulp, both signs, the specials, binades,
//! seeded); field pairs over a stride of field.rs's vectors and its edge vectors; tensors from those fields and
//! arbitrary ones; boost velocities at β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9} along x, y and the diagonal, at and past c;
//! four-vectors over the specials and physical events; Liénard–Wiechert geometries with the 1e-30 cuts of R and of
//! 1 − n̂⋅β, and β → 1.

use crate::grids::{Rng, pow2, specials};
use crate::m_field::{cat, edge_vectors, res_f, res_v, v3, vectors};
use crate::m_material::{dedup, ulps};
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::field::{FieldVector, SPEED_OF_LIGHT};
use bijli::relativity::{
    EmTensor, FourVector, beta, larmor_power, lienard_wiechert_b, lienard_wiechert_e, lienard_wiechert_scalar_potential,
    lienard_wiechert_vector_potential, lorentz_factor, lorentz_transform_fields, lorentz_transform_x,
    relativistic_larmor_power, retarded_scalar_potential, retarded_vector_potential,
};
use std::hint::black_box;

const C: f64 = SPEED_OF_LIGHT;

/// A tensor's tokens: its sixteen components, row-major.
fn tt(t: &EmTensor) -> Vec<Tok> {
    t.components.iter().flatten().map(|&x| f(x)).collect()
}

/// A four-vector's tokens: (ct, x, y, z).
fn fv4(v: FourVector) -> Vec<Tok> {
    vec![f(v.t), f(v.x), f(v.y), f(v.z)]
}

fn res_4(r: bijli::error::Result<FourVector>) -> Vec<Tok> {
    match r {
        Ok(v) => fv4(v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn res_pair(r: bijli::error::Result<(FieldVector, FieldVector)>) -> Vec<Tok> {
    match r {
        Ok((a, b)) => cat(&[v3(a), v3(b)]),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// Speeds: the specials, c and its neighbours, fractions of c, both signs, binades, seeded.
fn speeds(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    for b in [0.1, 0.5, 0.9, 0.99, 1.0 - 1e-9] {
        xs.push(b * C);
        xs.push(-b * C);
    }
    for k in -3..=3 {
        xs.push(ulps(C, k));
        xs.push(-ulps(C, k));
    }
    xs.extend([2.0 * C, 1e-300, 1e300, 1.0, 1000.0]);
    for k in (-1074..=1023).step_by(150) {
        xs.push(pow2(k));
    }
    for _ in 0..8 {
        xs.push(rng.uniform(-C, C));
    }
    dedup(xs)
}

/// Boost velocities: β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9} along x, y and the diagonal, at and past c, NaN, tiny, huge.
fn boosts() -> Vec<FieldVector> {
    let mut vs = Vec::new();
    for b in [0.0, 0.1, 0.9, 0.99, 1.0 - 1e-9, 1.0, 1.5] {
        vs.push(FieldVector::new(b * C, 0.0, 0.0));
        vs.push(FieldVector::new(0.0, -b * C, 0.0));
        let d = b * C / 3.0_f64.sqrt();
        vs.push(FieldVector::new(d, d, d));
    }
    vs.push(FieldVector::new(f64::NAN, 0.0, 0.0));
    vs.push(FieldVector::new(1e-300, 0.0, 0.0));
    vs.push(FieldVector::new(1e-200, 1e-200, 0.0));
    vs.push(FieldVector::new(1e300, 0.0, 0.0));
    vs.push(FieldVector::new(f64::INFINITY, 0.0, 0.0));
    vs
}

/// Four-vectors: specials in each slot, events, null and timelike vectors.
fn four_vectors(rng: &mut Rng) -> Vec<FourVector> {
    let mut out = Vec::new();
    for &x in &specials() {
        out.push(FourVector::new(x, 1.0, 2.0, 3.0));
        out.push(FourVector::new(1.0, x, 0.5, -0.25));
    }
    out.extend([
        FourVector::new(0.0, 0.0, 0.0, 0.0),
        FourVector::new(C, C, 0.0, 0.0),
        FourVector::new(5.0, 3.0, 4.0, 0.0),
        FourVector::new(1e200, 1e200, 0.0, 0.0),
        FourVector::new(1e-200, 0.0, 1e-200, 0.0),
        FourVector::new(1e300, -1e300, 1e300, 0.0),
    ]);
    for _ in 0..12 {
        out.push(FourVector::new(rng.uniform(-1e9, 1e9), rng.uniform(-10.0, 10.0), rng.uniform(-10.0, 10.0), rng.uniform(-10.0, 10.0)));
    }
    out
}

pub fn emit() -> Module {
    let mut m = Module::new("relativity");
    let mut rng = Rng::new(0x5e1_a71);
    let sp = speeds(&mut rng);

    for &v in &sp {
        m.row("relativity::lorentz_factor", Class::A, &[f(v)], &res_f(lorentz_factor(black_box(v))));
        m.row("relativity::beta", Class::A, &[f(v)], &[f(beta(black_box(v)))]);
    }

    // Fields: a stride of field.rs's vectors and its edge vectors, crossed in pairs.
    let mut vs: Vec<FieldVector> = vectors().into_iter().step_by(61).collect();
    vs.extend(edge_vectors());
    vs.extend([FieldVector::new(1000.0, 0.0, 0.0), FieldVector::new(0.0, 0.0, 1e-6), FieldVector::new(0.0, 1e5, 0.0)]);
    let mut tensors = Vec::new();
    for &e in &vs {
        for &b in &vs {
            let ins = cat(&[v3(e), v3(b)]);
            let t = EmTensor::from_fields(black_box(&e), &b);
            m.row("EmTensor::from_fields", Class::A, &ins, &tt(&t));
            tensors.push(t);
        }
    }
    // Arbitrary tensors: specials in one slot, seeded components.
    for k in 0..16 {
        for &x in &[f64::NAN, f64::INFINITY, 1e300, -0.0] {
            let mut c = [[0.5, -1.0, 2.0, 0.25]; 4];
            c[k / 4][k % 4] = x;
            tensors.push(EmTensor { components: c });
        }
    }
    for _ in 0..40 {
        let mut c = [[0.0; 4]; 4];
        for row in &mut c {
            for x in row.iter_mut() {
                *x = rng.signed_log_uniform(-6.0, 6.0);
            }
        }
        tensors.push(EmTensor { components: c });
    }
    let us = [
        FourVector::new(C, 0.0, 0.0, 0.0),
        FourVector::new(C / 0.6, 0.8 * C / 0.6, 0.0, 0.0),
        FourVector::new(1.0, -2.0, 0.5, 3.0),
        FourVector::new(f64::NAN, 0.0, 0.0, 0.0),
        FourVector::new(1e300, 1e300, 0.0, 0.0),
    ];
    for (j, t) in tensors.iter().enumerate() {
        let ins = tt(t);
        m.row("EmTensor::electric_field", Class::A, &ins, &v3(black_box(t).electric_field()));
        m.row("EmTensor::magnetic_field", Class::A, &ins, &v3(black_box(t).magnetic_field()));
        m.row("EmTensor::first_invariant", Class::A, &ins, &[f(black_box(t).first_invariant())]);
        m.row("EmTensor::second_invariant", Class::A, &ins, &[f(black_box(t).second_invariant())]);
        m.row("EmTensor::dual", Class::A, &ins, &tt(&black_box(t).dual()));
        if j % 3 == 0 {
            for (q, u) in [(1.0, us[j % 5]), (-1.6e-19, us[(j + 1) % 5])] {
                let ins4 = cat(&[ins.clone(), vec![f(q)], fv4(u)]);
                m.row("EmTensor::four_force", Class::A, &ins4, &fv4(black_box(t).four_force(q, &u)));
            }
        }
    }

    // The field transforms: fields from a short list, boosts.
    let fs: Vec<FieldVector> = vs.iter().step_by(3).copied().collect();
    let bs = boosts();
    for &e in &fs {
        for &b in fs.iter().step_by(2) {
            for &v in &bs {
                let ins = cat(&[v3(e), v3(b), v3(v)]);
                m.row(
                    "relativity::lorentz_transform_fields",
                    Class::A,
                    &ins,
                    &res_pair(lorentz_transform_fields(black_box(&e), &b, &v)),
                );
            }
            for &u in sp.iter().step_by(3) {
                let ins = cat(&[v3(e), v3(b), vec![f(u)]]);
                m.row("relativity::lorentz_transform_x", Class::A, &ins, &res_pair(lorentz_transform_x(black_box(&e), &b, u)));
            }
        }
    }

    // Four-vectors.
    let fvs = four_vectors(&mut rng);
    for &a in &fvs {
        let ins = fv4(a);
        m.row("FourVector::new", Class::A, &ins, &fv4(FourVector::new(black_box(a.t), a.x, a.y, a.z)));
        m.row("FourVector::interval_sq", Class::A, &ins, &[f(black_box(a).interval_sq())]);
        m.row("FourVector::spatial_magnitude", Class::A, &ins, &[f(black_box(a).spatial_magnitude())]);
        for &b in fvs.iter().step_by(4) {
            m.row("FourVector::dot", Class::A, &cat(&[ins.clone(), fv4(b)]), &[f(black_box(a).dot(&b))]);
        }
        for &u in sp.iter().step_by(4) {
            m.row("FourVector::boost_x", Class::A, &cat(&[ins.clone(), vec![f(u)]]), &res_4(black_box(a).boost_x(u)));
        }
    }
    for &t in &sp {
        for p in [[1.0, 2.0, 3.0], [f64::NAN, 0.0, -0.0]] {
            let ins = [f(t), f(p[0]), f(p[1]), f(p[2])];
            m.row("FourVector::from_time_position", Class::A, &ins, &fv4(FourVector::from_time_position(black_box(t), p)));
        }
    }

    // The retarded potentials over (q, d) and a velocity.
    let mut ds = specials();
    ds.extend([2.0, 1e-30, ulps(1e-30, -1), ulps(1e-30, 1), 1e-300, 1e300, -2.0]);
    let qs = [1.0, -1.602_176_634e-19, 0.0, f64::NAN, 1e300];
    for &q in &qs {
        for &d in &ds {
            m.row("relativity::retarded_scalar_potential", Class::A, &[f(q), f(d)], &res_f(retarded_scalar_potential(black_box(q), d)));
            for &v in bs.iter().step_by(5) {
                let ins = cat(&[vec![f(q)], v3(v), vec![f(d)]]);
                m.row(
                    "relativity::retarded_vector_potential",
                    Class::A,
                    &ins,
                    &res_v(retarded_vector_potential(black_box(q), &v, d)),
                );
            }
        }
    }

    // Liénard–Wiechert: R over the field vectors, β over fractions of c in a few directions, β̇.
    let mut betas: Vec<FieldVector> = bs.iter().map(|v| v.scale(1.0 / C)).collect();
    betas.push(FieldVector::new(1.0 - f64::EPSILON / 2.0, 0.0, 0.0));
    betas.push(FieldVector::new(ulps(1.0, -3), 0.0, 0.0));
    let mut rs: Vec<FieldVector> = vs.iter().step_by(2).copied().collect();
    rs.extend([
        FieldVector::new(1.0, 0.0, 0.0),
        FieldVector::new(1e-30, 0.0, 0.0),
        FieldVector::new(ulps(1e-30, -1), 0.0, 0.0),
        FieldVector::new(0.0, 0.0, 0.0),
        FieldVector::new(f64::NAN, 1.0, 0.0),
    ]);
    let bdots = [FieldVector::new(0.0, 0.0, 0.0), FieldVector::new(0.0, 1e8, 0.0), FieldVector::new(1.0, -2.0, 0.5)];
    for &q in &qs[0..4] {
        for &r in &rs {
            for &bv in &betas {
                let ins = cat(&[vec![f(q)], v3(r), v3(bv)]);
                m.row(
                    "relativity::lienard_wiechert_scalar_potential",
                    Class::A,
                    &ins,
                    &res_f(lienard_wiechert_scalar_potential(black_box(q), &r, &bv)),
                );
                m.row(
                    "relativity::lienard_wiechert_vector_potential",
                    Class::A,
                    &ins,
                    &res_v(lienard_wiechert_vector_potential(black_box(q), &r, &bv)),
                );
                for &bd in &bdots {
                    let ins4 = cat(&[ins.clone(), v3(bd)]);
                    m.row("relativity::lienard_wiechert_e", Class::A, &ins4, &res_v(lienard_wiechert_e(black_box(q), &r, &bv, &bd)));
                }
            }
        }
    }
    for &n in &vs {
        for &e in vs.iter().step_by(2) {
            m.row("relativity::lienard_wiechert_b", Class::A, &cat(&[v3(n), v3(e)]), &v3(lienard_wiechert_b(black_box(&n), &e)));
        }
    }

    // Larmor: pairs over (q, a); the relativistic form over velocities and accelerations.
    let mut la = specials();
    la.extend([1.602_176_634e-19, 1e10, 1e20, 1e-300, 1e300, -3.0]);
    for &q in &la {
        for &a in &la {
            m.row("relativity::larmor_power", Class::A, &[f(q), f(a)], &[f(larmor_power(black_box(q), a))]);
        }
    }
    for &q in &qs {
        for &v in &bs {
            for &a in &vs {
                let ins = cat(&[vec![f(q)], v3(v), v3(a)]);
                m.row(
                    "relativity::relativistic_larmor_power",
                    Class::A,
                    &ins,
                    &res_f(relativistic_larmor_power(black_box(q), &v, &a)),
                );
            }
        }
    }
    m
}
