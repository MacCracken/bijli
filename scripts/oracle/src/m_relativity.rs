//! `relativity.tsv` — rust-old/src/relativity.rs (M2 part 1d-1): the Lorentz factor, the EM tensor, the field
//! transforms, four-vectors, the retarded and Liénard–Wiechert potentials and fields, and the Larmor powers; and
//! (part 1d-2) the retarded-time solver and `lienard_wiechert_fields` over named trajectories.
//!
//! Every row is class A (+ − × ÷ and sqrt in Rust's order) but the solver's over the circular orbit and the
//! oscillation, whose trajectories call sin and cos: class C. A `Result` is its value, or `E:<Variant>` and the error's
//! Display text (ADR 0002 row 4). A vector is three `f:` tokens, a four-vector four (ct, x, y, z), an `EmTensor` its
//! sixteen components row-major, a pair of vectors six.
//!
//! The grids: speeds through c (0.1c, 0.9c, 0.99c, (1 − 1e-9)c and c ulp by ulp, both signs, the specials, binades,
//! seeded); field pairs over a stride of field.rs's vectors and its edge vectors; tensors from those fields and
//! arbitrary ones; boost velocities at β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9} along x, y and the diagonal, at and past c;
//! four-vectors over the specials and physical events; Liénard–Wiechert geometries with the 1e-30 cuts of R and of
//! 1 − n̂⋅β, and β → 1; and (M2 part 2c-4) the geometries at κ = 1 − n̂⋅β's own zero (`emit_kappa_edges`).

use crate::grids::{Rng, knife_edge, pow2, specials};
use crate::m_field::{cat, edge_vectors, res_f, res_v, v3, vectors};
use crate::m_material::{dedup, ulps};
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::field::{FieldVector, SPEED_OF_LIGHT};
use bijli::relativity::{
    EmTensor, FourVector, beta, larmor_power, lienard_wiechert_b, lienard_wiechert_e, lienard_wiechert_fields,
    lienard_wiechert_scalar_potential, lienard_wiechert_vector_potential, lorentz_factor, lorentz_transform_fields,
    lorentz_transform_x, relativistic_larmor_power, retarded_scalar_potential, retarded_vector_potential,
    solve_retarded_time,
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

/// The solver's trajectories, by number, from their parameters; validate.py and the harness build the same ones.
/// 0 a static charge at (x, y, z); 1 uniform motion r₀ + v(t − t₀) (r₀, v, t₀); 2 a circular orbit in the xy-plane,
/// radius R and angular frequency ω, (R cos ωt, R sin ωt, 0); 3 hyperbolic motion along x at proper acceleration α,
/// x = (c²/α)(√(1 + (αt/c)²) − 1); 4 an oscillation along z, A sin ωt; 5 a static charge at (x, y, z) with a given
/// acceleration (aₓ, a_y, a_z). Each returns (position, velocity, acceleration).
fn trajectory(id: i64, p: &[f64], t: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    match id {
        0 => ([p[0], p[1], p[2]], [0.0; 3], [0.0; 3]),
        1 => {
            let dt = t - p[6];
            ([p[0] + p[3] * dt, p[1] + p[4] * dt, p[2] + p[5] * dt], [p[3], p[4], p[5]], [0.0; 3])
        }
        2 => {
            let (r, w) = (p[0], p[1]);
            let (s, c) = (w * t).sin_cos();
            let a2 = r * w * w;
            ([r * c, r * s, 0.0], [-r * w * s, r * w * c, 0.0], [-a2 * c, -a2 * s, 0.0])
        }
        3 => {
            let a = p[0];
            let u = a * t / C;
            let g = (1.0 + u * u).sqrt();
            ([C * C / a * (g - 1.0), 0.0, 0.0], [a * t / g, 0.0, 0.0], [a / (g * g * g), 0.0, 0.0])
        }
        4 => {
            let (amp, w) = (p[0], p[1]);
            let (s, c) = (w * t).sin_cos();
            ([0.0, 0.0, amp * s], [0.0, 0.0, amp * w * c], [0.0, 0.0, -(amp * w * w) * s])
        }
        _ => ([p[0], p[1], p[2]], [0.0; 3], [p[3], p[4], p[5]]),
    }
}

/// One solver case: the trajectory, its parameters, the observer, the observation time and max_iterations.
type SolveCase = (i64, Vec<f64>, [f64; 3], f64, usize);

/// One fields case: the charge, the trajectory, its parameters, the observer and the observation time.
type FieldsCase = (f64, i64, Vec<f64>, [f64; 3], f64);

fn traj_class(id: i64) -> Class {
    if id == 2 || id == 4 { Class::C } else { Class::A }
}

/// Part 1d-2: `solve_retarded_time` and `lienard_wiechert_fields` over the trajectories. A solver row's inputs are the
/// trajectory number, the observer, the time, max_iterations as its 64-bit word (usize::MAX is `i:-1`) and the
/// parameters (`v:<k>`); a fields row's the charge, the trajectory number, the observer, the time and the parameters.
fn emit_solver(m: &mut Module, rng: &mut Rng) {
    let mut cases: Vec<SolveCase> = Vec::new();
    let origin = vec![0.0, 0.0, 0.0];
    for t in [1e-8, 1e-3, 1.0, 10.0, 1e4, 0.0, -5.0, 1e300, f64::NAN, f64::INFINITY] {
        cases.push((0, origin.clone(), [1.0, 0.0, 0.0], t, 50));
    }
    for (obs, k) in [([0.0, 1.0, 0.0], 50), ([3.0, -4.0, 12.0], 50), ([0.0, 0.0, 0.0], 50), ([f64::NAN, 0.0, 0.0], 50)] {
        cases.push((0, origin.clone(), obs, 1e-8, k));
    }
    for k in [0, 1, 2, 3, usize::MAX] {
        cases.push((0, vec![0.25, 0.0, 0.0], [1.0, 0.0, 0.0], 1e-6, k));
    }
    cases.push((0, vec![f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 100));
    // Uniform motion: the Rust tests' cases, speeds up to 1 − 1e-9 both ways, superluminal ones.
    let vt = [0.8 * C, -0.3 * C, 0.1 * C];
    for (t0, t) in [(0.0, 2e-9), (1.0, 1.0 + 3e-9), (-5.0, -5.0 + 1e-8)] {
        cases.push((1, vec![0.1, -0.2, 0.05, vt[0], vt[1], vt[2], t0], [1.0, 0.5, 0.2], t, 50));
    }
    for b in [0.1, 0.5, 0.9, 0.99, 1.0 - 1e-9, -0.5, -0.99, -(1.0 - 1e-9)] {
        for t in [1e-8, 1.0] {
            cases.push((1, vec![0.0, 0.0, 0.0, b * C, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0], t, 50));
            cases.push((1, vec![0.0, 0.0, 0.0, b * C, 0.0, 0.0, t], [0.3, 0.8, -0.4], t, 50));
        }
    }
    for k in [0, 1, 50] {
        cases.push((1, vec![0.0, 0.0, 0.0, 0.5 * C, 0.0, 0.0, 0.0], [1.0, 0.3, 0.0], 1e-8, k));
    }
    cases.push((1, vec![0.0, 0.0, 0.0, 0.5 * C, 0.0, 0.0, 1.0], [1.0, 1.0, 0.0], 1.0, 100));
    cases.push((1, vec![0.0, 0.0, 0.0, -2.0 * C, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 100));
    cases.push((1, vec![0.0, 0.0, 0.0, 1.5 * C, 0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 1e-8, 100));
    // The β = 0.99 circular orbit (the Rust tests' three observers and the bench's), slower orbits, seeded observers.
    let w99 = 0.99 * C;
    for (obs, t) in [([0.48, -1.02, 0.0], 8.37e-8), ([1.23, -0.59, 0.01], 2e-10), ([0.89, 0.805, -0.081], 6.93e-8)] {
        for k in [50, 5] {
            cases.push((2, vec![1.0, w99], obs, t, k));
        }
    }
    for b in [0.5, 0.9] {
        cases.push((2, vec![1.0, b * C], [0.5, 2.0, 0.3], 1e-8, 50));
    }
    for _ in 0..12 {
        let obs = [rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0), rng.uniform(-1.0, 1.0)];
        cases.push((2, vec![1.0, w99], obs, rng.uniform(0.0, 1e-7), 50));
    }
    // Hyperbolic motion and the oscillation.
    for a in [9.8, 1e10, 1e15] {
        for t in [0.0, 1e-6, 1.0, 1e3] {
            cases.push((3, vec![a], [1.0, 2.0, 0.0], t, 50));
        }
    }
    for (amp, w) in [(1e-3, 1e9), (0.1, 2e9), (1e-6, 1e12)] {
        for t in [0.0, 1e-9, 3.7e-8] {
            cases.push((4, vec![amp, w], [0.5, 0.0, 0.2], t, 50));
        }
    }
    for (id, ps, obs, t, k) in &cases {
        let word = i64::from_ne_bytes(k.to_ne_bytes());
        let mut ins = vec![crate::tsv::i(*id), f(obs[0]), f(obs[1]), f(obs[2]), f(*t), crate::tsv::i(word)];
        ins.push(Tok::V(ps.len()));
        ins.extend(ps.iter().map(|&x| f(x)));
        let pp = ps.clone();
        let idv = *id;
        let r = solve_retarded_time(black_box(*obs), *t, move |tt| {
            let (a, b, _) = trajectory(idv, &pp, tt);
            (a, b)
        }, *k);
        m.row("relativity::solve_retarded_time", traj_class(*id), &ins, &res_f(r));
    }
    // lienard_wiechert_fields over the converging cases, the singular observer and an infinite acceleration.
    let mut fields: Vec<FieldsCase> = Vec::new();
    for (id, ps, obs, t, k) in &cases {
        if *k == 50 && t.is_finite() && t.abs() < 1e5 {
            fields.push((1.602_176_634e-19, *id, ps.clone(), *obs, *t));
        }
    }
    fields.push((1e-6, 0, origin.clone(), [0.0, 1.0, 0.0], 1e-8));
    fields.push((1e-6, 0, origin.clone(), [0.0, 0.0, 0.0], 1e-8));
    fields.push((1e-9, 1, vec![0.0, 0.0, 0.0, 0.9 * C, 0.0, 0.0, 2.0], [0.3, 0.8, -0.4], 2.0));
    fields.push((1.0, 5, vec![0.0, 0.0, 0.0, f64::INFINITY, 0.0, 0.0], [1.0, 0.0, 0.0], 1e-8));
    fields.push((1.0, 5, vec![0.0, 0.0, 0.0, 1e20, -1e20, 0.0], [1.0, 2.0, 3.0], 1e-8));
    for (q, id, ps, obs, t) in fields {
        let mut ins = vec![f(q), crate::tsv::i(id), f(obs[0]), f(obs[1]), f(obs[2]), f(t)];
        ins.push(Tok::V(ps.len()));
        ins.extend(ps.iter().map(|&x| f(x)));
        let r = lienard_wiechert_fields(black_box(q), obs, t, move |tt| trajectory(id, &ps, tt));
        m.row("relativity::lienard_wiechert_fields", traj_class(id), &ins, &res_pair(r));
    }
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
    // The |R| < 1e-30 cut's knife edge (bijli ADR 0002 row 17's rule, M2 part 2c-1): grids.rs's points whose exact |R|²
    // lies within two ulps of 10^-60, through the potentials and E.
    for (src, at) in knife_edge(true) {
        if src != [0.0, 0.0, 0.0] {
            continue;
        }
        let r = FieldVector::new(at[0], at[1], at[2]);
        let bv = FieldVector::new(0.1, 0.0, 0.0);
        let ins = cat(&[vec![f(1.0)], v3(r), v3(bv)]);
        m.row(
            "relativity::lienard_wiechert_scalar_potential",
            Class::A,
            &ins,
            &res_f(lienard_wiechert_scalar_potential(black_box(1.0), &r, &bv)),
        );
        m.row(
            "relativity::lienard_wiechert_vector_potential",
            Class::A,
            &ins,
            &res_v(lienard_wiechert_vector_potential(black_box(1.0), &r, &bv)),
        );
        let bd = FieldVector::new(0.0, 0.0, 0.0);
        let ins4 = cat(&[ins.clone(), v3(bd)]);
        m.row("relativity::lienard_wiechert_e", Class::A, &ins4, &res_v(lienard_wiechert_e(black_box(1.0), &r, &bv, &bd)));
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
    emit_solver(&mut m, &mut rng);
    emit_kappa_edges(&mut m);
    emit_nb_edges(&mut m);
    emit_fv_edges(&mut m);
    emit_transform_edges(&mut m);
    emit_lw_inf_edges(&mut m);
    emit_e_tail_edges(&mut m);
    m
}

/// A unit direction, uniform on the sphere (rejection in the cube; no transcendental).
fn unit3(rng: &mut Rng) -> [f64; 3] {
    loop {
        let v = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
        let n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if n2 > 0.01 && n2 <= 1.0 {
            let n = n2.sqrt();
            return [v[0] / n, v[1] / n, v[2] / n];
        }
    }
}

/// x = m·2^e exactly, for a positive finite x (m below 2^53).
fn sig_exp(x: f64) -> (i128, i32) {
    let b = x.to_bits();
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    let ex = ((b >> 52) & 0x7ff) as i32;
    let mant = i128::from(b & 0xf_ffff_ffff_ffff);
    if ex == 0 { (mant, -1074) } else { (mant | (1 << 52), ex - 1075) }
}

/// The two potentials (when `pots`) and E at one geometry.
fn lw_rows(m: &mut Module, q: f64, r: FieldVector, bv: FieldVector, bd: FieldVector, pots: bool) {
    let ins = cat(&[vec![f(q)], v3(r), v3(bv)]);
    if pots {
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
    }
    let ins4 = cat(&[ins, v3(bd)]);
    m.row("relativity::lienard_wiechert_e", Class::A, &ins4, &res_v(lienard_wiechert_e(black_box(q), &r, &bv, &bd)));
}

/// An exact midpoint of κ (`kappa`) or of R_eff/s along the Pythagorean direction (1, b, c) of norm d: β₂ and β₃ on a
/// 2^-10 grid, β₁ solving for the midpoint μ above y in units of 2^-64 (None unless β₁ is a double and |β| < 0.99).
fn lw_tie(rng: &mut Rng, kappa: bool) -> Option<(FieldVector, FieldVector)> {
    const DIRS: [(i128, i128, i128); 5] = [(2, 2, 3), (4, 8, 9), (12, 12, 17), (6, 18, 19), (8, 4, 9)];
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let (bq, cq, d) = DIRS[rng.int(0, 4) as usize];
    #[allow(clippy::cast_precision_loss)]
    let df = d as f64;
    let y = if kappa { rng.uniform(0.05, 1.5) } else { rng.uniform(0.05, 1.5) * df };
    let one: i128 = 1 << 64;
    let (ym, ye) = sig_exp(y);
    let mu = (ym << (64 + ye)) + (1 << (63 + ye));  // (y + ulp/2)·2^64: y ≥ 2^-5, so ye ≥ −57
    let b2i = i128::from(rng.int(-64, 64));
    let b3i = i128::from(rng.int(-64, 64));
    let target = if kappa { d * (one - mu) } else { d * one - mu };
    let b1i = target - bq * (b2i << 54) - cq * (b3i << 54);
    #[allow(clippy::cast_precision_loss)]
    let b1 = b1i as f64 * pow2(-64);
    #[allow(clippy::cast_possible_truncation)]
    let back = (b1 * pow2(64)) as i128;
    #[allow(clippy::cast_precision_loss)]
    let (b2, b3) = (b2i as f64 * pow2(-10), b3i as f64 * pow2(-10));
    if back != b1i || b1 * b1 + b2 * b2 + b3 * b3 >= 0.99 {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    let sc = pow2(rng.int(-20, 20) as i32);
    #[allow(clippy::cast_precision_loss)]
    let r = FieldVector::new(sc, bq as f64 * sc, cq as f64 * sc);
    Some((r, FieldVector::new(b1, b2, b3)))
}

/// M2 part 2c-4 (bijli ADR 0002 row 17): Liénard–Wiechert geometries at κ = 1 − n̂⋅β's own zero, through both
/// potentials and E. R along β within small angles as |β| → 1; an ultra-relativistic charge seen near its axis (κ
/// about 6e-17, which Rust's 1 − n̂⋅β rounds to 0); κ's 1e-30 cut over a sweep of β₃'s ulps across 10^-30 (E alone;
/// D = 1 − |β|² near 2·10^-30 built in integers); exact midpoints of κ and R_eff along Pythagorean directions; the
/// near/far edge n̂⋅β = ½; extreme scales and spreads. Its own seed, after every older row.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn emit_kappa_edges(m: &mut Module) {
    let mut rng = Rng::new(0x2C4_4A99_5EED);
    let bdots = [FieldVector::new(0.0, 0.0, 0.0), FieldVector::new(1.0, -2.0, 0.5), FieldVector::new(0.0, 1e8, 0.0)];
    let qs = [1.0, -1.602_176_634e-19];
    // Near the cone: R along n̂ at 2^±90, β = |β|·(n̂ + θe)/|n̂ + θe|, |β| = 1 − u·2^-k.
    for i in 0..400 {
        let n = unit3(&mut rng);
        let sc = pow2(rng.int(-90, 90) as i32);
        let r = FieldVector::new(n[0] * sc, n[1] * sc, n[2] * sc);
        let e = unit3(&mut rng);
        let th = if i % 5 == 0 { 0.0 } else { pow2(-(rng.int(0, 60) as i32)) * rng.uniform(1.0, 2.0) };
        let sp = 1.0 - pow2(-(rng.int(1, 52) as i32)) * rng.uniform(1.0, 2.0) * 0.5;
        let w = [n[0] + th * e[0], n[1] + th * e[1], n[2] + th * e[2]];
        let wn = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        let bv = FieldVector::new(w[0] * sp / wn, w[1] * sp / wn, w[2] * sp / wn);
        lw_rows(m, qs[i % 2], r, bv, bdots[i % 3], true);
    }
    // An ultra-relativistic charge (|β| = 1 − m·2^-53) seen just off its axis, the parts rotated.
    for i in 0..100 {
        let x = pow2(-(rng.int(20, 40) as i32)) * rng.uniform(1.0, 2.0);
        let y = pow2(-(rng.int(24, 29) as i32)) * rng.uniform(1.0, 2.0);
        let a = [1.0, 3.0, 1e-20, 1e20][i % 4];
        let b1 = 1.0 - pow2(-53) * rng.int(1, 8) as f64;
        let (rv, bv) = ([a, x * a, 0.0], [b1, y, 0.0]);
        let k = i % 3;
        let r = FieldVector::new(rv[k], rv[(k + 1) % 3], rv[(k + 2) % 3]);
        let bvv = FieldVector::new(bv[k], bv[(k + 1) % 3], bv[(k + 2) % 3]);
        lw_rows(m, qs[i % 2], r, bvv, bdots[i % 3], true);
    }
    // κ's 1e-30 cut: β₁ = 1 − m·2^-53 (T₁ = 1 − β₁² = (m·2^54 − m²)·2^-106), β₂ leaving T₂ = T₁ − β₂² a little above
    // 2·10^-30 (exactly, in integers), β₃ = √(T₂ − 2·10^-30) within a few ulps of the crossing, swept ±32 ulps;
    // R = β·2^s.
    for j in 0..4 {
        let m1 = rng.int(1, 4);
        let b1 = 1.0 - pow2(-53) * m1 as f64;
        let mi = i128::from(m1);
        let t1 = (mi << 54) - mi * mi;
        let gap = rng.uniform(3e-31, 9e-31);
        let b2 = (t1 as f64 * pow2(-106) - 2e-30 - gap).sqrt();
        let (m2, e2) = sig_exp(b2);
        let t2 = (t1 << (-106 - 2 * e2)) - m2 * m2;
        let b3 = (t2 as f64 * pow2(2 * e2) - 2e-30).sqrt();
        let s = pow2(rng.int(-60, 60) as i32);
        for k in -32..=32 {
            let b3k = ulps(b3, k);
            let bb = [[b1, b2, b3k], [b2, b3k, b1], [b3k, b1, b2], [b1, b3k, b2]][j];
            let bv = FieldVector::new(bb[0], bb[1], bb[2]);
            let r = FieldVector::new(bb[0] * s, bb[1] * s, bb[2] * s);
            lw_rows(m, 1.0, r, bv, bdots[1], false);
        }
    }
    // Exact midpoints: 60 of κ, 60 of R_eff.
    for kappa in [true, false] {
        let mut got = 0;
        while got < 60 {
            if let Some((r, bv)) = lw_tie(&mut rng, kappa) {
                lw_rows(m, 1.0, r, bv, bdots[got % 3], true);
                got += 1;
            }
        }
    }
    // The near/far edge: n̂⋅β = ½(1 ± 2^-k).
    for i in 0..100 {
        let r = [rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0)];
        let n = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        let sg = if i % 2 == 0 { 1.0 } else { -1.0 };
        let fr = 0.5 * (1.0 + sg * pow2(-(rng.int(1, 60) as i32)));
        let bv = FieldVector::new(r[0] / n * fr, r[1] / n * fr, r[2] / n * fr);
        lw_rows(m, 1.0, FieldVector::new(r[0], r[1], r[2]), bv, bdots[i % 3], true);
    }
    // Extremes: R past DBL_MAX's root, near 10^-30, wide exponent spreads, subnormal parts; β zero, subnormal, at
    // 1 − ε.
    let rx = [
        [1e308, 1e308, -1e308],
        [1.7e308, 1.7e308, 1.7e308],
        [1e-30, 0.0, 0.0],
        [6e-31, 6e-31, 6e-31],
        [1e300, 1e-300, 5e-324],
        [5e-324, 1.0, 0.0],
        [1.0, 5e-324, -5e-324],
        [pow2(-100), pow2(1000), 1.0],
    ];
    let bx = [
        [0.0, 0.0, 0.0],
        [5e-324, 0.0, -5e-324],
        [0.999, 1e-200, 0.0],
        [1.0 - pow2(-53), 0.0, 0.0],
        [-0.6, 0.7, 0.3],
        [0.0, 1.0 - pow2(-53), pow2(-27)],
    ];
    for (i, r) in rx.iter().enumerate() {
        for (j, b) in bx.iter().enumerate() {
            let rv = FieldVector::new(r[0], r[1], r[2]);
            lw_rows(m, 1.0, rv, FieldVector::new(b[0], b[1], b[2]), bdots[(i + j) % 3], true);
        }
    }
}

/// M2 part 2d-1 (bijli ADR 0002 row 17, the user's answer of 2026-11-13): `lienard_wiechert_e` where n̂ − β cancels
/// part by part: R along Pythagorean directions (a, b, c)/d at scales 2^s, β = f·RN(n̂) for f from ½ to 1 − 2^-52, the
/// parts rotated, so each n̂ᵢ − βᵢ lands within ulps of 0 for f → 1. After every older row.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn emit_nb_edges(m: &mut Module) {
    let dirs = [
        (3.0, 4.0, 0.0, 5.0),
        (1.0, 2.0, 2.0, 3.0),
        (2.0, 3.0, 6.0, 7.0),
        (1.0, 4.0, 8.0, 9.0),
        (4.0, 4.0, 7.0, 9.0),
    ];
    let bd = FieldVector::new(1.0, -2.0, 0.5);
    for (j, &(a, b, c, d)) in dirs.iter().enumerate() {
        for s in [-40, 0, 7, 300] {
            let sc = pow2(s);
            for fr in [0.5, 0.9, 1.0 - pow2(-20), 1.0 - pow2(-52)] {
                let n = [a / d, b / d, c / d];
                let k = j % 3;
                let p = [n[k], n[(k + 1) % 3], n[(k + 2) % 3]];
                let r = FieldVector::new(p[0] * d * sc, p[1] * d * sc, p[2] * d * sc);
                let bv = FieldVector::new(p[0] * fr, p[1] * fr, p[2] * fr);
                lw_rows(m, 1.0, r, bv, bd, false);
            }
        }
    }
}

/// M2 part 2d-3c (bijli ADR 0002 row 18's rule): the four-vector and tensor chains at their range edges. Four-vectors
/// whose products pass DBL_MAX and cancel (light-like at 10^300, DBL_MAX parts), whose squares overflow alone
/// (10^155, 10^200) or fall below 2^-1022 (10^-160, 10^-300), through `interval_sq`, `spatial_magnitude`, `dot` and
/// `boost_x` (with speeds whose β is subnormal: 10^-310, 2^-1074, 10^-300). Tensors from fields whose invariants'
/// squares and products pass DBL_MAX and cancel (a radiation field at 10^200, E ∥ B at 10^150 and 10^154/10^155), whose
/// E/c is subnormal under a vast B, whose B² is subnormal, and DBL_MAX parts, through both invariants and the
/// four-force (four-velocities at 10^300 and 10^10, the latter under q = 10^-300). After every older row.
fn emit_fv_edges(m: &mut Module) {
    let big = f64::MAX;
    let fvs = [
        FourVector::new(1e300, 1e300, 0.0, 0.0),
        FourVector::new(1e300, 0.6e300, 0.8e300, 0.0),
        FourVector::new(big, big, 0.0, 0.0),
        FourVector::new(-big, 0.5 * big, big, 0.25 * big),
        FourVector::new(1e200, 1e200, 1e200, 1e200),
        FourVector::new(1e155, 1e155, 0.0, 0.0),
        FourVector::new(1e-160, 1e-161, 0.0, 0.0),
        FourVector::new(3e-160, 1e-160, 2e-160, 2e-160),
        FourVector::new(1e-300, 1e-300, 1e-300, 0.0),
        FourVector::new(1.0, 1e-170, 1e-170, 1e-170),
    ];
    for &a in &fvs {
        let ins = fv4(a);
        m.row("FourVector::interval_sq", Class::A, &ins, &[f(black_box(a).interval_sq())]);
        m.row("FourVector::spatial_magnitude", Class::A, &ins, &[f(black_box(a).spatial_magnitude())]);
        for &b in &fvs {
            m.row("FourVector::dot", Class::A, &cat(&[ins.clone(), fv4(b)]), &[f(black_box(a).dot(&b))]);
        }
        for u in [1e-310, pow2(-1074), 1e-300, 0.5 * C, -0.9 * C, 0.0] {
            m.row("FourVector::boost_x", Class::A, &cat(&[ins.clone(), vec![f(u)]]), &res_4(black_box(a).boost_x(u)));
        }
    }
    let fields = [
        (FieldVector::new(0.0, C * 1e200, 0.0), FieldVector::new(0.0, 0.0, 1e200)),
        (FieldVector::new(C * 1e150, 0.0, 0.0), FieldVector::new(1e150, 0.0, 0.0)),
        (FieldVector::new(1e154, 0.0, 0.0), FieldVector::new(1e155, 0.0, 0.0)),
        (FieldVector::new(1e-303, 2e-303, 0.0), FieldVector::new(1e300, -1e300, 5e299)),
        (FieldVector::new(C * 1e-160, 0.0, 0.0), FieldVector::new(2e-160, 0.0, 0.0)),
        (FieldVector::new(big, big, 0.0), FieldVector::new(big, 0.0, big)),
        (FieldVector::new(1e300, -1e300, 1e300), FieldVector::new(1e300, 1e300, 0.0)),
    ];
    let us = [
        (1.0, FourVector::new(1e300, 1e300, 0.0, 0.0)),
        (1e-300, FourVector::new(1e10, 1e10, 1e10, 1e10)),
        (-1.6e-19, FourVector::new(C / 0.6, 0.8 * C / 0.6, 0.0, 0.0)),
    ];
    for (e, b) in fields {
        let t = EmTensor::from_fields(black_box(&e), &b);
        let ins = tt(&t);
        m.row("EmTensor::first_invariant", Class::A, &ins, &[f(black_box(&t).first_invariant())]);
        m.row("EmTensor::second_invariant", Class::A, &ins, &[f(black_box(&t).second_invariant())]);
        for (q, u) in us {
            let ins4 = cat(&[ins.clone(), vec![f(q)], fv4(u)]);
            m.row("EmTensor::four_force", Class::A, &ins4, &fv4(black_box(&t).four_force(q, &u)));
        }
    }
}

/// M2 part 2d-3d (bijli ADR 0002 row 18's rule): the field transforms and the relativistic Larmor power at their range
/// edges. Boosts along x at γ = 2, 1.25 and near c, diagonal, and with subnormal or tiny parts; fields along the boost
/// at 10^308 and DBL_MAX (E′∥ = E, where Rust's γE overflowed), past DBL_MAX/c across it, subnormal, tiny, ordinary
/// and zero. lorentz_transform_x at speeds whose v·E/c² passes DBL_MAX or falls below 2^-1022. The Larmor power for
/// charges from 0 and 10^-300 to 10^300 under accelerations whose a² passes DBL_MAX (10^170) or falls below 2^-1022.
/// After every older row.
fn emit_transform_edges(m: &mut Module) {
    let big = f64::MAX;
    let vs = [
        FieldVector::new(0.866_025_403_784_438_6 * C, 0.0, 0.0),
        FieldVector::new(0.6 * C, 0.0, 0.0),
        FieldVector::new(0.5 * C, 0.5 * C, 0.0),
        FieldVector::new(1e-310, 0.0, 0.0),
        FieldVector::new(1e-100, 1e-100, 0.0),
        FieldVector::new(C * (1.0 - 1e-15), 0.0, 0.0),
    ];
    let es = [
        FieldVector::new(1e308, 0.0, 0.0),
        FieldVector::new(big, 0.0, 0.0),
        FieldVector::new(0.0, 1e301, 0.0),
        FieldVector::new(1e-310, 5e-324, 0.0),
        FieldVector::new(1e-160, 1e-160, 1e-160),
        FieldVector::new(1.0, 2.0, 3.0),
    ];
    let bs = [
        FieldVector::new(0.0, 0.0, 1e301),
        FieldVector::new(1e-310, 0.0, 0.0),
        FieldVector::new(1e300, 1e300, 0.0),
        FieldVector::new(0.0, 0.0, 0.0),
        FieldVector::new(1e-160, 0.0, 1e-160),
    ];
    for &e in &es {
        for &b in &bs {
            for &v in &vs {
                let ins = cat(&[v3(e), v3(b), v3(v)]);
                m.row(
                    "relativity::lorentz_transform_fields",
                    Class::A,
                    &ins,
                    &res_pair(lorentz_transform_fields(black_box(&e), &b, &v)),
                );
            }
            for u in [0.8 * C, 0.5 * C, 1e-310, 1e-150, 2e8] {
                let ins = cat(&[v3(e), v3(b), vec![f(u)]]);
                m.row("relativity::lorentz_transform_x", Class::A, &ins, &res_pair(lorentz_transform_x(black_box(&e), &b, u)));
            }
        }
    }
    let accs = [
        FieldVector::new(1e170, 0.0, 0.0),
        FieldVector::new(1e150, 1e150, 0.0),
        FieldVector::new(1e-200, 0.0, 0.0),
        FieldVector::new(5e-324, 0.0, 1e-310),
        FieldVector::new(1e300, 1e-300, 0.0),
    ];
    for q in [1.6e-19, 1e-200, 1e300, 0.0, 1e-300] {
        for &v in &vs {
            for &a in &accs {
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
}

/// M2 part 2d-3e (bijli ADR 0002 rows 27, amended, and 18): the Liénard–Wiechert functions on an R with an infinite part
/// (now refused, as a NaN R is; a NaN beside the ∞ keeps the NaN text) and the potentials where R_eff passes DBL_MAX
/// (R near DBL_MAX against β, φ and A from the 53-bit R_eff) or μ₀qc leaves the range (q = DBL_MAX, subnormal q).
/// After every older row.
fn emit_lw_inf_edges(m: &mut Module) {
    let inf = f64::INFINITY;
    let big = f64::MAX;
    let zero = FieldVector::new(0.0, 0.0, 0.0);
    let betas = [zero, FieldVector::new(0.5, 0.0, 0.0), FieldVector::new(0.0, 0.9, 0.0)];
    let infs = [
        FieldVector::new(inf, 0.0, 0.0),
        FieldVector::new(-inf, 1.0, 0.0),
        FieldVector::new(inf, inf, 0.0),
        FieldVector::new(1.0, -inf, 2.0),
        FieldVector::new(inf, f64::NAN, 0.0),
    ];
    for &r in &infs {
        for &bv in &betas {
            lw_rows(m, 1.0, r, bv, zero, true);
            lw_rows(m, 1.0, r, bv, FieldVector::new(1e10, 0.0, 0.0), false);
        }
    }
    let bigs = [FieldVector::new(big, big, 0.0), FieldVector::new(1e308, 1e308, 1e308), FieldVector::new(-big, 0.0, 0.0)];
    let bbs = [FieldVector::new(-0.9, 0.0, 0.0), FieldVector::new(0.0, 0.5, 0.0)];
    for &r in &bigs {
        for &bv in &bbs {
            for q in [1.0, 1e300, 1e-300, 0.0] {
                lw_rows(m, q, r, bv, zero, true);
            }
        }
    }
    for q in [big, 1e-310, pow2(-1074)] {
        for r in [FieldVector::new(1.0, 0.0, 0.0), FieldVector::new(1e-20, 0.0, 0.0), FieldVector::new(1e300, 0.0, 0.0)] {
            for &bv in &bbs {
                lw_rows(m, q, r, bv, zero, true);
            }
        }
    }
}

/// M2 part 2d-3f (bijli ADR 0002 row 18's rule): lienard_wiechert_e's tail at its range edges, outside its box: charges
/// of 10^300, 2^-1074 and 10^-200; R of 10^200 (R² past DBL_MAX) and just above 10^-30; β̇ of 10^300 and 10^-300 beside
/// zero. After every older row.
fn emit_e_tail_edges(m: &mut Module) {
    let zero = FieldVector::new(0.0, 0.0, 0.0);
    let rs = [
        FieldVector::new(1e200, 0.0, 0.0),
        FieldVector::new(1e200, 1e200, 0.0),
        FieldVector::new(2e-30, 0.0, 0.0),
        FieldVector::new(3e-30, 4e-30, 0.0),
    ];
    let betas = [zero, FieldVector::new(0.5, 0.0, 0.0), FieldVector::new(0.0, 0.9, 0.0)];
    let bds = [zero, FieldVector::new(0.0, 1e300, 0.0), FieldVector::new(1e-300, 0.0, 0.0)];
    for q in [1e300, pow2(-1074), 1e-200] {
        for &r in &rs {
            for &bv in &betas {
                for &bd in &bds {
                    lw_rows(m, q, r, bv, bd, false);
                }
            }
        }
    }
}
