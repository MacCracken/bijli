//! `polarization.tsv` — rust-old/src/polarization.rs (M3 part 1a-1): Jones vectors and matrices, Stokes vectors and
//! Mueller matrices.
//!
//! Classes: the functions that call `cos`, `sin`, `atan2` or `hypot` are class C (`validate.py` writes the correctly
//! rounded readings to `polarization.cr.tsv`); the rest are class A (+ − × ÷ and sqrt in Rust's order, through hisab's
//! `Complex`). A Jones vector is its four parts (x.re, x.im, y.re, y.im), a Jones matrix its eight (row-major, re
//! before im), a Stokes vector its four and a Mueller matrix its sixteen (row-major). A `Result` is its value, or
//! `E:<Variant>` and the error's Display text (ADR 0002 row 4).
//!
//! The grids: complex amplitudes (the specials' pairs, the unit circle, seeded uniform and log-uniform values, the
//! scales where |z|² under- and overflows) paired into Jones vectors; the named states and elements at seeded angles,
//! seeded complex matrices and matrices at the range edges; angles (the specials, multiples of π/8, seeded, tiny and
//! huge) and retardances; Stokes vectors from the specials, seeded values, the named states and the scales where
//! S₁² + S₂² + S₃² under- and overflows; Mueller matrices from the named elements and seeded entries. Each value appears
//! once per function: a row's (function, inputs) is its key in the mapping files.
//!
//! M3 part 2a-1 (ADR 0002 row 33) appends the waveplates' edge grid after every earlier row: fast axes within two ulps
//! of the π/4 family (where c² − s² cancels), tiny ones, and at DBL_MAX/2 and 2^1023 (where 2θ overflows); retardances
//! across cos δ = 1/2 (where 1 − cos δ starts to cancel) and down to 1e-300.

use crate::grids::{Rng, pow2, specials};
use crate::m_field::res_f;
use crate::tsv::{Class, Module, Tok, err, f, s};
use bijli::polarization::{Complex, JonesMatrix, JonesVector, MuellerMatrix, StokesVector};
use std::collections::BTreeSet;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
use std::hint::black_box;

fn jv(v: JonesVector) -> Vec<Tok> {
    vec![f(v.x.re), f(v.x.im), f(v.y.re), f(v.y.im)]
}

fn jm(m: JonesMatrix) -> Vec<Tok> {
    let mut out = Vec::new();
    for row in &m.m {
        for z in row {
            out.push(f(z.re));
            out.push(f(z.im));
        }
    }
    out
}

fn sv(v: StokesVector) -> Vec<Tok> {
    v.s.iter().map(|&x| f(x)).collect()
}

fn mm(m: MuellerMatrix) -> Vec<Tok> {
    m.m.iter().flat_map(|row| row.iter().map(|&x| f(x))).collect()
}

fn cx(z: Complex) -> Vec<Tok> {
    vec![f(z.re), f(z.im)]
}

fn res_jv(r: bijli::error::Result<JonesVector>) -> Vec<Tok> {
    match r {
        Ok(v) => jv(v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// The bits of every f64 in a token list, for keeping each row's inputs once.
fn key(name: &str, toks: &[Tok]) -> String {
    let mut k = String::from(name);
    for t in toks {
        if let Tok::F(x) = t {
            k.push_str(&format!(" {:016x}", x.to_bits()));
        }
    }
    k
}

/// A module whose rows are kept once per (function, inputs).
struct Rows {
    m: Module,
    seen: BTreeSet<String>,
}

impl Rows {
    fn row(&mut self, name: &'static str, class: Class, ins: &[Tok], outs: &[Tok]) {
        if self.seen.insert(key(name, ins)) {
            self.m.row(name, class, ins, outs);
        }
    }
}

/// Complex amplitudes: the specials' pairs, the unit circle, seeded values, and the scales where |z|² leaves the range.
fn complexes(rng: &mut Rng) -> Vec<Complex> {
    let sp = specials();
    let mut zs = Vec::new();
    for &re in &sp {
        for &im in &sp {
            zs.push(Complex::new(re, im));
        }
    }
    for k in 0..16 {
        let t = f64::from(k) * (PI / 8.0);
        zs.push(Complex::new(black_box(t).cos(), black_box(t).sin()));
    }
    for e in [-320.0, -300.0, -200.0, -160.0, -155.0, 150.0, 155.0, 160.0, 200.0, 300.0, 308.0] {
        let x: f64 = 10.0_f64.powf(e);
        zs.push(Complex::new(x, 0.0));
        zs.push(Complex::new(0.0, -x));
        zs.push(Complex::new(x, x));
    }
    for _ in 0..24 {
        zs.push(Complex::new(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0)));
    }
    for _ in 0..12 {
        zs.push(Complex::new(rng.signed_log_uniform(-30.0, 30.0), rng.signed_log_uniform(-30.0, 30.0)));
    }
    zs
}

/// Jones vectors: the named states, and pairs of complex amplitudes (a short list for the pairs of vectors).
fn vectors(rng: &mut Rng) -> (Vec<JonesVector>, Vec<JonesVector>) {
    let zs = complexes(rng);
    let named = vec![
        JonesVector::horizontal(),
        JonesVector::vertical(),
        JonesVector::diagonal(),
        JonesVector::anti_diagonal(),
        JonesVector::right_circular(),
        JonesVector::left_circular(),
        JonesVector::linear(0.3),
    ];
    let mut all = named.clone();
    let pick: Vec<Complex> = zs.iter().step_by(3).copied().collect();
    for &a in &pick {
        for &b in pick.iter().step_by(2) {
            all.push(JonesVector::new(a, b));
        }
    }
    for &a in zs.iter().skip(100) {
        all.push(JonesVector::new(a, zs[(zs.len() / 2 + all.len()) % zs.len()]));
    }
    let mut short = named;
    for _ in 0..20 {
        let a = Complex::new(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0));
        let b = Complex::new(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0));
        short.push(JonesVector::new(a, b));
    }
    for e in [-300.0, -200.0, -160.0, 160.0, 200.0, 300.0] {
        let x: f64 = 10.0_f64.powf(e);
        short.push(JonesVector::new(Complex::new(x, 0.5 * x), Complex::new(-0.25 * x, x)));
    }
    short.push(JonesVector::new(Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)));
    short.push(JonesVector::new(Complex::new(f64::NAN, 0.0), Complex::new(1.0, 0.0)));
    short.push(JonesVector::new(Complex::new(f64::INFINITY, 0.0), Complex::new(0.0, 0.0)));
    (all, short)
}

/// Angles: the specials, multiples of π/8, seeded, tiny and huge.
fn angles(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    for k in -16..=16 {
        xs.push(f64::from(k) * (PI / 8.0));
    }
    xs.extend([0.3, 0.37, 0.5, 0.7, -1.2, 1e-10, 1e-300, 1e10, 1e20, 1e300, FRAC_PI_2, FRAC_PI_4]);
    for _ in 0..24 {
        xs.push(rng.uniform(-7.0, 7.0));
    }
    xs
}

/// Retardances: the specials, the named plates' and seeded values.
fn retardances(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([FRAC_PI_2, PI, 2.0 * PI, 0.4, 2.9, -1.0, 1e-8, 1e10]);
    for _ in 0..8 {
        xs.push(rng.uniform(-7.0, 7.0));
    }
    xs
}

/// Jones matrices: the named elements at a few angles, seeded complex matrices and matrices at the range edges.
fn matrices(rng: &mut Rng) -> Vec<JonesMatrix> {
    let mut ms = vec![
        JonesMatrix::identity(),
        JonesMatrix::horizontal_polarizer(),
        JonesMatrix::vertical_polarizer(),
        JonesMatrix::new(Complex::new(-0.0, 0.0), Complex::new(0.0, -0.0), Complex::new(-0.0, -0.0), Complex::new(0.0, 0.0)),
    ];
    for t in [0.0, 0.3, FRAC_PI_4, -1.2] {
        ms.push(JonesMatrix::linear_polarizer(t));
        ms.push(JonesMatrix::quarter_wave_plate(t));
        ms.push(JonesMatrix::half_wave_plate(t));
        ms.push(JonesMatrix::waveplate(t, 0.4));
        ms.push(JonesMatrix::rotator(t));
    }
    for _ in 0..24 {
        let mut c = || Complex::new(rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0));
        ms.push(JonesMatrix::new(c(), c(), c(), c()));
    }
    for e in [-200.0, -160.0, 160.0, 200.0] {
        let x: f64 = 10.0_f64.powf(e);
        ms.push(JonesMatrix::new(Complex::new(x, x), Complex::new(0.0, -x), Complex::new(x, 0.5), Complex::new(-x, x)));
    }
    ms.push(JonesMatrix::new(
        Complex::new(f64::NAN, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(f64::INFINITY, 0.0),
        Complex::new(0.0, 1.0),
    ));
    ms
}

/// Reals: intensities and transmissions.
fn reals(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([0.5, 2.0, -3.0, 1e-300, 1e300, 1e-40, 0.25]);
    for _ in 0..12 {
        xs.push(rng.signed_log_uniform(-8.0, 8.0));
    }
    xs
}

/// Stokes vectors: the named states, the specials' combinations, seeded and the range edges.
fn stokes(rng: &mut Rng) -> Vec<StokesVector> {
    let mut out = vec![
        StokesVector::unpolarized(1.0),
        StokesVector::horizontal(1.0),
        StokesVector::vertical(2.0),
        StokesVector::right_circular(1.0),
        StokesVector::left_circular(0.5),
        StokesVector::default(),
        StokesVector::new(1e-40, 0.5e-40, 0.0, 0.0),
        StokesVector::new(1e-170, 0.6e-170, 0.0, 0.8e-170),
        StokesVector::new(1e170, 0.6e170, 0.0, -0.8e170),
        StokesVector::new(-1.0, 0.0, 0.0, 0.0),
        StokesVector::new(1.0, -0.0, 0.0, -0.0),
        StokesVector::new(2.0, 1.0, 1.0, 1.0),
    ];
    let sp = specials();
    for &a in &sp {
        out.push(StokesVector::new(1.0, a, 0.5, -0.25));
        out.push(StokesVector::new(a, 0.5, -0.5, 0.25));
        out.push(StokesVector::new(2.0, 0.5, a, 0.0));
        out.push(StokesVector::new(2.0, 0.0, -0.5, a));
    }
    for _ in 0..40 {
        let s0 = rng.uniform(0.0, 3.0);
        out.push(StokesVector::new(s0, rng.uniform(-s0, s0), rng.uniform(-s0, s0), rng.uniform(-s0, s0)));
    }
    for _ in 0..12 {
        out.push(StokesVector::new(
            rng.signed_log_uniform(-300.0, 300.0).abs(),
            rng.signed_log_uniform(-300.0, 300.0),
            rng.signed_log_uniform(-300.0, 300.0),
            rng.signed_log_uniform(-300.0, 300.0),
        ));
    }
    out
}

/// Mueller matrices: the named elements and seeded ones.
fn muellers(rng: &mut Rng) -> Vec<MuellerMatrix> {
    let mut ms = vec![
        MuellerMatrix::identity(),
        MuellerMatrix::horizontal_polarizer(),
        MuellerMatrix::vertical_polarizer(),
        MuellerMatrix::depolarizer(),
        MuellerMatrix::attenuator(0.5),
        MuellerMatrix::right_circular_polarizer(),
        MuellerMatrix::left_circular_polarizer(),
    ];
    for t in [0.0, 0.3, FRAC_PI_4, -1.2] {
        ms.push(MuellerMatrix::linear_polarizer(t));
        ms.push(MuellerMatrix::quarter_wave_plate(t));
        ms.push(MuellerMatrix::half_wave_plate(t));
        ms.push(MuellerMatrix::waveplate(t, 0.4));
        ms.push(MuellerMatrix::rotation(t));
    }
    for _ in 0..16 {
        let mut a = [[0.0; 4]; 4];
        for row in &mut a {
            for x in row.iter_mut() {
                *x = rng.uniform(-1.0, 1.0);
            }
        }
        ms.push(MuellerMatrix::new(a));
    }
    let mut z = [[-0.0; 4]; 4];
    z[0][0] = 1.0;
    ms.push(MuellerMatrix::new(z));
    ms
}

pub fn emit() -> Module {
    let mut rows = Rows { m: Module::new("polarization"), seen: BTreeSet::new() };
    let m = &mut rows;
    let mut rng = Rng::new(0x3a_1a1);
    let (vs, short) = vectors(&mut rng);
    let ths = angles(&mut rng);
    let ds = retardances(&mut rng);
    let ms = matrices(&mut rng);
    let rs = reals(&mut rng);
    let ss = stokes(&mut rng);
    let mus = muellers(&mut rng);

    // Jones vectors.
    for (name, v) in [
        ("JonesVector::horizontal", JonesVector::horizontal()),
        ("JonesVector::vertical", JonesVector::vertical()),
        ("JonesVector::diagonal", JonesVector::diagonal()),
        ("JonesVector::anti_diagonal", JonesVector::anti_diagonal()),
        ("JonesVector::right_circular", JonesVector::right_circular()),
        ("JonesVector::left_circular", JonesVector::left_circular()),
        ("JonesVector::default", JonesVector::default()),
    ] {
        m.row(name, Class::A, &[], &jv(v));
    }
    for &t in &ths {
        m.row("JonesVector::linear", Class::C, &[f(t)], &jv(JonesVector::linear(black_box(t))));
    }
    for &v in &vs {
        m.row("JonesVector::new", Class::A, &jv(v), &jv(JonesVector::new(black_box(v.x), v.y)));
        m.row("JonesVector::intensity", Class::A, &jv(v), &[f(black_box(v).intensity())]);
        m.row("JonesVector::normalized", Class::A, &jv(v), &res_jv(black_box(v).normalized()));
    }
    for &a in &short {
        for &b in &short {
            let ins = [jv(a), jv(b)].concat();
            m.row("JonesVector::inner", Class::A, &ins, &cx(black_box(a).inner(&b)));
            m.row("JonesVector::overlap", Class::A, &ins, &res_f(black_box(a).overlap(&b)));
        }
    }

    // Jones matrices.
    for (name, j) in [
        ("JonesMatrix::identity", JonesMatrix::identity()),
        ("JonesMatrix::horizontal_polarizer", JonesMatrix::horizontal_polarizer()),
        ("JonesMatrix::vertical_polarizer", JonesMatrix::vertical_polarizer()),
    ] {
        m.row(name, Class::A, &[], &jm(j));
    }
    for &t in &ths {
        m.row("JonesMatrix::linear_polarizer", Class::C, &[f(t)], &jm(JonesMatrix::linear_polarizer(black_box(t))));
        m.row("JonesMatrix::half_wave_plate", Class::C, &[f(t)], &jm(JonesMatrix::half_wave_plate(black_box(t))));
        m.row("JonesMatrix::quarter_wave_plate", Class::C, &[f(t)], &jm(JonesMatrix::quarter_wave_plate(black_box(t))));
        m.row("JonesMatrix::rotator", Class::C, &[f(t)], &jm(JonesMatrix::rotator(black_box(t))));
        for &d in &ds {
            m.row("JonesMatrix::waveplate", Class::C, &[f(t), f(d)], &jm(JonesMatrix::waveplate(black_box(t), d)));
        }
    }
    for &j in &ms {
        let [[a, b], [c, d]] = j.m;
        m.row("JonesMatrix::new", Class::A, &jm(j), &jm(JonesMatrix::new(black_box(a), b, c, d)));
        m.row("JonesMatrix::determinant", Class::A, &jm(j), &cx(black_box(j).determinant()));
        for &v in &short {
            m.row("JonesMatrix::apply", Class::A, &[jm(j), jv(v)].concat(), &jv(black_box(j).apply(&v)));
        }
        for &k in ms.iter().step_by(2) {
            m.row("JonesMatrix::compose", Class::A, &[jm(j), jm(k)].concat(), &jm(black_box(j).compose(&k)));
        }
    }

    // Stokes vectors.
    m.row("StokesVector::default", Class::A, &[], &sv(StokesVector::default()));
    for &x in &rs {
        m.row("StokesVector::unpolarized", Class::A, &[f(x)], &sv(StokesVector::unpolarized(black_box(x))));
        m.row("StokesVector::horizontal", Class::A, &[f(x)], &sv(StokesVector::horizontal(black_box(x))));
        m.row("StokesVector::vertical", Class::A, &[f(x)], &sv(StokesVector::vertical(black_box(x))));
        m.row("StokesVector::right_circular", Class::A, &[f(x)], &sv(StokesVector::right_circular(black_box(x))));
        m.row("StokesVector::left_circular", Class::A, &[f(x)], &sv(StokesVector::left_circular(black_box(x))));
    }
    for &v in &vs {
        m.row("StokesVector::from_jones", Class::A, &jv(v), &sv(StokesVector::from_jones(black_box(&v))));
    }
    for &st in &ss {
        let [a, b, c, d] = st.s;
        m.row("StokesVector::new", Class::A, &sv(st), &sv(StokesVector::new(black_box(a), b, c, d)));
        m.row("StokesVector::intensity", Class::A, &sv(st), &[f(black_box(st).intensity())]);
        m.row("StokesVector::degree_of_polarization", Class::C, &sv(st), &res_f(black_box(st).degree_of_polarization()));
        m.row("StokesVector::ellipticity_angle", Class::C, &sv(st), &[f(black_box(st).ellipticity_angle())]);
        m.row("StokesVector::orientation_angle", Class::C, &sv(st), &[f(black_box(st).orientation_angle())]);
    }
    for &a in ss.iter().step_by(3) {
        for &b in ss.iter().step_by(5) {
            m.row("StokesVector::add", Class::A, &[sv(a), sv(b)].concat(), &sv(black_box(a) + b));
        }
    }

    // Mueller matrices.
    for (name, x) in [
        ("MuellerMatrix::identity", MuellerMatrix::identity()),
        ("MuellerMatrix::horizontal_polarizer", MuellerMatrix::horizontal_polarizer()),
        ("MuellerMatrix::vertical_polarizer", MuellerMatrix::vertical_polarizer()),
        ("MuellerMatrix::depolarizer", MuellerMatrix::depolarizer()),
        ("MuellerMatrix::right_circular_polarizer", MuellerMatrix::right_circular_polarizer()),
        ("MuellerMatrix::left_circular_polarizer", MuellerMatrix::left_circular_polarizer()),
    ] {
        m.row(name, Class::A, &[], &mm(x));
    }
    for &x in &rs {
        m.row("MuellerMatrix::attenuator", Class::A, &[f(x)], &mm(MuellerMatrix::attenuator(black_box(x))));
    }
    for &t in &ths {
        m.row("MuellerMatrix::linear_polarizer", Class::C, &[f(t)], &mm(MuellerMatrix::linear_polarizer(black_box(t))));
        m.row("MuellerMatrix::quarter_wave_plate", Class::C, &[f(t)], &mm(MuellerMatrix::quarter_wave_plate(black_box(t))));
        m.row("MuellerMatrix::half_wave_plate", Class::C, &[f(t)], &mm(MuellerMatrix::half_wave_plate(black_box(t))));
        m.row("MuellerMatrix::rotation", Class::C, &[f(t)], &mm(MuellerMatrix::rotation(black_box(t))));
        for &d in &ds {
            m.row("MuellerMatrix::waveplate", Class::C, &[f(t), f(d)], &mm(MuellerMatrix::waveplate(black_box(t), d)));
        }
    }
    for &x in &mus {
        m.row("MuellerMatrix::new", Class::A, &mm(x), &mm(MuellerMatrix::new(black_box(x.m))));
        for &st in ss.iter().step_by(2) {
            m.row("MuellerMatrix::apply", Class::A, &[mm(x), sv(st)].concat(), &sv(black_box(x).apply(&st)));
        }
        for &y in mus.iter().step_by(2) {
            m.row("MuellerMatrix::compose", Class::A, &[mm(x), mm(y)].concat(), &mm(black_box(x).compose(&y)));
        }
    }
    waveplate_edges(m);
    rows.m
}

/// M3 part 2a-1 (ADR 0002 row 33, the user's answer of 2026-11-14): the waveplates at their cancellation edges.
fn waveplate_edges(m: &mut Rows) {
    let mut ths = vec![0.0, -0.0, 5e-324, 1e-200, 1e-8, f64::MAX / 2.0, pow2(1023), 1.5 * pow2(1023)];
    ths.extend([-f64::MAX, f64::MAX]);
    for base in [FRAC_PI_4, 3.0 * FRAC_PI_4, -FRAC_PI_4, 5.0 * FRAC_PI_4, PI / 8.0] {
        let b = base.to_bits();
        for k in -2_i64..=2 {
            ths.push(f64::from_bits(b.checked_add_signed(k).expect("a neighbour of a finite angle")));
        }
    }
    let third = PI / 3.0;
    let mut ds = vec![1e-300, 1e-160, 1e-20, 1e-10, 1e-5, 1e-3, 0.1, 1.0, -1e-8, -0.4, 6.0, 4.0 * PI];
    ds.extend([FRAC_PI_2, PI, -PI]);
    ds.extend([third.next_down(), third, third.next_up(), 2.0 * PI - 1e-8]);
    for &t in &ths {
        let jhw = JonesMatrix::half_wave_plate(black_box(t));
        let jqw = JonesMatrix::quarter_wave_plate(black_box(t));
        let mhw = MuellerMatrix::half_wave_plate(black_box(t));
        let mqw = MuellerMatrix::quarter_wave_plate(black_box(t));
        m.row("JonesMatrix::half_wave_plate", Class::C, &[f(t)], &jm(jhw));
        m.row("JonesMatrix::quarter_wave_plate", Class::C, &[f(t)], &jm(jqw));
        m.row("MuellerMatrix::half_wave_plate", Class::C, &[f(t)], &mm(mhw));
        m.row("MuellerMatrix::quarter_wave_plate", Class::C, &[f(t)], &mm(mqw));
        for &d in &ds {
            m.row("JonesMatrix::waveplate", Class::C, &[f(t), f(d)], &jm(JonesMatrix::waveplate(black_box(t), d)));
            m.row("MuellerMatrix::waveplate", Class::C, &[f(t), f(d)], &mm(MuellerMatrix::waveplate(black_box(t), d)));
        }
    }
}
