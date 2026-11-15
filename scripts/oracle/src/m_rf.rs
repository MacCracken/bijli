//! `rf.tsv` — rust-old/src/rf.rs (M4 part 1a): transmission lines, the reflection coefficient, VSWR and the mismatch
//! losses, the Smith chart, matching, and the antenna array factors.
//!
//! Classes: a function that calls atan2, sin, cos, atan, log10 or acos is class C (`characteristic_impedance_lossy`,
//! `propagation_constant`, `input_impedance_lossless_complex`, `return_loss`, `insertion_loss`, `mismatch_loss`, the
//! two stub matches and the three array factors; `validate.py` writes the correctly rounded readings to `rf.cr.tsv`);
//! the rest are class A (+ − × ÷ and sqrt in Rust's order, through hisab's `Complex` and the module's own Smith
//! division). A complex value is its two parts; an element count (`usize`) is its 64-bit word (usize::MAX is `i:-1`).
//! A `Result` is its value, or `E:<Variant>` and the error's Display text (ADR 0002 row 4).
//!
//! The grids: each argument crossed with the specials, the Rust tests' cases, typical lines (coax, twin lead) from DC
//! to 10 GHz, loads across the Smith chart and its poles, every rejection, each threshold straddled ulp by ulp (the
//! 10⁻³⁰ and 10⁻⁶⁰ divisors, the unbounded input impedance's 10⁻¹⁵ relative test, the stub match's |R_L − Z₀| ≤
//! 10⁻¹²·Z₀, the array factor's |sin(ψ/2)| < 10⁻¹⁵), element counts to usize::MAX, and seeded values over wide ranges.
//! Each value appears once per function: a row's (function, inputs) is its key in the mapping files.
//!
//! M4 part 1b appends the S-parameter matrices: `SMatrix::new`, `two_port`, `get`, `set` (its result, then the
//! matrix after it), `cascade` and `s_to_z`. A matrix is `i:<ports>`, `v:<2n²>` and its entries row-major as (re,
//! im), `v:<n>` and its reference impedances; only consistent shapes, since the port's BjSMatrix is fixed at
//! construction (ADR 0002 row 36). Rows where a second port's Z₀ is ±∞ are mapped by ADR 0002 row 37.
//!
//! M4 part 1c-1 appends `parse_touchstone`: its inputs are the content (`s:`) and the port count, its outputs the
//! data (`i:<ports>`, `f:<z0>`, `v:<points>` and the frequencies, then each matrix as above) or the error. It is class
//! C: MA and DB records take a sine and cosine of the angle, and DB 10^(dB/20).
//!
//! M4 part 1c-2 appends `try_write_touchstone`: its inputs are the data (as `parse_touchstone`'s outputs) and the
//! format (`i:` 0 MagnitudeAngle, 1 RealImag, 2 DbAngle, the port's enum), its output the text (`s:`) or the error.
//! Class C: MA and DB write √(re² + im²), atan2 in degrees and, for DB, log10.
//!
//! M4 part 2a-1 appends the stub match's cancellation edges (ADR 0002 row 39): R_L approaching Z₀ from 10⁻¹ to the
//! 10⁻¹²·Z₀ band on both sides, against reactances from 10⁻¹¹·Z₀ to 10⁶·Z₀ of both signs; the loads where
//! R_L(Z₀ − R_L) = X_L², ulp by ulp and at relative offsets, with their exact zeros; the scaled guard's 2^-480 on
//! either side; and seeded loads in both regimes.

use crate::grids::{Rng, specials};
use crate::m_field::res_f;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::polarization::Complex;
use bijli::rf::{
    SMatrix, TouchstoneData, TouchstoneFormat, characteristic_impedance_lossless, characteristic_impedance_lossy,
    circular_array_factor, denormalize_impedance, gamma_from_z, input_impedance_lossless_complex, insertion_loss,
    linear_array_factor, mismatch_loss, normalize_impedance, parse_touchstone, phase_velocity, planar_array_factor,
    propagation_constant, quarter_wave_transformer, reflection_coefficient, reflection_coefficient_complex,
    return_loss, s_to_z, single_stub_match, single_stub_match_complex, try_write_touchstone, vswr, z_from_gamma,
};
use std::collections::BTreeSet;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};
use std::hint::black_box;

fn cx(z: Complex) -> Vec<Tok> {
    vec![f(z.re), f(z.im)]
}

fn res<T>(r: bijli::error::Result<T>, ok: impl FnOnce(&T) -> Vec<Tok>) -> Vec<Tok> {
    match r {
        Ok(v) => ok(&v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// A `usize` as its 64-bit word.
fn word(n: usize) -> i64 {
    i64::from_ne_bytes(n.to_ne_bytes())
}

/// The bits of every f64 and the value of every integer in a token list, for keeping each row's inputs once.
fn key(name: &str, toks: &[Tok]) -> String {
    let mut k = String::from(name);
    for t in toks {
        match t {
            Tok::F(x) => k.push_str(&format!(" {:016x}", x.to_bits())),
            Tok::I(n) => k.push_str(&format!(" i{n}")),
            Tok::S(x) => k.push_str(&format!(" s{x:?}")),
            _ => {}
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
    fn row(&mut self, name: &'static str, class: Class, ins: &[Tok], outs: impl FnOnce() -> Vec<Tok>) {
        if self.seen.insert(key(name, ins)) {
            let o = outs();
            self.m.row(name, class, ins, &o);
        }
    }

    fn lossless(&mut self, l: f64, c: f64) {
        let ins = [f(l), f(c)];
        self.row("rf::characteristic_impedance_lossless", Class::A, &ins, || {
            res_f(characteristic_impedance_lossless(black_box(l), black_box(c)))
        });
        self.row("rf::phase_velocity", Class::A, &ins, || res_f(phase_velocity(black_box(l), black_box(c))));
    }

    fn lossy(&mut self, v: [f64; 5]) {
        let ins: Vec<Tok> = v.iter().map(|&x| f(x)).collect();
        self.row("rf::characteristic_impedance_lossy", Class::C, &ins, || {
            cx(characteristic_impedance_lossy(black_box(v[0]), v[1], v[2], v[3], black_box(v[4])))
        });
    }

    fn propagation(&mut self, v: [f64; 5]) {
        let ins: Vec<Tok> = v.iter().map(|&x| f(x)).collect();
        self.row("rf::propagation_constant", Class::C, &ins, || {
            cx(propagation_constant(black_box(v[0]), v[1], v[2], v[3], black_box(v[4])))
        });
    }

    fn zin(&mut self, z0: f64, zl: (f64, f64), beta: f64, length: f64) {
        let ins = [f(z0), f(zl.0), f(zl.1), f(beta), f(length)];
        self.row("rf::input_impedance_lossless_complex", Class::C, &ins, || {
            res(
                input_impedance_lossless_complex(
                    black_box(z0),
                    black_box(Complex::new(zl.0, zl.1)),
                    black_box(beta),
                    black_box(length),
                ),
                |z| cx(*z),
            )
        });
    }

    fn gamma_real(&mut self, zl: f64, z0: f64) {
        let ins = [f(zl), f(z0)];
        self.row("rf::reflection_coefficient", Class::A, &ins, || {
            res_f(reflection_coefficient(black_box(zl), black_box(z0)))
        });
        self.row("rf::quarter_wave_transformer", Class::A, &[f(z0), f(zl)], || {
            res_f(quarter_wave_transformer(black_box(z0), black_box(zl)))
        });
    }

    fn gamma_complex(&mut self, zl: (f64, f64), z0: (f64, f64)) {
        let ins = [f(zl.0), f(zl.1), f(z0.0), f(z0.1)];
        self.row("rf::reflection_coefficient_complex", Class::A, &ins, || {
            let (a, b) = (black_box(Complex::new(zl.0, zl.1)), black_box(Complex::new(z0.0, z0.1)));
            res(reflection_coefficient_complex(a, b), |z| cx(*z))
        });
    }

    fn losses(&mut self, g: f64) {
        let ins = [f(g)];
        self.row("rf::vswr", Class::A, &ins, || res_f(vswr(black_box(g))));
        self.row("rf::return_loss", Class::C, &ins, || res_f(return_loss(black_box(g))));
        self.row("rf::insertion_loss", Class::C, &ins, || res_f(insertion_loss(black_box(g))));
        self.row("rf::mismatch_loss", Class::C, &ins, || res_f(mismatch_loss(black_box(g))));
    }

    fn norm(&mut self, z: (f64, f64), z0: f64) {
        let ins = [f(z.0), f(z.1), f(z0)];
        self.row("rf::normalize_impedance", Class::A, &ins, || {
            res(normalize_impedance(black_box(Complex::new(z.0, z.1)), black_box(z0)), |w| cx(*w))
        });
        self.row("rf::denormalize_impedance", Class::A, &ins, || {
            cx(denormalize_impedance(black_box(Complex::new(z.0, z.1)), black_box(z0)))
        });
    }

    fn smith(&mut self, z: (f64, f64)) {
        let ins = [f(z.0), f(z.1)];
        self.row("rf::gamma_from_z", Class::A, &ins, || cx(gamma_from_z(black_box(Complex::new(z.0, z.1)))));
        self.row("rf::z_from_gamma", Class::A, &ins, || cx(z_from_gamma(black_box(Complex::new(z.0, z.1)))));
    }

    fn stub(&mut self, zl: f64, z0: f64) {
        let ins = [f(zl), f(z0)];
        self.row("rf::single_stub_match", Class::C, &ins, || {
            res(single_stub_match(black_box(zl), black_box(z0)), |&(d, l)| vec![f(d), f(l)])
        });
    }

    fn stub_complex(&mut self, zl: (f64, f64), z0: f64) {
        let ins = [f(zl.0), f(zl.1), f(z0)];
        self.row("rf::single_stub_match_complex", Class::C, &ins, || {
            res(single_stub_match_complex(black_box(Complex::new(zl.0, zl.1)), black_box(z0)), |s| {
                vec![f(s[0].0), f(s[0].1), f(s[1].0), f(s[1].1)]
            })
        });
    }

    fn linear(&mut self, n: usize, d: f64, k: f64, beta: f64, theta: f64) {
        let ins = [i(word(n)), f(d), f(k), f(beta), f(theta)];
        self.row("rf::linear_array_factor", Class::C, &ins, || {
            vec![f(linear_array_factor(black_box(n), black_box(d), black_box(k), black_box(beta), black_box(theta)))]
        });
    }

    fn planar(&mut self, n: (usize, usize), v: [f64; 7]) {
        let mut ins = vec![i(word(n.0)), i(word(n.1))];
        ins.extend(v.iter().map(|&x| f(x)));
        self.row("rf::planar_array_factor", Class::C, &ins, || {
            vec![f(planar_array_factor(
                black_box(n.0),
                black_box(n.1),
                black_box(v[0]),
                v[1],
                v[2],
                v[3],
                v[4],
                black_box(v[5]),
                black_box(v[6]),
            ))]
        });
    }

    fn circular(&mut self, n: usize, v: [f64; 5]) {
        let mut ins = vec![i(word(n))];
        ins.extend(v.iter().map(|&x| f(x)));
        self.row("rf::circular_array_factor", Class::C, &ins, || {
            vec![f(circular_array_factor(black_box(n), black_box(v[0]), v[1], v[2], v[3], black_box(v[4])))]
        });
    }
}

fn up(x: f64, k: u64) -> f64 {
    let mut y = x;
    for _ in 0..k {
        y = y.next_up();
    }
    y
}

fn down(x: f64, k: u64) -> f64 {
    let mut y = x;
    for _ in 0..k {
        y = y.next_down();
    }
    y
}

/// x and its neighbours k ulps away on both sides, k = 1…n.
fn around(x: f64, n: u64) -> Vec<f64> {
    let mut v = vec![x];
    for k in 1..=n {
        v.push(up(x, k));
        v.push(down(x, k));
    }
    v
}

/// Per-unit-length values and impedances: typical lines and the range ends where a product or quotient over- or
/// underflows.
fn line_vals() -> Vec<f64> {
    let mut v = specials();
    v.extend([250e-9, 100e-12, 4e-7, 1e-10, 2.0, 0.5, 1e-300, 1e300, 1e-160, 1e160, 5e-310, -250e-9]);
    v
}

/// Impedances for the real reflection coefficient and the quarter-wave transformer.
fn imp_vals() -> Vec<f64> {
    let mut v = specials();
    v.extend([50.0, 75.0, 100.0, 25.0, 1e15, -50.0, 1e-30, -1e-30, 5e-31, 1e300, 1e-300, 1e-160, 1e160]);
    v
}

/// Complex impedances and reflection coefficients: across the Smith chart, its poles, the specials in either part.
fn cx_vals() -> Vec<(f64, f64)> {
    let parts = [0.0, -0.0, 1.0, -1.0, 2.0, 0.5, 1e300, -1e300, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 1e-300];
    let mut v = Vec::new();
    for &a in &parts {
        for &b in &parts {
            v.push((a, b));
        }
    }
    v.extend([(5e-324, 0.0), (f64::MAX, f64::MAX), (-1.0, 1e-300), (1.0, -1e-300), (-1.0, 5e-324), (0.3, 0.4)]);
    v.extend([(-0.6, 0.8), (0.0, 1.0), (0.0, -1.0), (3.0, 4.0), (1e154, 1e154), (1e-160, 1e-160)]);
    v
}

/// Lines for the lossy impedance and the propagation constant: (R′, L′, G′, C′, ω).
fn line_bases() -> Vec<[f64; 5]> {
    vec![
        [0.001, 250e-9, 1e-6, 100e-12, 2e9 * PI],
        [0.0, 250e-9, 0.0, 100e-12, 2e9 * PI],
        [5.0, 400e-9, 0.01, 67e-12, 2e6 * PI],
        [1.0, 1e-6, 1e-3, 1e-11, 1e3],
        [0.0, 250e-9, 0.0, 100e-12, 0.0],
        [1.0, 250e-9, 0.0, 100e-12, 0.0],
        [4.0, 250e-9, 1.0, 100e-12, 0.0],
    ]
}

pub fn emit() -> Module {
    let mut rows = Rows { m: Module::new("rf"), seen: BTreeSet::new() };
    let r = &mut rows;
    let mut rng = Rng::new(0x00f1_e1d5);

    // ── transmission lines ─────────────────────────────────────────────────────────────
    let lv = line_vals();
    for &l in &lv {
        for &c in &lv {
            r.lossless(l, c);
        }
    }
    for base in line_bases() {
        r.lossy(base);
        r.propagation(base);
        for slot in 0..5 {
            for &x in &specials() {
                let mut v = base;
                v[slot] = x;
                r.lossy(v);
                r.propagation(v);
            }
        }
    }
    // The lossy impedance's DC branch: ω = ±0 and G′ + jωC′ = 0, with every sign of R′, L′ and C′.
    for &rr in &[0.0, -0.0, 1.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for &l in &[250e-9, -250e-9, 0.0, f64::NAN, f64::INFINITY] {
            for &g in &[0.0, -0.0, 1.0, f64::NAN] {
                for &c in &[100e-12, 0.0, -0.0, -1e-10, f64::INFINITY, f64::NAN] {
                    for &w in &[0.0, -0.0] {
                        r.lossy([rr, l, g, c, w]);
                    }
                }
            }
        }
    }
    for _ in 0..200 {
        let v = [
            rng.log_uniform(-6.0, 3.0),
            rng.log_uniform(-9.0, -5.0),
            rng.log_uniform(-9.0, -1.0),
            rng.log_uniform(-13.0, -9.0),
            rng.log_uniform(0.0, 12.0),
        ];
        r.lossless(v[1], v[3]);
        r.lossy(v);
        r.propagation(v);
    }
    for _ in 0..100 {
        let v = [
            rng.signed_log_uniform(-300.0, 300.0),
            rng.signed_log_uniform(-300.0, 300.0),
            rng.signed_log_uniform(-300.0, 300.0),
            rng.signed_log_uniform(-300.0, 300.0),
            rng.log_uniform(-300.0, 300.0),
        ];
        r.lossless(v[1].abs(), v[3].abs());
        r.lossy(v);
        r.propagation(v);
    }

    // ── input impedance ────────────────────────────────────────────────────────────────
    // The Rust tests' cases.
    r.zin(50.0, (100.0, 0.0), 1.0, FRAC_PI_4);
    r.zin(50.0, (30.0, 40.0), TAU, 0.3);
    r.zin(50.0, (0.0, 0.0), TAU, 0.25);
    r.zin(0.0, (50.0, 0.0), 1.0, 1.0);
    r.zin(f64::NAN, (50.0, 0.0), 1.0, 1.0);
    r.zin(50.0, (f64::NAN, 0.0), 1.0, 1.0);
    r.zin(50.0, (50.0, 0.0), f64::INFINITY, 1.0);
    r.zin(50.0, (100.0, -25.0), TAU, 0.125);
    let loads = [
        (0.0, 0.0),
        (50.0, 0.0),
        (100.0, 0.0),
        (25.0, 0.0),
        (30.0, 40.0),
        (60.0, -80.0),
        (1e-300, 0.0),
        (1e300, 1e300),
        (-50.0, 0.0),
        (0.0, 50.0),
        (0.0, -50.0),
        (1e15, 0.0),
        (-0.0, -0.0),
    ];
    for &zl in &loads {
        for &len in &[0.0, -0.0, 0.125, 0.25, 0.5, 0.3, 1.0, 0.75, 100.25, 1e10, -0.125] {
            r.zin(50.0, zl, TAU, len);
        }
        r.zin(75.0, zl, 1.0, FRAC_PI_2);
        r.zin(1e300, zl, TAU, 0.125);
        r.zin(1e-300, zl, TAU, 0.125);
    }
    let base_zin = [50.0, 30.0, 40.0, TAU, 0.3];
    for slot in 0..5 {
        for &x in &specials() {
            let mut v = base_zin;
            v[slot] = x;
            r.zin(v[0], (v[1], v[2]), v[3], v[4]);
        }
    }
    // The unbounded test |den| ≤ 10⁻¹⁵·(Z₀ + |Z_L|): a short and tiny loads near βl = π/2, ulp by ulp.
    for &len in &around(0.25, 8) {
        r.zin(50.0, (0.0, 0.0), TAU, len);
        r.zin(50.0, (1e-14, 0.0), TAU, len);
    }
    for &len in &around(FRAC_PI_2, 8) {
        r.zin(50.0, (0.0, 0.0), 1.0, len);
        r.zin(1.0, (0.0, 1e-16), 1.0, len);
    }
    for _ in 0..200 {
        let z0 = rng.log_uniform(0.0, 3.0);
        let zl = (rng.log_uniform(-2.0, 4.0), rng.signed_log_uniform(-2.0, 4.0));
        r.zin(z0, zl, TAU, rng.uniform(0.0, 2.0));
    }
    for _ in 0..60 {
        let z0 = rng.log_uniform(-300.0, 300.0);
        let zl = (rng.signed_log_uniform(-300.0, 300.0), rng.signed_log_uniform(-300.0, 300.0));
        r.zin(z0, zl, rng.signed_log_uniform(-3.0, 3.0), rng.signed_log_uniform(-3.0, 3.0));
    }

    // ── reflection, VSWR and the losses ────────────────────────────────────────────────
    let iv = imp_vals();
    for &zl in &iv {
        for &z0 in &iv {
            r.gamma_real(zl, z0);
        }
    }
    // |Z_L + Z₀| < 10⁻³⁰: the sum exact at 10⁻³⁰ and its neighbours, and rounded (2·10⁻³⁰ − 10⁻³⁰).
    for &x in &around(1e-30, 4) {
        r.gamma_real(x, 0.0);
        r.gamma_real(-x, 0.0);
        r.gamma_real(2.0 * x, -x);
        r.gamma_real(x, 1e-46);
    }
    let cv = [
        (0.0, 0.0),
        (50.0, 0.0),
        (75.0, 25.0),
        (30.0, -40.0),
        (-50.0, 0.0),
        (0.0, 50.0),
        (1e-30, 0.0),
        (1e300, 1e300),
        (f64::NAN, 0.0),
        (0.0, f64::NAN),
        (f64::INFINITY, 0.0),
        (0.0, f64::NEG_INFINITY),
        (-0.0, -0.0),
        (1e-160, 1e-160),
        (50.0, -0.0),
    ];
    for &zl in &cv {
        for &z0 in &cv {
            r.gamma_complex(zl, z0);
        }
    }
    // |Z_L + Z₀|² < 10⁻⁶⁰: a real sum at 10⁻³⁰ ± ulps, and an equal-part one at 10⁻³⁰/√2.
    for &x in &around(1e-30, 4) {
        r.gamma_complex((x, 0.0), (0.0, 0.0));
        r.gamma_complex((0.0, x), (0.0, -0.0));
    }
    for &x in &around(7.071_067_811_865_475e-31, 4) {
        r.gamma_complex((x, x), (0.0, 0.0));
    }
    for _ in 0..100 {
        let zl = (rng.signed_log_uniform(-3.0, 4.0), rng.signed_log_uniform(-3.0, 4.0));
        let z0 = (rng.log_uniform(0.0, 3.0), rng.signed_log_uniform(-3.0, 1.0));
        r.gamma_complex(zl, z0);
    }
    let mut gv = specials();
    gv.extend([0.1, 1.0 / 3.0, 0.333, 0.5, 0.9, 0.99, 0.999_999, 1e-300, 1e-8, 1.5, -0.1, 2.0, 5e-324, 1e-160]);
    gv.extend(around(1.0, 3));
    gv.extend(around(0.0, 2));
    for &g in &gv {
        r.losses(g);
    }
    for _ in 0..60 {
        let g = rng.unit();
        r.losses(g);
    }

    // ── the Smith chart ────────────────────────────────────────────────────────────────
    let mut z0s = specials();
    z0s.extend([50.0, 75.0, -50.0, 1e-30, 9.999_999_999_999_999e-31, 1e300, 1e-300, 3.0]);
    for &z in &cv {
        for &z0 in &z0s {
            r.norm(z, z0);
        }
    }
    for &z in &cx_vals() {
        r.smith(z);
    }
    for _ in 0..100 {
        let z = (rng.signed_log_uniform(-3.0, 3.0), rng.signed_log_uniform(-3.0, 3.0));
        r.smith(z);
        r.norm(z, rng.log_uniform(0.0, 3.0));
    }
    // The poles' neighbourhoods: z = −1 and Γ = 1 one ulp off in either part.
    for &x in &around(-1.0, 2) {
        for &y in &[0.0, -0.0, 5e-324, -5e-324, 1e-300] {
            r.smith((x, y));
            r.smith((-x, y));
        }
    }

    // ── matching ───────────────────────────────────────────────────────────────────────
    let mut sv = specials();
    sv.extend([50.0, 100.0, 25.0, 75.0, 1e-300, 1e300, 1e-160, 1e160, -100.0]);
    for &zl in &sv {
        for &z0 in &sv {
            r.stub(zl, z0);
        }
    }
    // The matched test |R_L − Z₀| ≤ 10⁻¹²·Z₀ on both sides.
    for &zl in &around(50.0 * (1.0 + 1e-12), 4) {
        r.stub(zl, 50.0);
    }
    for &zl in &around(50.0 * (1.0 - 1e-12), 4) {
        r.stub(zl, 50.0);
    }
    for &zl in &around(50.0, 3) {
        r.stub(zl, 50.0);
    }
    let zls = [
        (100.0, 0.0),
        (25.0, 0.0),
        (60.0, -80.0),
        (50.0, 50.0),
        (10.0, 120.0),
        (50.0, 0.0),
        (50.0, -30.0),
        (0.0, 30.0),
        (-10.0, 5.0),
        (1e-300, 1e300),
        (1e300, -1e300),
        (f64::NAN, 0.0),
        (50.0, f64::INFINITY),
        (5e-324, 0.0),
        (-0.0, 1.0),
    ];
    for &zl in &zls {
        for &z0 in &[50.0, 75.0, 1.0, 0.0, -50.0, f64::NAN, f64::INFINITY, 1e300, 1e-300] {
            r.stub_complex(zl, z0);
        }
    }
    for &x in &around(50.0 * (1.0 + 1e-12), 3) {
        r.stub_complex((x, 20.0), 50.0);
        r.stub_complex((x, 0.0), 50.0);
    }
    for _ in 0..150 {
        let z0 = rng.log_uniform(0.0, 3.0);
        let zl = (rng.log_uniform(-1.0, 4.0), rng.signed_log_uniform(-1.0, 4.0));
        r.stub_complex(zl, z0);
        r.stub(zl.0, z0);
    }

    // ── the array factors ──────────────────────────────────────────────────────────────
    // Past 2^53 the count rounds as it converts: above 2^63 a tie to even (…C00), a sticky bit past the tie (…C01).
    let counts = [
        0,
        1,
        2,
        3,
        4,
        8,
        16,
        100,
        1000,
        (1 << 53) + 1,
        (1 << 63) + 1,
        (1 << 63) + 0xc00,
        (1 << 63) + 0xc01,
        (1 << 63) + 0x400,
        (1 << 63) + 0x401,
        usize::MAX - 0x3ff,
        usize::MAX,
    ];
    for &n in &counts {
        for &d in &[0.5, 0.25, 1.0] {
            for &beta in &[0.0, FRAC_PI_2, -PI] {
                for &theta in &[0.0, PI / 6.0, FRAC_PI_4, FRAC_PI_2, 0.7, PI, 2.0] {
                    r.linear(n, d, TAU, beta, theta);
                }
            }
        }
    }
    r.linear(4, 0.5, TAU, 0.0, FRAC_PI_2);
    r.linear(1, 0.5, 10.0, 0.0, 0.5);
    r.linear(0, 0.5, 10.0, 0.0, 0.5);
    for k in 0..100_u32 {
        r.linear(8, 0.5, TAU, 0.0, f64::from(k) * PI / 100.0);
    }
    // |sin(ψ/2)| < 10⁻¹⁵ around the main beam and the grating lobe ψ/2 = π, ulp by ulp.
    for j in -6..=6_i32 {
        r.linear(8, 0.5, TAU, 2e-15 * (1.0 + f64::from(j) / 20.0), FRAC_PI_2);
    }
    for &b in &around(TAU, 6) {
        r.linear(8, 0.5, TAU, b, FRAC_PI_2);
    }
    for slot in 0..4 {
        for &x in &specials() {
            let mut v = [0.5, TAU, 0.0, 0.7];
            v[slot] = x;
            r.linear(8, v[0], v[1], v[2], v[3]);
        }
    }
    for _ in 0..100 {
        let n = usize::try_from(rng.int(2, 64)).unwrap_or(2);
        r.linear(n, rng.uniform(0.1, 2.0), TAU, rng.uniform(-PI, PI), rng.uniform(0.0, PI));
    }
    r.planar((4, 4), [0.5, 0.5, TAU, 0.0, 0.0, 0.3, 0.5]);
    for &(nx, ny) in &[(4, 4), (1, 1), (0, 3), (8, 2), (16, 16)] {
        for &theta in &[0.0, 0.3, FRAC_PI_4, FRAC_PI_2, PI, -0.5] {
            for &phi in &[0.0, 0.5, FRAC_PI_2, PI, -1.0] {
                r.planar((nx, ny), [0.5, 0.5, TAU, 0.0, 0.0, theta, phi]);
                r.planar((nx, ny), [0.5, 0.25, TAU, FRAC_PI_4, -FRAC_PI_2, theta, phi]);
            }
        }
    }
    for slot in 0..7 {
        for &x in &specials() {
            let mut v = [0.5, 0.5, TAU, 0.0, 0.0, 0.3, 0.5];
            v[slot] = x;
            r.planar((4, 4), v);
        }
    }
    for _ in 0..60 {
        let n = (usize::try_from(rng.int(1, 16)).unwrap_or(1), usize::try_from(rng.int(1, 16)).unwrap_or(1));
        let v = [
            rng.uniform(0.1, 1.0),
            rng.uniform(0.1, 1.0),
            TAU,
            rng.uniform(-PI, PI),
            rng.uniform(-PI, PI),
            rng.uniform(0.0, PI),
            rng.uniform(-PI, PI),
        ];
        r.planar(n, v);
    }
    r.circular(1, [0.5, 10.0, 0.5, 0.0, 0.0]);
    r.circular(8, [0.25, TAU, 0.0, 0.0, 0.0]);
    for &n in &[0, 1, 2, 3, 8, 16, 64] {
        for &a in &[0.25, 0.5, 1.0] {
            for &theta in &[0.0, PI / 6.0, FRAC_PI_2, 1.0] {
                for &phi in &[0.0, FRAC_PI_4, PI] {
                    for &beta in &[0.0, FRAC_PI_4] {
                        r.circular(n, [a, TAU, theta, phi, beta]);
                    }
                }
            }
        }
    }
    r.circular(1000, [1.0, TAU, 1.0, 0.5, 0.01]);
    for slot in 0..5 {
        for &x in &specials() {
            let mut v = [0.5, TAU, 1.0, 0.5, 0.25];
            v[slot] = x;
            r.circular(8, v);
        }
    }
    for _ in 0..60 {
        let n = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        let v = [rng.uniform(0.1, 2.0), TAU, rng.uniform(0.0, PI), rng.uniform(-PI, PI), rng.uniform(-PI, PI)];
        r.circular(n, v);
    }

    emit_1b(r, &mut rng);
    emit_1c1(r, &mut rng);
    emit_1c2(r, &mut rng);
    emit_2a1(r, &mut rng);

    rows.m
}

/// M4 part 2a-1: the loads where the stub match's root (X_L ∓ disc)/(R_L − Z₀) cancels, so ADR 0002 row 39 maps
/// them: R_L → Z₀ (disc → |X_L|) and R_L(Z₀ − R_L) → X_L² (the root → 0).
fn emit_2a1(r: &mut Rows, rng: &mut Rng) {
    let z0s = [50.0, 75.0, 1.0, 1e-3, 1e6, 3.0 * 2f64.powi(-400), 2f64.powi(400), 1e150, 1e-150];
    let rels = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 3e-12, 1.5e-12];
    let xs = [1e-11, 1e-9, 1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6];
    for &z0 in &z0s {
        for &rel in &rels {
            for sr in [1.0, -1.0] {
                let rl = z0 * (1.0 + sr * rel);
                for &x in &xs {
                    r.stub_complex((rl, z0 * x), z0);
                    r.stub_complex((rl, -z0 * x), z0);
                }
            }
        }
    }
    // R_L(Z₀ − R_L) = X_L²: exact zeros, then ulp neighbours and relative offsets of the root's X_L.
    let zeros = [(25.0, 25.0, 50.0), (10.0, 20.0, 50.0), (40.0, 20.0, 50.0), (0.5, 0.5, 1.0), (18.0, 24.0, 50.0)];
    for &(rl, x, z0) in &zeros {
        for &xv in &around(x, 4) {
            r.stub_complex((rl, xv), z0);
            r.stub_complex((rl, -xv), z0);
        }
    }
    for &z0 in &[50.0_f64, 1.0, 1e5, 1e-4] {
        for &frac in &[0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99] {
            let rl = z0 * frac;
            let x0 = (rl * (z0 - rl)).sqrt();
            for &xv in &around(x0, 2) {
                r.stub_complex((rl, xv), z0);
                r.stub_complex((rl, -xv), z0);
            }
            for k in [-1.0, -3.0, -6.0, -9.0, -12.0, -15.0] {
                for s in [1.0, -1.0] {
                    let xv = x0 * (1.0 + s * 10f64.powf(k));
                    r.stub_complex((rl, xv), z0);
                    r.stub_complex((rl, -xv), z0);
                }
            }
        }
    }
    // The scaled guard: R_L/Z₀ at 2^-479 (the cancellation-free root) and 2^-481 (Rust's chain), X_L = √(R_L·Z₀).
    for e in [-479, -480, -481] {
        let rl = 2f64.powi(e);
        r.stub_complex((rl, rl.sqrt()), 1.0);
        r.stub_complex((rl, -rl.sqrt()), 1.0);
        r.stub_complex((rl * 50.0, (rl * 50.0 * 50.0).sqrt()), 50.0);
    }
    for _ in 0..400 {
        let z0 = rng.log_uniform(-3.0, 6.0);
        let rl = z0 * (1.0 + rng.signed_log_uniform(-11.9, -0.3));
        let x = z0 * rng.signed_log_uniform(-6.0, 4.0);
        r.stub_complex((rl, x), z0);
    }
    for _ in 0..200 {
        let z0 = rng.log_uniform(-3.0, 6.0);
        let rl = z0 * rng.uniform(0.001, 0.999);
        let x = (rl * (z0 - rl)).sqrt() * (1.0 + rng.signed_log_uniform(-15.0, -1.0));
        r.stub_complex((rl, x), z0);
        r.stub_complex((rl, -x), z0);
    }
}

// ── M4 part 1b: the S-parameter matrices ─────────────────────────────────────────────────────────────────────────

/// A matrix as tokens: `i:<ports>`, `v:<2n²>` and the entries row-major (re, im), `v:<n>` and the reference
/// impedances. Only consistent shapes appear: the port's BjSMatrix is fixed at construction (ADR 0002 row 36).
fn sm_toks(m: &SMatrix) -> Vec<Tok> {
    let mut v = vec![i(word(m.ports)), Tok::V(2 * m.ports * m.ports)];
    for row in &m.s {
        for c in row {
            v.push(f(c.re));
            v.push(f(c.im));
        }
    }
    v.push(Tok::V(m.z0.len()));
    v.extend(m.z0.iter().map(|&x| f(x)));
    v
}

/// An n-port from its entries, row-major as (re, im), and its reference impedances.
fn sm_of(ports: usize, s: &[(f64, f64)], z0: &[f64]) -> SMatrix {
    assert_eq!(s.len(), ports * ports, "a consistent matrix");
    assert_eq!(z0.len(), ports, "one Z₀ per port");
    let row = |r: usize| (0..ports).map(|c| Complex::new(s[r * ports + c].0, s[r * ports + c].1)).collect();
    SMatrix { ports, s: (0..ports).map(row).collect(), z0: z0.to_vec() }
}

/// A 2-port [S11, S12, S21, S22] with reference impedances z0.
fn tp(s: [(f64, f64); 4], z0: [f64; 2]) -> SMatrix {
    sm_of(2, &s, &z0)
}

impl Rows {
    fn sm_new(&mut self, ports: usize, z0: f64) {
        let ins = [i(word(ports)), f(z0)];
        self.row("SMatrix::new", Class::A, &ins, || sm_toks(&SMatrix::new(black_box(ports), black_box(z0))));
    }

    fn sm_two_port(&mut self, s: [(f64, f64); 4], z0: f64) {
        let mut ins: Vec<Tok> = s.iter().flat_map(|&(a, b)| [f(a), f(b)]).collect();
        ins.push(f(z0));
        self.row("SMatrix::two_port", Class::A, &ins, || {
            let c = |k: usize| black_box(Complex::new(s[k].0, s[k].1));
            sm_toks(&SMatrix::two_port(c(0), c(1), c(2), c(3), black_box(z0)))
        });
    }

    fn sm_get(&mut self, m: &SMatrix, a: usize, b: usize) {
        let mut ins = sm_toks(m);
        ins.extend([i(word(a)), i(word(b))]);
        self.row("SMatrix::get", Class::A, &ins, || res(black_box(m).get(black_box(a), black_box(b)), |z| cx(*z)));
    }

    fn sm_set(&mut self, m: &SMatrix, a: usize, b: usize, v: (f64, f64)) {
        let mut ins = sm_toks(m);
        ins.extend([i(word(a)), i(word(b)), f(v.0), f(v.1)]);
        self.row("SMatrix::set", Class::A, &ins, || {
            let mut w = m.clone();
            let mut out = res(w.set(black_box(a), black_box(b), black_box(Complex::new(v.0, v.1))), |_| vec![i(0)]);
            out.extend(sm_toks(&w));
            out
        });
    }

    fn cascade(&mut self, a: &SMatrix, b: &SMatrix) {
        let mut ins = sm_toks(a);
        ins.extend(sm_toks(b));
        self.row("SMatrix::cascade", Class::A, &ins, || res(black_box(a).cascade(black_box(b)), sm_toks));
    }

    fn s_to_z(&mut self, m: &SMatrix) {
        let ins = sm_toks(m);
        self.row("rf::s_to_z", Class::A, &ins, || {
            res(s_to_z(black_box(m)), |z| {
                let mut v = Vec::new();
                for row in z {
                    for c in row {
                        v.extend([f(c.re), f(c.im)]);
                    }
                }
                v
            })
        });
    }

    /// cascade(a, b), cascade(b, a) and s_to_z of each: the grids feed the three together.
    fn two_ports(&mut self, a: &SMatrix, b: &SMatrix) {
        self.cascade(a, b);
        self.cascade(b, a);
        self.s_to_z(a);
        self.s_to_z(b);
    }
}

/// Part 1b's grids: `new` and `two_port` over port counts and impedances; `get` and `set` on every index from 0 to
/// one past the last and usize::MAX; the cascade and S → Z of the Rust tests' networks, of non-2-ports, of every
/// reference impedance rule (equal, mixed, within and just outside z0_equal's 10⁻¹², zero, negative, NaN, ±∞ on
/// either port), of S21's, T22's and det's 10⁻⁶⁰ tests straddled ulp by ulp, of the specials in each entry, and of
/// seeded networks over wide ranges.
fn emit_1b(r: &mut Rows, rng: &mut Rng) {
    for &n in &[0, 1, 2, 3, 4, 8] {
        for &z0 in &[50.0, 75.0, 0.0, -50.0, f64::NAN, f64::INFINITY, 5e-324, 1e300] {
            r.sm_new(n, z0);
        }
    }
    let star_a = [(0.1, 0.2), (0.3, -0.1), (0.8, 0.05), (-0.2, 0.1)];
    let star_b = [(-0.15, 0.05), (0.6, 0.2), (0.7, -0.3), (0.25, -0.05)];
    let thru = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 0.0)];
    for &s in &[star_a, star_b, thru, [(0.1, 0.0), (0.9, 0.0), (0.9, 0.0), (0.1, 0.0)]] {
        for &z0 in &[50.0, f64::NAN, -0.0] {
            r.sm_two_port(s, z0);
        }
    }
    for &x in &specials() {
        r.sm_two_port([(x, -x), (0.5, x), (x, 0.25), (-x, x)], 50.0);
    }

    // get and set on every index from 0 to one past the last, and usize::MAX and 2^63.
    let mut three = Vec::new();
    for k in 0..9_u32 {
        three.push((f64::from(k) / 10.0, -f64::from(k) / 100.0));
    }
    let mats = [
        SMatrix::new(0, 50.0),
        sm_of(1, &[(0.5, 0.1)], &[75.0]),
        tp(star_a, [50.0, 50.0]),
        sm_of(3, &three, &[50.0, 50.0, 50.0]),
    ];
    let idx = [0, 1, 2, 3, 4, usize::MAX, 1 << 63];
    for m in &mats {
        for &a in &idx {
            for &b in &idx {
                r.sm_get(m, a, b);
                r.sm_set(m, a, b, (-0.5, 0.125));
            }
        }
    }
    r.sm_set(&mats[2], 1, 2, (f64::NAN, f64::INFINITY));
    r.sm_set(&mats[2], 2, 1, (-0.0, -0.0));

    // The Rust tests' networks: thru, the star product's pair, the attenuator, the shunt resistor.
    let z50 = [50.0, 50.0];
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let atten = [(0.0, 0.0), (h, 0.0), (h, 0.0), (0.0, 0.0)];
    let shunt = [(-0.2, 0.0), (0.8, 0.0), (0.8, 0.0), (-0.2, 0.0)];
    let bench = [(0.1, 0.0), (0.9, 0.0), (0.9, 0.0), (0.1, 0.0)];
    r.two_ports(&tp(thru, z50), &tp(thru, z50));
    r.two_ports(&tp(star_a, z50), &tp(star_b, z50));
    r.two_ports(&tp(atten, z50), &tp(shunt, z50));
    r.two_ports(&tp(bench, z50), &tp(bench, z50));
    r.two_ports(&tp(atten, z50), &tp(thru, z50));
    // Non-2-ports on either side.
    for m in [&mats[0], &mats[1], &mats[3], &SMatrix::new(3, 50.0), &SMatrix::new(4, 50.0)] {
        r.two_ports(m, &tp(thru, z50));
        r.cascade(m, m);
    }
    // The reference impedance rules: equal, mixed, within and just past z0_equal's 10⁻¹², invalid, ±∞ on the second
    // port (Rust's ∞ ≤ 10⁻¹²·∞ accepts it: ADR 0002 row 37).
    let lim = 50.0 * (1.0 + 1e-12);
    let mut z0s = vec![
        [50.0, 50.0],
        [50.0, 75.0],
        [75.0, 75.0],
        [75.0, 50.0],
        [0.0, 0.0],
        [-50.0, -50.0],
        [f64::NAN, 50.0],
        [50.0, f64::NAN],
        [f64::INFINITY, 50.0],
        [50.0, f64::INFINITY],
        [50.0, f64::NEG_INFINITY],
        [1e300, f64::INFINITY],
        [5e-324, 5e-324],
        [1e300, 1e300],
        [50.0, -50.0],
        [-0.0, 50.0],
    ];
    for &x in &around(lim, 3) {
        z0s.push([50.0, x]);
        z0s.push([x, 50.0]);
    }
    for z in &z0s {
        r.two_ports(&tp(star_a, *z), &tp(star_b, z50));
        r.cascade(&tp(star_a, z50), &tp(star_b, *z));
    }
    for &x in &around(lim, 3) {
        r.cascade(&tp(star_a, z50), &tp(star_b, [x, x]));
        r.cascade(&tp(star_a, [x, x]), &tp(star_b, z50));
    }
    // S21's |S21|² < 10⁻⁶⁰ on either network: a real S21 at 10⁻³⁰ ± ulps, an equal-part one at 10⁻³⁰/√2.
    for &x in &around(1e-30, 3) {
        r.two_ports(&tp([(0.1, 0.0), (0.5, 0.0), (x, 0.0), (0.2, 0.0)], z50), &tp(star_b, z50));
        r.cascade(&tp(star_a, z50), &tp([(0.1, 0.0), (0.5, 0.0), (0.0, -x), (0.2, 0.0)], z50));
    }
    for &x in &around(7.071_067_811_865_475e-31, 3) {
        r.cascade(&tp([(0.1, 0.0), (0.5, 0.0), (x, x), (0.2, 0.0)], z50), &tp(star_b, z50));
    }
    // T22 = 0: S22 of the first times S11 of the second is 1 (a lossless loop), and nearly so.
    let loop_pair = |s22: (f64, f64), s11: (f64, f64)| {
        (tp([(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), s22], z50), tp([s11, (1.0, 0.0), (1.0, 0.0), (0.0, 0.0)], z50))
    };
    for &(s22, s11) in &[(1.0, 1.0), (-1.0, -1.0), (0.5, 2.0), (0.0, 1.0), (1.0, 1.0 + f64::EPSILON)] {
        let (a, b) = loop_pair((s22, 0.0), (s11, 0.0));
        r.cascade(&a, &b);
    }
    let (a, b) = loop_pair((0.0, 1.0), (0.0, -1.0));
    r.cascade(&a, &b);
    // S → Z's |det(I − S)|² < 10⁻⁶⁰: S11 = 1 − 2⁻⁵⁰ and S22 = 1 − k·2⁻⁵³, det = k·2⁻¹⁰³ across 10⁻³⁰.
    for k in 6..=14_u32 {
        let s22 = 1.0 - f64::from(k) * f64::EPSILON / 2.0;
        r.s_to_z(&tp([(1.0 - 4.0 * f64::EPSILON, 0.0), (0.0, 0.0), (0.0, 0.0), (s22, 0.0)], z50));
        r.s_to_z(&tp([(1.0, 1e-15), (0.0, 0.0), (0.0, 0.0), (s22, 0.0)], z50));
    }
    // The specials in each entry, against the star product's pair.
    for slot in 0..8 {
        for &x in &specials() {
            let mut s = star_a;
            if slot % 2 == 0 {
                s[slot / 2].0 = x;
            } else {
                s[slot / 2].1 = x;
            }
            r.two_ports(&tp(s, z50), &tp(star_b, z50));
        }
    }
    // Seeded networks: passive ones, and parts over wide ranges.
    for _ in 0..150 {
        let mut s = [(0.0, 0.0); 4];
        for e in &mut s {
            *e = (rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0));
        }
        let mut t = [(0.0, 0.0); 4];
        for e in &mut t {
            *e = (rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0));
        }
        let z = rng.log_uniform(0.0, 3.0);
        r.two_ports(&tp(s, [z, z]), &tp(t, [z, z]));
    }
    for _ in 0..60 {
        let mut s = [(0.0, 0.0); 4];
        for e in &mut s {
            *e = (rng.signed_log_uniform(-300.0, 300.0), rng.signed_log_uniform(-300.0, 300.0));
        }
        r.two_ports(&tp(s, z50), &tp(star_b, z50));
    }
}

// ── M4 part 1c-1: reading Touchstone ─────────────────────────────────────────────────────────────────────────────

/// A parse's tokens: `i:<ports>`, `f:<z0>`, `v:<points>` and the frequencies, then each matrix (`sm_toks`).
fn ts_toks(d: &TouchstoneData) -> Vec<Tok> {
    let mut v = vec![i(word(d.ports)), f(d.z0), Tok::V(d.frequencies.len())];
    v.extend(d.frequencies.iter().map(|&x| f(x)));
    for m in &d.s_matrices {
        v.extend(sm_toks(m));
    }
    v
}

impl Rows {
    fn parse(&mut self, content: &str, ports: usize) {
        let ins = [s(content), i(word(ports))];
        self.row("rf::parse_touchstone", Class::C, &ins, || {
            res(parse_touchstone(black_box(content), black_box(ports)), ts_toks)
        });
    }
}

/// A token for the reader's number grammar, in an RI record of a 1-port (its value lands in S11's real part as
/// parsed) and as a frequency in Hz.
fn tok_rows(r: &mut Rows, tok: &str) {
    r.parse(&format!("# Hz S RI R 50\n1 {tok} 0\n"), 1);
    r.parse(&format!("# Hz S RI R 50\n{tok} 0.5 0\n"), 1);
}

/// Part 1c-1's grids: the Rust tests' files; every option (units, parameter types, formats, R and its values, case,
/// order, unknown tokens, a second option line, one after the data or the noise block); Rust's float grammar token by
/// token, with the long, tied and range-edge decimals that need a correctly rounded parse; comments, blank and
/// whitespace-only lines, CRLF and the Unicode White_Space characters; records split across lines, partial and
/// unrepresentable ones; negative, zero and overflowing frequencies; the 2-port noise block; the port limits; error
/// snippets past 32 characters and with characters `{:?}` escapes; and seeded files in every format and unit.
fn emit_1c1(r: &mut Rows, rng: &mut Rng) {
    // The Rust tests' files.
    let rust_files: [(&str, usize); 18] = [
        ("# GHz S RI R 50\n1.0  0.5  0.1\n2.0  0.3  -0.2\n", 1),
        ("! comment only\n", 1),
        ("! 2-port\n# GHZ S RI R 50\n1.0  0.1 0.01  0.9 0.09  0.2 0.02  0.3 0.03\n", 2),
        (
            "# GHz S RI R 50\n1.0  0.11 -0.11  0.12 -0.12  0.13 -0.13\n     0.21 -0.21  0.22 -0.22  0.23 -0.23\n     \
             0.31 -0.31  0.32 -0.32  0.33 -0.33 ! row 3\n2.0  0.11 0 0.12 0 0.13 0\n 0.21 0 0.22 0 0.23 0\n 0.31 0 \
             0.32 0 0.33 0\n",
            3,
        ),
        ("# GHz S RI R 50 ! options\n1.0 0.5 0.1 ! first\n", 1),
        ("# GHz S RI R 50\n1.0 0.5 abc 0.1\n", 1),
        ("# GHz S RI R 50\n1.0 nan 0.1\n", 1),
        ("# GHz S RI R 50\n1.0 inf 0.1\n", 1),
        ("# GHz S RI R 50\n-1.0 0.5 0.1\n", 1),
        ("# GHz S DB R 50\n1.0 7000 0\n", 1),
        ("# GHz S RI R 50\n1.0 0.5 0.1\n2.0 0.3\n", 1),
        ("# GHz S RI R 50\n1.0 0.1 0 0.9 0 0.2 0\n", 2),
        ("# GHz S RI R 50\n# MHz S MA R 75\n1.0 0.5 0.1\n", 1),
        ("# GHz S RI XYZ R 50\n1.0 0.5 0.1\n", 1),
        ("1.0 0.5 0.1\n# GHz S RI R 50\n2.0 0.5 0.1\n", 1),
        ("[Version] 2.0\n# GHz S RI R 50\n1.0 0.5 0.1\n", 1),
        ("1.0 0.5 90\n", 1),
        ("# MHZ S RI R 75\n100.0  0.1  0.2\n", 1),
    ];
    for (c, n) in rust_files {
        r.parse(c, n);
    }
    for p in ["Y", "Z", "H", "G", "y", "s", "S"] {
        r.parse(&format!("# GHz {p} RI R 50\n1.0 0.5 0.1\n"), 1);
    }
    for v in ["", "abc", "-50", "0", "nan", "inf", "-0", "1e400", "5e-324", "50.0", "+50", "5e1", "1_0", "50 R 75"] {
        r.parse(&format!("# GHz S RI R {v}\n1.0 0.5 0.1\n"), 1);
    }
    for n in [0, 100_000, usize::MAX, 257, 1 << 63] {
        r.parse("1.0 0.5 0.1\n", n);
    }
    r.parse("1.0 0.5 0.1\n", 256);
    r.parse(
        "# GHZ S MA R 50\n1.0 0.5 10 0.9 -20 0.1 30 0.4 40\n2.0 0.6 11 0.8 -21 0.1 31 0.3 41\n! noise parameters\n\
         1.0 1.5 0.3 45 0.2\n2.0 1.8 0.35 50 0.25\n",
        2,
    );
    // The five-port wrapped file of test_parse_touchstone_five_port_wrapped, as Rust builds it.
    let mut five = String::from("# MHz S MA R 75\n");
    for (fi, fr) in [100.0_f64, 200.0].iter().enumerate() {
        for row in 0..5_usize {
            let pairs: Vec<String> = (0..5_usize)
                .map(|j| format!("{} {}", f64::from(u32::try_from(10 * (row + 1) + j + 1).unwrap_or(0)) / 100.0, fi))
                .collect();
            let head = if row == 0 { format!("{fr}") } else { String::new() };
            five.push_str(&format!("{head} {}\n", pairs[..4].join(" ")));
            five.push_str(&format!("   {}\n", pairs[4]));
        }
    }
    r.parse(&five, 5);
    // The bench's 201-point 2-port file.
    let mut s2p = String::from("! bench\n# GHz S MA R 50\n");
    for k in 0..201_u32 {
        let fr = 1.0 + f64::from(k) * 0.01;
        s2p.push_str(&format!("{fr} 0.1 -30 0.9 -45 0.01 60 0.2 15\n"));
    }
    r.parse(&s2p, 2);

    // Every option: units, formats, case and order, the defaults, S repeated, R anywhere.
    for unit in ["Hz", "kHz", "MHz", "GHz", "hz", "KHZ", "mhz", "gHz"] {
        for fmt in ["RI", "MA", "DB", "ri", "Ma", "db"] {
            r.parse(&format!("# {unit} S {fmt} R 75\n1.5 0.25 -45\n2.5 -0.125 30\n"), 1);
        }
    }
    for opts in [
        "#", "# ", "#GHz S RI R 50", "# R 25 RI S MHz", "# S S S", "# RI MA DB", "# GHz GHz Hz", "# R 50 R 75",
        "# r 75", "# R", "# S2P", "# GHz S RI R 50 extra", "#\tGHz\tS\tRI\tR\t50", "# GHz S RI R\u{a0}50",
        "# GHz S RI R 50\u{3000}", "#! all comment", "# GHz # S RI", "# \u{1F600}", "# MHz\u{0}",
    ] {
        r.parse(&format!("{opts}\n1.0 0.5 0.1\n"), 1);
    }
    // Option lines after data, after a partial record, after the noise block, and repeated.
    r.parse("1.0 0.5\n# GHz S RI R 50\n0.1\n", 1);
    r.parse("# GHz S RI R 50\n1.0 0.5 0.1\n# MHz\n", 1);
    r.parse("# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n1 0 0 0 0\n# MHz\n", 2);
    r.parse("# GHz S RI R 50\n# bogus option\n1.0 0.5 0.1\n", 1);
    r.parse("# GHz S RI R 50\n1.0 0.5 0.1\n[Network Data]\n", 1);
    r.parse("   [Version] 2.0\n1.0 0.5 0.1\n", 1);
    r.parse("! [Version] 2.0\n1.0 0.5 0.1\n", 1);
    r.parse("1.0 0.5 0.1 [x]\n", 1);

    // Rust's float grammar, token by token.
    for tok in [
        "1", "1.", ".5", "+1", "-0", "1.e5", "1E5", "1e+5", "-.5e-3", "00001.2500", "+.5",
        "1e0000000000000000000000000001",
        "1.0e-99999999999999999999", "0e0", "-0.0e-0", "9007199254740993", "9007199254740993.0000000000000000000001",
        "1.7976931348623157e308", "1.7976931348623158e308", "1.7976931348623159e308", "2.2250738585072011e-308",
        "4.9406564584124654e-324", "2.4703282292062328e-324", "2.4703282292062327e-324", "1e-325", "1e400", "-1e400",
        ".", "e5", "1e", "1e+", "+-1", ".e5", "inf", "Infinity", "+inf", "-NaN", "nan", "infinit", "1_0", "0x10", "١",
        "1,5", "1.5.5", "--1", "1e5.5", "+", "-", "1d5", "1e-", "0.5e", "\u{ff11}", "1\u{301}", "\u{feff}1",
    ] {
        tok_rows(r, tok);
    }
    // Long decimals: half of the smallest subnormal written out (a tie, to 0) and one digit past it (to 5e-324); a
    // tie between DBL_MAX and the next binade; 2^53 + 1 with long tails; a leading-zero run against a large exponent.
    let half_min = "\
        2.470328229206232720882843964341106861825299013071623822127928412503377536351043759326499181808179961898982823\
        47722858865463328355177969898199387398005390939063150356595155702263922908583924491051844359318028499365361525\
        00319370457678249219365623669863658480757001585769269903706311928279558551332927834338409351978015531246597263\
        57957462276646527282722005637400648549997709659947045402082816622623785739345073633900796776193057750674017632\
        46736009689513405355374585166611342237666786041621596804619144672918403005300575308490487653917113865916462395\
        24912623653881879636239373280423891018672348497668235089863388587925628302755995657524455507255189313690836254\
        7791869486679949683240497058210285131854513962138377228261454376934125320985913276672363281255e-324";
    tok_rows(r, half_min);
    tok_rows(r, &half_min.replace("1255e-324", "125e-324"));
    tok_rows(r, &format!("0.{}1e360", "0".repeat(350)));
    tok_rows(r, &format!("{}e-400", "1".repeat(400)));
    tok_rows(r, &format!("1{}e-300", "0".repeat(350)));
    tok_rows(r, &format!("9007199254740992.{}1", "0".repeat(60)));
    let max_tie = "\
        17976931348623158079372897140530341507993413271003782693617377898044496829276475094664901797758720709633028641\
        66928879109465555478519404026306574886715058206819089020007083836762738548458177115317644757302700698555713669\
        59622842914819860834936475292719074168444365510704342711559699508093042880177904174497791";
    tok_rows(r, &format!("{max_tie}.{}", "9".repeat(20)));
    // Comments, blank lines, line endings and the Unicode White_Space characters around and between tokens.
    for ws in [
        "\t", "\u{b}", "\u{c}", "\r", "\u{85}", "\u{a0}", "\u{1680}", "\u{2000}", "\u{200a}", "\u{2028}",
        "\u{2029}", "\u{202f}", "\u{205f}", "\u{3000}", "\u{200b}", "\u{feff}", "\u{180e}",
    ] {
        r.parse(&format!("#{ws}GHz S RI R 50{ws}\n{ws}1.0{ws}0.5{ws}0.1{ws}\n"), 1);
    }
    for c in [
        "", "\n", "\n\n\n", "   \n\t\n", "!\n!!\n", "# GHz S RI R 50\n", "# GHz S RI R 50\r\n1.0 0.5 0.1\r\n",
        "# GHz S RI R 50\r\n1.0 0.5 0.1\r", "1.0 0.5 0.1", "1.0 0.5 0.1\r\r\n", "1.0!0.5 0.1\n0.5 0.1\n",
        "1.0 0.5!\n0.1\n", "\u{2028}1.0 0.5 0.1\n", "1.0\u{85}0.5\u{2029}0.1\n",
    ] {
        r.parse(c, 1);
    }
    // Records across lines, partial, too many values on the last line, unrepresentable entries.
    for c in [
        "# GHz S RI R 50\n1.0\n0.5\n0.1\n2.0 0.3 -0.2 3.0\n0.25 0.5\n",
        "# GHz S RI R 50\n1.0 0.5 0.1 2.0\n",
        "# GHz S DB R 50\n1.0 6200 0\n",
        "# GHz S DB R 50\n1.0 6164 0\n2.0 6165 0\n",
        "# GHz S DB R 50\n1.0 -6500 0\n",
        "# GHz S DB R 50\n1.0 -400 45\n",
        "# GHz S MA R 50\n1.0 1.7976931348623157e308 45\n",
        "# GHz S MA R 50\n1.0 1.7976931348623157e308 1e300\n",
        "# GHz S MA R 50\n1.0 -1e308 180\n",
        "# GHz S MA R 50\n1.0 0.5 1e22\n",
        "# GHz S RI R 50\n1.0 1e308 -1e308\n",
    ] {
        r.parse(c, 1);
    }
    r.parse("# GHz S DB R 50\n1 0 0 6200 0\n0 0 0 0\n", 2);
    r.parse("# GHz S RI R 50\n1 0 0 0 0\n0 0 0 0\n1 0 0\n", 2);
    // Frequencies: negative, −0, zero, overflowing through the unit, decreasing (1-port and 3-port: accepted).
    for c in [
        "# GHz S RI R 50\n-0.0 0.5 0.1\n", "# GHz S RI R 50\n0 0.5 0.1\n", "# GHz S RI R 50\n1e300 0.5 0.1\n",
        "# Hz S RI R 50\n1.7976931348623157e308 0.5 0.1\n", "# kHz S RI R 50\n1.7976931348623157e305 0.5 0.1\n",
        "# GHz S RI R 50\n2.0 0.5 0.1\n1.0 0.3 0.2\n1.0 0.1 0.1\n", "# GHz S RI R 50\n1.0 0.5 0.1\n-0.5 0.5 0.1\n",
    ] {
        r.parse(c, 1);
    }
    r.parse("# GHz S RI R 50\n2 1 0 2 0 3 0 4 0 5 0 6 0 7 0 8 0 9 0\n1 1 0 2 0 3 0 4 0 5 0 6 0 7 0 8 0 9 0\n", 3);
    // The 2-port noise block: an equal frequency starts it, its tokens still parse, the record left open by it.
    for c in [
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n1 0.1 0 0.9 0 0.2 0 0.3 0\n",
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n0.5 1 2 3\nabc\n",
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n0.5 1 2 3 inf\n",
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n0.5 1 2\n[Noise]\n",
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n2 0.1 0 0.9\n0 0.2 0 0.3 0\n1.5 0 0 0 0 0 0 0 0\n",
        "# GHz S RI R 50\n-0.0 0.1 0 0.9 0 0.2 0 0.3 0\n0 0.1 0 0.9 0 0.2 0 0.3 0\n",
        "# GHz S RI R 50\n1 0.1 0 0.9 0 0.2 0 0.3 0\n0.5 -1\n",
    ] {
        r.parse(c, 2);
    }
    // Error snippets: tokens past 32 characters (ASCII and multi-byte), and characters `{:?}` escapes.
    for tok in [
        "x".repeat(32), "x".repeat(33), "é".repeat(40), format!("{}\u{301}", "a".repeat(31)), "a\"b\\c".to_string(),
        "\u{0}\u{1}\u{7f}\u{80}".to_string(), "\u{301}x".to_string(), "\u{200b}\u{e000}\u{10ffff}".to_string(),
        "'quote'".to_string(), "\u{1F600}\u{1F1FA}\u{1F1F8}".to_string(), "\u{ad}\u{378}\u{fff9}".to_string(),
    ] {
        r.parse(&format!("1.0 {tok} 0.1\n"), 1);
        r.parse(&format!("# GHz S RI R {tok}\n1.0 0.5 0.1\n"), 1);
        r.parse(&format!("# GHz {tok}\n1.0 0.5 0.1\n"), 1);
    }
    // Seeded files: 1- to 4-port, every format and unit, values in Rust's three float forms.
    for _ in 0..80 {
        let n = usize::try_from(rng.int(1, 4)).unwrap_or(1);
        let fmt = ["RI", "MA", "DB"][usize::try_from(rng.int(0, 2)).unwrap_or(0)];
        let unit = ["Hz", "kHz", "MHz", "GHz"][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
        let z0 = rng.log_uniform(0.0, 3.0);
        let mut c = format!("! seeded\n# {unit} S {fmt} R {z0}\n");
        let mut fr = 0.0;
        for _ in 0..usize::try_from(rng.int(1, 4)).unwrap_or(1) {
            fr += rng.uniform(0.0, 10.0);
            let mut vals = vec![fr];
            for _ in 0..n * n {
                match fmt {
                    "RI" => vals.extend([rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)]),
                    "MA" => vals.extend([rng.uniform(0.0, 1.0), rng.uniform(-180.0, 180.0)]),
                    _ => vals.extend([rng.uniform(-80.0, 0.0), rng.uniform(-180.0, 180.0)]),
                }
            }
            let words: Vec<String> = vals
                .iter()
                .map(|&x| match rng.int(0, 2) {
                    0 => format!("{x}"),
                    1 => format!("{x:e}"),
                    _ => format!("{x:?}"),
                })
                .collect();
            for (k, w) in words.iter().enumerate() {
                c.push_str(w);
                c.push(if k % 9 == 8 { '\n' } else { ' ' });
            }
            c.push('\n');
        }
        r.parse(&c, n);
    }
}

// ── M4 part 1c-2: writing Touchstone ─────────────────────────────────────────────────────────────────────────────

/// The formats in the port's enum order (BjTouchstoneFormat: MagnitudeAngle 0, RealImag 1, DbAngle 2).
const FORMATS: [TouchstoneFormat; 3] =
    [TouchstoneFormat::MagnitudeAngle, TouchstoneFormat::RealImag, TouchstoneFormat::DbAngle];

/// Data of `ports` ports at the given frequencies, each matrix's entries `ent(point, row-major index)` and impedances
/// `zs(point, port)`.
fn ts_data(
    ports: usize,
    z0: f64,
    freqs: &[f64],
    ent: impl Fn(usize, usize) -> (f64, f64),
    zs: impl Fn(usize, usize) -> f64,
) -> TouchstoneData {
    let s_matrices = (0..freqs.len())
        .map(|k| {
            let s: Vec<(f64, f64)> = (0..ports * ports).map(|q| ent(k, q)).collect();
            let z: Vec<f64> = (0..ports).map(|p| zs(k, p)).collect();
            sm_of(ports, &s, &z)
        })
        .collect();
    TouchstoneData { ports, z0, frequencies: freqs.to_vec(), s_matrices }
}

impl Rows {
    fn write(&mut self, d: &TouchstoneData, fmt: usize) {
        let mut ins = ts_toks(d);
        ins.push(i(i64::try_from(fmt).unwrap_or(0)));
        self.row("rf::try_write_touchstone", Class::C, &ins, || {
            res(try_write_touchstone(black_box(d), black_box(FORMATS[fmt])), |t| vec![s(t.clone())])
        });
    }

    fn write_all(&mut self, d: &TouchstoneData) {
        for fmt in 0..3 {
            self.write(d, fmt);
        }
    }
}

/// Part 1c-2's grids: the Rust tests' data (the one-port pair, the 2-port order, 3-, 4-, 5- and 9-ports' layout, the
/// precision of low frequencies and tiny entries, a zero entry in DB), in each format; every rule of the validation
/// (the port count, R, each frequency, a matrix's impedances against R within and past 10⁻¹², ±∞ and NaN among them,
/// a non-finite entry); entries at the edges of each format (zeros of either sign and the quadrants' angles,
/// subnormals and DB's MIN_POSITIVE clamp, |S|² past DBL_MAX); R and frequencies whose Display and `{:e}` forms are
/// long; and seeded data.
fn emit_1c2(r: &mut Rows, rng: &mut Rng) {
    let one = |e: [(f64, f64); 2]| ts_data(1, 50.0, &[1e9, 2e9], move |k, _| e[k], |_, _| 50.0);
    r.write_all(&one([(0.5, 0.1), (0.3, -0.2)]));
    let two = ts_data(2, 50.0, &[1e9], |_, q| [(0.1, 0.0), (0.2, 0.0), (0.9, 0.0), (0.3, 0.0)][q], |_, _| 50.0);
    r.write_all(&two);
    for n in [3_usize, 4, 5, 9, 8, 16] {
        let d = ts_data(
            n,
            50.0,
            &[1.5e9],
            move |_, q| {
                let (row, col) = (q / n, q % n);
                let x = |v: usize| f64::from(u32::try_from(v).unwrap_or(0));
                (0.01 * x(row * n + col), -0.003 * x(col))
            },
            |_, _| 50.0,
        );
        r.write_all(&d);
    }
    r.write_all(&ts_data(1, 50.0, &[1234.5, 1e6 + 0.25], |_, _| (1e-9, -2e-12), |_, _| 50.0));
    r.write_all(&ts_data(1, 50.0, &[1e9], |_, _| (0.0, 0.0), |_, _| 50.0));
    // The validation: the port count, R, the frequencies, a matrix's impedances, its entries.
    r.write(&ts_data(0, 50.0, &[1e9], |_, _| (0.5, 0.0), |_, _| 50.0), 1);
    r.write(&ts_data(257, 50.0, &[], |_, _| (0.5, 0.0), |_, _| 50.0), 1);
    r.write(&ts_data(4, 50.0, &[], |_, _| (0.5, 0.0), |_, _| 50.0), 0);
    for z0 in [0.0, -0.0, -50.0, f64::NAN, f64::INFINITY, 5e-324, 1e300, 1e-3, 1.0 / 3.0, 75.0, f64::MAX] {
        r.write_all(&ts_data(1, z0, &[1e9], |_, _| (0.5, 0.1), move |_, _| z0));
    }
    for fr in [
        vec![0.0], vec![-0.0], vec![-1.0], vec![f64::NAN], vec![f64::INFINITY], vec![2e9, 1e9], vec![1e9, 1e9],
        vec![0.0, -0.0], vec![1e300], vec![5e-324, 1e-300], vec![1e9, 2e9, 2e9], vec![1.0, 2.0, 3.0, f64::NAN],
        vec![123_456_789.012_345_6, 9_007_199_254_740_993.0], vec![],
    ] {
        r.write(&ts_data(1, 50.0, &fr, |_, _| (0.5, 0.1), |_, _| 50.0), 1);
        r.write(&ts_data(2, 50.0, &fr, |_, q| (0.1 * f64::from(u32::try_from(q).unwrap_or(0)), 0.0), |_, _| 50.0), 0);
    }
    let lim = 50.0 * (1.0 + 1e-12);
    let mut zv = vec![75.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 0.0, -50.0, 1e300];
    zv.extend(around(lim, 3));
    zv.extend(around(50.0 * (1.0 - 1e-12), 2));
    for &z in &zv {
        let second = move |k: usize, p: usize| if k == 1 && p == 1 { z } else { 50.0 };
        r.write(&ts_data(2, 50.0, &[1e9, 2e9], |_, _| (0.5, 0.1), second), 1);
        r.write(&ts_data(1, 50.0, &[1e9], |_, _| (0.5, 0.1), move |_, _| z), 0);
    }
    for &x in &[f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let ent = move |k: usize, q: usize| if k == 1 && q == 7 { (0.0, x) } else { (0.1, 0.1) };
        r.write(&ts_data(3, 50.0, &[1e9, 2e9], ent, |_, _| 50.0), 2);
    }
    // Entries at each format's edges.
    let edges = [
        (0.0, 0.0), (-0.0, 0.0), (0.0, -0.0), (-0.0, -0.0), (-1.0, 0.0), (-1.0, -0.0), (1.0, 0.0), (0.0, 1.0),
        (0.0, -1.0), (-1.0, 1.0), (-1.0, -1.0), (5e-324, 0.0), (0.0, -5e-324), (2.2250738585072014e-308, 0.0),
        (1e-160, 1e-160), (1e-170, 0.0), (1e154, 1e154), (1.3407807929942596e154, 0.0), (1e200, 1e200), (1e300, 0.0),
        (f64::MAX, 0.0), (-f64::MAX, f64::MAX), (0.1, 0.2), (3.0, 4.0), (1e-9, -2e-12), (0.999_999_999_999_999_9, 0.0),
    ];
    for &e in &edges {
        r.write_all(&ts_data(1, 50.0, &[1e9], move |_, _| e, |_, _| 50.0));
    }
    r.write_all(&ts_data(2, 50.0, &[1e9], |_, q| edges[q * 5], |_, _| 50.0));
    // Seeded data: 1- to 5-port, passive entries, increasing frequencies over wide units.
    for _ in 0..60 {
        let n = usize::try_from(rng.int(1, 5)).unwrap_or(1);
        let points = usize::try_from(rng.int(1, 4)).unwrap_or(1);
        let z0 = rng.log_uniform(0.0, 3.0);
        let mut fr = Vec::new();
        let mut acc = rng.log_uniform(0.0, 12.0);
        for _ in 0..points {
            fr.push(acc);
            acc += rng.log_uniform(0.0, 10.0);
        }
        let vals: Vec<(f64, f64)> =
            (0..points * n * n).map(|_| (rng.uniform(-1.0, 1.0), rng.signed_log_uniform(-12.0, 0.0))).collect();
        let d = ts_data(n, z0, &fr, move |k, q| vals[k * n * n + q], move |_, _| z0);
        r.write(&d, usize::try_from(rng.int(0, 2)).unwrap_or(0));
    }
}
