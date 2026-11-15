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

use crate::grids::{Rng, specials};
use crate::m_field::res_f;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::polarization::Complex;
use bijli::rf::{
    characteristic_impedance_lossless, characteristic_impedance_lossy, circular_array_factor, denormalize_impedance,
    gamma_from_z, input_impedance_lossless_complex, insertion_loss, linear_array_factor, mismatch_loss,
    normalize_impedance, phase_velocity, planar_array_factor, propagation_constant, quarter_wave_transformer,
    reflection_coefficient, reflection_coefficient_complex, return_loss, single_stub_match, single_stub_match_complex,
    vswr, z_from_gamma,
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

    rows.m
}
