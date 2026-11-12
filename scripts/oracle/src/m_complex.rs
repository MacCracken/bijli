//! `complex.tsv` — Rust hisab 1.4.0's `Complex`, as bijli 1.2.0 re-exports it (`bijli::polarization::Complex`),
//! and polarization.rs's `complex_real` / `complex_zero` / `complex_from_polar` and `ComplexExt::norm`.
//!
//! Each row is one operation on run-time values. A complex value is two `f:` tokens, the real part first. The pure
//! IEEE sequences (construction, conj, norm_sq, abs, the four operators and their f64 forms, neg, inv, is_zero) are
//! class A: the port evaluates the same operations in the same order, so the bits are equal. Rows that call a
//! transcendental (arg, from_polar, exp, ln, sqrt, powf, sin, cos) are class C: 0 ulp against the corrected
//! oracle, `complex.cr.tsv`, where glibc misrounded one of the row's calls (sin and cos call sinh and cosh too: M1
//! part 2d). The grid: the specials crossed with each other, values over every decade from 1e-310 to 1e308
//! with both signs (the range edges of re² + im² among them), the unit circle, and seeded uniform values.

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, b, f, s};
use bijli::polarization::{Complex, ComplexExt, complex_from_polar, complex_real, complex_zero};
use std::hint::black_box;

fn c(z: Complex) -> [crate::tsv::Tok; 2] {
    [f(z.re), f(z.im)]
}

/// The complex grid.
fn values() -> Vec<Complex> {
    let sp = specials();
    let mut zs = Vec::new();
    for &re in &sp {
        for &im in &sp {
            zs.push(Complex::new(re, im));
        }
    }
    let mut rng = Rng::new(0xc0_9c1e);
    for e in -310_i32..=308 {
        if e % 7 != 0 && e.abs() < 290 {
            continue;
        }
        let x: f64 = format!("1e{e}").parse().unwrap_or(0.0);
        let y = rng
            .signed_log_uniform(f64::from(e) - 2.0, f64::from(e) + 2.0)
            .clamp(-f64::MAX, f64::MAX);
        zs.push(Complex::new(x, y));
        zs.push(Complex::new(-y, x));
    }
    for k in 0..48 {
        let t = f64::from(k) * (std::f64::consts::PI / 24.0);
        // Points on (and near) the unit circle, built from the oracle's own cos and sin of exact angles: inputs
        // need no transcendental of their own to be reproducible (grids.rs), and these are recorded as bits.
        zs.push(Complex::new(black_box(t).cos(), black_box(t).sin()));
    }
    for _ in 0..400 {
        zs.push(Complex::new(
            rng.uniform(-10.0, 10.0),
            rng.uniform(-10.0, 10.0),
        ));
    }
    for _ in 0..200 {
        zs.push(Complex::new(
            rng.signed_log_uniform(-30.0, 30.0),
            rng.signed_log_uniform(-30.0, 30.0),
        ));
    }
    zs
}

/// Real operands: the specials, and seeded values.
fn reals() -> Vec<f64> {
    let mut xs = specials();
    xs.extend([0.5, 2.0, -3.0, 1.0 / 3.0, 1e-300, 1e300, 2.5, -0.75]);
    let mut rng = Rng::new(0x5ca1e);
    for _ in 0..40 {
        xs.push(rng.signed_log_uniform(-8.0, 8.0));
    }
    xs
}

pub fn emit() -> Module {
    let mut m = Module::new("complex");
    let zs = values();
    let xs = reals();
    let bb = |z: Complex| black_box(z);
    for &x in &xs {
        m.row(
            "Complex::from_real",
            Class::A,
            &[f(x)],
            &c(Complex::from_real(black_box(x))),
        );
        m.row(
            "polarization::complex_real",
            Class::A,
            &[f(x)],
            &c(complex_real(black_box(x))),
        );
        for &t in xs.iter().step_by(3) {
            m.row(
                "Complex::from_polar",
                Class::C,
                &[f(x), f(t)],
                &c(Complex::from_polar(black_box(x), black_box(t))),
            );
            m.row(
                "polarization::complex_from_polar",
                Class::C,
                &[f(x), f(t)],
                &c(complex_from_polar(black_box(x), black_box(t))),
            );
        }
    }
    m.row("Complex::default", Class::A, &[], &c(Complex::default()));
    m.row(
        "polarization::complex_zero",
        Class::A,
        &[],
        &c(complex_zero()),
    );
    for &z in &zs {
        let ins = c(z);
        m.row(
            "Complex::new",
            Class::A,
            &ins,
            &c(Complex::new(black_box(z.re), black_box(z.im))),
        );
        m.row("Complex::abs", Class::A, &ins, &[f(bb(z).abs())]);
        m.row(
            "ComplexExt::norm",
            Class::A,
            &ins,
            &[f(ComplexExt::norm(bb(z)))],
        );
        m.row("Complex::norm_sq", Class::A, &ins, &[f(bb(z).norm_sq())]);
        m.row("Complex::conj", Class::A, &ins, &c(bb(z).conj()));
        m.row("Complex::neg", Class::A, &ins, &c(-bb(z)));
        m.row("Complex::inv", Class::A, &ins, &c(bb(z).inv()));
        m.row("Complex::arg", Class::C, &ins, &[f(bb(z).arg())]);
        m.row("Complex::exp", Class::C, &ins, &c(bb(z).exp()));
        m.row("Complex::ln", Class::C, &ins, &c(bb(z).ln()));
        m.row("Complex::sqrt", Class::C, &ins, &c(bb(z).sqrt()));
        m.row("Complex::sin", Class::C, &ins, &c(bb(z).sin()));
        m.row("Complex::cos", Class::C, &ins, &c(bb(z).cos()));
        m.row("Complex::fmt", Class::A, &ins, &[s(format!("{}", bb(z)))]);
    }
    let mut rng = Rng::new(0xb1_7a2);
    let picks: Vec<Complex> = zs.iter().step_by(2).copied().collect();
    for (k, &a) in picks.iter().enumerate() {
        for j in 0..6 {
            let bz = picks[(k * 7 + j * 13 + 1) % picks.len()];
            let ins = [f(a.re), f(a.im), f(bz.re), f(bz.im)];
            m.row("Complex::add", Class::A, &ins, &c(bb(a) + bb(bz)));
            m.row("Complex::sub", Class::A, &ins, &c(bb(a) - bb(bz)));
            m.row("Complex::mul", Class::A, &ins, &c(bb(a) * bb(bz)));
            m.row("Complex::div", Class::A, &ins, &c(bb(a) / bb(bz)));
        }
        let x = xs[k % xs.len()];
        let ins = [f(a.re), f(a.im), f(x)];
        m.row("Complex::mul_f64", Class::A, &ins, &c(bb(a) * black_box(x)));
        m.row(
            "f64::mul_complex",
            Class::A,
            &[f(x), f(a.re), f(a.im)],
            &c(black_box(x) * bb(a)),
        );
        m.row("Complex::div_f64", Class::A, &ins, &c(bb(a) / black_box(x)));
        m.row(
            "Complex::is_zero",
            Class::A,
            &ins,
            &[b(bb(a).is_zero(black_box(x)))],
        );
        let n = [0.5, 2.0, -1.0, 3.7, -0.25, 1e-3, rng.uniform(-4.0, 4.0)][k % 7];
        m.row(
            "Complex::powf",
            Class::C,
            &[f(a.re), f(a.im), f(n)],
            &c(bb(a).powf(black_box(n))),
        );
    }
    m
}
