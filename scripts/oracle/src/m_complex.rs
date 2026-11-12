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

use crate::grids::{Rng, pow2, specials};
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
    edges(&mut m);
    m
}

/// The transcendental range edges (bijli ADR 0002 row 16, M1 part 5c-2), appended so the rows above keep their
/// places: e^re past DBL_MAX or below 2^-1022 (exp), cosh and sinh past DBL_MAX (sin, cos), ln and sqrt of an |z|
/// past DBL_MAX or below 2^-1022, sqrt's angle arg z / 2 below 2^-1022, and powf's |z|, |z|^n, arg z and n·arg z
/// past either end, with exact 53-bit midpoints of |z|^n. Every input is exact (literals, products by powers of two).
fn edges(m: &mut Module) {
    use std::f64::consts::{FRAC_PI_2, PI};
    let bb = |z: Complex| black_box(z);
    let tiny = f64::from_bits(1);
    let min = f64::MIN_POSITIVE;
    let max = f64::MAX;
    let exp_re = [
        709.0, 709.5, 709.78, 709.8, 710.0, 745.0, 1000.0, 1499.0, 1501.0, 1e10, -708.0, -708.5, -720.0, -745.0,
        -745.2, -800.0, -1000.0, -1501.0, -1e10,
    ];
    let exp_im = [0.0, -0.0, tiny, 1e-300, min, 1.0, -1.0, FRAC_PI_2, PI, 1e10, 1e300, -2.5];
    for &re in &exp_re {
        for &im in &exp_im {
            let z = Complex::new(re, im);
            m.row("Complex::exp", Class::C, &c(z), &c(bb(z).exp()));
        }
    }
    let hyp_im = [710.0, 710.4, 710.5, 711.0, 800.0, 1200.0, 1499.0, 1501.0, 1e5, 1e300];
    let hyp_re = [0.0, -0.0, tiny, 1e-300, min, 1.0, -2.0, PI, FRAC_PI_2, 1e10, 1e300];
    for &im in &hyp_im {
        for sg in [1.0, -1.0] {
            for &re in &hyp_re {
                let z = Complex::new(re, sg * im);
                m.row("Complex::sin", Class::C, &c(z), &c(bb(z).sin()));
                m.row("Complex::cos", Class::C, &c(z), &c(bb(z).cos()));
            }
        }
    }
    let k = 13.0 * pow2(1018); // (3k, 4k): |z| = 5k = 1.015625·2^1024, past DBL_MAX, exactly
    // (MAX, MAX), (MAX, 1), (2^-1074, 2^-1074) and (2^-1022, 0) are the specials grid's already.
    let mods = [
        (max, -max),
        (-max, max),
        (1e308, 1e308),
        (-1e308, 1e308),
        (1.5e308, -1.2e308),
        (3.0 * k, 4.0 * k),
        (-4.0 * k, -3.0 * k),
        (1e-310, 1e-310),
        (pow2(-1023), pow2(-1023)),
        (-1e-320, 3e-321),
        (3.0 * tiny, 4.0 * tiny),
        (1e-308, -1e-309),
        // sqrt's angle below 2^-1021: arg z subnormal, underflowed, or just above 2^-1022 with components in range
        (1.0, 1e-310),
        (4.0, 1e-310),
        (3.0, tiny),
        (1e300, 1e-300),
        (1.0, -1e-320),
        (1e200, 1e-200),
        (f64::from_bits(0x5FDF_FFFF_FFFF_FFFF), pow2(-511)),
        (1.0, pow2(-1021)),
    ];
    for &(re, im) in &mods {
        let z = Complex::new(re, im);
        m.row("Complex::ln", Class::C, &c(z), &c(bb(z).ln()));
        m.row("Complex::sqrt", Class::C, &c(z), &c(bb(z).sqrt()));
    }
    // powf: |z| past either end, |z|^n past either end with |z| normal, both angles, and exact midpoints of |z|^n:
    // (2^18 − 1)³ and (2^27 − 1)² have 54 bits, odd.
    let a3 = 262_143.0 * 262_143.0 * pow2(700);
    let a2 = 134_217_727.0 * pow2(500);
    let pzs = [
        (max, max),
        (1e308, 1e308),
        (3.0 * k, 4.0 * k),
        (tiny, tiny),
        (1e-310, 1e-310),
        (3.0 * tiny, 4.0 * tiny),
        (1e200, 0.0),
        (1e200, 1e199),
        (-1e200, 0.0),
        (1e-200, 0.0),
        (0.0, 1e-200),
        (1e154, 1e154),
        (1.0, 1e-310),
        (1e300, 1e-300),
    ];
    let ns = [0.5, 1.0, 2.0, -1.0, -2.0, 0.25, 1.5, 3.0, -0.5, 1e-20, -1.6, 0.0];
    let mut pow_rows: Vec<(f64, f64, f64)> = Vec::new();
    for &(re, im) in &pzs {
        for &n in &ns {
            pow_rows.push((re, im, n));
        }
    }
    pow_rows.extend([(a3, pow2(-1000), 1.5), (a3, 0.0, 1.5), (a2, pow2(-1000), 2.0), (a2, 0.0, 2.0)]);
    for &(re, im) in &[(-1.0, 0.0), (-2.0, 0.0), (-0.5, 0.0), (0.0, 1.0), (-1.0, 1e-300)] {
        for n in [1e308, -1e308, max, 6e307] {
            pow_rows.push((re, im, n));
        }
    }
    for &(re, im) in &[(1.0, 1.0), (3.0, -4.0)] {
        for n in [tiny, 1e-310, -1e-315, min, 4e-308] {
            pow_rows.push((re, im, n));
        }
    }
    for &(re, im, n) in &pow_rows {
        let z = Complex::new(re, im);
        m.row(
            "Complex::powf",
            Class::C,
            &[f(re), f(im), f(n)],
            &c(bb(z).powf(black_box(n))),
        );
    }
}
