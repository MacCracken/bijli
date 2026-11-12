//! `math.tsv` — what Rust **std** returns for every f64 method bijli 1.2.0's library code calls:
//! the reference for the M0 step 8 baselines (port plan §6 step 5, §9.1).
//!
//! The calls, counted in `rust-old/src` outside test code (2026-11-11): is_finite 101, abs 81,
//! sqrt 72, is_nan 62, min 40, sin 27, cos 21, max 18, exp 11, sin_cos 8, clamp 8, atan2 5,
//! log10 5, powf 4, hypot 4, atan 4, asin 4, acos 3, powi 3 (every exponent a run-time value),
//! ceil 2, to_radians 2, to_degrees 2, cbrt 2, tan 2, copysign 1; and sinh and cosh, which bijli reaches only
//! through the hisab `Complex` it re-exports (its `sin` and `cos`).
//!
//! Where each result comes from (rustc 1.99.1, x86_64-unknown-linux-gnu; port plan §5 P3, from
//! a disassembly of the release binaries):
//! - **glibc's libm** (dynamic): sin, cos, sincos (LLVM merges a `sin` and a `cos` of one
//!   argument into it; glibc 2.44's `sincos` gives the bits of the separate calls), exp, pow
//!   (`10.0.powf(y)` stays `pow`: there is no `exp10` rewrite), atan, atan2, asin, acos, tan,
//!   hypot, log10, sinh, cosh. These rows are glibc's bits (`glibc.rs`): the header records the glibc, and
//!   `--check` refuses another.
//! - **compiler-builtins' `cbrt`** (linked in, a local `t cbrt` symbol): Rust std's `f64::cbrt`
//!   resolves to the `libm` crate's port of CORE-MATH's correctly rounded cbrt, not glibc's.
//!   `validate.py` found it (glibc disagrees on 504 of the 1,084 rows); MPFR 4.2.2 agrees on all.
//! - **compiler-rt** (`__powidf2`, linked in): powi with a run-time exponent. bijli writes every
//!   exponent as a run-time value (`l`, `i`, `li`), so the exponent is black-boxed here too.
//! - **inline code**: sqrt (`sqrtsd`), abs, ceil (`roundsd`), min, max, clamp, is_nan, is_finite,
//!   copysign, to_radians and to_degrees (one multiply by the folded constant).
//!
//! Every input passes through `black_box`, so each row is a run-time call, never a compile-time
//! fold (P3: LLVM folds asin, acos, tan and cbrt of constants).
//!
//! Classes:
//! - **A**: exact on both sides — sqrt, abs, ceil, min, max, clamp, is_nan, is_finite, copysign.
//! - **B**: arithmetic — to_radians, to_degrees, powi.
//! - **C**: glibc — sin_cos, sin, cos, exp, atan2, powf, hypot, log10, atan, asin, acos, tan,
//!   sinh, cosh, cbrt. The port matches the corrected oracle (MPFR), not glibc (port plan §0.1).
//!
//! Inputs are the special values, the ranges bijli feeds each function, and seeded sweeps.

use std::hint::black_box;

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, b, f, i};

const PI: f64 = core::f64::consts::PI;

fn unary(m: &mut Module, name: &str, class: Class, fun: fn(f64) -> f64, xs: &[f64]) {
    for &x in xs {
        let y = fun(black_box(x));
        m.row(name, class, &[f(x)], &[f(y)]);
    }
}

fn binary(m: &mut Module, name: &str, class: Class, fun: fn(f64, f64) -> f64, xy: &[(f64, f64)]) {
    for &(x, y) in xy {
        let r = fun(black_box(x), black_box(y));
        m.row(name, class, &[f(x), f(y)], &[f(r)]);
    }
}

/// Every ordered pair of `xs`.
fn grid(xs: &[f64]) -> Vec<(f64, f64)> {
    xs.iter()
        .flat_map(|&x| xs.iter().map(move |&y| (x, y)))
        .collect()
}

/// Specials, ±2^k at every 16th binade, and a wide signed log-uniform sweep.
fn wide(seed: u64, n: usize) -> Vec<f64> {
    let mut v = specials();
    for k in (-1074..=1023).step_by(16) {
        let x = pow2(k);
        v.push(x);
        v.push(-x);
    }
    let mut rng = Rng::new(seed);
    v.extend((0..n).map(|_| rng.signed_log_uniform(-300.0, 300.0)));
    v
}

fn sweep(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

fn log_sweep(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| rng.log_uniform(lo_exp10, hi_exp10))
        .collect()
}

fn pairs(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<(f64, f64)> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let x = rng.signed_log_uniform(lo_exp10, hi_exp10);
            (x, rng.signed_log_uniform(lo_exp10, hi_exp10))
        })
        .collect()
}

fn uniform_pairs(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<(f64, f64)> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let x = rng.uniform(lo, hi);
            (x, rng.uniform(lo, hi))
        })
        .collect()
}

/// x, its two neighbours, and the three negated.
fn around(x: f64) -> [f64; 6] {
    [
        x,
        x.next_up(),
        x.next_down(),
        -x,
        -x.next_up(),
        -x.next_down(),
    ]
}

/// The correctly rounded sin and cos's thresholds (impetus ADR 0011, its `src/crmath.cyr`, which
/// bijli reuses), each straddled: the tiny-argument cut 2^-27, the slow path's no-reduction cut
/// 0.78, the fast path's 2^20, the cancellation cut on a reduced 2^-24 (the doubles nearest kπ/2
/// for k = 1 … 8 reduce below it), and the double closest to a multiple of π/2
/// (6381956970095103·2^797, 2^-60.9 away).
fn sin_cos_edges() -> Vec<f64> {
    let mut v = Vec::new();
    let pio2 = f64::from_bits(0x3FF9_21FB_5444_2D18);
    for x in [
        pow2(-27),
        0.78,
        pow2(20),
        6_381_956_970_095_103.0 * pow2(797),
    ] {
        v.extend(around(x));
    }
    for k in 1..=8_i32 {
        v.extend(around(f64::from(k) * pio2));
    }
    v
}

/// The correctly rounded exp's thresholds, each straddled: the 1.0 cut 2^-55, the fast path's
/// −708 and 709.7, the overflow cut 710 and the underflow cut −746.
fn exp_edges() -> Vec<f64> {
    let mut v = Vec::new();
    for x in [pow2(-55), 708.0, 709.7, 710.0, 746.0] {
        v.extend(around(x));
    }
    v
}

/// `x.powi(n)` with both operands opaque: what bijli's run-time exponents reach (`__powidf2`).
fn powi_rt(x: f64, n: i32) -> f64 {
    black_box(x).powi(black_box(n))
}

#[allow(clippy::too_many_lines)] // one grid per function, in the module comment's order
pub fn emit() -> Module {
    let mut m = Module::new("math");

    // ── exact (class A) ──
    // sqrt: field magnitudes √(x²+y²+z²) at atomic to astronomical scales, impedances, skin
    // depths, the Courant limit.
    let mut xs = wide(0x5157, 300);
    xs.extend(log_sweep(0x5158, 300, -60.0, 60.0));
    xs.extend((1..=100).map(|k| f64::from(k * k)));
    xs.extend(sweep(0x5159, 200, 0.0, 1e4));
    unary(&mut m, "f64::sqrt", Class::A, f64::sqrt, &xs);
    unary(&mut m, "f64::abs", Class::A, f64::abs, &wide(0xab5, 100));

    // ceil: Mie's n_max = ceil(x + 4·cbrt(x) + 2) and n_start. Ties, the doubles either side of
    // an integer, then sweeps.
    let mut halves: Vec<f64> = (-20..=20).map(|k| f64::from(k) + 0.5).collect();
    halves.extend((-20..=20).map(|k| f64::from(k) - f64::EPSILON * 8.0));
    halves.extend((-20..=20).map(f64::from));
    halves.extend([
        0.499_999_999_999_999_94,
        -0.499_999_999_999_999_94,
        4_503_599_627_370_495.5,
        -4_503_599_627_370_495.5,
        4_503_599_627_370_497.0,
    ]);
    halves.extend(sweep(0x7a11, 300, -1e4, 1e4));
    halves.extend(log_sweep(0x7a13, 200, -3.0, 8.0));
    halves.extend(wide(0x7a12, 100));
    unary(&mut m, "f64::ceil", Class::A, f64::ceil, &halves);

    // min / max: every special pair (NaN on either side, ±0 both ways, ±∞), equal pairs, a value
    // against 0.0 and MIN_POSITIVE both ways (bijli writes `x.max(f64::MIN_POSITIVE)`), and
    // random pairs.
    let mut mm = grid(&specials());
    let mut rng = Rng::new(0x3a3a);
    for _ in 0..60 {
        let x = rng.signed_log_uniform(-30.0, 30.0);
        mm.push((x, x));
        mm.push((x, 0.0));
        mm.push((0.0, x));
        mm.push((x, f64::MIN_POSITIVE));
        mm.push((x, f64::NAN));
    }
    mm.extend(pairs(0x3a3b, 300, -300.0, 300.0));
    binary(&mut m, "f64::max", Class::A, f64::max, &mm);
    binary(&mut m, "f64::min", Class::A, f64::min, &mm);

    // clamp: valid bounds only (min ≤ max, neither NaN; Rust panics otherwise), with the probes
    // at and around each bound, ±0, NaN and a random value. bijli clamps cosines to [−1, 1].
    let mut rng = Rng::new(0xc1a2);
    let mut bounds = vec![
        (0.0, 0.0),
        (-0.0, 0.0),
        (-0.0, -0.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (0.0, f64::INFINITY),
        (f64::NEG_INFINITY, 0.0),
    ];
    for _ in 0..80 {
        let a = rng.signed_log_uniform(-10.0, 10.0);
        let c = rng.signed_log_uniform(-10.0, 10.0);
        bounds.push(if a <= c { (a, c) } else { (c, a) });
    }
    for &(lo, hi) in &bounds {
        let probes = [
            lo,
            hi,
            0.0,
            -0.0,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            rng.signed_log_uniform(-10.0, 10.0),
            lo.next_down(),
            hi.next_up(),
        ];
        for &x in &probes {
            let r = black_box(x).clamp(black_box(lo), black_box(hi));
            m.row("f64::clamp", Class::A, &[f(x), f(lo), f(hi)], &[f(r)]);
        }
    }

    let mut preds = specials();
    preds.extend([
        f64::from_bits(0xfff8_0000_0000_0000), // -NaN
        f64::from_bits(0x7ff0_0000_0000_0001), // signalling-pattern NaN
        -f64::MAX,
        -f64::from_bits(1),
    ]);
    preds.extend(log_sweep(0xf1f1, 20, -300.0, 300.0));
    for &x in &preds {
        m.row(
            "f64::is_finite",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_finite())],
        );
        m.row(
            "f64::is_nan",
            Class::A,
            &[f(x)],
            &[b(black_box(x).is_nan())],
        );
    }

    // copysign: the special grid in both slots (±0, ±∞ and NaN of each sign as the sign source).
    let mut cs = specials();
    cs.extend([
        f64::from_bits(0xfff8_0000_0000_0000), // -NaN
        2.5,
        -2.5,
        1e-300,
        -1e300,
    ]);
    binary(&mut m, "f64::copysign", Class::A, f64::copysign, &grid(&cs));

    // ── arithmetic (class B) ──
    // to_radians / to_degrees: antenna-pattern and Smith-chart angles.
    let mut deg = specials();
    deg.extend((-720..=720).step_by(15).map(f64::from));
    deg.extend(sweep(0xde61, 300, -720.0, 720.0));
    deg.extend(wide(0xde62, 100));
    unary(&mut m, "f64::to_radians", Class::B, f64::to_radians, &deg);
    let mut rad = specials();
    rad.extend((-16..=16).map(|k| f64::from(k) * (PI / 4.0)));
    rad.extend(sweep(0xde63, 300, -4.0 * PI, 4.0 * PI));
    rad.extend(wide(0xde64, 100));
    unary(&mut m, "f64::to_degrees", Class::B, f64::to_degrees, &rad);

    // powi with a run-time exponent: Laguerre–Gauss ρ^|l| (l ≤ 20), a non-uniform grid's
    // ratio^i (i up to the cell count), the multipoles' r^l. Negative exponents for the
    // reciprocal path.
    let mut pw = specials();
    pw.extend(sweep(0x9101, 60, -10.0, 10.0));
    pw.extend(sweep(0x9104, 40, 0.9, 1.1));
    pw.extend(log_sweep(0x9102, 40, -20.0, 20.0));
    let mut exps: Vec<i32> = (-12..=24).collect();
    exps.extend([
        31,
        32,
        33,
        63,
        64,
        100,
        255,
        256,
        999,
        1000,
        4095,
        -64,
        -100,
        i32::MIN,
        i32::MAX,
    ]);
    for &x in &pw {
        for &n in &exps {
            m.row("f64::powi", Class::B, &[f(x), i(n)], &[f(powi_rt(x, n))]);
        }
    }

    // ── transcendental (class C) ──
    // sin_cos / sin / cos: dipole and array-factor angles, NF2FF phases, the modulated pulse's
    // 2π·f·(t − t₀) (thousands of radians), multipole φ, Gouy phase. Two output tokens for
    // sin_cos, sin then cos.
    let mut ang = specials();
    ang.extend(sweep(0x5101, 1000, -2.0 * PI, 2.0 * PI));
    ang.extend(sweep(0x5102, 300, -1e4, 1e4));
    ang.extend(log_sweep(0x5103, 200, -12.0, 0.0));
    ang.extend(log_sweep(0x5104, 200, 0.0, 15.0));
    ang.extend(wide(0x5105, 100));
    ang.extend(sin_cos_edges());
    for &x in &ang {
        let (s, c) = black_box(x).sin_cos();
        m.row("f64::sin_cos", Class::C, &[f(x)], &[f(s), f(c)]);
    }
    unary(&mut m, "f64::sin", Class::C, f64::sin, &ang);
    unary(&mut m, "f64::cos", Class::C, f64::cos, &ang);

    // exp: Gaussian pulses exp(−a²), RC/RL decays exp(−t/τ), the CPML's b = exp(−(σ/κ+α)·dt/ε₀),
    // beam envelopes: arguments ≤ 0 above all, down to the underflow; then the full range.
    let mut ex = specials();
    ex.extend(sweep(0xe001, 600, -50.0, 0.0));
    ex.extend(sweep(0xe002, 300, -745.2, 0.0));
    ex.extend(sweep(0xe003, 200, -1e-3, 0.0));
    ex.extend(sweep(0xe004, 200, -745.2, 709.8));
    ex.extend((-745..=709).step_by(7).map(f64::from));
    ex.extend([
        709.782_712_893_384,
        709.782_712_893_384_1,
        -708.396_418_532_264_1,
        -745.133_219_101_941_1,
        -745.133_219_101_941_2,
        1e-300,
        -1e-300,
    ]);
    ex.extend(exp_edges());
    unary(&mut m, "f64::exp", Class::C, f64::exp, &ex);

    // atan2(y, x): complex arguments (Jones phases, S-parameter angles, Mie coefficients), the
    // multipoles' φ. Every quadrant, the special grid and wide pairs.
    let mut at = grid(&specials());
    at.extend(uniform_pairs(0xa7a1, 600, -1.0, 1.0));
    at.extend(pairs(0xa7a2, 300, -10.0, 10.0));
    at.extend(pairs(0xa7a3, 200, -300.0, 300.0));
    binary(&mut m, "f64::atan2", Class::C, f64::atan2, &at);

    // hypot: Mie's |m|, complex magnitudes. Pythagorean triples (exact answers), unit and
    // physical-scale sweeps, the special grid, and pairs spanning the range.
    let mut hy = grid(&specials());
    for (a, c) in [
        (3.0, 4.0),
        (5.0, 12.0),
        (8.0, 15.0),
        (7.0, 24.0),
        (20.0, 21.0),
    ] {
        hy.push((a, c));
        hy.push((-a, c));
        hy.push((c, -a));
    }
    hy.extend(uniform_pairs(0x4791, 600, -1.0, 1.0));
    hy.extend(uniform_pairs(0x4792, 400, -1e3, 1e3));
    hy.extend(pairs(0x4793, 400, -300.0, 300.0));
    hy.extend(pairs(0x4794, 200, -160.0, -150.0));
    binary(&mut m, "f64::hypot", Class::C, f64::hypot, &hy);

    // powf: the CPML grading ρ^m (ρ ∈ [0, 1], m ∈ [1, 6]), a non-uniform grid's
    // (d_max/d_min)^(1/(n−1)), the multipoles' r^l (l a whole number), Touchstone dB 10^(v/20);
    // then random positive pairs, the special grid, and negative bases with integer exponents.
    let mut pf: Vec<(f64, f64)> = Vec::new();
    let mut rng = Rng::new(0x9001);
    for _ in 0..400 {
        pf.push((rng.unit(), rng.uniform(1.0, 6.0)));
    }
    for _ in 0..200 {
        let ratio = 1.0 + rng.log_uniform(-6.0, 3.0);
        #[allow(clippy::cast_precision_loss)]
        let n = rng.int(2, 4096) as f64;
        pf.push((ratio, 1.0 / (n - 1.0)));
    }
    for l in 0..=20 {
        for _ in 0..10 {
            pf.push((rng.log_uniform(-12.0, 12.0), f64::from(l)));
        }
    }
    for _ in 0..300 {
        pf.push((10.0, rng.uniform(-15.0, 15.0)));
    }
    for _ in 0..300 {
        let x = rng.log_uniform(-30.0, 30.0);
        pf.push((x, rng.uniform(-10.0, 10.0)));
    }
    pf.extend(grid(&specials()));
    for _ in 0..50 {
        let x = -rng.uniform(0.1, 10.0);
        #[allow(clippy::cast_precision_loss)]
        let y = rng.int(-5, 5) as f64;
        pf.push((x, y));
    }
    binary(&mut m, "f64::powf", Class::C, f64::powf, &pf);

    // log10: return loss −20·log10|Γ| (|Γ| ∈ (0, 1]), mismatch loss −10·log10(1 − |Γ|²), the
    // Touchstone dB writer (|S| ≥ MIN_POSITIVE).
    let mut lg = specials();
    lg.extend(sweep(0x1001, 400, 0.0, 1.0));
    lg.extend(log_sweep(0x1002, 400, -307.0, 3.0));
    lg.extend((-20..=20).map(pow2));
    lg.extend([10.0, 100.0, 1000.0, 0.1, 0.01, 1e-10, 1e22, 1e23]);
    lg.extend(wide(0x1003, 100));
    unary(&mut m, "f64::log10", Class::C, f64::log10, &lg);

    // atan: the Gouy phase atan(z/z_R), the stub length atan(t)/2π, Brewster's atan(n₂/n₁).
    let mut aa = specials();
    aa.extend(sweep(0xa001, 400, -10.0, 10.0));
    aa.extend(log_sweep(0xa002, 300, -20.0, 20.0));
    aa.extend(log_sweep(0xa003, 100, -20.0, 20.0).into_iter().map(|x| -x));
    aa.extend(wide(0xa004, 100));
    unary(&mut m, "f64::atan", Class::C, f64::atan, &aa);

    // asin / acos: Snell's refraction angle and the critical angle (asin), the resonator's
    // Gouy-phase acos and the antenna-array acos of products of sines and cosines: [−1, 1],
    // dense near ±1 and 0; then values outside the domain (NaN) and the specials.
    let mut unit = specials();
    unit.extend(sweep(0xa501, 600, -1.0, 1.0));
    for x in [1.0, 0.5, 2f64.sqrt() / 2.0, 3f64.sqrt() / 2.0] {
        unit.extend(around(x));
    }
    unit.extend(log_sweep(0xa502, 200, -20.0, 0.0));
    unit.extend((1..=200).map(|k| 1.0 - f64::from(k) * f64::EPSILON));
    unit.extend(sweep(0xa503, 50, -3.0, 3.0));
    unary(&mut m, "f64::asin", Class::C, f64::asin, &unit);
    unary(&mut m, "f64::acos", Class::C, f64::acos, &unit);

    // tan: the same angle grids as sin and cos.
    unary(&mut m, "f64::tan", Class::C, f64::tan, &ang);

    // cbrt: Mie's 4·cbrt(x) for x up to MIE_MAX_SIZE_PARAMETER and |m|·x up to MIE_MAX_MX; then
    // perfect cubes (exact answers) and the full range.
    let mut cb = specials();
    cb.extend(log_sweep(0xcb01, 500, -6.0, 8.0));
    cb.extend((1..=100).map(|k| f64::from(k * k * k)));
    cb.extend(wide(0xcb02, 200));
    unary(&mut m, "f64::cbrt", Class::C, f64::cbrt, &cb);

    // sinh and cosh: bijli's code calls neither, but the Complex it re-exports does (hisab's `sin` and `cos`
    // take `im.cosh()` and `im.sinh()`, glibc's), so they are pinned here too (M1 part 2d). Small arguments, the
    // port's thresholds (2^-26, 1/2, e^-|x|'s last term at 708, exp's 709.78), the overflow edge, and the range.
    let mut hy = specials();
    hy.extend(sweep(0x51a1, 400, -3.0, 3.0));
    hy.extend(log_sweep(0x51a2, 200, -12.0, 2.85));
    for x in [
        pow2(-26),
        0.5,
        22.0,
        708.0,
        709.78,
        710.475_860_073_943_9,
        711.0,
    ] {
        hy.extend(around(x));
        hy.extend(around(-x));
    }
    hy.extend(sweep(0x51a3, 200, -720.0, 720.0));
    unary(&mut m, "f64::sinh", Class::C, f64::sinh, &hy);
    unary(&mut m, "f64::cosh", Class::C, f64::cosh, &hy);

    m
}
