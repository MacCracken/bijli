//! `scattering.tsv` — rust-old/src/scattering.rs (M3 part 1b-1): Mie scattering by a homogeneous sphere and the
//! Rayleigh limit.
//!
//! Classes: `mie` is class C (sin and cos of x, hypot of m, cbrt; `validate.py` writes the correctly rounded readings
//! to `scattering.cr.tsv`); the Rayleigh functions and `cross_section_from_efficiency` are class A (+ − × ÷ in Rust's
//! order, through hisab's `Complex` and the module's own Smith division). A complex index is its two parts. A `Result`
//! is its value, or `E:<Variant>` and the error's Display text (ADR 0002 row 4).
//!
//! A `mie` row's outputs are Q_sca, Q_ext, Q_abs, g, the coefficient count n (`i:`), then (aₖ.re, aₖ.im, bₖ.re, bₖ.im)
//! at each sampled k: every k when n ≤ 128, else the 128 indices k_j = j·(n − 1)/127 (integer division, j = 0…127)
//! followed by `i:` the FNV-1a fold of all 4n words (`digest`). The sample keeps a 10⁵-term row a few kilobytes; the
//! digest still pins every coefficient bit.
//!
//! The grids: a cross of size parameters (the ψ₁ series threshold and its neighbours, perfect cubes where Wiscombe's
//! n_max is an integer, Bohren & Huffman's sample, the Rust tests' cases) with indices (real, low and high contrast,
//! absorbing, metal-like, matched, below 1, negative, pure imaginary, signed zeros); large spheres to x = 10⁵ and the
//! |m|·x limit; the tiny-x and tiny-m range edges where the series over- or underflows or mx underflows to zero; every
//! rejection (non-finite, zero, gain, the limits); and seeded spheres. Rayleigh: radii, wavelengths and indices with
//! the specials and the range edges where m² over- or underflows. Each value appears once per function: a row's
//! (function, inputs) is its key in the mapping files.

use crate::grids::{Rng, specials};
use crate::m_field::res_f;
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::polarization::Complex;
use bijli::scattering::{
    MIE_MAX_MX, MIE_MAX_SIZE_PARAMETER, MieResult, cross_section_from_efficiency, mie, rayleigh_angular_intensity,
    rayleigh_cross_section, rayleigh_efficiency, rayleigh_phase_function,
};
use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::hint::black_box;

/// Coefficient rows written whole; above this, a sample of this many and a digest.
const SAMPLE: usize = 128;

/// The indices a `mie` row samples from n coefficients.
fn sample(n: usize) -> Vec<usize> {
    if n <= SAMPLE {
        return (0..n).collect();
    }
    (0..SAMPLE).map(|j| j * (n - 1) / (SAMPLE - 1)).collect()
}

/// FNV-1a over 64-bit words: (h ^ w) · 0x100000001b3 from 0xcbf29ce484222325, over (aₖ.re, aₖ.im, bₖ.re, bₖ.im) for
/// k = 0…n − 1.
fn digest(r: &MieResult) -> i64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for (a, b) in r.a.iter().zip(&r.b) {
        for w in [a.re, a.im, b.re, b.im] {
            h = (h ^ w.to_bits()).wrapping_mul(0x0100_0000_01b3);
        }
    }
    i64::from_ne_bytes(h.to_ne_bytes())
}

fn mie_toks(r: bijli::error::Result<MieResult>) -> Vec<Tok> {
    match r {
        Ok(m) => {
            assert_eq!(m.a.len(), m.b.len(), "mie returns as many b as a");
            let n = m.a.len();
            let mut out = vec![f(m.q_sca), f(m.q_ext), f(m.q_abs), f(m.g), i(i64::try_from(n).expect("n fits"))];
            for k in sample(n) {
                out.extend([f(m.a[k].re), f(m.a[k].im), f(m.b[k].re), f(m.b[k].im)]);
            }
            if n > SAMPLE {
                out.push(i(digest(&m)));
            }
            out
        }
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
    fn row(&mut self, name: &'static str, class: Class, ins: &[Tok], outs: impl FnOnce() -> Vec<Tok>) {
        if self.seen.insert(key(name, ins)) {
            let o = outs();
            self.m.row(name, class, ins, &o);
        }
    }

    fn mie(&mut self, x: f64, m: (f64, f64)) {
        let ins = [f(x), f(m.0), f(m.1)];
        self.row("scattering::mie", Class::C, &ins, || {
            mie_toks(mie(black_box(x), black_box(Complex::new(m.0, m.1))))
        });
    }
}

fn up(x: f64) -> f64 {
    x.next_up()
}

fn down(x: f64) -> f64 {
    x.next_down()
}

/// Size parameters of the cross: the ψ₁ series threshold and its neighbours, perfect cubes (n_max an integer before its
/// ceiling), Bohren & Huffman's sample, the Rust tests' cases.
fn sizes() -> Vec<f64> {
    vec![
        1e-8,
        1e-5,
        1e-3,
        0.01,
        0.05,
        0.0999,
        down(0.1),
        0.1,
        up(0.1),
        0.1001,
        0.2,
        0.5,
        1.0,
        1.5,
        2.0,
        3.0,
        5.0,
        8.0,
        2.0 * PI * 0.525 / 0.6328,
        7.5,
        10.0,
        2.0 * PI * 1.0 / 0.5,
        20.0,
        27.0,
        33.3,
        50.0,
    ]
}

/// Indices of the cross: real (low and high contrast, matched, below 1, negative), absorbing, metal-like, pure
/// imaginary, signed zeros in Im(m), tiny and large.
fn indices() -> Vec<(f64, f64)> {
    vec![
        (1.5, 0.0),
        (1.5, -0.0),
        (1.33, 0.0),
        (1.05, 0.0),
        (1.02, 0.0),
        (1.55, 0.0),
        (2.0, 0.0),
        (3.5, 0.0),
        (4.0, 0.0),
        (1.0, 0.0),
        (0.5, 0.0),
        (0.1, 0.0),
        (1e-3, 0.0),
        (10.0, 0.0),
        (100.0, 0.0),
        (-1.5, 0.0),
        (-1.5, -0.1),
        (1.5, -0.1),
        (1.33, -0.1),
        (1.5, -0.01),
        (1.5, -1e-10),
        (0.2, -3.0),
        (0.05, -2.0),
        (1.0, -1.0),
        (0.0, -1.0),
        (-0.0, -2.5),
    ]
}

pub fn emit() -> Module {
    let mut rows = Rows { m: Module::new("scattering"), seen: BTreeSet::new() };
    let r = &mut rows;
    let mut rng = Rng::new(0x5ca7_7e11);

    // ── mie ───────────────────────────────────────────────────────────────────────────
    for &x in &sizes() {
        for &m in &indices() {
            r.mie(x, m);
        }
    }
    // Large spheres, low and high contrast, absorbing.
    for &x in &[100.0, 200.0, 500.0, 1000.0] {
        for &m in &[(1.5, 0.0), (1.05, 0.0), (1.02, 0.0), (1.33, -0.1), (4.0, 0.0), (0.2, -3.0)] {
            r.mie(x, m);
        }
    }
    r.mie(1e4, (1.05, 0.0));
    r.mie(1e4, (1.33, -0.01));
    r.mie(MIE_MAX_SIZE_PARAMETER, (1.0001, 0.0));
    r.mie(up(MIE_MAX_SIZE_PARAMETER), (1.0001, 0.0));
    r.mie(1e10, (1.5, 0.0));
    // The |m|·x limit: exactly 10⁷ is accepted (the D_n start at about 10⁷ + 877), the next double above refused.
    r.mie(10.0, (MIE_MAX_MX / 10.0, 0.0));
    r.mie(10.0, (up(MIE_MAX_MX / 10.0), 0.0));
    r.mie(10.0, (2e6, 0.0));
    r.mie(1.0, (up(MIE_MAX_MX), 0.0));
    r.mie(1.0, (1e308, -1e308));
    r.mie(10.0, (1e308, -1e308));
    // Tiny spheres: the series under- and overflows (2/x² past DBL_MAX), mx underflows to zero.
    for &x in &[1e-320, 5e-324, 1e-310, 1e-300, 1e-200, 1e-160, 1e-155, 1e-154, 1e-150, 1e-100, 1e-50, 1e-20, 1e-12] {
        r.mie(x, (1.5, 0.0));
        r.mie(x, (1.5, -0.01));
    }
    for &m in &[(1e-320, 0.0), (5e-324, 0.0), (1e-300, 0.0), (1e-160, -1e-160), (0.0, -5e-324), (1e-20, 0.0)] {
        r.mie(1e-10, m);
        r.mie(1.0, m);
    }
    // Rejections: every special x, then every special part of m.
    for &x in &specials() {
        r.mie(x, (1.5, 0.0));
    }
    r.mie(-1.0, (1.5, 0.0));
    for &a in &specials() {
        r.mie(1.0, (a, 0.0));
        r.mie(1.0, (1.5, a));
        r.mie(1.0, (a, -0.1));
    }
    for &m in &[(0.0, 0.0), (-0.0, 0.0), (0.0, -0.0), (-0.0, -0.0), (1.5, 0.1), (1.5, 5e-324), (0.0, 1.0)] {
        r.mie(1.0, m);
    }
    // Seeded spheres: x in [10⁻⁴, 10²], Re(m) in [0.1, 10] (one in ten negative), Im(m) zero or −[10⁻⁶, 10].
    for _ in 0..240 {
        let x = rng.log_uniform(-4.0, 2.0);
        let mut re = rng.log_uniform(-1.0, 1.0);
        if rng.int(0, 9) == 0 {
            re = -re;
        }
        let im = if rng.int(0, 2) == 0 { 0.0 } else { -rng.log_uniform(-6.0, 1.0) };
        r.mie(x, (re, im));
    }

    // ── cross_section_from_efficiency ───────────────────────────────────────────────────
    let mut qs = specials();
    qs.extend([2.0, 3.105_425_531_465_88, 1e-33, 0.5]);
    let mut radii = specials();
    radii.extend([1e-6, 1e-8, 2e-8, 1e-160, 1e160, 1e-300, 1e300]);
    for &q in &qs {
        for &a in &radii {
            r.row("scattering::cross_section_from_efficiency", Class::A, &[f(q), f(a)], || {
                vec![f(cross_section_from_efficiency(black_box(q), black_box(a)))]
            });
        }
    }

    // ── Rayleigh ─────────────────────────────────────────────────────────────────────────
    let mut ms = indices();
    ms.extend([
        (1e200, 1e200),
        (1e160, 0.0),
        (1e155, 1e155),
        (1e-200, 1e-200),
        (3e-170, 0.0),
        (1e-160, -1e-160),
        (0.0, std::f64::consts::SQRT_2),
        (0.0, 0.0),
        (f64::NAN, 0.0),
        (1.5, f64::NAN),
        (f64::INFINITY, 0.0),
        (1.5, f64::NEG_INFINITY),
        (1.5, 0.1),
    ]);
    let mut lambdas = vec![500e-9, 250e-9, 1e-6];
    lambdas.extend(specials());
    lambdas.push(-500e-9);
    let mut radii_r = vec![1e-8, 2e-8, 1e-9, 1e-7];
    radii_r.extend(specials());
    radii_r.extend([-1e-8, 1e-160, 1e160, 1e-300, 1e300]);
    for &a in &radii_r {
        for &lam in &lambdas {
            for &m in &[(1.5, 0.0), (1.33, -0.1)] {
                let ins = [f(a), f(lam), f(m.0), f(m.1)];
                r.row("scattering::rayleigh_cross_section", Class::A, &ins, || {
                    res_f(rayleigh_cross_section(black_box(a), black_box(lam), black_box(Complex::new(m.0, m.1))))
                });
            }
        }
    }
    for &m in &ms {
        for &(a, lam) in &[(1e-8, 500e-9), (1e-9, 1e-6), (1e-160, 1e-150), (1e100, 1e101)] {
            let ins = [f(a), f(lam), f(m.0), f(m.1)];
            r.row("scattering::rayleigh_cross_section", Class::A, &ins, || {
                res_f(rayleigh_cross_section(black_box(a), black_box(lam), black_box(Complex::new(m.0, m.1))))
            });
        }
    }
    let mut xs = specials();
    xs.extend([0.1, 0.01, 1e-8, 1e-3, 1.0, 1e-80, 1e-78, 1e80, 1e77, -0.1]);
    for &x in &xs {
        for &m in &ms {
            let ins = [f(x), f(m.0), f(m.1)];
            r.row("scattering::rayleigh_efficiency", Class::A, &ins, || {
                res_f(rayleigh_efficiency(black_box(x), black_box(Complex::new(m.0, m.1))))
            });
        }
    }
    let mut coss = specials();
    coss.extend([0.5, -0.5, 0.3, 1e-8, 1e160, -1e160, 1e155]);
    for _ in 0..60 {
        coss.push(rng.uniform(-1.0, 1.0));
    }
    for &c in &coss {
        r.row("scattering::rayleigh_phase_function", Class::A, &[f(c)], || {
            vec![f(rayleigh_phase_function(black_box(c)))]
        });
    }
    let mut sigmas = specials();
    sigmas.extend([1e-20, 3.3e-17]);
    let mut dists = specials();
    dists.extend([2.0, 1e-30, down(1e-30), -1e-30, up(-1e-30), 1e-20, 1e160, 1e-160, -3.0]);
    for &sg in &sigmas {
        for &d in &dists {
            for &c in &[0.5, 0.3, 1.0, f64::NAN, f64::INFINITY] {
                let ins = [f(sg), f(d), f(c)];
                r.row("scattering::rayleigh_angular_intensity", Class::A, &ins, || {
                    res_f(rayleigh_angular_intensity(black_box(sg), black_box(d), black_box(c)))
                });
            }
        }
    }
    rows.m
}
