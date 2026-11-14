//! `beam.tsv` — rust-old/src/beam.rs (M3 part 1c-1): Gaussian beams, ABCD matrices, resonator stability and the
//! Hermite–Gaussian and Laguerre–Gaussian modes.
//!
//! Classes: the functions that call atan, acos, exp, cos or sin (`gouy_phase`, `transverse_mode_spacing`, the three
//! mode functions) are class C (`validate.py` writes the correctly rounded readings to `beam.cr.tsv`); the rest are
//! class A (+ − × ÷ and sqrt in Rust's order). A beam is its three public fields (wavelength, waist, Rayleigh range),
//! so a method row may hold fields `GaussianBeam::new` would refuse, as Rust's public fields may; an ABCD matrix is its
//! four (A, B, C, D); a resonator stability is (g₁, g₂, g₁g₂, `i:` stable). A `Result` is its value, or
//! `E:<Variant>` and the error's Display text (ADR 0002 row 4).
//!
//! The grids: beams from `new` over laser wavelengths and waists, and at the range edges where z_R under- or
//! overflows; distances at the waist, the 10⁻³⁰ test on both sides, ±z_R, the far field and the specials; mode orders
//! to 10 (and 50 to 200, where H overflows) at transverse points scaled by the waist, LG azimuthal orders including
//! i32::MIN and i32::MAX; ABCD elements from every constructor, seeded and special matrices, the 10⁻⁶⁰ denominator
//! test and the |det − 1| ≤ 10⁻¹² wavelength rule on both sides; cavities across the stability diagram's interior,
//! edges and outside, flat and zero-radius mirrors. Each value appears once per function: a row's (function, inputs)
//! is its key in the mapping files.
//!
//! M3 part 2a-2 (ADR 0002 row 34) appends the resonators' edge grid after every earlier row: cavities within ulps of
//! the stability diagram's edges (confocal, half-confocal, concentric with L at R₁ + R₂ and its neighbours, the
//! convex–concave edge L = R₁ + R₂ with R₁R₂ < 0), near-planar ones whose g rounds to 1, flat mirrors, lengths and
//! radii past 2^1020 and a subnormal length, seeded stable cavities, and mode-spacing g pairs whose product rounds to
//! 1 or underflows.
//!
//! M3 part 2a-3 (ADR 0002 row 35) appends the ABCD determinant's edge grid: determinants exactly 1, within an ulp of 1
//! (rounding to 1) and within 10⁻¹² of it, exactly ±2⁻¹⁰⁴ where the rounded AD − BC is 0, past DBL_MAX and below
//! 2^-1074, and seeded near-singular matrices, each through `transform_beam` with beams whose λ puts λ·det in range.
//!
//! M3 part 2b (ADR 0002 row 13, amended) appends the round-once grid: `new` where πw₀² over- or underflows though z_R
//! is representable, `divergence` where π·w₀ does, the confocal waist √(Lλ/2π) where Lλ does, and seeded laser beams
//! and confocal cavities.

use crate::grids::{Rng, specials};
use crate::m_field::res_f;
use crate::tsv::{Class, Module, Tok, b, err, f, i, s};
use bijli::beam::{
    AbcdMatrix, GaussianBeam, ResonatorStability, hermite_gaussian, hermite_gaussian_field,
    laguerre_gaussian_intensity, resonator_beam_waist, resonator_stability, transverse_mode_spacing,
};
use std::collections::BTreeSet;
use std::hint::black_box;

fn gb(g: &GaussianBeam) -> Vec<Tok> {
    vec![f(g.wavelength), f(g.waist), f(g.rayleigh_range)]
}

fn am(m: &AbcdMatrix) -> Vec<Tok> {
    vec![f(m.a), f(m.b), f(m.c), f(m.d)]
}

fn rs(r: &ResonatorStability) -> Vec<Tok> {
    vec![f(r.g1), f(r.g2), f(r.g_product), b(r.is_stable)]
}

fn res<T>(r: bijli::error::Result<T>, ok: impl FnOnce(&T) -> Vec<Tok>) -> Vec<Tok> {
    match r {
        Ok(v) => ok(&v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
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
}

/// A beam from its three fields, as Rust's public fields allow.
fn beam_of(v: [f64; 3]) -> GaussianBeam {
    GaussianBeam { wavelength: black_box(v[0]), waist: black_box(v[1]), rayleigh_range: black_box(v[2]) }
}

fn mat_of(v: [f64; 4]) -> AbcdMatrix {
    AbcdMatrix::new(black_box(v[0]), black_box(v[1]), black_box(v[2]), black_box(v[3]))
}

fn up(x: f64) -> f64 {
    x.next_up()
}

fn down(x: f64) -> f64 {
    x.next_down()
}

pub fn emit() -> Module {
    let mut rows = Rows { m: Module::new("beam"), seen: BTreeSet::new() };
    let r = &mut rows;
    let mut rng = Rng::new(0xbea3_0001);

    // ── GaussianBeam::new and the beams ──────────────────────────────────────────────────────
    let lams = [1064e-9, 633e-9, 10.6e-6, 1e-6, 1.0, 1e-3, 1e-300, 1e300, 1e-160, 1e160];
    let waists = [100e-6, 1e-3, 1.0, 1e-5, 1e-200, 1e-160, 1e160, 1e-300, 1e300];
    let mut lv = lams.to_vec();
    lv.extend(specials());
    let mut wv = waists.to_vec();
    wv.extend(specials());
    let mut beams: Vec<[f64; 3]> = Vec::new();
    for &lam in &lv {
        for &w0 in &wv {
            r.row("GaussianBeam::new", Class::A, &[f(lam), f(w0)], || {
                res(GaussianBeam::new(black_box(lam), black_box(w0)), gb)
            });
            if let Ok(g) = GaussianBeam::new(lam, w0) {
                beams.push([g.wavelength, g.waist, g.rayleigh_range]);
            }
        }
    }
    // Beams the constructor refuses, from their public fields: a zero, negative, NaN or infinite field.
    let lab = [1064e-9, 100e-6, std::f64::consts::PI * 100e-6 * 100e-6 / 1064e-9];
    let mut odd: Vec<[f64; 3]> = Vec::new();
    for k in 0..3 {
        for x in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut v = lab;
            v[k] = x;
            odd.push(v);
        }
    }
    // The beams the methods run on: the laboratory beam, a unit beam, beams at the range edges, the odd ones.
    let mut mb: Vec<[f64; 3]> = vec![lab];
    if let Ok(g) = GaussianBeam::new(1.0, 1.0) {
        mb.push([g.wavelength, g.waist, g.rayleigh_range]);
    }
    for v in &beams {
        if v[2] < 1e-280 || v[2] > 1e280 || v[0] == 633e-9 {
            mb.push(*v);
        }
    }
    mb.extend(odd.iter().copied());
    let zr0 = lab[2];
    let mut zs = specials();
    zs.extend([
        1e-30,
        -1e-30,
        down(1e-30),
        up(-1e-30),
        5e-324,
        0.01,
        -0.01,
        zr0,
        -zr0,
        10.0 * zr0,
        100.0 * zr0,
        1e10 * zr0,
        -1e10 * zr0,
        1e300,
        std::f64::consts::PI,
        2.3,
    ]);
    for _ in 0..30 {
        zs.push(rng.signed_log_uniform(-6.0, 3.0));
    }
    for v in &mb {
        let g = beam_of(*v);
        let bt = [f(v[0]), f(v[1]), f(v[2])];
        for &z in &zs {
            let mut ins = bt.to_vec();
            ins.push(f(z));
            r.row("GaussianBeam::q", Class::A, &ins, || {
                let q = g.q(black_box(z));
                vec![f(q[0]), f(q[1])]
            });
            r.row("GaussianBeam::spot_size", Class::A, &ins, || vec![f(g.spot_size(black_box(z)))]);
            r.row("GaussianBeam::radius_of_curvature", Class::A, &ins, || {
                res_f(g.radius_of_curvature(black_box(z)))
            });
            r.row("GaussianBeam::gouy_phase", Class::C, &ins, || vec![f(g.gouy_phase(black_box(z)))]);
            r.row("GaussianBeam::peak_intensity_factor", Class::A, &ins, || {
                vec![f(g.peak_intensity_factor(black_box(z)))]
            });
        }
        r.row("GaussianBeam::divergence", Class::A, &bt, || vec![f(g.divergence())]);
        r.row("GaussianBeam::confocal_parameter", Class::A, &bt, || vec![f(g.confocal_parameter())]);
        r.row("GaussianBeam::depth_of_focus", Class::A, &bt, || vec![f(g.depth_of_focus())]);
    }

    // ── ABCD matrices ────────────────────────────────────────────────────────────────────────
    r.row("AbcdMatrix::identity", Class::A, &[], || am(&AbcdMatrix::identity()));
    r.row("AbcdMatrix::default", Class::A, &[], || am(&AbcdMatrix::default()));
    let mut ds = specials();
    ds.extend([0.1, 0.25, 2.0, -0.5, 1e-30, down(1e-30), -1e-30, up(-1e-30), 1e-300, 1e300, 0.5, 1.5, 3.0]);
    let mut mats: Vec<[f64; 4]> = Vec::new();
    for &d in &ds {
        r.row("AbcdMatrix::free_space", Class::A, &[f(d)], || am(&AbcdMatrix::free_space(black_box(d))));
        mats.push([1.0, d, 0.0, 1.0]);
        r.row("AbcdMatrix::thin_lens", Class::A, &[f(d)], || res(AbcdMatrix::thin_lens(black_box(d)), am));
        r.row("AbcdMatrix::curved_mirror", Class::A, &[f(d)], || {
            res(AbcdMatrix::curved_mirror(black_box(d)), am)
        });
        if let Ok(m) = AbcdMatrix::thin_lens(d) {
            mats.push([m.a, m.b, m.c, m.d]);
        }
    }
    let mut ns = specials();
    ns.extend([1.5, 1.33, 2.4, 1e-30, down(1e-30), -1e-30, 1.0 + 1e-12, 1.0 + 2e-12, 1.0 - 1e-12, 3.0, 1e300]);
    for &n1 in &ns {
        for &n2 in &ns {
            r.row("AbcdMatrix::flat_interface", Class::A, &[f(n1), f(n2)], || {
                res(AbcdMatrix::flat_interface(black_box(n1), black_box(n2)), am)
            });
            if let Ok(m) = AbcdMatrix::flat_interface(n1, n2) {
                mats.push([m.a, m.b, m.c, m.d]);
            }
        }
    }
    let rads = [0.1, -0.1, 1e-30, down(1e-30), f64::NAN, f64::INFINITY, 0.0, -0.0, 2.0, 1e-300];
    for &n1 in &[1.0, 1.5, f64::NAN, 0.0] {
        for &n2 in &[1.5, 1.0, 0.0, 1e-30, down(1e-30), f64::NAN, f64::INFINITY] {
            for &rad in &rads {
                r.row("AbcdMatrix::curved_interface", Class::A, &[f(n1), f(n2), f(rad)], || {
                    res(AbcdMatrix::curved_interface(black_box(n1), black_box(n2), black_box(rad)), am)
                });
                if let Ok(m) = AbcdMatrix::curved_interface(n1, n2, rad) {
                    mats.push([m.a, m.b, m.c, m.d]);
                }
            }
        }
    }
    // Composite systems (imaging at 2f), seeded matrices, special entries, the denominator and determinant edges.
    let lens = AbcdMatrix::thin_lens(0.1).expect("a lens");
    let f2 = AbcdMatrix::free_space(0.2).compose(&lens).compose(&AbcdMatrix::free_space(0.2));
    mats.push([f2.a, f2.b, f2.c, f2.d]);
    for _ in 0..80 {
        mats.push([rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0)]);
    }
    for &x in &specials() {
        mats.push([x, 0.0, 0.0, 1.0]);
        mats.push([1.0, x, 0.0, 1.0]);
        mats.push([1.0, 0.0, x, 1.0]);
        mats.push([1.0, 0.0, 0.0, x]);
    }
    mats.extend([
        [1.0, 0.0, 0.0, 1e-30],
        [1.0, 0.0, 0.0, down(1e-30)],
        [1.0, 0.0, 0.0, up(1e-30)],
        [1.0, 0.0, 1e-31, 0.0],
        [1.0, 0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0, 1.0 + 1e-12],
        [1.0, 0.0, 0.0, up(1.0 + 1e-12)],
        [1.0, 0.0, 0.0, 1.0 - 1e-12],
        [1.0, 0.0, 0.0, down(1.0 - 1e-12)],
        [1e200, 1e200, 1e200, 1e200],
        [1e-200, 0.0, 0.0, 1e-200],
        [0.5, 0.0, 0.0, 2.0],
        [2.0, 0.0, 0.0, 0.5],
    ]);
    let mut qs: Vec<(f64, f64)> = vec![(0.0, zr0), (0.01, zr0), (-0.3, 1.0), (1.0, 0.0), (0.0, 0.0), (1e300, 1e300)];
    for &x in &specials() {
        qs.push((x, zr0));
        qs.push((0.0, x));
    }
    let mut seen_m = BTreeSet::new();
    let mats: Vec<[f64; 4]> = mats.into_iter().filter(|m| seen_m.insert(m.map(f64::to_bits))).collect();
    for (j, v) in mats.iter().enumerate() {
        let m = mat_of(*v);
        let mt = am(&m);
        r.row("AbcdMatrix::new", Class::A, &mt, || am(&m));
        r.row("AbcdMatrix::determinant", Class::A, &mt, || vec![f(m.determinant())]);
        for &(y, th) in &[(1.0, 0.1), (0.0, 0.3), (1.0, 0.0), (f64::NAN, 1.0), (1e300, 1e300)] {
            let mut ins = mt.clone();
            ins.extend([f(y), f(th)]);
            r.row("AbcdMatrix::apply_ray", Class::A, &ins, || {
                let (a, c) = m.apply_ray(black_box(y), black_box(th));
                vec![f(a), f(c)]
            });
        }
        for &(q0, q1) in &qs {
            let mut ins = mt.clone();
            ins.extend([f(q0), f(q1)]);
            r.row("AbcdMatrix::propagate_beam", Class::A, &ins, || {
                res(m.propagate_beam(black_box([q0, q1])), |q| vec![f(q[0]), f(q[1])])
            });
        }
        for w in [mats[(j * 7 + 3) % mats.len()], mats[(j * 13 + 1) % mats.len()], [1.0, 0.25, -10.0, 1.0]] {
            let o = mat_of(w);
            let mut ins = mt.clone();
            ins.extend(am(&o));
            r.row("AbcdMatrix::compose", Class::A, &ins, || am(&m.compose(black_box(&o))));
        }
        let picks = [mb[0], mb[1], mb[(j + 2) % mb.len()]];
        for bv in picks {
            let g = beam_of(bv);
            let mut ins = mt.clone();
            ins.extend(gb(&g));
            r.row("AbcdMatrix::transform_beam", Class::A, &ins, || res(m.transform_beam(black_box(&g)), gb));
            r.row("AbcdMatrix::transform_beam_with_waist_position", Class::A, &ins, || {
                res(m.transform_beam_with_waist_position(black_box(&g)), |(o, d)| {
                    let mut t = gb(o);
                    t.push(f(*d));
                    t
                })
            });
        }
    }

    // ── resonators ───────────────────────────────────────────────────────────────────────────
    let lens_l = [0.1, 0.5, 1.0, 0.0, -1.0, f64::NAN, f64::INFINITY, 5e-324, 1e300];
    let mut radii = vec![0.05, 0.2, 0.4, 0.5, 1.0, 2.0, -0.5, -1.0, 1e10, 1e-30, down(1e-30), -1e-30, 0.25, 0.1];
    radii.extend(specials());
    radii.push(0.5 * (1.0 + 1e-9));
    for &l in &lens_l {
        for &r1 in &radii {
            for &r2 in &radii {
                if l != 0.5 && l != 0.1 && r1 != r2 && r2 != 0.2 && r2 != 0.5 {
                    continue;
                }
                let ins = [f(l), f(r1), f(r2)];
                r.row("beam::resonator_stability", Class::A, &ins, || {
                    res(resonator_stability(black_box(l), black_box(r1), black_box(r2)), rs)
                });
                for &lam in &[1064e-9, 0.0, -1.0, f64::NAN, f64::INFINITY] {
                    if lam != 1064e-9 && (r1 != 0.2 || r2 != 0.2) {
                        continue;
                    }
                    let mut ins4 = ins.to_vec();
                    ins4.push(f(lam));
                    r.row("beam::resonator_beam_waist", Class::A, &ins4, || {
                        res_f(resonator_beam_waist(black_box(l), black_box(r1), black_box(r2), black_box(lam)))
                    });
                }
            }
        }
    }
    let mut gs = vec![-2.0, -1.0, -0.5, 0.0, -0.0, 0.5, 1.0, 1e-300, -1e-300, 0.999_999_999, -0.999_999_999, 2.0];
    gs.extend(specials());
    for &fsr in &[1.5e9, 1.0, 0.0, -1.0, f64::NAN, f64::INFINITY, 5e-324] {
        for &g1 in &gs {
            for &g2 in &gs {
                if fsr != 1.5e9 && (g1 != 0.5 || g2 != 0.5) {
                    continue;
                }
                r.row("beam::transverse_mode_spacing", Class::C, &[f(fsr), f(g1), f(g2)], || {
                    res_f(transverse_mode_spacing(black_box(fsr), black_box(g1), black_box(g2)))
                });
            }
        }
    }

    // ── the modes ────────────────────────────────────────────────────────────────────────────
    let mode_beams = [lab, mb[1]];
    let orders: [u32; 6] = [0, 1, 2, 3, 5, 10];
    for bv in mode_beams {
        let g = beam_of(bv);
        let w0 = bv[1];
        let zr = bv[2];
        let pts = [
            (0.0, 0.0),
            (0.5 * w0, 0.0),
            (w0, 0.0),
            (-w0, 0.5 * w0),
            (0.3 * w0, -0.2 * w0),
            (2.0 * w0, w0),
            (0.4 * w0, 0.7 * w0),
            (-50e-6 * w0 / 100e-6, 0.0),
            (3.0 * w0, -2.5 * w0),
        ];
        for &z in &[0.0, zr, -zr, 10.0 * zr, 2.3 * zr / std::f64::consts::PI] {
            for &m in &orders {
                for &n in &orders {
                    for &(x, y) in &pts {
                        let ins = [f(bv[0]), f(bv[1]), f(bv[2]), i(m), i(n), f(x), f(y), f(z)];
                        let (mb_, nb_) = (black_box(m), black_box(n));
                        let (xb, yb, zb) = (black_box(x), black_box(y), black_box(z));
                        r.row("beam::hermite_gaussian", Class::C, &ins, || {
                            vec![f(hermite_gaussian(&g, mb_, nb_, xb, yb, zb))]
                        });
                        r.row("beam::hermite_gaussian_field", Class::C, &ins, || {
                            let e = hermite_gaussian_field(&g, mb_, nb_, xb, yb, zb);
                            vec![f(e.re), f(e.im)]
                        });
                    }
                }
            }
            for &p in &orders {
                for &l in &[0i32, 1, -1, 2, -2, 3, 5, -5, 10, i32::MIN, i32::MAX] {
                    for &rr in &[0.0, 0.5 * w0, w0, 0.7 * w0, 2.0 * w0, 5.0 * w0, 1.0] {
                        let ins = [f(bv[0]), f(bv[1]), f(bv[2]), i(p), i(l), f(rr), f(z)];
                        let (pb, lb, rb, zb) = (black_box(p), black_box(l), black_box(rr), black_box(z));
                        r.row("beam::laguerre_gaussian_intensity", Class::C, &ins, || {
                            vec![f(laguerre_gaussian_intensity(&g, pb, lb, rb, zb))]
                        });
                    }
                }
            }
        }
        // High orders, where H_m overflows, and the specials in every coordinate.
        for &m in &[50u32, 100, 200] {
            for &x in &[w0, 10.0 * w0, 30.0 * w0] {
                let ins = [f(bv[0]), f(bv[1]), f(bv[2]), i(m), i(1), f(x), f(0.0), f(0.0)];
                r.row("beam::hermite_gaussian", Class::C, &ins, || {
                    vec![f(hermite_gaussian(&g, black_box(m), 1, black_box(x), 0.0, 0.0))]
                });
                r.row("beam::hermite_gaussian_field", Class::C, &ins, || {
                    let e = hermite_gaussian_field(&g, black_box(m), 1, black_box(x), 0.0, 0.0);
                    vec![f(e.re), f(e.im)]
                });
            }
        }
        for &sp in &specials() {
            for (x, y, z) in [(sp, 0.0, 0.0), (0.0, sp, 0.0), (0.3 * w0, 0.1 * w0, sp)] {
                let ins = [f(bv[0]), f(bv[1]), f(bv[2]), i(1), i(2), f(x), f(y), f(z)];
                r.row("beam::hermite_gaussian", Class::C, &ins, || {
                    vec![f(hermite_gaussian(&g, 1, 2, black_box(x), black_box(y), black_box(z)))]
                });
                r.row("beam::hermite_gaussian_field", Class::C, &ins, || {
                    let e = hermite_gaussian_field(&g, 1, 2, black_box(x), black_box(y), black_box(z));
                    vec![f(e.re), f(e.im)]
                });
            }
            for (rr, z) in [(sp, 0.0), (0.5 * w0, sp)] {
                let ins = [f(bv[0]), f(bv[1]), f(bv[2]), i(1), i(2), f(rr), f(z)];
                r.row("beam::laguerre_gaussian_intensity", Class::C, &ins, || {
                    vec![f(laguerre_gaussian_intensity(&g, 1, 2, black_box(rr), black_box(z)))]
                });
            }
        }
    }
    // The odd beams through the modes, once each.
    for v in &odd {
        let g = beam_of(*v);
        let ins = [f(v[0]), f(v[1]), f(v[2]), i(1), i(1), f(50e-6), f(-20e-6), f(0.01)];
        r.row("beam::hermite_gaussian", Class::C, &ins, || vec![f(hermite_gaussian(&g, 1, 1, 50e-6, -20e-6, 0.01))]);
        r.row("beam::hermite_gaussian_field", Class::C, &ins, || {
            let e = hermite_gaussian_field(&g, 1, 1, 50e-6, -20e-6, 0.01);
            vec![f(e.re), f(e.im)]
        });
        let ins2 = [f(v[0]), f(v[1]), f(v[2]), i(1), i(2), f(50e-6), f(0.01)];
        r.row("beam::laguerre_gaussian_intensity", Class::C, &ins2, || {
            vec![f(laguerre_gaussian_intensity(&g, 1, 2, 50e-6, 0.01))]
        });
    }
    resonator_edges(r);
    abcd_det_edges(r);
    round_once_edges(r);
    rows.m
}

/// M3 part 2b (ADR 0002 row 13, amended): the single-product functions at their range edges, and seeded values.
fn round_once_edges(r: &mut Rows) {
    let mut rng = Rng::new(0xbea3_02b0);
    let mut pairs: Vec<(f64, f64)> =
        vec![(1e10, 1e155), (1e-100, 1e-170), (1e300, 1e300), (5e-324, 1e-162), (1e-300, 1e-160), (1e-10, 1e154)];
    for _ in 0..200 {
        pairs.push((rng.log_uniform(-7.0, -4.0), rng.log_uniform(-6.0, -2.0)));
    }
    for (lam, w0) in pairs {
        r.row("GaussianBeam::new", Class::A, &[f(lam), f(w0)], || {
            res(GaussianBeam::new(black_box(lam), black_box(w0)), gb)
        });
    }
    let mut fields: Vec<[f64; 3]> =
        vec![[1e-10, 1e308, 1.0], [5e-324, 1e-300, 1.0], [1e300, 1e-10, 1.0], [-1e-6, 1e-4, 1.0], [1e-6, -1e-4, 1.0]];
    for _ in 0..200 {
        fields.push([rng.log_uniform(-7.0, -4.0), rng.log_uniform(-6.0, -2.0), 1.0]);
    }
    for v in &fields {
        let g = beam_of(*v);
        r.row("GaussianBeam::divergence", Class::A, &gb(&g), || vec![f(g.divergence())]);
    }
    let mut conf: Vec<(f64, f64)> =
        vec![(1e200, 1e200), (1e-200, 1e-200), (1e300, 1e-300), (1e-300, 1e300), (0.5, 1e300)];
    for _ in 0..100 {
        conf.push((rng.log_uniform(-3.0, 1.0), rng.log_uniform(-7.0, -4.0)));
    }
    for (l, lam) in conf {
        r.row("beam::resonator_beam_waist", Class::A, &[f(l), f(l), f(l), f(lam)], || {
            res_f(resonator_beam_waist(black_box(l), black_box(l), black_box(l), black_box(lam)))
        });
    }
}

/// M3 part 2a-3 (ADR 0002 row 35, the user's answer of 2026-11-14): `transform_beam` on the exact determinant.
fn abcd_det_edges(r: &mut Rows) {
    let e52 = 2f64.powi(-52);
    let e53 = 2f64.powi(-53);
    let mut mats: Vec<[f64; 4]> = vec![
        [2.0, 3.0, 1.0, 2.0],
        [4.0, 1.0, 3.0, 1.0],
        [3.0, 5.0, 1.0, 2.0],
        [0.5, 0.0, 0.0, 2.0],
        [1.0 + e52, 0.0, 0.0, 1.0 - e53],
        [1.0 - e53, 0.0, 0.0, 1.0 + e52],
        [1.0 + e52, 0.0, 0.0, 1.0 + e52],
        [1.0, 0.0, 0.0, 1.0 + 1e-15],
        [1.0, 0.0, 0.0, 1.0 - 1e-14],
        [1.0, 0.0, 0.0, 1.0 + 1e-13],
        [1.0 + 1e-12, 0.0, 0.0, 1.0],
        [1.0 + e52, 1.0 + 2.0 * e52, 1.0, 1.0 + e52],
        [1.0 + 2.0 * e52, 1.0 + e52, 1.0 + e52, 1.0],
        [1.0, 1.0 + e52, 1.0 + e52, 1.0 + 2.0 * e52],
        [1e200, 0.0, 0.0, 1e200],
        [1e200, 1.0, 0.0, 1e200],
        [1e-200, 0.0, 0.0, 1e-200],
        [1e-170, 0.0, 0.0, 1e-170],
        [1.0, 0.0, 0.0, 5e-324],
    ];
    let mut rng = Rng::new(0xbea3_02a3);
    for _ in 0..60 {
        let a = rng.signed_log_uniform(-2.0, 2.0);
        let b = rng.signed_log_uniform(-2.0, 2.0);
        let c = rng.signed_log_uniform(-2.0, 2.0);
        let d = b * c / a;
        mats.push([a, b, c, d]);
        mats.push([a, b, c, up(d)]);
    }
    let beams: [[f64; 3]; 3] = [[1064e-9, 100e-6, std::f64::consts::PI * 100e-6 * 100e-6 / 1064e-9], [1e-300, 1.0, 1.0],
        [1e300, 1.0, 1.0]];
    for v in &mats {
        let m = mat_of(*v);
        for bv in beams {
            let g = beam_of(bv);
            let mut ins = am(&m);
            ins.extend(gb(&g));
            r.row("AbcdMatrix::transform_beam", Class::A, &ins, || res(m.transform_beam(black_box(&g)), gb));
            r.row("AbcdMatrix::transform_beam_with_waist_position", Class::A, &ins, || {
                res(m.transform_beam_with_waist_position(black_box(&g)), |(o, d)| {
                    let mut t = gb(o);
                    t.push(f(*d));
                    t
                })
            });
        }
    }
}

/// M3 part 2a-2 (ADR 0002 row 34, the user's answer of 2026-11-14): the resonators at their exact edges.
fn resonator_edges(r: &mut Rows) {
    let big = f64::MAX;
    let mut cav: Vec<[f64; 3]> = vec![
        [0.5, 0.5, 0.5],
        [0.5, up(0.5), 0.5],
        [0.5, up(0.5), up(0.5)],
        [0.5, down(0.5), down(0.5)],
        [0.5, down(0.5), up(0.5)],
        [0.3, 1.0, 0.3],
        [0.3, -1.0, 0.3],
        [0.3, f64::INFINITY, 0.3],
        [0.3, 1e-30, 0.3],
        [0.3, f64::INFINITY, 1.0],
        [0.3, f64::NEG_INFINITY, 0.5],
        [0.3, f64::INFINITY, f64::INFINITY],
        [0.3, f64::INFINITY, -1.0],
        [0.3, f64::INFINITY, down(0.3)],
        [0.3, f64::INFINITY, up(0.3)],
        [1e307, 1e308, 1e308],
        [1e300, big, big],
        [big, big, big],
        [1e308, -big, big],
        [big / 2.0, big, f64::INFINITY],
        [1e300, 3e300, 4e300],
        [5e-324, 1.0, 1.0],
        [5e-324, 1e300, 2e300],
        [1e-310, 1e-30, 1e-30],
        [1e-31, 1e-30, 1e-30],
        [2e-30, 1e-30, 1e-30],
    ];
    // Concentric: L at R₁ + R₂ (rounded) and its neighbours, and the convex–concave edge L = R₁ + R₂ with R₁R₂ < 0.
    let pairs = [(0.1, 0.2), (0.5, 0.5), (0.25, 0.75), (1.0, 2f64.powi(-60)), (1e-20, 1e-20), (3.0, 7.0), (-0.5, 2.0)];
    for (r1, r2) in pairs {
        let s = r1 + r2;
        for l in [s, up(s), down(s), up(up(s)), down(down(s))] {
            cav.push([l, r1, r2]);
        }
    }
    for l in [2.0, up(2.0), down(2.0)] {
        cav.push([l, -0.5, 2.0]);
    }
    // Near-planar: g rounds to 1 (L/R below 2^-54), concave, convex and mixed.
    for rr in [1e6, 1e10, 1e16, 1e17, 1e30, 1e100, 1e300, big] {
        cav.push([0.1, rr, rr]);
        cav.push([0.1, rr, -rr]);
        cav.push([0.1, -rr, -rr]);
        cav.push([0.1, rr, 0.2]);
    }
    let mut rng = Rng::new(0xbea3_02a2);
    for _ in 0..120 {
        let l = rng.log_uniform(-3.0, 1.0);
        let s1 = if rng.int(0, 3) == 0 { -1.0 } else { 1.0 };
        let s2 = if rng.int(0, 3) == 0 { -1.0 } else { 1.0 };
        cav.push([l, s1 * l * rng.log_uniform(-0.5, 4.0), s2 * l * rng.log_uniform(-0.5, 4.0)]);
    }
    for c in &cav {
        let [l, r1, r2] = *c;
        r.row("beam::resonator_stability", Class::A, &[f(l), f(r1), f(r2)], || {
            res(resonator_stability(black_box(l), black_box(r1), black_box(r2)), rs)
        });
        for lam in [1064e-9, 5e-324, 1e300] {
            r.row("beam::resonator_beam_waist", Class::A, &[f(l), f(r1), f(r2), f(lam)], || {
                res_f(resonator_beam_waist(black_box(l), black_box(r1), black_box(r2), black_box(lam)))
            });
        }
    }
    let third = 1.0 / 3.0;
    let gs: [(f64, f64); 16] = [
        (1e-200, -1e-200),
        (1e-200, 1e-200),
        (-1e-200, -1e-200),
        (5e-324, 5e-324),
        (-5e-324, 5e-324),
        (-0.0, 0.5),
        (up(1.0), down(1.0)),
        (down(1.0), up(1.0)),
        (up(1.0), up(1.0)),
        (2.0, 0.5),
        (3.0, third),
        (3.0, up(third)),
        (3.0, down(third)),
        (-3.0, -third),
        (1e200, 1e-200),
        (1e300, 1e-300),
    ];
    for (g1, g2) in gs {
        r.row("beam::transverse_mode_spacing", Class::C, &[f(1.5e9), f(g1), f(g2)], || {
            res_f(transverse_mode_spacing(black_box(1.5e9), black_box(g1), black_box(g2)))
        });
    }
}
