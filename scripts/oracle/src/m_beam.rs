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
//!
//! M3 part 2c-2 (ADR 0002 row 17, amended) appends `propagate_beam`'s grid: |C·q + D|² ulp by ulp around 10⁻⁶⁰,
//! C·q₀ + D cancelling to its rounding error, a lens and the distance to its focus, near-singular matrices, exact zeros,
//! the range ends and scaled systems; and the transforms at the threshold, at a focus and through seeded systems.
//!
//! M3 part 2d-3a (ADR 0002 row 18, amended) appends the range edges of the beam's, the ABCD matrices' and the
//! resonators' chains: spot sizes and radii of curvature whose square passes DBL_MAX, Gouy phases whose z/z_R falls
//! below 2^-1022 (exact midpoints of the subnormal grid among them), peak intensity factors where π·w·w leaves the
//! range, curved interfaces whose n₁ − n₂ or n₂·R passes DBL_MAX, two-term sums whose products over- or underflow
//! (determinants, rays and products of matrices), transforms whose z_R′·λ′ leaves the range, cavities whose L/R passes
//! DBL_MAX, and mode spacings for a free spectral range below π·2^-1022.
//!
//! M3 part 2d-3b-1 (ADR 0002 row 18, amended) appends the Hermite–Gaussian amplitude's range edges, each through
//! `hermite_gaussian` and `hermite_gaussian_field`: high orders whose H_m·H_n passes DBL_MAX against a small envelope,
//! envelopes below 2^-1022 (and past e^-2300) lifted back by a large H, √2x/w and the partial products below 2^-1022,
//! w·w and x² + y² past either end, the far field where w passes DBL_MAX, and on-axis points with exact zeros.
//!
//! M3 part 2d-3b-2 (ADR 0002 row 18, amended) appends the Laguerre–Gaussian intensity's range edges: azimuthal orders
//! ±(2³¹ − 1), i32::MIN and ±1000 at the radius where ρ²^|l| and the envelope balance (its envelope e^-5·10¹⁰), radial
//! orders to 1000 whose L·L passes DBL_MAX or meets an envelope below 2^-1022, w₀·w₀, w·w and r·r past either end, the
//! axis, the far stand-in, and one row for each of beam.cyr's tests failing alone.

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
    q_out_edges(r);
    chain_edges(r);
    hg_edges(r);
    lg_edges(r);
    rows.m
}

/// A `laguerre_gaussian_intensity` row.
fn lg_row(r: &mut Rows, v: [f64; 3], p: u32, l: i32, rr: f64, z: f64) {
    let g = beam_of(v);
    let ins = [f(v[0]), f(v[1]), f(v[2]), i(p), i(l), f(rr), f(z)];
    let (pb, lb, rb, zb) = (black_box(p), black_box(l), black_box(rr), black_box(z));
    r.row("beam::laguerre_gaussian_intensity", Class::C, &ins, || {
        vec![f(laguerre_gaussian_intensity(&g, pb, lb, rb, zb))]
    });
}

/// M3 part 2d-3b-2: the Laguerre–Gaussian intensity's range edges (ADR 0002 row 18).
fn lg_edges(r: &mut Rows) {
    let unit = [1.0, 1.0, std::f64::consts::PI];
    let rt = |x: f64| (x / 2.0).sqrt();
    // ρ²^|l| against e^{−ρ²} at their balance point x = |l|·ln x (the unit beam at its waist: ρ² = 2r²), and around it.
    for rr in [162834.63426861665, 162834.63380802667, 162834.63472920662] {
        lg_row(r, unit, 0, i32::MAX, rr, 0.0);
        lg_row(r, unit, 2, i32::MAX, rr, 0.0);
    }
    for (p, l) in [(0, i32::MIN), (2, i32::MIN)] {
        lg_row(r, unit, p, l, 162834.63430812964, 0.0);
    }
    for (p, l) in [(0, 1000), (2, 1000), (2, -1000), (5, 1000)] {
        lg_row(r, unit, p, l, 67.52039125480071, 0.0);
    }
    // L·L or ρ²^|l| past DBL_MAX, an envelope below 2^-1022, at radial orders to 1000.
    for (p, l, x) in [(100, 0, 1000.0), (0, 150, 400.0), (3, 150, 400.0), (500, 3, 1800.0), (60, 100, 450.0)] {
        lg_row(r, unit, p, l, rt(x), 0.0);
    }
    // w₀·w₀, w·w and r·r past DBL_MAX; the axis; the far stand-in.
    lg_row(r, [1.0, 1e160, 1.0], 2, 3, 1e160, 0.0);
    lg_row(r, [1.0, 1e160, 1.0], 2, 3, 1e160 * 2f64.powi(-23), 0.0);
    lg_row(r, [1.0, 1e-170, 1.0], 2, 3, 1e-170, 0.0);
    lg_row(r, unit, 4, 7, 0.0, 0.0);
    lg_row(r, unit, 0, 3, 0.0, 0.0);
    lg_row(r, unit, 3, 2, 1e10, 0.0);
    // Each of beam.cyr's tests the only one to fail: the envelope (|l| = 50 at ρ² = 720), w₀·w₀ (a waist below 2^-511
    // in the far field), w₀²/w² (w/w₀ = 2^512·0.97), (−2r)·r below 2^-1022, the power (L(0.01) lifting a ρ²^200 that
    // underflowed to 0), (w₀²/w²)·ρ²^|l| below 2^-1022 lifted back by L, and the last product past DBL_MAX.
    lg_row(r, unit, 0, 50, rt(720.0), 0.0);
    let w0 = 1.3 * 2f64.powi(-520);
    lg_row(r, [1.0, w0, 1.3 * 2f64.powi(-530)], 1, 1, 1.3 * 2f64.powi(-510), 1000.0 * 1.3 * 2f64.powi(-530));
    let w1 = 1.3 * 2f64.powi(-500);
    lg_row(r, [1.0, w1, 1.0], 0, 30, w1 * 2f64.powi(512) * 0.97, 2f64.powi(512) * 0.97);
    lg_row(r, [1.0, 1.5 * 2f64.powi(-511), 1.0], 1, 1, 1.3 * 2f64.powi(-515), 0.0);
    lg_row(r, unit, 1000, 200, rt(0.01), 0.0);
    lg_row(r, [1.0, 1.0, 1.0], 400, 300, 2f64.powi(300) * rt(2f64.powf(-500.0 / 300.0)), 2f64.powi(300));
    lg_row(r, unit, 200, 30, rt(600.0), 0.0);
}

/// A `hermite_gaussian` and a `hermite_gaussian_field` row.
fn hg_rows(r: &mut Rows, v: [f64; 3], m: u32, n: u32, x: f64, y: f64, z: f64) {
    let g = beam_of(v);
    let ins = [f(v[0]), f(v[1]), f(v[2]), i(m), i(n), f(x), f(y), f(z)];
    let (mb, nb, xb, yb, zb) = (black_box(m), black_box(n), black_box(x), black_box(y), black_box(z));
    r.row("beam::hermite_gaussian", Class::C, &ins, || vec![f(hermite_gaussian(&g, mb, nb, xb, yb, zb))]);
    r.row("beam::hermite_gaussian_field", Class::C, &ins, || {
        let e = hermite_gaussian_field(&g, mb, nb, xb, yb, zb);
        vec![f(e.re), f(e.im)]
    });
}

/// M3 part 2d-3b-1: the Hermite–Gaussian amplitude's range edges (ADR 0002 row 18).
fn hg_edges(r: &mut Rows) {
    let unit = [1.0, 1.0, std::f64::consts::PI];
    let r2 = std::f64::consts::SQRT_2;
    // H_m·H_n past DBL_MAX against the envelope (the unit beam at its waist: √2x/w = s), on both sides of the range.
    for (m, n, s) in [(150, 150, 2.0), (150, 150, 3.0), (150, 150, 5.0), (150, 150, 8.0), (150, 150, 12.0)] {
        hg_rows(r, unit, m, n, s / r2, s / r2, 0.0);
    }
    for (m, n, s) in [(160, 140, 4.0), (170, 170, 3.0), (200, 200, 3.0), (300, 300, 10.0), (120, 180, 6.0)] {
        hg_rows(r, unit, m, n, s / r2, -s / r2, 0.0);
    }
    // A large H against an envelope below 2^-1022 and past e^-2300 (exp's own reduction), one order zero.
    for (m, n, s) in [(200, 0, 30.0), (200, 0, 40.0), (200, 0, 50.0), (500, 0, 40.0), (500, 0, 70.0), (1000, 1, 60.0)] {
        hg_rows(r, unit, m, n, s / r2, 0.25, 0.0);
    }
    for (m, s) in [(4000, 220.0), (3000, 200.0), (2500, 180.0)] {
        hg_rows(r, unit, m, 0, s / r2, 0.0, 0.0);
    }
    // √2x below 2^-1022 or √2x/w below it, w·w past either end, the far field (w past DBL_MAX), the partial product
    // (w₀/w)·H_m below 2^-1022 lifted back by H_n.
    hg_rows(r, [1e-6, 1e-300, 1.0], 3, 3, 1e-310, 1e-310, 0.0);
    hg_rows(r, [1e-6, 1e-300, 1.0], 3, 1, 5e-324, 2e-300, 0.0);
    hg_rows(r, [1e-6, 1e200, 1.0], 1, 1, 1e-200, 1e-200, 0.0);
    hg_rows(r, [1.0, 1e160, 1.0], 2, 2, 1e160, 1e160, 0.0);
    hg_rows(r, [1.0, 1e-170, 1.0], 2, 2, 1e-170, 1e-170, 0.0);
    hg_rows(r, [1.0, 1e-170, 1.0], 3, 0, 1e-170, 0.0, 0.0);
    hg_rows(r, [1e-6, 1e10, 1e-10], 1, 1, 1e300, 1e300, 1e300);
    hg_rows(r, [1e-6, 1e10, 1e-10], 0, 0, 0.0, 0.0, 1e300);
    hg_rows(r, [1e-6, 1e-100, 1e-300], 1, 150, 1e-110, 3e100, 1e-100);
    hg_rows(r, [1e-6, 1e-100, 1e-300], 1, 60, 1e-110, 1e100, 1e-100);
    // Each of beam.cyr's tests the only one to fail: the envelope alone below 2^-1022 (m = 100 at √2x/w = 38), w·w
    // alone past DBL_MAX (Rust's e = 1 where q = 0.0098, and where q = 2^-45), √2x alone subnormal (√2x/w normal),
    // √2x/w alone just below 2^-1022 at m = 1, and x² + y² alone below 2^-1022 with q near 1.
    hg_rows(r, unit, 100, 0, 38.0 / r2, 0.0, 0.0);
    hg_rows(r, [1.0, 1e155, 1.0], 1, 1, 7e153, 7e153, 0.0);
    hg_rows(r, [1.0, 1e155, 1.0], 1, 1, 1e155 * 2f64.powi(-23), 1e155 * 2f64.powi(-23), 0.0);
    hg_rows(r, [1.0, 2f64.powi(-200), 1.0], 1, 0, 1.3 * 2f64.powi(-1030), 2f64.powi(-200), 0.0);
    hg_rows(r, [1.0, 1.3 * 2f64.powi(123), 1.0], 1, 0, (1.0 + 1.0 / 512.0) * 2f64.powi(-900), 2f64.powi(-400), 0.0);
    let xs = f64::from_bits(0x3ff1_b083_c6c0_b46a) * 2f64.powi(-512);
    hg_rows(r, [1.0, 1.5 * 2f64.powi(-511), 1.0], 1, 1, xs, xs, 0.0);
    // On axis: exact zeros of the odd orders (Rust's chain kept), and envelopes far past the stand-in's 2^47.
    hg_rows(r, unit, 3, 2, 0.0, 0.5, 0.0);
    hg_rows(r, unit, 2, 5, -0.5, 0.0, 0.0);
    hg_rows(r, unit, 1, 1, 1e10, 1e10, 0.0);
    hg_rows(r, unit, 40, 40, 1e8, -1e8, 0.0);
}

/// The Gaussian beam's methods at z (one row each).
fn beam_rows(r: &mut Rows, v: [f64; 3], z: f64) {
    let g = beam_of(v);
    let ins = [f(v[0]), f(v[1]), f(v[2]), f(z)];
    r.row("GaussianBeam::spot_size", Class::A, &ins, || vec![f(g.spot_size(black_box(z)))]);
    r.row("GaussianBeam::radius_of_curvature", Class::A, &ins, || res_f(g.radius_of_curvature(black_box(z))));
    r.row("GaussianBeam::gouy_phase", Class::C, &ins, || vec![f(g.gouy_phase(black_box(z)))]);
    r.row("GaussianBeam::peak_intensity_factor", Class::A, &ins, || {
        vec![f(g.peak_intensity_factor(black_box(z)))]
    });
}

/// M3 part 2d-3a: the chains' range edges (ADR 0002 row 18).
fn chain_edges(r: &mut Rows) {
    let p2 = |k: i32| 2f64.powi(k);
    // Spot sizes and radii of curvature: z/z_R around 2^512, where its square passes DBL_MAX, far past it, with a zero
    // or negative waist; z_R/z the same for the radius, with |z| small enough that R is representable.
    for (w0, zr) in [(1e-6, 1e-10), (0.0, 1e-10), (-1e-6, 1e-10), (1e-150, 1.0), (1.0, 1e-300), (5e-324, 1.0)] {
        for z in [1e200, -1e200, p2(511) * zr, p2(512) * zr, up(p2(512)) * zr, 1e154 * zr, 1e300, f64::MAX] {
            beam_rows(r, [1e-6, w0, zr], z);
        }
    }
    for zr in [1e140, 1e160, 1e300, f64::MAX, p2(512), up(p2(512))] {
        for z in [1e-20, -1e-20, 1e-30, 2e-15, 1.0, 1e5, 1e-10] {
            beam_rows(r, [1e-6, 1e-3, zr], z);
        }
    }
    // Gouy phases below 2^-1022: exact midpoints (2j + 1)·2^-1075 of the subnormal grid, either sign, and quotients
    // near and off them, 2^-1022 itself and its neighbours.
    let u = 1.375;
    for j in [0.0, 1.0, 2.0, 3.0, 7.0, 1000.0, 1001.0] {
        let z = (2.0 * j + 1.0) * p2(-75) * u;
        for zz in [z, -z, up(z), down(z)] {
            beam_rows(r, [1e-6, 1e-3, p2(1000) * u], zz);
        }
    }
    for (z, zr) in [(1e-10, 1e300), (3e-15, 1e300), (5e-324, 1.0), (1e-20, 1e300), (p2(-22), p2(1000)), (1.0, 3e307)] {
        beam_rows(r, [1e-6, 1e-3, zr], z);
        beam_rows(r, [1e-6, 1e-3, zr], -z);
    }
    for z in [p2(-22), up(p2(-22)), down(p2(-22))] {
        beam_rows(r, [1e-6, 1e-3, p2(1000)], z);
    }
    // Peak intensity factors: w past DBL_MAX's root (π·w·w past DBL_MAX, 2/(πw²) subnormal), w where π·w·w is
    // subnormal and 2/(πw²) normal, and w from a ratio past 2^512.
    for w0 in [7e153, 7.6e153, 8e153, 1e154, 1e155, 1e160, 3e161, 4e161, 1e162, 5.9e-155, 6e-155, 6.2e-155, 7e-155] {
        beam_rows(r, [1e-6, w0, 1.0], 0.0);
        beam_rows(r, [1e-6, w0, 1.0], 0.5);
    }
    for z in [1e144, 1e150, 1e154, 1e158, 1e165] {
        beam_rows(r, [1e-6, 1.0, 1e-10], z);
    }
    // π·w·w rounds to 2^-1022 from a 53-bit 2^-1022 − 2^-1075, where 2/(πw²) rounds up past 2^1023.
    beam_rows(r, [1e-6, f64::from_bits(0x1ff2_0dd7_5042_9b6d), 1.0], 0.0);
    // Curved interfaces: n₁ − n₂ or n₂·R past DBL_MAX.
    let ci = [
        (1e300, 1e200, 1e200),
        (1e300, -1e200, 1e200),
        (1e308, -1e308, 1.0),
        (-1.7e308, 1.7e308, 2.0),
        (f64::MAX, -f64::MAX, -3.0),
        (1e200, 1e300, 1e10),
        (1.5, 1e300, 1e10),
        (1.5, 1e-30, 1e-30),
        (1e300, 1e300, 1e300),
        (1e-300, 1e300, 1e300),
        (2.0, f64::MAX, f64::MAX),
    ];
    for (n1, n2, rad) in ci {
        r.row("AbcdMatrix::curved_interface", Class::A, &[f(n1), f(n2), f(rad)], || {
            res(AbcdMatrix::curved_interface(black_box(n1), black_box(n2), black_box(rad)), am)
        });
    }
    // Two-term sums: products past DBL_MAX that cancel, products below 2^-1022 whose sum is subnormal or normal, and
    // mixed scales, through the determinant, a ray and a product of matrices.
    let ms: [[f64; 4]; 12] = [
        [1e200, 1e200, 1e200, 1e200],
        [1e200, -1e200, -1e200, 1e200],
        [p2(600), p2(600), p2(600), p2(600) * 1.5],
        [1e-155, 3e-160, 7e-156, 1.5e-155],
        [1.5e-155, 1e-155, 1e-155, 1.25e-155],
        [3e-160, 1e-160, 2e-160, 1e-160],
        [1e-155, -1e-155, 1e-155, 1e-155],
        [1e-200, 1e-150, 1e150, 1e-200],
        [1e300, 1e-300, 1e-300, 1e300],
        [1e160, 1e-160, 1e160, 1e-160],
        [1.0, 1e-160, 1e-160, 1.0],
        [f64::MAX, 2.0, 2.0, f64::MAX],
    ];
    let rays = [(1e-155, 1e-155), (1e200, 1e200), (1e-160, -3e-161), (1e150, 1e-150), (1.5e-308, 2e-309)];
    for (j, v) in ms.iter().enumerate() {
        let m = mat_of(*v);
        let mt = am(&m);
        r.row("AbcdMatrix::determinant", Class::A, &mt, || vec![f(m.determinant())]);
        for &(y, th) in &rays {
            let mut ins = mt.clone();
            ins.extend([f(y), f(th)]);
            r.row("AbcdMatrix::apply_ray", Class::A, &ins, || {
                let (a, c) = m.apply_ray(black_box(y), black_box(th));
                vec![f(a), f(c)]
            });
        }
        for w in [ms[(j + 3) % ms.len()], ms[(j + 7) % ms.len()], *v] {
            let o = mat_of(w);
            let mut ins = mt.clone();
            ins.extend(am(&o));
            r.row("AbcdMatrix::compose", Class::A, &ins, || am(&m.compose(black_box(&o))));
        }
    }
    // A product past DBL_MAX meeting −DBL_MAX (Rust's ∞, for 2^1024 − DBL_MAX = 2^971), and two subnormal products
    // that round alike and cancel ((5 ± 7/16)·2^-1074: Rust's 0, for 2^-1074), through each two-term sum.
    let big = f64::MAX;
    let (x1, x2, t) = (5.4375 * p2(-537), 4.5625 * p2(-537), p2(-537));
    let pairs2: [([f64; 4], [f64; 4]); 2] =
        [([p2(512), -big, 1.0, 1.0], [p2(512), 0.0, 1.0, 1.0]), ([x1, -x2, 1.0, 1.0], [1.0, t, 0.0, t])];
    for (sv, ov) in pairs2 {
        let (m, o) = (mat_of(sv), mat_of(ov));
        let mut ins = am(&m);
        ins.extend(am(&o));
        r.row("AbcdMatrix::compose", Class::A, &ins, || am(&m.compose(black_box(&o))));
    }
    for (v, y, th) in [([p2(512), -big, 0.0, 1.0], p2(512), 1.0), ([x1, -x2, 1.0, 1.0], t, t)] {
        let m = mat_of(v);
        let mut ins = am(&m);
        ins.extend([f(y), f(th)]);
        r.row("AbcdMatrix::apply_ray", Class::A, &ins, || {
            let (a, c) = m.apply_ray(black_box(y), black_box(th));
            vec![f(a), f(c)]
        });
    }
    for v in [[p2(512), big, 1.0, p2(512)], [x1, x2, t, t]] {
        let m = mat_of(v);
        r.row("AbcdMatrix::determinant", Class::A, &am(&m), || vec![f(m.determinant())]);
    }
    // Transforms: z_R′·λ′ past DBL_MAX, below 2^-1022 or 0, and its quotient by π below 2^-1022, through the identity,
    // free space and a lens; and a determinant row 18 changes, quoted by the error.
    let tb: [[f64; 3]; 7] = [
        [1e200, 5.6e199, 1e200],
        [1e-200, 5.6e-201, 1e-200],
        [1e-160, 1e-155, 1e-150],
        [3e-158, 1e-154, 1e-150],
        [p2(-511), p2(-511), p2(-511)],
        [f64::MAX, 1e154, 1e10],
        [1e-300, 1e-160, 1e-20],
    ];
    let tm: [[f64; 4]; 5] = [
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.25, 0.0, 1.0],
        [1.0, 0.0, -1.0, 1.0],
        [1e200, 1e200, 1e200, 1e200],
        [1e-155, 0.0, 0.0, 1e-155],
    ];
    for v in tm {
        let m = mat_of(v);
        for bv in tb {
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
    // Cavities whose L/R passes DBL_MAX: g = ±∞ against a g₂ of 0, an ulp of 1 off 0, ½, 1 (flat) and −1.
    let l = 1e300;
    for r1 in [1e-10, -1e-10, 1e-8] {
        for r2 in [l, up(l), down(l), 2.0 * l, f64::INFINITY, -l, 1e-10] {
            r.row("beam::resonator_stability", Class::A, &[f(l), f(r1), f(r2)], || {
                res(resonator_stability(black_box(l), black_box(r1), black_box(r2)), rs)
            });
            r.row("beam::resonator_beam_waist", Class::A, &[f(l), f(r1), f(r2), f(1064e-9)], || {
                res_f(resonator_beam_waist(black_box(l), black_box(r1), black_box(r2), black_box(1064e-9)))
            });
        }
    }
    // Mode spacings for Δν_FSR/π below 2^-1022 (Δν_FSR below π·2^-1022 ≈ 7.0e-308) and around it, and g₁g₂ below
    // 2^-1022 with an ordinary Δν_FSR.
    let pi_min = std::f64::consts::PI * p2(-1022);
    for fsr in [5e-324, 1e-323, 1e-310, 3.3e-309, 2.2e-308, pi_min, up(pi_min), down(pi_min), 7e-308] {
        for (g1, g2) in [(0.5, 0.5), (0.0, 0.0), (-0.5, -0.5), (1.0, 1.0), (0.3, 0.9), (-1.0, -1.0)] {
            r.row("beam::transverse_mode_spacing", Class::C, &[f(fsr), f(g1), f(g2)], || {
                res_f(transverse_mode_spacing(black_box(fsr), black_box(g1), black_box(g2)))
            });
        }
    }
    for (g1, g2) in [(1e-160, 1e-160), (-1e-160, -1e-160), (3e-170, 0.5), (1e-300, 1e-20)] {
        r.row("beam::transverse_mode_spacing", Class::C, &[f(1.5e9), f(g1), f(g2)], || {
            res_f(transverse_mode_spacing(black_box(1.5e9), black_box(g1), black_box(g2)))
        });
    }
}

/// A `propagate_beam` row.
fn prop_row(r: &mut Rows, v: [f64; 4], q0: f64, q1: f64) {
    let m = mat_of(v);
    let mut ins = am(&m);
    ins.extend([f(q0), f(q1)]);
    r.row("AbcdMatrix::propagate_beam", Class::A, &ins, || {
        res(m.propagate_beam(black_box([q0, q1])), |q| vec![f(q[0]), f(q[1])])
    });
}

/// The two transform rows of a matrix and a beam.
fn transform_rows(r: &mut Rows, v: [f64; 4], bv: [f64; 3]) {
    let m = mat_of(v);
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

/// M3 part 2c-2 (ADR 0002 row 17, amended; the user's answer of 2026-11-14): `propagate_beam` on the exact
/// |C·q + D|² against 10⁻⁶⁰ with q_out rounded once from the exact quotient, and the transforms through it.
fn q_out_edges(r: &mut Rows) {
    let mut rng = Rng::new(0xbea3_02c2);
    let e52 = 2f64.powi(-52);
    // The threshold without cancellation: X = D alone, then X and Y together, ulp by ulp around 10⁻⁶⁰.
    for d in [1e-30, up(1e-30), down(1e-30), up(up(1e-30)), down(down(1e-30)), -1e-30, -down(1e-30)] {
        prop_row(r, [1.0, 0.5, 0.0, d], 0.25, 2.0);
        prop_row(r, [2.0, -1.0, 1.0, d], 0.0, 0.0);
    }
    for (x, y) in [(6e-31, 8e-31), (up(6e-31), 8e-31), (6e-31, down(8e-31)), (down(6e-31), down(8e-31))] {
        prop_row(r, [1.0, 1.0, 1.0, x], 0.0, y);
        prop_row(r, [1.0, 1.0, 2.0, x], 0.0, 0.5 * y);
    }
    // C·q₀ + D cancelling: Rust's rounded denominator 0 where the exact one is 2^(k − 104) (k = 4: below 10⁻⁶⁰, refused
    // by both; k = 5 and up: above it), and C = 0.1, q₀ = 3, D = −RN(0.3·…): −2.8·10⁻¹⁷.
    for k in 0..12 {
        let c = (1.0 + e52) * 2f64.powi(k / 2);
        let q0 = (1.0 + e52) * 2f64.powi(k - k / 2);
        for q1 in [0.0, 1e-40, 1e-31] {
            prop_row(r, [1.0, 0.25, c, -(c * q0)], q0, q1);
            prop_row(r, [-3.0, 0.5, -c, c * q0], q0, q1);
        }
    }
    for q1 in [0.0, 1e-20, 2.0] {
        prop_row(r, [1.0, 0.0, 0.1, -(0.1 * 3.0)], 3.0, q1);
        prop_row(r, [2.0, -0.5, 0.1, -(0.1 * 3.0)], 3.0, q1);
    }
    for _ in 0..120 {
        let c = rng.signed_log_uniform(-3.0, 3.0);
        let q0 = rng.signed_log_uniform(-3.0, 3.0);
        let k = rng.int(-3, 3);
        let p = c * q0;
        let d = if k == 0 { -p } else { f64::from_bits(((-p).to_bits() as i64 + k) as u64) };
        let q1 = if rng.int(0, 2) == 0 { 0.0 } else { rng.log_uniform(-40.0, -10.0) };
        prop_row(r, [rng.signed_log_uniform(-1.0, 1.0), rng.signed_log_uniform(-1.0, 1.0), c, d], q0, q1);
    }
    // U·X + V·Y cancelling: a lens, then the distance to its focus (Re q_out near 0); near-singular matrices
    // (AD − BC a few ulps of AD); exact zeros and their signs; and the range ends, where Rust's chain over- or
    // underflows while q_out is a double.
    for _ in 0..80 {
        let f0 = rng.signed_log_uniform(-2.0, 0.0);
        let q0 = if rng.int(0, 1) == 0 { 0.0 } else { rng.uniform(-0.5, 0.5) };
        let q1 = rng.log_uniform(-4.0, 0.0);
        let lens = AbcdMatrix::new(1.0, 0.0, -1.0 / f0, 1.0);
        let qq = lens.propagate_beam([q0, q1]).unwrap_or([0.0, 0.0]);
        let dist = -qq[0];
        let m = AbcdMatrix::free_space(dist).compose(&lens);
        prop_row(r, [m.a, m.b, m.c, m.d], q0, q1);
    }
    for _ in 0..80 {
        let a = rng.signed_log_uniform(-2.0, 2.0);
        let b = rng.signed_log_uniform(-2.0, 2.0);
        let c = rng.signed_log_uniform(-2.0, 2.0);
        let d = b * c / a;
        let k = rng.int(-3, 3);
        let d = if d == 0.0 || k == 0 { d } else { f64::from_bits((d.to_bits() as i64 + k) as u64) };
        prop_row(r, [a, b, c, d], rng.uniform(-1.0, 1.0), rng.log_uniform(-3.0, 0.0));
    }
    for (m, q0, q1) in [
        ([2.0, 6.0, 1.0, 1.0], -3.0, 0.0),
        ([2.0, -0.0, 1.0, 1.0], -0.0, 0.0),
        ([2.0, -0.0, 1.0, 1.0], -0.0, -0.0),
        ([-2.0, 0.0, 1.0, 1.0], 0.0, -0.0),
        ([1.0, 2.0, 2.0, 4.0], 1.0, 3.0),
        ([1.0, 2.0, 2.0, 4.0], 1.0, -0.0),
        ([0.0, 0.0, 1.0, 1.0], 1.0, 1.0),
        ([-0.0, -0.0, 1.0, 1.0], 1.0, -1.0),
    ] {
        prop_row(r, m, q0, q1);
    }
    for (m, q0, q1) in [
        ([1e200, 1e200, 1e200, 1e200], 1e200, 1e200),
        ([1e300, 1e-300, 1e300, 1e300], 1e10, 1e10),
        ([1.0, 0.0, 1e-300, 1e-300], 1e-20, 1e-20),
        ([1e-300, 1e-300, 1e-300, 1e-300], 1.0, 1.0),
        ([5e-324, 5e-324, 1.0, 1e-300], 5e-324, 1e-200),
        ([f64::MAX, f64::MAX, f64::MAX, f64::MAX], f64::MAX, f64::MAX),
        ([1.0, f64::MAX, 1e-300, 1.0], 1.0, 1.0),
        ([1.0, 0.5, 2f64.powi(-200), 2f64.powi(200)], 2f64.powi(201), 2f64.powi(-201)),
    ] {
        prop_row(r, m, q0, q1);
    }
    for _ in 0..60 {
        let s = rng.int(-1000, 1000) as i32;
        let m = [rng.uniform(0.5, 2.0), rng.uniform(-1.0, 1.0) * 2f64.powi(s), rng.uniform(-1.0, 1.0) * 2f64.powi(-s),
            rng.uniform(0.5, 2.0)];
        prop_row(r, m, rng.uniform(-1.0, 1.0) * 2f64.powi(s), rng.uniform(0.0, 1.0) * 2f64.powi(s));
    }
    // The transforms (q = j·z_R): |C·z_R|² around 10⁻⁶⁰ with D = 0 (det = −BC), a lens and the distance to its focus,
    // and seeded systems of free space, lenses and interfaces, on the lab beam.
    let lab = [1064e-9, 100e-6, std::f64::consts::PI * 100e-6 * 100e-6 / 1064e-9];
    let c0 = 1e-30 / lab[2];
    for c in [c0, up(c0), down(c0), up(up(c0)), down(down(c0))] {
        transform_rows(r, [1.0, -1.0 / c, c, 0.0], lab);
    }
    for _ in 0..40 {
        let f0 = rng.log_uniform(-2.0, 0.0);
        let lens = AbcdMatrix::thin_lens(f0).unwrap();
        let qq = lens.propagate_beam([0.0, lab[2]]).unwrap_or([0.0, 0.0]);
        let m = AbcdMatrix::free_space(-qq[0]).compose(&lens);
        transform_rows(r, [m.a, m.b, m.c, m.d], lab);
    }
    for _ in 0..120 {
        let mut m = AbcdMatrix::identity();
        for _ in 0..rng.int(1, 4) {
            let e = match rng.int(0, 2) {
                0 => AbcdMatrix::free_space(rng.log_uniform(-3.0, 0.5)),
                1 => AbcdMatrix::thin_lens(rng.signed_log_uniform(-2.0, 0.0)).unwrap(),
                _ => AbcdMatrix::flat_interface(1.0, rng.uniform(1.0, 2.5)).unwrap(),
            };
            m = e.compose(&m);
        }
        transform_rows(r, [m.a, m.b, m.c, m.d], lab);
    }
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
