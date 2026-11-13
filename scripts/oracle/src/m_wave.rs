//! `wave.tsv` — rust-old/src/wave.rs (M2 part 1b): the Poynting vector and the plane-wave relations, Snell's law and
//! the angle-form Fresnel coefficients, the trig-free Fresnel family (`*_direct`, unpolarized, Schlick), the
//! material-based interface, the waveguides and the antennas.
//!
//! Classes: the functions that call `sin`, `cos`, `asin` or `atan` are class C (their primitives are glibc's here, and
//! `validate.py` writes the correctly rounded readings to `wave.cr.tsv`); the rest are class A (+ − × ÷ and sqrt in
//! Rust's order). A `Result` is its value, or `E:<Variant>` and the error's Display text (ADR 0002 row 4); a
//! `Result<Option<f64>>` writes `n:` for `None`. A vector is three `f:` tokens, a `Material` its three fields, a `u32`
//! mode index an `i:` token.
//!
//! The grids: wave values (the specials, zero of both signs, the thresholds 1e-30 and 1 with their neighbours, SI
//! values, every 131st binade, seeded) crossed in pairs; angles over [0, π] with the specials, huge arguments and the
//! dipole's 1e-15 axis cut; the trig-free family over valid index pairs crossed with cosines through the clamp edges
//! (−1e-6, 1 + 1e-6) and each pair's critical cosine approached ulp by ulp, and invalid indices on a short list; the
//! materials of material.tsv's reflectance pairs; waveguide dimensions over the binades and the u32 mode edge; and the
//! `{ratio:.4}` of `guide_wavelength`'s error over exact ties (odd multiples of 1/32) and 300-digit ratios.

use crate::grids::{Edge2, Rng, knife_edge2, pow2, specials};
use crate::m_field::{cat, edge_vectors, res_f, v3, vectors};
use crate::m_material::{dedup, ulps};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::field::{FieldVector, SPEED_OF_LIGHT};
use bijli::material::Material;
use bijli::wave::{
    angular_frequency, b_from_e, brewster_angle, critical_angle, cylindrical_waveguide_cutoff, e_from_b,
    effective_aperture, fresnel_rp, fresnel_rp_direct, fresnel_rs, fresnel_rs_direct, friis_transmission,
    guide_wavelength, half_wave_dipole_directivity, half_wave_dipole_pattern, half_wave_dipole_radiation_resistance,
    hertzian_dipole_directivity, hertzian_dipole_pattern, hertzian_dipole_radiation_resistance,
    material_reflection_coefficient, momentum_density, plane_wave_e, plane_wave_intensity, poynting_vector,
    radiation_pressure_absorbed, radiation_pressure_reflected, rectangular_waveguide_cutoff,
    rectangular_waveguide_cutoff_wavelength, reflectance_at_interface, reflectance_normal, reflectance_normal_materials,
    reflectance_p, reflectance_p_direct, reflectance_s, reflectance_s_direct, reflectance_unpolarized,
    schlick_reflectance, snell_cos_theta_t, snell_refraction_angle, transmittance_normal, transmittance_unpolarized,
    wave_number,
};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
use std::hint::black_box;

/// A `Result<Option<f64>>`'s tokens: the value, `n:` for `None`, or the error.
fn res_o(r: bijli::error::Result<Option<f64>>) -> Vec<Tok> {
    match r {
        Ok(Some(x)) => vec![f(x)],
        Ok(None) => vec![Tok::None],
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// The wave values: specials, the thresholds with their neighbours, SI values, binades, seeded.
fn scalars(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 1.5, 2.0, 2.25, 4.0, 60.0, 1000.0, 0.03, 0.05, 5e14, 1e-300, 1e300, SPEED_OF_LIGHT, 3e8]);
    for t in [1e-30, 1.0] {
        for k in -1..=1 {
            xs.push(ulps(t, k));
            xs.push(-ulps(t, k));
        }
    }
    for k in (-1074..=1023).step_by(131) {
        xs.push(pow2(k));
    }
    for _ in 0..12 {
        xs.push(rng.signed_log_uniform(-15.0, 15.0));
    }
    dedup(xs)
}

/// The values the two-argument functions cross in pairs: the specials, the thresholds' neighbours, a few binades.
fn pairs(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 1.5, 2.25, 0.03, 0.05, 1e-300, 1e300, SPEED_OF_LIGHT]);
    for t in [1e-30, 1.0] {
        for k in -1..=1 {
            xs.push(ulps(t, k));
            xs.push(-ulps(t, k));
        }
    }
    for k in (-1074..=1023).step_by(300) {
        xs.push(pow2(k));
    }
    for _ in 0..4 {
        xs.push(rng.signed_log_uniform(-15.0, 15.0));
    }
    dedup(xs)
}

/// A shorter list, for triples and quadruples.
fn short() -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 2.25, 4.0, 1e-300, 1e300]);
    xs
}

/// Angles: the specials, the round angles, the axis cut, huge arguments, seeded.
fn angles(rng: &mut Rng) -> Vec<f64> {
    let mut ts = specials();
    ts.extend([
        PI / 6.0,
        FRAC_PI_4,
        PI / 3.0,
        FRAC_PI_2,
        PI,
        2.0 * PI,
        -FRAC_PI_4,
        0.3,
        0.5,
        0.76,
        1.2,
        1.5,
        3.0,
        100.0,
        1e-8,
        1e-15,
        1e22,
        1e300,
    ]);
    for _ in 0..20 {
        ts.push(rng.uniform(0.0, FRAC_PI_2));
    }
    dedup(ts)
}

/// A `Material`'s tokens: its three fields.
fn mt(m: Material) -> Vec<Tok> {
    vec![f(m.eps_r), f(m.mu_r), f(m.conductivity)]
}

/// The trig-free family's eight functions on one (n₁, n₂, cos θ_i).
fn trig_free(m: &mut Module, n1: f64, n2: f64, c: f64) {
    let ins = [f(n1), f(n2), f(c)];
    m.row("wave::snell_cos_theta_t", Class::A, &ins, &res_o(snell_cos_theta_t(black_box(n1), n2, c)));
    m.row("wave::fresnel_rs_direct", Class::A, &ins, &res_o(fresnel_rs_direct(black_box(n1), n2, c)));
    m.row("wave::fresnel_rp_direct", Class::A, &ins, &res_o(fresnel_rp_direct(black_box(n1), n2, c)));
    m.row("wave::reflectance_s_direct", Class::A, &ins, &res_f(reflectance_s_direct(black_box(n1), n2, c)));
    m.row("wave::reflectance_p_direct", Class::A, &ins, &res_f(reflectance_p_direct(black_box(n1), n2, c)));
    m.row("wave::reflectance_unpolarized", Class::A, &ins, &res_f(reflectance_unpolarized(black_box(n1), n2, c)));
    m.row("wave::transmittance_unpolarized", Class::A, &ins, &res_f(transmittance_unpolarized(black_box(n1), n2, c)));
    m.row("wave::schlick_reflectance", Class::A, &ins, &res_f(schlick_reflectance(black_box(n1), n2, c)));
}

/// The four angle-form Fresnel functions on one (n₁, θ_i, n₂, θ_t).
fn fresnel_angles(m: &mut Module, n1: f64, ti: f64, n2: f64, tt: f64) {
    let ins = [f(n1), f(ti), f(n2), f(tt)];
    m.row("wave::fresnel_rs", Class::C, &ins, &res_f(fresnel_rs(black_box(n1), ti, n2, tt)));
    m.row("wave::fresnel_rp", Class::C, &ins, &res_f(fresnel_rp(black_box(n1), ti, n2, tt)));
    m.row("wave::reflectance_s", Class::C, &ins, &res_f(reflectance_s(black_box(n1), ti, n2, tt)));
    m.row("wave::reflectance_p", Class::C, &ins, &res_f(reflectance_p(black_box(n1), ti, n2, tt)));
}

fn emit_plane(m: &mut Module, rng: &mut Rng, xs: &[f64], ps: &[f64]) {
    for &x in xs {
        let one = [f(x)];
        m.row("wave::plane_wave_intensity", Class::A, &one, &[f(plane_wave_intensity(black_box(x)))]);
        m.row("wave::radiation_pressure_absorbed", Class::A, &one, &[f(radiation_pressure_absorbed(black_box(x)))]);
        m.row("wave::radiation_pressure_reflected", Class::A, &one, &[f(radiation_pressure_reflected(black_box(x)))]);
        m.row("wave::e_from_b", Class::A, &one, &[f(e_from_b(black_box(x)))]);
        m.row("wave::b_from_e", Class::A, &one, &[f(b_from_e(black_box(x)))]);
        m.row("wave::angular_frequency", Class::A, &one, &[f(angular_frequency(black_box(x)))]);
    }
    for &a in ps {
        for &bb in ps {
            m.row("wave::wave_number", Class::A, &[f(a), f(bb)], &res_f(wave_number(black_box(a), bb)));
        }
    }
    // E(x, t) = E₀ sin(kx − ωt + φ): seeded physical calls, each argument through the specials, and phases past 2^60.
    let base = [100.0, 1.0, 0.5, 1.0, 0.2, 0.1];
    let mut calls: Vec<[f64; 6]> = Vec::new();
    for slot in 0..6 {
        for &x in &specials() {
            let mut c = base;
            c[slot] = x;
            calls.push(c);
        }
    }
    for _ in 0..400 {
        calls.push([
            rng.signed_log_uniform(-3.0, 6.0),
            rng.log_uniform(-2.0, 8.0),
            rng.uniform(-10.0, 10.0),
            rng.log_uniform(6.0, 16.0),
            rng.uniform(0.0, 1e-6),
            rng.uniform(-PI, PI),
        ]);
    }
    for big in [1e20, pow2(60), 1e300] {
        calls.push([1.0, big, 1.0, 0.0, 0.0, 0.0]);
        calls.push([1.0, 1.0, 0.0, big, 1.0, 0.25]);
    }
    for c in calls {
        let ins: Vec<Tok> = c.iter().map(|&x| f(x)).collect();
        m.row("wave::plane_wave_e", Class::C, &ins, &[f(plane_wave_e(black_box(c[0]), c[1], c[2], c[3], c[4], c[5]))]);
    }
    // The Poynting vector and the momentum density over pairs of vectors, and the edge vectors.
    let mut vs: Vec<FieldVector> = vectors().into_iter().step_by(43).collect();
    vs.extend(edge_vectors());
    for &e in &vs {
        for &bv in &vs {
            let ins = cat(&[v3(e), v3(bv)]);
            m.row("wave::poynting_vector", Class::A, &ins, &v3(poynting_vector(black_box(&e), &bv)));
            m.row("wave::momentum_density", Class::A, &ins, &v3(momentum_density(black_box(&e), &bv)));
        }
    }
}

fn emit_angles(m: &mut Module, rng: &mut Rng, xs: &[f64]) {
    let ts = angles(rng);
    // Snell's law, the critical and Brewster angles over index pairs (the angle forms take any index).
    let ns = [1.0, 1.5, 2.42, 0.0, -0.0, -1.5, f64::NAN, f64::INFINITY, 1e-30, ulps(1e-30, -1), 1e300];
    let ts_short: Vec<f64> = ts.iter().step_by(2).copied().collect();
    for &n1 in &ns {
        for &n2 in &ns {
            for &t in &ts_short {
                let ins = [f(n1), f(n2), f(t)];
                m.row("wave::snell_refraction_angle", Class::C, &ins, &res_o(snell_refraction_angle(black_box(n1), n2, t)));
            }
        }
    }
    // |sin θ_t| against 1, ulp by ulp around the critical angle of glass → air.
    let theta_c = (1.0_f64 / 1.5).asin();
    for k in -4..=4 {
        let t = ulps(theta_c, k);
        let ins = [f(1.5), f(1.0), f(t)];
        m.row("wave::snell_refraction_angle", Class::C, &ins, &res_o(snell_refraction_angle(black_box(1.5), 1.0, t)));
    }
    for &a in xs {
        for &bb in xs {
            let two = [f(a), f(bb)];
            m.row("wave::critical_angle", Class::C, &two, &res_o(critical_angle(black_box(a), bb)));
            m.row("wave::brewster_angle", Class::C, &two, &res_f(brewster_angle(black_box(a), bb)));
        }
    }
    // The ratio n₂/n₁ against 1, ulp by ulp.
    for k in -3..=3 {
        let n2 = ulps(1.5, k);
        m.row("wave::critical_angle", Class::C, &[f(1.5), f(n2)], &res_o(critical_angle(black_box(1.5), n2)));
    }

    // The angle-form Fresnel coefficients: index pairs × angle pairs, Snell-consistent physical calls, and the
    // |den| < 1e-30 cut (normal incidence, n₁ + n₂ about 1e-30).
    let nf = [1.0, 1.5, 0.0, -1.0, f64::NAN, 1e-31];
    let tf = [0.0, PI / 6.0, FRAC_PI_4, FRAC_PI_2, PI, f64::NAN, 1e300];
    for &n1 in &nf {
        for &n2 in &nf {
            for &ti in &tf {
                for &tt in &tf {
                    fresnel_angles(m, n1, ti, n2, tt);
                }
            }
        }
    }
    for _ in 0..150 {
        let n1 = rng.uniform(1.0, 3.0);
        let n2 = rng.uniform(1.0, 3.0);
        let ti = rng.uniform(0.0, FRAC_PI_2);
        if let Ok(Some(tt)) = snell_refraction_angle(n1, n2, ti) {
            fresnel_angles(m, n1, ti, n2, tt);
        }
    }
    for k in -3..=3 {
        let n = ulps(5e-31, k);
        fresnel_angles(m, n, 0.0, n, 0.0);
        fresnel_angles(m, n, 0.0, -ulps(1e-30, k) - n, 0.0);
    }

    // The dipole patterns over [0, π] in 720 steps, the angles, the 1e-15 axis cut ulp by ulp at 0 and π, seeded.
    let mut th: Vec<f64> = (0..=720).map(|j| f64::from(j) * PI / 720.0).collect();
    th.extend(ts.iter().copied());
    for k in -3..=3 {
        th.push(ulps(1e-15, k));
        th.push(-ulps(1e-15, k));
        th.push(PI - ulps(1e-15, k));
    }
    th.extend([pow2(-30), pow2(-60), pow2(-1074), f64::MAX, -PI, 4.0 * PI]);
    for _ in 0..200 {
        th.push(rng.uniform(-10.0, 10.0));
    }
    for &t in &dedup(th) {
        m.row("wave::hertzian_dipole_pattern", Class::C, &[f(t)], &[f(hertzian_dipole_pattern(black_box(t)))]);
        m.row("wave::half_wave_dipole_pattern", Class::C, &[f(t)], &[f(half_wave_dipole_pattern(black_box(t)))]);
    }
}

fn emit_trig_free(m: &mut Module, rng: &mut Rng) {
    let good = [1.0, 1.5, 2.42, 1.0 + f64::EPSILON, 1.0 - f64::EPSILON / 2.0, 1e-31, f64::from_bits(1), 1e154, f64::MAX];
    let bad = [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    let mut cs = specials();
    cs.extend((0..=8).map(|j| f64::from(j) / 8.0));
    cs.extend([-1e-6, ulps(1e-6, 1), ulps(1e-6, -1)].map(|x| if x > 0.0 { -x } else { x }));
    cs.extend([1.0 + 1e-6, ulps(1.0 + 1e-6, 1), ulps(1.0 + 1e-6, -1)]);
    cs.extend([-1e-9, 1.0 + 1e-9, 0.3, 0.76, 0.99, FRAC_PI_4.cos(), 1e-300, 1.0 - f64::EPSILON / 2.0, 2.0]);
    let cs = dedup(cs);
    for &n1 in &good {
        for &n2 in &good {
            for &c in &cs {
                trig_free(m, n1, n2, c);
            }
        }
    }
    let all: Vec<f64> = good.iter().chain(bad.iter()).copied().collect();
    for &n1 in &all {
        for &n2 in &all {
            if bad.iter().any(|b| b.to_bits() == n1.to_bits()) || bad.iter().any(|b| b.to_bits() == n2.to_bits()) {
                for c in [0.5, 1.0, f64::NAN, -1.0] {
                    trig_free(m, n1, n2, c);
                }
            }
        }
    }
    // Each dense → rare pair's critical cosine √(1 − (n₂/n₁)²), ulp by ulp: sin²θ_t against 1.
    for (n1, n2) in [(1.5, 1.0), (2.42, 1.0), (1.33, 1.0), (3.5, 1.5), (1.0 + f64::EPSILON, 1.0), (1e200, 1e-100)] {
        let r: f64 = n2 / n1;
        let cc = (1.0 - r * r).sqrt();
        for k in -4..=4 {
            trig_free(m, n1, n2, ulps(cc, k));
        }
    }
    for _ in 0..300 {
        let n1 = rng.uniform(1.0, 3.0);
        let n2 = rng.uniform(1.0, 3.0);
        trig_free(m, n1, n2, rng.unit());
    }
    // M2 part 2c-3 (ADR 0002 row 17: TIR decided exactly, cos θ_t rounded once): the critical cosine ±8 ulps for 24
    // more dense → rare pairs (ratios of 1 to 3, within 2^-50 of 1, and up to 2^26; some scaled by 2^±600); n₁ within a
    // few ulps of n₂ with small and tiny cos θ_i (Rust's 1 − (1 − c²) cancels to 0); c = 1 with a ratio past DBL_MAX
    // (Rust's ∞·0 = NaN); and ratios at the shortcuts' exponent differences. A seed of its own keeps the rest as it was.
    let mut trng = Rng::new(0x0C05_72C3);
    for j in 0..24 {
        let n2 = trng.uniform(1.0, 3.0);
        let q = match j % 3 {
            0 => trng.uniform(1.0, 3.0),
            1 => 1.0 + pow2(-i32::try_from(trng.int(1, 50)).unwrap_or(1)),
            _ => pow2(i32::try_from(trng.int(1, 25)).unwrap_or(1)) * trng.uniform(1.0, 2.0),
        };
        let n1 = n2 * q;
        let sc = match j % 4 {
            2 => pow2(-600),
            3 => pow2(600),
            _ => 1.0,
        };
        let r = n2 / n1;
        let cc = (1.0 - r * r).sqrt();
        for k in (-8..=8).filter(|_| cc > 0.0) {
            trig_free(m, n1 * sc, n2 * sc, ulps(cc, k).min(1.0));
        }
    }
    for k in [-3, -1, 1, 3] {
        for c in [pow2(-10), pow2(-27), pow2(-40), pow2(-200), pow2(-700), 1e-9, 1e-3] {
            trig_free(m, ulps(1.5, k), 1.5, c);
        }
    }
    for (n1, n2) in [(1e300, 1e-300), (f64::MAX, f64::from_bits(1)), (1e-300, 1e300)] {
        for c in [1.0, 1.0 - f64::EPSILON / 2.0, 0.5] {
            trig_free(m, n1, n2, c);
        }
    }
    for d in [26, 27, 28, -60, -61, -62] {
        for c in [1.0 - f64::EPSILON / 2.0, 1.0 - pow2(-40), 0.5] {
            trig_free(m, 1.25 * pow2(d), 1.0, c);
        }
    }
}

fn emit_materials(m: &mut Module, rng: &mut Rng) {
    let mut rs = vec![
        Material::vacuum(),
        Material::dielectric(2.25),
        Material::dielectric(4.0),
        Material::new(4.0, 4.0, 0.0),
        Material::new(16.0, 1.0, 0.0),
        Material::new(1.0, 2.25, 0.0),
        Material::new(80.0, 1.0, 0.0),
        Material::conductor(5.96e7),
        Material::new(2.0, 1.0, f64::NAN),
        Material::new(0.0, 1.0, 0.0),
        Material::new(-2.0, 1.0, 0.0),
        Material::new(f64::INFINITY, 1.0, 0.0),
        Material::new(1.0, -1.0, 0.0),
        Material::new(1.0, f64::NAN, 0.0),
        Material::new(f64::MAX, f64::MAX, 0.0),
        Material::new(1e300, 1e-300, 0.0),
        Material::new(1e-200, 1e-200, 0.0),
        Material::new(1e200, 1e-200, 0.0),
    ];
    for _ in 0..4 {
        rs.push(Material::new(rng.log_uniform(0.0, 2.0), rng.log_uniform(0.0, 1.0), 0.0));
    }
    let cs = [0.0, 0.3, 0.8, 1.0, -1e-9, -0.5, FRAC_PI_4.cos()];
    for &a in &rs {
        for &bm in &rs {
            let two = cat(&[mt(a), mt(bm)]);
            m.row(
                "wave::reflectance_normal_materials",
                Class::A,
                &two,
                &res_f(reflectance_normal_materials(black_box(&a), &bm)),
            );
            m.row(
                "wave::material_reflection_coefficient",
                Class::A,
                &two,
                &res_f(material_reflection_coefficient(black_box(&a), &bm)),
            );
            for &c in &cs {
                let ins = cat(&[mt(a), mt(bm), vec![f(c)]]);
                m.row("wave::reflectance_at_interface", Class::A, &ins, &res_f(reflectance_at_interface(black_box(&a), &bm, c)));
            }
        }
    }
    // M2 part 2c-3: each dense → rare dielectric pair's critical cosine ±8 ulps, through the interface.
    for (e1, e2) in [(2.25, 1.0), (4.0, 1.0), (1.7, 1.1), (5.3, 2.25), (1.0 + 4.0 * f64::EPSILON, 1.0), (80.0, 1.0)] {
        let (a, bm) = (Material::dielectric(e1), Material::dielectric(e2));
        let r = e2.sqrt() / e1.sqrt();
        let cc = (1.0 - r * r).sqrt();
        for k in (-8..=8).filter(|_| cc > 0.0) {
            let c = ulps(cc, k).min(1.0);
            let ins = cat(&[mt(a), mt(bm), vec![f(c)]]);
            m.row("wave::reflectance_at_interface", Class::A, &ins, &res_f(reflectance_at_interface(black_box(&a), &bm, c)));
        }
    }
}

fn emit_normal_antennas(m: &mut Module, xs: &[f64], sh: &[f64]) {
    for &a in xs {
        for &bb in xs {
            let two = [f(a), f(bb)];
            m.row("wave::reflectance_normal", Class::A, &two, &res_f(reflectance_normal(black_box(a), bb)));
            m.row("wave::transmittance_normal", Class::A, &two, &res_f(transmittance_normal(black_box(a), bb)));
            m.row(
                "wave::hertzian_dipole_radiation_resistance",
                Class::A,
                &two,
                &res_f(hertzian_dipole_radiation_resistance(black_box(a), bb)),
            );
            m.row("wave::effective_aperture", Class::A, &two, &res_f(effective_aperture(black_box(a), bb)));
        }
    }
    // The |n₁ + n₂| < 1e-30 cut, ulp by ulp.
    for k in -3..=3 {
        let n = ulps(5e-31, k);
        for (a, bb) in [(n, n), (-n, -n), (n, ulps(5e-31, -k))] {
            let two = [f(a), f(bb)];
            m.row("wave::reflectance_normal", Class::A, &two, &res_f(reflectance_normal(black_box(a), bb)));
            m.row("wave::transmittance_normal", Class::A, &two, &res_f(transmittance_normal(black_box(a), bb)));
        }
    }
    // Its knife edge (bijli ADR 0002 row 17's rule, M2 part 2c-1): n₁ + n₂ within an ulp of ±10^-30, the sum exact
    // (two halves) or rounded (a tiny n₂ against n₁ near 1e-30, or both near 1e-20 with cancellation).
    for (a, bb) in knife_edge2(Edge2::Sum, &[5e-31, 3e-31, 1e-30, 1e-20, 1e-31], 6, 24) {
        let two = [f(a), f(bb)];
        m.row("wave::reflectance_normal", Class::A, &two, &res_f(reflectance_normal(black_box(a), bb)));
        m.row("wave::transmittance_normal", Class::A, &two, &res_f(transmittance_normal(black_box(a), bb)));
    }
    m.row("wave::hertzian_dipole_directivity", Class::A, &[], &[f(black_box(hertzian_dipole_directivity()))]);
    m.row("wave::half_wave_dipole_directivity", Class::A, &[], &[f(black_box(half_wave_dipole_directivity()))]);
    m.row(
        "wave::half_wave_dipole_radiation_resistance",
        Class::A,
        &[],
        &[f(black_box(half_wave_dipole_radiation_resistance()))],
    );
    // Friis: gain pairs × (λ, d) over the short list, and the distance cut.
    let gains = [(1.0, 1.0), (1.5, 1.64), (-1.0, 1.0), (f64::NAN, 1.0), (1e300, 1e300), (1e-300, 1e-300)];
    for (gt, gr) in gains {
        for &l in sh {
            for &d in sh {
                let ins = [f(gt), f(gr), f(l), f(d)];
                m.row("wave::friis_transmission", Class::A, &ins, &res_f(friis_transmission(black_box(gt), gr, l, d)));
            }
        }
    }
    for k in -1..=1 {
        for d in [ulps(1e-30, k), -ulps(1e-30, k)] {
            let ins = [f(1.0), f(1.0), f(1.0), f(d)];
            m.row("wave::friis_transmission", Class::A, &ins, &res_f(friis_transmission(black_box(1.0), 1.0, 1.0, d)));
        }
    }
}

fn emit_waveguides(m: &mut Module, rng: &mut Rng, xs: &[f64], sh: &[f64]) {
    let dims = [22.86e-3, 10.16e-3, 0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300, 1e300, f64::MAX, 1e-160];
    let modes: [(u32, u32); 7] = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (3, 7), (u32::MAX, 0)];
    let vels = [SPEED_OF_LIGHT, 0.0, f64::NAN, f64::INFINITY, 1e300];
    for &a in &dims {
        for &bb in &dims {
            for (mm, nn) in modes {
                let four = [f(a), f(bb), i(mm), i(nn)];
                m.row(
                    "wave::rectangular_waveguide_cutoff_wavelength",
                    Class::A,
                    &four,
                    &res_f(rectangular_waveguide_cutoff_wavelength(black_box(a), bb, mm, nn)),
                );
                for &v in &vels {
                    let five = [f(a), f(bb), i(mm), i(nn), f(v)];
                    m.row(
                        "wave::rectangular_waveguide_cutoff",
                        Class::A,
                        &five,
                        &res_f(rectangular_waveguide_cutoff(black_box(a), bb, mm, nn, v)),
                    );
                }
            }
        }
    }
    // guide_wavelength: pairs, the ratio against 1 ulp by ulp, ties of its {ratio:.4} (odd multiples of 1/32), long
    // ratios, and the TEM line (λ_c = ∞).
    for &a in xs {
        for &bb in xs {
            m.row("wave::guide_wavelength", Class::A, &[f(a), f(bb)], &res_f(guide_wavelength(black_box(a), bb)));
        }
    }
    let mut gw: Vec<(f64, f64)> = Vec::new();
    for k in -3..=3 {
        gw.push((ulps(0.03, k), 0.03));
        gw.push((0.03, ulps(0.03, k)));
    }
    for j in [33, 35, 37, 63, 65, 97, 161, 3201] {
        gw.push((f64::from(j) / 32.0, 1.0));
    }
    for (l, lc) in [(1e300, 1e-5), (f64::MAX, 1.0), (1.0, f64::from_bits(1)), (123_456.789_05, 1.0), (9.99995, 1.0)] {
        gw.push((l, lc));
    }
    for _ in 0..100 {
        gw.push((rng.uniform(0.01, 0.06), rng.uniform(0.02, 0.05)));
    }
    for (l, lc) in gw {
        m.row("wave::guide_wavelength", Class::A, &[f(l), f(lc)], &res_f(guide_wavelength(black_box(l), lc)));
    }
    // The circular waveguide: radius × Bessel root × velocity.
    let roots = [1.8412, 2.4048, 0.0, -1.8412, f64::NAN, f64::INFINITY, 1e300, 1e-300];
    let rvels = [SPEED_OF_LIGHT, 0.0, f64::NAN, f64::INFINITY, 1e300, 1e-300];
    let mut radii: Vec<f64> = sh.to_vec();
    radii.extend([10e-3, 22.86e-3]);
    for &r in &radii {
        for &x in &roots {
            for &v in &rvels {
                let ins = [f(r), f(x), f(v)];
                m.row(
                    "wave::cylindrical_waveguide_cutoff",
                    Class::A,
                    &ins,
                    &res_f(cylindrical_waveguide_cutoff(black_box(r), x, v)),
                );
            }
        }
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("wave");
    let mut rng = Rng::new(0x3a7_e1b);
    let xs = scalars(&mut rng);
    let ps = pairs(&mut rng);
    let sh = short();
    emit_plane(&mut m, &mut rng, &xs, &ps);
    emit_angles(&mut m, &mut rng, &ps);
    emit_trig_free(&mut m, &mut rng);
    emit_materials(&mut m, &mut rng);
    emit_normal_antennas(&mut m, &ps, &sh);
    emit_waveguides(&mut m, &mut rng, &ps, &sh);
    emit_eta_guard(&mut m);
    emit_ratio_edges(&mut m);
    emit_flux_edges(&mut m);
    m
}

/// M2 part 2d-1 (bijli ADR 0002 row 25 amended): `material_reflection_coefficient` between materials of tiny relative
/// impedance η = √(μ_r/ε_r) (ε_r = 1, μ_r = η²), across Rust's cut (its absolute sum 376.7·(η₁ + η₂) < 1e-30) and the
/// port's former one (η₁ + η₂ < 1e-30). After every older row.
fn emit_eta_guard(m: &mut Module) {
    let etas = [1e-34, 1e-33, 3e-33, 1e-31, 4e-31, 6e-31, 1e-29];
    for &e1 in &etas {
        for &e2 in &etas[2..] {
            let a = Material::new(1.0, e1 * e1, 0.0);
            let bm = Material::new(1.0, e2 * e2, 0.0);
            let two = cat(&[mt(a), mt(bm)]);
            m.row(
                "wave::material_reflection_coefficient",
                Class::A,
                &two,
                &res_f(material_reflection_coefficient(black_box(&a), &bm)),
            );
        }
    }
}


/// M2 part 2d-2a (bijli ADR 0002 row 18's rule): the ratio chains (a·x − b·y)/(a·x + b·y) at their range edges, where
/// a sum passes DBL_MAX or a product falls below 2^-1022: `reflectance_normal` and `transmittance_normal`, the
/// trig-free family and Schlick, the angle forms. Indices at DBL_MAX and its fractions, at 10^300, subnormal and
/// tiny; cosines subnormal, tiny and ordinary. After every older row.
fn emit_ratio_edges(m: &mut Module) {
    let big = [f64::MAX, f64::MAX / 2.0, f64::MAX / 3.0, 1e308, 1.5e308, 9e307];
    let small = [5e-324, 1e-323, 2.5e-322, 1e-310, 3e-310, 2.2250738585072014e-308, 1e-300, 3e-300];
    for &a in &big {
        for &b in &big {
            let two = [f(a), f(b)];
            m.row("wave::reflectance_normal", Class::A, &two, &res_f(reflectance_normal(black_box(a), b)));
            m.row("wave::transmittance_normal", Class::A, &two, &res_f(transmittance_normal(black_box(a), b)));
            let nb = -b * 0.5;
            let mixed = [f(a), f(nb)];
            m.row("wave::reflectance_normal", Class::A, &mixed, &res_f(reflectance_normal(black_box(a), nb)));
            for c in [1.0, 0.5, 1e-310] {
                trig_free(m, a, b, c);
            }
        }
    }
    for &a in &small {
        for &b in &small {
            let two = [f(a), f(b)];
            m.row("wave::reflectance_normal", Class::A, &two, &res_f(reflectance_normal(black_box(a), b)));
            for c in [1.0, 0.3, 1e-20] {
                trig_free(m, a, b, c);
            }
        }
    }
    for (a, b) in [(1.0, 1.5), (1.5, 1.0), (3.0, 1.0), (1e-300, 2e-300), (1e300, 7e299)] {
        for c in [5e-324, 1e-323, 1e-310, 2.2250738585072014e-308, 1e-160] {
            trig_free(m, a, b, c);
        }
    }
    for &(a, b) in &[(f64::MAX, f64::MAX / 2.0), (1e308, 1e308), (1e308, -1e308), (1e-310, 3e-310), (5e-324, 1e-323)] {
        for &(ti, tt) in &[(0.0, 0.0), (0.3, 0.5), (1.2, 0.2), (2.5, 0.4), (FRAC_PI_2, 0.1)] {
            fresnel_angles(m, a, ti, b, tt);
        }
    }
}

/// M2 part 2d-2b (bijli ADR 0002 row 18's rule): the plane waves' chains at their range edges. The Poynting vector and
/// the momentum density with a cross below 2^-1022 (S normal), with S past DBL_MAX and with the cross past it where g
/// is representable; `plane_wave_e` with kx and ωt past DBL_MAX (cancelling, or not: a sine past DBL_MAX) and below
/// 2^-1022 under a large E₀; Snell's n₁ sin θ_i below 2^-1022 over a small n₂. After every older row.
fn emit_flux_edges(m: &mut Module) {
    let pairs = [
        ([pow2(-515), 0.0, 0.0], [0.0, 3.0 * pow2(-515), 0.0]),
        ([1e-155, 2e-155, 0.0], [0.0, 3e-156, 5e-156]),
        ([1.5e-160, -2.5e-160, 1e-161], [7e-161, 3e-162, -4e-160]),
        ([pow2(-540), pow2(-541), 0.0], [pow2(-500), -pow2(-499), pow2(-520)]),
        ([1e152, 0.0, 0.0], [0.0, 1e152, 0.0]),
        ([1e159, 2e158, 0.0], [-3e159, 1e160, 7e158]),
        ([1e159, 0.0, 0.0], [0.0, 1e160, 0.0]),
        ([1e300, 0.0, 0.0], [0.0, 1e300, 0.0]),
        ([f64::MAX, 0.0, 0.0], [0.0, 1e-10, 0.0]),
    ];
    for (e, b) in pairs {
        for flip in [1.0, -1.0] {
            let ev = FieldVector::new(e[0], e[1] * flip, e[2]);
            let bv = FieldVector::new(b[0], b[1], b[2] * flip);
            let ins = cat(&[v3(ev), v3(bv)]);
            m.row("wave::poynting_vector", Class::A, &ins, &v3(poynting_vector(black_box(&ev), &bv)));
            m.row("wave::momentum_density", Class::A, &ins, &v3(momentum_density(black_box(&ev), &bv)));
        }
    }
    let waves = [
        [1.0, 1e200, 1e200, 1e200, 1e200, 0.5],
        [1.0, 1e200, 1e200, 0.0, 0.0, 0.0],
        [2.0, -1e200, 3e200, 1e100, 1e100, 0.25],
        [1.0, 1e200, 1e200, 1e200, 9e199, 0.0],
        [1e300, 1e-200, 1e-200, 0.0, 0.0, 0.0],
        [1e300, 1e-200, 3e-201, 1e-200, 1e-200, 0.0],
        [1e250, 1e-170, 1e-150, 0.0, 1.0, 0.0],
        [5.0, 1e-200, 1e-200, 2.0, 0.5, 1.0],
        [1.0, 1e308, 10.0, 1e308, 10.0, 0.0],
    ];
    for c in waves {
        let ins = [f(c[0]), f(c[1]), f(c[2]), f(c[3]), f(c[4]), f(c[5])];
        m.row("wave::plane_wave_e", Class::C, &ins, &[f(plane_wave_e(black_box(c[0]), c[1], c[2], c[3], c[4], c[5]))]);
    }
    for (n1, n2, t) in [(1e-300, 1e-29, 1e-20), (3e-300, 2e-29, 1e-15), (1e-310, 1e-28, 0.5), (-2e-305, 1e-29, 1e-4)] {
        let ins = [f(n1), f(n2), f(t)];
        m.row("wave::snell_refraction_angle", Class::C, &ins, &res_o(snell_refraction_angle(black_box(n1), n2, t)));
    }
}
