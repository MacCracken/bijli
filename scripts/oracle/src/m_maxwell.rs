//! `maxwell.tsv` — rust-old/src/maxwell.rs (M1 part 3d): the divergences, the displacement current, wave speed,
//! refractive index, impedance, wavelength and frequency, and skin depth.
//!
//! Every row is class A: + − × ÷ and sqrt in Rust's order, which the port repeats. A `Result` is its value on `Ok`,
//! or `E:<Variant>` and the error's Display text on `Err` (ADR 0002 row 4). A vector is three `f:` tokens.
//!
//! The grids: material and wave values (the specials, zero of both signs, negatives, ε₀, μ₀, c, copper's σ, every
//! 131st binade, seeded SI-scale values) crossed in pairs; skin depth's triples crossed over a shorter list (the
//! paired negatives, which a product alone would let through, among them) and seeded conductor cases.

use crate::grids::{Rng, pow2, specials};
use crate::m_field::{res_f, v3, vectors};
use crate::tsv::{Class, Module, f};
use bijli::field::{EPSILON_0, MU_0, SPEED_OF_LIGHT};
use bijli::maxwell::{
    displacement_current, free_space_impedance, frequency, gauss_electric_divergence,
    gauss_magnetic_divergence, impedance, refractive_index, skin_depth, wave_speed, wavelength,
};
use std::f64::consts::PI;
use std::hint::black_box;

/// Material and wave values: the specials, edges, SI constants, binades, seeded.
fn scalars(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([
        -2.0,
        -1e-300,
        1e-300,
        0.5,
        2.0,
        2.25,
        80.0,
        EPSILON_0,
        MU_0,
        SPEED_OF_LIGHT,
        5.96e7,
        1e300,
    ]);
    for k in (-1074..=1023).step_by(131) {
        xs.push(pow2(k));
    }
    for _ in 0..16 {
        xs.push(rng.signed_log_uniform(-15.0, 15.0));
    }
    xs
}

pub fn emit() -> Module {
    let mut m = Module::new("maxwell");
    let mut rng = Rng::new(0x3a8_3e11);
    let xs = scalars(&mut rng);
    for &x in &xs {
        m.row(
            "maxwell::gauss_electric_divergence",
            Class::A,
            &[f(x)],
            &[f(gauss_electric_divergence(black_box(x)))],
        );
    }
    m.row(
        "maxwell::gauss_magnetic_divergence",
        Class::A,
        &[],
        &[f(black_box(gauss_magnetic_divergence()))],
    );
    m.row(
        "maxwell::free_space_impedance",
        Class::A,
        &[],
        &[f(black_box(free_space_impedance()))],
    );
    for v in vectors() {
        m.row(
            "maxwell::displacement_current",
            Class::A,
            &v3(v),
            &v3(displacement_current(&black_box(v))),
        );
    }
    for &a in &xs {
        for &b in &xs {
            let ins = [f(a), f(b)];
            m.row(
                "maxwell::wave_speed",
                Class::A,
                &ins,
                &res_f(wave_speed(black_box(a), b)),
            );
            m.row(
                "maxwell::refractive_index",
                Class::A,
                &ins,
                &res_f(refractive_index(black_box(a), b)),
            );
            m.row(
                "maxwell::impedance",
                Class::A,
                &ins,
                &res_f(impedance(black_box(a), b)),
            );
            m.row(
                "maxwell::wavelength",
                Class::A,
                &ins,
                &res_f(wavelength(black_box(a), b)),
            );
            m.row(
                "maxwell::frequency",
                Class::A,
                &ins,
                &res_f(frequency(black_box(a), b)),
            );
        }
    }
    let mut ts = specials();
    ts.extend([
        2.0 * PI * 60.0,
        2.0 * PI * 1e6,
        MU_0,
        -MU_0,
        5.96e7,
        -5.96e7,
    ]);
    for &w in &ts {
        for &mu in &ts {
            for &sigma in &ts {
                m.row(
                    "maxwell::skin_depth",
                    Class::A,
                    &[f(w), f(mu), f(sigma)],
                    &res_f(skin_depth(black_box(w), mu, sigma)),
                );
            }
        }
    }
    for _ in 0..300 {
        let w = rng.log_uniform(0.0, 12.0);
        let mu = MU_0 * rng.log_uniform(0.0, 4.0);
        let sigma = rng.log_uniform(-2.0, 8.0);
        m.row(
            "maxwell::skin_depth",
            Class::A,
            &[f(w), f(mu), f(sigma)],
            &res_f(skin_depth(black_box(w), mu, sigma)),
        );
    }
    m
}
