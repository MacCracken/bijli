//! `material.tsv` — rust-old/src/material.rs (M2 part 1a): the dielectric and magnetic relations, Clausius–Mossotti,
//! the Curie laws, `classify_magnetic` and the `Material` struct.
//!
//! Every row is class A: + − × ÷ and sqrt in Rust's order, which the port repeats. A `Result` is its value on `Ok`,
//! or `E:<Variant>` and the error's Display text on `Err` (ADR 0002 row 4). A vector is three `f:` tokens, a
//! `Material` its three fields in declaration order (eps_r, mu_r, conductivity), a `MagneticType` an `e:` token.
//!
//! The grids: material values (the specials, zero of both signs, negatives, ε₀, μ₀, the thresholds 1e-15, 1e-30 and
//! 1 with their neighbours, every 131st binade, seeded SI-scale values) crossed in pairs for the two-argument
//! functions; a shorter list crossed in triples for `curie_weiss`, `bound_volume_charge` and the `Material`s; scale
//! factors over a stride of field.rs's vectors; and the Clausius–Mossotti and permeability thresholds approached ulp
//! by ulp.

use crate::grids::{Edge2, Rng, knife_edge2, pow2, specials};
use crate::m_field::{cat, res_f, res_v, v3, vectors};
use crate::tsv::{Class, Module, Tok, b, e, f};
use bijli::field::{EPSILON_0, FieldVector, MU_0};
use bijli::material::{
    Material, absolute_permeability, absolute_permittivity, b_field_from_h, bound_surface_charge,
    bound_surface_current, bound_volume_charge, clausius_mossotti, classify_magnetic, curie_law, curie_weiss,
    dielectric_energy_density, displacement_field, electric_susceptibility, h_field_from_b, magnetic_energy_density_material,
    magnetic_susceptibility, magnetization, polarization, relative_permeability, relative_permittivity,
};
use std::hint::black_box;

/// `x` moved `k` ulps (k > 0 up, k < 0 down), for a positive finite x.
pub(crate) fn ulps(x: f64, k: i64) -> f64 {
    let bits = x.to_bits();
    let moved = if k >= 0 { bits + k.unsigned_abs() } else { bits - k.unsigned_abs() };
    f64::from_bits(moved)
}

/// `xs` with each value (by bits) kept once, in order: a row's (function, inputs) is its key in the mapping files.
pub(crate) fn dedup(xs: Vec<f64>) -> Vec<f64> {
    let mut seen = std::collections::BTreeSet::new();
    xs.into_iter().filter(|x| seen.insert(x.to_bits())).collect()
}

/// The material values: specials, edges and thresholds with their neighbours, SI constants, binades, seeded.
fn scalars(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 2.0, 2.25, 4.0, 80.0, 1000.0, 1.001, 0.001, 300.0, 770.0, 1e-300, 1e300, EPSILON_0, MU_0]);
    for t in [1e-15, 1e-30, 1.0] {
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

/// A shorter list, for triples.
fn short() -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 2.25, 4.0, 1e-300, 1e300]);
    xs
}

/// A `Material`'s tokens: its three fields.
fn mt(m: Material) -> Vec<Tok> {
    vec![f(m.eps_r), f(m.mu_r), f(m.conductivity)]
}

pub fn emit() -> Module {
    let mut m = Module::new("material");
    let mut rng = Rng::new(0x3a7_e41a);
    let xs = scalars(&mut rng);
    let sh = short();

    // One-argument relations and the classification.
    for &x in &xs {
        let one = [f(x)];
        m.row("material::electric_susceptibility", Class::A, &one, &[f(electric_susceptibility(black_box(x)))]);
        m.row("material::relative_permittivity", Class::A, &one, &[f(relative_permittivity(black_box(x)))]);
        m.row("material::absolute_permittivity", Class::A, &one, &[f(absolute_permittivity(black_box(x)))]);
        m.row("material::magnetic_susceptibility", Class::A, &one, &[f(magnetic_susceptibility(black_box(x)))]);
        m.row("material::relative_permeability", Class::A, &one, &[f(relative_permeability(black_box(x)))]);
        m.row("material::absolute_permeability", Class::A, &one, &[f(absolute_permeability(black_box(x)))]);
        m.row("material::classify_magnetic", Class::A, &one, &[e("MagneticType", classify_magnetic(black_box(x)))]);
    }

    // Pairs: the energy densities, Clausius–Mossotti and Curie's law.
    for &a in &xs {
        for &bb in &xs {
            let two = [f(a), f(bb)];
            m.row(
                "material::dielectric_energy_density",
                Class::A,
                &two,
                &[f(dielectric_energy_density(black_box(a), bb))],
            );
            m.row("material::clausius_mossotti", Class::A, &two, &res_f(clausius_mossotti(black_box(a), bb)));
            m.row(
                "material::magnetic_energy_density_material",
                Class::A,
                &two,
                &res_f(magnetic_energy_density_material(black_box(a), bb)),
            );
            m.row("material::curie_law", Class::A, &two, &res_f(curie_law(black_box(a), bb)));
        }
    }
    // Clausius–Mossotti at its two edges, X = Nα/(3ε₀) → 1 (the catastrophe) and X → −1/2 (ε_r → 0), ulp by ulp; and
    // seeded physical cases (α ~ 1e-40 F·m², N ~ 1e28 m⁻³).
    for n in [1e28, 2.5e25, 1.0] {
        for k in -3..=3 {
            let up = ulps(3.0 * EPSILON_0 / n, k);
            let dn = -ulps(1.5 * EPSILON_0 / n, k);
            for alpha in [up, dn] {
                m.row(
                    "material::clausius_mossotti",
                    Class::A,
                    &[f(alpha), f(n)],
                    &res_f(clausius_mossotti(black_box(alpha), n)),
                );
            }
        }
    }
    for _ in 0..200 {
        let alpha = rng.signed_log_uniform(-42.0, -38.0);
        let n = rng.log_uniform(24.0, 29.0);
        m.row("material::clausius_mossotti", Class::A, &[f(alpha), f(n)], &res_f(clausius_mossotti(black_box(alpha), n)));
    }
    // The permeability threshold |μ_r·μ₀| < 1e-30, ulp by ulp.
    for k in -3..=3 {
        for mu in [ulps(1e-30 / MU_0, k), -ulps(1e-30 / MU_0, k)] {
            m.row(
                "material::magnetic_energy_density_material",
                Class::A,
                &[f(mu), f(1.5)],
                &res_f(magnetic_energy_density_material(black_box(mu), 1.5)),
            );
        }
    }
    // Its knife edge (bijli ADR 0002 row 17's rule, M2 part 2c-1): μ_r·μ₀ within an ulp of 10^-30.
    for (_mu0, mu) in knife_edge2(Edge2::Prod, &[MU_0], 0, 400) {
        m.row(
            "material::magnetic_energy_density_material",
            Class::A,
            &[f(mu), f(1.5)],
            &res_f(magnetic_energy_density_material(black_box(mu), 1.5)),
        );
    }
    // classify_magnetic about −1e-15, ulp by ulp (the exact test χ < −10^-15).
    for k in -3..=3 {
        let x = -ulps(1e-15, k);
        m.row("material::classify_magnetic", Class::A, &[f(x)], &[e("MagneticType", classify_magnetic(black_box(x)))]);
    }

    // Triples: Curie–Weiss and the bound volume charge.
    for &c in &sh {
        for &t in &sh {
            for &tc in &sh {
                let three = [f(c), f(t), f(tc)];
                m.row("material::curie_weiss", Class::A, &three, &res_f(curie_weiss(black_box(c), t, tc)));
                m.row("material::bound_volume_charge", Class::A, &three, &[f(bound_volume_charge(black_box(c), t, tc))]);
            }
        }
    }
    for (c, t, tc) in [(1.0, 1000.0, 770.0), (1.0, 770.0, 770.0), (2.0, 300.0, 0.0), (1.0, ulps(770.0, 1), 770.0)] {
        m.row("material::curie_weiss", Class::A, &[f(c), f(t), f(tc)], &res_f(curie_weiss(black_box(c), t, tc)));
    }

    // Scale factors over a stride of field.rs's vectors.
    let vs: Vec<FieldVector> = vectors().into_iter().step_by(13).collect();
    let mut scales = specials();
    scales.extend([-2.0, 0.5, 4.0, 1000.0, 1e-30, 1e-31, 1e300, 1e-300]);
    for &k in &scales {
        for &v in &vs {
            let ins = cat(&[vec![f(k)], v3(v)]);
            m.row("material::polarization", Class::A, &ins, &v3(polarization(black_box(k), &v)));
            m.row("material::displacement_field", Class::A, &ins, &v3(displacement_field(black_box(k), &v)));
            m.row("material::magnetization", Class::A, &ins, &v3(magnetization(black_box(k), &v)));
            m.row("material::b_field_from_h", Class::A, &ins, &v3(b_field_from_h(black_box(k), &v)));
            m.row("material::h_field_from_b", Class::A, &ins, &res_v(h_field_from_b(black_box(k), &v)));
        }
    }
    // The bound surface charge and current over pairs of them.
    let ws: Vec<FieldVector> = vectors().into_iter().step_by(29).collect();
    for &p in &ws {
        for &nv in &ws {
            let ins = cat(&[v3(p), v3(nv)]);
            m.row("material::bound_surface_charge", Class::A, &ins, &[f(bound_surface_charge(black_box(&p), &nv))]);
            m.row("material::bound_surface_current", Class::A, &ins, &v3(bound_surface_current(black_box(&p), &nv)));
        }
    }

    // Material: the constructors, then every method over (ε_r, μ_r, σ) triples.
    m.row("Material::vacuum", Class::A, &[], &mt(black_box(Material::vacuum())));
    m.row("Material::default", Class::A, &[], &mt(black_box(Material::default())));
    for &x in &xs {
        m.row("Material::dielectric", Class::A, &[f(x)], &mt(Material::dielectric(black_box(x))));
        m.row("Material::conductor", Class::A, &[f(x)], &mt(Material::conductor(black_box(x))));
    }
    let mut em = sh.clone();
    em.extend([80.0, 1000.0]);
    let sigmas = [0.0, -0.0, 1.0, -1.0, 5.96e7, f64::NAN, f64::INFINITY];
    let mut mats = Vec::new();
    for &ep in &em {
        for &mu in &em {
            for &sg in &sigmas {
                mats.push(Material::new(ep, mu, sg));
            }
        }
    }
    for &mat in &mats {
        let ins = mt(mat);
        m.row("Material::new", Class::A, &ins, &mt(Material::new(black_box(mat.eps_r), mat.mu_r, mat.conductivity)));
        m.row("Material::refractive_index", Class::A, &ins, &[f(black_box(mat).refractive_index())]);
        m.row("Material::permittivity", Class::A, &ins, &[f(black_box(mat).permittivity())]);
        m.row("Material::permeability", Class::A, &ins, &[f(black_box(mat).permeability())]);
        m.row("Material::impedance", Class::A, &ins, &res_f(black_box(mat).impedance()));
        m.row("Material::wave_speed", Class::A, &ins, &res_f(black_box(mat).wave_speed()));
        m.row("Material::is_lossy", Class::A, &ins, &[b(black_box(mat).is_lossy())]);
    }
    // try_reflectance_at over pairs: lossless dielectrics and magnetic media, conductors and invalid materials.
    let mut rs = vec![
        Material::vacuum(),
        Material::dielectric(2.25),
        Material::dielectric(4.0),
        Material::new(4.0, 4.0, 0.0),
        Material::new(1.0, 1000.0, 0.0),
        Material::new(80.0, 1.0, 0.0),
        Material::new(2.25, 1.0, -0.0),
        Material::conductor(5.96e7),
        Material::new(2.0, 1.0, f64::NAN),
        Material::new(0.0, 1.0, 0.0),
        Material::new(-1.0, 1.0, 0.0),
        Material::new(f64::NAN, 1.0, 0.0),
        Material::new(f64::INFINITY, 1.0, 0.0),
        Material::new(1.0, 0.0, 0.0),
        Material::new(1.0, f64::NAN, 0.0),
        Material::new(1.0, f64::INFINITY, 0.0),
        Material::new(f64::MIN_POSITIVE, 1.0, 0.0),
        Material::new(f64::MAX, f64::MAX, 0.0),
        Material::new(1e300, 1e-300, 0.0),
    ];
    for _ in 0..12 {
        rs.push(Material::new(rng.log_uniform(0.0, 2.0), rng.log_uniform(0.0, 1.0), 0.0));
    }
    for &a in &rs {
        for &bm in &rs {
            m.row(
                "Material::try_reflectance_at",
                Class::A,
                &cat(&[mt(a), mt(bm)]),
                &res_f(black_box(a).try_reflectance_at(&bm)),
            );
        }
    }
    m
}
