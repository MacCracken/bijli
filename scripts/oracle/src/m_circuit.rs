//! `circuit.tsv` — rust-old/src/circuit.rs (M2 part 1c): Ohm's law, series and parallel combinations, capacitance
//! and inductance, the RC and RL transients and the RLC relations.
//!
//! Classes: the four transients call `exp` and are class C (`validate.py` writes the correctly rounded readings to
//! `circuit.cr.tsv`); the rest are class A (+ − × ÷ and sqrt in Rust's order). A `Result` is its value, or
//! `E:<Variant>` and the error's Display text (ADR 0002 row 4). A slice is `v:<n>` then its n values.
//!
//! The grids: circuit values (the specials, the 1e-30 threshold and 1 with their neighbours, SI-scale values, binades,
//! seeded) crossed in pairs; a short list crossed in triples and an eight-value list in quadruples; slices of length
//! 0 to 3, with the parallel combinations' 1/Σ(1/x) cut at 1e-30 approached ulp by ulp and the transients' τ = RC and
//! τ = L/R cut likewise; seeded physical circuits; and (M2 part 2c-5) the inverse combinations at Σ 1/x = 0 and at
//! the 1e-30 cut (`emit_inverse_edges`).

use crate::grids::{Edge2, Rng, knife_band2, knife_edge2, pow2, specials};
use crate::m_field::res_f;
use crate::m_material::{dedup, ulps};
use crate::tsv::{Class, Module, Tok, f};
use bijli::circuit::{
    capacitance_parallel, capacitance_series, capacitor_charge, capacitor_energy, current, damping_ratio,
    inductance_parallel, inductance_series, inductor_energy, parallel_plate_capacitance, power, quality_factor,
    rc_charging_voltage, rc_discharging_voltage, rc_time_constant, resistance, resistance_from_geometry,
    resistance_parallel, resistance_series, resonant_angular_frequency, resonant_frequency, rl_current_decay,
    rl_current_rise, rl_time_constant, rlc_impedance, solenoid_inductance, voltage,
};
use std::hint::black_box;

/// The values crossed in pairs: specials, the thresholds' neighbours, SI values, a few binades, seeded.
fn pairs(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 1.5, 100.0, 1e-6, 1e-3, 4.7e-9, 1e-300, 1e300]);
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
        xs.push(rng.signed_log_uniform(-12.0, 12.0));
    }
    dedup(xs)
}

/// A shorter list, for triples.
fn short() -> Vec<f64> {
    let mut xs = specials();
    xs.extend([-2.0, 0.5, 2.25, 1e-6, 1e-300, 1e300]);
    xs
}

/// Eight values, for quadruples.
fn eight() -> Vec<f64> {
    vec![0.0, -1.0, 1.0, 1e-6, 1000.0, f64::NAN, f64::INFINITY, 1e300]
}

/// A slice's tokens: `v:<n>` then its values.
fn sl(xs: &[f64]) -> Vec<Tok> {
    let mut t = vec![Tok::V(xs.len())];
    t.extend(xs.iter().map(|&x| f(x)));
    t
}

fn slice_row(m: &mut Module, xs: &[f64]) {
    let ins = sl(xs);
    let x = black_box(xs);
    m.row("circuit::resistance_series", Class::A, &ins, &res_f(resistance_series(x)));
    m.row("circuit::resistance_parallel", Class::A, &ins, &res_f(resistance_parallel(x)));
    m.row("circuit::capacitance_series", Class::A, &ins, &res_f(capacitance_series(x)));
    m.row("circuit::capacitance_parallel", Class::A, &ins, &res_f(capacitance_parallel(x)));
    m.row("circuit::inductance_series", Class::A, &ins, &res_f(inductance_series(x)));
    m.row("circuit::inductance_parallel", Class::A, &ins, &res_f(inductance_parallel(x)));
}

fn transient_row(m: &mut Module, a: f64, r: f64, c: f64, t: f64) {
    let ins = [f(a), f(r), f(c), f(t)];
    m.row("circuit::rc_charging_voltage", Class::C, &ins, &[f(rc_charging_voltage(black_box(a), r, c, t))]);
    m.row("circuit::rc_discharging_voltage", Class::C, &ins, &[f(rc_discharging_voltage(black_box(a), r, c, t))]);
    m.row("circuit::rl_current_rise", Class::C, &ins, &res_f(rl_current_rise(black_box(a), r, c, t)));
    m.row("circuit::rl_current_decay", Class::C, &ins, &res_f(rl_current_decay(black_box(a), r, c, t)));
}

pub fn emit() -> Module {
    let mut m = Module::new("circuit");
    let mut rng = Rng::new(0xc1_2c17);
    let ps = pairs(&mut rng);
    let sh = short();
    let e8 = eight();

    // Pairs: Ohm's law, power, τ = RC, the capacitor's charge and energy, the inductor's energy, τ = L/R and the
    // resonances.
    for &a in &ps {
        for &b in &ps {
            let two = [f(a), f(b)];
            m.row("circuit::voltage", Class::A, &two, &[f(voltage(black_box(a), b))]);
            m.row("circuit::current", Class::A, &two, &res_f(current(black_box(a), b)));
            m.row("circuit::resistance", Class::A, &two, &res_f(resistance(black_box(a), b)));
            m.row("circuit::power", Class::A, &two, &[f(power(black_box(a), b))]);
            m.row("circuit::capacitor_energy", Class::A, &two, &[f(capacitor_energy(black_box(a), b))]);
            m.row("circuit::capacitor_charge", Class::A, &two, &[f(capacitor_charge(black_box(a), b))]);
            m.row("circuit::inductor_energy", Class::A, &two, &[f(inductor_energy(black_box(a), b))]);
            m.row("circuit::rc_time_constant", Class::A, &two, &[f(rc_time_constant(black_box(a), b))]);
            m.row("circuit::rl_time_constant", Class::A, &two, &res_f(rl_time_constant(black_box(a), b)));
            m.row("circuit::resonant_frequency", Class::A, &two, &res_f(resonant_frequency(black_box(a), b)));
            m.row(
                "circuit::resonant_angular_frequency",
                Class::A,
                &two,
                &res_f(resonant_angular_frequency(black_box(a), b)),
            );
        }
    }

    // Triples: resistance from geometry, the parallel plate, Q and ζ.
    for &a in &sh {
        for &b in &sh {
            for &c in &sh {
                let three = [f(a), f(b), f(c)];
                m.row(
                    "circuit::resistance_from_geometry",
                    Class::A,
                    &three,
                    &res_f(resistance_from_geometry(black_box(a), b, c)),
                );
                m.row(
                    "circuit::parallel_plate_capacitance",
                    Class::A,
                    &three,
                    &res_f(parallel_plate_capacitance(black_box(a), b, c)),
                );
                m.row("circuit::quality_factor", Class::A, &three, &res_f(quality_factor(black_box(a), b, c)));
                m.row("circuit::damping_ratio", Class::A, &three, &res_f(damping_ratio(black_box(a), b, c)));
            }
        }
    }
    // Q's R < 1e-30 cut, ulp by ulp; seeded physical RLC circuits.
    for k in -2..=2 {
        let r = ulps(1e-30, k);
        m.row("circuit::quality_factor", Class::A, &[f(r), f(1e-3), f(1e-6)], &res_f(quality_factor(black_box(r), 1e-3, 1e-6)));
    }
    for _ in 0..200 {
        let (r, l, c) = (rng.log_uniform(-1.0, 4.0), rng.log_uniform(-9.0, -1.0), rng.log_uniform(-12.0, -3.0));
        let three = [f(r), f(l), f(c)];
        m.row("circuit::quality_factor", Class::A, &three, &res_f(quality_factor(black_box(r), l, c)));
        m.row("circuit::damping_ratio", Class::A, &three, &res_f(damping_ratio(black_box(r), l, c)));
        let w = rng.log_uniform(2.0, 9.0);
        let four = [f(r), f(l), f(c), f(w)];
        m.row("circuit::rlc_impedance", Class::A, &four, &res_f(rlc_impedance(black_box(r), l, c, w)));
        m.row("circuit::resonant_frequency", Class::A, &[f(l), f(c)], &res_f(resonant_frequency(black_box(l), c)));
    }

    // Quadruples: the solenoid and the RLC impedance over eight values, the impedance's two zero cuts.
    for &a in &e8 {
        for &b in &e8 {
            for &c in &e8 {
                for &d in &e8 {
                    let four = [f(a), f(b), f(c), f(d)];
                    m.row("circuit::solenoid_inductance", Class::A, &four, &[f(solenoid_inductance(black_box(a), b, c, d))]);
                    m.row("circuit::rlc_impedance", Class::A, &four, &res_f(rlc_impedance(black_box(a), b, c, d)));
                }
            }
        }
    }
    for k in -1..=1 {
        for x in [ulps(1e-30, k), -ulps(1e-30, k)] {
            let w = [f(10.0), f(1e-3), f(1e-6), f(x)];
            m.row("circuit::rlc_impedance", Class::A, &w, &res_f(rlc_impedance(black_box(10.0), 1e-3, 1e-6, x)));
            let c = [f(10.0), f(1e-3), f(x), f(1000.0)];
            m.row("circuit::rlc_impedance", Class::A, &c, &res_f(rlc_impedance(black_box(10.0), 1e-3, x, 1000.0)));
        }
    }

    // The transients over eight-value quadruples, the τ cuts ulp by ulp, and seeded physical circuits.
    for &a in &[1.0, -2.5, 0.0, f64::NAN, f64::INFINITY] {
        for &r in &e8 {
            for &c in &e8 {
                for &t in &e8 {
                    transient_row(&mut m, a, r, c, t);
                }
            }
        }
    }
    for k in -2..=2 {
        // τ = R·C about 1e-30 (R = 1e-15, C near 1e-15), and τ = L/R about 1e-30 (L near 1e-30, R = 1).
        transient_row(&mut m, 5.0, 1e-15, ulps(1e-15, k), 1e-31);
        transient_row(&mut m, 5.0, 1.0, ulps(1e-30, k), 1e-31);
    }
    // The knife edges of the exact τ tests (bijli ADR 0002 row 17's rule, M2 part 2c-1): R·C and L/R within an ulp of
    // 10^-30, where Rust's rounded τ and the exact one can fall on either side.
    for (r, c) in knife_edge2(Edge2::Prod, &[1e-15, 3e-16, 7e-11, 1.0, 1e10], 6, 24) {
        transient_row(&mut m, 5.0, r, c, 1e-31);
    }
    for (l, r) in knife_edge2(Edge2::Quot, &[1.0, 3.0, 7.0, 0.1, 1e10], 0, 48) {
        transient_row(&mut m, 5.0, r, l, 1e-31);
    }
    // And the band where they disagree, 0.024 ulp wide (grids.rs's knife_band2).
    for (r, c) in knife_band2(Edge2::Prod, &[1.7e-15, 2.9e-11], 1500) {
        transient_row(&mut m, 5.0, r, c, 1e-31);
    }
    for (l, r) in knife_band2(Edge2::Quot, &[3.3, 0.71], 1500) {
        transient_row(&mut m, 5.0, r, l, 1e-31);
    }
    for _ in 0..300 {
        let (r, c) = (rng.log_uniform(0.0, 6.0), rng.log_uniform(-12.0, -3.0));
        let t = rng.log_uniform(-12.0, 1.0);
        transient_row(&mut m, rng.uniform(-10.0, 10.0), r, c, t);
    }

    // Slices: empty, every pair value alone, short-list pairs, a few triples, and the 1/Σ(1/x) cut.
    slice_row(&mut m, &[]);
    for &x in &ps {
        slice_row(&mut m, &[x]);
    }
    for &a in &sh {
        for &b in &sh {
            slice_row(&mut m, &[a, b]);
        }
    }
    for xs in [[1.0, 2.0, 3.0], [100.0, 220.0, 470.0], [1.0, -1.0, 2.0], [1e-6, 2.2e-6, 4.7e-6], [-0.0, 0.0, -0.0]] {
        slice_row(&mut m, &xs);
    }
    slice_row(&mut m, &[-0.0]);
    slice_row(&mut m, &[1.0, -1.0]);
    for k in -3..=3 {
        slice_row(&mut m, &[ulps(1e30, k)]);
        slice_row(&mut m, &[ulps(2e30, k), ulps(2e30, -k)]);
    }
    for _ in 0..100 {
        let n = usize::try_from(rng.int(1, 6)).unwrap_or(1);
        let xs: Vec<f64> = (0..n).map(|_| rng.log_uniform(-12.0, 6.0)).collect();
        slice_row(&mut m, &xs);
    }
    emit_inverse_edges(&mut m);
    emit_fold_edges(&mut m);
    emit_transient_edges(&mut m);
    m
}

/// One list through the three inverse combinations (1/Σ 1/x).
fn inverse_rows(m: &mut Module, xs: &[f64]) {
    let ins = sl(xs);
    let x = black_box(xs);
    m.row("circuit::resistance_parallel", Class::A, &ins, &res_f(resistance_parallel(x)));
    m.row("circuit::capacitance_series", Class::A, &ins, &res_f(capacitance_series(x)));
    m.row("circuit::inductance_parallel", Class::A, &ins, &res_f(inductance_parallel(x)));
}

/// `x` moved `k` ulps away from zero (toward it for k < 0), either sign.
fn ulps_signed(x: f64, k: i64) -> f64 {
    if x < 0.0 { -ulps(-x, k) } else { ulps(x, k) }
}

/// M2 part 2c-5 (bijli ADR 0002 row 17): the inverse combinations at 1/Σ(1/x)'s own pole, Σ 1/x = 0, and at its cut.
/// |Σ 1/x| within ulps of 10^-30 (one value near ±10^30, k equal values near k·10^30, a pair a, −(a + j ulps) whose
/// difference lands on 10^-30·a²); sums that cancel (pairs ulps apart, triples closing on a third value); exact
/// midpoints of 1/S (M₁·2^s and (2^28 − M₁)·2^s for odd M₁ near 2^27); mixed signs; ±DBL_MAX, 1e-30 and ±∞ among
/// them; long lists. Its own seed, after every older row.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn emit_inverse_edges(m: &mut Module) {
    let mut rng = Rng::new(0x2C5_1A7E);
    for k in -8..=8 {
        inverse_rows(m, &[ulps(1e30, k)]);
        inverse_rows(m, &[-ulps(1e30, k)]);
    }
    for c in [2.0, 3.0, 7.0, 15.0] {
        for j in -6..=6 {
            let v = ulps(c * 1e30, j);
            inverse_rows(m, &vec![v; c as usize]);
        }
    }
    for a in [3e14, 7e14, 1.5e15, 4e15] {
        let u = ulps(a, 1) - a;
        let j0 = (1e-30 * a * a / u).round() as i64;
        for j in (j0 - 6).max(1)..=j0 + 6 {
            inverse_rows(m, &[a, -ulps(a, j)]);
        }
    }
    for i in 0..100 {
        let a = rng.uniform(1.0, 2.0) * pow2(rng.int(-60, 60) as i32);
        let k = [1, 2, 3, 4, -1, -2][i % 6];
        inverse_rows(m, &[a, -ulps(a, k)]);
    }
    for i in 0..100 {
        let a = rng.uniform(1.0, 2.0) * pow2(rng.int(-40, 40) as i32);
        let b = rng.uniform(1.0, 2.0) * pow2(rng.int(-40, 40) as i32);
        let t = -1.0 / (1.0 / a + 1.0 / b);
        inverse_rows(m, &[a, b, ulps_signed(t, [0, 1, -1, 2, -2][i % 5])]);
    }
    for i in 0..60 {
        let m1 = 2 * rng.int((1 << 25) + (1 << 24), (1 << 26) + (1 << 25) - 1) + 1;
        let sc = pow2(rng.int(-60, 60) as i32);
        let (a, b) = (m1 as f64 * sc, ((1_i64 << 28) - m1) as f64 * sc);
        match i % 3 {
            0 => inverse_rows(m, &[a, b]),
            1 => inverse_rows(m, &[b, a]),
            _ => inverse_rows(m, &[-a, -b]),
        }
    }
    for _ in 0..200 {
        let n = rng.int(2, 15) as usize;
        let xs: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect();
        inverse_rows(m, &xs);
    }
    let pool = [f64::MAX, -f64::MAX, 1e-30, -1e-30, 1e300, f64::INFINITY, f64::NEG_INFINITY, 1.0, 3.0, -3.0];
    for i in 0..40 {
        let n = rng.int(1, 6) as usize;
        let xs: Vec<f64> = if i % 10 == 0 {
            vec![f64::INFINITY; n]
        } else {
            (0..n).map(|_| pool[rng.int(0, 9) as usize]).collect()
        };
        inverse_rows(m, &xs);
    }
    for n in [64_usize, 200] {
        let xs: Vec<f64> = (0..n).map(|_| rng.log_uniform(0.0, 4.0)).collect();
        inverse_rows(m, &xs);
        let alt: Vec<f64> = (0..n).map(|j| if j % 2 == 0 { 3.0 } else { -3.0 }).collect();
        inverse_rows(m, &alt);
        let mut near = alt.clone();
        near[n - 1] = -ulps(3.0, 1);
        inverse_rows(m, &near);
    }
}

/// The three series folds on one slice.
fn series_rows(m: &mut Module, xs: &[f64]) {
    let ins = sl(xs);
    let x = black_box(xs);
    m.row("circuit::resistance_series", Class::A, &ins, &res_f(resistance_series(x)));
    m.row("circuit::capacitance_parallel", Class::A, &ins, &res_f(capacitance_parallel(x)));
    m.row("circuit::inductance_series", Class::A, &ins, &res_f(inductance_series(x)));
}

/// M2 part 2d-3a (bijli ADR 0002 row 18's rule): the series folds whose partial sums pass DBL_MAX (built, and drawn
/// from a pool of values near DBL_MAX of both signs), with absorption and non-finite controls; and the impedance's
/// chain at its range edges: a grid of R, L, C and ω from 0 and 2^-1074 to DBL_MAX, and exact resonances ωL = 1/(ωC)
/// at scales 2^±100 to 2^±580 under tiny and vast R. Its own seed, after every older row.
#[allow(clippy::cast_possible_truncation)]
fn emit_fold_edges(m: &mut Module) {
    let mut rng = Rng::new(0x0002_D3A5);
    let big = f64::MAX;
    let built = [
        vec![big, big, -big],
        vec![big, big, -big, -big],
        vec![-big, -big, big],
        vec![1e308, 1e308, -1e308, 5.0],
        vec![big, big / 2.0, -big],
        vec![big, big, -big, 1.0, -big / 2.0],
        vec![1.5e308, 1.5e308, -1e308, -1e308, -1e308],
        vec![big, 1.0, -big],
        vec![big, big, f64::INFINITY],
        vec![big, big, -big, f64::NAN],
    ];
    for xs in &built {
        series_rows(m, xs);
    }
    let pool = [big, -big, big / 2.0, -big / 2.0, 1e308, -1e308, 1.5e308, -1.5e308, 1.0, -0.0, pow2(-1074)];
    for _ in 0..60 {
        let n = rng.int(2, 7) as usize;
        let xs: Vec<f64> = (0..n).map(|_| pool[rng.int(0, 10) as usize]).collect();
        series_rows(m, &xs);
    }
    let rs = [0.0, pow2(-1074), 1e-200, 1e-160, 10.0, 1e160, 1e200, big];
    let ls = [0.0, 1e-300, 1e-3, 1e300, big];
    let cs = [2e-30, 1e-6, 1e300, big];
    let ws = [2e-30, 1e3, 1e160, 1e300, big, -1e3];
    for &r in &rs {
        for &l in &ls {
            for &c in &cs {
                for &w in &ws {
                    let four = [f(r), f(l), f(c), f(w)];
                    m.row("circuit::rlc_impedance", Class::A, &four, &res_f(rlc_impedance(black_box(r), l, c, w)));
                }
            }
        }
    }
    for k in [100, 300, 500, 580] {
        let (l, c, w) = (pow2(90 - 2 * k), pow2(-90), pow2(k));
        for r in [0.0, 1e-200, 1e-160, 10.0, 1e160, 1e200] {
            let four = [f(r), f(l), f(c), f(w)];
            m.row("circuit::rlc_impedance", Class::A, &four, &res_f(rlc_impedance(black_box(r), l, c, w)));
        }
        let ln = ulps(l, 1);
        let four = [f(0.0), f(ln), f(c), f(w)];
        m.row("circuit::rlc_impedance", Class::A, &four, &res_f(rlc_impedance(black_box(0.0), ln, c, w)));
    }
}

/// M2 part 2d-3b (bijli ADR 0002 row 18's rule): the four transients at their range edges. Each (a, R, C or L, t) feeds
/// all four, so the pairs set both τ's: τ = RC past DBL_MAX (10^200 twice, 2^1000 twice) with τ = L/R = 1; τ = L/R
/// past DBL_MAX with V/R past it too (R = 10^-29); τ = 10^10, where t = 10^-300 gives a subnormal t/τ; and 1. The
/// times put e^(−t/τ) subnormal (720), below it (800), past DBL_MAX (−750, −800), past e^1500 (−2000: with V/R =
/// 2^-2074 the product is 2^811), and t/τ past DBL_MAX (±10^300); the amplitudes span ±10^300, 10^-300, 2^-1074 and 0.
/// After every older row.
fn emit_transient_edges(m: &mut Module) {
    let amps = [1.0, 1e300, -1e300, 1e-300, 0.0, pow2(-1074)];
    let pairs = [(1.0, 1.0), (1e200, 1e200), (1e-29, 1e300), (1e10, 1e10), (pow2(1000), pow2(1000)), (1.0, 1e10)];
    let ts = [1.0, -1.0, 1e-300, 720.0, 800.0, -750.0, -800.0, -2000.0, 1e300, -1e300, 1e-30];
    for &a in &amps {
        for &(r, c) in &pairs {
            for &t in &ts {
                transient_row(m, a, r, c, t);
            }
        }
    }
}
