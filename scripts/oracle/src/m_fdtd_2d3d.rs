//! `fdtd_2d3d.tsv` — rust-old/src/fdtd.rs's 2D and 3D grids and their CPML (M6). `fdtd.tsv` holds M5's FDTD core and
//! is near the corpus's per-file cap, so M6's rows start this module.
//!
//! M6 part 1a: the CPML layer's grading and profiles.
//! - `CpmlParams::new` (class A): `i:<thickness>`; its five fields.
//! - `Cpml::with_params` (class C: powf and exp): the thickness, grading order, sigma_ratio, kappa_max, alpha_max, dx
//!   and dt; its thickness, `grading_order()`, `sigma_max()`, `kappa_max()`, `alpha_max()`, then the six profiles
//!   `b_node`, `c_node`, `km1_node`, `b_half`, `c_half`, `km1_half` (each `v:<len>` and its values, read from the
//!   layer's `Debug` text), or the error.
//! - `Cpml::new` (class C): the thickness, dx and dt; the same outputs (an invalid input gives the inert layer: σ_max 0
//!   and empty profiles).
//! - `Cpml::attach` (class A): the thickness, dx and dt of a `Cpml::new` layer, then `nx`, `ny`; the attached grid's
//!   nx and ny and the four Ψ fields' lengths (`psi_fx_y`, `psi_fz_y`, `psi_fy_x`, `psi_fz_x`), from `Debug`. Sizes
//!   are kept small, or past 2^64 (where Rust's `checked_mul` leaves the layer unattached): an allocation between
//!   depends on the machine.
//!
//! The grid: the validation of each parameter over the specials and its boundary; thicknesses 2^60 + 1 and usize::MAX,
//! whose profiles no allocation can hold (smaller huge ones depend on the machine); the profiles around the classic
//! grading, each parameter swept alone, at 2D and 3D Courant steps; the extremes of dx, dt, κ_max, α_max and m; and
//! seeded gradings.

use std::collections::BTreeSet;
use std::hint::black_box;

use crate::grids::{Rng, specials};
use crate::m_fdtd::{dbg_field, dbg_f64, dbg_vec, res_toks, seq, u};
use crate::tsv::{Class, Module, Tok, f};
use bijli::fdtd::{Cpml, CpmlParams};

/// The module's rows, each (function, inputs) once.
struct Rows {
    m: Module,
    keys: BTreeSet<String>,
}

impl Rows {
    /// A row, skipped when its (function, inputs) is already in.
    fn emit(&mut self, func: &str, cls: Class, ins: &[Tok], outs: &[Tok]) {
        if self.keys.insert(format!("{func} {ins:?}")) {
            self.m.row(func, cls, ins, outs);
        }
    }
}

/// A grading's tokens: thickness, m, sigma_ratio, kappa_max, alpha_max.
fn params_toks(p: &CpmlParams) -> Vec<Tok> {
    vec![u(p.thickness), f(p.grading_order), f(p.sigma_ratio), f(p.kappa_max), f(p.alpha_max)]
}

/// A grading of the five values.
fn params(t: usize, m: f64, sr: f64, km: f64, am: f64) -> CpmlParams {
    let mut p = CpmlParams::new(t);
    p.grading_order = m;
    p.sigma_ratio = sr;
    p.kappa_max = km;
    p.alpha_max = am;
    p
}

/// A layer's outputs (module comment).
fn layer_toks(c: &Cpml) -> Vec<Tok> {
    let d = format!("{c:?}");
    let mut out = vec![u(c.thickness), f(c.grading_order()), f(c.sigma_max()), f(c.kappa_max()), f(c.alpha_max())];
    assert_eq!(dbg_f64(&d, "sigma_max").to_bits(), c.sigma_max().to_bits(), "Debug's sigma_max is the accessor's");
    for name in ["b_node", "c_node", "km1_node", "b_half", "c_half", "km1_half"] {
        out.extend(seq(&dbg_vec(&d, name)));
    }
    out
}

impl Rows {
    fn with_params(&mut self, p: &CpmlParams, dx: f64, dt: f64) {
        let mut ins = params_toks(p);
        ins.extend([f(dx), f(dt)]);
        let outs = res_toks(Cpml::with_params(black_box(p), black_box(dx), black_box(dt)), |c| layer_toks(&c));
        self.emit("Cpml::with_params", Class::C, &ins, &outs);
    }

    fn new_layer(&mut self, t: usize, dx: f64, dt: f64) {
        let c = Cpml::new(black_box(t), black_box(dx), black_box(dt));
        self.emit("Cpml::new", Class::C, &[u(t), f(dx), f(dt)], &layer_toks(&c));
    }

    fn attach(&mut self, t: usize, dx: f64, dt: f64, nx: usize, ny: usize) {
        let mut c = Cpml::new(t, dx, dt);
        c.attach(black_box(nx), black_box(ny));
        let d = format!("{c:?}");
        let mut outs = vec![u(dbg_field(&d, "nx").parse().expect("nx")), u(dbg_field(&d, "ny").parse().expect("ny"))];
        for name in ["psi_fx_y", "psi_fz_y", "psi_fy_x", "psi_fz_x"] {
            outs.push(u(dbg_vec(&d, name).len()));
        }
        self.emit("Cpml::attach", Class::A, &[u(t), f(dx), f(dt), u(nx), u(ny)], &outs);
    }
}

/// The 2D grid's step for spacing dx, dx/(2√2·c), as `Fdtd2d::new` computes it.
fn dt_2d(dx: f64) -> f64 {
    dx / (2.0 * std::f64::consts::SQRT_2 * bijli::field::SPEED_OF_LIGHT)
}

/// The 3D grid's step for spacing dx, dx/(2√3·c), as `Fdtd3d::new` computes it.
fn dt_3d(dx: f64) -> f64 {
    dx / (2.0 * 3.0_f64.sqrt() * bijli::field::SPEED_OF_LIGHT)
}

fn emit_1a(r: &mut Rows) {
    for t in [0, 1, 8, 12, usize::MAX] {
        let p = CpmlParams::new(black_box(t));
        r.emit("CpmlParams::new", Class::A, &[u(t)], &params_toks(&p));
    }
    let (dx, dt) = (1e-4, dt_2d(1e-4));
    // The validation: each parameter over the specials and either side of its bound, then dx and dt; Rust's order
    // (the parameters, then dx, then dt, then the profiles' allocation, then their finiteness).
    let ulp_below = |x: f64| f64::from_bits(x.to_bits() - 1);
    for v in specials().into_iter().chain([ulp_below(1.0), 0.5, 2.0, 1e300]) {
        r.with_params(&params(8, v, 1.0, 1.0, 0.0), dx, dt);
        r.with_params(&params(8, 3.0, v, 1.0, 0.0), dx, dt);
        r.with_params(&params(8, 3.0, 1.0, v, 0.0), dx, dt);
        r.with_params(&params(8, 3.0, 1.0, 1.0, v), dx, dt);
        r.with_params(&CpmlParams::new(8), v, dt);
        r.with_params(&CpmlParams::new(8), dx, v);
        r.with_params(&params(0, v, v, v, v), v, v);
    }
    for t in [0, 1, (1 << 60) + 1, usize::MAX] {
        r.with_params(&CpmlParams::new(t), dx, dt);
        r.with_params(&CpmlParams::new(t), f64::NAN, dt);
    }
    r.with_params(&params(3, 3.0, f64::MAX, 1.0, 0.0), dx, dt);
    r.with_params(&params(3, 3.0, 1.0, 1.0, 0.0), f64::from_bits(1), dt);
    // The profiles: the classic grading at every thickness, then each parameter swept alone, at the 2D and 3D steps.
    for dxs in [1e-4, 1e-3, 1e-2] {
        for dts in [dt_2d(dxs), dt_3d(dxs)] {
            for t in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 32] {
                r.with_params(&CpmlParams::new(t), dxs, dts);
            }
            for m in [1.0, 2.0, 2.5, 3.5, 4.0, 5.0] {
                r.with_params(&params(10, m, 1.0, 1.0, 0.0), dxs, dts);
            }
            for sr in [0.5, 0.7, 1.2, 1.5, 2.0] {
                r.with_params(&params(10, 3.0, sr, 1.0, 0.0), dxs, dts);
            }
            for km in [1.5, 5.0, 7.0, 11.0, 15.0] {
                r.with_params(&params(10, 3.0, 1.0, km, 0.0), dxs, dts);
            }
            for am in [1e-6, 1e-3, 0.01, 0.05, 1.0] {
                r.with_params(&params(10, 3.0, 1.0, 1.0, am), dxs, dts);
                r.with_params(&params(10, 4.0, 1.0, 5.0, am), dxs, dts);
            }
        }
    }
    r.with_params(&params(6, 3.0, 1.0, 7.0, 0.0), 1e-4, 1e-13);
    // The extremes: dx, dt, κ_max, α_max and m far from their use (Rust's values; the range edges are a later part).
    for (dxe, dte) in [(1e-300, 1e-310), (1e300, 1e290), (1e-12, 1e-21), (1e10, 1e-13), (1e-4, 1e-300), (1e-4, 1e300),
                       (1e-4, f64::MAX), (f64::MAX, f64::MAX), (f64::MIN_POSITIVE, f64::MIN_POSITIVE)] {
        r.with_params(&CpmlParams::new(10), dxe, dte);
        r.with_params(&params(10, 4.0, 1.3, 9.0, 0.02), dxe, dte);
    }
    for km in [1e10, 1e300, f64::MAX] {
        r.with_params(&params(5, 3.0, 1.0, km, 0.0), dx, dt);
        r.with_params(&params(5, 3.0, 1.0, km, 0.5), dx, dt);
    }
    for am in [1e10, 1e300, f64::MAX] {
        r.with_params(&params(5, 3.0, 1.0, 1.0, am), dx, dt);
        r.with_params(&params(5, 3.0, 1.0, 7.0, am), dx, dt);
    }
    for m in [10.0, 100.0, 1e3, 1e10, 1e300] {
        r.with_params(&params(7, m, 1.0, 1.0, 0.0), dx, dt);
        r.with_params(&params(7, m, 1.0, 3.0, 0.1), dx, dt);
    }
    for sr in [1e-300, 1e300] {
        r.with_params(&params(7, 3.0, sr, 1.0, 0.0), dx, dt);
    }
    let mut rng = Rng::new(0x6c97_1a00);
    for _ in 0..120 {
        let t = usize::try_from(rng.int(1, 24)).unwrap_or(1);
        let m = 1.0 + 4.0 * rng.unit();
        let sr = rng.log_uniform(-1.0, 0.5);
        let km = if rng.unit() < 0.5 { 1.0 } else { 1.0 + 14.0 * rng.unit() };
        let am = if rng.unit() < 0.4 { 0.0 } else { rng.log_uniform(-6.0, 0.0) };
        let dxs = rng.log_uniform(-6.0, -1.0);
        let dts = if rng.unit() < 0.5 { dt_2d(dxs) } else { dt_3d(dxs) * (0.5 + rng.unit()) };
        r.with_params(&params(t, m, sr, km, am), dxs, dts);
    }
    // The infallible constructor: valid input, and each invalid input's inert layer.
    for (t, dxn, dtn) in [(10, 1e-4, dt_2d(1e-4)), (1, 1e-3, dt_3d(1e-3)), (32, 1e-2, 1e-13), (0, 1e-4, 1e-13),
                          (8, -1.0, 1e-13), (8, f64::NAN, 1e-13), (8, 0.0, 1e-13), (8, 1e-4, 0.0),
                          (8, 1e-4, f64::INFINITY), (usize::MAX, 1e-4, 1e-13), (3, f64::from_bits(1), 1e-13)] {
        r.new_layer(t, dxn, dtn);
    }
    // Attaching: small grids, a grid too small to use, an empty one, and sizes past 2^64.
    for (t, nx, ny) in [(1, 5, 5), (4, 11, 20), (10, 60, 60), (10, 81, 81), (4, 10, 20), (2, 3, 3), (1, 1, 1000),
                        (3, 0, 0), (3, 0, 7), (13, 30, 40), (2, usize::MAX, 5), (2, 5, usize::MAX),
                        (4, 1 << 62, 1), (1, 1 << 63, 1), (usize::MAX, 4, 4)] {
        r.attach(t, 1e-4, dt_2d(1e-4), nx, ny);
    }
    r.attach(8, -1.0, 1e-13, 60, 60);
    r.attach(0, 1e-4, 1e-13, 9, 9);
}

pub fn emit() -> Module {
    let mut r = Rows { m: Module::new("fdtd_2d3d"), keys: BTreeSet::new() };
    emit_1a(&mut r);
    r.m
}
