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
//! M6 part 1b-1: `Fdtd2d` without its CPML, as scripts (class T, port plan §6 step 5's trajectory format): the
//! grid's `i:<nx> i:<ny> f:<dx> i:<mode>` (0 TM, 1 TE), then its operations in order, each an opcode and its arguments:
//! - `i:1|2|3 i:<x0> i:<y0> i:<x1> i:<y1> f:<v>`: `set_permittivity`, `set_permeability`, `set_conductivity`;
//! - `i:4 i:<x0> i:<y0> i:<x1> i:<y1> f:<eps_r> f:<mu_r> f:<sigma>`: `set_material`;
//! - `i:5 i:<x> i:<y> f:<v>`: `add_source`;
//! - `i:6 i:<k>`: `step_once`, k times;
//! - `i:7 i:<steps> i:<x> i:<y>|n: f:<frequency> f:<amplitude>`: `run`, `n:` for `None`;
//! - `i:8 i:<field> f:<v>`: every value of `fz` (0), `fx` (1) or `fy` (2) set to v (Rust's public vectors);
//! - `i:9 i:<thickness>`: `try_enable_cpml`; `i:10 i:<thickness> f:<m> f:<sigma_ratio> f:<kappa_max> f:<alpha_max>`:
//!   `try_enable_cpml_with`; `i:11 i:<thickness>`: `enable_cpml`; `i:12`: `disable_cpml` (part 1b-2).
//!
//! The outputs: `new`'s error, or each failing operation's error in order, then the state: `i:<step>`, `f:<dt>`,
//! `f:<time()>`, `f:<total_energy()>`, the digests (`i:`, as `fdtd.tsv`'s) of `fz`, `fx`, `fy`, `permittivity`,
//! `permeability`, `ce` and `ch`, then `ca`'s and `conductivity`'s lengths and digests (empty while lossless), then
//! `f:` `fz`, `fx` and `fy` at the probes (nx/2, ny/2), (1, 1), (nx − 2, ny − 2), (0, ny/2) and (nx/2, 1), then,
//! while a CPML is attached, its thickness and the digests of its four Ψ fields (`psi_fx_y`, `psi_fz_y`, `psi_fy_x`,
//! `psi_fz_x`). The private vectors are read from the grid's and the layer's `Debug` text.
//!
//! M6 part 1c: `nf2ff_2d` (class C: sin_cos and hypot): `i:<box_margin> i:<n_angles> f:<frequency>`, then a 2D script
//! as above building the snapshot (its operations' errors aren't kept); its outputs the pattern and the angles (each
//! `v:<len>` and its values), or the error.
//!
//! The grid: the validation of each parameter over the specials and its boundary; thicknesses 2^60 + 1 and usize::MAX,
//! whose profiles no allocation can hold (smaller huge ones depend on the machine); the profiles around the classic
//! grading, each parameter swept alone, at 2D and 3D Courant steps; the extremes of dx, dt, κ_max, α_max and m; and
//! seeded gradings.

use std::collections::BTreeSet;
use std::hint::black_box;

use crate::grids::{Rng, specials};
use crate::m_fdtd::{dbg_field, dbg_f64, dbg_vec, digest, res_toks, seq, u};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::fdtd::{Cpml, CpmlParams, Fdtd2d, Mode2d, nf2ff_2d};
use bijli::material::Material;

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

// ── part 1b-1: Fdtd2d without its CPML ─────────────────────────────────────────────────────────────────────────

/// One operation of a 2D script.
#[derive(Clone, Copy)]
enum Op2 {
    Eps(usize, usize, usize, usize, f64),
    Mu(usize, usize, usize, usize, f64),
    Sigma(usize, usize, usize, usize, f64),
    Mat(usize, usize, usize, usize, f64, f64, f64),
    Source(usize, usize, f64),
    Step(usize),
    Run(usize, Option<(usize, usize)>, f64, f64),
    Fill(usize, f64),
    Cpml(usize),
    CpmlWith(usize, f64, f64, f64, f64),
    EnableCpml(usize),
    DisableCpml,
}

impl Op2 {
    fn toks(self) -> Vec<Tok> {
        match self {
            Self::Eps(a, b, c, d, x) => vec![i(1), u(a), u(b), u(c), u(d), f(x)],
            Self::Mu(a, b, c, d, x) => vec![i(2), u(a), u(b), u(c), u(d), f(x)],
            Self::Sigma(a, b, c, d, x) => vec![i(3), u(a), u(b), u(c), u(d), f(x)],
            Self::Mat(a, b, c, d, e, m, sg) => vec![i(4), u(a), u(b), u(c), u(d), f(e), f(m), f(sg)],
            Self::Source(x, y, v) => vec![i(5), u(x), u(y), f(v)],
            Self::Step(k) => vec![i(6), u(k)],
            Self::Run(k, src, fr, am) => {
                let mut t = vec![i(7), u(k)];
                match src {
                    Some((x, y)) => t.extend([u(x), u(y)]),
                    None => t.push(Tok::None),
                }
                t.extend([f(fr), f(am)]);
                t
            }
            Self::Fill(w, v) => vec![i(8), u(w), f(v)],
            Self::Cpml(t) => vec![i(9), u(t)],
            Self::CpmlWith(t, m, sr, km, am) => vec![i(10), u(t), f(m), f(sr), f(km), f(am)],
            Self::EnableCpml(t) => vec![i(11), u(t)],
            Self::DisableCpml => vec![i(12)],
        }
    }
}

/// A 2D grid's state tokens (module comment).
fn state2(sim: &Fdtd2d) -> Vec<Tok> {
    let d = format!("{sim:?}");
    let (nx, ny) = (sim.nx, sim.ny);
    let mut out = vec![u(sim.step), f(sim.dt), f(sim.time()), f(sim.total_energy()), i(digest(&sim.fz)),
                       i(digest(&sim.fx)), i(digest(&sim.fy))];
    for name in ["permittivity", "permeability", "ce", "ch"] {
        let v = dbg_vec(&d, name);
        assert_eq!(v.len(), nx * ny, "{name} has nx·ny values");
        out.push(i(digest(&v)));
    }
    for name in ["ca", "conductivity"] {
        let v = dbg_vec(&d, name);
        out.extend([u(v.len()), i(digest(&v))]);
    }
    for (x, y) in [(nx / 2, ny / 2), (1, 1), (nx - 2, ny - 2), (0, ny / 2), (nx / 2, 1)] {
        let k = y * nx + x;
        out.extend([f(sim.fz[k]), f(sim.fx[k]), f(sim.fy[k])]);
    }
    if let Some(c) = sim.cpml() {
        let dc = format!("{c:?}");
        out.push(u(c.thickness));
        for name in ["psi_fx_y", "psi_fz_y", "psi_fy_x", "psi_fz_x"] {
            out.push(i(digest(&dbg_vec(&dc, name))));
        }
    }
    out
}

/// A 2D script being built: its grid (the error when `new` failed), inputs and the operations' errors so far.
struct Script2 {
    sim: Result<Fdtd2d, Vec<Tok>>,
    ins: Vec<Tok>,
    errs: Vec<Tok>,
}

impl Script2 {
    fn new(nx: usize, ny: usize, dx: f64, mode: Mode2d) -> Self {
        let sim = Fdtd2d::new(black_box(nx), black_box(ny), black_box(dx), mode)
            .map_err(|e| vec![err(&e), s(e.to_string())]);
        let m = if mode == Mode2d::Te { 1 } else { 0 };
        Self { sim, ins: vec![u(nx), u(ny), f(dx), i(m)], errs: Vec::new() }
    }

    fn apply(&mut self, op: Op2) {
        self.ins.extend(op.toks());
        let sim = self.sim.as_mut().expect("operations only on a grid that was made");
        let r = match op {
            Op2::Eps(a, b, c, d, x) => sim.set_permittivity(a, b, c, d, x),
            Op2::Mu(a, b, c, d, x) => sim.set_permeability(a, b, c, d, x),
            Op2::Sigma(a, b, c, d, x) => sim.set_conductivity(a, b, c, d, x),
            Op2::Mat(a, b, c, d, e, m, sg) => sim.set_material(a, b, c, d, &Material::new(e, m, sg)),
            Op2::Source(x, y, v) => {
                sim.add_source(x, y, v);
                Ok(())
            }
            Op2::Step(k) => {
                for _ in 0..k {
                    sim.step_once();
                }
                Ok(())
            }
            Op2::Run(k, src, fr, am) => {
                sim.run(k, src, fr, am);
                Ok(())
            }
            Op2::Fill(w, v) => {
                match w {
                    0 => sim.fz.fill(v),
                    1 => sim.fx.fill(v),
                    _ => sim.fy.fill(v),
                }
                Ok(())
            }
            Op2::Cpml(t) => sim.try_enable_cpml(t),
            Op2::CpmlWith(t, m, sr, km, am) => sim.try_enable_cpml_with(&params(t, m, sr, km, am)),
            Op2::EnableCpml(t) => {
                sim.enable_cpml(t);
                Ok(())
            }
            Op2::DisableCpml => {
                sim.disable_cpml();
                Ok(())
            }
        };
        if let Err(e) = r {
            self.errs.extend([err(&e), s(e.to_string())]);
        }
    }

    fn outs(&self) -> Vec<Tok> {
        match &self.sim {
            Err(e) => e.clone(),
            Ok(sim) => {
                let mut out = self.errs.clone();
                out.extend(state2(sim));
                out
            }
        }
    }
}

impl Rows {
    /// A 2D scenario: `setup`, then a row at each frame (cumulative steps), `advance` bringing the grid there.
    #[allow(clippy::too_many_arguments)]
    fn frames2(&mut self, nx: usize, ny: usize, dx: f64, mode: Mode2d, setup: &[Op2], frames: &[usize],
               advance: impl Fn(usize) -> Op2) {
        let mut sc = Script2::new(nx, ny, dx, mode);
        if sc.sim.is_err() {
            self.emit("Fdtd2d::script", Class::T, &sc.ins, &sc.outs());
            return;
        }
        for &op in setup {
            sc.apply(op);
        }
        if frames.is_empty() {
            self.emit("Fdtd2d::script", Class::T, &sc.ins, &sc.outs());
        }
        let mut at = 0;
        for &fr in frames {
            sc.apply(advance(fr - at));
            at = fr;
            self.emit("Fdtd2d::script", Class::T, &sc.ins, &sc.outs());
        }
    }
}

fn emit_1b1(r: &mut Rows) {
    let modes = [Mode2d::Tm, Mode2d::Te];
    let step = Op2::Step;
    // Construction's errors, in Rust's order (the size, dx, the cell count, the allocation).
    for (nx, ny, dx) in [(2, 50, 1e-3), (50, 2, 1e-3), (0, 0, 1e-3), (2, 2, f64::NAN), (1, usize::MAX, 1e-3),
                         (usize::MAX, 3, 1e-3), (3, usize::MAX, 1e-3), (1 << 32, 1 << 30, 1e-3),
                         (1 << 33, 1 << 31, 1e-3)] {
        r.frames2(nx, ny, dx, Mode2d::Tm, &[], &[], step);
    }
    for &v in &specials() {
        r.frames2(10, 10, v, Mode2d::Tm, &[], &[], step);
    }
    for &mode in &modes {
        // A point source in vacuum on the bench's grid, and the smallest grids.
        r.frames2(100, 100, 1e-4, mode, &[Op2::Source(50, 50, 1.0)], &[1, 2, 10, 40, 100], step);
        for (nx, ny) in [(3, 3), (3, 5), (5, 3), (4, 4), (7, 5), (16, 9), (9, 16)] {
            let src = [Op2::Source(nx / 2, ny / 2, 1.0), Op2::Source(1, 1, -0.5)];
            r.frames2(nx, ny, 1e-3, mode, &src, &[1, 3, 10, 30], step);
        }
        // A source on the wall is cleared by the next step; one outside is ignored.
        let walls = [Op2::Source(0, 5, 1.0), Op2::Source(11, 9, 2.0), Op2::Source(12, 5, 1.0), Op2::Source(5, 10, 1.0),
                     Op2::Source(6, 4, 0.25)];
        r.frames2(12, 10, 1e-3, mode, &walls, &[1, 5], step);
        // Dielectric and magnetic regions.
        r.frames2(60, 40, 1e-4, mode, &[Op2::Eps(10, 5, 30, 35, 4.0), Op2::Mu(35, 0, 50, 40, 2.0),
                                          Op2::Source(20, 20, 1.0)], &[1, 20, 80], step);
        r.frames2(30, 30, 1e-3, mode, &[Op2::Mat(5, 5, 25, 25, 2.25, 1.5, 0.0), Op2::Source(15, 15, 1.0)], &[1, 30],
                  step);
        // Lossy regions: the bench's slab, a material with σ, then ε over it (σ kept), then σ = 0 (lossless again).
        r.frames2(100, 100, 1e-4, mode, &[Op2::Sigma(60, 0, 100, 100, 1.0), Op2::Source(30, 50, 1.0)], &[1, 20, 60],
                  step);
        r.frames2(40, 30, 1e-3, mode, &[Op2::Mat(10, 10, 20, 20, 2.0, 1.0, 0.3), Op2::Eps(10, 10, 20, 20, 4.0),
                                          Op2::Sigma(15, 10, 20, 20, 0.0), Op2::Source(12, 15, 1.0)], &[1, 25], step);
        r.frames2(40, 30, 1e-3, mode, &[Op2::Sigma(0, 0, 40, 30, 50.0), Op2::Source(20, 15, 1.0)], &[1, 10, 40], step);
        // A conductor's uniform decay (fdtd.rs's test): σΔt/ε = 0.05 on the whole grid, E set to 1.
        let sigma = 0.05 * 2.0 * bijli::field::EPSILON_0 / dt_2d(1e-3);
        let fill = if mode == Mode2d::Tm {
            vec![Op2::Fill(0, 1.0)]
        } else {
            vec![Op2::Fill(1, 1.0), Op2::Fill(2, 1.0)]
        };
        let mut setup = vec![Op2::Mat(0, 0, 41, 41, 2.0, 1.0, sigma)];
        setup.extend(fill);
        r.frames2(41, 41, 1e-3, mode, &setup, &[1, 19], step);
        // run with a sinusoidal source (the bench's), none, and one outside the grid.
        r.frames2(80, 80, 1e-4, mode, &[], &[10, 50], |k| Op2::Run(k, Some((40, 40)), 1e9, 1.0));
        r.frames2(40, 40, 1e-3, mode, &[Op2::Source(20, 20, 1.0)], &[5, 20], |k| Op2::Run(k, None, 1e9, 1.0));
        r.frames2(40, 40, 1e-3, mode, &[], &[5], |k| Op2::Run(k, Some((999, 0)), 1e9, 1.0));
        r.frames2(100, 100, 1e-4, mode, &[Op2::Eps(0, 0, 100, 100, 4.0)], &[50],
                  |k| Op2::Run(k, Some((50, 50)), 3e9, 1.0));
        r.frames2(30, 30, 1e-3, mode, &[], &[7, 30], |k| Op2::Run(k, Some((15, 15)), 3e10, 0.5));
    }
    // The regions' errors, each followed by the unchanged state.
    let bad_regions = [(10, 0, 12, 5), (0, 10, 5, 12), (5, 0, 4, 5), (0, 5, 5, 4), (usize::MAX, 0, 1, 1)];
    for (a, b, c, d) in bad_regions {
        r.frames2(10, 10, 1e-3, Mode2d::Tm, &[Op2::Eps(a, b, c, d, 2.0), Op2::Mu(a, b, c, d, 2.0),
                                                Op2::Sigma(a, b, c, d, 2.0), Op2::Mat(a, b, c, d, 2.0, 2.0, 2.0)],
                  &[], step);
    }
    r.frames2(10, 10, 1e-3, Mode2d::Tm, &[Op2::Eps(8, 8, 100, 100, 2.0), Op2::Eps(3, 3, 3, 3, 5.0),
                                            Op2::Sigma(0, 0, usize::MAX, usize::MAX, 0.0), Op2::Mu(9, 9, 10, 10, 3.0)],
              &[1], step);
    for &v in specials().iter().chain(&[-1.0, 1e-300, 1e300]) {
        let ops = [Op2::Eps(0, 0, 5, 5, v), Op2::Mu(0, 0, 5, 5, v), Op2::Sigma(0, 0, 5, 5, v), Op2::Source(4, 4, 1.0)];
        r.frames2(10, 10, 1e-3, Mode2d::Te, &ops, &[3], step);
    }
    for (e, m, sg) in [(4.0, 0.0, 0.0), (4.0, f64::NAN, 0.0), (f64::NAN, 1.0, 0.0), (4.0, 1.0, -1.0),
                       (4.0, 1.0, f64::INFINITY), (0.0, 0.0, -1.0), (2.0, 3.0, 0.0), (2.0, 3.0, 1e-3)] {
        let ops = [Op2::Mat(0, 0, 20, 20, e, m, sg), Op2::Mat(25, 0, 30, 20, 2.0, 1.0, 0.0), Op2::Source(10, 10, 1.0)];
        r.frames2(20, 20, 1e-3, Mode2d::Tm, &ops, &[2], step);
    }
    // Seeded scripts: grids from 3 to 30 a side, either mode, random regions, sources, fills and steps.
    let mut rng = Rng::new(0x1b1_2d00);
    for _ in 0..40 {
        let nx = usize::try_from(rng.int(3, 30)).unwrap_or(3);
        let ny = usize::try_from(rng.int(3, 30)).unwrap_or(3);
        let mode = if rng.unit() < 0.5 { Mode2d::Tm } else { Mode2d::Te };
        let dx = rng.log_uniform(-5.0, -1.0);
        let mut ops = Vec::new();
        let pick = |n: usize, rng: &mut Rng| usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(1))).unwrap_or(0);
        for _ in 0..usize::try_from(rng.int(1, 6)).unwrap_or(1) {
            let (x0, y0) = (pick(nx - 1, &mut rng), pick(ny - 1, &mut rng));
            let (x1, y1) = (x0 + pick(nx, &mut rng), y0 + pick(ny, &mut rng));
            let v = 1.0 + 9.0 * rng.unit();
            ops.push(match rng.int(0, 4) {
                0 => Op2::Eps(x0, y0, x1, y1, v),
                1 => Op2::Mu(x0, y0, x1, y1, v),
                2 => Op2::Sigma(x0, y0, x1, y1, rng.log_uniform(-3.0, 2.0)),
                3 => Op2::Mat(x0, y0, x1, y1, v, 1.0 + rng.unit(), if rng.unit() < 0.5 { 0.0 } else { rng.unit() }),
                _ => Op2::Fill(pick(2, &mut rng), rng.signed_log_uniform(-3.0, 1.0)),
            });
        }
        for _ in 0..usize::try_from(rng.int(1, 4)).unwrap_or(1) {
            ops.push(Op2::Source(pick(nx, &mut rng), pick(ny, &mut rng), rng.signed_log_uniform(-2.0, 1.0)));
        }
        let k = usize::try_from(rng.int(1, 40)).unwrap_or(1);
        r.frames2(nx, ny, dx, mode, &ops, &[k], step);
    }
}

// ── part 1b-2: the CPML attached to Fdtd2d ──────────────────────────────────────────────────────────────────────

fn emit_1b2(r: &mut Rows) {
    let step = Op2::Step;
    let cfs = |t| Op2::CpmlWith(t, 3.0, 1.0, 5.0, 2.0 * std::f64::consts::PI * bijli::field::EPSILON_0 * 1e9);
    for &mode in &[Mode2d::Tm, Mode2d::Te] {
        // The classic layer around a point source, the CFS layer around a driven one, Rust's bench grid.
        r.frames2(60, 60, 1e-4, mode, &[Op2::Cpml(10), Op2::Source(30, 30, 1.0)], &[1, 10, 50, 150], step);
        r.frames2(61, 61, 1e-4, mode, &[cfs(10)], &[10, 80], |k| Op2::Run(k, Some((30, 30)), 1e10, 1.0));
        r.frames2(100, 100, 1e-4, mode, &[Op2::Cpml(10), Op2::Source(50, 50, 1.0)], &[1, 20], step);
        // Graded otherwise: m, σ ratio, κ_max and α_max, on thin layers.
        for op in [Op2::CpmlWith(4, 2.0, 0.7, 1.0, 0.0), Op2::CpmlWith(6, 4.0, 1.5, 11.0, 0.05),
                   Op2::CpmlWith(3, 1.0, 1.0, 1.0, 1.0), Op2::CpmlWith(5, 3.5, 1.2, 3.0, 1e-3)] {
            r.frames2(24, 20, 1e-3, mode, &[op, Op2::Source(12, 10, 1.0), Op2::Source(3, 3, -0.5)], &[1, 7, 40], step);
        }
        // A dielectric and lossy block overlapping the strips (fdtd.rs's long-run test, shortened).
        r.frames2(60, 60, 1e-4, mode, &[Op2::Cpml(8), Op2::Mat(0, 20, 30, 40, 4.0, 1.0, 0.5)], &[50, 300],
                  |k| Op2::Run(k, Some((40, 30)), 3e10, 1.0));
        // The thinnest grids a layer fits (2t + 3 = min(nx, ny)) and a source inside a strip.
        r.frames2(7, 7, 1e-3, mode, &[Op2::Cpml(2), Op2::Source(3, 3, 1.0)], &[1, 4, 20], step);
        r.frames2(9, 13, 1e-3, mode, &[Op2::Cpml(3), Op2::Source(1, 6, 1.0), Op2::Source(4, 11, 0.5)], &[1, 12], step);
        r.frames2(8, 30, 1e-3, mode, &[Op2::Cpml(2), Op2::Source(4, 15, 1.0)], &[3, 25], step);
        // Enable, disable, replace: Ψ restarts at zero, a refused request keeps the layer.
        let ops = [Op2::Source(15, 15, 1.0), Op2::Step(5), Op2::Cpml(4), Op2::Step(5), Op2::Cpml(100),
                   Op2::EnableCpml(100), Op2::Step(3), Op2::Cpml(6), Op2::Step(4), Op2::DisableCpml, Op2::Step(2),
                   Op2::EnableCpml(3), Op2::EnableCpml(0)];
        r.frames2(30, 30, 1e-3, mode, &ops, &[1, 6], step);
        r.frames2(30, 30, 1e-3, mode, &[Op2::Cpml(5), Op2::Step(4), Op2::DisableCpml, Op2::Source(15, 15, 1.0)],
                  &[3], step);
        // Lossy cells inside and outside the strips.
        r.frames2(40, 40, 1e-3, mode, &[Op2::Cpml(6), Op2::Sigma(0, 0, 40, 10, 2.0), Op2::Sigma(15, 15, 25, 25, 50.0),
                                          Op2::Source(20, 30, 1.0)], &[1, 30], step);
    }
    // The thickness's and the grading's refusals, in Rust's order (the fit first, then the parameters).
    for t in [0, 1, 13, 14, 15, 100, usize::MAX] {
        r.frames2(30, 40, 1e-4, Mode2d::Te, &[Op2::Cpml(t), Op2::EnableCpml(t)], &[1], step);
    }
    for op in [Op2::CpmlWith(4, 0.5, 1.0, 1.0, 0.0), Op2::CpmlWith(4, 3.0, 0.0, 1.0, 0.0),
               Op2::CpmlWith(4, 3.0, 1.0, 0.5, 0.0), Op2::CpmlWith(4, 3.0, 1.0, 1.0, -1.0),
               Op2::CpmlWith(4, f64::NAN, f64::NAN, 1.0, 0.0), Op2::CpmlWith(0, 0.5, 0.0, 0.5, -1.0),
               Op2::CpmlWith(20, 0.5, 1.0, 1.0, 0.0), Op2::CpmlWith(4, 3.0, f64::MAX, 1.0, 0.0),
               Op2::CpmlWith(4, 3.0, 1.0, f64::INFINITY, 0.0)] {
        r.frames2(30, 30, 1e-3, Mode2d::Tm, &[op, Op2::Source(15, 15, 1.0)], &[2], step);
    }
    // Seeded scripts: a fitting layer of random grading on grids from 7 to 30 a side, with random materials and
    // sources.
    let mut rng = Rng::new(0x1b2_c9b1);
    for _ in 0..30 {
        let nx = usize::try_from(rng.int(7, 30)).unwrap_or(7);
        let ny = usize::try_from(rng.int(7, 30)).unwrap_or(7);
        let mode = if rng.unit() < 0.5 { Mode2d::Tm } else { Mode2d::Te };
        let tmax = i64::try_from((nx.min(ny) - 3) / 2).unwrap_or(2);
        let t = usize::try_from(rng.int(2, tmax.max(2))).unwrap_or(2);
        let m = 1.0 + 3.0 * rng.unit();
        let km = if rng.unit() < 0.5 { 1.0 } else { 1.0 + 10.0 * rng.unit() };
        let am = if rng.unit() < 0.5 { 0.0 } else { rng.log_uniform(-4.0, 0.0) };
        let mut ops = vec![Op2::CpmlWith(t, m, 0.5 + rng.unit(), km, am)];
        if rng.unit() < 0.5 {
            ops.push(Op2::Mat(0, 0, nx / 2, ny, 1.0 + 3.0 * rng.unit(), 1.0, if rng.unit() < 0.5 { 0.0 } else { 1.0 }));
        }
        let pick =
            |n: usize, rng: &mut Rng| usize::try_from(rng.int(0, i64::try_from(n).unwrap_or(1) - 1)).unwrap_or(0);
        for _ in 0..usize::try_from(rng.int(1, 3)).unwrap_or(1) {
            ops.push(Op2::Source(pick(nx, &mut rng), pick(ny, &mut rng), rng.signed_log_uniform(-2.0, 1.0)));
        }
        let k = usize::try_from(rng.int(1, 60)).unwrap_or(1);
        r.frames2(nx, ny, rng.log_uniform(-5.0, -2.0), mode, &ops, &[k], step);
    }
}

// ── part 1c: nf2ff_2d, the snapshot NF2FF ──────────────────────────────────────────────────────────────────────

impl Rows {
    /// An `nf2ff_2d` row: the grid of the script (nx, ny, dx, mode, ops), then the transform.
    #[allow(clippy::too_many_arguments)]
    fn nf2ff2(&mut self, margin: usize, na: usize, freq: f64, nx: usize, ny: usize, dx: f64, mode: Mode2d,
              ops: &[Op2]) {
        let mut sc = Script2::new(nx, ny, dx, mode);
        for &op in ops {
            sc.apply(op);
        }
        let mut ins = vec![u(margin), u(na), f(freq)];
        ins.extend(sc.ins.clone());
        let sim = sc.sim.as_ref().expect("a grid that was made");
        let outs = res_toks(nf2ff_2d(black_box(sim), margin, na, black_box(freq)), |res| {
            let mut o = seq(&res.pattern);
            o.extend(seq(&res.angles));
            o
        });
        self.emit("fdtd::nf2ff_2d", Class::C, &ins, &outs);
    }
}

fn emit_1c(r: &mut Rows) {
    let tm = Mode2d::Tm;
    // fdtd.rs's tests and its bench.
    let run60 = [Op2::Run(100, Some((30, 30)), 3e9, 1.0)];
    r.nf2ff2(5, 36, 3e9, 60, 60, 1e-4, tm, &run60);
    r.nf2ff2(6, 36, 1e9, 10, 10, 1e-4, tm, &[]);
    r.nf2ff2(2, 0, 1e9, 10, 10, 1e-4, tm, &[]);
    r.nf2ff2(3, 8, 1e9, 30, 30, 1e-4, Mode2d::Te, &[]);
    for (m, na, fr) in [(0, 8, 1e9), (3, 8, 0.0), (3, 8, f64::NAN), (3, 4, 1e9)] {
        r.nf2ff2(m, na, fr, 30, 30, 1e-4, tm, &[]);
    }
    let bench = [Op2::Run(60, Some((40, 40)), 3e10, 1.0)];
    r.nf2ff2(5, 180, 3e10, 80, 80, 1e-3, tm, &bench);
    // The validation: the frequency's specials, the margin at 0, 1, its last valid value and nx/2 or ny/2 on a
    // non-square grid, angle counts from 0 to 360.
    let run40 = [Op2::Source(20, 12, 1.0), Op2::Step(30)];
    for &v in &specials() {
        r.nf2ff2(4, 12, v, 40, 24, 1e-4, tm, &run40);
    }
    for m in [0, 1, 2, 11, 12, 19, 20] {
        r.nf2ff2(m, 12, 1e9, 40, 24, 1e-4, tm, &run40);
    }
    for na in [0, 1, 2, 3, 7, 36, 90, 360] {
        r.nf2ff2(3, na, 1e9, 40, 24, 1e-4, tm, &run40);
    }
    r.nf2ff2(4, 12, 1e9, 40, 24, 1e-4, Mode2d::Te, &run40);
    // Transforms: a driven point source with and without the CPML, a dielectric, the smallest contours, frequencies
    // from 1 MHz to 300 GHz.
    for (nx, ny, m) in [(7, 7, 1), (8, 9, 1), (9, 9, 2), (20, 20, 3), (33, 21, 4), (48, 48, 6)] {
        let ops = [Op2::Source(nx / 2, ny / 2, 1.0), Op2::Step(12), Op2::Source(nx / 3, ny / 3, -0.5), Op2::Step(5)];
        for na in [1, 4, 24] {
            r.nf2ff2(m, na, 2e9, nx, ny, 1e-4, tm, &ops);
        }
    }
    for fr in [1e6, 1e8, 1e9, 3e9, 1e10, 3e11] {
        r.nf2ff2(5, 36, fr, 60, 60, 1e-4, tm, &run60);
    }
    let cpml = [Op2::Cpml(8), Op2::Run(120, Some((30, 30)), 3e9, 1.0)];
    r.nf2ff2(10, 72, 3e9, 60, 60, 1e-4, tm, &cpml);
    let die = [Op2::Eps(20, 20, 40, 40, 4.0), Op2::Run(80, Some((30, 30)), 3e9, 1.0)];
    r.nf2ff2(5, 36, 3e9, 60, 60, 1e-4, tm, &die);
    // At the core's range edges (ADR 0002 row 18, M5 part 2d-6): k below 2^-1022 and near DBL_MAX, fields of 10^±300
    // and past DBL_MAX in their sums, a non-finite field, dx of 10^±300.
    for fr in [f64::from_bits(1), 1e-310, 1e-300, f64::MIN_POSITIVE, 1e300, f64::MAX] {
        r.nf2ff2(3, 8, fr, 20, 20, 1e-4, tm, &[Op2::Source(10, 10, 1.0), Op2::Step(15)]);
    }
    for v in [1e300, 1e-300, 1e-310, f64::MAX] {
        let ops = [Op2::Source(10, 10, v), Op2::Step(15), Op2::Source(4, 4, v)];
        r.nf2ff2(3, 8, 1e9, 20, 20, 1e-4, tm, &ops);
    }
    for v in [f64::MAX, 1e-305, f64::INFINITY, f64::NAN] {
        let ops = [Op2::Fill(0, 1.0), Op2::Source(3, 3, v), Op2::Fill(2, 0.25)];
        r.nf2ff2(3, 8, 1e9, 20, 20, 1e-4, tm, &ops);
    }
    for dx in [1e-300, 1e300, 1e-310] {
        r.nf2ff2(3, 8, 1e9, 20, 20, dx, tm, &[Op2::Source(10, 10, 1.0), Op2::Step(15)]);
        r.nf2ff2(3, 8, 1e-300, 20, 20, dx, tm, &[Op2::Source(10, 10, 1.0), Op2::Step(15)]);
    }
    // Seeded transforms on grids from 7 to 40 a side.
    let mut rng = Rng::new(0x1c_2ff2);
    for _ in 0..30 {
        let nx = usize::try_from(rng.int(7, 40)).unwrap_or(7);
        let ny = usize::try_from(rng.int(7, 40)).unwrap_or(7);
        let mmax = i64::try_from(nx.min(ny) / 2).unwrap_or(2) - 1;
        let m = usize::try_from(rng.int(1, mmax.max(1))).unwrap_or(1);
        let na = usize::try_from(rng.int(1, 48)).unwrap_or(1);
        let fr = rng.log_uniform(6.0, 11.0);
        let k = usize::try_from(rng.int(1, 40)).unwrap_or(1);
        let ops = [Op2::Source(nx / 2, ny / 2, rng.signed_log_uniform(-2.0, 1.0)), Op2::Step(k)];
        r.nf2ff2(m, na, fr, nx, ny, rng.log_uniform(-5.0, -2.0), tm, &ops);
    }
}

pub fn emit() -> Module {
    let mut r = Rows { m: Module::new("fdtd_2d3d"), keys: BTreeSet::new() };
    emit_1a(&mut r);
    emit_1b1(&mut r);
    emit_1b2(&mut r);
    emit_1c(&mut r);
    r.m
}
