//! `fdtd.tsv` — rust-old/src/fdtd.rs's FDTD core (M5). Part 1a: the 1D grid `Fdtd1d`.
//!
//! A row is a **script** (port plan §6 step 5's trajectory format): its inputs are the grid, `i:<num_cells>` and
//! `f:<dx>`, then the operations applied to it in order, each an opcode and its arguments:
//! - `i:1 i:<start> i:<end> f:<eps_r>`: `set_permittivity`;
//! - `i:2 i:<start> i:<end> f:<mu_r>`: `set_permeability`;
//! - `i:3 i:<cell> f:<value>`: `add_source`;
//! - `i:4 i:<k>`: `step_once`, k times;
//! - `i:5 i:<steps> i:<cell>|n: f:<frequency> f:<amplitude>`: `run`, `n:` for `None`.
//!
//! A usize is its 64-bit word (usize::MAX is `i:-1`). The outputs: `E:<Variant>` and the Display text when `new`
//! fails; otherwise the `E:` and text of each operation that failed, in order (an `Ok(())` adds nothing), then the
//! state: `i:<step>`, `f:<dt>`, `f:<time()>`, `f:<total_energy()>`, `f:<abc_left>`, `f:<abc_right>`, the digests
//! (`i:`) of `e_field`, `h_field`, `permittivity`, `permeability`, `e_coeff` and `h_coeff`, and the probe cells
//! `f:` E then H at 0, 1, n/2, n − 2 and n − 1. The private fields are read from the grid's `Debug` text, whose
//! shortest round-trip digits parse back to the same bits (a NaN's sign and payload aside, which the digest folds to
//! one NaN anyway). A digest is the harness's cross-target digest (digest v2, port-parity.md §9, `_par_dg_step`):
//! y = (h ⊕ w)·0x100000001b3, h = y ⊕ (y >> 29) with an arithmetic shift, from 0xcbf29ce484222325, every NaN as
//! 0x7ff8….
//!
//! Classes: T, every row (port plan §7): the state is compared exactly, frame by frame. A `run` that injects calls
//! sin, so `validate.py` writes the correctly rounded readings of those rows to `fdtd.cr.tsv`.
//!
//! The grids: `new` over cell counts and every special spacing; `set_permittivity` and `set_permeability` over the
//! specials and every kind of range (empty, reversed, clamped, past the grid, usize::MAX); `add_source` over the
//! specials, at every edge cell and outside; trajectories at frames {1, 2, 10, 100, 1000} in vacuum, ε_r = 4, μ_r = 9,
//! layered media, two- and three-cell grids, non-finite and huge sources, signed zeros; `run` with the Rust tests' and
//! benches' sources, None and out-of-grid cells, special and huge frequencies and amplitudes; and seeded scenarios.
//!
//! M5 part 1b appends the sources. The pulses (`fdtd::gaussian_pulse`, `ricker_wavelet`, `differentiated_gaussian`,
//! `modulated_gaussian_pulse`) are ordinary class-C rows: every special for t, t₀ and τ (and for the frequency on
//! typical pulses), the Rust tests' points, seeded pulses from 1 fs to 1 µs, seeded wide values, and the range's edges
//! ((t − t₀)/τ past DBL_MAX or subnormal, u² past either end, exp's underflow). The 1D scripts gain opcode 6, k steps
//! each after adding gaussian_pulse(time(), t₀, τ) at a cell: the Rust tests' Mur absorption in vacuum and ε_r = 4,
//! a pulse in vacuum, special widths and delays. A `TfsfSource` script is `i:<nx> i:<boundary_x> f:<dx> f:<dt>`, built
//! by the row's constructor (`TfsfSource::new` rows by `new`, the others by `try_new`), then `i:1 f:<v>` step(v),
//! `i:2 i:<k> i:<s0> f:<dt> f:<t0> f:<tau>` k steps driven by gaussian_pulse((s₀ + j)·dt, t₀, τ), `i:3 i:<b>` and
//! `i:4 i:<c>` the public `boundary_x` and `source_cell` written. Its state: nx, the boundary, the source cell,
//! `incident_e()`, `incident_h()`, the three private coefficients (from `Debug`), both arrays' digests and the probes.
//! The grid: every constructor case (nx below and at 3, the boundary at and past its range's ends, the specials for dx
//! and dt, dt straddling dx·(1 + 10⁻¹²)/c ulp by ulp), the Rust tests' drives, a unit impulse and non-finite sources,
//! the written fields, and seeded sources. Rust's `new` clamps what `try_new` refuses; ADR 0002 row 3 maps those rows.
//!
//! M5 part 1c appends the media and the subpixel averages, ordinary class-A rows (+ − × ÷ in Rust's order). A
//! material's inputs are its fields (then dt, and dx for the lossy one), its outputs its fields, its error, or its
//! coefficient tuple. The grid: `new` and `try_new` over every special crossed (Lorentz's `new` at γ = 1), the Rust
//! tests' cases and Debye's ε_s ≥ ε_∞ straddled ulp by ulp; the three presets; the coefficients over every special
//! for each argument (Lorentz's unread ε_∞ at 2), the Rust tests' cases, the presets across FDTD time steps, 1,000
//! seeded typical media and 500 seeded wide ones; the subpixel averages over every special crossed, fill fractions at
//! and one ulp past [0, 1]'s ends and at −0, and seeded media.

use std::collections::BTreeSet;

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::fdtd::{
    DebyeMaterial, DrudeMaterial, Fdtd1d, LorentzMaterial, LossyMaterial, TfsfSource, differentiated_gaussian,
    gaussian_pulse, modulated_gaussian_pulse, ricker_wavelet, subpixel_permittivity_parallel,
    subpixel_permittivity_perpendicular, subpixel_permittivity_perpendicular_checked,
};

/// One operation of a script.
#[derive(Clone, Copy)]
enum Op {
    Eps(usize, usize, f64),
    Mu(usize, usize, f64),
    Source(usize, f64),
    Step(usize),
    Run(usize, Option<usize>, f64, f64),
    Pulse(usize, usize, f64, f64),
}

fn u(x: usize) -> Tok {
    i(i64::from_ne_bytes(x.to_ne_bytes()))
}

impl Op {
    fn toks(self) -> Vec<Tok> {
        match self {
            Self::Eps(a, b, x) => vec![i(1), u(a), u(b), f(x)],
            Self::Mu(a, b, x) => vec![i(2), u(a), u(b), f(x)],
            Self::Source(c, x) => vec![i(3), u(c), f(x)],
            Self::Step(k) => vec![i(4), u(k)],
            Self::Run(k, c, fr, am) => vec![i(5), u(k), c.map_or(Tok::None, u), f(fr), f(am)],
            Self::Pulse(k, c, t0, tau) => vec![i(6), u(k), u(c), f(t0), f(tau)],
        }
    }
}

/// FNV-1a with the harness's v2 mix over the words of xs, every NaN as 0x7ff8000000000000: y = (h ⊕ w)·0x100000001b3,
/// then h = y ⊕ (y >> 29) with the shift arithmetic (sign-filling), as Cyrius's `>>>` in `_par_dg_step`.
fn digest(xs: &[f64]) -> i64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &x in xs {
        let w = if x.is_nan() { 0x7ff8_0000_0000_0000 } else { x.to_bits() };
        let y = (h ^ w).wrapping_mul(0x0100_0000_01b3);
        let sh = i64::from_ne_bytes(y.to_ne_bytes()) >> 29;
        h = y ^ u64::from_ne_bytes(sh.to_ne_bytes());
    }
    i64::from_ne_bytes(h.to_ne_bytes())
}

/// The text of field `name` in a `Debug` rendering: a scalar up to the next `, ` or ` }`, or a `[..]` list.
fn dbg_field<'a>(d: &'a str, name: &str) -> &'a str {
    let key = format!(" {name}: ");
    let at = d.find(&key).unwrap_or_else(|| panic!("no field {name} in the Debug text")) + key.len();
    let rest = &d[at..];
    if rest.starts_with('[') {
        let end = rest.find(']').expect("a closing bracket");
        &rest[1..end]
    } else {
        let end = rest.find([',', ' ']).expect("a field end");
        &rest[..end]
    }
}

fn dbg_f64(d: &str, name: &str) -> f64 {
    dbg_field(d, name).parse().expect("a Debug f64")
}

fn dbg_vec(d: &str, name: &str) -> Vec<f64> {
    let body = dbg_field(d, name);
    if body.is_empty() {
        return Vec::new();
    }
    body.split(", ").map(|x| x.parse().expect("a Debug f64")).collect()
}

/// The state tokens of a grid (module comment).
fn state(sim: &Fdtd1d) -> Vec<Tok> {
    let d = format!("{sim:?}");
    let n = sim.num_cells;
    let mut out = vec![
        u(sim.step),
        f(sim.dt),
        f(sim.time()),
        f(sim.total_energy()),
        f(dbg_f64(&d, "abc_left")),
        f(dbg_f64(&d, "abc_right")),
        i(digest(&sim.e_field)),
        i(digest(&sim.h_field)),
    ];
    for name in ["permittivity", "permeability", "e_coeff", "h_coeff"] {
        let v = dbg_vec(&d, name);
        assert_eq!(v.len(), n, "{name} has num_cells values");
        out.push(i(digest(&v)));
    }
    let probes = [0, 1, n / 2, n - 2, n - 1];
    for p in probes {
        out.push(f(sim.e_field[p]));
    }
    for p in probes {
        out.push(f(sim.h_field[p]));
    }
    out
}

/// A script being built: its grid (None when `new` failed), its input tokens and its operations' errors so far.
struct Script {
    sim: Result<Fdtd1d, Vec<Tok>>,
    ins: Vec<Tok>,
    errs: Vec<Tok>,
}

impl Script {
    fn new(n: usize, dx: f64) -> Self {
        let sim = Fdtd1d::new(n, dx).map_err(|e| vec![err(&e), s(e.to_string())]);
        Self { sim, ins: vec![u(n), f(dx)], errs: Vec::new() }
    }

    fn apply(&mut self, op: Op) {
        self.ins.extend(op.toks());
        let sim = self.sim.as_mut().expect("operations only on a grid that was made");
        let r = match op {
            Op::Eps(a, b, x) => sim.set_permittivity(a, b, x),
            Op::Mu(a, b, x) => sim.set_permeability(a, b, x),
            Op::Source(c, x) => {
                sim.add_source(c, x);
                Ok(())
            }
            Op::Step(k) => {
                for _ in 0..k {
                    sim.step_once();
                }
                Ok(())
            }
            Op::Run(k, c, fr, am) => {
                sim.run(k, c, fr, am);
                Ok(())
            }
            Op::Pulse(k, c, t0, tau) => {
                for _ in 0..k {
                    let v = gaussian_pulse(sim.time(), t0, tau);
                    sim.add_source(c, v);
                    sim.step_once();
                }
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
                out.extend(state(sim));
                out
            }
        }
    }
}

/// The module's rows, each (function, inputs) once.
struct Rows {
    m: Module,
    keys: BTreeSet<String>,
}

impl Rows {
    fn emit(&mut self, func: &str, sc: &Script) {
        self.emit_toks(func, Class::T, &sc.ins, &sc.outs());
    }

    /// A row of any shape, its (function, inputs) checked unique.
    fn emit_toks(&mut self, func: &str, cls: Class, ins: &[Tok], outs: &[Tok]) {
        let key = format!("{func} {ins:?}");
        assert!(self.keys.insert(key), "a repeated {func} row");
        self.m.row(func, cls, ins, outs);
    }

    /// A script of `setup`, then frames: at each frame, `advance` brings the grid from the previous frame to it.
    fn frames(&mut self, func: &str, n: usize, dx: f64, setup: &[Op], frames: &[usize], advance: impl Fn(usize) -> Op) {
        let mut sc = Script::new(n, dx);
        for &op in setup {
            sc.apply(op);
        }
        let mut at = 0;
        for &fr in frames {
            sc.apply(advance(fr - at));
            at = fr;
            self.emit(func, &sc);
        }
    }

    /// One script with no frames: `setup` applied, one row.
    fn one(&mut self, func: &str, n: usize, dx: f64, setup: &[Op]) {
        let mut sc = Script::new(n, dx);
        for &op in setup {
            sc.apply(op);
        }
        self.emit(func, &sc);
    }
}

const FRAMES: [usize; 5] = [1, 2, 10, 100, 1000];

/// 1D spacings: the specials, typical grids, and the range's ends.
fn spacings() -> Vec<f64> {
    let mut v = specials();
    v.extend([1e-3, 1e-4, 2.5e-2, 1e-300, 1e300, f64::MIN_POSITIVE * 3.0, f64::from_bits(0x000f_ffff_ffff_ffff)]);
    v
}

/// Material values: the specials, typical media, and wide values.
fn media() -> Vec<f64> {
    let mut v = specials();
    v.extend([4.0, 2.25, 0.5, 9.0, 81.0, 1e-300, 1e300, 1e-30, 1e30]);
    v
}

fn emit_1a(r: &mut Rows) {
    // new: every cell count from 0 to a large grid, against every spacing.
    for n in [0_usize, 1, 2, 3, 4, 7, 100, 1000, 65_537] {
        for dx in spacings() {
            r.one("Fdtd1d::new", n, dx, &[]);
        }
    }

    // set_permittivity / set_permeability: every value on an inner range, every range with ε_r = 4, then each kind
    // of range touching the Mur coefficients' nodes.
    let ranges: [(usize, usize); 13] = [
        (0, 10),
        (3, 7),
        (0, 1),
        (1, 2),
        (0, 2),
        (8, 9),
        (8, 10),
        (9, 10),
        (5, 5),
        (7, 3),
        (0, usize::MAX),
        (usize::MAX, 5),
        (20, 30),
    ];
    for (func, mk) in [
        ("Fdtd1d::set_permittivity", Op::Eps as fn(usize, usize, f64) -> Op),
        ("Fdtd1d::set_permeability", Op::Mu as fn(usize, usize, f64) -> Op),
    ] {
        for x in media() {
            r.one(func, 10, 1e-3, &[mk(3, 7, x)]);
            r.one(func, 10, 1e-3, &[mk(0, 10, x)]);
        }
        for &(a, b) in &ranges {
            if (a, b) != (3, 7) && (a, b) != (0, 10) {
                r.one(func, 10, 1e-3, &[mk(a, b, 4.0)]);
            }
            r.one(func, 10, 1e-3, &[mk(a, b, 2.25), mk(0, 1, f64::NAN)]);
        }
        // Two cells: both ends' coefficients read node 0 or 1.
        for (a, b) in [(0, 1), (1, 2), (0, 2)] {
            r.one(func, 2, 1e-3, &[mk(a, b, 9.0)]);
        }
        // Repeated and overlapping settings, and the last one refused.
        r.one(func, 10, 1e-3, &[mk(0, 10, 4.0), mk(2, 6, 2.0), mk(5, 9, 0.5), mk(0, 10, -1.0)]);
        // At the range's ends: a spacing whose coefficients under- or overflow.
        for dx in [1e-300, 1e300, f64::MIN_POSITIVE, f64::MAX] {
            r.one(func, 10, dx, &[mk(0, 10, 1e300)]);
            r.one(func, 10, dx, &[mk(0, 10, 1e-300)]);
        }
    }

    // add_source: every value at every edge cell, outside the grid, and twice.
    for x in media() {
        for c in [0_usize, 1, 5, 9] {
            r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(c, x)]);
        }
    }
    for c in [10_usize, 11, 1 << 40, usize::MAX - 1, usize::MAX] {
        r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(c, 1.0)]);
    }
    r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(5, 1.0), Op::Source(5, 1.0), Op::Source(5, -2.0)]);
    r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(5, f64::MAX), Op::Source(5, f64::MAX)]);
    r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(5, f64::INFINITY), Op::Source(5, f64::NEG_INFINITY)]);
    r.one("Fdtd1d::add_source", 10, 1e-3, &[Op::Source(5, -0.0), Op::Source(6, 0.0)]);

    // step_once trajectories.
    let step = |k| Op::Step(k);
    let st = "Fdtd1d::step_once";
    r.frames(st, 200, 1e-3, &[Op::Source(100, 1.0)], &FRAMES, step);
    r.frames(st, 200, 1e-3, &[Op::Eps(0, 200, 4.0), Op::Source(100, 1.0)], &FRAMES, step);
    r.frames(st, 150, 1e-3, &[Op::Mu(0, 150, 9.0), Op::Source(75, 1.0)], &FRAMES, step);
    r.frames(st, 200, 1e-3, &[Op::Eps(120, 200, 4.0), Op::Mu(50, 80, 2.0), Op::Source(30, 1.0)], &FRAMES, step);
    r.frames(st, 200, 1e-4, &[Op::Eps(0, 100, 2.25), Op::Eps(100, 200, 81.0), Op::Source(60, -3.5)], &FRAMES, step);
    r.frames(st, 100, 1e-3, &[Op::Source(50, 1.0), Op::Source(51, -1.0)], &[1, 10, 100, 1000, 5000], step);
    r.frames(st, 2, 1e-3, &[Op::Source(1, 1.0)], &[1, 2, 3, 10, 100], step);
    r.frames(st, 2, 1e-3, &[Op::Eps(0, 2, 4.0), Op::Source(0, 1.0)], &[1, 2, 10], step);
    r.frames(st, 3, 1e-3, &[Op::Source(1, 1.0)], &[1, 2, 3, 10, 100], step);
    r.frames(st, 4, 1e-3, &[Op::Mu(0, 4, 4.0), Op::Source(2, 1.0)], &[1, 2, 10, 100], step);
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX, 1e300, 1e-300, f64::from_bits(1), -0.0] {
        r.frames(st, 20, 1e-3, &[Op::Source(10, x)], &[1, 2, 10, 100], step);
    }
    r.frames(st, 20, 1e-3, &[Op::Source(0, 1.0)], &[1, 2, 10], step);
    r.frames(st, 20, 1e-3, &[Op::Source(19, 1.0)], &[1, 2, 10], step);
    r.frames(st, 10, 1e-300, &[Op::Source(5, 1.0)], &[1, 2, 10], step);
    r.frames(st, 10, 1e300, &[Op::Source(5, 1.0)], &[1, 2, 10], step);
    r.frames(st, 10, f64::MIN_POSITIVE, &[Op::Source(5, 1.0)], &[1, 2, 10], step);
    r.frames(st, 10, 1e-3, &[Op::Step(0)], &[0, 1], step);

    // run: the Rust tests' and benches' sources, None and out-of-grid cells, special frequencies and amplitudes.
    let rn = "Fdtd1d::run";
    let run_at = |c: Option<usize>, fr: f64, am: f64| move |k| Op::Run(k, c, fr, am);
    r.frames(rn, 200, 1e-3, &[], &FRAMES, run_at(Some(100), 1e9, 1.0));
    r.frames(rn, 500, 1e-4, &[], &[1, 10, 100, 200], run_at(Some(50), 3e9, 1.0));
    r.frames(rn, 500, 1e-4, &[Op::Eps(0, 500, 4.0)], &[1, 10, 100, 200], run_at(Some(50), 3e9, 1.0));
    r.frames(rn, 500, 1e-4, &[], &[1, 10, 100], run_at(Some(250), 1e9, 1.0));
    r.frames(rn, 40, 1e-3, &[Op::Source(20, 1.0)], &[1, 10, 30], run_at(None, 1e9, 1.0));
    r.frames(rn, 40, 1e-3, &[Op::Source(20, 1.0)], &[1, 10, 30], run_at(Some(40), 1e9, 1.0));
    r.frames(rn, 40, 1e-3, &[Op::Source(20, 1.0)], &[1, 10], run_at(Some(usize::MAX), 1e9, 1.0));
    r.frames(rn, 40, 1e-3, &[], &[0, 1, 2], run_at(Some(0), 1e10, 2.0));
    r.frames(rn, 40, 1e-3, &[], &[1, 2, 10], run_at(Some(39), 1e10, 2.0));
    for fr in [0.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e9, 1e20, 1e300, f64::MAX, 1e-300] {
        r.frames(rn, 30, 1e-3, &[], &[1, 2, 10], run_at(Some(15), fr, 1.0));
    }
    for am in [0.0, -0.0, f64::NAN, f64::INFINITY, -1.0, f64::MAX, 1e-300] {
        r.frames(rn, 30, 1e-3, &[], &[1, 2, 10], run_at(Some(15), 1e9, am));
    }

    // Seeded scenarios: random media, sources and drives.
    let mut rng = Rng::new(0x0fd7_d1a0);
    for k in 0..16 {
        let n = usize::try_from(rng.int(2, 300)).unwrap_or(2);
        let dx = rng.log_uniform(-6.0, 0.0);
        let mut setup = Vec::new();
        for _ in 0..rng.int(0, 3) {
            let a = usize::try_from(rng.int(0, 300)).unwrap_or(0);
            let b = usize::try_from(rng.int(0, 320)).unwrap_or(0);
            let x = rng.uniform(1.0, 16.0);
            setup.push(if rng.int(0, 1) == 0 { Op::Eps(a, b, x) } else { Op::Mu(a, b, x) });
        }
        for _ in 0..rng.int(1, 3) {
            let c = usize::try_from(rng.int(0, 310)).unwrap_or(0);
            setup.push(Op::Source(c, rng.signed_log_uniform(-3.0, 3.0)));
        }
        if k % 2 == 0 {
            r.frames(st, n, dx, &setup, &[1, 10, 100, 1000], step);
        } else {
            let cell = usize::try_from(rng.int(0, 300)).unwrap_or(0);
            let fr = rng.log_uniform(6.0, 12.0);
            let am = rng.signed_log_uniform(-2.0, 2.0);
            r.frames(rn, n, dx, &setup, &[1, 10, 100, 1000], run_at(Some(cell), fr, am));
        }
    }
}

// ── part 1b: the pulses and the TF/SF source ──────────────────────────────────────────────────────────────────────

/// One operation of a TF/SF script.
#[derive(Clone, Copy)]
enum TfOp {
    Step(f64),
    Pulse(usize, usize, f64, f64, f64),
    SetBoundary(usize),
    SetSource(usize),
}

impl TfOp {
    fn toks(self) -> Vec<Tok> {
        match self {
            Self::Step(v) => vec![i(1), f(v)],
            Self::Pulse(k, s0, dt, t0, tau) => vec![i(2), u(k), u(s0), f(dt), f(t0), f(tau)],
            Self::SetBoundary(b) => vec![i(3), u(b)],
            Self::SetSource(c) => vec![i(4), u(c)],
        }
    }
}

/// A TF/SF source's state tokens (module comment).
fn tf_state(src: &TfsfSource) -> Vec<Tok> {
    let d = format!("{src:?}");
    let n = src.inc_e.len();
    assert_eq!(src.inc_h.len(), n, "inc_e and inc_h have nx values");
    let mut out = vec![
        u(n),
        u(src.boundary_x),
        u(src.source_cell),
        f(src.incident_e()),
        f(src.incident_h()),
        f(dbg_f64(&d, "e_coeff_1d")),
        f(dbg_f64(&d, "h_coeff_1d")),
        f(dbg_f64(&d, "abc_coeff")),
        i(digest(&src.inc_e)),
        i(digest(&src.inc_h)),
    ];
    let probes = [0, 1, n / 2, n - 2, n - 1];
    for p in probes {
        out.push(f(src.inc_e[p]));
    }
    for p in probes {
        out.push(f(src.inc_h[p]));
    }
    out
}

/// A TF/SF script: the constructor (`try_new` when checked, else `new`), then operations.
struct TfScript {
    src: Result<TfsfSource, Vec<Tok>>,
    ins: Vec<Tok>,
}

impl TfScript {
    fn new(checked: bool, nx: usize, bx: usize, dx: f64, dt: f64) -> Self {
        let src = if checked {
            TfsfSource::try_new(nx, bx, dx, dt).map_err(|e| vec![err(&e), s(e.to_string())])
        } else {
            Ok(TfsfSource::new(nx, bx, dx, dt))
        };
        Self { src, ins: vec![u(nx), u(bx), f(dx), f(dt)] }
    }

    fn apply(&mut self, op: TfOp) {
        self.ins.extend(op.toks());
        let src = self.src.as_mut().expect("operations only on a source that was made");
        match op {
            TfOp::Step(v) => src.step(v),
            TfOp::Pulse(k, s0, dt, t0, tau) => {
                for j in 0..k {
                    src.step(gaussian_pulse((s0 + j) as f64 * dt, t0, tau));
                }
            }
            TfOp::SetBoundary(b) => src.boundary_x = b,
            TfOp::SetSource(c) => src.source_cell = c,
        }
    }

    fn outs(&self) -> Vec<Tok> {
        match &self.src {
            Err(e) => e.clone(),
            Ok(src) => tf_state(src),
        }
    }
}

impl Rows {
    fn tf(&mut self, func: &str, sc: &TfScript) {
        self.emit_toks(func, Class::T, &sc.ins, &sc.outs());
    }

    /// A Gaussian-driven TF/SF scenario (try_new): a frame after each of `frames` steps, the drive's time s·dt.
    fn tf_frames(&mut self, func: &str, g: (usize, usize, f64, f64), setup: &[TfOp], drive: (f64, f64),
                 frames: &[usize]) {
        let (nx, bx, dx, dt) = g;
        let mut sc = TfScript::new(true, nx, bx, dx, dt);
        for &op in setup {
            sc.apply(op);
        }
        let mut at = 0;
        for &fr in frames {
            sc.apply(TfOp::Pulse(fr - at, at, dt, drive.0, drive.1));
            at = fr;
            self.tf(func, &sc);
        }
    }

    /// A pulse row: `fdtd::<name>` of the inputs, class C (exp and sin).
    fn pulse(&mut self, name: &str, ins: &[f64]) {
        let v = match (name, ins) {
            ("gaussian_pulse", &[t, t0, tau]) => gaussian_pulse(t, t0, tau),
            ("ricker_wavelet", &[t, t0, tau]) => ricker_wavelet(t, t0, tau),
            ("differentiated_gaussian", &[t, t0, tau]) => differentiated_gaussian(t, t0, tau),
            ("modulated_gaussian_pulse", &[t, t0, tau, fr]) => modulated_gaussian_pulse(t, t0, tau, fr),
            _ => panic!("no pulse {name} of {} inputs", ins.len()),
        };
        let toks: Vec<Tok> = ins.iter().map(|&x| f(x)).collect();
        self.emit_toks(&format!("fdtd::{name}"), Class::C, &toks, &[f(v)]);
    }
}

const PULSES: [&str; 3] = ["gaussian_pulse", "ricker_wavelet", "differentiated_gaussian"];

fn emit_1b_pulses(r: &mut Rows) {
    let sp = specials();
    let fr0 = 1e9;
    // Every special for t, t0 and τ; the modulated pulse at 1 GHz, then every special frequency on a typical pulse.
    for &t in &sp {
        for &t0 in &sp {
            for &tau in &sp {
                for name in PULSES {
                    r.pulse(name, &[t, t0, tau]);
                }
                r.pulse("modulated_gaussian_pulse", &[t, t0, tau, fr0]);
            }
        }
    }
    for &fr in &sp {
        for (t, t0, tau) in [(1.5e-9, 1e-9, 1e-10), (1e-9, 1e-9, 1e-10), (0.0, 3e-9, 5e-10)] {
            if fr != fr0 {
                r.pulse("modulated_gaussian_pulse", &[t, t0, tau, fr]);
            }
        }
    }
    // The Rust tests' points.
    let (t0, tau) = (1e-9, 1e-10);
    let mut pts = vec![(t0 - 0.5e-10, t0, tau), (t0 + 0.5e-10, t0, tau), (t0 + 5.0 * tau, t0, tau), (0.0, 0.0, 1e-20)];
    for k in -40..=40 {
        pts.push((2e-9 + f64::from(k) * 0.1 * tau, 2e-9, tau));
    }
    pts.push((t0 - tau / 2f64.sqrt(), t0, tau));
    pts.push((t0 + tau / 2f64.sqrt(), t0, tau));
    pts.push((t0 + tau * 1.5f64.sqrt(), t0, tau));
    for (t, t0, tau) in pts {
        for name in PULSES {
            r.pulse(name, &[t, t0, tau]);
        }
        r.pulse("modulated_gaussian_pulse", &[t, t0, tau, 1e9]);
    }
    // Seeded typical pulses: τ from 1 fs to 1 µs, t₀ a few τ, t within ±8τ of it; carriers to 1 THz.
    let mut rng = Rng::new(0x0b1b_5ace);
    for _ in 0..2000 {
        let tau = rng.log_uniform(-15.0, -6.0);
        let t0 = tau * rng.uniform(0.0, 8.0);
        let t = t0 + tau * rng.uniform(-8.0, 8.0);
        for name in PULSES {
            r.pulse(name, &[t, t0, tau]);
        }
        r.pulse("modulated_gaussian_pulse", &[t, t0, tau, rng.log_uniform(6.0, 12.0)]);
    }
    // Seeded wide values, and the edges: (t − t₀)/τ past DBL_MAX or under 2^-1022, u² past either end, exp's
    // underflow near u² = 745.
    for _ in 0..1000 {
        let t = rng.signed_log_uniform(-300.0, 300.0);
        let t0 = rng.signed_log_uniform(-300.0, 300.0);
        let tau = rng.signed_log_uniform(-300.0, 300.0);
        for name in PULSES {
            r.pulse(name, &[t, t0, tau]);
        }
        r.pulse("modulated_gaussian_pulse", &[t, t0, tau, rng.signed_log_uniform(-300.0, 300.0)]);
    }
    let edges = [
        (1e300, -1e300, 1e-10),
        (1e308, 0.0, 0.5),
        (f64::MAX, -f64::MAX, f64::MAX),
        (1e-300, 0.0, 1e300),
        (f64::MIN_POSITIVE * 3.0, 0.0, 1.0),
        (1e-160, 0.0, 1.0),
        (1e160, 0.0, 1.0),
        (1e155, 0.0, 1e-1),
        (27.29, 0.0, 1.0),
        (27.3, 0.0, 1.0),
        (38.6, 0.0, 1.0),
        (1.0, 0.5, f64::from_bits(1)),
        (5e-324, -5e-324, 5e-324),
    ];
    for (t, t0, tau) in edges {
        for name in PULSES {
            r.pulse(name, &[t, t0, tau]);
        }
        r.pulse("modulated_gaussian_pulse", &[t, t0, tau, 1e9]);
        r.pulse("modulated_gaussian_pulse", &[t, t0, tau, 1e300]);
    }
}

fn emit_1b(r: &mut Rows) {
    emit_1b_pulses(r);

    // Fdtd1d driven by a Gaussian pulse (opcode 6): the Rust tests' Mur absorption in vacuum and ε_r = 4, a pulse
    // in vacuum at the frames, and special widths and delays.
    let st = "Fdtd1d::step_once";
    for eps in [1.0_f64, 4.0] {
        let dt = Fdtd1d::new(600, 1e-4).expect("a grid").dt;
        let (tau, t0) = (40.0 * dt, 240.0 * dt);
        let travel = (300.0 / (0.5 / eps.sqrt())) as usize;
        let mut sc = Script::new(600, 1e-4);
        sc.apply(Op::Eps(0, 600, eps));
        let mut at = 0;
        for fr in [1, 10, 100, 480] {
            sc.apply(Op::Pulse(fr - at, 300, t0, tau));
            at = fr;
            r.emit(st, &sc);
        }
        for fr in [481, 480 + travel / 2, 480 + travel, 480 + travel + travel / 3] {
            sc.apply(Op::Step(fr - at));
            at = fr;
            r.emit(st, &sc);
        }
    }
    let dt = Fdtd1d::new(200, 1e-3).expect("a grid").dt;
    r.frames(st, 200, 1e-3, &[], &FRAMES, |k| Op::Pulse(k, 100, 40.0 * dt, 10.0 * dt));
    for (t0, tau) in [(0.0, 0.0), (0.0, f64::NAN), (f64::NAN, dt), (f64::INFINITY, dt), (0.0, 1e-300), (0.0, 1e300)] {
        r.frames(st, 20, 1e-3, &[], &[1, 2, 10], |k| Op::Pulse(k, 10, t0, tau));
    }

    // TF/SF: try_new and new over every kind of argument: nx below and at 3, boundary_x at and past its range's ends,
    // the specials for dx and dt, and dt straddling the 1D Courant limit dx·(1 + 10⁻¹²)/c ulp by ulp.
    let c = bijli::field::SPEED_OF_LIGHT;
    let dx = 1e-4;
    let half = dx / (2.0 * c);
    let lim = dx * (1.0 + 1e-12) / c;
    let mut dts = specials();
    dts.extend([half, dx / c, 1e-13, 1e-11, 1.01 * dx / c]);
    for k in -3_i64..=3 {
        let bits = lim.to_bits().wrapping_add_signed(k);
        dts.push(f64::from_bits(bits));
    }
    let mut cases: Vec<(usize, usize, f64, f64)> = Vec::new();
    for nx in [0_usize, 1, 2, 3, 4, 10, 100] {
        for bx in [0_usize, 1, 2, 5, 50, 99, 100, 101, 500, usize::MAX] {
            cases.push((nx, bx, dx, half));
        }
    }
    for &dxx in &specials() {
        cases.push((100, 50, dxx, half));
        cases.push((100, 50, dxx, 1e-13));
        cases.push((100, 50, dxx, f64::NAN));
    }
    for &dtt in &dts {
        cases.push((100, 50, dx, dtt));
    }
    cases.extend([
        (0, 0, 1e-4, 1e-13),
        (1, 0, 1e-4, 1e-13),
        (100, 500, 1e-4, 1e-13),
        (100, 50, f64::NAN, 1e-13),
        (100, 50, 1e-4, f64::NAN),
        (100, 50, 1e-4, 1e-11),
        (100, 50, -1.0, -1.0),
        (100, 0, 1e-4, 1e-13),
        (100, 50, 1e-4, 1e-13),
        (200, 100, 1e-4, 1e-13),
        (50, 20, 1e-4, 1e-13),
        (3, 1, 1e-300, 1e-310),
        (3, 2, 1e300, 1e290),
        (4, 3, f64::MAX, 1e290),
    ]);
    let mut seen = BTreeSet::new();
    cases.retain(|&(nx, bx, dxx, dtt)| seen.insert((nx, bx, dxx.to_bits(), dtt.to_bits())));
    for &(nx, bx, dxx, dtt) in &cases {
        r.tf("TfsfSource::try_new", &TfScript::new(true, nx, bx, dxx, dtt));
        r.tf("TfsfSource::new", &TfScript::new(false, nx, bx, dxx, dtt));
    }
    // try_new refuses a huge nx's other arguments before allocating (Rust would abort on the allocation).
    for (bx, dxx, dtt) in [(0, 1e-4, 1e-13), (usize::MAX, 1e-4, 1e-13), (5, f64::NAN, 1e-13), (5, 1e-4, 1.0)] {
        r.tf("TfsfSource::try_new", &TfScript::new(true, usize::MAX, bx, dxx, dtt));
    }
    // The Rust test that stepped each clamped grid 2,000 times with gaussian_pulse(step, 30, 10).
    for (nx, bx, dxx, dtt) in [(0, 0, 1e-4, 1e-13), (100, 500, 1e-4, 1e-13), (100, 50, -1.0, -1.0)] {
        let mut sc = TfScript::new(false, nx, bx, dxx, dtt);
        sc.apply(TfOp::Pulse(2000, 0, 1.0, 30.0, 10.0));
        r.tf("TfsfSource::new", &sc);
    }

    // TF/SF steps: the Rust tests' drives at their frames, a unit impulse, every boundary from 1 to 3, the magic time
    // step (S = 1, the Mur coefficient 0), the smallest grids, non-finite and huge sources.
    let ts = "TfsfSource::step";
    let dt13 = 1e-13;
    r.tf_frames(ts, (200, 100, 1e-4, dt13), &[], (20.0 * dt13, 5.0 * dt13), &[1, 2, 10, 50, 100]);
    for bx in 1..4 {
        r.tf_frames(ts, (50, bx, dx, half), &[], (20.0 * half, 6.0 * half), &[1, 10, 60]);
    }
    r.tf_frames(ts, (200, 100, dx, half), &[], (120.0 * half, 20.0 * half), &[1, 10, 100, 1000, 2000]);
    r.tf_frames(ts, (300, 120, dx, half), &[], (90.0 * half, 15.0 * half), &[1, 100, 900]);
    r.tf_frames(ts, (100, 50, dx, dx / c), &[], (30.0 * dx / c, 8.0 * dx / c), &[1, 10, 100, 300]);
    r.tf_frames(ts, (3, 1, dx, half), &[], (2.0 * half, half), &[1, 2, 10]);
    r.tf_frames(ts, (3, 2, dx, half), &[], (2.0 * half, half), &[1, 2, 10]);
    r.tf_frames(ts, (4, 3, dx, half), &[], (2.0 * half, half), &[1, 2, 10]);
    for v in [1.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX, 1e300, 1e-300, f64::from_bits(1)] {
        let mut sc = TfScript::new(true, 100, 50, dx, dt13);
        for k in [1, 1, 8, 90] {
            for _ in 0..k {
                sc.apply(TfOp::Step(if k == 1 { v } else { 0.0 }));
            }
            r.tf(ts, &sc);
        }
    }

    // The public fields written: the boundary and the source cell inside, at the edges, outside and at usize::MAX
    // (Rust's test of mutated fields among them).
    let ti = "TfsfSource::incident_h";
    for b in [0_usize, 1, 2, 48, 49, 50, 51, 1000, usize::MAX] {
        let mut sc = TfScript::new(true, 50, 20, 1e-4, dt13);
        sc.apply(TfOp::SetBoundary(b));
        sc.apply(TfOp::Step(1.0));
        r.tf(ti, &sc);
        sc.apply(TfOp::Pulse(20, 1, dt13, 10.0 * dt13, 3.0 * dt13));
        r.tf(ti, &sc);
    }
    for c in [0_usize, 1, 2, 47, 48, 49, 1000, usize::MAX] {
        let mut sc = TfScript::new(true, 50, 20, 1e-4, dt13);
        sc.apply(TfOp::SetSource(c));
        sc.apply(TfOp::Step(1.0));
        r.tf(ti, &sc);
    }
    let mut sc = TfScript::new(true, 50, 20, 1e-4, dt13);
    sc.apply(TfOp::SetBoundary(1000));
    sc.apply(TfOp::SetSource(1000));
    sc.apply(TfOp::Step(1.0));
    sc.apply(TfOp::SetBoundary(0));
    r.tf(ti, &sc);

    // Seeded sources: random sizes, boundaries, Courant numbers and drives.
    let mut rng = Rng::new(0x7f5f_0b1b);
    for _ in 0..12 {
        let nx = usize::try_from(rng.int(3, 400)).unwrap_or(3);
        let bx = usize::try_from(rng.int(1, i64::try_from(nx).unwrap_or(3) - 1)).unwrap_or(1);
        let dxx = rng.log_uniform(-6.0, 0.0);
        let dtt = dxx / c * rng.uniform(0.05, 1.0);
        let tau = dtt * rng.uniform(2.0, 30.0);
        r.tf_frames(ts, (nx, bx, dxx, dtt), &[], (tau * rng.uniform(3.0, 8.0), tau), &[1, 10, 100, 1000]);
    }
}

// ── part 1c: the media and the subpixel averages ──────────────────────────────────────────────────────────────────

/// A `Result` of fields, or its error.
fn res_fields(r: bijli::error::Result<Vec<f64>>) -> Vec<Tok> {
    match r {
        Ok(v) => v.into_iter().map(f).collect(),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

fn floats(xs: &[f64]) -> Vec<Tok> {
    xs.iter().map(|&x| f(x)).collect()
}

/// Every combination of `k` values from `vals`, in lexicographic order.
fn cross(vals: &[f64], k: usize) -> Vec<Vec<f64>> {
    let mut out = vec![Vec::new()];
    for _ in 0..k {
        out = out.into_iter().flat_map(|p| vals.iter().map(move |&v| [p.clone(), vec![v]].concat())).collect();
    }
    out
}

impl Rows {
    /// A class-A row; a hand-picked case the crossed specials already hold is skipped (the same inputs give the same
    /// outputs).
    fn a(&mut self, func: &str, ins: &[f64], outs: &[Tok]) {
        let toks = floats(ins);
        if self.keys.contains(&format!("{func} {toks:?}")) {
            return;
        }
        self.emit_toks(func, Class::A, &toks, outs);
    }

    fn lossy(&mut self, v: &[f64]) {
        let m = LossyMaterial::new(v[0], v[1]);
        let (ca, cb) = m.update_coefficients(v[2], v[3]);
        self.a("LossyMaterial::update_coefficients", v, &[f(ca), f(cb)]);
    }

    fn drude(&mut self, v: &[f64]) {
        let (a, b, k) = DrudeMaterial::new(v[0], v[1], v[2]).ade_coefficients(v[3]);
        self.a("DrudeMaterial::ade_coefficients", v, &[f(a), f(b), f(k)]);
    }

    fn lorentz(&mut self, v: &[f64]) {
        let (a, b, g) = LorentzMaterial::new(v[0], v[1], v[2], v[3]).ade_coefficients(v[4]);
        self.a("LorentzMaterial::ade_coefficients", v, &[f(a), f(b), f(g)]);
    }

    fn debye(&mut self, v: &[f64]) {
        let (a, b) = DebyeMaterial::new(v[0], v[1], v[2]).ade_coefficients(v[3]);
        self.a("DebyeMaterial::ade_coefficients", v, &[f(a), f(b)]);
    }

    fn subpixel(&mut self, v: &[f64]) {
        let (e1, e2, fr) = (v[0], v[1], v[2]);
        let par = subpixel_permittivity_parallel(e1, e2, fr).map(|x| vec![x]);
        self.a("fdtd::subpixel_permittivity_parallel", v, &res_fields(par));
        self.a("fdtd::subpixel_permittivity_perpendicular", v, &[f(subpixel_permittivity_perpendicular(e1, e2, fr))]);
        let chk = subpixel_permittivity_perpendicular_checked(e1, e2, fr).map(|x| vec![x]);
        self.a("fdtd::subpixel_permittivity_perpendicular_checked", v, &res_fields(chk));
    }
}

fn emit_1c(r: &mut Rows) {
    let sp = specials();
    // The constructors: `new` stores, `try_new` validates (its checks in order: every special crossed).
    for v in cross(&sp, 2) {
        let m = LossyMaterial::new(v[0], v[1]);
        r.a("LossyMaterial::new", &v, &floats(&[m.eps_r, m.conductivity]));
        let t = LossyMaterial::try_new(v[0], v[1]).map(|m| vec![m.eps_r, m.conductivity]);
        r.a("LossyMaterial::try_new", &v, &res_fields(t));
    }
    for v in cross(&sp, 3) {
        let m = DrudeMaterial::new(v[0], v[1], v[2]);
        r.a("DrudeMaterial::new", &v, &floats(&[m.eps_inf, m.omega_p, m.gamma]));
        let t = DrudeMaterial::try_new(v[0], v[1], v[2]).map(|m| vec![m.eps_inf, m.omega_p, m.gamma]);
        r.a("DrudeMaterial::try_new", &v, &res_fields(t));
        let m = DebyeMaterial::new(v[0], v[1], v[2]);
        r.a("DebyeMaterial::new", &v, &floats(&[m.eps_inf, m.eps_s, m.tau]));
        let t = DebyeMaterial::try_new(v[0], v[1], v[2]).map(|m| vec![m.eps_inf, m.eps_s, m.tau]);
        r.a("DebyeMaterial::try_new", &v, &res_fields(t));
    }
    for v in cross(&sp, 4) {
        let t = LorentzMaterial::try_new(v[0], v[1], v[2], v[3])
            .map(|m| vec![m.eps_inf, m.delta_eps, m.omega_0, m.gamma]);
        r.a("LorentzMaterial::try_new", &v, &res_fields(t));
        if v[3].to_bits() == 1.0_f64.to_bits() {
            let m = LorentzMaterial::new(v[0], v[1], v[2], v[3]);
            r.a("LorentzMaterial::new", &v, &floats(&[m.eps_inf, m.delta_eps, m.omega_0, m.gamma]));
        }
    }
    // The Rust tests' constructor cases, the presets, and Debye's ε_s ≥ ε_∞ straddled.
    for v in [[1.0, 1.37e16, 4.05e13], [1.0, 1e16, 1e13], [0.0, 1e16, 1e13], [f64::NAN, 1e16, 1e13], [1.0, -1.0, 1e13],
              [1.0, 1e16, f64::INFINITY], [1.0, 1e15, 1e14], [1.0, -0.0, 0.0]] {
        let t = DrudeMaterial::try_new(v[0], v[1], v[2]).map(|m| vec![m.eps_inf, m.omega_p, m.gamma]);
        r.a("DrudeMaterial::try_new", &v, &res_fields(t));
    }
    for v in [[1.0, 2.0, 1e15, 1e13], [0.0, 2.0, 1e15, 1e13], [1.0, -2.0, 1e15, 1e13], [1.0, 2.0, 0.0, 1e13],
              [1.0, 2.0, 1e15, f64::NAN], [2.0, 2.0, 1e15, 1e14], [1.0, 2.0, -0.0, 1.0], [-1.0, -2.0, 0.0, f64::NAN]] {
        let t = LorentzMaterial::try_new(v[0], v[1], v[2], v[3])
            .map(|m| vec![m.eps_inf, m.delta_eps, m.omega_0, m.gamma]);
        r.a("LorentzMaterial::try_new", &v, &res_fields(t));
    }
    let below = |x: f64| f64::from_bits(x.to_bits() - 1);
    let above = |x: f64| f64::from_bits(x.to_bits() + 1);
    for v in [[4.9, 78.4, 8.27e-12], [0.0, 80.0, 1e-11], [80.0, 2.0, 1e-11], [2.0, 80.0, 0.0], [2.0, f64::NAN, 1e-11],
              [4.9, 78.4, 1e-11], [2.0, 80.0, 1e-11], [80.0, 2.5, 0.0], [2.0, 2.0, f64::INFINITY], [2.0, 2.0, 1e-11],
              [2.0, below(2.0), 1e-11], [2.0, above(2.0), 1e-11], [below(2.0), 2.0, 1e-11]] {
        let t = DebyeMaterial::try_new(v[0], v[1], v[2]).map(|m| vec![m.eps_inf, m.eps_s, m.tau]);
        r.a("DebyeMaterial::try_new", &v, &res_fields(t));
    }
    for v in [[4.0, 0.5], [-1.0, 0.5], [4.0, -0.5], [4.0, f64::NAN], [1.0, 1e10], [0.0, 1.0], [4.0, -0.0]] {
        let t = LossyMaterial::try_new(v[0], v[1]).map(|m| vec![m.eps_r, m.conductivity]);
        r.a("LossyMaterial::try_new", &v, &res_fields(t));
    }
    let presets = [("DrudeMaterial::gold", DrudeMaterial::gold()), ("DrudeMaterial::silver", DrudeMaterial::silver())];
    for (name, m) in presets {
        r.a(name, &[], &floats(&[m.eps_inf, m.omega_p, m.gamma]));
    }
    let w = DebyeMaterial::water();
    r.a("DebyeMaterial::water", &[], &floats(&[w.eps_inf, w.eps_s, w.tau]));

    // The coefficients: every special for each argument (Lorentz's ε_∞, which its coefficients don't read, at 2).
    for v in cross(&sp, 4) {
        r.lossy(&v);
        r.drude(&v);
        r.lorentz(&[2.0, v[0], v[1], v[2], v[3]]);
        r.debye(&v);
    }
    // The Rust tests' cases and the presets across FDTD time steps.
    for v in [[1.0, 0.0, 1e-13, 1e-4], [1.0, 1.0, 1e-13, 1e-4], [1.0, 1e10, 1e-13, 1e-4], [4.0, 0.5, 1e-13, 1e-4]] {
        r.lossy(&v);
    }
    let g = DrudeMaterial::gold();
    let s = DrudeMaterial::silver();
    for dt in [1e-19, 1e-18, 1e-17, 1e-16, 1e-15] {
        r.drude(&[g.eps_inf, g.omega_p, g.gamma, dt]);
        r.drude(&[s.eps_inf, s.omega_p, s.gamma, dt]);
        r.debye(&[w.eps_inf, w.eps_s, w.tau, dt * 1000.0]);
    }
    r.drude(&[1.0, 1e16, 1e13, 1e-16]);
    r.drude(&[1.0, 1e15, 1e14, 1e-17]);
    r.lorentz(&[1.0, 2.0, 1e15, 1e13, 1e-16]);
    r.lorentz(&[2.0, 2.0, 1e15, 1e14, 1e-17]);
    r.debye(&[2.0, 80.0, 1e-11, 1e-13]);
    r.debye(&[4.9, 78.4, 1e-11, 1e-14]);
    // Seeded typical media and FDTD steps, then seeded wide values.
    let mut rng = Rng::new(0x3d1a_0c0c);
    for _ in 0..1000 {
        r.lossy(&[rng.uniform(1.0, 20.0), rng.log_uniform(-6.0, 8.0), rng.log_uniform(-16.0, -10.0),
                  rng.log_uniform(-6.0, -2.0)]);
        r.drude(&[rng.uniform(1.0, 10.0), rng.log_uniform(14.0, 17.0), rng.log_uniform(11.0, 15.0),
                  rng.log_uniform(-19.0, -15.0)]);
        r.lorentz(&[rng.uniform(1.0, 10.0), rng.uniform(0.0, 10.0), rng.log_uniform(13.0, 16.0),
                    rng.log_uniform(11.0, 15.0), rng.log_uniform(-19.0, -15.0)]);
        let ei = rng.uniform(1.0, 10.0);
        r.debye(&[ei, rng.uniform(ei, 100.0), rng.log_uniform(-13.0, -9.0), rng.log_uniform(-16.0, -12.0)]);
    }
    for _ in 0..500 {
        let mut wide = || rng.signed_log_uniform(-300.0, 300.0);
        let v = [wide(), wide(), wide(), wide(), wide()];
        r.lossy(&v[..4]);
        r.drude(&v[..4]);
        r.lorentz(&v);
        r.debye(&v[..4]);
    }

    // The subpixel averages: every special crossed, the fill fraction at and one ulp past [0, 1]'s ends and at −0,
    // the Rust tests' cases, seeded media.
    for v in cross(&sp, 3) {
        r.subpixel(&v);
    }
    let fills = [0.0, -0.0, 1.0, below(1.0), above(1.0), f64::from_bits(1), -f64::from_bits(1), 0.5, 0.3, 4.0 / 3.0,
                 -0.1, 2.0, -3.0, f64::NAN, f64::INFINITY];
    for fr in fills {
        for (e1, e2) in [(1.0, 4.0), (2.0, 9.0), (2.0, 8.0), (1e-300, 1e300), (f64::MAX, f64::MIN_POSITIVE)] {
            r.subpixel(&[e1, e2, fr]);
        }
    }
    for v in [[0.0, 4.0, 0.5], [1.0, 0.0, 0.5], [f64::NAN, 4.0, 0.5], [1.0, -4.0, 0.5]] {
        r.subpixel(&v);
    }
    for _ in 0..1000 {
        r.subpixel(&[rng.log_uniform(-1.0, 2.0), rng.log_uniform(-1.0, 2.0), rng.unit()]);
    }
    for _ in 0..300 {
        r.subpixel(&[rng.log_uniform(-300.0, 300.0), rng.log_uniform(-300.0, 300.0), rng.unit()]);
    }
}

pub fn emit() -> Module {
    let mut r = Rows { m: Module::new("fdtd"), keys: BTreeSet::new() };
    emit_1a(&mut r);
    emit_1b(&mut r);
    emit_1c(&mut r);
    r.m
}
