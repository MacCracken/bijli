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

use std::collections::BTreeSet;

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::fdtd::Fdtd1d;

/// One operation of a script.
#[derive(Clone, Copy)]
enum Op {
    Eps(usize, usize, f64),
    Mu(usize, usize, f64),
    Source(usize, f64),
    Step(usize),
    Run(usize, Option<usize>, f64, f64),
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
        let key = format!("{func} {:?}", sc.ins);
        assert!(self.keys.insert(key), "a repeated {func} script");
        self.m.row(func, Class::T, &sc.ins, &sc.outs());
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

pub fn emit() -> Module {
    let mut r = Rows { m: Module::new("fdtd"), keys: BTreeSet::new() };
    emit_1a(&mut r);
    r.m
}
