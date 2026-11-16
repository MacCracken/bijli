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
//!
//! M5 part 1d appends the utilities. `NonUniformGrid` rows take a grid as three sequences (`v:<len>` and its widths);
//! `uniform` and `try_uniform` over sizes and the specials, `graded` (class C: powf) over cell counts, every special
//! pair, wide ratios and seeded grades, and each query (`validate`, `min_spacing`, `courant_dt`, `try_courant_dt`,
//! `extents`) over uniform, empty-axis, special-holding, graded, extreme and seeded grids. `fourth_order_derivative`
//! takes a sequence and dx (x, x², x³, specials inserted, seeded fields). A `CwAccumulator` script is
//! `i:<n> f:<frequency> f:<dt>` (built by `new` for `CwAccumulator::new` rows, by `try_new` otherwise), then
//! `i:1 i:<step> v:<len> f:…` accumulate, `i:2 …` try_accumulate, `i:3 i:<k> i:<s0> f:<a> f:<b>` k accumulations of the
//! exact field a·((step + c) mod 11) − b; its state the samples, ω and dt (from `Debug`), n, the digests of real, imag,
//! the amplitudes and the phases, and each at 0, n/2, n − 1. `find_eigenmodes` (class C) takes a signal given or
//! generated by exact recurrences x_{k+1} = a·x_k + b·x_{k−1} (`i:<n> v:<g>` and g of `f:<x0> f:<x1> f:<a> f:<b>`),
//! then dt, the range and n_freq; its outputs the frequencies, the Qs and the amplitudes, each a sequence.
//!
//! M5 part 1e appends `nf2ff_2d_phasor` (class C: sin_cos and hypot). Its inputs are three accumulators, each `i:<n>
//! f:<frequency> f:<dt>`, an exact generator `i:<p> i:<q> i:<m> f:<a> f:<b> f:<c> f:<d>` (re[i] = a·((i·p + q) mod m)
//! − b, im[i] = c·((i·p + q + 7) mod m) − d) and `v:<k>` overrides `i:<index> f:<re> f:<im>`, then nx, ny, dx, the
//! margin and the angle count; its outputs the pattern and the angles, each a sequence, or the error. The grid: Rust's
//! test's validation; the accumulators' 10⁻¹² agreement straddled ulp by ulp; every special frequency, step and dx; the
//! margin at 0, nx/2 and ny/2; nx·ny past 2^64; short accumulators; grid shapes from 4 × 4 to the bench's 80 × 80 with
//! 1 to 180 angles; a zero field; specials and huge values on each side of the contour; wide frequencies, steps and
//! spacings; and seeded transforms.

use std::collections::BTreeSet;
use std::hint::black_box;

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, Tok, err, f, i, s};
use bijli::fdtd::{
    CwAccumulator, DebyeMaterial, DrudeMaterial, Fdtd1d, LorentzMaterial, LossyMaterial, NonUniformGrid, TfsfSource,
    differentiated_gaussian, find_eigenmodes, fourth_order_derivative, gaussian_pulse, modulated_gaussian_pulse,
    nf2ff_2d_phasor, ricker_wavelet, subpixel_permittivity_parallel, subpixel_permittivity_perpendicular,
    subpixel_permittivity_perpendicular_checked,
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
        let bx: Vec<f64> = ins.iter().map(|&x| black_box(x)).collect();
        let v = match (name, bx.as_slice()) {
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

// ── part 1d: the non-uniform grid, the fourth-order derivative, the CW accumulator, the eigenmodes ────────────────

/// A sequence: `v:<len>`, then its values.
fn seq(xs: &[f64]) -> Vec<Tok> {
    let mut out = vec![Tok::V(xs.len())];
    out.extend(xs.iter().map(|&x| f(x)));
    out
}

/// A non-uniform grid's tokens: its three axes as sequences.
fn grid_toks(g: &NonUniformGrid) -> Vec<Tok> {
    let mut out = seq(&g.dx);
    out.extend(seq(&g.dy));
    out.extend(seq(&g.dz));
    out
}

/// A `Result` as its value's tokens, or its error.
fn res_toks<T>(r: bijli::error::Result<T>, ok: impl Fn(T) -> Vec<Tok>) -> Vec<Tok> {
    match r {
        Ok(v) => ok(v),
        Err(e) => vec![err(&e), s(e.to_string())],
    }
}

/// One operation of a CW accumulator script.
#[derive(Clone)]
enum CwOp {
    Accumulate(usize, Vec<f64>),
    TryAccumulate(usize, Vec<f64>),
    /// k accumulations from step s₀ of the field g[c] = a·((step + c) mod 11) − b, the accumulator's length.
    Generated(usize, usize, f64, f64),
}

/// The generated field of `CwOp::Generated` at `step` over n points: exact operations only.
fn cw_field(n: usize, step: usize, a: f64, b: f64) -> Vec<f64> {
    (0..n).map(|c| a * ((step + c) % 11) as f64 - b).collect()
}

/// A CW accumulator script: its constructor (`try_new` when checked), then operations.
struct CwScript {
    acc: Result<CwAccumulator, Vec<Tok>>,
    ins: Vec<Tok>,
    errs: Vec<Tok>,
}

impl CwScript {
    fn new(checked: bool, n: usize, freq: f64, dt: f64) -> Self {
        let acc = if checked {
            CwAccumulator::try_new(n, black_box(freq), black_box(dt)).map_err(|e| vec![err(&e), s(e.to_string())])
        } else {
            Ok(CwAccumulator::new(n, black_box(freq), black_box(dt)))
        };
        Self { acc, ins: vec![u(n), f(freq), f(dt)], errs: Vec::new() }
    }

    fn apply(&mut self, op: CwOp) {
        let acc = self.acc.as_mut().expect("operations only on an accumulator that was made");
        match op {
            CwOp::Accumulate(step, field) => {
                self.ins.extend([i(1), u(step)]);
                self.ins.extend(seq(&field));
                acc.accumulate(black_box(&field), step);
            }
            CwOp::TryAccumulate(step, field) => {
                self.ins.extend([i(2), u(step)]);
                self.ins.extend(seq(&field));
                if let Err(e) = acc.try_accumulate(black_box(&field), step) {
                    self.errs.extend([err(&e), s(e.to_string())]);
                }
            }
            CwOp::Generated(k, s0, a, b) => {
                self.ins.extend([i(3), u(k), u(s0), f(a), f(b)]);
                let n = acc.real.len();
                for j in 0..k {
                    acc.accumulate(black_box(&cw_field(n, s0 + j, a, b)), s0 + j);
                }
            }
        }
    }

    fn outs(&self) -> Vec<Tok> {
        let acc = match &self.acc {
            Err(e) => return e.clone(),
            Ok(a) => a,
        };
        let d = format!("{acc:?}");
        let (amp, ph) = (acc.amplitude(), acc.phase());
        let n = acc.real.len();
        let mut out = self.errs.clone();
        out.extend([
            u(acc.n_samples),
            f(dbg_f64(&d, "omega")),
            f(dbg_f64(&d, "dt")),
            u(n),
            i(digest(&acc.real)),
            i(digest(&acc.imag)),
            i(digest(&amp)),
            i(digest(&ph)),
        ]);
        if n > 0 {
            for p in [0, n / 2, n - 1] {
                out.extend([f(acc.real[p]), f(acc.imag[p]), f(amp[p]), f(ph[p])]);
            }
        }
        out
    }
}

/// An eigenmode signal: given (`v:<len>` and its samples) or generated, `i:<n>`, `v:<g>`, then g second-order
/// recurrences (x₀, x₁, a, b), x_{k+1} = a·x_k + b·x_{k−1} (a damped sinusoid for a² < −4b), summed sample by sample
/// in order: exact operations only.
enum Signal {
    Given(Vec<f64>),
    Generated(usize, Vec<[f64; 4]>),
}

impl Signal {
    fn samples(&self) -> Vec<f64> {
        match self {
            Self::Given(v) => v.clone(),
            Self::Generated(n, gens) => {
                let mut out = vec![0.0; *n];
                for (gi, g) in gens.iter().enumerate() {
                    let (mut prev, mut cur) = (g[0], g[1]);
                    for (k, o) in out.iter_mut().enumerate() {
                        let x = match k {
                            0 => g[0],
                            1 => g[1],
                            _ => {
                                let next = g[2] * cur + g[3] * prev;
                                prev = cur;
                                cur = next;
                                next
                            }
                        };
                        *o = if gi == 0 { x } else { *o + x };
                    }
                }
                out
            }
        }
    }

    fn toks(&self) -> Vec<Tok> {
        match self {
            Self::Given(v) => seq(v),
            Self::Generated(n, gens) => {
                let mut out = vec![u(*n), Tok::V(gens.len())];
                for g in gens {
                    out.extend(g.iter().map(|&x| f(x)));
                }
                out
            }
        }
    }
}

impl Rows {
    fn grid_rows(&mut self, g: &NonUniformGrid) {
        let ins = grid_toks(g);
        if self.keys.contains(&format!("NonUniformGrid::validate {ins:?}")) {
            return;
        }
        let gb = black_box(g);
        self.emit_toks("NonUniformGrid::validate", Class::A, &ins, &res_toks(gb.validate(), |()| vec![i(0)]));
        self.emit_toks("NonUniformGrid::min_spacing", Class::A, &ins, &[f(gb.min_spacing())]);
        self.emit_toks("NonUniformGrid::courant_dt", Class::A, &ins, &[f(gb.courant_dt())]);
        let tc = res_toks(gb.try_courant_dt(), |x| vec![f(x)]);
        self.emit_toks("NonUniformGrid::try_courant_dt", Class::A, &ins, &tc);
        let e = gb.extents();
        self.emit_toks("NonUniformGrid::extents", Class::A, &ins, &floats(&e));
    }

    fn fourth(&mut self, field: &[f64], dx: f64) {
        let mut ins = seq(field);
        ins.push(f(dx));
        if self.keys.contains(&format!("fdtd::fourth_order_derivative {ins:?}")) {
            return;
        }
        let d = fourth_order_derivative(black_box(field), black_box(dx));
        self.emit_toks("fdtd::fourth_order_derivative", Class::A, &ins, &seq(&d));
    }

    fn graded(&mut self, n: usize, d_min: f64, d_max: f64) {
        let ins = [u(n), f(d_min), f(d_max)];
        if self.keys.contains(&format!("NonUniformGrid::graded {ins:?}")) {
            return;
        }
        let r = NonUniformGrid::graded(black_box(n), black_box(d_min), black_box(d_max));
        self.emit_toks("NonUniformGrid::graded", Class::C, &ins, &res_toks(r, |v| seq(&v)));
    }

    fn cw(&mut self, func: &str, sc: &CwScript) {
        self.emit_toks(func, Class::T, &sc.ins, &sc.outs());
    }

    fn eigen(&mut self, sig: &Signal, dt: f64, fmin: f64, fmax: f64, n_freq: usize) {
        let mut ins = sig.toks();
        ins.extend([f(dt), f(fmin), f(fmax), u(n_freq)]);
        if self.keys.contains(&format!("fdtd::find_eigenmodes {ins:?}")) {
            return;
        }
        let samples = sig.samples();
        let (bdt, bmin, bmax) = (black_box(dt), black_box(fmin), black_box(fmax));
        let r = find_eigenmodes(black_box(&samples), bdt, bmin, bmax, black_box(n_freq));
        let outs = res_toks(r, |m| {
            let mut o = seq(&m.frequencies);
            o.extend(seq(&m.quality_factors));
            o.extend(seq(&m.amplitudes));
            o
        });
        self.emit_toks("fdtd::find_eigenmodes", Class::C, &ins, &outs);
    }
}

fn emit_1d_grids(r: &mut Rows) {
    let mut sp = specials();
    sp.extend([1e-3, 2.5e-4]);
    // uniform and try_uniform over sizes and every special spacing.
    for nx in [0_usize, 1, 2, 5] {
        for ny in [0_usize, 1, 2, 5] {
            for nz in [0_usize, 1, 2, 5] {
                for &s in &sp {
                    let ins = [u(nx), u(ny), u(nz), f(s)];
                    let g = NonUniformGrid::uniform(black_box(nx), ny, nz, black_box(s));
                    r.emit_toks("NonUniformGrid::uniform", Class::A, &ins, &grid_toks(&g));
                    let t = NonUniformGrid::try_uniform(black_box(nx), ny, nz, black_box(s));
                    r.emit_toks("NonUniformGrid::try_uniform", Class::A, &ins, &res_toks(t, |g| grid_toks(&g)));
                }
            }
        }
    }
    // graded: cell counts against every special pair, the Rust tests' cases, ratios far from 1, seeded grades.
    for n in [0_usize, 1, 2, 3, 10, 64] {
        for &a in &sp {
            for &b in &sp {
                r.graded(n, a, b);
            }
        }
    }
    let cases = [(10, 1e-4, 1e-3), (1, 1e-4, 1e-3), (10, 0.0, 1e-3), (10, f64::NAN, 1e-3), (10, 1e-4, f64::INFINITY),
                 (7, 1e-3, 1e-5), (1000, 1e-6, 1e-2), (100, 1e-300, 1e300), (100, 1e300, 1e-300),
                 (3, f64::MIN_POSITIVE, f64::MAX), (40, 1.0, 1.0), (5, 5e-324, 1.0)];
    for (n, a, b) in cases {
        r.graded(n, a, b);
    }
    let mut rng = Rng::new(0x9ad3_d0c0);
    for _ in 0..300 {
        let n = usize::try_from(rng.int(2, 100)).unwrap_or(2);
        r.graded(n, rng.log_uniform(-7.0, -2.0), rng.log_uniform(-6.0, -1.0));
    }

    // The grid queries over uniform grids, an empty axis each, a special at each axis's ends and middle, the Rust
    // tests' grids, graded axes, extremes and seeded grids.
    for n in [0_usize, 1, 3] {
        for &s in &sp {
            r.grid_rows(&NonUniformGrid::uniform(n, n, n, s));
        }
    }
    for (nx, ny, nz) in [(0, 2, 2), (2, 0, 2), (2, 2, 0), (0, 0, 2), (0, 0, 0)] {
        r.grid_rows(&NonUniformGrid::uniform(nx, ny, nz, 1e-3));
    }
    for &v in &sp {
        for axis in 0..3 {
            for pos in [0_usize, 1, 2] {
                let mut g = NonUniformGrid::uniform(3, 4, 5, 1e-3);
                let a = match axis {
                    0 => &mut g.dx,
                    1 => &mut g.dy,
                    _ => &mut g.dz,
                };
                let k = pos * (a.len() - 1) / 2;
                a[k] = v;
                r.grid_rows(&g);
            }
        }
    }
    let mut g = NonUniformGrid::uniform(5, 5, 5, 1e-3);
    g.dz[3] = 2.5e-4;
    r.grid_rows(&g);
    let mut g = NonUniformGrid::uniform(4, 4, 4, 1e-3);
    g.dy[2] = f64::NAN;
    r.grid_rows(&g);
    let graded = |n, a, b| NonUniformGrid::graded(n, a, b).expect("a valid grading");
    r.grid_rows(&NonUniformGrid { dx: graded(10, 1e-4, 1e-3), dy: graded(7, 1e-3, 1e-5), dz: vec![2e-4; 3] });
    r.grid_rows(&NonUniformGrid { dx: vec![1e300; 4], dy: vec![f64::MAX; 2], dz: vec![1.0] });
    r.grid_rows(&NonUniformGrid { dx: vec![5e-324; 3], dy: vec![f64::MIN_POSITIVE], dz: vec![1e-300, 2e-300] });
    r.grid_rows(&NonUniformGrid { dx: vec![-0.0, 1.0], dy: vec![1.0], dz: vec![1.0] });
    for _ in 0..100 {
        let mut axis = || {
            let n = usize::try_from(rng.int(1, 12)).unwrap_or(1);
            (0..n).map(|_| rng.log_uniform(-6.0, -2.0)).collect::<Vec<f64>>()
        };
        let mut g = NonUniformGrid { dx: axis(), dy: axis(), dz: axis() };
        if rng.int(0, 4) == 0 {
            let bad = [0.0, -1e-3, f64::NAN, f64::INFINITY][usize::try_from(rng.int(0, 3)).unwrap_or(0)];
            let k = usize::try_from(rng.int(0, i64::try_from(g.dy.len()).unwrap_or(1) - 1)).unwrap_or(0);
            g.dy[k] = bad;
        }
        r.grid_rows(&g);
    }

    // fourth_order_derivative: lengths 0 to 8, 16 and 33 over x, x², x³, specials inserted and seeded fields, against
    // the special and typical spacings.
    let mut dxs = specials();
    dxs.extend([0.1, 0.05, 3.0, 1e-300, 1e300]);
    for n in [0_usize, 1, 2, 3, 4, 5, 6, 7, 8, 16, 33] {
        let lin: Vec<f64> = (0..n).map(|k| k as f64 * 0.1).collect();
        let quad: Vec<f64> = lin.iter().map(|x| x * x).collect();
        let cube: Vec<f64> = (0..n).map(|k| k as f64 * 0.05).map(|x| x * x * x).collect();
        for &dx in &dxs {
            r.fourth(&lin, dx);
            r.fourth(&quad, dx);
            r.fourth(&cube, dx);
        }
        for &v in &specials() {
            if n > 0 {
                for pos in [0, n / 2, n - 1] {
                    let mut fld = quad.clone();
                    fld[pos] = v;
                    r.fourth(&fld, 0.1);
                }
            }
        }
    }
    for _ in 0..400 {
        let n = usize::try_from(rng.int(0, 40)).unwrap_or(0);
        let fld: Vec<f64> = (0..n).map(|_| rng.signed_log_uniform(-3.0, 3.0)).collect();
        r.fourth(&fld, rng.log_uniform(-6.0, 0.0));
    }
    let wide: Vec<f64> = (0..12).map(|_| rng.signed_log_uniform(-300.0, 300.0)).collect();
    r.fourth(&wide, 1e-300);
    r.fourth(&wide, 1e300);
}

fn emit_1d_cw(r: &mut Rows) {
    // new and try_new over sizes and every special frequency and step.
    let sp = specials();
    for n in [0_usize, 1, 3] {
        for &fr in &sp {
            for &dt in &sp {
                r.cw("CwAccumulator::new", &CwScript::new(false, n, fr, dt));
                r.cw("CwAccumulator::try_new", &CwScript::new(true, n, fr, dt));
            }
        }
    }
    for (n, fr, dt) in [(4, 1e9, 1e-12), (4, f64::NAN, 1e-12), (4, -1e9, 1e-12), (4, 1e9, 0.0), (4, 1e9, f64::INFINITY),
                        (4, -0.0, 1e-12), (2, 1e9, 1e-11)] {
        r.cw("CwAccumulator::try_new", &CwScript::new(true, n, fr, dt));
    }
    // accumulate: the Rust test's mismatches, then given and generated fields at the sizes the packed kernel splits
    // (a multiple of 4 and its tail), the time step at its ends, special samples, frequencies and steps.
    let ra = "CwAccumulator::accumulate";
    let mut sc = CwScript::new(true, 3, 1e9, 1e-12);
    for op in [CwOp::Accumulate(0, vec![1.0; 5]), CwOp::Accumulate(0, vec![1.0]),
               CwOp::TryAccumulate(1, vec![1.0; 4]), CwOp::TryAccumulate(1, vec![1.0; 3])] {
        sc.apply(op);
        r.cw("CwAccumulator::try_accumulate", &sc);
    }
    for n in [1_usize, 2, 3, 4, 5, 7, 8, 9, 100] {
        let mut sc = CwScript::new(true, n, 1e9, 1e-12);
        let mut at = 0;
        for fr in [1, 10, 100, 1000] {
            sc.apply(CwOp::Generated(fr - at, at, 0.25, 1.375));
            at = fr;
            r.cw(ra, &sc);
        }
    }
    for step in [0_usize, 1, 7, 1 << 20, 1 << 40, usize::MAX - 1, usize::MAX] {
        let mut sc = CwScript::new(true, 4, 1e9, 1e-12);
        sc.apply(CwOp::Accumulate(step, vec![1.0, -2.0, 0.5, 3.0]));
        r.cw(ra, &sc);
    }
    for &v in &sp {
        let mut sc = CwScript::new(true, 5, 3e9, 1e-11);
        sc.apply(CwOp::Accumulate(3, vec![1.0, v, -0.5, v, 2.0]));
        r.cw(ra, &sc);
    }
    // An unchecked accumulator (`CwAccumulator::new` rows are built by `new`, the others by `try_new`) sampled at
    // every special frequency and step.
    for &fr in &sp {
        for &dt in &sp {
            let mut sc = CwScript::new(false, 2, fr, dt);
            sc.apply(CwOp::Accumulate(5, vec![1.0, -1.0]));
            r.cw("CwAccumulator::new", &sc);
        }
    }
    for fr in [1e-300, 1e15, 1e300] {
        let mut sc = CwScript::new(true, 6, fr, 1e-12);
        sc.apply(CwOp::Generated(50, 0, 1.0, 5.0));
        r.cw(ra, &sc);
    }
    for (len, step) in [(0, 2), (2, 2), (6, 2), (1, 0)] {
        let mut sc = CwScript::new(true, 4, 1e9, 1e-12);
        sc.apply(CwOp::Generated(3, 0, 1.0, 0.5));
        sc.apply(CwOp::TryAccumulate(step, vec![0.75; len]));
        r.cw("CwAccumulator::try_accumulate", &sc);
    }
    let mut rng = Rng::new(0x0c3a_cc00);
    for _ in 0..10 {
        let n = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        let (fr, dt) = (rng.log_uniform(6.0, 12.0), rng.log_uniform(-15.0, -9.0));
        let (a, b) = (rng.signed_log_uniform(-3.0, 3.0), rng.signed_log_uniform(-3.0, 3.0));
        let mut sc = CwScript::new(true, n, fr, dt);
        let mut at = 0;
        for fr2 in [1, 10, 100, 1000] {
            sc.apply(CwOp::Generated(fr2 - at, at, a, b));
            at = fr2;
            r.cw(ra, &sc);
        }
    }
}

fn emit_1d_eigen(r: &mut Rows) {
    // The validation: the Rust tests' cases, every special for dt and the range, short and non-finite signals.
    let s3 = Signal::Given(vec![1.0, 0.5, 0.25]);
    let cases = [(1e-12, 1e9, 10e9, 1), (0.0, 1e9, 10e9, 10), (f64::NAN, 1e9, 10e9, 10), (1e-12, f64::NAN, 10e9, 10),
                 (1e-12, 1e9, f64::NAN, 10), (1e-12, 1e9, f64::INFINITY, 10), (1e-12, -1e9, 10e9, 10),
                 (1e-12, 10e9, 1e9, 100), (1e-12, 1e9, 10e9, 0), (1e-12, 1e9, 1e9, 10)];
    for (dt, a, b, n) in cases {
        r.eigen(&s3, dt, a, b, n);
    }
    r.eigen(&Signal::Given(vec![]), 1e-12, 1e9, 10e9, 100);
    r.eigen(&Signal::Given(vec![1.0, f64::NAN]), 1e-12, 1e9, 10e9, 10);
    r.eigen(&Signal::Given(vec![1.0]), 1e-12, 10e9, 1e9, 100);
    for &v in &specials() {
        r.eigen(&Signal::Given(vec![0.5, v, -0.25, 1.0]), 1.0, 0.0, 0.4, 9);
        r.eigen(&s3, v, 0.0, 0.4, 9);
        for &w in &specials() {
            r.eigen(&s3, 1.0, v, w, 9);
        }
    }
    // Damped modes from exact recurrences, at the bench's size and smaller: one mode, two, merged, at the band's
    // edges, past Nyquist; flat and constant signals.
    let one = |a: f64, b: f64| [0.0, 1.0, a, b];
    let s30 = Signal::Generated(4096, vec![one(1.9626, -0.998)]);
    r.eigen(&s30, 1e-12, 1e10, 5e10, 256);
    r.eigen(&Signal::Generated(4096, vec![one(1.9626, -0.998), one(1.7, -0.995)]), 1e-12, 1e10, 8e10, 256);
    for (n, nf) in [(256, 16), (1000, 64), (2048, 128), (1500, 2), (1500, 3)] {
        r.eigen(&Signal::Generated(n, vec![one(1.95, -0.99)]), 1e-12, 1e9, 1e11, nf);
        r.eigen(&Signal::Generated(n, vec![one(1.95, -0.99), one(1.6, -0.993)]), 1e-12, 1e9, 1.6e11, nf);
    }
    r.eigen(&Signal::Generated(1024, vec![one(1.95, -0.99), one(1.949, -0.99)]), 1e-12, 1e9, 1e11, 64);
    r.eigen(&Signal::Generated(1024, vec![one(1.95, -0.99)]), 1e-12, 0.0, 4e9, 40);
    r.eigen(&Signal::Generated(1024, vec![one(1.95, -0.99)]), 1e-12, 4e9, 6e11, 64);
    r.eigen(&Signal::Generated(800, vec![one(-1.9, -0.98)]), 1e-12, 1e9, 5e11, 50);
    r.eigen(&Signal::Generated(2100, vec![one(1.999, -0.9999)]), 1e-12, 0.0, 1e10, 101);
    r.eigen(&Signal::Given(vec![0.0; 8]), 1.0, 0.0, 0.25, 5);
    r.eigen(&Signal::Given(vec![1.0; 64]), 1.0, 0.0, 0.5, 33);
    r.eigen(&Signal::Given((0..64).map(|k| if k % 2 == 0 { 1.0 } else { -1.0 }).collect()), 1.0, 0.0, 0.5, 33);
    let mut rng = Rng::new(0x3e16_e000);
    for _ in 0..20 {
        let n = usize::try_from(rng.int(64, 1500)).unwrap_or(64);
        let gens = (0..rng.int(1, 3)).map(|_| one(rng.uniform(-1.98, 1.98), -rng.uniform(0.95, 0.999))).collect();
        let fmax = rng.uniform(0.05, 0.6);
        let nf = usize::try_from(rng.int(2, 120)).unwrap_or(2);
        r.eigen(&Signal::Generated(n, gens), 1.0, 0.0, fmax, nf);
    }
}

// ── part 1e: NF2FF's phasor form ───────────────────────────────────────────────────────────────────────────────────

/// A phasor accumulator for the NF2FF rows: `n` points at `freq` and `dt` (built by `CwAccumulator::new`), its parts
/// generated exactly, re[i] = a·((i·p + q) mod m) − b and im[i] = c·((i·p + q + 7) mod m) − d, then `overrides`
/// (index, re, im) written over them.
#[derive(Clone)]
struct Phasors {
    n: usize,
    freq: f64,
    dt: f64,
    spec: (usize, usize, usize, f64, f64, f64, f64),
    overrides: Vec<(usize, f64, f64)>,
}

impl Phasors {
    fn new(n: usize, freq: f64, dt: f64, spec: (usize, usize, usize, f64, f64, f64, f64)) -> Self {
        Self { n, freq, dt, spec, overrides: Vec::new() }
    }

    fn build(&self) -> CwAccumulator {
        let mut acc = CwAccumulator::new(self.n, self.freq, self.dt);
        let (p, q, m, a, b, c, d) = self.spec;
        for i in 0..self.n {
            acc.real[i] = a * ((i * p + q) % m) as f64 - b;
            acc.imag[i] = c * ((i * p + q + 7) % m) as f64 - d;
        }
        for &(k, re, im) in &self.overrides {
            acc.real[k] = re;
            acc.imag[k] = im;
        }
        acc
    }

    fn toks(&self) -> Vec<Tok> {
        let (p, q, m, a, b, c, d) = self.spec;
        let mut out = vec![u(self.n), f(self.freq), f(self.dt), u(p), u(q), u(m), f(a), f(b), f(c), f(d)];
        out.push(Tok::V(self.overrides.len()));
        for &(k, re, im) in &self.overrides {
            out.extend([u(k), f(re), f(im)]);
        }
        out
    }
}

impl Rows {
    /// An `nf2ff_2d_phasor` row (class C: sin_cos and hypot): the three accumulators, then nx, ny, dx, the margin and
    /// the angle count; its outputs the pattern and the angles, or the error.
    #[allow(clippy::too_many_arguments)]
    fn nf2ff(&mut self, ez: &Phasors, hx: &Phasors, hy: &Phasors, nx: usize, ny: usize, dx: f64, margin: usize,
             n_angles: usize) {
        let mut ins = ez.toks();
        ins.extend(hx.toks());
        ins.extend(hy.toks());
        ins.extend([u(nx), u(ny), f(dx), u(margin), u(n_angles)]);
        if self.keys.contains(&format!("fdtd::nf2ff_2d_phasor {ins:?}")) {
            return;
        }
        let (a, b, c) = (ez.build(), hx.build(), hy.build());
        let r = nf2ff_2d_phasor(black_box(&a), black_box(&b), black_box(&c), nx, ny, black_box(dx), margin, n_angles);
        let outs = res_toks(r, |res| {
            let mut o = seq(&res.pattern);
            o.extend(seq(&res.angles));
            o
        });
        self.emit_toks("fdtd::nf2ff_2d_phasor", Class::C, &ins, &outs);
    }
}

fn emit_1e(r: &mut Rows) {
    let g1 = (3, 1, 13, 0.25, 1.5, 0.125, 0.75);
    let g2 = (5, 2, 17, -0.5, 2.0, 0.375, -1.25);
    let g3 = (7, 4, 11, 1.0, 5.0, -0.25, 1.0);
    let ph = |n, fr, dt, g| Phasors::new(n, fr, dt, g);
    // The validation: Rust's test's cases, the accumulators' agreement straddled at 10⁻¹², every special frequency
    // and step, the margin at 0, nx/2 and ny/2, every special dx, no angle, short accumulators, nx·ny past 2^64.
    let (fr, dt) = (1e9, 1e-12);
    let a = ph(400, fr, dt, g1);
    r.nf2ff(&a, &a, &a, 20, 20, 1e-3, 3, 8);
    r.nf2ff(&a, &ph(400, 2e9, dt, g1), &a, 20, 20, 1e-3, 3, 8);
    r.nf2ff(&a, &a, &ph(10, fr, dt, g1), 20, 20, 1e-3, 3, 8);
    let z = ph(400, 0.0, dt, g1);
    r.nf2ff(&z, &z, &z, 20, 20, 1e-3, 3, 8);
    r.nf2ff(&a, &a, &a, 20, 20, f64::NAN, 3, 8);
    r.nf2ff(&a, &a, &a, 20, 20, 1e-3, 10, 8);
    r.nf2ff(&a, &a, &a, 20, 20, 1e-3, 3, 0);
    for k in [-3_i64, -2, -1, 0, 1, 2, 3] {
        let edge = fr * (1.0 + 1e-12);
        let fk = f64::from_bits(edge.to_bits().wrapping_add_signed(k));
        r.nf2ff(&a, &ph(400, fk, dt, g1), &a, 20, 20, 1e-3, 3, 8);
        r.nf2ff(&a, &a, &ph(400, fr, f64::from_bits((dt * (1.0 + 1e-12)).to_bits().wrapping_add_signed(k)), g1),
                20, 20, 1e-3, 3, 8);
    }
    for &v in &specials() {
        let s = ph(400, v, dt, g1);
        r.nf2ff(&s, &s, &s, 20, 20, 1e-3, 3, 8);
        let s = ph(400, fr, v, g1);
        r.nf2ff(&s, &s, &s, 20, 20, 1e-3, 3, 8);
        r.nf2ff(&a, &a, &a, 20, 20, v, 3, 8);
    }
    for (nx, ny, m) in [(20, 20, 0), (20, 20, 9), (20, 20, 10), (20, 21, 10), (21, 20, 10), (4, 4, 1), (4, 4, 2),
                        (3, 20, 1), (usize::MAX, 20, 3), (20, usize::MAX, 3), (1 << 33, 1 << 33, 3)] {
        r.nf2ff(&a, &a, &a, nx, ny, 1e-3, m, 8);
    }
    for (ne, nh, ny2) in [(399, 400, 400), (400, 399, 400), (400, 400, 399), (1000, 1000, 1000)] {
        r.nf2ff(&ph(ne, fr, dt, g1), &ph(nh, fr, dt, g1), &ph(ny2, fr, dt, g1), 20, 20, 1e-3, 3, 8);
    }
    // Transforms: grid shapes from the smallest contour to the bench's 80 × 80, angle counts from 1 to 180, a zero
    // field, specials and huge values on the contour, wide frequencies and steps.
    for (nx, ny, m) in [(4, 4, 1), (5, 7, 1), (6, 6, 2), (20, 20, 3), (30, 24, 4), (48, 48, 5), (80, 80, 5)] {
        for na in [1_usize, 2, 4, 36, 180] {
            let n = nx * ny;
            r.nf2ff(&ph(n, 3e9, 1e-12, g1), &ph(n, 3e9, 1e-12, g2), &ph(n, 3e9, 1e-12, g3), nx, ny, 1e-3, m, na);
        }
    }
    let n = 24 * 24;
    let zero = ph(n, 1e9, 1e-12, (1, 0, 1, 0.0, 0.0, 0.0, 0.0));
    r.nf2ff(&zero, &zero, &zero, 24, 24, 1e-3, 3, 36);
    let idx = |x: usize, y: usize| y * 24 + x;
    for &v in &[f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e300, f64::MAX, 5e-324, -0.0] {
        for (which, at) in [(0, idx(3, 3)), (0, idx(12, 3)), (1, idx(12, 2)), (2, idx(2, 12)), (2, idx(20, 12))] {
            let mut fs = [ph(n, 1e9, 1e-12, g1), ph(n, 1e9, 1e-12, g2), ph(n, 1e9, 1e-12, g3)];
            fs[which].overrides.push((at, v, 0.5));
            r.nf2ff(&fs[0], &fs[1], &fs[2], 24, 24, 1e-3, 3, 12);
        }
    }
    for (fr2, dt2, dx2) in [(1e6, 1e-9, 1.0), (1e12, 1e-16, 1e-6), (1e300, 1e-300, 1e-3), (1e-300, 1e300, 1e300),
                            (3e9, 1e-12, 1e300), (3e9, 1e-12, 1e-300), (1.5e9, 2.5e-11, 0.0125)] {
        let n = 20 * 20;
        r.nf2ff(&ph(n, fr2, dt2, g1), &ph(n, fr2, dt2, g2), &ph(n, fr2, dt2, g3), 20, 20, dx2, 3, 36);
    }
    let mut rng = Rng::new(0x2f2f_ef00);
    for _ in 0..20 {
        let nx = usize::try_from(rng.int(4, 60)).unwrap_or(4);
        let ny = usize::try_from(rng.int(4, 60)).unwrap_or(4);
        let m = usize::try_from(rng.int(1, i64::try_from(nx.min(ny) / 2).unwrap_or(2) - 1)).unwrap_or(1);
        let (fr3, dt3, dx3) = (rng.log_uniform(7.0, 11.0), rng.log_uniform(-14.0, -10.0), rng.log_uniform(-5.0, -2.0));
        let mut spec = || {
            let p = usize::try_from(rng.int(1, 40)).unwrap_or(1);
            let q = usize::try_from(rng.int(0, 40)).unwrap_or(0);
            let mm = usize::try_from(rng.int(2, 50)).unwrap_or(2);
            (p, q, mm, rng.signed_log_uniform(-3.0, 3.0), rng.signed_log_uniform(-3.0, 3.0),
             rng.signed_log_uniform(-3.0, 3.0), rng.signed_log_uniform(-3.0, 3.0))
        };
        let n = nx * ny;
        let (ga, gb, gc) = (spec(), spec(), spec());
        let na = usize::try_from(rng.int(1, 200)).unwrap_or(1);
        r.nf2ff(&ph(n, fr3, dt3, ga), &ph(n, fr3, dt3, gb), &ph(n, fr3, dt3, gc), nx, ny, dx3, m, na);
    }
}

// ── part 2a-1a: the media's coefficients at their cancellation (ADR 0002 row 40) ─────────────────────────────────

/// `n` values from `x` walking up one ulp at a time, starting `n/2` ulps below it.
fn ulp_walk(x: f64, n: usize) -> Vec<f64> {
    let mut y = x;
    for _ in 0..n / 2 {
        y = y.next_down();
    }
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push(y);
        y = y.next_up();
    }
    out
}

/// The lossy, Drude, Lorentz and Debye coefficients where Rust's (1 ∓ s), (2 ∓ γΔt), 2 − ω₀²Δt², 2τ/Δt ∓ 1 and
/// ε_s − ε_∞ cancel (the port rounds each once from its exact value, ADR 0002 row 40): relative offsets 10⁻¹ to 10⁻¹⁶
/// either side, ulp walks, exact zeros, the fast path's box edges at 2^±120 (src/crmath_edge.cyr), subnormal and huge
/// operands, and seeded media near the cancellation.
fn emit_2a1a(r: &mut Rows) {
    let e0 = bijli::field::EPSILON_0;
    let offs: Vec<f64> = (1..=16).flat_map(|k| [10f64.powi(-k), -10f64.powi(-k)]).collect();
    let p2 = |k: i32| 2f64.powi(k);
    // The lossy C_a: σΔt = 2ε₀ε_r(1 + ρ).
    for eps_r in [1.0, 2.25, 11.7] {
        for dt in [1e-12, 1e-15, p2(-40)] {
            let s0 = 2.0 * e0 * eps_r / dt;
            for &rho in &offs {
                r.lossy(&[eps_r, s0 * (1.0 + rho), dt, 1e-4]);
            }
            for s in ulp_walk(s0, 16) {
                r.lossy(&[eps_r, s, dt, 1e-4]);
            }
        }
    }
    r.lossy(&[1.0, 2.0 * e0 * p2(40), p2(-40), 1e-4]);
    r.lossy(&[4.0, 8.0 * e0 * p2(40), p2(-40), 1e-4]);
    // Drude's and Lorentz's γΔt = 2(1 + ρ), Lorentz's ω₀Δt = √2(1 + ρ).
    for dt in [1e-16, 1e-17, p2(-50)] {
        let g0 = 2.0 / dt;
        let w0 = std::f64::consts::SQRT_2 / dt;
        for &rho in &offs {
            r.drude(&[1.0, 1.37e16, g0 * (1.0 + rho), dt]);
            r.lorentz(&[2.0, 2.0, w0 * (1.0 + rho), 1e13, dt]);
            r.lorentz(&[2.0, 2.0, 1e15, g0 * (1.0 + rho), dt]);
        }
        for g in ulp_walk(g0, 16) {
            r.drude(&[1.0, 1.37e16, g, dt]);
            r.lorentz(&[2.0, 2.0, 1e15, g, dt]);
        }
        for w in ulp_walk(w0, 16) {
            r.lorentz(&[2.0, 2.0, w, 1e13, dt]);
        }
    }
    r.drude(&[1.0, 1e16, p2(51), p2(-50)]);
    r.lorentz(&[2.0, 2.0, 1e15, p2(51), p2(-50)]);
    // Debye's 2τ = Δt(1 + ρ) and ε_s = ε_∞(1 + ρ).
    for dt in [1e-13, 1e-14, p2(-44)] {
        for &rho in &offs {
            r.debye(&[4.9, 78.4, dt / 2.0 * (1.0 + rho), dt]);
            r.debye(&[4.9, 4.9 * (1.0 + rho), 1e-11, dt]);
        }
        for tau in ulp_walk(dt / 2.0, 16) {
            r.debye(&[4.9, 78.4, tau, dt]);
        }
        for es in ulp_walk(4.9, 16) {
            r.debye(&[4.9, es, 1e-11, dt]);
        }
    }
    r.debye(&[4.9, 78.4, p2(-45), p2(-44)]);
    // The fast path's box edges, each operand in turn, from a typical medium.
    let lossy = [4.0, 0.5, 1e-13, 1e-4];
    let drude = [1.0, 1.37e16, 4.05e13, 1e-16];
    let lorentz = [2.0, 2.0, 1e15, 1e13, 1e-16];
    let debye = [4.9, 78.4, 8.27e-12, 1e-13];
    let edges: Vec<f64> =
        [120, -120, 121, -121, 119, -119].iter().flat_map(|&k| [p2(k), -p2(k), 1.5 * p2(k)]).collect();
    let tiny = [f64::from_bits(1), 1e-310, f64::MIN_POSITIVE, 1e-300, 1e300, f64::MAX];
    for &x in edges.iter().chain(tiny.iter()) {
        for j in 0..4 {
            let mut v = lossy;
            v[j] = x;
            r.lossy(&v);
            let mut v = drude;
            v[j] = x;
            r.drude(&v);
            let mut v = debye;
            v[j] = x;
            r.debye(&v);
        }
        for j in 1..5 {
            let mut v = lorentz;
            v[j] = x;
            r.lorentz(&v);
        }
    }
    // Seeded media near their cancellation.
    let mut rng = Rng::new(0x2a1a_0040);
    for _ in 0..400 {
        let rho = if rng.unit() < 0.5 { 1.0 } else { -1.0 } * rng.log_uniform(-16.0, -1.0);
        let dt = rng.log_uniform(-17.0, -11.0);
        let eps_r = rng.uniform(1.0, 20.0);
        r.lossy(&[eps_r, 2.0 * e0 * eps_r / dt * (1.0 + rho), dt, rng.log_uniform(-6.0, -2.0)]);
        r.drude(&[rng.uniform(1.0, 10.0), rng.log_uniform(14.0, 17.0), 2.0 / dt * (1.0 + rho), dt]);
        r.lorentz(&[2.0, rng.uniform(0.0, 10.0), std::f64::consts::SQRT_2 / dt * (1.0 + rho),
                    rng.log_uniform(11.0, 15.0), dt]);
        let ei = rng.uniform(1.0, 10.0);
        r.debye(&[ei, rng.uniform(ei, 100.0), dt / 2.0 * (1.0 + rho), dt]);
    }
}

// ── part 2a-1b: the Mur coefficients where they cancel (ADR 0002 row 40) ──────────────────────────────────────────

/// The Mur coefficients where Rust's cancel (the port rounds the TF/SF one once and gives the 1D one a
/// cancellation-free form, ADR 0002 row 40): the TF/SF source at the magic time step S = c·dt/dx = 1 (dt = dx/c ulp by
/// ulp and S = 1 − 10⁻ᵏ), and the 1D grid where v·dt = c·dt/√(ε_rμ_r) meets dx, ε_rμ_r = 1/4 at dt = dx/(2c) (relative
/// offsets either side and ulp by ulp, through ε_r, μ_r or both), and at its band's edges ε_rμ_r = 1/16 and 9/16; each
/// as constructor or setup rows and as Gaussian-driven scenarios.
impl Rows {
    /// As `tf`, skipping a (function, inputs) an earlier part emitted.
    fn tf_new(&mut self, func: &str, sc: &TfScript) {
        if !self.keys.contains(&format!("{func} {:?}", sc.ins)) {
            self.tf(func, sc);
        }
    }

    /// As `one`, skipping a (function, inputs) already emitted.
    fn one_new(&mut self, func: &str, n: usize, dx: f64, setup: &[Op]) {
        let mut sc = Script::new(n, dx);
        for &op in setup {
            sc.apply(op);
        }
        if !self.keys.contains(&format!("{func} {:?}", sc.ins)) {
            self.emit(func, &sc);
        }
    }
}

fn emit_2a1b(r: &mut Rows) {
    let c = bijli::field::SPEED_OF_LIGHT;
    let offs: Vec<f64> = (1..=15).map(|k| 10f64.powi(-k)).collect();
    for dx in [1e-4, 1e-3, 3e-7, 0.01, 2f64.powi(-20)] {
        let m = dx / c;
        let mut dts = ulp_walk(m, 8);
        dts.extend(offs.iter().map(|&o| m * (1.0 - o)));
        for dt in dts {
            r.tf_new("TfsfSource::try_new", &TfScript::new(true, 100, 50, dx, dt));
        }
        r.tf_frames("TfsfSource::step", (200, 100, dx, m), &[], (40.0 * m, 10.0 * m), &[1, 60, 150, 400]);
    }
    let st = "Fdtd1d::step_once";
    let sf = "Fdtd1d::set_permittivity";
    let quarter: Vec<f64> = offs.iter().flat_map(|&o| [0.25 * (1.0 + o), 0.25 * (1.0 - o)]).chain(ulp_walk(0.25, 8))
        .collect();
    for &q in &quarter {
        r.one_new(sf, 20, 1e-3, &[Op::Eps(0, 20, q)]);
        r.one_new("Fdtd1d::set_permeability", 20, 1e-3, &[Op::Mu(0, 20, q)]);
        r.one_new(sf, 20, 1e-3, &[Op::Mu(0, 20, 0.5), Op::Eps(0, 20, 2.0 * q)]);
    }
    for q in [0.25, 1.0 / 16.0, 9.0 / 16.0] {
        for e in ulp_walk(q, 4) {
            r.one_new(sf, 20, 1e-3, &[Op::Eps(0, 20, e)]);
        }
    }
    for (dx, q) in [(1e-3, 0.25), (1e-4, 0.25 * (1.0 + 1e-9)), (2f64.powi(-12), 0.3), (0.05, 0.2)] {
        let dt = Fdtd1d::new(200, dx).expect("a grid").dt;
        r.frames(st, 200, dx, &[Op::Eps(0, 200, q)], &[1, 60, 200, 500], |k| Op::Pulse(k, 100, 40.0 * dt, 10.0 * dt));
    }
}

// ── part 2a-2: the subpixel means rounded once (ADR 0002 row 41) ─────────────────────────────────────────────────

/// The subpixel means where rounding once shows (ADR 0002 row 41): uniform cells (ε₁ = ε₂) of common media at fills
/// in tenths, thirds and quarters; adjacent doubles; fills one ulp or a subnormal from 0 and 1; media from 10^±300 to
/// DBL_MAX and the subnormals; seeded media and seeded uniform cells.
fn emit_2a2(r: &mut Rows) {
    let media = [1.0, 1.5, 2.1, 2.25, 3.0, 4.0, 9.8, 11.7, 12.25, 13.0];
    let mut fills: Vec<f64> = (0..=10).map(|k| f64::from(k) / 10.0).collect();
    fills.extend([1.0 / 3.0, 2.0 / 3.0, 0.25, 0.75]);
    for &e in &media {
        for &fr in &fills {
            r.subpixel(&[e, e, fr]);
            r.subpixel(&[e, e.next_up(), fr]);
        }
    }
    let edges = [f64::from_bits(1), 1e-300, 1e-17, 1.0f64.next_down(), 1.0 - 1e-16, 0.5f64.next_up()];
    for fr in edges {
        for (e1, e2) in [(2.0, 9.0), (1.0, 4.0), (3.0, 3.0), (9.0, 2.0)] {
            r.subpixel(&[e1, e2, fr]);
        }
    }
    let big = [f64::MAX, f64::MAX / 2.0, 1e300, 1e-300, f64::MIN_POSITIVE, f64::from_bits(1), 1e-310];
    for &e1 in &big {
        for &e2 in &big {
            for fr in [0.5, 0.3, 0.9] {
                r.subpixel(&[e1, e2, fr]);
            }
        }
    }
    let mut rng = Rng::new(0x2a2_0041);
    for _ in 0..500 {
        r.subpixel(&[rng.log_uniform(-1.0, 2.0), rng.log_uniform(-1.0, 2.0), rng.unit()]);
        let e = rng.log_uniform(-1.0, 2.0);
        r.subpixel(&[e, e, rng.unit()]);
    }
}

pub fn emit() -> Module {
    let mut r = Rows { m: Module::new("fdtd"), keys: BTreeSet::new() };
    emit_1a(&mut r);
    emit_1b(&mut r);
    emit_1c(&mut r);
    emit_1d_grids(&mut r);
    emit_1d_cw(&mut r);
    emit_1d_eigen(&mut r);
    emit_1e(&mut r);
    emit_2a1a(&mut r);
    emit_2a1b(&mut r);
    emit_2a2(&mut r);
    r.m
}
