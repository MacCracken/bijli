#!/usr/bin/env python3
"""Golden rows written from the reference model (ADR 0003): the corpus's emitter since the Rust tree was retired.

Each module's file keeps its frozen rows (the retired Rust generator's, FROZEN of them) and is followed by the rows the
parts below emit, in order: their inputs from deterministic grids, their outputs from validate.py's transcription of
Rust 1.2.0 on glibc's primitives (the oracle's convention: `validate.py --write` then puts MPFR's correctly rounded
values in `<module>.cr.tsv` where they differ), their tokens from validate.py's `wire()` (format 1, README.md). A
part's rows are appended once; a later part is appended after them, so every earlier row stays byte for byte.

Usage, from the repo root:
    python3 -I -B scripts/oracle/emit.py .            # rewrite the emitted rows of each module
    python3 -I -B scripts/oracle/emit.py . --check    # exit 1 if a module's file differs from what the parts emit
"""
import math
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate as V  # noqa: E402  (the reference model; its main runs only as a script)

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
check = "--check" in sys.argv[2:]

# The rows of each module written by the retired Rust generator; never rewritten. A module absent here has none: its
# file is all emitted rows, under the header new_header() writes.
FROZEN = {"fdtd_2d3d": 1107}


# ── tokens and rows ─────────────────────────────────────────────────────────────────────────────────────────────────

def u(x: int):
    """A usize as its row token (its 64-bit word)."""
    return ("i", x - (1 << 64) if x >= 1 << 63 else x)


def i(x: int):
    return ("i", x)


def f(x: float):
    return ("f", float(x))


NONE = ("n", "")


def human(tok) -> str:
    """A token for the note, as the frozen rows render it (Rust's `{:?}`)."""
    k, b = tok
    if k == "f":
        return V.rust_debug(b)
    if k == "i":
        return str(b)
    if k == "s":
        return '"' + b.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if k == "n":
        return "None"
    if k == "E":
        return f"Err({b})"
    return f"{k}:{b}"


def wire(toks) -> str:
    return V.wire(toks) if toks else "-"


class Rows:
    """One module's emitted rows, each (function, inputs) once, none repeating a frozen row."""

    def __init__(self, frozen_keys):
        self.lines, self.keys = [], set(frozen_keys)

    def emit(self, fn, cls, ins, outs):
        key = (fn, wire(ins))
        if key in self.keys:
            return
        self.keys.add(key)
        note = f"{fn}({', '.join(human(t) for t in ins)}) -> {', '.join(human(t) for t in outs)}"
        self.lines.append(f"{fn}\t{cls}\t{wire(ins)}\t{wire(outs)}\t{note}")


class Rng:
    """A seeded xorshift64*: draws for the grids (integers and plain floats only)."""

    def __init__(self, seed: int):
        self.s = seed or 0x9E3779B97F4A7C15

    def bits(self) -> int:
        x = self.s
        x ^= (x >> 12)
        x ^= (x << 25) & V.FD_U64
        x ^= (x >> 27)
        self.s = x
        return (x * 0x2545F4914F6CDD1D) & V.FD_U64

    def unit(self) -> float:
        return (self.bits() >> 11) / float(1 << 53)

    def int(self, lo: int, hi: int) -> int:
        """An integer in [lo, hi]."""
        return lo + self.bits() % (hi - lo + 1)

    def log_uniform(self, a: float, b: float) -> float:
        """10^e for e uniform in [a, b], through a decimal string (no transcendental)."""
        e = a + (b - a) * self.unit()
        k = math.floor(e)
        m = 1.0 + 9.0 * (e - k)
        return float(f"{m:.6f}e{k}")


# ── 3D scripts (M6 part 2b on) ──────────────────────────────────────────────────────────────────────────────────────

class Script3:
    """A 3D script being built (m_fdtd_2d3d.rs's form): the grid, its operations' tokens, a row at each frame; the
    outputs from validate.py's Fd3Scripts on glibc's primitives."""

    def __init__(self, rows, model, dims, dx):
        self.rows, self.model = rows, model
        self.ins = [u(dims[0]), u(dims[1]), u(dims[2]), f(dx)]

    def op(self, *toks):
        self.ins += list(toks)

    def row(self):
        outs = self.model.eval(V.wire(self.ins))
        self.rows.emit("Fdtd3d::script", "T", self.ins, outs)


def step(k):
    return (i(6), u(k))


def source(x, y, z, v):
    return (i(5), u(x), u(y), u(z), f(v))


def cpml(t):
    return (i(9), u(t))


def cpml_with(t, m, sr, km, am):
    return (i(10), u(t), f(m), f(sr), f(km), f(am))


def disable():
    return (i(12),)


def run(k, src, fr, am):
    head = (i(7), u(k)) + ((u(src[0]), u(src[1]), u(src[2])) if src is not None else (NONE,))
    return head + (f(fr), f(am))


def box(op, c, *vals):
    return (i(op),) + tuple(u(q) for q in c) + tuple(f(v) for v in vals)


def frames3(rows, model, dims, dx, setup, frames, advance=step):
    """A 3D scenario: the setup operations, then a row at each frame (cumulative steps)."""
    sc = Script3(rows, model, dims, dx)
    for op in setup:
        sc.op(*op)
    if not frames:
        sc.row()
    at = 0
    for fr in frames:
        sc.op(*advance(fr - at))
        at = fr
        sc.row()


def emit_2b(rows, model):
    """M6 part 2b: `Fdtd3d`'s CPML. The classic layer and the CFS layer around point and driven sources, the
    gradings swept, a dielectric and lossy block across the slabs, the thinnest grids a layer fits on each axis,
    sources inside the slabs, enable / disable / replace, the thickness's and the grading's refusals in Rust's order,
    and seeded scripts."""
    two_pi_e0 = 2.0 * math.pi * model.C["field::EPSILON_0"]
    cfs = lambda t, km: cpml_with(t, 3.0, 1.0, km, two_pi_e0 * 1e9)  # noqa: E731
    # The bench's grid (8 cells) a few steps; fdtd.rs's reflection and energy grids (10 and 6 cells), shortened.
    frames3(rows, model, (30, 30, 30), 1e-3, [cpml(8), source(15, 15, 15, 1.0)], [1, 2])
    frames3(rows, model, (17, 17, 17), 1e-3, [cpml(4), source(8, 8, 8, 1.0)], [1, 5, 20, 40])
    frames3(rows, model, (15, 15, 15), 1e-3, [cfs(4, 3.0)], [6, 24], lambda k: run(k, (7, 7, 7), 1e10, 1.0))
    # Gradings: m, σ ratio, κ_max and α_max, on thin layers, each axis a different length.
    for op in (cpml_with(2, 2.0, 0.7, 1.0, 0.0), cpml_with(3, 4.0, 1.5, 11.0, 0.05), cpml_with(2, 1.0, 1.0, 1.0, 1.0),
               cpml_with(2, 3.5, 1.2, 3.0, 1e-3)):
        frames3(rows, model, (12, 11, 10), 1e-3, [op, source(6, 5, 5, 1.0), source(2, 2, 2, -0.5)], [1, 6, 18])
    # A dielectric and lossy block overlapping the x-low slab (fdtd.rs's long-run test, smaller and shorter).
    frames3(rows, model, (13, 13, 13), 1e-3, [cfs(3, 5.0), box(4, (0, 4, 4, 6, 9, 9), 4.0, 2.0, 0.5)], [10, 40],
            lambda k: run(k, (8, 6, 6), 3e10, 1.0))
    # The thinnest grids a layer fits (2t + 3 on one axis, then on each), sources inside the slabs and their edges.
    frames3(rows, model, (7, 7, 7), 1e-3, [cpml(2), source(3, 3, 3, 1.0)], [1, 4, 15])
    for dims in ((7, 12, 11), (12, 7, 11), (12, 11, 7), (9, 9, 13)):
        ops = [cpml(2), source(dims[0] // 2, dims[1] // 2, dims[2] // 2, 1.0), source(1, 1, 1, 0.5),
               source(dims[0] - 2, 1, dims[2] - 3, -0.25)]
        frames3(rows, model, dims, 1e-3, ops, [1, 8, 20])
    # Every field filled, then stepped with a layer: each slab term sees a uniform source field.
    for w in range(6):
        frames3(rows, model, (9, 8, 7), 1e-3, [cpml(2), (i(8), u(w), f(1.0 + w))], [1, 3])
    # Enable, disable, replace: Ψ restarts at zero, a refused request keeps the layer.
    ops = [source(7, 7, 7, 1.0), step(3), cpml(3), step(4), cpml(100), step(2), cpml(4), step(3), disable(),
           step(2), cfs(2, 7.0), cpml(0)]
    frames3(rows, model, (15, 14, 13), 1e-3, ops, [1, 5])
    frames3(rows, model, (11, 11, 11), 1e-3, [cpml(2), step(4), disable(), source(5, 5, 5, 1.0)], [3])
    # Lossy cells inside and outside the slabs.
    frames3(rows, model, (14, 13, 12), 1e-3, [cpml(3), box(3, (0, 0, 0, 14, 4, 12), 2.0),
                                              box(3, (5, 5, 4, 9, 8, 8), 50.0), source(7, 9, 6, 1.0)], [1, 12])
    # The thickness's refusals (the fit first), then the grading's, in Rust's order.
    for t in (0, 1, 4, 5, 6, 100, (1 << 64) - 1):
        frames3(rows, model, (11, 13, 15), 1e-3, [cpml(t)], [1])
    for op in (cpml_with(2, 0.5, 1.0, 1.0, 0.0), cpml_with(2, 3.0, 0.0, 1.0, 0.0), cpml_with(2, 3.0, 1.0, 0.5, 0.0),
               cpml_with(2, 3.0, 1.0, 1.0, -1.0), cpml_with(2, math.nan, math.nan, 1.0, 0.0),
               cpml_with(0, 0.5, 0.0, 0.5, -1.0), cpml_with(9, 0.5, 1.0, 1.0, 0.0),
               cpml_with(2, 3.0, sys.float_info.max, 1.0, 0.0), cpml_with(2, 3.0, 1.0, math.inf, 0.0)):
        frames3(rows, model, (9, 9, 9), 1e-3, [op, source(4, 4, 4, 1.0)], [2])
    # Seeded scripts: a fitting layer of random grading on grids from 7 to 14 a side, random materials and sources.
    rng = Rng(0x2B_3D_C9)
    for _ in range(24):
        dims = tuple(rng.int(7, 14) for _ in range(3))
        tmax = max(2, (min(dims) - 3) // 2)
        t = rng.int(2, tmax)
        m = 1.0 + 3.0 * rng.unit()
        km = 1.0 if rng.unit() < 0.5 else 1.0 + 10.0 * rng.unit()
        am = 0.0 if rng.unit() < 0.5 else rng.log_uniform(-4.0, 0.0)
        ops = [cpml_with(t, m, 0.5 + rng.unit(), km, am)]
        if rng.unit() < 0.5:
            sg = 0.0 if rng.unit() < 0.5 else 1.0
            ops.append(box(4, (0, 0, 0, dims[0] // 2, dims[1], dims[2]), 1.0 + 3.0 * rng.unit(), 1.0, sg))
        for _ in range(rng.int(1, 3)):
            at = tuple(rng.int(0, d - 1) for d in dims)
            ops.append(source(*at, (-1.0 if rng.unit() < 0.5 else 1.0) * rng.log_uniform(-2.0, 1.0)))
        frames3(rows, model, dims, rng.log_uniform(-5.0, -2.0), ops, [rng.int(1, 30)])


# ── compute.rs (M6 Part 3) ──────────────────────────────────────────────────────────────────────────────────────

class Script2:
    """A `GpuFdtd2d::script` being built: `i:<nx> i:<ny> f:<dx>`, then the 2D scripts' add_source (5), step_once (6)
    and run (7); the outputs from validate.py's Gpu2Scripts (Fd2's TM physics, which GpuFdtd2d on CpuBackend is)."""

    def __init__(self, rows, model, nx, ny, dx):
        self.rows, self.model = rows, model
        self.ins = [u(nx), u(ny), f(dx)]

    def op(self, *toks):
        self.ins += list(toks)

    def row(self):
        self.rows.emit("GpuFdtd2d::script", "T", self.ins, self.model.eval(V.wire(self.ins)))


def source2(x, y, v):
    return (i(5), u(x), u(y), f(v))


def run2(k, src, fr, am):
    head = (i(7), u(k)) + ((u(src[0]), u(src[1])) if src is not None else (NONE,))
    return head + (f(fr), f(am))


def frames_gpu(rows, model, nx, ny, dx, setup, frames, advance=step):
    sc = Script2(rows, model, nx, ny, dx)
    for op in setup:
        sc.op(*op)
    if not frames:
        sc.row()
    at = 0
    for fr in frames:
        sc.op(*advance(fr - at))
        at = fr
        sc.row()


def kernel(rows, op, nx, ny, arg, lens, seed):
    """A `CpuBackend::kernel` row (class A): `i:<op> i:<nx> i:<ny> i:<axis|index> i:<len0> .. i:<len3> f:<seed>`, four
    arrays filled by validate.py's cpu_fill (j, len_j, seed), then the operation: 1 fdtd_update_h(b0, b1, b2, nx, ny,
    axis), 2 fdtd_update_e(b0, b1, b2, b3, nx, ny), 3 apply_zero_boundary(b0, nx, ny), 4 read_buffer(b0) into the host
    array 1, 5 write_buffer(b0) from the host array 1, 6 CpuBackend's add_to_element(b0, index, seed), 7 the trait's
    default add_to_element through the host array 1 as scratch; its outputs the error, or the four arrays' digests."""
    ins = [i(op), u(nx), u(ny), u(arg)] + [u(q) for q in lens] + [f(seed)]
    rows.emit("CpuBackend::kernel", "A", ins, V.cpu_kernel([t[1] for t in ins]))


def emit_3(rows, model):
    """M6 Part 3: GpuFdtd2d on CpuBackend (compute.rs's tests' grids, Rust's 2D grid's shapes, sources inside, on the
    wall and outside, run with and without a source, seeded scripts) and the CPU backend's kernels on every shape and
    every refusal, in Rust's order."""
    mx = (1 << 64) - 1
    for nx, ny, dx in ((2, 50, 1e-3), (50, 2, 1e-3), (0, 0, 1e-3), (2, 2, math.nan), (1, mx, 1e-3), (mx, 3, 1e-3),
                       (1 << 32, 1 << 32, 1e-3), (10, 10, math.nan), (10, 10, math.inf), (10, 10, -1.0),
                       (10, 10, 0.0), (10, 10, -0.0)):
        frames_gpu(rows, model, nx, ny, dx, [], [])
    frames_gpu(rows, model, 50, 50, 1e-3, [], [])
    frames_gpu(rows, model, 30, 30, 1e-4, [source2(15, 15, 1.0)], [1, 5, 20])
    frames_gpu(rows, model, 40, 40, 1e-4, [], [10, 20, 30], lambda k: run2(k, (20, 20), 1e9, 1.0))
    frames_gpu(rows, model, 30, 30, 1e-4, [], [25], lambda k: run2(k, (12, 17), 2e9, 1.0))
    for src in ((20, 5), (0, 20), (mx, 0)):
        frames_gpu(rows, model, 20, 20, 1e-4, [], [10], lambda k, src=src: run2(k, src, 1e9, 1.0))
    frames_gpu(rows, model, 10, 10, 1e-3, [source2(3, 4, 2.5), source2(3, 4, 0.5), source2(10, 4, 1.0)], [0, 1])
    for nx, ny in ((3, 3), (3, 5), (5, 3), (4, 4), (7, 5), (16, 9), (9, 16)):
        ops = [source2(nx // 2, ny // 2, 1.0), source2(1, 1, -0.5), source2(0, ny // 2, 2.0)]
        frames_gpu(rows, model, nx, ny, 1e-3, ops, [1, 3, 12])
    frames_gpu(rows, model, 24, 24, 1e-3, [source2(12, 12, 1.0)], [4, 15], lambda k: run2(k, None, 1e9, 1.0))
    frames_gpu(rows, model, 100, 100, 1e-4, [source2(50, 50, 1.0)], [1, 3])
    rng = Rng(0x3C_0FF)
    for _ in range(16):
        nx, ny = rng.int(3, 24), rng.int(3, 24)
        ops = [source2(rng.int(0, nx), rng.int(0, ny), (-1.0 if rng.unit() < 0.5 else 1.0) * rng.log_uniform(-2.0, 1.0))
               for _ in range(rng.int(1, 3))]
        frames_gpu(rows, model, nx, ny, rng.log_uniform(-5.0, -2.0), ops, [rng.int(1, 25)],
                   (lambda k, s=(rng.int(0, nx - 1), rng.int(0, ny - 1)), fr=rng.log_uniform(8.0, 11.0):
                    run2(k, s, fr, 1.0)) if rng.unit() < 0.5 else step)
    # The kernels on every small shape (1 or 2 wide included), each with its buffers nx·ny long.
    shapes = ((3, 3), (4, 5), (5, 4), (1, 4), (4, 1), (2, 2), (1, 1), (2, 7), (17, 9), (40, 3))
    for nx, ny in shapes:
        n = nx * ny
        for seed in (1.0, 0.37, -2.5):
            kernel(rows, 1, nx, ny, 0, (n, n, n, n), seed)
            kernel(rows, 1, nx, ny, 1, (n, n, n, n), seed)
            kernel(rows, 2, nx, ny, 0, (n, n, n, n), seed)
            kernel(rows, 3, nx, ny, 0, (n, n, n, n), seed)
    # The refusals, in Rust's order: the axis first, then the grid (a zero or overflowing nx·ny), then the first
    # buffer of the wrong length.
    for axis in (2, 3, 255):
        kernel(rows, 1, 3, 3, axis, (9, 9, 9, 9), 1.0)
        kernel(rows, 1, 0, 3, axis, (9, 9, 9, 9), 1.0)
    for nx, ny in ((3, 0), (0, 3), (0, 0), (mx, 2), (2, mx), (1 << 32, 1 << 32)):
        for op in (1, 2, 3):
            kernel(rows, op, nx, ny, 0, (9, 9, 9, 9), 1.0)
    for lens in ((8, 9, 9, 9), (9, 8, 9, 9), (9, 9, 8, 9), (9, 9, 9, 8), (16, 16, 16, 16), (8, 8, 8, 8)):
        for op in (1, 2, 3):
            kernel(rows, op, 3, 3, 1, lens, 0.5)
    # The transfers and add_to_element: equal and unequal lengths, indices to and past the end.
    for lens in ((0, 0, 0, 0), (3, 3, 0, 0), (9, 9, 0, 0), (3, 2, 0, 0), (3, 4, 0, 0), (0, 1, 0, 0)):
        kernel(rows, 4, 1, 1, 0, lens, 0.25)
        kernel(rows, 5, 1, 1, 0, lens, 0.25)
    for n, idx in ((4, 0), (4, 3), (4, 4), (4, 5), (4, mx), (1, 0), (0, 0), (9, 8)):
        kernel(rows, 6, 1, 1, idx, (n, n, 0, 0), 1.5)
        kernel(rows, 7, 1, 1, idx, (n, n, 0, 0), 1.5)
    for idx in (0, 3, 4):
        kernel(rows, 7, 1, 1, idx, (4, 3, 0, 0), 1.5)
        kernel(rows, 7, 1, 1, idx, (4, 5, 0, 0), 1.5)


# ── M6 Part 4: the range-edge review ────────────────────────────────────────────────────────────────────────────────

class Script2d:
    """An `Fdtd2d::script` being built: `i:<nx> i:<ny> f:<dx> i:<mode>`, then the 2D scripts' operations; the outputs
    from validate.py's Fd2Scripts on glibc's primitives."""

    def __init__(self, rows, model, nx, ny, dx, mode):
        self.rows, self.model = rows, model
        self.ins = [u(nx), u(ny), f(dx), i(mode)]

    def op(self, *toks):
        self.ins += list(toks)

    def row(self):
        self.rows.emit("Fdtd2d::script", "T", self.ins, self.model.eval(V.wire(self.ins)))


def region2(op, c, *vals):
    return (i(op),) + tuple(u(q) for q in c) + tuple(f(v) for v in vals)


def source2d(x, y, v):
    return (i(5), u(x), u(y), f(v))


def frames2(rows, model, nx, ny, dx, mode, setup, frames, advance=step):
    """A 2D scenario: the setup operations, then a row at each frame (cumulative steps)."""
    sc = Script2d(rows, model, nx, ny, dx, mode)
    for op in setup:
        sc.op(*op)
    if not frames:
        sc.row()
    at = 0
    for fr in frames:
        sc.op(*advance(fr - at))
        at = fr
        sc.row()


# The spacings at the range's ends: past them Rust's ε₀·dx, (2√2)·c·… or dt itself left the normal range.
EDGE_DX = (sys.float_info.max, 1e300, 1e200, 1.0, 0.032, 1e-200, 1e-300, 1e-305, 2.2250738585072014e-308, 1e-310,
           5e-324)


def emit_4a(rows, m2, m3, mg):
    """M6 part 4a (ADR 0002 rows 13 and 40): the grids' dt, update coefficients and lossy C_a and C_b rounded once.
    Seeded 2D (both modes), 3D and GpuFdtd2d grids over ten decades of dx, each through seeded materials; lossy cells
    where C_a's numerator cancels (σΔt/2ε within 10⁻¹⁵ to 10⁻³ of 1); the spacings at the range's ends; and CPML
    layers on the rounded dt."""
    c, eps0 = m2.K[0], m2.K[1]
    # The spacings at the range's ends, a source in each grid.
    for dx in EDGE_DX:
        for mode in (0, 1):
            frames2(rows, m2, 7, 6, dx, mode, [source2d(3, 3, 1.0), source2d(2, 4, -0.5)], [1, 3])
        frames3(rows, m3, (5, 6, 5), dx, [source(2, 3, 2, 1.0)], [1, 3])
        frames_gpu(rows, mg, 7, 6, dx, [source2(3, 3, 1.0)], [1, 3])
    # Lossy cells where C_a = (1 − s)/(1 + s) cancels: σ = (1 + δ)·2ε₀ε_r/dt, so s = σΔt/2ε is 1 + δ to an ulp or two.
    for delta in (0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-9, 1e-6, -1e-3):
        for eps_r in (1.0, 2.5):
            dt2 = V.fdiv(1e-3, 2.0 * math.sqrt(2.0) * c)
            sg = (1.0 + delta) * 2.0 * eps0 * eps_r / dt2
            ops = [region2(4, (2, 2, 6, 5), eps_r, 1.0, sg), source2d(4, 3, 1.0)]
            frames2(rows, m2, 9, 8, 1e-3, 0 if eps_r == 1.0 else 1, ops, [1, 4])
            dt3 = V.fdiv(1e-3, 2.0 * math.sqrt(3.0) * c)
            sg3 = (1.0 + delta) * 2.0 * eps0 * eps_r / dt3
            frames3(rows, m3, (6, 6, 6), 1e-3, [box(4, (1, 1, 1, 5, 5, 4), eps_r, 1.0, sg3), source(3, 3, 2, 1.0)],
                    [1, 3])
    # Seeded 2D grids: dx over ten decades, both modes, up to three material regions, sources.
    rng = Rng(0x4A_2D)
    for _ in range(64):
        nx, ny, mode = rng.int(3, 12), rng.int(3, 12), rng.int(0, 1)
        ops = []
        for _ in range(rng.int(0, 3)):
            x0, y0 = rng.int(0, nx - 1), rng.int(0, ny - 1)
            cc = (x0, y0, rng.int(x0, nx + 1), rng.int(y0, ny + 1))
            kind = rng.int(1, 4)
            eps, mu = rng.log_uniform(0.0, 2.0), rng.log_uniform(0.0, 1.0)
            sg = 0.0 if rng.unit() < 0.3 else rng.log_uniform(-3.0, 3.0)
            ops.append(region2(kind, cc, *{1: (eps,), 2: (mu,), 3: (sg,), 4: (eps, mu, sg)}[kind]))
        for _ in range(rng.int(1, 2)):
            ops.append(source2d(rng.int(0, nx - 1), rng.int(0, ny - 1), rng.log_uniform(-2.0, 1.0)))
        frames2(rows, m2, nx, ny, rng.log_uniform(-9.0, 1.0), mode, ops, [rng.int(1, 6)])
    # Seeded 3D grids, as the 2D ones.
    rng = Rng(0x4A_3D)
    for _ in range(32):
        dims = tuple(rng.int(3, 7) for _ in range(3))
        ops = []
        for _ in range(rng.int(0, 3)):
            lo = tuple(rng.int(0, d - 1) for d in dims)
            cc = lo + tuple(rng.int(lo[k], dims[k] + 1) for k in range(3))
            kind = rng.int(1, 4)
            eps, mu = rng.log_uniform(0.0, 2.0), rng.log_uniform(0.0, 1.0)
            sg = 0.0 if rng.unit() < 0.3 else rng.log_uniform(-3.0, 3.0)
            ops.append(box(kind, cc, *{1: (eps,), 2: (mu,), 3: (sg,), 4: (eps, mu, sg)}[kind]))
        ops.append(source(*(rng.int(0, d - 1) for d in dims), rng.log_uniform(-2.0, 1.0)))
        frames3(rows, m3, dims, rng.log_uniform(-9.0, 1.0), ops, [rng.int(1, 4)])
    # Seeded GpuFdtd2d grids.
    rng = Rng(0x4A_6C)
    for _ in range(24):
        nx, ny = rng.int(3, 16), rng.int(3, 16)
        ops = [source2(rng.int(0, nx - 1), rng.int(0, ny - 1), rng.log_uniform(-2.0, 1.0))]
        frames_gpu(rows, mg, nx, ny, rng.log_uniform(-9.0, 1.0), ops, [rng.int(1, 8)])
    # CPML layers on the rounded dt: the profiles read the grid's dt.
    rng = Rng(0x4A_C9)
    for _ in range(6):
        dx = rng.log_uniform(-6.0, -1.0)
        frames2(rows, m2, 11, 10, dx, rng.int(0, 1), [cpml(2), source2d(5, 5, 1.0)], [3])
    for _ in range(4):
        frames3(rows, m3, (9, 8, 9), rng.log_uniform(-6.0, -1.0), [cpml(2), source(4, 4, 4, 1.0)], [2])


def layer(rows, C, t, m, sr, km, am, dx, dt):
    """A `Cpml::with_params` row (class C): `i:<t> f:<m> f:<sigma_ratio> f:<kappa_max> f:<alpha_max> f:<dx> f:<dt>`."""
    ins = [u(t), f(m), f(sr), f(km), f(am), f(dx), f(dt)]
    rows.emit("Cpml::with_params", "C", ins, V.fd2_eval("Cpml::with_params", V.wire(ins), V.GLIBC_P, C))


def emit_4b(rows, m2, m3):
    """M6 part 4b (ADR 0002 rows 40 and 18): the CPML profile. Thick layers and high orders on the 2D and 3D steps,
    where c's b − 1 and 1/κ − 1 cancel in the inner layers (κ_max near 1 for the second); σ_max past either end of the
    range (sigma_ratio and dx from 2^-1074 to DBL_MAX, m past it); grades below 2^-1022 (m to DBL_MAX); α_max, κ_max
    and dt at their extremes; seeded gradings over the whole range; and grids whose layers κ_max near 1 grades."""
    C = m2.C
    mx, tiny = sys.float_info.max, 5e-324
    c0 = C["field::SPEED_OF_LIGHT"]
    dt2 = V.rn_mono([1e-3], [2.0 * math.sqrt(2.0), c0])
    dt3 = V.rn_mono([1e-3], [2.0 * math.sqrt(3.0), c0])
    # The cancellation: thick layers, high orders, CFS κ_max near 1.
    for dt in (dt2, dt3):
        for t, m in ((16, 3.0), (20, 4.0), (32, 4.0), (40, 6.0), (64, 3.0)):
            layer(rows, C, t, m, 1.0, 1.0, 0.0, 1e-3, dt)
        for km in (1.0 + 1e-12, 1.0 + 1e-8, 1.001, 1.5, 3.0, 15.0):
            for am in (0.0, 1e-3, 0.05):
                layer(rows, C, 12, 4.0, 1.0, km, am, 1e-3, dt)
    # σ_max past either end: sigma_ratio and dx over the range, m past DBL_MAX/2.
    for sr in (tiny, 1e-310, 2.2250738585072014e-308, 1.0, 1e300, mx):
        for dx in (tiny, 1e-310, 1e-300, 1e-3, 1e300, mx):
            layer(rows, C, 3, 3.0, sr, 1.0, 0.0, dx, 1e-12)
    for m in (1e300, mx):
        for sr in (1e-300, 1.0, 1e10):
            layer(rows, C, 3, m, sr, 1.0, 0.0, 1e-3, dt2)
    # Grades below 2^-1022: high orders on a fine dx (σ_max large enough that σ stays normal), CFS and classic.
    for m in (200.0, 300.0, 500.0, 1000.0, 1e5, 1e300):
        for t in (4, 10):
            layer(rows, C, t, m, 1.0, 1.0, 0.0, 1e-300, 1e-310)
            layer(rows, C, t, m, 1e-200, 1e100, 0.01, 1e-3, dt2)
    # α_max, κ_max and dt at their extremes; σ/κ + α past DBL_MAX on a dt that brings it back.
    for am in (tiny, 1e-310, 1e300, mx):
        layer(rows, C, 4, 3.0, 1.0, 1.0, am, 1e-3, dt2)
        layer(rows, C, 4, 3.0, 1.0, 5.0, am, 1e-3, tiny)
    for km in (1e300, mx):
        layer(rows, C, 4, 3.0, 1.0, km, 0.0, 1e-3, dt2)
        layer(rows, C, 4, 3.0, 1.0, km, 1.0, 1e-3, dt2)
    for dt in (tiny, 1e-310, 1e-300, 1e300, mx):
        layer(rows, C, 4, 3.0, 1.0, 2.0, 0.01, 1e-3, dt)
    for dt in (tiny, 1e-320, 1e-300):
        layer(rows, C, 3, 3.0, mx, 1.0, mx, tiny, dt)
        layer(rows, C, 3, 3.0, 1e300, 1.0, mx, 1e-20, dt)
    # Seeded gradings over the whole range.
    rng = Rng(0x4B_C9)
    for _ in range(64):
        t = rng.int(1, 12)
        m = 1.0 + rng.log_uniform(-2.0, 3.0)
        km = 1.0 if rng.unit() < 0.3 else 1.0 + rng.log_uniform(-15.0, 300.0)
        am = 0.0 if rng.unit() < 0.3 else rng.log_uniform(-320.0, 300.0)
        layer(rows, C, t, m, rng.log_uniform(-300.0, 300.0), km, am, rng.log_uniform(-300.0, 300.0),
              rng.log_uniform(-320.0, 300.0))
    # Grids whose layers κ_max near 1 grades: 1/κ − 1 cancels in every layer.
    for km in (1.0 + 1e-9, 1.001):
        frames2(rows, m2, 15, 14, 1e-3, 0, [cpml_with(4, 4.0, 1.0, km, 0.0), source2d(7, 7, 1.0)], [6])
        frames3(rows, m3, (11, 11, 11), 1e-3, [cpml_with(3, 4.0, 1.0, km, 0.0), source(5, 5, 5, 1.0)], [4])


def fill2(which, v):
    return (i(8), i(which), f(v))


def emit_4c(rows, m2, mg):
    """M6 part 4c (ADR 0002 row 18): the 2D grid's and the CPU backend's per-cell updates at their range edges, in both
    modes: neighbours of ±DBL_MAX (their differences past it), of 2^-1074 to 2^-1022 (products below the normal range)
    and of opposite signs (products that round to ±0 on ±0 fields), lossy cells among them, seeded scenarios, the same
    through GpuFdtd2d, and the CPU kernels on buffers near DBL_MAX."""
    mx, tiny, dmin = sys.float_info.max, 5e-324, 2.2250738585072014e-308
    big = ((mx, -mx), (mx, mx), (1e308, -1e308), (8.98846567431158e307, -8.98846567431158e307), (mx, 1e300))
    small = ((tiny, -tiny), (tiny, 0.0), (1e-310, -1e-310), (dmin, -dmin), (-1e-320, 3e-321), (dmin, -tiny))
    for mode in (0, 1):
        for a, b in big + small:
            frames2(rows, m2, 9, 8, 1e-3, mode, [source2d(4, 4, a), source2d(5, 4, b), source2d(4, 5, b)], [1, 2, 4])
        # The ±0 fields: a product of a subnormal difference rounds to ∓0 against them.
        for z in (0.0, -0.0):
            for w in (0, 1, 2):
                frames2(rows, m2, 8, 7, 1e-3, mode, [fill2(w, z), source2d(3, 3, -tiny), source2d(4, 3, tiny)], [1, 3])
        # Lossy cells beside the edges: C_a·E first, then the curl term.
        ops = [region2(3, (2, 2, 7, 6), 50.0), source2d(4, 4, mx), source2d(5, 4, -mx), source2d(3, 3, tiny)]
        frames2(rows, m2, 9, 8, 1e-3, mode, ops, [1, 3])
    # Seeded scenarios: sources over every binade of either sign.
    rng = Rng(0x4C_2D)
    for _ in range(40):
        nx, ny, mode = rng.int(3, 10), rng.int(3, 10), rng.int(0, 1)
        ops = []
        for _ in range(rng.int(2, 5)):
            v = rng.log_uniform(-323.0, 308.0)
            ops.append(source2d(rng.int(0, nx - 1), rng.int(0, ny - 1), -v if rng.unit() < 0.5 else v))
        frames2(rows, m2, nx, ny, rng.log_uniform(-6.0, -1.0), mode, ops, [rng.int(1, 4)])
    # GpuFdtd2d: the same kernels through the trait.
    for a, b in big + small:
        frames_gpu(rows, mg, 9, 8, 1e-3, [source2(4, 4, a), source2(5, 4, b), source2(4, 5, b)], [1, 3])
    # The CPU kernels on buffers near DBL_MAX: every difference, curl and product at its edge.
    for seed in (1.7e307, 8.98846567431158e307 / 5, 1e307, 3e306):
        for nx, ny in ((3, 3), (5, 4), (9, 7)):
            n = nx * ny
            kernel(rows, 1, nx, ny, 0, (n, n, n, n), seed)
            kernel(rows, 1, nx, ny, 1, (n, n, n, n), seed)
            kernel(rows, 2, nx, ny, 0, (n, n, n, n), seed)
            kernel(rows, 2, nx, ny, 0, (n, n, n, n), -seed)


def emit_4d(rows, m3):
    """M6 part 4d (ADR 0002 row 18): the 3D grid's six curl updates at their range edges: neighbours of ±DBL_MAX, of
    2^-1074 to 2^-1022 and of opposite signs on E_z (H_x and H_y from it, then E_x, E_y and E_z from H, then H_z), on
    a grid whose rows have a packed part and a scalar tail; ±0 fields under subnormal sources; lossy cells among the
    edges; materials whose coefficients let an E curl pass DBL_MAX; and seeded scenarios. Field fills take Fd3.FIELDS's
    index (E_x, E_y, E_z, H_x, H_y, H_z)."""
    mx, tiny, dmin = sys.float_info.max, 5e-324, 2.2250738585072014e-308
    big = ((mx, -mx), (mx, mx), (1e308, -1e308), (8.98846567431158e307, -8.98846567431158e307), (mx, 1e300))
    small = ((tiny, -tiny), (tiny, 0.0), (1e-310, -1e-310), (dmin, -dmin), (-1e-320, 3e-321), (dmin, -tiny))
    dims = (9, 7, 7)
    for a, b in big + small:
        ops = [source(4, 3, 3, a), source(5, 3, 3, b), source(4, 4, 3, b), source(4, 3, 4, b)]
        frames3(rows, m3, dims, 1e-3, ops, [1, 2, 4])
    # The ±0 fields: a product of a subnormal difference rounds to ∓0 against them.
    for z in (0.0, -0.0):
        for w in range(6):
            frames3(rows, m3, (8, 6, 6), 1e-3, [fill2(w, z), source(3, 3, 3, -tiny), source(4, 3, 3, tiny)], [1, 3])
    # Lossy cells beside the edges: C_a·E first, then the curl term.
    ops = [box(3, (2, 2, 2, 7, 6, 6), 50.0), source(4, 3, 3, mx), source(5, 3, 3, -mx), source(3, 3, 3, tiny)]
    frames3(rows, m3, dims, 1e-3, ops, [1, 3])
    # Coefficients that let a curl pass DBL_MAX where its update doesn't (μ_r ≈ 10⁻³: C_h ≈ 0.6-0.77; ε_r = 1000:
    # C_e ≈ 0.11), E_z pairs at ±DBL_MAX/2 in the E rows' packed part (x = 3, 4) and their scalar tail (x = 6, 7).
    h = 8.98846567431158e307
    for mu in (1.28e-3, 1e-3):
        ops = [box(1, (0, 0, 0, 9, 7, 7), 1000.0), box(2, (0, 0, 0, 9, 7, 7), mu), source(3, 3, 3, h),
               source(4, 3, 3, -h), source(6, 3, 3, h), source(7, 3, 3, -h)]
        frames3(rows, m3, dims, 1e-3, ops, [1, 2, 3])
    # Seeded scenarios: sources over every binade of either sign.
    rng = Rng(0x4D_3D)
    for _ in range(30):
        d3 = (rng.int(3, 9), rng.int(3, 8), rng.int(3, 8))
        ops = []
        for _ in range(rng.int(2, 5)):
            v = rng.log_uniform(-323.0, 308.0)
            ops.append(source(rng.int(0, d3[0] - 1), rng.int(0, d3[1] - 1), rng.int(0, d3[2] - 1),
                              -v if rng.unit() < 0.5 else v))
        frames3(rows, m3, d3, rng.log_uniform(-6.0, -1.0), ops, [rng.int(1, 4)])


def emit_4e(rows, m2, m3):
    """M6 part 4e (ADR 0002 row 18): the CPML corrections' two statements at their range edges. 2D (both modes, t = 3
    on 12 × 11): pairs across an x-strip (scalar), a y-strip (packed, with the E strips' scalar tail) and a corner; 3D
    (t = 2 on 9³): pairs across the x-slabs (scalar) and the y-slabs (packed), frames on to the z-slabs. Each with the
    classic layer and a CFS one (κ_max = 5, so 1/κ − 1 ≠ 0, and α_max = 0.05): neighbours of ±DBL_MAX, of 2^-1074 to
    2^-1022, ±0 fields under subnormal sources; then layers whose coefficients reach the range's ends (σ ratios of
    10^±300, a grade below 2^-1022 that makes c and 1/κ − 1 subnormal, κ_max and α_max of 10^300) under unit sources;
    then seeded scenarios. Field fills take each grid's field index."""
    mx, tiny, dmin = sys.float_info.max, 5e-324, 2.2250738585072014e-308
    big = ((mx, -mx), (1e308, -1e308), (8.98846567431158e307, -8.98846567431158e307))
    small = ((tiny, -tiny), (1e-310, -1e-310), (dmin, -dmin), (-1e-320, 3e-321))
    for mode in (0, 1):
        for lay in (cpml(3), cpml_with(3, 4.0, 1.0, 5.0, 0.05)):
            for a, b in big + small:
                ops = [lay, source2d(1, 5, a), source2d(2, 5, b), source2d(6, 1, a), source2d(6, 2, b),
                       source2d(1, 1, a), source2d(9, 8, b)]
                frames2(rows, m2, 12, 11, 1e-3, mode, ops, [1, 2, 4])
            for z in (0.0, -0.0):
                for w in (0, 1, 2):
                    ops = [lay, fill2(w, z), source2d(1, 5, -tiny), source2d(6, 1, tiny), source2d(9, 9, -tiny)]
                    frames2(rows, m2, 12, 11, 1e-3, mode, ops, [1, 3])
    ends = ((3, 3.0, 1e300, 1.0, 0.0), (3, 3.0, 1e-300, 1.0, 0.0), (3, 3.0, 1e-305, 5.0, 0.05),
            (4, 350.0, 1.0, 2.0, 0.0), (4, 350.0, 1.0, 1.0, 0.0), (3, 3.0, 1.0, 1e300, 0.0), (3, 3.0, 1.0, 1.0, 1e300))
    for t, m, sr, km, am in ends:
        for mode in (0, 1):
            ops = [cpml_with(t, m, sr, km, am), source2d(2, 5, 1.0), source2d(6, 2, -1.0), source2d(7, 5, 1e-300)]
            frames2(rows, m2, 14, 13, 1e-3, mode, ops, [1, 3, 6])
    # 3D: E_z pairs across an x-slab (x = 1, 2) and a y-slab (y = 1, 2), one near the z faces.
    for lay in (cpml(2), cpml_with(2, 4.0, 1.0, 5.0, 0.05)):
        for a, b in big + small:
            ops = [lay, source(1, 4, 4, a), source(2, 4, 4, b), source(4, 1, 4, a), source(4, 2, 4, b),
                   source(4, 4, 1, a), source(5, 4, 1, b)]
            frames3(rows, m3, (9, 9, 9), 1e-3, ops, [1, 2, 4])
        for z in (0.0, -0.0):
            for w in (0, 3):
                ops = [lay, fill2(w, z), source(1, 4, 4, -tiny), source(4, 1, 4, tiny)]
                frames3(rows, m3, (9, 9, 9), 1e-3, ops, [1, 3])
    for t, m, sr, km, am in ((2, 3.0, 1e300, 1.0, 0.0), (2, 3.0, 1e-300, 5.0, 0.05), (3, 400.0, 1.0, 2.0, 0.0),
                             (2, 3.0, 1.0, 1e300, 1e300)):
        ops = [cpml_with(t, m, sr, km, am), source(2, 4, 4, 1.0), source(4, 2, 4, -1.0), source(4, 4, 2, 1e-300)]
        frames3(rows, m3, (9, 9, 9), 1e-3, ops, [1, 3])
    # Seeded scenarios: classic and CFS layers, sources over every binade of either sign, most inside the layers.
    rng = Rng(0x4E_CB)
    for _ in range(20):
        nx, ny, mode, t = rng.int(9, 14), rng.int(9, 14), rng.int(0, 1), rng.int(1, 3)
        if rng.unit() < 0.5:
            lay = cpml(t)
        else:
            lay = cpml_with(t, 1.0 + rng.log_uniform(-1.0, 1.0), rng.log_uniform(-2.0, 2.0),
                            1.0 + rng.log_uniform(-3.0, 1.0), rng.log_uniform(-4.0, 0.0))
        ops = [lay]
        for _ in range(rng.int(2, 5)):
            v = rng.log_uniform(-323.0, 308.0)
            x = rng.int(0, t) if rng.unit() < 0.5 else rng.int(nx - 1 - t, nx - 1)
            ops.append(source2d(x, rng.int(0, ny - 1), -v if rng.unit() < 0.5 else v))
        frames2(rows, m2, nx, ny, rng.log_uniform(-6.0, -1.0), mode, ops, [rng.int(1, 4)])
    for _ in range(10):
        d3, t = (rng.int(7, 9), rng.int(7, 9), rng.int(7, 9)), rng.int(1, 2)
        ops = [cpml(t) if rng.unit() < 0.5 else cpml_with(t, 3.0, 1.0, 5.0, 0.05)]
        for _ in range(rng.int(2, 4)):
            v = rng.log_uniform(-323.0, 308.0)
            x, y, z = rng.int(0, t), rng.int(0, d3[1] - 1), rng.int(0, d3[2] - 1)
            ops.append(source(x, y, z, -v if rng.unit() < 0.5 else v))
        frames3(rows, m3, d3, rng.log_uniform(-6.0, -1.0), ops, [rng.int(1, 3)])


def emit_4f(rows, m2, m3, mg):
    """M6 part 4f (ADR 0002 row 18): `total_energy` and `run`'s source at their range edges. The energies, 2D in both
    modes and 3D: fields of 10^±160 on dx of 10^∓150 (f² past DBL_MAX, or below the range, where dx² or dx³ brings the
    energy back), dx = 10⁻¹⁶⁰ (dx² subnormal, dx³ zero), media of 10⁻³⁰⁰ (½·m·k subnormal), an in-plane field of
    10⁻¹⁷⁰ alone (its square 0); the source on all three grids: f = DBL_MAX (2πf past it, the phase ∞·0 at t = 0),
    f = 10³⁰⁰ (the phase past DBL_MAX within a few steps), dx = 10⁻³⁰⁰ (dt and t subnormal). Field fills take each
    grid's field index."""
    for mode in (0, 1):
        for w, v, dx in ((0, 1e160, 1e-150), (1, 1e160, 1e-150), (0, 1e-160, 1e150), (2, 1e-160, 1e150),
                         (1, 1e-170, 1e150), (0, 1.0, 1e-160), (0, 3e-155, 1e150)):
            frames2(rows, m2, 8, 7, dx, mode, [fill2(w, v)], [])
            frames2(rows, m2, 8, 7, dx, mode, [fill2(w, v), source2d(3, 3, -v)], [1, 2])
        frames2(rows, m2, 8, 7, 1e-3, mode, [region2(1, (0, 0, 8, 7), 1e-300), region2(2, (0, 0, 8, 7), 1e-300),
                                             fill2(0, 1.0), fill2(1, 1.0)], [])
    for w, v, dx in ((2, 1e160, 1e-150), (3, 1e160, 1e-150), (0, 1e-160, 1e150), (5, 1e-160, 1e150),
                     (1, 1e-170, 1e150), (2, 1.0, 1e-160), (2, 1.0, 1e-110)):
        frames3(rows, m3, (6, 6, 6), dx, [fill2(w, v)], [])
        frames3(rows, m3, (6, 6, 6), dx, [fill2(w, v), source(2, 2, 2, -v)], [1])
    frames3(rows, m3, (6, 6, 6), 1e-3, [box(1, (0, 0, 0, 6, 6, 6), 1e-300), box(2, (0, 0, 0, 6, 6, 6), 1e-300),
                                        fill2(2, 1.0), fill2(3, 1.0)], [])
    mx = sys.float_info.max
    for fr, dx in ((mx, 1e-3), (1e300, 1e-3), (1e9, 1e-300), (3e307, 1e-3)):
        for mode in (0, 1):
            frames2(rows, m2, 8, 7, dx, mode, [], [1, 3], advance=lambda k, fr=fr: run2(k, (3, 3), fr, 1.0))
        frames3(rows, m3, (6, 6, 6), dx, [], [1, 3], advance=lambda k, fr=fr: run(k, (2, 2, 2), fr, 1.0))
        frames_gpu(rows, mg, 8, 7, dx, [], [1, 3], advance=lambda k, fr=fr: run2(k, (3, 3), fr, 1.0))


# ── M7: bridge.rs (part 1a) ─────────────────────────────────────────────────────────────────────────────────────────

def emit_1a(rows):
    """M7 part 1a: bridge.rs's eight functions (class A), their outputs validate.py's bridge_eval. Rust's tests' inputs;
    each argument through ±0, ±2^-1074, ±2^-1022, ±10^±300, ±1, ±DBL_MAX, ±∞ and NaN beside ordinary values; the
    ranges the functions are documented on (BK7's 300-2500 nm, the Lyman and Balmer lines, the light nuclei); the range
    edges a later Part reviews (squares and products past DBL_MAX or below 2^-1022 where the result is representable,
    the Sellmeier poles at √C₁, √C₂ and √C₃ µm); and seeded values over every binade of either sign."""
    inf, nan = math.inf, math.nan
    mx, tiny, dmin = sys.float_info.max, 5e-324, 2.2250738585072014e-308
    sp = (0.0, -0.0, tiny, -tiny, dmin, -dmin, 1e-300, -1e-300, 1.0, -1.0, 1e300, -1e300, mx, -mx, inf, -inf, nan)

    def put(name, *v):
        fn = "bridge::" + name
        ins = tuple(i(x) if isinstance(x, int) else f(x) for x in v)
        rows.emit(fn, "A", ins, V.bridge_eval(fn, list(v)))

    rng = Rng(0x7A_B1)

    def draw(lo=-320.0, hi=308.0):
        x = rng.log_uniform(lo, hi)
        return -x if rng.unit() < 0.5 else x

    # ushma: J²·ρ, then p/(ρ·c_p) with its sentinel.
    put("joule_heating_power_density", 1.0, 1e-7)
    put("joule_heating_power_density", 1e6, 1.68e-8)
    for a in sp:
        put("joule_heating_power_density", a, 1.68e-8)
        put("joule_heating_power_density", 1e6, a)
    for j, r in ((1e155, 1e-10), (1e200, 1e-300), (1e-160, 1e100), (1e-200, 1e300), (-1e160, -1e-20)):
        put("joule_heating_power_density", j, r)
    for _ in range(40):
        put("joule_heating_power_density", draw(), draw())
    p = V.bridge_eval("bridge::joule_heating_power_density", [1e6, 1.68e-8])[0][1]
    for v in ((1000.0, 1000.0, 4186.0), (p, 8960.0, 385.0), (1000.0, 0.0, 4186.0), (1000.0, nan, 4186.0)):
        put("em_energy_to_heating_rate", *v)
    for a in sp:
        put("em_energy_to_heating_rate", a, 8960.0, 385.0)
        put("em_energy_to_heating_rate", 1.68e4, a, 385.0)
        put("em_energy_to_heating_rate", 1.68e4, 8960.0, a)
    for v in ((1.0, 1e-200, 1e-200), (1e300, 1e200, 1e200), (1e-300, 1e150, 1e160), (1e308, 1e-160, 1e-160)):
        put("em_energy_to_heating_rate", *v)
    for _ in range(40):
        put("em_energy_to_heating_rate", draw(), draw(), draw())
    # prakash: BK7's Sellmeier fit over its range, beyond it, at its poles; the plane wave's intensity.
    for nm in range(300, 2501, 50):
        put("wavelength_to_refractive_index_bk7", float(nm))
    for nm in (550.0, 450.0, 650.0, 587.5618, 486.1327, 656.2725, 1064.0, 1550.0, 77.0, 141.0, 10200.0, 1e6):
        put("wavelength_to_refractive_index_bk7", nm)
    for a in sp:
        put("wavelength_to_refractive_index_bk7", a)
    for c in ("BK7_C1", "BK7_C2", "BK7_C3"):
        pole = math.sqrt(V.bridge_consts()[c]) * 1000.0
        for k in range(-3, 4):
            put("wavelength_to_refractive_index_bk7", pole + k * pole * 2.0 ** -50)
    for nm in (1e-150, 1e-160, 1e150, 1e160, 3e156):
        put("wavelength_to_refractive_index_bk7", nm)
    for _ in range(30):
        put("wavelength_to_refractive_index_bk7", draw(-3.0, 8.0))
    for a in sp + (100.0, 27.5, 1e154, 1.3e154, 1e-160):
        put("e_field_to_optical_intensity", a)
    for _ in range(30):
        put("e_field_to_optical_intensity", draw())
    # tanmatra: hc/E, then k·Z·e/r.
    for a in sp + (10.2, 1.0, 1.89, 13.6, 2.55, 1e-310, 1e-306):
        put("energy_level_to_wavelength_nm", a)
    for _ in range(30):
        put("energy_level_to_wavelength_nm", draw())
    for z in (0, 1, 2, 26, 79, 92, 118, (1 << 31) - 1, 1 << 31, (1 << 32) - 1):
        for r in (5.29e-11, 1e-15, 1.0, 0.0, -0.0, -1.0, nan, inf, tiny, dmin, mx):
            put("nuclear_charge_potential", z, r)
    for _ in range(20):
        put("nuclear_charge_potential", rng.int(0, (1 << 32) - 1), draw())
    # dravya: e·E, then (3/2)·λ_s·(B/B_sat)² clamped.
    put("e_field_to_piezo_stress", 1000.0, 15.0)
    for a in sp:
        put("e_field_to_piezo_stress", a, 15.0)
        put("e_field_to_piezo_stress", 1000.0, a)
    for _ in range(30):
        put("e_field_to_piezo_stress", draw(), draw())
    for b in (1.0, 2.15, 3.0, 10.0, -3.0, -1.0, 0.5):
        put("magnetic_flux_to_magnetostrictive_strain", b, -7e-6, 2.15)
    put("magnetic_flux_to_magnetostrictive_strain", 1.0, 35e-6, 0.0)
    for a in sp:
        put("magnetic_flux_to_magnetostrictive_strain", a, -7e-6, 2.15)
        put("magnetic_flux_to_magnetostrictive_strain", 1.0, a, 2.15)
        put("magnetic_flux_to_magnetostrictive_strain", 1.0, -7e-6, a)
    for v in ((1e-160, -7e-6, 1e-150), (1e-300, 1e-10, 1e300), (1e300, -7e-6, 1e-300), (1.0, 1e300, 1.0)):
        put("magnetic_flux_to_magnetostrictive_strain", *v)
    for _ in range(30):
        sat = abs(draw(-10.0, 10.0))
        put("magnetic_flux_to_magnetostrictive_strain", sat * (rng.unit() * 6.0 - 3.0), draw(-12.0, -2.0), sat)


# ── modules ─────────────────────────────────────────────────────────────────────────────────────────────────────────

def new_header(name):
    """The header of a module with no frozen rows (M7 on): the platform line validate.py's check_platform reads."""
    return [
        f"# bijli golden corpus, format 1: {name}",
        f"# oracle: bijli 1.2.0 via validate.py; libm: glibc {V.glibc_version()} {platform.machine()}. Written by "
        "scripts/oracle/emit.py; do not edit.",
        "# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note (tokens: scripts/oracle/README.md)",
        "# rows: 0",
    ]


def module(name, parts):
    """The module file's text: its header (the row count updated), its frozen rows, then the parts' rows."""
    path = root / "tests/golden" / f"{name}.tsv"
    lines = path.read_text().split("\n") if path.exists() else new_header(name)
    header = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if ln and not ln.startswith("#")]
    frozen = body[:FROZEN.get(name, 0)]
    keys = [(ln.split("\t")[0], ln.split("\t")[2]) for ln in frozen]
    rows = Rows(keys)
    for part in parts:
        part(rows)
    total = len(frozen) + len(rows.lines)
    header = [f"# rows: {total}" if h.startswith("# rows: ") else h for h in header]
    text = "\n".join(header + frozen + rows.lines) + "\n"
    return path, text, len(rows.lines)


def main() -> int:
    model2 = V.Fd2Scripts(root, V.GLIBC_P)
    model3 = V.Fd3Scripts(root, V.GLIBC_P)
    model_gpu = V.Gpu2Scripts(root, V.GLIBC_P)
    problems = 0
    parts = [lambda r: emit_2b(r, model3), lambda r: emit_3(r, model_gpu),
             lambda r: emit_4a(r, model2, model3, model_gpu), lambda r: emit_4b(r, model2, model3),
             lambda r: emit_4c(r, model2, model_gpu), lambda r: emit_4d(r, model3),
             lambda r: emit_4e(r, model2, model3), lambda r: emit_4f(r, model2, model3, model_gpu)]
    for name, parts in (("fdtd_2d3d", parts), ("bridge", [emit_1a])):
        path, text, n = module(name, parts)
        big = V.too_big(f"tests/golden/{name}.tsv", text)
        if big:
            print("  " + big[0])
            problems += 1
            continue
        if check:
            if path.read_text() != text:
                print(f"BAD  {name}.tsv differs from the emitted rows (rewrite with scripts/oracle/emit.py .)")
                problems += 1
            else:
                print(f"ok   {name}.tsv: {FROZEN.get(name, 0)} frozen rows, {n} emitted")
        else:
            path.write_text(text)
            print(f"wrote tests/golden/{name}.tsv ({FROZEN.get(name, 0)} frozen rows, {n} emitted)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
