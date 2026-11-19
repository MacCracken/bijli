# Getting started with bijli

bijli is a Cyrius library: electromagnetism for AGNOS, from point charges and Maxwell's equations to waveguides,
Gaussian beams, RF networks and FDTD solvers in one, two and three dimensions. Every result is bit-identical on
x86_64 (ELF and PE) and aarch64, and every transcendental it calls is correctly rounded. Coming from the Rust crate
1.2? Read the [migration guide](migration-2.0.md).

## Use it from your project

bijli ships as a bundle, one Cyrius file. In your `cyrius.cyml`, take it as a dependency and list the stdlib leaves it
needs (its sidecar, `dist/bijli.deps`, names them; `cyrius deps` reads it):

```toml
[deps]
stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita", "tagged", "result", "fnptr",
          "freelist", "mmap", "atomic", "thread_local", "sakshi", "bayan"]

[deps.bijli]
git = "https://github.com/MacCracken/bijli.git"
tag = "2.0.1"
modules = ["dist/bijli.cyr"]
```

`cyrius deps` vendors it as `lib/bijli.cyr`. A whole program: a 2D TM grid, a block of glass, a 3 GHz source, and the
field energy after 200 steps.

```cyrius
include "lib/bijli.cyr"

# A 120 × 120 TM grid (1 mm cells), a glass block in its middle, a 3 GHz soft source, 200 steps.
fn main(): i64 {
    alloc_init();
    var sim: BjFdtd2d;
    var t, v = BjFdtd2d_new(120, 120, 0.001, BJ_MODE2D_TM, &sim);
    if (is_err_result(t) == 1) { return 1; }          # bj_last_error_display() says why
    var glass: BjMaterial = BjMaterial_dielectric(2.25);
    var tm, vm = sim.set_material(70, 40, 90, 80, &glass);
    if (is_err_result(tm) == 1) { return 2; }
    sim.run(200, 30, 60, 3000000000.0, 1.0);
    var u: f64 = sim.total_energy();
    BjFdtd2d_free(&sim);
    if (u > 0.0) { return 0; }
    return 3;
}
```

Every fallible function returns a value-form `Result` (`var t, v = f(..)`); on an error `v` is a `BJ_ERR_*` code and
`bj_last_error_display()` is the error's text. A struct that owns memory (the FDTD grids, an S-matrix, a field
slice) is given back with its `_free`; a constructor that fails leaves it empty, so a `_free` after the failure (a
`defer`, say) does nothing. Formatting writes into a text over a buffer you own:

```cyrius
var z: BjComplex = BjComplex_new(1.0, -2.0);
var w: BjComplex = z.exp();
var buf[64];
var t[24];
bj_text_init(&t, &buf, 64);
w.display(&t);                                       # Rust's Display of the Complex
if (bj_text_len(&t) == 0) { return 1; }
```

Logging is off until you ask: `bj_log_init()` reads `BIJLI_LOG` (`debug`, or `bijli=trace`) and sends bijli's events
through sakshi.

## The ai profile

The daimon/hoosh client (registering bijli as a daimon agent) is in a separate bundle, `dist/bijli-ai.cyr`, so the
default one never carries TLS. Take `modules = ["dist/bijli-ai.cyr"]` instead, add `"sandhi"` to your stdlib list and
include `lib/bijli-ai.cyr`. sandhi's stack is large: with bijli it passes the compiler's 128 KiB static-data advisory,
which the compiler reports as a warning on every build of a program that includes it.

## Build bijli itself

```sh
cyrius deps                              # resolve the stdlib (vendors lib/)
cyrius build src/main.cyr build/bijli    # compile
cyrius test                              # the unit tests and the golden parity harness
cyrius bench                             # the benchmarks
cyrius distlib --all                     # dist/bijli.cyr and dist/bijli-ai.cyr
```

## Layout

- `src/` — one file per Rust module (or part of one), each `private` with its API marked `public`; the bundle's
  order is `cyrius.cyml`'s `[lib]` list. `src/main.cyr` is the entry point the build compiles.
- `tests/*.tcyr` — every Rust 1.2.0 test ported by name, the Cyrius-only checks (`cyrius_*` groups) and
  `tests/parity.tcyr`, the harness that replays the golden corpus.
- `tests/golden/` — the corpus: frozen rows from Rust 1.2.0 and rows emitted by `scripts/oracle/emit.py`.
  `scripts/oracle/validate.py` is its reference model (Python); Rust bijli 1.2.0 itself is read from the `1.2.0` tag
  (ADR 0003).
- `docs/` — the ADRs (`docs/adr/`), the port's records (`docs/development/`) and these guides.

## Adding a feature

1. Read the Rust it ports (`git show 1.2.0:src/<module>.rs`) and the port plan's notes on it.
2. Write the Cyrius in `src/`, in Rust's order of operations; a fix to Rust's behaviour is an ADR 0002 row.
3. Add golden rows: transcribe the function in `scripts/oracle/validate.py`, emit rows with `scripts/oracle/emit.py`,
   and dispatch them in `tests/parity.tcyr`.
4. Port the Rust tests by name into `tests/<module>.tcyr`, and add a bench twin for each Rust bench.
5. Run the cleanliness check in `CLAUDE.md`, then the cross-target runs, one at a time.

See [`../adr/template.md`](../adr/template.md) when a non-trivial design choice deserves an ADR.
