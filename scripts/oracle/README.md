# scripts/oracle: what the port derives from the Rust oracle

`oracle` drives the frozen Rust bijli 1.2.0 in [`rust-old/`](../../rust-old/) and writes two kinds
of output:
- **The golden corpus, `tests/golden/<module>.tsv`** (port plan §6 step 5). It has one row per
  call, with every float as its IEEE-754 bit pattern. The Cyrius parity harness
  (`tests/parity.tcyr`, M0 step 7) replays every row against the port and checks it within its
  tolerance class.
- **The generated Cyrius sources**: `src/constants.cyr` (the public constants), plus data
  modules and private constants as milestones need them (step 6). No value the port reads is
  typed by hand.

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml              # write
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # verify; exit 1 on any difference
python3 -I -B scripts/oracle/validate.py .                                          # independent re-check
python3 -I -B scripts/oracle/mpfr_ref.py                                            # the MPFR helper's self-checks
```

The format, the comparison rules and the generator's structure are impetus's
(`/home/macro/Repos/impetus/scripts/oracle/README.md`), which are kimiya's and tanmatra's: every
port's harness reads the same files. This README records what is specific to bijli.

## Rules

- **`rust-old/` is never built in place or edited.** This crate depends on it by path, with
  `features = ["full"]` (every module but the `compute-wgpu` scaffold). `Cargo.lock` was seeded
  from `rust-old/Cargo.lock`: the 193 shared crates are the versions the oracle was tested with,
  Rust hisab 1.4.0 (the oracle's `Complex`) among them.
- **The output is deterministic.** Two runs give byte-identical files, and `--check` passes on
  any machine with the recorded glibc. Grids use only + − × ÷, `floor`, decimal parsing, bit
  construction and seeded xorshift (`src/grids.rs`); never a transcendental std method, never a
  timestamp or a hash-map order.
- **glibc, mostly.** bijli calls the std f64 methods. On x86_64 Linux those reach the host glibc
  for sin, cos, sincos, tan, asin, acos, atan, atan2, exp, pow, hypot and log10, and through the re-exported
  `Complex` for ln, sinh and cosh (`src/glibc.rs`; port plan §5 P3). Every TSV header records the glibc and the CPU class
  (`glibc 2.44 x86_64 fma`), and `--check` refuses a corpus from another one. **`cbrt` is the
  exception**: Rust std links compiler-builtins' own, correctly rounded `cbrt` (the `libm`
  crate's port of CORE-MATH), so its rows are MPFR's bits. `powi` with a run-time exponent is
  compiler-rt's `__powidf2`.
- **Every input is a run-time value** (`std::hint::black_box`), so no row is a compile-time fold:
  LLVM folds asin, acos, tan and cbrt of constants (P3).
- **The corpus is regenerated only from the oracle.** A wrong oracle value is fixed in Rust as
  1.2.x, then regenerated. It is never patched in a `.tsv`.
- **Modules arrive with their milestone.** Each one starts by adding its golden emitter
  (`src/m_<module>.rs`, the grid from the port plan §8) and, where it has data, its Cyrius
  emitter (`src/t_<module>.rs`), registered in `src/main.rs`. `validate.py` gets a check for the
  module at the same time; it reports a listed module with no check as a problem.

| Module | Milestone | Golden rows | Golden emitter | Generated Cyrius |
|---|---|---|---|---|
| `complex` | M1 | 30,767 (+ 8 corrected) | `src/m_complex.rs` | — |
| `constants` | M0 | 12 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`, with the completeness proof) |
| `error` | M1 | 1,635 | `src/m_error.rs` | — |
| `field` | M1 (parts 3a–3b) | 22,325 | `src/m_field.rs` | — |
| `fmtf` | M1 | 37,576 | `src/m_fmtf.rs` (impetus's grid) | — |
| `math` | M0 (+ sinh, cosh and 125 corrected, M1) | 36,522 | `src/m_math.rs` | — |
| (correctly rounded tables) | M1 | — | — | `src/crmath_data.cyr` (`crtables.py`, written and checked by `validate.py`) |
| (correctly rounded tables, kimiya's) | M1 | — | — | `src/crmath_k_data.cyr` (`crtables_k.py`, written and checked by `validate.py`) |

Beyond the corpus, `cr_verify.py` checks each correctly rounded primitive against MPFR on seeded argument sets,
natively and under qemu-aarch64 (`python3 -I -B scripts/oracle/cr_verify.py . --n 100000 --aarch64`, after
`cyrius distlib`). `tan_slow` calls tan's fixed-point path directly, so the path the rounding test rarely reaches is
checked on arguments of its own.

## Format 1

Five header lines, then one row per call:

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **Header**: the module name; `# oracle: <provenance>; libm: <platform>` (the `--check`
  platform test reads it); the regenerate command; the column legend; `# rows: N`.
- **Tokens** (space-separated; `-` for an empty list):
  - `f:<16 lowercase hex>`: an f64 as its bit pattern;
  - `i:<decimal>`: an integer, a boolean as 0 / 1;
  - `e:<Enum>.<Variant>`: a Rust enum variant by name;
  - `s:<percent-encoded UTF-8>`: a string (bytes outside printable ASCII, and `%`, as `%XX`);
  - `n:`: `None`;
  - `E:<Variant>`: `Err(BijliError::<Variant>)`; the error rows follow it with the variant's fields and its Display
    text (`s:`), which the port reads back from the last-error record (ADR 0002 row 4);
  - `v:<len>`: a sequence, its elements' tokens following.
- **A complex value is two `f:` tokens**, the real part first.
- The **note** column is a human-readable rendering for review; the harness ignores it.

## Classes and comparison (port plan §7)

| Class | Rule |
|---|---|
| A | bit-identical (any NaN matches any NaN), or equal integers / strings / names |
| B | ≤ 2 ulp, expected 0 |
| C | 0 ulp against the corrected oracle: the row's value, or its `<module>.cr.tsv` correction where glibc misrounded |
| D | relative ≤ 1e-12 (expected to stay empty) |
| T | an FDTD trajectory frame: each field component's digest equal, probe cells by their path's class (M5) |

`<module>.cr.tsv` (the corrected oracle, kimiya ADR 0006) holds, for a row whose glibc calls
misrounded, the oracle's sequence of operations evaluated with correctly rounded primitives
(`mpfr_ref.py`). `validate.py --write` produces it; the generator only checks it has a module.

## bijli specifics

- `const` rows hold f64 bit patterns: every public constant is an `f64`.
- `f64::sin_cos` has two outputs, sin then cos (glibc's `sincos`, which equals the separate calls).
- `f64::powi` rows take the exponent as an `i:` input: bijli's exponents are run-time values,
  which reach `__powidf2`.
