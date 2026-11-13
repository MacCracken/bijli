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
| `charge` | M1 (parts 3d, 5d-1 to 5d-3) | 12,769 | `src/m_charge.rs` (field.rs's grids) | — |
| `circuit` | M2 (part 1c) | 56,106 | `src/m_circuit.rs` | — |
| `complex` | M1 (+ 5c-2's edge grid) | 31,677 (+ 8 corrected) | `src/m_complex.rs` | — |
| `constants` | M0 | 12 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`, with the completeness proof) |
| `error` | M1 | 1,635 | `src/m_error.rs` | — |
| `field` | M1 (parts 3a–3c, 5d-1 to 5d-5) | 26,770 | `src/m_field.rs` | — |
| `fmtf` | M1 (+ Debug, part 3c; + fixed precision, M2 part 1b) | 70,740 | `src/m_fmtf.rs` (impetus's grid) | — |
| `json` | M1 (part 4), M2 (part 1a: `Material`) | 19,489 (+ 8,899 corrected) | `src/m_json.rs` (serde_json 1.0.154 on `Complex`, `FieldVector`, `PointCharge`, `Material`) | — |
| `material` | M2 (part 1a) | 51,814 | `src/m_material.rs` (field.rs's vectors) | — |
| `math` | M0 (+ sinh, cosh and 125 corrected, M1) | 36,522 | `src/m_math.rs` | — |
| `maxwell` | M1 (part 3d) | 20,349 | `src/m_maxwell.rs` | — |
| `relativity` | M2 (parts 1d-1, 1d-2: the solver over six trajectories) | 26,244 | `src/m_relativity.rs` (field.rs's vectors) | — |
| `wave` | M2 (part 1b) | 66,690 (+ 3 corrected) | `src/m_wave.rs` (field.rs's vectors, material.tsv's materials) | `src/wave_data.cyr`, the dipole constants of ADR 0002 row 22 (`validate.py`'s `render_wave_data`: Cin(2π) by its series, π by Machin's formula checked against MPFR, η₀ from BJ_MU_0 and BJ_EPSILON_0, in rational arithmetic) |
| (correctly rounded tables) | M1 | — | — | `src/crmath_data.cyr` (`crtables.py`, written and checked by `validate.py`) |
| (2/π for the angle past DBL_MAX, bijli's) | M1 (part 5d-3) | — | — | `src/crmath_ux_data.cyr` (`crtables.py`'s `render_ux`, from π to 2,800 bits checked against MPFR) |
| (correctly rounded tables, kimiya's) | M1 | — | — | `src/crmath_k_data.cyr` (`crtables_k.py`, written and checked by `validate.py`) |

Beyond the corpus, `cr_verify.py` checks each correctly rounded primitive (kimiya's expm1 among them since M2 part
2a-1) against MPFR on seeded argument sets,
natively and under qemu-aarch64 (`python3 -I -B scripts/oracle/cr_verify.py . --n 100000 --aarch64`, after
`cyrius distlib`). `tan_slow` calls tan's fixed-point path directly, so the path the rounding test rarely reaches is
checked on arguments of its own. Its exact families (`cr_verify.py . exact`) check src/crmath_kx.cyr's copy of kimiya's
correctly rounded monomials, square roots and unbounded-exponent operations against exact values (Python's Fraction and
an integer square root), not MPFR. Its `r2_below` family (an exact family, M1 part 5d-1) checks src/crmath_cmp.cyr's exact Σ (aᵢ − bᵢ)² < 10^-p
against Fractions, on points crowded around 10^-30 and 10^-60. Its gamma families (exact families, M2 part 2a-2) check
src/crmath_gamma.cyr's correctly rounded γ and γ² (`gamma`, `gamma2`), its exact path alone (`gamma_slow`,
`gamma2_slow`: the rounding test sends a value there only within about 2^-100 of a midpoint) and its validity test
(`gamma_ok`) against Fractions and an integer square root, over c = the speed of light, 1, 2^-30, 2^30 and four random
c between, with c outside that range refused. Its prod5 families (exact families, M2 part 2b-3) check
src/crmath_mono.cyr's product of five full-precision doubles (`prod5`, and `prod5_acc` its exact path alone) against
Fractions, on any binade, products steered to the subnormal boundary, exact ties and signed zeros. Its mono_wide
families (M2 part 2b-4) check src/crmath_mono.cyr's monomial of up to six factors a side (`mono_wide`, and
`mono_wide_acc` its exact path alone), within the 253-bit capacity, the Larmor power's shape among them. `cr_driver.cyr`
reads up to sixteen words a call, and `cr_verify.py` refuses a longer one. Its ux families (`cr_verify.py . ux`, M1 part 5c-2) check src/crmath_ux.cyr's
e^x, cosh, sinh, ln(s·2^E), (s·2^E)^y and sin and cos past DBL_MAX against `mpfr_ref.wide`: MPFR at 53 bits with its
widest exponent range, where nothing over- or underflows. Each (s, e) result is read in two calls, its significand
and its exponent.

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

`<module>.adr.tsv` (M1 part 5) holds the ADR 0002 mappings: for a golden row whose port result an ADR 0002 row
changes, the port's outputs as that row gives them, written by `validate.py --write` from exact values
(`exact.py`: Fraction and an integer square root; `mpfr_ref.wide` for the transcendental steps of ADR 0002 row 16),
never from the port. Its columns are fn, inputs (the key, byte for byte), outputs, the ADR 0002 row number, and a note;
the harness loads it beside `<module>.cr.tsv`, compares the port against it under the row's class, and reports each
row number's count. A golden row never has both. The generator only checks that it has a module; it has no platform
header, since exact values are the same everywhere.

`<module>.cr.tsv` (the corrected oracle, kimiya ADR 0006) holds, for a row whose glibc calls
misrounded, the oracle's sequence of operations evaluated with correctly rounded primitives
(`mpfr_ref.py`). `validate.py --write` produces it; the generator only checks it has a module. `json.cr.tsv` holds, for a reader row serde_json misread
(its default float parse is best-effort: ADR 0002 row 11), the document read with every number correctly rounded
(Python's `float`).

## bijli specifics

- `const` rows hold f64 bit patterns: every public constant is an `f64`.
- `f64::sin_cos` has two outputs, sin then cos (glibc's `sincos`, which equals the separate calls).
- `f64::powi` rows take the exponent as an `i:` input: bijli's exponents are run-time values,
  which reach `__powidf2`.
- **`json` rows** (`m_json.rs`). A writer row (`json::<Type>::to_string`, class A) gives serde_json's text as parsed
  tokens: each number outside strings an `f:` token (its correctly rounded value), each run between numbers an `s:`
  token, so the port's text matches it whatever notation each printer chose (D5). A reader row
  (`json::<Type>::from_str`, class C at 0 ulp) takes the document (`s:`) and `i:<flags>`, what serde_json reads
  without checking (1 a repeated key, 2 a lone surrogate escape, 4 nesting past 128, set only on a document it
  reads: ADR 0002 row 12), and gives the f64 in field order, or `E:Deserialize` with serde_json's message in the
  note. Each (function, inputs) appears once, since a correction keys on it.
