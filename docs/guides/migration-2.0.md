# Migrating from bijli 1.2 (Rust) to 2.0 (Cyrius)

bijli 2.0.0 is the Cyrius port of the Rust crate 1.2.0 ([ADR 0001](../adr/0001-cyrius-port-with-rust-oracle.md)). It
computes the same physics with the same results, errors and branches, checked row by row against a golden corpus of
770,000 rows. Where the Rust crate was wrong, inconsistent or incomplete, 2.0 fixes it, and every such change is a row
of [ADR 0002](../adr/0002-port-api-divergences.md). This guide covers both: how Rust's idioms map to Cyrius's, and
each change a 1.2 caller will notice.

## Using the library

bijli ships as two bundles, built by `cyrius distlib --all`:

| Bundle | Holds | Extra stdlib |
|---|---|---|
| `dist/bijli.cyr` | every module but the ai client | none past the sidecar's (`dist/bijli.deps`) |
| `dist/bijli-ai.cyr` | the same, plus the daimon/hoosh client | sandhi, the stdlib's HTTP/TLS stack |

A consumer vendors one bundle into `lib/` with `cyrius deps`, includes it, and lists the stdlib leaves its sidecar
(`dist/bijli.deps` or `dist/bijli-ai.deps`) names: atomic, thread_local, syscalls, tagged, string, freelist, io,
sakshi, result, ganita, math, str, alloc, mmap, bayan and fnptr (and sandhi for the ai bundle). Rust's cargo features
have no counterpart: every module is in the bundle, and the ai client is the only one in a profile of its own, so the
default bundle never carries TLS. [Getting started](getting-started.md) has a whole program.

## Rust idioms in Cyrius

| Rust 1.2 | Cyrius 2.0 |
|---|---|
| `bijli::field::electric_field_point_charge(..)` | `bj_electric_field_point_charge(..)`: free functions are `bj_<name>`, one namespace |
| `FieldVector`, `PointCharge`, `Fdtd2d` | `BjFieldVector`, `BjPointCharge`, `BjFdtd2d` |
| `v.magnitude()`, `sim.step_once()` | the same: methods keep their Rust names |
| `FieldVector::new(..)`, `Fdtd2d::new(..)` | `BjFieldVector_new(..)`, `BjFdtd2d_new(..)` |
| `field::EPSILON_0` | `BJ_EPSILON_0` |
| `Mode2d::Tm`, `BijliError::InvalidParameter` | `BJ_MODE2D_TM`, `BJ_ERR_INVALID_PARAMETER`: append-only numeric values |
| `Result<f64>` | a value-form `Result`: `var t, v = f(..)`; `is_err_result(t)`; on Ok bind `v` to an f64 (`var x: f64 = v;`), on Err `v` is the `BJ_ERR_*` code |
| `Result<FieldVector>`, `Result<Vec<_>>` | the value goes to an `out` pointer the caller supplies; the Result carries Ok(0) or the code |
| `Option<f64>` | Ok(`Some`) with the value in `out`, or Ok(`None`) |
| `[f64; 3]` arguments | a pointer to three f64 (`var p: f64[3];`) |
| `[f64; 3]` struct fields | tuple fields (`pc.position.0`) |
| `&str` input (Touchstone text, a soorat component) | bytes and a length, checked to be UTF-8 |
| `Vec` inside a struct | a block the struct owns, with its length; give it back with `BjX_free` |
| `impl Display` / `Debug` | `x.display(t)` / `x.debug(t)` into a bounded text over your buffer (`bj_text_init(&t, &buf, cap)`, `bj_text_ptr`, `bj_text_len`) |
| `serde_json::to_string(&x)` | `BjX_to_json(&x, sb)` into a `str_builder` |
| `serde_json::from_str::<X>(s)` | `BjX_from_json(buf, len, &out)` |
| `logging::init()`, `BIJLI_LOG` | `bj_log_init()`, the same variable |
| closures (`trace_field_line`'s field, `solve_retarded_time`'s trajectory) | Cyrius closures or function pointers |
| `hisab::num::Complex` (re-exported) | `BjComplex`, by value, with hisab's whole surface |

A fallible call, its error and the error's text:

```cyrius
var t, v = bj_wave_speed(BJ_EPSILON_0, BJ_MU_0);
if (is_err_result(t) == 1) { return 1; }
var c: f64 = v;                                     # the Ok payload is the f64: bind it to a typed local
if (c < 299792457.0 || c > 299792459.0) { return 2; }
var tn, vn = bj_wave_speed(-1.0, BJ_MU_0);
if (is_ok(tn) == 1) { return 3; }
var code = bj_last_error();                         # the same code as vn
if (code != vn) { return 4; }
var why = bj_last_error_display();                  # Rust's Display text of the error, byte for byte
if (strlen(why) == 0) { return 5; }
```

Structs come back through out pointers, and arrays are pointers:

```cyrius
var origin: f64[3];
origin[0] = 0.0;
origin[1] = 0.0;
origin[2] = 0.0;
var p: BjPointCharge = BjPointCharge_proton(&origin);
var at: f64[3];
at[0] = 1.0;
at[1] = 0.0;
at[2] = 0.0;
var e: BjFieldVector;
var t, v = p.electric_field_at(&at, &e);
if (is_err_result(t) == 1) { return 1; }
var m: f64 = e.magnitude();                         # about 1.44e-9 V/m, a proton at 1 m
if (m <= 0.0) { return 2; }
var buf[128];
var tx[24];
bj_text_init(&tx, &buf, 128);
e.display(&tx);                                      # Rust's Display: "(x, y, z)"
if (bj_text_len(&tx) == 0 || bj_text_overflowed(&tx) == 1) { return 3; }
```

An `Option` is `Some` or `None` beside the value:

```cyrius
var theta: f64[1];
var t, v = bj_snell_refraction_angle(1.0, 1.5, 0.5, &theta);
if (is_err_result(t) == 1 || v != Some) { return 1; }
var tir, none = bj_snell_refraction_angle(1.5, 1.0, 1.2, &theta);
if (none != None) { return 2; }                     # total internal reflection: Rust's None
```

JSON reads and writes as serde_json did:

```cyrius
var z: BjComplex = BjComplex_new(3.0, -4.0);
var sb = str_builder_new();
var tw, vw = BjComplex_to_json(&z, sb);
var doc = str_builder_build(sb);                     # {"re":3.0,"im":-4.0}
var back: BjComplex;
var tr, vr = BjComplex_from_json(str_data(doc), str_len(doc), &back);
if (is_err_result(tw) == 1 || is_err_result(tr) == 1) { return 1; }
var im: f64 = back.im;
if (im != -4.0) { return 2; }
```

## What changed (ADR 0002)

Each item names its ADR 0002 row; the row has the details and how the golden corpus maps it.

### Removed and reshaped

- **The 1.2 deprecations are gone** (row 1): `rf::input_impedance_lossless`, `rf::write_touchstone` and
  `Material::reflectance_at`. Use `bj_input_impedance_lossless_complex`, `bj_try_write_touchstone` and
  `try_reflectance_at`.
- **Shapes are fixed at construction** (rows 2 and 36). The FDTD grids' arrays, `BjSMatrix` and `BjTouchstoneData`
  can't be resized, so Rust's per-step layout checks and their errors are gone. Read a grid's fields through its
  struct (`sim.e_field`, `sim.num_cells`) or its accessors; a matrix's entries through `get`, `set` and `s_at`.
- **The sanitizing constructors are fallible** (row 3): `BjTfsfSource_try_new` and `BjCpml_new` refuse what Rust's
  `new` clamped with a warning, and `Fdtd2d::enable_cpml` is `try_enable_cpml` (2.0.1: before it, `BjCpml_new` and
  `enable_cpml` still sanitized).
- **`GpuFdtd2d<B>`'s methods are generic functions** (row 43): `BjGpuFdtd2d_step_once(&sim)` for `sim.step_once()`.
  `BjComputeBackend` is a Cyrius trait with an explicit `free_buffer`, and `BjGpuFdtd2d_free` gives every buffer back.

### Errors

- **Errors are codes** (row 4): `BJ_ERR_*` in a value-form `Result`, the variant's fields and Rust's Display text in
  the calling thread's last-error record (`bj_last_error`, `bj_last_error_display`, `bj_last_error_detail`,
  `bj_last_error_value`, …). Nothing is logged.
- **Three codes are new**: `BJ_ERR_TRANSPORT` (row 7) for the ai client's network failures, `BJ_ERR_ALLOCATION`
  (row 9) where Rust aborted on a failed allocation, and `BJ_ERR_DESERIALIZE` (row 10) for a refused JSON document.

### One time convention

- **e^{+jωt}, IEEE handedness, in every module** (row 5). Rust 1.2 used e^{−iωt} in `polarization`, `beam`, the
  field's Green's function and the FDTD media, and e^{+jωt} elsewhere. In those four, 2.0's complex amplitudes are
  1.2's conjugates and its phases their negations: an absorbing index is n − jκ, a wave toward +z is A·e^{−jkz}, and
  right-hand circular polarization is (1, −j)/√2. The convention is stated once, in `src/complex.cyr`.

### Exact and correctly rounded results

The physics is the same; the last bits may differ, and the range edges no longer over- or underflow.

- **Bit-identical on every target** (row 6): x86_64 ELF, x86_64 PE and aarch64 give the same bits, with correctly
  rounded transcendentals ([ADR 0001](../adr/0001-cyrius-port-with-rust-oracle.md)). 1.2 matched only on one platform.
- **Round once** (row 13): every function that is one product, quotient or root of its inputs and constants is
  correctly rounded, where Rust rounded two to four times (one ulp off on about a quarter of inputs).
- **No spurious overflow or underflow** (rows 14–18): at a range edge, where an intermediate left binary64's range
  but the result is representable, the chain runs with an unbounded exponent. `abs(1e160 + j1e160)` is 1.41e160, not
  ∞; `powf` has no cutoff below |z|² = 1e-300 (row 15). The Singularity tests compare the exact r² against 10⁻³⁰
  (row 17).
- **Cancellation fixed**: the half-wave dipole near its axis (row 28), the RC and RL rise (row 29, as expm1), the
  Lorentz factor near c (row 30), the waveplates (row 33), the resonators near the planar edge (row 34), the beam's
  ABCD determinant (row 35), the complex stub match (row 39), the FDTD media coefficients and subpixel means (rows 40
  and 41), and the graded grid's last cell, exactly d_max (row 42).
- **Exact constants**: the dipoles use bijli's own η₀ = √(μ₀/ε₀), not 120π (row 22).
- **max and min order −0 < +0** (row 47): IEEE 754-2019's rule in either order, where Rust's x86-64 code returned its
  second operand on a tie.

### Validation

- **NaN no longer slips through** a `<= 0` test: the Curie laws (row 19), `classify_magnetic`, now a Result (row 21),
  `effective_aperture` (row 23), `resistance_from_geometry` and `parallel_plate_capacitance` (row 26), and
  `lienard_wiechert_e` (row 27) refuse it with Rust's own messages. `Material::wave_speed` validates as `impedance`
  does (row 20).
- **Fixed inconsistencies**: `schlick_reflectance` takes R₀ without the sum guard (row 24); the material reflection
  coefficient is Γ on the relative impedances (row 25); an infinite Z₀ no longer equals a finite one (row 37); a
  superluminal trajectory reports the bracket error (row 31); `laguerre_gaussian_intensity` takes |l| = 2³¹ for
  `i32::MIN` (row 32).

### Text, JSON and logging

- **Text is bytes** (rows 38 and 44): `bj_parse_touchstone` and the soorat field-slice constructors take an address
  and a length, and refuse bytes that are not UTF-8 (the lenient soorat forms fall back as for an unknown name).
- **JSON numbers are correctly rounded** (row 11), where serde_json's default parse misread 6,342 of the corpus's
  documents; an f32 field is rounded once from the decimal (row 45), and a number is refused only when its rounded
  value overflows.
- **The whole JSON document is checked** (row 12): a repeated key, a lone surrogate or nesting past 128 levels anywhere
  refuses it, where serde_json skipped unknown values unchecked. A usize count above 2⁶³ − 1 is refused (row 36).
- **Logging is bijli's own** (row 8): `bj_log_init` sets bijli's gate from `BIJLI_LOG` (a level or `bijli=<level>`)
  and only raises sakshi's process-wide level, never lowering another library's.

### The ai client

- **Its own profile**: the daimon/hoosh client is `dist/bijli-ai.cyr` only (Rust's `ai` feature).
- **daimon's route** (row 46): `register_agent` posts to `{endpoint}/v1/agents` and returns daimon's integer `id` as
  text. Rust's `/v1/agents/register` and its `agent_id` field exist on no daimon.
- **Its failures are Transport** (row 7), with Rust's `daimon:` reasons.

## `BjComplex` for hisab's `Complex`

1.2 re-exported hisab 1.4's `Complex`; 2.0 has its own `BjComplex`, a 16-byte value with fields `re` and `im` and
hisab's surface in hisab's order of operations: `BjComplex_new`, `BjComplex_from_polar`, the arithmetic
(`BjComplex_add`, `_sub`, `_mul`, `_div`, `_neg`), and the methods `conj`, `abs`, `norm`, `norm_sq`, `arg`, `inv`,
`exp`, `ln`, `sqrt`, `powf`, `sin`, `cos` and `is_zero`. Its JSON is hisab's `{"re":…,"im":…}`. Only the range edges
differ (rows 14–16), and its transcendentals are correctly rounded.

## Where to look

- [ADR 0002](../adr/0002-port-api-divergences.md): every change, with its golden comparison.
- [`docs/development/port-parity.md`](../development/port-parity.md): how each module was ported and measured.
- `tests/*.tcyr`: every Rust test ported by name, a worked example of each function's Cyrius form.
