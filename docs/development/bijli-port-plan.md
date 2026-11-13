# bijli Rust → Cyrius port plan

> **Plan of record** for porting bijli to Cyrius. Rust **1.2.0** (tag `1.2.0`, commit `c95fd30`) is
> the final Rust release and the **reference implementation** (the oracle); the Cyrius port ships as **2.0.0**.
> The method is tanmatra's ([its ADR 0006](../../../tanmatra/docs/adr/0006-cyrius-port-with-rust-oracle.md),
> [its port plan](../../../tanmatra/docs/development/tanmatra-port-plan.md)), as executed by impetus
> ([its port plan](../../../impetus/docs/development/impetus-port-plan.md), finished 2026-11-10) and kimiya
> ([its port plan](../../../kimiya/docs/development/kimiya-port-plan.md), in progress). Where this plan
> says "as impetus" or "as kimiya", follow that document. §2 is the review of the Cyrius
> language that every scope decision rests on; §3 lists what is different from the siblings.
>
> ▶ **Prep first (§5), then `cyrius port`.** Start on Cyrius 6.7.18. The oracle is tagged (§5 P1), and
> nothing in `src/` changes after the tag.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | bijli **1.2.0**, tag `1.2.0` at `c95fd30` (2026-11-11), which includes the 3D CPML that the tagged CHANGELOG still lists under `[Unreleased]`. Never published: crates.io stops at 1.1.0. 8,589 code lines in `src/` (no tests, no comments, no blank lines), 19 files, 430 `pub fn`s, 827 `#[test]`s, 80 criterion benches, 1 example | `Cargo.toml`, `VERSION`, a line count of `src/` |
| Rust health, 2026-11-11 on rustc 1.99.1 | ✅ `fmt --check`, clippy `-D warnings` (all features, all targets), rustdoc `-D warnings`, `cargo audit`, `cargo deny`. Tests: 830 unit, 26 integration and 2 doc tests pass in 6.2 s | local run |
| Correctness baseline | The 1.2.0 audit of every module (CHANGELOG 1.2.0, *Breaking (behavioral)* and *Fixed*), plus the 3D CPML measured at −80 dB | `CHANGELOG.md` |
| Cyrius toolchain | **6.7.18** (installed; still the newest) | `~/.cyrius/current` |
| hisab | Rust **1.4.0** gives bijli exactly one thing, `hisab::num::Complex`. Cyrius hisab **3.4.1** isn't needed in production (§0.2) | `grep hisab src/`, `Cargo.lock` |
| ganita | **1.2.15**, bundled in the 6.7.18 stdlib | `~/.cyrius/versions/6.7.18/lib/ganita.cyr` |
| Oracle's math | Rust **std** calls the host **glibc 2.44** for `sin` (27 calls), `cos` (21), `exp` (11), `sin_cos` (8), `atan2` (5), `log10` (5), `powf` (4), `hypot` (4), `atan` (4), `asin` (4), `acos` (3), `tan` (2) and `cbrt` (2). Complex values use Rust hisab 1.4.0's naive formulas (§9.2) | `grep` of `src/` outside tests |
| Cyrius consumers | **prakash 2.15.2** has primitive-value bijli bridges of its own and its own Gaussian-beam module. **tanmatra 2.0.x** has deferred its bijli bridge until bijli ports. Still Rust: soorat 1.0.0 (feature `em`) and ushma 1.4.0. No local repo: kiran, joshua, dravya | `~/Repos` |
| Target | bijli **2.0.0**, Cyrius | this document |

## 0. Decisions

### 0.1 Carried from the earlier ports (no new decision needed)

- **The oracle is frozen.** Rust 1.2.0 moves to `rust-old/` and is never edited. A generator
  crate in `scripts/oracle/` drives it by path, as tanmatra and impetus do. A wrong oracle value
  is fixed in Rust as 1.2.x and the corpus is regenerated.
- **Correct rounding is the accuracy target**, from the first milestone. This follows
  [kimiya ADR 0006](../../../kimiya/docs/adr/0006-correct-rounding-is-the-accuracy-target.md) and
  [impetus ADR 0011](../../../impetus/docs/adr/0011-correctly-rounded-exp-sin-cos.md); the user's ruling
  was "we want to be better than glibc". Every transcendental the port calls is correctly
  rounded and gives the same bits on every target. A class-C golden row compares against the
  **corrected oracle**: the Rust sequence of operations, transcribed in `validate.py` and
  evaluated with MPFR-correct primitives. Where glibc misrounded a call, the port must match the
  corrected value, never glibc's.
- **Inexact formulas and range edges are fixed in the port, within a set scope.** This follows
  kimiya's ADR 0002 rows #12 onward and its handoff memory `accuracy-over-glibc-parity.md`. The
  port fixes three kinds of case:
  - a simple exact constant that Rust misrepresents;
  - a spurious overflow or underflow at a range edge, where the true result is representable.
    The fix is a power-of-two prescale applied only when the edge is hit, so in-range rows keep
    Rust's bits;
  - a physically exact value that comes out wrong.

  Ordinary chains of IEEE operations stay as Rust wrote them. Each fix is an ADR 0002 row, with
  its oracle rows mapped. Anything debatable goes to the user.

  **The user's decisions for bijli (2026-11-12, M1 part 5), each the recommended option, as kimiya's
  were:**
  - **range edges**: where an intermediate leaves the normal range and the result is representable,
    Rust's chain is evaluated with an unbounded exponent and its last step rounded once (kimiya ADR
    0002 #20); in-range rows keep Rust's bits;
  - **round once in range** for every formula that is one product, quotient or root of inputs and
    constants (kimiya's overflow-audit answer): correctly rounded everywhere, not only at the edges.
    Formulas with sums inside get the edge fix only;
  - **exact thresholds**: the Singularity tests compare the exact r² (or r, or distance) against
    10⁻³⁰, not the rounded value against the double 1e-30;
  - **hisab's `powf` cutoff is dropped**: `|z|² < 1e-300` no longer returns 0 whatever the exponent,
    so z⁰ = 1 and z⁻¹ = 1/z for every z ≠ 0; z = 0 follows pow's rules.

  The machinery is kimiya's, copied by `scripts/port/kimiya_closure.py` (src/crmath_kx.cyr: its
  correctly rounded monomials, √(x/y), w·√((x/y)·2^e), the multi-word integers and binary64 with
  an unbounded exponent). M1's Part 5 runs as 5a (the machinery, verified), 5b (round once: the
  single-product functions of maxwell, charge and field), 5c-1 and 5c-2 (complex: the algebraic edges, then
  the transcendental ones) and 5d (field's chains and the thresholds), each its own commit. 5d, about 25
  functions with sums inside, runs in bites:
  - 5d-1, the exact Singularity threshold (ADR 0002 row 17);
  - 5d-2, BjFieldVector's magnitude, magnitude_sq, dot and cross, the Lorentz force, and the plane and dielectric
    images;
  - 5d-3, the point charges (E, V, B, superposition) and `coulomb_force`, and the Green's functions;
  - 5d-4, the dipoles, the quadrupole, the ring and disk axes, and the sphere image;
  - 5d-5, the stress tensor and its force, the Poynting flux, the multipole moments, and the field line's step.

  All five bites are done (2026-11-12; port-parity §20–§24).

  kimiya ran this as a 12-batch, library-wide audit after M4. Here it is **one Part at the end of
  each milestone**. bijli's obvious candidates:
  - `re² + im²` in complex `abs`/`norm_sq`/division;
  - `x² + y² + z²` in field magnitudes;
  - `1/r²` near the r² < 1e-30 m² cutoff (a threshold, so it goes to the user);
  - the products of physical constants in skin depth and Poynting.
- **Performance is measured, not gated.** Every Rust bench gets a Cyrius twin. Both columns are
  recorded per milestone in `docs/benchmarks-rust-v-cyrius.md`. impetus's median ratio was
  10× in 2D and 13× in 3D. A slow row is a finding to file with its numbers, not a reason to
  change the port.
- **No Rust 2.0.** The roadmap's "2.0 cleanup (breaking)" items land in the Cyrius API only,
  accepted by ADR (§4.3). New physics stays out of the port.
- **The port removes Rust; it doesn't carry Rust's faults.** Parity with the oracle is
  functional: the same physics, results, errors and branches. It is not a reason to keep an
  inconsistency or a gap the Rust crate has. The time-convention split (§4.4) and the missing
  transport error (§4.3 row 7) are fixed in the port, and each fix is an ADR 0002 row with its
  oracle rows mapped exactly. "As Rust" is never a justification on its own; the oracle's
  operation order is kept only where it is what makes a golden row comparable.
- **Bit-identical across targets** for everything bijli computes: x86_64 ELF, x86_64 PE (wine)
  and aarch64 (qemu), as impetus.

### 0.2 New for bijli (proposed here; each is confirmed in ADR 0001 or 0002)

- **No production hisab dependency.** bijli uses only `hisab::num::Complex` from hisab. Cyrius
  hisab's `cx_*` doesn't reproduce it: it uses Smith division, a hypot-rescued `abs` and a C99
  `sqrt`, and it allocates 16 B on the never-freed bump heap for every result. bijli therefore
  owns a **by-value complex type**: a transcription of Rust hisab 1.4.0's `Complex`, operation
  for operation (§9.2). hisab may join at `scope = "test"`, as prakash and impetus have it, if a
  cross-check needs it. This is impetus's [ADR 0012](../../../impetus/docs/adr/0012-impetus-owned-3d-math-stopgap.md)
  and its second ADR 0007 amendment, applied from day one.
- **Reuse the correctly rounded primitives; don't write a third copy.** impetus has
  `src/crmath.cyr` (exp, sin, cos, sin_cos, atan2, hypot, pow), and kimiya has `src/mathx.cyr`
  (exp, expm1, ln, log10, pow, sin, cos, sincos, atan2, acos). bijli copies the impetus file,
  renamed to `_bj_cr_*`, regenerates its tables with impetus's `crtables.py`, and adds kimiya's
  acos and log10.
  - atan comes free: `atan(y)` is `atan2(y, 1)`, and both are the correctly rounded value of the
    same real number.
  - cbrt is already correctly rounded in ganita.
  - That leaves **asin and tan** as the only new primitives. sinh appears only in tests. *(M1 part 2c/2d: the
    `Complex` bijli re-exports is bijli's API too, and its `sin`, `cos` and `ln` reach sinh, cosh and ln, so
    **sinh and cosh** are new as well, and ln is kimiya's.)*

  File the shared module upstream (§10 G1), so a fourth port can depend on it instead of copying.
- **The FDTD solvers are compared over whole trajectories** (§3 item 2). Their updates are
  linear and use only `+ − × ÷`, so a whole run is class B.

### 0.3 Decisions for the maintainer before §5 (all answered 2026-11-11)

| # | Question | Recommendation | Alternative |
|---|---|---|---|
| D1 | Which Rust version is the oracle? | ✅ **Settled 2026-11-11: the `1.2.0` tag at `c95fd30`.** It includes the 3D CPML, so no 1.3.0 is cut | — |
| D2 | ~~`bridge.rs`~~ | **Settled by §2: ported, all 8 (M7).** They are dependency-free functions, and two partners (prakash, tanmatra) are already in Cyrius | — |
| D3 | Rust mixes two time conventions (§4.4). Which **one** does the port use everywhere? | ✅ **Accepted 2026-11-11: e^{+jωt} (engineering) with IEEE handedness, in every module, in 2.0.0** (§4.4). Keeping Rust's split was wrong: the port removes Rust's faults instead of carrying them | e^{−iωt} everywhere: rejected (§4.4) |
| D4 | ~~`ai`~~ | **Settled by §2: ported, in a `[lib.ai]` profile (M7)**, on prakash's model (sandhi + bayan) | — |
| D5 | JSON for the 36 serde types (§2 rows *serde*) | ✅ **Accepted 2026-11-11:** **`#derive(Serialize)` writers for every type, and hand-written readers on bayan's tree** (impetus's `json.cyr` pattern) until the three open reader defects are fixed upstream. Then switch the readers to the derive's scanner. The writer emits serde's shape with the same digits; JSON golden rows compare parsed bits, because the exponent notation differs (`1e+16` against `10000000000000000.0`). `BjComplex` serializes as `{"re":..,"im":..}`, as Rust hisab's does. **Amended 2026-11-12 (M1 part 4):** the writers are hand-written too, calling the derive's own printer (`bayan_f64_to_json`), so the text is the derive's byte for byte: on 6.7.18 a derive in a second `private` file fails to compile (its shared helper is private to the first; probe p16b), and in a library file the derive exports its zero-filling, bump-allocating decoders as public API. The readers apply serde's semantics with ADR 0002 rows 10–12 | (a) derive for both directions now, accepting the reader defects on untrusted input; or (b) impetus ADR 0010 in full: hand-written, text-identical serde codecs in both directions |
---

## 1. Porting models and references

- **impetus 2.0.0** (`/home/macro/Repos/impetus`): the closest model and the newest finished port.
  - The `scripts/oracle/` generator crate, the corrected oracle (`mpfr_ref.py`) and
    `validate.py`'s transcription style;
  - `src/crmath.cyr` + `src/crmath_data.cyr` (correctly rounded primitives), `src/mathx.cyr` (the
    Rust-semantics helpers), `src/fmtf.cyr` (shortest round-trip digits) and `src/json.cyr`
    (serde-compatible JSON);
  - `scripts/rust_test_ledger.py`, `scripts/surface_parity.py`, `scripts/bench_record.py` and
    `scripts/doc_examples.py`;
  - the M0 probe set and the CLAUDE.md translation of the Rust principles;
  - the per-milestone *Parts*, each closing with green gates.
- **tanmatra 2.0.x**: the template that impetus copied: TSV format 1, `tests/parity.tcyr`,
  `tests/parity_core.cyr`, `tests/ulp.cyr`, `docs/development/port-parity.md` and the class-C rule.
- **kimiya** (in progress): ADRs [0005](../../../kimiya/docs/adr/0005-correctly-rounded-primitives-for-ill-conditioned-functions.md)
  and [0006](../../../kimiya/docs/adr/0006-correct-rounding-is-the-accuracy-target.md), `_km_cr_acos` and
  `_km_cr_log10`, and the `.cr.tsv` corrected-row files; also `handoff/HANDOFF.md`.
- **prakash 2.15.2**: the `[lib.ai]` profile and `src/ai.cyr` (bayan and sandhi). Its
  `src/wave_beam.cyr` (Gaussian beams, e^{+iωt}) and `src/bridge.cyr` (bijli bridges) are the
  Cyrius consumer side of bijli's beam and bridge modules.
- **Language and toolchain** (`/home/macro/Repos/cyrius`):
  - `docs/guides/cyrius-guide.md`, `docs/stdlib-reference.md` and
    `docs/development/ecosystem-migration.md` (the 6.7.x entries);
  - `cyrius help manifest`;
  - the open issues bijli will hit: `2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`,
    `2026-11-08-small-struct-result-payload-read-as-address.md`, `2026-11-09-field-store-of-struct-call-copies-a-field-argument.md`,
    `2026-11-09-no-scoped-allocator-for-alloc.md`, `2026-11-06-f64v-builtins-no-scalar-tail.md`,
    `2026-11-02-memcpy-memset-byte-loops-dominate-frame-copies.md` and `2026-11-09-exp-polyfill-double-rounds-subnormal-results.md`.

---

## 2. Language review: Cyrius 6.7.18, checked before scoping

Every scope and mapping decision below rests on this section. Each row was read in
`cyrius/docs/guides/cyrius-guide.md` or `docs/stdlib-reference.md`, and the ✅/⚠ rows were run as
programs on the installed 6.7.18 compiler (`scripts/probes/`, 2026-11-11). Numbers
were checked bit for bit against Python's IEEE arithmetic. **The sibling ports' limitation lists
were not taken on trust**: two of their claims are refined here (marked "corrects").

| Rust need in bijli | Cyrius 6.7.18 | Evidence | Port consequence |
|---|---|---|---|
| `Complex`, `FieldVector`, Jones/Stokes types with `+ − × ÷` | struct values (16 B in two registers, 24 B by hidden pointer) with operator fns `T_add` / `T_sub` / `T_mul` / `T_div` (exactly two parameters; struct or `f64` right operand); typed `f64` fields | guide *Structs*, *Returning a struct by value*; **probe p1 ✅** `a + b`, `a * b`, `a / b`, `(v + w) * 2.0`, neg: all bit-exact | `BjComplex`, `BjFieldVector` by value with operators, as in Rust |
| Arrays of complex values (Mie `a`/`b`, `SMatrix`) | `var a: C[N]`, `*C` pointers stepping `sizeof(C)` | guide *Arrays*, *Pointers*; **p1 ✅** | stride buffers of `BjComplex` |
| `&[f64]` / `&mut [f64]` grid rows | typed `*f64` with `p[i]` (float semantics, 6.7.8); `slice<f64>` (bounds-checked, fn-local) | guide *Subscripts*, *Slices*; **p2 ✅** | the FDTD kernels index `*f64` directly |
| f64 compares (`<`, `==`, NaN) | on typed `f64` values: IEEE (NaN ≠ NaN, −0 == +0). On untyped words and builtin results: integer compares; use `f64_lt` etc. | guide *Number Literals*; **p1 ✅** | type every float binding |
| `trait ComputeBackend` + `GpuFdtd2d<B: ComputeBackend>` | **traits with checked `impl`s and defaults (6.7.0); generic structs with trait bounds (6.7.1); `dyn Tr` vtables (6.7.8)**; at most two type parameters; no associated types | guide *Traits*, *Trait objects*, *Generic structs*, *Trait bounds*; **p2 ✅** bounded `Grid<B: Backend>` dispatching `g.be.update_h(..)`, `dyn` dispatch 79 = 7·10 + 9 | `trait BjComputeBackend`, `struct BjGpuFdtd2d<B: BjComputeBackend>`. The associated `Buffer` type becomes an opaque handle word |
| `ComplexExt` (a trait on hisab's type) | bijli owns `BjComplex`, so these are inherent `impl BjComplex` methods | guide *Traits* | `z.norm_sq()`, `z.norm()`, `z.arg()` keep their Rust names |
| `impl Fn` / `F: Fn` arguments (`trace_field_line`, `solve_retarded_time`, `lienard_wiechert_fields`, NF2FF) | **closures `\|a, b\| { .. }` with by-value capture** (6.3.8). The value is a fn pointer or an opaque handle, called with `callptr`; ≤ 5 parameters on Linux/macOS (3 on Windows); parameters are untyped words; a callee held in a variable returns ONE word | guide *Closures*, *Function Pointers*; **p3 ✅** a captured point-charge field traced 1,000 steps, bit-exact (`0x4025ffffffffff95`); the environment allocates 40 B once, the calls 0 B | a field callback is `\|pos, out\| -> status` (pointers in, a `*BjFieldVector` out, a status word back). No separate context word is needed: capture carries it |
| — | ⚠ **a captured `f64` (a local or a parameter) is an untyped word inside the closure**: as a left operand it does integer arithmetic (`a * 2.0` gave 0) | **p6 ⚠** (not filed upstream yet, §10) | rebind every captured float as `var a2: f64 = a;` inside the body, or capture a struct with `f64` fields (**p3b ✅**) |
| `Result<f64>` / `?` | value-form `Result` (tag, one-word payload), zero allocation, `?` propagates | guide *Sum Types*, *`?` Propagation*; **p6 ✅** `Ok(a / b)?` then `* 2.0` | f64 results as `Result` directly |
| `Result<Complex>` / `Result<FieldVector>` | one-word payload; a struct over 8 B can't be one | guide *Sum Types* | code + out-pointer, or a declared `(f64, f64)` multi-return |
| `(f64, f64)` returns (stub match, CPML coefficients, complex parts) | declared multi-return `fn f(): (f64, f64)`, bindings typed by the declaration; arity 2 or 3 | guide *Multi-Return*; **p6 ✅** | tuples port as multi-returns; 4+ values go through an out-struct |
| `BijliError` with payloads (`DimensionMismatch { expected, got }`, `InvalidPermittivity { value }`, `String` reasons) | a `: stack` enum carries ONE word per variant; a boxed enum allocates on the bump heap at every error | guide *Sum Types* | `BJ_ERR_*` codes in `Result`. The payload goes to a per-thread last-error record (`thread_local.cyr`) read by `bj_last_error_*`, and the reason to sakshi. Nothing is lost and nothing allocates |
| `Vec<f64>` grids, freed on drop | **`freelist.cyr`: individual `fl_free`; blocks over 4,096 B are mmap'd and munmapped on free**; arenas; the `Allocator` vtable; `defer` | stdlib ref *freelist*, *alloc*; guide *Defer*; **p4 ✅** 1,000 cycles of 6 × 216 KB grids: 0 failed frees, 0 B global-heap growth, flat RSS | grid storage is `fl_calloc` / `fl_free` behind `Bj*_free`. No custom allocator. `defer { BjFdtd3d_free(g); }` in tests |
| serde derives (36 types) | **`#derive(Serialize)`**: `T_to_json`, `T_from_json_str` (scanner), `T_from_json(pairs)`; f64/f32/int/bool/Str/cstring/`Vec<f64>`/`Vec<derived>`/nested derived fields. ⚠ **in two `private` files the second derive fails to compile** (`'_cy_json_skipval' is private to its file`), and a `public` struct's derive exports its decoders | guide *Derive Serialize*; **p5 ✅/⚠** see the next three rows; **p16b ⚠** | D5 (amended: hand-written writers on bayan's printer) |
| — | ✅ writer: serde's shape (`{"re":0.1,"im":1e-12}`, `1.0`, `-0.0`, nested objects, arrays) with the same shortest round-trip digits | **p5** against serde_json 1.x | — |
| — | ⚠ number **notation** differs for some exponents: bayan `10000000000000000.0` / `0.000001` where serde_json writes `1e+16` / `1e-6` (each parses the other's) | **p5b** | JSON golden rows compare parsed bits, not text, unless a serde-identical printer is written |
| — | ✅ the scanner decoder reads `Vec<Complex>` back; ⚠ the pairs decoder `T_from_json(bayan_json_parse(..))` returns an empty `Vec<struct>` (index out of bounds); three reader defects are open upstream (2026-11-07: an over-read on an unterminated key, a hang on a stray `}` in `Vec<struct>`, a memset of a refused allocation) | **p5**; cyrius issues | D5 |
| — | bayan's value tree against serde_json: ⚠ `-0` is the integer 0 (an f64 field loses the sign), an integer outside i64 refuses the document, `1e400` is +∞ with no error, a raw control character or bytes that are not UTF-8 in a string are taken; a repeated key, a lone surrogate escape and nesting past 128 refuse the document wherever they are (serde_json: only where it reads) | **p16** | src/json.cyr's pre-pass rewrites the two integer forms and refuses the string defects; ADR 0002 rows 11 and 12 |
| Touchstone text: parse, write `{:e}` | `f64_parse_n` (correctly rounded, span-bounded); bayan's `bayan_f64_to_json` / `bayan_f64_parse` | stdlib ref *math*; bayan header | parsing is direct. Rust's `{:e}` output (`1e-1`, `-2.99792458e8`, `1.2345678901234568e17`) needs a small exponent-form printer over the shortest digits (impetus's `fmtf.cyr` has the digit generation) |
| `tracing` (67 call sites), `logging::init` | sakshi in the stdlib: `sakshi_error/warn/info/debug/trace`, `sakshi_log_kv`, spans, `sakshi_set_level`, outputs (fd, file, ring, UDP) | `lib/sakshi.cyr` | `log.cyr` wraps sakshi; `init` / `init_with_level` port as level setup |
| `ai` (reqwest + tokio + serde_json) | sandhi (HTTP/2, JSON-RPC, TLS) and bayan in the stdlib; prakash ported its daimon/hoosh client the same way | stdlib ref *Folded sibling distfiles*; `prakash/src/ai.cyr` | **ported**, in a `[lib.ai]` profile (TLS stays out of the default bundle) |
| `compute-wgpu` `WgpuBackend` | Rust's is a no-op scaffold: `new()` returns `Ok`, and it implements no trait. mabda 4.2.0 (stdlib) offers wgpu-native, native AMD (with an f64 SPIR-V compiler) and native NVIDIA compute | `src/compute.rs:459-491`; mabda README | **the scaffold is ported 1:1** (`BjWgpuBackend_new`). A working mabda backend is new work in either language (§11) |
| `soorat-compat` (8 structs, 14 tests) | **no soorat dependency at all**: bijli-owned data structs plus serde. goonj (Cyrius 2.3.1) ported its own `integration/soorat.rs` as `src/soorat.cyr` | `src/integration/soorat.rs`; `goonj/src/soorat.cyr` | **ported** (M7). Fixed-size array fields (`[f64; 3]`) become typed fields or a tuple field, since struct fields can't be `T[N]` |
| `#[inline]` (401 sites) | `#inline`, **ignored on a fn of more than two parameters**; small generic bodies inline themselves | **p2 ⚠** the compiler says *#inline ignored: fn takes more than 2 parameters* | a performance hint only; no semantic change |
| `const ADD: bool` generics (`row_diff`, `row_curl`) | no const generics; at most two type parameters | guide *Generic Functions* | two functions each (`_add`, `_sub`) |
| `pub` / private helpers | **`private` files with `public` items** (6.5.0): a hard error across files; private fns leave the symbol table and `cyrius api-surface` | guide *Visibility* | each `src/*.cyr` is `private`, with its Rust `pub` items `public`. The `_bj_` prefix stays, against collisions inside one compilation unit |
| scientific literals (`8.854_187_812_8e-12`) | no exponent form; decimal literals are correctly rounded; `const` takes no annotation (`const EPS0 = 0.0000000000088541878128;` is an f64) | guide *Number Literals*, *Constants*; **p6 ✅** bit-exact | generated `const`s, as decimals |
| `#[cfg(feature)]` gates | `[features]` (optional deps only) + `[lib.<profile>]` bundles; `#define` capped at 16 names | guide *Manifest keys*, *Preprocessor* | profiles: `[lib]`, `[lib.ai]` |
| tests over the library | `[test] modules` (prepended to every test, bench and fuzz unit); `[test.embed]` (embed golden files); `cyrius coverage --min N`; `cyrius api-surface`; `#bounds` | guide *Tests*, *Manifest keys* | ~~corrects impetus's `tests/lib_*.cyr` mirrors: list the modules in `[test] modules`~~ **Reversed in M1 part 2a** (port-parity.md §8): a prepended module becomes part of the test's own file, so a `private` src/ file makes the test file private and `cyrius coverage --exec` credits no src/ fn. Tests include `tests/lib.cyr`, the `[lib]` list as files of their own (impetus's pattern), checked against `[lib]`. Golden TSVs ride `[test.embed]`. CLAUDE.md's 80 % coverage target becomes `cyrius coverage --exec --min 80` |

Struct fields can't be fixed arrays (`T[N]`); an inline tuple field `(f64, f64, f64)` or named
fields replace them.

---

## 3. What is different from the earlier ports

1. **Complex arithmetic is everywhere.** It appears in `rf` (106 mentions), `polarization` (46),
   `scattering` (30), `beam` and NF2FF (a private `Cplx` in `fdtd.rs`). The earlier ports had
   almost none.
   - `BjComplex` is a by-value struct with operator fns (§2).
   - impetus measured struct-value kernels at 2.5× slower than `f64` locals. So the API and the
     cold paths use struct values, while the hot loops (Mie recurrences, NF2FF accumulation,
     the FDTD kernels) work on `f64` locals. A probe (§6 step 9) re-measures this for complex.
   - Never write `s.f = g(t.e)`: the open field-store defect silently stores the argument instead
     of the result. Bind the result to a local first. impetus's `surface_parity.py` scans for the
     pattern.
   - Mie (`scattering`) runs complex recurrences (the downward log-derivative D_n and
     Riccati–Bessel ψ/χ). Its public `a`/`b` coefficient vectors are golden outputs too.
2. **FDTD is a linear, non-chaotic time stepper.** impetus couldn't compare long trajectories,
   because a 1-ulp change grows in a pile of bodies. A Yee update uses only `+ − × ÷` in a fixed
   order, and Rust performs no FMA contraction (bijli calls no `mul_add`), so the same operation
   order gives the same bits at every step.
   - **Trajectory rows (class T, §7)** compare the whole field state at chosen frames of a run
     that starts at t = 0. No mid-run state has to be loaded, so no oracle hook is needed for the
     private CPML Ψ arrays.
   - A frame's state is large: a 30³ grid holds 162,000 values in its six components. A row
     therefore stores a **digest** of each component's bit patterns, plus a few probe cells in
     full. A digest mismatch reruns the generator in `--dump` mode on that frame, to find the
     first differing cell.
3. **Grids are large and long-lived.** They use freelist storage (§2: mmap'd above 4 KB, munmapped
   on free) and are released by `BjFdtd1d_free` / `BjFdtd2d_free` / `BjFdtd3d_free`. The Rust code
   allocates CPML Ψ slabs and ADE arrays only when they are enabled, and so does the port.
   **Never `alloc_reset`**: it would wipe the consumer's heap.
4. **Closures carry the callbacks** (§2). The one trap is float captures: rebind them as typed
   locals.
5. **Text I/O.** Touchstone parsing uses `f64_parse_n`. Writing needs Rust's `{:e}` form, and
   malformed-file rejection becomes golden rows.
6. **More transcendentals.** asin and tan aren't in either sibling's correctly rounded set (§0.2).
7. **A wide surface in one profile.** bijli has 430 `pub fn`s. tanmatra's `parity_core`
   per-function table overflows past 96 functions, so size it for 512 from the start (§6 step 7).

---

## 4. Scope

**Everything in the Rust crate is ported**, except the three functions its own roadmap deprecates
(§4.2). Nothing is deferred for want of a language feature: §2 shows each one exists.

### 4.1 In scope for 2.0.0

| Milestone | Rust module | Code lines | `pub fn` | Tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M1 | `error.rs` | 21 | — | 10 | `src/error.cyr` (`BJ_ERR_*`, the last-error record) |
| M1 | (constants in `field.rs`, `fdtd.rs`) | — | — | — | `src/constants.cyr` (generated and proven, as impetus) |
| M1 | (math) | — | — | — | `src/mathx.cyr`, `src/crmath.cyr`, `src/crmath_data.cyr` |
| M1 | `logging.rs` + every `tracing` site | 14 | 2 | 1 | `src/log.cyr` (sakshi; `init`, `init_with_level`, the `_bj_warn` family used from M1 on) |
| M1 | `polarization.rs` `Complex` helpers and the hisab 1.4.0 `Complex` | — | — | — | `src/complex.cyr` (`BjComplex` by value, `impl` methods, operators) |
| M1 | `field.rs` | 786 | 38 | 105 | `src/field.cyr` (`BjFieldVector` by value with operators, point fields, dipoles, multipoles, images, Green's functions, stress tensor, field lines with closure callbacks) |
| M1 | `maxwell.rs` | 108 | 10 | 18 | `src/maxwell.cyr` |
| M1 | `charge.rs` | 172 | 14 | 26 | `src/charge.cyr` |
| M2 | `material.rs` | 283 | 33 | 42 | `src/material.cyr` |
| M2 | `wave.rs` | 464 | 42 | 63 | `src/wave.cyr` |
| M2 | `circuit.rs` | 297 | 27 | 51 | `src/circuit.cyr` |
| M2 | `relativity.rs` | 517 | 27 | 48 | `src/relativity.cyr` (trajectory closures) |
| M3 | `polarization.rs` | 565 | 53 | 69 | `src/polarization.cyr` |
| M3 | `scattering.rs` | 250 | 6 | 32 | `src/scattering.cyr` |
| M3 | `beam.rs` | 411 | 28 | 60 | `src/beam.cyr` |
| M4 | `rf.rs` | 984 | 31 | 80 | `src/rf.cyr`, `src/touchstone.cyr`, `src/fmtf.cyr` |
| M5 | `fdtd.rs`: 1D, pulses, TF/SF, lossy and dispersive materials, NF2FF, subpixel, non-uniform grid, 4th order, CW, eigenmodes | 1,209 | ≈ 45 | ≈ 80 | `src/fdtd1d.cyr`, `src/fdtd_source.cyr`, `src/fdtd_media.cyr`, `src/nf2ff.cyr`, `src/fdtd_util.cyr` |
| M6 | `fdtd.rs`: CPML 2D and its helpers, `Fdtd2d`, CPML 3D, `Fdtd3d` | 1,691 | ≈ 47 | ≈ 85 | `src/cpml.cyr`, `src/fdtd2d.cyr`, `src/fdtd3d.cyr` |
| M6 | `compute.rs`: `ComputeBackend` (a trait), `CpuBackend`, `GpuFdtd2d<B>` (a bounded generic struct), the `WgpuBackend` scaffold | 313 | 9 | 18 | `src/compute.cyr` |
| M7 | `integration/soorat.rs` | 271 | 8 | 14 | `src/soorat.cyr` |
| M7 | `bridge.rs` (all 8: ushma, prakash, tanmatra, dravya) | 78 | 8 | 18 | `src/bridge.cyr` |
| M7 | `ai.rs` | 117 | 2 | 7 | `src/ai.cyr`, in the `[lib.ai]` profile only |

fdtd.rs's test count is split by section. Its `cpml_tests` and `cpml3d_tests` modules go to M6.
Surface parity (§8 gates) covers every `pub fn` and every `pub` method.

### 4.2 Not carried

| Rust surface | Why |
|---|---|
| `rf::input_impedance_lossless`, `rf::write_touchstone`, `Material::reflectance_at` | deprecated in 1.2.0, and bijli's roadmap ("2.0 cleanup") removes them. Their replacements are ported, and the Rust tests that call them are n/a in the ledger |

### 4.3 API changes in the Cyrius version (accepted divergences, ADR 0002)

These come from the roadmap's "2.0 cleanup (breaking)". Each row says how the harness still
compares it, so no function leaves parity entirely.

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | three deprecated functions | removed | through their replacements |
| 2 | public FDTD grid fields (`ex`, `ez`, `nx`, …), with a per-step layout check because callers can resize them | the buffers are fixed at construction and read through accessors (pointer + length); nothing can resize them, so the per-step layout check and its `tracing::error` path go | class T on the accessors. The Rust layout-mismatch tests are n/a |
| 3 | `TfsfSource::new` and `Cpml::new` sanitize with a warning; `Fdtd2d::enable_cpml` warns and ignores | fallible only: the constructors return `Err`, and `enable_cpml` is the `try_` form | oracle sanitize-or-ignore rows ↔ `Err`, mapped by name |
| 4 | `BijliError` variants carry data and `String`s | `BJ_ERR_*` codes (append-only) in `Result`; the payload in the per-thread last-error record (§2) | class A on the code **and** on the recorded payload |
| 5 | two time conventions (polarization, beam, field and the FDTD media in e^{−iωt}; scattering, rf and the FDTD phasors in e^{+jωt}) | **e^{+jωt} everywhere**, IEEE handedness (D3, §4.4) | amplitudes ↔ `conj(oracle)`, phases ↔ `−oracle`, observables unchanged: all exact |
| 6 | results bit-identical only on the same platform | **bit-identical across targets** | measured on three targets, as impetus |
| 7 | no transport error: the `ai` client maps network failures to `InvalidParameter` (roadmap backlog) | `BJ_ERR_TRANSPORT`, appended | oracle `InvalidParameter` rows from the `ai` client ↔ `BJ_ERR_TRANSPORT` |

### 4.4 One time convention (D3)

Rust 1.2.0 mixes the two harmonic time conventions, module by module, and in `fdtd.rs` within one
file:

| Rust module | Convention | Where it shows |
|---|---|---|
| `polarization` | e^{−iωt}, Born & Wolf handedness | Jones vectors and matrices; the waveplate's fast axis at e^{−iδ/2}; `right_circular` = [1, −i]/√2; S₃ = 2 Im(E_x E_y*) |
| `beam` | e^{−iωt} (field ∝ e^{+ikz}) | `hermite_gaussian_field`'s complex field and its Gouy phase |
| `field` | e^{−iωt} | the Green's function e^{ikR}/(4πR) |
| `fdtd` media | e^{−iωt} | Drude, Lorentz and Debye ε(ω) (the ADE update coefficients themselves are real) |
| `fdtd` CW / NF2FF | e^{+jωt} | `CwAccumulator` (`imag −= f·sin ωt`), `nf2ff_2d_phasor`, the e^{+jω·dt/2} H rotation |
| `scattering` | e^{+iωt} | Mie `a`/`b`; absorbing Im(m) < 0 |
| `rf` | e^{+jωt} | +jωL, Γ, Smith chart, S-parameters, Touchstone |

**The Cyrius 2.0.0 uses e^{+jωt} in every module** (written with `j`), with **IEEE handedness**
(IEEE Std 145), so that one consistent engineering system runs through the library. The
maintainer chose it on 2026-11-11 (D3). The reasons:
- It is the convention of the literature bijli implements: FDTD (Taflove & Hagness), RF (Pozar),
  antennas (Balanis) and Gaussian beams (Siegman). prakash's `wave_beam.cyr` uses it too.
- `rf`, Mie and the FDTD phasors already use it. Touchstone files and Smith-chart coordinates need
  no conversion: both are defined in e^{+jωt}.
- The modules that change are `polarization`, `beam`, the Green's function and the media's ε(ω).
  Absorbing indices are written n − jκ, as bijli's Mie already requires.

**What changes and what doesn't**, with how the harness compares each (ADR 0002 row 5). Every
mapping is exact, so no tolerance widens:

| Quantity | Change | Golden comparison |
|---|---|---|
| complex amplitudes: Jones vectors and matrices, waveplates, the complex beam field, the Green's function, complex ε(ω) | conjugated | port = `conj(oracle)`: the imaginary part's sign bit flips |
| phases (`arg` of an amplitude, the field's Gouy phase term) | negated | port = `−oracle` |
| circular-state names (`right_circular`, `left_circular`, the circular polarizers) | named by IEEE handedness | mapped by physical state: each name is pinned by a convention-free test that samples the real field E(t) at two instants and checks its rotation sense against the propagation direction |
| S₃ and the ellipticity angle's sign | keep their physical meaning under the IEEE names | the formula from a Jones vector changes sign; the rows map with the names |
| intensities, \|Γ\|, VSWR, return, insertion and mismatch loss, Mie Q_ext/Q_sca/Q_abs and g, cross sections, Rayleigh, S₀–S₂, Mueller matrices, degree of polarization, orientation angle, array-factor magnitudes, NF2FF patterns, every real FDTD field | unchanged | class as before, against the oracle directly |

**Rejected: e^{−iωt} everywhere** (Jackson, Born & Wolf, Bohren & Huffman). It would keep
`polarization` and `beam` and let tabulated n + iκ go into Mie without conjugation. But `rf`,
Mie and the FDTD phasors would change instead. Touchstone files would be conjugated at the file
boundary, because their spec is e^{+jωt}, and Smith-chart coordinates would come out conjugated
relative to every printed chart.

The convention is stated once, in `src/complex.cyr`'s header and the README, and every module's doc
points there instead of carrying its own. `docs/guides/migration-2.0.md` (M8) gives consumers the
conversion: conjugate amplitudes from bijli 1.x's `polarization`, `beam`, Green's function and
ε(ω). Circular-state names now follow IEEE.

---

## 5. Prep: everything before `cyrius port`

These steps run on the Rust tree. They end in a Rust-side **"port plan" commit**, which is what
tanmatra (`2febf11`), impetus (`a26296b`) and kimiya (`6281f19`) each landed before their
initial-port commit.

- **P0. Read this plan and settle D1, D3 and D5** (§0.3; D2 and D4 are settled by §2). Done 2026-11-11: D1 the 1.2.0 tag, D3 e^{+jωt}, D5 derive writers with hand-written readers. Record the answers in ADR 0001 / 0002 drafts.
- **P1. Tag the oracle.** ✅ *Done 2026-11-11.* There had been no 1.2.0 tag, so the maintainer
  tagged `c95fd30` (HEAD, `origin/main`) as `1.2.0`. `VERSION` and `Cargo.toml` both read 1.2.0,
  which is what the release workflow checks.
  - **The full P(-1) gate is green on the tagged tree.** `src/`, `tests/`, `benches/` and the
    manifests are byte-identical to `1e0a889`, where `cargo fmt --check`, clippy
    `--all-features --all-targets -D warnings`, `cargo test --all-features` (830 unit, 26
    integration, 2 doc), `RUSTDOCFLAGS="-D warnings" cargo doc`, `cargo audit` and
    `cargo deny check` passed. `make features` (each of the 16 features alone, plus
    no-default) passed on `c95fd30` itself.
  - **The 3D CPML is in 1.2.0.** The tagged CHANGELOG still lists it under `[Unreleased]`. On
    `main` the entry is folded into `[1.2.0]` (dated 2026-11-11, the tag's date); the tag itself
    keeps the old text unless it is re-tagged.
  - crates.io stops at 1.1.0, so the tag is the only 1.2.0 there is.
- **P2. Decide on oracle hooks: none planned.** impetus needed a hidden `oracle` feature to see
  state its API hid. bijli's API already shows almost everything:
  - the FDTD field arrays are public;
  - Mie's `a`/`b` coefficients are public;
  - class T starts every run at t = 0 (§3 item 2), so no mid-run state has to be loaded.

  Add a `#[doc(hidden)]` `oracle` feature in 1.2.1 only if an M0 probe shows a need. The likely
  candidates are the CPML σ/κ/α profile coefficients (`cpml_coefficients`) and the private Mur
  coefficient.
- **P3. Record the oracle's environment** in the plan's header table and in ADR 0001:
  - rustc 1.99.1, glibc 2.44 (Rust std's libm), Rust hisab 1.4.0 (`Complex`);
  - the `Cargo.lock` at the tag.
  - **Which calls actually reach glibc** in the generator's build profile. ✅ *Done 2026-11-11.*
    Method: every non-test call site in `src/`, then a disassembly of the release unit-test and
    bench binaries (`objdump -d`), with each libm call attributed to its calling function.
    - **glibc symbols reached:** `sin`, `cos`, `sincos`, `exp`, `pow`, `atan`, `atan2`, `asin`,
      `acos`, `tan`, `hypot`, `log10`. (`sinh` is reached only from tests; `log` in the bench
      binary is criterion's.)
    - **`cbrt` is not glibc's** (corrected at M0 step 5, when `validate.py` disagreed on 504 of
      1,084 `cbrt` rows). Rust std links compiler-builtins' own `cbrt`, a local symbol: the
      `libm` crate's port of CORE-MATH, which is correctly rounded. MPFR 4.2.2 agrees on all
      1,084 rows. So the oracle's `cbrt` rows are already the correctly rounded target.
    - **`sin` + `cos` of one argument become one `sincos` call** (`nf2ff_2d_phasor`, `mie`,
      `parse_touchstone`, the multipoles, `single_stub_match_complex`). glibc 2.44's `sincos`
      equals separate `sin` and `cos` bit for bit: 0 mismatches over 399,944 arguments in four
      ranges, checked through `ctypes`. So the corrected oracle treats it as two calls.
    - **`10.0_f64.powf(v1 / 20.0)` (`rf.rs:1247`, Touchstone dB) stays `pow`.** `parse_touchstone`
      calls `pow@GLIBC_2.29`; there is no `exp10` rewrite. bijli has no `2.0.powf` site.
    - **`powi` with a runtime exponent** (`beam.rs:643`, `fdtd.rs:4664`, `field.rs:934`) reaches
      compiler-rt's `__powidf2`, not glibc. Constant exponents expand inline instead. Both are
      square-and-multiply in the same order, and `_bj_powi` copies that order (class B).
    - **Constant folding:** with constant arguments LLVM folds `asin`, `acos`, `tan` and `cbrt`
      at compile time (none of their calls survive in the test binary). **The generator must feed
      every input at run time** (from its grids, through `std::hint::black_box`), so each golden
      row is what glibc computes at run time.
- **P4. Take a quiet Rust bench baseline at the tag.** Run `./scripts/bench-history.sh` on an
  idle machine. That run is the Rust column of `docs/benchmarks-rust-v-cyrius.md`: 80 benches,
  including `fdtd/3d_step_30x30x30_cpml8`.
  - ✅ **Done 2026-11-11 (maintainer's choice): the `2026-11-12T03:15:21Z` row at `c95fd30` is the
    Rust column.** It ran under load: kimiya's `cr_verify.py` and wine gate jobs, at a load
    average of about 2.6. Median CI width 0.8 %, p90 4.2 %, maximum 13 %;
    `fdtd/3d_step_30x30x30_cpml8` read 495 µs against 452 µs quiet.
  - The Rust-vs-Cyrius ratios carry that caveat. `docs/benchmarks-rust-v-cyrius.md` states it in
    its header, and a ratio within the Rust row's own CI width is read as parity, not as a
    finding.
  - Before Cyrius bench runs, check `ps` for other repos' gates and benches, not only the load
    average.
- **P5. Tidy what will move.** ✅ *Done 2026-11-11.*
  - `target/` (5.8 GiB) was removed with `cargo clean`, so `cyrius port` has no build cache to
    move. The generator builds `rust-old/` into its own target directory at M0.
  - No tracked build output, unlike impetus's `fuzz/target/`.
  - `.claude/` is ignored. Nothing else needs untracking.
- **P6. Pre-stage what the scaffold will overwrite.** `cyrius port` 6.7.18 still overwrites
  `CLAUDE.md`, `docs/development/roadmap.md`, `docs/development/state.md` and both workflows. It
  also pins `version = "0.1.0"`, leaves `bench` out of `[test] stdlib`, and its stubs fail
  `cyrius audit` (the open cyrius issue in §1). Everything is in git, so restore from `HEAD`
  after the scaffold.
  - ✅ *CLAUDE.md drafted 2026-11-11*: [`CLAUDE.cyrius.md`](CLAUDE.cyrius.md) merges bijli's
    CLAUDE.md, the scaffold's `claude-md-port` template, impetus's translation and this plan's
    rules. At M0 step 2 it is copied over `CLAUDE.md`, and the draft is deleted.
- **P7. Draft ADRs 0001 and 0002.**
  - bijli has no `docs/adr/` yet, so start directly in the scaffold's four-digit form.
  - [0001 — Port to Cyrius with Rust 1.2.0 frozen as the oracle](../adr/0001-cyrius-port-with-rust-oracle.md)
    is drafted with this plan.
  - 0002, the §4.3 API divergences, is written at M0 step 2.
- **P8. Land the "port plan" commit:** this file, the ADR 0001 draft, and the roadmap's pointer
  to this plan. No `src/` change.

**Gate for P:** 1.2.0 tagged and green ✅; D1, D3 and D5 answered; the bench baseline recorded; the port
plan commit landed.

---

## 6. M0 — Initialize with `cyrius port`

1. **Scaffold.** Run `cyrius port --dry-run .`, then `cyrius port .`.
   - The dry run reports 22,696 Rust lines.
   - Expect these to move to `rust-old/`: `Cargo.toml`, `Cargo.lock`, `src/`, `tests/`,
     `benches/`, `examples/`, `Makefile`, `deny.toml`, `codecov.yml`, `rust-toolchain.toml`
     and `target/`.
   - The dry run doesn't list them all (cyrius issue item 5), so check with `git status` after.
2. **Repair what the scaffold overwrites** (as impetus §4 step 2):
   - **CLAUDE.md**: `cp docs/development/CLAUDE.cyrius.md CLAUDE.md`, drop its draft note and
     delete the draft (P6). It already keeps the identity, process, DO NOT list and documentation
     structure, and translates the principles as impetus did:
     - `#[non_exhaustive]` → append-only enum values;
     - `#[must_use]` → `#must_use` / `#pure`;
     - no `unwrap`/`panic!` → no `+?`, `-?`, `*?`, `assert_fatal`, `result_unwrap` or `panic`
       in `src/`;
     - "tracing on all operations" → sakshi, behind a level;
     - "Feature-gate optional deps" → `[lib.<profile>]` bundles;
     - **never `hashmap`** where order matters (its seed is random per process);
     - **never `alloc_reset`**.

     Replace the cargo commands in the cleanliness check with `cyrius build`, `test`,
     `test --bounds`, `bench`, `vet`, `lint`, `fmt --check`, `distlib --check`, `deps --locked`
     and `audit`.
   - **Roadmap**: restore from git, then add the port section. **`state.md`**: write it for
     1.2.0 → 2.0.0-dev.
   - **Workflows**: keep the scaffold's `release.yml`, which is identical in all three siblings.
     Take impetus's `ci.yml`, with `test --bounds`, `deps --locked`, `audit`, the ledger,
     bench-twin and doc-example checks on from the start, plus an **`[lib.ai]` profile job**
     once M7 lands. Move the Rust workflows to `rust-old/.github/workflows/`. Append
     `/rust-old/target/` (not `/rust-old/`) and `/build/` to `.gitignore`.
   - **Move Rust-only material** to `rust-old/`: `scripts/bench-history.sh`,
     `scripts/version-bump.sh`, `bench-history.csv`, `benchmarks.md` and
     `docs/architecture/overview.md`. Rewrite the overview for Cyrius in M8.
   - **ADRs**: 0001 (the port) and 0002 (§4.3), with `docs/adr/README.md` and `template.md` from
     the scaffold.
   - `cyrius.cyml` reads `version = "${file:VERSION}"`, with `bench` in `[test] stdlib`.
     Document the six stub functions so `cyrius audit` passes.
3. **`cyrius.cyml`.** Confirm the keys with `cyrius help manifest`.
   ```toml
   [package]
   name = "bijli"
   version = "${file:VERSION}"
   license = "GPL-3.0-only"
   language = "cyrius"
   cyrius = "6.7.18"

   [deps]
   stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita",
             "tagged", "result", "freelist", "fnptr", "sakshi", "overflow", "bayan"]

   [test]
   stdlib  = ["assert", "bench"]
   modules = ["src/error.cyr", ...]   # the [lib] list, so every test, bench and fuzz unit sees the library
   timeout = 300

   [test.embed]                       # golden corpora compiled into the parity units
   GOLDEN_FIELD = "tests/golden/field.tsv"

   [lib]            # the default bundle: everything except ai
   modules = ["src/error.cyr", "src/constants.cyr", "src/mathx.cyr", "src/crmath_data.cyr",
              "src/crmath.cyr", "src/log.cyr", "src/complex.cyr", "src/field.cyr", ...]

   [lib.ai]         # the default modules plus src/ai.cyr (sandhi, bayan): dist/bijli-ai.cyr
   modules = [...]
   ```
   - `[deps] stdlib` also takes `thread_local` (the last-error record) and `mmap` (freelist's
     large blocks); `[lib.ai]` needs `sandhi` and its companions, as prakash lists them.
   - **Module order** follows Rust's feature graph:
     - error → constants → mathx → crmath → complex → log;
     - then field → maxwell → charge;
     - then material → wave → circuit → relativity;
     - then polarization → scattering → beam → rf;
     - then the fdtd files → compute → soorat → bridge.
   - A profile repeats the whole list, and the order is load-bearing (prakash's comment in its
     `cyrius.cyml`).
   - No `[deps.hisab]` (§0.2).
4. **Toolchain and math check.** Confirm the Cyrius counterpart of every Rust f64 method bijli
   calls (§9.1). Port impetus's `mathx.cyr` helpers as `_bj_*`: Rust-semantics
   min/max/clamp, `is_nan`, `is_finite`, `copysign`, `to_radians`/`to_degrees` and the literal
   `powi`.
5. **Golden-value generator** in `scripts/oracle/`: copy impetus's crate and adapt it.
   - TSV format 1, with f64 values as bit patterns; deterministic; a `--check` mode; glibc in
     the header.
   - `mpfr_ref.py` and the `.cr.tsv` corrected rows (kimiya ADR 0006).
   - **Complex values** take two bit-pattern columns.
   - **Trajectory format** (new): a row holds a scenario name, a frame index and a component.
     Its output is a 64-bit digest of the component's bit patterns in index order, plus the bits
     of the probe cells. The same digest is in `validate.py` and `tests/parity_core.cyr`.
     `--dump <scenario> <frame>` writes the full component for diffing.
   - **Closure inputs**: the scenarios name the field and trajectory functions, which exist on
     both sides: a point charge, a dipole, uniform motion, a circular orbit at β = 0.99 and a
     hyperbolic trajectory.
6. **Table and constant generator.**
   - The constants: CODATA 2018 `EPSILON_0`, `MU_0`, `COULOMB_K`, `SPEED_OF_LIGHT`(`_SQ`) and
     `ETA_0`; `FOURTH_ORDER_C1`/`C2`; and the limits `MIE_MAX_SIZE_PARAMETER`, `MIE_MAX_MX`,
     `CPML_MIN_INTERIOR` and `TRACE_INITIAL_CAPACITY`.
   - The material presets.
   - They are extracted and proven as impetus's `constants_private.cyr`.
7. **Parity harness and port tooling.** Copy impetus's `tests/parity_core.cyr` (with its
   `_par_error_name` hook) and `tests/ulp.cyr`, and kimiya's `_par_note_p` / `_par_force_map` for
   the ADR mapping rows. Size the per-function table for 512 functions.
   - Add a complex comparator. Each part is classed separately, so a row reports which part
     failed.
   - Add the trajectory runner and the digest comparator.
   - Add the divergence table that reports each ADR 0002 mapping by name.
   - ~~Tests see the library through `[test] modules`~~ (§2). **Reversed in M1 part 2a**: tests include
     `tests/lib.cyr`, the `[lib]` modules as files of their own (impetus's `tests/lib_*.cyr` pattern), and a
     check keeps it equal to `[lib]` (port-parity.md §8).
   - The harness reads `tests/golden/` at run time, as impetus's does (`file_read_whole`).
     `[test.embed]` stays available if a unit ever has to run without the files.
   - **Coverage**: `cyrius coverage --exec --min 80`, CLAUDE.md's target, from M1 on. It measures
     public fns only, and M0's `src/` has none (the `_bj_` helpers and the constants).
   - **API surface**: `cyrius api-surface --update` writes the public-API snapshot that the
     surface-parity script compares with Rust's `pub` items.
   - **Bring M8's tooling forward**, so it runs at every milestone, not just at release:
     - `scripts/surface_parity.py`, with the field-store scan;
     - `scripts/rust_test_ledger.py --check`;
     - `scripts/bench_record.py`, which maps the bench twins;
     - `scripts/doc_examples.py`;
     - a bundle smoke test (`tests/bijli.tcyr`, as kimiya);
     - a fuzz harness for the Touchstone parser.
   - Probe programs go in `scripts/probes/` (as impetus).
8. **Baseline measurements** go in `docs/development/port-parity.md`.
   - The math sweep runs every counterpart against glibc and against MPFR on x86_64 ELF, x86_64 PE
     and aarch64.
   - The new correctly rounded primitives (asin, tan) are gated at 0 ulp against MPFR on
     ≥ 10⁶ pairs, as impetus measured its pow and atan2. ganita's `ganita_f64_asin` /
     `ganita_f64_tan` are measured beside them, as the fallback's cost.
9. **Probes.** The 2026-11-11 language review (§2, `scripts/probes/`) has
   already answered six. Move them to `scripts/probes/` and re-run them at the M0 pin:
   - ✅ complex by value with operators, arrays and pointer stepping (p1);
   - ✅ traits, a bounded generic struct and `dyn` (p2);
   - ✅ closure callbacks, and the float-capture trap with its two workarounds (p3, p3b, p6);
   - ✅ grid memory through freelist: flat over 1,000 cycles (p4);
   - ✅/⚠ `#derive(Serialize)` against serde_json (p5, p5b);
   - ✅ `Result` with an f64 payload, `(f64, f64)` multi-returns, decimal `const`s (p6).

   ✅ **Answered 2026-11-11** (M0 Part 3, [`port-parity.md`](port-parity.md) §4, probes p7–p14). The
   questions were:
   1. **Grid memory at full size.** p4 with real `BjFdtd3d` grids (CPML and dispersive arrays
      on), and freelist against a single mmap per grid for 53³.
   2. **Struct values against `f64` locals in the hot loops**: a Mie recurrence step and an
      NF2FF accumulation, both ways (impetus measured 2.5× for vectors).
   3. **Callback cost.** `trace_field_line`'s 10⁴-step loop through `callptr` on a capturing
      closure, against a direct call.
   4. **The FDTD inner loop.** `row_diff`/`row_curl` written four ways: `*f64` indexing, pointer
      stepping, `f64v_*` (mind the open scalar-tail bug) and unrolled. Gate: bits equal to the
      Rust kernel. The fastest form is recorded for M6.
   5. **Rust iterator semantics.** 13 `.sum()` sites: Rust's float `Sum` starts from −0.0, so
      check the empty and all-negative-zero cases. Also `.min_by`/`partial_cmp` tie order (the
      eigenmode peak search) and `chunks`/`windows` bounds.
   6. **Touchstone text.** `f64_parse_n` against Rust's `parse::<f64>` on the 1.2.0 test
      corpus; Rust's `{:e}` form on every value in the 201-point bench file.
   7. **Test runtime.** 3D CPML tests are costly: a 53³ reference run for 300 steps is about 45M
      cell-steps. Measure the Cyrius cost per cell-step, and size the class-T scenarios to fit
      `cyrius test` in CI. Bigger scenarios become class R in a separate slow file.
   8. **The last-error record.** `thread_local.cyr` holding `BJ_ERR_*` payloads, read after an
      `Err` on two threads.

**Gate:** `cyrius build` and `cyrius bench` pass; `cyrius test` runs the harness (0 modules)
with no warnings; `cyrius lint` and `fmt --check` are clean; `cyrius audit` exits 0. Committed:
the generator, `tests/ulp.cyr`, ADRs 0001 and 0002, the probe results in `port-parity.md` and the
upstream issues (§10).

---

## 7. Parity policy

tanmatra's classes apply, with impetus's and kimiya's corrected oracle, plus **T** for the FDTD:

| Class | Tolerance | bijli examples |
|---|---|---|
| **A — exact** | bit-identical / equal integers | error codes, validation verdicts, constants, material presets, Touchstone parse results and the bytes written, the Smith-chart and stub-match branch choices, `n_max` of Mie, eigenmode peak indices |
| **B — arithmetic** | ≤ 2 ulp, expected 0 | `+ − × ÷ √` paths: point fields, superposition, Coulomb, Lorentz force, Poynting, impedances, circuits, ABCD composition, Jones/Mueller products, complex `+ − × ÷` and `conj`, the subpixel averages, `fourth_order_derivative` |
| **C — transcendental** | **0 ulp against the corrected oracle** (kimiya ADR 0006) | anything that calls sin, cos, exp, atan2, hypot, pow, log10, atan, asin, acos, tan or cbrt: dipole patterns, Fresnel with angles, skin depth, Gaussian pulses, Ricker wavelets, beam Gouy phase, LG/HG modes, Mie, array factors, complex `sqrt`/`arg`/`from_polar` |
| **D — iterative** | relative ≤ 1e-12, or a per-function bound set before its first comparison | expected to be **empty**. The retarded-time solver (bracket doubling, then Newton with bisection) and the field-line tracer are deterministic sequences of the classes above, so they are compared as B/C, step for step. Class D is used only if a hisab or ganita routine replaces bijli code |
| **T — trajectory** | the digest of every field component equals the oracle's at each listed frame; probe cells per their path's class (B, or C where a source term uses exp) | `Fdtd1d`, `Fdtd2d` (TM, TE), `Fdtd3d`, each with PEC/PMC walls, Mur, CPML (classic and CFS) and lossy, Drude, Lorentz and Debye media; TF/SF; `GpuFdtd2d` on `CpuBackend`. Frames {1, 2, 10, 100, last} |
| **R — reference** | the physical tolerance stated in the test | CPML reflection (−88 dB 2D, −80 dB 3D), dispersive ε(ω) against the analytic form (< 1e-5), Mie against Bohren & Huffman tables, the Lorentz/energy invariants, the Rust unit tests that assert physics |

Rules:
- **Never widen a class to make a row pass** (as tanmatra). A failure is one of three things:
  - a port defect;
  - a Cyrius or ganita defect, which is filed upstream;
  - an oracle defect, which is fixed in Rust 1.2.x and regenerated.

  The only exception is a divergence accepted by ADR, and the harness still reports it by name.
- **A class T digest mismatch is never "close enough".** Bisect to the first frame and cell that
  differ. Then find the operation order that differs (§9.3 FDTD kernels).
- **Branch points on both sides:**
  - Fresnel total internal reflection and the `cos θ_i` 1e-6 clamp;
  - the Mie limits (`MIE_MAX_SIZE_PARAMETER`, `MIE_MAX_MX`) and the gain-index rejection (Im(m) > 0; absorbing indices are n − jκ);
  - the r² < 1e-30 m² singular cutoff;
  - the CPML fit check `2·t + 3 ≤ min(n)`;
  - `dt` at and just above the Courant limit;
  - the matched-load stub `(0.0, 0.25)`;
  - the damping regimes in `circuit`;
  - β → 1 in `relativity`.

---

## 8. Milestones

Every milestone, in Parts where it is large, ends with impetus's gates:
- `cyrius build`; `cyrius test` and `test --bounds` **with no warnings**;
- `cyrius bench`, recorded against the Rust numbers;
- `vet`, `lint` and `fmt --check` on the milestone's files;
- `distlib --all`, then `distlib --check`, then `deps --locked`;
- **surface parity**: every in-scope Rust `pub fn` has a Cyrius counterpart or a written reason;
- the grid soak (§6 step 9.1) from M5 on;
- the parity digests equal on x86_64 ELF, x86_64 PE and aarch64;
- the CHANGELOG and `state.md`.

Gate mechanics the siblings learned the hard way:
- Run the targets **one at a time**. wine under parallel load truncates output and still exits
  0, so check that the `=== parity summary` line is present.
- Never pipe `cyrius audit` through `tail`, because the exit status you see is then `tail`'s.
- `cyrius audit` fails on comment lines over 120 characters.
- One compile is capped at 24 MB of source, and `cyrius bench` stops a file after 300 s. Split
  the corpora and the `.bcyr` files.
- distlib's sidecar misses stdlib leaves, so `src/` files that need one `include` it
  explicitly.
- Run Python as `python3 -I -B`.
- The full gate takes about 25 minutes, run sequentially.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation** | error, constants, mathx, crmath, complex, log, field, maxwell, charge | `BijliError` → `BJ_ERR_*`. Rust hisab 1.4.0's `Complex` transcribed by value (§9.2), with `complex_real`/`zero`/`from_polar` and `norm`/`norm_sq`/`arg`. The correctly rounded primitives, including the new atan/asin/tan. `BjFieldVector` with operators. Point, line, ring, disk and sphere fields; dipoles; superposition; the `associated_legendre` multipoles; images; Green's functions; the stress tensor; `trace_field_line` with callbacks. Maxwell and charge functions. The JSON (D5) of `Complex`, `FieldVector` and `PointCharge` | every function over seeded grids including 0, ±tiny, ±huge, NaN and ±∞; positions near the singular cutoff; multipoles l ≤ 10; field lines for point, dipole and quadrupole fields, with singular and escaping endings |
| **M2 Media, waves, circuits, relativity** | material, wave, circuit, relativity | `MagneticType` as an enum. Clausius–Mossotti, Curie–Weiss and the Fresnel family (the trig-free `*_direct` forms, Schlick, TIR). RC/RL/RLC with the damping regimes. Lorentz factor, four-vectors, the EM tensor, boosts, Liénard–Wiechert, and `solve_retarded_time` with callbacks | each Fresnel function over θ ∈ [0, π/2] × n₁/n₂ including TIR and the clamp edge; the circuit regimes; β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9}; retarded time for the five scenario trajectories |
| **M3 Optics** | polarization, scattering, beam | Jones vectors and matrices, Stokes, Mueller, converted to the one convention with IEEE handedness names (§4.4), each circular state pinned by a rotation-sense test. Mie with its downward D_n recurrence and complex ψ/χ; its `a`/`b` vectors are golden too. Rayleigh. Beams: q-parameter, ABCD, resonator stability, `hermite_polynomial`/`laguerre_polynomial`, HG/LG fields, the complex field in the one convention | Mie over x ∈ {1e-3 … 1e5} × m including high index, absorbing and the rejection limits; every Jones element over seeded angles; HG/LG up to order 10 at z ∈ {0, ±z_R, 10 z_R} |
| **M4 RF** | rf, touchstone, fmtf | Transmission lines, VSWR, Smith chart, S-parameters and cascades, stub matching, the array factor, `input_impedance_complex`. Touchstone read and write in the spec's `N11 N21 N12 N22` order: the rejection of malformed files and the scientific output, byte for byte | the 1.2.0 Touchstone tests as files; seeded round trips; every malformed case; stub matching including the matched load |
| **M5 FDTD core** | fdtd1d, fdtd_source, fdtd_media, nf2ff, fdtd_util | `Fdtd1d` with the per-end Mur coefficient. Gaussian, Ricker, differentiated Gaussian and modulated pulses. `TfsfSource` (with Mur at both ends; fallible, §4.3 #3). Lossy, Drude, Lorentz and Debye ADE. Subpixel averages. NF2FF (J and M, complex accumulation) on `BjComplex`. `NonUniformGrid`, `fourth_order_derivative`, `CwAccumulator`, `find_eigenmodes` (DFT and half-power Q) | class T for 1D: vacuum, ε_r = 4, lossy, each dispersive model, at frames {1, 10, 100, 1000}; pulses over seeded t; NF2FF on 80 × 80 × 180 angles; eigenmodes on 4096 × 256 |
| **M6 FDTD 2D/3D** | cpml, fdtd2d, fdtd3d, compute | Freelist grid storage (§2) and the kernels from probe 4. Roden–Gedney CPML in 2D (strips) and 3D (12 Ψ slabs, κ, b/c at the Yee positions) as post-loop corrections. TM/TE, conductivity C_a/C_b, `set_material` validation before any write, `total_energy`. `ComputeBackend` as a trait, `CpuBackend` as its impl, `GpuFdtd2d<B>` as a bounded generic struct, the `WgpuBackend` scaffold | class T: 2D TM and TE at 100 × 100, PEC/PMC, CPML 10 (classic, CFS), lossy; 3D at 30³, PEC and CPML 8; `GpuFdtd2d` against `Fdtd2d`. Class R: the −88 dB and −80 dB reflection tests and the energy decay to 2e-8, sized by probe 8 |
| **M7 Integration surfaces** | soorat, bridge, ai | The soorat visualization structs (`FieldSlice2D`/`3D`, `FieldLineVisualization`, `ChargeVisualization`, `RadiationPattern`, `VectorFieldSample`) with their constructors from FDTD grids, field lines, charges and NF2FF, and their JSON (D5), as goonj ported its own. The 8 bridges. The daimon/hoosh client in `[lib.ai]`, on prakash's `ai.cyr`: config JSON, a redacting debug print, client construction | each soorat constructor from M1/M5/M6 scenarios, and its JSON compared by parsed bits against serde_json's; the bridge functions over seeded grids; the config JSON round trips; the `[lib.ai]` profile builds and its tests run in CI |
| **M8 Release 2.0.0** | all | Close the Rust-test ledger, which has run since M0: each of the 830 unit, 26 integration and 2 doc tests is ported, folded into golden rows, or marked n/a with its reason (the deprecated functions and the per-step layout check, §4.3). Nothing is "deferred". The full surface-parity diff. `docs/benchmarks-rust-v-cyrius.md` with all 80 benches in both columns. `docs/guides/migration-2.0.md`: Rust idioms, the ADR 0002 changes, and `BjComplex` replacing hisab's `Complex`. README, overview, getting-started, SECURITY. `cyrius distlib --all`, `cyrius audit`, VERSION 2.0.0 | the full corpus green, with equal digests on x86_64 ELF, x86_64 PE and aarch64 |

Order: P → M0 → M1 → … → M8. The reasons for this order:
- M1–M4 are pure functions, the shape of tanmatra and kimiya, and they build the complex type and
  the correctly rounded primitives that M5–M6 rely on.
- The FDTD comes last because it is the biggest piece (2,900 Rust lines). It is also the one
  that needs the memory and kernel probes settled.
- Split M5 and M6 into Parts of about 600 Rust lines: M6 Part 1 CPML 2D and `Fdtd2d`; Part 2
  CPML 3D and `Fdtd3d`; Part 3 compute.

**Pacing**, from the siblings:
- impetus ported 10.4k lines in 24 Parts of about 400–800 Rust lines, 20–50 minutes and one
  commit each.
- tanmatra (4.2k) took one commit per milestone.
- kimiya's time went mostly to retrofitting correct rounding and its overflow audit. §0.1 builds
  both in from the start to avoid that.

At impetus's granularity, bijli's 8.6k lines come to about 16–20 Parts:

| Milestone | Parts |
|---|---|
| M1 | 3 |
| M2 | 2 |
| M3 | 2 |
| M4 | 2 |
| M5 | 2 |
| M6 | 3 |
| M7 | 1 |
| M8 | 3 |

Each Part ends with:
- the gates;
- the range-edge review of §0.1, at a milestone's last Part;
- an independent Python model of the new functions, written from `rust-old/` by a background
  agent, as kimiya's `handoff/verify/`. All of kimiya's model mismatches but one were model
  mistakes, and that one was found by targeted inputs.

The maintainer does all git. A Part ends by asking whether to commit it.

---

## 9. Rust → Cyrius mapping (bijli specifics)

### 9.1 Math functions

Every std f64 method in `src/` outside tests, at the `1.2.0` tag:

| Rust (calls) | Cyrius | Parity |
|---|---|---|
| `sqrt` (72), `abs` (81), `ceil` (2) | `f64_sqrt`, `f64_abs`, `f64_ceil` (builtins) | exact on every target |
| `sin` (27), `cos` (21), `sin_cos` (8) | `_bj_cr_sin`, `_bj_cr_cos`, `_bj_cr_sin_cos` (from impetus `crmath.cyr`) | class C, 0 ulp against the corrected oracle |
| `exp` (11) | `_bj_cr_exp` (impetus) | class C |
| `atan2` (5) | `_bj_cr_atan2` (impetus) | class C |
| `powf` (4) | `_bj_cr_pow` (impetus) | class C |
| `hypot` (4) | `_bj_cr_hypot` (impetus) | class C |
| `log10` (5) | `_bj_k_cr_log10` (from kimiya) | class C |
| `acos` (3) | `_bj_k_cr_acos` (from kimiya) | class C |
| `atan` (4) | `_bj_cr_atan2(y, 1.0)` | class C; exact, since both are the rounding of the same real |
| `asin` (4), `tan` (2) | **new**: `_bj_k_cr_asin` (on kimiya's atan core), `_bj_cr_tan` (on impetus's sin/cos phases, with a fixed-point slow path) | class C; gated on ≥ 10⁶ MPFR arguments each (M1 parts 1b, 1c) |
| `cbrt` (2) | `ganita_f64_cbrt` (correctly rounded) | class C |
| `sinh`, `cosh` (through hisab's complex `sin` / `cos` only) | **new** (M1 part 2d): `_bj_cr_sinh`, `_bj_cr_cosh` on impetus's exp phases, sinh's small arguments by series | class C; 10⁶ MPFR arguments each |
| `ln` (through hisab's complex `ln` only) | `_bj_k_cr_ln` (from kimiya, M1 part 2c) | class C |
| `powi` (3) | explicit products in LLVM's square-and-multiply order (as tanmatra) | class B |
| `min` (40), `max` (18), `clamp` (8) | `_bj_min` / `_bj_max` / `_bj_clamp`, Rust's NaN-ignoring semantics | exact |
| `is_nan` (62), `is_finite` (101) | `x != x`; `x - x == 0` | exact |
| `to_radians` (2), `to_degrees` (2) | `x * (PI / 180.0)` and `x * (180.0 / PI)` in Rust's order | class B |
| `copysign` (1) | a sign-bit merge on the bits | exact |

**Typing traps** (as impetus §7.1, probed on 6.7.18):
- An operator is typed by its left operand. Bind every builtin and ganita result to a `: f64`
  local first; `(0 - x)` is an integer subtraction, so write `-x`.
- There is no exponent notation, so `8.854_187_812_8e-12` becomes a generated hex constant.
- A hex pattern in a `const` is an integer; `F64_*` stdlib constants are untyped `var`s.

### 9.2 Complex: Rust hisab 1.4.0 → `BjComplex`

`src/complex.cyr` transcribes `hisab-1.4.0/src/num/complex.rs`, operation for operation, so
that class B/C rows can be compared:

| Rust hisab 1.4.0 | Formula (keep exactly) | Note |
|---|---|---|
| `abs` / bijli `norm` | `sqrt(re·re + im·im)` | no hypot rescue; Cyrius hisab's `cx_abs` has one, so it would differ. In range: Rust's bits. Beyond about 2^±511: a power-of-two prescale, as an ADR 0002 row (§0.1) |
| `norm_sq`, `conj`, `+ − ×`, `× f64`, `neg` | the textbook forms, Rust's order | class A (M1 part 2c: the same IEEE operations in the same order, so bit for bit) |
| `/` | `((a·c + b·d)/den, (b·c − a·d)/den)`, `den = c² + d²` | naive. bijli's own Smith `cdiv`/`crecip` (rf) and `complex_div` (scattering) are ported as written too |
| `/ f64`, `inv` | componentwise | class A |
| `arg` | `atan2(im, re)` | class C |
| `from_polar(r, θ)` | `(r·cos θ, r·sin θ)` | class C |
| `exp` | `e^re · (cos im, sin im)` | class C |
| `sqrt` | `from_polar(sqrt(abs), arg/2)` | class C. Not the C99 branch of Cyrius hisab's `cx_sqrt` |
| `ln`, `powf`, `sin`, `cos`, `is_zero` | as written | ~~ported only if a non-test call exists~~ **all ported**: bijli 1.2.0 re-exports `Complex` (`pub use hisab::num::Complex` in polarization.rs), so its whole surface is bijli's API. `ln` on kimiya's correctly rounded ln, `powf` and `is_zero` (M1 part 2c); `sin` and `cos` on bijli's new correctly rounded `sinh` and `cosh` (M1 part 2d). Display (`"{re}+{im}i"`) through src/fmtf.cyr |

### 9.3 Types and idioms

| Rust | Cyrius |
|---|---|
| `FieldVector` (`Add`, `Sub`, `Mul<f64>`, `Neg`, `*Assign`, `f64 * v`) | `BjFieldVector` struct value (24 B) with operator overloads and explicit `_scale_left` |
| `Complex` | `BjComplex` struct value (16 B, two registers), §9.2 |
| `JonesVector`, `JonesMatrix`, `StokesVector`, `MuellerMatrix` | struct values. Matrices are flat fields (`m00 … m33`) because struct fields can't be fixed arrays |
| `Vec<f64>` grid fields | **one mmap per grid**, sliced into its `*f64` arrays (zeroed by the kernel; 10× faster than `fl_calloc` per array, p7), unmapped by `Bj*_free`; accessors return pointer + length; typed through `_bj_alloc_f64(n): *f64` (an untyped word into a `*T` warns, a typed return does not) |
| `Vec<Complex>` (Mie `a`/`b`), `Vec<Vec<Complex>>` (`SMatrix`) | one stride buffer of `BjComplex` with n (or n × n × freq) and an index helper |
| `Result<T, BijliError>` | value-form `Result` + `BJ_ERR_*` code, with `?` (p6); the error's payload in the last-error record. A result over one word (complex, field vector) is the code plus an out-pointer. No `Result` through `callptr` |
| `Option<f64>` (`half_power_offset`), tuples (`(f64, f64)` stub match, `(f64, f64, f64)` CPML coefficients) | `Option` (value form) for one word; a declared `: (f64, f64)` / `: (f64, f64, f64)` multi-return (p6) for the tuples; an out-struct past three values |
| `impl Fn` / generic `F: Fn` | a closure (or `&fn`) called with `callptr`: `\|pos, out\| -> status`, its float captures rebound as typed locals (§2, p3/p3b/p6) |
| `ComputeBackend` trait + `GpuFdtd2d<B>` | `trait BjComputeBackend` with `impl BjComputeBackend for BjCpuBackend`; `struct BjGpuFdtd2d<B: BjComputeBackend>` (p2). The associated `Buffer` type is an opaque handle word. `BjWgpuBackend` is Rust's scaffold, ported 1:1 |
| `row_diff::<ADD>` / `row_curl::<ADD>` | two functions each (`_add`, `_sub`), since there are no const generics, on `f64v_*` with a scalar tail for the open lane-multiple bug: bit-identical to Rust and 4.5× faster than the scalar forms (p10) |
| `try_filled` (`try_reserve`) | `fl_calloc` returns 0 → `BJ_ERR_ALLOCATION` (a new code, appended) |
| `checked_mul` in `grid_cells` | `overflow.cyr` |
| `tracing::{warn, error, debug}` (67 sites) | `_bj_warn` / `_bj_error` / `_bj_debug` over sakshi, behind a level |
| iterator chains, `.sum()`, `.zip()`, `.windows()` | loops in Rust's order. `.sum()` is a left fold from −0.0 (p11) |
| `#[non_exhaustive]` / `#[must_use]` | append-only enum values / `#must_use` |
| `#[derive(Serialize, Deserialize)]` (36 types) | `#derive(Serialize)` writers; readers per D5 |
| `pub` / private | `private` file + `public` items (§2) |
| `Drop` of a grid | an explicit `Bj*_free`; `defer` in tests and examples |

**Typing traps beyond §9.1** (from the siblings' `port-parity.md`):
- An **undeclared** multi-return (`return (a, b);` with no `: (T1, T2)`) gives untyped bindings,
  so `-x` on one is integer negation. Declare `: (f64, f64)` and the bindings are f64 (p6; this
  corrects the siblings' general claim). Arity is 2 or 3; four values don't parse.
- A captured `f64` is untyped inside a closure (§2, p6): rebind it. A closure's parameters are untyped words:
  bind them to typed pointers through a typed-return helper (`_bj_as_vec(p): *BjFieldVector`), since a direct
  `var q: *T = param;` warns (p9).
- `var x = f()` with `fn f(): *f64` leaves `x` untyped: write `var x: *f64 = f();` (p13).
- A `var u = x` copy of an f64 is f64-typed, which can turn a NaN self-compare into a trap.
- `&c.field` is refused.
- `#inline` takes at most two parameters (p2: the compiler says so and ignores it).
- `f64_min`/`f64_max` are NaN-asymmetric, `f64_round` ties to even, and `f64_sign` isn't Rust's
  `signum`.
- `vec_sort_by` is not stable, so sort on (key, index). bijli has one sort site, the eigenmode
  peak search.

### 9.4 Names

- **Public structs** `Bj<Name>` with `impl` blocks, so methods keep their Rust names:
  `g.step_tm()`, `g.try_enable_cpml(t)`. Constructors are `BjFdtd2d_new`, `BjPointCharge_new` and
  so on. Free functions keep the Rust name with the `bj_` prefix: `bj_electric_field`,
  `bj_mie`, `bj_fresnel_rs`.
- **Globals and enum members** carry `BJ_`; private helpers carry `_bj_`.
- The prefix scan is clean today: no `bj_`, `BJ_`, `_bj_` or `Bj[A-Z]` in the 6.7.18 stdlib,
  hisab, prakash or tanmatra bundles. Re-scan at M1 for unprefixed collisions such as `vec_*`,
  `Ok`/`Err` and ganita's `mat_*`. "Last definition wins" across files is only a warning.

### 9.5 Memory

Settled by the language review (§2, p4); M0 probe 1 confirms it at full grid size.
- **Grid lifetime**: a grid's field, coefficient, material, Ψ and ADE arrays are one anonymous mmap,
  sliced (p7: flat over 200 cycles of a 19.2 MB grid; 0.85 ms per grid against 8.7 ms for
  `fl_calloc` per array, whose zeroing is a byte loop), released by one munmap in `Bj*_free`.
  Smaller per-call buffers use freelist.
- **Call lifetime**: NF2FF accumulators, eigenmode spectra and Touchstone buffers are allocated
  per call, and their caller frees them. Field-line polylines go into a caller-owned growable
  buffer.
- **No per-step allocation** in any `step_*`: the kernels work in place, as Rust's do.
- **Never `alloc_reset`.**

---

## 10. Gaps to file upstream (M0, then as found)

Each one gets an issue in the owning repo, with the probe or benchmark that shows it.

**Shared, filed by bijli as the third port to need it**
- **G1. Correctly rounded f64 primitives as a library.** impetus (`_im_cr_*`, 72 functions) and
  kimiya (`_km_cr_*`, 175 functions) each wrote their own, and bijli needs a third set. Proposal:
  ganita, or a shared bundle, provides correctly rounded exp, ln, log10, pow, sin, cos, sincos,
  tan, atan, atan2, asin, acos, hypot, expm1 and cbrt. Then bijli and the next ports depend on it
  instead of copying.

**ganita / hisab** (special functions bijli owns today; filed for the ecosystem, not blocking)
- **G2.** Spherical and Riccati–Bessel functions and Hankel functions, with complex arguments
  (Mie); cylindrical J_n/Y_n and their zeros (`cylindrical_waveguide_cutoff`).
- **G3.** Associated Legendre P_l^m and spherical harmonics (`field` multipoles).
- **G4.** Hermite H_n and generalized Laguerre L_p^l polynomials (`beam`).
- **G5.** Complex values without a heap allocation per result (hisab D011; impetus's H1/H9).

**Cyrius stdlib / compiler** (all open; add any new ones the probes find)
- `cyrius port` overwriting docs (open, 2026-11-08). Add a ✚ to the issue for bijli.
- The small-struct `Result` payload bug (open, 2026-11-08): it hits every complex-valued
  `Result`.
- Field-store of a struct call copies a field argument (open, 2026-11-09).
- No scoped allocator, S1 (open, 2026-11-09).
- `f64v_*` writes past `dst` when n isn't a lane multiple (open, 2026-11-06). It matters for the
  FDTD kernels.
- `memcpy`/`memset` are byte loops (open, 2026-11-02). It matters for grid copies and resets.
- Struct fields can't be fixed arrays. There are no exponent literals. `#define` is capped at 16
  names.
- distlib's sidecar omits stdlib leaves (noted by kimiya, not yet filed). If bijli hits it,
  file it.
- **The 6.7.18 stdlib has no `log10`** (only `f64_ln` and `f64_log2`), so the `log10` rows stay
  pending until bijli's correctly rounded copy of kimiya's `_km_cr_log10` lands (M1). It belongs in
  G1's shared module.
- **New from bijli's language review (2026-11-11), not yet filed:**
  - a captured `f64` (local or parameter) is an untyped word inside a closure, so `a * 2.0` there
    is integer arithmetic and gives 0 (p6). The guide's 6.7.16 note covers a captured `f64` only
    for generic binding;
  - `#derive(Serialize)`'s pairs decoder `T_from_json(bayan_json_parse(..))` returns an empty
    `Vec<struct>` field, while the scanner reads it (p5);
  - (a note, maybe not a defect) bayan's JSON number notation differs from serde_json's/ryu's
    for exponents 16–20 and for 1e-6 (p5b).
- **Filed 2026-11-12 (M1 part 4):** `#derive(Serialize)` in a second `private` file fails to
  compile (p16b; cyrius `2026-11-12-derive-serialize-in-a-second-private-file-fails-to-compile.md`),
  and bayan's JSON tree turns `1e400` into +∞ without an error (bayan
  `2026-11-12-json-number-past-f64-range-reads-as-infinity-with-no-error.md`). bayan's other readings
  that bijli corrects (p16) were filed by sibling ports, and bijli's note is added to each: `-0` as
  the integer 0 (pravash), integers past i64 refused, raw control bytes and non-UTF-8 in strings
  (ushma), and no opt-in for repeated ignored keys. bijli works around each in src/json.cyr.

**Sibling notes**
- tanmatra: its `parity_core` per-function table overflows past 96 functions. kimiya noted
  this; bijli's sizing fix goes back to tanmatra as a note.

---

## 11. After 2.0.0

1. **The roadmap backlog**, now in Cyrius, with reference tests only:
   - CPML and materials for `GpuFdtd2d`;
   - TF/SF on the `Fdtd2d` contour;
   - frequency-aware lossy Fresnel;
   - sub-femtometre point-charge separations.
2. **A working mabda GPU backend**: `impl BjComputeBackend for BjMabdaBackend` behind the ported
   `BjWgpuBackend` scaffold's API, in a `[lib.gpu]` profile. It is new code in either language,
   so its parity is against `BjCpuBackend`, bit for bit where the device's f64 allows.
3. **Speed work with the ecosystem.** The benchmark columns and the §10 issues feed compiler,
   ganita and hisab tuning.
4. **Removing `rust-old/`** once the test ledger is closed and the corpus and generators are
   archived. Nothing waits on a consumer, because nothing is deferred.

---

## 12. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| P Prep | ✅ | | | 80 (Rust, 2026-11-12T03:15, under load) | D1–D5 answered; 1.2.0 tagged and green; P3, P4, P5 done (2026-11-11); P6 is done at the scaffold; the commit remains |
| M0 Initialize | ✅ | 34,746 (constants 12, math 34,734): 15,507 ok, 0 failed, 18,006 measured, 1,233 pending, on all three targets | | 80 Rust (pending twins) | done 2026-11-11: Parts 1–3. Baselines and probes in `port-parity.md` |
| M1 Foundation | ✅ | 9 modules, 189,015 ok, 11,841 accepted (ADR 01: 6,505; ADR 0002: 5,336), 0 measured | error.rs 10/10, logging.rs 1/1, field.rs 105/105, maxwell.rs 18/18, charge.rs 26/26 | 17 of 80 twinned; 11 Cyrius-only | Part 1 (the math) done 2026-11-11: 1a (crmath, corrected oracle, cr_verify), 1b (acos, log10, asin), 1c (tan). Part 2 done 2026-11-11: 2a (fmtf, log, error), 2c (complex, ln, digest v2), 2d (sinh, cosh, complex sin and cos). Part 3 done 2026-11-11: 3a–3c (field: vectors, point charges, dipoles; Gauss, multipoles, images, Green's, stress; trace_field_line), 3d (maxwell, charge). Part 4 done 2026-11-12: the JSON (D5, amended) of `BjComplex`, `BjFieldVector` and `BjPointCharge`, writers on bayan's printer and readers on bayan's tree with serde's semantics; ADR 0002 rows 10–12; `json.cr.tsv` (6,390 rows serde_json misreads). Part 5 (the range-edge review, §0.1, with the user's 2026-11-12 decisions) runs as 5a–5d. 5a done 2026-11-12: kimiya's correctly rounded monomials, square roots of quotients and unbounded-exponent arithmetic, copied (src/crmath_kx.cyr, tests/crmath_kx.tcyr) and verified against exact values (`cr_verify.py`'s exact families: 1,700,000 calls on x86_64 and aarch64, 0 mismatches). 5b done 2026-11-12: round once (ADR 0002 row 13) for the 14 single-product functions of maxwell, charge and field (2,031 mapped rows in `<module>.adr.tsv`), with bijli's correctly rounded √ of a monomial. 5c-1 done 2026-11-12: BjComplex's algebraic range edges (ADR 0002 row 14) and powf without its cutoff (row 15), 1,592 mapped rows. 5c-2 done 2026-11-12: the transcendental edges (ADR 0002 row 16: e^re, cosh and sinh, ln and pow of an |z| or a power past either end, sqrt's and powf's angles) on src/crmath_ux.cyr, verified against MPFR's widest exponent range; a 910-row edge grid; 1,987 mapped complex rows. 5d-1 done 2026-11-12: the exact Singularity threshold (ADR 0002 row 17) on src/crmath_cmp.cyr, a 3,826-row knife-edge grid, 302 mapped rows. 5d-2 done 2026-11-12: the vector chains' range edges (ADR 0002 row 18) on src/crmath_vx.cyr, 294 mapped rows. 5d-3 done 2026-11-12: the point charges, coulomb_force and the Green's functions (row 18, amended), 251 more. 5d-4 done 2026-11-12: the dipoles, the quadrupole, the axes and the sphere image, 251 more. 5d-5 done 2026-11-12: the stress tensor and its force, the Poynting flux, the multipole moments and the field line's step, 166 more. M1 done 2026-11-12 |
| M2 Media, waves, circuits, relativity | ◐ | 13 modules, 406,282 ok, 46,212 accepted (ADR 01: 8,998; ADR 0002: 37,214), 0 measured | material.rs 39/42 (3 n/a, ADR 0002 row 1), wave.rs 63/63, circuit.rs 51/51, relativity.rs 48/48 | 46 of 80 twinned | Part 1a done 2026-11-12: material.rs and `Material`'s JSON; ADR 0002 rows 19–21 (the Curie laws refuse NaN, `Material::wave_speed` validates as `impedance`, `classify_magnetic` returns a Result). Part 1b done 2026-11-12: wave.rs and Rust's `{x:.p}`; ADR 0002 rows 22–25 (the dipoles on bijli's η₀ with Cin(2π) exact, `effective_aperture` refuses NaN, Schlick's R₀ without its guard, Γ on the relative impedances). Part 1c done 2026-11-12: circuit.rs; ADR 0002 row 26 (a NaN area or plate separation refused). Part 1d-1 done 2026-11-12: relativity.rs without the callbacks; ADR 0002 row 27 (`lienard_wiechert_e` refuses a NaN R). Part 1d-2 done 2026-11-12: the retarded-time solver and `lienard_wiechert_fields` with trajectory closures over six scenarios; M2's four modules ported. Part 2 (the range-edge review) runs as 2a–2d: 2a-1 done 2026-11-12, ADR 0002 rows 28 (the half-wave pattern's axis) and 29 (the transients' −expm1); 2a-2 done 2026-11-12, ADR 0002 row 30 (a correctly rounded γ and γ², validity exact); 2b (2b-1–2b-4) done 2026-11-12, row 13 amended for material, wave, circuit and relativity; 2c-1 done 2026-11-12, row 17 amended for M2; 2c-2 done 2026-11-12, Clausius–Mossotti's tests exact and its ε_r rounded once; 2c-3 done 2026-11-12, the trig-free Fresnel TIR exact and its cos θ_t rounded once; 2c-4 done 2026-11-13, Liénard–Wiechert's κ < 10^-30 exact and its κ and R_eff rounded once; 2c-5 done 2026-11-13, the inverse combinations' |Σ 1/xᵢ| < 10^-30 exact and 1/Σ 1/xᵢ rounded once. **2c is done.** 2d-1 done 2026-11-13, the user's answers of 2026-11-13: the η-sum guard dropped, Liénard–Wiechert's n̂ − β rounded once; 2d-2a done 2026-11-13, wave's ratio chains at their range edges (row 18 amended); 2d-2b done 2026-11-13, the plane waves' and Snell's chains; 2d-2c done 2026-11-13, the material interface's roots and Γ and the waveguides' root; 2d-3a done 2026-11-13, circuit's series folds and impedance; 2d-3b done 2026-11-13, the transients' edges; 2d-3c done 2026-11-13, FourVector's chains, `boost_x`, the four-force and the invariants. Next: 2d-3d, the field transforms, the relativistic Larmor power and an ∞ R in E, then the solver's |
| M3 Optics | ☐ | | | | |
| M4 RF | ☐ | | | | |
| M5 FDTD core | ☐ | | | | |
| M6 FDTD 2D/3D | ☐ | | | | |
| M7 Integration surfaces | ☐ | | | | |
| M8 Release 2.0.0 | ☐ | | | | |
