# bijli Rust → Cyrius port plan

> **Plan of record** for porting bijli to Cyrius. Rust **1.3.0** will be the final Rust feature
> release and the **reference implementation** (the oracle); the Cyrius port ships as **2.0.0**.
> The method is tanmatra's ([its ADR 0006](../../../tanmatra/docs/adr/0006-cyrius-port-with-rust-oracle.md),
> [its port plan](../../../tanmatra/docs/development/tanmatra-port-plan.md)), as executed by impetus
> ([its port plan](../../../impetus/docs/development/impetus-port-plan.md), finished 2026-11-10) and kimiya
> ([its port plan](../../../kimiya/docs/development/kimiya-port-plan.md), in progress). Where this plan
> says "as impetus" or "as kimiya", follow that document. §2 lists what is different here and why.
>
> ▶ **Prep first (§4), then `cyrius port`.** Start on Cyrius 6.7.18, after 1.3.0 is committed and
> tagged (§4 P1). Nothing in `src/` changes after the tag.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | bijli **1.3.0** (to be cut: 1.2.0 + the `[Unreleased]` 3D CPML). 8,589 code lines in `src/` (no tests, no comments, no blank lines), 19 files, 430 `pub fn`s, 827 `#[test]`s, 80 criterion benches, 1 example | `Cargo.toml`, `VERSION`, a line count of `src/` |
| Rust health, 2026-11-11 on rustc 1.99.1 | ✅ `fmt --check`, clippy `-D warnings` (all features, all targets), rustdoc `-D warnings`, `cargo audit`, `cargo deny`. Tests: 830 unit, 26 integration and 2 doc tests pass in 6.2 s | local run |
| Correctness baseline | The 1.2.0 audit of every module (CHANGELOG 1.2.0, *Breaking (behavioral)* and *Fixed*), plus the 3D CPML measured at −80 dB | `CHANGELOG.md` |
| Cyrius toolchain | **6.7.18** (installed; still the newest) | `~/.cyrius/current` |
| hisab | Rust **1.4.0** gives bijli exactly one thing, `hisab::num::Complex`. Cyrius hisab **3.4.1** isn't needed in production (§2 item 1) | `grep hisab src/`, `Cargo.lock` |
| ganita | **1.2.15**, bundled in the 6.7.18 stdlib | `~/.cyrius/versions/6.7.18/lib/ganita.cyr` |
| Oracle's math | Rust **std** calls the host **glibc 2.44** for `sin` (27 calls), `cos` (21), `exp` (11), `sin_cos` (8), `atan2` (5), `log10` (5), `powf` (4), `hypot` (4), `atan` (4), `asin` (4), `acos` (3), `tan` (2) and `cbrt` (2). Complex values use Rust hisab 1.4.0's naive formulas (§8.2) | `grep` of `src/` outside tests |
| Cyrius consumers | **prakash 2.15.2** has primitive-value bijli bridges of its own and its own Gaussian-beam module. **tanmatra 2.0.x** has deferred its bijli bridge until bijli ports. Still Rust: soorat 1.0.0 (feature `em`) and ushma 1.4.0. No local repo: kiran, joshua, dravya | `~/Repos` |
| Target | bijli **2.0.0**, Cyrius | this document |

## 0. Decisions

### 0.1 Carried from the earlier ports (no new decision needed)

- **The oracle is frozen.** Rust 1.3.0 moves to `rust-old/` and is never edited. A generator
  crate in `scripts/oracle/` drives it by path, as tanmatra and impetus do. A wrong oracle value
  is fixed in Rust as 1.3.x and the corpus is regenerated.
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
- **No Rust 2.0.** The roadmap's "2.0 cleanup (breaking)" items that change the API's shape
  land in the Cyrius API only, accepted by ADR (§3.3). New physics stays out of the port.
- **Bit-identical across targets** for everything bijli computes: x86_64 ELF, x86_64 PE (wine)
  and aarch64 (qemu), as impetus.

### 0.2 New for bijli (proposed here; each is confirmed in ADR 0001 or 0002)

- **No production hisab dependency.** bijli uses only `hisab::num::Complex` from hisab. Cyrius
  hisab's `cx_*` doesn't reproduce it: it uses Smith division, a hypot-rescued `abs` and a C99
  `sqrt`, and it allocates 16 B on the never-freed bump heap for every result. bijli therefore
  owns a **by-value complex type**: a transcription of Rust hisab 1.4.0's `Complex`, operation
  for operation (§8.2). hisab may join at `scope = "test"`, as prakash and impetus have it, if a
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
  - That leaves **asin and tan** as the only new primitives. sinh appears only in tests.

  File the shared module upstream (§9 G1), so a fourth port can depend on it instead of copying.
- **The FDTD solvers are compared over whole trajectories** (§2 item 2). Their updates are
  linear and use only `+ − × ÷`, so a whole run is class B.

### 0.3 Decisions for the maintainer before §4 (recommendation first)

| # | Question | Recommendation | Alternative |
|---|---|---|---|
| D1 | Which Rust version is the oracle? | **Cut 1.3.0**: the `[Unreleased]` 3D CPML adds API (`try_enable_cpml*`, `disable_cpml`, `cpml`), so it is a minor release | Port from 1.2.0 and carry the 3D CPML as Cyrius-only work. That loses the oracle for the newest and hardest code, so it is not recommended |
| D2 | `bridge.rs` (8 fns, 18 tests) | **Port all 8 in M7.** They are primitive-value functions with no dependencies, and two partners (prakash, tanmatra) are already in Cyrius | Follow the sibling rule: port only the prakash and tanmatra bridges, and defer ushma and dravya |
| D3 | Time conventions: polarization uses e^{−iωt}, scattering uses e^{+iωt} | **Carry both as in Rust** in 2.0.0. Unifying them changes physics, not API shape, and the earlier ports kept physics out. Track the unification as a 2.x item: new convention-explicit functions, then deprecations | Unify in 2.0.0 (the one breaking release). The harness maps each row by `conj`, which is exact; prakash's conventions have to be checked first |
| D4 | `ai` (daimon/hoosh client, 117 lines) | **Port in M7 as a `[lib.ai]` profile**, on prakash's model: bayan for JSON and sandhi for HTTP, kept out of the default bundle | Defer until a consumer asks |
| D5 | serde (36 types derive it) | **Serde-compatible JSON for every derived type**, as [impetus ADR 0010](../../../impetus/docs/adr/0010-serde-compatible-json.md). That reversed impetus's own plan, and it lets Rust's serde round-trip tests port 1:1. Use hand-written codecs on bayan, because of the open `#derive(Serialize)` defects. Each type's codec lands in the milestone that ports the type | JSON only for the exchange types (`JonesVector`, `MieResult`, `SMatrix`, the `ai` configs), and the rest on demand |

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

## 2. What is different from the earlier ports

1. **Complex arithmetic is everywhere.** It appears in `rf` (106 mentions), `polarization` (46),
   `scattering` (30), `beam` and NF2FF (a private `Cplx` in `fdtd.rs`). The earlier ports had
   almost none.
   - A complex value is 16 bytes: a struct value passed in two registers, with operator
     overloads. Probe 2 confirms this on 6.7.18.
   - A value-form `Result` holds one word, so a complex or vector result is a code plus an
     out-pointer. A struct of 8 bytes or less in a payload hits the open small-struct bug.
   - impetus measured struct-value kernels at 2.5× slower than `f64` locals. So the API and
     the cold paths use struct values, while the hot loops (Mie recurrences, NF2FF accumulation,
     the FDTD kernels) work on `f64` locals.
   - Never write `s.f = g(t.e)`: the open field-store defect silently stores the argument instead
     of the result. Bind the result to a local first. impetus's `surface_parity.py` scans for
     the pattern.
   - Mie (`scattering`) runs complex recurrences (the downward log-derivative D_n and
     Riccati–Bessel ψ/χ). Its `a`/`b` coefficient vectors are public, so they are golden
     outputs too.
2. **FDTD is a linear, non-chaotic time stepper.** impetus couldn't compare long trajectories,
   because a 1-ulp change grows in a pile of bodies. A Yee update uses only `+ − × ÷` in a fixed
   order, and Rust performs no FMA contraction (bijli calls no `mul_add`). So the same operation
   order gives the same bits at every step.
   - **Trajectory rows (class T, §6)** compare the whole field state at chosen frames of a run
     that starts at t = 0. No mid-run state has to be loaded, so no oracle hook is needed for the
     private CPML Ψ arrays.
   - A frame's state is large: a 30³ grid holds 162,000 values in its six components. A row
     therefore stores a **digest** of each component's bit patterns, plus a few probe cells in
     full. A digest mismatch reruns the generator in `--dump` mode on that frame, to find the
     first differing cell.
3. **Grid memory is large and long-lived.**
   - The stdlib `alloc` is a bump allocator that never frees, and one test file builds hundreds
     of grids. Grids live in freelist- or mmap-backed storage and are released by
     `Fdtd1d_free`, `Fdtd2d_free` and `Fdtd3d_free`.
   - CPML Ψ slabs and dispersive-material ADE arrays are allocated only when they are enabled,
     as in Rust.
   - **Never `alloc_reset`**: it would wipe the consumer's heap. This is M0 probe 1.
4. **Functions that take closures.**
   - Public: `trace_field_line<F>` (field.rs), `solve_retarded_time<F>` and
     `lienard_wiechert_fields<F>` (relativity.rs) take a field or trajectory function.
   - Private: `nf2ff_tm_core` takes three `impl Fn` arguments.
   - In Cyrius each becomes a function pointer plus a context word, called through `fncall*`.
   - Rust reports field-function errors through `Result`, and `Result` can't go through
     `callptr`, so the callback writes a status word instead (probe 4).
5. **Text I/O.** `rf` reads and writes Touchstone files.
   - Parsing: `f64_parse_n` is correctly rounded, so it matches Rust's `str::parse::<f64>`.
   - Writing: Rust prints scientific notation with `{:e}`, the shortest round-trip digits. Reuse
     impetus's `src/fmtf.cyr` and add the exponent form.
   - Malformed-file rejection becomes golden rows too.
6. **A trait-generic compute backend.** `GpuFdtd2d<B: ComputeBackend>` becomes a struct of
   function pointers (a vtable) plus a backend context. `CpuBackend` is ported.
   `WgpuBackend` is a scaffold in Rust and isn't carried. A real GPU backend comes after 2.0, on
   **mabda** (the Cyrius GPU layer, 4.2.0, already in the stdlib).
7. **Error payloads.** Several `BijliError` variants carry data: `DimensionMismatch { expected,
   got }`, `InsufficientResolution { cells, domain_size }`, `InvalidPermittivity { value }` and
   the `String` reasons. They become append-only `BJ_ERR_*` codes with no detail string, as
   tanmatra. The reason text goes to a sakshi warning at the failure site.
8. **More transcendentals.** asin and tan aren't in either sibling's correctly rounded set
   (§0.2).
9. **A wide surface in one profile.** bijli has 430 `pub fn`s, all in the default bundle.
   tanmatra's `parity_core` per-function table overflows past 96 functions, so size it for 512
   from the start (M0 step 7).

---

## 3. Scope

### 3.1 In scope for 2.0.0

| Milestone | Rust module | Code lines | `pub fn` | Tests | Cyrius module(s) |
|---|---|---|---|---|---|
| M1 | `error.rs` | 21 | — | 10 | `src/error.cyr` |
| M1 | (constants in `field.rs`, `fdtd.rs`) | — | — | — | `src/constants.cyr` (generated and proven, as impetus) |
| M1 | (math) | — | — | — | `src/mathx.cyr`, `src/crmath.cyr`, `src/crmath_data.cyr` |
| M1 | `polarization.rs` `Complex` helpers and the hisab 1.4.0 `Complex` | — | — | — | `src/complex.cyr` (`BjComplex` by value) |
| M1 | `field.rs` | 786 | 38 | 105 | `src/field.cyr` (`BjFieldVector` by value, point fields, dipoles, multipoles, images, Green's functions, stress tensor, field lines) |
| M1 | `maxwell.rs` | 108 | 10 | 18 | `src/maxwell.cyr` |
| M1 | `charge.rs` | 172 | 14 | 26 | `src/charge.cyr` |
| M2 | `material.rs` | 283 | 33 | 42 | `src/material.cyr` |
| M2 | `wave.rs` | 464 | 42 | 63 | `src/wave.cyr` |
| M2 | `circuit.rs` | 297 | 27 | 51 | `src/circuit.cyr` |
| M2 | `relativity.rs` | 517 | 27 | 48 | `src/relativity.cyr` |
| M3 | `polarization.rs` | 565 | 53 | 69 | `src/polarization.cyr` |
| M3 | `scattering.rs` | 250 | 6 | 32 | `src/scattering.cyr` |
| M3 | `beam.rs` | 411 | 28 | 60 | `src/beam.cyr` |
| M4 | `rf.rs` | 984 | 31 | 80 | `src/rf.cyr`, `src/touchstone.cyr`, `src/fmtf.cyr` |
| M5 | `fdtd.rs`: 1D, pulses, TF/SF, lossy and dispersive materials, NF2FF, subpixel, non-uniform grid, 4th order, CW, eigenmodes | 1,209 | ≈ 45 | ≈ 80 | `src/fdtd1d.cyr`, `src/fdtd_source.cyr`, `src/fdtd_media.cyr`, `src/nf2ff.cyr`, `src/fdtd_util.cyr` |
| M6 | `fdtd.rs`: CPML 2D and its helpers, `Fdtd2d`, CPML 3D, `Fdtd3d` | 1,691 | ≈ 47 | ≈ 85 | `src/grid.cyr` (storage), `src/cpml.cyr`, `src/fdtd2d.cyr`, `src/fdtd3d.cyr` |
| M6 | `compute.rs` (`ComputeBackend`, `CpuBackend`, `GpuFdtd2d`) | 313 | 9 | 18 | `src/compute.cyr` |
| M7 | `bridge.rs` (D2) | 78 | 8 | 18 | `src/bridge.cyr` |
| M7 | `logging.rs` → sakshi | 14 | 2 | 1 | `src/log.cyr` (level setup, the `_bj_warn` helpers that every module uses from M1) |
| M7 | `ai.rs` (D4) | 117 | 2 | 7 | `src/ai.cyr`, in the `[lib.ai]` profile only |

fdtd.rs's test count is split by section. Its `cpml_tests` and `cpml3d_tests` modules go to M6.
Surface parity (§7 gates) covers every `pub fn` and every `pub` method.

### 3.2 Deferred or not carried

| Rust surface | Decision | Why |
|---|---|---|
| `integration/soorat.rs` (feature `soorat-compat`, 8 structs, 14 tests) | **deferred** until soorat ports | soorat is Rust 1.0.0 (its `em` feature). This is the same rule as tanmatra's and kimiya's |
| `compute::WgpuBackend` (feature `compute-wgpu`) | **not carried** | a scaffold in Rust. A mabda backend is a 2.x item (§10) |
| serde derives (36 types) | **serde-compatible JSON per D5** (recommended: every derived type, hand-written on bayan, as impetus ADR 0010) | the derive-serialize defects are open (bayan S3/S4); the soorat structs wait with soorat |
| the three 1.2.0 deprecations: `rf::input_impedance_lossless`, `rf::write_touchstone`, `Material::reflectance_at` | **not carried** (§3.3 #1) | their replacements are ported. Their Rust tests are n/a in the ledger |
| `ComplexExt` trait | its three methods become plain functions on `BjComplex` | Cyrius has no traits on foreign types |
| `Send`/`Sync`, `#[inline]`, `const ADD: bool` generics | **not carried** | Rust-only mechanics. Each `row_diff`/`row_curl` instance becomes its own function |

### 3.3 API changes in the Cyrius version (accepted divergences, ADR 0002)

These come from the roadmap's "2.0 cleanup (breaking)". Each row says how the harness still
compares it, so no function leaves parity entirely.

| # | Rust 1.3.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | three deprecated functions | removed | through their replacements |
| 2 | public FDTD grid fields (`ex`, `ez`, `nx`, …), with a per-step layout check because callers can resize them | buffers private behind accessors (pointer and length). Nothing can resize them, so the per-step layout check and its `tracing::error` path go | class T on the accessors. The Rust layout-mismatch tests are n/a |
| 3 | `TfsfSource::new` and `Cpml::new` sanitize with a warning; `Fdtd2d::enable_cpml` warns and ignores | fallible only: `BjTfsfSource_new` and the CPML constructors return `Err`, and `enable_cpml` is the `try_` form | oracle sanitize-or-ignore rows ↔ `Err`, mapped by name |
| 4 | `BijliError` variants carry data and `String`s | `BJ_ERR_*` codes, append-only, with no payload | class A on the code. The payload is logged |
| 5 | the polarization/scattering time-convention split | **as Rust** (D3) | — |
| 6 | results bit-identical only on the same platform | **bit-identical across targets** | measured on three targets, as impetus |

---

## 4. Prep: everything before `cyrius port`

These steps run on the Rust tree. They end in a Rust-side **"port plan" commit**, which is what
tanmatra (`2febf11`), impetus (`a26296b`) and kimiya (`6281f19`) each landed before their
initial-port commit.

- **P0. Read this plan and settle D1–D5** (§0.3). Record the answers in ADR 0001 / 0002 drafts.
- **P1. Cut the oracle release.** *Maintainer action.*
  - Move `[Unreleased]` to `[1.3.0] — <date>`. Bump `VERSION` and `Cargo.toml` with
    `scripts/version-bump.sh`.
  - Run the full P(-1) gate on that commit: `cargo fmt --check`; clippy
    `--all-features --all-targets -D warnings`; `cargo test --all-features`;
    `RUSTDOCFLAGS="-D warnings" cargo doc`; `cargo audit`; `cargo deny check`; and
    **`make features`**, which builds every feature on its own (1.2.0's per-feature isolation).
  - Commit and tag `1.3.0`. The release workflow requires the tag to equal `VERSION`.
  - **Also check that the 1.2.0 tag exists on GitHub.** The local clone has `0.24.3`, `1.0.0`,
    `1.0.1` and `1.1.0`, but not `1.2.0`.
- **P2. Decide on oracle hooks: none planned.** impetus needed a hidden `oracle` feature to see
  state its API hid. bijli's API already shows almost everything:
  - the FDTD field arrays are public;
  - Mie's `a`/`b` coefficients are public;
  - class T starts every run at t = 0 (§2 item 2), so no mid-run state has to be loaded.

  Add a `#[doc(hidden)]` `oracle` feature in 1.3.1 only if an M0 probe shows a need. The likely
  candidates are the CPML σ/κ/α profile coefficients (`cpml_coefficients`) and the private Mur
  coefficient.
- **P3. Record the oracle's environment** in the plan's header table and in ADR 0001:
  - rustc 1.99.1, glibc 2.44 (Rust std's libm), Rust hisab 1.4.0 (`Complex`);
  - the `Cargo.lock` at the tag.
  - **Which calls actually reach glibc** in the generator's build profile. Release LLVM rewrites
    some calls: `2.0.powf(y)` becomes `exp2`, as kimiya's `port-parity.md` records, and
    constant arguments fold at compile time. List each `powf`/`powi`/`exp` site with what the
    optimized oracle executes.
- **P4. Take a quiet Rust bench baseline at the tag.** Run `./scripts/bench-history.sh` on an
  idle machine; the 2026-11-10T04:57 row ran under load. That run is the Rust column of
  `docs/benchmarks-rust-v-cyrius.md`: 80 benches, including `fdtd/3d_step_30x30x30_cpml8`.
- **P5. Tidy what will move.**
  - `target/` is 5.5 GB, ignored and untracked. `cyrius port` moves it into `rust-old/target/`.
    Either `cargo clean` first, or accept the move and ignore `/rust-old/target/`.
  - No tracked build output, unlike impetus's `fuzz/target/`.
  - `.claude/` is ignored. Nothing else needs untracking.
- **P6. Pre-stage what the scaffold will overwrite.** `cyrius port` 6.7.18 still overwrites
  `CLAUDE.md`, `docs/development/roadmap.md`, `docs/development/state.md` and both workflows. It
  also pins `version = "0.1.0"`, leaves `bench` out of `[test] stdlib`, and its stubs fail
  `cyrius audit` (the open cyrius issue in §1). Everything is in git, so restore from `HEAD`
  after the scaffold. Draft the merged CLAUDE.md now (§5 step 2) so the repair is a paste.
- **P7. Draft ADRs 0001 and 0002.**
  - bijli has no `docs/adr/` yet, so start directly in the scaffold's four-digit form.
  - [0001 — Port to Cyrius with Rust 1.3.0 frozen as the oracle](../adr/0001-cyrius-port-with-rust-oracle.md)
    is drafted with this plan.
  - 0002, the §3.3 API divergences, is written at M0 step 2.
- **P8. Check upstream readiness.**
  - Is a Cyrius newer than 6.7.18 out? If so, read its `ecosystem-migration.md` entries and
    decide the pin.
  - Re-read the open issues in §1. A fixed small-struct `Result` bug or scoped allocator (S1)
    changes §8.3 and §8.5.
- **P9. Land the "port plan" commit:** this file, the ADR 0001 draft, and the roadmap's pointer
  to this plan. No `src/` change.

**Gate for P:** 1.3.0 tagged and green; D1–D5 answered; the bench baseline recorded; the port
plan commit landed.

---

## 5. M0 — Initialize with `cyrius port`

1. **Scaffold.** Run `cyrius port --dry-run .`, then `cyrius port .`.
   - The dry run reports 22,696 Rust lines.
   - Expect these to move to `rust-old/`: `Cargo.toml`, `Cargo.lock`, `src/`, `tests/`,
     `benches/`, `examples/`, `Makefile`, `deny.toml`, `codecov.yml`, `rust-toolchain.toml`
     and `target/`.
   - The dry run doesn't list them all (cyrius issue item 5), so check with `git status` after.
2. **Repair what the scaffold overwrites** (as impetus §4 step 2):
   - **CLAUDE.md**: merge, keeping the identity, process, DO NOT list and documentation
     structure. Translate the principles as impetus did:
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
     1.3.0 → 2.0.0-dev.
   - **Workflows**: keep the scaffold's `release.yml`, which is identical in all three siblings.
     Take impetus's `ci.yml`, with `test --bounds`, `deps --locked`, `audit`, the ledger,
     bench-twin and doc-example checks on from the start, plus an **`[lib.ai]` profile job**
     once M7 lands. Move the Rust workflows to `rust-old/.github/workflows/`. Append
     `/rust-old/target/` (not `/rust-old/`) and `/build/` to `.gitignore`.
   - **Move Rust-only material** to `rust-old/`: `scripts/bench-history.sh`,
     `scripts/version-bump.sh`, `bench-history.csv`, `benchmarks.md` and
     `docs/architecture/overview.md`. Rewrite the overview for Cyrius in M8.
   - **ADRs**: 0001 (the port) and 0002 (§3.3), with `docs/adr/README.md` and `template.md` from
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
   stdlib = ["assert", "bench"]

   [lib]            # the default bundle: everything except ai
   modules = ["src/error.cyr", "src/constants.cyr", "src/mathx.cyr", "src/crmath_data.cyr",
              "src/crmath.cyr", "src/complex.cyr", "src/log.cyr", "src/field.cyr", ...]

   [lib.ai]         # the default modules plus src/ai.cyr (bayan, sandhi): dist/bijli-ai.cyr
   modules = [...]
   ```
   - **Module order** follows Rust's feature graph:
     - error → constants → mathx → crmath → complex → log;
     - then field → maxwell → charge;
     - then material → wave → circuit → relativity;
     - then polarization → scattering → beam → rf;
     - then the fdtd files → compute → bridge.
   - A profile repeats the whole list, and the order is load-bearing (prakash's comment in its
     `cyrius.cyml`).
   - No `[deps.hisab]` (§0.2).
4. **Toolchain and math check.** Confirm the Cyrius counterpart of every Rust f64 method bijli
   calls (§8.1). Port impetus's `mathx.cyr` helpers as `_bj_*`: Rust-semantics
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
   - **`[lib]` is never prepended to tests.** Mirror the module list in `tests/lib_bijli.cyr`,
     with a drift check against `cyrius.cyml`, as impetus.
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
   - New correctly rounded primitives (atan, asin, tan) are gated at 0 ulp against MPFR on
     ≥ 10⁶ pairs, as impetus measured its pow and atan2.
9. **Probes**, each answered with a small program and numbers before M1 writes code:
   1. **Grid memory.** Build and free 1,000 `Fdtd3d` grids of 30³ in a loop, with and without
      CPML and dispersive materials. Gate: the heap high-water mark is flat after the first
      grid. Also time freelist against mmap-backed storage for a 53³ grid.
   2. **Complex by value.**
      - Returns of 16 bytes in two registers; operator overloads (`BjComplex_add`, …);
      - complex values as struct fields (`BjJonesVector` holds two);
      - arrays of complex values as a stride buffer (Mie's `a`/`b`, `SMatrix`);
      - the field-store defect (`2026-11-09-field-store-of-struct-call-copies-a-field-argument.md`),
        which impetus had to work around.
   3. **FieldVector by value.** 24-byte returns through the hidden pointer, operator overloads
      for `+ − × neg`, and the compound-assignment forms.
   4. **Callbacks.** A field function through `fncall*` with a context word and a status
      out-pointer. Time it in `trace_field_line`'s loop of 10⁴ steps.
   5. **The FDTD inner loop.** `row_diff`/`row_curl` written four ways: indexed loads, pointer
      stepping, `f64v_*` (mind the scalar-tail bug) and unrolled. Gate: bits equal to the Rust
      kernel. The fastest form is recorded for M6.
   6. **Rust iterator semantics.** 13 `.sum()` sites: Rust's float `Sum` starts from −0.0, so
      check the empty and all-negative-zero cases. Also `.min_by`/`partial_cmp` tie order (the
      eigenmode peak search) and `chunks`/`windows` bounds.
   7. **Touchstone text.** `f64_parse_n` against Rust's `parse::<f64>` on the 1.2.0 test
      corpus; `{:e}` formatting through `fmtf.cyr`, on every value in the 201-point bench file.
   8. **Test runtime.** 3D CPML tests are costly: a 53³ reference run for 300 steps is
      ≈ 45M cell-steps. Measure the Cyrius cost per cell-step, and size the class-T scenarios
      to fit `cyrius test` in CI. Bigger scenarios become class R in a separate slow file.

**Gate:** `cyrius build` and `cyrius bench` pass; `cyrius test` runs the harness (0 modules)
with no warnings; `cyrius lint` and `fmt --check` are clean; `cyrius audit` exits 0. Committed:
the generator, `tests/ulp.cyr`, ADRs 0001 and 0002, the probe results in `port-parity.md` and the
upstream issues (§9).

---

## 6. Parity policy

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
  - an oracle defect, which is fixed in Rust 1.3.x and regenerated.

  The only exception is a divergence accepted by ADR, and the harness still reports it by name.
- **A class T digest mismatch is never "close enough".** Bisect to the first frame and cell that
  differ. Then find the operation order that differs (§8.3 FDTD kernels).
- **Branch points on both sides:**
  - Fresnel total internal reflection and the `cos θ_i` 1e-6 clamp;
  - the Mie limits (`MIE_MAX_SIZE_PARAMETER`, `MIE_MAX_MX`) and the Im(m) > 0 rejection;
  - the r² < 1e-30 m² singular cutoff;
  - the CPML fit check `2·t + 3 ≤ min(n)`;
  - `dt` at and just above the Courant limit;
  - the matched-load stub `(0.0, 0.25)`;
  - the damping regimes in `circuit`;
  - β → 1 in `relativity`.

---

## 7. Milestones

Every milestone, in Parts where it is large, ends with impetus's gates:
- `cyrius build`; `cyrius test` and `test --bounds` **with no warnings**;
- `cyrius bench`, recorded against the Rust numbers;
- `vet`, `lint` and `fmt --check` on the milestone's files;
- `distlib --all`, then `distlib --check`, then `deps --locked`;
- **surface parity**: every in-scope Rust `pub fn` has a Cyrius counterpart or a written reason;
- the grid soak (§5 step 9.1) from M5 on;
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
| **M1 Foundation** | error, constants, mathx, crmath, complex, log, field, maxwell, charge | `BijliError` → `BJ_ERR_*`. Rust hisab 1.4.0's `Complex` transcribed by value (§8.2), with `complex_real`/`zero`/`from_polar` and `norm`/`norm_sq`/`arg`. The correctly rounded primitives, including the new atan/asin/tan. `BjFieldVector` with operators. Point, line, ring, disk and sphere fields; dipoles; superposition; the `associated_legendre` multipoles; images; Green's functions; the stress tensor; `trace_field_line` with callbacks. Maxwell and charge functions | every function over seeded grids including 0, ±tiny, ±huge, NaN and ±∞; positions near the singular cutoff; multipoles l ≤ 10; field lines for point, dipole and quadrupole fields, with singular and escaping endings |
| **M2 Media, waves, circuits, relativity** | material, wave, circuit, relativity | `MagneticType` as an enum. Clausius–Mossotti, Curie–Weiss and the Fresnel family (the trig-free `*_direct` forms, Schlick, TIR). RC/RL/RLC with the damping regimes. Lorentz factor, four-vectors, the EM tensor, boosts, Liénard–Wiechert, and `solve_retarded_time` with callbacks | each Fresnel function over θ ∈ [0, π/2] × n₁/n₂ including TIR and the clamp edge; the circuit regimes; β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9}; retarded time for the five scenario trajectories |
| **M3 Optics** | polarization, scattering, beam | Jones vectors and matrices, Stokes, Mueller (the e^{−iωt} convention, D3). Mie with its downward D_n recurrence and complex ψ/χ (e^{+iωt}); its `a`/`b` vectors are golden too. Rayleigh. Beams: q-parameter, ABCD, resonator stability, `hermite_polynomial`/`laguerre_polynomial`, HG/LG fields | Mie over x ∈ {1e-3 … 1e5} × m including high index, absorbing and the rejection limits; every Jones element over seeded angles; HG/LG up to order 10 at z ∈ {0, ±z_R, 10 z_R} |
| **M4 RF** | rf, touchstone, fmtf | Transmission lines, VSWR, Smith chart, S-parameters and cascades, stub matching, the array factor, `input_impedance_complex`. Touchstone read and write in the spec's `N11 N21 N12 N22` order: the rejection of malformed files and the scientific output, byte for byte | the 1.2.0 Touchstone tests as files; seeded round trips; every malformed case; stub matching including the matched load |
| **M5 FDTD core** | fdtd1d, fdtd_source, fdtd_media, nf2ff, fdtd_util | `Fdtd1d` with the per-end Mur coefficient. Gaussian, Ricker, differentiated Gaussian and modulated pulses. `TfsfSource` (with Mur at both ends; fallible, §3.3 #3). Lossy, Drude, Lorentz and Debye ADE. Subpixel averages. NF2FF (J and M, complex accumulation) on `BjComplex`. `NonUniformGrid`, `fourth_order_derivative`, `CwAccumulator`, `find_eigenmodes` (DFT and half-power Q) | class T for 1D: vacuum, ε_r = 4, lossy, each dispersive model, at frames {1, 10, 100, 1000}; pulses over seeded t; NF2FF on 80 × 80 × 180 angles; eigenmodes on 4096 × 256 |
| **M6 FDTD 2D/3D** | grid, cpml, fdtd2d, fdtd3d, compute | The grid storage from probe 1 and the kernels from probe 5. Roden–Gedney CPML in 2D (strips) and 3D (12 Ψ slabs, κ, b/c at the Yee positions) as post-loop corrections. TM/TE, conductivity C_a/C_b, `set_material` validation before any write, `total_energy`. `ComputeBackend` as a vtable, `CpuBackend`, `GpuFdtd2d` | class T: 2D TM and TE at 100 × 100, PEC/PMC, CPML 10 (classic, CFS), lossy; 3D at 30³, PEC and CPML 8; `GpuFdtd2d` against `Fdtd2d`. Class R: the −88 dB and −80 dB reflection tests and the energy decay to 2e-8, sized by probe 8 |
| **M7 Integration surfaces** | bridge, log, ai | The 8 bridges (D2). Sakshi level setup (`init`, `init_with_level`). The daimon/hoosh client in `[lib.ai]` (D4), on prakash's `ai.cyr`: config JSON, a redacting debug print, client construction | the bridge functions over seeded grids; the config JSON round trips; the `[lib.ai]` profile builds and its tests run in CI |
| **M8 Release 2.0.0** | all | Close the Rust-test ledger, which has run since M0: each of the 830 unit, 26 integration and 2 doc tests is ported, folded into golden rows, or marked n/a or deferred. The full surface-parity diff. `docs/benchmarks-rust-v-cyrius.md` with all 80 benches in both columns. `docs/guides/migration-2.0.md`: Rust idioms, the ADR 0002 changes, and `BjComplex` replacing hisab's `Complex`. README, overview, getting-started, SECURITY. `cyrius distlib --all`, `cyrius audit`, VERSION 2.0.0 | the full corpus green, with equal digests on x86_64 ELF, x86_64 PE and aarch64 |

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

## 8. Rust → Cyrius mapping (bijli specifics)

### 8.1 Math functions

Every std f64 method in `src/` outside tests, at 1.2.0 + `[Unreleased]`:

| Rust (calls) | Cyrius | Parity |
|---|---|---|
| `sqrt` (72), `abs` (81), `ceil` (2) | `f64_sqrt`, `f64_abs`, `f64_ceil` (builtins) | exact on every target |
| `sin` (27), `cos` (21), `sin_cos` (8) | `_bj_cr_sin`, `_bj_cr_cos`, `_bj_cr_sin_cos` (from impetus `crmath.cyr`) | class C, 0 ulp against the corrected oracle |
| `exp` (11) | `_bj_cr_exp` (impetus) | class C |
| `atan2` (5) | `_bj_cr_atan2` (impetus) | class C |
| `powf` (4) | `_bj_cr_pow` (impetus) | class C |
| `hypot` (4) | `_bj_cr_hypot` (impetus) | class C |
| `log10` (5) | `_bj_cr_log10` (from kimiya) | class C |
| `acos` (3) | `_bj_cr_acos` (from kimiya) | class C |
| `atan` (4) | `_bj_cr_atan2(y, 1.0)` | class C; exact, since both are the rounding of the same real |
| `asin` (4), `tan` (2) | **new**: `_bj_cr_asin`, `_bj_cr_tan`, Ziv two-phase as kimiya's | class C; gated on ≥ 10⁶ MPFR pairs at M0 step 8 |
| `cbrt` (2) | `ganita_f64_cbrt` (correctly rounded) | class C |
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

### 8.2 Complex: Rust hisab 1.4.0 → `BjComplex`

`src/complex.cyr` transcribes `hisab-1.4.0/src/num/complex.rs`, operation for operation, so
that class B/C rows can be compared:

| Rust hisab 1.4.0 | Formula (keep exactly) | Note |
|---|---|---|
| `abs` / bijli `norm` | `sqrt(re·re + im·im)` | no hypot rescue; Cyrius hisab's `cx_abs` has one, so it would differ. In range: Rust's bits. Beyond about 2^±511: a power-of-two prescale, as an ADR 0002 row (§0.1) |
| `norm_sq`, `conj`, `+ − ×`, `× f64`, `neg` | the textbook forms, Rust's order | class B |
| `/` | `((a·c + b·d)/den, (b·c − a·d)/den)`, `den = c² + d²` | naive. bijli's own Smith `cdiv`/`crecip` (rf) and `complex_div` (scattering) are ported as written too |
| `/ f64`, `inv` | componentwise | class B |
| `arg` | `atan2(im, re)` | class C |
| `from_polar(r, θ)` | `(r·cos θ, r·sin θ)` | class C |
| `exp` | `e^re · (cos im, sin im)` | class C |
| `sqrt` | `from_polar(sqrt(abs), arg/2)` | class C. Not the C99 branch of Cyrius hisab's `cx_sqrt` |
| `ln`, `powf`, `sin`, `cos`, `is_zero` | as written | ported only if a non-test call exists. Today there is none, except `sqrt` (rf.rs:145) |

### 8.3 Types and idioms

| Rust | Cyrius |
|---|---|
| `FieldVector` (`Add`, `Sub`, `Mul<f64>`, `Neg`, `*Assign`, `f64 * v`) | `BjFieldVector` struct value (24 B) with operator overloads and explicit `_scale_left` |
| `Complex` | `BjComplex` struct value (16 B, two registers), §8.2 |
| `JonesVector`, `JonesMatrix`, `StokesVector`, `MuellerMatrix` | struct values. Matrices are flat fields (`m00 … m33`) because struct fields can't be fixed arrays |
| `Vec<f64>` grid fields | typed stride buffers from the grid storage (probe 1) with a length word; accessors return pointer + length |
| `Vec<Complex>` (Mie `a`/`b`), `Vec<Vec<Complex>>` (`SMatrix`) | one stride buffer of `BjComplex` with n (or n × n × freq) and an index helper |
| `Result<T, BijliError>` | value-form `Result` + `BJ_ERR_*` code, with `?`. For a payload over one word (complex, field vector), use the code plus an out-pointer (the small-struct bug). No `Result` through `callptr` |
| `Option<f64>` (`half_power_offset`), tuples (`(f64, f64)` stub match, CPML coefficients) | an out-pointer + 1/0, or several out-pointers |
| `impl Fn` / generic `F: Fn` | a function pointer plus a context word, called with `fncall*` (probe 4) |
| `ComputeBackend` trait + `GpuFdtd2d<B>` | `BjComputeBackend` struct of function pointers + a context; `BjCpuBackend_new()` |
| `row_diff::<ADD>` / `row_curl::<ADD>` | two functions each (`_add`, `_sub`), in probe 5's fastest bit-exact form |
| `try_filled` (`try_reserve`) | the grid allocator returns 0 → `BJ_ERR_ALLOCATION` (a new code, appended) |
| `checked_mul` in `grid_cells` | `overflow.cyr` |
| `tracing::{warn, error, debug}` (67 sites) | `_bj_warn` / `_bj_error` / `_bj_debug` over sakshi, behind a level |
| iterator chains, `.sum()`, `.zip()`, `.windows()` | loops in Rust's order. `.sum()` starts from −0.0 (probe 6) |
| `#[non_exhaustive]` / `#[must_use]` | append-only enum values / `#must_use` |

**Typing traps beyond §8.1** (from the siblings' `port-parity.md`):
- Multi-return values are untyped words: `-x` on one is integer negation, with no diagnostic.
  Returns of four values don't parse.
- A `var u = x` copy of an f64 is f64-typed, which can turn a NaN self-compare into a trap.
- `&c.field` is refused.
- `#inline` takes at most two parameters.
- `f64_min`/`f64_max` are NaN-asymmetric, `f64_round` ties to even, and `f64_sign` isn't Rust's
  `signum`.
- `vec_sort_by` is not stable, so sort on (key, index). bijli has one sort site, the eigenmode
  peak search.

### 8.4 Names

- **Public structs** `Bj<Name>` with `impl` blocks, so methods keep their Rust names:
  `g.step_tm()`, `g.try_enable_cpml(t)`. Constructors are `BjFdtd2d_new`, `BjPointCharge_new` and
  so on. Free functions keep the Rust name with the `bj_` prefix: `bj_electric_field`,
  `bj_mie`, `bj_fresnel_rs`.
- **Globals and enum members** carry `BJ_`; private helpers carry `_bj_`.
- The prefix scan is clean today: no `bj_`, `BJ_`, `_bj_` or `Bj[A-Z]` in the 6.7.18 stdlib,
  hisab, prakash or tanmatra bundles. Re-scan at M1 for unprefixed collisions such as `vec_*`,
  `Ok`/`Err` and ganita's `mat_*`. "Last definition wins" across files is only a warning.

### 8.5 Memory

Settled by probe 1. The starting design:
- **Grid lifetime**: field, coefficient, Ψ and ADE arrays live in storage the grid owns
  (freelist, or mmap for large grids), released by `*_free`.
- **Call lifetime**: NF2FF accumulators, eigenmode spectra and Touchstone buffers are allocated
  per call, and their caller frees them. Field-line polylines go into a caller-owned growable
  buffer.
- **No per-step allocation** in any `step_*`: the kernels work in place, as Rust's do.
- **Never `alloc_reset`.**

---

## 9. Gaps to file upstream (M0, then as found)

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

**Cyrius stdlib / compiler** (re-check each at P8; add any new ones the probes find)
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

**Sibling notes**
- tanmatra: its `parity_core` per-function table overflows past 96 functions. kimiya noted
  this; bijli's sizing fix goes back to tanmatra as a note.

---

## 10. After 2.0.0

1. **soorat integration** (`src/soorat.cyr`), once soorat ports. Parity comes from the Rust
   tests in `rust-old/`.
2. **The roadmap backlog**, now in Cyrius, with reference tests only:
   - CPML and materials for `GpuFdtd2d`;
   - TF/SF on the `Fdtd2d` contour;
   - frequency-aware lossy Fresnel;
   - sub-femtometre point-charge separations.
3. **A mabda GPU backend** for `BjComputeBackend` (replacing the never-real `WgpuBackend`),
   in a `[lib.gpu]` profile.
4. **Unifying the time conventions** (D3), if it is not taken in 2.0.0.
5. **Speed work with the ecosystem.** The benchmark columns and the §9 issues feed compiler,
   ganita and hisab tuning.
6. **Removing `rust-old/`** after soorat is ported or dropped, the corpus and generators are
   archived, and the test ledger is closed.

---

## 11. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| P Prep | ☐ | | | | D1–D5; 1.3.0 tag |
| M0 Initialize | ☐ | | | | |
| M1 Foundation | ☐ | | | | |
| M2 Media, waves, circuits, relativity | ☐ | | | | |
| M3 Optics | ☐ | | | | |
| M4 RF | ☐ | | | | |
| M5 FDTD core | ☐ | | | | |
| M6 FDTD 2D/3D | ☐ | | | | |
| M7 Integration surfaces | ☐ | | | | |
| M8 Release 2.0.0 | ☐ | | | | |
