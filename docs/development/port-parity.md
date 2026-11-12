# bijli port parity — measured baselines

Measurements that set the parity classes and the M1+ implementation choices (port plan §7). The oracle is Rust
bijli 1.2.0 (tag `1.2.0`, `c95fd30`) on **glibc 2.44**, x86_64 with FMA: Rust std's f64 methods call the host
libm, except `cbrt` (compiler-builtins) and `powi` (compiler-rt). The port is Cyrius 6.7.18 with ganita 1.2.15
and no hisab. References: MPFR 4.2.2 for correctly rounded values, Python's IEEE doubles for the exact functions.

## 1. Semantics (M0 steps 4–5, 2026-11-11, x86_64 Linux ELF)

The language itself is reviewed in the port plan's §2, with the probes `scripts/probes/p1`–`p6`. What the
corpus adds, per function:

| Rust (std) | Cyrius | Result | Consequence |
|---|---|---|---|
| `f64::min` / `max` ignore a NaN operand on either side | `f64_min` / `f64_max` | differ from Rust on **69 of 700** rows each (a NaN second operand gives NaN) | `_bj_min` / `_bj_max` / `_bj_clamp` (`src/mathx.cyr`): 700 / 700 / 870 rows exact |
| `is_nan`, `is_finite`, `copysign` | — (no public stdlib predicate) | `x != x`; `x − x == 0`; a sign-bit merge | `src/mathx.cyr`, exact |
| `to_radians`, `to_degrees` | — | one multiply by Rust's folded π/180 and 180/π | exact (781 and 717 rows) |
| `powi` (run-time exponent: compiler-rt `__powidf2`) | — | square-and-multiply, then 1/x for n < 0 | `_bj_powi`, exact on all 7,800 rows (exponents −12 … 4,095 and ±i32 limits) |
| `cbrt` (compiler-builtins: CORE-MATH, correctly rounded) | `ganita_f64_cbrt` (correctly rounded) | **0 ulp on all 1,084 rows** | gated at 0 now |
| `log10` | — | **the 6.7.18 stdlib has no log10** (only `f64_ln`, `f64_log2`) | 1,233 rows pending until M1's correctly rounded log10 (from kimiya's `_km_cr_log10`) |
| `exp`, `atan` | `f64_exp`, `f64_atan` | x87 on x86_64, the polyfills on aarch64: bits differ by target (§2) | not used; M1's correctly rounded primitives replace them |
| `asin`, `acos`, `atan2` | ganita's, built on `f64_atan` | inherit `f64_atan`'s per-target bits | as above |
| `powf(−1, DBL_MAX)` = 1.0 (an even integer exponent) | `ganita_f64_pow` | **+∞** (`0x7ff0…` against `0x3ff0…`): ganita's known \|x\| == 1 defect (tanmatra ADR 0007, impetus ADR 0009) | M1's correctly rounded pow replaces it; no ADR needed |

## 2. Corpus baselines (M0 step 8, 2026-11-11)

`tests/golden/math.tsv` (34,734 rows, 25 functions) and `constants.tsv` (12 rows), replayed by `tests/parity.tcyr` on
**x86_64 ELF, x86_64 PE (wine 11) and aarch64 (qemu-aarch64)**, one target at a time, each run's
`=== parity summary` line checked. Every target: **2 modules, 0 failing; 15,507 ok, 0 failed, 18,006 measured,
1,233 pending.**

The digest is the harness's 64-bit FNV-1a fold of every output word the port returned for the function, in row
order: equal digests mean bit-identical results on every row.

| Function | x86_64 ELF | x86_64 PE | aarch64 |
|---|---|---|---|
| const | `12f052aaad529019` | same | same |
| f64::sqrt | `cd983c72c4b5292e` | same | same |
| f64::abs | `b964a4a06222be77` | same | same |
| f64::ceil | `5a389d80a7c26c2d` | same | same |
| f64::max | `eef666276eded487` | same | same |
| f64::min | `431a92b5449ab1bc` | same | same |
| f64::clamp | `8d3643ed5a793d9d` | same | same |
| f64::is_finite | `87f6e55c27f7530c` | same | same |
| f64::is_nan | `4a02e9ec41493c60` | same | same |
| f64::copysign | `11283fab954a6cc4` | same | same |
| f64::to_radians | `d1f7aa9b4495223d` | same | same |
| f64::to_degrees | `f174e398f581d951` | same | same |
| f64::powi | `c3233c8f6c27c4ec` | same | same |
| f64::sin_cos | `ed65a58db9156cea` | same | same |
| f64::sin | `5c9237a620fe4908` | same | same |
| f64::cos | `a208610badde85c3` | same | same |
| f64::exp | `3d4d27ae6189313e` | same | `1b5f4c9d75f5d327` |
| f64::atan2 | `437185336b2f7645` | same | `7d234dec79df8a94` |
| f64::hypot | `d6d60db0f1471e6a` | same | same |
| f64::powf | `b375cbbe99e0707b` | same | same |
| f64::atan | `c025934e4462a8a2` | same | `6a77a6dd170f8556` |
| f64::asin | `20ed726e1bcba944` | same | `db64c2ede37a8959` |
| f64::acos | `341de438f1489f1f` | same | `15e546e4b24ca78b` |
| f64::tan | `666d07d3c49475f3` | same | same |
| f64::cbrt | `617466f8d5175026` | same | same |

**20 of 25 functions are bit-identical on all three targets**, every exact and arithmetic one among them. The five
that differ on aarch64 are exactly the x87 ones (§1). Cross-checks: the `atan2` digest on x86_64 and the `hypot`
digest everywhere equal impetus's (the same ganita functions over the same seeded grids), and the aarch64 `atan2`
digest equals impetus's polyfill-form digest.

The measured class-C rows, ulp distance from glibc (PE identical to ELF):

| Function → Cyrius candidate | Rows | x86_64 | aarch64 |
|---|---|---|---|
| sin_cos → `f64_sin` + `f64_cos` | 2,156 | 0: 4,222 · 1: 90 | same |
| sin → `f64_sin` | 2,156 | 0: 2,115 · 1: 41 | same |
| cos → `f64_cos` | 2,156 | 0: 2,107 · 1: 49 | same |
| tan → `ganita_f64_tan` | 2,156 | 0: 2,091 · 1: 65 | same |
| exp → `f64_exp` | 1,555 | 0: 1,513 · 1: 42 | 0: 1,438 · 1: 117 (the polyfill) |
| atan2 → `ganita_f64_atan2` | 1,200 | 0: 969 · 1: 231 | 0: 962 · 1: 238 |
| hypot → `ganita_f64_hypot` | 1,715 | 0: 1,301 · 1: 411 · 2: 3 | same |
| powf → `ganita_f64_pow` | 1,560 | 0: 1,430 · 1: 129 · ∞: 1 (§1) | same |
| atan → `f64_atan` | 1,184 | 0: 1,184 | 0: 1,156 · 1: 28 (the polyfill) |
| asin → `ganita_f64_asin` | 1,084 | 0: 891 · 1: 192 · 2: 1 | 0: 886 · 1: 197 · 2: 1 |
| acos → `ganita_f64_acos` | 1,084 | 0: 821 · 1: 263 | 0: 820 · 1: 263 · 2: 1 |

**The corrected oracle's size.** glibc misrounds **125 of the 20,323 class-C rows** against MPFR 4.2.2: log10 92 (7.5 %),
hypot 12, tan 8, sin_cos 5, sin 3, cos 2, asin 1, atan2 1, exp 1; acos, atan, powf and cbrt none. These rows go into
`tests/golden/math.cr.tsv` (`validate.py --write`) when M1's correctly rounded primitives land, and the port must
match the corrected values (kimiya ADR 0006).

**Decisions.**
- **Class C gates at 0 ulp against the corrected oracle** (ADR 0001), function by function, as src/crmath.cyr
  dispatches each in M1. Until then the rows are measured. `cbrt` gates at 0 now.
- No native x87 primitive (`f64_exp`, `f64_atan`) and no ganita function built on one is used in `src/`.

## 3. Oracle provenance findings (M0 steps 5–8)

- **glibc for the transcendentals, except `cbrt`.** sin, cos, sincos, tan, asin, acos, atan, atan2, exp, pow, hypot
  and log10 reach glibc. Rust std links compiler-builtins' own `cbrt`, the `libm` crate's port of CORE-MATH, which
  is correctly rounded. `validate.py` found this: glibc disagreed on 504 of 1,084 rows, MPFR on none.
- **`sin` + `cos` of one argument is one `sincos` call** in the optimized binaries, and glibc 2.44's `sincos` equals
  the separate calls bit for bit (0 of 399,944 arguments differ).
- **`10.0_f64.powf(y)` stays `pow`**; there is no `exp10` rewrite. bijli has no `2.0.powf` site.
- **`powi` with a run-time exponent is compiler-rt's `__powidf2`**, square-and-multiply in a fixed order; constant
  exponents expand inline in the same order.
- **Rust's float `Sum` is a left fold from −0.0** (probe p11): the empty sum is −0.0, a sum of −0.0 is −0.0,
  `[1e16, 1, −1e16]` sums to 0. bijli's 13 `.sum()` sites port as loops from `-0.0`.
- **Constant folding**: LLVM folds asin, acos, tan and cbrt of constants, so the generator black-boxes every input.

## 4. Probes (M0 step 9, 2026-11-11, x86_64 Linux ELF)

`scripts/probes/` (README there). Timings ran beside kimiya's gates (load average about 1), so they are approximate.

| Probe | Question | Result | Decision |
|---|---|---|---|
| p1–p6 | the language review (plan §2) | complex and vector structs with operators, traits + bounded generics + `dyn`, closures, freelist, derive Serialize, `Result` with f64, `(f64, f64)` returns | plan §2, §9.3 |
| p7 | a full-size 53³ grid: six fields, three coefficient and three material arrays, twelve 10-cell CPML Ψ slabs (19.2 MB), 200 alloc/free cycles | 0 failed frees, 0 B bump growth, peak RSS **19,376 kB for a 19,224 kB grid** (flat). `fl_calloc` per array: **8.7 ms/grid**; one `mmap` per grid: **0.85 ms** | a grid is **one mmap**, sliced into its arrays (zeroed by the kernel, no byte-loop memset), unmapped by `Bj*_free` |
| p8 | struct values vs `f64` locals in the Mie recurrence and an NF2FF-style accumulation | bit-identical; structs **2.2×** (Mie) and **2.7×** (accumulation) slower | hot loops on `f64` locals, the API and cold paths on `BjComplex` |
| p9 | a field callback per step (trace_field_line) | direct 50.1 ns/step, `callptr` on a fn pointer 49.8 ns, a capturing closure 61.7 ns (+23 %); identical bits | closures for the public callbacks; the cost is small against a step |
| p10 | `row_diff` four ways against Rust's | all four **bit-identical to Rust** (digest `bb2c6a49c09f2a7b`). Rust 0.84 ns/cell; Cyrius `*f64` index 5.7, pointer step 5.5, unrolled 5.1, **`f64v_*` 1.26 ns** | FDTD kernels on `f64v_*` (with a scalar tail: the open lane-multiple bug), a scalar reference kept for tests |
| p11 | Rust's `.sum()` | a left fold from −0.0 (§3) | sums start from `-0.0` |
| p12 | Touchstone text: `f64_parse_n` and bayan's digits on 120,018 strings (repr, %E, %g, integers, ±0, subnormals, overflow, 30-digit strings) | `f64_parse_n` **= Rust's parse on all 120,018**. bayan's output always reads back, but is **longer than the shortest digits on 104 values** (`3.2557220000000003e+77` for `3.255722e77`) | parse with `f64_parse_n`; write Rust's `{:e}` / `{}` with impetus's `src/fmtf.cyr` (shortest digits, proven against Rust's Display), not bayan |
| p13 | one 3D Yee cell-step, scalar, 30³ | **1,282 µs/step** (47.5 ns per cell-step, 4.3× Rust's 298 µs) | Rust's largest 3D CPML test (about 45M cell-steps) takes about 2 s: full-size class-T scenarios fit `cyrius test` |
| p14 | the last-error record (ADR 0002 row 4) on two threads | each thread reads back only its own payload (100,000 sets each); the main thread's record is separate | `thread_local_alloc()` slot → a per-thread 32-byte record |

Language facts the probes found, beyond plan §2:
- Storing an untyped word (`fl_alloc`'s result, a closure parameter) in a `*T` warns; returning it from a fn
  declared `: *T` does not. The port uses typed helpers (`_bj_alloc_f64(n): *f64`, `_bj_as_vec(p): *BjFieldVector`).
- `var x = f()` with `fn f(): *f64` leaves `x` untyped: declare `var x: *f64 = f();`.

## 5. M1 part 1a: correctly rounded exp, sin, cos, hypot, pow, atan2 and atan (2026-11-11)

- **`src/crmath.cyr` and `src/crmath_data.cyr`** are impetus's (its ADR 0011 and amendments), with the prefixes
  renamed. `scripts/oracle/crtables.py` regenerates the tables byte for byte, and `validate.py` checks them: π and
  ln 2 to 2,048 bits and every table's hi equal MPFR 4.2.2's. impetus's 78 Cyrius-only tests
  (`tests/crmath.tcyr`: the IEEE specials, fast and slow paths agreeing, phase 1 deciding) pass.
- **`src/mathx.cyr`** wraps them as `_bj_exp`, `_bj_sin`, `_bj_cos`, `_bj_sin_cos`, `_bj_hypot`, `_bj_powf`,
  `_bj_atan2`, and **`_bj_atan` = `atan2(x, 1)`**, which is exact: both are the correctly rounded value of the same
  real number, and the specials agree.
- **The corrected oracle**, `tests/golden/math.cr.tsv`: the 125 rows of §2 where glibc misrounds. The harness
  compares them against MPFR's value: **all 24 rows of the dispatched functions equal the corrected value, none
  glibc's**. The rest belong to tan, asin and acos (measured) and log10 (pending).
- **Beyond the corpus** (`scripts/oracle/cr_verify.py`, the driver `scripts/oracle/cr_driver.cyr` built against
  `dist/bijli.cyr`): 100,000 arguments per function (the whole domain, every binade, the implementation's
  thresholds, the special grid), **700,000 calls on x86_64 and on aarch64 (qemu): every result equals MPFR's**.
- **Bounds: 0** for exp, sin, cos, sin_cos, hypot, powf, atan2, atan and cbrt.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): 2 modules, 0 failing; **29,165 ok, 0 failed, 4,324
  measured (tan, asin, acos), 1,233 pending (log10), 24 accepted (ADR 01)**. Digests: 23 of 25 equal on all three
  targets; asin and acos (ganita's, on the native atan) still differ on aarch64 until part 1b/1c replaces them.

| Function | Digest, every target |
|---|---|
| f64::exp | `9763e0c4a6dc5471` (= impetus's: the same grid) |
| f64::sin | `cf422e1f146f3036` |
| f64::cos | `ccbaeb8b2a98fec2` |
| f64::sin_cos | `cca2d7bb6f9b190b` |
| f64::atan2 | `6b403ab79185f15f` |
| f64::hypot | `b904b4aff7f5703d` |
| f64::powf | `e25863c3e07f8bf6` |
| f64::atan | `c025934e4462a8a2` |

## 6. M1 part 1b: correctly rounded acos, log10 and asin (2026-11-11)

- **`src/crmath_k.cyr` and `src/crmath_k_data.cyr`**: kimiya's acos and log10 (its ADR 0006), copied function for
  function with the prefixes renamed (`_km_` → `_bj_k_`): the 24 functions their call closure reaches, the atan and
  ln cores (a double-double fast path with a proven bound and a rounding test, then a triple-double accurate path)
  and the triple-double arithmetic under them. `scripts/oracle/crtables_k.py` (kimiya's `mathx_tables.py`, its ln
  and atan tables) regenerates the data, and `validate.py` checks it is that output.
- **asin is bijli's**, on the same atan core as kimiya's acos: asin |x| = atan2(|x|, √(1 − x²)). 1 − x² is a
  double-double ((1 − |x|)·(1 + |x|) by Dekker's product for |x| ≥ ½, 1 − x² below), its square root a double-double
  by one Newton step, and the core takes the pair in the octant that keeps its quotient ≤ 1, so its rounding test
  and accurate path decide the last bit. |x| < 2⁻²⁶ returns x (asin x − x < x³/6 is below half an ulp), ±1 returns ±π/2, |x| > 1
  and NaN give NaN, and a negative x is the negated result.
- **Beyond the corpus** (`cr_verify.py`): asin on **1,000,000 arguments** (the plan's §6 step 8 floor), acos and
  log10 on 300,000 each, drawn over the whole domain, the binades near 0 and 1, the last 2,000 doubles below ±1,
  the thresholds, every power of ten from 10⁻³⁰ to 10³⁰ and the subnormals: **every result equals MPFR's, on
  x86_64 and on aarch64 (qemu)**.
- **Bounds: 0** for asin, acos and log10, against the corrected oracle: log10 1,141 ok and 92 accepted (ADR 01), asin
  1,083 ok and 1 accepted, acos 1,084 ok. Every one of those 93 rows equals MPFR's value, none glibc's.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): 2 modules, 0 failing; **32,473 ok, 0 failed, 2,156
  measured (tan), 0 pending, 117 accepted (ADR 01)**. **Digests: all 26 equal on all three targets.** tan is still
  ganita's, measured (2,097 rows at 0 ulp, 59 at 1 ulp), until part 1c.

| Function | Digest, every target |
|---|---|
| f64::log10 | `401a912eb68a4da1` |
| f64::asin | `eb26b0ebcaa77e11` |
| f64::acos | `35c6927605dbf37c` |

## 7. M1 part 1c: correctly rounded tan (2026-11-11) — the math primitives complete

- **tan is bijli's**, in `src/crmath.cyr`'s tan section, on impetus's sin and cos phases and reduction. Dividing
  two rounded values is not correctly rounded, so sin and cos stay unrounded, are divided, and the quotient is
  rounded once:
  - |x| < 2⁻²⁷ returns x (tan x − x ≈ x³/3 is below half an ulp); NaN and ±∞ give NaN;
  - phase 1 (|x| ≤ 2²⁰): the i/64 table gives sin x and cos x within 2⁻⁶² each, and their double-double quotient
    is within 2⁻⁶⁰ (2⁻⁶² + 2⁻⁶² + the division's 2⁻¹⁰⁰, rounded up);
  - phase 2: the double-double series, the quotient within es + ec + 2⁻⁹⁸ (es and ec are 2⁻⁷⁴ and 2⁻⁸⁸, swapped by
    quadrant);
  - the slow path, in 364-bit fixed point: Payne–Hanek reduction to [0, π/2), both series, tan y = S/C or
    −C/S by quadrant, the divisor scaled into [2, 4) for a Newton reciprocal, one rounding. No double is within
    2⁻⁶¹ of a multiple of π/2, so numerator and divisor are at least that, and the quotient stays within 2⁻²⁹⁰.
- **Measured** on 100,000 arguments in [−10⁵, 10⁵] and 100,000 in [−1.5, 1.5]: phase 1 decides 98.9 %, phase 2
  the rest, and none reaches the slow path. One call costs about 275 ns against `ganita_f64_tan`'s 93 ns, and the
  slow path about 40 µs (a 100,000-call loop on `dist/bijli.cyr`, with kimiya's bench and a qemu run on the
  machine). bijli calls tan twice, both in `rf.rs`'s `single_stub_match`, so no hot path.
- **Beyond the corpus** (`cr_verify.py`): tan on **1,000,000 arguments** (the plan's §6 step 8 floor: the whole
  domain, the doubles around k·π/2 for k to 64 and at random to 700,000, the double closest to a multiple of π/2,
  the thresholds 2⁻²⁷, 2⁻²⁴, 0.78 and 2²⁰), and the slow path called directly on 157,132 more (`tan_slow`): **every
  result equals MPFR's, on x86_64 and on aarch64 (qemu)**.
- **Cyrius tests**: `tests/crmath_bijli.tcyr`, 33 assertions over bijli's own tan and asin. They cover the IEEE
  specials and thresholds, odd symmetry, MPFR's values (tan at fl(π/4), 1, fl(π/2), fl(π), 355, 2²⁰, 2¹⁰²³ and the
  double closest to a multiple of π/2; asin at ±1, 2⁻²⁶, 0.5, 0.7 and −0.99), the fast path equal to the slow one
  on 500 tan arguments (200 next to k·π/2) and to the accurate one on 2,000 asin arguments, and phase 1 deciding
  ordinary arguments. They pass on all three targets.
- **Bound: 0** for tan, against the corrected oracle: 2,148 ok and 8 accepted (ADR 01), all 8 equal to MPFR's
  value, none to glibc's.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): 2 modules, 0 failing; **34,621 ok, 0 failed, 0
  measured, 0 pending, 125 accepted (ADR 01)**. Every row of the corpus is gated, and all 26 digests are equal on
  all three targets. The oracle's 13 class-C functions are all correctly rounded: sin, cos, sin_cos, exp, atan2,
  hypot, powf, log10, atan, asin, acos, tan and cbrt.

| Function | Digest, every target |
|---|---|
| f64::tan | `b33eeff0806f9a96` |

## 8. M1 part 2a: fmtf, log and error (2026-11-11)

- **`src/fmtf.cyr`**: impetus's f64 printer (its ADR 0010), the bignum and Burger & Dybvig's digit generator line for
  line. Its output goes to a bounded `BjText` buffer instead of a growing string builder, so the error record and
  log lines format into fixed storage. `tests/golden/fmtf.tsv` (impetus's grid: every binade, the powers of ten, the
  2^53 edge, 300 ties and a seeded sweep, both signs): **37,576 rows, Rust's Display and serde_json's form, every
  one byte for byte**. `validate.py` transcribes both printers independently, Python's shortest `repr` with Rust's
  round-up on an exact tie (751 Display rows are such ties), and matches every row.
- **`src/error.cyr`**: the `BJ_ERR_*` codes (Rust's 8 variants in order, plus `BJ_ERR_TRANSPORT`, ADR 0002 row 7) and
  the per-thread **last-error record** (ADR 0002 row 4, amended): the code, the variant's fields, and Rust's Display
  text of the error byte for byte, the reason or context included. `tests/golden/error.tsv`: **1,635 rows**, every
  variant over its payload grid (usize fields to i64::MAX, f64 fields over the specials, binades, powers of ten and
  seeded bits, strings with UTF-8, `%`, quotes and spaces, and reasons built `format!("{prefix}{x}")` as bijli's
  are). Each row checks the code, each field through `bj_last_error_*`, and the text: all exact.
  - **Where the record lives.** Cyrius has no thread-local variables. `lib/thread_local.cyr`'s slots need a
    per-thread block, which a Linux main thread has only after `thread_local_init`. So the Linux main thread uses a
    static record; other threads (and, off Linux, every thread) keep one from `fl_alloc` behind a slot claimed
    once. bijli never calls `thread_local_init`, so it never replaces another library's block. Until
    `thread_create` has run (`_threads_active`), a Linux caller can only be the main thread, so making an error
    costs no syscall. **The PE run caught a defect in the first version**: off Linux, that shortcut moved the main
    thread from the static record to its slot's record once workers started, losing its last error. The shortcut is
    now Linux-only, and `tests/error.tcyr`'s `last_error_per_thread` (two workers, 20,000 errors each, then the main
    thread's own) passes on all three targets.
- **`src/log.cyr`**: Rust's `logging::init` / `init_with_level` as `bj_log_init` / `bj_log_init_with_level`, and the
  `_bj_error` / `_bj_warn` / `_bj_debug` / `_bj_trace` site helpers on sakshi. bijli's gate starts off (Rust with no
  subscriber prints nothing), and init raises sakshi's process-wide level but never lowers it (ADR 0002 row 8).
- **The visibility finding.** A `private` marker applies to the whole file it sits in, and **a bundle is one file**,
  so `dist/bijli.cyr` is private as a whole: a consumer sees exactly the `public` items (17 `bj_*` fns, the
  constants, now `public const`, and `BjError`'s members), as a probe confirmed. **`[test] modules` was the wrong
  way to give tests the library** (plan §2, reversed): a prepended module becomes part of the test's own file
  (`<source>` in a coverage dump), so the test file turned private (`parity_core.cyr` could not call its hooks), and
  `cyrius coverage --exec` could never credit a src/ fn (0 of 17 executed though every one ran). Tests now include
  `tests/lib.cyr`, the `[lib]` list as files of their own, as impetus's and kimiya's do; `surface_parity.py` checks
  it against `[lib]`. Files of internals only (the math, fmtf) carry no marker, so `cr_driver.cyr` includes the math
  sources directly.
- **Coverage**: `cyrius coverage --exec --min 80` joins the gate: **17 of 17 public fns executed (100 %)**. `cyrius
  api-surface` has its snapshot (`docs/api-surface.snapshot`; the tool's only scope lists the vendored stdlib too).
- **Ledger and surface**: error.rs's 10 tests and logging.rs's 1 ported (`tests/error.tcyr`, `tests/log.tcyr`; Rust's
  `test_error_is_debug` checks the variant's name through `bj_error_name`); surface parity 11 items, 0 problems.
- **Costs** (`cyrius bench`, Cyrius-only; another repo's tests ran beside it, load 1.9): a log site with the gate off
  3 ns (the loop's own cost); making an error 122 ns (960 ns in the first version: two syscalls per error to tell the
  main thread, then 311 ns before `_bj_tx_cstr` became one copy loop); an error with an f64 field 2.9 µs; one f64 in
  Display form 1.9 µs (impetus's exact bignum printer, never on a step path).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **4 modules, 0 failing; 73,832 ok, 0 failed, 0
  measured, 0 pending, 125 accepted (ADR 01)**; all 37 digests equal on all three targets. The digest now folds the
  bytes of a text token (its word is an address).

| Function | Digest, every target |
|---|---|
| BijliError::DimensionMismatch | `7ce0df9f10df6a39` |
| BijliError::ZeroCharge | `1cb0386e31fb2aaf` |
| BijliError::Singularity | `4aeb42d525e34032` |
| BijliError::InsufficientResolution | `f59bbe8c7ab5b6c7` |
| BijliError::InvalidPermittivity | `4b2990a05ea1e62c` |
| BijliError::InvalidPermeability | `7a37a32db98ea3b4` |
| BijliError::InvalidParameter/f64 | `dcb4da1e5054485f` |
| BijliError::DivisionByZero | `36f278fcb8b1b69b` |
| BijliError::InvalidParameter | `187d253e26632343` |
| f64::fmt_display | `737b6e256576067a` |
| serde_json::f64 | `b64474e0f3090f73` |

## 9. M1 part 2c: BjComplex and ln (2026-11-11), and digest v2

(Part 2b, the logging, went into 2a.)

- **`src/complex.cyr`**: `BjComplex`, a 16-byte struct value, the port of Rust hisab 1.4.0's `num::Complex`.
  bijli 1.2.0 re-exports that type (`pub use hisab::num::Complex` in polarization.rs), so **its whole surface is
  bijli's API** and is ported, not only the methods bijli calls. The plan's §9.2 first said "ported only if a non-test
  call exists", which would have dropped re-exported surface without an ADR. The file holds:
  - the constructors, `from_polar`, Default and `From`;
  - `abs`, `norm`, `conj`, `norm_sq`, `arg`, `exp`, `ln`, `sqrt`, `powf`, `inv` and `is_zero`;
  - the four operators, `neg`, `scale` (Rust's `z * s` and `s * z`) and `div_f64`;
  - Display;
  - polarization.rs's `complex_real`, `complex_zero` and `complex_from_polar`.

  Complex `sin` and `cos` need correctly rounded `sinh` and `cosh`, new primitives, so they arrive in part 2d. The
  header states **the one time convention** (e^{+jωt}, IEEE handedness) for the whole library.
- **The language review** was p1 (`scripts/probes/p1_complex.cyr`): by-value `T_add` / `T_mul` / `T_div` operators,
  a scalar right operand and struct arrays, bit-exact against Python. One operator fn per name means
  `Complex * f64` is the method `scale`. Cyrius refuses `+=` on a struct value, so hisab's compound assignments are
  written `a = a + b` (the same operations).
- **`ln`**: kimiya's correctly rounded ln (`_bj_k_cr_ln`, with its fast path and accurate entry; the core was already
  in `src/crmath_k.cyr` for log10). `cr_verify.py`: 300,000 arguments, **all equal to MPFR on x86_64 and aarch64**.
- **`tests/golden/complex.tsv`**: **28,779 rows**. The grid crosses the specials with each other, spans every
  seventh decade from 1e-310 to 1e308 with both signs, takes points on the unit circle, and adds seeded values;
  binary operations take six partners each. The rows are classed by what they compute:
  - pure IEEE sequences are **class A**, bit for bit, which is stricter than the plan's class B, and they all pass;
  - rows that call a transcendental are **class C**, at 0 ulp against `complex.cr.tsv`. glibc misrounded inside
    5 of them (2 `exp`, 2 `sqrt`, 1 `powf`), and in each of those the port equals MPFR's value, not glibc's.

  `validate.py` transcribes every operation in hisab's order twice: with glibc's primitives it reproduces all
  28,779 oracle rows exactly, and with MPFR's it writes the corrected oracle.
- **`tests/complex.tcyr`**, 23 assertions on all three targets:
  - value semantics: copies are independent, and a method leaves its receiver alone;
  - the operators equal their fns inside a chained expression;
  - hisab's formulas at the specials: 1/0 is NaN + jNaN by the textbook division, `abs(∞ + jNaN)` is NaN with no
    hypot rescue, `abs(1e160)` overflows (a range-edge candidate), and `powf` under 1e-300 gives 0;
  - Display's edges: `1+-0i`, `1NaNi`, `inf+infi`.
- **Coverage** 29 of 29 public fns executed. **Costs** (Cyrius-only benches, with kimiya's and agnosai's benches
  running beside them, load 2.2): a complex multiply and divide costs 37 ns, and `exp` 403 ns (three correctly
  rounded calls).
- **Digest v2.** The complex run gave equal rows but four unequal digests (`sub`, `div_f64`, `sqrt`, `powf`), and
  they differed **in bit 63 alone**. A generated NaN (∞ − ∞, 0/0, √−1) is `0xFFF8…` on x86 and `0x7FF8…` on ARM:
  IEEE 754 leaves a NaN's sign and payload open, and a probe on both targets confirmed it. The v1 fold
  `h = (h ^ w)·p` carries a difference only upward, so a sign-bit difference reached nothing but the digest's own top
  bit, and **an even number of them cancelled**. That is why `f64::sqrt` (√−1) and `f64::sin` (sin ±∞) had shown
  equal digests in §2 and §5–§8. The guarantee itself was never in doubt: every non-NaN output is checked exactly
  against the x86 oracle on each target. But the digests were weaker than those sections said. v2 folds any NaN as
  one quiet NaN (the class rules already match any NaN with any NaN; ADR 0002 row 6, clarified) and adds a
  right-shift mix to each step, so every bit reaches the whole digest. **With v2, all 64 digests are equal on
  x86_64 ELF, x86_64 PE and aarch64.** The digests in §2 and §5–§8 are v1's; v2's, for every function so far, are
  below.
- **Result, every target**: **5 modules, 0 failing; 102,606 ok, 0 failed, 0 measured, 0 pending, 130 accepted
  (ADR 01)**.

| Function | Digest v2, every target |
|---|---|
| Complex::from_real | `93890af3323383b` |
| polarization::complex_real | `93890af3323383b` |
| Complex::from_polar | `1a5bcded97aba4fc` |
| polarization::complex_from_polar | `1a5bcded97aba4fc` |
| Complex::default | `15e30f71b82b3f49` |
| polarization::complex_zero | `15e30f71b82b3f49` |
| Complex::new | `138ee5f603ef66ac` |
| Complex::abs | `6f941a7f2f954e45` |
| ComplexExt::norm | `6f941a7f2f954e45` |
| Complex::norm_sq | `69cf04fa49c4eda4` |
| Complex::conj | `263b60fdb4ff4662` |
| Complex::neg | `592be8942100dbdd` |
| Complex::inv | `1b1f1d8dc4717f15` |
| Complex::arg | `650917640d35cc4d` |
| Complex::exp | `7d71c1b844cea38f` |
| Complex::ln | `16a502fce112f802` |
| Complex::sqrt | `6378c83897ce0cc0` |
| Complex::fmt | `4f5314e9e82d1a16` |
| Complex::add | `733a991f7bf422ad` |
| Complex::sub | `1557d6338798be51` |
| Complex::mul | `7382450d17a73c68` |
| Complex::div | `8632bf6e16883e9` |
| Complex::mul_f64 | `29e8507f5c559dc8` |
| f64::mul_complex | `29e8507f5c559dc8` |
| Complex::div_f64 | `34631a46c487397d` |
| Complex::is_zero | `715afc37cdf20bd4` |
| Complex::powf | `14a53d59bf180a53` |
| const | `4a4ee3d6e449876a` |
| BijliError::DimensionMismatch | `668b75c621ba8449` |
| BijliError::ZeroCharge | `35eda1ee9785b119` |
| BijliError::Singularity | `4f2dcf083f999d83` |
| BijliError::InsufficientResolution | `1b01215ab3ab9bfb` |
| BijliError::InvalidPermittivity | `4408414f52c5e4fc` |
| BijliError::InvalidPermeability | `6ec818c146d4da14` |
| BijliError::InvalidParameter/f64 | `7936ed192f15f08b` |
| BijliError::DivisionByZero | `14045824a4fbbd7` |
| BijliError::InvalidParameter | `65fa3d84c5daf63f` |
| f64::fmt_display | `562cda75c1763c4c` |
| serde_json::f64 | `d4f03a226e67586` |
| f64::sqrt | `3023e740efc97dc0` |
| f64::abs | `2c74b9f1ad2c373c` |
| f64::ceil | `5071c4cbb2ba4d80` |
| f64::max | `912c1a3db38cfcf` |
| f64::min | `4c864c68de7d2954` |
| f64::clamp | `28c00100ce689ae1` |
| f64::is_finite | `76f0ad36418ef63a` |
| f64::is_nan | `29b13373c47492e6` |
| f64::copysign | `17f3cefae04fdbe9` |
| f64::to_radians | `637906467d68ba9e` |
| f64::to_degrees | `2e02ab9265252d52` |
| f64::powi | `7e6904ca7ebd9dfc` |
| f64::sin_cos | `2e0ed56d67abe021` |
| f64::sin | `cec83064268b1ed` |
| f64::cos | `511ffb33517201e1` |
| f64::exp | `2ff1a66e408d4d69` |
| f64::atan2 | `588c44b505f60fb1` |
| f64::hypot | `74b7edba7f47975` |
| f64::powf | `34a482f213b211ed` |
| f64::log10 | `651a02ab4d58fb7d` |
| f64::atan | `34c1cb536a8ce5d0` |
| f64::asin | `35d87844c33daff4` |
| f64::acos | `478d08451f218cc5` |
| f64::tan | `6b5c3a573a69716c` |
| f64::cbrt | `7fbcb357e22caa6f` |

## 10. M1 part 2d: correctly rounded sinh and cosh, complex sin and cos (2026-11-11)

- **sinh and cosh are bijli's**, in `src/crmath.cyr`'s second bijli section, on impetus's exp phases and slow path:
  - **cosh** |x| = (e^|x| + e^−|x|)/2 and **sinh** |x| = (e^|x| − e^−|x|)/2. Both terms come from exp's phase 1
    (2^-62) or phase 2 (2^-88) as 2^k·(h + l), are summed as a double-double at e^|x|'s scale, and are scaled by
    2^(k−1) once. So nothing overflows before the result does: the phases hold to |x| = 711, and only exp's own
    final scaling stops at 709.7. Each phase's error bound is absolute, twice both terms' errors, so sinh's
    cancellation lowers only how often the rounding test passes, never its correctness.
  - **sinh below |x| = ½** is a 16-term Taylor series in double-double (bound 2^-96), where the difference would
    cancel more than a bit.
  - **The slow path**, in 364-bit fixed point: both exponentials (exp's reduction and series, the rounding left
    out), summed or subtracted at e^|x|'s scale and rounded once; below ½, sinh's series.
  - **Shortcuts:** |x| < 2^-26 gives 1 for cosh and x for sinh. At 2^-26 itself cosh is 1 + 2^-52, not 1: x²/2 is
    exactly half an ulp there, and x⁴/24 tips it up, so the shortcut is strictly below. |x| ≥ 711 gives ±∞, and
    cosh overflows from 710.475860073944 on.
- **Measured**: on 100,000 arguments in [−50, 50], phase 1 decides 99.4 % of the 200,000 sinh and cosh calls,
  phase 2 the rest, and none reaches the slow path. The series decides every one of 100,000 sinh arguments below ½.
  `cyrius bench` (machine quiet, load 1.1): a sinh plus a cosh costs 580 ns over [0, 70), and a complex sin 833 ns.
- **Beyond the corpus** (`cr_verify.py`): sinh and cosh on **1,000,000 arguments each**, and their slow paths called
  directly on 114,877 and 114,843 more. The sets cover the whole domain, the thresholds 2^-26, ½, 708, 709.7,
  709.78 and 711, both sides of the overflow edge, and dense [700, 712]. **Every result equals MPFR's, on x86_64
  and on aarch64 (qemu).**
- **The corpus**:
  - `tests/golden/math.tsv` gains `f64::sinh` and `f64::cosh`, 894 rows each. bijli calls neither itself, but the
    `Complex` it re-exports does: hisab's `sin` and `cos` take glibc's `sinh` and `cosh`. glibc rounds every one of
    those rows correctly, so `math.cr.tsv` is unchanged.
  - `complex.tsv` gains `Complex::sin` and `Complex::cos`, 994 rows each. glibc misrounded inside 3 of them; the
    port equals MPFR's value in each (`complex.cr.tsv`: 8 rows now).
  - `validate.py` transcribes all four with glibc (exact) and MPFR (the corrections).
- **Cyrius tests**:
  - `tests/crmath_bijli.tcyr`, 24 new assertions: the specials, the 2^-26 edge on both sides, MPFR's values at ½,
    1, 3, 22, 710 and the largest finite cosh, fast equal to slow on 940 arguments (40 of them near overflow), and
    phase 1 deciding;
  - `tests/complex.tcyr`, 3 new: sin and cos on the axes, the sign of zero included.

  All pass on all three targets.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **5 modules, 0 failing; 106,379 ok, 0 failed, 0
  measured, 0 pending, 133 accepted (ADR 01)**. **All 68 digests (v2) are equal** on all three targets. The 64 of §9
  are unchanged, and the four new ones follow.

| Function | Digest v2, every target |
|---|---|
| f64::sinh | `2fdbec357845cfe2` |
| f64::cosh | `5e635b625102bd04` |
| Complex::sin | `3a216378d8e33b5c` |
| Complex::cos | `5b3bc04c1482e544` |

## 11. M1 part 3a: field — BjFieldVector, point charges, superposition, energy densities, dipoles (2026-11-11)

- **`src/field.cyr`** ports field.rs's first part, the rest coming in 3b and 3c:
  - `BjFieldVector` by value: `new`, `zero`, `default`, `magnitude`, `magnitude_sq`, `normalized` (with its
    rescaled path for huge or tiny components), `dot`, `cross`, `scale`, `eq`, `to_array`, `from_array` and
    Display;
  - the operators `+ − ×f64` and `neg`, and `scale_left` for `f64 * v`;
  - point-charge E, V and moving-charge B, the superposition, the two energy densities, and the electric dipole's
    field and potential and the magnetic dipole's field.

  Shapes (plan §9.3): a Rust `[f64; 3]` position is a `*f64`. A typed-pointer subscript is a typed f64, and a
  probe confirmed it bit for bit. A `Result<FieldVector>` writes through `out: *BjFieldVector` and returns
  `Ok(0)` or `Err(code)`, leaving `out` untouched on an error; a `Result<f64>` keeps its f64 payload. Rust's `+=`,
  `-=` and `*=` are written `a = a + b` and so on, since Cyrius refuses compound assignment on a struct value.
- **`tests/golden/field.tsv`**: **12,975 rows, every one class A** (+ − × ÷ and sqrt in Rust's order). The grids:
  - vectors: the specials crossed, every 37th binade with both signs, and seeded;
  - source and field-point pairs: a lattice, random points, offsets either side of the singular cutoff r² = 1e-30,
    and the specials;
  - charges: SI-scale, zero, the specials, and superpositions of 0–5 charges, one of them singular.

  Every error row carries Rust's Display text. `validate.py` transcribes each function independently and matches
  all 12,975.
- **Rust tests**: field.rs's 42 tests on these functions are ported to `tests/field.tcyr` under their own names, and
  Cyrius-only checks join them: the cutoff is 1e-30 to the bit, `out` is untouched on an error, and the normalize
  message is Rust's text. **58 assertions, passing on all three targets.** The ledger now counts 53 of 858 Rust
  tests as ported.
- **Surface parity** gains `IN_PROGRESS`, for a module ported Part by Part. field.rs's 3a items are checked now:
  48 matched, 0 problems. Its 21 functions for 3b and 3c count as pending. The six trait impls with a different
  Cyrius form (`From`, the three compound assignments, `f64 * v`, Display) are EXCEPTIONS entries, each naming its
  form.
- **Benchmarks, the first Rust twins** (`docs/benchmarks-rust-v-cyrius.md`): 8 of field's benches. Each Cyrius
  call costs **14–25×** Rust's time:

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | point-charge E | 1.52 ns | 31 ns | 20× |
  | potential | 1.23 ns | 18 ns | 15× |
  | moving-charge B | 1.73 ns | 43 ns | 25× |
  | superposition of 10 | 39.8 ns | 552 ns | 14× |
  | magnitude | 0.50 ns | 9 ns | 18× |
  | cross | 0.99 ns | 18 ns | 18× |
  | each dipole field | 5.5 ns | 84 ns | 15× |

  impetus measured 10–13×. Rust's benches pass constant inputs and black-box only the result, so LLVM may fold
  part of each call, and the Rust column is a floor. The Cyrius run overlapped another project's `cargo bench`
  (load 1.3); it was started without waiting, which it should not have been.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **6 modules, 0 failing; 119,354 ok, 0 failed, 0
  measured, 0 pending, 133 accepted (ADR 01)**; **all 94 digests equal** on all three targets (the field ones below;
  every earlier one unchanged).

| Function | Digest v2, every target |
|---|---|
| FieldVector::new | `72b01b0fd53fd139` |
| FieldVector::magnitude | `2f60dfe27dbbf625` |
| FieldVector::magnitude_sq | `530e8372112c0470` |
| FieldVector::normalized | `696d8d70ff04a1b2` |
| FieldVector::dot | `70dd4f6547171b2d` |
| FieldVector::cross | `74b338e0fd4c3e88` |
| FieldVector::add | `7098164f178bb9a4` |
| FieldVector::sub | `264f93f12420552f` |
| FieldVector::neg | `7dd25cf7f9444002` |
| FieldVector::scale | `744f9269ea73faa8` |
| FieldVector::mul | `744f9269ea73faa8` |
| f64::mul_field_vector | `744f9269ea73faa8` |
| FieldVector::fmt | `22a7ee8ed0fb7fd` |
| FieldVector::from_array | `24050fe88a1107e6` |
| FieldVector::into_array | `24050fe88a1107e6` |
| FieldVector::zero | `5c12873911e8b0d4` |
| FieldVector::default | `5c12873911e8b0d4` |
| field::electric_field_point_charge | `7d35fbcfd36153d9` |
| field::electric_potential_point_charge | `3a8e93094703bbfe` |
| field::magnetic_field_moving_charge | `41f275f887ae2310` |
| field::electric_dipole_field | `7c00e1c229105b7b` |
| field::electric_dipole_potential | `2061d141d5d37a67` |
| field::magnetic_dipole_field | `3fdfb0298d42e286` |
| field::electric_field_superposition | `3768b8946570930b` |
| field::electric_energy_density | `368e46cc416720b3` |
| field::magnetic_energy_density | `5b959635d3878461` |
