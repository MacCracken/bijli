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

## 12. M1 part 3b: field — the Gauss's-law fields, multipoles, images, Green's functions, stress tensor (2026-11-11)

- **`src/field.cyr`** gains field.rs's other 20 non-callback functions:
  - the Gauss's-law fields (plane, sphere outside and inside, line, cylinder outside and inside), the ring and disk
    on their axes (the disk in Rust's cancellation-free form, with its one-sided limits at ±0), and the internal
    `check_radial_distance`, as `_bj_check_radial_distance`;
  - the quadrupole potential;
  - the real and the normalized complex multipole moments, with both associated-Legendre recurrences and Rust's
    pole-accurate spherical coordinates;
  - the plane, sphere and dielectric image charges;
  - the scalar and static Green's functions;
  - the Maxwell stress tensor, its force, and the Poynting flux through a sampled surface.

  Shapes: a `[[f64; 3]; 3]` is 9 f64, row-major. A `(f64, [f64; 3])` image result returns the charge (in a `Result`
  where Rust's is fallible) and writes the position. A complex `(f64, f64)` is a `BjComplex` out. l and m stay
  Rust's u32 and i32: outside those types the port refuses by name, a Cyrius-only error, since Rust's types make the
  case unrepresentable.
- **The first ADR 0002 row 5 mapping.** Rust's `greens_function_scalar` is e^{+ikR}/(4πR), the outgoing wave in
  e^{−iωt}. The port's is **e^{−jkR}/(4πR)**, the outgoing wave in bijli's e^{+jωt}. The harness conjugates the
  port's value back and compares exactly, so a wrong value still fails, and it reports the count after the module:
  *"ADR 0002 row 5: 320 field::greens_function_scalar rows compared as port = conj(oracle)"*. Every other function
  so far is convention-free.
- **`tests/golden/field.tsv`: 22,325 rows (9,350 new).**
  - Class A for the Gauss's-law fields, images, quadrupole, static Green's function, stress tensor and flux, with
    distances either side of the 1e-30 cutoff, negative and NaN.
  - Class C for the multipole moments (atan2, cos, sin; l ≤ 10 for every m, |m| = l + 1 refused, and l up to 200)
    and the scalar Green's function. **All at 0 ulp, and no row needed a correction**: glibc rounded every call
    in them correctly, so there is no `field.cr.tsv`.
  - `validate.py` transcribes all 20 functions independently, with glibc's primitives (exact) and MPFR's (the
    corrections, none), and matches every row.
- **Rust tests**: field.rs's 56 tests on these functions are ported, the `direct_potential` helper included, plus
  Cyrius-only checks of the l/m type messages. `tests/field.tcyr` holds **123 assertions, passing on all three
  targets**. The ledger counts 109 of 858 Rust tests as ported, 102 of field.rs's 105 (its 7
  `trace_field_line` tests come with 3c), and surface parity matches 68 items, with only `trace_field_line` pending.
- **Benchmarks** (machine quiet, load 0.95): charged sphere 14 ns, ring axis 24 ns, disk axis 38 ns. Rust's three
  read 1.23 ns, the `black_box` floor: LLVM folds those constant-input calls, so the ratios (11–31×) compare with
  no work at all. The 3a numbers reproduced on the quiet machine. 11 of 80 benches are twinned.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **6 modules, 0 failing; 128,704 ok, 0 failed, 0
  measured, 0 pending, 133 accepted (ADR 01)**; **all 113 digests equal** on all three targets. The 94 of §11 are
  unchanged, and the new ones follow.

| Function | Digest v2, every target |
|---|---|
| field::electric_field_infinite_plane | `61ba4ecfd28100d4` |
| field::electric_field_charged_sphere | `522c3aa7c06fcd4c` |
| field::electric_field_infinite_line | `7a0b43b1b59fc4fd` |
| field::electric_field_outside_cylinder | `7a0b43b1b59fc4fd` |
| field::electric_field_inside_cylinder | `597828e09ab62715` |
| field::electric_field_inside_sphere | `146c581b3c6f19e7` |
| field::electric_field_ring_axis | `4a0593674452e0df` |
| field::electric_field_disk_axis | `31cc0e6af1e0623c` |
| field::electric_quadrupole_potential | `707d3cf6285b2c13` |
| field::multipole_moment | `714bd0e852839289` |
| field::spherical_multipole_moment | `5ae627701c328e31` |
| field::image_charge_plane | `24f3f6041cf03fc3` |
| field::image_charge_sphere | `6d957079cb835d4b` |
| field::image_charge_dielectric | `3f83c9af53116711` |
| field::greens_function_scalar | `1271e84da12f6cce` |
| field::greens_function_static | `17d9959d2c9868ff` |
| field::maxwell_stress_tensor | `15043abe365d8997` |
| field::stress_tensor_force | `6de79c1edb44fddc` |
| field::poynting_flux_surface | `b2072256d29a879` |

## 13. M1 part 3c: field — trace_field_line (2026-11-11); field.rs complete

- **`bj_trace_field_line(start, step, max_steps, field_fn, out)`**, Rust's forward-Euler tracer:
  - **The callback.** Rust's `Fn([f64; 3]) -> Result<FieldVector>` is a closure or `&fn` called through
    `callptr` as `field_fn(pos, f)`. It returns `BJ_ERR_NONE` with the field in `f`, or an error code, its
    error recorded: a Result can't pass through `callptr` (plan §9.3). Probe p9 measured that call at 62 ns.
  - **The output.** Rust's `Vec<[f64; 3]>` is a **`BjPolyline`**, a caller-owned buffer that grows by `fl_alloc`,
    is reused across calls (`count`, `at`) and is freed by `BjPolyline_free`. On an error it is left empty.
  - **`max_steps`** is Rust's usize, so its 64 bits count unsigned: the loop runs to `i != max_steps`, so `-1` is
    `usize::MAX`, as Rust's huge-max-steps test needs. The first reservation is Rust's min(max_steps + 1, 4096).
  - **A null point** is tested before `normalized`, so a successful trace leaves the last-error record as the
    field function left it, as Rust's discarded `Err` has no effect.
  - **A buffer that can't grow** returns the new `BJ_ERR_ALLOCATION`, where Rust aborts: **ADR 0002 row 9**.
  - The two `tracing::debug!` sites log through `_bj_log_on` and a stack BjText.
- **`_bj_fmt_f64_debug`**, Rust's `{:?}`, for the positions in two of the error messages: exponential below 1e-4
  and from 1e16, otherwise positional with at least one fractional digit. **18,788 new `f64::fmt_debug` rows** in
  `fmtf.tsv`, on impetus's grid, all exact, and `validate.py` transcribes it independently.
- **`tests/golden/field.tsv` +36 `trace_field_line` rows** (2,265 coordinates). Each row names a field-function
  scenario that the oracle (a Rust closure) and the harness (a Cyrius `&fn`) both build:
  - fields: uniform (a 1e-35 field, zero and −0 fields, 1e300 components), a ±point charge (the line escaping, or
    falling into a sink), electric and magnetic dipoles, a linear quadrupole;
  - ends: a sink reporting Singularity, a propagated ZeroCharge, NaN and infinite fields;
  - rejections: steps of −1, 0, −0, NaN and ∞, and non-finite starts;
  - `usize::MAX` steps.

  All exact, error texts with Rust's Debug positions included (`[1e-5, -0.0, 1e16]`). `validate.py` replays every
  scenario with its own transcription of the 3a functions.
- **Rust tests**: field.rs's 7 trace tests are ported with closures, as Rust writes them, and Cyrius-only checks
  join them: the polyline growing past 4,096 points and being reused, an error emptying it, and a null point
  leaving the record alone. **`tests/field.tcyr`: 140 assertions; field.rs's 105 tests are all ported.** Lint
  compiles a file standalone, and a capturing closure there needs `include "lib/alloc.cyr"`, which the compiler
  asks for by name; the test and bench files with closures now carry it.
- **Surface parity**: field.rs is a ported module, 69 items matched, 0 problems.
- **Benchmark**: `field/trace_field_line_10` costs 1.09 µs against Rust's 318 ns (3.4×). Rust allocates a Vec per
  call, while the port reuses its polyline. 12 of 80 Rust benches are twinned.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **6 modules, 0 failing; 147,528 ok, 0 failed, 0
  measured, 0 pending, 133 accepted (ADR 01)**; **all 115 digests equal**, the 113 of §12 unchanged.

| Function | Digest v2, every target |
|---|---|
| field::trace_field_line | `382668389d61f92e` |
| f64::fmt_debug | `708dfd83398d149f` |

## 14. M1 part 3d: maxwell and charge (2026-11-11)

- **`src/maxwell.cyr`**: maxwell.rs's ten functions, free and by value. The material checks run in Rust's order, ε
  before μ, and refuse NaN and every value ≤ 0. Skin depth checks ω, μ and σ one at a time, so a pair of negatives
  whose product is positive is still refused. The error texts are Rust's (`"{name} must be positive for skin depth,
  got {value}"`).
- **`src/charge.cyr`**: `BjPointCharge` and charge.rs's eleven functions.
  - **The struct keeps Rust's field names.** Rust's `position: [f64; 3]` and `velocity: [f64; 3]` are
    `(f64, f64, f64)` tuples stored inline, so `q.position[0]` is `q.position.0`. Probe **p15** checks the
    typed element reads, a literal store, a struct copy, by-value passing and `&x + offset` as a `*f64` (the
    language review: a field can't be a fixed array, and `&x.f` is refused by name). The struct is 64 bytes.
  - **Methods.** The methods copy the tuple into a local `f64[3]` for the field functions, which take a `*f64`.
  - **Constructors.** `BjPointCharge_new(charge, mass, position, velocity, out)` writes `out` only on Ok;
    `electron` and `proton` return the struct.
  - **Errors.** charge.rs's private `check_mass` keeps its two kinds: DivisionByZero for ±0, InvalidParameter for
    NaN or a negative mass. `PointCharge::new` refuses a zero mass as InvalidParameter, as Rust does: it divides by
    nothing, so the two kinds are not an inconsistency to fix.
- **New golden files, class A throughout:**
  - **`maxwell.tsv`, 20,349 rows.** The material and wave values (specials, both zeros, negatives, ε₀, μ₀, c, copper's
    σ, every 131st binade, seeded) are crossed in pairs for the five two-argument functions. Skin depth is crossed
    over a 16-value list (16³ triples, the paired negatives among them) plus 300 seeded conductor cases.
  - **`charge.tsv`, 12,277 rows.** It reuses field.rs's grids (`m_field.rs`'s are now `pub(crate)`): its (source,
    field point) pairs, the r² = 1e-30 cutoff on both sides among them, for the methods and Coulomb's force; its
    vectors for the Lorentz force; its distances for the Coulomb energy. The masses cross field.rs's charges for
    `new`, and the cyclotron and Larmor arguments cross the specials and SI scales (4,800 Larmor rows).
  - **`validate.py`.** It transcribes both modules independently. The point-charge field, potential and
    moving-charge B are now module-level helpers, which the field and charge checks share.
- **Rust tests**: maxwell.rs's 18 (`tests/maxwell.tcyr`, 39 assertions) and charge.rs's 26 (`tests/charge.tcyr`, 62
  assertions) are all ported. Cyrius-only checks join them: the recorded error codes and Display texts, the tuple
  fields read back, and `out` untouched on an error. The ledger stands at 160 of 858.
- **Surface parity**: 5 modules ported, 101 items matched (32 new), 0 problems; the struct fields match by name.
- **Benchmarks** (17 of 80 Rust benches now twinned; the part 3c rows are within 3 % of §13's run):

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `maxwell/wave_speed` | 1.24 ns | 15 ns | 12.1× |
  | `maxwell/impedance` | 1.24 ns | 13 ns | 10.5× |
  | `maxwell/skin_depth` | 1.29 ns | 27 ns | 21.0× |
  | `charge/coulomb_force` | 5.54 ns | 43 ns | 7.8× |
  | `charge/lorentz_force` | 1.51 ns | 61 ns | 40.5× |

  Rust's bench arguments are constants, so LLVM inlines the whole call and folds most of it. The port pays a call
  and the Result. `lorentz_force` composes three by-value `BjFieldVector` operations (cross, add, scale), each
  copying 24-byte structs: probe p8's struct cost (2.2–2.7× over f64 locals), three times over. It is a finding
  with its numbers, not a change (speed is measured, not gated).
- **For the JSON Part**: `#derive(Serialize)` takes neither a tuple field nor a fixed array (guide, *Derive Serialize
  on a struct*), and serde writes `"position":[x,y,z]`. So `BjPointCharge`'s writer will write the two arrays by
  hand. A flat `px, py, pz` would have needed that too.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **8 modules, 0 failing; 180,154 ok, 0 failed, 0
  measured, 0 pending, 133 accepted (ADR 01)**; **all 139 digests equal**, the 115 of §13 unchanged.

| Function | Digest v2, every target |
|---|---|
| PointCharge::new | `5f66f64b93dd8d79` |
| PointCharge::electron | `7f928a29355c790e` |
| PointCharge::proton | `33a0a40bc641ec6c` |
| PointCharge::electric_field_at | `55439647304ccf15` |
| PointCharge::electric_potential_at | `5d0051c79b612c87` |
| PointCharge::magnetic_field_at | `30d35891a8cab153` |
| PointCharge::lorentz_force_in | `3be703a9d6e419cf` |
| charge::coulomb_force | `3ffc7fb063fb41f6` |
| charge::lorentz_force | `150d626b3dbf856a` |
| charge::dipole_moment | `4015cfd9c418ec4` |
| charge::potential_energy | `4015cfd9c418ec4` (the same products as dipole_moment over the same grid) |
| charge::coulomb_potential_energy | `282e990e5a1d44fa` |
| charge::cyclotron_frequency | `70307729c72e9591` |
| charge::larmor_radius | `275b8f2d129dd7d6` |
| maxwell::gauss_electric_divergence | `86dcc5146508918` |
| maxwell::gauss_magnetic_divergence | `509c42b1fd1c5dbb` |
| maxwell::displacement_current | `6cc785941eed9fa3` |
| maxwell::wave_speed | `22eae1326a0af53f` |
| maxwell::refractive_index | `2f12f45c2d468080` |
| maxwell::impedance | `33ab6cf945d39235` |
| maxwell::free_space_impedance | `53335ce0be4f5134` |
| maxwell::wavelength | `2a51a212b9db5c1` |
| maxwell::frequency | `4d7f9c82d5504df7` |
| maxwell::skin_depth | `46c2da50c032dc3d` |

## 15. M1 part 4: the JSON of BjComplex, BjFieldVector and BjPointCharge (2026-11-12)

- **`src/json.cyr`** (D5, amended): serde_json 1.0.154's documents for M1's three serde types, written and read.
  - **Writers, by hand on the derive's printer.** `BjComplex_to_json`, `BjFieldVector_to_json` and
    `BjPointCharge_to_json` append serde's shape (fields in Rust's order, no spaces, `[f64; 3]` as an array) to a
    `str_builder`, every f64 through `bayan_f64_to_json`, the printer `#derive(Serialize)` uses: the same text the
    derive would write, `null` for NaN and ±∞ as serde_json writes. The document is built on the stack (the longest
    is 251 bytes) and appended in one grow, so a refused allocation returns `BJ_ERR_ALLOCATION` with the builder as
    it was. The derive itself isn't used. **Probe p16b** shows why: on 6.7.18 a derive in a second `private` file
    fails to compile (`'_cy_json_skipval' is private to its file`; bijli's two types live in two private files),
    and in a library the derive exports its zero-filling, bump-allocating decoders as public API.
  - **Readers, on bayan's tree with serde's semantics.** `T_from_json(buf, len, out)`: Ok(0), or
    `BJ_ERR_DESERIALIZE` (ADR 0002 row 10) with the reason in the record and `out` untouched. A struct reads from an
    object by name (unknown keys skipped) or from an array of exactly its fields; a `[f64; 3]` from an array of
    three numbers; an f64 from any number, never `null`. The reasons are serde's words where they apply
    (`missing field `z``, `invalid type: null, expected f64`, `invalid length 2, expected struct FieldVector with 3
    elements`) and bayan's for syntax, at the byte offset in the caller's text.
  - **Where bayan reads otherwise than serde_json** (**probe p16**): one read-only pass before bayan puts it
    right. The integer `-0` (bayan's integer 0, so the sign would be lost) and integers outside i64 (bayan refuses
    them) get `.0` appended, which bayan reads as a correctly rounded float. A raw control character in a string
    (bayan takes it, serde_json and RFC 8259 §7 refuse it) and bytes that are not UTF-8 in a string (a Rust `&str`
    can't hold them) refuse the document. bayan's own error wins when it comes earlier in the text, and offsets
    after a rewrite are mapped back. A document with none of these is parsed in place. `1e400` (bayan: +∞, no
    error) is refused as serde_json refuses it.
  - **Memory.** Every node comes from a decode region on the reader's stack: freelist chunks of one 4,096-byte
    class, chained, and given back when the read ends (a request larger than a chunk gets a block of its own). A
    read leaves both heaps as it found it, and two threads reading at once share nothing (impetus's reader keeps one
    module-level arena and isn't thread-safe).
- **What serde_json does that the port doesn't** (ADR 0002 rows 11 and 12):
  - **Its default float parse is best-effort** (`float_roundtrip` is off in bijli's build): it misreads 6,342 of
    the corpus's documents by 1–2 ulp, its own output among them. The elementary charge as serde_json writes it,
    `1.602176634e-19`, reads back as `0x3c07a4da290c1654`; correctly rounded (Rust's `str::parse`, Python, the
    port) it is `…653`. It also refuses 48 numbers that round correctly to `f64::MAX`
    (`1.7976931348623158e308`, `17976931348623157` followed by 292 zeros) as out of range. The port reads every
    number correctly rounded, through ADR 0001's corrected oracle: `tests/golden/json.cr.tsv`, 6,390 rows, written by
    `validate.py` from Python's `float`.
  - **It checks only what it deserializes.** An unknown field's value is skipped with no check of repeated keys,
    lone surrogate escapes or nesting (serde_json's recursion limit doesn't reach a skipped value: 200 levels read),
    while the same defects in the struct's own keys refuse the document. The port reads the whole document as one
    tree and refuses any of them anywhere: 54 rows, which the oracle flags (`i:` 1 a repeated key, 2 a lone
    surrogate, 4 nesting past 128) and the harness reports as row 12's.
- **`json.tsv`, 14,758 rows** (`m_json.rs`):
  - **Writer rows, class A** (3,558): 1,186 value sets per type (the specials, the notation edges 1e15–1e23 and
    1e-5–1e-7, SI and bijli constants, every 131st binade of both signs, 600 seeded bit patterns of every class,
    NaN payloads among them, 200 seeded SI-scale values). The output is serde's text *as parsed tokens*: each
    number an `f:` token (its correctly rounded value) and each run between numbers an `s:` token. The port's text,
    tokenized the same way (`_par_got_json`), must equal it: the same punctuation, keys and `null`s byte for byte,
    and the same bits for every number in whichever notation.
  - **Reader rows, class C at 0 ulp** (11,200): serde's own text for every writer value; the same values in bayan's
    notation, in Rust's `{}` (every digit, no exponent: 300-digit integers, 320-digit fractions) and `{:e}`; 120
    number spellings in each field position (signs, `-0`, integers past i64 and u64, the overflow, underflow and
    half-way edges of f64, 800-digit significands, `1e00000000000000000000001`, and refused forms: `01`, `.1`, `1.`,
    `+1`, `NaN`, `Infinity`, `0x10`, a fullwidth digit); 1,200 seeded decimals per type of up to 40 digits with
    exponents across the range; and the shapes: field order, JSON and non-JSON whitespace (form feed, NBSP,
    U+2028, BOM), unknown keys of every kind (deep, 20 KB strings, 2,000-element arrays, 100-key objects, escapes,
    lone and paired surrogates), escaped field names, missing and repeated fields, wrong types in each field, the
    sequence form (short, long, nested), every prefix of a document, trailing data, comments, nesting at 120–200.
  - **`validate.py`** recomputes every row from Python's json and a transcription of serde's derive semantics,
    numbers correctly rounded: the writer rows' shape and values, the readers' outcomes, and the flags. Where its
    reading and the oracle's differ, only a misrounding or an overflow judgment may explain it, and that row is a
    correction; anything else is a problem. It agreed on every row the first time it ran.
- **Tests**: `tests/json.tcyr`, 60 assertions in 7 groups, all Cyrius-only (rust-old tests no serde of these types):
  the exact text, round trips of 3,000 seeded f64 bit for bit, the reasons and byte offsets, `out` untouched, 3,000
  reads leaving the freelist's arena and the bump heap unchanged, and four threads reading 16,000 documents at once,
  each its own answer and reason. The ledger stays at 160 of 858.
- **Surface parity** now checks the serde derives: every `#[derive(Serialize)]` / `#[derive(Deserialize)]` of a
  ported module needs `Bj<T>_to_json` / `Bj<T>_from_json`, and hisab's re-exported `Complex` its pair. 107 items
  matched (6 new); the pending modules' 30 serde types count among their items.
- **Benchmarks** (Cyrius-only: rust-old benches no serde; `cyrius bench` on a quiet machine, ushma's gate done):

  | Bench | Cyrius |
  |---|---|
  | `json/field_vector_to_json` | 2.87 µs |
  | `json/field_vector_from_json` | 1.94 µs |
  | `json/point_charge_from_json` | 4.59 µs |

  The writer's time is bayan's printer: `bayan_f64_to_json` alone takes 710–800 ns per f64 (1.5, c, e), so three
  numbers are about 2.3 µs of the 2.87, and the one builder append 161 ns. src/fmtf.cyr's `_bj_fmt_f64_json`,
  which writes serde_json's own text (zmij's notation, `1e16`), takes 0.42–1.92 µs per f64 over the same values, so
  a writer whose text equals serde_json's byte for byte would cost about twice as much. A reader's time is bayan's
  tree (a call per byte through its bounds-checked slices) and the decode region's two freelist calls. Rust has no
  twin to set beside these numbers. They are findings with their numbers, not changes.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **9 modules, 0 failing; 188,486 ok, 0 failed, 0
  measured, 0 pending, 6,559 accepted (ADR 01: 6,505; ADR 0002 row 12: 54)**; **all 145 digests equal**, the 139 of
  §14 unchanged. `tests/json.tcyr`'s 60 assertions pass on all three, its four threads included.

| Function | Digest v2, every target |
|---|---|
| json::Complex::to_string | `4d9ed1ed8014f65d` |
| json::Complex::from_str | `1c5a58b59a93cae1` |
| json::FieldVector::to_string | `16d67f22edc24710` |
| json::FieldVector::from_str | `2d9a9d3270002358` |
| json::PointCharge::to_string | `4429d64e05c2e255` |
| json::PointCharge::from_str | `2b9687ddb5ad7846` |

## 16. M1 part 5a: the correctly rounded machinery for the range-edge review (2026-11-12)

- **The review's scope** (plan §0.1, the user's decisions of 2026-11-12): the survey of M1's chains found four
  families. (1) Raw sums of squares and what uses them: `magnitude`, complex `abs`, `ln`, `sqrt`, `inv`, `/` and
  `is_zero`, and every field, Green's, image and Coulomb chain on r². These go to 0, ±∞ or NaN, or lose most of their
  bits, past about 1e±154; their cubes (`inv_r³`) fail near 1e±102, and the quadrupole's r⁵ near 1e61. The
  point-charge field of q = 1 is exact at r = 1e100, 2 ulp off at 1e105, 3×10¹¹ ulp off at 1e110 and 0 from about
  1e120, where the true field is normal. (2) Single products, quotients and roots of inputs and constants:
  `wave_speed(1e-160, 1e-160)` is 3.6×10¹⁰ ulp off and `(1e160, 1e160)` gives 0; in range `wave_speed` isn't
  correctly rounded on 28 % of seeded inputs and `skin_depth` on 22 %. (3) The `r² < 1e-30` thresholds. (4) hisab's
  `powf`, 0 for every |z| < 1e-150 whatever the exponent. The user chose, each the recommended option: Rust's chain
  with an unbounded exponent at the edges, round once in range for the single products, exact thresholds, and no
  `powf` cutoff.
- **`src/crmath_kx.cyr`, 78 items, generated** by `scripts/port/kimiya_closure.py` from kimiya's `src/mathx.cyr`:
  the call closure of `_km_cr_mono` and its nine shapes (x·y/z, x·y·z, x/(y·z), x·y·z/w, five factors, x·y/(z·w),
  x·y·z/(v·w), four factors, x/(y·z·w)), `_km_cr_sqrt_div` (√(x/y)), `_km_cr_mul_sqrt_ratio` (w·√((x/y)·2^e)), the
  multi-word integers and `_km_cr_ratio` (RN(p/q·2^e)), and binary64 with an unbounded exponent (`_km_ux_*`, without
  its atan2: bijli's correctly rounded atan2 is scale-free). The prefixes are renamed as src/crmath_k.cyr's copy was,
  and what that copy already holds is used from it (checked identical). Comment lines the rename pushed past 120
  characters are rewrapped, and the file goes through `cyrius fmt`. `--check` exits 1 when either output differs
  from what kimiya's current text gives.
- **`tests/crmath_kx.tcyr`, 148 assertions in 11 groups**, from kimiya's `tests/mathx.tcyr` by the same script (its
  ux, multi-word, ratio, monomial, B2 and B3 groups). The script drops six exact passages, the ux atan2 checks, and
  refuses to run if kimiya's text no longer has them. Every expected value is kimiya's (MPFR 4.2.2 or Python's
  Fraction). They pass on x86_64 ELF, x86_64 PE and aarch64.
- **`scripts/oracle/cr_verify.py`'s exact families**: the generic monomial and its nine shapes, √(x/y),
  w·√((x/y)·2^e) and the five unbounded final operations (×, ÷, +, −, √), against exact values. The monomials and
  operations use Python's Fraction, rounded by CPython's int division, which rounds once, ties to even, onto the
  subnormal grid. The roots use an integer square root and an exact midpoint comparison. The arguments cover every
  binade, subnormal factors, products that leave the range both ways, midpoints (a 54-bit 3y, w·3), exact squares,
  cancellation, and exponents to ±4,300. `cr_driver.cyr` now takes up to twelve words per call. **1,700,000 calls on
  x86_64 and aarch64: 0 mismatches**, and the MPFR families still pass after the driver change (319,071 calls).
  `mul5`'s arguments keep kimiya's documented domain: the exact path's capacity is four doubles and a 41-bit integer
  per side, and past it the result is NaN by design. No bijli formula has five doubles on a side.
- **The corpus is unchanged**: 9 modules, 188,486 ok, 6,559 accepted on all three targets, all 145 digests as §15.

## 17. M1 part 5b: round once (2026-11-12)

- **ADR 0002 row 13**: every formula that is one product, quotient or root of inputs and constants is correctly
  rounded for finite operands (the user's "round once" decision). Fourteen functions change:
  - maxwell: `wave_speed` (1/√(ε·μ)), `refractive_index` (√(ε_r·μ_r)), `impedance` (√(μ/ε)), `skin_depth`
    (√(2/(ω·μ·σ)));
  - charge: `coulomb_potential_energy` (kq₁q₂/r), `cyclotron_frequency` (|q|B/m), `larmor_radius` (m|v⊥|/(|q|B));
  - field: `electric_energy_density` (ε₀E²/2), `magnetic_energy_density` (B²/(2μ₀)), and the `charged_sphere`
    (kQ/r²), `inside_sphere` (kQr/R³), `infinite_line` and `outside_cylinder` (λ/(2πε₀r)) and `inside_cylinder`
    (ρr/(2ε₀)) fields.

  `free_space_impedance` is rounded once too; Rust's value already was. A physical constant enters as its double
  and π as Rust's. An infinite or NaN operand keeps Rust's chain (the callers test finiteness first), and so does
  every error. A single IEEE operation (`wavelength`, `frequency`, the divergences, `dipole_moment`,
  `potential_energy`, the plane field) was already rounded once and is unchanged. The functions with a sum inside
  (ring, disk, point charges, dipoles, Green's) get the edge treatment in 5d.
- **`src/crmath_mono.cyr`**: bijli's RN(√(Π num / Π den)), up to 3 over 4 positive finite factors, for
  `wave_speed`, `refractive_index` and `skin_depth`. kimiya's shapes stop at √(x/y) and w·√((x/y)·2^e), and no
  exact rewrite takes 1/√(ε·μ) or √(2/(ω·μ·σ)) to one of them.
  - **Fast path**: kimiya's significand split (`_bj_k_mono_sig`) and double-double products (`_bj_k_mono_prod`), a
    double-double quotient and one Newton step. The value is within 13u² of the exact root, inside kimiya's
    rounding test's 32u² budget (`_bj_k_dd_round_ends`), and `_bj_k_mono_fin` scales it.
  - **Exact path**: R = ⌊√(P·2^S/Q)⌋ to 56 or 57 bits by exact comparisons of R²·Q with P·2^S, a sticky bit, and
    kimiya's `_bj_k_cr_from_int`, which rounds once onto the binade's grid or the subnormal one.
  - Its multi-word integers are its own, 12 limbs of 32 bits: R²·Q reaches 2^327, past kimiya's 2^256. Cyrius's `>>`
    is the logical shift (`>>>` the arithmetic one), so a limb product's carry is `t >> 32`.
  - **Checks**: `cr_verify.py`'s `sqrt_mono` and `sqrt_mono_acc` families, 1,000,000 calls on x86_64 and aarch64 (the
    full function and the exact path alone), all equal to an exact integer square root. The arguments cover every
    size in the domain, every binade, subnormal factors, roots past DBL_MAX and on the subnormal grid, exact squares,
    and both maxwell shapes.
- **The impedances use kimiya's `_bj_k_cr_sqrt_div`**, which gives the same bits as the general root for one factor
  over one (72 ns against 123 ns).
- **The mappings**: `tests/golden/<module>.adr.tsv`, written by `validate.py --write` from exact values
  (`scripts/oracle/exact.py`: Fraction, an integer square root, IEEE's signed zero for a zero factor), never from
  the port. Each row is fn, inputs, the port's outputs, the ADR 0002 row and a note. The harness
  (`tests/parity_core.cyr`) loads it beside `<module>.cr.tsv` into one table, tags each row with its ADR 0002 row,
  compares the port against it, and reports each row number's count (`ADR 0002 mappings, …: row 13 N accepted`).
  The summary now reads `accepted (ADR 01: …, ADR 0002: …)`.
  - **2,031 mapped rows**: maxwell 1,613 (wave_speed 642, impedance 333, refractive_index 322, skin_depth 316),
    charge 248 (larmor 137, cyclotron 67, Coulomb energy 44), field 170 (inside_sphere 70, magnetic energy 32,
    electric energy 19, line and outside cylinder 15 each, charged sphere 11, inside cylinder 8).
  - Those counts are exactly the rows the harness failed before the mappings existed: the exact reference changes
    precisely the rows the port changed.
  - The rest of these functions' golden rows (in range, where Rust's chain happened to round once, and every error
    row) stay Rust's bits.
- **Tests**: a `cyrius_round_once (ADR 0002 row 13)` group in `tests/maxwell.tcyr` (9 assertions), `tests/charge.tcyr`
  (4) and `tests/field.tcyr` (5) pins the edges Rust fails outright, among them:
  - 1/√(1e-160·1e-160) = 1e160, where Rust gives 1.0000055664551363e160;
  - √(2/1e600) = 1.4142135623730952e-300, where Rust gives 0;
  - kQr/R³ at 1e300 = 8.987551792261171e-291, where Rust gives NaN;
  - B²/(2μ₀) at B = 1e-160, exact on the subnormal grid, where Rust loses 15 bits;
  - infinite operands keeping Rust's results.
- **Result, every target**: 9 modules, 0 failing; **186,455 ok, 0 failed, 0 measured, 0 pending, 8,590 accepted (ADR
  01: 6,505; ADR 0002: 2,085)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64. Fourteen changed,
  for the fourteen functions above (`free_space_impedance`'s is unchanged); the other 131 are §16's.

| Function | Digest v2, every target |
|---|---|
| charge::coulomb_potential_energy | `5ce7145fa9a402d6` |
| charge::cyclotron_frequency | `568aae27e1dcff05` |
| charge::larmor_radius | `1d41c82248b61a02` |
| field::electric_energy_density | `49ecf4e86ae15ee1` |
| field::magnetic_energy_density | `48ab3b2a35bf0eed` |
| field::electric_field_charged_sphere | `42f9335f2ee931e5` |
| field::electric_field_inside_sphere | `4054ea92538a2564` |
| field::electric_field_infinite_line | `22b383a69be51ff` |
| field::electric_field_outside_cylinder | `22b383a69be51ff` (the line's) |
| field::electric_field_inside_cylinder | `5185a3a13729d380` |
| maxwell::wave_speed | `203d2973e8936288` |
| maxwell::refractive_index | `5b31a22f587bb545` |
| maxwell::impedance | `70b4787eec185117` |
| maxwell::skin_depth | `6f7fbf4ff1f1cf3` |

- **Benchmarks** (`cyrius bench` once agnosai's benches were done; the other rows are within 3 % of §15's):

  | Bench | Rust | Cyrius (5a) | Cyrius (5b) | Ratio to Rust |
  |---|---|---|---|---|
  | `maxwell/wave_speed` | 1.24 ns | 15 ns | 141 ns | 114× |
  | `maxwell/impedance` | 1.24 ns | 13 ns | 84 ns | 68× |
  | `maxwell/skin_depth` | 1.29 ns | 27 ns | 186 ns | 144× |

  Correct rounding costs a double-double root and its rounding test: about 70 ns for kimiya's √(x/y), and 120 ns for
  the general monomial root (its significand split and products through memory; Cyrius inlines none of it). Rust's
  numbers are LLVM folding constant arguments. These are findings with their numbers (speed is measured, not
  gated). A shape-specific fast path for 1/√(x·y) and √(2/(x·y·z)), like kimiya's √(x/y), would close part of the
  gap.

## 18. M1 part 5c-1: BjComplex's range edges and powf's cutoff (2026-11-12)

- **ADR 0002 row 14**: hisab squares and multiplies raw parts, so past [2^-511, 2^511) a square, a cross product or
  their sum leaves the normal range although the result is representable. At such an edge (every part finite, one
  outside the range) `src/complex.cyr` evaluates hisab's chain with an unbounded exponent (kimiya's `_bj_k_ux_*`): every
  operation is rounded to 53 bits and the last one once into binary64. This covers `abs`/`norm`, `norm_sq`, the product,
  `inv`, the quotient, `is_zero`'s comparison, and `sqrt`'s two roots and its polar products. `ln` and `powf` take the
  unbounded |z| where it is a normal double; past DBL_MAX or below 2^-1022 their transcendental of an unbounded |z|
  waits for 5c-2, and those rows keep hisab's value.
  - **The guard** (`_bj_cx_plain`, `_bj_cx_plain4` for the binary operations) is all or nothing. Any non-finite
    part keeps hisab's chain. With every part zero or in [2^-511, 2^511), every square and product is normal and
    every sum below 2^1023, so the unbounded evaluation IS hisab's in binary64 and the rows keep Rust's bits.
  - A first version tested each operand on its own, so (1 − j∞)/(1e182 − j6.5e183) took the unbounded path with an
    infinite part. The harness's 10 failures showed it.
- **ADR 0002 row 15**: `powf` has no cutoff. hisab returned 0 for every |z|² < 1e-300, whatever the exponent. Now
  every z ≠ 0 takes from_polar(|z|ⁿ, n·arg z), and z = 0 follows pow's rules: 0, 1, +∞ (im +0) and NaN for n > 0,
  n = 0, n < 0 and a NaN n.
- **The model** (`validate.py`'s `cx_model`; `exact.py` gains `rnd53` and `rnd53_sqrt`, 53-bit rounding with an
  unbounded exponent) evaluates each in-scope function exactly, with IEEE's signed zeros through the products and
  sums and MPFR's correctly rounded atan2, cos, sin, ln and pow. Where the port runs the unbounded chain, a row whose
  model value differs from the oracle's goes to `tests/golden/complex.adr.tsv`, and the corrected oracle
  (`complex.cr.tsv`) skips it. Where the port keeps hisab's chain, the model must agree with that chain, on every row:
  the guard's proof.
  - **1,592 mapped rows**: div 514, sqrt 181, abs and norm 180 each, inv 178, ln 169, mul 117, powf 31 (row 14) and 39
    (row 15), norm_sq 2, is_zero 1, every one passed after the guard fix; before it, the port's ten mixed rows (a
    non-finite part in one operand, an edge in the other) failed against hisab's chain.
- **Tests**: `tests/complex.tcyr`'s `cyrius_range_edges (ADR 0002 rows 14 and 15)`, 13 assertions on the model's
  values:
  - |1e160 + j1e160| = 1.414213562373095e160;
  - (1e200 + j1e200)² = (0, ∞), where hisab gives (NaN, ∞);
  - 1/1e160 rounded once, and 1/(1e200 + j1e200) = 5e-201 − j5e-201;
  - is_zero at 1e200 against 1e201;
  - √(1e300 + j1e300) and ln|1e200 + j1e200| = 460.8635921890891;
  - (1e-160)⁻¹ = 1e160, 0⁰ = 1, 0⁻¹ = +∞, (1e-151)² = 9.999999999999998e-303;
  - a non-finite part keeping hisab's chain.

  Two older assertions that pinned hisab's edge values (abs(1e160) = ∞, powf of |z|² below 1e-300 = 0) now pin the
  fixed values.
- **Result, every target**: 9 modules, 0 failing; **184,863 ok, 0 failed, 10,182 accepted (ADR 01: 6,505; ADR 0002:
  3,677)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64. Ten changed, those of the functions above;
  `ComplexExt::norm`'s equals `abs`'s.

| Function | Digest v2, every target |
|---|---|
| Complex::abs, ComplexExt::norm | `3dc86167dcf02563` |
| Complex::norm_sq | `5d2fbc215ae7cb6c` |
| Complex::mul | `47b801ceee8d67de` |
| Complex::div | `510e32b88696b639` |
| Complex::inv | `3de67cb534ca7167` |
| Complex::is_zero | `1e4707ccf49d3a28` |
| Complex::sqrt | `1443f1f8570fefbd` |
| Complex::ln | `33c41c5da9de3624` |
| Complex::powf | `78c102e1f284216d` |

- **Benchmarks** (Cyrius-only; `cyrius bench` with no other bench running): `complex_mul_div` (a product and a
  quotient by the operators) went from 36 ns to 84 ns, and `complex_exp` (386 ns) and `complex_sin` (837 ns) are
  unchanged.
  - **The guard costs about 9 ns per operation.** It reads the four biased exponents and folds the eight bound tests
    into one sign test of an OR. Without the operators the same product and quotient take 53 ns, against 34 ns
    without the guards (a scratch bench).
  - **The remaining 31 ns are an open finding.** The operator expression `a * b / b` costs 90 ns where the
    equivalent calls cost 53 ns in the same binary, a gap the old bodies didn't have. A two-field struct whose
    operator fn has the same branch and side call shows no gap (11 ns either way), so the cause isn't yet isolated
    enough to file.
  - **Two attempts changed nothing:** moving the edge paths out of the operator fns, and reading the parts in place.
  - Speed is measured, not gated.

## 19. M1 part 5c-2: BjComplex's transcendental range edges (2026-11-12)

- **ADR 0002 row 16** carries row 14's rule into the chain's transcendental steps. Where e^re (`exp`), cosh and
  sinh of im (`sin`, `cos`), ln or pow of |z|, the power |z|ⁿ, or an angle (`sqrt`'s arg z / 2, `powf`'s n·arg z)
  is not a normal double, it is taken to 53 bits with an unbounded exponent and each final product rounded once.
  - **The guards**: `exp` off re ∈ [−708, 709]; `sin` and `cos` past |im| = 710; `ln` at row 14's edge where the
    unbounded |z| is outside the normal range; `sqrt` where im ≠ 0 and |arg z| < 2^-1020; `powf` wherever |z|, |z|ⁿ,
    arg z or n·arg z is not a normal double above 2^-1022 (a binary64 value of exactly 2^-1022 may be the rounding of
    a smaller 53-bit one, so it takes the edge too). A non-finite part or exponent keeps hisab's chain.
  - **A fix to 5c-1**: `sqrt`'s edge path used the binary64 arg z / 2, which loses bits below 2^-1022. √(DBL_MAX + j)
    came out as 2^-513 less an ulp, where 1/(2√DBL_MAX) rounds to 2^-513 exactly.
- **`src/crmath_ux.cyr`** (new; in the bundle after crmath_mono) computes the steps on crmath's 364-bit fixed point
  (impetus's slow paths, error below 2^-290):
  - e^x by `_bj_cr_exp_fx_parts`;
  - cosh and sinh past 710 as e^|x|/2, since e^-|x| is below 2^-2048 of it;
  - ln(s·2^E) by `_bj_cr_log_fx` plus E·ln 2;
  - (s·2^E)^y as crmath's pow slow path, with its exactness test given an unbounded exponent;
  - sin and cos of a 53-bit angle past DBL_MAX by Payne–Hanek against the same 1,568 bits of 2/π;
  - `_bj_cr_round_ux` for the 53-bit rounding, with no exponent range.

  Two of these are bijli-side copies of crmath functions taking an unbounded input: the exactness test
  (`_bj_ux_pow_exact`) and the reduction (`_bj_ux_reduce`). src/crmath.cyr is impetus's copy, kept in step with it.
  - **Midpoints**: an exact dyadic power can be a 53-bit midpoint. For example, (2^18 − 1)²·2^700 to the power 1.5
    is (2^18 − 1)³·2^1050, and (2^18 − 1)³ is an odd 54-bit integer. With the exactness test disabled, `cr_verify`'s
    pow family shows 65 mismatches in 2,000 calls; with it, none.
  - **The stand-in**: far past the range a value is (1, ±2^20) in the exponent, whose products overflow or vanish
    exactly as the true value's do. For e^x, cosh and sinh that is past |x| = 1500, since their factors are doubles.
    For a power it is past |y·ln x| = 10000. Its factor can be the sine of an angle far below 2^-1074, but never
    below 2^-2098·|n|.
- **The model** (`validate.py`'s `cx_model`) is the whole unbounded chain. Its transcendental steps are MPFR at 53
  bits with its widest exponent range (`mpfr_ref.wide`, new), with Fractions throughout. `cx_governed` mirrors the
  port's guards, and on every row the port keeps hisab's chain for, the model must equal that chain. That is the
  guard's proof, now for exp, sin and cos too.
- **The edge grid**: 910 rows appended to `complex.tsv` (`m_complex.rs`'s `edges`):
  - exp at re ∈ ±[708, 1e10] across a dozen angles;
  - sin and cos at |im| ∈ [710, 1e300];
  - ln and sqrt of 20 moduli past either end or with a tiny angle;
  - 202 powf rows: moduli and powers past either end, the two 54-bit midpoints, n·arg z past DBL_MAX, and n below
    2^-1022.
- **`complex.adr.tsv`: 1,987 rows**, all passed. Row 14 has 1,596, row 15 has 92, and row 16 has 299 (sin 82, cos 82,
  powf 54, exp 51, ln 22, sqrt 8).
- **Tests**: `tests/complex.tcyr`'s `cyrius_transcendental_edges (ADR 0002 row 16)`, 13 assertions on the model's
  values:
  - e^710·(1, 0) = (∞, +0), and e^1000·sin 2^-1074 = 9.733444573000164e110;
  - e^-740·cos π on the subnormal grid, rounded once;
  - sin(j800) = (0, ∞), and sin(−0 + j1e5) = (−0, ∞) past the stand-in;
  - cos(π/2 + j720) = (1.5065301609522464e296, −∞);
  - ln|DBL_MAX·(1 + j)| = 710.1292864836639, and ln|2^-1074·(1 + j)| = −744.0934983311013;
  - √(DBL_MAX + j) with im 2^-513;
  - (1e200)² = (∞, +0), and the 53-bit midpoint (2^18 − 1)³·2^1050 broken to even;
  - (−1)^1e308 finite, and (2^-1074·(3 + j4))^0.5 = 2^-536 + j2^-537.
- **`cr_verify.py`'s ux families** (`cr_verify.py . --n 2000000 --aarch64 ux`): 1,087,540 calls on x86_64 and on
  aarch64, 0 mismatches against MPFR at 53 bits with its widest exponent range. They cover 687,540 arguments: e^x,
  cosh, sinh, ln and pow 100,000 each, sin and cos 93,770 each, and each (s, e) result is read in two calls. The
  arguments include binary64's ends, the stand-ins, powers of two and near-one values for ln, and for pow exact
  squares and fourth powers, both 54-bit midpoint families, near-one powers, and powers up to |y·ln x| = 9900. The
  angles include ones of the form RN53(arg·n) past DBL_MAX.
- **Result, every target**: 9 modules, 0 failing; **185,378 ok, 0 failed, 10,577 accepted (ADR 01: 6,505; ADR 0002:
  4,072)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and `tests/complex.tcyr`'s 52 assertions
  pass on all three. Six digests changed, those of the functions row 16 touches, whose rows the edge grid also
  extends.

| Function | Digest v2, every target |
|---|---|
| Complex::exp | `65d6dc5819b417b8` |
| Complex::ln | `47ef68bf6eba60fd` |
| Complex::sqrt | `07394dbe48c66d2d` |
| Complex::sin | `188c17d4248e51c1` |
| Complex::cos | `3427714b3af83662` |
| Complex::powf | `4f6fcb559a507331` |
- **Benchmarks** (Cyrius-only). The machine is never free of other repos' gates now, so this is an interleaved A/B:
  HEAD (an export of 6685842) then this tree, twice over, at a load average of 2.9 to 4.1.
  - `complex_exp` costs about 10 ns more, 395/399 ns against 405/408 ns (≈ 2.5 %): the guard's two compares.
  - `complex_sin` is unchanged (867 ns against 867/869 ns; HEAD's second round read 959 ns, a load spike), and so is
    `complex_mul_div` (86 ns against 87 ns).
  - The edge paths are new benches: `complex_exp_edge` (exp(1000 + j0.5)) takes 30.7 µs, and `complex_powf_edge`
    (|z|^2.5 past DBL_MAX, on the fixed point) 66–69 µs. Only the edges reach them.
  - The twins' record (`docs/benchmarks-rust-v-cyrius.md`) keeps its earlier column. This run read every twin 4–8 %
    slower, and HEAD read the same under the same load, while nothing it covers changed.
  - Speed is measured, not gated.

## 20. M1 part 5d-1: the exact Singularity threshold (2026-11-12)

- **ADR 0002 row 17**, the user's "exact thresholds" decision. Rust tests the rounded r² against the double nearest
  10^-30, which lies 0.37 ulp above 10^-30, and the Green's functions test the rounded R the same way. The port tests
  the exact Σ (aᵢ − bᵢ)² of the inputs against 10^-30, or against 10^-60 for R. That covers `electric_field_point_charge`,
  `electric_potential_point_charge`, `magnetic_field_moving_charge` (and through them superposition and the
  `PointCharge` methods), the three dipole functions, the quadrupole, both Green's functions and `coulomb_force`.
  - **Radial distances need nothing.** A distance is a double, and the double 1e-30's predecessor lies below 10^-30,
    so `distance < 1e-30` is exact; `validate.py` checks the claim.
- **`src/crmath_cmp.cyr`** (new; an internals file, in the bundle after crmath_ux):
  - **The band test** (`_bj_x_band`): the rounded value decides when it is more than 64 ulps from the cutoff (its
    error is a few ulps), and a NaN keeps Rust's false.
  - **The exact test** (`_bj_x_r2_below`), inside the band only: each difference is s + t by TwoSum, so r² is nine
    products of two doubles. Times 10^p their sum is compared with 1 in 2,560-bit integers (80 limbs of 32 bits),
    which hold the widest spread the band allows (about 2,300 bits, from a subnormal t's square).
  - A first version named its integers `_bj_bn_*`, which `fmtf.cyr` already defines. The later definitions won,
    with a warning, and the harness segfaulted. They are `_bj_xn_*` now.
- **The knife-edge grid** (`grids.rs`'s `knife_edge`, selected in double-double with 10^-30 from 5^30 in a u128):
  field points along four directions stepped by ulps, against three sources (the origin; −3·10^-32·s, and 2^-1070,
  whose differences round). It keeps the 421 pairs whose exact r² is within an ulp of 10^-30, and 781 within two ulps
  of 10^-60 for R. A Python search over 3,185 candidates found 13 flips, all inside a one-ulp window.
  - The rows: 3,405 in `field.tsv`, through every threshold function (the dipoles and the quadrupole only on the
    145 pairs about the origin), and 421 `coulomb_force` rows in `charge.tsv`.
- **The model**: `validate.py`'s transcriptions test through `cut_r2` and `cut_r`. With `EXACT_THR` off they are Rust's
  comparisons, which `validate_field` and `validate_charge` check against the oracle. `adr17_rows` turns it on, maps
  each row whose outcome changes, and writes the port's outputs (MPFR's primitives for the Green's function's cos and
  sin). The corrected oracle skips any row that row 17 makes an error.
  - **302 mapped rows**, all passed: field 263 (point E, V and B 39 each, the dipoles and the quadrupole 13 each, the
    two Green's functions 47 each) and charge 39 (`coulomb_force`). 202 go from Ok to Singularity and 100 from
    Singularity to Ok.
- **`cr_verify.py`'s `r2_below` family** (an exact family, against Fractions): 200,000 calls on x86_64 and on aarch64,
  0 mismatches. 55 % lie within 2^-46 of the cutoff, and the points include rounding differences with tiny and
  subnormal partners, nearly equal large points, and p = 30 and 60.
- **Tests**: `tests/field.tcyr`'s `cyrius_exact_threshold (ADR 0002 row 17)`, 7 assertions on golden knife-edge points:
  - an r² that rounds to the double 1e-30 but lies below 10^-30 is Singularity;
  - a source 16·2^-1074 from the origin decides through the exact difference;
  - an r² that rounds below 1e-30 but lies above 10^-30 gives Rust's value, for both the potential and the field;
  - R just above 10^-30 gives Rust's 1/(4πR);
  - r = 0 stays Singularity, and a NaN point keeps Rust's test.
- **Result, every target**: 9 modules, 0 failing; **188,902 ok, 0 failed, 10,879 accepted (ADR 01: 6,505; ADR 0002:
  4,374)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and field's 152 and charge's 66 assertions
  pass on all three. The ten threshold functions' digests changed with their new rows.
- **Benchmarks** (interleaved A/B against HEAD fbfe263, twice, load average 2.4 to 3.1):
  - `field/electric_field_point` +2–3 ns (30/32 → 33/34 ns) and `electric_potential_point` +4 ns (18/19 → 22/23);
  - the two dipole fields +4 ns (83/86 → 87/90);
  - `magnetic_field_moving`, `coulomb_force` and `superposition_10` within their noise (+0–3 ns, +1–12 ns over 10
    charges).

  The cost is the band test's call: its body is over `#inline`'s 32 tokens. A first version called it through a
  five-argument wrapper and cost 5–9 ns. Speed is measured, not gated.

## 21. M1 part 5d-2: the vector chains' range edges (2026-11-12)

- **ADR 0002 row 18**: row 14's rule for field.rs's and charge.rs's vector chains: `FieldVector::magnitude`,
  `magnitude_sq`, `dot`, `cross` and `normalized`, the Lorentz force (and `PointCharge::lorentz_force_in`), and the
  plane and dielectric images. At an edge the chain runs with an unbounded exponent, and its last operation rounds once.
  - **`normalized`** at an edge is its first branch, |v| to 53 bits and then each quotient, where Rust rescaled by
    its largest part with two more roundings. The result is exactly scale-invariant: 2^1074·v and v give the same
    bits wherever the parts stay exact. Rust's rescale branch is gone; only the zero vector and a non-finite part
    reach the code after the fast test.
  - **The images**: the plane image's 2·z_p is exact unless |z_p| ≥ 2^1023. The dielectric image checks its product
    −q·(ε₂ − ε₁) and its sum ε₁ + ε₂ after Rust computes them.
- **The guard** (`src/crmath_vx.cyr`) tests what Rust's chain computed. A product below 2^-1022 differs from its
  53-bit value, but it lies under half an ulp of any operand of 2^-968 or more, and a sum of such small values is
  absorbed in turn by any term of 2^-913 or more. So a final sum or difference S with |S| in [2^-900, DBL_MAX] is
  the unbounded chain's, an overflow shows as a non-finite S, and an exact 0 is the chain's when every product is
  exact (a zero factor, or normal). In code, `cross` and the Lorentz force fold the three results' biased exponents
  into one sign test ([123, 2046] is the band), then test components that are 0 through zero factors. A non-finite
  part, anywhere in a call, keeps Rust's chain.
  - **Two guards were tried first.** The parts in [2^-511, 2^511), as row 14's, called through a helper: magnitude
    9 → 19 ns, cross 19 → 44 ns, the Lorentz force 61 → 216 ns. Inline and with zeros folded in, it was still 18 ns.
    A scratch bench put the cost in the guard's integer work itself (Cyrius keeps every local in memory), and two
    compares on the sum Rust computes anyway cost 2 ns.
  - **A per-component non-finite test** was also wrong at first. A NaN in E's y left the z component on the
    unbounded path, and the harness's 13 failures showed it. The rule is per call, as for complex.
  - The dipoles, the stress tensor and the Poynting flux call Rust's raw dot and cross (`_bj_fv_dot_rust`,
    `_bj_fv_cross_rust`) until their own parts of 5d. The dipoles had picked up the new dot at once (4 harness
    failures).
- **The model** (`validate.py`'s `fv_model`, `fv_governed`, `fv_rows`): the exact chain rounded to 53 bits (`rnd53`),
  the last operation once. `fv_governed` mirrors the port's guards on Rust's binary64 values, and every row they leave
  on Rust's chain must equal Rust's transcription: the guard's proof, which held for the band argument on every row.
- **The edge grid**: 455 field rows over `m_field.rs`'s 13 `edge_vectors` (cancellations past DBL_MAX, the subnormal
  grid, three squares of 2.25·2^-1074 that Rust rounds to 2 each), the plane image at |z_p| ≥ 2^1023, and the
  dielectric image where its product or sum leaves the range; 65 Lorentz rows in `charge.tsv`.
- **294 mapped rows**, all passed: magnitude 145, cross 55, dot 28, normalized 11, the plane and dielectric images 6
  each, magnitude_sq 1 (7·2^-1074 for Rust's 6), the Lorentz force 40 and `lorentz_force_in` 2.
- **Tests**: `tests/field.tcyr`'s `cyrius_vector_edges (ADR 0002 row 18)`, 8 assertions:
  - |2^600·(3, 4, 12)| = 13·2^600;
  - the three subnormal squares summed exactly;
  - a cross component 2^1100 − 2^1100 = 0, and a dot product past DBL_MAX still ∞;
  - normalized on the subnormal grid, and its scale invariance against Rust's in-range bits;
  - the plane image at z_p = DBL_MAX, and the dielectric image with ε₁ + ε₂ past DBL_MAX.
- **Result, every target**: 9 modules, 0 failing; **189,128 ok, 0 failed, 11,173 accepted (ADR 01: 6,505; ADR 0002:
  4,668)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and field's 160 and charge's 66 assertions
  pass on all three. The nine functions' digests changed.
- **Benchmarks** (interleaved A/B against HEAD 16d9c20, twice, load average 3.6 to 5.3):
  - `field/vector_magnitude` +1–2 ns (10/11 → 12 ns) and `vector_cross` +4–7 ns (19/22 → 26/27 ns);
  - `charge/lorentz_force` −10 to −14 ns (65/66 → 52/55 ns): the inline chain no longer calls three struct-valued
    methods (cross, add, scale);
  - the dipoles and `trace_field_line_10` within the load's noise.

## 22. M1 part 5d-3: the point charges' and Green's functions' range edges (2026-11-12)

- **ADR 0002 row 18, amended.** The point charge's E and V, the moving charge's B, superposition, `coulomb_force`, the
  `PointCharge` methods through them, and both Green's functions now run with an unbounded exponent at a range edge,
  each output's last product rounding once. Typical edges:
  - |E| = kq/r² is representable while the factor kq/r³ under- or overflows (q = 1.1e90 at r = 1e200 gives
    9.886306971487286e-301, where Rust's factor underflowed);
  - k·q₁·q₂ overflows while the force does not (1e200·1e200 at 1e100 gives 8.987551792261174e209);
  - two fields past DBL_MAX cancel in a superposition (0, where Rust has ∞ − ∞ = NaN);
  - 4πR passes DBL_MAX (1/(4π·1e308) = 7.95774715459477e-310, where Rust has 1/∞ = 0);
  - kR passes DBL_MAX (kR = 1e600 gives a finite e^{−jkR}/(4πR), where Rust has NaN).
- **The guards** test Rust's own values:
  - Past the Singularity cutoff r² ≥ 10^-30, so r² finite keeps the squares' sum exact by crmath_vx's band argument.
  - A chain of products by 1/r is the unbounded one when its first and last products are normal (biased exponent in
    [2, 2046]), since the ones between lie between them; an overflowed offset makes the factor 0, so E and
    `coulomb_force` need no r² test.
  - B adds its unit-vector parts and cross components, folded branch-free with exact zeros masked in.
  - Superposition needs each field component normal (or 0 with a zero offset part) and the totals finite.
  - The Green's functions test 4πR, its reciprocal and kR.
  - **A zero charge** keeps Rust's chain only where r² is finite. Otherwise Rust's 1/r is 0, its r̂ and cross
    components are signed zeros, and its 0·∞ is NaN. `validate.py`'s guard proof caught that on 2 rows of a moving
    zero charge at DBL_MAX (−0 where the chain gives +0).
- **The shared helpers** (`_bj_vx_offset_ux`, `_bj_vx_coulomb_ux`, `_bj_vx_coulomb_edge`, `_bj_vx_coulomb_plain`,
  `_bj_vx_finite_n`) moved into `src/crmath_vx.cyr`. charge.cyr can't reach a private field.cyr function (the
  compiler says so). The 53-bit angle helpers moved there from complex.cyr too (`_bj_vx_cos53`, `_bj_vx_sin53`).
- **A longer 2/π.** kR reaches 2^2050, past 5c-2's reduction, which reads impetus's 56 chunks of 2/π and stops at
  2^1026. The first run gave NaN on the 4 corpus rows with k = 1e300 and R ≥ 1e300.
  - `src/crmath_ux_data.cyr` (generated: `crtables.py`'s `render_ux`, from π to 2,800 bits checked against MPFR) holds
    96 chunks, whose first 56 equal impetus's, and `_bj_ux_reduce` reads them up to index 88.
  - `cr_verify.py`'s angle families now reach te = 2052 and RN53(k·R), and the whole ux set passes again: 1,093,984
    calls on x86_64 and on aarch64, 0 mismatches.
- **The model** (`validate.py`'s `pc_plain`, `pc_model`): the exact threshold, then the exact chain rounded to 53 bits
  (MPFR's widest range for the Green's cosine and sine), mirrored guards, and the guard's proof on every row the guards
  leave on Rust's chain. That proof compares against Rust's chain with the exact threshold and MPFR's primitives, and
  skips the rows row 17 maps. The corrected oracle skips rows row 18 governs.
- **The edge grid** (`m_field.rs`'s `emit_5d3`, `m_charge.rs`), 263 rows:
  - 8 charges over 9 point pairs through E, V and B;
  - 5 superpositions;
  - 6 Green's positions with 5 wavenumbers;
  - 6 `coulomb_force` pairs.
- **251 more mapped rows** (row 18: 545), all passed: point E 63 and `electric_field_at` 24, moving B 50 and
  `magnetic_field_at` 4, V 32 and `electric_potential_at` 4, `coulomb_force` 36, the scalar Green's function 24, the
  static 9, superposition 5.
- **Tests**: `tests/field.tcyr`'s `cyrius_point_charge_edges (ADR 0002 row 18)`, 8 assertions: the six edges above
  with their bits, and the two charges' construction.
- **Result, every target**: 9 modules, 0 failing; **189,140 ok, 0 failed, 11,424 accepted (ADR 01: 6,505; ADR 0002:
  4,919)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and field's 167 and charge's 66 assertions
  pass on all three.
- **Benchmarks** (interleaved A/B against HEAD 72e1a56, twice, load average 3.5 to 4.2):
  - `field/electric_field_point` +4 ns (35/36 → 39/40 ns) and `electric_potential_point` +3 ns (24 → 27);
  - `charge/coulomb_force` +6–8 ns (46/48 → 54);
  - `superposition_10` +24–44 ns (569/583 → 607/613);
  - `trace_field_line_10` within noise;
  - **`magnetic_field_moving` +38–39 ns (45/46 → 84/85).** Its bench geometry is axis-aligned, so two r̂ parts and two
    cross components are exact zeros. The first guard sent those through a second tier of compares and calls (163 ns,
    then 98 ns inline); the branch-free fold with masked zeros gets 85. Cyrius spends about a nanosecond per integer
    operation here (every local is in memory), and the guard needs about 40. This is a finding with its numbers; speed
    is measured, not gated.

## 23. M1 part 5d-4: the dipoles', the quadrupole's, the axes' and the sphere image's range edges (2026-11-12)

- **ADR 0002 row 18, amended**: the electric and magnetic dipole fields, the dipole potential, the quadrupole potential,
  the ring and disk axial fields and the sphere image run with an unbounded exponent at a range edge. Typical edges:
  - a DBL_MAX dipole's E_z is −∞ where Rust had ∞ − ∞ = NaN;
  - a subnormal dipole's potential is −1.2e-312 for Rust's −0;
  - the ring axis at z = 2^361 is kQz/z³ where Rust's z³ overflowed to give 0;
  - the disk axis at R = z = 2^-1022 is finite where Rust's was ∞;
  - a quadrupole whose terms overflow and cancel is +∞ where Rust's was NaN.
- **The guards** test Rust's own values.
  - **The dipoles**: r̂, m·r̂ (in the band), (1/r)³, the factor, 3(m·r̂), k and each output component (in the band) in
    one branch-free fold. An exact 0 passes: r̂ᵢ for pᵢ = 0, m·r̂, 3(m·r̂) and k when every product has a zero factor,
    and an output when both its products do.
  - **The quadrupole**: r², r⁵ and 8πε₀r⁵ (r⁴ lies between), each Qᵢⱼxᵢ and each term normal or an exact 0, and the
    sum finite.
  - **The ring**: z² + R² in the band, then (z² + R²)^(3/2), k·Q and k·Q·z.
  - **The disk**: z² + R² in the band, r + |z| finite, then R/r, R/(r + |z|), their product and σ/(2ε₀).
  - **The sphere image**: |r|² in the band (which makes Rust's |r|, and so its `|r| > R` decision, the unbounded
    chain's), then R/|r| and its square. Outside the band the comparison itself runs on the unbounded |r|, and its
    error reports |r| rounded once.
- **Two slips on the way.**
  - The first dipole guard read its values as `load64(h + 8)` on a typed `*f64`, which advances 8 elements, not 8
    bytes. Every dipole took the unbounded path (819 ns to 1.2 µs on the bench) until a probe printed the fold.
  - The guard then cost 143–150 ns as 18 small calls, and 118 ns inline in the shared dipole function.
  - `fv_rows` also compared the model's decoded error text with the oracle's still-encoded `s:` token and wrote one
    identical mapping; it decodes now.
- **The edge grid** (`m_field.rs`'s `emit_5d4`), 88 rows:
  - 4 moments over 5 field points through the three dipole functions;
  - 3 quadrupole tensors over 4 points;
  - 5 axis triples;
  - 6 sphere images, among them a charge at 10^-200 outside a sphere of 10^-210, whose |r|² underflowed in Rust to
    0, "inside".
- **251 more mapped rows** (row 18: 796), all passed: disk axis 124, ring axis 45, the electric dipole field 29, the
  magnetic 28, the dipole potential 13, the sphere image 6 (one Ok where Rust refused, two refused where it
  accepted), the quadrupole 6.
- **Tests**: `tests/field.tcyr`'s `cyrius_dipole_axis_image_edges (ADR 0002 row 18)`, 6 assertions on mapped rows'
  values.
- **Result, every target**: 9 modules, 0 failing; **188,977 ok, 0 failed, 11,675 accepted (ADR 01: 6,505; ADR 0002:
  5,170)**. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and field's 173 assertions pass on all
  three.
- **Benchmarks** (interleaved A/B against HEAD 1301d63, twice, load average 4.5 to 5.0):
  - both dipole fields +18–20 ns (96/98 → 114/118 ns);
  - `field/ring_axis` +8–9 ns (28 → 36/37) and `field/disk_axis` +10–14 ns (42/44 → 54/56);
  - `trace_field_line_10` within noise.

  The dipole guard reads eleven values and six zero bits, about 70 integer operations. Speed is measured, not gated.

## 24. M1 part 5d-5: the stress tensor's, the force's, the Poynting flux's, the multipoles' and the field line's range edges (2026-11-12)

- **ADR 0002 row 18, amended** (the last bite: Part 5 and M1 are complete): `maxwell_stress_tensor`,
  `stress_tensor_force`, `poynting_flux_surface`, both multipole moments and `trace_field_line`'s step run with an
  unbounded exponent at a range edge. Typical edges:
  - the stress tensor of E = (1e155, 1e155, 0), B = x̂ is (3.98e5, 8.85e298, 0, …) where Rust's E² overflowed and
    eight elements were NaN (the off-diagonal ½·0·∞ among them);
  - a force row (1e200, −1e200, 0)·(1e200, 1e200, 0) is +0 for NaN;
  - E × B = 10^400 ẑ through n̂ = 10^-200 ẑ and ΔA = 10^-200 gives 7.96e5 W for ∞;
  - a 1e300 C charge at z = 1e-170 has M₁₀ = 1e130: Rust's z² underflowed, so r = 0 put it at the origin;
  - P₁₇₀¹⁷⁰ past DBL_MAX times rˡ below 2^-1074 is 7.0e-103 for NaN;
  - a field line past a 1e290 C charge traces its 4 steps, where Rust's kq/r³ = ∞ made ∞·0 = NaN and refused it.
- **The guards** test Rust's own values.
  - **The stress tensor**: a box. Every part of E and B ±0 or with its biased exponent in [773, 1273) (2^±250) keeps
    each product, E², B², their halves, each difference (exact when it is that small) and its products by ε₀ and 1/μ₀
    normal or an exact 0. Branch-free: per part, a word that is negative exactly in the box or at 0.
  - **The force and the flux**: sums of products, so src/crmath_vx.cyr's band applies at each level. A force component,
    each cross component, each dot product and the flux's total in [2^-900, DBL_MAX] is the unbounded chain's, as is
    an exact 0 whose products are exact. A non-finite part anywhere in the call keeps Rust's chain.
  - **The multipoles**: no sum to absorb an underflow, so the checked twins (`_bj_spherical_coords`,
    `_bj_assoc_legendre_k`, `_bj_norm_legendre_k`) clear a flag unless every product and quotient is normal above
    2^-1022 or an exact 0 from a zero operand. rˡ needs only its own test: with r normal, powi's intermediates lie
    between 1 and rˡ. Rust's origin test must agree (r = 0 only when every coordinate is 0), and φ = atan2(y, x)
    below 2^-1021 with y ≠ 0 fails (the unbounded chain takes y/x to 53 bits there, which atan2 rounds alike).
  - **The field line's step**: |F|² in the band (normalized's first branch), h < 2^1023, and per part F_i = 0 or
    v_i normal, v_i = h·F_i/|F| for h ≤ 1 (|F_i/|F|| ≥ |v_i| then) and F_i/|F| for h > 1 (then |h·F_i/|F|| ≥ it, and
    |F_i/|F|| ≤ 1 + 2^-51 keeps it finite). The trace's field function decides its own chain.
- **The unbounded twins**: `_bj_ux_sph` (r, cos θ, sin θ and φ with Rust's clamps), `_bj_ux_radial` (powi's
  square-and-multiply to 53 bits, `_bj_ux_pow` past 2^31 − 1), `_bj_ux_assoc_legendre`, `_bj_ux_norm_legendre`, and
  each moment's sum with its last addition pending (`_bj_ux_sum_add`, `_bj_ux_sum_fin`), so it rounds once.
- **validate.py**: `m5_plain` mirrors each guard on Rust's values and `m5_model` is the whole unbounded chain. The
  trace runs in a port mode: each step's field is the field function's row-18 model where its own guard leaves Rust's
  chain, each step the unbounded one, and on every step the port keeps, the unbounded step must equal Rust's (the
  guard's proof, step by step).
- **The edge grid** (`m_field.rs`'s `emit_5d5`), 204 rows: 8 (E, B) pairs; 4 tensors over 4 normals; 6 flux sample
  sets; 8 charge sets over 10 (l, m) through both moments; 14 field lines (a uniform field on the subnormal grid or past
  either end, steps of 5e-324 to 1e308, and the point charge, dipoles and quadrupole far out).
- **166 more mapped rows** (row 18: 962), all passed: the stress tensor 70, the real multipole moment 54, the
  spherical 26, the force 7, the Poynting flux 5, the field line 4. 97 of them are existing rows of `field.tsv`
  (stress 64, multipole 30, force 3), now read at the edge.
- **Tests**: `tests/field.tcyr`'s `cyrius_stress_multipole_trace_edges (ADR 0002 row 18)`, 9 assertions on mapped
  rows' values.
- **Result, every target**: 9 modules, 0 failing; **189,015 ok, 0 failed, 11,841 accepted (ADR 01: 6,505; ADR 0002:
  5,336)** over 200,856 rows. All 145 digests are equal on x86_64 ELF, x86_64 PE and aarch64, and field's 182
  assertions pass on all three.
- **Benchmarks** (interleaved A/B against HEAD dacd8b9, twice, load average 1.9 to 2.3, another repo's validate.py
  on one core):
  - `trace_field_line_10` +5 % (1.128/1.131 → 1.182/1.194 µs), the step guard's ~5.5 ns a step. Its first version
    (both F/|F| and h·F/|F| tested per part, ~60 integer operations) cost +13 %. A third pair at load 3.1–3.5 read
    +7 % and +11 %, with untouched twins (the dipole field, the disk axis) moving up to +21 % in one round of it;
  - every other twin within ±5 % at load 2.2, with no consistent sign.

  Ad hoc, with no Rust twin (`tests/zz_mp.bcyr` in the scratch trees, not committed): the stress tensor 142 → 170 ns,
  the force 60 → 72 ns, the flux over 3 samples 229 → 294 ns, a real l = 4 moment of 5 charges 2.22 → 2.75 µs (+24 %),
  the spherical 2.63 → 3.24 µs. The force's and the flux's first guards re-ran every product through a check
  (172 and 645 ns); the band made them a test per sum. The multipoles' cost is the per-product check, about 20 a
  charge. Speed is measured, not gated.

## 25. M2 part 1a: material (2026-11-12)

- **`src/material.cyr`**: material.rs's 21 free functions, `BjMagneticType` (append-only values 0–2) and `BjMaterial`
  (24 bytes, by value: `BjMaterial_new` / `vacuum` / `dielectric` / `conductor` / `default`, and the methods with
  Rust's names). Every formula is Rust's in Rust's order. Rust's `pub(crate) lossless_optical_constants` is the
  public internal `_bj_material_lossless_optical_constants`, for wave.cyr's Fresnel functions (part 1b): a file's
  non-public functions are private to it in the tests' build. The deprecated `reflectance_at` is not carried (ADR
  0002 row 1).
- **ADR 0002 rows 19–21**, the user's answers of 2026-11-12:
  - **19**: the Curie laws refuse a NaN temperature (`!(T > 0)`, `!(T − T_c > 0)`); Rust's `<= 0` let NaN through as
    `Ok(NaN)`. 597 rows mapped.
  - **20**: `BjMaterial.wave_speed` validates as `impedance` does: InvalidPermittivity / InvalidPermeability with the
    value, for a non-finite or non-positive ε_r or μ_r (Rust: InvalidParameter, and ∞ accepted). 1,421 rows mapped.
  - **21**: `bj_classify_magnetic` returns a Result: InvalidParameter for NaN (Rust: Paramagnetic). 1 row mapped.

  The fourth question, the two zero-permeability guards (|μ_r| < 1e-30 for H, |μ_r·μ₀| < 1e-30 for the energy), stays
  as Rust's: each guards its own division.
- **Row 18 through FieldVector**: the bound surface charge and current are `dot` and `cross`, so their range edges
  are 5d-2's: 29 rows mapped from `fv_model`.
- **The JSON**: `BjMaterial_to_json` / `BjMaterial_from_json` (D5) on json.cyr's shared reader. `json.tsv` grows by
  4,731 Material rows; `json.cr.tsv` by 2,509 rows serde_json's best-effort parse misreads (gated at 0 ulp, as the
  other readers).
- **New golden file `material.tsv`, 51,814 rows, class A**: the one-argument relations over 69 values (specials, the
  1e-15, 1e-30 and 1 thresholds with their neighbours, ε₀, μ₀, binades, seeded), crossed in pairs for the energy
  densities, Clausius–Mossotti and Curie's law; triples for Curie–Weiss and ρ_b; the Clausius–Mossotti edges X → 1 and
  X → −1/2 and the permeability threshold ulp by ulp; scale factors over a stride of field.rs's vectors; 2,268
  materials through every method; 31² material pairs for `try_reflectance_at`. Each value appears once: a row's
  (function, inputs) is its key in the mapping files. `validate.py`'s independent transcription agreed with every row
  on its first run.
- **Rust tests**: 39 of material.rs's 42 ported (`tests/material.tcyr`, 75 assertions with the Cyrius-only groups).
  The three `reflectance_at` tests are n/a in the ledger (ADR 0002 row 1); `cyrius_reflectance_via_try` checks their
  values through the replacement. The ledger: 199 ported, 3 n/a, 656 pending of 858.
- **Surface parity**: 6 modules ported, 151 items matched (44 new), 0 problems. The API snapshot gains 27 public fns.
- **Benchmarks** (the five material twins; 22 of 80 now twinned; load average 2.5–3.5, each the lower of two runs):

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `material/polarization` | 1.00 ns | 25 ns | 25.1× |
  | `material/displacement_field` | 0.99 ns | 23 ns | 23.2× |
  | `material/h_field_from_b` | 5.51 ns | 22 ns | 4.0× |
  | `material/clausius_mossotti` | 1.24 ns | 26 ns | 20.9× |
  | `material/curie_weiss` | 1.24 ns | 8 ns | 6.5× |

  Rust's arguments are constants that LLVM folds; the port pays a call, a by-value struct and the Result.
- **For M2's range-edge Part** (§0.1), from the audit:
  - `classify_magnetic`'s `x < −1e-15`: the double nearest 1e-15 lies 0.39 ulp above 10^-15, so x = −D is
    paramagnetic where the exact test says diamagnetic (an exact threshold, as row 17). `|μ_r| < 1e-30` on an input
    agrees with the exact test (D(1e-30) lies above 10^-30), but `|μ_r·μ₀| < 1e-30` tests a rounded product;
  - Clausius–Mossotti: N·α past DBL_MAX (X = −∞ gives Rust's ε_r = NaN where the limit is −2), and `X ≥ 1`, `ε_r ≤ 0`
    on the rounded X;
  - one-product chains to round once (row 13): the dielectric energy density, B²/(2μ), the scale factors times E or H,
    and BjMaterial's √(ε_r μ_r), impedance and wave speed (maxwell's twins round once since M1 part 5b).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **10 modules, 0 failing; 241,003 ok, 0 failed, 0
  measured, 0 pending, 16,398 accepted (ADR 01: 8,996; ADR 0002: 7,402)** over 257,401 rows; **all 179 digests
  equal**, the 145 of M1 unchanged.


## 26. M2 part 1b: wave (2026-11-12)

- **`src/wave.cyr`**: wave.rs's 42 public functions: the Poynting vector and the plane-wave relations, Snell's law and
  the angle-form Fresnel coefficients, the trig-free family (`*_direct`, unpolarized, Schlick), the material-based
  interface, the rectangular and circular waveguides and the antennas. Every formula is Rust's in Rust's order. The
  sin, cos, asin and atan are bijli's correctly rounded ones. Rust's `Result<Option<f64>>` (Snell, the critical angle,
  cos θ_t, the direct r_s and r_p) is `Ok(Some)` with the value in `*out`, or `Ok(None)` under total internal
  reflection, out untouched as on an error. A waveguide's `u32` mode indices are i64s with a Cyrius-only
  InvalidParameter outside a u32 ("m and n must be u32s, got m=-1, n=0"), field.cyr's precedent for the multipoles.
- **ADR 0002 rows 22–25**, the user's answers of 2026-11-12:
  - **22**: the dipoles on bijli's η₀ = √(μ₀/ε₀) with Cin(2π) exact (Rust's 1.6409, 73.1 Ω and 80π² assume
    η₀ = 120π): D = 1.640922376984585, R = 73.07901028545056 Ω, and (2π/3)·η₀ = 789.0221238669274 Ω in the
    Hertzian resistance. The constants are generated (`src/wave_data.cyr`, by `validate.py` in rational arithmetic).
    627 rows mapped.
  - **23**: `bj_effective_aperture` refuses a NaN wavelength (`!(λ ≥ 0)`), as row 19. 39 rows mapped.
  - **24**: `bj_schlick_reflectance` computes R₀ with no |n₁ + n₂| < 1e-30 guard: the index check already makes the sum
    positive. 76 rows mapped.
  - **25**: `bj_material_reflection_coefficient` takes Γ on the relative impedances √(μ_r/ε_r), so Γ² is
    `bj_reflectance_normal_materials`' R bit for bit (Rust's absolute impedances missed it on about 68 % of pairs). 74
    rows mapped.

  Reviewed and kept as Rust's: the zero guards on an input (`wave_number`, Snell, the critical and Brewster angles,
  the Hertzian resistance, Friis), as the user's material.rs answer; the angle forms' acceptance of any index (a
  negative index refracts negatively, which is physical); check_indices' order (n₂ = 0 is DivisionByZero before any
  other index check).
- **Row 18 through FieldVector**: the Poynting vector and the momentum density are `cross` before their scale
  factors: 180 rows mapped from `fv_model`.
- **`{ratio:.4}`**: `guide_wavelength`'s below-cutoff error prints the ratio as Rust's fixed precision does: the exact
  value rounded half to even (`1.03125` → `1.0312`), `-0.0000` for a negative zero, every integer digit of a ratio
  past 10^300. `src/fmtf.cyr` gains `_bj_fmt_f64_fixed` (0 ≤ p ≤ 17) on its bignum, with shift-right, bit and
  divide-by-small helpers; beam.rs's `{g:.4}` (M3) uses it too. `fmtf.tsv` gains 14,376 `f64::fmt_fixed` rows
  (precisions 0, 1 and 4 over a stride of the grid, the ties at each precision and the carries).
- **New golden file `wave.tsv`, 66,690 rows**: class C for the ten functions that call a transcendental (the plane
  wave, Snell, the critical and Brewster angles, the four angle-form Fresnel functions, the two patterns), class A for
  the rest. Pairs over 39 values (specials, the 1e-30 and 1 thresholds' neighbours, SI values, binades); Snell over
  index pairs and angles; the Fresnel angle forms over index × angle pairs and Snell-consistent physical calls; the
  trig-free family over valid index pairs × cosines through the clamp edges (−1e-6, 1 + 1e-6, each ± an ulp), invalid
  indices on a short list, and each dense → rare pair's critical cosine ulp by ulp; the patterns over [0, π] in 720
  steps with the 1e-15 axis cut ulp by ulp; the materials of material.tsv's reflectance pairs; the waveguides over
  dimensions, the u32 edge and the ties of `{ratio:.4}`. A first cut was 51 MB, mostly long error texts in crossed
  grids, so the grids are thinned to 19 MB with every edge kept. `wave.cr.tsv` holds the 3 rows where glibc
  misrounds (one asin, two cos). `validate.py`'s independent transcription agreed with every row on its first run.
- **Rust tests**: all 63 of wave.rs's tests ported (`tests/wave.tcyr`, 111 assertions with the Cyrius-only groups).
  `test_hertzian_radiation_resistance` compares against (2π/3)·η₀·10⁻⁴ (row 22). The ledger: 262 ported, 3 n/a, 593
  pending of 858.
- **Surface parity**: 7 modules ported, 193 items matched (42 new), 0 problems. The API snapshot gains 42 public fns.
- **Benchmarks** (the eleven wave twins; 33 of 80 now twinned; load average 2.3–3.3, each the lower of two runs, an
  A/B against HEAD showing no change in the 22 earlier rows):

  | Bench | Rust | Cyrius |
  |---|---|---|
  | `wave/poynting_vector` | 0.98 ns | 62 ns |
  | `wave/plane_wave_intensity` | 0.52 ns | 6 ns |
  | `wave/plane_wave_e` | 0.51 ns | **20.99 µs** |
  | `wave/fresnel_rs` | 1.25 ns | 275 ns |
  | `wave/snell_refraction` | 1.75 ns | 284 ns |
  | `wave/fresnel_rs_direct` | 1.73 ns | 58 ns |
  | `wave/reflectance_unpolarized` | 1.24 ns | 68 ns |
  | `wave/reflectance_at_interface` | 7.78 ns | 117 ns |
  | `wave/schlick_reflectance` | 1.24 ns | 44 ns |
  | `wave/half_wave_dipole_pattern` | 0.50 ns | 483 ns |
  | `wave/rect_waveguide_cutoff` | 1.24 ns | 47 ns |

  Rust's arguments are constants that LLVM folds, so its column is near a clock tick. **Finding**: `plane_wave_e`'s
  argument kx − ωt is about 5·10⁶ (k = 10⁷, x = 0.5), and crmath's sin and cos keep their fast path only for
  |x| ≤ 2^20: past it every call takes the multi-precision slow path, about 21 µs. An optical plane wave lands there
  routinely. Filed in the roadmap as a medium-range reduction (Cody–Waite with a longer π/2, or crmath_ux's
  Payne–Hanek, then the fast path's polynomials). Speed is measured, not gated: the port is unchanged. The other
  trigonometric rows (150–200 ns a call) are crmath's fast path.
- **For M2's range-edge Part** (§0.1), from the audit:
  - the trig-free family at extreme index ratios: (n₁/n₂)² past DBL_MAX gives sin²θ_t = ∞·0 = NaN at cos θ_i = 1, so
    r_s = NaN for valid indices (`fresnel_rs_direct(10^200, 10^-100, 1)` is Some(NaN), where r_s ≈ 1); n₁ cos θ_i ±
    n₂ cos θ_t past DBL_MAX gives ∞/∞. The ratios are scale-invariant, so a power-of-two scaling fixes both;
  - `sin²θ_t > 1`, `|n₁ sin θ_i/n₂| > 1`, `n₂/n₁ ≥ 1` and `1 − (λ/λ_c)² ≤ 0` decided on rounded values (TIR or
    cutoff within an ulp of the boundary; continuous in R, so low priority); the |den| < 1e-30 guards of the angle
    forms and `reflectance_normal` on rounded sums;
  - `reflectance_normal` and Schlick's R₀ with n₁ + n₂ past DBL_MAX (DBL_MAX and DBL_MAX/2 give R = 0 for 1/9);
  - the material interface: √(ε_r μ_r) past DBL_MAX or below 2^-1074 turns valid materials into an index error
    ("got n₁=1, n₂=inf"), and √(μ_r/ε_r) under- or overflows;
  - the waveguides: (m/a)² + (n/b)² under- or overflows (sides of 10^300 m are "not representable" where
    λ_c = 1.41·10^300 is; the cutoff frequency is 0 there), a hypot-style scaling; x·v and 2πa in the circular cutoff;
  - round-once candidates (row 13): the Poynting vector's and momentum density's scale factors after the cross
    (E × B past DBL_MAX makes g = ∞ where S/c² is representable), the plane-wave intensity, k = 2πf/v, the
    effective aperture, the Hertzian resistance and Friis's chains;
  - accuracy, not range: `half_wave_dipole_pattern` near the axis computes cos((π/2) cos θ) by cancellation, a relative
    error of about 2^-53/(πθ²/4) (3 % at θ = 10^-7, O(1) below 10^-8); the identity cos((π/2) cos θ) =
    sin(π sin²(θ/2)) has none. A question for the user in that Part.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **11 modules, 0 failing; 321,070 ok, 0 failed, 0
  measured, 0 pending, 17,397 accepted (ADR 01: 8,999; ADR 0002: 8,398)** over 338,467 rows; **all 222 digests
  equal**, the 179 before unchanged.

## 27. M2 part 1c: circuit (2026-11-12)

- **`src/circuit.cyr`**: circuit.rs's 27 functions: Ohm's law, the series and parallel combinations, capacitance and
  inductance, the RC and RL transients and the series RLC relations. Every formula is Rust's in Rust's order; the
  transients take bijli's correctly rounded exp. Rust's `&[f64]` is a pointer and a count, and a series sum is Rust's
  `Sum`, a left fold from −0.0 (Σ[−0] = −0, probe p11).
- **ADR 0002 row 26**, rows 19 and 23's rule: `bj_resistance_from_geometry` and `bj_parallel_plate_capacitance` refuse
  a NaN area or separation, which Rust's `<= 0` let through. 512 rows mapped. Reviewed and kept as Rust's: the zero
  guards on inputs (the user's material.rs answer), and `rlc_impedance`'s lighter validation, since L = 0 (a series RC) and
  R = 0 (an LC) are circuits it must accept, where `validate_lc` would refuse them.
- **New golden file `circuit.tsv`, 56,106 rows (13 MB)**: class C for the four transients, class A for the rest; pairs
  over 39 values, triples over 16, quadruples over 8 (the solenoid, the impedance, the transients with five source
  values), slices of length 0 to 6, the 1e-30 cuts of Q's R, the impedance's ω and C, the transients' τ and the
  parallel combinations' Σ 1/x approached ulp by ulp, and seeded physical circuits. glibc's exp is correctly rounded
  on every class-C row, so there is no `circuit.cr.tsv`. `validate.py`'s independent transcription agreed with every
  row on its first run.
- **Rust tests**: all 51 of circuit.rs's tests ported (`tests/circuit.tcyr`, 69 assertions with the Cyrius-only
  groups). The ledger: 313 ported, 3 n/a, 542 pending of 858.
- **Surface parity**: 8 modules ported, 220 items matched (27 new), 0 problems. The API snapshot gains 27 public fns.
- **Benchmarks** (the five circuit twins; 38 of 80 now twinned; load average 3.9–5.0 and rising, each the lower of
  two runs; the A/B against HEAD shows the earlier rows unchanged within that drift, `error_record` the most, 128–130 →
  141–142 ns with its code untouched):

  | Bench | Rust | Cyrius |
  |---|---|---|
  | `circuit/rc_charging` | 0.50 ns | 131 ns |
  | `circuit/rl_current_rise` | 2.03 ns | 125 ns |
  | `circuit/rlc_impedance` | 1.53 ns | 17 ns |
  | `circuit/resonant_frequency` | 1.31 ns | 30 ns |
  | `circuit/resistance_parallel` | 7.09 ns | 34 ns |

  Rust's arguments are constants that LLVM folds; the transients' cost is crmath's exp.
- **For M2's range-edge Part** (§0.1), from the audit:
  - accuracy, a question for the user with the half-wave pattern's (§26): `rc_charging_voltage` and `rl_current_rise`
    compute 1 − e^(−t/τ) by cancellation, a relative error of about 2^-53·τ/t (8e-8 at t/τ = 10^-10), and 0 for
    t/τ below 1.1·10^-16; −expm1(−t/τ), correctly rounded, has none;
  - the |τ| < 1e-30 tests on the rounded RC and L/R, and |Σ 1/xᵢ| < 1e-30 on a rounded sum;
  - chains that over- or underflow where the result is representable: √(LC) in the resonances (L = C = 10^-200 gives
    f₀ = ∞), ωL − 1/(ωC) and R² + X² in the impedance (a hypot-style scaling), 0.5·C·V², μn²Al, ρL/A, εA/d.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **12 modules, 0 failing; 376,664 ok, 0 failed, 0
  measured, 0 pending, 17,909 accepted (ADR 01: 8,999; ADR 0002: 8,910)** over 394,573 rows; **all 249 digests
  equal**, the 222 before unchanged.

## 28. M2 part 1d-1: relativity, the closed forms (2026-11-12)

relativity.rs is two bites: 1d-1, everything but the callbacks, and 1d-2, `solve_retarded_time` and
`lienard_wiechert_fields` over the plan's trajectories.

- **`src/relativity.cyr`**: the Lorentz factor and β; `BjEmTensor` (`from_fields`, `electric_field`, `magnetic_field`,
  both invariants, `dual`, `four_force`, `eq`) with Rust's `components: [[f64; 4]; 4]` as a 16-element tuple field,
  row-major, and `BjEmTensor_from_components` for Rust's struct literal; the two field transforms (two out pointers);
  `BjFourVector` (`new`, `from_time_position`, `dot`, `interval_sq`, `spatial_magnitude`, `boost_x`, `eq`);
  `BjFourPotential` (a nested `a: BjFieldVector`); the retarded and Liénard–Wiechert potentials, `lienard_wiechert_e`
  and `_b`, and both Larmor powers. Every formula is Rust's in Rust's order; the FieldVector operations are
  BjFieldVector's.
- **Probe p17** (`scripts/probes/p17_tensor_struct.cyr`): the 16-element tuple field, the 128-byte struct returned by
  value and copied, and the nested struct field all work. A method's untyped `self` is a `*T`, so the struct's address
  is `self` itself; `&self` is the slot that holds it, as the guide says of any pointer.
- **ADR 0002 row 27**, row 20's rule: `bj_lienard_wiechert_e` refuses a NaN R as the potentials do, where Rust went on
  to NaN. 228 rows mapped.
- **Row 18 through FieldVector**: the chains on magnitude, magnitude_sq, dot and cross take 5d-2's unbounded chain at a
  range edge. `validate.py` transcribes each function once on a set of vector operations: Rust's plain ones, which
  reproduce the oracle, and the port's (`fv_governed` / `fv_model` per operation), which give the port's outputs. 720
  rows mapped.
- **New golden file `relativity.tsv`, 26,052 rows, class A**: speeds through c ulp by ulp; 841 field pairs and their
  tensors plus arbitrary ones (specials in each slot, seeded); boosts at β ∈ {0, 0.1, 0.9, 0.99, 1 − 1e-9} along x, y
  and the diagonal, at and past c; four-vectors; Liénard–Wiechert geometries with the 1e-30 cuts of R and β → 1.
  `validate.py`'s transcription agreed with every row on its first run.
- **Rust tests**: 39 of relativity.rs's 48 ported (`tests/relativity.tcyr`, 64 assertions with the Cyrius-only
  groups); the other nine test the solver (1d-2). The ledger: 352 ported, 3 n/a, 503 pending of 858.
- **Surface parity**: relativity in progress, 255 items matched (35 new), 0 problems; the solver's two functions
  pending. The API snapshot gains 16 public fns.
- **Benchmarks** (six of relativity's eight twins; 44 of 80 now twinned; load average 4.3–10.8, each the lower of two
  runs; the first round ran at load 9–11 and doubled every row, the second agrees with HEAD):

  | Bench | Rust | Cyrius |
  |---|---|---|
  | `relativity/lorentz_factor` | 1.24 ns | 13 ns |
  | `relativity/em_tensor_from_fields` | 8.11 ns | 36 ns |
  | `relativity/lorentz_transform_x` | 5.91 ns | 42 ns |
  | `relativity/lienard_wiechert_e` | 21.64 ns | 227 ns |
  | `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | 64 ns |
  | `relativity/four_vector_boost` | 5.45 ns | 33 ns |

- **For M2's range-edge Part** (§0.1), from the audit:
  - accuracy, a question for the user with §26's and §27's: γ = 1/√(1 − v²/c²) loses digits near c (1 − β² is a
    cancellation, a relative error of about 2^-53/(1 − β²)), where 1/√((1 − β)(1 + β)) does not;
  - FourVector's dot, interval_sq and spatial_magnitude and `four_force` are raw chains (no FieldVector), so they over-
    and underflow where the result is representable, as FieldVector's did before row 18;
  - `dual` round-trips B·c/c and E·(−1/c) (an ulp, and ∞ past DBL_MAX/c); an ∞ R in `lienard_wiechert_e` makes
    n̂ = R·0 = NaN; the Larmor chains q²a².
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 401,768 ok, 0 failed, 0
  measured, 0 pending, 18,857 accepted (ADR 01: 8,999; ADR 0002: 9,858)** over 420,625 rows; **all 274 digests
  equal**, the 249 before unchanged.

## 29. M2 part 1d-2: relativity, the retarded-time solver (2026-11-12)

- **`bj_solve_retarded_time` and `bj_lienard_wiechert_fields`** in `src/relativity.cyr`: Rust's bracket-doubling and
  safeguarded Newton/bisection on the delay τ, its stopping tests (4ε·max(|t|, τ), and |g| at its own rounding level
  with one last Newton correction) and its errors, step for step. The trajectory is a closure called through `callptr`
  as `trajectory_fn(t, out)`: t an f64 in a word, position and velocity (and for the fields the acceleration) written
  through `out`, the return value ignored since Rust's trajectory can't fail. `max_iterations` is Rust's usize as a
  word, −1 being usize::MAX. The fields call one solver core with the user's nine-value closure, so no wrapping
  closure is allocated. The four `tracing::debug!` events and the fields' `trace!` are `_bj_log_text` lines in
  tracing's `message key=value` form.
- **Probe p18** (`scripts/probes/p18_trajectory_closure.cyr`): an f64 time passed through `callptr` and rebound inside
  the closure, six f64 written through `out`.
- **Scenarios** (the plan's closure inputs), the same in `m_relativity.rs`, `validate.py` and the harness: 0 a static
  charge, 1 uniform motion, 2 a circular orbit (sin, cos: class C), 3 hyperbolic motion, 4 an oscillation along z
  (class C), 5 a static charge with a given acceleration. **192 rows** in `relativity.tsv` (now 26,244): 102 solver
  rows (the Rust tests' cases, speeds to (1 − 1e-9)c both ways, superluminal ones, the β = 0.99 orbit's observers and
  seeded ones, budgets 0, 1, 2, 3 and usize::MAX, non-finite observers and states) and 90 fields rows. `validate.py`
  transcribes the solver on the primitives P; glibc's sin and cos are correctly rounded on every trajectory call, so
  there is no `relativity.cr.tsv`. Its transcription agreed with every row on its first run.
- **Rust tests**: relativity.rs is complete, all 48 ported (`tests/relativity.tcyr`, 79 assertions). The ledger: 361
  ported, 3 n/a, 494 pending of 858.
- **Surface parity**: relativity ported, 257 items matched, 0 problems. **M2's four modules are ported.**
- **Benchmarks** (the solver's two twins; 46 of 80 now twinned; load average 2.8–3.5, each the lower of two runs, the
  A/B against HEAD unchanged):

  | Bench | Rust | Cyrius |
  |---|---|---|
  | `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | 839 ns |
  | `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | 5.44 µs |

  The orbit's cost is its trajectory's sin and cos at every evaluation (crmath's fast path, about 200 ns each).
- **For M2's range-edge Part** (§0.1), from the audit:
  - a superluminal trajectory is reported as "trajectory returned a non-finite state" (its position overflows while τ
    doubles), not as the bracket error Rust's docs give for it ("could not bracket the light-cone root (is the
    trajectory faster than light?)"); with finite positions that error is unreachable, since cτ passes any distance
    first;
  - the light-cone function's |R| is a raw √(dx² + dy² + dz²): past DBL_MAX it is ∞ and g = ∞ − cτ turns NaN when cτ
    overflows, which ends the bracket loop without an error.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 401,960 ok, 0 failed, 0
  measured, 0 pending, 18,857 accepted (ADR 01: 8,999; ADR 0002: 9,858)** over 420,817 rows; **all 276 digests
  equal**, the 274 before unchanged.

## 30. M2 part 2a-1: the pattern's and the transients' cancellations (2026-11-12)

M2's range-edge Part (§0.1) runs as 2a–2d on the user's answers of 2026-11-12: 2a the three accuracy rows (2a-1 the
dipole pattern and the transients, 2a-2 a correctly rounded γ), then M1's three rules over M2's candidates (§25–§29):
2b round once, 2c exact thresholds, 2d an unbounded exponent at the edges.

- **ADR 0002 row 28**: `bj_half_wave_dipole_pattern`'s numerator is sin(π·sin²(θ/2)), or sin(π·cos²(θ/2)) where
  cos θ < 0: no cancellation near either axis. sin θ and cos θ come from one `_bj_sin_cos`, which keeps the bench at
  HEAD's speed (504 and 483 ns against 508 and 511; the four separate calls cost 610). 661 rows mapped; the two
  corrected-oracle rows of the old formula are among them, so `wave.cr.tsv` keeps 1 row.
- **ADR 0002 row 29**: the charging transients take −expm1(−t/τ). kimiya's correctly rounded expm1 joins
  `src/crmath_kx.cyr`: `scripts/port/kimiya_closure.py` gains the root, and its copy of kimiya's tests now stops at
  kimiya's B10a group (the unbounded exp family, which bijli doesn't carry: kimiya's tests had moved on, and the
  closure's check had gone stale). `scripts/oracle/crtables_k.py` renders kimiya's exp tables and limits into
  `src/crmath_k_data.cyr`, and `cr_verify.py` gains an expm1 family: 400,000 arguments through its thresholds, all
  equal to MPFR on x86_64 and aarch64. `_bj_expm1` in `src/mathx.cyr` stands for `f64::exp_m1`. 657 rows mapped.
- **Tests**: `cyrius_adr_row_28` checks the pattern at θ = 10^-8, 10^-12 and π − 10^-10 to 10^-15 of 80-digit
  references; `cyrius_adr_row_29` the transients at t/τ = 10^-20 and 10^-10.
- **Benchmarks** (an A/B against HEAD, load 3.4–8, each the lower of two runs): `circuit/rc_charging` 123 → 70 ns
  (−43 %), `circuit/rl_current_rise` 121 → 71 ns (−41 %): kimiya's expm1 is cheaper than crmath's exp; the other rows
  unchanged.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 400,644 ok, 0 failed, 0
  measured, 0 pending, 20,173 accepted (ADR 01: 8,997; ADR 0002: 11,176)**; **all 276 digests equal** (the pattern's
  and the two transients' are row 28's and 29's).

## 31. M2 part 2a-2: a correctly rounded γ (2026-11-12)

- **ADR 0002 row 30** (the user's answer of 2026-11-12): `src/crmath_gamma.cyr`'s `_bj_cr_lorentz` gives
  γ = c/√(c² − Σvᵢ²) or γ² = c²/(c² − Σvᵢ²), rounded once, from the exact difference D of the squares, and
  `_bj_cr_lorentz_ok` its validity test alone (D > 0 exactly). c lies in [2^-30, 2^30]: bijli passes the speed of
  light, or 1 for a β, and any other c is refused (the slow path's buffers and the fast paths' ranges assume it).
  - One nonzero component a: D = (c − a)(c + a) as a double-double (two TwoSums and one product), within 2^-103.
  - Several: D as Shewchuk's nonoverlapping expansion of c² and the −vᵢ² (exact), compressed to its top two
    components (within 2^-104). Components below 2^-450 are left out and flagged.
  - Both go through kimiya's rounding test (`_bj_k_dd_div_cr`) on c/√D (a double-double square root) or c²/D,
    unscaled: in the c range every quantity is inside the test's ranges with no term subnormal, so a scaling by 4^j
    would change no bit.
  - The slow path, when the test fails (a value within about 2^-100 of a midpoint) or a flagged component leaves D
    under 2^-790: D exact as a bignum in units of 2^K, the candidate the lower of the test's two ends
    (`_bj_k_dd_div_ends`), then the exact signs of value − μ at the neighbouring midpoints, ties to even.
- **Callers**: `bj_lorentz_factor` (and through it `bj_lorentz_transform_x` and `BjFourVector.boost_x`),
  `bj_lorentz_transform_fields` and `bj_relativistic_larmor_power` take γ with c the speed of light;
  `bj_lienard_wiechert_e` (and through it `bj_lienard_wiechert_fields`) takes γ² with c = 1; the Liénard–Wiechert
  potentials test their β with `_bj_cr_lorentz_ok`. A β whose exact Σβ² lies just below 1 is accepted
  ((0.2903717016735131, 0.9569139328420442, 0), Σβ² = 1 − 6.2·10^-17, which Rust refused); (0.6, 0.8), exactly just
  above 1, is refused.
- **Golden rows**: `validate.py`'s transcription takes the exact γ (`rel_exact_gamma`, on Fractions) where the port
  does. `relativity.adr.tsv` maps 1,308 row-30 rows; row 18's 720 became 717, its three rows that γ also changes
  being row 30's. The class-C rows that γ changes are row 30's too, never a correction and a mapping for one row.
- **Verification**: `cr_verify.py`'s gamma families check `_bj_cr_lorentz` (both paths, and the exact path alone) and
  `_bj_cr_lorentz_ok` against Fractions and an integer square root, over c = the speed of light, 1, 2^-30, 2^30 and
  four random c between, with c outside the range, ∞ and NaN refused: **1,500,020 arguments, all equal to the exact
  value on x86_64 and aarch64**. Widening the families to the whole c range found two defects before this Part
  landed: the first γ² family scaled its β by 1/c twice, so every case fell under the 2^-28·c shortcut and the γ²
  path went unchecked; and the slow path took its candidate from the rounding test's verdict, 0 exactly when that path
  is needed, so the values just above a midpoint came out near 2^-1020 (β = 1 − 2^-52 gave `lienard_wiechert_e` a
  γ² of 8.9·10^-308). The candidate is now the test's lower end; `cyrius_adr_row_30` keeps four such values.
- **Tests**: `cyrius_adr_row_30` (relativity.tcyr, 96 assertions): γ at β = 1 − 10^-9 and one ulp below c against
  exact references, |v| = c refused, a vector at 0.991c, the two potentials' β above, the c range's edges (2^±31 and
  NaN refused), and four values just above a midpoint, on the fast path and on the slow path alone.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds, load 3.4–3.9 for the first and 8–10 for the third;
  each the lowest of the three):

  | Bench | HEAD | 2a-2 |
  |---|---|---|
  | `relativity/lorentz_factor` | 12 ns | 153 ns |
  | `relativity/lorentz_transform_x` | 39 ns | 177 ns |
  | `relativity/four_vector_boost` | 30 ns | 168 ns |
  | `relativity/lienard_wiechert_e` | 208 ns | 343 ns |
  | `relativity/lienard_wiechert_scalar_potential` | 58 ns | 84 ns |
  | `relativity/lienard_wiechert_fields_orbit_b099` | 5.56 µs | 6.26 µs |

  `em_tensor_from_fields` and the solver are unchanged. The first draft cost 330 ns for the factor; scaling D to
  [1, 4) and back took 25 ns of it and is gone (exact in the c range, so no bit changed). **Finding**: the rest is the
  arithmetic itself. A dependent f64 operation costs about 1 ns in cc5's code (every operand through the stack and a
  general register) and a call about 3 ns; the one-component path is about 90 operations, two divisions and a
  square root: D in double-double 39 ns, the square root's correction 54 ns, kimiya's quotient and test 52 ns. A
  vector of two or three nonzero components builds and compresses an expansion first (about 115 ns more). A scheme
  with one division (the residuals of c/√D taken together) would need a new error proof, not taken for about 20 %.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 399,339 ok, 0 failed, 0
  measured, 0 pending, 21,478 accepted (ADR 01: 8,997; ADR 0002: 12,481)** over 420,817 rows; **all 276 digests
  equal**. Seven changed: `lorentz_factor`, `lorentz_transform_x`, `lorentz_transform_fields`, `FourVector::boost_x`,
  `relativistic_larmor_power`, `lienard_wiechert_e` and `lienard_wiechert_fields`, row 30's.

## 32. M2 part 2b-1: material rounded once (2026-11-12)

2b applies M1's round-once rule (ADR 0002 row 13) to M2's candidates (§25–§29), one module at a time; material first.

- **ADR 0002 row 13, amended**: material's single products, quotients and roots of inputs and constants are correctly
  rounded for finite operands, on kimiya's monomials (`src/crmath_kx.cyr`) and bijli's √ of one
  (`src/crmath_mono.cyr`):
  - `bj_dielectric_energy_density`, ε_r·ε₀·E²/2 (`mul5`, with 0.5 one of its factors: one bit of the exact path's
    253);
  - `bj_magnetic_energy_density_material`, B²/(2μ_rμ₀) (`mul2_div2` on 2μ₀, an exact doubling), after Rust's test of
    the rounded |μ_r·μ₀| < 1e-30 (a threshold, 2c's);
  - `bj_displacement_field`, `bj_b_field_from_h` and `bj_h_field_from_b`, part by part: ε_r·ε₀·Eᵢ and μ_r·μ₀·Hᵢ
    (`mul3`), Bᵢ/(μ_r·μ₀) (`div2`). A non-finite ε_r or μ_r keeps Rust's whole chain, and a non-finite part its own
    product by Rust's factor;
  - `BjMaterial`'s `refractive_index` (√(ε_rμ_r)), `impedance` (√(μ_rμ₀/(ε_rε₀))), `wave_speed`
    (1/√(ε_rε₀μ_rμ₀)) and the lossless optical constants n = √(ε_rμ_r) and η = √(μ_r/ε_r) (kimiya's √(x/y)), which
    the Material-based Fresnel and reflectance functions use. `refractive_index` takes any operands: a zero or
    non-finite one keeps Rust's √(ε_r·μ_r), both negative give √(|ε_r||μ_r|), and opposite signs give NaN, the exact
    product being negative (Rust's −0 where it underflows).
  - Out of row 13's scope, as M1's precedent (§17: a sum inside gets the edge treatment): `polarization`
    ((ε_r − 1)·ε₀·E), `magnetization` ((μ_r − 1)·H), Clausius–Mossotti and Curie–Weiss.
- **Mappings** (`validate.py`'s `material_port13`, from exact values): `material.adr.tsv` maps 3,414 row-13 rows
  (dielectric energy 903, magnetic energy 690, wave speed 525, impedance 441, H 209, the index 189, B 168, D 166, the
  reflectance 123) and `wave.adr.tsv` 509: `reflectance_at_interface` 493, `reflectance_normal_materials` 14 and
  `material_reflection_coefficient`'s 2 that only the rounding moves; its 74 row-25 rows keep row 25 (the relative η
  alone moves them), now with η rounded once. 3,931 rows failed before the mappings existed, and exactly those are
  mapped. The harness passed every one on its first run.
- **Tests**: `cyrius_round_once (ADR 0002 row 13)` in `tests/material.tcyr` (11 assertions) against exact
  references: the energy densities at 1e-300 and 1e200 (Rust: …397905e-12, ∞), D and B with a subnormal ε_r·ε₀ or
  μ_r (a NaN and an ∞ part kept), H one ulp off in Rust, n = 1e200 and DBL_MAX (Rust: ∞), the NaN of opposite signs,
  v = 2.997924580000065e208 at 1e-200 (Rust: ∞), η = 3.77e-298 (Rust: 0); and in `tests/wave.tcyr`
  `reflectance_at_interface` into a DBL_MAX/DBL_MAX material, R = 0 where Rust's n₂ = ∞ was an error.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds at load 1.6–2.0, the lowest of each):
  `material/displacement_field` 23 → 231 ns and `material/h_field_from_b` 22 → 256 ns (three correctly rounded parts
  each, about 70 ns a monomial), `wave/reflectance_at_interface` 115 → 478 ns (two materials' n and η); the rest
  unchanged. **Finding**: the three parts share x·c, so one exact TwoProduct and a rounding test per part could serve
  them, but that is a new proof, not taken here.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 395,416 ok, 0 failed, 0
  measured, 0 pending, 25,401 accepted (ADR 01: 8,997; ADR 0002: 16,404)** over 420,817 rows; **all 276 digests
  equal**. Twelve changed: the nine above and wave's three Material-based functions.

## 33. M2 part 2b-2: wave rounded once (2026-11-12)

- **ADR 0002 row 13, amended**: wave's products and quotients of inputs and constants are correctly rounded for
  finite operands (kimiya's monomials, `src/crmath_kx.cyr`); a constant enters as its double, 2π and 4π as Rust's π
  times 2 and 4 (exact), and a non-finite operand keeps Rust's chain:
  - `bj_plane_wave_intensity`, c·ε₀·E₀²/2 (`mul5`, 0.5 a factor), Rust's four roundings an ulp off at 1 kV/m;
  - `bj_radiation_pressure_reflected`, 2I/c: Rust's 2I is exact up to DBL_MAX/2, where its result already rounds
    once; past it 2I overflowed, and RN(I/c)·2 is RN(2I/c) exactly (I/c normal, the doubling exact), so no monomial
    is needed;
  - `bj_wave_number`, 2πf/v (`mul_div`), after Rust's |v| < 1e-30 test;
  - `bj_cylindrical_waveguide_cutoff`, x·v/(2πa) (`mul2_div2`): Rust's x·v and 2πa over- or underflowed;
  - `bj_hertzian_dipole_radiation_resistance`, K·(dl/λ)² on row 22's double K = RN((2π/3)·η₀) (`mul3_div2`);
  - `bj_effective_aperture`, λ²D/(4π) (`mul3_div`), after row 23's test;
  - `bj_friis_transmission`, G_tG_rλ²/((4π)²d²) (kimiya's general monomial, four over four).
  - Already rounded once: `radiation_pressure_absorbed`, `e_from_b`, `b_from_e`, `angular_frequency`. Out of row 13's
    scope: the Poynting vector and momentum density (a cross product inside; their range edge, E × B past DBL_MAX
    where S/c² is representable, is 2d's), the Fresnel, Snell and waveguide formulas (sums), and the patterns.
- **Mappings** (`validate.py`'s `wave_port13`): `wave.adr.tsv` maps 459 more row-13 rows (friis 175, aperture 129,
  wave number 99, the circular cutoff 40, intensity 15, the reflected pressure 1); the Hertzian resistance's 625 rows
  keep row 22 (bijli's K alone moves them), their outputs now the product rounded once. Every one passed on the
  harness's first run.
- **Tests**: `cyrius_round_once (ADR 0002 row 13)` in `tests/wave.tcyr`, 7 more assertions against exact references:
  the intensity at 1 kV/m (Rust: …4783), 2I/c at DBL_MAX, 2πf/v at 1e308, x·v = 1e400, λ² = 1e400 and G_tG_r = 1e400
  (Rust: ∞ each), and K·(dl/λ)² at 0.01 (the chain: …275).
- **Benchmarks** (an interleaved A/B against HEAD, three rounds at load 4–10, the lowest of each):
  `wave/plane_wave_intensity` 7 → 148 ns (kimiya's five-factor product); the other wave rows within the load's noise
  (±8 %). The other six functions have no Rust bench.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 394,957 ok, 0 failed, 0
  measured, 0 pending, 25,860 accepted (ADR 01: 8,997; ADR 0002: 16,863)** over 420,817 rows; **all 276 digests
  equal**. Seven changed: the seven above.

## 34. M2 part 2b-3: circuit rounded once (2026-11-12)

- **ADR 0002 row 13, amended**: circuit's products, quotients and roots of inputs and constants are correctly rounded
  for finite operands; 2π is Rust's π doubled, and a non-finite operand keeps Rust's chain:
  - `bj_resistance_from_geometry` (ρL/A) and `bj_parallel_plate_capacitance` (εA/d), `mul_div`, after row 26's tests;
  - `bj_capacitor_energy` and `bj_inductor_energy` (CV²/2, LI²/2), `mul4`;
  - `bj_solenoid_inductance` (μn²Al), five full-precision doubles: 265 significant bits, past the 253 of kimiya's
    exact path (four doubles and a 41-bit integer per side, NaN beyond), which a subnormal result or a value within
    2^-100 of a midpoint reaches. **`_bj_cr_mul5`** (`src/crmath_mono.cyr`) takes kimiya's fast path
    (`_bj_k_mul5_fast`, whose error budget doesn't depend on the factors' sizes), then **`_bj_cr_prod_acc`**: the
    integer significands' product in bijli's limbs, its top 60 bits and a sticky bit, rounded once by
    `_bj_k_cr_from_int`. `cr_verify.py`'s new `prod5` and `prod5_acc` families (the full function and the exact path
    alone: any binade, near one, products steered to the subnormal boundary, exact ties, signed zeros): **600,000
    calls, all equal to the exact value on x86_64 and aarch64**. (The `mul5` calls of 2b-1 and 2b-2 carry 0.5, one bit,
    so they stay inside kimiya's domain: 213 and 188 bits.)
  - `bj_resonant_frequency` (√(1/((2π)²LC))), `bj_resonant_angular_frequency` (1/√(LC)) and `bj_quality_factor`
    (√(L/(C·R²))), bijli's √ of a monomial;
  - `bj_damping_ratio` (R·√(C/L)/2), kimiya's RN(w·√((x/y)·2^e)) for R > 0; a zero R keeps Rust's chain (±0).
  - Already rounded once: Ohm's law, the power, the charge and the time constants. Out of row 13's scope (sums
    inside): the combinations, the transients and `rlc_impedance`.
- **Mappings** (`validate.py`'s `circuit_port13`): `circuit.adr.tsv` maps 1,738 row-13 rows (resonant frequency 398,
  angular 292, damping 249, Q 192, the energies 135 each, the solenoid 133, ρL/A and εA/d 102 each), every one
  passing on the harness's first run. `rel_exact_gamma` and `cr_verify.py`'s γ reference no longer copy the port's
  |vᵢ| < 2^-28·c shortcut: they compute those values exactly too (the mappings are unchanged: the shortcut's bound
  holds).
- **Tests**: `cyrius_round_once (ADR 0002 row 13)` in `tests/circuit.tcyr` (11 assertions): ρL = 1e310 and
  εA = 1e-400 (Rust: ∞, 0), the energies an ulp off in Rust at 4.7 nF and 230 V, the solenoid an ulp off in Rust and
  on the subnormal grid (the exact path), the resonances at L = C = 1e-200 and 1e200 (Rust: ∞, 0), Q with L/C = 1e600
  and ζ with C/L = 1e-600 (Rust: ∞, 0), and a −0 R keeping −0.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds at load 1.9–3.8, the lowest of each):
  `circuit/resonant_frequency` 29 → 217 ns (the √ of a four-factor monomial); the other circuit rows unchanged.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 393,219 ok, 0 failed, 0
  measured, 0 pending, 27,598 accepted (ADR 01: 8,997; ADR 0002: 18,601)** over 420,817 rows; **all 276 digests
  equal**. Nine changed: the nine above.

## 35. M2 part 2b-4: relativity rounded once (2026-11-12)

- **ADR 0002 row 13, amended**: relativity's products and quotients of inputs and constants are correctly rounded for
  finite operands; a non-finite operand keeps Rust's chain:
  - `bj_retarded_scalar_potential`, q/(4πε₀R) (`div3`), and `bj_retarded_vector_potential`, μ₀qvᵢ/(4πR) part by part
    (`mul3_div2`), after Rust's radial-distance checks;
  - `bj_larmor_power`, q²a²/(6πε₀c³) with c³ of c's double (c = 299792458 is exact; Rust's C3 = c·RN(c²) rounded
    twice) and 6π on Rust's π. Six factors below, and none merge exactly to fewer than five (6c is a double; c², 6c²
    and πε₀ are not), past kimiya's four: **`_bj_cr_mono_wide`** (`src/crmath_mono.cyr`) takes kimiya's steps with
    six factors a side. Its fast path's two double-doubles carry at most 12.04u² each (3.01u² per significand past the
    second), inside the rounding test's 32u²; its exact path is kimiya's (the integer significands' quotient rounded
    by `_bj_k_cr_ratio`), with kimiya's 253-bit capacity a side (the Larmor power needs 212 over 192).
    `cr_verify.py`'s `mono_wide` and `mono_wide_acc` families (up to six a side within the capacity, the Larmor shape
    at SI scales and at the range ends, ties): **511,324 calls, all equal to the exact value on x86_64 and aarch64**.
    They need fourteen words a call, so `cr_driver.cyr` now reads sixteen (it read twelve and silently dropped the
    rest, which the first run caught as 114 mismatches), and `cr_verify.py` refuses a call longer than that. Every
    exact family passes on the new driver (574,084 calls).
  - `BjEmTensor.dual`: its components are (cB)/c = B and −(cF^{i0})/c = −F^{i0}, F's own exactly. With all six
    sources finite the dual is built from them with no arithmetic (60 ns, Rust's chain 122); a non-finite source
    keeps Rust's chain for its components.
  - Out of row 13's scope (sums inside): the invariants, the four-vector dot, interval and magnitude, `four_force`,
    the transforms, the Liénard–Wiechert potentials and fields, `lienard_wiechert_b` and the relativistic Larmor power.
- **Mappings** (`validate.py`'s `REL13_FNS`, `rel_dual13`): `relativity.adr.tsv` maps 320 row-13 rows (the dual 273,
  the vector potential 25, Larmor 20, the scalar potential 2), every one passing on the harness's first run.
- **Tests**: `cyrius_round_once (ADR 0002 row 13)` in `tests/relativity.tcyr` (6 assertions): the dual of B = 1e305
  (Rust: ∞) and of a B whose (cB)/c Rust rounds an ulp off, a NaN part kept; q/(4πε₀R) an ulp off in Rust;
  μ₀qv/(4πR) on the subnormal grid with an ∞ part kept; Larmor with q² = 1e320 (Rust: ∞) and at the electron's
  scale (Rust: an ulp off).
- **Benchmarks**: no Rust bench covers the four, and the relativity rows of an interleaved A/B (load 4–5.6) are
  unchanged. A probe, both trees in a row: `larmor_power` 9 → 359 ns (the six-factor monomial),
  `retarded_vector_potential` 41 → 478 ns (three monomials), `dual` 122 → 60 ns.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 392,899 ok, 0 failed, 0
  measured, 0 pending, 27,918 accepted (ADR 01: 8,997; ADR 0002: 18,921)** over 420,817 rows; **all 276 digests
  equal**. Four changed: the four above. **2b is done**: row 13 now covers M2's four modules.

## 36. M2 part 2c-1: exact thresholds (2026-11-12)

2c applies M1's exact-threshold rule (ADR 0002 row 17) to M2. The inventory (every comparison in the four modules) and
the user's answers of 2026-11-12 split it five ways: 2c-1 the tests on one rounded product, quotient, sum or |R| and the
literal readings; 2c-2 to 2c-5 the four tests at their value's own pole or zero, where an exact test needs an accurate
value beside it (Clausius–Mossotti, the trig-free TIR boundary, Liénard–Wiechert's κ, the parallel combinations). A
test on a sin or cos (Snell's |n₁ sin θᵢ/n₂| > 1, the angle-form Fresnel denominators, the half-wave pattern's
|sin θ| < 1e-15) stays on bijli's correctly rounded value (the user's answer): its exact decision has no bounded
precision. The comparisons of an input with a literal are exact already (no double lies between 10^-30 and the double
1e-30), and so are critical_angle's n₂/n₁ ≥ 1, guide_wavelength's 1 − (λ/λ_c)² ≤ 0 (λ ≥ λ_c on doubles: RN(r²) < 1
whenever r < 1) and Curie–Weiss's T − T_c ≤ 0.

- **ADR 0002 row 17, amended**:
  - the transients' |R·C| < 10^-30 (rc_charging and rc_discharging) and |L/R| < 10^-30 (rl_current_rise and decay),
    `magnetic_energy_density_material`'s |μ_r·μ₀| < 10^-30 and `reflectance_normal`'s |n₁ + n₂| < 10^-30 (and so
    `transmittance_normal`'s) on the exact value. The rounded value is within half an ulp, so it decides outside
    `_bj_x_band`'s 64 ulps of the double 1e-30 (an inline pre-filter, 5e-31 and 2e-30, decides most inputs with no
    call), and inside it `src/crmath_cmp.cyr`'s new `_bj_x_sign` does: the sign of a sum of ±c·m₁·m₂·2^k·10^q terms in
    integers, as `_bj_x_r2_below` does for r², with `_bj_x_prod_below`, `_bj_x_quot_below` and `_bj_x_sum_below` (a
    TwoSum first) on it. They hold for any finite operands (the widest alignment, a subnormal product against 10^60,
    stays under about 2,250 of the 2,560 bits);
  - the Liénard–Wiechert potentials' and E's |R| < 10^-30 (Σ Rᵢ² < 10^-60, M1's `_bj_x_r2_below`), and the fields' on
    the observer and the trajectory's position, the inputs (their difference rounds);
  - classify_magnetic's χ < −10^-15 is χ ≤ −1e-15 on doubles: the double 1e-15 lies 0.476 ulp above 10^-15 and no
    double falls between, so χ = −1e-15 itself is diamagnetic (Rust: paramagnetic). The Fresnel family's cos θᵢ band
    [−10^-6, 1 + 10^-6] already agrees with its doubles (the double 1e-6 lies below 10^-6, so c ≥ −1e-6 is c ≥ −10^-6;
    RN(1 + 1e-6) is the double nearest 1.000001, below it, and the next lies above), contrary to the estimate in the
    question put to the user, which said the band's edges would move by an ulp.
- **Knife edges** (`scripts/oracle/src/grids.rs`): the double 1e-30 lies 0.476 ulp above 10^-30, so a rounded product
  or quotient decides differently only when the exact value lies in [1e-30 − ulp/2, 10^-30), a band 0.024 ulp wide,
  which the corpus's ulp-stepped rows never reached. `knife_edge2` (pairs within an ulp of ±10^-30: R·C, L/R, μ_r·μ₀,
  n₁ + n₂, sums exact, rounded or cancelling) and `knife_band2` (a factor swept ±1,500 ulps, kept within
  [−0.05, 0.02] ulp: the band and its two sides) add 9,996 rows: circuit 8,480 (2,120 a transient), wave 720 (360
  a function), relativity 783 (261 a function, M1's `knife_edge(true)` points as R), material 13. μ_r·μ₀'s fixed μ₀
  leaves only a few lattice points near 10^-30, none in the band, so its rows all agree.
- **Mappings** (`validate.py`, Fraction): 1,420 row-17 rows: circuit 1,292 (`circuit_flip17`: rc 320 each, rl 326
  each; a flipped transient is row 17's, not row 29's or the corrected oracle's), wave 84, relativity 42 (classified
  before rows 13, 18 and 30), material 2 (classify at −1e-15). All passed on the harness's first run.
- **`cr_verify.py`'s below_prod, below_quot and below_sum** (exact families): pairs within a few ulps of 10^-p (p = 30
  mostly, every p to 60), sums that round or cancel, far ones: **895,656 calls, all equal to Fraction's verdict on
  x86_64 and aarch64**, and r2_below's 300,000 still.
- **Tests**: `cyrius_exact_thresholds (ADR 0002 row 17)` groups in tests/circuit.tcyr (R·C and L/R below 10^-30
  exactly, rounded to 1e-30: instant), tests/wave.tcyr (n₁ + n₂ refused) and tests/relativity.tcyr (|R| at least
  10^-30 exactly, computed where Rust gave Singularity, the potential and E); tests/material.tcyr's classify of
  −1e-15 (diamagnetic) and of the double next to it toward zero (paramagnetic).
- **Benchmarks** (an interleaved A/B against HEAD, load 2.2–2.9): the first draft called the band test on every input,
  `circuit/rc_charging` 67 → 88 ns, `circuit/rl_current_rise` 66 → 85 ns, the LW potential 79 → 92 ns; with the
  pre-filter, 69, 69 and 82 ns (HEAD 67–69, 66–69, 80–90 in the same runs): within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 399,374 ok, 0 failed, 0
  measured, 0 pending, 31,439 accepted (ADR 01: 8,998; ADR 0002: 22,441)** over 430,813 rows; **all 276 digests
  equal**. Eleven changed, the eleven functions with knife-edge rows.

## 37. M2 part 2c-2: Clausius–Mossotti, exact tests and an accurate ε_r (2026-11-12)

The first of 2c's four boundary values (the user's answer of 2026-11-12: exact tests, accurate values).

- **ADR 0002 row 17, amended**: `clausius_mossotti` decides X ≥ 1 and ε_r ≤ 0 on the exact P = Nα (3ε₀ − P ≤ 0, and
  below it 3ε₀ + 2P ≤ 0, ε₀ its double), returns ε_r = (3ε₀ + 2P)/(3ε₀ − P) rounded once, and reports the
  catastrophe's X = Nα/(3ε₀) rounded once (kimiya's `_bj_k_cr_mul2_div2`). Rust's χ rounds twice, and 1 + 2χ and
  1 − χ cancel at the edges: ε_r is many ulps off near the pole (578 at an exact tie below it) and 0 or negative near
  the zero. ε₀'s significand is 4·1,370,119,201,102,221 and three times the odd part fits 53 bits, so 3ε₀ and 1.5ε₀ are
  doubles and Rust's rounded tests are monotone against them: they only ever refuse inputs the exact tests accept (X
  within half an ulp below 1, where χ rounds to 1; ε_r up to about 2^-53, where 1 + 2χ rounds to 0), never the
  reverse. A negative Nα past DBL_MAX made χ = −∞ and ε_r = NaN where the value is −2.
- **`src/crmath_edge.cyr`** (new: the boundary values, 2c-3 to 2c-5's to join it): the rough χ sends |χ| < 2^-60 to
  ε_r = 1 and |χ| ≥ 2^136 to X ≥ 1 or ε_r = −2 (each inside its rounding interval); between, 3ε₀ = t_h + t_l and
  P = b_h + b_l exactly (TwoProduct, α and N scaled to [1, 2) first), every part normal. For χ in [−0.24, 0.49]
  (X in [−¼, ½]) num and den are double-doubles within 20.1u² together; elsewhere they are Shewchuk expansions of their
  four parts, whose top components' signs decide both tests exactly, compressed to within 2^-104. kimiya's rounding
  test rounds |num|/den; when it fails (within about 2^-99 of a midpoint), the candidate's two neighbouring midpoints
  are checked exactly: |num| − y·den − h·den is sixteen doubles (TwoProducts, and the power-of-two h exactly), summed by
  kimiya's `_bj_k_exp_sign`, ties to even. A nonzero ε_r lies in [2^-109, 2^110], so nothing under- or overflows.
  Exact ties exist, and only here: a midpoint needs the odd part of its denominator to divide 9e (e = ε₀'s odd
  significand = 3·71·6,432,484,512,217), which leaves ε_r = 9e·2^-j − 2 (j = 0 … 50, X just below 1) and
  ε_r = −2 + 3f·2^-53 (f dividing e, X far below −½).
- **Rows** (`m_material.rs`, on a seed of its own so the other rows stay put): 2,046 more: both edges ±6 ulps for 48 N
  of every binade (1,248), the exact ties at three scales (177), the exact pole and zero (9), Nα past DBL_MAX, on the
  subnormal grid and at the extremes (12), X uniform in (−4, 1) over N of every binade (600); 7,049 in all.
- **Mappings** (`validate.py`'s `material_cm_port`, Fraction): 1,809 row-17 rows: 1,117 values an ulp or more off
  (X = 1/3 exactly: 2.5, Rust 2.4999999999999996), 481 InvalidPermittivity values (X = −2/3: −0.2, Rust
  −0.19999999999999998), 146 catastrophe texts (X rounded once), 29 NaN payloads now −2, and 36 refusals now accepted
  (18 at the pole, ε_r up to 6.3e16; 18 at the zero, ε_r down to 3.2e-17). All passed on the harness's first run.
- **`cr_verify.py`'s cm and cm_slow** (exact families; cm_slow skips the rounding test, so the exact check always
  runs): physical pairs, X uniform, both edges ±6 ulps, any finite pair, the exact ties (41,301 of the 300,000
  cases a midpoint), the region and shortcut edges, the exact pole and zero: **600,000 calls, all equal to Fraction's value on
  x86_64 and aarch64**.
- **Tests**: tests/material.tcyr's `cyrius_exact_boundaries (ADR 0002 row 17)`, 9 assertions: the 2.5 and −0.2, the
  two refusals now accepted, the exact pole's text, the exact zero, Nα past −DBL_MAX, and a tie on each side.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 1.0–2.3): `material/clausius_mossotti`
  (α = 1e-40, N = 1e28: the fast region) 26 → 152 ns as first written. Forming 3ε₀'s parts by FastTwoSum(2ε₀, ε₀)
  instead of a TwoProduct call, and Nα's TwoProduct inline on α and N themselves when both exponents lie in
  [−500, 500] (the scaled path stays for the rest; the same parts, so no result moves: every digest held), brings it
  to 128–134 ns (HEAD 26–27 in the same runs). The rest is the double-double num and den and kimiya's rounding test,
  about 80 dependent f64 operations on this compiler: a finding, recorded (Rust 1.24 ns), not a reason to round twice.
  The other material rows held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 399,611 ok, 0 failed, 0
  measured, 0 pending, 33,248 accepted (ADR 01: 8,998; ADR 0002: 24,250)** over 432,859 rows; **all 276 digests
  equal**. One changed: `material::clausius_mossotti`'s.

## 38. M2 part 2c-3: the trig-free Fresnel TIR, exact, with an accurate cos θ_t (2026-11-12)

The second of 2c's boundary values (the user's answer of 2026-11-12: exact tests, accurate values).

- **ADR 0002 row 17, amended**: the trig-free family's cosines (wave.rs's `fresnel_cosines`, under
  `snell_cos_theta_t`, the four `*_direct` functions, `reflectance_unpolarized`, `transmittance_unpolarized`,
  `schlick_reflectance` and `reflectance_at_interface`) decide total internal reflection as n₁²(1 − c²) > n₂² exactly
  and return cos θ_t = √(n₂² − n₁²(1 − c²))/n₂ rounded once; each function's own formula on top stays Rust's. Rust's
  1 − (n₁/n₂)²(1 − c²) rounds four times and cancels near the critical angle: just past it cos θ_t came out 0 (accepted)
  and just short of it TIR (cos θ_t = 4.3e-7 refused); n₁ = n₂ near grazing incidence gave 0 for cos θ_t = c
  (cos θ_i = 1e-200, even 0.3 came out 0.29999999999999993); and normal incidence with n₁/n₂ past DBL_MAX gave NaN
  (∞·0).
- **`src/crmath_edge.cyr`**: n₁ and n₂ scaled by n₂'s exponent (exact); shortcuts for n₁ = n₂ (|c|), c = 1 (1), a
  ratio past 2^27 (TIR, since 1 − c² ≥ 2^-52) and below 2^-61 (1). Where sin²θ_t ≤ ½, D = b² − a²(1 − c²) is a
  double-double within 18.6u² (TwoProducts, FastTwoSum, kimiya's `_bj_k_dd_mul`), √D within 15.4u², then kimiya's
  rounding test on √D/b. Near the critical angle D is the Shewchuk expansion of its ten parts, whose sign decides TIR;
  compressed and rooted within 8.06u². When the test fails, the exact check compares D with (μb)² at the candidate's two
  neighbouring midpoints: twenty-two doubles summed exactly. A c below 2^-300 is taken as 0, which provably moves
  nothing: no dyadic μ makes D − μ²b² vanish at c = 0 (it would need a Pythagorean triple with an odd hypotenuse
  dividing a leg), so D − μ²b² stays a multiple of 2^-442 while (ac)² < 2^-542. A nonzero D lies in [2^-620, 4] and
  cos θ_t ≥ 2^-311, inside every helper's range. `_bj_g_compress`'s scratch grew from 8 components to 16 (D has ten).
- **Rows** (`m_wave.rs`, on a seed of their own): 3,806 more: 463 for each of the eight trig-free functions (the
  critical cosine ±8 ulps for 24 pairs, ratios of 1 to 3, within 2^-50 of 1 and up to 2^26, a third scaled by
  2^±600; n₁ within 3 ulps of n₂ at seven grazing cosines down to 2^-700; ratios past DBL_MAX at three cosines; the
  shortcuts' exponent differences) and 102 for `reflectance_at_interface` (six dielectric pairs' critical cosines
  ±8 ulps).
- **Mappings** (`validate.py`'s `wave_cosines(…, exact=True)`, Fraction and an integer square root): 3,420 more row-17
  rows (`snell_cos_theta_t` 574, among them 50 values now TIR and 1 TIR now a value; `fresnel_rs_direct` 454,
  `fresnel_rp_direct` 467, `reflectance_s_direct` 405, `reflectance_p_direct` 418, `reflectance_unpolarized` 417,
  `transmittance_unpolarized` 330, `schlick_reflectance` 260, `reflectance_at_interface` 95), 68 more row-24 rows
  (Schlick's R₀ on the indices scaled by 2^-600, whose sum Rust's guard refused), and 4 interface rows relabelled from
  13 to 17; one interface row the port now matches Rust on (the exact cosine undoes row 13's ulp). The labels: Schlick
  is row 24 where R₀'s guard alone moves the row, else 17; the interface is 17 where the exact cosines alone move it,
  else 13. All passed once the n₁ = n₂ shortcut returned |c| (cos θ_i = −0 gives +0, √(c²)).
- **`cr_verify.py`'s cos_t and cos_t_slow** (exact families; cos_t_slow skips the rounding tests): uniform indices and
  angles, the critical cosine ±8 ulps at scales to 2^±990, n₁ within ulps of n₂ at small and tiny c, any positive
  finite indices, the shortcuts' edges, the fast region's edge a²(1 − c²) = b²/2, c = 0, c = 1 and n₁ = n₂.
- **Tests**: tests/wave.tcyr's `cyrius_exact_thresholds (ADR 0002 row 17)` gains 7 assertions: cos θ_t = c for n₁ = n₂
  (0.3 and 2^-1074), the two flipped rows, R = 1 past the critical angle, normal incidence past DBL_MAX, and R on the
  rounded-once cos θ_t.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 3.2–5.3: another repository's bench ran
  beside it, so both columns are inflated alike): `wave/fresnel_rs_direct` 61 → 258 ns, `wave/reflectance_unpolarized`
  72 → 276 ns, `wave/schlick_reflectance` 47 → 244 ns, `wave/reflectance_at_interface` 507 → 714 ns (HEAD's figures
  from the same runs). Each bench is in the fast region (sin²θ_t ≤ ½), so the cost is the double-double D, its root and
  kimiya's rounding test: about 130 dependent f64 operations and eight calls on this compiler, some 200 ns. A finding,
  recorded (Rust 1.24–7.78 ns), as Clausius–Mossotti's was; the other wave rows held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 399,933 ok, 0 failed, 0
  measured, 0 pending, 36,732 accepted (ADR 01: 8,998; ADR 0002: 27,734)** over 436,665 rows; **all 276 digests
  equal**. Nine changed: the nine functions on the cosines.

## 39. M2 part 2c-4: Liénard–Wiechert's κ and R_eff, exact tests and accurate values (2026-11-13)

The third of 2c's boundary values (the user's answer of 2026-11-12: exact tests, accurate values).

- **ADR 0002 row 17, amended**: the effective distance R_eff = |R| − R⋅β and κ = 1 − n̂⋅β = R_eff/|R| are each
  rounded once from the exact S = ΣRᵢ² and P = R⋅β of a finite R, and `lienard_wiechert_e`'s DivisionByZero is
  κ < 10^-30 exactly (the decimal, the user's literal reading). Both potentials take R_eff, E takes κ and the fields
  take E's; each function's formula around them stays Rust's, and an ∞ part keeps Rust's chain. Rust's two differences
  cancel as n̂⋅β → 1, a relative error of about 2^-53/κ:
  - on the axis at β = 1 − 2^-53, R_eff = 1000·2^-53 rounded to an ulp of 1000, so φ came out 2.3 % low;
  - just off the axis of an ultra-relativistic charge (R = (1, 2^-27, 0), β = (1 − 2^-53, 2^-26.5, 0)) κ = 6.0e-17
    rounded to 0, a DivisionByZero, and nearby R_eff to 0 (φ = ∞);
  - near the cone, E was off by up to a factor of 10.8 over the new rows;
  - past DBL_MAX's root, Rust's |R| overflowed and φ came out 0 or NaN (the port's |R| was row 18's already; R_eff
    now comes from the scaled S and P).
- **`src/crmath_edge.cyr`**: D = 1 − |β|² ≥ 2^-266 for any doubles with |β| < 1 (a lattice argument in the section
  comment), so κ lies in (2^-267, 2) and R_eff above 2^-367. R is scaled by its largest component's exponent; each part
  is split once and every product is a TwoProduct on the halves. Far from the cone's edge (P ≤ |R|/2) R_eff = √S − P
  in double-double (AccurateDWPlusDW, Joldes–Muller–Popescu 2017), within 51.3u²; near it Lagrange's identity gives
  S − P² = S(1 − |β|²) + |R × β|², a sum of non-negative terms built from D's exact expansion and the exact cross
  components, so R_eff = (S − P²)/(√S + P) is within 46.2u² and κ within 66.3u² however small. The rounding test is
  kimiya's with ε = 2^-98 (operands to 239u²). When it fails, the candidate's neighbouring midpoints decide exactly on
  bignums (S against (P + μ)², or (1 − μ)²S against P² by the signs), with src/fmtf.cyr's bignum and a product here, in
  buffers of 216 words (values below 2^6,560). κ's cut is decided by the test's ends where they clear 10^-30 by 2^-40,
  else as (10^30 − 1)²S < 10^60·P². Exact ties exist (κ = ½ + 2^-54 along (1, 2, 2)) and go to even.
- **Rows** (`m_relativity.rs`'s `emit_kappa_edges`, on a seed of its own, after every older row): 2,564 more:
  400 geometries along β within 2^-60…1 rad as |β| → 1 at scales 2^±90, 100 ultra-relativistic charges seen just off
  their axis, 120 exact midpoints of κ and R_eff along five Pythagorean directions, 100 at the far/near edge
  n̂⋅β = ½(1 ± 2^-k), and 48 extremes (R past DBL_MAX's root, near 10^-30, exponent spreads to 2^-1074; β zero,
  subnormal, at 1 − 2^-53), each through both potentials and E; and κ's cut, four sweeps of β₃ over ±32 ulps across
  10^-30 (D built in integers), through E. `validate.py`'s transcription agreed with every row on its first run.
- **Mappings** (`validate.py`'s `rn_lw` and `lw_below`, scripts/oracle/exact.py: an integer square root to 2^-700,
  then the exact signs at the neighbouring midpoints): `relativity.adr.tsv` grows from 2,615 rows to 4,718, row 17
  labelling 2,140 (E 809, the scalar potential 673, the vector potential 638, the fields 20; 2c-1's were 42). 191 moved
  from rows 18 (148) and 30 (43), because κ or R_eff changes them too; 205 older corpus rows are newly mapped. Of the
  sweep's 260 E rows Rust refused all (InvalidParameter: its |β|² rounded to 1; row 30 makes them valid): 133 are
  now computed and 127 refused as DivisionByZero, by the exact κ. All passed on the harness's first run.
- **`cr_verify.py`'s lw_reff, lw_kappa and their *_slow** (exact families; the slow ones skip the rounding test and
  always decide κ's cut exactly): generic directions and speeds, R along β as |β| → 1, the near-axis ultra-relativistic
  charge, κ's knife edge (β₃ bisected to the last lattice points on either side of 10^-30, each step about 2^-54 of κ),
  constructed ties, the far/near edge, extreme scales and spreads, any finite parts. **600,000 calls, all equal to the
  exact value on x86_64 and aarch64.**
- **Tests**: tests/relativity.tcyr's `cyrius_exact_thresholds (ADR 0002 row 17)` gains 9 assertions: E at κ = 6.0e-17
  (Rust: DivisionByZero), the two sides of κ's cut an ulp of β₃ apart, φ and A on the axis at β = 1 − 2^-53, φ finite
  where Rust's was ∞, φ past DBL_MAX's root, and `_bj_cr_lw` on a tie of κ and one of R_eff by both paths.
- **For M2's range-edge review** (a finding, beside the η-sum guard question): E's velocity term multiplies 1/κ³ by
  n̂ − β, which cancels in the same limit (|n̂ − β|² = 2κ − (1 − |β|²), between κ² and 2κ), and n̂ = R·(1/|R|) carries
  Rust's roundings: off the axis its absolute error of about 2^-53 is a relative one of 2^-53/|n̂ − β|, up to the
  2^-53/κ that κ's was. The agreed scope was κ; n̂ − β stays Rust's.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 1.5–2.1, another repository's verifier
  running beside it; each figure the lower of the three): `relativity/lienard_wiechert_e` 327 → 582 ns,
  `relativity/lienard_wiechert_scalar_potential` 80 → 292 ns, `relativity/lienard_wiechert_fields_orbit_b099` 5.61 →
  5.83 µs (its trajectory's sin and cos dominate). Both twins sit far from the cone's edge (E on the axis at β = 0.1,
  φ at n̂⋅β = 0.45). A probe near it (R = (1, 0.5, 0.25), β = (0.8, 0.41, 0.2), n̂⋅β = 0.89): E 575 → 1,402 ns, φ
  338 → 1,130 ns, the cost of D's exact expansion, the three cross components and five double-double products. The
  first draft cost more (E 693 ns, φ 447 ns): splitting each part once, forming the twelve products inline on the
  halves, and summing S and P with plain TwoSums took `_bj_cr_lw`'s far path from 285 to 208 ns, with every result
  unchanged. A finding, recorded (Rust 1.74–21.64 ns); the other rows held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 400,394 ok, 0 failed, 0
  measured, 0 pending, 38,835 accepted (ADR 01: 8,998; ADR 0002: 29,837)** over 439,229 rows; **all 276 digests
  equal**. Four changed: both potentials, E and the fields.

## 40. M2 part 2c-5: the inverse combinations, an exact zero test and an accurate 1/Σ(1/xᵢ) (2026-11-13)

The last of 2c's boundary values (the user's answer of 2026-11-12: exact tests, accurate values). **2c is done.**

- **ADR 0002 row 17, amended**: `resistance_parallel`, `capacitance_series` and `inductance_parallel` return
  1/Σ 1/xᵢ from the exact sum, rounded once, and refuse |Σ 1/xᵢ| < 10^-30 exactly (DivisionByZero "sum of inverse
  …s is zero"). Rust's element test (|xᵢ| < 1e-30, already exact on doubles) and its order stay; an ±∞ value adds
  1/±∞ = ±0 exactly, as IEEE's; a NaN keeps Rust's chain. Rust's 1/(Σ RN(1/xᵢ)) rounds n + 1 times and cancels when
  the signs differ:
  - two values an ulp apart, 15.887776697373296 and −15.887776697373294, summed to 0 and were refused, where
    1/S = −1.42e17;
  - 1e-30, 3 and −1e-30: 1e30 absorbed the 1/3 and the sum came out 0 (refused), where the value is 3;
  - two values near 2e30 gave 1e30, where Σ 1/x lies below 10^-30;
  - over the new rows, per function, 28 refusals are now values, 7 values are now refused, 3 signs are corrected, and
    189 values were more than 1,000 ulps off (the worst −1.34e8 for +7.06e8).
- **`src/crmath_edge.cyr`**, three tiers, each with a rigorous bound E on |S − S_dd|, the first that decides
  answering. Each 1/xᵢ is q + l, q = RN(1/xᵢ) and l = RN(r·q) from the exact remainder r = 1 − q·xᵢ (TwoProduct;
  2.02u²|q|, exact when r = 0).
  - T1, a compensated sum (TwoSums of the heads, their errors and the lows into one correction): E/|S| grows as
    n²u²/2 for like-signed values.
  - T2, when T1 cannot decide: the exact Shewchuk expansion of every q and l (at most 252 components), about 6u² of S
    whatever n is, and exact when every reciprocal is (sums that cancel to 0 included).
  - The verdict: |S| below or above 10^-30 by the bound against the doubles either side of it, then, with
    E/|S| ≤ 2^-60, kimiya's quotient steps for 1/S with the rounding test's ε widened by E/|S|.
  - T3, the exact path, for what neither decides (a sum near 10^-30, a cancellation below the bounds, a value within
    about 2^-101 of a midpoint, ties among them): each finite xᵢ = ±Mᵢ·2^Eᵢ with Mᵢ odd, heapsorted by M, each group's
    dyadic sum in integers, the groups folded into one fraction N/D on bignums, D the product of the distinct M
    (so [3, −3, …] stays small). |S| < 10^-30 is X·5^30·2^30 < D·2^G; the candidate for 1/|S| = D·2^G/X comes from the
    integers' top 62 bits, and the neighbouring midpoints decide, ties to even (they exist: M₁ and 2^28 − M₁ for odd
    M₁ near 2^27 put 1/S on one). The buffers, sized from the group count, come from fl_calloc and go back by fl_free;
    a refused allocation is the new Allocation error ("the exact … combination's buffers").
  - **A finding, recorded**: T3 is quadratic in the number of distinct significands (D grows by 53 bits a group). Only
    sums the bounds cannot decide reach it: knife edges, deep cancellation, or a value within about 2^-101 of a
    midpoint, which T2 makes about a 2^-48 event for like-signed values at any n.
- **Rows** (`m_circuit.rs`'s `emit_inverse_edges`, on a seed of its own, after every older row): 635 lists through the
  three functions, 1,905 rows: one value near ±10^30 (±8 ulps), k equal values near k·10^30 (k = 2, 3, 7, 15), pairs
  a and −(a + j ulps) whose 1/S crosses 10^-30 (a from 3e14 to 4e15), pairs cancelling to 1–4 ulps at scales 2^±60,
  triples closing on a third value, 60 exact midpoints, 200 mixed-sign lists, 40 lists over ±DBL_MAX, ±1e-30, 1e300
  and ±∞ (all ∞ among them), and lists of 64 and 200 values (positive, alternating ±3, and alternating with the last
  an ulp off). `validate.py`'s transcription agreed with every row on its first run.
- **The parity harness** read a slice into an `f64[8]`: the older rows held at most three values, and the new
  64- and 200-value lists overran it (1,119 mapped rows failed and one went pending on the first run). It now takes
  up to 512 and leaves a longer slice pending rather than read past its buffer.
- **Mappings** (`validate.py`'s `circuit_inverse17`, Fraction): `circuit.adr.tsv` maps 1,293 more row-17 rows (431 a
  function, 150 of them older rows), all passing once the harness was fixed.
- **`cr_verify.py`'s par, par_x, par_slow** (the three tiers; from the expansion on; the exact path alone) over
  lists of up to 15 values (physical, mixed signs, cancelling pairs and triples, knife edges at 10^-30, constructed
  ties, extremes with ±∞, repeated values of alternating sign, any finite values), and **par_rep, par_rep_x,
  par_rep_slow** over short patterns repeated up to 20,000 times (`cr_driver.cyr` builds the list on the heap):
  **452,250 calls, all equal to the exact value on x86_64 and aarch64.** cr_driver now includes lib/freelist.cyr.
- **Tests**: tests/circuit.tcyr's `cyrius_exact_thresholds (ADR 0002 row 17)` gains 9 assertions: the pair an ulp
  apart (Rust: DivisionByZero), the two values near 2e30 (refused), ±1e-30 around 3, an ∞ resistance, fifty 3 and −3
  (Σ = 0 exactly) and the same with the last an ulp off (Rust 11 % low), 64 values 1 + k/7 and the same alternating
  in sign (Rust an ulp low), and a tie by all three tiers.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 2.9–5.5, another repository's benches
  running beside it; each figure the lower of the three): `circuit/resistance_parallel` (five positive values, T1) 30
  → 384 ns. Two probes: 64 values, 241 → 3,453 ns (54 ns a value against Rust's chain's 3.8); a pair cancelling to an
  ulp ([3, −3 an ulp away from zero], decided only by T3), 25 → 3,554 ns. T1 costs about 38 f64 operations a value
  (two TwoSums, a TwoProduct, two divisions' worth), at about 1.5 ns each on this compiler; T3's cost is its bignums
  and buffers. A finding, recorded (Rust 7.09 ns); the other rows held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 401,006 ok, 0 failed, 0
  measured, 0 pending, 40,128 accepted (ADR 01: 8,998; ADR 0002: 31,130)** over 441,134 rows; **all 276 digests
  equal**. Three changed: the three inverse combinations.
- **The user's answers of 2026-11-13** (the end of 2c, the two findings held for it): `material_reflection_coefficient`
  drops its |η₁ + η₂| < 1e-30 guard, since positive relative impedances never sum to zero and Γ depends only on their
  ratio (row 24's rule, as for Schlick's R₀); and `lienard_wiechert_e`'s n̂ − β is rounded once from the exact parts,
  as κ is (§39's finding). Both open 2d.

## 41. M2 part 2d-1: the η-sum guard dropped, Liénard–Wiechert's n̂ − β rounded once (2026-11-13)

2d opens with the user's answers of 2026-11-13 on the two findings held for the end of 2c.

- **ADR 0002 row 25, amended**: `material_reflection_coefficient` no longer refuses |η₁ + η₂| < 1e-30. Both relative
  impedances are positive (RN(√(μ_r/ε_r)) of positive finite values never rounds to 0), so their sum never vanishes,
  and Γ = (η₂ − η₁)/(η₂ + η₁) depends only on their ratio: the guard refused only pairs of tiny impedances whose Γ is
  well defined (row 24's rule, as for Schlick's R₀). Rust's guard sat on the absolute impedances' sum, 376.7 times the
  relative one, so the port's had refused pairs Rust computed (η₁ + η₂ from about 2.7·10^-33 to 10^-30); now
  neither refuses. Γ = 3/5 for η = 2^-112 and 2^-110 (Rust: DivisionByZero).
- **ADR 0002 row 17, amended**: `lienard_wiechert_e`'s n̂ − β is rounded once part by part
  (`src/crmath_edge.cyr`'s `_bj_cr_lw` mode 2, beside κ, sharing its scaling, splits and √S'). It vanishes with κ
  (|n̂ − β|² = 2κ − (1 − |β|²)), and Rust's n̂ = R·RN(1/|R|) carried two roundings into the difference:
  - along (3, 4, 0) with β = 0.6 x̂, n̂ₓ − βₓ = 3/5 − RN(0.6) = 2.2e-17 came out 1.1e-16, Eₓ five times too large;
  - just off the ultra-relativistic axis of §39 (κ = 6.0e-17), n̂ₓ − βₓ is 0.75·2^-53, not 2^-53: Eₓ was a third
    high even with κ exact;
  - of the 903 E rows 2c-4 had mapped and this moves again, 118 move by more than 10^-6 relative, 32 by more than 1 %.
- **The parts** (the section comment): dᵢ = R'ᵢ/√S' − βᵢ by kimiya's quotient steps on √S' 's double-word, βᵢ taken off
  by TwoSum, within 27.4u²|n̂ᵢ| + 1.01u²|dᵢ| + 2^-1060 (an underflowing step); equal nonzero ends of the widened test
  are RN(dᵢ). Otherwise the exact sign of sg·dᵢ − μ is that of R̃ − B√S' (R̃ = sg·R'ᵢ, B = sg·βᵢ + μ), by their signs or
  R̃² − B²S' on bignums; the candidate is the fast estimate when its ends are adjacent, else N/(√S'(|R'ᵢ| + |βᵢ|√S'))
  with N = R'ᵢ² − βᵢ²S' exact (R'ᵢ and βᵢ of one sign, the only cancelling case). Signs and scales come from R
  itself: the first run returned +0 for 80 parts whose scaled R'ᵢ had underflowed to −0 (a tiny negative part, so
  RN(dᵢ) = −0); fixed before any further run. A zero Rᵢ keeps Rust's ±0 − βᵢ, exact. Ties need R'ᵢ/√S' dyadic, which
  only an axis gives.
- **Rows**: `m_wave.rs`'s `emit_eta_guard`, 35 (η from 10^-34 to 10^-29 across both cuts), and `m_relativity.rs`'s
  `emit_nb_edges`, 80 (five Pythagorean directions at scales 2^-40 to 2^300, β = f·RN(n̂) for f = ½, 0.9, 1 − 2^-20 and
  1 − 2^-52, the parts rotated), each after every older row. `validate.py`'s transcription agreed with every row.
- **Mappings**: `wave.adr.tsv` maps 22 more row-25 rows (the 4 older rows Rust refused, η = 10^-300 and 10^-200, and
  18 new ones). `relativity.adr.tsv` grows from 4,718 rows to 5,131: row 17 labels 1,379 E rows (809 under 2c-4; 89
  moved from row 18 and 83 from row 30, since n̂ − β changes them too) and 40 fields rows. All pass.
- **`cr_verify.py`'s lw_nb and lw_nb_slow** (exact families; each part of n̂ − β, the slow one the exact check alone)
  on lw_cases' pairs and on pairs built to cancel (β along n̂ to ulps, Pythagorean R with β at its own n̂'s doubles,
  zero and signed-zero parts): **374,982 calls, all equal to the exact value on x86_64 and aarch64**, and the four lw
  families rerun on the shared code, 600,000 calls, all equal.
- **Tests**: tests/wave.tcyr's row-25 group gains Γ = 3/5 for η = 2^-112 and 2^-110 (Rust: DivisionByZero);
  tests/relativity.tcyr's row-17 group gains E along (3, 4, 0) and the parts of n̂ − β by both paths, and §39's
  near-axis E takes its new Eₓ.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 4.3–5.8: this bite's verifier and another
  repository's tests ran beside it, so both columns are inflated alike; each figure the lower of the three):
  `relativity/lienard_wiechert_e` 639 → 763 ns, the near-cone E probe (§39's) 1,521 → 1,822 ns: the three parts' fast
  paths, about 40 dependent f64 operations each. The potentials and the other rows held within the noise
  (`lienard_wiechert_fields_orbit_b099` 6.30 → 6.75 µs, its rounds spread over 6.3–7.3). A finding, recorded.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 400,686 ok, 0 failed, 0
  measured, 0 pending, 40,563 accepted (ADR 01: 8,998; ADR 0002: 31,565)** over 441,249 rows; **all 276 digests
  equal**. Three changed: E, the fields and the reflection coefficient.

## 42. M2 part 2d-2a: wave's ratio chains at their range edges (2026-11-13)

2d's rule is ADR 0002 row 18's (M1 part 5d): where Rust's chain over- or underflows an intermediate but the result is
representable, the chain runs with an unbounded exponent and its last operation rounds once; a guard keeps Rust's
chain, bit for bit, everywhere else. 2d-2 takes wave's chains in three bites: 2d-2a the ratio chains, 2d-2b the plane
waves (Poynting, momentum density, `plane_wave_e`) and Snell's angle, 2d-2c the material interface, Γ and the
waveguides.

- **ADR 0002 row 18, amended**: every ratio (a·x − b·y)/(a·x + b·y) in wave: `reflectance_normal` and
  `transmittance_normal` ((n₁ − n₂)/(n₁ + n₂)), Schlick's R₀, the four angle-form Fresnel functions, and the
  trig-free family's ratios (`fresnel_rs_direct`, `fresnel_rp_direct`, the two reflectances, the unpolarized pair, and
  through it `reflectance_at_interface`).
  - **The guard** (`src/wave.cyr`'s `_bj_wave_ratio_of`): both products' magnitudes above 2^-1022 and both sums at
    most DBL_MAX, inline; failing that, src/crmath_vx.cyr's full test (a zero factor makes a product exact; a
    non-finite operand keeps Rust's chain). A sum or difference of normal values that lands below 2^-1022 is exact, so
    then Rust's chain is the unbounded one.
  - **The edge** (`_bj_vx_ratio_ux`): the products, the difference and the sum to 53 bits with an unbounded exponent,
    the quotient rounded once; a zero difference gives IEEE's ±0. A ratio's square needs nothing more: a quotient
    rounded once is its 53-bit value unless it lies below 2^-1022, where the square underflows to 0 either way. And
    r_s² + r_p² never nears the underflow: r_s vanishes only for n₁ = n₂, a shortcut, so |r_s| ≳ 2^-53 for distinct
    doubles.
  - **The tests**: the angle forms' |den| < 1e-30 (on the chain's value, the user's 2c answer for trig tests) reads
    the unbounded denominator at an edge, the trig-free family's den = 0 likewise; `reflectance_normal`'s exact
    |n₁ + n₂| < 10^-30 (2c-1) comes first.
  - Rust: DBL_MAX and DBL_MAX/2 gave r = 0 for 1/3 (the sum past DBL_MAX: num/∞) and R = 0 for 1/9; Schlick's R₀
    came out 0 and its R 0.0005 for 0.111; at indices of 10^-310 both products are subnormal, and r_s was 227 ulps
    off, R 309.
- **Rows** (`m_wave.rs`'s `emit_ratio_edges`, after every older row): 2,872 more: `reflectance_normal` and
  `transmittance_normal` over indices at DBL_MAX, its half and third, 10^308, 1.5·10^308 and 9·10^307, a mixed-sign
  half, and over eight subnormal and tiny indices; the trig-free family at those indices with cosines 1, ½ and 10^-310
  (big), 1, 0.3 and 10^-20 (small), and five ordinary pairs with five subnormal and tiny cosines; the angle forms at
  five index pairs (DBL_MAX and its half, ±10^308, subnormal) and five angle pairs. `validate.py`'s transcription
  agreed with every row.
- **Mappings** (`validate.py`'s `wave_ratio_of`: the guard, then exact.py's `rnd53` at each step): `wave.adr.tsv`
  grows from 6,215 rows to 6,876: row 18 labels 463 more (28 older rows; the rest new), and the new rows add 139
  row-24 and 56 row-17 mappings (Schlick's guard, the exact cosines). Row 18 goes to a row only where the edge's chain
  moves it; an angle-form row it moves leaves the corrected oracle for the mappings. All passed on the harness's first
  run.
- **Tests**: tests/wave.tcyr's `cyrius_adr_rows_22_to_25` group gains 6 assertions: R, the angle form's r_s, the
  trig-free r_s and Schlick's R at DBL_MAX and DBL_MAX/2, and the trig-free r_s and R at subnormal indices.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 2.3–3.9, another repository's verifier
  beside it; each figure the lower of the three): `wave/fresnel_rs` 283 → 304 ns, `wave/fresnel_rs_direct` 255 →
  277 ns, `wave/reflectance_unpolarized` 264 → 301 ns, `wave/reflectance_at_interface` 705 → 747 ns,
  `wave/schlick_reflectance` 243 → 265 ns: about 22 ns a ratio, the helper's call and its Result. The first draft
  called crmath_vx's full test (three more calls) on every ratio, 45 to 95 ns; the inline common case took it down.
  A finding, recorded; the other rows held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 402,897 ok, 0 failed, 0
  measured, 0 pending, 41,224 accepted (ADR 01: 8,998; ADR 0002: 32,226)** over 444,121 rows; **all 276 digests
  equal**. Fourteen changed: the functions with new rows.

## 43. M2 part 2d-2b: the plane waves' and Snell's chains at their range edges (2026-11-13)

- **ADR 0002 row 18, amended** (the rule of §42):
  - **The Poynting vector and the momentum density** (`src/wave.cyr`'s `_bj_wave_flux_part`, per part of the
    cross): Rust's cᵢ·RN(1/μ₀) (then ·RN(1/c²)) where the part BjFieldVector's cross gives is normal (it is then the
    chain's 53-bit value), an exact 0 (zero factors, or two normal products), or past DBL_MAX for S (k > 1: S is past
    it too); and, for g, where S is finite. Otherwise the part runs from the products with an unbounded exponent and
    its last product rounds once. A part below 2^-1022 had given S from a degraded value; a cross or S past DBL_MAX had
    made g = ∞ where (E × B)ε₀ is representable: E = 10^159 x̂ and B = 10^160 ŷ give g_z = 8.85·10^307.
  - **`plane_wave_e`**: the argument kx − ωt + φ with an unbounded exponent where a product leaves the normal range
    or a sum passes DBL_MAX, its sine to 53 bits (`_bj_vx_sin53`: t itself below 2^-1022, Payne–Hanek past DBL_MAX),
    and the product by E₀ rounded once. Rust's sin(±∞) was NaN where kx and ωt cancel exactly (10^200·10^200 twice and
    φ = ½: sin ½) or pass DBL_MAX together (sin 10^400 = 0.18651); an underflowed kx left E₀ = 10^300 times sin 0,
    where the value is 10^-100.
  - **Snell's angle**: n₁ sin θ_i below 2^-1022 (or 0 by underflow) runs unbounded and its quotient by n₂ rounds once;
    θ_t at n₁ = 10^-300, n₂ = 10^-29 and θ_i = 10^-20 had been off in its fifth digit. The |sin θ_t| > 1 test reads
    that quotient (the user's 2c answer for trig tests: on the computed value).
- **Rows** (`m_wave.rs`'s `emit_flux_edges`, after every older row): 49 more: 18 Poynting and momentum pairs (crosses
  below 2^-1022 at two scales and mixed signs, S past DBL_MAX, the cross past DBL_MAX with g representable, DBL_MAX
  against 10^-10), nine plane waves (kx and ωt past DBL_MAX cancelling or not, below 2^-1022 under large amplitudes)
  and four Snell rows. `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `wave_flux_part`, and `ux` through the plane wave and Snell, on `rnd53` and
  `mpfr_ref.wide`'s sine): `wave.adr.tsv` maps 88 more row-18 rows (the momentum density 48, the Poynting vector 26,
  `plane_wave_e` 8, Snell 6; 65 older rows among them), all passing on the harness's first run. A class-C row the edge
  moves leaves the corrected oracle for the mappings, as §42's angle forms.
- **Tests**: tests/wave.tcyr's `cyrius_adr_rows_22_to_25` group gains 5 assertions: g past DBL_MAX's S, the three
  plane waves above, and Snell's θ_t.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 2.6–3.1; each figure the lower of the
  three): `wave/poynting_vector` 65 → 95 ns, the guard's call and its tests (a bench vector with zero parts: the
  first draft's inline test took only normal parts, so it called the per-part function three times, 122 ns);
  `wave/plane_wave_e` and `wave/snell_refraction` held within the noise, as every other row. A finding, recorded.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 402,858 ok, 0 failed, 0
  measured, 0 pending, 41,312 accepted (ADR 01: 8,998; ADR 0002: 32,314)** over 444,170 rows; **all 276 digests
  equal**. Four changed: the Poynting vector, the momentum density, the plane wave and Snell's angle.

## 44. M2 part 2d-2c: the material interface's roots and Γ, the waveguides' root (2026-11-13)

- **ADR 0002 row 18, amended** (the rule of §42):
  - **Γ** (`src/crmath_vx.cyr`'s `_bj_vx_gamma`, under `material_reflection_coefficient`,
    `BjMaterial.try_reflectance_at` and so `reflectance_normal_materials`, and the interface's equal-index branch):
    Rust's (η₂ − η₁)/(η₂ + η₁) where both correctly rounded η's lie in (2^-1022, DBL_MAX] and their sum is finite
    (each is then its 53-bit value, and a difference of normal values below 2^-1022 is exact). Otherwise the η's are
    taken to 53 bits with an unbounded exponent (`_bj_vx_sqrt_ux`: the correctly rounded root of the significands'
    quotient, the exponents' difference made even first) and the chain's quotient rounds once. The lower bound is
    strict: a root just below 2^-1022 rounds up to it in binary64 but not in 53 bits, and the corpus has that pair
    (Γ = 2^-54 where the f64 η's are equal). Rust: an η of ∞ (ε_r = 2^-1074 under μ_r = DBL_MAX) made Γ NaN (∞ − ∞);
    η's of 2^1023.5 and 2^1023 summed past DBL_MAX for Γ = −0 where it is −0.17157; subnormal η's were refused by
    Rust's absolute guard, and the port's Γ on them had kept only their few bits.
  - **The interface** (`bj_reflectance_at_interface`): an index at or below 2^-1022 sends the TIR test, cos θ_t and
    the equal-index test to the 53-bit indices (`_bj_vx_sqrt_ux`'s √(ε_r μ_r)), scaled together into the normal range
    by the larger's binade, the smaller clamped at 2^-200 below it. Those three depend only on n₁/n₂, and past a ratio
    of 2^27 both decide alike (TIR below grazing incidence, or cos θ_t = 1). An η outside (2^-1022, DBL_MAX] sends
    r_s and r_p to `_bj_vx_ratio_ux2`, the ratio chain on the 53-bit η's. Rust's n = √(ε_r·μ_r) underflowed to 0 for
    every such pair and refused it (n₂ = 0). The port's correctly rounded n's (row 13) were equal for ε_r μ_r of
    2^-2148 and 2^-2147, giving Γ² = 0.0294 at cos θ_i = 0.5 where R = 0.0745. An η of ∞ gave NaN where R = 1.
  - **The rectangular waveguides** (`_bj_wave_rect_root`, which now also reports whether Rust's chain is the
    unbounded one): each square above 2^-1022 or an exact 0, their sum finite and, for f_c, v/2 above 2^-1022.
    Otherwise √((m/a)² + (n/b)²) runs with an unbounded exponent, f_c = (v/2)·r and λ_c = 2/r round once, and a λ_c
    past DBL_MAX is refused with Rust's "not representable" text. Rust: b = 10^-300 gave f_c = ∞ for 1.5·10^308 and
    λ_c = 0 for 2·10^-300; sides of 2^600 m (λ_c = 2^601) and 8·10^307 m (1.6·10^308) were refused; and
    v = 2^-1021 − 2^-1074 had its half tie up to 2^-1022, moving f_c by an ulp (157 rows).
- **Rows**: 5,327 more. `m_material.rs`'s `edge_materials` has 19 lossless materials: η past DBL_MAX three ways,
  η's of 2^1023.5 and 2^1023, subnormal η's three ways, an η of 2^-1022 rounded up and one exactly, n subnormal four
  ways, an n of 2^-1022 rounded up and one exactly, and vacuum, glass and a tiny normal pair as controls. Every ordered
  pair gives:
  - `try_reflectance_at`: 361 rows in material.tsv;
  - `emit_material_edges` (m_wave.rs): 2,166 rows of Γ, R at normal incidence, and R at the interface at
    cos θ_i = 1, 0.7, 0.2 and 0.

  `emit_guide_edges` adds 2,800 waveguide rows:
  - ten sides: 2^511 and its neighbours; 2^-512; 0.625·2^-511, whose squares' sum overflows; 0.75·2^-511, whose sum
    stays finite; 8·10^307 and 10^308; and two ordinary sides;
  - four modes;
  - six velocities, among them 2^-1021 − 2^-1074 (its half ties to 2^-1022) and 2^-1021.

  `validate.py`'s transcription agreed with every row.
- **Mappings**: `validate.py`'s `wave_gamma_port`, `wave_scaled_indices`, `wave_ratio_ux` (wave_ratio_of's unbounded
  chain, now on 53-bit values) and the waveguides' `ux` chain, all on exact.py's `rnd53` and `rnd53_sqrt`.
  `wave.adr.tsv` maps 3,117 more rows:
  - row 18, 2,000: the cutoff 1,200, λ_c 252, the interface 314, Γ 117, R at normal incidence 117;
  - row 13, 1,012: the interface 881, R at normal incidence 73, Γ 58;
  - row 25, 105.

  `material.adr.tsv` maps 190 more (`try_reflectance_at`: row 18 117, row 13 73). 299 of the 3,307 are older
  waveguide rows (the cutoff 168, λ_c 131). All pass. Labels: row 18 where the edge alone moves the row, checked
  first, as in §42.
- **cr_verify** (`vx_sqrt`, `vx_sqrt_e` and `vx_gamma`; `cr_driver.cyr` now includes src/crmath_vx.cyr): 750,000
  calls equal to the exact value on x86_64 and aarch64. `vx_gamma`'s expected value is the unbounded chain alone, so
  its 150,000 calls check the guard's claim wherever Rust's quotient is taken.
- **Tests**:
  - tests/wave.tcyr's `cyrius_adr_rows_22_to_25` gains 9 assertions: Γ past DBL_MAX and between subnormal η's;
    Γ = 2^-54 at the rounded-up 2^-1022; R = 1 from an η past DBL_MAX; the subnormal indices' R; the rounded-up n at
    grazing incidence; λ_c = 2^601 and 1.6·10^308; and f_c where Rust gave ∞. Its 2^600-side refusal became a
    10^308-side one (λ_c = 2·10^308).
  - tests/material.tcyr gains 1: `try_reflectance_at` with η's summing past DBL_MAX.
- **Benchmarks**: an interleaved A/B against HEAD, three rounds each, load 1.9–3.4, each figure the lower of the
  three.
  - `wave/rect_waveguide_cutoff`: 46 → 51 ns, steady across the rounds. The cost is the guard's three tests and its
    flag through a pointer (the guard shares the root's quotients and squares).
  - `wave/reflectance_at_interface`: 715 → 755 ns. HEAD's three rounds spread over 715–755, so a focused A/B of the
    two rows ran too (five repetitions a round, three rounds, load 2.7–3.0). It puts the cost at about 740 → 757 ns,
    the six range tests on the indices and η's.
  - Every other row held within the noise: untouched rows drifted 2–4 % with the load. A finding, recorded.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 404,878 ok, 0 failed, 0
  measured, 0 pending, 44,619 accepted (ADR 01: 8,998; ADR 0002: 35,621)** over 449,497 rows; **all 276 digests
  equal**. Six changed: Γ, both normal-incidence reflectances, the interface and the two waveguide functions.

## 45. M2 part 2d-3a: circuit's series folds and impedance at their range edges (2026-11-13)

2d-3 takes §27–§29's remaining range-edge candidates in bites. The resonances, Q, ζ, ½CV², ½LI², μn²Al, ρL/A and
εA/d were already rounded once on exact monomials (row 13, §34), which reach both ends of the range: the corpus maps
L = C = 2^-1022 to f₀ = 7.2·10^306 where Rust gave ∞. The audit's other circuit items:

- **ADR 0002 row 18, amended** (the rule of §42):
  - **The series folds** (`_bj_circuit_sum`, under `resistance_series`, `capacitance_parallel` and
    `inductance_series`): where partial sums of finite values pass DBL_MAX, Rust's left fold ends ±∞, or NaN from
    ∞ − ∞, though the sum is representable. The fold then runs with an unbounded exponent, each step to 53 bits and
    the last rounded once. The guard is Rust's own result: a finite fold had no partial past DBL_MAX, and a sum of
    doubles loses nothing to the underflow (a sum below 2^-1022 is exact). Absorption stays as Rust's
    (DBL_MAX + 1 − DBL_MAX = 0 in both chains), and a non-finite value keeps Rust's fold. Rust: DBL_MAX + DBL_MAX −
    DBL_MAX = ∞ for DBL_MAX; 1.5·10^308 twice less 10^308 three times ∞ for +0.
  - **The impedance** (`bj_rlc_impedance`): Rust's √(R² + (ωL − 1/(ωC))²) where its sum R² + X² lies in
    [2^-900, DBL_MAX], src/crmath_vx.cyr's band. A product or quotient below 2^-1022 is then absorbed by its partner,
    and X is the chain's wherever X² is not negligible. Otherwise, for finite operands, ωL, ωC, 1/(ωC), X, R², X² and
    their sum run to 53 bits with an unbounded exponent and √ rounds once (`_bj_k_ux_sqrt_f64`). Rust: R = 10^200 and
    ωL = 10^300 gave |Z| = ∞ (a square past DBL_MAX); R = 10^-200 at an exact resonance gave 0 (R² underflowed);
    ωC past DBL_MAX left 1/(ωC) = 0, which the band absorbs wherever it does not matter.
- **Rows** (`m_circuit.rs`'s `emit_fold_edges`, its own seed, after every older row): 1,198 more.
  - 70 slices for the three series folds: ten built (cancelling partials past DBL_MAX, an absorption control, ∞ and
    NaN controls) and 60 drawn from values near ±DBL_MAX, 1, −0 and 2^-1074.
  - A 960-row impedance grid: R from 0 and 2^-1074 to DBL_MAX; L from 0 to DBL_MAX; C from 2·10^-30 to DBL_MAX;
    ω from 2·10^-30 to DBL_MAX, and −10^3.
  - 28 exact resonances: L = 2^(90−2k), C = 2^-90 and ω = 2^k for k = 100, 300, 500 and 580, under six R from 0 to
    10^200, and L an ulp off.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `circuit_combo` and `circuit_rust` with `ux`, on exact.py's `rnd53`):
  `circuit.adr.tsv` maps 886 row-18 rows: the impedance 844 (353 of them older rows) and each series fold 14. All
  pass on the harness's first run.
- **Guard audit** (a sampling check in Python, not a gate): on random operands of every binade, wherever a guard
  keeps Rust's chain its result equals the unbounded chain's. That held for 44,613 impedances (half with ordinary R
  and extreme L, C and ω) and 132,940 folds.
- **Tests**: tests/circuit.tcyr gains a `cyrius_range_edges (ADR 0002 row 18)` group of 4 assertions: the two folds
  above, R = 10^-200 at resonance, and R = 10^200.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 1.8–2.2, each figure the lower of the
  three): `circuit/rlc_impedance` 16 → 17 ns, the band test on Rust's sum (the record's 17.00 ns stands). No bench
  covers a series fold, whose guard adds one subtraction and one test. Every other row held within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,190 ok, 0 failed, 0
  measured, 0 pending, 45,505 accepted (ADR 01: 8,998; ADR 0002: 36,507)** over 450,695 rows; **all 276 digests
  equal**. Four changed: the three series folds and the impedance.
- **Next (2d-3b)**: the transients' edges, a finding of this bite's review. τ = RC or L/R past DBL_MAX (Rust's t/τ is
  then 0), t/τ below 2^-1022, e^x or expm1 past either end under an amplitude that brings the product back, and
  V/R subnormal or past DBL_MAX.

## 46. M2 part 2d-3b: the transients at their range edges (2026-11-13)

A finding of §45's review: the four transients compute amp·w with x = −t/τ, w = −expm1(x) for a charge or rise and e^x
for a discharge or decay, and each step can leave the range while the result stays representable.

- **ADR 0002 row 18, amended** (the rule of §42):
  - **The guard** (`_bj_circuit_plain`, a few tests on Rust's own values): τ finite (past the instant test it is above
    10^-30), x in (2^-1022, DBL_MAX] or t = 0, w finite (for e^x, in (2^-1022, DBL_MAX]), and for the rise V/R in
    (2^-1022, DBL_MAX] or V = 0. Each step is then its 53-bit value and the product rounds once in both chains. A
    non-finite operand keeps Rust's chain.
  - **The edge** (`_bj_circuit_edge`): τ = RC or L/R, x, V/R and w to 53 bits with an unbounded exponent, the product
    rounded once.
    - |x| below 2^-60: w = −x or 1.
    - |x| from 4,096: a charge's w = 1 for x < 0. Otherwise the product overflows, or (a decay with x < 0) rounds to
      zero, since the amplitude lies between 2^-2098 (V/R) and 2^1024 (an input).
    - Between, x is a double: bijli's correctly rounded expm1 or exp where w is normal, else e^x by `_bj_ux_exp_lim`
      (past 709, −expm1 and −e^x round alike).
    - A zero amplitude gives IEEE's signed zero against w's sign (Rust's 0·∞ was NaN). t = 0 under a τ past DBL_MAX
      gives the signed zero, or the amplitude.
  - **`_bj_ux_exp_lim`** (src/crmath_ux.cyr): `_bj_ux_exp` with its stand-in beyond a chosen |x| (`_bj_ux_exp` keeps
    1500). The transients use 2300, since V/R = 2^-2098 times e^x is representable up to x ≈ 2147. The reduction holds
    there: k ≤ 3319 times the fixed-point ln 2 errs below 2^-352.
  - **Rust's results at these edges**:
    - RC = 10^400 gave V = 0 for 10^200 (τ = ∞, t/τ = 0).
    - t/RC = 10^-310 left V a part in 10^13 off (a subnormal x).
    - 10^300·e^-720 kept only the subnormal e^x's 36 bits.
    - 10^300·e^-800 = 3.7·10^-48 underflowed to 0, and 10^-300·e^750 = 5.3·10^25 overflowed to ∞.
    - V/R = 2^-2074 against e^2000 gave NaN (0·∞) for −1.8·10^244.
    - V₀ = 0 against e^1000 gave NaN for −0.
- **Rows** (`m_circuit.rs`'s `emit_transient_edges`, after every older row): 1,584 more, through all four transients.
  The grid takes six amplitudes (±10^300, 1, 10^-300, 2^-1074, 0) and six (R, C or L) pairs:
  - τ = RC past DBL_MAX with τ = L/R = 1 (10^200 twice, 2^1000 twice);
  - τ = L/R and V/R past DBL_MAX (R = 10^-29);
  - τ = 10^10;
  - and 1.

  It crosses them with eleven times, which put e^(−t/τ) subnormal, below it, past DBL_MAX and past e^1500, and t/τ
  subnormal or past DBL_MAX. `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `circuit_plain` and `circuit_edge`: −expm1 and e^x from MPFR's widest range at every x,
  not the port's shortcuts; past |x| = 10^5 the outcome is certain): `circuit.adr.tsv` maps 504 more rows.
  - Row 18, 277: the rise 101, the charge 80, the decay 64, the discharge 32. 70 of them are older rows, mostly zero
    amplitudes where Rust's 0·∞ was NaN; four of those had been row-29 rows.
  - Row 29, 227 among the new (−expm1 alone moves them).

  All pass on the harness's first run with the final model. Labels: row 18 where the edge moves the row, checked
  before rows 17 and 29.
- **Guard audit** (a sampling check in Python, not a gate): 102,528 transients on random operands of every binade,
  98,678 of them on Rust's chain. The port's result equalled the pure unbounded chain on every one.
- **cr_verify** (`ux_exp_w`, new): `_bj_ux_exp_lim` at 2300 on 200,000 arguments, two thirds between 1500 and 2300,
  each read as its significand and exponent; with `ux_exp`, 800,000 calls equal to MPFR (widest range) on x86_64 and
  aarch64.
- **Tests**: tests/circuit.tcyr's `cyrius_range_edges` group gains 7 assertions, the cases above.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 0.8–1.0, each figure the lower of the
  three): `circuit/rc_charging` 65 → 67 ns, `circuit/rl_current_rise` 65 → 72 ns. The other rows held, among them
  `complex_exp_edge`, whose `_bj_ux_exp` now goes through `_bj_ux_exp_lim`.
  - **A finding** (the pinned compiler): a function pays for each f64 local it declares, executed or not, about
    0.6 ns a call. HEAD's rise with 16 dead f64 locals in an unreachable branch ran 65 → 75 ns; three dead tuple-binding
    calls cost nothing.
  - The first draft carried the guard and the edge inline and ran the rise at 85 ns. The four functions now keep only
    the common test (three |·| locals against the bounds) and call one out-of-line slow path,
    `_bj_circuit_slow`, which recomputes Rust's chain, applies the full test and then the edge.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 406,270 ok, 0 failed, 0
  measured, 0 pending, 46,009 accepted (ADR 01: 8,998; ADR 0002: 37,011)** over 452,279 rows; **all 276 digests
  equal**. Four changed: the four transients.

## 47. M2 part 2d-3c: relativity's four-vector and tensor chains at their range edges (2026-11-13)

§28's audit named FourVector's dot, `interval_sq`, `spatial_magnitude` and `four_force` as raw chains. Reviewing them
turned up three more of the same kind: `boost_x` and the two invariants. The field transforms, the relativistic Larmor
power and an ∞ R in E are 2d-3d's.

- **ADR 0002 row 18, amended** (the rule of §42; each guard tests Rust's own values, the edge is out of line):
  - **`dot` and `interval_sq`**: Rust's −a₀b₀ + a₁b₁ + a₂b₂ + a₃b₃ in the band, or an exact 0 of exact products, or a
    non-finite operand. Otherwise `_bj_rel_dot4_edge`: the products and sums to 53 bits with an unbounded exponent,
    the last sum rounded once.
  - **`spatial_magnitude`**: BjFieldVector `magnitude`'s guard and edge (`_bj_vx_sumsq`, √ rounded once).
  - **`boost_x`** (`_bj_rel_boost_edge`, per output): γ(p − β·q) keeps Rust's chain where β = v/c is normal or 0 and
    β·q is exact or the difference lies in the band. Otherwise β, β·q and the difference run unbounded and the product
    by γ rounds once. A subnormal β is coarse, and its product with a vast q is not negligible, so β must be exact.
  - **`four_force`** (`_bj_rel_ff_edge`, per row): each row's sum in the band, else the products and sums unbounded,
    the product by q rounded once.
  - **The first invariant** 2(B² − E²/c²): Rust's difference in the band. A degraded E (from F·c) or B² is then
    negligible: its square lies far below the difference. Otherwise E = c·F, E², B², E²/c² and the difference run
    unbounded and the doubling rounds once (`_bj_rel_inv1_edge`).
  - **The second invariant** −E⋅B/c: Rust's where each E = c·F is exact and E⋅B is normal. A degraded E is not
    negligible here, since B multiplies it. Otherwise E⋅B runs unbounded and the quotient by c rounds once
    (`_bj_rel_inv2_edge`). An exact 0 keeps IEEE's signs where every product is 0, or gives a cancellation's +0.
  - **Rust's results at these edges**:
    - a light-like four-vector at 10^300 gave s² = NaN (∞ − ∞) for 0;
    - |(10^200, 10^200, 10^200)| = ∞ for 1.73·10^200;
    - E⋅B = 10^309 gave −∞ for −3.3·10^300;
    - a radiation field at B = 10^200 gave the first invariant NaN for 0;
    - q = 10^-300 against four-force sums past DBL_MAX gave (0, ∞, NaN, −∞) for (0, 1.5·10^10, 5·10^9, −2·10^10).
- **Rows** (`m_relativity.rs`'s `emit_fv_edges`, after every older row): 215 more.
  - Ten four-vectors through `interval_sq`, `spatial_magnitude`, `dot` and `boost_x`: products past DBL_MAX that
    cancel (light-like at 10^300, DBL_MAX parts), squares that overflow alone (10^155, 10^200), squares below 2^-1022
    (10^-160, 10^-300). The boosts use speeds whose β is subnormal (10^-310, 2^-1074, 10^-300) as well as ±0.5c, −0.9c
    and 0.
  - Seven tensors through both invariants and the four-force (four-velocities at 10^300 and 10^10, the latter under
    q = 10^-300): a radiation field at 10^200, E ∥ B at 10^150 and at 10^154/10^155, E/c subnormal under a vast B, a
    subnormal B², DBL_MAX parts.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `rel_dot4`, `rel_boost_edge`, `rel_ff_row` and the invariants' chains under
  `rel_ops`'s new `ux`, on exact.py's `rnd53`): `relativity.adr.tsv` maps 203 more rows net.
  - Row 18, 194 net: the first invariant 88, the four-force 31, the second invariant 28, `dot` 28,
    `spatial_magnitude` 13, `interval_sq` 8, `boost_x` 6. 144 are older rows, eight of them second-invariant rows
    whose outputs changed: E⋅B past DBL_MAX had been ∞, and is now finite after the quotient.
  - Row 30, 9 among the new (`boost_x`, its exact γ alone).

  All pass. Labels: as §28's precedence, 27 → 17 → 13 → 18 → 30.
- **Guard audit** (a sampling check in Python, not a gate): on random operands of every binade, wherever a guard keeps
  Rust's chain its result equalled the pure unbounded chain's. Samples kept: dot 84,761, `boost_x` 80,946, the
  four-force 77,782, the first invariant 18,190, the second 80,053.
- **Tests**: tests/relativity.tcyr gains a `cyrius_range_edges (ADR 0002 row 18)` group of 5 assertions, the cases
  above.
- **A pointer note** (the pinned compiler): arithmetic on a typed pointer is scaled by the element size (a `*f64`
  plus 4 is the next row of four), as field.cyr's `charges + i * 4` relies on.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 1.1–1.9, each figure the lower of the
  three): `relativity/four_vector_boost` 155 → 161 ns, the guard's three tests (the record's 168 ns, from a noisier
  run, becomes 161). `relativity/em_tensor_from_fields`, whose code is unchanged, read 31 → 34 ns within HEAD's own
  31–35. Every other row held, `lienard_wiechert_e` among them. No bench covers `dot`, `interval_sq`, the invariants
  or the four-force, whose guards are one band test on a value Rust already computed.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 406,282 ok, 0 failed, 0
  measured, 0 pending, 46,212 accepted (ADR 01: 8,998; ADR 0002: 37,214)** over 452,494 rows; **all 276 digests
  equal**. Seven changed: `dot`, `interval_sq`, `spatial_magnitude`, `boost_x`, the four-force and both invariants.

## 48. M2 part 2d-3d: the field transforms and the relativistic Larmor power at their range edges (2026-11-13)

- **ADR 0002 row 18, amended** (the rule of §42), with a box guard instead of the band:
  - **Why a box.** These chains multiply intermediates by γ (up to 2^27, or far more from a three-part v near c) and
    by v (up to 2^28.2): (E + v × B)·γ, k·(v⋅E)·v, (B − (v × E)/c²)·γ, pre·γ⁶·(a² − |v × a|²/c²). A subnormal
    intermediate is then amplified, not absorbed, so Rust's own result can sit in the band and still be wrong. The
    guard bounds the inputs instead, as the stress tensor's does (§24).
  - **The transforms** (`lorentz_transform_fields`, `lorentz_transform_x`): every nonzero E, B and v part (the speed)
    has a biased exponent in [723, 1923], i.e. lies in [2^-300, 2^901), and γ < 2^27. Then every product and quotient
    of the chain is normal or an exact 0: the smallest, v·k·(v⋅E), is at least 2^(-300-57.3-652) = 2^-1009.3; the
    largest, γ(E + v × B), is at most 2^(27+930.3) = 2^957.3. A cancellation is an exact 0, or at least an ulp of its
    operands (2^-952 for the smallest).
  - **The Larmor power** (`relativistic_larmor_power`): q in [2^-200, 2^401), a in [2^-200, 2^481), the nonzero v parts
    at least 2^-200, and γ < 2^27. The tightest are |v × a|² ≤ 2^1022 and (q²/(6πε₀c³))·γ⁶ ≤ 2^912.
  - **Outside the box**, for finite operands, the whole chain runs with an unbounded exponent: γ + 1, (γ + 1)c², γ², k
    and γ⁶ included, each output's last operation rounded once (`_bj_rel_ltf_edge`, `_bj_rel_ltx_edge`,
    `_bj_rel_larmor_edge`). It agrees with Rust's chain wherever that is exact, so the box need only be sufficient.
    A non-finite operand keeps Rust's chain.
  - **Rust's results at these edges**:
    - E = 10^308 along a γ = 2 boost gave E′ = NaN: γE and the γ²(v⋅E)v/((γ + 1)c²) term were both ∞, while E′∥ = E.
      The port gives 10^308 less an ulp, each step rounded.
    - vE_z past DBL_MAX gave B′ = ∞ where γ(B_y + vE_z/c²) = 3·10^292.
    - a = 10^170 gave P = ∞ for 2.2·10^287 (a² past DBL_MAX).
    - q = 0 against an a² past DBL_MAX gave NaN for 0.
    - Many rows of a² − |v × a|²/c² with both terms ∞ gave NaN, or −∞, where P is +∞ (it is at least a²(1 − β²)·…).
- **An ∞ R in `lienard_wiechert_e`** (§28's audit): the sibling potentials give NaN too (|R| − R⋅β is ∞ − ∞, or R·0),
  and Rust refuses only a NaN R (`r_mag.is_nan()`). A non-finite operand keeps Rust's chain, so this bite leaves it.
  The user's answer of 2026-11-13: all three refuse an ∞ R (InvalidParameter, a text of its own), in 2d-3e. With it,
  the retarded-time solver reports its documented bracket error for a non-finite state met while the bracket expands
  (a superluminal trajectory's position overflowing), not "trajectory returned a non-finite state".
- **Rows** (`m_relativity.rs`'s `emit_transform_edges`, after every older row): 480 more.
  - 180 `lorentz_transform_fields` rows: six boosts (γ = 2, 1.25, diagonal, near c, a subnormal v and tiny parts),
    six E (10^308 and DBL_MAX along the boost, 10^301 across it, subnormal, tiny, ordinary) and five B.
  - 150 `lorentz_transform_x` rows: the same fields at five speeds (0.8c, 0.5c, 2·10^8, 10^-150, 10^-310).
  - 150 Larmor rows: five charges (1.6·10^-19, 10^-200, 10^300, 0, 10^-300), six velocities and five accelerations
    (10^170, 10^150 twice, 10^-200, subnormal, mixed).

  `validate.py`'s transcription agreed with every row.
- **Mappings**: `validate.py`'s new `uz_*` chain keeps IEEE's signed zeros as kimiya's operations do: a product's
  zero takes the xor of the signs, two zeros add as IEEE's, an exact cancellation is +0. `relativity.adr.tsv` maps:
  - row 18, 706 more: the Larmor power 551, `lorentz_transform_fields` 107, `lorentz_transform_x` 48. 538 are older
    rows, mostly the Larmor power's NaN and wrong-signed ∞.
  - row 30, 152 more net: the γ-sensitive rows (port ≠ the vector-operations-only model), 58 of them older rows whose
    outputs changed.

  All pass.
- **Box audit** (a sampling check in Python, not a gate): inside each box, weighted to its edges, Rust's chain equalled
  the pure unbounded chain on every sample. `lorentz_transform_fields` 40,000, `lorentz_transform_x` 36,573, the
  Larmor power 40,000.
- **Tests**: tests/relativity.tcyr's `cyrius_range_edges` group gains 4 assertions, the cases above.
- **Benchmarks** (interleaved A/B against HEAD, three rounds each, load 0.7–1.9, each figure the lower):
  - `relativity/lorentz_transform_x`: 164 → 186 ns. A focused A/B of both transforms (five repetitions a round)
    also gave `lorentz_transform_fields` 367 → 407 ns. It has no bench twin; the focused bench gave it the same fields
    and a 0.5c boost along x.
  - The first box, three looping calls and a one-element array on the hot path, cost 166 → 214 ns in the full A/B.
    One straight-line integer test per caller, on untyped struct addresses (no pointer-conversion calls), with the
    speed and the charge compared as f64s, brought it to 186. The remaining cost is the box's nine or five tests.
  - Every other row held. No bench covers the Larmor power.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,904 ok, 0 failed, 0
  measured, 0 pending, 47,070 accepted (ADR 01: 8,998; ADR 0002: 38,072)** over 452,974 rows; **all 276 digests
  equal**. Three changed: both transforms and the Larmor power.

## 49. M2 part 2d-3e: an ∞ R refused, the solver's bracket error, the potentials' tails (2026-11-13)

On the user's answers of 2026-11-13 (asked at the end of §48), and the potentials' range edges this bite's review found.
E's own tail (γ²κ³R², pre/dv, κ³Rc) is 2d-3f's, with the solver's light-cone sample and `lienard_wiechert_fields`'
R = r_obs − r(t_r) past DBL_MAX: their chains amplify an intermediate as the transforms' do.

- **ADR 0002 row 27, amended** (the user's answer): the three Liénard–Wiechert functions refuse an R with an infinite
  part, InvalidParameter "R vector must be finite", right after the NaN test (a NaN beside the ∞ keeps the NaN text).
  Rust's |R| − R⋅β was ∞ − ∞ and R·(1/|R|) ∞·0: NaN, silently. With NaN and ∞ refused, R is always finite past the
  checks, so the Rust-chain fallbacks for a non-finite R (E's κ = 1 − n̂⋅β with its |κ| < 10^-30 test, the effective
  distance's |R| − R⋅β) are gone.
- **ADR 0002 row 31, new** (the user's answer): `solve_retarded_time` reports its documented "could not bracket the
  light-cone root (is the trajectory faster than light?)" for a non-finite state met while the bracket still expands.
  The doubling drove the trajectory past the range, as a superluminal one does. The first two samples and the
  refinement keep "trajectory returned a non-finite state at t = …", and `lienard_wiechert_fields` follows.
- **ADR 0002 row 18, amended** (the rule of §42):
  - **φ** keeps Rust's chain where R_eff (row 17's, rounded once) is finite: 4πε₀R_eff is then normal (R_eff ≥ 2^-367)
    and the quotient rounds once.
  - **R_eff past DBL_MAX** (|R| near DBL_MAX) is taken to 53 bits with an unbounded exponent: `_bj_cr_lw`'s new mode 3
    returns its scaled R_eff′ = R_eff·2^-k (normal) and k, the value before mode 0's final `_bj_g_scale`. 4πε₀R_eff
    then runs unbounded too and the quotient rounds once.
  - **A** keeps Rust's chain where μ₀q, μ₀qc and their quotient by 4πR_eff are normal (or q = 0) and R_eff is finite;
    otherwise the chain runs unbounded and each part's product by βᵢ rounds once.
  - **Rust's results at these edges**:
    - |R| = √3·DBL_MAX gave φ = 0 for 2.9·10^-299;
    - q = 10^300 at |R| ≈ 2.5·10^308 gave 0 for 54.7;
    - q = DBL_MAX gave A = (−∞, NaN, NaN) (0·∞ across β) for (−∞, 0, 0).
- **Rows** (`m_relativity.rs`'s `emit_lw_inf_edges`, after every older row): 186 more.
  - 60 rows of five R with an infinite part (one beside a NaN) against three β, through all three functions.
  - 72 rows of R near DBL_MAX against two β under four charges.
  - 54 rows of q = DBL_MAX, 10^-310 and 2^-1074 at three distances.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`exact.py`'s new `rnd53_lw`, rn_lw's method on the 53-bit grid with ties to even; the potentials'
  tails on the `uz_*` chain): `relativity.adr.tsv` maps 160 more rows net.
  - Row 27, 54: E 30, each potential 12.
  - Row 18, 106: the scalar potential 52, the vector potential 48, E 6. 63 of them are older rows; 14 older row-17
    rows are now labelled 18.
  - Row 31, 2: the superluminal trajectories.

  All pass.
- **cr_verify** (`lw_reff_ux`, new: mode 3 as its significand and exponent, on lw_cases' pairs and on them scaled by
  2^k toward DBL_MAX): 415,724 calls (207,862 values, each its significand and exponent) equal to the exact value
  on x86_64 and aarch64. With them every older lw family ran again on the changed function (lw_reff, lw_kappa, lw_nb
  and their slow paths), all equal: 1,390,706 calls on each target.
- **Guard audit** (a sampling check in Python, not a gate): the vector potential's guard kept Rust's chain on
  161,745 samples, all equal to the pure unbounded chain.
- **Tests**: tests/relativity.tcyr gains 6 assertions:
  - the two refusals of an ∞ R, and the NaN text beside an ∞;
  - φ at √3·DBL_MAX and under q = 10^300;
  - A under q = DBL_MAX.

  The solver's superluminal case now asserts the bracket text.
- **Benchmarks** (an interleaved A/B against HEAD, three rounds each, load 2.4–4.0 with this bite's cr_verify running
  in the first two, each figure the lower of the three):
  - `relativity/lienard_wiechert_e`: 705 → 662 ns. E no longer computes Rust's κ = 1 − n̂⋅β and n̂ − β only to
    overwrite them with row 17's exact values; the ∞ R test costs less than they did. A focused A/B at load 1.0 (five
    repetitions a round, three rounds) gave 679 → 653 ns (−4 %). The record's 763 ns, from 2d-1, becomes 662.
  - `lienard_wiechert_scalar_potential` 303 → 297 ns, `solve_retarded_time_uniform_t1s` 812 → 808 ns,
    `lienard_wiechert_fields_orbit_b099` 6,187 → 6,014 ns: within the noise. Every other row held too.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,930 ok, 0 failed, 0
  measured, 0 pending, 47,230 accepted (ADR 01: 8,998; ADR 0002: 38,232)** over 453,160 rows; **all 276 digests
  equal**. Four changed: E, both potentials and the solver.

## 50. M2 part 2d-3f: lienard_wiechert_e's tail at its range edges (2026-11-13)

§49 left E's own tail: dv = γ²κ³R², fv = (q/(4πε₀))/dv, the velocity term (n̂ − β)·fv, the cross products
n̂ × ((n̂ − β) × β̇), da = κ³Rc, fa = (q/(4πε₀))/da and the acceleration term. The solver's light-cone sample and
`lienard_wiechert_fields`' R = r_obs − r(t_r) past DBL_MAX are 2d-3g's.

- **ADR 0002 row 18, amended** (the rule of §42, with a box as §48's):
  - **Why a box.** The chain multiplies by up to 1/(γ²κ³R²) (to 2^498) and by β̇ (to 2^201), so a subnormal
    intermediate is amplified, not absorbed.
  - **The box**: every nonzero q, R part and β̇ part has a biased exponent in [923, 1223], i.e. lies in
    [2^-100, 2^201). Every nonzero n̂ − β part (row 17's, bounded by 2) is at least 2^-100, and γ² < 2^100. κ needs
    nothing: past its 10^-30 test, κ³ ∈ [2^-299, 8].
  - **Within the box** every intermediate is normal or an exact 0:
    - n̂ = R·(1/|R|) ≥ 2^-301.8;
    - (q/(4πε₀)) ∈ [2^-67, 2^234.1];
    - the velocity factor ∈ [2^-573.6, 2^732.5];
    - (n̂ − β) × β̇, a difference of normal products, is 0 or at least 2^-252, and n̂ × that at least 2^-605.8;
    - the acceleration factor ∈ [2^-300, 2^604.6] and its term ≥ 2^-905.8.

    Each term is then its 53-bit value and each part's sum rounds once in both chains.
  - **Outside it**, for finite q and β̇, `_bj_rel_e_edge` runs the whole tail with an unbounded exponent: |R| (the
    sum of squares and its root), 1/|R|, n̂, κ³, q/(4πε₀), γ²κ³R², both factors, both cross products (Rust's order)
    and both terms. Each part's sum rounds once. κ, n̂ − β (row 17) and γ² (row 30) enter as the doubles their rows
    define. A non-finite q or β̇ keeps Rust's chain.
  - **Rust's results at these edges**:
    - q = 10^300 at R = 10^200 gave E = NaN for 9·10^-91 (R² past DBL_MAX: dv = ∞, and pre past DBL_MAX too);
    - |R| = √3·DBL_MAX gave 1/|R| = 0 and E = 0, where the radiation term is 3.2·10^-300.
- **Rows** (`m_relativity.rs`'s `emit_e_tail_edges`, after every older row): 108 more. Charges of 10^300, 2^-1074 and
  10^-200; R of 10^200 along an axis and diagonal, 2·10^-30 and (3, 4, 0)·10^-30; three β; β̇ of 0, 10^300 and
  10^-300. `validate.py`'s transcription agreed with every row.
- **Mappings** (the tail on the `uz_*` chain and `rnd53_sqrt`): `relativity.adr.tsv` maps 206 more rows net.
  - Row 17, 121, and row 18, 83. The corpus's older E rows at extreme scales (lw_cases' R past 10^154 and near 10^-30,
    subnormal parts) now take the unbounded tail; 88 of their outputs changed. Rows whose κ moves too keep row 17 by
    precedence.
  - Row 30, 2.

  All pass.
- **Box audit** (a sampling check in Python, not a gate): 16,364 samples inside the box, weighted to its edges, with
  κ, n̂ − β and γ² from the exact references. Rust's chain equalled the pure unbounded chain on every one.
- **Tests**: tests/relativity.tcyr gains 2 assertions, the two cases above.
- **Benchmarks** (a focused interleaved A/B against HEAD, five repetitions a round, three rounds, load 1.9):
  `relativity/lienard_wiechert_e` 672 → 711 ns, the box's nine bit-pattern tests and two compares. With §49's
  −26 ns that is about +12 ns against 2d-3d's version. The record's 662 ns becomes 710.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,832 ok, 0 failed, 0
  measured, 0 pending, 47,436 accepted (ADR 01: 8,998; ADR 0002: 38,438)** over 453,268 rows; **all 276 digests
  equal**. One changed: E.

## 51. M2 part 2d-3g: the retarded-time solver's light-cone sample at its range edges (2026-11-13)

§50 left the solver's light-cone sample: R = r_obs − r(t_r), |R|, the light-travel term c·(t − t_r),
g = |R| − c(t − t_r), dg = R̂⋅v − c and the noise 2ε·max(|R|, c(t − t_r), ‖r_obs‖∞, ‖r(t_r)‖∞, ‖v‖∞|t_r|).
`lienard_wiechert_fields`' R past DBL_MAX, which now meets row 27's refusal of an infinite R, is 2d-3h's.

- **ADR 0002 row 18, amended** (the rule of §42; a guard on Rust's own values, not a box: the sample's chain
  absorbs a degraded part rather than amplifying it):
  - **The guard.** Rust's chain is the unbounded one when:
    - dx, dy and dz are finite;
    - dx² + dy² + dz² is in the band [2^-900, DBL_MAX], or every d is 0, so |R| is its 53-bit value;
    - c·(t − t_r) is normal or 0;
    - R⋅v is in the band, or 0 from exact products (or |R| = 0, where dg is −c), so R̂⋅v − c absorbs a degraded part;
    - ‖v‖∞|t_r| is finite and the noise's maximum is at least 2^-900 (or every term is 0).

    `_bj_rel_sample` tests the common case inline: |R|² and |R⋅v| in the band, c(t − t_r) 0 or normal, and ‖v‖∞|t_r|
    finite. The noise's maximum is then at least |R| ≥ 2^-450. `_bj_rel_sample_edge` runs the full test out of line.
  - **Outside it** `_bj_rel_sample_ux` runs the sample with an unbounded exponent: R part by part, |R|, c(t − t_r),
    g (kept as its 53-bit pair and rounded once), R⋅v, dg and the noise (each rounded once).
  - **The solver** reads g's pair for:
    - g's sign (g = 0, the bracket's expansion, the lo/hi update);
    - |g_lo| < |g_hi|, |g| ≤ the noise and |2g| ≤ |last_step·dg|;
    - τ₀ = g/c, and the Newton point τ − g/dg with its step |g/dg|.

    Rust's doubles are the chain's where g is normal with 2g finite and g/dg is normal. Then |2g| ≤ |last_step·dg|
    is decided alike when the product underflows (both false: |2g| > 2^-1021) or overflows (both true). That is one
    test per iteration; `_bj_rel_step_ux` and `_bj_rel_gabs_lt` run off it.
  - **Rust's results at these edges**:
    - a static charge 2·10^308 from the observer gave t_r = −∞ (|R| = ∞, so g and τ₀ = ∞); the port gives
      −2·10^308/c = −6.67·10^299;
    - one 10^-310 away gave t_r = t, "observer on the trajectory" (|R|² underflowed to 0); the port gives
      t − |R|/c = −3.3·10^-319 at t = 0;
    - uniform motion across 2·10^308 reported a non-finite state at t = −∞; the port reports it at the finite t where
      the trajectory itself overflows (−6.67·10^299). That is the first refinement sample, not a doubling, so row 31's
      bracket error does not apply.
- **Rows** (`m_relativity.rs`'s `emit_solver_edges`, after every older row): 10 more solver rows.
  - Static charges 2·10^308 and 2·DBL_MAX from the observer, at t = 0 and 1, and a diagonal at 10^300.
  - Separations of 10^-310, 10^-160 and 5·10^-170, whose squares underflow.
  - Three uniform motions across 2·10^308.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `rel_sample_plain`, `rel_sample_ux` and `rel_solve` on the `uz_*` chain, with exact
  comparisons of g's pair): `relativity.adr.tsv` maps 9 more row-18 rows. 10^-160 apart, Rust's t_r = 1 is already
  the unbounded chain's. No older row changed. All pass.
- **Guard audit** (a sampling check in Python, not a gate): 34,946 samples kept on Rust's chain, weighted to the
  guard's edges. Rust's chain equalled the pure unbounded chain on every one.
- **Tests**: tests/relativity.tcyr gains 2 assertions, the static 2·10^308 and 10^-310 cases above.
- **Benchmarks** (a focused interleaved A/B against HEAD, five repetitions a round, three rounds, load 1.2–1.4):
  - **The first cut cost.** With a helper per pair operation it cost `relativity/solve_retarded_time_uniform_t1s`
    812 → about 1,150 ns. Both benches run 5 and 6 samples a solve and never leave the common case.
  - **The cost model.** A probe on the pinned compiler (20 million calls a construct, scratch only) put the cost at:
    - each local, integer or f64, about 0.75 ns a call;
    - each compare about 1 ns, integer and f64 alike (bit-pattern tests save nothing here);
    - a call that never runs 2–3 ns.

    HEAD's sample made four `_bj_max` calls and three `_bj_rel_max_abs3` calls (ten `_bj_max` in all), and six
    `_bj_is_finite` calls.
  - **The final shape.** The sample takes its maxima with compares (every operand is finite or +∞ and none is −0 or
    NaN, so f64::max's fold is order-free), tests the state with one sum of x − x, and keeps its error and edge paths
    in one few-argument call each. The loop keeps one test an iteration.
  - **Final**: `relativity/solve_retarded_time_uniform_t1s` 804 → 612 ns (−24 %) and
    `relativity/lienard_wiechert_fields_orbit_b099` 6.06 → 5.87 µs (−3 %). The record's 839 ns and 6.26 µs become
    612 ns and 5.87 µs.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,833 ok, 0 failed, 0
  measured, 0 pending, 47,445 accepted (ADR 01: 8,998; ADR 0002: 38,447)** over 453,278 rows; **all 276 digests
  equal**. One changed: `solve_retarded_time`.

## 52. M2 part 2d-3h: lienard_wiechert_b at its range edges (2026-11-13)

`lienard_wiechert_b` (B = n̂ × E/c) had no range-edge review: §35 left it outside row 13's scope, and the cross
product's own row-18 guard (M1 part 5d-2) rounds each component once, which is the first of B's two roundings.
`lienard_wiechert_fields`' chain is split in two bites:
- **2d-3i**: per-part exponent offsets in the exact kernels `_bj_cr_lw` and `_bj_cr_lorentz`. An R past DBL_MAX
  needs one more bit of exponent range than a double holds, and halving R would round an odd subnormal companion part.
- **2d-3j**: the fields' chain itself: R past DBL_MAX; 1/|R| and n̂ for |R| > 2^1022; β = v/c and β̇ = a/c
  subnormal; E and B from 53-bit values.

- **ADR 0002 row 18, amended** (the rule of §42; a guard on Rust's own values):
  - **The guard.** A cross component that is normal or ±0 is its 53-bit value: the cross rounds it once, correctly,
    and a ±0 carries the sign of a value too small to matter. B's product with 1/c then rounds once. That is three
    bit-pattern tests.
  - **Outside it** (a component past DBL_MAX or subnormal, for finite n̂ and E), `_bj_rel_b_edge` takes that
    component's n̂ⱼEₖ − n̂ₖEⱼ to 53 bits with an unbounded exponent (`_bj_vx_cross1`). Times 1/c, it rounds once. A
    non-finite part keeps Rust's chain.
  - **Rust's results at these edges**: n̂ = (0, 0.8, 0.6) against E = (0, −DBL_MAX, DBL_MAX) gave B_x = ∞ for
    1.4·DBL_MAX/c = 8.4·10^299. The corpus's E = DBL_MAX·(1, 1, 1) against n̂ of norm about 12 gave ±∞ components for
    about 10^300. A subnormal cross rounded twice; that differs only within about 10^-332 of a midpoint.
- **Rows** (`m_relativity.rs`'s `emit_b_edges`, after every older row): 49 more, seven n̂ against seven E:
  - crosses past DBL_MAX, and products past it;
  - subnormal and cancelling crosses;
  - a subnormal n̂ part against E near DBL_MAX.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s B on `uz_cross`): `relativity.adr.tsv` maps 31 more row-18 rows, all among the new
  rows (the other 18 are Rust's). Six older rows already on row 18 changed outputs, the E = DBL_MAX·(1, 1, 1) rows
  above. All pass.
- **Guard audit** (a sampling check in Python, not a gate): 60,000 random (n̂, E) pairs weighted to the edges.
  155,104 components were kept on Rust's product, every one equal to the pure unbounded chain; 24,896 took the edge.
- **Tests**: tests/relativity.tcyr gains 1 assertion, the 1.4·DBL_MAX cross above, with Rust's signed zeros.
- **Benchmarks** (a focused interleaved A/B against HEAD, five repetitions a round, three rounds, load 1.6–1.8):
  `relativity/lienard_wiechert_fields_orbit_b099` 5,966 → 5,956 ns and `relativity/solve_retarded_time_uniform_t1s`
  620 → 620 ns. The guard costs nothing measurable, and the record stands.
- **Toolchain**: stable rustc moved to 1.100.0 on 2026-11-13 (the corpus was built on 1.99.1). `oracle --check`
  rebuilds the corpus byte for byte under it. Its clippy flags five new lints in `m_json.rs` and `m_material.rs` under
  `-D warnings`: a hex literal's digit groups, a collapsible `if let`, a counter equal to the loop's index, and two
  `is_multiple_of`. They are fixed without changing a row.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,851 ok, 0 failed, 0
  measured, 0 pending, 47,476 accepted (ADR 01: 8,998; ADR 0002: 38,478)** over 453,327 rows; **all 276 digests
  equal**. One changed: `lienard_wiechert_b`.

## 53. M2 part 2d-3i: the exact kernels' offset forms (2026-11-13)

`lienard_wiechert_fields`' chain (2d-3j) hands E an R and a β that doubles cannot hold:
- R = r_obs − r(t_r) to 53 bits can pass DBL_MAX, up to 2·DBL_MAX. That is one bit more exponent range than a double
  has. Halving R would round an odd subnormal part beside it, so R/2 is not exact.
- β = v/c on the 53-bit grid can fall below 2^-1022, to 2^-1103.

Row 17's κ and n̂ − β (`_bj_cr_lw`) and row 30's γ² (`_bj_cr_lorentz`) must stay exact on those values. This bite gives
both kernels an offset form; 2d-3j uses it.

- **The offset form.**
  - `_bj_cr_lw(r, b, xo, …)` takes Rᵢ = r[i]·2^xo[i] and βᵢ = b[i]·2^xo[3 + i]. An offset part's double is normal,
    |Rᵢ| ≤ 2·DBL_MAX, and an offset βᵢ lies in [2^-1110, 2^-1022).
  - `_bj_cr_lorentz(v, vo, …)` takes vᵢ = v[i]·2^vo[i] for such tiny parts.
  - Every existing caller passes 0, the old form.
- **`_bj_cr_lw`, how it works.**
  - k is the largest ⌊log₂|Rᵢ|⌋, now up to 1023, so the largest |R'ᵢ| stays in [1, 4).
  - The fast path takes R' = R·2^-k and β rounded (`_bj_g_scale`). An underflowing part's absolute error stays below
    2^-1074, inside the section comment's 2^-1060.
  - The exact checks take each part's significand, and its exponent with its offset: `_bj_e_lw_parts`, `_bj_e_nb_s`,
    `_bj_e_nb_sign`, `_bj_e_nb_cand`.
  - The bignums widen: R' down to 2^-2097 and β down to 2^-1162 in their last bit take P'² about 180 bits higher, below
    2^6,740 against the buffers' 6,880 bits.
  - A zero Rᵢ's n̂ᵢ − βᵢ is now −RN(βᵢ) for a nonzero βᵢ, so an offset βᵢ that rounds to 0 keeps −βᵢ's sign. For
    doubles this is Rust's ±0 − βᵢ.
- **`_bj_cr_lorentz`.** An offset part is below every cut (2^-28·c, 2^-450), so only the slow path sees it. There 2E
  reaches −2324, D_u stays at most 2^2384 and the midpoint products below 2^2496, inside the 81 limbs (2,592 bits).
- **A port defect found on the way** (2d-1's n̂ − β, on ordinary doubles):
  - **Where.** R'ᵢ near 2^-2096 beside R near DBL_MAX, with a subnormal βᵢ of Rᵢ's sign.
  - **What went wrong.** The fast ends straddle many subnormal ulps there, and the exact path's candidate puts |βᵢ| on
    R'ᵢ's scale, |βᵢ|·2^2096, which overflowed. The candidate came out 0, and the ten-step settle loop stopped ten ulps
    above it.
  - **Example.** R = (1.5·2^1023, 2^-1074, 0) with β = (0.5, 10^-310, 0): n̂ − β's y part came out −4.9·10^-323 for
    −10^-310. With β̇ = (10^300, 0, 0), E_y was −8.7·10^-29 for −1.78·10^-16.
  - **The fix.** When |βᵢ| > 2^59·|R'ᵢ| nothing cancels, so RN(|βᵢ|) − RN(|n̂ᵢ|) is within an ulp of |dᵢ| and is the
    candidate.
  - **How it was found.** The offset families caught it first: they reach the same spread through an offset R.
    `cr_verify.py`'s `nb_families` now carries the doubles' case as a regression.
- **Verification** (`cr_verify.py`, exact references from Python's Fraction):
  - `lwo_reff`, `lwo_kappa`, `lwo_nb` (each with its `*_slow`) and `lwo_ux`, on 5,000 cases: R with one to three parts
    past DBL_MAX beside odd subnormals; `lw_cases`' pairs (their knife edges and ties) scaled past DBL_MAX; ordinary
    pairs in the offset form; β parts on the 53-bit grid down to 2^-1110; the widest spreads.
  - `gammao`, `gammao2` and their `*_slow`, on 10,000 cases (c = 1 and the speed of light).
  - The whole exact suite: **986,385 calls on x86_64 and aarch64, all equal to the exact value**.
- **Rows** (`m_relativity.rs`'s `emit_nb_spread_edges`, after every older row): 48 more E rows over three spread R,
  four β (subnormal and tiny parts of R's sign or not, and a normal control), two β̇ and two charges.
  `validate.py`'s transcription agreed with every row.
- **Mappings**: `relativity.adr.tsv` maps all 48: row 17 24 (n̂ − β) and row 18 24 (the tail). No older row changed.
  All pass. The 2d-3h sources fail 20 of them (errors past 2^57 ulps).
- **Tests**: tests/relativity.tcyr gains 1 assertion, the E_y case above (Rust: NaN).
- **Benchmarks**: (a focused interleaved A/B against HEAD over every relativity row, five repetitions a
  round, three rounds, load 1.5–2.1; medians):
  - `lienard_wiechert_e` 705 → 704 ns, `lorentz_factor` 145 → 146 and `four_vector_boost` 165 → 167;
  - `lorentz_transform_x` 188 → 192 and `lienard_wiechert_scalar_potential` 290 → 300, which run one more argument
    and one test, likely the kernels' code layout;
  - `lienard_wiechert_fields_orbit_b099` 5,922 → 5,974 and the solver 616 → 610.

  A first cut, testing the offsets part by part in `_bj_cr_lorentz`'s loop and copying β in mode 2, cost E 14 ns.
  The offset form now has its own function (`_bj_g_offset`), entered on one test. The record stands.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,851 ok, 0 failed, 0
  measured, 0 pending, 47,524 accepted (ADR 01: 8,998; ADR 0002: 38,526)** over 453,375 rows; **all 276 digests
  equal**. One changed: `lienard_wiechert_e`.

## 54. M2 part 2d-3j: lienard_wiechert_fields' chain at its range edges (2026-11-13)

The last of relativity's range edges. 2d-3i readied the exact kernels for values doubles cannot hold.
`lienard_wiechert_fields`' own chain is R = r_obs − r(t_r), |R|, 1/|R|, n̂ = R·(1/|R|), β = v/c, β̇ = a/c, then
E and B = n̂ × E/c.

- **ADR 0002 row 18, amended** (the rule of §42, a guard on Rust's own values):
  - **The guard** (`_bj_rel_fields_plain`). Rust's chain up to E is the unbounded one when:
    - R is finite (a difference of finite doubles is otherwise its 53-bit value, and a subnormal one is exact);
    - |R| ≤ 2^1022, so 1/|R| is normal;
    - each part of n̂, β and β̇ is normal, or the exact zero of a zero part of R, v or a.
  - **E's parts.** Each must be normal, or a zero of Rust's chain inside E's box (`_bj_rel_e_exact3`). Inside the box
    every intermediate is normal or an exact 0, so a zero there is exact; outside it a zero may be an underflow.
    `_bj_rel_e_impl`, E's body, reports which. The public E is now a tail call into it.
  - **Within the guard** E's doubles are its 53-bit values and n̂'s are too, so B keeps Rust's chain (with its own
    guard, §52).
  - **Outside it** (a finite charge), `_bj_rel_fields_edge`:
    - takes R, β and β̇ to 53 bits with an unbounded exponent, with an R past DBL_MAX and a β below 2^-1022 in 2d-3i's
      offset forms;
    - runs E's checks on them: |β| < 1 and γ² by `_bj_cr_lorentz`, κ < 10^-30 and n̂ − β by `_bj_cr_lw`;
    - runs E's tail on the pairs (`_bj_rel_e_tail`), which keeps E's and n̂'s 53-bit values beside E's doubles;
    - computes B from those, rounded once.
  - **A non-finite charge** keeps Rust's chain. Where R passes DBL_MAX that gives NaN fields: Rust's ∞ R makes n̂ and κ
    NaN.
  - **E's tail moved onto pairs.** The public E converts its doubles first, so its results are unchanged: every E
    digest stayed.
- **Rust's results at these edges** (q = 10^300 unless named):
  - a static charge 2·10^308 away: NaN fields, where E = 2.25·10^-307;
  - the same with an acceleration of 10^300: NaN, where E_y = −5·10^284 and B_z = −1.67·10^276;
  - q = DBL_MAX at 4000 m with β = 0.5: NaN, where E_y = ∞ (correctly rounded) and B_z = 1.9·10^302 from E's 53-bit
    value;
  - a β̇ below 2^-1022: E_x = −10^-7 (Rust: NaN);
  - a β below 2^-1022: E_x = 0, where the exact n̂ − β gives −4.4·10^-24.
- **Rows** (`m_relativity.rs`'s `emit_fields_edges`, after every older row): 13 more fields rows.
  - R past DBL_MAX: along an axis, a diagonal, with an acceleration, and beside a 2^-1074 part.
  - |R| past 2^1022, and an n̂ part subnormal under a large E.
  - β and β̇ subnormal, and E past DBL_MAX.
  - Two non-finite charges across DBL_MAX, and a control.

  `validate.py`'s transcription agreed with every row.
- **Mappings** (`validate.py`'s `rel_fields_plain` and `rel_fields_edge`, the latter on exact values and the `uz_*`
  chain; E's model now `rel_lw_e`, returning its box): `relativity.adr.tsv` maps 10 more row-18 rows, all passing.
  The other 3 are Rust's: the two non-finite charges, and the subnormal n̂ part, whose effect stays below B's last
  bit. No older row changed.
- **Guard audit** (a sampling check in Python, not a gate): 8,648 cases kept on Rust's chain over two seeds, weighted
  to the edges. Every one equalled the pure unbounded chain; 27,027 took the edge.
- **Tests**: tests/relativity.tcyr gains 2 assertions, the static 2·10^308 case and the DBL_MAX charge above.
- **Benchmarks** (a focused interleaved A/B against HEAD over every relativity row, three rounds, load 1.0–1.7):
  - `relativity/lienard_wiechert_e` 707 → 711 ns. A first cut that unpacked and re-wrapped the inner Result cost
    14 ns.
  - `relativity/lienard_wiechert_fields_orbit_b099` 5,968–5,972 → 6,028–6,075 ns (+1–2 %): the guard, E's exactness
    test and the box flag. The record's 5.87 µs becomes 6.05.
  - The other rows within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **13 modules, 0 failing; 405,854 ok, 0 failed, 0
  measured, 0 pending, 47,534 accepted (ADR 01: 8,998; ADR 0002: 38,536)** over 453,388 rows; **all 276 digests
  equal**. One changed: `lienard_wiechert_fields`.
- **Relativity's range-edge review is complete**, every public function and method checked (`electric_field` is one
  product, `magnetic_field` copies). With circuit's 2d-3a and 2d-3b, **2d-3 is done, so 2d is done: Part 2 is complete,
  and with it M2.**

## 55. M3 part 1a-1: polarization (2026-11-13)

- **`src/polarization.cyr`**: polarization.rs's Jones vectors and matrices, Stokes vectors and Mueller matrices, with
  all their constructors and methods under Rust's names.
  - `BjJonesVector` holds two BjComplex. `BjJonesMatrix` holds Rust's `m: [[Complex; 2]; 2]` as four fields `m00` to
    `m11`: a tuple's slots are 8 bytes, so it can't hold a BjComplex.
  - `BjStokesVector` holds `s` as a 4-tuple. `BjMuellerMatrix` holds `m` as a row-major 16-tuple, as BjEmTensor does.
  - The Default and Add traits are carried.
  - Every formula is Rust's, in Rust's order, and complex arithmetic is BjComplex's (row 14 at a range edge).
  - Mueller's row sums start from −0.0, as Rust's `f64` `Sum` does (checked on rustc 1.99.1).
  - The four structs' JSON (their serde derives) is part 1a-2's. `test_mie_result_serde` waits for scattering
    (part 1b).
- **The one convention in polarization (ADR 0002 row 5, D3, amended)**:
  - **The elements.** The port works in e^{+jωt} with IEEE handedness, so for one physical state or element its
    complex amplitudes are Rust's conjugated. Only the retarders differ in their numbers:
    - `BjJonesMatrix_waveplate` is Rust's matrix conjugated: the fast axis takes e^{+jδ/2}, the slow e^{−jδ/2}.
    - `BjMuellerMatrix_waveplate` is D·M·D, D = diag(1, 1, 1, −1): its four S₃-coupling entries are the negatives
      of Rust's.

    Each changed entry is computed as the exact negative of Rust's expression.
  - **The names.** `right_circular` is IEEE right-hand circular, (1, −j)/√2, turning from x̂ toward ŷ about +z.
    It is Rust's `left_circular` (Born & Wolf's "left") conjugated, which has the numbers of Rust's
    `right_circular`. So the circular Jones states, the circular Mueller polarizers and the circular Stokes states
    carry Rust's numbers under the IEEE names.
  - **S₃.** S₃ = 2 Im(E_x E_y*) keeps its meaning, +1 for right-hand (IEEE). For one physical state it is Rust's S₃
    negated, and so is the ellipticity angle.
  - **How the harness compares.**
    - Every function but the retarders is the same formula in both conventions, so its rows compare on the
      oracle's numbers directly, bit for bit. The named states compare by name, since their numbers coincide.
    - The 4,104 retarder rows compare exactly: as port = conj(oracle) for Jones and D·oracle·D for Mueller.
    - The D·M·D mapping negates only the entries a retarder computes. Its m[0][3] and m[3][0] are zeros in both
      conventions; a first form negated them too, which mapped every Mueller retarder row for a zero's sign.
    - A first form also conjugated every amplitude both ways. That is exact on values but not on IEEE's zero signs:
      conj(1 + 0j) is 1 − 0j while the port's horizontal is 1 + 0j, and (+0) + (−0) = +0. It would have mapped most
      Jones rows for a zero's sign alone.
  - **Pinned by tests** (`tests/polarization.tcyr`):
    - the rotation-sense group samples Re{J·e^{+jωt}} at two instants for each named state and for a QWP's output,
      and checks the turn against +z;
    - each circular polarizer passes its own state and blocks the other;
    - S₃ and χ have the IEEE signs;
    - the Jones and Mueller waveplates agree on nine (θ, δ) pairs.
- **New golden file `polarization.tsv`, 23,135 rows over 53 functions.** Class C, 4,779 rows, is the cos, sin, atan2
  and hypot chains; the rest is class A. The grids:
  - complex amplitudes (the specials' pairs, the unit circle, seeded values, and the scales where |z|² under- and
    overflows), paired into Jones vectors;
  - the named elements at four angles, plus seeded and range-edge matrices;
  - 76 angles and 25 retardances;
  - Stokes vectors from the named states, the specials and seeded values, at the range edges;
  - Mueller matrices from the named elements and seeded entries.

  `validate.py`'s independent transcription (`pol_eval`) agreed with every row on its first run. Glibc's cos, sin,
  atan2 and hypot equal MPFR's on every class-C row, so there is no `polarization.cr.tsv`. The class-C rows gate at
  0 ulp.
- **Mappings**: `polarization.adr.tsv` has 61 rows, all row 14 (BjComplex's products and |z|² at a range edge):
  `from_jones` 40, `inner` 9, `apply` 6, `compose` 4, `determinant` 2.
- **Rust tests**: 68 of polarization.rs's 69 are ported (`tests/polarization.tcyr`, 166 assertions with the
  Cyrius-only groups). `test_waveplate_slow_axis_lags` reads the same light under the IEEE name: phase −π/2, the
  port's `right_circular`, S₃ = +1. The ledger: 429 ported, 3 n/a, 426 pending of 858.
- **Surface parity**: polarization is in progress, 66 items matched (323 in all), with the eight serde trait items
  pending for part 1a-2. `surface_parity.py`'s in-progress list now takes owner-qualified trait keys. `ComplexExt`
  is an exception: its methods are BjComplex's own `norm_sq`, `norm` and `arg`. The API snapshot gains 40 public fns.
- **Benchmarks** (the four polarization twins, so 50 of 80 are twinned; load 1.1–1.3, medians of 30):

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `polarization/jones_apply_polarizer` | 0.99 ns | 213 ns | 215× |
  | `polarization/jones_qwp_apply` | 0.99 ns | 216 ns | 218× |
  | `polarization/mueller_apply` | 0.99 ns | 86 ns | 87× |
  | `polarization/stokes_from_jones` | 5.44 ns | 117 ns | 21.5× |

  Rust's arguments are constants that LLVM folds. The port pays the by-value structs (a Jones matrix is 64 bytes) and
  a call per complex product (`BjComplex_mul` and its range test). A finding: Jones `apply`'s four products are most
  of its 213 ns.
- **For M3's range-edge Part** (§0.1), candidates from writing the port:
  - `from_jones`'s S₀ = |E_x|² + |E_y|², a sum of two rounded norms;
  - the ellipticity angle's √(S₁² + S₂²), where S₁² overflows while χ is representable;
  - the degree of polarization's two hypot roundings;
  - the waveplate's c² − s² and 2cs chains near a fast axis at π/4 or 0;
  - the products past DBL_MAX in Jones `apply`, `compose` and `determinant` (row 14 rounds each product once, not
    the sums).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **14 modules, 0 failing; 428,928 ok, 0 failed, 0
  measured, 0 pending, 47,595 accepted (ADR 01: 8,998; ADR 0002: 38,597)** over 476,523 rows; **all 329 digests
  equal**, the 276 of M1 and M2 unchanged.

## 56. M3 part 1a-2: polarization's JSON (2026-11-13)

- **`src/json.cyr`**: `BjJonesVector`, `BjJonesMatrix`, `BjStokesVector` and `BjMuellerMatrix` get `_to_json` and
  `_from_json` (D5, amended), on the shared writer and reader. serde's shapes:
  - a Jones vector is `{"x":<Complex>,"y":<Complex>}`;
  - a Jones matrix is `{"m":[[<Complex>,<Complex>],[<Complex>,<Complex>]]}`, Rust's `m[r][c]` being the port's field
    m<r><c>;
  - a Stokes vector is `{"s":[S₀,S₁,S₂,S₃]}`;
  - a Mueller matrix is `{"m":[[..],[..],[..],[..]]}`, row by row.

  All four lie in memory in serde's field order. Each writer copies its struct to an f64 array and writes from it
  (buffers of 160, 320, 128 and 512 bytes, sized for the longest numbers bayan prints). Each reader decodes into an
  array and copies it to `out` only on success.
- **The reader** applies serde's semantics at every level. An inner struct (a Complex) reads from an object or by
  position. A Rust `[T; n]` reads from an array of exactly n T: `_bj_jr_f64x3` became `_bj_jr_array` and `_bj_jr_f64s`,
  and PointCharge's arrays read through them. The whole-document checks (ADR 0002 row 12) are unchanged.
- **The convention (ADR 0002 row 5, amended)**: a document holds the stored numbers, and the conjugation and the
  renaming cancel on them. So a document bijli 1.x wrote keeps its handedness names: a named state, a polarizer and
  a Stokes vector read back under the same names. A retarder's numbers are odd in δ where they differ, so a matrix
  1.x wrote as waveplate(θ, δ) reads as 2.0's waveplate(θ, −δ). This holds in value on all 584 (θ, ±δ) waveplate
  pairs of `polarization.tsv`, Jones and Mueller. The harness compares the JSON rows on the oracle's numbers
  directly: no field changes its meaning in a document.
- **`json.tsv` grows by 22,199 rows to 41,688.** The new rows:
  - writers: 1,186 per type;
  - readers: 3,798 JonesVector, 4,436 JonesMatrix, 3,695 StokesVector and 5,526 MuellerMatrix.

  The generator appends the four kinds after Material, drawing a set as wide as the type (16 for a Mueller matrix).
  The 19,489 earlier rows are an exact prefix of the new file, and their 329 digests are unchanged.
  `nested` adds each inner Complex and array in every form serde reads or refuses:
  - by position, reordered;
  - a field missing, repeated or unknown, an escaped key;
  - a lone surrogate in a key or an unknown string;
  - short, long and wrong-typed arrays, an object in an array's place;
  - mixed forms at each level.

  `validate.py`'s model became a type grammar (`f64`, a struct, `("arr", n, T)`) and agreed with every new row on its
  first run. `json.cr.tsv` grows by 12,157 rows, to 21,056: the readings serde_json's best-effort parse gets wrong,
  gated at 0 ulp. 82 more rows are ADR 0002 row 12's (154 in all).
- **Tests**: `cyrius_polarization_json` in `tests/polarization.tcyr` (19 assertions, 185 in the file):
  - serde's exact text for one named element of each type;
  - bit-for-bit round trips of a QWP and a waveplate, and of the longest numbers in every slot;
  - both of serde's forms at every level, and any numbers read as written;
  - the reasons (`invalid length 3, expected an array of length 4`, `missing field `im``), with out untouched.

  `test_mie_result_serde` stays pending for scattering (part 1b). The ledger is unchanged: 429 ported, 3 n/a, 426
  pending.
- **Surface parity**: polarization is a ported module, 10 in all; 331 items matched, 0 in progress. The API snapshot
  gains the 8 functions (2,002 public).
- **Benchmarks** (Cyrius-only: rust-old benches no serde), interleaved against HEAD, three rounds at load 1.1–1.3:
  - `json/mueller_matrix_to_json`: 9.21 µs;
  - `json/mueller_matrix_from_json`: 8.78 µs (about 4.6 µs over a FieldVector's 1.88 µs, for 16 numbers in five
    containers against 3 in one).

  The shared-reader refactor costs nothing measurable: `json/point_charge_from_json` is 4.48 µs before and 4.47 µs
  after (medians). The one row that moved, `relativity/em_tensor_from_fields`, is bimodal at 31 and 34 ns in both
  builds.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **14 modules, 0 failing; 438,934 ok, 0 failed, 0
  measured, 0 pending, 59,788 accepted (ADR 01: 21,109; ADR 0002: 38,679)** over 498,722 rows; **all 337 digests
  equal**, the 329 before unchanged. `polarization.tcyr`, `json.tcyr` and `material.tcyr` pass on all three.

## 57. M3 part 1b-1: scattering (2026-11-13)

- **`src/scattering.cyr`**: scattering.rs's Mie and Rayleigh functions under Rust's names (`bj_mie`,
  `bj_cross_section_from_efficiency`, `bj_rayleigh_*`), with the limits `BJ_MIE_MAX_SIZE_PARAMETER` and
  `BJ_MIE_MAX_MX` (already in the generated `src/constants.cyr`).
  - **`BjMieResult`** holds Q_sca, Q_ext, Q_abs and g, and owns `a` and `b` (Rust's `Vec<Complex>`) as fl_alloc
    buffers with their lengths (`a_len`, `b_len`; equal from `bj_mie`, apart only for a document part 1b-2 reads).
    `a_at(k)` and `b_at(k)` read them. `bj_mie` writes new buffers into `out` and never reads it, so a stale `out`
    can't be freed by mistake. `BjMieResult_free` gives them back. On an error, `out` is untouched and the call's
    buffers are freed.
  - **The convention**: scattering.rs already used e^{+jωt} with absorbing Im(m) < 0, so ADR 0002 row 5 changes no
    number here and the harness compares every row on the oracle's numbers. The gain-index reason keeps Rust's
    words ("e^{+iωt}", "n − iκ"): error rows carry the Display text (ADR 0002 row 4).
  - **Mie's recurrences** (the downward Dₙ(mx), the upward ψₙ and χₙ, the coefficient loop and its three sums) run
    on f64 locals, as probe p8 found (a struct-valued step is 2.2× slower). They follow Rust's order of operations,
    with hisab's complex products and sums written out. The module's own `complex_recip` (one division when |z|² is
    normal) and Smith `complex_div` are ported as written, the zero test at each call. A `+ 0.0` that Rust's
    `complex_real` adds to an imaginary part is kept (it turns −0 into +0); the `− 0.0` it subtracts is the identity
    and is dropped. sin and cos are src/crmath.cyr's, hypot is bijli's, and cbrt is ganita's, correctly rounded as
    Rust's is (CORE-MATH's, port plan P3). It is the bundle's first direct ganita call, so `dist/bijli.deps` now
    lists ganita, a dependency `cyrius.cyml` already declares.
  - **Rayleigh** is cold, so it uses BjComplex: `BjComplex_mul` for m·m and `norm_sq` for |K|², with row 14 at a
    range edge. 138 rows reach the m² edge (m = 10²⁰⁰ + 10²⁰⁰j, 10¹⁶⁰, 10⁻²⁰⁰ + 10⁻²⁰⁰j and the like), and every
    one ends in the same NaN under both chains, so no mapping is needed.
- **New golden file `scattering.tsv`, 3,907 rows.** `mie` (1,027 rows) is class C; the rest is class A. A `mie`
  row's outputs are Q_sca, Q_ext, Q_abs, g and the count n, then (aₖ, bₖ) at every k when n ≤ 128. Above that, the
  row takes 128 indices k_j = j·(n − 1)/127 and the FNV-1a fold of all 4n words, which keeps the 100,188-term row
  small and still pins every coefficient bit. The grids:
  - **Mie's cross**: 26 sizes (the ψ₁ threshold and its neighbours, the perfect cubes where Wiscombe's count is an
    integer, Bohren & Huffman's sample, the Rust tests' cases) × 26 indices (real low and high contrast, matched,
    below 1, negative, absorbing, metal-like, pure imaginary, signed zeros);
  - **Mie's limits**:
    - large spheres to x = 10⁵;
    - |m|·x = 10⁷ accepted (the D_n start at about 10⁷ + 877) and its next double refused;
    - tiny x, where the series overflows somewhere between x = 10⁻⁵⁰ (computed) and 10⁻¹⁰⁰ (refused);
    - tiny m, and mx underflowing to zero (DivisionByZero);
    - every rejection, and 240 seeded spheres;
  - **Rayleigh**: radii, wavelengths and indices with the specials and the range edges of m², the 10⁻³⁰ singularity
    test on both sides, and the phase function.

  `validate.py`'s independent transcription (`sc_eval`) agreed with every row on its first run. One row is past its
  2·10⁶-step limit (|m|·x = 10⁷); validate checks it through its primitives, whose agreement with MPFR's makes the
  oracle's reading the corrected one. `scattering.cr.tsv` holds 6 rows, every x = 200 sphere: **glibc 2.44's
  cos(200) misrounds**, and the port gives the correctly rounded readings (0 ulp). No ADR 0002 mapping is needed.
- **Rust tests**: all 32 of scattering.rs's are ported (`tests/scattering.tcyr`, 90 assertions with the Cyrius-only
  groups). The Cyrius-only groups cover:
  - every reason's exact text in Rust's Display forms;
  - out untouched by each refusal;
  - 2,000 calls leaving the freelist arena and the bump heap as they were;
  - |m|·x = 10⁷ and x = 10⁵ accepted.

  The ledger: 461 ported, 3 n/a, 394 pending of 858.
- **Surface parity**: scattering is in progress, 15 items matched (346 in all); MieResult's two serde items are part
  1b-2's. The API snapshot gains 8 public functions (2,010).
- **Benchmarks** (the 5 scattering twins: 55 of 80 twinned). The recorded run was at load 1.06; an interleaved A/B
  against HEAD, three rounds at load 1.5–1.9 beside another repo's `cr_verify.py`, gave medians within 2 % of it,
  and left every other row within noise:

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `scattering/mie_x1_real` | 423.87 ns | 1.81 µs | 4.3× |
  | `scattering/mie_x10_real` | 859.83 ns | 3.62 µs | 4.2× |
  | `scattering/mie_x50_m3_5` | 3.44 µs | 11.74 µs | 3.4× |
  | `scattering/rayleigh_cross_section` | 1.74 ns | 129 ns | 74× |
  | `scattering/rayleigh_phase_function` | 0.50 ns | 4 ns | 8× |

  Mie is 3.4–4.3× Rust's: two Smith divisions a term out of line, and three buffers allocated and freed a call, as
  Rust's Vecs are. Rayleigh's constant arguments fold in Rust; the port pays BjComplex_mul and norm_sq with their
  range tests, the Smith division and the Result.
- **For M3's range-edge Part** (§0.1), candidates from writing the port:
  - Mie's complex products in the recurrences (hisab's plain chain on f64 locals), reachable with a tiny or huge m;
  - the tiny-x refusals: between x = 10⁻⁵⁰ and 10⁻¹⁰⁰, χₙ ≈ (2n − 1)!!/xⁿ overflows, so the series is refused
    although Q is a representable (or zero) value;
  - |m|·x > 10⁷ decided on the rounded product hypot(m)·x (row 17's exact-threshold question);
  - the products of many factors in Rayleigh's σ and Q, `cross_section_from_efficiency` and
    `rayleigh_angular_intensity` (row 13);
  - x⁴ and m² past either end of the range in Rayleigh. For |m| past about 10¹⁵⁴, m² overflows and K = (m² − 1)/(m² + 2)
    is NaN, though K → 1 and σ is finite; row 14's unbounded m² still ends in ∞/∞.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **15 modules, 0 failing; 442,835 ok, 0 failed, 0
  measured, 0 pending, 59,794 accepted (ADR 01: 21,115; ADR 0002: 38,679)** over 502,629 rows; **all 343 digests
  equal**, the 337 before unchanged. `scattering.tcyr` passes on all three.

## 58. M3 part 1b-2: MieResult's JSON (2026-11-13)

- **`src/json.cyr`**: `BjMieResult_to_json` and `BjMieResult_from_json` (D5, amended), on the shared writer and
  reader. serde's shape is `{"q_sca":..,"q_ext":..,"q_abs":..,"g":..,"a":[<Complex>,..],"b":[<Complex>,..]}`.
  - **The writer** builds the document in one fl_alloc block: 160 bytes and 64 a coefficient, since an f64 is at
    most 25 bytes and a Complex object with its comma at most 64. It appends the document in one grow, so a refused
    allocation leaves the builder as it was.
  - **The reader** takes `a` and `b` as arrays of any length (`_bj_jr_cvec`: serde's `Vec`, "invalid type: map,
    expected a sequence" for anything else). Each element is a Complex in either of serde's forms, so the two
    lengths may differ, as Rust's public fields may. Each vector's block is stored before its elements are read: a
    read that fails later gives back what it took, and `out` is untouched. A successful read writes new blocks into
    `out` and never reads it, as `bj_mie` does.
- **`json.tsv` grows by 5,298 rows, to 46,986**: 1,029 writer rows and 4,269 reader rows. `m_json.rs` emits them
  after the other types, so the 41,688 earlier rows are an exact prefix, and their 343 digests are unchanged. A row's
  value is its four f64, then `i:` a's length and parts, then `i:` b's length and parts. The values:
  - 21 real Mie results, to x = 200;
  - the corpus's values with every pair of lengths up to 7, and 120 seeded bit patterns;
  - two 500-element vectors.

  The documents:
  - serde's text, bayan's notation, Rust's `{}` and `{:e}`;
  - every spelling in each number position of a two-and-two document, and 600 seeded decimals;
  - the shapes: both forms at every level, unknown keys of every kind, each field missing, repeated or of the wrong
    type, a `Vec` as an object or holding non-Complex elements, every inner Complex form, vectors of 3 to 300
    elements, a deep unknown value inside an element at serde's nesting limit, every prefix, trailing data.

  `validate.py`'s JSON grammar gained `("vec", T)` (a length, then the elements) and agreed with every row on its
  first run. `json.cr.tsv` grows by 3,060 rows to 24,116: the readings serde_json's best-effort parse gets wrong. 15
  more rows are ADR 0002 row 12's (169 in all).
- **Tests**:
  - **`test_mie_result_serde`**, polarization.rs's last test, is ported (`tests/polarization.tcyr`, 186 assertions),
    so polarization.rs's 69 are all ported.
  - **`cyrius_mie_result_json`** (`tests/scattering.tcyr`, 102 assertions) covers:
    - serde's exact text, an empty `b` included;
    - a real result's bit-for-bit round trip, all 42 coefficients;
    - both forms at every level, and unequal lengths;
    - the reasons, with out untouched;
    - 2,000 reads, half of them refused after a vector was read, leaving the freelist arena and the bump heap as
      they were.

  The ledger: 462 ported, 3 n/a, 393 pending of 858.
- **Surface parity**: scattering is a ported module (11 in all, 348 items matched, 0 in progress). The API snapshot
  gains the 2 functions (2,012 public).
- **Benchmarks** (Cyrius-only; rust-old benches no serde), on a Mie result with 21 + 21 coefficients, interleaved
  against HEAD over three rounds at load 1.3–1.7:
  - `json/mie_result_to_json`: 86.5 µs;
  - `json/mie_result_from_json`: 86.1 µs, its blocks freed each time through.

  That is 88 numbers at about 1 µs each, in line with FieldVector's 0.77 µs a number. Every other row is within
  noise (±3.7 %).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **15 modules, 0 failing; 445,069 ok, 0 failed, 0
  measured, 0 pending, 62,858 accepted (ADR 01: 24,164; ADR 0002: 38,694)** over 507,927 rows; **all 345 digests
  equal**, the 343 before unchanged. `scattering.tcyr`, `polarization.tcyr` and `json.tcyr` pass on all three.

## 59. M3 part 1c-1: beam (2026-11-13)

- **`src/beam.cyr`**: beam.rs under Rust's names.
  - `BjGaussianBeam` (wavelength, waist, Rayleigh range) and its methods; `q` returns (Re q, Im q).
  - `BjAbcdMatrix`: its constructors, the refusing ones into `out`; `apply_ray`, `compose`, `determinant`,
    `propagate_beam` (into `out[0..2]`) and the two beam transforms (into `out` and the waist distance).
  - `BjResonatorStability` (`is_stable` 1 or 0), `bj_resonator_stability`, `bj_resonator_beam_waist` and
    `bj_transverse_mode_spacing`, whose unstable-cavity reasons keep Rust's `{:.4}` (`_bj_fmt_f64_fixed`).
  - `bj_hermite_gaussian`, `bj_hermite_gaussian_field` and `bj_laguerre_gaussian_intensity`.
  - The private polynomials are `public fn _bj_hermite_polynomial` and `_bj_laguerre_polynomial`, as
    `_bj_associated_legendre` is, because their Rust tests reach them.

  The mode indices are Rust's u32 and i32, as documented, and every formula runs in Rust's order of operations.
- **The convention (ADR 0002 row 5, amended)**:
  - **The field**: only `hermite_gaussian_field` carries a phase. The port's is Rust's conjugated part for part,
    (u·cos φ, −(u·sin φ)) from Rust's own phase, a field ∝ e^{−jkz} as Siegman writes it.
  - **q keeps its numbers**: q = z + j·z_R is already Siegman's e^{+jωt} form (1/q = 1/R − jλ/(πw²)). Rust's ABCD
    law used it beside an e^{−iωt} field, so the port makes the two agree with no change.
  - **Everything else is convention-free**: the Gouy phase ζ = arctan(z/z_R) is a geometric angle, and every other
    value is real.
  - **The harness** compares the 3,325 field rows as port = conj(oracle), exactly, and every other row on the
    oracle's numbers.
  - **Pinned by tests**: `cyrius_beam_convention` checks the phase falling along +z (arg E = −(kz − ζ)).
    `test_hermite_gaussian_field_closed_form` reads its imaginary parts with the opposite sign, and says so.
  - **One sincos**: Rust's `Complex::from_polar` compiles to one sincos, and the port calls `_bj_sin_cos` once (the
    correctly rounded bits of separate calls; the field's digest is unchanged by the switch). That made the field
    bench 669 ns, from 789 with `BjComplex_from_polar`'s two reductions.
- **ADR 0002 row 32, new**: `laguerre_gaussian_intensity` raised 2r²/w² to `l.unsigned_abs() as i32`. For
  l = i32::MIN the cast wraps |l| = 2³¹ to −2³¹, so Rust returned the reciprocal power (∞ where the true power
  underflows to 0). The port raises it to 2³¹ through the same `__powidf2` loop. `beam.adr.tsv` maps 366 rows (every
  LG row at l = i32::MIN that the wrap changes), all passing. `cyrius_lg_i32_min` pins 0 and ∞.
- **New golden file `beam.tsv`, 40,086 rows over 29 functions.** Class C, 13,516 rows, is `gouy_phase`,
  `transverse_mode_spacing` and the three mode functions; the rest is class A. A beam row holds the beam's three
  public fields, so the methods also see fields `new` would refuse (zero, negative, NaN, ∞). The grids:
  - laser wavelengths and waists, and the range edges where z_R under- or overflows;
  - distances at the waist, the 10⁻³⁰ test on both sides, ±z_R, the far field and the specials;
  - mode orders to 10 at nine transverse points and five planes, orders 50–200 where H overflows, and l in ±1, ±2,
    ±5, 3, 10, i32::MIN and i32::MAX;
  - every ABCD constructor over its grid, 2f imaging, seeded and special matrices, the 10⁻⁶⁰ denominator test, and
    the |det − 1| ≤ 10⁻¹² wavelength rule on both sides;
  - cavities across the stability diagram's interior, edges and outside, with flat and zero-radius mirrors.

  `validate.py`'s independent transcription (`bm_eval`) agreed with every row on its first run, and a one-bit
  perturbation of every row's first output was caught. Glibc's atan, acos, exp, cos and sin equal MPFR's on every
  class-C row, so there is no `beam.cr.tsv`.
- **Rust tests**: all 60 of beam.rs's are ported (`tests/beam.tcyr`, 90 assertions with the Cyrius-only groups:
  the convention, the reasons' exact text, out untouched by a refusal, row 32). The ledger: 522 ported, 3 n/a, 333
  pending of 858.
- **Surface parity**: beam is in progress, 43 items matched (391 in all), with the three structs' six serde items
  pending for part 1c-2. The API snapshot gains 15 public functions (2,027).
- **Benchmarks** (the 6 beam twins, so 61 of 80 are twinned). The recorded run was at load 1.13, and an interleaved
  A/B against HEAD (three rounds, load 1.4–1.9) agreed within 1 ns and left every other row within noise (±3.6 %):

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `beam/spot_size` | 0.49 ns | 10 ns | 20× |
  | `beam/abcd_propagate_beam` | 1.79 ns | 19 ns | 11× |
  | `beam/abcd_compose` | 0.99 ns | 20 ns | 20× |
  | `beam/hermite_gaussian_hg33` | 4.46 ns | 160 ns | 36× |
  | `beam/resonator_stability` | 6.48 ns | 27 ns | 4.2× |
  | `beam/hermite_gaussian_field_hg33` | 30.33 ns | 669 ns | 22× |

  Rust's constant arguments fold, the HG ones nearly whole. The port pays the calls (exp, atan, sincos, the
  recurrences) and by-value matrices. A finding: the field costs about 510 ns over the envelope's 160 ns, in the
  Gouy atan and the sincos of a phase near 5.9·10⁴ rad (kz at z = 1 cm); the split between the two is not measured.
- **For M3's range-edge Part** (§0.1), candidates from writing the port:
  - `propagate_beam`'s |C·q + D|² < 10⁻⁶⁰, on a rounded sum of squares (row 17's exact-threshold question);
  - `transform_beam`'s |det − 1| ≤ 10⁻¹² wavelength rule, a tolerance Rust chose (it keeps λ for a det just off 1
    that two real media could give);
  - `resonator_beam_waist`'s exact s = 0 test and its g₁g₂(1 − g₁g₂) near the stability edges;
  - the products of many factors (z_R = πw₀²/λ, the envelope's scale, the LG power and its prefactors) past either
    end of the range (rows 13 and 18);
  - the phase kz + k(x² + y²)/(2R) − (m + n + 1)ζ, which loses its digits where kz is large.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 484,789 ok, 0 failed, 0
  measured, 0 pending, 63,224 accepted (ADR 01: 24,164; ADR 0002: 39,060)** over 548,013 rows; **all 374 digests
  equal**, the 345 before unchanged. `beam.tcyr` passes on all three.

## 60. M3 part 1c-2: beam's JSON (2026-11-14)

- **`src/json.cyr`**: `BjGaussianBeam`, `BjAbcdMatrix` and `BjResonatorStability` get `_to_json` and `_from_json`
  (D5, amended), on the shared writer and reader. serde's shapes:
  - a beam is `{"wavelength":..,"waist":..,"rayleigh_range":..}`;
  - an ABCD matrix is `{"a":..,"b":..,"c":..,"d":..}`;
  - a stability is `{"g1":..,"g2":..,"g_product":..,"is_stable":true}`.

  Each writer builds its document in a 128-byte buffer: at most 122 bytes with the longest numbers bayan prints
  (25 bytes, `-0.0000012345678901234567`). Each reader decodes into words and copies them to `out` only on success.
- **The bool**: `is_stable` is the port's first JSON bool.
  - The writer writes `true` for any value but 0, as `bj_resonator_beam_waist` reads it, and `false` for 0.
  - The reader (`_bj_jr_bool`) takes only JSON's two literals, as serde's `deserialize_bool` does: `true` is 1,
    `false` is 0. The oracle's rows show serde_json refusing everything else in that position: `1`, `0`, `-0`,
    `1.0`, `"true"`, `null`, `[true]`, `{}`, and every near miss of the literals (`True`, `tru`, `true1`, `trUe`).
    The port refuses each, with serde's wording (`invalid type: integer, expected a boolean`).
- **The convention**: the three structs hold real numbers (q is computed, not stored), so ADR 0002 row 5 changes
  nothing in their documents, and the harness compares every row on the oracle's numbers. No ADR 0002 row is
  needed: the readers and writers are serde's, bit for bit.
- **Read as written**: serde's derive calls no constructor, so a beam `BjGaussianBeam_new` would refuse (a zero
  waist), a z_R that isn't πw₀²/λ, and an `is_stable` that disagrees with g₁g₂ are all read as the document holds
  them, as Rust's public fields are.
- **`json.tsv` grows by 14,416 rows, to 61,402**: 1,186 writer rows per type, and 3,549 GaussianBeam, 3,702
  AbcdMatrix and 3,607 ResonatorStability reader rows.
  - `m_json.rs`'s per-type loop became `emit_kind`, and the three new `Kind`s (`Kind::BEAM`) are emitted after
    MieResult. So the 46,986 earlier rows are an exact prefix of the new file, and their 374 digests are unchanged.
  - A `Kind` value may now be a bool (`is_bool`): 0.0 or 1.0 in a value, `i:0` or `i:1` in a row (as `beam.tsv`'s
    `resonator_stability` rows), `false` or `true` in a document. It alternates over the value sets, and its position
    takes its own 48 spellings (`bool_spellings`); the number spellings and the seeded decimals go to the f64
    positions only.

  `validate.py`'s JSON grammar gained `"bool"` (`JsonBool`, an `i:` token) and agreed with every new row on its first
  run. Flipping each new row's first output (a bit of its first f64, or its outcome) was caught on all 14,416 rows,
  and flipping `is_stable` on all 4,380 rows that carry it. `json.cr.tsv` grows by 5,295 rows, to 29,411: the
  readings serde_json's best-effort parse gets wrong, gated at 0 ulp. 54 more rows are ADR 0002 row 12's (223 in
  all).
- **Tests**: `cyrius_beam_json` in `tests/beam.tcyr` (23 assertions, 113 in the file) covers:
  - serde's exact text for a beam (z_R = π), a thin lens and a stable and an unstable cavity, and `true` for any
    `is_stable` but 0;
  - bit-for-bit round trips with the longest number in every slot;
  - both of serde's forms, fields in any order, and values read as written;
  - the bool's refusals (`1`, `"true"`, `null`, `True`), the reasons, and out untouched;
  - 2,000 reads, half of them refused, leaving the freelist arena and the bump heap as they were.

  beam.rs has no serde test of its own, so the ledger is unchanged: 522 ported, 3 n/a, 333 pending of 858.
- **Surface parity**: beam is a ported module (12 in all), 397 items matched, 0 in progress. The API snapshot gains
  the 6 functions (2,033 public).
- **Benchmarks** (Cyrius-only; rust-old benches no serde), medians of five rounds interleaved against HEAD at load
  0.8–1.9:
  - `json/resonator_stability_to_json`: 2.28 µs;
  - `json/resonator_stability_from_json`: 2.45 µs, a FieldVector's 1.86 µs and one more field.

  Every other row is within noise (±2.9 %) but two readers whose code is unchanged: `json/point_charge_from_json`
  4.56 → 4.74 µs (+4.0 %) and `json/field_vector_from_json` 1.857 → 1.905 µs (+2.6 %). It is code layout. A variant
  with the ten new functions appended at `src/json.cyr`'s end, so that no existing function moves, read 4.40–4.59 and
  1.86–1.87 µs, HEAD's numbers, in three rounds interleaved with this build's 4.70–4.77 and 1.90–1.92 µs; the
  Mueller and Mie readers were equal in both. The file keeps its sections (writers, decoders, readers): any later
  edit would move the readers again. `material/curie_weiss` sits at 7 or 8 ns in both builds.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,892 ok, 0 failed, 0
  measured, 0 pending, 68,537 accepted (ADR 01: 29,423; ADR 0002: 39,114)** over 562,429 rows; **all 380 digests
  equal**, the 374 before unchanged. `beam.tcyr` and `json.tcyr` pass on all three. M3 Part 1 is done: polarization.rs,
  scattering.rs and beam.rs are ported.

## 61. M3 part 2a-1: the waveplates cancellation-free (2026-11-14)

- **The user's answers of 2026-11-14** open M3's range-edge Part (plan §0.1):
  - the waveplates are cancellation-free (this bite);
  - the resonators decide stability on the exact (1 − L/R₁)(1 − L/R₂) and compute w₀ without cancellation (2a-2);
  - `transform_beam` decides on the exact determinant (2a-3);
  - Mie keeps its documented limit below x ≈ 10⁻⁵⁰.

  The Part runs as 2a–2d (plan §0.1). The HG field's phase keeps Rust's chain: rounding kz costs about ulp(kz),
  which the input z's own rounding already costs.
- **ADR 0002 row 33, new**:
  - `BjJonesMatrix_waveplate`, and so the quarter- and half-wave plates, takes cos 2θ and sin 2θ, each correctly
    rounded (`_bj_pol_cos_sin_2t`): one `_bj_sin_cos` of the exact double 2θ. Where |θ| ≥ 2^1023, 2θ passes
    DBL_MAX, and they come from `_bj_ux_sin_cos` of 2θ = (θ/2^1023)·2^1024, since Rust's c² − s² was finite there.
    Rust's c² − s² from cos θ and sin θ lost every digit near a fast axis at 45°: a quarter-wave plate at the double
    π/4 had 2.2e-16 for cos(π/2)'s 6.1e-17, its m₀₀.im 1.57e-16 for 4.33e-17. sin 2θ replaces 2cs, which was
    accurate, so the two come from one call.
  - `BjMuellerMatrix_waveplate` takes 1 − cos δ as 2·sin²(δ/2) where cos δ > 1/2, and keeps Rust's 1 − cos δ
    elsewhere. 2·sin²(δ/2) is within about 1.5 ulps everywhere. Rust's difference is within about an ulp where
    cos δ ≤ 1/2 (the quarter- and half-wave plates sit there and keep Rust's bits), and loses 2^-53/(δ²/2)
    relative above. At δ = 1e-8, Rust's CS(1 − cos δ) is 0 and the port's 2.33e-17.
  - Kept: the Mueller form's C² + S²cos δ and S² + C²cos δ (an absolute error of about 2^-53, their zeros on
    curves of (θ, δ)), and its angle 2θ past DBL_MAX, NaN in Rust and the port alike (2d-1's).
- **The oracle**: `m_polarization.rs` appends a 1,380-row edge grid (`waveplate_edges`) after every earlier row, so
  the 23,135 earlier rows are an exact prefix of the new file. The grid takes the Jones and Mueller waveplates over
  35 fast axes × 19 retardances, and the quarter- and half-wave plates at each axis:
  - fast axes within two ulps of the π/4 family (π/4, 3π/4, −π/4, 5π/4) and of π/8, tiny ones, and DBL_MAX/2 and
    2^1023 on either side of the overflow;
  - retardances across cos δ = 1/2 (π/3 and its neighbours) and down to 1e-300.
- **validate.py**:
  - `pol_eval` takes the port's forms under `r33`, cos 2θ and sin 2θ through `pol_cos_sin_2t` (MPFR with an
    unbounded exponent past DBL_MAX/2);
  - `polarization_adr_rows` attributes row 33 where the forms alone explain the difference;
  - `polarization_corrections` skips a row an ADR 0002 mapping holds, as wave's does, so no row sits in both files
    (the harness's rule).

  The model agreed with every new row on its first run, and a one-bit flip of each new row's first output was caught
  on all 1,380.
- **glibc 2.44 misrounds cos(π/2 + 2 ulps)**: the new `polarization.cr.tsv` holds 18 rows, the Mueller retarders at
  θ = ±(π/4 + 1 ulp), where 2θ = 1.5707963267948968. 20 more such rows are row 33's mappings. The port gives the
  correctly rounded readings, gated at 0 ulp.
- **Mappings**: `polarization.adr.tsv` holds 2,340 rows: row 33's 2,279 (Jones waveplate 1,393, half-wave plate 70,
  quarter-wave plate 67, Mueller waveplate 749) and row 14's 61.
- **Tests**: `cyrius_waveplate_row33` in `tests/polarization.tcyr` (4 assertions, 190 in the file) pins MPFR's exact
  bits:
  - a quarter-wave plate at 45°;
  - a waveplate at θ = DBL_MAX;
  - the Mueller entry at δ = 1e-8;
  - a Mueller quarter-wave plate keeping Rust's bits.

  `cyrius_waveplate_forms_agree` (Jones against Mueller on nine (θ, δ) pairs) passes unchanged.
- **A trap** (cyrius 6.7.18): a hex f64 global whose exponent field is 0 (a subnormal, `0x0008…`) warns "an integer
  stored into an f64 slot". The port divides by 2^1023 instead of multiplying by 2^-1023.
- **Benchmarks**: Rust's benches build their waveplates outside the loop, and so do the port's twins. A one-off
  bench of the construction, interleaved against 1c-2's tree over three rounds at load 0.9–1.5:
  - `BjJonesMatrix_waveplate(0.3, 1.1)`: 535 → 452 ns (−16 %), one sin_cos of 2θ for Rust's cos θ and sin θ;
  - `BjMuellerMatrix_waveplate(0.3, 0.4)`, where cos δ > 1/2: 571 → 697 ns (+22 %), the one sine more;
  - `BjMuellerMatrix_waveplate(0.3, 1.1)`: 625 → 627 ns, unchanged.

  Every recorded row is within noise. The nanosecond rows move by one tick either way, and
  `json/mueller_matrix_from_json` by +4.2 % (8.77 → 9.14 µs), with no change to its code: a code-placement shift of
  the kind §60 measured.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 492,975 ok, 0 failed, 0
  measured, 0 pending, 70,834 accepted (ADR 01: 29,441; ADR 0002: 41,393)** over 563,809 rows; **all 380 digests
  equal**. The 374 outside the six waveplate functions are unchanged. `polarization.tcyr` passes on all three.

## 62. M3 part 2a-2: the resonators' exact stability and their waist without cancellation (2026-11-14)

- **ADR 0002 row 34, new** (the user's answer of 2026-11-14). beam.rs works on g = 1 − L/R rounded, which loses
  L/R's digits near the planar edge, and then on the rounded g₁g₂.
  - **Stability** (`_bj_beam_stable34`) is decided on the exact (1 − L/R₁)(1 − L/R₂):
    - its sign comes from each factor's, the sign of R − L, or + for a convex (R < 0) or flat (±∞) mirror;
    - ≤ 1 is ab ≤ a + b with a = L/R₁, b = L/R₂: always for two flat mirrors, R > 0 with one, and otherwise
      L ≤ R₁ + R₂ (L ≥ R₁ + R₂ when R₁R₂ < 0). L is compared with RN(R₁ + R₂), and at equality with TwoSum's exact
      error (`_bj_beam_cmp_sum`).

    The returned g₁, g₂ and g₁g₂ stay Rust's rounded values, so `is_stable` can disagree with g₁g₂'s rounding at an
    edge: two convex mirrors at L = 5e-324 have g = 1 + 1e-323 exactly, and are unstable though their g rounds to 1.
  - **The mode spacing** decides on the exact product of its own g's (`_bj_beam_gg_ok`): the sign from the factors'
    (a product that underflows to −0 is still negative), RN's side of 1, and at a rounded 1 the exact sign of
    g₁g₂ − 1 (src/crmath_cmp.cyr's `_bj_x_sign`). Its value keeps Rust's chain.
  - **The waist** (`_bj_beam_waist34`) comes from the Kogelnik form,
    w₀² = (λ/π)·√(L(R₁ − L)(R₂ − L)(R₁ + R₂ − L))/|R₁ + R₂ − 2L|. That is Rust's g₁g₂(1 − g₁g₂) and g₁ + g₂ − 2g₁g₂
    written in the inputs; with one flat mirror it is (λ/π)·√(L(R − L)), the limit as the other R → ∞.
    - The factors are differences and compensated sums of the inputs (`_bj_beam_sum3`: two TwoSums, so a sum that
      cancels is the exact sum rounded once), so none cancels.
    - Past 2^1020 they are summed scaled by 2^-4, so no sum overflows; the form is homogeneous, and the scale comes
      back exactly in the exponent.
    - The products run with an unbounded exponent (kimiya's `_bj_k_ux_*`), and w₀ is rounded once.
    - A fast path runs the same operations on doubles. Where every result is normal (kimiya's `_bj_k_ux_normal`, which
      also excludes 2^-1022 itself), each is the same 53-bit rounding of the same exact value, so its bits are the
      unbounded chain's; only an edge takes the slow path. The golden rows prove it: the digests are those of the
      unbounded chain alone.
    - The confocal limit √(Lλ/2π) at R₁ = R₂ = L and the planar refusal (two flat mirrors) keep Rust's branches. A
      near-planar cavity of finite radii, which Rust refused as planar when its g rounded to 1, gets its finite waist.
- **Accuracy**: on every golden row where the port computes a w₀, it is within one ulp of the exact value (538
  correctly rounded, 185 one ulp off). `validate.py` checks a bound of four on every row, against the exact fourth
  root (`bm_rn_root4`). Rust's w₀ is more than 1,000 ulps off on 170 rows. At L = 0.5, R₁ = 0.4, R₂ = 0.1 the rounded
  g₁g₂ is exactly 1, so Rust's w₀ is 0 for 2.2e-8 m (the exact product is 0.99999999999999965).
- **The oracle**: `m_beam.rs` appends an 861-row edge grid (`resonator_edges`) after every earlier row, so the 40,086
  earlier rows are an exact prefix of the new file. The cavities:
  - the confocal, half-confocal and concentric edges within ulps (L at R₁ + R₂ and its neighbours), and the
    convex–concave edge L = R₁ + R₂ with R₁R₂ < 0;
  - near-planar cavities of radii 10⁶ to DBL_MAX, concave, convex and mixed;
  - flat mirrors;
  - lengths and radii past 2^1020, and subnormal lengths;
  - 120 seeded cavities;

  each with three wavelengths, plus 16 mode-spacing g pairs (products that round to 1 or underflow).
- **validate.py**:
  - `bm_stable34` decides stability in Fractions, independently of the port's method;
  - `bm_waist34` transcribes the port's operations (TwoSums, the unbounded-exponent products, `rn_sqrt` for the
    last rounding);
  - `beam_adr_rows` attributes row 34;
  - `beam_corrections` skips a row a mapping holds.

  The model agreed with every new row on its first run, and a one-bit flip of each new row's first output was caught
  on all 861.
- **Mappings**: `beam.adr.tsv` holds 949 rows: row 32's 366 and row 34's 583. Of row 34's 583:
  - 432 waist values change;
  - 99 near-planar cavities get a waist where Rust refused them as planar, and 3 where it called them unstable;
  - 23 change their error (a convex–convex near-planar cavity is unstable, not planar);
  - 11 stability verdicts become unstable and 1 stable;
  - 14 mode spacings are refused (a product above 1 that rounds to 1, a negative one that underflows).
- **Tests**: `cyrius_resonator_row34` in `tests/beam.tcyr` (8 assertions, 121 in the file) pins the convex–convex
  verdict, a near-planar waist's correctly rounded bits, Rust's bits where they were right, the scaled path, the
  planar and confocal branches, and both mode-spacing edges with their reasons.
- **Benchmarks**: interleaved against 2a-1's tree, three rounds at load 0.5–1.2.
  - `beam/resonator_stability`, the twin: 26 → 39 ns (Rust 6.48 ns, so 4.2× → 6.0×), for the exact test's comparisons
    and TwoSum.
  - A one-off bench of the rest (Rust benches neither): `bj_resonator_beam_waist(0.1, 0.2, 0.25, λ)` 60 → 159 ns; with
    one flat mirror 59 → 96 ns; `bj_transverse_mode_spacing` 170 → 182 ns.
  - The waist's first form ran the unbounded chain on every call (324 ns). The bit-identical fast path took it to
    221 ns, and inlining the test and moving the slow chain into a helper of its own (the cost model: locals and calls
    cost even when unused) to 159.

  Every other row is within noise. Two A/B runs along the way ran beside another repo's benchmarks (load 3.4–8.7);
  these numbers are from the quiet one.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,253 ok, 0 failed, 0
  measured, 0 pending, 71,417 accepted (ADR 01: 29,441; ADR 0002: 41,976)** over 564,670 rows; **all 380 digests
  equal**. The 377 outside the three resonator functions are unchanged. `beam.tcyr` passes on all three.

## 63. M3 part 2a-3: transform_beam on the exact determinant (2026-11-14)

- **ADR 0002 row 35, new** (the user's answer of 2026-11-14). beam.rs's `transform_beam` and
  `transform_beam_with_waist_position` tested the rounded det = AD − BC for det > 0, and kept λ where
  |det − 1| ≤ 10⁻¹², a tolerance Rust chose to absorb rounding noise.
  - **The exact determinant** (`_bj_beam_det35`): X = AD and Y = BC are integers (two significands' product, at most
    106 bits, times a power of two) in src/fmtf.cyr's bignums (41 words, 1,280 bits).
    - A product more than 220 bits below the other only decides roundings as a sticky bit. It becomes a unit 230
      bits below the other, so every shift stays under 330 bits and every value under 2^500.
    - X − Y's sign decides det > 0.
    - det = 1 exactly when its magnitude is 2^−K (one bit, no bit below it), K the alignment's exponent.
  - **λ′**: λ is kept only where det is exactly 1. Otherwise λ′ = λ·det, rounded once: λ's significand times the
    bignum, its top 60 bits and a sticky bit for the rest, rounded by kimiya's `_bj_k_cr_from_int` (the subnormal
    grid, ties to even, +∞ past DBL_MAX). A zero, infinite or NaN λ (the beam's fields are public), which a positive
    det leaves as it is, stays.
  - **The fast path** (`_bj_beam_det35_fast`) answers first where it can, with the same results.
    - It takes entries that are 0 or of magnitude in [2^-450, 2^450], where TwoProduct is exact. AD − BC is then
      Shewchuk's exact expansion of p₁ + e₁ − p₂ − e₂ (crmath_gamma's `_bj_g_grow`): its sign is its largest
      component's, and it is 1 when it less 1 is empty.
    - λ·det comes from the compressed top two components (within 2^-104) times |λ| (`_bj_k_dd_mul_d`, 3.01u²), and
      kimiya's rounding test (`_bj_k_dd_round_ends`, a 32u² budget) proves the correctly rounded value.
    - It declines (to the bignums) for an entry outside the range, |λ·det| outside [2^-900, 2^900], or a value within
      about 2^-100 of a midpoint.
  - **Kept**: the reason quotes Rust's rounded det; `determinant` itself returns Rust's value; a non-finite entry
    keeps Rust's test.
- **The oracle**: `m_beam.rs` appends an 828-row edge grid (`abcd_det_edges`) after every earlier row, so the 40,947
  earlier rows are an exact prefix of the new file. It has 139 matrices through both transforms, each with three
  beams (the lab beam, and λ = 1e-300 and 1e300, so λ·det lands in range; 414 new rows a function, three pairs being
  in the grid already):
  - dets exactly 1 from integers and powers of two;
  - within an ulp of 1, rounding to 1, and within 10⁻¹² of it;
  - ±2^-104 exactly where the rounded det is 0;
  - past DBL_MAX and below 2^-1074;
  - 120 seeded near-singular matrices (D = RN(BC/A) and its next double).
- **validate.py**: `bm_transform(port=True)` takes the exact determinant in Fractions and λ′ by `rn_frac` (Python's
  correctly rounded integer division), independently of the port's bignums. `beam_adr_rows` attributes row 35. The
  model agreed with every new row on its first run, and a one-bit flip of each new row's first output was caught on
  all 828.
- **Mappings**: `beam.adr.tsv` holds 1,423 rows: row 32's 366, row 34's 583 and row 35's 474 (237 for each function).
  Of a function's 237:
  - 90 have λ·det rounded once, an ulp from Rust's double rounding;
  - 27 have a det within 10⁻¹² of 1 (but not 1) that now scales λ;
  - 119 Rust refused at the determinant. Of those, 103 are now refused at z_R′, because `propagate_beam`'s
    Im q_out = z·(AD − BC)/|Cq + D|² cancels in Rust's chain as the det did; 8 at its zero denominator; and 8 give a
    beam;
  - 1 is a beam with λ = DBL_MAX through det = 1 + 10⁻¹²: Rust's snap kept λ, and the port's λ·det overflows, which
    `BjGaussianBeam_new` refuses.

  The 8 new beams come from seeded near-singular matrices, and their z_R′ still comes from `propagate_beam`'s chain.
  Its cancellation in z·(AD − BC) is `propagate_beam`'s own to settle, with its threshold and range edges (2c, 2d).
- **Tests**: `cyrius_abcd_row35` in `tests/beam.tcyr` (4 assertions, 125 in the file) pins:
  - det = 1 + 1e-13 scaling λ, its bits from the model;
  - det = 1 from integers keeping λ;
  - det = −2^-104 refused with the rounded det's reason;
  - det = +2^-104 passing the determinant and refused at z_R′, out untouched.
- **Benchmarks** (Rust benches no transform; a one-off bench, interleaved against 2a-2's tree, three rounds at load
  1.1–2.0):
  - `transform_beam` with det = 1.2666…: 96 → 300 ns;
  - with a thin lens (det exactly 1): 94 → 221 ns.

  The bignums alone took the first to 781 ns. The fast path, which gives the same bits (the digests are unchanged by
  it), brought it to 300. Every other row is within this run's noise (median 2.6 %; the nanosecond rows move a tick).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,607 ok, 0 failed, 0
  measured, 0 pending, 71,891 accepted (ADR 01: 29,441; ADR 0002: 42,450)** over 565,498 rows; **all 380 digests
  equal**. The 378 outside the two transforms are unchanged. `beam.tcyr` passes on all three. **M3 Part 2a is done.**

## 64. M3 part 2b: round once (2026-11-14)

- **ADR 0002 row 13, amended for M3** (the user's "round once" decision of 2026-11-12: a formula that is one product,
  quotient or root of inputs and constants is correctly rounded). M3 has four:
  - beam's z_R = πw₀²/λ in `BjGaussianBeam_new` (`_bj_cr_mono_wide`). Its representability test now reads the
    rounded-once value: λ = 1e10, w₀ = 1e155 was refused as ∞ for z_R = 3e300, and λ = 5e-324, w₀ = 1e-162 gave 1.0
    for 0.636 (w₀² went subnormal);
  - `divergence`, λ/(πw₀), for finite fields and a nonzero waist (`_bj_cr_mono_wide`): λ = 1e-10, w₀ = 1e308 gave 0
    for 3.18e-319;
  - the confocal waist √(Lλ/2π) in `bj_resonator_beam_waist` (`_bj_cr_sqrt_mono`): L = λ = 1e200 gave ∞ for 4e199;
  - scattering's `bj_cross_section_from_efficiency`, Q·π·r², for finite operands (`_bj_cr_mul5` with a factor 1):
    Q = 1e308, r = 0.1 gave ∞ for 3.1e306, and Q = 1e-320, r = 1e10 lost four digits to a subnormal Q·π.

  An infinite or NaN operand keeps Rust's chain. `transform_beam` builds its output beam with `BjGaussianBeam_new`,
  so its z_R rounds once too. M3's other formulas have a sum or a computed value inside (the spot size, the curvature,
  the Rayleigh functions' |K|², the waveplates' and the Stokes sums) and are left to the range edges (2d).
- **The oracle**: each grid gets a round-once edge grid after every earlier row, the old files exact prefixes.
  - `beam.tsv` gains 513 rows: `new` where πw₀² over- or underflows though z_R is representable, `divergence` where
    π·w₀ does, the confocal waist where Lλ does, and 200 seeded laser beams and 100 confocal cavities, where the
    change is an ulp.
  - `scattering.tsv` gains 205 rows: Q·π·r² where Q·π or Q·π·r over- or underflows, and 200 seeded efficiencies and
    radii.
- **validate.py**:
  - `bm_new`, `divergence`, the confocal branch and `sc_eval` take the rounded-once forms under `port` (`rn_mono`,
    `rn_sqrt`: exact rationals rounded once);
  - `beam_adr_rows` attributes row 13 to them, and to a transform whose only change is its z_R (`r35=False` isolates
    it from row 35);
  - `scattering_adr_rows` attributes row 13 to the cross-section.

  The model agreed with every new row on its first run, and a one-bit flip of each new row's first output was caught
  on all 513 and all 205.
- **Mappings**:
  - `beam.adr.tsv` holds 1,968 rows. Row 13's 545 are `new` 115, `divergence` 60, the confocal waist 14 and the
    transforms 356, whose only change is z_R. Rows 32, 34 and 35 keep their 366, 583 and 474.
  - `scattering.adr.tsv` maps 93 rows, all row 13's.
- **Tests**:
  - `cyrius_beam_row13` in `tests/beam.tcyr` (5 assertions, 130 in the file) pins both `new` edges, the lab beam's z_R
    (correctly rounded, where Rust's was an ulp low), the divergence's subnormal and the confocal waist past Lλ's
    overflow.
  - `cyrius_cross_section_row13` in `tests/scattering.tcyr` (3 assertions, 105 in the file) pins both edges and ∞·0's
    NaN.
  - `cyrius_abcd_row35`'s transform now expects the waist that the lab beam's rounded-once z_R gives (…432D for
    …432C).
- **Benchmarks** (none of the four has a Rust bench; a one-off bench, interleaved against 2a-3's tree, three rounds at
  load 1.0–1.3):
  - `cross_section_from_efficiency` 5 → 128 ns;
  - `divergence` 5 → 110 ns;
  - `BjGaussianBeam_new` 33 → 179 ns;
  - the confocal waist 66 → 177 ns.

  That is the monomial machinery's cost (significand splitting, double-double products, the rounding test), and M1's
  and M2's round-once functions pay the same: `maxwell/wave_speed` 141 ns, `material/displacement_field` 231 ns.
  Every other row is within noise (median 0.2 %).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,687 ok, 0 failed, 0
  measured, 0 pending, 72,529 accepted (ADR 01: 29,441; ADR 0002: 43,088)** over 566,216 rows; **all 380 digests
  equal**. The 374 outside the six changed functions (`new`, `divergence`, the waist, the two transforms, the
  cross-section) are unchanged. `beam.tcyr` and `scattering.tcyr` pass on all three.
- **A slip, caught by the audit**: the scattering test group went in between `cyrius_json`'s doc comment and its
  function, leaving `cyrius_json` undocumented; the comment is back on its function.

## 65. M3 part 2c-1: Mie's exact |m|·x limit (2026-11-14)

- **ADR 0002 row 17, amended for M3** (the user's "exact thresholds" decision of 2026-11-12). `mie` refused
  |m|·x > 10⁷ on the rounded hypot(Re m, Im m)·x: two roundings, so the verdict followed them within an ulp of the
  limit. m = (60 + ulp) − 80i at x = 10⁵ has |m|·x = 10⁷ + 4·10⁻¹⁰, and Rust's hypot rounds it back to 100, so it was
  accepted (with a 10⁷-step recurrence). x = 4834.4, m = 367.9 − 2035.5i has (|m|·x)² 0.0019 below 10¹⁴, and its
  rounded |m|·x is 10⁷ + 1 ulp, so it was refused. The port decides (Re m·x)² + (Im m·x)² > 10¹⁴ exactly; the
  refusal still quotes the rounded |m|·x, which can read 10000000.
- **src/crmath_cmp.cyr's sum-of-squares sign**, `_bj_x_sumsq_sign`: the sign of (a·b + c)² + (e·f)² − 10^p for finite
  doubles. Bounds on |X| and |Y| from the rounded X̃ = a·b + c and Ỹ = e·f (each within 2^-53 of its operands, plus
  2^-1022) decide outside a 2^-40 band around the double nearest 10^p; inside it, or after a cancellation in X̃,
  `_bj_x_sumsq_exact` compares X², Y² and 5^p·2^p as integers (src/fmtf.cyr's bignums, 160-word buffers). It has
  propagate_beam's shape too (X = C·q₀ + D, Y = C·q₁), which 2c-2 takes up.
- **cr_verify.py's sumsq and sumsq_exact families** (the band then the exact test, and the exact test alone): sums
  within ulps of 10^p for p = −60, 14 and any p in [−60, 30], X cancelling to its product's rounding error, subnormal
  and tiny partners, Pythagorean ties (10^p = X² + Y² exactly) and their one-ulp neighbours, Mie- and
  propagate-shaped operands, and far ones from any binade. **200,000 calls equal to the exact sign on both targets**;
  two seeded mutations (the 10^p alignment off by one bit, the exact test replaced by the bounds' midpoint) were
  caught on 30 % and 13 % of calls.
- **The oracle**: `scattering.tsv` gains 77 rows after every earlier row (the old file an exact prefix): the knife edge
  (60 − 80i)·10⁵, |m|·x = 10⁷ exactly (accepted; the port runs it in full), and the one-ulp neighbours above the
  Pythagorean edges (60 − 80i)·10⁵, (6·10⁶ − 8·10⁶i)·1, (6·10⁵ − 8·10⁵i)·10 and (6·10³ − 8·10³i)·10³, then 60 seeded
  spheres whose index lies two to six ulps outside |m| = 10⁷/x. A sphere the exact limit accepts and Rust's refused
  runs a 10⁷-step recurrence that validate.py cannot transcribe, so the grid holds none (`scattering_adr_rows` stops
  on one); `cyrius_mie_row17` pins one in the tests.
- **validate.py**: `sc_mie_checks(port=True)` decides on Fractions; `sc_row17` names the rows the two verdicts split;
  `scattering_adr_rows` maps them to the port's refusal, and `scattering_corrections` skips them. The model agreed with
  every new row on its first run. A one-bit flip of each new row's first output was caught on the 68 refusals; the 9
  value rows are past SC_PY_START, so validate checks their primitives, and parity checks the tie against the port's
  own run.
- **Mappings**: `scattering.adr.tsv` holds 101 rows: row 13's 93 and row 17's 8, the spheres Rust accepted an ulp
  past the limit: (±(60 + ulp) − 80i)·10⁵; ±(6·10⁶ + ulp) − 8·10⁶i and 6·10⁶ − (8·10⁶ + ulp)i at x = 1; and the same
  three scaled to x = 10³.
- **Tests**: `cyrius_mie_row17` in `tests/scattering.tcyr` (2 assertions, 107 in the file): (60 + ulp − 80i)·10⁵
  refused with the text quoting 10000000, and the 4834.4 sphere accepted with its 4,905 coefficients (its D_n
  recurrence from 10⁷ costs about a second in the debug test).
- **propagate_beam is not in this bite.** The exact test makes |C·q + D|² ≥ 10⁻⁶⁰ pass where Rust's rounded
  denominator cancels to 0 (C = 0.1, q = 3, D = −0.30000000000000004: the exact C·q + D is −2.8·10⁻¹⁷), and Rust's
  chain then gives 0/0 = NaN. Which value it should give goes to the user before 2c-2.
- **Benchmarks** (interleaved against 2b's tree, three rounds at load 2.3–3.7, another repository's tests running):
  `scattering/mie_x1_real` 1.98 → 1.99 µs, `mie_x10_real` 4.02 → 3.96 µs, `mie_x50_m3_5` 12.83 → 12.70 µs. The band
  test is a few ns inside a µs-scale call, so the record keeps its rows. The median over the 81 rows is +0.8 %; the
  largest moves (±5–8 %) are on code this bite does not touch (`field/electric_potential_point` 25 → 27 ns, the JSON
  Mie readers −6 %): the load, and the bundle's layout.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,756 ok, 0 failed, 0
  measured, 0 pending, 72,537 accepted (ADR 01: 29,441; ADR 0002: 43,096)** over 566,293 rows; **all 380 digests
  equal**. Only `scattering::mie`'s changed (its 77 new rows); the other 379 are unchanged. `scattering.tcyr` passes
  on all three.

## 66. M3 part 2c-2: propagate_beam on the exact quotient (2026-11-14)

- **ADR 0002 row 17, amended for `propagate_beam`** (the user's answer of 2026-11-14: accurate q_out). The exact
  threshold alone would have accepted inputs where Rust's rounded C·q₀ + D cancels to 0, and Rust's chain would then
  give 0/0. So `propagate_beam` decides |C·q + D|² < 10⁻⁶⁰ exactly and rounds each part of q_out once from the exact
  quotient: Re = (U·X + V·Y)/(X² + Y²) and Im = q₁·(AD − BC)/(X² + Y²), with U = A·q₀ + B, V = A·q₁, X = C·q₀ + D and
  Y = C·q₁. A part exactly 0 keeps Rust's zero where Rust's chain gives one (+0 otherwise); a non-finite entry or q
  keeps Rust's test and chain.
- **src/crmath_edge.cyr's `_bj_cr_abcd`**, three paths:
  - C = 0 (free space, interfaces): the threshold is |D| < 1e-30 (no double lies between 10⁻³⁰ and it), Re is B/D,
    RN(A·q₀ + B) by TwoProduct and kimiya's `_bj_k_rn3`, or a double-double over D; Im is q₁·A, q₁/D or
    TwoProduct(q₁, A) over D;
  - otherwise, for entries and q in [2^-200, 2^200] or 0: X, U and AD − BC as double-doubles with error bounds,
    kimiya's products and AccurateDWPlusDW sums, the bounds carried to each quotient, and `_bj_e_div_ends`'s test
    when they total at most 200u²; the threshold first, by 2c-1's `_bj_x_sumsq_sign`;
  - the exact path: signed bignums with an exponent (302 words: any finite input spans at most 8,393 bits), and a
    restoring division to 55 bits and a sticky bit, rounded by `_bj_k_cr_from_int` (subnormals and ±∞ included).

  `_bj_x_sumsq_exact` now squares with `_bj_e_bn_mul`, and 2c-1's `_bj_xq_mul`, its duplicate, is gone.
- **cr_verify.py's abcd and abcd_slow families** (every path, and the exact one alone; each case's Re, Im and flags):
  composed systems, free space and interfaces with ties, cancellations at every depth, the threshold, midpoints of
  Im q_out, exact zeros, a lens followed by the distance to its focus, near-singular matrices, the box edges, scaled
  systems and any finite values. **150,000 calls equal to the exact value on both targets**, and the sumsq families'
  200,000 again. Seeded mutations caught: the fast rounding test skipped (200 calls of 15,000), the budget
  unbounded or Re's error bound zeroed (48), the sticky bit dropped (822; 1,069 on the exact path alone), C = 0's
  `_bj_k_rn3` replaced by p + B (10) and its Im low part dropped (17). Zeroing the bound on X, U, AD − BC, Im or the
  denominator alone changed nothing: under a cancellation those double-doubles come out exact (Sterbenz), and without
  one their error stays far inside `_bj_e_div_ends`'s 239u², so those bounds are the proof's, not the cases'.
- **The oracle**: `beam.tsv` gains 786 rows after every earlier row (the old file an exact prefix): 456
  `propagate_beam` rows (|C·q + D|² ulp by ulp around 10⁻⁶⁰ without and with cancellation, C·q₀ + D cancelling to
  2^(k − 104), 120 seeded cancellations, 80 lenses followed by the distance to their focus, 80 near-singular matrices,
  exact zeros and their signs, the range ends, 60 scaled systems) and 165 of each transform (|C·z_R|² around 10⁻⁶⁰
  with D = 0, 40 lenses to their focus, 120 seeded systems on the lab beam).
- **validate.py**: `bm_propagate(port=True)` decides and divides in Fractions (`rn_frac`), with the zero rule;
  `bm_transform` takes it unless `r17=False`; `beam_adr_rows` labels `propagate_beam` row 17, and a transform the
  first of rows 13, 35 and 17 whose fix alone accounts for its change (row 35 where only two together do). The model
  agreed with every new row on its first run; a one-bit flip of each new row's first output was caught on 783 of 786.
  The 3 others are NaN outputs, whose bits are not compared (row 6).
- **Mappings**: `beam.adr.tsv` holds 5,028 rows: row 17's 2,987 (`propagate_beam` 2,589, `transform_beam` 154, the
  waist-position form 244), row 13's 564, row 35's 528, rows 32 and 34 as before. Of the 2,589: Rust was an ulp off
  on 1,154 and more on 255, had the wrong sign or a zero on 562 (Re q_out at a focus, subnormal parts), gave NaN or ±∞
  for a finite value on 550, had one part past DBL_MAX beside a misrounded or NaN other part on 20, and refused 48 the
  exact threshold accepts.
- **Tests**:
  - `cyrius_abcd_row17` in `tests/beam.tcyr` (5 assertions, 135 in the file): the cancelled denominator
    (−1.08·10¹⁷, where Rust refused), |D| an ulp below 1e-30 refused, D = 1e-30's Im q_out an ulp below Rust's, 1e200
    everywhere giving 1 (Rust: NaN), and exact zeros keeping Rust's signs.
  - `cyrius_abcd_row35`'s det = 2^-104 transform now passes with z_R′ = 1.45·10⁻³³ (λ′ and w₀′ pinned), where it was
    refused at Rust's z_R′ = 0.
- **Benchmarks** (interleaved against 2c-1's tree, three rounds; the load fell from 21 to 4 during the run, another
  repository's gate finishing, so rows that didn't change moved with it, the median −11 %):
  - `beam/abcd_propagate_beam` (free space, the C = 0 path) 19 → 89 ns, 49.7× Rust's 1.79 ns (the record updated);
  - a one-off bench at load 5.7: `propagate_beam` through a lens (the double-double path) 21 → 650 ns, and
    `transform_beam` through it 575 → 1,210 ns.

  The general path is the call overhead of six TwoProducts, five double-double products and two quotient tests
  (`_bj_k_two_prod` alone takes 17 ns here): a roadmap finding, with inlining as its first step. Mie's rows are
  unchanged (within 2 %).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 491,482 ok, 0 failed, 0
  measured, 0 pending, 75,597 accepted (ADR 01: 29,441; ADR 0002: 46,156)** over 567,079 rows; **all 380 digests
  equal**. `propagate_beam`'s and both transforms' changed; the other 377 are unchanged. `beam.tcyr` passes on all
  three.

## 67. M3 part 2d-1a: the Mueller matrices at their range edges (2026-11-14)

- **2d-1 runs as three bites** (plan §0.1): 2d-1a the Mueller matrices, 2d-1b the Jones vectors and matrices
  (`intensity`, `normalized`, `inner`, `overlap`, `apply`, `compose`, `determinant`), 2d-1c the Stokes vectors
  (`from_jones`, the degree of polarization, the ellipticity and orientation angles). The list is §55's candidates
  and a second reading of polarization.rs against the rule: where an intermediate leaves the normal range and the
  result is representable, Rust's chain runs with an unbounded exponent and its last step rounds once.
- **ADR 0002 row 18, amended for the Mueller matrices**:
  - `linear_polarizer`, `rotation` and `waveplate` (and through it the quarter- and half-wave plates) take cos 2θ
    and sin 2θ from row 33's `_bj_pol_cos_sin_2t`: one `_bj_sin_cos` of the double 2θ (the bits of the two separate
    calls), or the reduction past DBL_MAX where |θ| ≥ 2^1023. Rust's 2θ was ∞ there, and every entry it entered NaN.
  - `waveplate`'s CS·2·sin²(δ/2), row 33's form where cos δ > 1/2, rounds once from the 53-bit square where
    sin²(δ/2) falls below 2^-1022, |δ| below about 2^-510 (`_bj_pol_cs_om`). CS is always its 53-bit value: |C| and
    |S| stay above about 2^-62 for any double 2θ, except where 2θ is tiny, and there C = 1 and CS = S exactly.
  - `apply` and `compose` sum four products from −0.0, Rust's f64 Sum, inline and unrolled (−0.0 + p₀ is p₀), then
    test the sum against a band. Outside it, unless a row or column of zeros makes every product exact,
    `src/crmath_vx.cyr`'s `_bj_vx_dot4_edge` keeps a finite sum whose four products are exact, and otherwise
    `_bj_vx_dot4_ux` runs the chain with an unbounded exponent (kimiya's `_bj_k_ux_add`, the last sum by
    `_bj_k_ux_add_f64`). The helpers live in crmath_vx, a file of internals, so `cr_driver.cyr` can call them; its
    `_drv_dot4` is the inline chain and test.
- **The guard's argument** (crmath_vx's new section). A product below 2^-1022 is within 2^-1075 of its 53-bit value
  and at most 2^-1022 in magnitude. A partial sum or term of 2^-968 or more absorbs it in both chains (5d-2's
  argument). Otherwise the partial sums the two chains disagree on stay below 2^-967, then 2^-911, then 2^-855, as
  each later term fails to absorb them. So a finite sum of at least 2^-800 is the chain's: 5d-2's three-term bound,
  2^-900, is not enough for four terms. A finite sum whose products are all exact is the chain's at any size, since
  a sum of normal values below 2^-1022 is exact.
- **The oracle**: `polarization.tsv` gains 558 rows after every earlier row (`mueller_edges`; the old file an exact
  prefix):
  - the angle constructors at fast axes from DBL_MAX/2 to DBL_MAX of either sign (36 rows);
  - the waveplate at six fast axes and 13 retardances around sin²(δ/2) = 2^-1022 and 2^-1075, and at the large
    axes (124 in all);
  - 14 Mueller matrices against 15 Stokes vectors (`apply`, 210 rows) and against each other (`compose`, 188). They
    hold products past DBL_MAX that cancel, partial sums past it, products of 2.25·2^-1074, a product that rounds up
    to 2^-1022 (w = 1 − 2^-53 times 2^-1022), a last sum that rounding to 53 bits first would put on a midpoint
    (2^-1030 ± 2^-1075 ± 2^-1090), and small sums of exact products.
- **validate.py**: `pol_eval` takes `r18`, the port's forms under `wp`: the angles through `pol_cos_sin_2t`, the
  waveplate's square through `rnd53`, and every `apply` and `compose` sum with finite inputs through `pol_dot4`, the
  unbounded chain itself. So the harness checks the port's guard on every row. `polarization_adr_rows` labels a row
  18 where the range edges alone account for its change, before rows 33 and 14. The harness passed on its first run,
  and a one-bit flip of each new row's first output was caught on 547 of 558; the other 11 are NaN outputs, whose
  bits are not compared (row 6).
- **cr_verify.py's dot4 family**: `_drv_dot4` on 50,000 cases, with b contiguous (`apply`'s) and at a stride of 4
  (`compose`'s column), against exact.py's chain:
  - products anywhere in the range, with zeros of both signs;
  - products near 2^-1022, some below it, with normal terms among them;
  - products past DBL_MAX, cancelling or not;
  - sums cancelling to small values that subnormal products then meet;
  - products just below 2^-1022 that round up to it.

  **100,000 calls per target, 0 mismatches on x86_64 and aarch64.** 6,643 cases differ from Rust's sum, none of them
  inside the band.
- **Mutations**, each run through the parity harness, and the rows that caught it:
  - the band dropped: 24;
  - a product of 2^-1022 taken as exact: 1;
  - an overflowed sum kept: 24;
  - the exact-product test skipped: 24;
  - the square not doubled: 54;
  - the waveplate's edge switched off: 6;
  - the last sum rounded to 53 bits first: 1, the midpoint row, added for it after no row caught it;
  - the zero-row shortcut on two entries of the row: 5, and the zero-column one: 6.
- **Mappings**: `polarization.adr.tsv` holds 2,619 rows:
  - row 33's 2,321: 42 more Mueller waveplate rows from the new grid, 4 of which the 53-bit square changes
    (δ = 1e-155 at θ = 0.3: 2.330097714918e-311 for 2.3300977149184e-311);
  - row 18's 237: the waveplate 149 (all at |θ| ≥ 2^1023), `apply` 31, `compose` 17, and the linear polarizer, the
    rotation and the two plates 10 each. 107 of them are older rows, NaN in Rust and the port alike before;
  - row 14's 61.

  Rust gave NaN for every entry that 2θ enters on the 189 angle rows. On the 48 sums it gave ±∞ or NaN for a finite
  value (19 rows), NaN for ±∞ (7), the wrong zero sign (14), 0 for a nonzero value or the reverse (2), and a misrounded
  subnormal (9); a row can show several.
- **Tests**: `cyrius_mueller_row18` in `tests/polarization.tcyr` (7 assertions, 197 in the file):
  - the linear polarizer at θ = DBL_MAX and the rotation by −2^1023 (MPFR's bits);
  - the waveplate at δ = 1e-155, one unit below the twice-rounded square;
  - DBL_MAX + DBL_MAX − DBL_MAX = DBL_MAX;
  - three products of 2.25·2^-1074 summing to 7·2^-1074;
  - 10^200·10^200 − 10^200·10^200 = +0;
  - w·2^-1022 − 2^-1022 = −0.
- **Benchmarks.** A first form called a shared `_bj_vx_dot4` for each output. A full bench interleaved against 2c-2's
  tree (three rounds, load 4.8 to 7.3) put `polarization/mueller_apply` at 177–190 ns for 91–97 in its two quieter
  rounds: a 6-argument call per output, and the edge test's call for each of the polarizer's two zero rows. A
  micro-bench of three variants settled the form: `apply` and `compose` sum inline and unrolled, and test a row (or
  column) of zeros before calling the edge test. The final tree, interleaved with HEAD's three times (load 6.7):
  - `mueller_apply` on the bench's polarizer: 99–106 → 116–122 ns, the band and zero-row tests;
  - a dense `apply` (a waveplate on a general Stokes vector): 98–106 → 100–105 ns, unchanged;
  - a dense `compose`: 449–487 → 359–376 ns (−21 %), each row of the left matrix read once, not once per column.

  A second full bench, interleaved against 2c-2's tree as the load fell from 11 to 5.6, read `mueller_apply` 94 →
  107 ns in its quietest round. The record takes that 107 ns (108× Rust's 0.99 ns, whose arguments LLVM folds); every
  other row is within the load's drift (a median of +6 % over all 81, and within ±3 % in that round).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 491,761 ok, 0 failed, 0
  measured, 0 pending, 75,876 accepted (ADR 01: 29,441; ADR 0002: 46,435)** over 567,637 rows; **all 380 digests
  equal**. The seven Mueller constructors' and sums' digests changed (`linear_polarizer`, `rotation`, `waveplate`, the
  two plates, `apply`, `compose`); the other 373 are unchanged. `polarization.tcyr` passes on all three. The gate is
  clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 68. M3 part 2d-1b: the Jones vectors and matrices at their range edges (2026-11-14)

- **ADR 0002 row 18, amended for the Jones vectors and matrices.** Seven functions run hisab's complex chain in Rust's
  order: `intensity` ((re² + im²) + (re² + im²)), `inner`, `apply`, `compose` and `determinant` (two complex products,
  each four real products then a difference and a sum, added or subtracted in a pair), and `normalized` and `overlap`
  through them.
  - **The guard is row 14's box, on the inputs.** With every part ±0 or 2^-511 ≤ |x| < 2^511, every real product is
    normal or an exact 0, each pair's sum lies below 2^1023, and a sum below 2^-1022 is exact, so Rust's chain is the
    unbounded one up to its last operation (class 0). A finite part outside the box is an edge (class 1): the chain
    runs with an unbounded exponent and each output's last operation rounds once. A part that is not finite (class 2)
    leaves an output to BjComplex's row-14 products, as before, since the non-finite part decides its value.
  - **`normalized`** keeps Rust's first branch where the intensity i lies in [1e-290, DBL_MAX]: i is then the chain's
    (a part below 2^-511 adds a square under 2^-1022, which the larger square absorbs), and √i and 1/√i are normal.
    Elsewhere Rust rescaled by its largest part; the port runs the first branch with an unbounded exponent there (i,
    √i and 1/√i to 53 bits, each product rounded once), as row 18 did for `FieldVector::normalized` (§21). So the
    result is exactly scale-invariant: 2^-600·v gives v's bits.
  - **`overlap`** keeps Rust's first branch where both vectors are boxed, I₁·I₂ lies in [1e-290, DBL_MAX] and
    ⟨J₁|J₂⟩'s parts are boxed too (its squares normal or 0). Elsewhere, where a square under- or overflowed or Rust
    rescaled both vectors, the first branch runs with an unbounded exponent: the intensities, their product, ⟨J₁|J₂⟩
    and its |·|² to 53 bits, the quotient rounded once. The I₁·I₂ bound is slack: by Cauchy–Schwarz,
    |⟨J₁|J₂⟩|² ≤ I₁·I₂, so an I₁·I₂ below 2^-1022 forces ⟨J₁|J₂⟩ out of the box, and between 2^-1022 and 1e-290 Rust's
    first branch is the chain's.
  - Rust's errors are unchanged: a part not finite is InvalidParameter, a zero vector DivisionByZero, and `overlap`
    checks the first vector before the second.
- **src/crmath_vx.cyr's Jones section**: `_bj_vx_cls4` (the box, by exponent fields with one fast test),
  `_bj_vx_cdot`, `_bj_vx_inten`, `_bj_vx_jnorm` and `_bj_vx_overlap` on the parts' addresses, and the unbounded
  helpers `_bj_vx_cprod_ux`, `_bj_vx_inten_ux` and `_bj_vx_sq2`. polarization.cyr's methods call them; `apply`,
  `compose`, `inner` and `determinant` keep BjComplex's products for an output with a non-finite input.
  `_BJ_POL_SAFE_MIN` and `_bj_pol_rescaled` are gone.
- **The oracle**: `polarization.tsv` gains 1,631 rows after every earlier row (`jones_edges`; the old file an exact
  prefix): 24 Jones vectors and 12 matrices, each vector's intensity and normalization (22 rows each), every pair's
  inner product and overlap (576 each), each matrix's determinant (9), and its products with every vector (`apply`,
  288) and every matrix (`compose`, 138). They hold squares and products below 2^-1022 that round twice, products past
  DBL_MAX that cancel, intensities around 1e-290 and past DBL_MAX, subnormal vectors, a NaN beside an edge row, boxed
  pairs whose inner product lies below 2^-511, and a complex last sum that rounding to 53 bits first would put on a
  midpoint.
- **validate.py**: `pol_eval` under `r18` takes every Jones output with finite inputs through the unbounded chain
  (`cdot`, `inten_ux`), `normalized` and `overlap` through their unbounded first branches after Rust's errors
  (`jerr`), so the harness checks every box and fast path row by row. `polarization_adr_rows` labels a row 18 when
  turning row 18 off, or rows 14 and 18 together, restores Rust's value: the 21 rows row 14 held for `inner`,
  `apply`, `compose` and `determinant` are row 18's now. The harness passed on its first run.
- **cr_verify.py's jcdot, jint, jnorm and jover families**: the four kernels against exact.py's chain, on parts inside
  the box (7,438 of the 20,000 jcdot calls entirely), at its edges, subnormal and huge; products cancelling;
  intensities near 1e-290 and past DBL_MAX; nearly orthogonal pairs, 2,228 of the 10,000 overlaps on Rust's first
  branch. **80,000 calls per target, 0 mismatches on x86_64 and aarch64**, and dot4's 100,000 again.
- **Mutations**, each run through the parity harness:
  - the box's top at 2^1023: 22 failing rows; its bottom at 2^-1022: 6;
  - `normalized`'s fast path below 1e-290: 21;
  - `intensity` always on Rust's chain: 2;
  - 1/√i rounded into binary64 first: 317; |⟨J₁|J₂⟩|² rounded into binary64 first: 534;
  - the conjugate dropped from the unbounded inner product: 795;
  - `overlap`'s box on ⟨J₁|J₂⟩ dropped: 4, and `cdot`'s last sum rounded to 53 bits first: 1. No row caught these two
    at first; the boxed small-inner pairs and the complex midpoint row were added for them;
  - `overlap`'s I₁·I₂ bound lowered to 0: not caught, and not catchable (the slack bound above).
- **Mappings**: `polarization.adr.tsv` holds 3,186 rows: row 33's 2,321, row 18's 825 (the Mueller matrices' 237 and
  the Jones vectors' and matrices' 588), row 14's 40 (`from_jones`, 2d-1c's). The 588: `overlap` 422, `normalized`
  60, `inner` 47, `apply` 33, `compose` 19, `determinant` 5, `intensity` 2; 307 are older rows. What Rust gave on them
  (a row can show several):
  - NaN or ±∞ for a finite value: `overlap` 84, `inner` 23, `apply` 15, `compose` 14, `normalized` 5 (the rescale's
    1/5e-324 = ∞: (5e-324, 0) normalized to NaN for (1, 0)), `determinant` 2;
  - NaN for ±∞: 14 across the products;
  - the wrong zero sign: `inner` 16, `apply` 7, `compose` 4;
  - an ulp off: `overlap` 196, `normalized` 53, `inner` 6, `apply` 5, `intensity` 2, `compose` 1, `determinant` 1;
    2 to 16 ulps: `overlap` 134, `normalized` 2; more than 16: `overlap` 4; 0 for a nonzero value: `overlap` 4.
- **Tests**: `cyrius_jones_row18` in `tests/polarization.tcyr` (5 assertions, 202 in the file):
  - the intensity of three squares of 2.25·2^-1074: 7·2^-1074;
  - `normalized(2^-600·v)` = `normalized(v)` bit for bit, where Rust's rescale left a part an ulp off;
  - ⟨(1e200, 1e200)|(1e200, −1e200)⟩ and det [[1e200, 1e200], [1e200, 1e200]] = +0;
  - the overlap of a boxed pair whose inner product lies below 2^-511 (5.085658759397545e-36; Rust 5.0856527…e-36);
  - `apply`'s last sum 2^-1030 + 2^-1075 + 2^-1090 rounded once to 2^-1030 + 2^-1074.
- **Benchmarks.** The kernels replace BjComplex's per-product calls (four `BjComplex_mul` and two adds, each product
  with its own range test) with one box test per pair of complex inputs and the arithmetic inline:
  - a full bench interleaved against 2d-1a's tree (three rounds, load 3.6 to 5.4): `polarization/jones_apply_polarizer`
    229–242 → 173–193 ns (−25 %), `jones_qwp_apply` 229–248 → 179–196 ns (−22 %); every other row within the load's
    noise (median ratio 1.00 over 81). The record takes the medians, 175 and 183 ns;
  - a micro-bench, interleaved with HEAD twice: `inner` 123 → 79–85 ns, `normalized` 93 → 72–75, `overlap` 242–246 →
    201–213, `compose` 391–393 → 354–377, `intensity` unchanged (28 ns), and `determinant` 80 → 89–96 ns, one call of
    eight arguments for two struct-valued products.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 492,825 ok, 0 failed, 0
  measured, 0 pending, 76,443 accepted (ADR 01: 29,441; ADR 0002: 47,002)** over 569,268 rows; **all 380 digests
  equal**. The seven Jones functions' digests changed; the other 373 are unchanged. `polarization.tcyr` passes on all
  three. The gate is clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 69. M3 part 2d-1c: the Stokes vectors at their range edges (2026-11-14)

- **ADR 0002 row 18, amended for the Stokes vectors.** This closes 2d-1, so every polarization chain now follows the
  rule.
  - **`from_jones`** is hisab's chain on a Jones vector's parts: |x|², |y|² and x·y*, then S₀ = |x|² + |y|²,
    S₁ = |x|² − |y|², S₂ = 2 Re(x·y*) and S₃ = 2 Im(x·y*). It uses 2d-1b's box, and runs the chain with an unbounded
    exponent where a finite part leaves the box, each output rounded once; the doubling is exact there
    (`_bj_vx_stokes`). A part that is not finite keeps BjComplex's chain.
  - **The degree of polarization** keeps Rust's hypot(hypot(S₁, S₂), S₃)/S₀ where both hypots are their 53-bit values:
    an exact 0, or a finite value above 2^-1022 (`_bj_vx_sure`; a hypot rounded up to 2^-1022 may lie below it).
    Otherwise each hypot runs to 53 bits with an unbounded exponent and the quotient rounds once (`_bj_vx_dop`).
    `_bj_vx_hypot_ux` scales both operands by one power of two, the larger into [1, 2), and takes crmath's correctly
    rounded hypot, since rounding commutes with that scaling. An operand 2^600 times smaller leaves the larger: it adds
    less than 2^-1199 relative, and the larger is a 53-bit number.
  - **The ellipticity angle** keeps Rust's 0.5·atan2(S₃, √(S₁² + S₂²)) where S₁ and S₂ are boxed (their squares
    normal or 0, so the root is exact up to its rounding) and atan2's value is its 53-bit one. Otherwise S₁² + S₂², its
    root and atan2 run to 53 bits unbounded and the half rounds once (`_bj_vx_ellip`). **The orientation angle** does
    the same for an angle below 2^-1022 (`_bj_vx_orient`; S₁ > 0 there).
  - **`_bj_vx_atan2_ux`** gives atan2(y, x) to 53 bits for x ≥ 0 with an unbounded exponent:
    - where |y/x| < 2^-59, the 53-bit quotient y/x: atan2 = (y/x)(1 − δ) with δ < 2^-118, and a quotient of two 53-bit
      numbers that is not one itself lies at least 2^-106 (relative) from a midpoint of the 53-bit grid;
    - where |y/x| > 2^1000, ±π/2's;
    - between, both scaled into binary64's range, where atan2 is correctly rounded.
- **The oracle**: `polarization.tsv` gains 125 rows after every earlier row (`stokes_edges`; the old file an exact
  prefix). 2d-1b's Jones edge vectors moved into `jones_edge_vectors`, unchanged, so `stokes_edges` can take them:
  - `from_jones` on those 24 vectors (22 rows new);
  - 26 Stokes vectors for the degree of polarization and both angles: hypots past DBL_MAX in either hypot, hypots
    below 2^-1022, S₁² + S₂² past DBL_MAX and below 2^-1074, angles below 2^-1022, zeros and non-finite parts;
  - angles below 2^-1022 off the subnormal grid, S₂/S₁ and S₃/S₁ for twelve odd S₁ (24 rows);
  - a hypot whose operands lie 2^21.5 apart (S = (DBL_MAX, DBL_MAX, DBL_MAX, 2^-21·DBL_MAX)).
- **validate.py**: `pol_eval` under `r18` takes `from_jones` with finite parts through the unbounded chain, and the
  degree of polarization and the angles from MPFR at 53 bits with its widest exponent range (`mpfr_ref.wide`'s hypot
  and atan2, independent of the port's scaling), halved by `pol_half`. The harness passed on its first run.
- **cr_verify.py's sjones, sdop, sellip and sorient families**: each part of `from_jones` against exact.py's chain,
  and the degree of polarization and both angles against `mpfr_ref.wide`, on parts around the box, near DBL_MAX and
  subnormal, hypot operands 2^10 to 2^40 apart, and angles below 2^-1022. **70,000 calls per target, 0 mismatches on
  x86_64 and aarch64**; Rust's chain differs from the exact one on 4,233 of them. The run that checked them also
  checked dot4's and the Jones families' 180,000 again: 250,000 calls per target, 0 mismatches.
- **Mutations**, each run through the parity harness:
  - `from_jones` without its doubling: 575 failing rows; with S₁ as a sum: 493; without the conjugate: 991;
  - the degree of polarization's guard off: 7;
  - the ellipticity angle's box off: 19, and its 53-bit test off: 2;
  - the orientation angle's guard off: 2;
  - `_bj_vx_atan2_ux`'s small-quotient branch off: 34, and that quotient rounded into binary64 first: 4;
  - `_bj_vx_hypot_ux` scaling by one power of two too many: 8;
  - its far-ratio shortcut at 2^20 instead of 2^600: 1, the row added for it after no row caught it;
  - the first hypot's test dropped: not caught. A subnormal first hypot moves a normal second one by less than
    2^-1075, so only a near-tie could show it; the test stays, since the proof needs it.
- **Mappings**: `polarization.adr.tsv` holds 3,318 rows: row 33's 2,321 and row 18's 997. Row 14 maps no
  polarization row now: its 40 `from_jones` rows are row 18's. The Stokes 172 are `from_jones` 141, the ellipticity
  angle 21, the degree of polarization 8 and the orientation angle 2; 141 are older rows. What Rust gave on them (a row
  can show several):
  - `from_jones`: NaN for ±∞ on 46 rows, NaN or ±∞ for a finite value on 25, the wrong zero sign on 31, 0 for a
    nonzero value or the reverse on 4, an ulp off on 35 and 2 to 16 on 5;
  - the degree of polarization: ∞ for a finite value on 4 (√2 at S = (DBL_MAX, DBL_MAX, DBL_MAX, 0)), more than 16 ulps
    off on 3 (0.57724 for 1/√3 on subnormal parts), an ulp off on 1;
  - the ellipticity angle: 0 for a nonzero value on 12 (S₁² past DBL_MAX: 0 for 0.3077), more than 16 ulps off on 7
    (S₁² below 2^-1074: π/4 for π/8), an ulp off on 2;
  - the orientation angle: an ulp off on 2 (halves of angles below 2^-1022).
- **Tests**: `cyrius_stokes_row18` in `tests/polarization.tcyr` (6 assertions, 208 in the file):
  - `from_jones`' S₀ of three squares of 2.25·2^-1074: 7·2^-1074;
  - the degree of polarization √2 with hypot past DBL_MAX, and 1/√3 on subnormal parts;
  - the ellipticity angle 0.3077 with S₁² past DBL_MAX, and π/8 with S₁² below 2^-1074;
  - an orientation angle below 2^-1022 with its half rounded once.
- **Benchmarks.** `from_jones` replaces BjComplex's calls (two `norm_sq`, a `conj` and a `mul`, each with its range
  test) with the shared kernel's single box test:
  - a full bench interleaved against 2d-1b's tree (three rounds, load 2.7 to 3.3): `polarization/stokes_from_jones`
    120–147 → 56–59 ns (−55 %); every other row within the load's noise (median ratio 1.00 over 81). The record takes
    59 ns;
  - a micro-bench, interleaved with HEAD four times: the degree of polarization 345–394 → 362–399 ns and the
    orientation angle 319–364 → 327–445 ns (both within about 4 % once the load settled), the ellipticity angle
    291–329 → 320–341 ns (+5 to 10 %: its guard calls `_bj_vx_box1` twice and `_bj_vx_sure` once).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 492,818 ok, 0 failed, 0
  measured, 0 pending, 76,575 accepted (ADR 01: 29,441; ADR 0002: 47,134)** over 569,393 rows; **all 380 digests
  equal**. The four Stokes functions' digests changed; the other 376 are unchanged. `polarization.tcyr` passes on all
  three. The gate is clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 70. M3 part 2d-2: Rayleigh scattering at its range edges (2026-11-14)

- **ADR 0002 row 18, amended for Rayleigh.** σ = (8/3)·π·r·r·x⁴·|K|² (x = 2π·r/λ), Q = (8/3)·x·x·x·x·|K|² and
  K = (m² − 1)/(m² + 2) by scattering.rs's Smith division are long products whose factors over- and underflow where σ
  is representable. Rust's σ for r = 10⁻²⁰⁰ at λ = 10⁻³⁰⁰ was 0·∞ = NaN for 543.4 m², since (8/3)πr² underflowed
  while x⁴ overflowed, and |m| past 10¹⁵⁴ made K = ∞/∞ = NaN where K → 1.
  - **The guard is a box on the inputs**: r and λ in [2^-100, 2^100] for σ, x in [2^-200, 2^200] for Q, and each part
    of m ±0 or in [2^-40, 2^40].
    - x then lies in [2^-202.7, 2^202.7], and the product before |K|² in [2^-1007.5, 2^1013.7] (Q's in [2^-798.6,
      2^801.4]).
    - m², m² − 1 and m² + 2 are 0 or multiples of 2^-132. Smith's sums are 0 or at least 2^-345: a cancelling pair has
      terms of at least 2^-79, a lone term is at least 2^-345. Its denominator lies in [2^-132, 2^82], so |K| lies in
      [2^-427, 2^214] and |K|², the chain's last factor, is normal.
    - A bound of 2^-50 on m was not enough. In Smith's second branch a sum can fall to 2^-405, which leaves |K|² only
      eight binades above 2^-1022 in the worst case.
  - Inside the box Rust's chain is the unbounded one up to its last product. Outside it, for a finite m, the chain runs
    with an unbounded exponent (`_bj_vx_ray_k2_ux`: hisab's m·m, Smith's division on its branch decided by the exact
    |b.re| ≥ |b.im|, then re² + im²) and its last product rounds once (`_bj_vx_ray_sigma`, `_bj_vx_ray_q`). A non-finite
    m keeps Rust's chain: K is NaN either way. BjComplex's row-14 m·m no longer enters Rayleigh, and Smith's division
    moved into `src/crmath_vx.cyr` as `_bj_vx_cdiv`; Mie calls it unchanged.
  - **The phase function** 0.75·(1 + cos²θ) leaves the range only where cos²θ passes DBL_MAX, near |cos θ| = 2^512.
    Rust's value is kept when it is finite (every step is then its 53-bit value); otherwise, for |cos θ| ≥ 2^511, the
    sum runs to 53 bits and the product rounds once (`_bj_vx_ray_phase`). Rust gave ∞ for 1.47·10³⁰⁸ at cos θ =
    1.4·10¹⁵⁴.
  - **The angular intensity** σ/(d·d)·P/(4π) keeps Rust's chain where σ/d² and its product by P are normal (then d², P
    and every input are finite, and each step is its 53-bit value), for σ = 0 with d and P finite, and for a
    non-finite input. Otherwise it runs with an unbounded exponent, P's chain included (`_bj_vx_ray_ang`). A first
    guard also bounded |d| and |cos θ| by 2^511; both bounds were redundant, and dropping them took the cost from 44
    to 21 ns.
- **The oracle**: `scattering.tsv` gains 856 rows after every earlier row (`rayleigh_edges`; the old file an exact
  prefix):
  - σ over 8 radii, 6 wavelengths and 13 indices (624 rows), at the boxes' edges and far past them;
  - Q over 8 size parameters and the same indices (100);
  - the phase function at 8 cosines around 2^511 and 2^512;
  - the angular intensity over 6 cross-sections, 7 distances and 3 cosines, and 4 rows where σ/d² falls below 2^-1022
    and a large P lifts it back (124).

  The indices hold the box's edges (2^40 and just below, 2^-41, 3 − 2^-41 i), |m| past 10¹⁵⁴ and 10²⁰⁰, tiny ones,
  m = 1 + ε and 1 − 10⁻¹⁷⁰ i (K cancelling), and −√2 i.
- **validate.py**: `sc_eval(r18)` takes every Rayleigh row with finite inputs through the unbounded chain (`sc_k2_ux`,
  `sc_phase_ux`), so the harness checks the boxes row by row; `scattering_adr_rows` labels those rows 18. The harness
  passed on its first run.
- **cr_verify.py's rsig, rq, rphase and rang families**: σ, Q, the phase function and the angular intensity against
  exact.py's chain, on radii, wavelengths and size parameters around the boxes and far past them, indices around
  theirs, cosines around 2^511, and cross-sections and distances whose σ/d² and its product by P leave the range.
  **40,000 calls per target, 0 mismatches on x86_64 and aarch64**; Rust's chain differs from the exact one on 2,158.
- **Mutations**, each run through the parity harness. Seven of sixteen were caught:
  - σ's box on r and λ off: 9 failing rows; the box on m off: 50; Q's box off: 37;
  - the angular intensity's test on σ·P/d² off: 19; its test on σ/d² off: 3, the rows added for it after no row
    caught it; its σ = 0 case without the test on P: 14;
  - the phase function's fast path taken for an infinite P: 2.

  The nine others:
  - slack bounds: σ's box at 2^-101, m's at 2^±60, the phase function's at 2^512, and the angular intensity's
    distance bound;
  - the σ = 0 shortcut, which only saves time;
  - b.im = m².im ± 0, a zero's sign that |K|² squares away;
  - Smith's branch at an exact tie |b.re| = |b.im|;
  - the last product (σ's) and the phase function's last product rounded to 53 bits first. Those differ only where a
    53-bit product lands on a subnormal midpoint or the overflow boundary: a search of 20,000 random Q's found none,
    as expected for a chance near 2^-37 per case.
- **Mappings**: `scattering.adr.tsv` holds 426 rows: row 18's 325 (σ 144, Q 95, the angular intensity 84, the phase
  function 2; 106 are older rows), row 13's 93, row 17's 8. What Rust gave on them (a row can show several):
  - NaN or ±∞ for a finite value: σ 103, Q 67, the angular intensity 65, the phase function 2;
  - NaN for ±∞: σ 30, Q 12, the angular intensity 1;
  - 0 for a nonzero value or the reverse: the angular intensity 15, σ 11, Q 2;
  - misrounded: Q an ulp off on 9, 2 to 16 on 3, more on 2; the angular intensity 2 to 16 on 1, more on 2.
- **Tests**: `cyrius_rayleigh_row18` in `tests/scattering.tcyr` (5 assertions, 112 in the file): σ = 543.4 m² for
  r = 10⁻²⁰⁰ at λ = 10⁻³⁰⁰, σ for m = 10²⁰⁰ − 10²⁰⁰i, Q for m = 1.4·10¹⁵⁴ − 10¹⁵⁰i, P(1.4·10¹⁵⁴) = 1.47·10³⁰⁸, and
  the angular intensity through a subnormal σ/d².
- **Mie's tiny m**, found on the way and not changed here. The audit's candidate "Mie's complex products, reachable
  with a tiny or huge m" was not in 2d's plan. Reading the recurrences against Rust's chain:
  - at ordinary x, Rust refuses |m| below about 10⁻¹⁵⁵ ("Mie series over/underflowed": inv_m·Dₙ overflows), while the
    m → 0 limit is finite (Q_ext = 0.2769 at x = 1; m = 10⁻¹⁰⁰ still gives it);
  - a huge m is refused only below x = 10⁻⁵⁰, the limit the user kept on 2026-11-14.

  The tiny m is the user's call.
- **Benchmarks.** σ and Q no longer call BjComplex: their box test is a few compares on the inputs.
  - A full bench interleaved against 2d-1c's tree (three rounds, load 4.4 to 6.3): `scattering/rayleigh_cross_section`
    143–151 → 82–85 ns (−41 %), `rayleigh_phase_function` 4 → 6–13 ns (median 7), Mie's three rows within 2 %, and
    every other row within the load's noise (median ratio 1.00 over 81). The record takes 85 and 7 ns.
  - A micro-bench interleaved with HEAD twice: Q 90–91 → 67 ns, and the angular intensity 15–16 → 21 ns.
  - The phase function can no longer be inlined into its caller: its edge needs a branch, and `#inline` ignores a body
    with control flow. Its 2 to 3 ns more are that call.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,349 ok, 0 failed, 0
  measured, 0 pending, 76,900 accepted (ADR 01: 29,441; ADR 0002: 47,459)** over 570,249 rows; **all 380 digests
  equal**. The four Rayleigh functions' digests changed; the other 376 are unchanged. `scattering.tcyr` passes on all
  three. The gate is clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 71. M3 part 2d-2b: Mie's tiny m (2026-11-14)

- **The user's answer of 2026-11-14: fix it in M3.** §70 found that Rust's `mie` refuses an index |m| below about
  10⁻¹⁵⁵ at ordinary x, though the m → 0 limit is finite and Rust's own documentation says the series fails "only for
  vanishingly small x". Reading the chain closer found more:
  - mx underflows to 0 for a subnormal m and Rust divides by it (the old `cyrius_mie_errors` assertion);
  - just above the refusal, the series is silently wrong. At x = 1, m = 10⁻¹⁵³, Dₙ/m is finite but (Dₙ/m)·χ₄ is not,
    Smith's division by an infinite denominator gives a₄ = 0, and Q_sca comes out 0.2768511768752622 for
    0.2768511783189432.
- **ADR 0002 row 18, amended for Mie's m-dependent chain.**
  - **The unbounded path** (`_bj_sc_mie_ux`) runs every m-dependent operation to 53 bits with an unbounded exponent:
    m·x, 1/(mx), the Dₙ recurrence, 1/m, both factors, the coefficients' numerators, denominators and Smith
    divisions, and the three sums. Each coefficient part and each output rounds once. complex_recip takes its first
    branch, one division by |z|²: its Smith fallback only avoided the range's ends. The x-only ψₙ, χₙ, n/x, 2/x² and
    the n-only weights stay Rust's doubles.
  - **The cues**, read off Rust's fast path:
    - a nonzero part of m outside [2^-511, 2^511), or |m|² or |mx|² not normal (checked before the fast path, which
      is then skipped);
    - a reciprocal in the Dₙ recurrence falling back to Smith's division;
    - a denominator of aₙ or bₙ not finite: an infinite factor (Dₙ/m + n/x, m·Dₙ + n/x) or a numerator past DBL_MAX
      makes its denominator non-finite too, and for x ≥ 10⁻⁵⁰ the x-only values are always finite. A first cut also
      tested the factors, the numerators and the x-only values each step; dropping those took the cues' cost from
      11–15 % of a Mie call to about 4 %;
    - a coefficient part outside [2^-511, 2^511): its squares and products in the sums may leave the range (g's
      aₙ·bₙ underflows below x ≈ 10⁻³⁸).
  - **The boundary.** Below x = 10⁻⁵⁰ no cue counts and Rust's chain runs to its end, refusals included: the tiny-x
    limit the user kept, made exact here. The double 1e-50 is 10⁻⁵⁰'s nearest double above it, with none between,
    so `x >= 1e-50` is the exact test. At and above 10⁻⁵⁰ every x-only value is finite (χ₄ ≤ about 10²⁰²), so a
    flagged overflow is m-driven. A first cut flagged x-driven overflows at x = 10⁻¹⁰⁰ and would have changed
    Rust's tiny-x results; the boundary is what keeps them.
  - **A bug found by the boundary rows.** Below 10⁻⁵⁰ the Dₙ loop's fallback still stopped the series from running,
    and the flag was then cleared, so the port returned sums it had never computed (4 failing rows). The series now
    skips only when the unbounded path will run.
- **The model** (`validate.py`): `sc_mie_ux` is the unbounded path in exact Fractions with Rust's x-only doubles, and
  `sc_mie(r18)` mirrors the cues on Rust's doubles and takes it where they fire. `validate_scattering` now proves the
  cues: on every Mie row with x ≥ 10⁻⁵⁰ that they keep on Rust's chain and whose recurrence starts at most 300 terms up
  (`SC_UX_PROOF`), the unbounded path must give Rust's reading, bit for bit. A one-off run without the 300 cap checked
  all 957 such rows, 0 different. It first found two: x = 10⁻⁵⁰ with an ordinary m, where g's aₙ·bₙ underflowed. That
  gave the coefficient cue. Mie rows go to the mapping file or the corrected oracle, never both.
- **The oracle**: `scattering.tsv` gains 87 rows after every earlier row (`mie_tiny_m_edges`; the old file an exact
  prefix):
  - 13 indices at x = 10⁻¹⁰, 0.1, 1, 5 and 20: from 10⁻¹⁵⁰ down to 5·10⁻³²⁴, the silent zone (2·10⁻¹⁵⁴ to
    1.6·10⁻¹⁵³), tiny imaginary parts (1.33 − 10⁻³⁰⁰ i, 0 − 10⁻²⁰⁰ i);
  - 12 rows at x = 0.1, 0.01 and 0.001 where (Dₙ/m)·χ₄ passes DBL_MAX while the Dₙ recurrence stays in range, so
    only the denominator cue sees them (at x = 1 the Dₙ recurrence overflows too);
  - x = 10⁻⁵⁰, the double below it, 10⁻⁴⁰ and 10⁻⁶⁰ with tiny and ordinary indices.
- **Mappings**: `scattering.adr.tsv` holds 495 rows: row 18's 394 (Rayleigh's 325 and Mie's 69), row 13's 93, row
  17's 8. Of the 69 Mie rows Rust refused 40 ("over/underflowed"), divided by zero on 4, and computed 25 differently:
  the silent zone at every x, and the tiny-x underflow of g at x = 10⁻⁵⁰.
- **Mutations**, each run through the parity harness. Six of ten were caught:
  - the coefficient cue off: 3 failing rows; the input cues off: 4; the boundary moved: 4;
  - the denominator cue off: 11, once the rows only it sees were added (a first try at disabling it tested fd != fd,
    which an infinite denominator's NaN still trips);
  - g without its factor 2: 77; complex_recip's (−im) dropped: 21.

  The four others:
  - the Dₙ recurrence's fallback cue and m's box cue are redundant in practice: every row they would catch also
    trips |mx|², an overflow or the coefficient cue. They stay as cheap insurance;
  - Smith's branch differs from Rust's only at a tie |b.re| = |b.im|;
  - q_abs from the rounded doubles equals the unbounded difference above 10⁻⁵⁰, where Q is normal.
- **Tests**:
  - `cyrius_mie_errors`: m = 10⁻³²⁰ at x = 10⁻¹⁰ now computes, Q_sca at its Rayleigh limit (8/3)x⁴|K|² with K → −1/2,
    6.67·10⁻⁴¹ (Rust divided by zero); the leak check's refusal is x = 10⁻¹⁰⁰'s;
  - `cyrius_mie_row18` (4 assertions, 116 in the file): x = 1 at m = 10⁻²⁰⁰ (Rust refused) and at m = 10⁻¹⁵³ (Rust's
    silent a₄ = 0), x = 10⁻⁵⁰ at m = 10⁻²⁰⁰ computed, and the double below 10⁻⁵⁰ still refused.
- **Benchmarks.** The cues cost each series step a few compares: four denominators' finiteness and the coefficients'
  exponent fields.
  - A full bench interleaved against 2d-2's tree (three rounds, load 1.6 to 2.2): `scattering/mie_x1_real` 1.92 →
    2.06 µs (+7 %), `mie_x10_real` 3.78 → 4.11 µs (+9 %), `mie_x50_m3_5` 12.07 → 12.86 µs (+7 %); every other row
    within the load's noise (median ratio 1.00 over 81). The record takes the three new medians.
  - A first form, which also tested the factors, the numerators and the x-only values each step, cost 11 to 15 %.
  - The unbounded path itself is about 8 times the fast one on the rows that take it: 15.8 µs at x = 1,
    m = 10⁻²⁰⁰, for the fast path's 2.0 µs at m = 1.5.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,367 ok, 0 failed, 0
  measured, 0 pending, 76,969 accepted (ADR 01: 29,441; ADR 0002: 47,528)** over 570,336 rows; **all 380 digests
  equal**. Mie's digest changed; the other 379 are unchanged. `scattering.tcyr` passes on all three. The gate is clean
  (17 steps, no warnings; execution coverage 288 of 288 functions).

## 72. M3 part 2d-3a: the beam's, ABCD and resonator chains at their range edges (2026-11-14)

- **2d-3 runs in three bites.** beam's chains are more than one bite carries: 2d-3a (this one) the Gaussian beam's,
  the ABCD matrices' and the resonators' chains; 2d-3b the HG and LG amplitudes, whose Hermite and Laguerre recurrences,
  envelope and power over- and underflow where the mode's value is representable; 2d-3c the HG field's phase
  kz + k(x² + y²)/(2R) − (m + n + 1)ζ, which passes 2^2054 (kz alone reaches 2^2100 at λ = 5·10⁻³²⁴), where the 2,688
  bits of 2/π that `_bj_ux_reduce` reads end.
- **ADR 0002 row 18, amended** (`src/crmath_vx.cyr`'s `_bj_vx_bm_*`). Each function keeps Rust's chain and tests what
  it computed; where an intermediate left the normal range its chain runs with an unbounded exponent and the last
  operation rounds once. A non-finite input keeps Rust's chain.
  - **The spot size and the radius of curvature**, w₀·√(1 + (z/z_R)²) and z·(1 + (z_R/z)²): a square below 2^-1022 is
    rounded away against 1 by both chains, so they differ only where the square passes DBL_MAX, and Rust's result is
    then ±∞ or NaN. The test is one compare on the result. Rust: w = ∞ for 10²⁰⁴ at w₀ = 10⁻⁶, z = 10²⁰⁰,
    z_R = 10⁻¹⁰, and NaN for 0 at w₀ = 0; R = ∞ for 10³⁰⁰ at z = 10⁻²⁰, z_R = 10¹⁴⁰.
  - **The Gouy phase** arctan(z/z_R), where |z/z_R| is at most 2^-1022 from a nonzero z: the arctangent of the 53-bit
    quotient, rounded once. arctan t lies within |t|³ inside t, so it is t's rounding, with a midpoint of the subnormal
    grid broken toward 0 (a step of 2^-100 of t's significand into `_bj_k_ux_round`). Rust: ζ = 2^-1073 at
    z/z_R = 3·2^-1075 (2^-1074). Past DBL_MAX both chains give ±RN(π/2).
  - **The peak intensity factor** 2/(π·w·w), where π·w·w is not in (2^-1022, DBL_MAX], from the 53-bit w. Rust: 0
    at w = 10¹⁵⁵ (6.37·10⁻³¹¹); and where π·w·w rounds to 2^-1022 from the 53-bit 2^-1022 − 2^-1075, 2^1023 for
    2^1023 + an ulp.
  - **The curved interface's** C = (n₁ − n₂)/(n₂R), where n₁ − n₂ or n₂R passes DBL_MAX (|n₂R| ≥ 10⁻⁶⁰ is never
    below 2^-1022, and a difference below it is exact). Rust: −∞ for −2 at n₁ = 10³⁰⁸, n₂ = −10³⁰⁸, R = 1; 0 for
    10⁻¹⁰⁰ at n₁ = 10³⁰⁰, n₂ = R = 10²⁰⁰.
  - **`apply_ray`, `compose` and `determinant`**: each two-term sum outside [2^-900, DBL_MAX] that is not a sum of two
    zero products (a zero factor in each, common in ABCD products: a lens's B, free space's C) runs unbounded:
    `_bj_vx_bm_dot2` for a ray, beam.cyr's own tests and then `_bj_vx_bm_dot2_ux` for a product of matrices,
    `_bj_vx_cross1_f64` for the determinant (all three `_bj_vx_cross1_f64`'s result, bit for bit). Rust: AD − BC =
    ∞ − ∞ = NaN for 0 with every entry 10²⁰⁰; ∞ for 2^971 where a product past DBL_MAX meets −DBL_MAX; 0 for 2^-1074
    where (5 ± 7/16)·2^-1074, two subnormal products, round alike. A lone product below 2^-1022 rounds twice there,
    to 53 bits and then onto the subnormal grid, as in every row-18 chain: D = 0x0005555555555555 times θ = 0.1 gives
    …888 where Rust's single rounding gave …889.
  - **The transforms' w₀′** = √(z_R′λ′/π), where z_R′λ′ is not normal or its quotient by π falls below 2^-1022: the
    chain from the doubles z_R′ and λ′, each the exact value rounded once (rows 17 and 35) and the output beam's own,
    so its z_R stays consistent with them. Rust refused λ = z_R = 10²⁰⁰ through the identity (w₀′ = ∞, for
    5.64·10¹⁹⁹) and 10⁻²⁰⁰ (w₀′ = 0, for 5.64·10⁻²⁰¹).
  - **The stability's g₁g₂**, where a g is ±∞ (L/R past DBL_MAX): the 53-bit g's product, rounded once. g₁ and g₂ stay
    ±∞, their true values being unrepresentable, and `resonator_beam_waist`'s "unstable" reason quotes the new
    product. Rust: NaN for −0 at L = R₂ = 10³⁰⁰, R₁ = 10⁻¹⁰; −∞ for −1.1·10²⁹⁴ with R₂ an ulp above L.
  - **The mode spacing** (Δν_FSR/π)·arccos(±√(g₁g₂)), where Δν_FSR/π falls below 2^-1022 (Δν_FSR below π·2^-1022,
    7.0·10⁻³⁰⁸). The root's own underflow is invisible: arccos of anything below 2^-511 is RN(π/2), which lies more
    than 2^-55 from a midpoint. Rust: 0 for 5·10⁻³²⁴ at Δν_FSR = 5·10⁻³²⁴, g₁ = g₂ = 0.
- **The model** (`validate.py`'s `bm_eval(r18)`, with `bm_one_sq`, `bm_spot53`, `bm_gouy_ux`, `bm_dot2`, `bm_g53` and
  `bm_w0`): the unbounded chain itself on every row with finite inputs (and z_R ≠ 0 where it divides), so the
  harness checks the port's tests row by row. `beam_adr_rows` gives row 18 to every row whose reading row 18's chains
  change, whatever else changes it.
- **The oracle**: `beam.tsv` gains 936 rows after every earlier row (`chain_edges`; the old file an exact prefix):
  - ratios z/z_R around 2^512 and far past it, with zero, negative and subnormal waists; z_R/z past 2^512 where R is
    representable;
  - Gouy quotients below 2^-1022: the exact midpoints (2j + 1)·2^-1075 for seven j, either sign and their
    neighbours, and 2^-1022 itself;
  - waists around DBL_MAX's root and where π·w·w is subnormal, and the one found by a search where π·w·w rounds to
    2^-1022 from below;
  - curved interfaces past DBL_MAX; twelve matrices whose products pass DBL_MAX (cancelling or not) or fall below
    2^-1022, with five rays and three products each, and a product past DBL_MAX meeting −DBL_MAX and the two
    subnormal products that cancel, each through a product, a ray and a determinant;
  - seven beams whose z_R′λ′ leaves the range, through the identity, free space, a lens and two scaled matrices;
  - cavities with L = 10³⁰⁰ and R₁ = ±10⁻¹⁰ or 10⁻⁸ against seven R₂;
  - mode spacings for nine Δν_FSR around and below π·2^-1022, and four g pairs whose product is below 2^-1022.
- **Mappings**: `beam.adr.tsv` holds 6,152 rows, row 18's 1,260: `spot_size` 589, `peak_intensity_factor` 314, each
  transform 137, `compose` 19, `gouy_phase` 13, `transverse_mode_spacing` 13, `apply_ray` 11, `curved_interface` 9,
  `resonator_stability` 6, `determinant` 5, `resonator_beam_waist` 4, `radius_of_curvature` 3. 1,068 of them are older
  rows: the earlier grids' beams at the range edges (z_R near 10⁻³¹⁴ with z = 1 gives Rust's w = ∞ for 3.4·10¹⁵³), 141
  of them another row's mapping before. All passed on the harness's first run.
- **cr_verify** (`bspot`, `bpeak`, `broc`, `bgouy`, `bci`, `bdot2`, `bw0`, `bgp`, `bfsr`; `cr_driver.cyr` copies
  beam.cyr's few lines of test before each helper): 87,550 calls per target against exact Fractions, 0 mismatches on
  x86_64 and aarch64. The draws sit around the tests' bounds, so the calls kept on Rust's chain check that claim too.
- **Mutations**, each run through the parity harness: twenty of twenty-two caught. Three were caught only after
  rows were added for them: the peak factor's strict 2^-1022 bound (the searched waist), `compose`'s test against
  DBL_MAX (an ∞ entry passes the lower bounds) and its zero test (the cancelling subnormal products). The two left:
  - the z_R = 0 test in front of the unbounded chain: x/0 gives the same ±∞ or NaN in both chains, so it is cheap
    insurance;
  - the mode spacing's strict bound: π·2^-1022 is the only Δν_FSR whose quotient by π rounds to 2^-1022, and that
    quotient is exact.
- **Tests**: `cyrius_beam_row18` (17 assertions, 152 in tests/beam.tcyr): the examples above.
- **Benchmarks.** A full bench interleaved against 2d-2b's tree, with a one-off bench of the guarded functions the
  record has no twin for (three rounds, load 3.4 to 6.9; each figure the lowest of the three):
  - `beam/spot_size` 11 → 13 ns, `beam/abcd_compose` 20 → 31 ns, `beam/resonator_stability` 42 → 47 ns; the HG
    amplitude and field, `propagate_beam` and every other row within the load's noise (median ratio 1.00 over 89).
    The record takes the three new medians (13, 32 and 47 ns);
  - one-off: `radius_of_curvature` 8 → 12 ns, `determinant` 5 → 7, `curved_interface` 14 → 17,
    `peak_intensity_factor` 21 → 23; `gouy_phase`, `transform_beam` and the mode spacing unchanged; `apply_ray` 8 →
    46 ns on a ray the lens maps to an exact 0 by cancellation (−10·0.001 + 0.01), which takes the exactness test,
    and 9 → 14 ns otherwise (a focused probe);
  - a finding, recorded: in this compiler a band test costs about what the operations it guards do. `compose`'s
    four magnitudes, their sum and five compares cost 11 ns over 20, and an executed five-argument helper call
    about 8 ns. A first form called `_bj_vx_cross1_f64` per entry (`compose` 39 ns, a lens composed with itself
    65 ns for a structural zero); testing structural zeros in the rare path took the lens to 54 ns. Moving the tests
    inline into `compose`, or flagging entries for a single call, saved nothing: the dead code's locals cost the
    common case more.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,179 ok, 0 failed, 0
  measured, 0 pending, 78,093 accepted (ADR 01: 29,441; ADR 0002: 48,652)** over 571,272 rows; **all 380 digests
  equal**. The thirteen functions this bite guards changed digests; the other 367 are unchanged. `beam.tcyr` passes
  on all three. The gate is clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 73. M3 part 2d-3b-1: the Hermite–Gaussian amplitude at its range edges (2026-11-14)

- **2d-3b runs in two bites.** The HG amplitude needs an e^x whose argument reaches far past 2300, so it goes first
  (2d-3b-1); the LG intensity, its Laguerre recurrence and its power follow (2d-3b-2). The field's chain (its 53-bit
  amplitude and its phase) is 2d-3c; here the field only takes the new amplitude, rounded as before.
- **ADR 0002 row 18, amended.** Rust's ((w₀/w)·H_m(√2x/w))·H_n(√2y/w)·e^{−(x² + y²)/(w·w)}:
  - ∞ for 5.3·10³⁰⁶ at m = n = 150, √2x/w = √2y/w = 2: H_m·H_n passes DBL_MAX, while Cramér's bound
    |H_m(s)|·e^{−s²/2} ≤ 1.09·√(2^m·m!) keeps the mode finite;
  - ∞·0 = NaN for 1.85·10³⁰ at m = 200, √2x/w = 40, and for 2.55·10²⁶ at m = 4000, √2x/w = 220 (envelope e^-24200);
  - NaN for 4.87 at w₀ = x = y = 10¹⁶⁰ (w·w and x² + y² past DBL_MAX), and at 10⁻¹⁷⁰ (both 0);
  - −0 for 5.4·10⁻²⁵⁷ where (w₀/w)·H₁ fell below 2^-1022 before H₁₅₀ lifted it back;
  - −1.98·10⁻³¹⁹ for −2.79·10⁻³¹⁹ at x = 5·10⁻³²⁴, √2·x rounded onto the subnormal grid.
- **The guard** (beam.cyr): Rust's chain where it is provably the unbounded one, else `_bj_vx_bm_hg`:
  - each Hermite value from H₂ on at least 2^-900 for a nonzero argument (each step two products and a difference:
    the band argument), or exact integers and zeros on the axis; an ∞ or NaN reaches the partial products;
  - for an order from 1, √2x and √2x/w normal (tested inside the Hermite helper, so the TEM₀₀ case skips them);
  - w·w and the envelope normal, which makes w₀/w normal too (w = w₀·√(1 + (z/z_R)²) with the root below 2^512): the
    first draft's own test of w₀/w was implied, and so was a test of H's finiteness, and both went;
  - x² + y² in [2^-900, DBL_MAX], or x = y = 0;
  - both partial products normal or exact zeros.
- **The unbounded chain** (`src/crmath_vx.cyr`): w, √2x/w and √2y/w, both recurrences (`_bj_vx_bm_hermite_ux`),
  x² + y², w·w and their quotient q, w₀/w and the partial products, each to 53 bits; e^{−q} at 53 bits (1 below
  2^-60, which is under 2^-54); the last product rounded once.
- **e^x to |x| = 2^47** (`src/crmath_ux.cyr`'s `_bj_ux_exp_wide`). `_bj_ux_exp_lim` stopped at 2300. A product with the
  envelope can stay representable far below e^-2300 (the m = 4000 row), so the reduction is stretched:
  - the fixed point's integer part sits unmasked in limb 0, whose add, subtract and compare already carry into it;
  - k·ln 2 is formed by halves of 24 bits, so no limb product passes 2^52;
  - the error stays below 2^48·2^-364 = 2^-316, inside the slow path's 2^-290;
  - beyond 2^47 the stand-in (1, ±2^62) is exact: for m, n < 2^32 and |√2x/w| ≤ 2^2100 the factors before the
    envelope stay below 2^(2^44.1), so every product with e^{−2^47} rounds to 0.

  `mpfr_ref.py` gained `wide_split`, the widest-range value as (s, e): e^(−10¹⁴)'s Fraction would take 10¹⁴ bits.
- **The model** (`validate.py`'s `bm_hermite53` and `bm_hg_ux`): exact Fractions, the envelope from MPFR at 53 bits
  with its widest exponent range. A product whose binade lies more than 1200 below e^q's is ±0 without building the
  envelope. It is the unbounded chain on every row with finite inputs, so the harness checks the guard row by row.
- **The oracle**: `beam.tsv` gains 76 rows after every earlier row (`hg_edges`, each through `hermite_gaussian` and
  `hermite_gaussian_field`):
  - H_m·H_n past DBL_MAX at m = n = 150 and five more order pairs, on both sides of the range;
  - large H against envelopes below 2^-1022 and past e^-2300, to m = 4000;
  - √2x below 2^-1022, √2x/w below it, w·w past either end, the far field (w past DBL_MAX), the partial product below
    2^-1022 lifted back;
  - on-axis exact zeros, and envelopes past the stand-in's 2^47;
  - one row for each of beam.cyr's tests failing alone, the last two found by a search (subnormal roundings show only
    for some significands: with x = (1 + k/512)·2^-512, x² lands exactly on the subnormal grid).
- **Mappings**: `beam.adr.tsv` holds 6,225 rows, row 18's 1,333: 73 more, the amplitude 40 and the field 33 (12 of
  them older rows). All passed on the harness's first run.
- **cr_verify**: `ux_exp_big` (e^x to 2^47 by its own reduction, every binade from 2^11, integers and halves, around
  2300 and 2^47, and the stand-in; 10,000 calls per target) and `bhg` (`_bj_vx_bm_hg` itself on beams from 2^±600,
  orders to 400 and envelopes past the stand-in; 5,000 calls, about two thirds normal results): 0 mismatches on x86_64
  and aarch64.
- **Mutations**: eleven of fourteen caught. Six of them (the w·w, envelope, √2x, √2x/w and x² + y² tests, and the
  envelope's 2^-60 cutoff) were caught only once rows were added that fail that test alone. The three left:
  - the per-step Hermite test: for a tiny normal argument Rust's odd values are exact multiples of s and its even ones
    absorb s², so it fires without changing a result. It stays as the proof's;
  - the on-axis exemption, which only spares those points the unbounded path;
  - the zero waist's return: both chains give NaN.
- **Tests**: `cyrius_hg_row18` (7 assertions, 159 in tests/beam.tcyr).
- **Benchmarks.** A full bench interleaved against 2d-3a's tree with a one-off bench of the amplitude (three rounds,
  load 3.2 to 7.5, other repos' fuzzing, tests and benches beside it; each figure the lowest of three):
  - `beam/hermite_gaussian_hg33` 169 → 199 ns (+18 %); the record takes the new median, 204 ns, and the field's,
    753 ns (its lowest 706 → 716). Every other row within the load's noise (each round's median ratio 0.99 to 1.01).
  - one-off: TEM₀₀ on the axis 44 → 69 ns, HG₃₁ on the x axis 154 → 191, HG₁₀,₁₀ 226 → 277.
  - a finding, recorded: the tests cost about 30 ns. Removing them all left the amplitude within 3 ns of HEAD, and
    the per-step check alone costs about 0.6 ns a step. The partial products' tests are the dearest (about 10 ns),
    then x² + y², the axis tests, w·w and the envelope. A second return value for the Hermite flag cost more than
    passing it through memory.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,182 ok, 0 failed, 0
  measured, 0 pending, 78,166 accepted (ADR 01: 29,441; ADR 0002: 48,725)** over 571,348 rows; **all 380 digests
  equal**. The amplitude's and the field's changed; the other 378 are unchanged. `beam.tcyr` passes on all three. The
  gate is clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 74. M3 part 2d-3b-2: the Laguerre–Gaussian intensity at its range edges (2026-11-14)

- **ADR 0002 row 18, amended.** Rust's ((((w₀·w₀)/(w·w))·ρ²^|l|)·L·L)·e^{−2r·r/(w·w)}, ρ = √2r/w, L = L_p^|l|(ρ²):
  - NaN for 1.0000047 at l = 2³¹ − 1 where ρ²^|l| and the envelope balance (ρ² = |l|·ln ρ² ≈ 5.3·10¹⁰: the power
    past DBL_MAX, the envelope e^-5.3·10¹⁰ far below 2^-1074); 1.0000011 at l = i32::MIN, its power 2³¹ by row 32;
  - 0 for 1.08·10⁻¹⁶⁰ at p = 100, ρ² = 1000 (e^-1000 underflows, L² does not);
  - ∞ for 3.9·10²¹⁶ at l = 150, ρ² = 400 (400^150 overflows); NaN for 122297 at p = 500, ρ² = 1800;
  - 0 for 3.6·10⁶⁶ at p = 1000, l = 200, ρ² = 0.01 (ρ²^200 underflows to 0 against an L near 10²³³);
  - NaN for 4.33 at w₀ = r = 10¹⁶⁰ (w₀·w₀, w·w and r·r past DBL_MAX).
- **The guard** (beam.cyr): Rust's chain where it is provably the unbounded one, else `_bj_vx_bm_lg`:
  - each Laguerre value from L₂ on at least 2^-900: each step is two products and a difference, then a quotient by
    k ≥ 2, and L₁ and every coefficient (2k + |l| − 1) − ρ² is exact or normal for a normal ρ²;
  - for r ≠ 0, (−2r)·r normal, and ρ²^|l| normal for |l| ≥ 1: its squares fall below 2^-1022 only below 1, where the
    power falls too;
  - w₀·w₀, its quotient by w·w, and the envelope normal;
  - the first and last partial products normal or exact zeros.

  A first draft also tested w·w, √2r, ρ, ρ² and the middle product. Each is implied by the rest:
  - w ≥ w₀, so w·w past DBL_MAX makes w₀²/w² zero, and below 2^-1022 makes w₀·w₀ subnormal too;
  - a normal ρ² means a normal ρ, which comes from a normal √2r unless w·w is below 2^-1022;
  - a ρ² below 2^-1022 vanishes against every integer in L, and the power's test sees it there;
  - the middle product past DBL_MAX makes the last ∞ too; below 2^-1022 it makes the last smaller still, where both
    chains give the same signed 0.

  Those tests went.
- **The unbounded chain** (`src/crmath_vx.cyr`): w, ρ and ρ², the recurrence (`_bj_vx_bm_laguerre_ux`, its integer
  parts exact as in Rust), (2r)·r, w·w and their quotient q, w₀·w₀ and its quotient by w·w, ρ²^|l| in compiler-rt's
  __powidf2 order (`_bj_vx_bm_powi_ux`, the order `_bj_powi` already shares with Rust) and the partial products, each
  to 53 bits; e^{−q} from §73's `_bj_ux_exp_wide`; the last product rounded once. The stand-in past 2^47 is exact here
  too: for p < 2^32 and |l| ≤ 2^31, L·L and ρ²^|l| stay below 2^(2^45.4).
- **The model** (`validate.py`'s `bm_laguerre53`, `bm_me_mul`, `bm_powi_me` and `bm_lg_ux`): the power and the products
  run on (m, e) triples with the same 53-bit roundings, since ρ²^(2³¹) at ρ² = 2·10⁸ is a Fraction of 6·10¹⁰ bits; a
  result past ±1200 in binade is ±∞ or ±0. It is the unbounded chain on every row with finite inputs, so the harness
  checks the guard row by row.
- **The oracle**: `beam.tsv` gains 29 rows after every earlier row (`lg_edges`):
  - azimuthal orders ±(2³¹ − 1) and i32::MIN at the balance radius and two neighbours, ±1000 at theirs (the envelope
    e^-9118, past e^-2300);
  - radial orders to 1000 whose L·L or ρ²^|l| passes DBL_MAX or meets an envelope below 2^-1022;
  - w₀·w₀, w·w and r·r past either end (one with q = 2^-45 on the unbounded path), the axis, and the far stand-in;
  - one row for each of beam.cyr's seven remaining tests failing alone, found by a search that mirrors the guard on
    Rust's doubles. With l = 0 the last product never fails alone: |L_p(x)| ≤ e^{x/2}, so L² passes DBL_MAX only
    past x = 709, where the envelope's test fails too. The row uses l = 30, p = 200, ρ² = 600: Rust's ∞ for 3.2·10⁶⁶.
- **Mappings**: `beam.adr.tsv` holds 6,284 rows, row 18's 1,422: 89 more LG rows, 62 of them older rows. 30 of those
  were row 32's before: at l = i32::MIN the power 2³¹ now runs unbounded wherever it leaves the range. Row 32 keeps 336.
  All passed on the harness's first run.
- **cr_verify**: `blg` (`_bj_vx_bm_lg` itself on beams from 2^±600, radial orders to 300, powers to 2³¹, ρ² from
  2^-120 to 2^80; 3,000 calls, 56 % normal results): 0 mismatches on x86_64 and aarch64.
- **Mutations**: thirteen of fifteen caught, the envelope's 2^-60 cutoff after its row was added. The two left are the
  per-step Laguerre test and the zero waist's return, as for HG (§73).
- **Tests**: `cyrius_lg_row18` (8 assertions, 167 in tests/beam.tcyr).
- **Benchmarks.** Rust has no LG bench, so the record has no row. A one-off bench interleaved against 2d-3b-1's tree
  (three rounds, load 4.0 to 4.8):
  - LG₁,₂ 158 → 185 ns (+17 %), LG₅,₋₃ 202 → 238 (+18 %), LG₀₀ on the axis 51 → 73;
  - as for HG: the tests on the chain's values are what cost.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,152 ok, 0 failed, 0
  measured, 0 pending, 78,225 accepted (ADR 01: 29,441; ADR 0002: 48,784)** over 571,377 rows; **all 380 digests
  equal**. The LG intensity's changed; the other 379 are unchanged. `beam.tcyr` passes on all three. The gate is
  clean (17 steps, no warnings; execution coverage 288 of 288 functions).

## 75. M3 part 2d-3c: the Hermite–Gaussian field at its range edges; M3 done (2026-11-14)

- **ADR 0002 row 18, amended.** Rust's field (u·cos φ, u·sin φ), with φ = k·z + ((k/2)·(x² + y²))·z/(z² + z_R²) −
  (m + n + 1)·ζ, k = 2π/λ and ζ = arctan(z/z_R); the port gives the conjugate:
  - NaN where φ passes DBL_MAX: at λ = 10⁻⁹, z = 10³⁰⁰, k·z = 6.3·10³⁰⁹, where the field is (3.5·10⁻¹¹, −9.2·10⁻¹¹);
  - NaN where k does (λ = 10⁻³¹⁰), for (0.078, −1.5·10⁻⁵);
  - the curvature term lost where z² passes DBL_MAX (z = 10¹⁶⁰: Rust's 1/R = 0 for 10⁻¹⁶⁰, a field of 9.3·10⁻¹⁶¹
    given as 0);
  - subnormal intermediates (k/2, z² + z_R², 1/R, ζ, the amplitude) rounded before use.
- **The unbounded chain** (`src/crmath_vx.cyr`'s `_bj_vx_bm_hgf`):
  - the amplitude's 53-bit value (`_bj_vx_bm_hg_ux`, split out of the HG chain), not its binary64 rounding;
  - k, every term of φ, and ζ to 53 bits (`_bj_vx_bm_atan53`: t itself below 2^-1022, ±RN(π/2) past DBL_MAX);
  - sin and cos of the 53-bit φ (`_bj_vx_sin53` and `_bj_vx_cos53`, Payne–Hanek past DBL_MAX);
  - the two last products rounded once.

  |φ| stays below 2^4199, since k ≤ 2^1076.7, x² + y² ≤ 2^2049 and z/(z² + z_R²) ≤ 2^1074 (1/(2|z_R|) at most, 1/z when
  z_R = 0).
- **The longer 2/π.** `src/crmath_ux_data.cyr`'s table grows from 96 chunks of 28 bits to 176 (from π to 5,100 bits,
  checked against MPFR), and `_bj_ux_reduce` takes angles below 2^4454, its reads ending at chunk 173. A first form
  kept M1's loader, one `store64` per chunk. Each such statement cost a test binary 64 bytes of static data, so the
  table took two of them past the compiler's 128 KiB advisory and `cyrius test` warned. Bisecting by reverting one file
  at a time found the cause. The table is now an initialized `i64` array, and the binaries are 4,752 bytes below
  HEAD's.
- **The guard** (beam.cyr): Rust's chain where it is provably the unbounded one:
  - the amplitude 0 or normal (a 0 that underflowed gives the same signed zeros);
  - k/2 normal;
  - for z ≠ 0: z² + z_R² in [2^-900, DBL_MAX], and k·z, 1/R and ζ normal;
  - for x or y ≠ 0: x² + y² in [2^-900, DBL_MAX], and (k/2)(x² + y²) and its product by 1/R normal;
  - φ finite.

  A sum of two doubles the chains share rounds alike in both (below 2^-1022 it is exact), and cos and sin of a finite
  double are normal. Elsewhere the unbounded chain, which declines (keeping Rust's) for a non-finite input or a zero λ,
  z_R or w₀.
- **The model** (`validate.py`'s `bm_hg_amp53` and `bm_hgf_ux`): the amplitude to 53 bits, its binade past ∓1300
  deciding ±0 or ±∞ alone; φ in exact Fractions; ζ, cos φ and sin φ from MPFR at 53 bits with its widest exponent
  range. A first model gave ±∞ for a huge amplitude times sin 0. A finite amplitude times 0 is ±0, which the four
  failing rows (orders 200 to 1000 at z = 0) showed.
- **The oracle**: `beam.tsv` gains 20 field rows after every earlier row (`hgf_edges`):
  - k·z past DBL_MAX (λ = 10⁻⁹) and past 2^2054 (λ = 5·10⁻³²⁴), at z = ±10³⁰⁰;
  - φ near 2^4141 (z = z_R = 5·10⁻³²⁴, x = 10³⁰⁰), where the amplitude vanishes and the zeros' signs are cos φ's and
    sin φ's;
  - k past DBL_MAX, and z² past it;
  - seven rows that each fail one of beam.cyr's tests alone. A search found them by mirroring the guard on the port's
    doubles, against the chain the guard keeps (the port's amplitude with Rust's phase). A first search compared
    against Rust's whole chain, and its rows differed through the amplitude, not the test. The k/2 row needed its own
    search: k just below 2^-1021, x² + y² near DBL_MAX.
- **Mappings**: `beam.adr.tsv` holds 6,317 rows, row 18's 1,455: 66 field rows (33 more; 13 of the new ones older rows,
  and 3 older mappings with new values). All pass.
- **cr_verify**: the `ux_sin` and `ux_cos` families now draw angles to 2^4452 (about 9,750 calls each per target), and
  `bhgf` runs `_bj_vx_bm_hgf` itself (3,000 fields, both parts: 6,000 calls; 22 % of the phases past DBL_MAX, 3 %
  past 2^2054, 70 % of the parts nonzero). 0 mismatches on x86_64 and aarch64.
- **Mutations**: eleven of sixteen caught, four (k/2, z² + z_R², 1/R, ζ) once the searched rows were added. The five
  left:
  - k·z's and the curvature term's tests: their underflow moves φ by a few units of 2^-1075, visible only where φ
    cancels to near 2^-1074;
  - x² + y²'s and (k/2)(x² + y²)'s: their underflow shows only where z² + z_R²'s test fails too (60,000 draws found
    none alone);
  - a zero λ's return: k = ∞ gives NaN in both chains.
- **Tests**: `cyrius_hgf_row18` (5 assertions, 172 in tests/beam.tcyr).
- **Benchmarks.** A full bench interleaved against 2d-3b-2's tree with a one-off bench of the field (three rounds, load
  1.3 to 2.5; each figure the lowest of three):
  - `beam/hermite_gaussian_field_hg33` 713 → 734 ns (+3 %); the record takes the new median, 741 ns;
  - one-off: a field off the axis 693 → 723 ns, TEM₀₀ on the axis 126 → 140;
  - every other row unchanged (each round's median ratio 0.98 to 1.00).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **16 modules, 0 failing; 493,139 ok, 0 failed, 0
  measured, 0 pending, 78,258 accepted (ADR 01: 29,441; ADR 0002: 48,817)** over 571,397 rows; **all 380 digests
  equal**. The field's changed; the other 379 are unchanged, powf's and the Green's function's among them (they share
  the reduction). `beam.tcyr` passes on all three. The gate is clean (17 steps, no warnings; execution coverage 288 of
  288 functions).
- **M3 is done.** Part 1 ported polarization, scattering and beam in the one convention. Part 2, the range-edge
  review, took ADR 0002 rows 33–35 (2a), row 13 (2b), row 17 (2c) and row 18 (2d), in fifteen bites from 2a-1 to
  2d-3c.

## 76. M4 part 1a: rf's lines, reflection, Smith chart, matching and arrays (2026-11-15)

M4 ports rf.rs (1,378 lines with 80 tests; 2,345 with them) in Parts, as M3 did: Part 1 the API in bites (1a the
functions below; 1b the S-parameter matrices, their cascade and S → Z; 1c Touchstone read and write; 1d the JSON of
`SMatrix`, `TouchstoneFormat` and `TouchstoneData`), Part 2 the range-edge review.

- **`src/rf.cyr`**: rf.rs's free functions under Rust's names:
  - the lines: `bj_characteristic_impedance_lossless`, `bj_characteristic_impedance_lossy` (its DC branch: √(L′/C′)
    for R′ = 0, +∞ − j∞ for R′ > 0, NaN parts otherwise), `bj_propagation_constant`, `bj_phase_velocity` and
    `bj_input_impedance_lossless_complex` (into `out`);
  - reflection and the losses: `bj_reflection_coefficient`, `bj_reflection_coefficient_complex` (into `out`),
    `bj_vswr`, `bj_return_loss`, `bj_insertion_loss`, `bj_mismatch_loss`;
  - the Smith chart: `bj_normalize_impedance` (into `out`), `bj_denormalize_impedance`, `bj_gamma_from_z` and
    `bj_z_from_gamma` with their poles (+∞ + j0, and an infinite part to Γ = 1 or z = −1);
  - matching: `bj_quarter_wave_transformer`, `bj_single_stub_match` (into `out[0..2]`) and
    `bj_single_stub_match_complex` (both solutions into `out[0..4]`, sorted by d), with Rust's debug log line;
  - the arrays: `bj_linear_array_factor`, `bj_planar_array_factor`, `bj_circular_array_factor`.

  The deprecated `input_impedance_lossless` is not carried (ADR 0002 row 1).
- **No convention change**: rf.rs already used e^{+jωt} (+jωL, γ = jβ), so ADR 0002 row 5 changes no number, and the
  harness compares every row on the oracle's numbers.
- **Rust's chains, as written**:
  - Complex arithmetic is hisab's, written out on f64 locals in Rust's order: |z| = √(re² + im²) (hisab's `abs`, not
    hypot), arg z = atan2(im, re), from_polar as r·cos θ, r·sin θ with one sincos.
  - The division is rf.rs's Smith `cdiv`, the same as scattering.rs's (`_bj_vx_cdiv`); its `crecip` divides once when
    |z|² is a normal double. rf's `cdiv` returns NaN for a zero divisor where scattering's refuses, so each caller
    keeps its own test.
  - `complex_real`'s `+ 0.0` and `− 0.0` stay: `gamma_from_z` and `z_from_gamma` build z ± 1 from them, and the
    poles' tests and the signed zeros of Smith's quotient follow (z = 0 − j0 gives Γ = −1 + j0, as Rust's).
  - An element count is Rust's usize: the word read unsigned, `n as f64` rounded to nearest even above 2^63 (the
    halved word with a sticky bit, then doubled exactly). Rows at 2^63 + 1, + 0x400 and + 0xc00 (ties), + 0x401 and
    + 0xc01 (a sticky bit), usize::MAX − 0x3ff and usize::MAX pin it.
- **New golden file `rf.tsv`, 11,323 rows over 21 functions** (class C 6,940: every function that calls atan2, sin,
  cos, atan, log10 or acos; class A 4,383; 1,820 error rows). The grids:
  - each argument crossed with the specials; the Rust tests' cases; typical lines from DC to 10 GHz and the lossy
    impedance's DC branch over every sign of R′, L′, G′ and C′;
  - loads across the Smith chart and one ulp from its poles; every rejection;
  - each threshold straddled ulp by ulp: |Z_L + Z₀| < 10⁻³⁰, |Z_L + Z₀|² < 10⁻⁶⁰, |Z₀| < 10⁻³⁰, the unbounded input
    impedance's |den| ≤ 10⁻¹⁵·(Z₀ + |Z_L|) near βl = π/2, the stub match's |R_L − Z₀| ≤ 10⁻¹²·Z₀, the array factor's
    |sin(ψ/2)| < 10⁻¹⁵ at the main beam and the grating lobe;
  - element counts from 0 to usize::MAX, and seeded values over wide ranges (to 10^±300).

  `validate.py`'s independent transcription (`rf_eval`) agreed with every row on its first run, and a one-bit
  perturbation of every row's first output was caught (a NaN output matches any NaN, by class A's rule).
  **`rf.cr.tsv`, 62 corrected rows**: glibc's log10 misrounds in 43 loss rows (return, insertion, mismatch), its
  cos(π/6) in 12 linear-array rows, its sin or cos in 6 circular-array sums and its sincos in one input impedance.
- **Rust tests**: 45 of rf.rs's 80 are ported (`tests/rf.tcyr`, 90 assertions with the Cyrius-only groups: the
  reasons' exact text, out untouched by every refusal, the DC limits, the Smith chart's poles and signed zeros, the
  stub match's matched branch and its sort). The deprecated function's three are n/a (ADR 0002 row 1), and
  `cyrius_input_impedance_real_loads` checks their values through the complex form. The ledger: 567 ported, 6 n/a,
  285 pending of 858.
- **Surface parity**: rf is in progress, its 21 functions matched (420 items in all). `IN_PROGRESS` now names a
  struct or enum still to come (its fields or variants with it) and a method as `Owner::name`: `SMatrix`,
  `TouchstoneFormat`, `TouchstoneData`, their methods and serde impls, `s_to_z`, `parse_touchstone` and
  `try_write_touchstone` are pending. The deprecated `input_impedance_lossless` and `write_touchstone` are
  exceptions (ADR 0002 row 1). The API snapshot gains 21 public functions (2,054).
- **Benchmarks** (the 5 twins of part 1a, so 66 of 80 are twinned). An interleaved A/B against HEAD (three rounds,
  load 0.75 to 1.15) gave each twin's figure within 6 ns across the rounds, and left every other row unchanged
  (median ratio 1.00 over 81 rows; three 1 to 4 ns faster, the bimodal rows):

  | Bench | Rust | Cyrius | Ratio |
  |---|---|---|---|
  | `rf/reflection_coefficient` | 1.23 ns | 14 ns | 11× |
  | `rf/vswr` | 1.23 ns | 9 ns | 7.3× |
  | `rf/linear_array_factor_8` | 0.55 ns | 439 ns | 798× |
  | `rf/input_impedance_complex` | 7.43 ns | 221 ns | 30× |
  | `rf/single_stub_match` | 128.46 ns | 1.88 µs | 15× |

  The Rust inputs are constants, and LLVM folds the array factor's cos and two sins at compile time (0.55 ns is a
  load). A finding, with a one-off probe of the primitives on this tree (1M calls each, load about 1): atan 279 ns,
  atan2 339, sincos 189, sin 150, cos 124, log10 83. The stub match is two atans, two sincos and two atan2s, about
  1.6 µs of its 1.88; the input impedance is one sincos and two square roots. src/crmath.cyr's atan (as atan2(x, 1))
  and atan2 cost about twice its sin, so they are the first place a faster correctly rounded path would pay off here.
- **For M4's range-edge Part** (§0.1), candidates from writing the port:
  - the thresholds on rounded values: |Z_L + Z₀| < 10⁻³⁰ and |Z₀| < 10⁻³⁰ (an exact sum or input), |Z_L + Z₀|² <
    10⁻⁶⁰ (a rounded sum of squares, row 17's question), the unbounded input impedance's relative 10⁻¹⁵ test and the
    stub match's |R_L − Z₀| ≤ 10⁻¹²·Z₀ (tolerances Rust chose);
  - the complex chains past either end: |z|² in the lossy impedance, the propagation constant, the input impedance's
    scale and |den|, crecip's normal test; the stub match's (Z₀ − R_L)² + X_L²; `normalize_impedance`'s Z·(1/Z₀),
    rounded twice.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 504,400 ok, 0 failed, 0
  measured, 0 pending, 78,320 accepted (ADR 01: 29,503; ADR 0002: 48,817)** over 582,720 rows; **all 401 digests
  equal**, HEAD's 380 unchanged. `rf.tcyr` passes on all three. The gate is clean (17 steps, no warnings; execution
  coverage 309 of 309 functions).

## 77. M4 part 1b: the S-parameter matrices (2026-11-15)

- **`BjSMatrix`, fixed at construction (ADR 0002 row 36, the user's answer of 2026-11-15)**:
  - Rust's `SMatrix` had public `s: Vec<Vec<Complex>>` and `z0: Vec<f64>`. A caller could reshape them, so
    `check_shape`, `two_port_params` and `get`/`set` checked the shape on every call and could fail with "inconsistent
    SMatrix".
  - The port's matrix is `ports` plus one fl_alloc block: the ports² entries row-major (S_ij at (i − 1)·ports + j − 1,
    16 bytes each), then the ports reference impedances. `s` and `z0` are the two addresses, both 0 for a 0-port.
    Nothing can reshape it, so the shape errors are gone.
  - `get` and `set` keep Rust's 1-based ports and their reason ("port indices must be 1..={n}, got ({i}, {j})"),
    with usize indices: a word past 2^63 is refused and printed unsigned. `s_at(i, j)` (Rust's `s[i][j]`, 0-based)
    and `z0_at(k)` and `set_z0(k, v)` (Rust's `z0[k]`) are new accessors. `set_z0` refuses an index past the last,
    and leaves the value unchecked, as Rust's field is.
  - `BjSMatrix_new`, `BjSMatrix_two_port` and `cascade` return Allocation when their block can't be had (row 9,
    amended): Rust's `vec!` aborts. A refused `cascade` leaves `out` untouched and never reads it, and
    `BjSMatrix_free` gives the block back.
- **`cascade` and `bj_s_to_z`, in Rust's order**: the checks (2-ports, `uniform_z0` on each, the two equal), the
  transfer matrices with |S21|² < 10⁻⁶⁰ refused, T = T₁·T₂, T22's |T22|² < 10⁻⁶⁰ refused (a lossless loop,
  S22·S11′ = 1), and back to S; S → Z on Z₀(I + S)(I − S)⁻¹ with |det(I − S)|² < 10⁻⁶⁰ refused. hisab's products
  run on f64 locals; z0 times a matrix entry is the full complex product (z0 + j0)·w, as Rust's is, so a NaN or ∞ part
  meets 0·∞ exactly where Rust's did. The trace lines are Rust's ("cascading 2-port S-matrices z0=…").
- **ADR 0002 row 37, new**: rf.rs's `z0_equal(a, b)` is |a − b| ≤ 10⁻¹²·max(|a|, |b|), and for b = ±∞ that is
  ∞ ≤ ∞. So a 2-port whose second port's Z₀ was infinite passed `uniform_z0`, and its cascade and Z-parameters used
  the first port's 50 Ω. The port's `z0_equal` also requires both values finite, so such a matrix is refused with
  Rust's own reason ("all port reference impedances must be equal, got [50.0, inf]"). `rf.adr.tsv` maps the 12 rows
  (9 cascades, 3 conversions), all passing. The Touchstone writer (part 1c) uses the same test.
- **`rf.tsv` grows by 1,707 rows (13,030)**, appended after part 1a's (whose body is an exact prefix). A matrix is
  `i:<ports>`, `v:<2n²>` and its entries, `v:<n>` and its impedances, consistent shapes only; `set` rows hold its
  result and then the matrix after it. The grids:
  - `new` over 0 to 8 ports and eight impedances, `two_port` over the Rust tests' networks and the specials;
  - `get` and `set` on four matrices (0, 1, 2 and 3 ports) at every index from 0 to one past the last, 2^63 and
    usize::MAX;
  - the Rust tests' networks (the thru, the star product's pair, the attenuator, the shunt resistor, the bench's), and
    non-2-ports on either side;
  - every reference impedance rule: equal, mixed, within and just past z0_equal's 10⁻¹² on either side, zero,
    negative, NaN and ±∞ on either port;
  - |S21|² and |det(I − S)|² straddling 10⁻⁶⁰ ulp by ulp (det = k·2^-103 across 10⁻³⁰), and T22 = 0;
  - the specials in each of the eight entries, and 210 seeded networks, passive and over 10^±300.

  The transcription (`rf_sm_eval`) agreed with every row on its first run, and a one-bit (or one-integer) perturbation
  of every row's first output was caught.
- **Rust tests**: 12 more of rf.rs's are ported (57 of 80; `tests/rf.tcyr`, 132 assertions), and
  `test_smatrix_inconsistent_does_not_panic` is n/a (row 36). The Cyrius-only groups pin the layout, the accessors'
  refusals, Allocation at 2^40 ports, row 37, every reason's text, out untouched, and the heap (2,000 cascades, each
  freed, leave the arena as it was). The ledger: 579 ported, 7 n/a, 272 pending of 858.
- **Surface parity**: `SMatrix`, its fields and methods and `s_to_z` leave `IN_PROGRESS` (430 items matched). The API
  snapshot gains 4 public functions (2,058).
- **Benchmarks** (the `rf/s_cascade` twin, so 67 of 80 are twinned). An interleaved A/B against part 1a's tree (three
  rounds, load 0.8 to 1.3) gave 523 or 524 ns in each round, for Rust's 100.25 ns (5.2×), and left every other row
  unchanged (median ratio 1.00 over 86 rows; `relativity/em_tensor_from_fields` 31 → 34 ns, its known bimodal pair).
  A one-off probe split the cascade's 540 ns: building the result is 124 ns (`BjSMatrix_two_port` with four BjComplex
  by value, and freeing it; the block's `fl_calloc` and `fl_free` alone are 36 ns), and the arithmetic, 22 complex
  products through `_bj_rf_cmul` and three reciprocals, is the rest, at about 8 ns per executed helper call (§72).
  `bj_s_to_z` alone is 270 ns. The products written inline would trade the helper's call for
  locals; speed is measured, not gated, so this is a finding, not a change.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 506,095 ok, 0 failed, 0
  measured, 0 pending, 78,332 accepted (ADR 01: 29,503; ADR 0002: 48,829)** over 584,427 rows; **all 407 digests
  equal**, HEAD's 401 unchanged. `rf.tcyr` passes on all three. The gate is clean (17 steps, no warnings; execution
  coverage 313 of 313 functions).

## 78. M4 part 1c-1: reading Touchstone (2026-11-15)

Part 1c runs as two bites: 1c-1 the reader, 1c-2 the writer (Rust's `{:e}` form) with the round trips.

- **`src/touchstone.cyr`**: `bj_parse_touchstone(content, len, ports, out)` into a `BjTouchstoneData`, and
  `BjTouchstoneFormat` (MagnitudeAngle 0, RealImag 1, DbAngle 2).
  - **The data's shape is fixed** (ADR 0002 row 36, amended): one fl_alloc block of the `points` frequencies, the
    `points` matrix headers (each a `BjSMatrix`, so `get`, `cascade` and `bj_s_to_z` take it as it is) and their
    entries and impedances. `frequency_at`, `matrix_at` and `BjTouchstoneData_free` replace Rust's vectors.
  - **The reader buffers records** in a store that doubles (fl_alloc), then builds the block once at the end. A
    refused parse frees both and leaves `out` untouched.
  - **Rust's order and reasons**: the port count, the lines, the option line (only the first, before the data),
    v2 keywords, each token's grammar and finiteness, the frequency (finite, non-negative, through the unit), the
    2-port noise block, each pair's representability, a partial last record, no record at all. The three debug lines
    are Rust's.
- **Rust's text, measured first** (a probe on rustc 1.99.1, then the oracle):
  - `{:?}` of a str escapes a character the same wherever it sits: 954,590 scalar values, as `\u{…}` unless one of
    `\0 \t \n \r \" \\`. `char::is_whitespace` is White_Space's 25 code points.
  - `parse::<f64>` takes `1.`, `.5`, `+.5`, `1.e5` and an exponent with any number of leading zeros, and refuses `.`,
    `1e`, `e5`, `+-1`, `1_0`, `0x10`, `١` and a surrounding space. The stdlib's f64_parse_n gave Rust's value on every
    such form tried (and on probe p12's 120,018 strings): the reader checks Rust's grammar, then calls it. It is the
    bundle's first direct call into the stdlib's math, so `dist/bijli.deps` lists math, which `cyrius.cyml` declares.
  - **The tables are generated** by the oracle (`t_touchstone.rs`) into `src/touchstone_data.cyr`: 916 escape ranges
    (a binary search) and the 25 White_Space values. They are the compiler's own: rustc 1.99.1 has 910 ranges and
    1.100.0, the oracle's toolchain since 2026-11-13, has 916 (a Unicode update). Every earlier output is the same on
    both; this is the first that isn't. So `scripts/oracle/build.rs` records the compiler, the file's header names it,
    and `--check` fails if another toolchain moves the tables.
- **ADR 0002 row 38, new**: Rust's `&str` is UTF-8 by construction; the port takes bytes. Bytes that are not UTF-8
  (`str::from_utf8`'s rules: no overlong forms, surrogates or values past U+10FFFF) are refused with InvalidParameter
  "Touchstone content is not valid UTF-8 (an invalid sequence at byte {k})", after the port count. The oracle can't
  make such a row; `cyrius_touchstone_utf8` and the fuzz harness cover it.
- **`rf.tsv` grows by 417 rows (13,447)**, appended (the earlier body an exact prefix). A row's inputs are the content
  (`s:`) and the port count; its outputs the data (`i:<ports> f:<z0> v:<points>`, the frequencies, each matrix) or the
  error. Class C: MA and DB records take a sincos, DB a pow. The grids:
  - the Rust tests' files, the five-port wrapped file as Rust builds it, and the bench's 201-point file;
  - every option: the units, formats and parameter types in each case, R's values (empty, NaN, ∞, −0, 0, overflow,
    the forms of 50), S repeated, R twice, unknown tokens, a second option line, one after the data, the noise block
    or a partial record, v2 keywords after spaces and in comments;
  - 54 tokens through Rust's grammar, each as an S11 and as a frequency in Hz, then the decimals a correctly rounded
    parse needs: half the smallest subnormal written in full (a tie, to 0) and with one more digit (to 5·10⁻³²⁴), the
    DBL_MAX tie, 2^53 + 1 with long tails, 350 leading zeros against e360, 400 ones against e-400;
  - comments, blank lines, CRLF and a final bare '\r', and each White_Space character (with U+200B, U+FEFF and U+180E,
    which aren't) around and between tokens;
  - records across lines, partial ones, MA and DB pairs that overflow (dB = 6164 and 6165 straddle DBL_MAX) or
    underflow, angles of 10²²; negative, −0, zero and overflowing frequencies; decreasing ones (a 1-port, a 3-port);
  - the 2-port noise block: an equal frequency starts it, its tokens still parse, a record left open by it;
  - error snippets of 32 and 33 characters, of multi-byte characters, combining marks, quotes, NUL, controls,
    private-use and unassigned characters, emoji and regional indicators;
  - 80 seeded files (1 to 4 ports, every format and unit, values in Rust's Display, `{:e}` and Debug forms).

  The transcription (`rf_parse_touchstone`) restates Rust's text rules rather than borrowing Python's (`str.strip`
  treats U+001C–U+001F as space; `float()` takes `1_0` and Arabic digits), and escapes with the generated table. It
  agreed with every row on its first run; a perturbed first output was caught on every row, and a perturbed reason on
  all 152 error rows. `rf.cr.tsv` gains 3 rows where glibc's sin, cos or pow misrounds a pair.
- **Fuzzing** (`tests/bijli.fcyr`, the plan's §6 step 7 harness, until now a stub): six seed files mutated byte by
  byte (bytes of the grammar, controls, and ones that break UTF-8) and parsed as 1 to 4 ports, 20,000 mutants
  (about 1,070 accepted). Every parse returns, a refusal leaves `out` untouched, an accepted one has its ports, points
  and finite values; the same 20,000 replayed don't grow the freelist's arena or the bump heap.
- **Rust tests**: 12 more ported (69 of 80; `tests/rf.tcyr`, 182 assertions), the Cyrius-only groups pinning row 38,
  the text rules (White_Space, `{:?}`, the 32-character cut), every reason's text, the layout and the heap. The
  ledger: 591 ported, 7 n/a, 260 pending of 858.
- **Surface parity**: `TouchstoneFormat` (and its variants), `TouchstoneData` (and its fields) and
  `parse_touchstone` leave `IN_PROGRESS` (440 items matched); `try_write_touchstone` and the serde impls remain. The
  API snapshot gains 2 public functions (2,060).
- **Static data**: the escape table as an `i64[916]` cost every test binary 7.3 KB of static data and pushed
  `tests/crmath_kx.tcyr` (kimiya's generated vectors) 584 bytes past the compiler's 128 KiB advisory (131,656
  bytes). A probe showed a string literal isn't static data (12,000 characters of one added 16 bytes), so the table is
  one multi-line literal of fixed-width records (12 hex digits a range, nine to a line, read by offset in the binary
  search). The advisory is gone; the lookup runs only on a non-ASCII character in an error's quote.
- **Benchmarks** (the `rf/parse_touchstone_s2p_201` twin, so 68 of 80 are twinned). An interleaved A/B against part
  1b's tree (three rounds, load 0.2 to 1.1) gave 538, 538 and 546 µs, for Rust's 79.42 µs (6.8×), and left every other
  row unchanged (median ratio 1.00 over 87 rows). A probe split the first cut's 582 µs: the 804 sincos of the MA pairs
  are about 143 µs (the same file in RI took 438), and f64_parse_n on its 1,809 tokens about 99. The rest was the scan
  at about 42 ns a byte, a decode and a White_Space call per character; ASCII tested inline and one pass for '\n'
  and '!' brought it to 534 µs. What remains is the separate passes over each byte (UTF-8, the line, the token, the
  grammar) and a call per token: a finding, with those numbers.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 506,509 ok, 0 failed, 0
  measured, 0 pending, 78,335 accepted (ADR 01: 29,506; ADR 0002: 48,829)** over 584,844 rows; **all 408 digests
  equal**, HEAD's 407 unchanged. `rf.tcyr` passes on all three. The gate is clean (17 steps, no warnings; execution
  coverage 315 of 315 functions), and `cyrius fuzz` passes.

## 79. M4 part 1c-2: writing Touchstone (2026-11-15)

- **`src/fmtf.cyr` gains `_bj_fmt_f64_exp`**, Rust's `LowerExp` (`{x:e}`): the shortest digits Display prints (core's
  upper rounding of exact ties), always as `d.ddd` + `e` + the exponent with no `+` (`1e9`, `1.5e9`, `-2.5e-3`,
  `1.00000025e-3`), and `0e0`, `-0e0`, `inf`, `-inf`, `NaN`. A probe on rustc 1.100.0 fixed the form first.
  `fmtf.tsv` grows by 18,788 rows, the whole grid (every binade and power of ten with their neighbours, the 2^53
  edge, the 751 exact ties, seeded bit patterns) through `{x:e}`, appended after the fixed-precision rows; the
  transcription (`rust_exp`) agreed with each on its first run.
- **`bj_try_write_touchstone(data, format, sb)`**, Rust's `try_write_touchstone`:
  - the validation in Rust's order: the port count (1 to 256), R (finite, positive), then per point the frequency
    (finite, non-negative, strictly increasing), the matrix's port count, its impedances against R (row 37's
    `z0_equal`, now `public fn _bj_rf_z0_equal` beside `_bj_rf_require`, src/rf.cyr's helpers the writer shares) and
    its entries (finite). Rust's frequency and matrix counts can't disagree here (ADR 0002 row 36);
  - the text: the comment line, `# GHz S <fmt> R <z0>` (R in Display form), each record's frequency in GHz and its
    pairs as `  {v1:e}  {v2:e}`. 2-port records are `f N11 N21 N12 N22` on one line; for 3 ports and more each row
    starts a line and wraps after four pairs (`\n `). MA is (|S|, arg S in degrees), DB (20 log₁₀ max(|S|,
    MIN_POSITIVE), arg S), with hisab's |S| = √(re² + im²): an entry whose |S|² passes DBL_MAX can't be written in MA
    or DB, as in Rust (a range-edge candidate below);
  - the text is built in one fl_alloc buffer sized for the data (60 bytes a pair is more than the longest), then
    appended to the builder in one call, so a refusal (or a failed allocation) leaves the builder as it was.
- **`BjTouchstoneData_new(ports, z0, points, out)` and `set_frequency(k, f)`** build data to write, checking nothing
  (as Rust's struct literal), the writer checks it.
- **ADR 0002 row 37 in the writer**: Rust's `z0_equal` took a matrix impedance of ±∞ for the file's R and wrote the
  file; the port refuses it ("S-matrix #k reference impedances [50.0, inf] differ from the file's R 50"). `rf.adr.tsv`
  now maps 16 rows, the 4 writes among them (their port outputs taken on MPFR's primitives).
- **`rf.tsv` grows by 273 rows (13,720)**, appended. A row's inputs are the data, as `parse_touchstone`'s outputs, and
  the format (`i:` 0 MA, 1 RI, 2 DB, the port's enum); its output the text (`s:`) or the error. The grids:
  - the Rust tests' data (the one-port pair, the 2-port order, 3-, 4-, 5-, 8-, 9- and 16-port layouts, low
    frequencies and tiny entries, a zero entry in DB) in each format;
  - every rule of the validation: ports 0 and 257 (the latter with no points: a 257-port matrix would put 66,049
    entries in the row, 3 MB), R at 0, −0, negative, NaN, ∞, 5·10⁻³²⁴ and 10³⁰⁰ (long Display forms), frequencies
    at 0, −0, negative, NaN, ∞, 10³⁰⁰ and 5·10⁻³²⁴, equal and decreasing, a matrix's impedance at 75, ±∞, NaN, 0 and
    within and past 10⁻¹² of R, a non-finite entry;
  - entries at each format's edges: signed zeros in every combination and the four quadrants' angles, subnormals
    and DB's clamp at MIN_POSITIVE, |S|² just under and past DBL_MAX;
  - 60 seeded data sets (1 to 5 ports, 1 to 4 points, entries over 10^±12).

  The transcription (`rf_write_touchstone`) agreed with every row on its first run, and a perturbed first output was
  caught on each. `rf.cr.tsv` gains 21 rows where glibc's atan2 or log10 misrounds a written digit.
- **Rust tests**: the last 6 of rf.rs's are ported (`tests/rf.tcyr`, 205 assertions), the round trips among them:
  RI, MA and DB text read back by `bj_parse_touchstone`, the multi-port layouts bit for bit. In
  `test_write_touchstone_inconsistent_does_not_panic` two cases can't be built (a frequency count unlike the
  matrices', a matrix of another port count: row 36), and its deprecated `write_touchstone` isn't carried (row 1); the
  other four are ported. `test_write_touchstone_deprecated_matches_try` is n/a (row 1). rf.rs's 80 tests are all
  accounted for: 75 ported, 5 n/a. The ledger: 597 ported, 8 n/a, 253 pending of 858. The Cyrius-only groups pin
  the MA and RI text byte for byte, the header alone for no points, the reasons, the builder untouched by a refusal,
  and the freelist's arena over 2,000 writes.
- **Surface parity**: `try_write_touchstone` leaves `IN_PROGRESS` (441 items matched); only rf's six serde impls
  remain, for part 1d. The API snapshot gains 2 public functions (2,062).
- **For M4's range-edge Part**: the writer's MA and DB refuse an entry whose |S| is representable but whose |S|²
  isn't (hisab's √(re² + im²)), as Rust does: 10²⁰⁰ + j cannot be written in DB.
- **Benchmarks** (the `rf/write_touchstone_s2p_201` twin, so 69 of 80 are twinned). An interleaved A/B against part
  1c-1's tree (three rounds, load 0.7 to 1.1) gave 3.89, 3.86 and 3.86 ms, for Rust's 213.82 µs (18×), and left every
  other row unchanged (median ratio 1.00 over 88 rows; the outliers 1 to 3 ns on the smallest rows). The write prints
  1,809 numbers (the 201 frequencies and 804 pairs), and src/fmtf.cyr's shortest digits, Burger and Dybvig's
  algorithm on bignums, cost about 1.8 µs a number (the `fmt_f64_display` row, 1.82 µs): about 3.3 ms of the 3.86.
  The 804 atan2 for the angles are about 0.27. A finding: a faster shortest-digit generator (Ryu's or Grisu's with
  this one as the fallback) would serve every writer, the JSON's and the error texts' among them.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 525,545 ok, 0 failed, 0
  measured, 0 pending, 78,360 accepted (ADR 01: 29,527; ADR 0002: 48,833)** over 603,905 rows; **all 410 digests
  equal**, HEAD's 408 unchanged. `rf.tcyr` passes on all three. The gate is clean (17 steps, no warnings; execution
  coverage 317 of 317 functions), and `cyrius fuzz` passes.

## 80. M4 part 1d: rf's JSON; M4 Part 1 done (2026-11-15)

- **serde's forms, measured first** (a probe of serde_json 1.0.154 on the oracle):
  - `SMatrix` writes `{"ports":2,"s":[[{"re":..,"im":..},..],..],"z0":[..]}`, `TouchstoneData`
    `{"ports":..,"z0":..,"frequencies":[..],"s_matrices":[<SMatrix>,..]}`, `TouchstoneFormat` its variant's name;
  - a `usize` reads from an integer of ASCII digits up to u64::MAX: `-0` is refused (the float −0.0), as are `-1`,
    `2.0`, `2e0`, `18446744073709551616` (a float) and every non-number;
  - a unit variant reads from its name, escapes decoded (`"DbAngle"`), or from a one-member object mapping the
    name to `null` (`{"RealImag":null}`); a second member, a non-null value, an array or a wrong name is refused;
  - serde reads any lengths: a matrix with no rows, a ragged row, a short `z0`, more matrices than frequencies.
- **`src/json.cyr`** gains `BjSMatrix_to_json`/`_from_json`, `BjTouchstoneFormat_to_json`/`_from_json` and
  `BjTouchstoneData_to_json`/`_from_json`. The writers build in one fl_alloc block (64 bytes a Complex) and append
  once; a port count is written as its unsigned word. The readers use the decode region and the whole-document checks
  (rows 11 and 12) as the others do, then:
  - **ADR 0002 row 36, amended**: a shape the fixed structs can't hold is refused as Deserialize ("… (its shape is
    fixed: ADR 0002 row 36)"), checked before anything is allocated: `s` not ports rows of ports entries, `z0` not ports
    long, frequencies and matrices unequal in number, a matrix of another port count than the data's;
  - **a port count past 2^63 − 1** is refused too. The port's word is an i64, and bayan's tree keeps no exact
    integer past i64 (the pre-pass reads it as a float), so `9223372036854775808` and `9223372036854775808.0` can't
    be told apart. A count that large fits no matrix, so only data with no points is touched; row 36 records it;
  - a read matrix or data is a new block, freed on any later refusal, so `out` is untouched and the heap unchanged.
- **The format's writers refuse a value outside the three**: `BjTouchstoneFormat_to_json` and, now first in its
  checks, `bj_try_write_touchstone` (InvalidParameter "unknown Touchstone format 3"; Rust's enum has no other value).
  `_bj_ts_format_name` becomes public for src/json.cyr.
- **`json.tsv` grows by 5,218 rows (66,620)**, appended after beam's (the earlier body an exact prefix, the seeded
  stream unchanged): 756 + 2,867 for `SMatrix`, 443 + 1,098 for `TouchstoneData`, 3 + 51 for `TouchstoneFormat`. The
  layouts are validate.py's JSON_FIELDS ones (`i:` a usize's word or a length, then the elements; the enum `i:` its
  value), and a read may be ragged. The grids: n-ports to 4 and data to 3 points from the shared value list and seeded
  bits, serde's text of each read back in bayan's, Rust's `{}` and `{:e}` notations; every number spelling in each of
  a 2-port's eleven positions (`ports` among them) and in data's port count, R, a frequency and an entry; seeded
  decimals; the shapes (field order, whitespace, unknown keys of every kind, missing, repeated and wrong-typed fields
  at each level, Complex's forms, deep unknown values inside an element, every prefix, trailing data); row 36's
  shapes and port counts at usize's edges; 51 format documents.
  - validate.py's model gains a `usize` (ASCII digits, < 2^64) and JSON_ENUMS; it agreed with every new row on its
    first run. `json.cr.tsv` gains 1,331 rows serde_json's best-effort float parse misreads; `json.adr.tsv` (new) maps
    row 36's 36 (27 matrices, 9 data sets); row 12's harness hook now takes a reading that starts with `i:`, so 10
    more of its rows (233 in all) are counted.
- **Tests**: `tests/json.tcyr`'s `cyrius_rf_json` groups (19 assertions, 73 in all): the text byte for byte, both
  enum forms, the sequence forms, round trips bit for bit (−0 kept), row 36's refusals and reason, `out` untouched, a
  2^62 port count with no points read, and the heap over 1,000 reads and refusals.
- **Surface parity**: rf is ported (13 modules, 447 items matched, none in progress). The API snapshot gains 6 public
  functions (2,068).
- **M4 Part 1 is done**: rf.rs's API in five bites (1a lines to arrays, 1b the matrices, 1c-1 and 1c-2 Touchstone,
  1d the JSON), its 80 tests ported or recorded, 8 of its 8 benches twinned. Part 2 is the range-edge review. Its
  candidates: §76's (the thresholds on rounded values, the complex chains past either end), §79's (MA and DB refusing
  an |S| whose square overflows), and the S-matrices': the 10⁻⁶⁰ tests on rounded |S21|², |T22|² and |det(I − S)|²,
  and the cascade's and S → Z's complex products past either end.
- **Benchmarks**: two Cyrius-only JSON rows (rust-old benches no serde): `json/s_matrix_to_json` 6.18 µs and
  `json/s_matrix_from_json` 7.31 µs for a 2-port. An interleaved A/B against part 1c-2's tree (three rounds, load 0.8
  to 1.2) left every other row unchanged (median ratio 1.00 over 89 rows, all within ±5 %).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 529,386 ok, 0 failed, 0
  measured, 0 pending, 79,737 accepted (ADR 01: 30,858; ADR 0002: 48,879)** over 609,123 rows; **all 416 digests
  equal**, HEAD's 410 unchanged. `json.tcyr` and `rf.tcyr` pass on all three. The gate is clean (17 steps, no
  warnings; execution coverage 323 of 323 functions), and `cyrius fuzz` passes.

## 81. The golden corpus under GitHub's 50 MiB advisory (2026-11-15)

- **The finding**: pushing part 1d, GitHub warned that `tests/golden/json.tsv` was 50.53 MiB (52,982,029 bytes),
  above its 50 MiB recommendation (it refuses a file above 100 MiB). The note column was 25 MB of it, and it repeated
  what each row already holds: a writer row's note gave its inputs, its outputs and its whole document (10.3 MB over
  15,277 rows), a reader row's its document's preview and every value read (15.0 MB over 51,343 rows).
- **The notes are short now** (`scripts/oracle/src/m_json.rs`): a writer's note is a 160-character Debug-quoted
  preview of its document, with its length when longer, and a reader's `unchecked <d> -> Ok` or
  `unchecked <d> -> Err(<serde's message>)`. validate.py reads only that message (serde's "number out of range").
  `json.tsv` is **31.0 MB** (31,006,749 bytes); all 66,620 rows keep their function, class, inputs and outputs byte
  for byte, so `json.cr.tsv` and `json.adr.tsv` are unchanged.
- **A guard**: the generator (`main.rs`) and validate.py fail on any file they write past 48 MiB
  (`MAX_FILE_BYTES`), in both modes, naming the file. The largest now: `json.tsv` 31.0 MB, `json.cr.tsv` 22.5 MB,
  `wave.tsv` 22.0 MB.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): validate.py 609,123 rows, 0 problems; **17 modules,
  0 failing; 529,386 ok, 79,737 accepted**, all 416 digests equal on all three. The gate is clean.

## 82. M4 part 2a-1: the stub match's cancellation-free root (2026-11-15)

M4's Part 2, the range-edge review, runs on the user's answers of 2026-11-15 (plan §0.1): 2a the answers (2a-1 the
stub's root, 2a-2 the cascade, 2a-3 S → Z), 2b round once, 2c the exact thresholds, 2d the range edges.

- **The finding**: `single_stub_match_complex` takes t = (X_L ± disc)/(R_L − Z₀), disc = √(R_L((Z₀ − R_L)² +
  X_L²)/Z₀), and one root cancels wherever disc ≈ |X_L| (X_L − disc for X_L > 0, X_L + disc for X_L < 0). That is
  so in two places: as R_L → Z₀ (outside the 10⁻¹²·Z₀ band the root keeps about 2^-53·Z₀/|R_L − Z₀| of itself:
  2.4·10⁻⁵ relative at R_L = Z₀(1 + 2·10⁻¹²)), and where R_L(Z₀ − R_L) = X_L², where the root → 0 and Rust's came
  out with the wrong sign or orders of magnitude off.
- **ADR 0002 row 39** (`src/rf.cyr`'s `_bj_rf_stub_root`): where disc lies within a factor 2 of |X_L| (the
  difference's cancelling range; Rust's root is within an ulp or two outside it), the cancelling root is
  N/(Z₀(X_L ± disc)) with N = R_L(Z₀ − R_L) − X_L², the same value, since X_L² − disc² = −(R_L − Z₀)(X_L² +
  R_L(R_L − Z₀))/Z₀. Its denominator adds like signs. N would cancel too near the root's zero, so it is R_L·Z₀ − R_L²
  − X_L² from TwoProduct's exact parts (`_bj_k_two_prod`), grown into one expansion (`_bj_g_grow`), compressed and
  rounded (`_bj_g_compress`), on the inputs scaled by the power of two that brings max(Z₀, |X_L|) into [1, 2): exact,
  and the root is scale-free. In the region R_L ≤ 4·Z₀, so the scaled values lie in [2^-480, 8] (every TwoProduct
  exact) unless one falls below 2^-480; there, and for max(Z₀, |X_L|) outside [2^-1022, 2^1023), Rust's chain stays.
  disc itself is Rust's, so its range edges (inner over- or underflowing) are part 2d's.
- **The grid**: 6,179 rows appended to `rf.tsv` (19,899; the earlier body an exact prefix): R_L = Z₀(1 ± ρ) for ρ
  from 10⁻¹ to 1.5·10⁻¹² on nine Z₀ (50, 75, 1, 10⁻³, 10⁶, 3·2^-400, 2^400, 10¹⁵⁰, 10⁻¹⁵⁰) against X_L = ±Z₀·10^k,
  k from −11 to 6; the loads where R_L(Z₀ − R_L) = X_L² ulp by ulp and at relative offsets from 10⁻¹ to 10⁻¹⁵, with
  five exact zeros (N = 0, d = +0); R_L/Z₀ at 2^-479, 2^-480 and 2^-481 (the guard); 600 seeded loads.
- **validate.py** transcribes the port's chain (`rf_stub_root`: TwoProduct, Grow-Expansion and Compress operation for
  operation) and maps **2,911 rows** to row 39 in `rf.adr.tsv` (18 of them part 1a's seeded loads), on MPFR's atan
  and atan2. It also checks the root against the exact one (200 digits) on every row that takes it: **3,775 rows,
  within 2 ulps**; Rust's root was up to 8.7·10¹⁸ ulps off on them. `rf.cr.tsv` gains one row (glibc's atan2
  misrounding on a new row).
- **Tests**: `tests/rf.tcyr`'s `cyrius_stub_root` (5 assertions, 210 in all): both signs of X_L near R_L = Z₀, the root
  just past the band (−0.30000000000113747, Rust's −0.3), an exact zero of N, and the band's limit unchanged.
- **Benchmarks**: the twin `rf/single_stub_match` (a real load, which never takes the new root) 1.90 → 1.91 µs. An
  interleaved A/B against HEAD (three rounds, load 1.15 to 1.45, another repo's tests running) left every row unchanged
  (median ratio 1.00 over 91 rows, all within ±4 % but `material/curie_weiss`'s 7 → 8 ns). A one-off probe of a complex
  load in the region, Z_L = 50.0000005 + j30 on 50 Ω: 43.4 µs, HEAD's 43.1 µs, so the root costs about 0.3 µs.
- **A finding, for the record** (not this bite's): that load costs 43 µs because its d₁ ≈ λ/4 puts the input
  impedance's βl within 10⁻⁸ of π/2, where `src/crmath.cyr`'s sin and cos leave their fast path: about 20 µs each and
  41 µs for `_bj_sin_cos`, against 150 to 200 ns elsewhere. Measured on this tree: π/2 − 10⁻⁶ takes 188 ns, π/2 −
  10⁻⁸, π/2 − 10⁻¹⁰, the double π/2 and π − 10⁻⁸ take 41.2 µs, 10⁻⁸ takes 169 ns. RF hits that band often: every
  matched stub's d = λ/4 (the band's second solution) and every quarter-wave line. A faster correctly rounded path near
  the nonzero multiples of π/2 (the reduction there needs about 2^-60 more of π/2; CORE-MATH's sin and cos do it in
  double-double) would pay off across rf.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 532,653 ok, 0 failed, 0
  measured, 0 pending, 82,649 accepted (ADR 01: 30,859; ADR 0002: 51,790)** over 615,302 rows; **all 416 digests
  equal**, HEAD's unchanged but `rf::single_stub_match_complex`'s (its new rows). `rf.tcyr` passes on all three
  (210). The gate is clean (no warnings; execution coverage 323 of 323 functions); validate.py 615,302 rows, 0
  problems.

## 83. M4 part 2a-2: the cascade rounded once (2026-11-15)

- **The finding**: Rust's `cascade` goes through the transfer matrices, T = T_a·T_b, and back, S = [T12/T22, det T/T22,
  1/T22, −T21/T22], each step rounded. Its T22 = D/(A21·B21), D = 1 − A22·B11, cancels near a resonance: for a mirror
  pair with A22·B11 = 1 − 2^-40 its S came out 4.6·10⁻⁶ off; Rust's own bench thru (0.1, 0.9) with itself gave 2/11
  two ulps high and S12 ≠ S21; wide entries gave ±∞, NaN or the wrong sign for finite values; and its rounded
  |T22|² < 10⁻⁶⁰ refused inputs the exact test accepts.
- **ADR 0002 row 17, amended** (the user's answer of 2026-11-15; `src/crmath_edge.cyr`'s `_bj_cr_cascade`): for finite
  entries the cascade is Redheffer's star product, S11 = A11 + A12·A21·B11/D, S12 = A12·B12/D, S21 = A21·B21/D,
  S22 = B22 + A22·B12·B21/D, each part rounded once from the exact value, and the DivisionByZero is decided exactly
  as 10⁶⁰·|D|² < |A21|²·|B21|²:
  - the fast path (the 16 parts each 0 or in [2^±120]): D and the four complex products as exact expansions of
    TwoProduct's parts, compressed (so D, the resonance's cancellation, carries none); the rest as double-doubles with
    carried absolute bounds (`_bj_e_ed_*`: AccurateDWPlusDW, kimiya's dd products), the threshold decided where the
    difference exceeds its bound, each quotient by `_bj_e_abcd_quot`'s rounding test;
  - the exact path (outside the box, an undecided test): signed scaled bignums in 352-word buffers from fl_calloc (a
    degree-5 value spans at most 10,495 bits), the quotients by restoring division; a refused allocation is Allocation
    "the cascade's exact buffers".

  A part exactly 0 keeps Rust's zero where Rust's chain gives one (+0 otherwise), as `propagate_beam`'s; a
  non-finite entry keeps Rust's test and chain (`_bj_rf_cascade_chain`). The S21 tests stay Rust's until part 2c.
- **Verified**: `cr_verify.py`'s cascade and cascade_slow families (passive networks, resonances to 2^-70, the
  threshold ulp by ulp, zeros and thrus, a cancelling S11, the box edges, 10^±300 entries, subnormals; the driver
  now reads 24 words a call): **18,000 calls equal to the exact value (Python's Fraction) on x86_64 and aarch64**,
  the fast path and the exact path alone, with the exact-zero bits.
- **The grid**: 732 rows appended to `rf.tsv` (20,631; the earlier body an exact prefix): resonances A22·B11 =
  1 − 2^-k·(1, (1 + j)/2, −j) for k to 70 on four A22; D = 2^-k exactly against B21 within three ulps of
  2^-k·10³⁰ (real and imaginary); thrus and an exact 0 in every position; a cancelling S11; one entry at 2^±119,
  2^±120 and 2^±121 in each of the 16 places; 60 wide pairs, 16 subnormal entries, 200 seeded networks.
- **validate.py** models the port's cascade on Fractions (`rf_cascade_exact`, independent of the Cyrius code) and maps
  **1,068 rows** to row 17 in `rf.adr.tsv` (498 of them part 1b's): 1,025 values (of their 7,949 compared parts,
  1,938 equal, 2,299 an ulp off, 2,708 within 1,000 ulps, 1,004 past that, 113 of the wrong sign, 251 ±∞ or NaN for
  a finite value), 24 now refused and 19 refusals now computed.
- **Tests**: `tests/rf.tcyr`'s `cyrius_cascade_exact` (7 assertions, 217 in all): the bench thru (2/11 and 9/11, S12
  = S21, Rust's −0 kept), the resonance's S11 and S22, the threshold either side of 2^-40·10³⁰ with out untouched, and
  the heap over 20 exact paths.
- **Benchmarks**: `rf/s_cascade` (Rust's bench thru with itself) **523 ns → 3.47 µs** in a focused interleaved A/B
  against HEAD (three rounds, 3.47 µs each, load 0.75), Rust's 100 ns: the cost of rounding once, which the user's
  answer accepted ("a few µs"). A full interleaved A/B (three rounds, load 2.0 to 3.0, another repo's benches pinned
  beside it) left the other 90 rows unchanged (median ratio 0.99). One-off probes of `_bj_cr_cascade` alone: the fast
  path 2.84 µs (real entries) and 3.71 µs (complex), the exact path 53 µs and 73 µs; an "ed" product 47 ns, a 5-term
  expansion compressed 131 ns. The complex products and Im D first took expansions too (3.0 and 4.4 µs); an
  AccurateDWPlusDW of two exact TwoProducts is as cancellation-proof (its bound is relative to the exact sum) at half
  the cost, so only Re D, with its 1, keeps the expansion. The cascade also runs Rust's chain (about 0.5 µs) when a part
  is exactly 0, for that zero's sign.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 532,317 ok, 0 failed, 0
  measured, 0 pending, 83,717 accepted (ADR 01: 30,859; ADR 0002: 52,858)** over 616,034 rows; **all 416 digests
  equal**, HEAD's unchanged but `SMatrix::cascade`'s. `rf.tcyr` passes on all three (217). The gate is clean (no
  warnings; execution coverage 323 of 323 functions); validate.py 616,034 rows, 0 problems.

## 84. M4 part 2a-3: S → Z rounded once; 2a done (2026-11-15)

- **The finding**: Rust's `s_to_z` rounds det = (1 − S11)(1 − S22) − S12·S21 first, and it cancels toward an ideal
  thru. Rust's own bench network (S11 = S22 = 0.1, S12 = S21 = 0.9) was refused as singular, since (1 − 0.1)² and
  0.9² round alike; the doubles' exact det is −5.0·10⁻¹⁷, and Z = −1.8·10¹⁸ in all four entries. A near-thru at
  det ≈ 2^-31 came out 2.2·10⁻¹⁰ off, and wide entries gave ±∞ or NaN for finite values.
- **ADR 0002 row 17, amended** (the user's answer of 2026-11-15; `src/crmath_edge.cyr`'s `_bj_cr_s2z`): for finite
  entries each part of Z11 = Z₀((1 + S11)(1 − S22) + S12·S21)/det, Z12 = 2Z₀·S12/det, Z21 = 2Z₀·S21/det and
  Z22 = Z₀((1 − S11)(1 + S22) + S12·S21)/det is rounded once from the exact value, and |det|² < 10⁻⁶⁰ is decided
  exactly (10⁶⁰·|det|² < 1):
  - the fast path (the cascade's box): det, N11 = 1 + a − d − a·d + b·c and N22 = 1 − a + d − a·d + b·c as exact
    expansions of the 1, the parts and the eight TwoProducts of a·d and b·c (up to 11 doubles, `_bj_e_ed_exact` now
    takes 16), so none carries a cancellation; 2b and 2c exact; then the cascade's carried-bound steps
    (`_bj_e_s2z_out`: N·conj(det), ×Z₀, the rounding test);
  - the exact path: signed scaled bignums from fl_calloc; Allocation "the S→Z conversion's exact buffers".

  A part exactly 0 keeps Rust's zero where Rust's chain gives one (`_bj_rf_s2z_chain`), +0 otherwise; a non-finite
  entry keeps Rust's test and chain.
- **Verified**: `cr_verify.py`'s s2z and s2z_slow families (passive networks, near-thrus to 2^-70, |det| straddling
  10⁻³⁰ ulp by ulp, lossless networks, exact zeros, the box edges, 10^±300 entries, subnormals): **18,000 calls equal
  to the exact value (Fraction) on x86_64 and aarch64**, with the cascade's 18,000 rerun beside them.
- **The grid**: 492 rows appended to `rf.tsv` (21,123): near-thrus at 2^-k, k to 70; |det| within three ulps of
  10⁻³⁰ either side (S11 = S22 = 1 + k·2^-52, S12 = 1, S21 = k²·2^-104 ∓ 10⁻³⁰); 40 lossless networks (S from a
  reactance matrix in floats, so Re Z is tiny); a matched 2-port and an exact 0 in each place; one part or Z₀ at
  2^±119, 2^±120, 2^±121; 60 wide and 8 subnormal; 200 seeded networks.
- **validate.py** models the port on Fractions (`rf_s2z_exact`) and maps **875 rows** to row 17 (441 on earlier
  grids): 874 values (of their 6,514 compared parts, 2,217 equal, 2,366 an ulp off, 1,473 within 1,000 ulps, 458 past
  that, 35 of the wrong sign, 478 ±∞ or NaN for a finite value) and the bench network's refusal, now computed.
- **Tests**: `tests/rf.tcyr`'s `cyrius_s_to_z_exact` (6 assertions, 223 in all): the bench network, the near-thru's
  Z11 and Z22, an ideal thru refused with out untouched, a matched 2-port's Z₀·I with Rust's zeros, the heap over 20
  exact paths.
- **Benchmarks** (rf.rs has no S → Z bench): a focused interleaved A/B against HEAD (three rounds, load 0.75): a
  2-port attenuator 185 ns → 3.26 µs (the fast path and, for its exact zeros' signs, Rust's chain), complex entries
  181 ns → 4.97 µs.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 531,934 ok, 0 failed, 0
  measured, 0 pending, 84,592 accepted (ADR 01: 30,859; ADR 0002: 53,733)** over 616,526 rows; **all 416 digests
  equal**, HEAD's unchanged but `rf::s_to_z`'s. `rf.tcyr` passes on all three (223). The gate is clean (no warnings;
  execution coverage 323 of 323 functions); validate.py 616,526 rows, 0 problems.
- **2a is done**: the user's three answers. Next, 2b (round once: the lossless line's √(L′/C′), the phase velocity
  1/√(L′C′), the quarter-wave transformer's √(Z₀Z_L) and `normalize_impedance`'s Z/Z₀).

## 85. M4 part 2b: round once in rf (2026-11-15)

- **ADR 0002 row 13, amended** (the user's "round once" decision of 2026-11-12): rf's formulas that are one product,
  quotient or root of inputs are correctly rounded:
  - `characteristic_impedance_lossless` √(L′/C′) by kimiya's √(x/y) (`_bj_k_cr_sqrt_div`, as maxwell's `impedance`);
  - `phase_velocity` 1/√(L′C′) and `quarter_wave_transformer` √(Z₀·Z_L) by bijli's √ of a monomial
    (`_bj_cr_sqrt_mono`, as maxwell's `wave_speed` and `refractive_index`);
  - `normalize_impedance` Z/Z₀ by one IEEE division a part, where Rust took Z·(1/Z₀).

  The inputs are refused unless finite and positive (Z₀ finite and not below 10⁻³⁰ for Z/Z₀), so no ∞ reaches the
  roots; a non-finite part of Z divides as Rust's product did. `denormalize_impedance` (z·Z₀), the Smith chart, the
  losses and the array factors have a sum or a transcendental inside, or were one IEEE operation already.
- **The finding**: over 20,000 typical inputs each (L′ from 10⁻⁹ to 10⁻⁵ H/m, C′ from 10⁻¹³ to 10⁻⁹ F/m, impedances
  from 1 to 1,000 Ω), Rust's value was an ulp off on 12.8 % (√(L′/C′)), 29.2 % (1/√(L′C′)), 12.5 % (√(Z₀·Z_L)) and
  25.5 % (Z·(1/Z₀)). And its quotient or product left the range though the root is representable: √(10³⁰⁰/10⁻³⁰⁰)
  was ∞ for 10³⁰⁰, 1/√(10²⁰⁰·10²⁰⁰) was 0 for 10⁻²⁰⁰, √(10²⁰⁰·10²⁰⁰) was ∞ for 10²⁰⁰.
- **The grid**: 3,135 rows appended to `rf.tsv` (24,258): 300 typical seeded values of each kind (the line pair also
  feeds both `characteristic_impedance_lossless` and `phase_velocity`, the impedance pair `normalize_impedance` and
  `denormalize_impedance`), 200 over 10^±307, the ranges' edges (10^±154 to DBL_MAX and the subnormals crossed, Z₀ near
  DBL_MAX so that 1/Z₀ is subnormal), and exact roots (k²·10⁻⁷/10⁻⁷, k²/4, k²·9).
- **validate.py** rounds each exactly (`rn_sqrt` of the Fraction, Python's IEEE division) and maps **1,105 rows** to
  row 13 in `rf.adr.tsv` (440 on earlier grids): `characteristic_impedance_lossless` 247, `phase_velocity` 413,
  `quarter_wave_transformer` 198, `normalize_impedance` 247; 399 of them were ∞ or 0 for a finite value.
- **Tests**: `tests/rf.tcyr`'s `cyrius_round_once` (5 assertions, 228 in all).
- **Benchmarks** (rf.rs benches none of the four): a focused interleaved A/B against HEAD (three rounds, load 0.75):
  `characteristic_impedance_lossless` 25 → 88 ns, `phase_velocity` 28 → 140 ns, `quarter_wave_transformer`
  23 → 142 ns (maxwell's `wave_speed`, the same root, is 141 ns), `normalize_impedance` 17 → 16 ns.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 533,964 ok, 0 failed, 0
  measured, 0 pending, 85,697 accepted (ADR 01: 30,859; ADR 0002: 54,838)** over 619,661 rows; **all 416 digests
  equal**, HEAD's unchanged but the four functions' and `denormalize_impedance`'s (new rows). `rf.tcyr` passes on all
  three (228). The gate is clean (no warnings; execution coverage 323 of 323 functions); validate.py 619,661 rows, 0
  problems. Next, 2c: the exact thresholds.

## 86. M4 part 2c: rf's exact thresholds (2026-11-15)

- **ADR 0002 row 17, amended**: rf's remaining tests against a literal threshold are decided exactly on the inputs,
  the threshold as its decimal value:
  - `reflection_coefficient`'s |Z_L + Z₀| < 10⁻³⁰: `src/crmath_cmp.cyr`'s band (the rounded sum decides outside 64
    ulps of the double 1e-30), then `_bj_x_sum_below` (TwoSum's exact parts), as wave's |n₁ + n₂|;
  - `reflection_coefficient_complex`'s |Z_L + Z₀|² < 10⁻⁶⁰ and `s_to_t`'s |S21|² < 10⁻⁶⁰ for a finite S21 (the
    cascade's first step): the band on the rounded sum of squares, then `_bj_x_r2_below` on the parts, as
    relativity's |R|. A NaN or ∞ part keeps Rust's test (not refused);
  - the stub match's |R_L − Z₀| ≤ 10⁻¹²·Z₀ and `z0_equal`'s |a − b| ≤ 10⁻¹²·max(|a|, |b|) (the cascade, S → Z,
    `uniform_z0` and Touchstone's writer): `src/rf.cyr`'s new `_bj_rf_within12` scales m into [2^-900, 2^900] by
    2^±600, lets the rounded comparison decide outside ±0.1 %, and inside it takes the exact |a − b| (Sterbenz),
    |a − b|·10¹² = p + e by TwoProduct (10¹² is a double) and the sign of m − p − e by Grow-Expansion.

  The input impedance's relative 10⁻¹⁵ stays on its rounded den (the user's answer of 2026-11-15).
- **What Rust's tests got wrong**: the double 1e-30 lies 0.37 ulp above 10⁻³⁰, so a sum Z_L + Z₀ whose exact value is
  below 10⁻³⁰ but rounds to 1e-30 was accepted (Γ = −1.2·10⁻¹⁵ for a DivisionByZero); the squared tests likewise.
  For normal values the 10⁻¹² tests were already exact: |a − b| is a multiple of the smaller ulp, and RN(1e-12·m)
  meets such a multiple only where it equals it exactly (with m = k·10¹²·2^-47, 10⁻¹²·m is exactly k ulps, and both
  agree on the grid's rows). Below about 2^-970, 1e-12·m rounds to the subnormal grid: at m = (j·10¹² − r)·2^-1042
  Rust took Z₀ values 10⁻¹² apart as equal.
- **The grid**: 1,974 rows appended to `rf.tsv` (26,232): `knife_edge2`'s sums at 10⁻³⁰ (with 50 Ω among the
  bases; the pair feeds `reflection_coefficient` and `quarter_wave_transformer`); `knife_edge`'s planar pairs at
  |d| ≈ 10⁻³⁰ as Z_L + Z₀ and, from the origin, as S21 of either cascaded 2-port; the 10⁻¹² cut at m = k·10¹²·2^-47
  (normal: no disagreement) and at m = (j·10¹² − r)·2^-1042 (subnormal 1e-12·m), each through the stub matches, a
  cascade of two Z₀'s and S → Z of a mixed pair.
- **validate.py**: row 17's exact tests (and 2a's rounding once) now apply together while `RF17["on"]`; a mapping
  names the first of rows 13, 37 and 39 that accounts for its whole change, else 17. With the switch the old grids
  map exactly as before (the files unchanged before the new rows). `rf.adr.tsv` maps 6,598 rows: of this bite's,
  148 reflections now refused (42 real, 106 complex), 72 cascades refused for |S21|² < 10⁻⁶⁰, 12 cascades and 12
  S → Z conversions refused for unequal Z₀; the grid's other cascade and S → Z rows are 2a's (282), its quarter-wave
  rows 2b's (69), its stub rows 2a-1's (28).
- **The stub's band shows in no row**: where its exact decision differs (Z₀ ≈ 2^-990), the input impedance's
  |den|² underflows and refuses first. That range edge is part 2d's; the band's rows are in the grid for it.
- **Tests**: `tests/rf.tcyr`'s `cyrius_exact_thresholds` (4 assertions, 232 in all).
- **Benchmarks**: a full interleaved A/B against HEAD (three rounds, load 0.75 to 1.1) left every row unchanged
  (median ratio 1.00 over 91 rows) but `rf/reflection_coefficient`, 15 → 23 ns: the band helper ran on every call.
  A first compare now skips it for a rounded |Z_L + Z₀| at or above 2·10⁻³⁰ (|·|² at or above 2·10⁻⁶⁰), which cannot
  lie below the threshold, and the squared tests' exact path is one helper (`_bj_rf_sq_below60`), so its locals stay
  off the common path. A focused interleaved A/B (three rounds, load 0.6): `rf/reflection_coefficient` 14–16 → 14–15
  ns, `reflection_coefficient_complex` 50–51 → 50–51 ns, `rf/s_cascade` 3.48 → 3.52 µs (+1.2 %: five exact
  `z0_equal` calls).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 535,315 ok, 0 failed, 0
  measured, 0 pending, 86,320 accepted (ADR 01: 30,859; ADR 0002: 55,461)** over 621,635 rows; **all 416 digests
  equal**, HEAD's unchanged but those of the functions the grid grew (the reflections, the quarter-wave transformer,
  the stub matches, the cascade, S → Z). `rf.tcyr` passes on all three (232). The gate is clean (no warnings;
  execution coverage 323 of 323 functions); validate.py 621,635 rows, 0 problems. Next, 2d: the range edges.

## 87. M4 part 2d-1: the input impedance at its range edges (2026-11-15)

M4 part 2d, the range edges (ADR 0002 row 18), runs in five bites: 2d-1 the input impedance; 2d-2 the lossy Z₀ and
the propagation constant (their complex chains into hisab's square root, which `BjComplex_sqrt` already carries at its
edges); 2d-3 the reflection coefficients, the Smith chart's Γ and z, and the stub match's discriminant and roots; 2d-4
Touchstone's MA and DB magnitudes; 2d-5 the array factors' ψ past DBL_MAX.

- **The finding** (the user's answer of 2026-11-15 named the first): `input_impedance_lossless_complex` refused a load
  past |Z_L| = 1.3·10¹⁵⁴, since Z₀ + |Z_L| overflowed and every |den| fell under 10⁻¹⁵·∞; refused finite β and l
  whose product passed DBL_MAX ("βl must be finite"); and refused impedances near 10⁻²⁰⁰ Ω as unbounded, since |den|²
  underflowed to 0. Where |den|² overflowed it passed the 10⁻¹⁵ test by accident: a λ/4 line of 10³⁰⁰ Ω on a 1 Ω
  load (Z_in = 10⁶⁰⁰) was accepted.
- **ADR 0002 row 18, amended** (`src/rf.cyr`): Rust's chain stays where every intermediate is normal, which the guard
  checks cheaply: each nonzero input in [2^-100, 2^101) (βl, its sine and cosine, the products and |Z_L|² then lie
  within 2^±400), the two cancelling sums X·c + Z₀·s and Z₀·c − X·s exactly 0 or at least 2^-400 (|den|², Smith's
  ratio and its products then stay above 2^-910), and the quotient's parts 0 or normal. Elsewhere `_bj_rf_zin_ux`
  runs the chain with an unbounded exponent: kimiya's 53-bit `_bj_k_ux_*`, βl's sine and cosine by src/crmath_vx.cyr
  (Payne–Hanek past DBL_MAX), the 10⁻¹⁵ test on the 53-bit |den| and 10⁻¹⁵·scale, Smith's division, and the products
  by Z₀ rounded once. A non-finite β or l is still refused (Rust's reason); a finite pair no longer is. The stub match
  calls it for each solution, so its tiny-Z₀ rows now reach row 17's exact band (M4 part 2c): 48 of them.
- **The grid**: 739 rows appended to `rf.tsv` (26,971): thirteen magnitudes from 10⁻³⁰⁰ to 10³⁰⁰ (the box's 2^±100
  and 2^±101 among them) for the load, Z₀ or both, at βl = 0.7 and π/2; βl at 10⁴⁰⁰, DBL_MAX·2, 10⁻⁴⁰⁰, subnormal
  and the like; 400 seeded wide values; the stub matches at each magnitude.
- **validate.py** transcribes the unbounded chain (`rf_zin_ux`: rnd53 per operation, MPFR's widest range for the sine
  and cosine) under `RF18`, so the whole chain is the model and the guard is proved on every row Rust's chain keeps.
  `rf.adr.tsv` maps **303 input impedances** (222 refusals as unbounded and 37 for βl now computed, 44 values) and
  **241 stub matches** (215 refusals now computed, 26 now refused: Rust's |den|² had overflowed at d = λ/4 against
  Z₀ = 10³⁰⁰ or DBL_MAX), all passing; with 2c's 48, 592 rows.
- **Tests**: `tests/rf.tcyr`'s `cyrius_zin_edges` (5 assertions, 237 in all).
- **Benchmarks**: a focused interleaved A/B against HEAD (four rounds; the load rose from 1.1 to 2.6, other repos'
  benches beside it): `rf/input_impedance_complex` 236–264 → 265–297 ns, about +30 ns for the guard, after inlining
  its common case (six helper calls first cost +54 ns); `rf/single_stub_match` (two input impedances) about +0.1 µs.
- **ADR 0002's table**: row 18 carried four unescaped |x|, |y| pipes from M3 part 2d-1c, which split it into extra
  cells; escaped (every row now has its four cells).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 535,456 ok, 0 failed, 0
  measured, 0 pending, 86,918 accepted (ADR 01: 30,860; ADR 0002: 56,058)** over 622,374 rows; **all 416 digests
  equal**, HEAD's unchanged but the input impedance's and the stub matches'. `rf.tcyr` passes on all three (237). The
  gate is clean (no warnings; execution coverage 323 of 323 functions); validate.py 622,374 rows, 0 problems. Next,
  2d-2.

## 88. M4 part 2d-2: the lossy Z₀ and γ at their range edges (2026-11-15)

- **The finding**: `characteristic_impedance_lossy` and `propagation_constant` end in hisab's square root,
  |z| = √(re² + im²) and from_polar(√|z|, arg z / 2), after Smith's division of (R + jωL) by (G + jωC) or their
  product. Rust gave NaN where ωL passed DBL_MAX (the true Z₀ is 10¹⁰⁵ at R = 0.1, L = 10²⁰⁰, C = 10⁻¹⁰, ω = 10²⁰⁰),
  NaN where ωC underflowed to 0, which sent a lossy line into the DC branch and its "ω = 0" test (Z₀ = 5·10⁹⁶), ∞ +
  jNaN where the DC branch's L′/C′ overflowed (Z₀ = 10³⁰⁰), and NaN where γ's product overflowed (γ = 10²⁰⁰(1 + j)).
- **ADR 0002 row 18, amended** (`src/rf.cyr`): Rust's chain stays where each nonzero input lies in [2^-60, 2^61)
  (ωL and ωC within 2^±122, Smith's ratio and its products within 2^±370) and the quotient's or product's parts are 0
  or within [2^-500, 2^501) (so |z|², √|z|, arg z and its half are normal); a non-finite input keeps it too. Elsewhere
  the chain runs on 53-bit values with an unbounded exponent (`_bj_rf_lossy_ux`, `_bj_rf_gamma_ux`):
  - Smith's division, shared with 2d-1's input impedance (`_bj_rf_cdiv_ux`);
  - arg z to 53 bits (`_bj_rf_atan2_ux`): a zero part leaves only the signs; parts whose exponents lie within 1,000
    of each other scale by one power of two into doubles, whose correctly rounded atan2 is normal; further apart,
    y/x for x > 0 (atan t lies within t³/3 of t, and a quotient of two 53-bit values is never a 54-bit midpoint),
    ±π for x < 0, ±π/2 when |x| ≪ |y|;
  - hisab's square root (`_bj_rf_sqrt_polar_ux`): √|z| to 53 bits, the exact half-angle, its sine and cosine to 53
    bits (src/crmath_vx.cyr), the polar products rounded once.

  The DC branch keeps Rust's outcomes (√(L′/C′) through the same square root, +∞ − j∞, NaN) on the exact zero of
  ωC.
- **The grid**: 724 rows appended to `rf.tsv` (27,695): 22 hand cases (ωL and ωC past either end, the DC branch's
  L′/C′ both ways, |z|² past either end, arguments with |im/re| below 2^-1000 for either sign of re and the reverse,
  subnormal and DBL_MAX inputs, negative ω), one input at 2^±59 to 2^±62 in each place, 300 seeded wide 5-tuples; each
  through both functions.
- **validate.py** transcribes the chains (`rf_lossy_ux`, `rf_gamma_ux`, `rf_sqrt_polar_ux`, `ux_atan2` on MPFR's
  widest range) under `RF18`. `rf.adr.tsv` maps **294 lossy Z₀** (150 NaN or ∞ in Rust) and **321 γ** (260), all
  passing; 2d-1's 303 input impedances are unchanged by the shared Smith division.
- **Tests**: `tests/rf.tcyr`'s `cyrius_lossy_edges` (5 assertions, 242 in all).
- **Benchmarks** (rf.rs benches neither): a focused interleaved A/B against HEAD (three rounds, load 6.3 to 7.4: other
  repos' jobs): a typical lossy Z₀ 499–515 → 529–534 ns, about +33 ns for its two guards; γ 544–650 → 564–578 ns,
  within the noise.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 535,565 ok, 0 failed, 0
  measured, 0 pending, 87,533 accepted (ADR 01: 30,860; ADR 0002: 56,673)** over 623,098 rows; **all 416 digests
  equal**, HEAD's unchanged but the two functions'. `rf.tcyr` passes on all three (242). The gate is clean (no
  warnings; execution coverage 323 of 323 functions); validate.py 623,098 rows, 0 problems. Next, 2d-3.

## 89. M4 part 2d-3a: the reflection coefficients and the Smith chart at their range edges (2026-11-15)

2d-3 runs in two bites: 2d-3a the Smith-division family here, 2d-3b the stub match (its discriminant, roots, crecip
and stub length, with an unbounded exponent end to end).

- **The finding**: `reflection_coefficient` gave 0 where Z_L + Z₀ passed DBL_MAX (Γ(1.5·10³⁰⁸, 10³⁰⁸) = 0.2) and ∞
  where Z_L − Z₀ did (Γ(10³⁰⁸, −0.9·10³⁰⁸) = 19); `reflection_coefficient_complex` gave NaN for such sums (Γ = j for
  Z_L = 10³⁰⁸(1 + j), Z₀ = 10³⁰⁸(1 − j)); `gamma_from_z` gave NaN for z = DBL_MAX(1 + j), where Γ = 1, as Smith's
  d = b.re + b.im·r overflowed. Smith's division is otherwise robust: most of the new edge rows agree with Rust.
- **ADR 0002 row 18, amended** (`src/rf.cyr`): `_bj_rf_cdiv_guard` keeps rf.rs's Smith division where its four
  inputs are each 0 or within [2^-250, 2^251) and the quotient's parts 0 or normal (every intermediate is then normal:
  the ratio at least 2^-501, its products 2^-751, a cancelling numerator 2^-804, a nonzero quotient 2^-1056); else
  `_bj_rf_cdiv_edge` runs it on 53-bit values (2d-2's `_bj_rf_cdiv_ux`, now with its two last divisions rounded once).
  The complex reflection and the Smith chart call it; the Smith chart's z ± 1 and 1 ± Γ never overflow. The reflection
  coefficients take Z_L ± Z₀ to 53 bits where a binary64 sum is not finite (`_bj_rf_refl_ux`, `_bj_rf_refl_cx_ux`).
- **The grid**: 1,137 rows appended to `rf.tsv` (28,832): ten real pairs near DBL_MAX of both signs (each through the
  real reflection, three complex ones and the quarter-wave transformer), nine pairs of far-apart magnitudes (the box's
  2^±250 and 2^±251 among them) through the complex reflection and the Smith chart, the Smith chart near DBL_MAX and at
  its poles with subnormal imaginary parts, and 200 seeded wide values.
- **validate.py** models the four with 53-bit sums and Smith's division (`rf_cdiv_ux(..., fin=True)`) under `RF18`;
  `rf.adr.tsv` maps **66 rows** (8 real, 29 complex, 12 Γ, 17 z), all passing.
- **Tests**: `tests/rf.tcyr`'s `cyrius_smith_edges` (4 assertions, 246 in all).
- **Benchmarks** (a focused interleaved A/B against HEAD, three rounds at load 12.5 to 12.8, so rough):
  `rf/reflection_coefficient` 19–24 → 22–27 ns; `reflection_coefficient_complex` 68–80 → 99–105 ns and
  `gamma_from_z` 43–48 → 67–93 ns, about +20–25 ns: the guard is a function between the caller and Smith's division,
  with its own locals. Splitting its edge path (four scratch arrays) into a separate helper and inlining its quotient
  test took it from about +45 ns; inlining the guard into its three callers would likely save most of the rest (a
  finding, left for now).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 536,609 ok, 0 failed, 0
  measured, 0 pending, 87,626 accepted (ADR 01: 30,860; ADR 0002: 56,766)** over 624,235 rows; **all 416 digests
  equal**, HEAD's unchanged but the reflections', the Smith chart's and the quarter-wave transformer's (new rows).
  `rf.tcyr` passes on all three (246). The gate is clean (no warnings; execution coverage 323 of 323 functions);
  validate.py 624,235 rows, 0 problems. Next, 2d-3b: the stub match.

## 90. M4 part 2d-3b: the stub match end to end at its range edges (2026-11-15)

- **The finding**: `single_stub_match_complex` (and so `single_stub_match`) refused a load of 10³⁰⁰(1 + j) on
  10³⁰⁰ Ω, and the same load scaled to 2^-1074 Ω, as 1/Z₀ overflowed, though both have the unit case's solutions.
  Where its discriminant overflowed (X_L = −3·10¹⁵², Z₀ = 10¹⁵⁰), both roots came out ±∞ and both d = λ/4, so the two
  solutions collapsed into one; 2d-1's 26 refusals at such loads were the same overflowed discriminant.
- **ADR 0002 row 18, amended** (`src/rf.cyr`): Rust's chain stays where R_L, X_L and Z₀ lie in 2^±100 and each Z_in
  is normal with |Z_in|² within 2^±1000 (`_bj_rf_zin_norm`); elsewhere `_bj_rf_stub_ux` runs the whole chain on
  53-bit values with an unbounded exponent:
  - the discriminant and both roots, the cancelling one by row 39's cancellation-free form on the same scaled values
    (`_bj_rf_stub_root_ux`; row 39's root split into `_bj_rf_stub_root_core`, which both chains share);
  - arctan t and d (`_bj_rf_stub_d_ux`), d kept to 53 bits for the input impedance and rounded once as an output;
  - the input impedance's 53-bit quotient times Z₀ (`_bj_rf_zin_q53`, split out of 2d-1's `_bj_rf_zin_ux`);
  - crecip's 1/|Z_in|² branch, y₀ = 1/Z₀ and atan2(y₀, Im y), and ℓ rounded once (`_bj_rf_stub_len_ux`).
- **The grid**: 550 rows appended to `rf.tsv` (29,382): Z₀ subnormal (1/Z₀ past DBL_MAX), at 2^±100 and 2^±101 (the
  box's edges), 10^±300 and near DBL_MAX (2Z₀ past it), each with six loads scaled by Z₀ and the same six unscaled;
  loads whose |Z_in|² leaves the range; roots below 2^-1022 (d subnormal); 200 seeded wide loads through each stub.
- **validate.py** transcribes the chain (`rf_stub_ux`, `rf_stub_d_ux`, `rf_stub_len_ux`, with `rf_atan53`) under
  `RF18`. `rf.adr.tsv` maps **2,415 stub rows** to row 18 (256 real, 2,159 complex; 241 at HEAD), all passing, most
  on 2a-1's grid at Z₀ of 2^±400 and 10^±150; **405** of them (174 real, 231 complex) were refusals in Rust and are
  computed now. A labeling fix: a row that differs only with both the port's row (13, 37 or 39) and row 18 on is
  labeled 18, where it had fallen through to 17; the stubs keep 72 row-17 rows (24 real, 48 complex). One corrected
  oracle row (a complex stub whose square root rf.cr.tsv corrected) is now a row-18 mapping: `rf.cr.tsv` 87 rows,
  `rf.adr.tsv` 10,101.
- **Tests**: `tests/rf.tcyr`'s `cyrius_stub_edges` (3 assertions, 249 in all): the 10³⁰⁰ and 2^-1074 loads give the
  unit case's solutions bit for bit, and the overflowed discriminant's load two distinct ones.
- **Benchmarks** (rf.rs benches neither stub at these values): a focused interleaved A/B against HEAD (three rounds,
  load 8.1 to 8.6, so rough): `rf/single_stub_match` 2.50–2.76 → 2.56–2.83 µs and the complex one 2.43–2.62 →
  2.48–2.59 µs, both within the noise; the guard is a box test on the inputs and a range test on each Z_in.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 534,962 ok, 0 failed, 0
  measured, 0 pending, 89,823 accepted (ADR 01: 30,859; ADR 0002: 58,964)** over 624,785 rows; **all 416 digests
  equal**, HEAD's unchanged but the two stubs'. `rf.tcyr` passes on all three (249). The gate is clean (no warnings;
  execution coverage 323 of 323 functions); validate.py 624,785 rows, 0 problems; `cr_verify.py`'s cascade and S→Z
  families 36,000 calls on x86_64 and aarch64, 0 mismatches. Next, 2d-4: Touchstone's MA and DB magnitudes.

## 91. M4 part 2d-4a: the Touchstone reader at its range edges (2026-11-15)

2d-4 runs in two bites: 2d-4a the reader here, 2d-4b the writer (|S| through row 14's |z|, DB past DBL_MAX, angles
below 2^-1022).

- **The finding**: `parse_touchstone` reads an MA or DB pair as (r·cos θ, r·sin θ), r = v₁ or 10^(v₁/20), θ = v₂·π/180.
  It refused a DB pair whose magnitude passed DBL_MAX though both parts were representable (6166 dB at 45°: |S| =
  10^308.3, both parts 1.41·10³⁰⁸; a pair at 45° stays representable up to 20·log₁₀(√2·DBL_MAX) ≈ 6168.1 dB). It
  read 10³⁰⁰ at 2^-1074° as 10³⁰⁰ + j0, as θ underflowed (8.6·10⁻²⁶ j), and lost digits where θ or 10^(v₁/20) was
  subnormal (10¹⁰ at 10⁻³¹⁰°: im 1.7453292519952424·10⁻³⁰² for 1.745329251994324·10⁻³⁰²; −6200 dB at 30°).
- **ADR 0002 row 18, amended** (`src/touchstone.cyr`): Rust's chain stays where DB's 10^(v₁/20) is a normal double
  above 2^-1022 and θ is above 2^-1022 in magnitude or v₂ = 0 (a power rounded to 2^-1022 may have come from below
  it); elsewhere `_bj_ts_polar_ux` runs the pair with an unbounded exponent: r = v₁, or 10^(v₁/20) to 53 bits
  (src/crmath_ux.cyr's `_bj_ux_pow` of 1.25·2³, its stand-in past 10^4343 overflowing every product as the true value
  does); θ to 53 bits where it is below 2^-1022, its sine to 53 bits θ itself (within θ³/6) and its cosine 1,
  elsewhere the correctly rounded sine and cosine of Rust's θ; each product rounded once. A part past DBL_MAX is
  still refused ("is not representable").
- **The grid**: 904 rows appended to `rf.tsv` (30,286), all `parse_touchstone`: DB magnitudes from 6165 to 10⁶ dB and
  −6150 to −10⁶ dB, with 20·log₁₀ DBL_MAX and 20·log₁₀ 2^-1022 and three doubles each side, at eleven angles; angles
  from 10⁻³⁰⁰ down to 2^-1074 degrees of either sign, with v at θ's 2^-1022 and two doubles each side, under nine MA
  and five DB magnitudes; one 2-port and one 3-port file mixing them; 240 seeded 1-port files.
- **validate.py** models the pair's chain (`rf_ts_polar_ux`: MPFR's widest-range power, sine and cosine) under
  `RF18`, and `rf_adr_rows` now evaluates `parse_touchstone` too. `rf.adr.tsv` maps **194 parse rows** to row 18
  (10,295 in all), all passing, every one on the new grid; **45** were refusals in Rust and read now. `rf.cr.tsv`
  is unchanged (87).
- **Tests**: `tests/rf.tcyr`'s `cyrius_touchstone_read_edges` (5 assertions, 254 in all): the four pairs above to the
  bit, and 6168.2 dB at 45° still refused.
- **Benchmarks** (a focused interleaved A/B against HEAD, three rounds at load 3.6): `rf/parse_touchstone_s2p_201`
  (MA) 568–584 → 572–593 µs, and the same file in DB 1.17–1.19 → 1.18–1.21 ms, both within the noise; the guard is
  two compares per pair.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 535,672 ok, 0 failed, 0
  measured, 0 pending, 90,017 accepted (ADR 01: 30,859; ADR 0002: 59,158)** over 625,689 rows; **all 416 digests
  equal**, HEAD's unchanged but `parse_touchstone`'s. `rf.tcyr` passes on all three (254). The gate is clean (no
  warnings; execution coverage 323 of 323 functions); validate.py 625,689 rows, 0 problems. Next, 2d-4b: the writer.

## 92. M4 part 2d-4b: the Touchstone writer at its range edges (2026-11-15)

- **The finding**: `try_write_touchstone` writes MA's (|S|, arg S°) and DB's (20·log₁₀ max(|S|, 2^-1022), arg S°)
  with |S| = √(re² + im²) squared raw, which is BjComplex's abs without its row-14 edge. It refused 10¹⁶⁰ in MA
  (re² past DBL_MAX) and DBL_MAX(1 + j) in DB, where |S| = 2^1024.5 has a finite 6168.1 dB; it wrote 10⁻²⁰⁰ as 0 in
  MA and 10⁻³⁰⁰(1 + j) at DB's clamp (−6153.05 dB for −5997 dB); and an argument below 2^-1022 lost bits before its
  product by 180/π (10³⁰⁰ + j10⁻¹⁰: 5.729577951308214·10⁻³⁰⁹° for 5.729577951308234·10⁻³⁰⁹°).
- **ADR 0002 row 18, amended** (`src/touchstone.cyr`): Rust's chain stays inside `_bj_cx_plain`'s box (each part 0
  or within 2^±511); elsewhere `_bj_ts_mag_ux` takes |S| as BjComplex's abs at its edge (re² + im² to 53 bits with an
  unbounded exponent, the root rounded once). In DB an |S| whose 53-bit root reaches 2^1024 has its log₁₀ to 53 bits
  (`_bj_ux_log10`), then Rust's 20·; below DBL_MAX the rounded |S| takes Rust's clamp, log₁₀ and 20·. Where atan2(im,
  re) is at most 2^-1022 with im ≠ 0, `_bj_ts_deg_ux` keeps it to 53 bits (`_bj_rf_atan2_ux`: im/re there) and
  rounds its product by 180/π once. MA still refuses an |S| past DBL_MAX; DB now writes every finite entry.
  `_bj_cx_plain`, `_bj_cx_ux_sumsq` and `_bj_rf_atan2_ux` are public for the writer, as `_bj_rf_require` is.
- **log₁₀ with an unbounded exponent** (`src/crmath_ux.cyr`): `_bj_ux_log10(s, E)` multiplies `_bj_ux_log_fx`'s
  |ln(s·2^E)| (364-bit fixed point) by 1/ln 10 (`_bj_ux_ivln10`, 14 limbs in `src/crmath_ux_data.cyr`) and rounds
  once. The constant is generated: crtables.py computes ln 10 = 6·atanh(1/3) + 2·atanh(1/9) to 1,024 bits, validate.py
  checks it against MPFR's (`mpfr_ref.log_bits`), and render_ux writes the truncated reciprocal. `cr_verify.py`'s
  `ux_log10` family: 100,000 arguments (half of them in the writer's binade, s·2^1024 with s < √2) equal to MPFR's
  widest-range log10 on x86_64 and aarch64.
- **The grid**: 428 rows appended to `rf.tsv` (30,714), all `try_write_touchstone` in MA and DB: |S| past DBL_MAX in
  each quadrant, squares past DBL_MAX or below 2^-1022 with |S| in range, subnormal parts, the box's 2^±510 to 2^±513,
  arguments at most 2^-1022 with the knife edge of atan2 = 2^-1022 and of degrees = 2^-1022, a 2-port and a 3-port
  mixing them, and 150 seeded entries.
- **validate.py** models the writer's pair (`rf_ts_pair_ux`) under `RF18`. `rf.adr.tsv` maps **213 write rows** to
  row 18 (10,508 in all), all passing: **165** were refusals in Rust, 48 had their text change, 18 of them on 1c-2's
  edge grid. One corrected-oracle row, a 2-port in DB with an entry of 10⁻¹⁷⁰ (−3,400 dB, which Rust's underflowed
  square put at the clamp), is now a row-18 mapping: `rf.cr.tsv` 86 rows.
- **Tests**: `tests/rf.tcyr`'s `cyrius_touchstone_write_edges` (6 assertions) and `cyrius_touchstone_write_reasons`
  amended: 10²⁰⁰ + j writes 4,000 dB, where it pinned Rust's refusal, and the refusal text is pinned on an |S| past
  DBL_MAX in MA (261 in all).
- **Benchmarks** (a focused interleaved A/B against HEAD, six rounds at load 2.9 to 3.8 with other repos' jobs and a
  core dump running): `rf/write_touchstone_s2p_201` medians 4.22 → 4.24 ms in MA and 5.26 → 5.29 ms in DB, within
  the noise; the guard is one exponent test per entry and one compare per argument.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 535,888 ok, 0 failed, 0
  measured, 0 pending, 90,229 accepted (ADR 01: 30,858; ADR 0002: 59,371)** over 626,117 rows; **all 416 digests
  equal**, HEAD's unchanged but `try_write_touchstone`'s. `rf.tcyr` passes on all three (261). The gate is clean (no
  warnings; execution coverage 323 of 323 functions); validate.py 626,117 rows, 0 problems. **2d-4 is done.** Next,
  2d-5: the array factors' ψ past DBL_MAX.

## 93. M4 part 2d-5: the array factors at their range edges; M4 done (2026-11-15)

- **The finding**: the linear factor's ψ = kd·cos θ + β and Nψ/2, and the circular factor's ψₙ = ka·sin θ·cos(φ −
  φₙ) + β·n, overflow to ±∞ where kd, a product, the sum or β·n passes DBL_MAX, and Rust takes the sine of ∞: NaN,
  for finite inputs. k = 10³⁰⁰ and d = 10²⁰ at θ = π/2, where ψ = 6.1·10³⁰³ is itself finite, gave NaN for 0.2674;
  N = 10¹⁰ with ψ = 10³⁰⁰ (Nψ/2 past DBL_MAX) NaN for 1.13·10⁻¹⁰. The planar factor inherits the linear one's.
- **ADR 0002 row 18, amended** (`src/rf.cyr`): the guard is Rust's own NaN, tested once (sin(Nψ/2) in the linear
  factor, Σ cos + Σ sin after the circular loop), with every input finite (a sum of x·0); then `_bj_rf_laf_ux` and
  `_bj_rf_caf_ux` rerun the chain with an unbounded exponent: each product and sum to 53 bits, ψ/2 exact, the sines
  and cosines by src/crmath_vx.cyr's `_bj_vx_sin53` and `_bj_vx_cos53` (Payne–Hanek past DBL_MAX; |Nψ/2| < 2^2112
  and |ψₙ| < 2^2049, within the table's 2^4454), rounded into binary64, then Rust's 10⁻¹⁵ test, quotient, sums and
  root. Below the range nothing changes: an intermediate under 2^-1022 either makes |ψ/2| < 10⁻¹⁵ (the factor is 1)
  or lies under half an ulp of β, and a tiny ψₙ's sine is lost in the sum beside sr ≈ N. Non-finite inputs keep
  Rust's NaN; the planar factor's cos(acos x) keeps Rust's chain (§0.1).
- **The grid**: 1,502 rows appended to `rf.tsv` (32,216): the linear factor at N from 2 to usize::MAX over six
  (k, d) pairs whose kd, kd·cos θ or ψ passes DBL_MAX (kd past it with ψ finite near θ = π/2 among them), six angles
  and four β; Nψ/2 alone past it; non-finite inputs; the planar factor through it; the circular factor at N from 1
  to 16 where ka·sin θ·cos(φ − φₙ), β·n or their sum passes DBL_MAX; 300 seeded wide values.
- **validate.py** models both chains (`rf_linear` and the circular transcription, MPFR's widest-range sine and cosine)
  under `RF18` for finite inputs. `rf.adr.tsv` maps **1,079 rows** to row 18 (793 linear, 52 planar, 234 circular;
  11,587 in all), all passing and every one NaN in Rust; 9 of them on older grids. `rf.cr.tsv` gains 5 rows of the
  new grid where the correctly rounded primitives differ from glibc's (91).
- **Tests**: `tests/rf.tcyr`'s `cyrius_array_factor_edges` (5 assertions, 266 in all).
- **Benchmarks** (a focused interleaved A/B against HEAD, three rounds at load 2.8): `rf/linear_array_factor_8`
  459–464 → 460–471 ns and the planar factor (4 × 4) 1.66–1.69 → 1.66–1.71 µs, within the noise; the circular factor
  (8 elements) 3.89–3.90 → 3.93 µs, +0.8 % in every round (the loop's end test and its locals).
- **The record at M4's close**: the full suite ran at load 11.6, and the rf section three more times at load 17 to 24
  (other repos' benches, an AV1 encoder and four core dumps): 1.5–4× this Part's A/B figures. So
  `docs/benchmarks-rust-v-cyrius.md` keeps its Cyrius column, the rf rows from 1c-2. Part 2's A/Bs give their costs
  since: `rf/s_cascade` 523 ns → 3.5 µs (2a-2's star product rounded once, §83), `input_impedance_complex` about
  280 ns (§87), `reflection_coefficient` 22–27 ns (§89), `single_stub_match` 2.6–2.8 µs (§90), the parse 570–590 µs
  (§91), the write about 4.2 ms (§92), `linear_array_factor_8` 460–470 ns (here).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **17 modules, 0 failing; 536,306 ok, 0 failed, 0
  measured, 0 pending, 91,313 accepted (ADR 01: 30,863; ADR 0002: 60,450)** over 627,619 rows; **all 416 digests
  equal**, HEAD's unchanged but the three array factors'. `rf.tcyr` passes on all three (266). The gate is clean (no
  warnings; execution coverage 323 of 323 functions); validate.py 627,619 rows, 0 problems. **2d is done: Part 2 is
  complete, and with it M4.** Next: M5, the FDTD core.

## 94. M5 part 1a: the 1D grid and the trajectory rows (2026-11-15)

M5 ports fdtd.rs's FDTD core (1,209 of its code lines; its 101 tests outside the 2D, 3D and CPML ones) in Parts, as
M3 and M4 did. Part 1 is the API
in bites: 1a the 1D grid; 1b the pulses and TF/SF (`src/fdtd_source.cyr`), with the three Mur-absorption tests that
drive a Gaussian pulse; 1c the lossy and dispersive materials and the subpixel averages (`src/fdtd_media.cyr`); 1d
`NonUniformGrid`, `fourth_order_derivative`, `CwAccumulator` and `find_eigenmodes` (`src/fdtd_util.cyr`); 1e NF2FF's
phasor form (`src/nf2ff.cyr`); 1f the JSON of the serde types. Part 2 is the range-edge review. `nf2ff_2d` reads an
`Fdtd2d`, so it lands with the 2D grid in M6 Part 1, and so does the `fdtd/nf2ff_2d_80x80_180` bench (the plan's
"NF2FF on 80 × 80 × 180 angles" with it).

- **`src/fdtd1d.cyr`**: `BjFdtd1d`, with `BjFdtd1d_new`, `BjFdtd1d_free` and Rust's methods `set_permittivity`,
  `set_permeability`, `add_source`, `step_once`, `run`, `time` and `total_energy` (and the private
  `update_abc_coefficients`).
  - **Storage** (ADR 0002 row 2, amended): the six arrays (E, H, ε_r, μ_r and the two update coefficients) and a scratch
    row are one anonymous mmap (plan §9.5; VirtualAlloc on PE, zeroed by the kernel either way). The fields `e_field`
    and `h_field` are the arrays' addresses and `num_cells` their length; nothing can resize them, so Rust's
    `active_cells` and the common-length guards in `step_once` and `update_abc_coefficients` go, and
    `test_fdtd1d_inconsistent_public_fields_no_panic` is n/a.
  - **Allocation** (row 9, amended): a grid past 2^63 bytes (56 a cell), or one whose mapping can't be had, is
    `Err(Allocation)` "a 1D FDTD grid of {n} cells", where Rust's `vec!` aborts.
  - **usize**: start, end and the cell are read unsigned, a negative word being past 2^63. An end past the grid
    clamps, a start past it sets nothing, and a cell past it is ignored.
  - **`run`'s `Option<usize>`** is `BJ_FDTD_NO_SOURCE` (usize::MAX). Rust's `Some(c)` with c outside the grid injects
    nothing either, so it is `None` exactly; the port skips that sine, which no output reads.
  - **Rust's chains as written**:
    - dt = dx/(2c), and the coefficients dt/(ε₀·dx)/ε_r and dt/(μ₀·dx)/μ_r;
    - the Mur coefficient c·dt/√(ε_rμ_r), then (v·dt − dx)/(v·dt + dx);
    - the updates in Rust's order, where the right Mur update reads the E₀ just written (two cells pin it);
    - `run`'s sin((2π·f)·t), correctly rounded;
    - the energy summed from −0.0 (p11).
  - **The kernel**: H[i] += c_h·(E[i+1] − E[i]) and E[i] += c_e·(H[i] − H[i−1]) run as the plan's `row_diff` on the
    packed `f64v_sub`, `f64v_mul` and `f64v_add` (§9.3, p10), through the scratch row, then a scalar tail. Each lane
    rounds Rust's three operations, so the bits are Rust's.

    **Probe p19** (new): on 6.7.18 `f64v_sub` and `f64v_mul` still write one word past `dst` at every odd count (2
    lanes), on x86_64 and aarch64, while every word they compute equals the scalar loop's. `f64v_add` never wrote
    past. The packed verbs therefore take the largest multiple of 4, and plan §10's defect stays open.
- **The trajectory format** (plan §6 step 5, new):
  - **A row is a script**: the grid (`i:<n> f:<dx>`), then operations, each an opcode and its arguments
    (set_permittivity, set_permeability, add_source, k × step_once, run with `n:` for `None`).
  - **The outputs**: the constructor's error; or each failed operation's error, then the state. The state is the step,
    dt, time, energy and the Mur coefficients, the digests of all six arrays, and E and H at the probe cells 0, 1,
    n/2, n − 2 and n − 1.
  - **Frames**: a scenario's frames are scripts that each extend the previous one by an advancing operation. The
    harness replays every row on a new grid, and `validate.py` continues from the previous row's grid, so it runs
    each scenario once.
  - **Private state**: the oracle reads Rust's private arrays and Mur coefficients from the grid's `Debug` text, whose
    shortest digits parse back exactly.
  - **The digest** is the harness's cross-target step, `_par_dg_step`: FNV-1a, then y ⊕ (y >> 29), with every NaN
    folded to 0x7ff8….
  - **Class T** is compared as class A, every digest and probe exact. `parity_core.cyr` accepts it.

  The first port run failed every digest. Rust's and Python's digests shifted logically, while Cyrius's `>>>` is the
  arithmetic shift (its `>>` is the logical one, the reverse of Java's). The two references now take the harness's
  form, and HEAD's cross-target digests are unchanged.
- **New golden file `fdtd.tsv`, 627 rows over 6 functions**, every one class T:
  - `new` over cell counts 0 to 65,537 against every special spacing (subnormal, 10^±300, DBL_MAX among them);
  - `set_permittivity` and `set_permeability` over the specials and every kind of range (empty, reversed, clamped, past
    the grid, usize::MAX), including those touching the Mur coefficients' nodes, plus a failed setting after good
    ones;
  - `add_source` over the specials at every edge cell, outside the grid and repeated;
  - `step_once` trajectories at frames {1, 2, 10, 100, 1000} (one to 5,000) in vacuum, ε_r = 4, μ_r = 9 and layered
    media, on two- to four-cell grids, with non-finite, huge, subnormal and signed-zero impulses, and with spacings
    whose coefficients under- or overflow;
  - `run` with the Rust tests' and benches' drives, with `None` and out-of-grid cells, and with special and huge
    frequencies and amplitudes;
  - 16 seeded scenarios.

  `validate.py`'s independent transcription (`Fd1`, `FdScripts`) agreed with every row on its first run.
  **`fdtd.cr.tsv` has 3 rows** where glibc's sin misrounds somewhere in a drive and the whole state after it follows:
  the 1,000-step frame of the Rust test's 1 GHz drive, and a seeded scenario's 100- and 1,000-step frames.
- **Rust tests**: 21 of the 1D grid's 24 are ported. `tests/fdtd1d.tcyr` has 68 assertions with the Cyrius-only groups:
  - the error records' text, and the cell count checked before dx, as Rust does;
  - the allocation refusal at usize::MAX and 2^59 cells;
  - `BjFdtd1d_free` and a freed grid;
  - usize ranges;
  - the no-source sentinel against an out-of-grid cell and plain steps;
  - the two-cell grid.

  Two tests drive `gaussian_pulse` and follow in 1b; one is n/a. The ledger: 618 ported, 9 n/a, 231 pending of 858.
- **Surface parity**: fdtd is in progress. Fdtd1d's 15 items and the two fourth-order constants are matched (464 items
  in all), and the 114 still to come are listed by name. The API snapshot gains `BjFdtd1d_new` and `BjFdtd1d_free`
  (2,070).
- **Benchmarks** (the 3 twins of part 1a, so 72 of 80 are twinned; load 0.4 to 2.5, interleaved with HEAD, two rounds
  each). Every other row is unchanged, with a median ratio of 1.01 and 0.97 in the two rounds. The packed kernel
  against the scalar loops first written (load about 1):

  | Bench | Rust | Cyrius, scalar | Cyrius | Ratio |
  |---|---|---|---|---|
  | `fdtd/step_100_cells` | 64.09 ns | 1.07–1.10 µs | 249–251 ns | 3.9× |
  | `fdtd/step_1000_cells` | 623.03 ns | 10.4–10.8 µs | 2.11–2.25 µs | 3.4× |
  | `fdtd/run_100_steps_500_cells` | 32.23 µs | 562–579 µs | 146 µs | 4.5× |

  A step is now about 2.1 ns a cell, against p10's 1.26 ns for one row (two rows a step here). The run adds a sine a
  step and a mapping per iteration.

  `docs/benchmarks-rust-v-cyrius.md` takes the three twins' figures, and keeps the rest of its Cyrius column as at
  M4's close (§93).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 536,930 ok, 0 failed, 0
  measured, 0 pending, 91,316 accepted (ADR 01: 30,866; ADR 0002: 60,450)** over 628,246 rows; **all 422 digests
  equal**, HEAD's 416 unchanged and the six Fdtd1d ones the same under the scalar and the packed kernel. `fdtd1d.tcyr`
  passes on all three (68). The gate is clean (no warnings; execution coverage 325 of 325 functions); validate.py
  628,246 rows, 0 problems; the oracle's `--check` and clippy are clean. Next: 1b, the pulses and TF/SF.

## 95. M5 part 1b: the pulses and the TF/SF source (2026-11-15)

- **`src/fdtd_source.cyr`**:
  - **The pulses**: `bj_gaussian_pulse`, `bj_ricker_wavelet`, `bj_differentiated_gaussian` and
    `bj_modulated_gaussian_pulse`, Rust's chains on the correctly rounded exp and sin: u = (t − t₀)/τ, exp((−u)·u),
    (1 − 2u²)·exp(−u²), ((−2u)/τ)·exp((−u)·u), and the carrier's ((2π·f₀)·(t − t₀)). A τ of 0 (either sign) or NaN
    gives 0.0, as Rust's.
  - **`BjTfsfSource`**: `BjTfsfSource_try_new`, `BjTfsfSource_new`, `BjTfsfSource_free` and Rust's `step`,
    `incident_e` and `incident_h`. The auxiliary grid's `inc_e`, `inc_h` and a scratch row are one mapping, and `nx` is
    their length (ADR 0002 row 2, amended); `step`'s common-length guard goes. `boundary_x` and `source_cell` stay
    writable, as Rust's public fields, so their guards stay: the source cell is clamped into 1..nx − 2 (unsigned: 0 to
    1, usize::MAX to nx − 2), and the incident fields read 0.0 outside the grid.
  - **ADR 0002 row 3, amended**: Rust's `new` clamped exactly the inputs `try_new` refuses (nx < 3, boundary_x outside
    1..nx, a dx or dt not finite and positive, c·dt > dx·(1 + 10⁻¹²)). `BjTfsfSource_new` is therefore `try_new`: it
    refuses them with the same error, and Rust's four warnings go.
  - **usize**: nx and boundary_x are compared unsigned, as Rust's. A huge nx with a huge boundary below it passes the
    checks and then fails to allocate (`Err(Allocation)`, row 9 amended, where Rust's `vec!` aborts).
  - **The kernel**: the two row updates run on `f64v_sub`, `f64v_scale` and `f64v_add` through the scratch row, with a
    scalar tail. Probe p19 now covers `f64v_scale`: bit-identical to Rust's `k·(a − b)`, and one word past `dst` at
    every odd count, as `f64v_mul` (x86_64 and aarch64), so the verbs again take the largest multiple of 4.
- **The grid**: 16,901 rows appended to `fdtd.tsv` (17,528).
  - **The pulses**, 16,447 class-C rows: every special for t, t₀ and τ (and for the frequency on three typical
    pulses); the Rust tests' points; 2,000 seeded pulses (τ from 1 fs to 1 µs, t within ±8τ of t₀, carriers to
    1 THz); 1,000 seeded wide values; and the range's edges: u past DBL_MAX or subnormal, u² past either end, exp's
    underflow near u² = 745, and carriers of 10³⁰⁰.
  - **The 1D grid driven by a Gaussian** (opcode 6), 39 rows: the Rust tests' Mur absorption in vacuum and ε_r = 4 at
    frames through the source and the boundaries' reflection; a pulse in vacuum at the frames; τ = 0 and NaN; t₀ NaN
    and +∞; τ of 10^±300.
  - **TF/SF**, 415 rows over `try_new`, `new`, `step` and the written fields:
    - every constructor case: nx 0 to 100, the boundary 0, 1, at and past nx, 500 and usize::MAX, the specials for dx
      and dt, dt at dx/(2c), dx/c and 1.01·dx/c, and dt straddling dx·(1 + 10⁻¹²)/c three ulps each way;
    - a huge nx refused before any allocation;
    - the clamped grids Rust's test stepped 2,000 times;
    - the Rust tests' drives at their frames, boundaries 1 to 3, the magic time step (S = 1, a Mur coefficient of 0),
      three- and four-cell grids, unit, non-finite, huge and subnormal impulses;
    - the boundary and source cell written inside, at the edges, outside and at usize::MAX (Rust's mutated-fields
      test among them);
    - 12 seeded sources to 1,000 steps.
- **validate.py**: its transcription (`fd_pulse`, `Tf`, `tf_check`, `tf_sanitize` and the TF/SF scripts, cached per
  scenario like the 1D ones) agreed with every row on its first run.
  - `fdtd.cr.tsv` gains 8 pulse rows where glibc's exp or sin misrounds (11 in all). No driven trajectory needed one.
  - **`fdtd.adr.tsv`, new: 105 rows** map the `TfsfSource::new` rows Rust clamped to `try_new`'s error and text (ADR
    0002 row 3). The other 29 `new` rows are unchanged.
- **Rust tests**: all 12 pulse and all 10 TF/SF tests are ported in `tests/fdtd_source.tcyr`, 64 assertions with the
  Cyrius-only groups (the error records, the unsigned boundary at usize::MAX, `BjTfsfSource_free`, the source cell's
  clamps, the pulses' specials).
  - Two tests assert the refusal (row 3) where Rust asserted the clamped grid: `test_tfsf_new_sanitizes_without_panic`
    checks each of its seven cases is refused with `try_new`'s code, and `test_tfsf_boundary_at_zero` checks both
    constructors refuse boundary_x = 0.
  - `test_tfsf_mutated_public_fields_no_panic` keeps all but its `inc_h.truncate(2)` (row 2), and adds the clamped
    source cell's injection.
  - `tests/fdtd1d.tcyr` ports the two Mur absorption tests (70 assertions). §94 counted them as three: the third name
    was their helper, `fdtd1d_reflected_energy_fraction`.

  The ledger: 642 ported, 9 n/a, 207 pending of 858.
- **Surface parity**: the four pulses and `TfsfSource` with its fields and methods are matched (478 items), and 104
  items are still to come. The API snapshot gains seven functions (2,077).
- **Benchmarks**: no Rust bench covers the pulses or TF/SF, so there are two Cyrius-only ones (`EXTRA`): a Gaussian
  pulse's value, and a 200-cell TF/SF step under a Gaussian drive. HEAD's 94 rows were unchanged against the new tree,
  with median ratios of 1.00 and 1.00 over two interleaved rounds at load 2.1–2.7. Two rows were noise:
  `step_1000_cells` once at 3.37 µs at HEAD, and both 1D steps once at 476 ns and 3.68 µs at load 1.9.

  | Bench | Scalar loops | Packed kernel |
  |---|---|---|
  | `fdtd/gaussian_pulse` | — | 122–127 ns |
  | `fdtd/tfsf_step_200` | 2.37–2.45 µs | 846–852 ns (−65 %) |

  The TF/SF step's figure includes its Gaussian, about 120 ns. That leaves about 3.6 ns a cell, against the 1D
  grid's 2.1.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 553,718 ok, 0 failed, 0
  measured, 0 pending, 91,429 accepted (ADR 01: 30,874; ADR 0002: 60,555)** over 645,147 rows; **all 430 digests
  equal**. Against 1a's, only `Fdtd1d::step_once`'s changed (its 39 new rows), and eight are new. `fdtd_source.tcyr`
  (64) and `fdtd1d.tcyr` (70) pass on all three. The gate is clean (no warnings; execution coverage 332 of 332
  functions); validate.py 645,147 rows, 0 problems; the oracle's `--check` and clippy are clean. Next: 1c, the lossy
  and dispersive materials and the subpixel averages.

## 96. M5 part 1c: the FDTD media and the subpixel averages (2026-11-15)

- **`src/fdtd_media.cyr`**:
  - **The materials**: `BjLossyMaterial`, `BjDrudeMaterial`, `BjLorentzMaterial` and `BjDebyeMaterial`, small structs
    passed by value. Each has an unchecked `new` that stores what it is given, as Rust's does, and a validating
    `try_new` into an out-pointer, with Rust's checks in Rust's order:
    - ε_∞ (or ε_r) finite and positive;
    - ω_p, γ, Δε and σ finite and non-negative (−0 passes, as Rust's `>= 0.0`);
    - ω₀ and τ finite and positive;
    - Debye's ε_s ≥ ε_∞, checked after both permittivities and before τ.
  - **The coefficients**: `update_coefficients` (C_a, C_b) and the three `ade_coefficients` return Rust's tuples
    (`: (f64, f64)` and `: (f64, f64, f64)`), each a chain of `+ − × ÷` written in Rust's order. For example
    Lorentz's β is −(1 − γΔt/2) negated before the division, and Drude's β is ((ε₀·ω_p)·ω_p)·Δt/(2 + γΔt).
  - **The presets**: `gold`, `silver` and `water` take their fields from the new generated `src/fdtd_media_data.cyr`
    (`scripts/oracle/src/t_fdtd_media.rs`, the oracle's own values as positional decimals; the oracle's `--check`
    regenerates it).
  - **The subpixel averages**: the arithmetic mean (checked), the harmonic mean (unchecked: NaN for an invalid
    permittivity; f clamped into [0, 1], NaN taken as 0, −0 kept as Rust's `clamp` keeps it), and its checked form.
  - **ADR 0002 row 5, amended**: only Rust's docs stated the media's ε(ω), in e^{−iωt}. The port's docs state it in
    e^{+jωt}, and no number changes, since the coefficients are time-domain. The ADE response tests already drove in
    e^{+jωt}.
- **The grid**: 68,381 class-A rows appended to `fdtd.tsv` (85,909, 23 MB):
  - `new` and `try_new` over every special crossed (Lorentz's 10⁴ for `try_new`; its `new` at γ = 1), the Rust tests'
    cases, Debye's ε_s ≥ ε_∞ straddled one ulp each way, and −0 for a non-negative field;
  - the three presets;
  - the four coefficient sets over every special for each argument (10⁴ each), the Rust tests' cases, gold, silver
    and water across FDTD time steps, 1,000 seeded typical media and 500 seeded wide ones;
  - the three subpixel averages over every special crossed, fill fractions at 0, −0, 1 and one ulp past each end,
    the Rust tests' cases, and seeded media to 10^±300.

  `validate.py`'s transcription (`fm_eval`) agreed with every row on its first run, and so did the port. No row needs
  a correction or a mapping.
- **Rust tests**: 12 of the 13 material tests and all 7 subpixel tests are ported in `tests/fdtd_media.tcyr`: 60
  assertions with the Cyrius-only groups (the error records' text and order, the presets' fields, `new` storing
  unchecked values and its κ = ∞ at ε_∞ = 0, the subpixel's −0). `test_lossy_material_serde` follows with the JSON
  in part 1f. The ledger: 661 ported, 9 n/a, 188 pending of 858.
- **Surface parity**: the four materials with their fields and methods, and the three subpixel functions, are matched
  (512 items). Their serde impls stay listed for 1f. The API snapshot gains 14 functions (2,091).
- **Benchmarks**: no Rust bench covers the media. A Cyrius-only `fdtd/media_coefficients` (`EXTRA`) takes one set of
  each material's coefficients and a harmonic average in 83–84 ns. HEAD's 96 rows were unchanged against the new
  tree, with median ratios of 1.01 and 1.01 over two interleaved rounds at load 0.7–2.8.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 622,099 ok, 0 failed, 0
  measured, 0 pending, 91,429 accepted (ADR 01: 30,874; ADR 0002: 60,555)** over 713,528 rows; **all 448 digests
  equal**, 1b's 430 unchanged and 18 new. `fdtd_media.tcyr` passes on all three (60). The gate is clean (no
  warnings; execution coverage 346 of 346 functions); validate.py 713,528 rows, 0 problems; the oracle's `--check`
  and clippy are clean. Next: 1d, `NonUniformGrid`, `fourth_order_derivative`, `CwAccumulator` and `find_eigenmodes`.

## 97. M5 part 1d: the non-uniform grid, the fourth-order derivative, the CW accumulator, the eigenmodes (2026-11-15)

- **`src/fdtd_util.cyr`**:
  - **`BjNonUniformGrid`**: its axes are addresses and lengths in one freelist block (ADR 0002 row 2, amended).
    - `uniform`, unchecked as Rust's, and `try_uniform`;
    - `graded` into a caller buffer: Rust's (d_max/d_min)^(1/(n − 1)) on the correctly rounded pow, then `powi` with
      i as an i32, as Rust truncates it;
    - `validate`, axes in Rust's order;
    - `min_spacing`, folded from +∞ with Rust's NaN-ignoring min;
    - `try_courant_dt` and `courant_dt`: d_min/((2·√3)·c), 0.0 with a warning for an invalid grid;
    - `extents` as a tuple, each axis summed from −0.0, so an empty axis is −0.0, as Rust's `Sum`.

    Rust builds a grid from its public vectors (`graded`'s among them); the port's `BjNonUniformGrid_from_axes` copies
    a caller's arrays.
  - **`bj_fourth_order_derivative`** writes into the caller's `out`, which must not overlap the field. It is Rust's
    stencil: (1/dx)·(c₁Δ₁ + c₂Δ₃) inside, (f[i+1] − f[i])·(1/dx) at the edges, and 0.0 padding. Fewer than two
    samples, or an invalid dx (with a warning), give zeros.
  - **`BjCwAccumulator`**: real, imag and a scratch row are one mapping.
    - `new`, unchecked with Rust's warning, and `try_new`;
    - `accumulate`: the correctly rounded sin_cos of ω·(step·dt), then real += f·(cos·dt) and imag −= f·(sin·dt). It
      runs on `f64v_scale`, `f64v_add` and `f64v_sub` over a multiple of 4 (probe p19), with Rust's common prefix and
      its one-time warning on a length mismatch;
    - `try_accumulate`, which refuses the mismatch;
    - `amplitude` √(re² + im²) and `phase` atan2(im, re), into the caller's `out`.

    The phasors are e^{+jωt} already (ADR 0002 row 5).
  - **`bj_find_eigenmodes`** into a `BjEigenmodeResult`, three arrays in one block that the caller frees:
    - Rust's checks and texts in Rust's order;
    - a Nyquist warning;
    - `dft_magnitude`'s phasor recurrence, reseeded from the correctly rounded sin_cos every 1,024 samples (public as
      `_bj_fu_dft_magnitude` for Rust's test of it);
    - the strict local maxima above 10 % of the strongest, with `half_power_offset`'s interpolated Q.

    A failed allocation is `Err(Allocation)` (row 9, amended).
- **The grid**: 5,965 rows appended to `fdtd.tsv` (91,874, 27 MB):
  - `uniform` and `try_uniform` over 64 sizes and 12 spacings;
  - 1,173 gradings (class C): cell counts against every special pair, ratios of 10^±600, and seeded grades;
  - 235 grids, each queried five ways: uniform grids, an empty axis, a special at each axis's ends and middle, the
    Rust tests' grids, graded axes, extremes and seeded grids;
  - 1,097 fourth-order derivatives: x, x², x³, specials inserted, and seeded fields, against special and typical dx;
  - 811 CW scripts: the constructors over the specials, the Rust test's mismatches, generated fields at the packed
    kernel's split sizes to 1,000 samples, the step at usize::MAX, special samples, frequencies and steps, and seeded
    runs;
  - 173 eigenmode searches: the Rust tests' validation; special dt, ranges and samples; and damped modes from exact
    recurrences (x_{k+1} = a·x_k + b·x_{k−1}, so no input depends on a libm). The modes run at the bench's 4096 × 256,
    in pairs, merged, at the band's edges, past Nyquist, and seeded.

  The rows follow the oracle README's `black_box` rule. Wrapping 1b's pulse inputs in it changed no row.
- **validate.py**: its transcription (`fu_eval`, `Cw`, `fu_eigen`) agreed with every row once the CW rows' builder was
  fixed. The first emitter built the special-frequency `accumulate` scripts with the unchecked `new`, so those rows
  are now named `CwAccumulator::new`, the name telling the model which constructor built them. `fdtd.cr.tsv` gains 7
  CW rows where glibc's sin_cos or atan2 misrounds (18 in all). No grading or eigenmode search needed a correction.
- **The port** passed every row on its first run. The harness's two 64 KB static scratch rows tripped the compiler's
  static-data advisory (a warning), so they are allocated on first use instead.
- **Rust tests**: all 25 are ported in `tests/fdtd_util.tcyr` (7 grid, 6 fourth-order, 5 CW, 7 eigenmode), with 73
  assertions in all including the Cyrius-only groups: the error records, `from_axes`, the extents' −0.0, the free
  functions, and a result with no mode. `test_cw_accumulator_length_mismatch` passes lengths and keeps its own copy
  where Rust clones. The ledger: 686 ported, 9 n/a, 163 pending of 858.
- **Surface parity**: the three types and two functions are matched (540 items), and their serde impls stay listed for
  1f. The API snapshot gains 11 functions (2,102).
- **Benchmarks** (the twin of part 1d, so 73 of 80 are twinned): `fdtd/find_eigenmodes_4096x256` takes 8.37–8.48 ms at
  load 1.8, against Rust's 2.12 ms (3.9×). HEAD's 97 rows were unchanged, with median ratios of 1.00 and 1.00 over two
  interleaved rounds at load 0.5–2.1; the 1D step rows swung both ways at HEAD too.

  First written, the twin took 8.8–9.1 ms. Updating the phasor's imaginary part in place, one local fewer a sample,
  saves about 5 % with the same bits. The rest is the per-sample chain (two products into the sums, four into the
  phasor), about 2 ns a sample in Rust and about 8 here, which a two-lane form across frequencies could halve. That is
  a finding for M8's performance pass, not this bite.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 628,057 ok, 0 failed, 0
  measured, 0 pending, 91,436 accepted (ADR 01: 30,881; ADR 0002: 60,555)** over 719,493 rows; **all 462 digests
  equal**, 1c's 448 unchanged and 14 new. `fdtd_util.tcyr` passes on all three (73). The gate is clean (no warnings;
  execution coverage 357 of 357 functions); validate.py 719,493 rows, 0 problems; the oracle's `--check` and clippy
  are clean. Next: 1e, NF2FF's phasor form.

## 98. M5 part 1e: NF2FF's phasor form (2026-11-15)

- **`src/nf2ff.cyr`**: `bj_nf2ff_2d_phasor` into a `BjNf2ffResult` (pattern and angles in one freelist block, freed by
  `BjNf2ffResult_free`). It is Love's equivalence on a TM-mode Yee grid's rectangular contour through the E_z nodes
  `box_margin` cells inside:
  - J = n̂×H from the mean of the two Yee H samples beside the contour, rotated by e^{+jω·dt/2} onto E's time;
  - M = −n̂×E;
  - trapezoid weights, half at each side's corners;
  - e^{jkρ̂·ρ'} advanced by a constant phasor along each side;
  - the pattern √(k/8π)·|η·N_z − L_φ|, with L_φ = −L_x sin φ + L_y cos φ.

  fdtd.rs's own `Cplx` arithmetic is written out on f64 locals in Rust's order. sin_cos and hypot are correctly
  rounded. The phasors are e^{+jωt}, as fdtd.rs's were (ADR 0002 row 5).
  - **Checks**, Rust's texts in Rust's order:
    - the accumulators agree on ω and dt within 10⁻¹² (relative, with Rust's NaN-ignoring max);
    - at least one angle;
    - the margin in 1..min(nx, ny)/2, compared unsigned;
    - dx, then the frequency ω/2π, finite and positive;
    - each accumulator holds nx·ny samples, the product overflowing 2^64 refused as `checked_mul`'s None;
    - dt finite.

    The last is unreachable: an infinite or NaN dt already fails the agreement test, since inf − inf is NaN. It is
    ported as written.
  - **The core** takes the fields as real/imaginary address pairs, an imaginary part of 0 reading as zero. M6's
    `nf2ff_2d` on an `Fdtd2d` snapshot (the 6 remaining NF2FF tests and the `fdtd/nf2ff_2d_80x80_180` bench) will
    reuse it with rotation 1.
  - **A finding, kept bit-identical**: Rust's samples also carry the zero magnetic component of each side (M_y on
    horizontal sides, M_x on vertical) and multiply it by the phasor at every angle. Those products are ±0, or NaN
    only where the phasor is NaN, and then N_z is NaN as well. So they reach neither the pattern's value nor its
    NaN-ness: hypot ignores a zero's sign, and a NaN N_z makes both parts of η·N_z − L_φ NaN either way. The port skips
    them, a quarter of the inner loop's complex products. The 161 rows, specials and NaN cases among them, gave the
    same digest before and after.
- **The grid**: 161 class-C rows appended to `fdtd.tsv` (92,035). Each accumulator's phasors come from an exact
  generator, re[i] = a·((i·p + q) mod m) − b, with overrides at chosen indices, so no input depends on a libm. The
  rows:
  - Rust's test's validation;
  - the 10⁻¹² agreement straddled three ulps each way in ω and dt;
  - every special frequency, step and dx;
  - the margin at 0, nx/2 and ny/2;
  - nx·ny past 2^64;
  - short accumulators;
  - grid shapes from 4 × 4 to 80 × 80 with 1 to 180 angles;
  - a zero field;
  - NaN, ±∞, 10³⁰⁰, DBL_MAX, the smallest subnormal and −0 at a corner and on each side;
  - frequencies, steps and spacings of 10^±300;
  - 20 seeded transforms.

  `validate.py`'s transcription (`nf_eval`, `nf_core`) agreed with every row on its first run, once its primitive
  sets gained hypot. `fdtd.cr.tsv` gains 13 rows where glibc's sin_cos or hypot misrounds somewhere across the angles
  (31 in all). The port passed every row on its first run.
- **Rust tests**: 4 of the 10 NF2FF tests are ported in `tests/nf2ff.tcyr`, 16 assertions with the Cyrius-only groups
  (the error records, nx·ny past 2^64, the angles, a zero field, the free):
  - `test_nf2ff_phasor_line_source`: an off-centre line current's exact phasors, isotropic within 1 % at all 72
    angles;
  - `test_nf2ff_phasor_plane_wave_extinction`: under 2 % of one side's coherent scale;
  - `test_nf2ff_phasor_validation`;
  - `test_bessel_reference_helpers`: Rust's trapezoid J_n and Simpson Y_n against A&S Table 9.1 within 10⁻⁹.

  The six that drive an `Fdtd2d` follow in M6. The unit takes about 4 s, mostly the line source's Bessel evaluations
  near the contour. The ledger: 690 ported, 9 n/a, 159 pending of 858.
- **Surface parity**: `Nf2ffResult` and `nf2ff_2d_phasor` are matched (544 items); `nf2ff_2d` stays listed for M6. The
  API snapshot gains 2 functions (2,104).
- **Benchmarks**: Rust's NF2FF bench drives the snapshot form on an `Fdtd2d`, so it stays pending for M6. A Cyrius-only
  `fdtd/nf2ff_phasor_80x80_180` (`EXTRA`) sizes the phasor form the same way (an 80 × 80 grid, margin 5, 180 angles):
  1.87–1.90 ms as first written, 1.52–1.54 ms after the skip above (−19 %, load 2.2–2.8), about 30 ns a contour
  sample.

  HEAD's 98 rows gave median ratios of 1.05 and 0.98 over two interleaved rounds. The load climbed from 1.0 to 2.8
  during the first round; no existing path changed.
- **The gate**: its audit first failed on one undocumented function, the field helper in `tests/nf2ff.tcyr` (its
  comment sat above the globals before it); the comment now sits on the function.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 628,205 ok, 0 failed, 0
  measured, 0 pending, 91,449 accepted (ADR 01: 30,894; ADR 0002: 60,555)** over 719,654 rows; **all 463 digests
  equal**, 1d's 462 unchanged and one new. `nf2ff.tcyr` passes on all three (16). The gate is clean (no warnings;
  execution coverage 359 of 359 functions); validate.py 719,654 rows, 0 problems; the oracle's `--check` and clippy
  are clean. Next: 1f, the JSON of the FDTD types.

## 99. M5 part 1f: the FDTD types' JSON; M5 Part 1 done (2026-11-15)

- **`src/json.cyr`**: writers and readers for the FDTD's six serde types. `CpmlParams` and `Mode2d` come with the CPML
  in M6.
  - **The four materials** (`BjLossyMaterial`, `BjDrudeMaterial`, `BjLorentzMaterial`, `BjDebyeMaterial`) are all-f64
    structs, written and read as `Material` is: fields in declaration order, bayan's printer, serde's object and
    sequence forms, every field as written (a Debye ε_s below ε_∞ too, as serde's derive calls no constructor).
  - **`BjNonUniformGrid`** writes its three axes as arrays and reads arrays of any lengths into a new grid (its block
    allocated once the lengths are known, and given back on any failure).
  - **`BjEigenmodeResult`**: its three vectors share one count (part 1d), so a document whose vectors differ in length
    is refused as Deserialize: "an EigenmodeResult whose three arrays differ in length (its shape is fixed: ADR 0002
    row 36)". serde reads such a document; `json.adr.tsv` maps those rows to row 36 (amended).
- **The grid** (`m_json.rs`): 25,336 rows appended to `json.tsv` (91,956, 38.7 MB).
  - The four materials are new `Kind`s through `emit_kind`: seeded and special values, serde's text, Rust's `{}` and
    `{:e}` and bayan's notations, every number spelling in each field position, seeded decimals, the shapes and
    unknown values.
  - The two three-vector types (`emit_fdtd_vecs`) cover every triple of lengths to 3, seeded, special and 300-element
    vectors, the same notations and spellings, and their own shapes: the sequence form, a vector as an object, `null`,
    string and nested-array elements, missing, repeated and wrong-typed fields, and prefixes.
  - The first emitter wrote one ragged `EigenmodeResult` as a writer value, which Rust's public fields allow and the
    port can't build. Writer values are now equal-length, as row 36 keeps the corpus to shapes the port holds.
- **validate.py**: the six types join `JSON_FIELDS`, and `json_row36` learns the eigenmode result's equal lengths. The
  model agreed with every row. `json.cr.tsv` gains 9,515 rows where serde_json's best-effort decimal parse misrounds
  (40,257 in all); `json.adr.tsv` maps 60 ragged eigenmode results (96 in all). ADR 0002 row 12's count rises to 309
  documents the port refuses whole.
- **Tests**: `test_lossy_material_serde` is ported in `tests/fdtd_media.tcyr` (62 assertions), so all 13 material
  tests are ported. `tests/json.tcyr` gains `cyrius_fdtd_json` (86 assertions in all): each type's exact text, round
  trips, the grid's free lengths, the eigenmode result's refusal and its reason, and out left alone. The ledger: 691
  ported, 9 n/a, 158 pending of 858.
- **Surface parity**: the twelve serde impls are matched (556 items). What remains of fdtd.rs (`nf2ff_2d`, the CPML,
  `Fdtd2d`, `Fdtd3d`, `CpmlParams`, `Mode2d`) is M6's. The API snapshot gains 12 functions (2,116).
- **Benchmarks**: two Cyrius-only JSON benches (`EXTRA`, as earlier types had): `json/drude_material_to_json`
  2.46–2.79 µs and `json/nonuniform_grid_from_json` 9.07–11.2 µs (a 10 × 7 × 3 grid). HEAD's 99 rows gave median
  ratios of 1.02, 1.05 and 0.98 over three interleaved rounds. The load ran 1.1 to 4.4, with other repos' jobs
  (`drishti-conform`, a checker) running; no existing path changed.
- **A note for M6 and M7**: `json.tsv` is 38.7 MB against the 48 MiB per-file cap (§81). M6 (`CpmlParams`, `Mode2d`)
  and M7 (the soorat structs) will add JSON rows, so the module will likely split by milestone before then.

**M5 Part 1 is done**: the 1D grid, the sources, the media, the utilities, NF2FF's phasor form and their JSON, in six
bites (§94–§99). fdtd.rs's M5 surface is ported, except what needs an `Fdtd2d` (`nf2ff_2d`, its six tests and its
bench), which moves to M6. Next: Part 2, the range-edge review (§0.1) of the FDTD core.

## 100. M5 part 2a-1a: the media coefficients rounded once (2026-11-15)

M5's Part 2, the range-edge review, runs on the user's answers of 2026-11-15 (plan §0.1): 2a the answers (2a-1a the
media coefficients, 2a-1b the Mur coefficients, 2a-2 the subpixel means, 2a-3 the grading), 2b round once, 2c the
exact thresholds, 2d the range edges, the per-cell step updates among them (the user's answer).

- **The finding**: the lossy C_a = (1 − s)/(1 + s), s = σΔt/(2ε), the Drude and Lorentz α and β on 2 ∓ γΔt and
  2 − ω₀²Δt², and the Debye α on 2τ/Δt − 1 and β on ε_s − ε_∞ cancel where the two sides meet: C_a was 5.6·10⁻⁵ off
  (relative) at s = 1 + 10⁻¹², 3.1·10⁻⁸ at 1 + 10⁻⁹. Each chain also rounds three to six times where once is possible,
  and leaves the range at its edges (an ε_r of 10⁻³²⁰ made ε₀ε_r 0 and C_a NaN where it is −1).
- **ADR 0002 row 40** (`src/crmath_edge.cyr`'s new `_bj_cr_rat`, called from `src/fdtd_media.cyr`): each coefficient is
  rounded once from its exact value, a ratio N/D of sums of products of the inputs and the double ε₀, written as a
  spec string ("+2ea-1bc/+2ea+1bc" is C_a over (ε_r, σ, Δt, Δx, ε₀)). Three paths:
  - the first tier, for operands each 0 or within [2^-120, 2^120]: each term TwoProduct's exact pair (then
    `_bj_e_ed_muld` steps past two factors), the sum by `_bj_e_ed_add`, so a sum of two exact pairs is within 3.02u² of
    its own value however much it cancels; the quotient takes the cascade's rounding test (`_bj_e_ed_quot`);
  - where that test fails, the second tier: each term of at most four factors expanded exactly (up to eight
    TwoProduct components), the sum's components grown into one expansion and compressed (`_bj_e_ed_exact`), the same
    test;
  - then the exact path: signed scaled bignums (`_bj_e_q_*`), rounded by `_bj_e_q_signed_ratio`.

  An exact zero takes IEEE's signs (a one-product numerator's factors', a cancelled sum +0, then D's): Lorentz's
  β at γΔt = 2 is now +0, Rust's −(1 − 1)/2 −0. An output keeps Rust's chain where an operand it reads isn't finite
  or Rust's formula divides by an exact zero (ε_r = 0, Δx = 0, Debye's Δt = 0, a denominator such as 1 + s exactly 0).
- **cr_verify.py**: the rat (both tiers, then the exact path), rat_x (the second tier, then the exact path) and
  rat_slow (the exact path alone) families check the engine against Python's Fraction on the ten coefficients:
  typical media, each cancellation by relative offsets and ulp by ulp, exact zeros, the box edges, signed zeros, any
  finite values and subnormals: **300,000 calls on x86_64 and on aarch64, 0 mismatches**.
- **The grid**: 3,052 rows appended to `fdtd.tsv` (95,087; the earlier body an exact prefix): each cancellation at
  relative offsets 10⁻¹ to 10⁻¹⁶ either side and ulp by ulp around it on three time steps (and three ε_r for the lossy
  material), exact zeros (σΔt = 2ε, γΔt = 2 and 2τ = Δt on powers of two), every operand at the box edges 2^±119 to
  2^±121 and at 2^-1074, 10⁻³¹⁰, DBL_MIN, 10^±300 and DBL_MAX, and 400 seeded media near their cancellation.
- **validate.py**: `fm_port` computes each output from its formula in Fractions (independent of the spec strings) and
  maps **9,775 rows** to row 40 in `fdtd.adr.tsv` (2,819 lossy, 2,831 Drude, 2,224 Lorentz, 1,901 Debye; part 1c's
  rows among them, wherever rounding once changes a bit).
- **Tests**: `tests/fdtd_media.tcyr`'s `cyrius_media_rounded_once` (8 assertions, 70 in all): C_a at s = 1 + 10⁻¹²,
  a typical lossy medium, Drude's β and κ, the exact zeros and their signs, ε_r = 0's NaNs kept.
- **Benchmarks**: the Cyrius-only `fdtd/media_coefficients` (the ten coefficients and a subpixel mean a call) went from
  79 ns to 3.90–3.93 µs, about 0.4 µs a coefficient where Rust's chains take under 8 ns, and the new Cyrius-only
  `fdtd/lossy_coefficients_exact_path` (C_a and C_b at Δt = 10⁻⁴⁰, outside the first tier's box) 10.3–10.4 µs. A first
  version, the expansion tier alone, took 6.2 µs; the double-double first tier takes every typical medium (216–524 ns a
  coefficient, measured one by one). The interleaved A/B against HEAD (two rounds, load 1.1 to 1.8) left the other 100
  rows unchanged (median ratio 1.00 both rounds). The cost is paid at setup: M6's 2D and 3D grids compute C_a and C_b
  per cell when a lossy region is set (Rust's `e_update_coefficients`), so a 10⁶-cell region will take about a second
  there. Recorded as a finding; a per-coefficient fast path without the spec parsing would cut it further.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 637,207 ok, 0 failed, 0
  measured, 0 pending, 110,835 accepted (ADR 01: 40,369; ADR 0002: 70,466)** over 748,042 rows; **all 475 digests
  equal**, HEAD's unchanged but the four media functions' (their new rows and values). `fdtd_media.tcyr` passes on all
  three (70). The gate is clean (no warnings; execution coverage 371 of 371 functions); validate.py 748,042 rows, 0
  problems; the oracle's `--check` and clippy clean. Next: 2a-1b, the Mur coefficients.

## 101. M5 part 2a-1b: the Mur coefficients (2026-11-15)

- **The finding**: both first-order Mur coefficients are (a − b)/(a + b) and cancel where the wave's step meets the
  cell. `TfsfSource`'s (S − 1)/(S + 1), S = c·dt/dx rounded twice, gave 0 at the magic time step dt = dx/c (S = 1 in
  binary64), where the exact (c·dt − dx)/(c·dt + dx) of the given doubles is −7.8·10⁻¹⁸ at dx = 10⁻⁴; near it the
  value was rounding noise. `Fdtd1d`'s (v·dt − dx)/(v·dt + dx), v·dt = c·dt/√(ε_rμ_r), gave 0 at ε_rμ_r = 1/4 and
  dt = dx/(2c) for 2.09·10⁻¹⁷, and over seeded inputs near it was up to 2.6·10²⁰ ulps off.
- **ADR 0002 row 40, amended**:
  - the TF/SF coefficient is rounded once from (c·dt − dx)/(c·dt + dx) (`_bj_cr_rat`'s rational 10; the arguments are
    validated finite and positive, so its denominator is never 0);
  - the 1D coefficient has √(ε_rμ_r) inside, so, as the user's answer has it, where it cancels it takes the
    cancellation-free form (c²dt² − ε_rμ_r·dx²)/(c·dt + dx·√(ε_rμ_r))², the same value: the numerator rounded once
    from the exact one (rational 11), the denominator on Rust's rounded c·dt and √(ε_rμ_r), whose sum has no
    cancellation. "Where it cancels" is v·dt in (2dx/3, 2dx), where the difference at least doubles v·dt's rounding
    error (the stub root's factor 2 would have reached the vacuum default at v·dt = dx/2, whose Rust chain is
    accurate). Outside the band, and where c·dt or dx·√(ε_rμ_r) lies outside [2^-400, 2^400], Rust's chain stays;
    its range edges are part 2d's.
- **Accuracy**: validate.py measures the 1D form against the exact coefficient (80-digit Decimals) on every argument
  set the golden scripts take it on, 119: within 3.21 ulps (gated at 8). A one-off sweep of 171,386 seeded band
  inputs (ε_rμ_r near 1/4 by relative offsets 10⁻¹ to 10⁻¹⁶, across the band, at dt = dx/(2c) and other steps): within
  5.67 ulps, where Rust's chain was up to 2.6·10²⁰; outside the band Rust's chain was within 4.74.
- **cr_verify.py**: the rat families gain rationals 10 and 11 (twelve in all): 360,000 calls on x86_64 and on
  aarch64, 0 mismatches.
- **The grid**: 272 rows appended to `fdtd.tsv` (95,359): `TfsfSource::try_new` at dt = dx/c ulp by ulp and at
  S = 1 − 10⁻ᵏ on five dx, with a Gaussian-driven scenario at S = 1 for each; `Fdtd1d` with ε_r, μ_r or both setting
  ε_rμ_r = 1/4 at relative offsets 10⁻¹ to 10⁻¹⁵ either side and ulp by ulp, the band's edges 1/16 and 9/16 ulp by
  ulp, and four Gaussian-driven scenarios. The generator skips a script an earlier part already emitted.
- **validate.py**: `Fd1` and `Tf` take a `port` flag (the port's Mur coefficients); `fdtd_adr_rows` replays each 1D
  and TF/SF script with and without it on correctly rounded primitives, and maps a row whose state the flag changes
  to row 40 (the corrected oracle skips it): **311 rows** (78 `set_permittivity`, 40 `set_permeability`, 16
  `step_once`, 119 `try_new`, 6 `new`, 52 `step`), part 1's among them. `fdtd.adr.tsv` maps 10,191 rows.
- **Tests**: `tests/fdtd1d.tcyr`'s `cyrius_mur_cancellation_free` (3 assertions, 73 in all: the vacuum default kept,
  ε_r = 1/4 and 0.3) and `tests/fdtd_source.tcyr`'s `cyrius_tfsf_mur_rounded_once` (1, 65 in all).
- **Benchmarks**: the step loops don't read the coefficient, and the interleaved A/B against HEAD (two rounds, load 2.4
  to 2.5; HEAD predates part 2a-1a, so `fdtd/media_coefficients` shows that part's cost again) left the FDTD rows
  unchanged: `fdtd/tfsf_step_200` 862–948 → 871–872 ns, `fdtd/step_100_cells` 254–334 → 264–268 ns,
  `fdtd/run_100_steps_500_cells` 150–156 → 151–153 µs (median ratios 1.01 and 0.99 over 101 rows). The band form is paid
  at setup: a one-off probe of a 20-cell grid's `set_permittivity`, 108 ns on Rust's chain (ε_r = 4), 1.13 µs in the
  band (ε_r = 0.3: two coefficients, the expansion tier) and 2.90 µs at ε_r = 1/4 (the exact path, where the numerator's
  terms cancel past the tiers' budget).
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 637,168 ok, 0 failed, 0
  measured, 0 pending, 111,146 accepted (ADR 01: 40,369; ADR 0002: 70,777)** over 748,314 rows; **all 475 digests
  equal**, part 2a-1a's unchanged but the 1D `set_permittivity`, `set_permeability` and `step_once` and the TF/SF
  `try_new`, `new` and `step` (their new rows and values). `fdtd1d.tcyr` (73) and `fdtd_source.tcyr` (65) pass on all
  three. The gate is clean (no warnings; execution coverage 371 of 371 functions); validate.py 748,314 rows, 0
  problems; the oracle's `--check` and clippy clean. Next: 2a-2, the subpixel means.

## 102. M5 part 2a-2: the subpixel means rounded once (2026-11-15)

- **The finding**: the subpixel means round three and four times, so a uniform cell (ε₁ = ε₂ = ε), which must give
  ε, missed it: the parallel mean f·ε₂ + (1 − f)·ε₁ in 2.5 % of seeded uniform cells (ε = 3 at f = 0.2 gave
  3.0000000000000004), the perpendicular mean 1/(f/ε₂ + (1 − f)/ε₁) in 23 % (12.25 at f = 0.1 gave 12.249999999999998).
  Their intermediates also left the range: the perpendicular mean of DBL_MAX cells was ∞ (f/ε₂ and (1 − f)/ε₁ fell to
  the subnormals, 1/s overflowed), of the least subnormal 0, and the parallel mean of DBL_MAX cells at f = 0.3 was an
  ulp low.
- **ADR 0002 row 41** (the user's answer of 2026-11-15): each mean is rounded once from its exact value on part
  2a-1a's engine: the parallel (1 − f)ε₁ + f·ε₂ (rational 12, written ε₁ − f·ε₁ + f·ε₂ so the first tier's first
  difference stays within its own bound) and the perpendicular ε₁ε₂/((1 − f)ε₂ + f·ε₁) (rational 13). A uniform cell
  gives ε exactly, the mean lies in [ε₁, ε₂] (both are convex combinations, and rounding once is monotone), and no
  intermediate leaves the range. The checks and errors, and the perpendicular mean's NaNs for an invalid ε and its
  clamp of f (NaN as 0), are Rust's; a valid cell's denominator is positive, so neither rational meets a zero.
- **cr_verify.py**: rationals 12 and 13 join the rat families (fourteen in all): equal, adjacent and unrelated media,
  fills in tenths and thirds, at 0 and 1, a subnormal or an ulp from them, and media from 10^±300 to DBL_MAX and the
  subnormals. 420,000 calls on x86_64 and on aarch64, 0 mismatches.
- **The grid**: 4,371 rows appended to `fdtd.tsv` (99,730), each input through the three functions: uniform cells of
  ten common media (1 to 13) at fills in tenths, thirds and quarters, the same with ε₂ one ulp above ε₁, fills of the
  least subnormal, 10⁻³⁰⁰, 10⁻¹⁷ and an ulp or 10⁻¹⁶ from 1 on four media pairs, media of DBL_MAX, DBL_MAX/2, 10^±300,
  DBL_MIN, 10⁻³¹⁰ and the least subnormal crossed at three fills, and 500 seeded cells and 500 seeded uniform cells.
- **validate.py**: `fm_port41` computes each mean in Fractions and maps **2,371 rows** to row 41 (434 parallel, 992
  perpendicular, 945 checked; part 1c's rows among them).
- **Tests**: `tests/fdtd_media.tcyr`'s `cyrius_subpixel_rounded_once` (5 assertions, 75 in all): the two uniform
  cells above, DBL_MAX cells through both means, the least subnormal's.
- **Benchmarks**: the Cyrius-only `fdtd/media_coefficients` (ten coefficients and a perpendicular mean) 3.96 µs at part
  2a-1b → 4.27–4.32 µs, so the rounded-once perpendicular mean costs about 0.35 µs where Rust's chain took a few ns.
  Three interleaved rounds against HEAD (which predates part 2a-1a): the first at load 1.7 to 2.5 left the other rows
  within the noise (median ratio 0.99); from the second on another repo's fuzzer ran five processes (load 8 to 15), and
  the second's new run and the third's HEAD run are inflated (median ratios 1.00 and 0.56), so only
  `media_coefficients`, steady at 4.27–4.32 µs whenever the new side ran, is read from them.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 639,168 ok, 0 failed, 0
  measured, 0 pending, 113,517 accepted (ADR 01: 40,369; ADR 0002: 73,148)** over 752,685 rows; **all 475 digests
  equal**, part 2a-1b's unchanged but the three subpixel functions' (their new rows and values). `fdtd_media.tcyr`
  passes on all three (75). The gate is clean (no warnings; execution coverage 371 of 371 functions); validate.py
  752,685 rows, 0 problems; the oracle's `--check` and clippy clean. Next: 2a-3, the grading.

## 103. M5 part 2a-3: the grading's ends; 2a done (2026-11-16)

- **The finding**: `NonUniformGrid::graded` documents a geometric grading "computed to reach d_max at the end", but
  its last cell d_min·ratio^(n − 1) misses d_max: the ratio (d_max/d_min)^(1/(n − 1)) carries its rounding, and the
  power multiplies it about n − 1 times, so in 97 % of seeded gradings the last cell was off, by up to 270 ulps for
  n ≤ 200 (`graded(10, 1e-4, 1e-3)` ends at 0.0010000000000000005). Rust also takes the power as
  `ratio.powi(i as i32)`: past 2³¹ cells the exponent wraps, so cell 2³¹ + 5 is ratio^-(2³¹ − 5), below d_min.
- **ADR 0002 row 42** (`src/fdtd_util.cyr`): the last cell is exactly d_max (the user's answer of 2026-11-15; the first
  already is d_min), and the interior keeps Rust's chain but never passes d_max (the user's answer of 2026-11-16). The
  power takes the whole usize index on the same square-and-multiply loop (`_bj_powi`, compiler-rt's `__powidf2`), the
  same bits below 2³¹; a Rust fault, fixed whatever the answer. The cell is `_bj_fu_graded_cell`, public for the test
  (no golden row can hold 2³¹ cells: 16 GiB of output). The range edges (d_max/d_min or ratio^i out of range) are part
  2d's.
- **A question asked in the bite**: pinning the last cell alone broke monotonicity where the ratio lies within about
  (n − 1)²·2^-53 of 1: there the grading's steps are below the accumulated rounding of ratio^i, and Rust's interior
  runs past d_max (`graded(10, 1, 1.000000000000001)` ends at 1.000000000000002), so the pinned last cell stepped back.
  Over 20,000 seeded gradings of such ratios (1 ± 10⁻⁴ to 10⁻¹⁵, n = 10, 100, 1,000) 1,473 became non-monotone, where
  Rust's were all monotone; over 20,000 typical ones none. The user's answer of 2026-11-16 (the recommended option):
  hold a cell past d_max at d_max (min for an increasing grading, max for a decreasing one), so the grading is
  monotone (not strictly) and ends at d_max; typical gradings keep Rust's bits.
- **The grid**: 620 rows appended to `fdtd.tsv` (100,350; 35.0 MB, under the 48 MiB cap): twelve gradings
  (increasing and decreasing by a decade, a ratio of 1 and ones an ulp either side of it, 0.1 to 0.3, 2 to 1,024 and
  1 to 3, 10⁻⁹ to 1 and 1 to 10⁹, and two at 10^±300's ends) at n = 2, 3, 4, 7, 10, 33, 100, 257 and 1,000, ratios
  1 ± 10⁻ᵏ (k = 4 to 15) from three d_min at n = 10, 100 and 1,000, whose interior passes d_max, and 300 seeded
  gradings.
- **validate.py**: `fu_graded_port` replays a grading on correctly rounded primitives, holds its cells at d_max and
  pins its last; a row whose state that changes is mapped to row 42 and left out of the corrected oracle: **913
  rows**, part 1d's among them. `fdtd.cr.tsv` keeps its 31 rows (it held no grading).
- **Tests**: `tests/fdtd_util.tcyr`'s `cyrius_graded_rounded_ends` (4 assertions, 77 in all): both ends of the Rust
  test's grading exact, a decreasing grading's last cell and monotonicity, `graded(10, 1, 1.000000000000001)`'s cells
  held at d_max, and cell 2³¹ + 5 of a ratio 1 + 10⁻¹⁰, 1.2395499649960704 (Rust's wrap gave 0.8067).
- **Benchmarks**: no bench calls `graded`, and nothing else changed. The interleaved A/B against HEAD (two rounds at
  load 4 to 15, another repo's fuzzer running throughout) is too noisy to read the other rows (median ratios 0.50 and
  0.94, HEAD's runs the more loaded); `fdtd/media_coefficients` held at 4.27–4.38 µs, as at part 2a-2. The clamp is
  one compare a cell.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 638,875 ok, 0 failed, 0
  measured, 0 pending, 114,430 accepted (ADR 01: 40,369; ADR 0002: 74,061)** over 753,305 rows; **all 475 digests
  equal**, part 2a-2's unchanged but `NonUniformGrid::graded`'s (its new rows and values). `fdtd_util.tcyr` passes on
  all three (77). The gate is clean (no warnings; execution coverage 371 of 371 functions); validate.py 753,305 rows,
  0 problems; the oracle's `--check` and clippy clean.

**Part 2a, the user's answers, is done** (§100–§103): the media and Mur coefficients (row 40), the subpixel means
(row 41) and the grading (row 42). Next, 2b: round once (row 13) for the update coefficients dt/(ε₀ε_r·dx) and
dt/(μ₀μ_r·dx), `courant_dt`'s d/(2√3·c) and `Fdtd1d::time`'s step·dt.

## 104. M5 part 2b: round once (2026-11-16)

- **ADR 0002 row 13, amended**: the standing "round once" decision for M5's single-quotient formulas:
  - the update coefficients dt/(ε₀·dx) and dt/(μ₀·dx) in `Fdtd1d::new` and `TfsfSource::try_new`, and dt/(ε₀·dx·ε_r)
    and dt/(μ₀·dx·μ_r) in `set_permittivity` and `set_permeability` (Rust's (dt/(k₀·dx))/m rounded three times), on
    kimiya's monomial `_bj_k_cr_mono`. The 1D grid's reads its `dt` and `dx` fields (writable), so where either isn't
    finite or dx is 0 Rust's chain stays (`_bj_fd1_coeff`);
  - `NonUniformGrid::try_courant_dt` and `courant_dt`'s d/(2√3·c), √3 as Rust's correctly rounded double;
  - `Fdtd1d::time`'s step·dt past 2⁵³, where Rust's `step as f64` rounds first: the step split exactly as
    s_hi + s_lo (s_hi a multiple of 2¹¹), rounded once on part 2a-1a's engine (rational 14).

  Rust was an ulp off on many grids (dx = 0.032: 188.36515683343495 for 188.36515683343492), and every 1D and TF/SF
  trajectory reads the coefficients, so most change in their last bits. `dt = dx/(2c)` and the TF/SF's Courant
  message dx/c were already one operation. Formulas with a sum, a computed value or a transcendental inside keep
  their own treatment (the run's source t = step·dt and the CW phase are arguments of a sine, as M3's HG phase).
- **A bug the mapping caught**: the 1D coefficient's first guard tested dt·dx for finiteness, which overflows at
  dx = DBL_MAX and kept Rust's chain there; 11 rows failed against validate.py's exact values until each operand
  was tested on its own.
- **The grid**: 2,800 rows appended to `fdtd.tsv` (103,150; 36.6 MB): 200 seeded 1D grids (dx from 10⁻⁹ to 10), each
  also through seeded `set_permittivity` and `set_permeability` media and a seeded TF/SF source, and 200 seeded
  uniform and 200 seeded non-uniform grids through the Courant step and the grid queries.
- **validate.py**: `Fd1` and `Tf` take the port's rows as a set (`FD_PORT`, rows 13 and 40); `fdtd_adr_rows` maps a
  trajectory to row 40 when the Mur coefficient alone accounts for its change, else to row 13, and
  `fu_courant_port` the Courant rows. **1,327 rows** map to row 13 (987 1D scripts, 218 TF/SF scripts, 122 Courant
  steps); `fdtd.cr.tsv` loses 3 `Fdtd1d::run` rows to them (28; the corrected oracle 40,366 in all).
- **Tests**: `cyrius_round_once` in `tests/fdtd1d.tcyr` (3 assertions, 76 in all: the coefficients at dx = 0.032, ε_r
  = 6.63's, `time()` at a step of 2⁵³ + 1,120,445), `cyrius_tfsf_round_once` in `tests/fdtd_source.tcyr` (1, 66) and
  `cyrius_courant_round_once` in `tests/fdtd_util.tcyr` (1, 78). The two-cell test now reads the grid's c_h.
- **Benchmarks**: the step loops don't compute a coefficient. The interleaved A/B against HEAD (two rounds, load 2.9 to
  9.6) left the FDTD rows within the noise: `fdtd/step_100_cells` 260–576 → 264–290 ns (HEAD's 576 an outlier),
  `fdtd/run_100_steps_500_cells` 152–200 → 153–154 µs, `fdtd/tfsf_step_200` 855–1,079 → 888–929 ns (median ratios 1.06
  and 0.99 over 101 rows). A one-off probe, this tree against HEAD at load 2.6: a 20-cell grid's `set_permittivity` 98 →
  254 ns and a 2×2×2 grid's `courant_dt` 139 → 252 ns, the cost of the monomial rounding (about 0.1–0.16 µs), paid at
  setup.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 640,428 ok, 0 failed, 0
  measured, 0 pending, 115,677 accepted (ADR 01: 40,366; ADR 0002: 75,311)** over 756,105 rows; **all 475 digests
  equal**. `fdtd1d.tcyr` (76), `fdtd_source.tcyr` (66) and `fdtd_util.tcyr` (78) pass on all three. The gate is clean
  (no warnings; execution coverage 371 of 371 functions); validate.py 756,105 rows, 0 problems; the oracle's
  `--check` and clippy clean. Next: 2c, the exact thresholds.

## 105. M5 part 2c: the exact thresholds (2026-11-16)

- **ADR 0002 row 17, amended**: fdtd's tests against a literal threshold on the inputs are decided exactly, the
  threshold as its decimal value:
  - `TfsfSource::try_new`'s Courant test c·dt ≤ dx·(1 + 10⁻¹²): Rust's rounded c·dt against RN(dx·RN(1 + 1e-12)),
    where the double 1 + 1e-12 lies 8.9·10⁻¹⁷ above the decimal, so it accepted dt up to an ulp past the limit: dx =
    4.115·10⁻⁶ with dt = 1.3726162517417684·10⁻¹⁴ (over 200,000 seeded walks around the limit, 54,345 decisions
    disagreed). `_bj_fs_courant_over` lets the rounded values decide outside ±10⁻⁶ and otherwise takes the sign of
    the exact 10¹²·c·dt − (10¹² + 1)·dx (rational 15); `TfsfSource::new` refuses as `try_new` does (row 3);
  - `nf2ff_2d_phasor`'s agreement of the three accumulators' ω and dt, |a − b| ≤ 10⁻¹²·max(|a|, |b|), on z0_equal's
    helper for finite values, Rust's test otherwise (an ∞ is refused afterwards as an invalid frequency or dt). The
    helper moves from `src/rf.cyr` (whose file-private `_bj_rf_within12` NF2FF couldn't reach) into
    `src/crmath_edge.cyr` as `_bj_e_within12`, byte for byte, and rf calls it there. As M4 found, the rounded test is
    already exact for normal values, and parts only below about 2^-970, where 1e-12·m rounds to the subnormal grid;
  - `find_eigenmodes`' Nyquist warning freq_max > 1/(2·dt), as freq_max·dt > ½ (`_bj_fu_over_nyquist`: the rounded
    product outside ±2^-50, else the sign of the exact difference, rational 16). It only gates a log line, so no
    golden row shows it: 48,505 of 200,000 seeded boundary cases disagreed with Rust's freq_max > RN(0.5/dt).

  The tests on the computed spectrum (the 10 % floor, the half-power crossing) stay on the rounded values, as M4's
  input impedance kept its test on the computed den (row 17's rule for a test on a computed value: the DFT's own
  rounding moves the spectrum more than the test's).
- **The grid**: 576 rows appended to `fdtd.tsv` (103,726): the Courant limit RN(dx·(1 + 1e-12)/c) walked 24 ulps on
  eleven dx (10⁻⁴, 10⁻³, 3·10⁻⁷, 0.01, three from the seeded search, 10⁻³⁰⁰, 10⁻³¹⁰, 10³⁰⁰ and 0.5) through `try_new`
  and `new`, and NF2FF with hx's dt on M4's cut: dt = (j·10¹² − r)·2^-1042 against dt ± j·2^-1042 for four j and six
  r, 12 of the 48 pairs on the side where the rounded and exact tests part.
- **validate.py**: `tf_check` takes the exact test under row 17, and in port mode the 1D and TF/SF models refuse
  `TfsfSource::new` as `try_new`; a trajectory is named by the first port subset that accounts for its whole change
  (40; 13 with 40; else 17). `fu_nf2ff_port` decides the agreement in Fractions. **28 rows** map to row 17 (8
  `try_new`, 8 `new` and 12 NF2FF calls, each now refused); row 3 grows to 228 with the walks' clamped `new` rows.
- **Tests**: `cyrius_tfsf_courant_exact` in `tests/fdtd_source.tcyr` (2 assertions, 68 in all), `cyrius_nyquist_exact`
  in `tests/fdtd_util.tcyr` (2, 80) and `cyrius_nf2ff_agreement_exact` in `tests/nf2ff.tcyr` (2, 18). rf's 266 pass
  on the moved helper.
- **Benchmarks**: the interleaved A/B against HEAD (two rounds, load 1.7 to 9.2) left every row within the noise (median
  ratios 1.02 and 1.02): the paths through the moved helper, `rf/s_cascade` 3.69–6.86 → 3.75–3.77 µs (HEAD's first run
  under the load spike) and `rf/single_stub_match` 2.18–3.52 → 2.23–2.25 µs, and NF2FF's `fdtd/nf2ff_phasor_80x80_180`
  1.48–1.59 → 1.55 ms. The exact tests run their slow step only near a threshold.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 640,587 ok, 0 failed, 0
  measured, 0 pending, 116,094 accepted (ADR 01: 40,366; ADR 0002: 75,728)** over 756,681 rows; **all 475 digests
  equal**. `fdtd_source.tcyr` (68), `fdtd_util.tcyr` (80), `nf2ff.tcyr` (18) and `rf.tcyr` (266) pass on all three.
  The gate is clean (no warnings; execution coverage 371 of 371 functions); validate.py 756,681 rows, 0 problems; the
  oracle's `--check` and clippy clean. Next: 2d, the range edges.

## 106. M5 part 2d-1: the pulses' range edges (2026-11-16)

M5's 2d, the range edges (row 18), runs in bites: the pulses, the 1D grid (its step updates, the user's answer of
2026-11-15, the Mur ends and `total_energy`), the TF/SF grid, the CW accumulator and the eigenmodes, the utilities
(`graded`, `extents`, the fourth-order derivative) and NF2FF.

- **The finding**: the pulses' chains lose the value where an intermediate leaves the normal range for finite
  arguments:
  - the Ricker wavelet (1 − 2u²)·e^{−u²}: once u² passes DBL_MAX, −∞·0 = NaN where the value is −0; for u² in
    (708, 753), e^{−u²} is subnormal (or 0) and the product hundreds of subnormal ulps off (u² = 730:
    −1.3461193976·10⁻³¹⁴ for −1.3461191496·10⁻³¹⁴);
  - the differentiated Gaussian (−2u/τ)·e^{−u²}: where −2u/τ passes DBL_MAX, ±∞ or NaN for a representable value
    (τ = 10⁻³⁰⁸, u = 3: −∞ for −7.404588245200748·10³⁰⁴); a subnormal u carried its lost bits into −2u/τ;
  - the modulated pulse e^{−u²}·sin(2πf₀(t − t₀)): where the phase passes DBL_MAX, the sine of ∞ (f₀ = 10³⁰⁰, t − t₀ =
    10¹⁰: NaN for −0.40564476173363645), and e^{−u²} subnormal against the sine;
  - all four: t − t₀ past DBL_MAX against a τ that brings u back (t = 1.5·10³⁰⁸, t₀ = −t, τ = 10³⁰⁸: 0 for e^-9).
- **ADR 0002 row 18, amended** (`src/fdtd_source.cyr`): each pulse runs Rust's chain and, where an intermediate it
  reads later (t − t₀, u, u², −2u/τ, e^{−u²}, 2πf₀, the phase) is non-finite or below 2^-1022 for finite arguments,
  reruns it with an unbounded exponent (the Gaussian only where t − t₀ overflows: its other edges already give the
  unbounded value): kimiya's `_bj_k_ux_*` to 53 bits, e^x by `_bj_ux_exp_wide` (1 below 2^-300,
  the far stand-in past 2^47), the sine by `_bj_vx_sin53` (Payne–Hanek past DBL_MAX), and the last product rounded once
  (`_bj_k_ux_mul_f64`). The Gaussian's last step is its exp, so it stays the correctly rounded exp of the 53-bit
  argument. In-range rows keep Rust's bits; a non-finite argument keeps Rust's chain.
- **validate.py**: `fd_pulse_port` transcribes the chains on Fractions (53 bits, unbounded exponent) with e^x and the
  sine from MPFR at its widest exponent range. e^x stays split as (s, e): past |x| = 2300 its power of two would take
  trillions of bits, and a first version that built it was killed for memory, silently; a product whose binade lies
  past ∓1200 is ±0 or ±∞ without it. The 1D and TF/SF models take the port's Gaussian under row 18, and name a
  trajectory row 18 when only it accounts for a change (none does: the scripts' pulses stay in range).
- **The grid**: 2,238 rows appended to `fdtd.tsv` (105,964; 37.5 MB): t − t₀ past DBL_MAX on six τ from 3·10³⁰⁶ to
  DBL_MAX (and four f₀ for the modulated pulse); u² from 700 to 760 on five τ (10⁻³⁰⁰ to 10³⁰⁰); u from 10¹⁵⁵ to
  10³⁰⁰; −2u/τ past DBL_MAX on five τ down to the least subnormal; three subnormal u; the phase past DBL_MAX on three
  f₀; and 400 seeded argument sets from 10^±300 of either sign.
- **The mapping**: **1,848 rows** to row 18 (766 Ricker, 626 differentiated, 443 modulated, 13 Gaussian; part 1b's
  special-value rows among them), all passing.
- **Tests**: `tests/fdtd_source.tcyr`'s `cyrius_pulse_edges` (5 assertions, 73 in all), the five cases above.
- **Benchmarks**: the interleaved A/B against HEAD (two rounds, load 1.7 to 4.0) left the other rows within the noise
  (median ratios 1.00 and 1.06) but the Cyrius-only `fdtd/gaussian_pulse`, 119–129 → 135–148 ns: four range checks on
  every call. The Gaussian's need only one (the bullet above): a subnormal t − t₀ is exact, and a subnormal or infinite
  u gives Rust's e⁰ = 1 or 0 as the unbounded chain does, so it now reruns only where t − t₀ overflows (`fdtd.adr.tsv`
  byte for byte the same). A focused interleaved A/B of the bench's loop, two rounds at load 1.7: 118–121 → 122–123 ns.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 640,977 ok, 0 failed, 0
  measured, 0 pending, 117,942 accepted (ADR 01: 40,366; ADR 0002: 77,576)** over 758,919 rows; **all 475 digests
  equal**. `fdtd_source.tcyr` passes on all three (73). The gate is clean (no warnings; execution coverage 371 of 371
  functions); validate.py 758,919 rows, 0 problems; the oracle's `--check` and clippy clean. Next: 2d-2, the 1D grid
  (its step updates, the Mur ends, `total_energy`).

## 107. M5 part 2d-2a: the 1D grid's updates at their range edges (2026-11-16)

- **The user's answer of 2026-11-15**: the per-cell time-step updates get the range-edge fix (row 18), at the cost
  of a scan each step; M6's 2D and 3D grids and the CPML follow the same rule. This bite is the 1D grid's; the TF/SF
  grid's and the CW accumulator's follow in 2d-3 and 2d-4.
- **The finding**: a cell's update d + k·(a − b) runs Rust's three roundings at any range: opposite fields near
  DBL_MAX side by side overflow a − b, so H₉ = c_h·(E₁₀ − E₉) was −∞ for −3.98·10³⁰⁵ (E = ±1.5·10³⁰⁸); a field
  decaying through the subnormals rounds its product to the subnormal grid before the sum (a 10⁻³⁰⁵ pulse gave
  H₉ = 1.327209363996536·10⁻³⁰⁸ after one step for 1.3272093639965355·10⁻³⁰⁸).
- **ADR 0002 row 18, amended** (`src/fdtd1d.cyr`): a − b is exact when subnormal (Hauser) and overflows only into an
  infinite product, so the product k·(a − b) alone marks the cells whose chain differs from the unbounded one: ±∞,
  NaN or a nonzero subnormal. `_bj_fd1_fix` scans the product row between the packed multiply and add (an inline test
  on the word's exponent field), and for each such cell with finite operands writes d + k·(a − b) with a − b and the
  product to 53 bits with an unbounded exponent and the sum rounded once (`_bj_fd1_cell`, kimiya's `_bj_k_ux_*`),
  setting the product to −0 so the packed add leaves it (r + −0 = r, signed zeros included). The scalar tail and both
  Mur ends take the same test. Every other cell keeps Rust's bits.
- **validate.py**: `Fd1.step_once` under row 18 updates each cell by `fd_cell18`, the same rule on Fractions.
- **The grid**: 79 scenario rows appended to `fdtd.tsv` (106,043): opposite sources of 1.5·10³⁰⁸ and DBL_MAX side by
  side in the interior and at either Mur end, five subnormal or near-DBL_MIN sources (one cell, and both ends), a
  10⁻³⁰⁵ pulse in a 200-cell grid over 1,000 steps, and ε_r or μ_r of 10³⁰⁰. The generator skips frames an earlier
  part emitted (`frames_new`).
- **The mapping**: **60 trajectory rows** to row 18 (59 of the new scenarios and one existing `run`), all passing;
  others of the new rows map to row 13 or 40 where only those change them.
- **Tests**: `tests/fdtd1d.tcyr`'s `cyrius_update_edges` (3 assertions, 79 in all): the opposite sources' H₉, and the
  10⁻³⁰⁵ pulse's H₉ and E₉ after one step.
- **Benchmarks**: the scan is the cost the user's answer accepted, and it is more than my estimate of 1.5–2× when
  asking: the interleaved A/B against HEAD (two rounds, load 0.7 to 1.0) gave `fdtd/step_100_cells` 249 → 862–866 ns,
  `fdtd/step_1000_cells` 2.04–2.09 → 7.86–7.87 µs and `fdtd/run_100_steps_500_cells` 144–146 → 520–523 µs (3.5–3.9×),
  every other row unchanged (median ratios 1.00 and 1.00). The scan is about 2.7 ns a cell and row, against the packed
  verbs' 0.35 ns an element. Variants measured on a 1,000-element row: the indexed loop 2.95 µs, a pointer walk 2.63 µs
  (taken: a 1,000-cell step 7.36–7.53 µs against HEAD's 1.97–2.12 µs, two focused rounds at load 0.1 to 0.4), unrolled
  by four with branches or ORed flags 2.8 µs; a one-test unrolled loop ran 1.7 µs but can't exclude zeros, which fill a
  grid. No vector route exists on 6.7.18: `f64v_*` have no compare or classify, `iv_*` no shift, AND or compare, and a
  float transform needs five passes to tell a zero from a subnormal. Filed in plan §10 as a language gap; M6's 2D and 3D
  grids will pay the same until it closes.
- **Result, every target** (x86_64 ELF, x86_64 PE, aarch64): **18 modules, 0 failing; 640,977 ok, 0 failed, 0
  measured, 0 pending, 118,021 accepted (ADR 01: 40,366; ADR 0002: 77,655)** over 758,998 rows; **all 475 digests
  equal**. `fdtd1d.tcyr` passes on all three (79). The gate is clean (no warnings; execution coverage 371 of 371
  functions); validate.py 758,998 rows, 0 problems; the oracle's `--check` and clippy clean. Next: 2d-2b, the 1D grid's
  Mur coefficient outside its band, `total_energy` and `run`'s source.
