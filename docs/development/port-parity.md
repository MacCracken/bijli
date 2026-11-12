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
