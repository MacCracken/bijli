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
