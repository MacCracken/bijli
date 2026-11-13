# Probes (Cyrius 6.7.18, 2026-11-11; p16 2026-11-12)

The evidence behind the port plan's §2 (the language review, p1–p6) and its M0 step 9 (p7–p14), recorded in
[`port-parity.md`](../../docs/development/port-parity.md) §4. Each one is a small program that checks one claim
the port relies on. Every float result was compared bit for bit
against Python's IEEE arithmetic, and `p5_serde_ref.rs` produces the serde_json reference text.

Run a probe from the repo root with `cyrius run scripts/probes/<file>.cyr`, or build it
(`cyrius build scripts/probes/<file>.cyr build/<name> && ./build/<name>`) where it is timed or sampled. The `.rs`
files are Rust references (`rustc -O`). Re-run them all at every Cyrius pin bump.

| Probe | Claim | Expected output (6.7.18) |
|---|---|---|
| `p1_complex.cyr` | By-value `C { re: f64; im: f64 }` and `V3` with `T_add`/`T_mul`/`T_div` operators, struct arrays, `*C` stepping, IEEE compares on typed `f64` | `0x4000000000000000 0x3ffc000000000000 0x4023800000000000 0x4013800000000000 0xbfe03f03f03f03f0 0xbfdc0fc0fc0fc0fc 0x403e000000000000`, `negzero == zero: 1`, `nan==nan: 0`, `nan<1: 0`, `0xc019000000000000` |
| `p2_traits.cyr` | `trait` + `impl`, a default method, `struct Grid<B: Backend>` dispatching `g.be.update_h(..)`, `dyn Backend`; `#inline` on four parameters | `0 2 0x3ff0000000000000 79 10`, plus the warning *#inline ignored: fn takes more than 2 parameters* |
| `p3_closures.cyr` | ⚠ a closure capturing `f64` parameters, used directly: wrong (NaN) | `1000 0xfff8000000000000 40 0`: the environment is 40 B, and the calls allocate 0 |
| `p3b_closures.cyr` | the two workarounds: rebinding as `var c: f64 = qx;`, or capturing a struct with `f64` fields | `1000 0x4025ffffffffff95` twice (Python: `0x4025ffffffffff95`) |
| `p4_freelist.cyr` | 1,000 cycles of 6 × 216 KB grids through `fl_calloc`/`fl_free` | `0 0` (no failed frees, no global-heap growth); peak RSS stays flat |
| `p5_derive.cyr` | `#derive(Serialize)` on nested `Complex`, `JonesVector`, `Vec<Complex>`; scanner decode | `{"ex":{"re":0.1,"im":1e-12},"ey":{"re":1.0,"im":-299792458.0}}`, `{"q_sca":8.8541878128e-12,"q_ext":100000000000000000000.0,"a":[{"re":1.5,"im":-0.0}]}`, `0x3da37876f1206634`, `0x8000000000000000`. With `T_from_json(bayan_json_parse(..))` instead: `vec: index out of bounds` |
| `p5b_fmt.cyr` | bayan's number notation at the edges (compare `p5_serde_ref.rs`: serde_json writes `1e+16`, `1.2345678901234568e+17`, `1e-6`, `1e-5`) | `10000000000000000.0 1000000000000000.0 123456789012345680.0 0.000001 1e-7 5e-324 1.7976931348623157e+308 0.00001 0.0001` |
| `p6_misc.cyr` | a decimal `const`; ⚠ captured `f64` as a left operand (`a * 2.0`); `Result` with an f64 payload and `?`; a declared `(f64, f64)` multi-return | `0x3da37876f1206634 0x4008000000000000 0x4008000000000000 0x0 1 0x3ff8000000000000 0 3 0x4023800000000000 0x402d400000000000`: the fourth value, `0x0`, is the trap |
| `p7_grid_memory.cyr` | a 53³ grid with CPML slabs (19.2 MB), 200 alloc/free cycles, `fl_calloc` per array vs one mmap | `failures 0 \| bump growth 0 B \| fl_calloc ≈ 8,700 us/grid \| one mmap ≈ 850 us/grid`; peak RSS ≈ one grid (sample `VmHWM` from outside) |
| `p8_hot_loops.cyr` | complex struct values vs f64 locals: the Mie recurrence and an NF2FF-style accumulation | equal digests per pair (`98b386cd5b5d633`, `8a6749a2e812e8f1`); structs ≈ 2.2× and 2.7× slower |
| `p9_callback.cyr` | a field callback per step: direct, `callptr` on `&fn`, a capturing closure | equal digests (`3fff85b424670a80`); ≈ 50, 50 and 62 ns/step |
| `p10_fdtd_kernel.cyr` + `p10_kernel_ref.rs` | `row_diff` four ways (index, pointer step, `f64v_*`, unrolled) against Rust's | every Cyrius digest `bb2c6a49c09f2a7b`, Rust's the same; `f64v_*` ≈ 1.3 ns/cell, the scalar forms ≈ 5–6, Rust ≈ 0.8 |
| `p11_rust_sum.rs` | Rust's float `Sum` | `sum([]) 8000000000000000` (−0.0), a left fold |
| `p12_text.cyr` + `p12_text_gen.py` + `p12_check.py` | `f64_parse_n` vs Python's (Rust's) parse, bayan's digits vs the shortest | `120018 lines \| parse mismatches vs float(): 0 \| bayan round-trip failures: 0 \| bayan longer than repr: 104, shorter: 0` |
| `p13_cell_step.cyr` | one scalar 3D Yee cell-step, 30³ | ≈ 1,280 us/step (≈ 47 ns per cell-step; Rust 298 us/step) |
| `p14_last_error.cyr` | the per-thread last-error record on two threads | `threads ok: 2 of 2 \| main's own record code 33` |
| `p15_tuple_field.cyr` | `(f64, f64, f64)` tuple fields for Rust's `[f64; 3]` (`PointCharge`): typed element reads, a literal store, a struct copy, `&x + offset` as a `*f64`, by-value passing | `4008000000000000 4022000000000000 4018000000000000 3fec000000000000 40` |
| `p16_json.cyr` | bayan 1.11.1's value tree and printer against serde_json 1.0.154 (M1 part 4): ⚠ `-0` is the integer 0, an integer past i64 is refused, `1e400` is +∞ with no error, a raw control character and an overlong UTF-8 form in a string are taken; a lone surrogate, a repeated key and the 129th level are refused | `2 0 \| ERR integer is out of range for i64 \| 3 7ff0000000000000 \| 5 0 \| 5 0 \| ERR incomplete surrogate pair \| ERR duplicate key \| ERR nesting too deep \|`, then `10000000000000000.0 1e-7 null null -0.0` |
| `p16b_derive_private.cyr` (+ `p16b_private_a.cyr`, `p16b_private_b.cyr`) | ⚠ `#derive(Serialize)` in two `private` files | the build **fails**: `'_cy_json_skipval' is private to its file` (it builds with the first file's `private` removed) |
| `p17_tensor_struct.cyr` | relativity's shapes (M2 part 1d): `EmTensor`'s `[[f64; 4]; 4]` as a 16-element tuple field, a 128-byte struct returned by value and copied, read by `.k` and through a `*f64` at the address an untyped `self` is (⚠ `&self` is the slot holding that pointer, as the guide says of any `*T`), and `FourPotential`'s nested `a: V3` field | `4009000000000000 4009000000000000 1 bff0000000000000 bfc0000000000000 80 3fd0000000000000 3fc0000000000000 20` |
