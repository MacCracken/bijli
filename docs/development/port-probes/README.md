# Language-review probes (Cyrius 6.7.18, 2026-11-11)

These are the evidence behind [`bijli-port-plan.md`](../bijli-port-plan.md) §2. Each one is a small
program that checks one claim the port relies on. Every float result was compared bit for bit
against Python's IEEE arithmetic, and `p5_serde_ref.rs` produces the serde_json reference text.

Run a probe with `cyrius run <file>.cyr`. Run p4 with `cyrius build p4_freelist.cyr p4 && ./p4`.
Re-run them all at every Cyrius pin bump, and move them to `scripts/probes/` at M0.

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
