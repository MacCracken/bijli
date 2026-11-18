# Architecture

bijli 2.0 is one Cyrius library, the port of the Rust crate 1.2.0 (ADRs
[0001](../adr/0001-cyrius-port-with-rust-oracle.md) and [0003](../adr/0003-retire-the-rust-tree.md)). Its sources are
the modules below, concatenated in `cyrius.cyml`'s
`[lib]` order into one bundle, `dist/bijli.cyr`; `[lib.ai]` adds the ai client for `dist/bijli-ai.cyr`. No source
file includes another: each resolves against what precedes it in that order.

## Module map

| Layer | Files | What it holds |
|---|---|---|
| Correctly rounded math | `crmath.cyr` (+ `crmath_data.cyr`), `crmath_k.cyr` (+ `crmath_k_data.cyr`), `crmath_kx.cyr`, `crmath_mono.cyr`, `crmath_ux.cyr` (+ `crmath_ux_data.cyr`), `crmath_gamma.cyr` | exp, sin, cos, tan, sinh, cosh, hypot, pow, atan2, acos, asin, ln, log10, expm1, correctly rounded and target-independent (impetus's and kimiya's, with bijli's own); monomials rounded once; binary64 with an unbounded exponent for the range edges; the Lorentz factor |
| Exact tests | `crmath_cmp.cyr`, `crmath_vx.cyr`, `crmath_edge.cyr` | the exact Singularity threshold, 3-vector chains at a range edge, the per-cell rule of the FDTD updates and the rationals rounded once |
| Std semantics | `mathx.cyr`, `constants.cyr` | Rust's f64 methods (min and max with −0 < +0, clamp, the predicates, `powi`), the physical constants |
| Text, errors, logging | `fmtf.cyr`, `text.cyr`, `error.cyr`, `log.cyr` | Rust's Display, Debug and serde_json number forms; the bounded text consumers format into; the `BJ_ERR_*` codes and the per-thread last-error record; logging on sakshi |
| Fields and charges | `field.cyr`, `maxwell.cyr`, `charge.cyr`, `material.cyr` | `BjFieldVector`, point charges, dipoles, multipoles, Green's functions; Maxwell's equations; `BjPointCharge` and the Lorentz force; `BjMaterial` and bound charges |
| Waves and circuits | `wave.cyr` (+ `wave_data.cyr`), `circuit.cyr` | plane waves, Snell and Fresnel, waveguides, antennas; DC and transient circuits |
| Relativity and optics | `relativity.cyr`, `complex.cyr`, `polarization.cyr`, `scattering.cyr`, `beam.cyr` | the EM tensor, transforms, four-vectors, retarded fields; `BjComplex` and the time convention (e^{+jωt}, stated once here); Jones, Stokes and Mueller; Mie and Rayleigh; Gaussian beams, ABCD, resonators, HG and LG modes |
| RF | `rf.cyr`, `touchstone.cyr` (+ `touchstone_data.cyr`) | lines, matching, the Smith chart, arrays, `BjSMatrix`; Touchstone v1.1 read and write |
| FDTD | `fdtd1d.cyr`, `fdtd_source.cyr`, `fdtd_media.cyr` (+ `fdtd_media_data.cyr`), `fdtd_util.cyr`, `cpml.cyr`, `fdtd2d.cyr`, `fdtd3d.cyr`, `nf2ff.cyr`, `compute.cyr` | the 1D, 2D and 3D Yee grids; pulses and TF/SF; lossy, Drude, Lorentz and Debye media; the non-uniform grid, the CW accumulator, eigenmodes; the CPML; near-to-far-field; the compute-backend trait and `GpuFdtd2d` |
| Integration | `bridge.cyr` (+ `bridge_data.cyr`), `soorat.cyr`, `json.cyr`, `ai.cyr` (`[lib.ai]` only) | the cross-crate bridges; soorat's visualization structs; every serializable type's JSON on bayan; the daimon/hoosh client on sandhi |

Files named `*_data.cyr` are generated (by `scripts/oracle/`) or frozen from Rust 1.2.0, and checked by
`scripts/oracle/validate.py`; they are never edited by hand.

## Data flow

- **Values in, values out.** A function takes f64 values, pointers to arrays and structs by value, and returns a
  value-form `Result`: Ok with the value (or 0 when the value went to an `out` pointer), or Err with a `BJ_ERR_*`
  code. The error's fields and Rust's Display text go into the calling thread's last-error record, which
  `bj_last_error_*` read; nothing is logged for an error and nothing allocates.
- **Owned memory.** A grid, an S-matrix, Touchstone data, a field slice or a JSON reading owns its blocks (freelist
  or one mmap), allocated at construction and given back by its `_free`; nothing allocates inside a time step. A
  failed allocation is `BJ_ERR_ALLOCATION`, never an abort.
- **Arithmetic.** Each function follows Rust 1.2.0's order of operations, so its bits are Rust's, except where ADR
  0002 records a fix: products and quotients rounded once (row 13), the range edges on an unbounded exponent (rows
  14–18), the cancellation fixes. Every transcendental is correctly rounded, so the bits are the same on every target.
- **Text.** Rust's `&str` inputs are bytes and a length, checked for UTF-8; `display` and `debug` write into a text
  the caller owns (`bj_text_*`); JSON writers append to a `str_builder`, readers parse with bayan's tree into a
  decode region freed before they return.

## How it is checked

- `tests/golden/` holds the corpus: 770,127 rows over 21 modules, frozen from Rust 1.2.0's generator or emitted by
  `scripts/oracle/emit.py`, with the corrected oracle (`*.cr.tsv`, MPFR where glibc misrounded) and the ADR 0002
  mappings (`*.adr.tsv`). `scripts/oracle/validate.py` re-derives every row in Python.
- `tests/parity.tcyr` replays every row against the port and prints one digest per function; the digests are equal on
  x86_64 ELF, x86_64 PE and aarch64.
- `tests/*.tcyr` port all of Rust 1.2.0's tests by name (`docs/development/rust-test-ledger.md`); `tests/bijli.bcyr`
  twins its 80 benchmarks (`docs/benchmarks-rust-v-cyrius.md`).

## Dependencies and consumers

- **Stdlib only**: syscalls, string, alloc, str, fmt, io, vec, args, math, ganita, tagged, result, fnptr, freelist,
  mmap, overflow, atomic, thread_local, sakshi (logging) and bayan (JSON); sandhi (HTTP/TLS) for `[lib.ai]` alone,
  never in the default bundle. No hisab: `BjComplex` and the math are bijli's own.
- **Consumers**: prakash (optics) and tanmatra (atomic physics) use the bridges; kiran (the game engine) and joshua
  (simulation) the whole library; soorat renders the soorat structs. `docs/development/state.md` keeps the list.

## Notes

Non-obvious invariants that the code alone doesn't show are numbered notes beside this file (`NNN-title.md`). There
are none yet; decisions are ADRs ([`../adr/`](../adr/)), not notes.
