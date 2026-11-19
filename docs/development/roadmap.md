# Bijli — Development Roadmap

> Electromagnetism simulation for AGNOS

This file lists planned work only. Everything that shipped is in the [CHANGELOG](../../CHANGELOG.md), the 1.x crate
included. 2.0.0 (2026-11-18) is the Cyrius port of the Rust crate 1.2.0: how it was checked is in
[port-parity.md](port-parity.md), and every behaviour it changed from Rust is in
[ADR 0002](../adr/0002-port-api-divergences.md).

---

## 2.0.x — after the port

These are patch releases: fixes and hardening, no new features. A fix that changes a result is an ADR 0002 row, and
any golden rows it moves are mapped in `*.adr.tsv` (the user's answer of 2026-11-18).

### 2.0.1 — the audit and its repairs

The audit's findings H1–H34 are in [port-parity §142](port-parity.md). They are repaired in these bites, in order:

- [x] **1a** Memory safety with no result changes (H1–H11; port-parity §143).
- [x] **1b** Constructors leave `out` empty on Err; `BjCpml_new` becomes fallible and `enable_cpml` is dropped, as
  ADR row 3 says; the per-step CPML attach check is restored (H12–H15; port-parity §144).
- [x] **2a** Hangs: Mie, the ai client's deadlines and allocations, the mode orders and the GPU constructor's loop
  (H16, H19–H22; port-parity §145).
- [x] **2b-1** The f32 reader's halfway floats: a per-document number index (H17; port-parity §146).
- [x] **2b-2** The parallel combination: a bounded fixed-point tier ahead of the exact one, so that only a sum exactly
  at a tie of its rounding reaches the quadratic tier (H18; port-parity §147).
- [ ] **3** Threads: the math tables become data, plus the error slot's publication and the 3D energy's scratch
  (H23–H25).
- [ ] **4** Results that change, each an ADR 0002 row (H26–H31).
- [ ] **5** Rules and docs: `#must_use`, the stale claims, CLAUDE.md's cleanliness check and CI's missing gates
  (H32–H34).
- [ ] **6** A JSON fuzz harness: the readers rest on bayan's parser and on the golden corpus's documents, not on
  mutation (the [threat model](threat-model.md)'s open item). A harness like `tests/bijli.fcyr`, the Touchstone
  reader's, closes it for `_bj_jr_doc` and a few readers.
- [ ] **7** The release: the post-review bench against §142's baseline, the cross-targets, the docs.

### Hardening

- [ ] **Retire the port's handoff**: [handoff.md](handoff.md) has served the port. Delete it once no 2.0.x work needs
  it; git keeps it.

### Performance

None of these changes a result: the bits stay the same. The numbers go in
[benchmarks-rust-v-cyrius.md](../benchmarks-rust-v-cyrius.md).

- [ ] **crmath's sin and cos past 2^20**: a medium-range argument reduction (Cody–Waite with a longer π/2, or
  crmath_ux's Payne–Hanek) ahead of the fast path's polynomials. Today every |x| > 2^20 takes the multi-precision slow
  path, about 21 µs a call, and an optical plane wave reaches it routinely (`wave/plane_wave_e`, port-parity §26).
- [ ] **The FDTD updates' per-cell scan**: the range-edge rule's scan for products that leave the normal range costs
  about 2.7 ns a cell beside the packed update's 1 ns. Against their unchecked forms:
  - 1D: the step runs 3.5–3.8×, the TF/SF step about 2×, `CwAccumulator::accumulate` about 3× (port-parity §107, §109,
    §110).
  - 2D: the step runs 2.4–2.7× (about 2.3 ns a cell), the CPML step 1.4×, `total_energy` 2.3×.
  - 3D: the step runs 1.8× (3.4 ns a cell update), the CPML step 1.5×, `total_energy` 2.8× (port-parity §127–§130).

  A lane-wise exponent-field test in the stdlib's vector verbs (the filing below) would bring each scan to about one
  pass.
- [ ] **`propagate_beam`'s general fast path**: with C ≠ 0 (a lens in the system), q_out rounded once takes about
  650 ns where Rust took 1.8 ns. That is the call overhead of its six TwoProducts, five double-double products and two
  quotient tests on this toolchain. With C = 0 (free space, interfaces) it takes 94 ns, and `transform_beam` takes
  1.2 µs where Rust took 0.58. The first steps are inlining TwoProduct and deciding the threshold on the fast path's
  own denominator (port-parity §66).

### Upstream

- [ ] **File the gaps the plan's §10 marks as not yet filed**, each with its probe:
  - a captured `f64` is an untyped word inside a closure (p6);
  - `#derive(Serialize)`'s pairs decoder returns an empty `Vec<struct>` (p5);
  - bayan's exponent notation differs from serde_json's (p5b, a note);
  - no vector verb classifies a lane (the FDTD scan above).
- [ ] **File the PE syscall warnings** that sandhi's stack prints under `CYRIUS_TARGET_WIN`: `lib/dynlib.cyr` lines
  674 and 950, and `lib/fdlopen.cyr` lines 472, 530, 626 and 954. Once a fixed stdlib is pinned, drop the PE part of
  CLAUDE.md's exception for `tests/ai.tcyr`.
- [ ] **File sandhi's client gaps** (port-parity §142):
  - DNS and the TLS handshake have no deadline;
  - the resolver skips `/etc/hosts`, so the ai client's default `localhost` endpoint resolves through DNS;
  - the response's header parse allocates from the never-freed bump heap.

---

## 2.1.x – 2.x.x — open work

These are minor releases: new API is appended (enum values too), and each feature gets reference tests, plus golden
rows from validate.py's model where it has one. The order follows plan §11.

Every addition must keep the default bundle under the compiler's 128 KiB static-data advisory. bijli alone is just
under it (port-parity §135), so a large table is allocated on first use.

### The 1.2.0 audit's backlog (plan §11.1)

- [ ] **CPML and materials for `GpuFdtd2d`**: the compute-backend solver has PEC walls and vacuum only.
- [ ] **TF/SF on the `Fdtd2d` contour**: `BjTfsfSource` is a standalone 1D auxiliary grid. The work is to inject its
  plane wave on the 2D grid's TF/SF contour.
- [ ] **Frequency-aware lossy Fresnel**: a complex ε = ε − jσ/ω (in the e^{+jωt} convention), so that conductors are
  supported. Today the Fresnel functions return an `Err` for them.
- [ ] **Sub-femtometre point-charge separations**: the field, potential and Coulomb functions treat r² < 10^-30 m² as
  singular (an exact test since ADR 0002 row 17). Accepting closer points changes results that 2.0 refuses, so this
  item needs an ADR.

### A GPU backend (plan §11.2)

- [ ] **`BjMabdaBackend`**: `impl BjComputeBackend` on mabda (wgpu-native, native AMD with an f64 SPIR-V compiler,
  native NVIDIA compute), in a `[lib.gpu]` profile. It sits behind the API of the ported `BjWgpuBackend` scaffold.
  The backend is new code, so its parity is against `BjCpuBackend`: bit for bit where the device's f64 allows.

### Shared numerics (plan §10, G1–G4)

- [ ] **Depend on the ecosystem's math once it exists.** The candidates:
  - correctly rounded f64 primitives as a library (G1; bijli, impetus and kimiya each carry their own copy);
  - spherical, Riccati and cylindrical Bessel functions with complex arguments (G2);
  - associated Legendre functions and spherical harmonics (G3);
  - Hermite and generalized Laguerre polynomials (G4).

  bijli moves to one only when it gives bijli's results bit for bit. Until then, `src/crmath*.cyr` and the special
  functions stay bijli's own.

### Research watch (not scheduled)

> These are tracked until each reaches practical maturity.

- [ ] Physics-informed neural networks (PINNs) as surrogate EM solvers (17× speedup published 2025)
- [ ] Adaptive mesh refinement (AMR) for FDTD — h/p/r refinement, subgridding
- [ ] Spectral element methods for smooth geometries
- [ ] Discrete dipole approximation (DDA) for non-spherical scattering
- [ ] SR-FDTDNet — super-resolution reconstruction of fine-grid FDTD from coarse runs

---

## Boundaries

| Domain | Owned by | NOT bijli |
|--------|----------|-----------|
| General PDE solvers | hisab | bijli's solvers are its EM ones (FDTD, eigenmodes); no hisab dependency (ADR 0001) |
| Particle dynamics | impetus | bijli computes forces, impetus integrates |
| Optics (ray, wave) | prakash | bijli provides EM foundations; prakash builds ray tracing, lens design, etc. |
| Quantum EM (QED) | future | bijli is classical EM |
| GPU compute | mabda | bijli provides the `BjComputeBackend` trait; backends are pluggable (mabda's is 2.x work) |
| ML surrogates | daimon/hoosh | bijli defines solver interfaces; AI agents provide trained surrogates |
