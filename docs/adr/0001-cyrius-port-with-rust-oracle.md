# 0001 — Port to Cyrius, with Rust 1.2.0 frozen as the parity oracle

**Status**: Accepted (2026-11-11, at M0 step 2 of the [port plan](../development/bijli-port-plan.md)).
Plan §0.3's decisions: D1 the `1.2.0` tag at `c95fd30`, which includes the 3D CPML; D3 e^{+jωt}
with IEEE handedness; D5 derive writers with hand-written readers.
**Date**: 2026-11-11

## Context

AGNOS libraries are moving from Rust to Cyrius. Three siblings set the method:
- tanmatra ported first (2.0.0, 2026-11-08);
- impetus followed (2.0.0, 2026-11-10);
- kimiya is partway through.

Each freezes an audited Rust release as the oracle and compares the port against it value by
value. Two of bijli's partners are already in Cyrius: prakash, and tanmatra, which has deferred
its bijli bridge until bijli ports.

bijli 1.2.0 audited and repaired every module, and it includes the 3D CPML, measured at −80 dB.
It is tagged `1.2.0` at `c95fd30` (2026-11-11). The tagged tree is green on rustc 1.99.1: fmt,
clippy `-D warnings`, rustdoc, audit, deny and `make features` all pass, as do the 830 unit, 26
integration and 2 doc tests. So it can serve as the oracle.

bijli differs from the siblings in ways that shape the decision:
- **Its only hisab use is `hisab::num::Complex`**, and Cyrius hisab's complex type differs from
  it. hisab's `cx_*` uses Smith division, a hypot-rescued `abs` and a C99 `sqrt`, and it
  allocates 16 B on the never-freed heap for every result.
- **Its FDTD solvers are linear and use only `+ − × ÷`**, so whole trajectories can be compared
  bit for bit. impetus's chaotic piles could not be.
- **Its math reaches glibc 2.44** for sin, cos, exp, atan2, log10, pow, hypot, atan, asin, acos,
  tan and cbrt. glibc isn't correctly rounded, and the user's standing ruling (kimiya
  ADR 0006) is that ports are more accurate than glibc, not bit-compatible with it.

## Decision

bijli 2.0.0 is a Cyrius library. Rust **1.2.0** (tag `1.2.0`, `c95fd30`) is frozen in `rust-old/` as the parity oracle,
and nothing in it is edited.

- **Initialization.** Start with `cyrius port` on Cyrius 6.7.18, after the prep in plan §4.
  Repair the files the scaffold overwrites, as the siblings did.
- **Scope, milestones and parity policy** are those of
  [`bijli-port-plan.md`](../development/bijli-port-plan.md):
  - **the whole crate is ported**: soorat integration, bridges, the `ai` client (a `[lib.ai]`
    profile) and the `WgpuBackend` scaffold included;
  - only the three functions bijli's own roadmap deprecates are dropped;
  - nothing is deferred for want of a language feature. Plan §2 reviews Cyrius 6.7.18 against
    every Rust construct bijli uses, with probe programs in `docs/development/port-probes/`.
- **Parity classes:**
  - A: exact;
  - B: arithmetic;
  - C: transcendental, at 0 ulp against the corrected oracle;
  - D: iterative, expected empty;
  - T: an FDTD trajectory, compared by digest at chosen frames;
  - R: reference.
- **Accuracy.** Every transcendental is correctly rounded and target-independent. The primitives
  are reused from impetus's `crmath.cyr` and kimiya's `mathx.cyr`; asin and tan are new. Inexact
  formulas and range edges are fixed within plan §0.1's scope, one ADR 0002 row each.
- **One time convention.** Rust mixes e^{−iωt} (polarization, beam, field, the FDTD media) and
  e^{+jωt} (scattering, rf, the FDTD phasors). The port uses **e^{+jωt} with IEEE handedness** in
  every module (plan §4.4), the convention of the FDTD, RF, antenna and beam literature bijli
  implements. Amplitudes are compared as `conj(oracle)` and phases as `−oracle`; observables are
  compared directly. The port removes Rust's faults; it doesn't carry them.
- **No production hisab dependency.** bijli owns a by-value `BjComplex` that transcribes Rust
  hisab 1.4.0's `Complex` operation for operation.
- **Divergences.** A divergence from the oracle beyond its class is a defect, unless ADR 0002
  accepts it by name. ADR 0002 also holds the API changes of plan §4.3.

## Consequences

- **Positive**
  - Correctness is measured, not asserted: every in-scope public function is diffed against the
    audited Rust crate, bit pattern by bit pattern. The FDTD is diffed over whole runs.
  - The 1.2.0 audit's reference-value tests travel with the code.
  - bijli has no dependency beyond the stdlib, so a consumer's build stays small. impetus's build
    shrank from 3.68 MB to 1.70 MB when hisab left its production path.
- **Negative**
  - Two toolchains are needed until the port completes: cargo for the oracle and its generator,
    cyrius for the library.
  - `rust-old/` stays in the tree until the Rust-test ledger is closed and the corpus archived.
  - A third copy of the correctly rounded primitives exists until the shared module (plan §10 G1)
    lands upstream.
  - Large grids need explicit `*_free` calls: the Cyrius API isn't RAII. freelist returns their
    memory (blocks over 4 KB are munmapped), and `defer` keeps tests tidy.
  - Closure callbacks must rebind captured floats as typed locals (plan §2), a trap no sibling
    met.
- **Neutral**
  - The generator drives `rust-old/` by path from `scripts/oracle/`, so the oracle itself is never
    modified. A wrong oracle value is fixed in Rust 1.2.x and the corpus is regenerated.
  - Cargo features become `[lib]` plus a `[lib.ai]` profile.

## Alternatives considered

- **Port and fix together, with no frozen oracle.** Rejected, as in every sibling: a translation
  error and a deliberate fix look the same.
- **Defer the soorat integration and drop `WgpuBackend`, as the siblings deferred theirs.**
  Rejected after the language review: the soorat module has no soorat dependency (goonj ported
  its own), and the scaffold is one struct and one function.
- **Leave the 3D CPML out of the oracle and carry it as Cyrius-only work.** Rejected: it leaves the
  newest and hardest code (12 Ψ slabs, edges and corners) without an oracle.
- **Use Cyrius hisab's `cx_*` for complex values.** Rejected: its results differ from the
  oracle's formulas, and it allocates per operation inside Mie and NF2FF loops.
- **Keep bijli in Rust.** Rejected: its Cyrius consumers would need a foreign-function bridge
  for every call, against the AGNOS direction.
