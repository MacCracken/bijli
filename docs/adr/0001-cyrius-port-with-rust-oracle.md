# 0001 — Port to Cyrius, with Rust 1.3.0 frozen as the parity oracle

**Status**: Proposed. It is accepted at M0 step 2 of the [port plan](../development/bijli-port-plan.md),
once D1–D5 (plan §0.3) are answered.
**Date**: 2026-11-11

## Context

AGNOS libraries are moving from Rust to Cyrius. Three siblings set the method:
- tanmatra ported first (2.0.0, 2026-11-08);
- impetus followed (2.0.0, 2026-11-10);
- kimiya is partway through.

Each freezes an audited Rust release as the oracle and compares the port against it value by
value. Two of bijli's partners are already in Cyrius: prakash, and tanmatra, which has deferred
its bijli bridge until bijli ports.

bijli 1.2.0 audited and repaired every module. The `[Unreleased]` 3D CPML is measured at −80 dB.
The Rust tree is green on rustc 1.99.1: fmt, clippy `-D warnings`, rustdoc, audit and deny pass,
and the 830 unit, 26 integration and 2 doc tests pass. So once it is tagged, it can serve as the
oracle.

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

bijli 2.0.0 is a Cyrius library. Rust **1.3.0** is frozen in `rust-old/` as the parity oracle,
and nothing in it is edited.

- **Initialization.** Start with `cyrius port` on Cyrius 6.7.18, after the prep in plan §4.
  Repair the files the scaffold overwrites, as the siblings did.
- **Scope, milestones and parity policy** are those of
  [`bijli-port-plan.md`](../development/bijli-port-plan.md):
  - physics first;
  - soorat integration deferred until soorat ports;
  - `WgpuBackend` not carried.
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
- **No production hisab dependency.** bijli owns a by-value `BjComplex` that transcribes Rust
  hisab 1.4.0's `Complex` operation for operation.
- **Divergences.** A divergence from the oracle beyond its class is a defect, unless ADR 0002
  accepts it by name. ADR 0002 also holds the API changes of plan §3.3.

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
  - `rust-old/` stays in the tree until soorat ports or is dropped.
  - A third copy of the correctly rounded primitives exists until the shared module (plan §9 G1)
    lands upstream.
  - Large grids need explicit `*_free` calls: the Cyrius API isn't RAII.
- **Neutral**
  - The generator drives `rust-old/` by path from `scripts/oracle/`, so the oracle itself is never
    modified. A wrong oracle value is fixed in Rust 1.3.x and the corpus is regenerated.
  - Cargo features become `[lib]` plus a `[lib.ai]` profile.

## Alternatives considered

- **Port and fix together, with no frozen oracle.** Rejected, as in every sibling: a translation
  error and a deliberate fix look the same.
- **Port from 1.2.0 and carry the 3D CPML as Cyrius-only work.** Rejected: it leaves the newest
  and hardest code (12 Ψ slabs, edges and corners) without an oracle.
- **Use Cyrius hisab's `cx_*` for complex values.** Rejected: its results differ from the
  oracle's formulas, and it allocates per operation inside Mie and NF2FF loops.
- **Keep bijli in Rust.** Rejected: its Cyrius consumers would need a foreign-function bridge
  for every call, against the AGNOS direction.
