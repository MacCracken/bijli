# 0003 — Retire the Rust tree: the frozen corpus and a Python reference model

**Status**: Accepted (2026-11-16, the user's decision, after M6 part 2a). Supersedes ADR 0001's oracle clause: Rust
1.2.0 stays the specification the port follows, but no Rust is kept, built or run in this repository.
**Date**: 2026-11-16

## Context

The port exists to leave Rust. ADR 0001 kept the Rust 1.2.0 crate in `rust-old/` as the parity oracle until the
Rust-test ledger closed, and the golden corpus was written by a generator that ran it. That generator was itself Rust
(`scripts/oracle/src`, 16,629 lines by M6 part 2a), and it grew with every bite, so the repository carried two Rust
trees (23,752 lines in `rust-old/`) and a `cargo` step in every gate. The user's answer of 2026-11-16: retire all of
it, now.

What the corpus needs from Rust is already in it. Its 766,539 rows are plain data: `cyrius test` reads only the TSVs,
and `scripts/oracle/validate.py` re-derives every row from a Python transcription of Rust's semantics (0
disagreements), with MPFR's correctly rounded values for the corrected oracle. The 1.2.0 source is tagged (`1.2.0`,
`c95fd30`) and every `.rs` file of `rust-old/` is byte for byte the tag's.

## Decision

Remove every Rust file from the tree and every `cargo` step from the process; Rust 1.2.0 remains the specification,
read from the tag.

- **The corpus is frozen data.** `tests/golden/*.tsv` keep their rows; the header line that named the Rust
  generator's command now records that it is retired. The generated sources it wrote once (`src/constants.cyr`,
  `src/fdtd_media_data.cyr`, `src/touchstone_data.cyr`) are frozen too: validate.py checks the first two's every value
  against the corpus (`check_frozen_sources`), and rf.tsv's Touchstone rows exercise the third.
- **validate.py is the reference model.** It keeps checking every frozen row. New rows (M6 part 2b on) are written
  from it: a Python emitter in `scripts/oracle/` builds the inputs, takes the expected outputs from validate.py's
  transcription, and writes them with its `wire()` (format 1, unchanged). Each row then compares two
  implementations, the Python transcription and the Cyrius port, where it compared three (the real Rust too).
- **The port's tools read an inventory.** `scripts/rust_1_2_0.json` is the snapshot of what surface parity, the
  Rust-test ledger, the benchmark record and the constants' check read from `rust-old/`, made by their own parsers
  when the tree was retired: 20 modules and their 733 public items, the 858 tests, the 80 benches and their BASELINE
  medians, and the 12 public constants' literals.
- **Deleted**: `rust-old/`, `scripts/oracle/src/`, `Cargo.toml`, `Cargo.lock` and `build.rs`, the three Rust probe
  references (`scripts/probes/p5_serde_ref.rs`, `p10_kernel_ref.rs`, `p11_rust_sum.rs`; their results are recorded
  in `port-parity.md`) and `scripts/rustsrc.py`. The Rust source is read with `git show 1.2.0:src/fdtd.rs`.

## Consequences

- **Positive** — no Rust anywhere: no toolchain, no `cargo` step in the gate or CI, no second codebase growing beside
  the port. The corpus stays checked end to end, in Python and Cyrius.
- **Negative** — new rows lose the real Rust run. A transcription error shared by validate.py and the port (both
  written from the same source) would go unseen where the Rust run would have caught it; the frozen rows keep that
  check, and new rows lean on the review of the transcription against `git show 1.2.0:…`. The 1.2.0 tests that remain
  (95 pending) are ported from the tag's text.
- **Neutral** — the Python emitter for new rows is written with its first rows (M6 part 2b), not ahead of them.
  ADR 0002's mappings and the corrected oracle are unchanged: they are validate.py's already.

## Alternatives considered

- **Keep the Rust generator until the ledger closes** (ADR 0001 as written): rejected by the user. It kept a growing
  Rust codebase in a project whose purpose is to leave Rust.
- **Freeze the corpus and stop adding rows**: rejected. The rest of M6, M7's surfaces and the range-edge reviews
  need rows, and a row compared against a second implementation is the port's guarantee.
- **Keep `rust-old/` as read-only text**: unneeded. The tag holds it byte for byte, and the tools' needs fit in a
  108 KB inventory.
