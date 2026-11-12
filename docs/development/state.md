# bijli — Current State

> Refreshed every release and at every port milestone. `CLAUDE.md` is
> preferences/process/procedures (durable); this file is **state** (volatile).

## Version

`VERSION` reads **1.2.0**: the oracle's version, kept until the Cyrius library ships as **2.0.0**
(M8). `cyrius.cyml` reads it through `version = "${file:VERSION}"`, and the release workflow
checks both against the tag. The last Rust release is **1.2.0**, tag `1.2.0` at `c95fd30`
(2026-11-11), frozen as the oracle in `rust-old/`. It includes the 3D CPML. crates.io stops at
1.1.0.

## Port progress

Plan of record: [`bijli-port-plan.md`](bijli-port-plan.md). Decisions:
[ADR 0001](../adr/0001-cyrius-port-with-rust-oracle.md) (the port) and
[ADR 0002](../adr/0002-port-api-divergences.md) (the 2.0 API changes).

| Milestone | Status |
|---|---|
| P Prep | ✅ **done 2026-11-11.** D1–D5 answered (1.2.0 tag; one convention, e^{+jωt} with IEEE handedness; derive writers with hand-written readers). Oracle tagged and green. Rust bench baseline at `2026-11-12T03:15:21Z` (under load). glibc call map recorded (plan §5 P3). `target/` cleaned. Language review with 8 probes (`port-probes/`) |
| M0 Initialize | ◐ **Part 1 done 2026-11-11, steps 1–3.** `cyrius port` scaffold. `rust-old/` matches the tag file for file. Repairs: CLAUDE.md from the P6 draft; roadmap restored; this file; CI extended (`test --bounds`, `deps --locked`, `fuzz`, `audit`), with the Rust workflows in `rust-old/.github/workflows/`; Rust-only scripts and the Rust architecture overview moved to `rust-old/`; ADR 0001 accepted and ADR 0002 written; `cyrius.cyml` with `${file:VERSION}`, `bench` in `[test] stdlib`, M0's stdlib leaves and no hisab; the scaffold's six stub fns documented. Gate: `build`, `test` (0 warnings), `test --bounds`, `bench`, `fuzz`, `lint`, `fmt --check`, `vet`, `deps --locked` and `audit` all pass. **Part 2a done 2026-11-11, steps 5–6.** `scripts/oracle/`: the generator crate (TSV format 1, glibc provenance, seeded grids, a `--check` mode), its lock seeded from `rust-old/`'s (193 crates, the same versions), fmt- and clippy-clean. It writes `tests/golden/constants.tsv` (12 rows), `math.tsv` (34,734 rows over the 25 std f64 methods bijli calls) and the generated `src/constants.cyr` (all 12 constants read back bit-exact in Cyrius). `--check` and `--locked` reproduce it byte for byte. `validate.py`, an independent Python transcription with glibc through `ctypes` and MPFR (`mpfr_ref.py`): 34,746 rows, 0 problems. **It caught a P3 error**: Rust std's `cbrt` is compiler-builtins' correctly rounded port of CORE-MATH, not glibc's (glibc disagrees on 504 of 1,084 rows, MPFR on none). **Next: Part 2b** (step 4 `mathx.cyr`, step 7 the parity harness and tooling), then Part 3 (steps 8–9) |
| M1 Foundation | ☐ |
| M2 Media, waves, circuits, relativity | ☐ |
| M3 Optics | ☐ |
| M4 RF | ☐ |
| M5 FDTD core | ☐ |
| M6 FDTD 2D/3D | ☐ |
| M7 Integration surfaces | ☐ |
| M8 Release 2.0.0 | ☐ |

## Toolchain

- **Cyrius pin**: `6.7.18` (`cyrius.cyml` `[package].cyrius`)
- **ganita**: 1.2.15, from the 6.7.18 stdlib bundle
- **hisab**: none in production (ADR 0001)
- **Rust oracle**: bijli 1.2.0 at `rust-old/`, on rustc 1.99.1. Its transcendentals come from
  the host **glibc 2.44** (Rust std), except `cbrt` (compiler-builtins' correctly rounded port of
  CORE-MATH) and `powi` (compiler-rt's `__powidf2`); its complex arithmetic from Rust hisab
  1.4.0's `Complex`
- **MPFR**: 4.2.2 (system `libmpfr.so.6`), for `validate.py` and the corrected oracle
  (development only)

## Source

- **Rust reference**: `rust-old/`, 22,696 lines in all, of which 8,589 are code in `src/` (no
  tests, comments or blank lines). Frozen, never edited.
- **Cyrius port**: the `src/main.cyr` stub and the generated `src/constants.cyr` (12 public
  constants). The golden corpus is in `tests/golden/` (34,746 rows).

## Tests

- **Cyrius**: the scaffold's smoke group (2 assertions), one bench (`noop`) and one fuzz stub.
- **Rust oracle**: 830 unit, 26 integration and 2 doc tests, green at the tag. They are green from
  `rust-old/` too (`cargo test --manifest-path rust-old/Cargo.toml --all-features`, 2026-11-11),
  with no warnings.

## Benchmarks

The Rust column is the `2026-11-12T03:15:21Z` run at `c95fd30`: 80 benches, in
`rust-old/bench-history.csv` and `rust-old/benchmarks.md`. It ran under load (kimiya's gates):
median CI width 0.8 %, p90 4.2 %. A Rust-vs-Cyrius ratio within the Rust row's CI width is read
as parity.

## Dependencies

Declared in `cyrius.cyml`:
- stdlib: syscalls, string, alloc, str, fmt, io, vec, args, math, ganita, tagged, result, fnptr,
  freelist, mmap, overflow;
- test stdlib (`[test] stdlib`, test / bench / fuzz builds only): assert, bench.

## Consumers

- Cyrius: prakash 2.15.2 (its own bijli bridges) and tanmatra 2.0.x (its bijli bridge waits on
  this port).
- Rust: soorat 1.0.0 (the soorat integration structs) and ushma 1.4.0, which is planning its own
  port.
- No local repo: kiran, joshua, dravya.
