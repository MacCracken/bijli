# Bijli — Claude Code Instructions

## Project Identity

**Bijli** (Hindi: बिजली — electricity, lightning) — Electromagnetism simulation — fields, Maxwell's equations, charge dynamics, EM waves

- **Type**: Flat crate (library)
- **License**: GPL-3.0
- **Version**: 0.1.0

## Consumers

kiran (game engine), joshua (simulation), prakash (optics)

**Foundation**: hisab (PDE solvers), impetus (particle dynamics)

**Modules**: field, maxwell, charge, wave, error, ai, logging

## Development Process

### P(-1): Scaffold Hardening (before any new features)

1. Test + benchmark sweep of existing code
2. Cleanliness check: `cargo fmt --check`, `cargo clippy --all-features --all-targets -- -D warnings`, `cargo audit`, `cargo deny check`
3. Get baseline benchmarks (`./scripts/bench-history.sh`)
4. Initial refactor + audit (performance, memory, security, edge cases)
5. Cleanliness check — must be clean after audit
6. Additional tests/benchmarks from observations
7. Post-audit benchmarks — prove the wins
8. Repeat audit if heavy

### Development Loop (continuous)

1. Work phase — new features, roadmap items, bug fixes
2. Cleanliness check
3. Test + benchmark additions for new code
4. Run benchmarks (`./scripts/bench-history.sh`)
5. Audit phase
6. Cleanliness check — must be clean after audit
7. Deeper tests/benchmarks from audit observations
8. Run benchmarks again — prove the wins
9. If audit heavy → return to step 5
10. Documentation — update CHANGELOG, roadmap, docs
11. Return to step 1

### Key Principles

- **Never skip benchmarks.** Numbers don't lie.
- **Tests + benchmarks are the way.** Minimum 80%+ coverage target.
- **Own the stack.** Depend on AGNOS crates, not external libs directly.
- **`#[non_exhaustive]`** on all public enums.
- **`#[must_use]`** on all pure functions.
- **`#[inline]`** on hot-path functions.
- **Feature-gate optional deps.**
- **tracing on all operations.**

## DO NOT

- **Do not commit or push** — the user handles all git operations
- **NEVER use `gh` CLI** — use `curl` to GitHub API only
- Do not add unnecessary dependencies
- Do not `unwrap()` or `panic!()` in library code
- Do not skip benchmarks before claiming performance improvements
