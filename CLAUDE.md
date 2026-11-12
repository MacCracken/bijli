# Bijli — Claude Code Instructions

> **Core rule**: this file is **preferences, process, and procedures** —
> durable rules that change rarely. Volatile state (current version,
> module line counts, port progress, test counts, consumers) lives in
> [`docs/development/state.md`](docs/development/state.md).
> Do not inline state here.

## Project Identity

**Bijli** (Hindi: बिजली — electricity, lightning) — Electromagnetism simulation — fields,
Maxwell's equations, charge dynamics, EM waves

- **Type**: flat library (Cyrius), ported from the Rust crate 1.2.0; the port ships as 2.0.0
- **License**: GPL-3.0-only
- **Language**: Cyrius. The toolchain is pinned in `cyrius.cyml` `[package].cyrius`
- **Version**: `VERSION` at the project root is the source of truth. Don't inline the number here
- **Rust oracle**: bijli 1.2.0 (tag `1.2.0`), frozen at `rust-old/`. It stays the parity reference
  until the Rust-test ledger is closed and the corpus archived
- **Plan of record**: [`docs/development/bijli-port-plan.md`](docs/development/bijli-port-plan.md),
  with [ADR 0001](docs/adr/0001-cyrius-port-with-rust-oracle.md) (the port) and ADR 0002 (the 2.0
  API changes)
- **Standards**: [First-Party Standards](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-standards.md) · [First-Party Documentation](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-documentation.md)

## Goal

bijli owns classical electromagnetism for AGNOS: static and dynamic fields, Maxwell's equations,
charges, waves and materials, circuits, relativistic electrodynamics, polarization, scattering,
Gaussian beams, RF and microwave networks, and FDTD solvers in 1D, 2D and 3D. Results are
deterministic and bit-identical across targets. bijli uses one harmonic time convention:
**e^{+jωt}, with IEEE handedness**.

## Consumers

prakash (optics, Cyrius) and tanmatra (Cyrius; its bijli bridge waits on this port) use bijli's
bridges. kiran (game engine) and joshua (simulation) use the full library. soorat (Rust) renders
the soorat integration structs. The current list is in `state.md`.

**Dependencies**: the Cyrius stdlib only: ganita, sakshi, freelist, bayan, plus sandhi in the
`[lib.ai]` profile. There is no production hisab dependency (ADR 0001).

## Quick Start

```sh
cyrius deps                                    # resolve dependencies (vendors lib/)
cyrius build src/main.cyr build/bijli          # compile
cyrius test                                    # tests/*.tcyr: parity, trajectory, reference
cyrius bench                                   # tests/*.bcyr
cyrius distlib --all                           # dist/bijli.cyr and dist/bijli-ai.cyr
cargo test --manifest-path rust-old/Cargo.toml --all-features                       # the Rust oracle, read-only
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # the golden corpus is up to date
python3 -I -B scripts/oracle/validate.py .                                          # independent check of the corpus
```

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read the roadmap, the port plan, the CHANGELOG and the open issues: know what was intended
   before auditing what was built
1. Test + benchmark sweep of the existing code
2. Cleanliness check (below)
3. Get baseline benchmarks (`cyrius bench`, recorded in `docs/benchmarks-rust-v-cyrius.md`)
4. Internal deep review — gaps, optimizations, security, logging/errors, docs
5. External research — domain completeness, missing capabilities, best practices, world-class accuracy
6. Cleanliness check — must be clean after the review
7. Additional tests/benchmarks from the findings
8. Post-review benchmarks — prove the wins
9. Repeat if heavy

### Cleanliness check

```sh
cyrius fmt --check <files>
cyrius lint <files>
cyrius vet src/main.cyr
cyrius build src/main.cyr build/bijli
cyrius test                          # must be warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius coverage --exec --min 80      # the 80 % target, as functions entered by the tests
cyrius distlib --all && cyrius distlib --check
cyrius deps --locked
cyrius api-surface                   # the public API against its snapshot
cyrius audit                         # read its exit status; never pipe it through `tail`
python3 -I -B scripts/surface_parity.py .           # every ported Rust item has its Cyrius name
python3 -I -B scripts/rust_test_ledger.py . --check  # every Rust test ported, pending or recorded
python3 -I -B scripts/bench_record.py . --check      # every Rust bench has a twin, or is pending
python3 -I -B scripts/doc_examples.py .              # the README's and guides' cyrius blocks run
```

`cyrius test` must be **warning-free**. The f64 type diagnostics come only from the compiles that
`test` and `build` drive, never from `cyrius lint`. Each one is a bug:
- *integer arithmetic with an f64 right operand*
- *f64 arithmetic with a non-f64 right operand*
- *an integer stored into an f64 slot*
- *comparison mixes f64 and integer*

Run the cross-target checks (x86_64 PE under wine, aarch64 under qemu) **one at a time**, and
check that each prints its `=== parity summary` line: wine under parallel load truncates output
and still exits 0.

### Work Loop / Working Loop (continuous)

1. Work phase — one module or Part at a time, in the order the port plan gives
2. Cleanliness check
3. Test + benchmark additions: golden parity rows, trajectory rows, reference tests, bench twins
4. Run benchmarks, recorded against the Rust numbers. Before any benchmark run, check `ps` for
   other repos' gates and benches, not only the load average
5. Audit phase — correctness against `rust-old/`, the range-edge review (plan §0.1), memory (no
   per-step allocation), a surface-parity diff
6. Cleanliness check — must be clean after the audit
7. Deeper tests/benchmarks from the audit observations
8. Run benchmarks again — prove the wins
9. If the audit was heavy, return to step 5
10. Documentation — CHANGELOG, `state.md`, roadmap, the port plan's tracking table
11. Version check — `VERSION` and `cyrius.cyml` (`${file:VERSION}`) in sync
12. Return to step 1

### Task Sizing

- **Low/Medium effort**: Batch freely — multiple items per work loop cycle
- **Large effort**: Small bites only — break into sub-tasks, verify each before moving to the next. Never batch large items together
- **If unsure**: Treat it as large. Smaller bites are always safer than overcommitting

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks
- Never refactor speculatively. Wait for the third instance before extracting an abstraction
- Refactoring is part of the work loop, not a separate phase. If a review (step 5) reveals structural issues, refactor before moving to step 6
- Every refactor must pass the same cleanliness + benchmark gates as new code

## Key Principles

- **The port removes Rust; it doesn't carry Rust's faults.** Parity with `rust-old/` is
  functional: the same physics, results, errors and branches. An inconsistency or gap in the Rust
  crate is fixed in the port as an ADR 0002 row, with its oracle rows mapped exactly. "As Rust"
  is never a justification on its own.
- **Review the language before deferring anything.** Before a Rust surface is deferred or
  dropped, read the relevant sections of `cyrius/docs/guides/cyrius-guide.md` and
  `docs/stdlib-reference.md`. Then prove the claim with a probe program on the pinned compiler
  (`scripts/probes/`). Never take another port's limitation list on trust. The default is to
  port everything.
- **Cross-check against `rust-old/`** within the tolerance class the plan assigns (§7). A failure
  is a port defect, an upstream defect (file it), or an oracle defect (fix it in Rust 1.2.x and
  regenerate). **Never widen a class to make a row pass.**
- **Correct rounding is the accuracy target.** Every transcendental is correctly rounded and
  target-independent (`src/crmath.cyr`). Class-C rows compare against the corrected oracle
  (MPFR), never against glibc's misroundings.
- **One time convention: e^{+jωt}, IEEE handedness.** It is stated once, in `src/complex.cyr`'s
  header, and every module's doc points there.
- **Deterministic and bit-identical** on x86_64 ELF, x86_64 PE and aarch64. Never iterate a
  `hashmap` where order matters (its seed is random per process). Every sort key is total, or
  the sort is stable.
- **Never skip benchmarks.** Numbers don't lie. Speed is measured, not gated: a slow row is a
  finding to file with its numbers, not a reason to change the port.
- **Tests + benchmarks are the way.** Minimum 80 % coverage (`cyrius coverage --exec --min 80`).
- **Own the stack.** Depend on the Cyrius stdlib and AGNOS libraries, not on copies of them.
- **Enums.** Public enums have explicit, append-only numeric values that are never renumbered
  (the Cyrius form of `#[non_exhaustive]`).
- **`#must_use` on all pure functions**, plus `#pure` where it applies.
- **`#inline` on hot-path functions**, which applies only to functions of at most two
  parameters (the compiler ignores it otherwise, with a warning).
- **Optional dependencies go in profiles.** The default `[lib]` bundle never pulls sandhi/TLS;
  the `ai` client lives in `[lib.ai]`.
- **Logging on all operations**: sakshi, through the `_bj_warn` / `_bj_error` / `_bj_debug`
  helpers in `src/log.cyr`, behind a level.
- **No aborts in library code.** `src/` never uses `+?`, `-?` or `*?` (they exit the process),
  `assert_fatal`, `panic` or `result_unwrap`. Errors return a value-form `Result` with a
  `BJ_ERR_*` code, and the error's payload goes into the per-thread last-error record.
- **Memory.**
  - Grid and other large buffers come from `fl_calloc` and are released by `fl_free` inside
    `Bj*_free`.
  - No allocation inside any `step_*`.
  - **Never `alloc_reset`**: it wipes the consumer's heap.
- **f64 hygiene.**
  - Bind every call result to a `: f64` local before using it in an operator.
  - Compares on typed `f64` values are IEEE; compares on untyped words and builtin results are
    integer compares (use `f64_lt` etc.).
  - Never put hex in a `const`: it is an integer. Write `-x`, never `(0 - x)`, on an f64.
  - No exponent notation in literals: the lexer has none. Generated files write decimals.
  - Declare multi-returns `: (f64, f64)`. An undeclared `return (a, b)` gives untyped bindings.
  - **A captured `f64` is an untyped word inside a closure**: rebind it as `var a2: f64 = a;`
    in the body, or capture a struct with `f64` fields.
  - **Never write `s.f = g(t.e)`** when `g` returns a struct (an open compiler defect): bind the
    result to a local first.
  - Use the canonical `ganita_f64_*` names, not the deprecated short aliases.
- **Constants and tables are generated, never typed by hand.** `scripts/oracle/` writes them from
  the oracle. Regenerate them; never edit them.
- **Names.**
  - Public structs are `Bj<Name>`, with `impl` blocks, so methods keep their Rust names.
  - Free functions are `bj_<rust_name>`.
  - Globals and enum members carry `BJ_`; private helpers carry `_bj_`.
  - A `src/*.cyr` with public API is a `private` file that marks that API `public`. The bundle is one
    file, so the marker makes all of it private, and consumers see exactly the `public` items. A file
    of internals only (the math, fmtf) carries no marker.
  - A bundle shares one namespace with its consumers.
  - Tests, benches and fuzz units include `tests/lib.cyr`, never `[test] modules`: a prepended module
    joins the test's own file, which breaks `private` and execution coverage (port-parity.md §8).
- Build with `cyrius build`, not raw `cat file | cc5` — the manifest auto-resolves deps.
- Source files only need project includes — stdlib auto-resolves from `cyrius.cyml`.
- `var buf[N]` = N **bytes**, not N entries. Declare `var a: f64[N]` for elements.

## DO NOT

- **Do not commit or push** — the user handles all git operations (commit, push, tag)
- **NEVER use `gh` CLI** — use `curl` to the GitHub API only
- **Do not modify `rust-old/`** — it is the parity oracle. Tools that drive it (the golden and
  table generators) live in `scripts/oracle/` and depend on it by path
- Do not modify `lib/` files (vendored stdlib / dep symlinks)
- Do not hardcode toolchain versions in CI YAML — `cyrius = "X.Y.Z"` in `cyrius.cyml` is the source of truth
- Do not widen a parity tolerance class to make a row pass
- Do not defer or drop any Rust surface without the language review above and an ADR
- Do not keep a Rust inconsistency "as Rust"
- Do not add unnecessary dependencies
- Do not `panic` or abort in library code
- Do not `alloc_reset`
- Do not skip benchmarks before claiming performance improvements
- Do not skip tests before claiming changes work

## Documentation Structure

```
Root files (required):
  README.md          — quick start, features, dependency stack, consumers, license
  CHANGELOG.md       — per-version changes (Added/Changed/Fixed/Removed)
  CLAUDE.md          — this file (development process, principles, DO NOTs)
  CONTRIBUTING.md    — fork, branch, the cleanliness check, PR workflow
  SECURITY.md        — supported versions, scope, reporting
  CODE_OF_CONDUCT.md — Contributor Covenant
  LICENSE            — GPL-3.0-only

docs/ (required):
  architecture/README.md          — module map, data flow, consumers, dependency stack
  development/bijli-port-plan.md  — the port plan of record
  development/state.md            — live state
  development/roadmap.md          — planned work only; done work goes in the CHANGELOG
  development/port-parity.md      — measured parity baselines (M0 step 8)
  benchmarks-rust-v-cyrius.md     — the Rust vs Cyrius benchmark record
  adr/NNNN-title.md               — architecture decision records (template.md)

docs/ (when earned — not scaffolded empty):
  architecture/math.md            — mathematical reference for algorithms/formulas
  development/threat-model.md     — attack surface, mitigations (Touchstone and JSON readers)
  development/rust-test-ledger.md — every Rust test: ported, folded into golden rows, or n/a
  guides/getting-started.md, guides/migration-2.0.md — how-tos; the 2.0 migration (M8)

scripts/probes/   — the language-review probe programs, re-run at every pin bump
rust-old/docs/    — the Rust 1.2.0 architecture overview, kept for reference
```

ADR format: `# NNNN — Title`, then Status, Context, Decision, Consequences, Alternatives
considered (see `docs/adr/template.md`).

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/):

```markdown
# Changelog

## [Unreleased]
### Added — new features
### Changed — changes to existing features
### Fixed — bug fixes
### Removed — removed features
### Security — vulnerability fixes
### Performance — benchmark-proven improvements (include numbers)

## [X.Y.Z] - YYYY-MM-DD
### Added
- **module_name** — what was added and why
### Changed
- item: old behavior → new behavior
### Fixed
- issue description (root cause → fix)
### Performance
- benchmark_name: before → after (−XX%)
```

Rules:
- Every change that alters behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers (Rust vs Cyrius where relevant)
- Breaking changes get a **Breaking** section with a migration guide
- Group by module when multiple changes land in one release
- Link the ADR when a change was driven by an architectural decision
