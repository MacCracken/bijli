# Contributing to bijli

Thank you for considering a contribution! This document covers the development
workflow, coding standards, and review process.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/MacCracken/bijli.git
cd bijli

# The Cyrius toolchain is pinned in cyrius.cyml ([package].cyrius)
cyrius deps

# Build and test
cyrius build src/main.cyr build/bijli
cyrius test
```

## Pull Request Process

1. **Fork and branch** — create a feature branch from `main`.
2. **Keep commits focused** — one logical change per commit.
3. **Write tests** — new features require tests; bug fixes require a regression
   test.
4. **Run the cleanliness check** (CLAUDE.md) before pushing: `cyrius fmt --check`, `cyrius lint`, `cyrius vet`,
   `cyrius test` (warning-free), `cyrius test --bounds`, `cyrius coverage --exec --min 80`, `cyrius distlib --check`,
   `cyrius api-surface`, `cyrius audit`, and the port's Python checks (`scripts/surface_parity.py`,
   `scripts/rust_test_ledger.py --check`, `scripts/bench_record.py --check`, `scripts/oracle/validate.py`).
5. **Open a PR** against `main` with a clear description of the change.
6. **Address review feedback** — maintainers may request changes before merging.

## Code Style

- `cyrius fmt` and `cyrius lint` clean; lines at most 120 characters.
- Public API items have a doc comment, with units for physics.
- `#must_use` on pure functions, `#inline` on hot-path functions of at most two parameters.
- Public enums have explicit, append-only numeric values.
- No Rust: no `.rs` file, Cargo manifest or `cargo` step ([ADR 0003](docs/adr/0003-retire-the-rust-tree.md)).

## Testing Requirements

- All public API changes include unit tests (`tests/*.tcyr`) and golden rows (`tests/golden/`, written by a Python
  emitter in `scripts/oracle/` from validate.py's model).
- Benchmarks go in `tests/*.bcyr`; run an interleaved A/B against the base before and after performance-sensitive
  changes.
- Execution coverage stays at 80 % or above.

```bash
cyrius test                          # tests and the golden parity harness
cyrius bench                         # benchmarks
cyrius coverage --exec --min 80      # execution coverage
```

## Commit Messages

Use clear, imperative-mood commit messages:

```
add Sellmeier dispersion coefficients for sapphire

fix NaN handling in wavelength_to_rgb boundary check
```

## License

By contributing, you agree that your contributions will be licensed under the
same license as the project (see [LICENSE](LICENSE)).
