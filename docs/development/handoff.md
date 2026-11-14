# Handoff: the Cyrius port, mid-M3 (2026-11-13)

For the next agent picking up the port. **Read first:** `CLAUDE.md` (the rules), [`state.md`](state.md) (live state),
[`bijli-port-plan.md`](bijli-port-plan.md) (the plan; its tracking table near the end), the latest sections of
[`port-parity.md`](port-parity.md) (§55–59 are M3's), and [ADR 0002](../adr/0002-port-api-divergences.md) (every
place the port differs from Rust, rows 1–32). This file is the working knowledge those don't hold: how a bite is
done, the scripts, and the traps. Delete it once it has served.

## Where it stands

- **Done and committed** (last commit `f81cc26`, "M3 part 1c-1"): M0, M1 and M2. In M3:
  - 1a-1: polarization (§55);
  - 1a-2: its JSON (§56);
  - 1b-1: scattering (§57);
  - 1b-2: MieResult's JSON (§58);
  - 1c-1: beam (§59).
- **Numbers at that commit**: 16 modules, 484,789 ok, 0 failed, 63,224 accepted (ADR 01: 24,164; ADR 0002: 39,060)
  over 548,013 rows; 374 digests equal on x86_64 ELF, x86_64 PE and aarch64. Also: ledger 522 ported, 3 n/a, 333
  pending of 858; surface parity 391 items matched; API snapshot 2,027 public functions; benches 61 of 80 twinned;
  coverage 282 of 282 functions.
- **Next bite, 1c-2**: the JSON of beam's three structs (`GaussianBeam`, `AbcdMatrix`, `ResonatorStability`; their
  six serde items sit in `scripts/surface_parity.py`'s `IN_PROGRESS`). Then **M3 Part 2**, the range-edge review: its
  candidates are listed at the end of §55, §57 and §59. Then M3 is done, and M4 (rf.rs) follows.
- **1c-2's design question**: `ResonatorStability.is_stable` is a `bool`, serialized `true`/`false`.
  - Neither `src/json.cyr` (f64 fields only) nor `scripts/oracle/src/m_json.rs`'s fixed-width `Kind` framework (f64
    values in every position) has a bool yet. `validate.py`'s `JSON_FIELDS` grammar needs a `"bool"` type.
  - serde reads only JSON `true`/`false` as a bool, never 0 or 1. Let the oracle's rows show serde_json's exact
    behaviour before writing the reader.
  - MieResult's own section in `m_json.rs` (`emit_mie`) shows how a type outside the `Kind` framework was added.
  - `GaussianBeam` and `AbcdMatrix` are plain f64 structs and fit `Kind` as Material did. The readers take the fields
    as written (serde's derive never calls `new`).

## How a bite is done

Each bite is one module's API (or its JSON, or one range-edge group), in this order. Rows are kept once per
(function, inputs), since the mapping files key on that.

1. **Read** the Rust module and its tests. Split it the way M3 did: the API first, the JSON second, range edges in
   the milestone's Part 2. Decide the convention mapping (ADR 0002 row 5) and any Rust fault (a new ADR 0002 row).
2. **Oracle**: write `scripts/oracle/src/m_<module>.rs`, and register it in `main.rs` (`mod` and `emit()`,
   file-name order).
   - A function is class C if it calls a glibc transcendental, else class A. A `Result` is its value, or
     `err(&e), s(e.to_string())`: the Display text is compared too (row 4).
   - Write the corpus with `cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml`, and verify it
     with the same plus `-- --check`.
   - Lint with `cargo clippy --release --locked --all-targets`. **Never `cargo fmt`**: there is no rustfmt config,
     and it rewraps ten files. Keep lines ≤ 120 characters.
   - **New JSON kinds go last** in `m_json.rs` (`Kind::ALL`, or after it), so the earlier rows and the seeded
     stream are unchanged. Check that the old `json.tsv` is an exact prefix of the new one.
3. **validate.py**: an independent transcription (`<x>_eval(fn, v, P, port=False)`), plus `validate_<module>`,
   `<module>_corrections` and `<module>_adr_rows`.
   - Wire them into `checks`, `check_corrections`, the `check_adr` module tuple and `adr_rows`.
   - Use `fdiv`, `fsqrt`, `powidf2_ieee` and `rust_display` / `rust_fixed` / `cx_fmt`. Python raises on x/0, and
     those helpers give IEEE's results instead.
   - The primitives come from `GLIBC_P` (the oracle) and `MPFR_P` (the corrected oracle). `cbrt` is correctly
     rounded in both: it is CORE-MATH's.
   - `python3 -I -B scripts/oracle/validate.py . --write` writes `<m>.cr.tsv` and `<m>.adr.tsv`; a plain run checks.
   - After a first clean run, prove the model compares something: flip a bit of each row's first output, and the
     run must catch every one.
4. **Port**: write `src/<module>.cyr`, a `private` file with `public` API and `_bj_` helpers. Add it to `cyrius.cyml`
   `[lib] modules` and `tests/lib.cyr` in the same order, before `src/json.cyr`.
5. **Harness**: in `tests/parity.tcyr`, add `_par_<module>` and its dispatch line.
   - Add `_par_c_bound` returning 0 for the class-C names (0 ulp against the corrected oracle).
   - For a convention mapping, add a counter and a `_par_module_done` line (see polarization's and beam's).
6. **Tests**: write `tests/<module>.tcyr`. Name each Rust test exactly as Rust does (the ledger matches by name), and
   add `cyrius_*` groups for the reasons' exact text, out untouched by an error, the convention's sense, and the
   heap. Then `python3 -I -B scripts/rust_test_ledger.py .` and the same with `--check`.
7. **Surfaces**: update `scripts/surface_parity.py`. A module in progress is in `IN_PROGRESS`, not `MODULES`; its
   serde items are keyed `"impl Serialize for X"`. Then run `cyrius api-surface --update`.
8. **Benches**: add twins in `tests/bijli.bcyr`, named as Rust's criterion ids.
   - A Cyrius-only bench must be in `scripts/bench_record.py`'s `EXTRA`.
   - Fill `docs/benchmarks-rust-v-cyrius.md`'s rows **by hand**, and its header count. Passing a log rewrites
     every row.
   - Compile-check with one `cyrius bench` before the A/B (variables are function-wide; a duplicate name fails the
     compile).
   - Run the A/B against HEAD (below), interleaved, three rounds, and note the load. Never wait for a quiet machine.
9. **Gates** (below), then **cross-target**, one at a time:
   - `cyrius test --win tests/parity.tcyr`, then `cyrius test --aarch64 tests/parity.tcyr`.
   - Check each run's `=== parity summary` line, and the module's own test file on both.
   - Compare the digests: `grep digest | tr -d '\r' | sort`, equal across ELF, PE and aarch64, and HEAD's all
     unchanged.
10. **Docs**: in `CHANGELOG.md`, a Changed entry per bite, and a Fixed entry per new ADR 0002 row (before
    `### Performance`). Then:
    - `port-parity.md`: a new §N with the result line;
    - `state.md`'s milestone row;
    - the plan's tracking row;
    - `scripts/oracle/README.md`'s module table (alphabetical);
    - ADR 0002's rows and amendments.

    Escape `|` inside table cells as `\|`. Check each table's pipe count and the 120-character width (tables are
    exempt). The user commits each bite; never commit or push.

### The gate script

```bash
#!/bin/bash
# The cleanliness check: every step's tail, its rc, and its warning lines (lint's real count is its
# "N warnings total" line; this grep also counts each file's "0 warnings").
cd /home/macro/Repos/bijli
S=${S:?set S to a scratch directory}
F="src/*.cyr tests/*.tcyr tests/*.bcyr tests/*.fcyr tests/*.cyr scripts/oracle/cr_driver.cyr"
run() { echo "### $*"; eval "$@" > $S/g.out 2>&1; rc=$?; tail -n 4 $S/g.out; echo "rc=$rc"; }
run cyrius fmt --check $F
run cyrius lint $F
run cyrius vet src/main.cyr
run cyrius build src/main.cyr build/bijli
run cyrius test
run cyrius test --bounds
run cyrius coverage --exec --min 80
run "cyrius distlib --all && cyrius distlib --check"
run cyrius deps --locked
run cyrius api-surface
run python3 -I -B scripts/surface_parity.py .
run python3 -I -B scripts/rust_test_ledger.py . --check
run python3 -I -B scripts/doc_examples.py .
run python3 -I -B scripts/oracle/validate.py .
run "cyrius audit > $S/audit.txt 2>&1"
run python3 -I -B scripts/bench_record.py . --check
run cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check
```

Run it as a background job (about four minutes), and grep `warning:` in the output: `cyrius test` must be
warning-free. `cyrius build src/main.cyr` builds only the scaffold; the library compiles through the tests.

### The A/B bench

```bash
S=/path/to/scratch; rm -rf $S/head && mkdir -p $S/head && git archive HEAD | tar -x -C $S/head && cp -a lib $S/head/
for r in 1 2 3; do uptime > $S/ab_load_$r.txt
  (cd $S/head && cyrius bench > $S/ab_head_$r.txt 2>&1)
  (cd /home/macro/Repos/bijli && cyrius bench > $S/ab_new_$r.txt 2>&1); done
```

Compare the per-bench medians of the three rounds, and check `ps` for other repos' gates first. A row bimodal in both
builds (`relativity/em_tensor_from_fields` sits at 31 or 34 ns) is layout noise. Costs go in the CHANGELOG's Fixed
entries; the Performance section is for proven improvements only.

## Traps found in M3 (Cyrius 6.7.18)

- **f64 literals and locals**:
  - There are no exponent literals: write decimals (generate long ones by script; Python's `Decimal(repr(x))`).
  - A hex bit pattern may initialize a typed f64 **global**. Into an f64 **local** or an argument it warns ("an
    integer stored into an f64 slot"); write `_bj_from_bits(0x…)` there.
- **Variables**: a `var` is function-wide. A name repeated in another block or loop of the same function is a
  compile error, and a `var` declared in a loop body is not visible after the loop.
- **Syntax**:
  - There is no `!` on a comparison: use a flag or the opposite test.
  - `#inline` is honoured only for a single-expression body (no `var`, no control flow) of at most two parameters.
- **Visibility**:
  - A `private` file's helpers are invisible to other files.
  - A Rust-private function its tests reach becomes `public fn _bj_…` (precedent: `_bj_associated_legendre`,
    `_bj_hermite_polynomial`).
- **Structs and Results**:
  - A `Result` with a struct or pair payload writes through `out` pointers, and leaves `out` untouched on an error.
  - A pair returns `: (f64, f64)`.
  - Inside an `impl` method `self` is `*T`: pass it as `self`, copy it with `*self`.
  - An array of structs can't be subscripted: use addresses and `load64`/`store64`, or an f64 array and `memcpy` (all
    four polarization structs lie in serde's field order).
- **Owned buffers**: buffers (`BjMieResult`'s coefficients) come from `fl_alloc`. The producer writes new ones into
  `out` and never reads it; `Bj*_free` gives them back. Tests pin the heap with `_fl_arena_pos` and `alloc_used()`.
- **Exactness**: Rust's `Complex::from_polar` compiles to one sincos, so call `_bj_sin_cos` once. That gives the same
  correctly rounded bits as two calls, at half the cost.
- **Tools**:
  - The editor's diagnostics on a lone `src/*.cyr` ("unknown type BjComplex") are noise: the types come from the
    bundle.
  - `cyrius fmt --check` wants two spaces of continuation indent per open parenthesis.
  - `cyrius lint`'s rc is 0 with warnings: read its "warnings total" line.

## Standing preferences (also in the memory directory)

- Full rigor: no speed at the cost of bugs.
- Review the Cyrius language (its guide and a probe on the pinned compiler) before deferring anything; port
  everything.
- The port removes Rust's faults: each one becomes an ADR 0002 row with exact oracle mappings. Never keep one "as
  Rust", and never widen a tolerance.
- Never propose side tasks or separate sessions: findings go in the docs, and the port goes on.
- Run long jobs in the background with commands that exit. Never write a polling loop. Don't touch the tree
  while a background job uses it.
- Never `gh`, never commit or push: the user commits each bite, sometimes under a loose title.
