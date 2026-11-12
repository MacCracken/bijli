# 0002 — The Cyrius 2.0.0 API diverges from Rust 1.2.0 where the Rust crate was wrong or incomplete

**Status**: Accepted
**Date**: 2026-11-11

## Context

[ADR 0001](0001-cyrius-port-with-rust-oracle.md) freezes Rust 1.2.0 as the parity oracle. Parity is
*functional*: the same physics, results, errors and branches. The port removes Rust; it doesn't
carry the Rust crate's faults. bijli's roadmap already lists a "2.0 cleanup (breaking)" that was
never released in Rust, and planning the port found more:
- the crate mixes two harmonic time conventions;
- it has no transport error;
- errors carry payloads that a one-word `Result` can't hold.

Every change here keeps its golden comparison: each oracle row is mapped exactly to what the
Cyrius API returns, and the harness reports each mapping by this ADR's row number. No tolerance
class widens.

## Decision

The Cyrius 2.0.0 API makes these changes. The table is append-only: rows found during the port
(the range-edge and exact-value fixes of plan §0.1) are added below with the milestone that
found them.

| # | Rust 1.2.0 | Cyrius 2.0.0 | Golden comparison |
|---|---|---|---|
| 1 | `rf::input_impedance_lossless`, `rf::write_touchstone` and `Material::reflectance_at`, deprecated in 1.2.0 | removed | through their replacements (`input_impedance_lossless_complex`, `try_write_touchstone`, the fallible `Material` Fresnel functions). Their Rust tests are n/a in the ledger |
| 2 | public FDTD grid fields (`ex`, `ez`, `nx`, …), with a per-step layout check because callers can resize them | buffers fixed at construction, read through accessors (pointer + length). Nothing can resize them, so the per-step layout check and its error path go | class T on the accessors. The Rust layout-mismatch tests are n/a |
| 3 | `TfsfSource::new` and `Cpml::new` sanitize with a warning; `Fdtd2d::enable_cpml` warns and ignores | fallible only: the constructors return `Err`, and `enable_cpml` is the `try_` form | each oracle sanitize-or-ignore row ↔ the `Err` code, mapped by name |
| 4 | `BijliError` variants carry data (`DimensionMismatch { expected, got }`, `InsufficientResolution { cells, domain_size }`, `InvalidPermittivity { value }`, …) and `String` reasons | append-only `BJ_ERR_*` codes in a value-form `Result`. The payload goes into a per-thread last-error record (`bj_last_error_*`): the variant's fields, and Rust's Display text of the error byte for byte, its reason or context included (`bj_last_error_display`, `bj_last_error_detail`). *Amended 2026-11-11 (M1 part 2a)*: the record holds the text, so nothing is logged (the first draft sent the reason to a sakshi warning; Rust logs no error) | class A on the code, on each field and on the Display text (`tests/golden/error.tsv`) |
| 5 | two time conventions: `polarization`, `beam`, the field Green's function and the FDTD media use e^{−iωt} with Born & Wolf handedness; `scattering`, `rf` and the FDTD CW/NF2FF phasors use e^{+jωt} | **e^{+jωt} in every module, with IEEE handedness** (IEEE Std 145), stated once in `src/complex.cyr`. Absorbing indices are n − jκ | complex amplitudes ↔ `conj(oracle)`; phases ↔ `−oracle`; circular-state names by physical state, each pinned by a convention-free rotation-sense test; real observables against the oracle directly. All exact |
| 6 | results bit-identical only on the same platform | bit-identical across x86_64 ELF, x86_64 PE and aarch64, with correctly rounded transcendentals (ADR 0001). *Clarified 2026-11-11 (M1 part 2c)*: a NaN is the same NaN in the same places, its sign and payload the hardware's (IEEE 754 leaves them open: x86 generates 0xFFF8…, ARM 0x7FF8…; Rust guarantees no NaN bits either). Canonicalizing would cost every operation for no observable value: Display prints `NaN` and JSON `null` either way | every row exact on every target (any NaN matching any NaN), and the parity digests equal on all three, with each NaN folded as one (port-parity.md §9) |
| 7 | no transport error: the `ai` client reports network failures as `InvalidParameter` | `BJ_ERR_TRANSPORT`, appended | the `ai` client's oracle `InvalidParameter` rows for transport failures ↔ `BJ_ERR_TRANSPORT` |
| 8 | `logging::init` installs a process-wide `tracing` subscriber that filters every crate's events, set by `BIJLI_LOG` in `EnvFilter`'s grammar (M1 part 2a) | `bj_log_init` sets bijli's own gate and only raises sakshi's process-wide level, never lowers it: other libraries' logging is not bijli's to silence. `BIJLI_LOG` takes a level or `bijli=<level>` directives; other targets' directives are ignored, and a value with none that applies falls back to the default, as `EnvFilter::try_from_env` does on an invalid one. With no init, nothing is printed, as Rust with no subscriber | none: logging is outside the corpus. `tests/log.tcyr` pins the gate and the grammar |
| 9 | a failed allocation aborts the process: `trace_field_line`'s `Vec` growth, and the FDTD grids' `try_filled` (M1 part 3c) | `BJ_ERR_ALLOCATION`, appended (Display "allocation failed: {reason}"): `bj_trace_field_line` returns it when its `BjPolyline` can't grow, leaving the polyline empty; the FDTD grids follow in M5–M6. Library code never aborts (CLAUDE.md) | none: the corpus can't make an allocation fail. `tests/field.tcyr` covers the growth path; the error path has no reachable test (freelist refuses only past its 2 GiB request limit or an exhausted heap) |

## Consequences

- **Positive**
  - The Cyrius API is internally consistent. One convention means phases can move between
    modules without conjugation, and consumers in the engineering literature (FDTD, RF,
    antennas, beams) read the formulas as written.
  - Every error keeps its information without an allocation.
  - The deprecations and the sanitize-with-warning constructors are gone, as the roadmap intended.
- **Negative**
  - 2.0.0 is a breaking release for every 1.x consumer. `docs/guides/migration-2.0.md` (M8)
    gives the conversions: conjugate amplitudes from bijli 1.x's `polarization`, `beam`, Green's
    function and ε(ω); circular-state names now follow IEEE; errors are codes, with the payload
    read through `bj_last_error_*`.
  - The harness carries a mapping table, and every mapped row is reported by number.
- **Neutral**
  - Rows 1–3 come from the roadmap's "2.0 cleanup (breaking)" item, rows 5 and 7 from its
    backlog, row 4 from the Cyrius `Result` shape (plan §2), row 8 from sakshi's process-wide level
    (M1 part 2a), and row 9 from the rule that library code never aborts (M1 part 3c).

## Alternatives considered

- **Carry the Rust API unchanged and fix it in a Cyrius 3.0.** Rejected: it ships the faults
  into a new language and then breaks consumers twice.
- **Unify on e^{−iωt}** (Jackson, Born & Wolf, Bohren & Huffman). Rejected: `rf`, Mie and the
  FDTD phasors would change instead. Touchstone files would need conjugating at the file boundary,
  and Smith-chart coordinates would come out conjugated relative to every printed chart (plan §4.4).
- **Error payloads in a boxed enum.** Rejected: it allocates on the never-freed bump heap at every
  error.
