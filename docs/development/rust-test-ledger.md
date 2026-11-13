# Rust-test ledger

> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6 step 7, the
> mela pattern). Every test of the frozen Rust oracle (`rust-old/`, bijli 1.2.0) is ported to a Cyrius test
> group of the same name, recorded here with its reason, or pending until its milestone. The golden corpus
> (`tests/golden/`) checks every ported function on its grid besides; this ledger is about the oracle's own
> assertions.

**858 Rust tests:** 313 ported, 3 not applicable, 542 pending.

| Rust file | Milestone | Cyrius tests | Rust tests | Ported | N/A | Pending |
|---|---|---|---:|---:|---:|---:|
| `src/ai.rs` | M7 | `tests/ai.tcyr` | 10 | 0 | 0 | 10 |
| `src/beam.rs` | M3 | `tests/beam.tcyr` | 60 | 0 | 0 | 60 |
| `src/bridge.rs` | M7 | `tests/bridge.tcyr` | 18 | 0 | 0 | 18 |
| `src/charge.rs` | M1 | `tests/charge.tcyr` | 26 | 26 | 0 | 0 |
| `src/circuit.rs` | M2 | `tests/circuit.tcyr` | 51 | 51 | 0 | 0 |
| `src/compute.rs` | M6 | `tests/compute.tcyr` | 18 | 0 | 0 | 18 |
| `src/error.rs` | M1 | `tests/error.tcyr` | 10 | 10 | 0 | 0 |
| `src/fdtd.rs` | M6 | `tests/fdtd1d.tcyr`, `tests/fdtd_media.tcyr`, `tests/nf2ff.tcyr`, `tests/fdtd_util.tcyr`, `tests/cpml.tcyr`, `tests/fdtd2d.tcyr`, `tests/fdtd3d.tcyr` | 167 | 0 | 0 | 167 |
| `src/field.rs` | M1 | `tests/field.tcyr` | 105 | 105 | 0 | 0 |
| `src/integration/soorat.rs` | M7 | `tests/soorat.tcyr` | 14 | 0 | 0 | 14 |
| `src/logging.rs` | M1 | `tests/log.tcyr` | 1 | 1 | 0 | 0 |
| `src/material.rs` | M2 | `tests/material.tcyr` | 42 | 39 | 3 | 0 |
| `src/maxwell.rs` | M1 | `tests/maxwell.tcyr` | 18 | 18 | 0 | 0 |
| `src/polarization.rs` | M3 | `tests/polarization.tcyr` | 69 | 0 | 0 | 69 |
| `src/relativity.rs` | M2 | `tests/relativity.tcyr` | 48 | 0 | 0 | 48 |
| `src/rf.rs` | M4 | `tests/rf.tcyr` | 80 | 0 | 0 | 80 |
| `src/scattering.rs` | M3 | `tests/scattering.tcyr` | 32 | 0 | 0 | 32 |
| `src/wave.rs` | M2 | `tests/wave.tcyr` | 63 | 63 | 0 | 0 |
| `tests/integration.rs` | M8 | `tests/integration.tcyr` | 26 | 0 | 0 | 26 |

## Not applicable

| Rust test | Reason |
|---|---|
| `src/material.rs::test_material_reflectance_air_glass` | deprecated reflectance_at (ADR 0002 row 1); its 0.04 is checked through try_reflectance_at |
| `src/material.rs::test_material_reflectance_same` | Material::reflectance_at is deprecated in 1.2.0 and not carried (ADR 0002 row 1); tests/material.tcyr's cyrius_reflectance_via_try checks the same value through try_reflectance_at |
| `src/material.rs::test_reflectance_at_deprecated_never_wrong` | deprecated reflectance_at (ADR 0002 row 1): the port has no NaN-returning form; try_reflectance_at refuses the conductor (cyrius_reflectance_via_try) |
