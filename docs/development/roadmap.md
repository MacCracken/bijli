# Bijli — Development Roadmap

> Electromagnetism simulation for AGNOS

## Completed (through 1.2.0)

All milestones through V2.0 are shipped. See [CHANGELOG.md](/CHANGELOG.md) for details.

**V0.1–V0.5**: Foundation — fields, Maxwell's equations, charge dynamics, waves, circuits, materials, FDTD, relativistic electrodynamics.

**V1.0**: API stabilization, feature audit, performance optimization.

**V1.1**: Optics & scattering — trig-free Fresnel, polarization (Jones/Stokes/Mueller), Mie/Rayleigh scattering.

**V1.2**: Gaussian beams — complex beam parameter, ABCD matrices, resonator stability, HG/LG modes.

**V1.3**: FDTD maturity — 2D Yee grid (TM/TE), CPML, sources (Gaussian/Ricker/TFSF), dispersive materials (Drude/Lorentz/Debye), NF2FF, subpixel smoothing.

**V1.4**: RF & microwave — transmission lines, VSWR, Smith chart, S-parameters, Touchstone I/O, antenna arrays.

**V1.5**: Advanced fields — multipole expansion, method of images, Green's functions, Maxwell stress tensor, Poynting flux, retarded time solver.

**V1.6**: Cross-module integration — PointCharge convenience methods, unified Material struct, FDTD/wave/circuit bridges.

**V2.0**: 3D FDTD & GPU — full vector Yee grid, ComputeBackend trait (CPU/wgpu pipes), non-uniform grids, 4th-order FDTD, CW solver, eigenmode solver.

**1.2.0 (crate release)**: Toolchain refresh (MSRV 1.99, hisab 1.4, wgpu 30) and a full physics/robustness audit — working 2D CPML (−88 dB), FDTD conductivity, correct Mie for high-index spheres, robust retarded-time solver, corrected Fresnel/Schlick/Touchstone/stub-matching/beam formulas, non-panicking public API, per-feature build isolation.

---

## Next — Cyrius port (2.0.0)

Plan of record: [bijli-port-plan.md](bijli-port-plan.md). Decision: [ADR 0001](../adr/0001-cyrius-port-with-rust-oracle.md), proposed.

- [ ] **Prep (plan §4)**:
  - decisions: D1 (cut 1.3.0) and D5 (JSON) accepted, D2 and D4 settled (ported), D3 settled as e^{+jωt} with IEEE handedness in every module (plan §4.4);
  - cut and tag **1.3.0**, the last Rust feature release and the frozen oracle (with the 3D CPML);
  - take a quiet bench baseline;
  - land the port-plan commit.
- [ ] **M0**: run `cyrius port`, which moves the Rust crate to `rust-old/`. Repair the scaffold. Build the oracle generator, the parity harness and the probes.
- [ ] **M1–M8**:
  - foundation;
  - media, waves, circuits and relativity;
  - optics;
  - RF;
  - FDTD core;
  - FDTD 2D/3D;
  - integration surfaces;
  - release 2.0.0.
- The "2.0 cleanup (breaking)" backlog item below lands in the Cyrius API only, as ADR 0002 (plan §3.3). There is no Rust 2.0.

---

## Cross-Crate Bridges (shipped in 1.1.0)

- [x] **`bridge.rs` module** — primitive-value conversions for cross-crate electromagnetism
- [x] **ushma bridge**: current density (A/m²), resistance (Ω) → Joule heating power (W); EM field energy density → temperature rise rate
- [x] **prakash bridge**: wavelength (nm) → refractive index in medium; E-field amplitude → optical intensity (W/m²)
- [x] **tanmatra bridge**: electron energy level (eV) → emission wavelength (nm); nuclear charge → Coulomb potential scaling
- [x] **dravya bridge**: electric field strength (V/m) → piezoelectric stress; magnetic flux density (T) → magnetostrictive strain

---

## Soorat Integration (rendering visualization — shipped in 1.1.0)

- [x] **`integration/soorat.rs` module** — feature-gated `soorat-compat`, visualization data structures for soorat to consume
- [x] **FDTD field slices**: Expose `FieldSlice2D` / `FieldSlice3D` structs from `Fdtd2d`/`Fdtd3d` grids (field component, dimensions, spacing, origin) for heatmap/volumetric rendering
- [x] **Field line traces**: Expose `FieldLineVisualization` wrapping `trace_field_line()` output as `Vec<[f64; 3]>` polylines with field magnitude per point
- [x] **Point charge visualization**: Expose `ChargeVisualization` with charge positions, magnitudes, and sign for particle rendering with field halos
- [x] **Radiation pattern**: Expose `RadiationPattern` wrapping `Nf2ffResult` (angles + gain) for polar plot / 3D balloon rendering
- [x] **Vector field sampling**: Expose `VectorFieldSample` — grid of `FieldVector` values at regular positions for arrow/streamline rendering

---

## Backlog — follow-ups from the 1.2.0 audit

- [x] **CPML for `Fdtd3d`** — Roden–Gedney CPML on all six faces (`try_enable_cpml`/`try_enable_cpml_with`/`disable_cpml`/`cpml`): 12 Ψ fields stored only in the face slabs, −80 dB reflection measured for a 10-cell layer around a 3D point source
- [ ] **CPML + materials for `GpuFdtd2d`** — the compute-backend solver has PEC walls and vacuum only
- [ ] **TFSF integration into `Fdtd2d`** — `TfsfSource` is a standalone auxiliary grid; add plane-wave injection on the TF/SF contour
- [ ] **Frequency-aware lossy Fresnel** — complex ε = ε − jσ/ω so conductors are supported (today they return `Err`)
- [ ] **Transport error variant** — `BijliError` has no network/transport variant; the `ai` client maps to `InvalidParameter`
- [ ] **Sub-femtometre point-charge separations** — field/potential/Coulomb functions treat r² < 1e-30 m² as singular
- [ ] **2.0 cleanup (breaking)** — remove the 1.2.0 deprecations (`rf::input_impedance_lossless`, `rf::write_touchstone`, `Material::reflectance_at`), make FDTD grid fields private behind accessors, make `TfsfSource::new`/`Cpml::new` fallible, unify the polarization (e^{−iωt}) and scattering (e^{+iωt}) time conventions
- [ ] **Real `WgpuBackend`** — currently a scaffold

---

## Future — Research Watch

> Not scheduled — track for when they reach practical maturity

- [ ] Physics-informed neural networks (PINNs) as surrogate EM solvers (17× speedup published 2025)
- [ ] Adaptive mesh refinement (AMR) for FDTD — h/p/r refinement, subgridding
- [ ] Spectral element methods for smooth geometries
- [ ] Discrete dipole approximation (DDA) for non-spherical scattering
- [ ] SR-FDTDNet — super-resolution reconstruction of fine-grid FDTD from coarse runs

---

## Consumer Mapping

| Consumer | Uses | Status |
|----------|------|--------|
| kiran | EM effects in game physics, sky rendering, volumetric scattering, beam effects | Unlocked |
| joshua | Full EM simulation, RF circuits, antenna design, 3D FDTD | Unlocked |
| prakash | Optical systems, polarization, beam propagation, material optics | Unlocked |

## Boundaries

| Domain | Owned by | NOT bijli |
|--------|----------|-----------|
| PDE solvers | hisab | bijli uses hisab for numerical methods |
| Particle dynamics | impetus | bijli computes forces, impetus integrates |
| Optics (ray, wave) | prakash | bijli provides EM foundations; prakash builds ray tracing, lens design, etc. |
| Quantum EM (QED) | future | bijli is classical EM |
| GPU compute | wgpu (interim) | bijli provides the trait; backends are pluggable — wgpu will be replaced by a custom GPU backend |
| ML surrogates | daimon/hoosh | bijli defines solver interfaces; AI agents provide trained surrogates |
