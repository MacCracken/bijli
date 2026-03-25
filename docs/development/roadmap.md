# Bijli — Development Roadmap

> Electromagnetism simulation for AGNOS

## Completed (V1.0.0)

All milestones through V2.0 are shipped. See [CHANGELOG.md](/CHANGELOG.md) for details.

**V0.1–V0.5**: Foundation — fields, Maxwell's equations, charge dynamics, waves, circuits, materials, FDTD, relativistic electrodynamics.

**V1.0**: API stabilization, feature audit, performance optimization.

**V1.1**: Optics & scattering — trig-free Fresnel, polarization (Jones/Stokes/Mueller), Mie/Rayleigh scattering.

**V1.2**: Gaussian beams — complex beam parameter, ABCD matrices, resonator stability, HG/LG modes.

**V1.3**: FDTD maturity — 2D Yee grid (TM/TE), CPML, sources (Gaussian/Ricker/TFSF), dispersive materials (Drude/Lorentz/Debye), NF2FF, subpixel smoothing.

**V1.4**: RF & microwave — transmission lines, VSWR, Smith chart, S-parameters, Touchstone I/O, antenna arrays.

**V1.5**: Advanced fields — multipole expansion, method of images, Green's functions, Maxwell stress tensor, Poynting flux, retarded time solver.

**V1.6**: Cross-module integration — PointCharge convenience methods, unified Material struct, FDTD/wave/circuit bridges.

**V2.0**: 3D FDTD & GPU — full vector Yee grid, ComputeBackend trait (CPU/soorat/wgpu pipes), non-uniform grids, 4th-order FDTD, CW solver, eigenmode solver.

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
| GPU compute | soorat/wgpu | bijli provides the trait; backends are pluggable |
| ML surrogates | daimon/hoosh | bijli defines solver interfaces; AI agents provide trained surrogates |
