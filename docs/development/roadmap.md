# Bijli — Development Roadmap

> Electromagnetism simulation for AGNOS

## V0.1 — Foundation (current)

- [x] Electric field from point charges (Coulomb's law)
- [x] Electric potential
- [x] Magnetic field from moving charges (Biot-Savart)
- [x] Field superposition
- [x] Lorentz force
- [x] Maxwell's equations (Gauss, wave speed, impedance)
- [x] EM wave propagation (Poynting, intensity, radiation pressure)
- [x] Physical constants and charge dynamics
- [x] Error types, logging, AI integration scaffold

## V0.2 — Field Computation

- [x] Charge distributions (line, surface, volume)
- [x] Gauss's law applications (infinite plane, sphere, cylinder)
- [x] Electric dipole field (full 3D)
- [x] Magnetic dipole field
- [x] Field line tracing

## V0.3 — Circuits & Materials

- [x] Ohm's law, resistance, capacitance, inductance
- [x] RC, RL, RLC circuits
- [x] Dielectric materials (polarization, bound charges)
- [x] Magnetic materials (diamagnetic, paramagnetic, ferromagnetic)

## V0.4 — Wave Propagation

- [x] EM wave in media (refraction, reflection, transmission)
- [x] Waveguides (rectangular, cylindrical)
- [x] Antenna patterns (dipole, monopole)
- [x] FDTD (Finite-Difference Time-Domain) solver

## V0.5 — Relativistic Electrodynamics

- [x] Electromagnetic tensor
- [x] Lorentz transformations of fields
- [x] Retarded potentials
- [x] Liénard-Wiechert fields

## V1.0 — Stable Release (0.24.3)

- [x] API review and stabilization
- [x] Feature audit
- [x] Complete documentation with physics references
- [x] Performance optimization pass

---

## V1.1 — Optics & Scattering Foundations

> Unlocks downstream: prakash (optical systems), kiran (sky rendering, volumetric effects)

- [ ] Trig-free Fresnel equations (`fresnel_rs_direct`, `fresnel_rp_direct` taking cos θ_i only)
- [ ] Polarization formalism — Jones vectors and Jones matrices (coherent light)
- [ ] Stokes parameters and Mueller matrices (partially polarized light)
- [ ] Standard optical elements: polarizer, waveplate, retarder matrices
- [ ] Mie scattering — coefficients a_n/b_n, Q_sca, Q_ext, Q_abs, asymmetry parameter g
- [ ] Rayleigh scattering — small-particle limit, λ⁻⁴ dependence, cross-sections
- [ ] Scattering cross-section, absorption cross-section, phase function

## V1.2 — Gaussian Beams & Paraxial Optics

> Unlocks downstream: prakash (laser/beam optics, resonator design)

- [ ] Complex beam parameter q(z) — waist, Rayleigh range, divergence, spot size
- [ ] ABCD ray transfer matrices — free space, thin lens, curved mirror, flat interface
- [ ] ABCD propagation: q_out = (Aq_in + B)/(Cq_in + D)
- [ ] Resonator stability analysis (g-parameters)
- [ ] Higher-order Hermite-Gaussian and Laguerre-Gaussian modes

## V1.3 — FDTD Maturity

> Unlocks downstream: joshua (production EM simulation), kiran (real-time EM viz)

- [ ] 2D FDTD on Yee grid (TM and TE modes)
- [ ] CPML (Convolutional PML) absorbing boundaries — replace Mur ABC
- [ ] Gaussian pulse source (broadband excitation)
- [ ] Total-field/scattered-field (TFSF) plane wave source
- [ ] Lossy/conductive materials (conductivity term in update equations)
- [ ] Dispersive materials — Drude model (metals/plasmas)
- [ ] Dispersive materials — Lorentz oscillator model (dielectric resonances)
- [ ] Dispersive materials — Debye model (polar molecules)
- [ ] Auxiliary Differential Equation (ADE) method for dispersive FDTD
- [ ] Near-field to far-field (NF2FF) transform
- [ ] Subpixel smoothing at material boundaries

## V1.4 — RF & Microwave

> Unlocks downstream: joshua (RF circuit simulation)

- [ ] Transmission line theory — characteristic impedance, propagation constant, attenuation
- [ ] Telegrapher's equations
- [ ] VSWR, return loss, insertion loss from reflection coefficient
- [ ] Smith chart computations (impedance normalization, SWR circles, matching)
- [ ] Quarter-wave transformer, stub matching
- [ ] S-parameter (scattering) matrix for N-port networks
- [ ] Parameter conversions: S ↔ Z ↔ Y ↔ ABCD ↔ T
- [ ] Touchstone file format (.s1p, .s2p, .snp) read/write
- [ ] Cascade/de-embed operations
- [ ] Antenna array factor (linear, planar, circular arrays)

## V1.5 — Advanced Field Computation

> Unlocks downstream: joshua (boundary value problems, precision simulation)

- [ ] Multipole expansion — quadrupole (electric + magnetic) fields and potentials
- [ ] Octupole and general spherical harmonic expansion
- [ ] Method of images — point charge near conducting plane
- [ ] Method of images — point charge near conducting sphere
- [ ] Method of images — dielectric half-space
- [ ] Free-space scalar Green's function: G(r,r') = e^{ikR}/(4πR)
- [ ] Dyadic Green's function for vector fields
- [ ] Maxwell stress tensor — force/torque on objects from EM fields
- [ ] Poynting flux surface integration (total power through a surface)
- [ ] Retarded time solver (Newton-Raphson on light-cone equation)
- [ ] Convenience: `lienard_wiechert_fields()` combining retarded time + E + B

## V1.6 — Cross-Module Integration

> Improves ergonomics for all downstream consumers

- [ ] `PointCharge::electric_field_at(pos)` → delegates to field module
- [ ] `PointCharge::lienard_wiechert_field_at(pos, time)` → delegates to relativity
- [ ] `wave` Fresnel functions accepting material structs from `material` module
- [ ] `fdtd.set_material(region, material_params)` using `material` module types
- [ ] `circuit` ↔ `wave` bridge for transmission line reflection/Fresnel overlap
- [ ] Unified `Material` struct spanning `material` + `fdtd` + `wave`

## V2.0 — 3D FDTD & GPU

> Unlocks downstream: joshua (production 3D simulation), kiran (real-time EM)

- [ ] 3D FDTD on full vector Yee grid
- [ ] Trait-based solver backend (CPU/GPU agnostic)
- [ ] GPU compute backend via wgpu (cross-platform)
- [ ] Non-uniform grid spacing (first step toward AMR)
- [ ] Higher-order FDTD (4th-order spatial derivatives)
- [ ] Frequency-domain (CW) solver
- [ ] Eigenmode solver for resonant structures

## Future — Research Watch

> Not scheduled — track for when they reach practical maturity

- [ ] Physics-informed neural networks (PINNs) as surrogate EM solvers (17× speedup published 2025)
- [ ] Adaptive mesh refinement (AMR) for FDTD — h/p/r refinement, subgridding
- [ ] Spectral element methods for smooth geometries
- [ ] Discrete dipole approximation (DDA) for non-spherical scattering
- [ ] SR-FDTDNet — super-resolution reconstruction of fine-grid FDTD from coarse runs

---

## Consumer Mapping

| Consumer | Uses | Unlocked By |
|----------|------|-------------|
| kiran | EM effects in game physics, sky rendering, volumetric scattering, beam effects | V1.1 (scattering), V1.3 (2D FDTD), V2.0 (GPU) |
| joshua | Full EM simulation, RF circuits, antenna design | V1.3 (FDTD), V1.4 (RF), V1.5 (advanced fields), V2.0 (3D) |
| prakash | Optical systems, polarization, beam propagation, material optics | V1.1 (scattering/polarization), V1.2 (Gaussian beams) |

## Boundaries

| Domain | Owned by | NOT bijli |
|--------|----------|-----------|
| PDE solvers | hisab | bijli uses hisab for numerical methods |
| Particle dynamics | impetus | bijli computes forces, impetus integrates |
| Optics (ray, wave) | prakash | bijli provides EM foundations; prakash builds ray tracing, lens design, etc. |
| Quantum EM (QED) | future | bijli is classical EM |
| GPU compute | wgpu (external) | bijli provides the algorithms; GPU backend is a pluggable trait |
| ML surrogates | daimon/hoosh | bijli defines solver interfaces; AI agents provide trained surrogates |
