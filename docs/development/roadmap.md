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

- [ ] Electromagnetic tensor
- [ ] Lorentz transformations of fields
- [ ] Retarded potentials
- [ ] Liénard-Wiechert fields

## V1.0 — Stable Release

- [ ] API review and stabilization
- [ ] Feature audit
- [ ] Complete documentation with physics references
- [ ] Performance optimization pass

## Consumer Mapping

| Consumer | Uses |
|----------|------|
| kiran | EM effects in game physics |
| joshua | EM simulation mode |
| prakash | Optics builds on EM foundations |

## Boundaries

| Domain | Owned by | NOT bijli |
|--------|----------|-----------|
| PDE solvers | hisab | bijli uses hisab for numerical methods |
| Particle dynamics | impetus | bijli computes forces, impetus integrates |
| Optics (ray, wave) | prakash | bijli provides EM foundations |
| Quantum EM (QED) | future | bijli is classical EM |
