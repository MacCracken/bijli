# Changelog

All notable changes to bijli will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased] — V0.2

### Added

- **field**: Electric/magnetic dipole fields (full 3D), dipole potential
- **field**: Gauss's law applications — infinite plane, charged sphere (inside/outside), infinite line, cylinder (inside/outside)
- **field**: Charge distributions — ring (axial), disk (axial)
- **field**: Field line tracing with Euler integration
- **field**: `std::ops` traits (`Add`, `Sub`, `Mul<f64>`, `Neg`) and `PartialEq` on `FieldVector`

### Changed

- **charge**: `PointCharge::new` now validates mass > 0 (returns `Result`)
- **maxwell**: `refractive_index` now validates inputs (returns `Result`)
- **maxwell**: `wavelength`/`frequency` now validate velocity > 0
- **ai**: `register_agent` returns error instead of silent fallback on missing `agent_id`

### Added (V0.3)

- **circuit**: New module — Ohm's law (V/I/R/P), resistance from geometry, series/parallel combinations (R, C, L), parallel plate capacitance, solenoid inductance, capacitor/inductor energy
- **circuit**: RC circuits — time constant, charging/discharging voltage
- **circuit**: RL circuits — time constant, current rise/decay
- **circuit**: RLC circuits — resonant frequency, quality factor, damping ratio, impedance
- **material**: New module — electric susceptibility, relative/absolute permittivity, polarization, displacement field, bound surface/volume charges, dielectric energy density, Clausius-Mossotti relation
- **material**: Magnetic materials — susceptibility, relative/absolute permeability, magnetization, H/B field conversions, bound surface current, magnetic energy density
- **material**: Curie's law, Curie-Weiss law, magnetic type classification (dia/para/ferromagnetic)

### Fixed

- Removed unused `serde` imports from `wave` and `maxwell` modules
- Removed unused `SPEED_OF_LIGHT` import from `maxwell` module scope

### Performance

- `inv_r` multiply pattern replacing sqrt+divide — 10-29% faster on field/charge computations
- `#[inline]` on all public functions
- `superposition_10`: 64ns → 41ns (-36%)

## [0.1.0] — 2026-03-24

### Added

- **field**: Electric/magnetic vector fields, point charge fields, potentials, superposition, energy density, constants (ε₀, μ₀, c, k)
- **maxwell**: Gauss's laws, displacement current, wave speed, refractive index, impedance, wavelength/frequency, skin depth
- **charge**: Point charges (electron, proton), Coulomb force, Lorentz force, dipole moment, potential energy, cyclotron frequency, Larmor radius
- **wave**: Poynting vector, plane wave intensity, radiation pressure, E/B amplitude conversion, wave number, momentum density
- **error**: BijliError with domain-specific variants (DimensionMismatch, ZeroCharge, Singularity, InvalidPermittivity, InvalidPermeability)
- **ai**: Daimon/hoosh client integration (feature-gated)
- **logging**: Structured logging via BIJLI_LOG (feature-gated)
- Infrastructure: CI/CD, deny.toml, codecov, benchmarks, Makefile
