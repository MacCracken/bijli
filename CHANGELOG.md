# Changelog

All notable changes to bijli will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Changed

- **license**: GPL-3.0 → GPL-3.0-only (no "or later" clause)
- **deps**: hisab 1.1.0 → 1.3.0
- **polarization**: `Complex` type now re-exported from `hisab::num::Complex`; added `ComplexExt` trait and `complex_real`/`complex_zero`/`complex_from_polar` free functions
- **relativity**: `solve_retarded_time` now uses `hisab::num::newton_raphson` for root-finding

### Removed

- **deps**: Removed unused `impetus` dependency
- **polarization**: Removed custom `Complex` implementation (~185 lines) in favor of hisab's

## [1.0.0] — 2026-03-25

### Added

**New modules:**
- **polarization**: Jones vectors/matrices (coherent light), Stokes parameters, Mueller matrices (partially polarized light), optical elements (polarizer, waveplate, retarder, rotator, depolarizer), Complex number type
- **scattering**: Mie scattering (Bohren-Huffman algorithm, a_n/b_n coefficients, Q_sca/Q_ext/Q_abs, asymmetry parameter g), Rayleigh scattering (λ⁻⁴ cross-section, phase function)
- **beam**: Gaussian beam parameter q(z), ABCD ray transfer matrices, resonator stability (g-parameters), Hermite-Gaussian and Laguerre-Gaussian mode profiles
- **rf**: Transmission lines (lossless/lossy Z₀, propagation constant), VSWR, return/insertion loss, Smith chart, matching networks (quarter-wave, single-stub), S-parameters (N-port, cascade, S↔Z↔T conversions), Touchstone I/O, antenna arrays (linear, planar, circular)
- **compute**: `ComputeBackend` trait (vendor-agnostic GPU/CPU), `CpuBackend` reference impl, `GpuFdtd2d<B>` generic solver, `SooratBackend`/`WgpuBackend` scaffolds

**Existing module additions:**
- **wave**: Trig-free Fresnel equations (`fresnel_rs_direct`, `fresnel_rp_direct`, `snell_cos_theta_t`), Schlick approximation, material-based Fresnel interface
- **fdtd**: 2D FDTD (TM/TE modes), 3D FDTD (full vector Yee grid, 6 field components), CPML absorbing boundaries, source functions (Gaussian pulse, Ricker wavelet, modulated Gaussian, TFSF), dispersive materials (Drude/Lorentz/Debye with ADE coefficients, gold/silver/water presets), NF2FF transform, subpixel smoothing, `NonUniformGrid`, 4th-order spatial stencil, `CwAccumulator`, `find_eigenmodes`
- **field**: Multipole expansion (quadrupole potential, spherical harmonic moments, associated Legendre polynomials), method of images (conducting plane, grounded sphere, dielectric half-space), Green's functions (scalar complex + static), Maxwell stress tensor, Poynting flux surface integration
- **relativity**: Retarded time solver (Newton-Raphson), convenience `lienard_wiechert_fields()` combining retarded time + E + B
- **charge**: `PointCharge` convenience methods — `electric_field_at`, `electric_potential_at`, `magnetic_field_at`, `lorentz_force_in`
- **material**: Unified `Material` struct — vacuum/dielectric/conductor constructors, refractive index, impedance, wave speed, reflectance

### Changed

- **circuit**: `resistance_series`, `capacitance_parallel`, `inductance_series` now return `Result` and reject empty slices
- **fdtd**: Mur ABC coefficient precomputed in constructor; added tracing instrumentation
- **fdtd**: `subpixel_permittivity_parallel` now returns `Result` (validates permittivity > 0)
- **fdtd**: `Fdtd2d::set_material` and `Fdtd3d::set_material` accept unified `Material` struct
- **rf**: `input_impedance_lossless` now returns `Result` (catches zero denominator)
- **rf**: `normalize_impedance` now returns `Result` (validates Z₀ > 0)
- **field**: Added `Default`, `From<[f64;3]>`, `Into<[f64;3]>` for `FieldVector`; exposed `SPEED_OF_LIGHT_SQ`
- **field**: Removed redundant `FieldVector::add()` method (use `+` operator)
- **polarization**: Added `Default` for `Complex`, `JonesVector`, `StokesVector`; `From<f64>` and `Div` for `Complex`
- **polarization**: Removed redundant `StokesVector::add()` method (use `+` operator)

### Dependencies

- hisab 0.22.4 → 1.1.0
- impetus 0.23.3 → 1.0.0
- reqwest 0.12 → 0.13
- criterion 0.5 → 0.8
- Added: wgpu 25 (optional, `compute-wgpu` feature)
- Removed: soorat dependency

## [0.24.3] — 2026-03-24

### Added

- **field**: Electric/magnetic dipole fields (full 3D), dipole potential
- **field**: Gauss's law applications — infinite plane, charged sphere (inside/outside), infinite line, cylinder (inside/outside)
- **field**: Charge distributions — ring (axial), disk (axial)
- **field**: Field line tracing with Euler integration
- **field**: `std::ops` traits (`Add`, `Sub`, `Mul<f64>`, `Neg`) and `PartialEq` on `FieldVector`
- **relativity**: New module — Lorentz factor, beta, four-vectors, electromagnetic field tensor F^μν, Lorentz transformation of fields, retarded potentials, Liénard-Wiechert fields, Larmor radiation power
- **wave**: Snell's law, Fresnel equations, waveguides (rectangular/cylindrical), antenna patterns, Friis transmission
- **fdtd**: New module — 1D FDTD Yee algorithm with Mur ABC
- **circuit**: New module — Ohm's law, series/parallel combinations, RC/RL/RLC circuits
- **material**: New module — dielectric/magnetic materials, Clausius-Mossotti, Curie's law

### Changed

- **charge**: `PointCharge::new` now validates mass > 0
- **maxwell**: `refractive_index`, `wavelength`, `frequency` now validate inputs

### Performance

- `inv_r` multiply pattern — 10-29% faster field/charge computations
- `superposition_10`: 64ns → 41ns (-36%)

## [0.1.0] — 2026-03-24

### Added

- **field**: Electric/magnetic vector fields, point charge fields, potentials, superposition, energy density, constants (ε₀, μ₀, c, k)
- **maxwell**: Gauss's laws, displacement current, wave speed, refractive index, impedance, wavelength/frequency, skin depth
- **charge**: Point charges (electron, proton), Coulomb force, Lorentz force, dipole moment, potential energy, cyclotron frequency, Larmor radius
- **wave**: Poynting vector, plane wave intensity, radiation pressure, E/B amplitude conversion, wave number, momentum density
- **error**: BijliError with domain-specific variants
- **ai**: Daimon/hoosh client integration (feature-gated)
- **logging**: Structured logging via BIJLI_LOG (feature-gated)
