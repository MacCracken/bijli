# Changelog

All notable changes to bijli will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased] — V1.1–V2.0 Full Roadmap

### Added

- **wave**: Trig-free Fresnel equations — `fresnel_rs_direct`, `fresnel_rp_direct`, `snell_cos_theta_t`, reflectance/transmittance variants, Schlick approximation
- **polarization**: New module — Jones vectors (H, V, D, A, R, L, arbitrary linear), Jones matrices (polarizer, waveplate, QWP, HWP, rotator)
- **polarization**: Stokes parameters — construction from Jones, unpolarized/polarized states, DOP, ellipticity/orientation angles, incoherent addition
- **polarization**: Mueller matrices — polarizers (linear, circular), waveplates, rotator, depolarizer, attenuator, composition
- **polarization**: Complex number type for polarization math
- **scattering**: New module — Mie scattering (Bohren-Huffman algorithm), coefficients a_n/b_n, Q_sca/Q_ext/Q_abs, asymmetry parameter g
- **scattering**: Rayleigh scattering — cross-section (λ⁻⁴ dependence), efficiency, phase function, angular intensity
- **beam**: New module — Gaussian beam parameter q(z), waist, Rayleigh range, spot size, divergence, Gouy phase, peak intensity
- **beam**: ABCD ray transfer matrices — free space, thin lens, curved mirror, flat/curved interface, compose, beam propagation
- **beam**: Resonator stability — g-parameters, stability condition, beam waist in cavity, transverse mode spacing
- **beam**: Hermite-Gaussian (HG_mn) and Laguerre-Gaussian (LG_pl) mode profiles
- **fdtd**: 2D FDTD on Yee grid — TM and TE modes, source injection, dielectric regions, energy computation
- **fdtd**: CPML absorbing boundary infrastructure — polynomial-graded conductivity, auxiliary field allocation
- **fdtd**: Source functions — Gaussian pulse, Ricker wavelet, modulated Gaussian, TFSF plane wave source with 1D auxiliary grid
- **fdtd**: Dispersive material models — LossyMaterial (conductivity), DrudeMaterial (metals: gold/silver presets), LorentzMaterial (dielectric resonances), DebyeMaterial (polar molecules: water preset), all with ADE coefficients
- **fdtd**: NF2FF near-field to far-field transform for 2D TM mode
- **fdtd**: Subpixel smoothing — parallel (harmonic mean) and perpendicular (arithmetic mean) effective permittivity
- **rf**: New module — transmission lines (lossless/lossy Z₀, propagation constant, phase velocity, input impedance)
- **rf**: Reflection coefficient, VSWR, return loss, insertion loss, mismatch loss
- **rf**: Smith chart — impedance normalization, Γ↔Z conversion
- **rf**: Matching networks — quarter-wave transformer, single-stub matching
- **rf**: S-parameter matrix (N-port), cascade via T-parameters, S→Z conversion
- **rf**: Touchstone (.snp) file parse/write (RI, MA, DB formats)
- **rf**: Antenna arrays — linear, planar, circular array factors
- **field**: Multipole expansion — quadrupole potential, general spherical harmonic moments, associated Legendre polynomials
- **field**: Method of images — conducting plane, grounded sphere, dielectric half-space
- **field**: Green's functions — free-space scalar (complex, k-dependent) and static
- **field**: Maxwell stress tensor (3×3), stress-to-force computation, Poynting flux surface integration
- **relativity**: Retarded time solver (Newton-Raphson on light-cone equation)
- **relativity**: Convenience `lienard_wiechert_fields()` combining retarded time solver + E + B computation
- **charge**: `PointCharge` convenience methods — `electric_field_at`, `electric_potential_at`, `magnetic_field_at`, `lorentz_force_in`
- **material**: Unified `Material` struct — vacuum/dielectric/conductor constructors, refractive index, impedance, wave speed, reflectance, Default impl
- **fdtd**: `Fdtd2d::set_material` accepting `Material` struct for region setup
- **wave**: Material-based Fresnel — `reflectance_at_interface`, `reflectance_normal_materials`, `material_reflection_coefficient` bridging wave/circuit/material modules
- **compute**: New module — `ComputeBackend` trait (vendor-agnostic GPU/CPU), `CpuBackend` reference impl, `GpuFdtd2d<B>` generic solver, `SooratBackend`/`WgpuBackend` scaffolds
- **fdtd**: 3D FDTD — `Fdtd3d` with full vector Yee grid (Ex,Ey,Ez,Hx,Hy,Hz), 3D Courant condition, material setting
- **fdtd**: `NonUniformGrid` — per-cell spacing, geometric grading, Courant dt computation
- **fdtd**: Higher-order FDTD — 4th-order spatial derivative stencil (c₁=9/8, c₂=-1/24)
- **fdtd**: `CwAccumulator` — running DFT for frequency-domain steady-state extraction
- **fdtd**: `find_eigenmodes` — eigenmode solver via broadband excitation + DFT peak detection

### Changed

- **circuit**: `resistance_series`, `capacitance_parallel`, `inductance_series` now return `Result` and reject empty slices
- **fdtd**: Mur ABC coefficient precomputed in constructor (minor FDTD step performance improvement)
- **fdtd**: Added tracing instrumentation to `new`, `run`
- **field**: Added tracing to `trace_field_line`
- **field**: Added `Default`, `From<[f64;3]>`, `Into<[f64;3]>` for `FieldVector`; exposed `SPEED_OF_LIGHT_SQ` constant
- **field**: Removed redundant `FieldVector::add()` method (use `+` operator)
- **polarization**: Added `Default` for `Complex`, `JonesVector`, `StokesVector`; `From<f64>` and `Div` for `Complex`
- **polarization**: Removed redundant `StokesVector::add()` method (use `+` operator)

### Dependencies

- hisab 0.22.4 → 1.1.0
- impetus 0.23.3 → 1.0.0
- reqwest 0.12 → 0.13
- criterion 0.5 → 0.8

## [0.24.3] — 2026-03-24

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

### Added (V0.5)

- **relativity**: New module — Lorentz factor, beta, four-vectors (Minkowski dot, interval, boost)
- **relativity**: Electromagnetic field tensor F^μν — construction from E/B, field extraction, first/second Lorentz invariants, Hodge dual
- **relativity**: Lorentz transformation of E/B fields — general boost and x-axis special case
- **relativity**: Retarded scalar and vector potentials for point charges
- **relativity**: Liénard-Wiechert fields — electric (velocity + acceleration terms) and magnetic
- **relativity**: Larmor radiation power (non-relativistic and relativistic/Liénard extension)
- **relativity**: Four-potential A^μ

### Added (V0.4)

- **wave**: Snell's law refraction, critical angle, Brewster's angle
- **wave**: Fresnel equations (r_s, r_p), reflectance/transmittance (s, p, normal)
- **wave**: Rectangular waveguide cutoff frequency/wavelength, guide wavelength
- **wave**: Cylindrical waveguide cutoff
- **wave**: Hertzian and half-wave dipole antenna patterns, directivity, radiation resistance
- **wave**: Effective aperture, Friis transmission equation
- **fdtd**: New module — 1D FDTD Yee algorithm solver with Mur absorbing boundaries, dielectric/magnetic material support, soft source injection, energy tracking

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
