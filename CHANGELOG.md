# Changelog

All notable changes to bijli will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

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
