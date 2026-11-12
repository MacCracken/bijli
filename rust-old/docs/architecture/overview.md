# Bijli — Architecture Overview

> Classical electromagnetism for the AGNOS stack. Flat library crate, pure Rust, SI units throughout.

## Module map

| Module | Feature | Depends on | Role |
|--------|---------|------------|------|
| `error` | always | — | `BijliError` (`#[non_exhaustive]`), `Result<T>` alias |
| `bridge` | always | — | Primitive-value conversions for sibling crates (ushma, prakash, tanmatra, dravya); keeps its own CODATA 2018 constants so it builds without `field` |
| `field` | `field` | — | `FieldVector`, physical constants (CODATA 2018), point/dipole/Gauss-law fields, multipoles, images, Green's functions, stress tensor, field-line tracing |
| `maxwell` | `maxwell` | field | Wave speed, impedance, refractive index, skin depth |
| `charge` | `charge` | field | `PointCharge`, Coulomb/Lorentz force, cyclotron motion |
| `material` | `material` | field | Unified `Material` (ε_r, μ_r, σ), Clausius–Mossotti, Curie laws |
| `wave` | `wave` | field, material | Plane waves, Snell, Fresnel (angle and trig-free forms), Schlick, waveguides, antennas, Friis |
| `circuit` | `circuit` | — | Ohm, series/parallel, RC/RL/RLC |
| `relativity` | `relativity` | field | Lorentz transforms, four-vectors, `EmTensor`, retarded-time solver, Liénard–Wiechert |
| `polarization` | `polarization` | hisab `Complex` | Jones, Stokes, Mueller |
| `scattering` | `scattering` | polarization | Mie (Bohren–Huffman), Rayleigh |
| `beam` | `beam` | — | Gaussian beams, ABCD, resonators, HG/LG modes |
| `rf` | `rf` | polarization | Transmission lines, Smith chart, stub matching, 2-port S-parameters, Touchstone v1.1 I/O, arrays |
| `fdtd` | `fdtd` | field, material | 1D/2D/3D Yee solvers, 2D/3D CPML, lossy + Drude/Lorentz/Debye ADE, TFSF, NF2FF, CW accumulation, eigenmodes |
| `compute` | `compute` | fdtd, field | `ComputeBackend` trait, `CpuBackend`, `GpuFdtd2d<B>`; `WgpuBackend` (`compute-wgpu`) is a scaffold |
| `integration::soorat` | `soorat-compat` | (fdtd, charge for some constructors) | Visualization DTOs for the soorat renderer |
| `ai` | `ai` | reqwest, tokio, serde_json | Daimon/hoosh agent registration client |
| `logging` | `logging` | tracing-subscriber | `BIJLI_LOG` env-filter setup |

Every feature compiles on its own (`make features`, CI job *Feature Isolation*): inter-module
imports must be mirrored as feature dependencies in `Cargo.toml`.

## Data flow

```
Material ──► maxwell / wave (Fresnel, impedance) ──► kiran (rendering), prakash (optics)
   │
   └──► fdtd (set_material: ε, μ, σ) ──► Fdtd1d / Fdtd2d (+CPML) / Fdtd3d (+CPML)
                                              │
              sources (gaussian, ricker, TFSF)┘── CwAccumulator ──► nf2ff_2d_phasor ──► RadiationPattern (soorat)
PointCharge ──► field (E, φ, B) ──► charge (forces) ──► impetus (integration, outside bijli)
trajectory fn ──► relativity::solve_retarded_time ──► lienard_wiechert_fields
```

## Conventions

- **Units**: SI everywhere; angles in radians unless a name says `_deg`; RF impedances in ohms, frequencies in Hz.
- **Time convention**: `polarization` and `beam` use e^{−iωt}; `scattering::mie` uses e^{+iωt}, so absorbing media have Im(m) < 0 (`mie` rejects Im(m) > 0). Unifying these is a 2.0 item.
- **Relativity**: `FourVector::dot` uses the (−,+,+,+) metric; `EmTensor` stores F^μν in the (+,−,−,−) Jackson/Griffiths layout. Use `EmTensor::four_force` rather than contracting by hand.
- **Quadrupoles**: traceless Qᵢⱼ = Σq(3xᵢxⱼ − r²δᵢⱼ) (Jackson 4.10).
- **FDTD**: Yee staggering with E_z on integer nodes in 2D TM; default walls are PEC (TM, 3D) or PMC (TE) — enable CPML for open boundaries. dt defaults to half the Courant limit.
- **Touchstone**: v1.1; 2-port data ordered `N11 N21 N12 N22`, ≥ 3 ports row-major with line wrapping.

## Error-handling policy

- Library code never panics on user input: no `unwrap`/`expect`/`panic!`, bounds-checked indexing, checked size arithmetic for untrusted inputs (Touchstone port counts, grid sizes).
- Fallible operations return `Result<_, BijliError>`; guards are NaN-safe (`x.is_nan() || x <= 0.0`, `is_finite`).
- Legacy constructors that return `Self` sanitize and `tracing::warn!`; their `try_*` siblings return `Err`.
- Functions whose correct answer does not fit their signature are deprecated in favor of a correctly-typed replacement and return `Err` rather than a wrong value.

## Dependency stack

- **hisab** (AGNOS numerics) — complex arithmetic (`hisab::num::Complex`)
- **serde**, **thiserror**, **tracing** — always on
- Optional: **reqwest** + **tokio** + **serde_json** (`ai`), **tracing-subscriber** (`logging`), **wgpu** (`compute-wgpu`)

## Consumers

- **kiran** (game engine) — Fresnel/Schlick on hot paths, scattering, beams
- **joshua** (simulation) — FDTD, RF, antennas
- **prakash** (optics) — polarization, beam propagation, material optics
