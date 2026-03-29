# Bijli — Electromagnetism Simulation

> Hindi: बिजली (bijli) — electricity, lightning

Electromagnetism simulation for the AGNOS ecosystem. Fields, Maxwell's equations, charge dynamics, EM wave propagation, optics, scattering, FDTD solvers, RF/microwave, and GPU compute.

Built on [hisab](https://crates.io/crates/hisab) for PDE solvers and [impetus](https://crates.io/crates/impetus) for particle interactions.

## Quick Start

```rust
use bijli::field;
use bijli::charge::PointCharge;
use bijli::material::Material;

// Electric field from a point charge
let q = PointCharge::proton([0.0; 3]);
let e = q.electric_field_at([1.0, 0.0, 0.0]).unwrap();

// Reflectance at air-glass interface
let air = Material::vacuum();
let glass = Material::dielectric(2.25);
let r = air.reflectance_at(&glass); // 4%

// 2D FDTD simulation
use bijli::fdtd::{Fdtd2d, Mode2d};
let mut sim = Fdtd2d::new(200, 200, 1e-4, Mode2d::Tm).unwrap();
sim.set_material(80, 80, 120, 120, &glass).unwrap();
sim.run(500, Some((100, 100)), 3e9, 1.0);
```

## Modules

| Module | Default | Description |
|--------|---------|-------------|
| `field` | Yes | Electric/magnetic fields, potentials, dipoles, multipoles, Green's functions, stress tensor |
| `maxwell` | Yes | Maxwell's equations, wave speed, impedance, skin depth |
| `charge` | Yes | Point charges, Coulomb/Lorentz force, cyclotron frequency |
| `wave` | Yes | EM waves, Fresnel equations, waveguides, antenna patterns |
| `circuit` | Yes | Ohm's law, RC/RL/RLC circuits, series/parallel combinations |
| `material` | Yes | Dielectric/magnetic materials, unified `Material` struct |
| `fdtd` | Yes | 1D/2D/3D FDTD solvers, CPML, dispersive materials, NF2FF, eigenmode |
| `relativity` | Yes | EM tensor, Lorentz transforms, Liénard-Wiechert fields, retarded time |
| `polarization` | Yes | Jones/Stokes/Mueller formalism, optical elements |
| `scattering` | Yes | Mie and Rayleigh scattering |
| `beam` | Yes | Gaussian beams, ABCD matrices, resonator stability, HG/LG modes |
| `rf` | Yes | Transmission lines, S-parameters, Smith chart, Touchstone I/O, arrays |
| `compute` | No | GPU compute backend trait (CPU/wgpu) |
| `ai` | No | Daimon/hoosh AI integration |
| `logging` | No | Structured logging via `BIJLI_LOG` |

## Consumers

- **kiran** — EM effects in game physics, sky rendering, volumetric scattering
- **joshua** — Full EM simulation, RF circuits, antenna design, 3D FDTD
- **prakash** — Optical systems, polarization, beam propagation

## Building

```bash
cargo build                    # default features
cargo build --all-features     # everything including GPU backends
cargo test --all-features      # full test suite (589+ tests)
make check                     # fmt + clippy + test + audit
make bench                     # criterion benchmarks with history
```

## License

GPL-3.0-only — see [LICENSE](LICENSE).
