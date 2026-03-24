# Bijli — Electromagnetism Simulation

> Hindi: बिजली (bijli) — electricity, lightning

Electromagnetism simulation for the AGNOS ecosystem. Electric and magnetic fields, Maxwell's equations, charge dynamics, and EM wave propagation.

Built on [hisab](https://crates.io/crates/hisab) for PDE solvers and [impetus](https://crates.io/crates/impetus) for particle interactions.

## Quick Start

```rust
use bijli::field;
use bijli::charge::{self, PointCharge};
use bijli::maxwell;

// Electric field from a point charge
let e = field::electric_field_point_charge(1e-6, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();

// Coulomb force
let q1 = PointCharge::proton([0.0; 3]);
let q2 = PointCharge::electron([5.3e-11, 0.0, 0.0]);
let f = charge::coulomb_force(&q1, &q2).unwrap();

// Speed of light from Maxwell's equations
let c = maxwell::wave_speed(field::EPSILON_0, field::MU_0).unwrap();
```

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `field` | Yes | Electric/magnetic vector fields, potentials |
| `maxwell` | Yes | Maxwell's equations, wave speed, impedance |
| `charge` | Yes | Point charges, Coulomb, Lorentz force |
| `wave` | Yes | EM wave propagation, Poynting vector |
| `ai` | No | Daimon/hoosh AI integration |
| `logging` | No | Structured logging via `BIJLI_LOG` |

## Architecture

```
bijli
├── field     — FieldVector, point charge fields, potentials, energy density
├── maxwell   — Gauss's law, wave speed, impedance, skin depth
├── charge    — PointCharge, Coulomb force, Lorentz force, cyclotron
├── wave      — Poynting vector, plane waves, radiation pressure
├── error     — BijliError
├── ai        — DaimonClient (optional)
└── logging   — Structured logging (optional)
```

## Consumers

- **kiran** — EM effects in game physics
- **joshua** — EM simulation mode
- **prakash** — optics builds on EM foundations

## Building

```bash
cargo build                    # default features
cargo build --all-features     # everything
cargo test --all-features      # full test suite
make check                     # fmt + clippy + test + audit
make bench                     # criterion benchmarks with history
```

## License

GPL-3.0 — see [LICENSE](LICENSE).
