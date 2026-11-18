# Bijli — Electromagnetism Simulation

> Hindi: बिजली (bijli) — electricity, lightning

Electromagnetism simulation for the AGNOS ecosystem: fields, Maxwell's equations, charge dynamics, EM wave
propagation, optics, scattering, FDTD solvers, RF/microwave networks and a compute backend. Written in
[Cyrius](https://github.com/MacCracken/cyrius) on its stdlib alone.

bijli 2.0 is the Cyrius port of the Rust crate 1.2.0: the same physics, results, errors and branches, checked
against a golden corpus of 770,000 rows, with Rust's inconsistencies fixed and recorded
([ADR 0002](docs/adr/0002-port-api-divergences.md)). Every result is bit-identical on x86_64 (ELF and PE) and
aarch64, every transcendental is correctly rounded, and one time convention holds everywhere: e^{+jωt}, IEEE
handedness. Coming from 1.2? See the [migration guide](docs/guides/migration-2.0.md).

## Quick Start

```cyrius
# The electric field of a proton 1 m away.
var origin: f64[3];
origin[0] = 0.0;
origin[1] = 0.0;
origin[2] = 0.0;
var proton: BjPointCharge = BjPointCharge_proton(&origin);
var at: f64[3];
at[0] = 1.0;
at[1] = 0.0;
at[2] = 0.0;
var e: BjFieldVector;
var te, ve = proton.electric_field_at(&at, &e);
if (is_err_result(te) == 1) { return 1; }

# Reflectance at an air-glass interface: about 4 %.
var air: BjMaterial = BjMaterial_vacuum();
var glass: BjMaterial = BjMaterial_dielectric(2.25);
var tr, vr = air.try_reflectance_at(glass);
var r: f64 = vr;
if (is_err_result(tr) == 1 || r < 0.03 || r > 0.05) { return 2; }

# A 2D TM FDTD run: a glass block and a 3 GHz source.
var sim: BjFdtd2d;
var ts, vs = BjFdtd2d_new(200, 200, 0.0001, BJ_MODE2D_TM, &sim);
if (is_err_result(ts) == 1) { return 3; }
var tm, vm = sim.set_material(80, 80, 120, 120, &glass);
sim.run(500, 100, 100, 3000000000.0, 1.0);
BjFdtd2d_free(&sim);
```

[Getting started](docs/guides/getting-started.md) shows a whole program and how to take bijli as a dependency.

## Modules

| Module | Bundle | Description |
|--------|--------|-------------|
| `field` | default | Electric/magnetic fields, potentials, dipoles, multipoles, Green's functions, stress tensor |
| `maxwell` | default | Maxwell's equations, wave speed, impedance, skin depth |
| `charge` | default | Point charges, Coulomb/Lorentz force, cyclotron frequency |
| `wave` | default | EM waves, Fresnel equations, waveguides, antenna patterns |
| `circuit` | default | Ohm's law, RC/RL/RLC circuits, series/parallel combinations |
| `material` | default | Dielectric/magnetic materials, the unified `BjMaterial` |
| `fdtd` | default | 1D/2D/3D FDTD solvers, CPML (2D/3D), lossy/dispersive materials, TF/SF, NF2FF, eigenmodes |
| `relativity` | default | EM tensor, Lorentz transforms, Liénard–Wiechert fields, retarded time |
| `polarization` | default | Jones/Stokes/Mueller formalism, optical elements |
| `scattering` | default | Mie and Rayleigh scattering |
| `beam` | default | Gaussian beams, ABCD matrices, resonator stability, HG/LG modes |
| `rf` | default | Transmission lines, S-parameters, Smith chart, Touchstone I/O, arrays |
| `bridge` | default | Cross-crate bridges: primitive values to and from ushma, prakash, tanmatra and dravya |
| `soorat` | default | Visualization data for soorat: FDTD field slices, field lines, charges, radiation patterns, vector fields |
| `compute` | default | The compute-backend trait, the CPU reference backend, `GpuFdtd2d` (the wgpu backend is a scaffold) |
| `json` | default | JSON for every serializable type, as serde_json wrote it |
| `logging` | default | Structured logging on sakshi, off until `bj_log_init` (`BIJLI_LOG`) |
| `ai` | ai | The daimon/hoosh client: `dist/bijli-ai.cyr`, with sandhi's HTTP/TLS |

[`docs/architecture/README.md`](docs/architecture/README.md) maps the modules to their source files.

## Consumers

- **kiran** — EM effects in game physics, sky rendering, volumetric scattering
- **joshua** — Full EM simulation, RF circuits, antenna design, 3D FDTD
- **prakash** — Optical systems, polarization, beam propagation (through the bridges)
- **tanmatra** — Atomic physics (through the bridges)
- **soorat** — Renders the soorat visualization structs

## Building

```bash
cyrius deps                                   # resolve the stdlib (vendors lib/)
cyrius build src/main.cyr build/bijli         # compile
cyrius test                                   # the unit tests and the golden parity harness
cyrius bench                                  # the benchmarks
cyrius distlib --all                          # dist/bijli.cyr and dist/bijli-ai.cyr
python3 -I -B scripts/oracle/validate.py .    # the reference model checks every golden row
```

There is no Rust in the tree ([ADR 0003](docs/adr/0003-retire-the-rust-tree.md)): Rust bijli 1.2.0 is the
specification, read from the `1.2.0` tag.

## License

GPL-3.0-only — see [LICENSE](LICENSE).
