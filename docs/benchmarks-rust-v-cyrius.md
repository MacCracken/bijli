# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: rust-old/bench-history.csv, run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 61 twinned, 0 not applicable, 19 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 30.00 ns | 19.7× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 18.00 ns | 14.6× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 41.00 ns | 23.7× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 522.00 ns | 13.1× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 9.00 ns | 18.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 18.00 ns | 18.2× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 82.00 ns | 14.8× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 82.00 ns | 14.9× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 108.00 ns | 87.5× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 25.00 ns | 20.2× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 38.00 ns | 30.8× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.06 µs | 3.3× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 141.00 ns | 114.0× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 84.00 ns | 67.8× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 186.00 ns | 144.7× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 42.00 ns | 7.6× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 59.00 ns | 39.2× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 95.00 ns | 96.9× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 148.00 ns | 284.6× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 20.99 µs | 40928.2× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 304.00 ns | 243.2× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 284.00 ns | 162.6× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 277.00 ns | 160.1× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 301.00 ns | 242.7× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 755.00 ns | 97.0× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 265.00 ns | 213.7× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 483.00 ns | 962.5× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 51.00 ns | 41.1× |
| `fdtd/step_100_cells` | 64.09 ns | pending (M5) | — | — |
| `fdtd/step_1000_cells` | 623.03 ns | pending (M5) | — | — |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | pending (M5) | — | — |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | pending (M6) | — | — |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | pending (M6) | — | — |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | pending (M6) | — | — |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | pending (M6) | — | — |
| `fdtd/3d_step_30x30x30` | 319.98 µs | pending (M6) | — | — |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | pending (M6) | — | — |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | pending (M5) | — | — |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | pending (M5) | — | — |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 67.00 ns | 134.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 72.00 ns | 35.5× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 17.00 ns | 11.1× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 217.00 ns | 165.6× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 384.00 ns | 54.2× |
| `material/polarization` | 1.00 ns | `material/polarization` | 25.00 ns | 25.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 231.00 ns | 233.3× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 256.00 ns | 46.5× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 128.00 ns | 103.2× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 8.00 ns | 6.5× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 153.00 ns | 123.4× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 36.00 ns | 4.4× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 186.00 ns | 31.5× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 710.00 ns | 32.8× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 612.00 ns | 5.9× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 6.05 µs | 18.9× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 292.00 ns | 167.8× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 161.00 ns | 29.5× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 213.00 ns | 215.2× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 216.00 ns | 218.2× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 86.00 ns | 86.9× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 117.00 ns | 21.5× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 1.81 µs | 4.3× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 3.62 µs | 4.2× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 11.74 µs | 3.4× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 129.00 ns | 74.1× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 4.00 ns | 8.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 10.00 ns | 20.4× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 19.00 ns | 10.6× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 20.00 ns | 20.2× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 160.00 ns | 35.9× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 39.00 ns | 6.0× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 669.00 ns | 22.1× |
| `rf/reflection_coefficient` | 1.23 ns | pending (M4) | — | — |
| `rf/vswr` | 1.23 ns | pending (M4) | — | — |
| `rf/linear_array_factor_8` | 0.55 ns | pending (M4) | — | — |
| `rf/s_cascade` | 100.25 ns | pending (M4) | — | — |
| `rf/input_impedance_complex` | 7.43 ns | pending (M4) | — | — |
| `rf/single_stub_match` | 128.46 ns | pending (M4) | — | — |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | pending (M4) | — | — |
| `rf/write_touchstone_s2p_201` | 213.82 µs | pending (M4) | — | — |
