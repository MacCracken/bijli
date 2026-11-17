# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: rust-old/bench-history.csv, run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 78 twinned, 0 not applicable, 2 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 36.00 ns | 23.7× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 25.00 ns | 20.3× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 76.00 ns | 44.0× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 546.00 ns | 13.7× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 11.00 ns | 22.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 23.00 ns | 23.3× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 106.00 ns | 19.2× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 105.00 ns | 19.0× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 109.00 ns | 88.3× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 33.00 ns | 26.7× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 50.00 ns | 40.5× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.17 µs | 3.7× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 138.00 ns | 111.6× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 84.00 ns | 67.8× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 183.00 ns | 142.3× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 49.00 ns | 8.8× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 47.00 ns | 31.2× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 90.00 ns | 91.4× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 133.00 ns | 256.9× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 20.82 µs | 40598.7× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 283.00 ns | 227.1× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 276.00 ns | 158.0× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 262.00 ns | 151.9× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 282.00 ns | 227.7× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 714.00 ns | 91.8× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 253.00 ns | 204.6× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 454.00 ns | 904.7× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 49.00 ns | 39.6× |
| `fdtd/step_100_cells` | 64.09 ns | `fdtd/step_100_cells` | 711.00 ns | 11.1× |
| `fdtd/step_1000_cells` | 623.03 ns | `fdtd/step_1000_cells` | 5.91 µs | 9.5× |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | `fdtd/run_100_steps_500_cells` | 647.93 µs | 20.1× |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | `fdtd/2d_step_100x100_tm` | 46.37 µs | 4.1× |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | `fdtd/2d_run_50_steps_80x80` | 1.74 ms | 4.2× |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | `fdtd/2d_step_100x100_tm_cpml10` | 114.62 µs | 5.6× |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | `fdtd/2d_step_100x100_tm_lossy` | 49.25 µs | 3.4× |
| `fdtd/3d_step_30x30x30` | 319.98 µs | pending (M6) | — | — |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | pending (M6) | — | — |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | `fdtd/nf2ff_2d_80x80_180` | 1.57 ms | 8.1× |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | `fdtd/find_eigenmodes_4096x256` | 10.68 ms | 5.0× |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 67.00 ns | 135.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 72.00 ns | 35.4× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 17.00 ns | 11.1× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 205.00 ns | 156.3× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 354.00 ns | 49.9× |
| `material/polarization` | 1.00 ns | `material/polarization` | 24.00 ns | 24.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 230.00 ns | 231.9× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 248.00 ns | 45.0× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 125.00 ns | 100.7× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 8.00 ns | 6.5× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 144.00 ns | 116.5× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 33.00 ns | 4.1× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 186.00 ns | 31.5× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 701.00 ns | 32.4× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 604.00 ns | 5.8× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 5.97 µs | 18.7× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 295.00 ns | 169.9× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 163.00 ns | 29.9× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 158.00 ns | 159.4× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 162.00 ns | 164.1× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 101.00 ns | 101.6× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 56.00 ns | 10.3× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 1.97 µs | 4.7× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 3.93 µs | 4.6× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 12.51 µs | 3.6× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 75.00 ns | 43.2× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 6.00 ns | 12.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 12.00 ns | 24.3× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 83.00 ns | 46.2× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 30.00 ns | 30.3× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 193.00 ns | 43.3× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 45.00 ns | 6.9× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 740.00 ns | 24.4× |
| `rf/reflection_coefficient` | 1.23 ns | `rf/reflection_coefficient` | 17.00 ns | 13.9× |
| `rf/vswr` | 1.23 ns | `rf/vswr` | 9.00 ns | 7.3× |
| `rf/linear_array_factor_8` | 0.55 ns | `rf/linear_array_factor_8` | 446.00 ns | 807.5× |
| `rf/s_cascade` | 100.25 ns | `rf/s_cascade` | 3.56 µs | 35.5× |
| `rf/input_impedance_complex` | 7.43 ns | `rf/input_impedance_complex` | 251.00 ns | 33.8× |
| `rf/single_stub_match` | 128.46 ns | `rf/single_stub_match` | 2.14 µs | 16.7× |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | `rf/parse_touchstone_s2p_201` | 546.63 µs | 6.9× |
| `rf/write_touchstone_s2p_201` | 213.82 µs | `rf/write_touchstone_s2p_201` | 3.93 ms | 18.4× |
