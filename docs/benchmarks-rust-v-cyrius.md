# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: rust-old/bench-history.csv, run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 76 twinned, 0 not applicable, 4 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 37.00 ns | 24.3× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 25.00 ns | 20.3× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 79.00 ns | 45.7× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 565.00 ns | 14.2× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 11.00 ns | 22.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 24.00 ns | 24.3× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 107.00 ns | 19.4× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 107.00 ns | 19.4× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 110.00 ns | 89.1× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 33.00 ns | 26.7× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 50.00 ns | 40.5× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.18 µs | 3.7× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 140.00 ns | 113.2× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 85.00 ns | 68.6× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 187.00 ns | 145.4× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 50.00 ns | 9.0× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 48.00 ns | 31.9× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 92.00 ns | 93.5× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 134.00 ns | 258.8× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 21.20 µs | 41349.5× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 289.00 ns | 231.9× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 279.00 ns | 159.8× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 269.00 ns | 155.9× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 288.00 ns | 232.5× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 726.00 ns | 93.3× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 258.00 ns | 208.7× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 463.00 ns | 922.7× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 49.00 ns | 39.6× |
| `fdtd/step_100_cells` | 64.09 ns | `fdtd/step_100_cells` | 711.00 ns | 11.1× |
| `fdtd/step_1000_cells` | 623.03 ns | `fdtd/step_1000_cells` | 5.91 µs | 9.5× |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | `fdtd/run_100_steps_500_cells` | 647.93 µs | 20.1× |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | `fdtd/2d_step_100x100_tm` | 46.37 µs | 4.1× |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | `fdtd/2d_run_50_steps_80x80` | 1.74 ms | 4.2× |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | pending (M6) | — | — |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | `fdtd/2d_step_100x100_tm_lossy` | 49.25 µs | 3.4× |
| `fdtd/3d_step_30x30x30` | 319.98 µs | pending (M6) | — | — |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | pending (M6) | — | — |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | pending (M6) | — | — |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | `fdtd/find_eigenmodes_4096x256` | 10.91 ms | 5.2× |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 68.00 ns | 137.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 74.00 ns | 36.4× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 18.00 ns | 11.8× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 211.00 ns | 160.9× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 368.00 ns | 51.9× |
| `material/polarization` | 1.00 ns | `material/polarization` | 24.00 ns | 24.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 232.00 ns | 233.9× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 253.00 ns | 45.9× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 127.00 ns | 102.3× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 8.00 ns | 6.5× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 148.00 ns | 119.8× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 33.00 ns | 4.1× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 189.00 ns | 32.0× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 717.00 ns | 33.1× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 625.00 ns | 6.0× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 6.04 µs | 18.9× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 303.00 ns | 174.5× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 169.00 ns | 31.0× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 159.00 ns | 160.4× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 164.00 ns | 166.1× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 101.00 ns | 101.6× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 58.00 ns | 10.7× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 2.00 µs | 4.7× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 3.95 µs | 4.6× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 12.86 µs | 3.7× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 77.00 ns | 44.3× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 6.00 ns | 12.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 13.00 ns | 26.3× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 83.00 ns | 46.2× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 30.00 ns | 30.3× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 194.00 ns | 43.5× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 45.00 ns | 6.9× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 744.00 ns | 24.5× |
| `rf/reflection_coefficient` | 1.23 ns | `rf/reflection_coefficient` | 17.00 ns | 13.9× |
| `rf/vswr` | 1.23 ns | `rf/vswr` | 9.00 ns | 7.3× |
| `rf/linear_array_factor_8` | 0.55 ns | `rf/linear_array_factor_8` | 450.00 ns | 814.8× |
| `rf/s_cascade` | 100.25 ns | `rf/s_cascade` | 3.57 µs | 35.6× |
| `rf/input_impedance_complex` | 7.43 ns | `rf/input_impedance_complex` | 253.00 ns | 34.0× |
| `rf/single_stub_match` | 128.46 ns | `rf/single_stub_match` | 2.14 µs | 16.7× |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | `rf/parse_touchstone_s2p_201` | 552.05 µs | 7.0× |
| `rf/write_touchstone_s2p_201` | 213.82 µs | `rf/write_touchstone_s2p_201` | 3.93 ms | 18.4× |
