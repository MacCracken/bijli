# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: Rust 1.2.0's bench-history.csv (scripts/rust_1_2_0.json), run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 80 twinned, 0 not applicable, 0 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 36.00 ns | 23.7× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 25.00 ns | 20.3× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 75.00 ns | 43.4× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 541.00 ns | 13.6× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 10.00 ns | 20.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 23.00 ns | 23.3× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 103.00 ns | 18.6× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 102.00 ns | 18.5× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 108.00 ns | 87.5× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 32.00 ns | 25.9× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 50.00 ns | 40.5× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.15 µs | 3.6× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 135.00 ns | 109.2× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 82.00 ns | 66.1× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 181.00 ns | 140.8× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 48.00 ns | 8.7× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 47.00 ns | 31.2× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 93.00 ns | 94.5× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 128.00 ns | 247.2× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 21.53 µs | 41977.4× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 294.00 ns | 235.9× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 284.00 ns | 162.6× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 268.00 ns | 155.4× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 291.00 ns | 234.9× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 733.00 ns | 94.2× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 264.00 ns | 213.5× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 468.00 ns | 932.6× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 50.00 ns | 40.4× |
| `fdtd/step_100_cells` | 64.09 ns | `fdtd/step_100_cells` | 711.00 ns | 11.1× |
| `fdtd/step_1000_cells` | 623.03 ns | `fdtd/step_1000_cells` | 6.25 µs | 10.0× |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | `fdtd/run_100_steps_500_cells` | 658.23 µs | 20.4× |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | `fdtd/2d_step_100x100_tm` | 108.81 µs | 9.5× |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | `fdtd/2d_run_50_steps_80x80` | 5.47 ms | 13.1× |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | `fdtd/2d_step_100x100_tm_cpml10` | 236.91 µs | 11.6× |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | `fdtd/2d_step_100x100_tm_lossy` | 112.02 µs | 7.8× |
| `fdtd/3d_step_30x30x30` | 319.98 µs | `fdtd/3d_step_30x30x30` | 1.35 ms | 4.2× |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | `fdtd/3d_step_30x30x30_cpml8` | 4.24 ms | 8.6× |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | `fdtd/nf2ff_2d_80x80_180` | 1.52 ms | 7.9× |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | `fdtd/find_eigenmodes_4096x256` | 10.67 ms | 5.0× |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 70.00 ns | 141.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 76.00 ns | 37.4× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 18.00 ns | 11.8× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 213.00 ns | 162.4× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 354.00 ns | 49.9× |
| `material/polarization` | 1.00 ns | `material/polarization` | 23.00 ns | 23.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 224.00 ns | 225.8× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 245.00 ns | 44.5× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 122.00 ns | 98.2× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 7.00 ns | 5.7× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 143.00 ns | 115.7× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 35.00 ns | 4.3× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 193.00 ns | 32.6× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 738.00 ns | 34.1× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 629.00 ns | 6.1× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 6.20 µs | 19.4× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 297.00 ns | 171.0× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 175.00 ns | 32.1× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 163.00 ns | 164.4× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 167.00 ns | 169.1× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 101.00 ns | 101.6× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 57.00 ns | 10.5× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 2.03 µs | 4.8× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 4.05 µs | 4.7× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 12.99 µs | 3.8× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 77.00 ns | 44.3× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 6.00 ns | 12.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 13.00 ns | 26.3× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 84.00 ns | 46.8× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 30.00 ns | 30.3× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 198.00 ns | 44.4× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 50.00 ns | 7.7× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 759.00 ns | 25.0× |
| `rf/reflection_coefficient` | 1.23 ns | `rf/reflection_coefficient` | 17.00 ns | 13.9× |
| `rf/vswr` | 1.23 ns | `rf/vswr` | 9.00 ns | 7.3× |
| `rf/linear_array_factor_8` | 0.55 ns | `rf/linear_array_factor_8` | 463.00 ns | 838.3× |
| `rf/s_cascade` | 100.25 ns | `rf/s_cascade` | 3.65 µs | 36.4× |
| `rf/input_impedance_complex` | 7.43 ns | `rf/input_impedance_complex` | 258.00 ns | 34.7× |
| `rf/single_stub_match` | 128.46 ns | `rf/single_stub_match` | 2.19 µs | 17.0× |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | `rf/parse_touchstone_s2p_201` | 580.29 µs | 7.3× |
| `rf/write_touchstone_s2p_201` | 213.82 µs | `rf/write_touchstone_s2p_201` | 4.12 ms | 19.2× |
