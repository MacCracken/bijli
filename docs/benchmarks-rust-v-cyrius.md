# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: Rust 1.2.0's bench-history.csv (scripts/rust_1_2_0.json), run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 80 twinned, 0 not applicable, 0 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 36.00 ns | 23.7× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 24.00 ns | 19.5× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 74.00 ns | 42.8× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 526.00 ns | 13.2× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 10.00 ns | 20.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 22.00 ns | 22.3× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 100.00 ns | 18.1× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 101.00 ns | 18.3× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 103.00 ns | 83.4× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 32.00 ns | 25.9× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 48.00 ns | 38.9× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.15 µs | 3.6× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 132.00 ns | 106.8× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 80.00 ns | 64.5× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 176.00 ns | 136.9× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 47.00 ns | 8.5× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 45.00 ns | 29.9× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 85.00 ns | 86.4× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 126.00 ns | 243.3× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 19.80 µs | 38601.8× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 268.00 ns | 215.0× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 266.00 ns | 152.3× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 254.00 ns | 147.2× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 276.00 ns | 222.8× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 688.00 ns | 88.5× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 252.00 ns | 203.8× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 430.00 ns | 856.9× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 46.00 ns | 37.2× |
| `fdtd/step_100_cells` | 64.09 ns | `fdtd/step_100_cells` | 689.00 ns | 10.8× |
| `fdtd/step_1000_cells` | 623.03 ns | `fdtd/step_1000_cells` | 5.69 µs | 9.1× |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | `fdtd/run_100_steps_500_cells` | 624.82 µs | 19.4× |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | `fdtd/2d_step_100x100_tm` | 105.66 µs | 9.3× |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | `fdtd/2d_run_50_steps_80x80` | 5.30 ms | 12.7× |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | `fdtd/2d_step_100x100_tm_cpml10` | 230.53 µs | 11.3× |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | `fdtd/2d_step_100x100_tm_lossy` | 108.83 µs | 7.6× |
| `fdtd/3d_step_30x30x30` | 319.98 µs | `fdtd/3d_step_30x30x30` | 1.30 ms | 4.1× |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | `fdtd/3d_step_30x30x30_cpml8` | 4.07 ms | 8.2× |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | `fdtd/nf2ff_2d_80x80_180` | 1.49 ms | 7.7× |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | `fdtd/find_eigenmodes_4096x256` | 10.35 ms | 4.9× |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 64.00 ns | 129.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 70.00 ns | 34.5× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 16.00 ns | 10.5× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 198.00 ns | 151.0× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 343.00 ns | 48.4× |
| `material/polarization` | 1.00 ns | `material/polarization` | 23.00 ns | 23.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 218.00 ns | 219.8× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 238.00 ns | 43.2× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 120.00 ns | 96.6× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 7.00 ns | 5.7× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 140.00 ns | 113.3× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 33.00 ns | 4.1× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 180.00 ns | 30.4× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 690.00 ns | 31.9× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 588.00 ns | 5.7× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 5.76 µs | 18.0× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 280.00 ns | 161.2× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 160.00 ns | 29.3× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 152.00 ns | 153.3× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 157.00 ns | 159.0× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 95.00 ns | 95.6× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 53.00 ns | 9.8× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 1.91 µs | 4.5× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 3.78 µs | 4.4× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 12.02 µs | 3.5× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 76.00 ns | 43.8× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 6.00 ns | 12.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 12.00 ns | 24.3× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 80.00 ns | 44.6× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 28.00 ns | 28.3× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 186.00 ns | 41.7× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 43.00 ns | 6.6× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 715.00 ns | 23.6× |
| `rf/reflection_coefficient` | 1.23 ns | `rf/reflection_coefficient` | 16.00 ns | 13.0× |
| `rf/vswr` | 1.23 ns | `rf/vswr` | 8.00 ns | 6.5× |
| `rf/linear_array_factor_8` | 0.55 ns | `rf/linear_array_factor_8` | 423.00 ns | 765.9× |
| `rf/s_cascade` | 100.25 ns | `rf/s_cascade` | 3.46 µs | 34.5× |
| `rf/input_impedance_complex` | 7.43 ns | `rf/input_impedance_complex` | 240.00 ns | 32.3× |
| `rf/single_stub_match` | 128.46 ns | `rf/single_stub_match` | 2.04 µs | 15.8× |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | `rf/parse_touchstone_s2p_201` | 522.50 µs | 6.6× |
| `rf/write_touchstone_s2p_201` | 213.82 µs | `rf/write_touchstone_s2p_201` | 3.72 ms | 17.4× |
