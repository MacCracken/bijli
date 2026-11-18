# Benchmarks: Rust 1.2.0 vs Cyrius

> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius
> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan
> §0.1): a slow row is a finding to file with its numbers.

**Rust column**: Rust 1.2.0's bench-history.csv (scripts/rust_1_2_0.json), run `2026-11-12T03:15:21Z` at `c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own noise is read as parity (port plan §5 P4).

**80 Rust benches:** 80 twinned, 0 not applicable, 0 pending.

| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |
|---|---:|---|---:|---:|
| `field/electric_field_point` | 1.52 ns | `field/electric_field_point` | 35.00 ns | 23.0× |
| `field/electric_potential_point` | 1.23 ns | `field/electric_potential_point` | 24.00 ns | 19.5× |
| `field/magnetic_field_moving` | 1.73 ns | `field/magnetic_field_moving` | 74.00 ns | 42.8× |
| `field/superposition_10` | 39.79 ns | `field/superposition_10` | 536.00 ns | 13.5× |
| `field/vector_magnitude` | 0.50 ns | `field/vector_magnitude` | 10.00 ns | 20.0× |
| `field/vector_cross` | 0.99 ns | `field/vector_cross` | 22.00 ns | 22.3× |
| `field/electric_dipole_field` | 5.53 ns | `field/electric_dipole_field` | 100.00 ns | 18.1× |
| `field/magnetic_dipole_field` | 5.52 ns | `field/magnetic_dipole_field` | 102.00 ns | 18.5× |
| `field/charged_sphere` | 1.23 ns | `field/charged_sphere` | 104.00 ns | 84.2× |
| `field/ring_axis` | 1.23 ns | `field/ring_axis` | 31.00 ns | 25.1× |
| `field/disk_axis` | 1.23 ns | `field/disk_axis` | 48.00 ns | 38.9× |
| `field/trace_field_line_10` | 317.87 ns | `field/trace_field_line_10` | 1.12 µs | 3.5× |
| `maxwell/wave_speed` | 1.24 ns | `maxwell/wave_speed` | 133.00 ns | 107.6× |
| `maxwell/impedance` | 1.24 ns | `maxwell/impedance` | 80.00 ns | 64.5× |
| `maxwell/skin_depth` | 1.29 ns | `maxwell/skin_depth` | 177.00 ns | 137.7× |
| `charge/coulomb_force` | 5.54 ns | `charge/coulomb_force` | 48.00 ns | 8.7× |
| `charge/lorentz_force` | 1.51 ns | `charge/lorentz_force` | 46.00 ns | 30.5× |
| `wave/poynting_vector` | 0.98 ns | `wave/poynting_vector` | 86.00 ns | 87.4× |
| `wave/plane_wave_intensity` | 0.52 ns | `wave/plane_wave_intensity` | 127.00 ns | 245.3× |
| `wave/plane_wave_e` | 0.51 ns | `wave/plane_wave_e` | 20.05 µs | 39095.2× |
| `wave/fresnel_rs` | 1.25 ns | `wave/fresnel_rs` | 275.00 ns | 220.6× |
| `wave/snell_refraction` | 1.75 ns | `wave/snell_refraction` | 267.00 ns | 152.9× |
| `wave/fresnel_rs_direct` | 1.73 ns | `wave/fresnel_rs_direct` | 249.00 ns | 144.3× |
| `wave/reflectance_unpolarized` | 1.24 ns | `wave/reflectance_unpolarized` | 271.00 ns | 218.8× |
| `wave/reflectance_at_interface` | 7.78 ns | `wave/reflectance_at_interface` | 690.00 ns | 88.7× |
| `wave/schlick_reflectance` | 1.24 ns | `wave/schlick_reflectance` | 242.00 ns | 195.7× |
| `wave/half_wave_dipole_pattern` | 0.50 ns | `wave/half_wave_dipole_pattern` | 439.00 ns | 874.9× |
| `wave/rect_waveguide_cutoff` | 1.24 ns | `wave/rect_waveguide_cutoff` | 46.00 ns | 37.2× |
| `fdtd/step_100_cells` | 64.09 ns | `fdtd/step_100_cells` | 732.00 ns | 11.4× |
| `fdtd/step_1000_cells` | 623.03 ns | `fdtd/step_1000_cells` | 5.88 µs | 9.4× |
| `fdtd/run_100_steps_500_cells` | 32.23 µs | `fdtd/run_100_steps_500_cells` | 633.38 µs | 19.7× |
| `fdtd/2d_step_100x100_tm` | 11.42 µs | `fdtd/2d_step_100x100_tm` | 111.62 µs | 9.8× |
| `fdtd/2d_run_50_steps_80x80` | 416.78 µs | `fdtd/2d_run_50_steps_80x80` | 5.32 ms | 12.8× |
| `fdtd/2d_step_100x100_tm_cpml10` | 20.47 µs | `fdtd/2d_step_100x100_tm_cpml10` | 236.59 µs | 11.6× |
| `fdtd/2d_step_100x100_tm_lossy` | 14.31 µs | `fdtd/2d_step_100x100_tm_lossy` | 113.61 µs | 7.9× |
| `fdtd/3d_step_30x30x30` | 319.98 µs | `fdtd/3d_step_30x30x30` | 1.34 ms | 4.2× |
| `fdtd/3d_step_30x30x30_cpml8` | 494.96 µs | `fdtd/3d_step_30x30x30_cpml8` | 4.14 ms | 8.4× |
| `fdtd/nf2ff_2d_80x80_180` | 193.01 µs | `fdtd/nf2ff_2d_80x80_180` | 1.50 ms | 7.8× |
| `fdtd/find_eigenmodes_4096x256` | 2.12 ms | `fdtd/find_eigenmodes_4096x256` | 10.44 ms | 4.9× |
| `circuit/rc_charging` | 0.50 ns | `circuit/rc_charging` | 65.00 ns | 131.0× |
| `circuit/rl_current_rise` | 2.03 ns | `circuit/rl_current_rise` | 70.00 ns | 34.5× |
| `circuit/rlc_impedance` | 1.53 ns | `circuit/rlc_impedance` | 16.00 ns | 10.5× |
| `circuit/resonant_frequency` | 1.31 ns | `circuit/resonant_frequency` | 199.00 ns | 151.8× |
| `circuit/resistance_parallel` | 7.09 ns | `circuit/resistance_parallel` | 344.00 ns | 48.5× |
| `material/polarization` | 1.00 ns | `material/polarization` | 22.00 ns | 22.1× |
| `material/displacement_field` | 0.99 ns | `material/displacement_field` | 218.00 ns | 219.8× |
| `material/h_field_from_b` | 5.51 ns | `material/h_field_from_b` | 239.00 ns | 43.4× |
| `material/clausius_mossotti` | 1.24 ns | `material/clausius_mossotti` | 120.00 ns | 96.6× |
| `material/curie_weiss` | 1.24 ns | `material/curie_weiss` | 7.00 ns | 5.7× |
| `relativity/lorentz_factor` | 1.24 ns | `relativity/lorentz_factor` | 144.00 ns | 116.5× |
| `relativity/em_tensor_from_fields` | 8.11 ns | `relativity/em_tensor_from_fields` | 31.00 ns | 3.8× |
| `relativity/lorentz_transform_x` | 5.91 ns | `relativity/lorentz_transform_x` | 183.00 ns | 30.9× |
| `relativity/lienard_wiechert_e` | 21.64 ns | `relativity/lienard_wiechert_e` | 678.00 ns | 31.3× |
| `relativity/solve_retarded_time_uniform_t1s` | 103.48 ns | `relativity/solve_retarded_time_uniform_t1s` | 601.00 ns | 5.8× |
| `relativity/lienard_wiechert_fields_orbit_b099` | 319.92 ns | `relativity/lienard_wiechert_fields_orbit_b099` | 5.83 µs | 18.2× |
| `relativity/lienard_wiechert_scalar_potential` | 1.74 ns | `relativity/lienard_wiechert_scalar_potential` | 288.00 ns | 165.8× |
| `relativity/four_vector_boost` | 5.45 ns | `relativity/four_vector_boost` | 164.00 ns | 30.1× |
| `polarization/jones_apply_polarizer` | 0.99 ns | `polarization/jones_apply_polarizer` | 152.00 ns | 153.3× |
| `polarization/jones_qwp_apply` | 0.99 ns | `polarization/jones_qwp_apply` | 157.00 ns | 159.0× |
| `polarization/mueller_apply` | 0.99 ns | `polarization/mueller_apply` | 96.00 ns | 96.6× |
| `polarization/stokes_from_jones` | 5.44 ns | `polarization/stokes_from_jones` | 53.00 ns | 9.8× |
| `scattering/mie_x1_real` | 423.87 ns | `scattering/mie_x1_real` | 1.90 µs | 4.5× |
| `scattering/mie_x10_real` | 859.83 ns | `scattering/mie_x10_real` | 3.82 µs | 4.4× |
| `scattering/mie_x50_m3_5` | 3.44 µs | `scattering/mie_x50_m3_5` | 12.15 µs | 3.5× |
| `scattering/rayleigh_cross_section` | 1.74 ns | `scattering/rayleigh_cross_section` | 73.00 ns | 42.0× |
| `scattering/rayleigh_phase_function` | 0.50 ns | `scattering/rayleigh_phase_function` | 6.00 ns | 12.0× |
| `beam/spot_size` | 0.49 ns | `beam/spot_size` | 12.00 ns | 24.3× |
| `beam/abcd_propagate_beam` | 1.79 ns | `beam/abcd_propagate_beam` | 80.00 ns | 44.6× |
| `beam/abcd_compose` | 0.99 ns | `beam/abcd_compose` | 29.00 ns | 29.3× |
| `beam/hermite_gaussian_hg33` | 4.46 ns | `beam/hermite_gaussian_hg33` | 186.00 ns | 41.7× |
| `beam/resonator_stability` | 6.48 ns | `beam/resonator_stability` | 43.00 ns | 6.6× |
| `beam/hermite_gaussian_field_hg33` | 30.33 ns | `beam/hermite_gaussian_field_hg33` | 712.00 ns | 23.5× |
| `rf/reflection_coefficient` | 1.23 ns | `rf/reflection_coefficient` | 16.00 ns | 13.0× |
| `rf/vswr` | 1.23 ns | `rf/vswr` | 9.00 ns | 7.3× |
| `rf/linear_array_factor_8` | 0.55 ns | `rf/linear_array_factor_8` | 431.00 ns | 780.4× |
| `rf/s_cascade` | 100.25 ns | `rf/s_cascade` | 3.48 µs | 34.7× |
| `rf/input_impedance_complex` | 7.43 ns | `rf/input_impedance_complex` | 243.00 ns | 32.7× |
| `rf/single_stub_match` | 128.46 ns | `rf/single_stub_match` | 2.04 µs | 15.9× |
| `rf/parse_touchstone_s2p_201` | 79.42 µs | `rf/parse_touchstone_s2p_201` | 529.88 µs | 6.7× |
| `rf/write_touchstone_s2p_201` | 213.82 µs | `rf/write_touchstone_s2p_201` | 3.75 ms | 17.5× |
