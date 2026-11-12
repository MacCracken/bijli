# Benchmarks

Three-point tracking: **baseline** (first run) / **previous** / **latest**

| Point | Date | Commit |
|-------|------|--------|
| Baseline | 2026-03-25T12:34:25Z | `2df8cdf` |
| Previous | 2026-11-10T06:23:08Z | `678f256` |
| Latest | 2026-11-12T03:15:21Z | `c95fd30` |

## field

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `electric_field_point` | 2.07 ns | 1.70 ns | 1.52 ns |
| `electric_potential_point` | 1.55 ns | 1.31 ns | 1.23 ns |
| `magnetic_field_moving` | 1.82 ns | 2.01 ns | 1.73 ns |
| `superposition_10` | 40.11 ns | 44.61 ns | 39.79 ns |
| `vector_magnitude` | 543.30 ps | 648.20 ps | 500.20 ps |
| `vector_cross` | 1.46 ns | 1.16 ns | 987.00 ps |
| `electric_dipole_field` | 1.81 ns | 5.93 ns | 5.53 ns |
| `magnetic_dipole_field` | 2.08 ns | 5.94 ns | 5.52 ns |
| `charged_sphere` | 1.56 ns | 1.34 ns | 1.23 ns |
| `ring_axis` | 1.55 ns | 1.39 ns | 1.23 ns |
| `disk_axis` | 1.56 ns | 1.39 ns | 1.23 ns |
| `trace_field_line_10` | 335.97 ns | 335.91 ns | 317.87 ns |

## maxwell

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `wave_speed` | 1.56 ns | 1.39 ns | 1.24 ns |
| `impedance` | 1.59 ns | 1.39 ns | 1.24 ns |
| `skin_depth` | 1.56 ns | 1.34 ns | 1.29 ns |

## charge

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `coulomb_force` | 1.82 ns | 5.36 ns | 5.54 ns |
| `lorentz_force` | 1.47 ns | 1.47 ns | 1.51 ns |

## wave

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `poynting_vector` | 1.46 ns | 1.15 ns | 984.30 ps |
| `plane_wave_intensity` | 516.10 ps | 561.40 ps | 517.80 ps |
| `plane_wave_e` | 801.90 ps | 536.30 ps | 512.80 ps |
| `fresnel_rs` | 2.01 ns | 1.31 ns | 1.25 ns |
| `snell_refraction` | 1.90 ns | 1.86 ns | 1.75 ns |
| `fresnel_rs_direct` | — | 1.84 ns | 1.73 ns |
| `reflectance_unpolarized` | — | 1.36 ns | 1.24 ns |
| `reflectance_at_interface` | — | 9.19 ns | 7.78 ns |
| `schlick_reflectance` | — | 1.24 ns | 1.24 ns |
| `half_wave_dipole_pattern` | 537.00 ps | 485.40 ps | 501.80 ps |
| `rect_waveguide_cutoff` | 1.67 ns | 1.23 ns | 1.24 ns |

## fdtd

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `step_100_cells` | 76.69 ns | 63.94 ns | 64.09 ns |
| `step_1000_cells` | 659.12 ns | 605.25 ns | 623.03 ns |
| `run_100_steps_500_cells` | — | 33.73 µs | 32.23 µs |
| `2d_step_100x100_tm` | — | 11.73 µs | 11.42 µs |
| `2d_run_50_steps_80x80` | — | 412.34 µs | 416.78 µs |
| `2d_step_100x100_tm_cpml10` | — | 20.76 µs | 20.47 µs |
| `2d_step_100x100_tm_lossy` | — | 14.26 µs | 14.31 µs |
| `3d_step_30x30x30` | — | 299.15 µs | 319.98 µs |
| `3d_step_30x30x30_cpml8` | — | — | 494.96 µs |
| `nf2ff_2d_80x80_180` | — | 206.43 µs | 193.01 µs |
| `find_eigenmodes_4096x256` | — | 2.27 ms | 2.12 ms |

## circuit

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `rc_charging` | 527.50 ps | 525.30 ps | 496.30 ps |
| `rl_current_rise` | 1.60 ns | 1.23 ns | 2.03 ns |
| `rlc_impedance` | 1.55 ns | 1.23 ns | 1.53 ns |
| `resonant_frequency` | 1.55 ns | 1.23 ns | 1.31 ns |
| `resistance_parallel` | 6.97 ns | 6.63 ns | 7.09 ns |

## material

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `polarization` | 1.47 ns | 1.01 ns | 995.80 ps |
| `displacement_field` | 1.43 ns | 1.02 ns | 991.90 ps |
| `h_field_from_b` | 2.07 ns | 5.75 ns | 5.51 ns |
| `clausius_mossotti` | 1.57 ns | 1.23 ns | 1.24 ns |
| `curie_weiss` | 1.55 ns | 1.22 ns | 1.24 ns |

## relativity

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `lorentz_factor` | 1.55 ns | 1.24 ns | 1.24 ns |
| `em_tensor_from_fields` | 8.20 ns | 8.14 ns | 8.11 ns |
| `lorentz_transform_x` | 3.81 ns | 6.36 ns | 5.91 ns |
| `lienard_wiechert_e` | 24.58 ns | 21.54 ns | 21.64 ns |
| `solve_retarded_time_uniform_t1s` | — | 97.24 ns | 103.48 ns |
| `lienard_wiechert_fields_orbit_b099` | — | 318.20 ns | 319.92 ns |
| `lienard_wiechert_scalar_potential` | — | 1.72 ns | 1.74 ns |
| `four_vector_boost` | 2.09 ns | 5.45 ns | 5.45 ns |

## polarization

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `jones_apply_polarizer` | — | 1.04 ns | 991.50 ps |
| `jones_qwp_apply` | — | 1.03 ns | 987.40 ps |
| `mueller_apply` | — | 995.80 ps | 994.00 ps |
| `stokes_from_jones` | — | 5.51 ns | 5.44 ns |

## scattering

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `mie_x1_real` | — | 458.47 ns | 423.87 ns |
| `mie_x10_real` | — | 902.11 ns | 859.83 ns |
| `mie_x50_m3_5` | — | 3.65 µs | 3.44 µs |
| `rayleigh_cross_section` | — | 1.89 ns | 1.74 ns |
| `rayleigh_phase_function` | — | 533.10 ps | 499.00 ps |

## beam

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `spot_size` | — | 536.50 ps | 494.70 ps |
| `abcd_propagate_beam` | — | 1.91 ns | 1.79 ns |
| `abcd_compose` | — | 1.07 ns | 989.20 ps |
| `hermite_gaussian_hg33` | — | 5.12 ns | 4.46 ns |
| `resonator_stability` | — | 6.71 ns | 6.48 ns |
| `hermite_gaussian_field_hg33` | — | 30.88 ns | 30.33 ns |

## rf

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `reflection_coefficient` | — | 1.23 ns | 1.23 ns |
| `vswr` | — | 1.23 ns | 1.23 ns |
| `linear_array_factor_8` | — | 557.00 ps | 552.30 ps |
| `s_cascade` | — | 96.47 ns | 100.25 ns |
| `input_impedance_complex` | — | 7.75 ns | 7.43 ns |
| `single_stub_match` | — | 133.74 ns | 128.46 ns |
| `parse_touchstone_s2p_201` | — | 69.89 µs | 79.42 µs |
| `write_touchstone_s2p_201` | — | 229.47 µs | 213.82 µs |

---

Generated by `./scripts/bench-history.sh`. Full history in `bench-history.csv`.
