# Benchmarks

Three-point tracking: **baseline** (first run) / **previous** / **latest**

| Point | Date | Commit |
|-------|------|--------|
| Baseline | 2026-03-25T12:34:25Z | `2df8cdf` |
| Previous | 2026-11-10T05:51:14Z | `678f256` |
| Latest | 2026-11-10T06:23:08Z | `678f256` |

## field

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `electric_field_point` | 2.07 ns | 1.65 ns | 1.70 ns |
| `electric_potential_point` | 1.55 ns | 1.31 ns | 1.31 ns |
| `magnetic_field_moving` | 1.82 ns | 1.84 ns | 2.01 ns |
| `superposition_10` | 40.11 ns | 42.17 ns | 44.61 ns |
| `vector_magnitude` | 543.30 ps | 535.00 ps | 648.20 ps |
| `vector_cross` | 1.46 ns | 1.04 ns | 1.16 ns |
| `electric_dipole_field` | 1.81 ns | 5.85 ns | 5.93 ns |
| `magnetic_dipole_field` | 2.08 ns | 5.84 ns | 5.94 ns |
| `charged_sphere` | 1.56 ns | 1.27 ns | 1.34 ns |
| `ring_axis` | 1.55 ns | 1.27 ns | 1.39 ns |
| `disk_axis` | 1.56 ns | 1.42 ns | 1.39 ns |
| `trace_field_line_10` | 335.97 ns | 335.02 ns | 335.91 ns |

## maxwell

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `wave_speed` | 1.56 ns | 1.30 ns | 1.39 ns |
| `impedance` | 1.59 ns | 1.30 ns | 1.39 ns |
| `skin_depth` | 1.56 ns | 1.31 ns | 1.34 ns |

## charge

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `coulomb_force` | 1.82 ns | 5.86 ns | 5.36 ns |
| `lorentz_force` | 1.47 ns | 1.64 ns | 1.47 ns |

## wave

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `poynting_vector` | 1.46 ns | 1.07 ns | 1.15 ns |
| `plane_wave_intensity` | 516.10 ps | 547.70 ps | 561.40 ps |
| `plane_wave_e` | 801.90 ps | 564.10 ps | 536.30 ps |
| `fresnel_rs` | 2.01 ns | 1.42 ns | 1.31 ns |
| `snell_refraction` | 1.90 ns | 1.91 ns | 1.86 ns |
| `fresnel_rs_direct` | — | 1.92 ns | 1.84 ns |
| `reflectance_unpolarized` | — | 1.37 ns | 1.36 ns |
| `reflectance_at_interface` | — | 9.47 ns | 9.19 ns |
| `schlick_reflectance` | — | 1.37 ns | 1.24 ns |
| `half_wave_dipole_pattern` | 537.00 ps | 546.50 ps | 485.40 ps |
| `rect_waveguide_cutoff` | 1.67 ns | 1.38 ns | 1.23 ns |

## fdtd

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `step_100_cells` | 76.69 ns | 73.42 ns | 63.94 ns |
| `step_1000_cells` | 659.12 ns | 697.47 ns | 605.25 ns |
| `run_100_steps_500_cells` | — | 39.09 µs | 33.73 µs |
| `2d_step_100x100_tm` | — | 13.58 µs | 11.73 µs |
| `2d_run_50_steps_80x80` | — | 469.92 µs | 412.34 µs |
| `2d_step_100x100_tm_cpml10` | — | 23.67 µs | 20.76 µs |
| `2d_step_100x100_tm_lossy` | — | 16.10 µs | 14.26 µs |
| `3d_step_30x30x30` | — | 360.28 µs | 299.15 µs |
| `nf2ff_2d_80x80_180` | — | 217.93 µs | 206.43 µs |
| `find_eigenmodes_4096x256` | — | 2.37 ms | 2.27 ms |

## circuit

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `rc_charging` | 527.50 ps | 580.20 ps | 525.30 ps |
| `rl_current_rise` | 1.60 ns | 1.39 ns | 1.23 ns |
| `rlc_impedance` | 1.55 ns | 1.41 ns | 1.23 ns |
| `resonant_frequency` | 1.55 ns | 1.46 ns | 1.23 ns |
| `resistance_parallel` | 6.97 ns | 7.34 ns | 6.63 ns |

## material

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `polarization` | 1.47 ns | 1.07 ns | 1.01 ns |
| `displacement_field` | 1.43 ns | 1.10 ns | 1.02 ns |
| `h_field_from_b` | 2.07 ns | 6.02 ns | 5.75 ns |
| `clausius_mossotti` | 1.57 ns | 1.37 ns | 1.23 ns |
| `curie_weiss` | 1.55 ns | 1.34 ns | 1.22 ns |

## relativity

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `lorentz_factor` | 1.55 ns | 1.35 ns | 1.24 ns |
| `em_tensor_from_fields` | 8.20 ns | 8.99 ns | 8.14 ns |
| `lorentz_transform_x` | 3.81 ns | 6.40 ns | 6.36 ns |
| `lienard_wiechert_e` | 24.58 ns | 23.28 ns | 21.54 ns |
| `solve_retarded_time_uniform_t1s` | — | 105.18 ns | 97.24 ns |
| `lienard_wiechert_fields_orbit_b099` | — | 348.54 ns | 318.20 ns |
| `lienard_wiechert_scalar_potential` | — | 1.87 ns | 1.72 ns |
| `four_vector_boost` | 2.09 ns | 5.88 ns | 5.45 ns |

## polarization

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `jones_apply_polarizer` | — | 1.07 ns | 1.04 ns |
| `jones_qwp_apply` | — | 1.06 ns | 1.03 ns |
| `mueller_apply` | — | 1.07 ns | 995.80 ps |
| `stokes_from_jones` | — | 6.01 ns | 5.51 ns |

## scattering

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `mie_x1_real` | — | 612.94 ns | 458.47 ns |
| `mie_x10_real` | — | 1.24 µs | 902.11 ns |
| `mie_x50_m3_5` | — | 4.75 µs | 3.65 µs |
| `rayleigh_cross_section` | — | 1.85 ns | 1.89 ns |
| `rayleigh_phase_function` | — | 554.40 ps | 533.10 ps |

## beam

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `spot_size` | — | 561.70 ps | 536.50 ps |
| `abcd_propagate_beam` | — | 2.01 ns | 1.91 ns |
| `abcd_compose` | — | 1.13 ns | 1.07 ns |
| `hermite_gaussian_hg33` | — | 5.22 ns | 5.12 ns |
| `resonator_stability` | — | 6.24 ns | 6.71 ns |
| `hermite_gaussian_field_hg33` | — | 30.55 ns | 30.88 ns |

## rf

| Benchmark | Baseline | Previous | Latest |
|-----------|----------|----------|--------|
| `reflection_coefficient` | — | 1.24 ns | 1.23 ns |
| `vswr` | — | 1.38 ns | 1.23 ns |
| `linear_array_factor_8` | — | 578.30 ps | 557.00 ps |
| `s_cascade` | — | 96.75 ns | 96.47 ns |
| `input_impedance_complex` | — | 8.26 ns | 7.75 ns |
| `single_stub_match` | — | 151.44 ns | 133.74 ns |
| `parse_touchstone_s2p_201` | — | 85.30 µs | 69.89 µs |
| `write_touchstone_s2p_201` | — | 248.75 µs | 229.47 µs |

---

Generated by `./scripts/bench-history.sh`. Full history in `bench-history.csv`.
