# Surface parity: Rust 1.2.0 and Cyrius 2.0

> **Generated** by `python3 -I -B scripts/surface_parity.py . --write`; a plain run fails when this is stale.
> Every public item of Rust bijli 1.2.0 (outside its test modules, from `scripts/rust_1_2_0.json`) beside the
> Cyrius name the port's naming rules give it (port plan §9.4), or the reason it has none of that form.

**735 items**, every one ported or explained.

## `ai` (15)

| Rust | Kind | Cyrius |
|---|---|---|
| `DaimonConfig` | struct | `BjDaimonConfig` |
| `impl Deserialize for DaimonConfig` | trait | `BjDaimonConfig_from_json` |
| `impl Serialize for DaimonConfig` | trait | `BjDaimonConfig_to_json` |
| `DaimonConfig::endpoint` | field | `BjDaimonConfig.endpoint` |
| `DaimonConfig::api_key` | field | `BjDaimonConfig.api_key` |
| `impl Debug for DaimonConfig` | trait | the method debug(t), into a BjText, the key redacted (src/ai.cyr) |
| `impl Default for DaimonConfig` | trait | `BjDaimonConfig_default` |
| `HooshConfig` | struct | `BjHooshConfig` |
| `impl Deserialize for HooshConfig` | trait | `BjHooshConfig_from_json` |
| `impl Serialize for HooshConfig` | trait | `BjHooshConfig_to_json` |
| `HooshConfig::endpoint` | field | `BjHooshConfig.endpoint` |
| `impl Default for HooshConfig` | trait | `BjHooshConfig_default` |
| `DaimonClient` | struct | `BjDaimonClient` |
| `DaimonClient::config` | fn | `BjDaimonClient::config` |
| `DaimonClient::new` | fn | `BjDaimonClient_new` |

## `beam` (49)

| Rust | Kind | Cyrius |
|---|---|---|
| `GaussianBeam` | struct | `BjGaussianBeam` |
| `impl Deserialize for GaussianBeam` | trait | `BjGaussianBeam_from_json` |
| `impl Serialize for GaussianBeam` | trait | `BjGaussianBeam_to_json` |
| `GaussianBeam::wavelength` | field | `BjGaussianBeam.wavelength` |
| `GaussianBeam::waist` | field | `BjGaussianBeam.waist` |
| `GaussianBeam::rayleigh_range` | field | `BjGaussianBeam.rayleigh_range` |
| `GaussianBeam::new` | fn | `BjGaussianBeam_new` |
| `GaussianBeam::q` | fn | `BjGaussianBeam::q` |
| `GaussianBeam::spot_size` | fn | `BjGaussianBeam::spot_size` |
| `GaussianBeam::radius_of_curvature` | fn | `BjGaussianBeam::radius_of_curvature` |
| `GaussianBeam::gouy_phase` | fn | `BjGaussianBeam::gouy_phase` |
| `GaussianBeam::divergence` | fn | `BjGaussianBeam::divergence` |
| `GaussianBeam::peak_intensity_factor` | fn | `BjGaussianBeam::peak_intensity_factor` |
| `GaussianBeam::confocal_parameter` | fn | `BjGaussianBeam::confocal_parameter` |
| `GaussianBeam::depth_of_focus` | fn | `BjGaussianBeam::depth_of_focus` |
| `AbcdMatrix` | struct | `BjAbcdMatrix` |
| `impl Deserialize for AbcdMatrix` | trait | `BjAbcdMatrix_from_json` |
| `impl Serialize for AbcdMatrix` | trait | `BjAbcdMatrix_to_json` |
| `AbcdMatrix::a` | field | `BjAbcdMatrix.a` |
| `AbcdMatrix::b` | field | `BjAbcdMatrix.b` |
| `AbcdMatrix::c` | field | `BjAbcdMatrix.c` |
| `AbcdMatrix::d` | field | `BjAbcdMatrix.d` |
| `AbcdMatrix::new` | fn | `BjAbcdMatrix_new` |
| `AbcdMatrix::identity` | fn | `BjAbcdMatrix_identity` |
| `AbcdMatrix::free_space` | fn | `BjAbcdMatrix_free_space` |
| `AbcdMatrix::thin_lens` | fn | `BjAbcdMatrix_thin_lens` |
| `AbcdMatrix::curved_mirror` | fn | `BjAbcdMatrix_curved_mirror` |
| `AbcdMatrix::flat_interface` | fn | `BjAbcdMatrix_flat_interface` |
| `AbcdMatrix::curved_interface` | fn | `BjAbcdMatrix_curved_interface` |
| `AbcdMatrix::apply_ray` | fn | `BjAbcdMatrix::apply_ray` |
| `AbcdMatrix::compose` | fn | `BjAbcdMatrix::compose` |
| `AbcdMatrix::determinant` | fn | `BjAbcdMatrix::determinant` |
| `AbcdMatrix::propagate_beam` | fn | `BjAbcdMatrix::propagate_beam` |
| `AbcdMatrix::transform_beam` | fn | `BjAbcdMatrix::transform_beam` |
| `AbcdMatrix::transform_beam_with_waist_position` | fn | `BjAbcdMatrix::transform_beam_with_waist_position` |
| `impl Default for AbcdMatrix` | trait | `BjAbcdMatrix_default` |
| `ResonatorStability` | struct | `BjResonatorStability` |
| `impl Deserialize for ResonatorStability` | trait | `BjResonatorStability_from_json` |
| `impl Serialize for ResonatorStability` | trait | `BjResonatorStability_to_json` |
| `ResonatorStability::g1` | field | `BjResonatorStability.g1` |
| `ResonatorStability::g2` | field | `BjResonatorStability.g2` |
| `ResonatorStability::g_product` | field | `BjResonatorStability.g_product` |
| `ResonatorStability::is_stable` | field | `BjResonatorStability.is_stable` |
| `resonator_stability` | fn | `bj_resonator_stability` |
| `resonator_beam_waist` | fn | `bj_resonator_beam_waist` |
| `transverse_mode_spacing` | fn | `bj_transverse_mode_spacing` |
| `hermite_gaussian` | fn | `bj_hermite_gaussian` |
| `hermite_gaussian_field` | fn | `bj_hermite_gaussian_field` |
| `laguerre_gaussian_intensity` | fn | `bj_laguerre_gaussian_intensity` |

## `bridge` (8)

| Rust | Kind | Cyrius |
|---|---|---|
| `joule_heating_power_density` | fn | `bj_joule_heating_power_density` |
| `em_energy_to_heating_rate` | fn | `bj_em_energy_to_heating_rate` |
| `wavelength_to_refractive_index_bk7` | fn | `bj_wavelength_to_refractive_index_bk7` |
| `e_field_to_optical_intensity` | fn | `bj_e_field_to_optical_intensity` |
| `energy_level_to_wavelength_nm` | fn | `bj_energy_level_to_wavelength_nm` |
| `nuclear_charge_potential` | fn | `bj_nuclear_charge_potential` |
| `e_field_to_piezo_stress` | fn | `bj_e_field_to_piezo_stress` |
| `magnetic_flux_to_magnetostrictive_strain` | fn | `bj_magnetic_flux_to_magnetostrictive_strain` |

## `charge` (24)

| Rust | Kind | Cyrius |
|---|---|---|
| `ELEMENTARY_CHARGE` | const | `BJ_ELEMENTARY_CHARGE` |
| `ELECTRON_MASS` | const | `BJ_ELECTRON_MASS` |
| `PROTON_MASS` | const | `BJ_PROTON_MASS` |
| `PointCharge` | struct | `BjPointCharge` |
| `impl Deserialize for PointCharge` | trait | `BjPointCharge_from_json` |
| `impl Serialize for PointCharge` | trait | `BjPointCharge_to_json` |
| `PointCharge::charge` | field | `BjPointCharge.charge` |
| `PointCharge::mass` | field | `BjPointCharge.mass` |
| `PointCharge::position` | field | `BjPointCharge.position` |
| `PointCharge::velocity` | field | `BjPointCharge.velocity` |
| `PointCharge::new` | fn | `BjPointCharge_new` |
| `PointCharge::electron` | fn | `BjPointCharge_electron` |
| `PointCharge::proton` | fn | `BjPointCharge_proton` |
| `PointCharge::electric_field_at` | fn | `BjPointCharge::electric_field_at` |
| `PointCharge::electric_potential_at` | fn | `BjPointCharge::electric_potential_at` |
| `PointCharge::magnetic_field_at` | fn | `BjPointCharge::magnetic_field_at` |
| `PointCharge::lorentz_force_in` | fn | `BjPointCharge::lorentz_force_in` |
| `coulomb_force` | fn | `bj_coulomb_force` |
| `lorentz_force` | fn | `bj_lorentz_force` |
| `dipole_moment` | fn | `bj_dipole_moment` |
| `potential_energy` | fn | `bj_potential_energy` |
| `coulomb_potential_energy` | fn | `bj_coulomb_potential_energy` |
| `cyclotron_frequency` | fn | `bj_cyclotron_frequency` |
| `larmor_radius` | fn | `bj_larmor_radius` |

## `circuit` (27)

| Rust | Kind | Cyrius |
|---|---|---|
| `voltage` | fn | `bj_voltage` |
| `current` | fn | `bj_current` |
| `resistance` | fn | `bj_resistance` |
| `power` | fn | `bj_power` |
| `resistance_from_geometry` | fn | `bj_resistance_from_geometry` |
| `resistance_series` | fn | `bj_resistance_series` |
| `resistance_parallel` | fn | `bj_resistance_parallel` |
| `capacitance_series` | fn | `bj_capacitance_series` |
| `capacitance_parallel` | fn | `bj_capacitance_parallel` |
| `inductance_series` | fn | `bj_inductance_series` |
| `inductance_parallel` | fn | `bj_inductance_parallel` |
| `parallel_plate_capacitance` | fn | `bj_parallel_plate_capacitance` |
| `capacitor_energy` | fn | `bj_capacitor_energy` |
| `capacitor_charge` | fn | `bj_capacitor_charge` |
| `solenoid_inductance` | fn | `bj_solenoid_inductance` |
| `inductor_energy` | fn | `bj_inductor_energy` |
| `rc_time_constant` | fn | `bj_rc_time_constant` |
| `rc_charging_voltage` | fn | `bj_rc_charging_voltage` |
| `rc_discharging_voltage` | fn | `bj_rc_discharging_voltage` |
| `rl_time_constant` | fn | `bj_rl_time_constant` |
| `rl_current_rise` | fn | `bj_rl_current_rise` |
| `rl_current_decay` | fn | `bj_rl_current_decay` |
| `resonant_frequency` | fn | `bj_resonant_frequency` |
| `resonant_angular_frequency` | fn | `bj_resonant_angular_frequency` |
| `quality_factor` | fn | `bj_quality_factor` |
| `damping_ratio` | fn | `bj_damping_ratio` |
| `rlc_impedance` | fn | `bj_rlc_impedance` |

## `compute` (19)

| Rust | Kind | Cyrius |
|---|---|---|
| `CpuBackend` | struct | `BjCpuBackend` |
| `impl ComputeBackend for CpuBackend` | trait | `impl BjComputeBackend for BjCpuBackend`: a Cyrius trait impl, checked by the compiler against the trait (src/compute.cyr) |
| `GpuFdtd2d` | struct | `BjGpuFdtd2d` |
| `GpuFdtd2d::backend` | field | `BjGpuFdtd2d.backend` |
| `GpuFdtd2d::nx` | field | `BjGpuFdtd2d.nx` |
| `GpuFdtd2d::ny` | field | `BjGpuFdtd2d.ny` |
| `GpuFdtd2d::dx` | field | `BjGpuFdtd2d.dx` |
| `GpuFdtd2d::dt` | field | `BjGpuFdtd2d.dt` |
| `GpuFdtd2d::step` | field | `BjGpuFdtd2d.step` |
| `GpuFdtd2d::new` | fn | `BjGpuFdtd2d_new` |
| `GpuFdtd2d::step_once` | fn | `BjGpuFdtd2d_step_once` |
| `GpuFdtd2d::add_source` | fn | `BjGpuFdtd2d_add_source` |
| `GpuFdtd2d::run` | fn | `BjGpuFdtd2d_run` |
| `GpuFdtd2d::time` | fn | `BjGpuFdtd2d_time` |
| `GpuFdtd2d::read_ez` | fn | `BjGpuFdtd2d_read_ez` |
| `GpuFdtd2d::read_hx` | fn | `BjGpuFdtd2d_read_hx` |
| `GpuFdtd2d::read_hy` | fn | `BjGpuFdtd2d_read_hy` |
| `WgpuBackend` | struct | `BjWgpuBackend` |
| `WgpuBackend::new` | fn | `BjWgpuBackend_new` |

## `error` (9)

| Rust | Kind | Cyrius |
|---|---|---|
| `BijliError` | enum | `BjError` |
| `BijliError::DimensionMismatch` | variant | `BJ_ERR_DIMENSION_MISMATCH` |
| `BijliError::ZeroCharge` | variant | `BJ_ERR_ZERO_CHARGE` |
| `BijliError::Singularity` | variant | `BJ_ERR_SINGULARITY` |
| `BijliError::InsufficientResolution` | variant | `BJ_ERR_INSUFFICIENT_RESOLUTION` |
| `BijliError::InvalidPermittivity` | variant | `BJ_ERR_INVALID_PERMITTIVITY` |
| `BijliError::InvalidPermeability` | variant | `BJ_ERR_INVALID_PERMEABILITY` |
| `BijliError::DivisionByZero` | variant | `BJ_ERR_DIVISION_BY_ZERO` |
| `BijliError::InvalidParameter` | variant | `BJ_ERR_INVALID_PARAMETER` |

## `fdtd` (187)

| Rust | Kind | Cyrius |
|---|---|---|
| `Fdtd1d` | struct | `BjFdtd1d` |
| `Fdtd1d::e_field` | field | `BjFdtd1d.e_field` |
| `Fdtd1d::h_field` | field | `BjFdtd1d.h_field` |
| `Fdtd1d::dx` | field | `BjFdtd1d.dx` |
| `Fdtd1d::dt` | field | `BjFdtd1d.dt` |
| `Fdtd1d::step` | field | `BjFdtd1d.step` |
| `Fdtd1d::num_cells` | field | `BjFdtd1d.num_cells` |
| `Fdtd1d::new` | fn | `BjFdtd1d_new` |
| `Fdtd1d::set_permittivity` | fn | `BjFdtd1d::set_permittivity` |
| `Fdtd1d::set_permeability` | fn | `BjFdtd1d::set_permeability` |
| `Fdtd1d::add_source` | fn | `BjFdtd1d::add_source` |
| `Fdtd1d::step_once` | fn | `BjFdtd1d::step_once` |
| `Fdtd1d::run` | fn | `BjFdtd1d::run` |
| `Fdtd1d::time` | fn | `BjFdtd1d::time` |
| `Fdtd1d::total_energy` | fn | `BjFdtd1d::total_energy` |
| `gaussian_pulse` | fn | `bj_gaussian_pulse` |
| `ricker_wavelet` | fn | `bj_ricker_wavelet` |
| `differentiated_gaussian` | fn | `bj_differentiated_gaussian` |
| `modulated_gaussian_pulse` | fn | `bj_modulated_gaussian_pulse` |
| `TfsfSource` | struct | `BjTfsfSource` |
| `TfsfSource::inc_e` | field | `BjTfsfSource.inc_e` |
| `TfsfSource::inc_h` | field | `BjTfsfSource.inc_h` |
| `TfsfSource::boundary_x` | field | `BjTfsfSource.boundary_x` |
| `TfsfSource::source_cell` | field | `BjTfsfSource.source_cell` |
| `TfsfSource::new` | fn | `BjTfsfSource_new` |
| `TfsfSource::try_new` | fn | `BjTfsfSource_try_new` |
| `TfsfSource::step` | fn | `BjTfsfSource::step` |
| `TfsfSource::incident_e` | fn | `BjTfsfSource::incident_e` |
| `TfsfSource::incident_h` | fn | `BjTfsfSource::incident_h` |
| `LossyMaterial` | struct | `BjLossyMaterial` |
| `impl Deserialize for LossyMaterial` | trait | `BjLossyMaterial_from_json` |
| `impl Serialize for LossyMaterial` | trait | `BjLossyMaterial_to_json` |
| `LossyMaterial::eps_r` | field | `BjLossyMaterial.eps_r` |
| `LossyMaterial::conductivity` | field | `BjLossyMaterial.conductivity` |
| `LossyMaterial::new` | fn | `BjLossyMaterial_new` |
| `LossyMaterial::try_new` | fn | `BjLossyMaterial_try_new` |
| `LossyMaterial::update_coefficients` | fn | `BjLossyMaterial::update_coefficients` |
| `DrudeMaterial` | struct | `BjDrudeMaterial` |
| `impl Deserialize for DrudeMaterial` | trait | `BjDrudeMaterial_from_json` |
| `impl Serialize for DrudeMaterial` | trait | `BjDrudeMaterial_to_json` |
| `DrudeMaterial::eps_inf` | field | `BjDrudeMaterial.eps_inf` |
| `DrudeMaterial::omega_p` | field | `BjDrudeMaterial.omega_p` |
| `DrudeMaterial::gamma` | field | `BjDrudeMaterial.gamma` |
| `DrudeMaterial::new` | fn | `BjDrudeMaterial_new` |
| `DrudeMaterial::try_new` | fn | `BjDrudeMaterial_try_new` |
| `DrudeMaterial::gold` | fn | `BjDrudeMaterial_gold` |
| `DrudeMaterial::silver` | fn | `BjDrudeMaterial_silver` |
| `DrudeMaterial::ade_coefficients` | fn | `BjDrudeMaterial::ade_coefficients` |
| `LorentzMaterial` | struct | `BjLorentzMaterial` |
| `impl Deserialize for LorentzMaterial` | trait | `BjLorentzMaterial_from_json` |
| `impl Serialize for LorentzMaterial` | trait | `BjLorentzMaterial_to_json` |
| `LorentzMaterial::eps_inf` | field | `BjLorentzMaterial.eps_inf` |
| `LorentzMaterial::delta_eps` | field | `BjLorentzMaterial.delta_eps` |
| `LorentzMaterial::omega_0` | field | `BjLorentzMaterial.omega_0` |
| `LorentzMaterial::gamma` | field | `BjLorentzMaterial.gamma` |
| `LorentzMaterial::new` | fn | `BjLorentzMaterial_new` |
| `LorentzMaterial::try_new` | fn | `BjLorentzMaterial_try_new` |
| `LorentzMaterial::ade_coefficients` | fn | `BjLorentzMaterial::ade_coefficients` |
| `DebyeMaterial` | struct | `BjDebyeMaterial` |
| `impl Deserialize for DebyeMaterial` | trait | `BjDebyeMaterial_from_json` |
| `impl Serialize for DebyeMaterial` | trait | `BjDebyeMaterial_to_json` |
| `DebyeMaterial::eps_inf` | field | `BjDebyeMaterial.eps_inf` |
| `DebyeMaterial::eps_s` | field | `BjDebyeMaterial.eps_s` |
| `DebyeMaterial::tau` | field | `BjDebyeMaterial.tau` |
| `DebyeMaterial::new` | fn | `BjDebyeMaterial_new` |
| `DebyeMaterial::try_new` | fn | `BjDebyeMaterial_try_new` |
| `DebyeMaterial::water` | fn | `BjDebyeMaterial_water` |
| `DebyeMaterial::ade_coefficients` | fn | `BjDebyeMaterial::ade_coefficients` |
| `Nf2ffResult` | struct | `BjNf2ffResult` |
| `Nf2ffResult::pattern` | field | `BjNf2ffResult.pattern` |
| `Nf2ffResult::angles` | field | `BjNf2ffResult.angles` |
| `nf2ff_2d` | fn | `bj_nf2ff_2d` |
| `nf2ff_2d_phasor` | fn | `bj_nf2ff_2d_phasor` |
| `subpixel_permittivity_parallel` | fn | `bj_subpixel_permittivity_parallel` |
| `subpixel_permittivity_perpendicular` | fn | `bj_subpixel_permittivity_perpendicular` |
| `subpixel_permittivity_perpendicular_checked` | fn | `bj_subpixel_permittivity_perpendicular_checked` |
| `CpmlParams` | struct | `BjCpmlParams` |
| `impl Deserialize for CpmlParams` | trait | `BjCpmlParams_from_json` |
| `impl Serialize for CpmlParams` | trait | `BjCpmlParams_to_json` |
| `CpmlParams::thickness` | field | `BjCpmlParams.thickness` |
| `CpmlParams::grading_order` | field | `BjCpmlParams.grading_order` |
| `CpmlParams::sigma_ratio` | field | `BjCpmlParams.sigma_ratio` |
| `CpmlParams::kappa_max` | field | `BjCpmlParams.kappa_max` |
| `CpmlParams::alpha_max` | field | `BjCpmlParams.alpha_max` |
| `CpmlParams::new` | fn | `BjCpmlParams_new` |
| `Cpml` | struct | `BjCpml` |
| `Cpml::thickness` | field | `BjCpml.thickness` |
| `Cpml::new` | fn | `BjCpml_new` |
| `Cpml::with_params` | fn | `BjCpml_with_params` |
| `Cpml::attach` | fn | `BjCpml::attach` |
| `Cpml::sigma_max` | fn | `BjCpml::sigma_max` |
| `Cpml::kappa_max` | fn | `BjCpml::kappa_max` |
| `Cpml::alpha_max` | fn | `BjCpml::alpha_max` |
| `Cpml::grading_order` | fn | `BjCpml::grading_order` |
| `Mode2d` | enum | `BjMode2d` |
| `impl Deserialize for Mode2d` | trait | `BjMode2d_from_json` |
| `impl Serialize for Mode2d` | trait | `BjMode2d_to_json` |
| `Mode2d::Tm` | variant | `BJ_MODE2D_TM` |
| `Mode2d::Te` | variant | `BJ_MODE2D_TE` |
| `Fdtd2d` | struct | `BjFdtd2d` |
| `Fdtd2d::mode` | field | `BjFdtd2d.mode` |
| `Fdtd2d::nx` | field | `BjFdtd2d.nx` |
| `Fdtd2d::ny` | field | `BjFdtd2d.ny` |
| `Fdtd2d::dx` | field | `BjFdtd2d.dx` |
| `Fdtd2d::dt` | field | `BjFdtd2d.dt` |
| `Fdtd2d::step` | field | `BjFdtd2d.step` |
| `Fdtd2d::fz` | field | `BjFdtd2d.fz` |
| `Fdtd2d::fx` | field | `BjFdtd2d.fx` |
| `Fdtd2d::fy` | field | `BjFdtd2d.fy` |
| `Fdtd2d::new` | fn | `BjFdtd2d_new` |
| `Fdtd2d::idx` | fn | `BjFdtd2d::idx` |
| `Fdtd2d::set_permittivity` | fn | `BjFdtd2d::set_permittivity` |
| `Fdtd2d::set_permeability` | fn | `BjFdtd2d::set_permeability` |
| `Fdtd2d::set_conductivity` | fn | `BjFdtd2d::set_conductivity` |
| `Fdtd2d::set_material` | fn | `BjFdtd2d::set_material` |
| `Fdtd2d::add_source` | fn | `BjFdtd2d::add_source` |
| `Fdtd2d::try_enable_cpml` | fn | `BjFdtd2d::try_enable_cpml` |
| `Fdtd2d::try_enable_cpml_with` | fn | `BjFdtd2d::try_enable_cpml_with` |
| `Fdtd2d::enable_cpml` | fn | `BjFdtd2d::enable_cpml` |
| `Fdtd2d::disable_cpml` | fn | `BjFdtd2d::disable_cpml` |
| `Fdtd2d::cpml` | fn | `BjFdtd2d::cpml` |
| `Fdtd2d::step_once` | fn | `BjFdtd2d::step_once` |
| `Fdtd2d::run` | fn | `BjFdtd2d::run` |
| `Fdtd2d::time` | fn | `BjFdtd2d::time` |
| `Fdtd2d::total_energy` | fn | `BjFdtd2d::total_energy` |
| `Fdtd3d` | struct | `BjFdtd3d` |
| `Fdtd3d::nx` | field | `BjFdtd3d.nx` |
| `Fdtd3d::ny` | field | `BjFdtd3d.ny` |
| `Fdtd3d::nz` | field | `BjFdtd3d.nz` |
| `Fdtd3d::dx` | field | `BjFdtd3d.dx` |
| `Fdtd3d::dt` | field | `BjFdtd3d.dt` |
| `Fdtd3d::step` | field | `BjFdtd3d.step` |
| `Fdtd3d::ex` | field | `BjFdtd3d.ex` |
| `Fdtd3d::ey` | field | `BjFdtd3d.ey` |
| `Fdtd3d::ez` | field | `BjFdtd3d.ez` |
| `Fdtd3d::hx` | field | `BjFdtd3d.hx` |
| `Fdtd3d::hy` | field | `BjFdtd3d.hy` |
| `Fdtd3d::hz` | field | `BjFdtd3d.hz` |
| `Fdtd3d::new` | fn | `BjFdtd3d_new` |
| `Fdtd3d::idx` | fn | `BjFdtd3d::idx` |
| `Fdtd3d::set_permittivity` | fn | `BjFdtd3d::set_permittivity` |
| `Fdtd3d::set_permeability` | fn | `BjFdtd3d::set_permeability` |
| `Fdtd3d::set_conductivity` | fn | `BjFdtd3d::set_conductivity` |
| `Fdtd3d::set_material` | fn | `BjFdtd3d::set_material` |
| `Fdtd3d::add_source_ez` | fn | `BjFdtd3d::add_source_ez` |
| `Fdtd3d::try_enable_cpml` | fn | `BjFdtd3d::try_enable_cpml` |
| `Fdtd3d::try_enable_cpml_with` | fn | `BjFdtd3d::try_enable_cpml_with` |
| `Fdtd3d::disable_cpml` | fn | `BjFdtd3d::disable_cpml` |
| `Fdtd3d::cpml` | fn | `BjFdtd3d::cpml` |
| `Fdtd3d::step_once` | fn | `BjFdtd3d::step_once` |
| `Fdtd3d::run` | fn | `BjFdtd3d::run` |
| `Fdtd3d::time` | fn | `BjFdtd3d::time` |
| `Fdtd3d::total_energy` | fn | `BjFdtd3d::total_energy` |
| `NonUniformGrid` | struct | `BjNonUniformGrid` |
| `impl Deserialize for NonUniformGrid` | trait | `BjNonUniformGrid_from_json` |
| `impl Serialize for NonUniformGrid` | trait | `BjNonUniformGrid_to_json` |
| `NonUniformGrid::dx` | field | `BjNonUniformGrid.dx` |
| `NonUniformGrid::dy` | field | `BjNonUniformGrid.dy` |
| `NonUniformGrid::dz` | field | `BjNonUniformGrid.dz` |
| `NonUniformGrid::uniform` | fn | `BjNonUniformGrid_uniform` |
| `NonUniformGrid::try_uniform` | fn | `BjNonUniformGrid_try_uniform` |
| `NonUniformGrid::validate` | fn | `BjNonUniformGrid::validate` |
| `NonUniformGrid::graded` | fn | `BjNonUniformGrid_graded` |
| `NonUniformGrid::min_spacing` | fn | `BjNonUniformGrid::min_spacing` |
| `NonUniformGrid::courant_dt` | fn | `BjNonUniformGrid::courant_dt` |
| `NonUniformGrid::try_courant_dt` | fn | `BjNonUniformGrid::try_courant_dt` |
| `NonUniformGrid::extents` | fn | `BjNonUniformGrid::extents` |
| `FOURTH_ORDER_C1` | const | `BJ_FOURTH_ORDER_C1` |
| `FOURTH_ORDER_C2` | const | `BJ_FOURTH_ORDER_C2` |
| `fourth_order_derivative` | fn | `bj_fourth_order_derivative` |
| `CwAccumulator` | struct | `BjCwAccumulator` |
| `CwAccumulator::real` | field | `BjCwAccumulator.real` |
| `CwAccumulator::imag` | field | `BjCwAccumulator.imag` |
| `CwAccumulator::n_samples` | field | `BjCwAccumulator.n_samples` |
| `CwAccumulator::new` | fn | `BjCwAccumulator_new` |
| `CwAccumulator::try_new` | fn | `BjCwAccumulator_try_new` |
| `CwAccumulator::accumulate` | fn | `BjCwAccumulator::accumulate` |
| `CwAccumulator::try_accumulate` | fn | `BjCwAccumulator::try_accumulate` |
| `CwAccumulator::amplitude` | fn | `BjCwAccumulator::amplitude` |
| `CwAccumulator::phase` | fn | `BjCwAccumulator::phase` |
| `EigenmodeResult` | struct | `BjEigenmodeResult` |
| `impl Deserialize for EigenmodeResult` | trait | `BjEigenmodeResult_from_json` |
| `impl Serialize for EigenmodeResult` | trait | `BjEigenmodeResult_to_json` |
| `EigenmodeResult::frequencies` | field | `BjEigenmodeResult.frequencies` |
| `EigenmodeResult::quality_factors` | field | `BjEigenmodeResult.quality_factors` |
| `EigenmodeResult::amplitudes` | field | `BjEigenmodeResult.amplitudes` |
| `find_eigenmodes` | fn | `bj_find_eigenmodes` |

## `field` (60)

| Rust | Kind | Cyrius |
|---|---|---|
| `EPSILON_0` | const | `BJ_EPSILON_0` |
| `MU_0` | const | `BJ_MU_0` |
| `SPEED_OF_LIGHT` | const | `BJ_SPEED_OF_LIGHT` |
| `COULOMB_K` | const | `BJ_COULOMB_K` |
| `SPEED_OF_LIGHT_SQ` | const | `BJ_SPEED_OF_LIGHT_SQ` |
| `FieldVector` | struct | `BjFieldVector` |
| `impl Deserialize for FieldVector` | trait | `BjFieldVector_from_json` |
| `impl Serialize for FieldVector` | trait | `BjFieldVector_to_json` |
| `FieldVector::x` | field | `BjFieldVector.x` |
| `FieldVector::y` | field | `BjFieldVector.y` |
| `FieldVector::z` | field | `BjFieldVector.z` |
| `FieldVector::new` | fn | `BjFieldVector_new` |
| `FieldVector::zero` | fn | `BjFieldVector_zero` |
| `FieldVector::magnitude` | fn | `BjFieldVector::magnitude` |
| `FieldVector::magnitude_sq` | fn | `BjFieldVector::magnitude_sq` |
| `FieldVector::normalized` | fn | `BjFieldVector::normalized` |
| `FieldVector::dot` | fn | `BjFieldVector::dot` |
| `FieldVector::cross` | fn | `BjFieldVector::cross` |
| `FieldVector::scale` | fn | `BjFieldVector::scale` |
| `impl Default for FieldVector` | trait | `BjFieldVector_default` |
| `impl From for FieldVector` | trait | BjFieldVector_from_array(p); the other direction, From<FieldVector> for [f64; 3], is the method to_array(p) |
| `impl Add for FieldVector` | trait | `BjFieldVector_add` |
| `impl Sub for FieldVector` | trait | `BjFieldVector_sub` |
| `impl Mul for FieldVector` | trait | `BjFieldVector_mul` |
| `impl Neg for FieldVector` | trait | `BjFieldVector_neg` |
| `impl AddAssign for FieldVector` | trait | `a = a + b`: Cyrius refuses `+=` on a struct value (guide, Compound assignment); the same operations |
| `impl SubAssign for FieldVector` | trait | `a = a - b`, as AddAssign |
| `impl MulAssign for FieldVector` | trait | `a = a * s`, as AddAssign |
| `impl Mul for f64` | trait | BjFieldVector_scale_left(s, v): an operator fn belongs to its left operand's type |
| `impl Display for FieldVector` | trait | the method display(t), into a BjText (src/fmtf.cyr) |
| `electric_field_point_charge` | fn | `bj_electric_field_point_charge` |
| `electric_potential_point_charge` | fn | `bj_electric_potential_point_charge` |
| `magnetic_field_moving_charge` | fn | `bj_magnetic_field_moving_charge` |
| `electric_field_superposition` | fn | `bj_electric_field_superposition` |
| `electric_energy_density` | fn | `bj_electric_energy_density` |
| `magnetic_energy_density` | fn | `bj_magnetic_energy_density` |
| `electric_dipole_field` | fn | `bj_electric_dipole_field` |
| `electric_dipole_potential` | fn | `bj_electric_dipole_potential` |
| `magnetic_dipole_field` | fn | `bj_magnetic_dipole_field` |
| `electric_field_infinite_plane` | fn | `bj_electric_field_infinite_plane` |
| `check_radial_distance` | fn | `_bj_check_radial_distance` |
| `electric_field_charged_sphere` | fn | `bj_electric_field_charged_sphere` |
| `electric_field_inside_sphere` | fn | `bj_electric_field_inside_sphere` |
| `electric_field_infinite_line` | fn | `bj_electric_field_infinite_line` |
| `electric_field_outside_cylinder` | fn | `bj_electric_field_outside_cylinder` |
| `electric_field_inside_cylinder` | fn | `bj_electric_field_inside_cylinder` |
| `electric_field_ring_axis` | fn | `bj_electric_field_ring_axis` |
| `electric_field_disk_axis` | fn | `bj_electric_field_disk_axis` |
| `trace_field_line` | fn | `bj_trace_field_line` |
| `electric_quadrupole_potential` | fn | `bj_electric_quadrupole_potential` |
| `multipole_moment` | fn | `bj_multipole_moment` |
| `spherical_multipole_moment` | fn | `bj_spherical_multipole_moment` |
| `image_charge_plane` | fn | `bj_image_charge_plane` |
| `image_charge_sphere` | fn | `bj_image_charge_sphere` |
| `image_charge_dielectric` | fn | `bj_image_charge_dielectric` |
| `greens_function_scalar` | fn | `bj_greens_function_scalar` |
| `greens_function_static` | fn | `bj_greens_function_static` |
| `maxwell_stress_tensor` | fn | `bj_maxwell_stress_tensor` |
| `stress_tensor_force` | fn | `bj_stress_tensor_force` |
| `poynting_flux_surface` | fn | `bj_poynting_flux_surface` |

## `hisab (re-exported)` (2)

| Rust | Kind | Cyrius |
|---|---|---|
| `impl Serialize for Complex` | trait | `BjComplex_to_json` |
| `impl Deserialize for Complex` | trait | `BjComplex_from_json` |

## `integration/soorat` (59)

| Rust | Kind | Cyrius |
|---|---|---|
| `FieldSlice2D` | struct | `BjFieldSlice2D` |
| `impl Deserialize for FieldSlice2D` | trait | `BjFieldSlice2D_from_json` |
| `impl Serialize for FieldSlice2D` | trait | `BjFieldSlice2D_to_json` |
| `FieldSlice2D::values` | field | `BjFieldSlice2D.values` |
| `FieldSlice2D::dimensions` | field | `BjFieldSlice2D.dimensions` |
| `FieldSlice2D::origin` | field | `BjFieldSlice2D.origin` |
| `FieldSlice2D::spacing` | field | `BjFieldSlice2D.spacing` |
| `FieldSlice2D::component` | field | `BjFieldSlice2D.component` |
| `FieldSlice2D::step` | field | `BjFieldSlice2D.step` |
| `FieldSlice3D` | struct | `BjFieldSlice3D` |
| `impl Deserialize for FieldSlice3D` | trait | `BjFieldSlice3D_from_json` |
| `impl Serialize for FieldSlice3D` | trait | `BjFieldSlice3D_to_json` |
| `FieldSlice3D::values` | field | `BjFieldSlice3D.values` |
| `FieldSlice3D::dimensions` | field | `BjFieldSlice3D.dimensions` |
| `FieldSlice3D::origin` | field | `BjFieldSlice3D.origin` |
| `FieldSlice3D::spacing` | field | `BjFieldSlice3D.spacing` |
| `FieldSlice3D::component` | field | `BjFieldSlice3D.component` |
| `FieldSlice3D::step` | field | `BjFieldSlice3D.step` |
| `FieldSlice2D::from_fdtd_2d` | fn | `BjFieldSlice2D_from_fdtd_2d` |
| `FieldSlice2D::try_from_fdtd_2d` | fn | `BjFieldSlice2D_try_from_fdtd_2d` |
| `FieldSlice3D::from_fdtd_3d` | fn | `BjFieldSlice3D_from_fdtd_3d` |
| `FieldSlice3D::try_from_fdtd_3d` | fn | `BjFieldSlice3D_try_from_fdtd_3d` |
| `FieldLineVisualization` | struct | `BjFieldLineVisualization` |
| `impl Deserialize for FieldLineVisualization` | trait | `BjFieldLineVisualization_from_json` |
| `impl Serialize for FieldLineVisualization` | trait | `BjFieldLineVisualization_to_json` |
| `FieldLineVisualization::lines` | field | `BjFieldLineVisualization.lines` |
| `FieldLine` | struct | `BjFieldLine` |
| `impl Deserialize for FieldLine` | trait | `BjFieldLine_from_json` |
| `impl Serialize for FieldLine` | trait | `BjFieldLine_to_json` |
| `FieldLine::points` | field | `BjFieldLine.points` |
| `FieldLine::magnitudes` | field | `BjFieldLine.magnitudes` |
| `ChargeVisualization` | struct | `BjChargeVisualization` |
| `impl Deserialize for ChargeVisualization` | trait | `BjChargeVisualization_from_json` |
| `impl Serialize for ChargeVisualization` | trait | `BjChargeVisualization_to_json` |
| `ChargeVisualization::charges` | field | `BjChargeVisualization.charges` |
| `ChargePoint` | struct | `BjChargePoint` |
| `impl Deserialize for ChargePoint` | trait | `BjChargePoint_from_json` |
| `impl Serialize for ChargePoint` | trait | `BjChargePoint_to_json` |
| `ChargePoint::position` | field | `BjChargePoint.position_x, position_y, position_z` |
| `ChargePoint::charge` | field | `BjChargePoint.charge` |
| `ChargePoint::magnitude` | field | `BjChargePoint.magnitude` |
| `ChargeVisualization::from_point_charges` | fn | `BjChargeVisualization_from_point_charges` |
| `RadiationPattern` | struct | `BjRadiationPattern` |
| `impl Deserialize for RadiationPattern` | trait | `BjRadiationPattern_from_json` |
| `impl Serialize for RadiationPattern` | trait | `BjRadiationPattern_to_json` |
| `RadiationPattern::angles` | field | `BjRadiationPattern.angles` |
| `RadiationPattern::pattern` | field | `BjRadiationPattern.pattern` |
| `RadiationPattern::max_gain` | field | `BjRadiationPattern.max_gain` |
| `RadiationPattern::from_nf2ff` | fn | `BjRadiationPattern_from_nf2ff` |
| `VectorFieldSample` | struct | `BjVectorFieldSample` |
| `impl Deserialize for VectorFieldSample` | trait | `BjVectorFieldSample_from_json` |
| `impl Serialize for VectorFieldSample` | trait | `BjVectorFieldSample_to_json` |
| `VectorFieldSample::vectors` | field | `BjVectorFieldSample.vectors` |
| `VectorFieldSample::dimensions` | field | `BjVectorFieldSample.dimensions` |
| `VectorFieldSample::origin` | field | `BjVectorFieldSample.origin` |
| `VectorFieldSample::spacing` | field | `BjVectorFieldSample.spacing` |
| `VectorFieldSample::max_magnitude` | field | `BjVectorFieldSample.max_magnitude` |
| `VectorFieldSample::from_components` | fn | `BjVectorFieldSample_from_components` |
| `VectorFieldSample::try_from_components` | fn | `BjVectorFieldSample_try_from_components` |

## `logging` (2)

| Rust | Kind | Cyrius |
|---|---|---|
| `init` | fn | `bj_log_init` |
| `init_with_level` | fn | `bj_log_init_with_level` |

## `material` (44)

| Rust | Kind | Cyrius |
|---|---|---|
| `electric_susceptibility` | fn | `bj_electric_susceptibility` |
| `relative_permittivity` | fn | `bj_relative_permittivity` |
| `absolute_permittivity` | fn | `bj_absolute_permittivity` |
| `polarization` | fn | `bj_polarization` |
| `displacement_field` | fn | `bj_displacement_field` |
| `bound_surface_charge` | fn | `bj_bound_surface_charge` |
| `bound_volume_charge` | fn | `bj_bound_volume_charge` |
| `dielectric_energy_density` | fn | `bj_dielectric_energy_density` |
| `clausius_mossotti` | fn | `bj_clausius_mossotti` |
| `magnetic_susceptibility` | fn | `bj_magnetic_susceptibility` |
| `relative_permeability` | fn | `bj_relative_permeability` |
| `absolute_permeability` | fn | `bj_absolute_permeability` |
| `magnetization` | fn | `bj_magnetization` |
| `h_field_from_b` | fn | `bj_h_field_from_b` |
| `b_field_from_h` | fn | `bj_b_field_from_h` |
| `magnetic_energy_density_material` | fn | `bj_magnetic_energy_density_material` |
| `bound_surface_current` | fn | `bj_bound_surface_current` |
| `curie_law` | fn | `bj_curie_law` |
| `curie_weiss` | fn | `bj_curie_weiss` |
| `MagneticType` | enum | `BjMagneticType` |
| `MagneticType::Diamagnetic` | variant | `BJ_MAGNETIC_TYPE_DIAMAGNETIC` |
| `MagneticType::Paramagnetic` | variant | `BJ_MAGNETIC_TYPE_PARAMAGNETIC` |
| `MagneticType::Ferromagnetic` | variant | `BJ_MAGNETIC_TYPE_FERROMAGNETIC` |
| `classify_magnetic` | fn | `bj_classify_magnetic` |
| `Material` | struct | `BjMaterial` |
| `impl Deserialize for Material` | trait | `BjMaterial_from_json` |
| `impl Serialize for Material` | trait | `BjMaterial_to_json` |
| `Material::eps_r` | field | `BjMaterial.eps_r` |
| `Material::mu_r` | field | `BjMaterial.mu_r` |
| `Material::conductivity` | field | `BjMaterial.conductivity` |
| `Material::new` | fn | `BjMaterial_new` |
| `Material::vacuum` | fn | `BjMaterial_vacuum` |
| `Material::dielectric` | fn | `BjMaterial_dielectric` |
| `Material::conductor` | fn | `BjMaterial_conductor` |
| `Material::refractive_index` | fn | `BjMaterial::refractive_index` |
| `Material::permittivity` | fn | `BjMaterial::permittivity` |
| `Material::permeability` | fn | `BjMaterial::permeability` |
| `Material::impedance` | fn | `BjMaterial::impedance` |
| `Material::lossless_optical_constants` | fn | `_bj_material_lossless_optical_constants` |
| `Material::wave_speed` | fn | `BjMaterial::wave_speed` |
| `Material::is_lossy` | fn | `BjMaterial::is_lossy` |
| `Material::try_reflectance_at` | fn | `BjMaterial::try_reflectance_at` |
| `Material::reflectance_at` | fn | deprecated in 1.2.0, not carried (ADR 0002 row 1): try_reflectance_at |
| `impl Default for Material` | trait | `BjMaterial_default` |

## `maxwell` (10)

| Rust | Kind | Cyrius |
|---|---|---|
| `gauss_electric_divergence` | fn | `bj_gauss_electric_divergence` |
| `gauss_magnetic_divergence` | fn | `bj_gauss_magnetic_divergence` |
| `displacement_current` | fn | `bj_displacement_current` |
| `wave_speed` | fn | `bj_wave_speed` |
| `refractive_index` | fn | `bj_refractive_index` |
| `impedance` | fn | `bj_impedance` |
| `free_space_impedance` | fn | `bj_free_space_impedance` |
| `wavelength` | fn | `bj_wavelength` |
| `frequency` | fn | `bj_frequency` |
| `skin_depth` | fn | `bj_skin_depth` |

## `polarization` (74)

| Rust | Kind | Cyrius |
|---|---|---|
| `impl ComplexExt for Complex` | trait | BjComplex's own methods norm_sq, norm and arg (src/complex.cyr): hisab 1.4 made them inherent, and the trait only delegates |
| `complex_real` | fn | `bj_complex_real` |
| `complex_zero` | fn | `bj_complex_zero` |
| `complex_from_polar` | fn | `bj_complex_from_polar` |
| `JonesVector` | struct | `BjJonesVector` |
| `impl Deserialize for JonesVector` | trait | `BjJonesVector_from_json` |
| `impl Serialize for JonesVector` | trait | `BjJonesVector_to_json` |
| `JonesVector::x` | field | `BjJonesVector.x` |
| `JonesVector::y` | field | `BjJonesVector.y` |
| `impl Default for JonesVector` | trait | `BjJonesVector_default` |
| `JonesVector::new` | fn | `BjJonesVector_new` |
| `JonesVector::horizontal` | fn | `BjJonesVector_horizontal` |
| `JonesVector::vertical` | fn | `BjJonesVector_vertical` |
| `JonesVector::diagonal` | fn | `BjJonesVector_diagonal` |
| `JonesVector::anti_diagonal` | fn | `BjJonesVector_anti_diagonal` |
| `JonesVector::right_circular` | fn | `BjJonesVector_right_circular` |
| `JonesVector::left_circular` | fn | `BjJonesVector_left_circular` |
| `JonesVector::linear` | fn | `BjJonesVector_linear` |
| `JonesVector::intensity` | fn | `BjJonesVector::intensity` |
| `JonesVector::normalized` | fn | `BjJonesVector::normalized` |
| `JonesVector::inner` | fn | `BjJonesVector::inner` |
| `JonesVector::overlap` | fn | `BjJonesVector::overlap` |
| `JonesMatrix` | struct | `BjJonesMatrix` |
| `impl Deserialize for JonesMatrix` | trait | `BjJonesMatrix_from_json` |
| `impl Serialize for JonesMatrix` | trait | `BjJonesMatrix_to_json` |
| `JonesMatrix::m` | field | `BjJonesMatrix.m00, m01, m10, m11` |
| `JonesMatrix::new` | fn | `BjJonesMatrix_new` |
| `JonesMatrix::identity` | fn | `BjJonesMatrix_identity` |
| `JonesMatrix::apply` | fn | `BjJonesMatrix::apply` |
| `JonesMatrix::compose` | fn | `BjJonesMatrix::compose` |
| `JonesMatrix::horizontal_polarizer` | fn | `BjJonesMatrix_horizontal_polarizer` |
| `JonesMatrix::vertical_polarizer` | fn | `BjJonesMatrix_vertical_polarizer` |
| `JonesMatrix::linear_polarizer` | fn | `BjJonesMatrix_linear_polarizer` |
| `JonesMatrix::half_wave_plate` | fn | `BjJonesMatrix_half_wave_plate` |
| `JonesMatrix::quarter_wave_plate` | fn | `BjJonesMatrix_quarter_wave_plate` |
| `JonesMatrix::waveplate` | fn | `BjJonesMatrix_waveplate` |
| `JonesMatrix::rotator` | fn | `BjJonesMatrix_rotator` |
| `JonesMatrix::determinant` | fn | `BjJonesMatrix::determinant` |
| `StokesVector` | struct | `BjStokesVector` |
| `impl Deserialize for StokesVector` | trait | `BjStokesVector_from_json` |
| `impl Serialize for StokesVector` | trait | `BjStokesVector_to_json` |
| `StokesVector::s` | field | `BjStokesVector.s` |
| `impl Default for StokesVector` | trait | `BjStokesVector_default` |
| `StokesVector::new` | fn | `BjStokesVector_new` |
| `StokesVector::from_jones` | fn | `BjStokesVector_from_jones` |
| `StokesVector::unpolarized` | fn | `BjStokesVector_unpolarized` |
| `StokesVector::horizontal` | fn | `BjStokesVector_horizontal` |
| `StokesVector::vertical` | fn | `BjStokesVector_vertical` |
| `StokesVector::right_circular` | fn | `BjStokesVector_right_circular` |
| `StokesVector::left_circular` | fn | `BjStokesVector_left_circular` |
| `StokesVector::intensity` | fn | `BjStokesVector::intensity` |
| `StokesVector::degree_of_polarization` | fn | `BjStokesVector::degree_of_polarization` |
| `StokesVector::ellipticity_angle` | fn | `BjStokesVector::ellipticity_angle` |
| `StokesVector::orientation_angle` | fn | `BjStokesVector::orientation_angle` |
| `impl Add for StokesVector` | trait | `BjStokesVector_add` |
| `MuellerMatrix` | struct | `BjMuellerMatrix` |
| `impl Deserialize for MuellerMatrix` | trait | `BjMuellerMatrix_from_json` |
| `impl Serialize for MuellerMatrix` | trait | `BjMuellerMatrix_to_json` |
| `MuellerMatrix::m` | field | `BjMuellerMatrix.m` |
| `MuellerMatrix::new` | fn | `BjMuellerMatrix_new` |
| `MuellerMatrix::identity` | fn | `BjMuellerMatrix_identity` |
| `MuellerMatrix::apply` | fn | `BjMuellerMatrix::apply` |
| `MuellerMatrix::compose` | fn | `BjMuellerMatrix::compose` |
| `MuellerMatrix::horizontal_polarizer` | fn | `BjMuellerMatrix_horizontal_polarizer` |
| `MuellerMatrix::vertical_polarizer` | fn | `BjMuellerMatrix_vertical_polarizer` |
| `MuellerMatrix::linear_polarizer` | fn | `BjMuellerMatrix_linear_polarizer` |
| `MuellerMatrix::waveplate` | fn | `BjMuellerMatrix_waveplate` |
| `MuellerMatrix::quarter_wave_plate` | fn | `BjMuellerMatrix_quarter_wave_plate` |
| `MuellerMatrix::half_wave_plate` | fn | `BjMuellerMatrix_half_wave_plate` |
| `MuellerMatrix::rotation` | fn | `BjMuellerMatrix_rotation` |
| `MuellerMatrix::depolarizer` | fn | `BjMuellerMatrix_depolarizer` |
| `MuellerMatrix::attenuator` | fn | `BjMuellerMatrix_attenuator` |
| `MuellerMatrix::right_circular_polarizer` | fn | `BjMuellerMatrix_right_circular_polarizer` |
| `MuellerMatrix::left_circular_polarizer` | fn | `BjMuellerMatrix_left_circular_polarizer` |

## `relativity` (37)

| Rust | Kind | Cyrius |
|---|---|---|
| `lorentz_factor` | fn | `bj_lorentz_factor` |
| `beta` | fn | `bj_beta` |
| `EmTensor` | struct | `BjEmTensor` |
| `EmTensor::components` | field | `BjEmTensor.components` |
| `EmTensor::from_fields` | fn | `BjEmTensor_from_fields` |
| `EmTensor::electric_field` | fn | `BjEmTensor::electric_field` |
| `EmTensor::magnetic_field` | fn | `BjEmTensor::magnetic_field` |
| `EmTensor::first_invariant` | fn | `BjEmTensor::first_invariant` |
| `EmTensor::second_invariant` | fn | `BjEmTensor::second_invariant` |
| `EmTensor::dual` | fn | `BjEmTensor::dual` |
| `EmTensor::four_force` | fn | `BjEmTensor::four_force` |
| `lorentz_transform_fields` | fn | `bj_lorentz_transform_fields` |
| `lorentz_transform_x` | fn | `bj_lorentz_transform_x` |
| `FourVector` | struct | `BjFourVector` |
| `FourVector::t` | field | `BjFourVector.t` |
| `FourVector::x` | field | `BjFourVector.x` |
| `FourVector::y` | field | `BjFourVector.y` |
| `FourVector::z` | field | `BjFourVector.z` |
| `FourVector::new` | fn | `BjFourVector_new` |
| `FourVector::from_time_position` | fn | `BjFourVector_from_time_position` |
| `FourVector::dot` | fn | `BjFourVector::dot` |
| `FourVector::interval_sq` | fn | `BjFourVector::interval_sq` |
| `FourVector::spatial_magnitude` | fn | `BjFourVector::spatial_magnitude` |
| `FourVector::boost_x` | fn | `BjFourVector::boost_x` |
| `FourPotential` | struct | `BjFourPotential` |
| `FourPotential::phi_over_c` | field | `BjFourPotential.phi_over_c` |
| `FourPotential::a` | field | `BjFourPotential.a` |
| `retarded_scalar_potential` | fn | `bj_retarded_scalar_potential` |
| `retarded_vector_potential` | fn | `bj_retarded_vector_potential` |
| `lienard_wiechert_scalar_potential` | fn | `bj_lienard_wiechert_scalar_potential` |
| `lienard_wiechert_vector_potential` | fn | `bj_lienard_wiechert_vector_potential` |
| `lienard_wiechert_e` | fn | `bj_lienard_wiechert_e` |
| `lienard_wiechert_b` | fn | `bj_lienard_wiechert_b` |
| `larmor_power` | fn | `bj_larmor_power` |
| `relativistic_larmor_power` | fn | `bj_relativistic_larmor_power` |
| `solve_retarded_time` | fn | `bj_solve_retarded_time` |
| `lienard_wiechert_fields` | fn | `bj_lienard_wiechert_fields` |

## `rf` (50)

| Rust | Kind | Cyrius |
|---|---|---|
| `characteristic_impedance_lossless` | fn | `bj_characteristic_impedance_lossless` |
| `characteristic_impedance_lossy` | fn | `bj_characteristic_impedance_lossy` |
| `propagation_constant` | fn | `bj_propagation_constant` |
| `phase_velocity` | fn | `bj_phase_velocity` |
| `input_impedance_lossless_complex` | fn | `bj_input_impedance_lossless_complex` |
| `input_impedance_lossless` | fn | deprecated in 1.2.0, not carried (ADR 0002 row 1): input_impedance_lossless_complex |
| `reflection_coefficient` | fn | `bj_reflection_coefficient` |
| `reflection_coefficient_complex` | fn | `bj_reflection_coefficient_complex` |
| `vswr` | fn | `bj_vswr` |
| `return_loss` | fn | `bj_return_loss` |
| `insertion_loss` | fn | `bj_insertion_loss` |
| `mismatch_loss` | fn | `bj_mismatch_loss` |
| `normalize_impedance` | fn | `bj_normalize_impedance` |
| `denormalize_impedance` | fn | `bj_denormalize_impedance` |
| `gamma_from_z` | fn | `bj_gamma_from_z` |
| `z_from_gamma` | fn | `bj_z_from_gamma` |
| `quarter_wave_transformer` | fn | `bj_quarter_wave_transformer` |
| `single_stub_match` | fn | `bj_single_stub_match` |
| `single_stub_match_complex` | fn | `bj_single_stub_match_complex` |
| `SMatrix` | struct | `BjSMatrix` |
| `impl Deserialize for SMatrix` | trait | `BjSMatrix_from_json` |
| `impl Serialize for SMatrix` | trait | `BjSMatrix_to_json` |
| `SMatrix::ports` | field | `BjSMatrix.ports` |
| `SMatrix::s` | field | `BjSMatrix.s` |
| `SMatrix::z0` | field | `BjSMatrix.z0` |
| `SMatrix::new` | fn | `BjSMatrix_new` |
| `SMatrix::two_port` | fn | `BjSMatrix_two_port` |
| `SMatrix::get` | fn | `BjSMatrix::get` |
| `SMatrix::set` | fn | `BjSMatrix::set` |
| `SMatrix::cascade` | fn | `BjSMatrix::cascade` |
| `s_to_z` | fn | `bj_s_to_z` |
| `TouchstoneFormat` | enum | `BjTouchstoneFormat` |
| `impl Deserialize for TouchstoneFormat` | trait | `BjTouchstoneFormat_from_json` |
| `impl Serialize for TouchstoneFormat` | trait | `BjTouchstoneFormat_to_json` |
| `TouchstoneFormat::MagnitudeAngle` | variant | `BJ_TOUCHSTONE_FORMAT_MAGNITUDE_ANGLE` |
| `TouchstoneFormat::RealImag` | variant | `BJ_TOUCHSTONE_FORMAT_REAL_IMAG` |
| `TouchstoneFormat::DbAngle` | variant | `BJ_TOUCHSTONE_FORMAT_DB_ANGLE` |
| `TouchstoneData` | struct | `BjTouchstoneData` |
| `impl Deserialize for TouchstoneData` | trait | `BjTouchstoneData_from_json` |
| `impl Serialize for TouchstoneData` | trait | `BjTouchstoneData_to_json` |
| `TouchstoneData::ports` | field | `BjTouchstoneData.ports` |
| `TouchstoneData::z0` | field | `BjTouchstoneData.z0` |
| `TouchstoneData::frequencies` | field | `BjTouchstoneData.frequencies` |
| `TouchstoneData::s_matrices` | field | `BjTouchstoneData.s_matrices` |
| `parse_touchstone` | fn | `bj_parse_touchstone` |
| `try_write_touchstone` | fn | `bj_try_write_touchstone` |
| `write_touchstone` | fn | deprecated in 1.2.0, not carried (ADR 0002 row 1): try_write_touchstone |
| `linear_array_factor` | fn | `bj_linear_array_factor` |
| `planar_array_factor` | fn | `bj_planar_array_factor` |
| `circular_array_factor` | fn | `bj_circular_array_factor` |

## `scattering` (17)

| Rust | Kind | Cyrius |
|---|---|---|
| `MIE_MAX_SIZE_PARAMETER` | const | `BJ_MIE_MAX_SIZE_PARAMETER` |
| `MIE_MAX_MX` | const | `BJ_MIE_MAX_MX` |
| `MieResult` | struct | `BjMieResult` |
| `impl Deserialize for MieResult` | trait | `BjMieResult_from_json` |
| `impl Serialize for MieResult` | trait | `BjMieResult_to_json` |
| `MieResult::q_sca` | field | `BjMieResult.q_sca` |
| `MieResult::q_ext` | field | `BjMieResult.q_ext` |
| `MieResult::q_abs` | field | `BjMieResult.q_abs` |
| `MieResult::g` | field | `BjMieResult.g` |
| `MieResult::a` | field | `BjMieResult.a` |
| `MieResult::b` | field | `BjMieResult.b` |
| `mie` | fn | `bj_mie` |
| `cross_section_from_efficiency` | fn | `bj_cross_section_from_efficiency` |
| `rayleigh_cross_section` | fn | `bj_rayleigh_cross_section` |
| `rayleigh_efficiency` | fn | `bj_rayleigh_efficiency` |
| `rayleigh_phase_function` | fn | `bj_rayleigh_phase_function` |
| `rayleigh_angular_intensity` | fn | `bj_rayleigh_angular_intensity` |

## `wave` (42)

| Rust | Kind | Cyrius |
|---|---|---|
| `poynting_vector` | fn | `bj_poynting_vector` |
| `plane_wave_intensity` | fn | `bj_plane_wave_intensity` |
| `radiation_pressure_absorbed` | fn | `bj_radiation_pressure_absorbed` |
| `radiation_pressure_reflected` | fn | `bj_radiation_pressure_reflected` |
| `e_from_b` | fn | `bj_e_from_b` |
| `b_from_e` | fn | `bj_b_from_e` |
| `plane_wave_e` | fn | `bj_plane_wave_e` |
| `wave_number` | fn | `bj_wave_number` |
| `angular_frequency` | fn | `bj_angular_frequency` |
| `momentum_density` | fn | `bj_momentum_density` |
| `snell_refraction_angle` | fn | `bj_snell_refraction_angle` |
| `critical_angle` | fn | `bj_critical_angle` |
| `brewster_angle` | fn | `bj_brewster_angle` |
| `fresnel_rs` | fn | `bj_fresnel_rs` |
| `fresnel_rp` | fn | `bj_fresnel_rp` |
| `reflectance_s` | fn | `bj_reflectance_s` |
| `reflectance_p` | fn | `bj_reflectance_p` |
| `reflectance_normal` | fn | `bj_reflectance_normal` |
| `transmittance_normal` | fn | `bj_transmittance_normal` |
| `snell_cos_theta_t` | fn | `bj_snell_cos_theta_t` |
| `fresnel_rs_direct` | fn | `bj_fresnel_rs_direct` |
| `fresnel_rp_direct` | fn | `bj_fresnel_rp_direct` |
| `reflectance_s_direct` | fn | `bj_reflectance_s_direct` |
| `reflectance_p_direct` | fn | `bj_reflectance_p_direct` |
| `reflectance_unpolarized` | fn | `bj_reflectance_unpolarized` |
| `transmittance_unpolarized` | fn | `bj_transmittance_unpolarized` |
| `schlick_reflectance` | fn | `bj_schlick_reflectance` |
| `reflectance_at_interface` | fn | `bj_reflectance_at_interface` |
| `reflectance_normal_materials` | fn | `bj_reflectance_normal_materials` |
| `material_reflection_coefficient` | fn | `bj_material_reflection_coefficient` |
| `rectangular_waveguide_cutoff` | fn | `bj_rectangular_waveguide_cutoff` |
| `rectangular_waveguide_cutoff_wavelength` | fn | `bj_rectangular_waveguide_cutoff_wavelength` |
| `guide_wavelength` | fn | `bj_guide_wavelength` |
| `cylindrical_waveguide_cutoff` | fn | `bj_cylindrical_waveguide_cutoff` |
| `hertzian_dipole_pattern` | fn | `bj_hertzian_dipole_pattern` |
| `half_wave_dipole_pattern` | fn | `bj_half_wave_dipole_pattern` |
| `hertzian_dipole_directivity` | fn | `bj_hertzian_dipole_directivity` |
| `half_wave_dipole_directivity` | fn | `bj_half_wave_dipole_directivity` |
| `hertzian_dipole_radiation_resistance` | fn | `bj_hertzian_dipole_radiation_resistance` |
| `half_wave_dipole_radiation_resistance` | fn | `bj_half_wave_dipole_radiation_resistance` |
| `effective_aperture` | fn | `bj_effective_aperture` |
| `friis_transmission` | fn | `bj_friis_transmission` |
