//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into bijli electromagnetism parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! ushma   (thermodynamics)   ──┐
//! prakash (optics)             ┼──> bridge ──> bijli EM parameters
//! tanmatra (atomic physics)    ┤
//! dravya  (material science)  ┘
//! ```

// ── Ushma bridges (thermodynamics) ─────────────────────────────────────────

/// Joule heating power per unit volume (W/m³) from current density and
/// resistivity.
///
/// p = J² ρ = J² / σ, with `current_density_a_m2` the current density
/// magnitude J (A/m²) and `resistivity_ohm_m` the bulk resistivity ρ (Ω·m).
/// The result can be fed to [`em_energy_to_heating_rate`].
#[must_use]
#[inline]
pub fn joule_heating_power_density(current_density_a_m2: f64, resistivity_ohm_m: f64) -> f64 {
    current_density_a_m2 * current_density_a_m2 * resistivity_ohm_m
}

/// Temperature rise rate (K/s) from an absorbed EM power density (W/m³),
/// assuming adiabatic heating (no conduction, convection or radiation).
///
/// dT/dt = p / (ρ c_p), where p is the volumetric heating rate (W/m³, e.g.
/// from [`joule_heating_power_density`] or ½σ|E|²), ρ the material mass
/// density (kg/m³) and c_p the specific heat capacity (J/(kg·K)).
///
/// Returns 0.0 if ρ·c_p is not positive (including NaN).
#[must_use]
#[inline]
pub fn em_energy_to_heating_rate(
    power_density_w_m3: f64,
    material_density_kg_m3: f64,
    specific_heat_j_per_kg_k: f64,
) -> f64 {
    let denom = material_density_kg_m3 * specific_heat_j_per_kg_k;
    if denom > 0.0 {
        power_density_w_m3 / denom
    } else {
        0.0
    }
}

// ── Prakash bridges (optics) ───────────────────────────────────────────────

/// Refractive index of Schott N-BK7 glass from vacuum wavelength (nm), via
/// the Sellmeier equation (λ in µm):
///
/// n² = 1 + B₁λ²/(λ²−C₁) + B₂λ²/(λ²−C₂) + B₃λ²/(λ²−C₃)
///
/// with B = (1.03961212, 0.231792344, 1.01046945) and
/// C = (0.00600069867, 0.0200179144, 103.560653) µm² (n_d = 1.5168 at
/// 587.56 nm).
///
/// **Validity:** the fit is specified for 300–2500 nm. Outside that range
/// the value is an extrapolation and should not be trusted (the formula has
/// poles near 77 nm, 141 nm and 10.2 µm). As this function cannot return an
/// error, non-positive (or NaN) wavelengths and wavelengths where n² ≤ 0
/// return the sentinel 1.0; callers needing range checking must validate
/// the wavelength themselves.
#[must_use]
pub fn wavelength_to_refractive_index_bk7(wavelength_nm: f64) -> f64 {
    if wavelength_nm.is_nan() || wavelength_nm <= 0.0 {
        return 1.0;
    }
    let l2 = (wavelength_nm * 1e-3) * (wavelength_nm * 1e-3); // μm²
    // BK7 Sellmeier coefficients
    let b1 = 1.039_612_12;
    let b2 = 0.231_792_344;
    let b3 = 1.010_469_45;
    let c1 = 0.006_000_698_67;
    let c2 = 0.020_017_914_4;
    let c3 = 103.560_653;

    let n2 = 1.0 + b1 * l2 / (l2 - c1) + b2 * l2 / (l2 - c2) + b3 * l2 / (l2 - c3);
    if n2 > 0.0 { n2.sqrt() } else { 1.0 }
}

/// Convert E-field amplitude (V/m) to optical intensity (W/m²).
///
/// I = ε₀ × c × E² / 2 (time-averaged, plane wave in vacuum; E is the peak
/// amplitude).
#[must_use]
#[inline]
pub fn e_field_to_optical_intensity(e_amplitude_v_m: f64) -> f64 {
    // CODATA 2018 (kept local: `bridge` compiles without the `field` feature).
    const EPSILON_0: f64 = 8.854_187_812_8e-12;
    const C: f64 = 299_792_458.0;
    0.5 * EPSILON_0 * C * e_amplitude_v_m * e_amplitude_v_m
}

// ── Tanmatra bridges (atomic/nuclear physics) ──────────────────────────────

/// Convert electron energy level transition (eV) to photon emission
/// wavelength (nm, vacuum).
///
/// λ = hc / E with hc = 1239.841984… eV·nm (exact SI h = 6.62607015e-34 J·s,
/// c = 299 792 458 m/s, e = 1.602176634e-19 C). Returns 0.0 for
/// non-positive (or NaN) energies.
#[must_use]
#[inline]
pub fn energy_level_to_wavelength_nm(energy_ev: f64) -> f64 {
    /// hc in eV·nm.
    const HC_EV_NM: f64 = 1_239.841_984_332_002_5;
    if energy_ev.is_nan() || energy_ev <= 0.0 {
        return 0.0;
    }
    HC_EV_NM / energy_ev
}

/// Convert nuclear charge number (Z) to Coulomb potential at distance r (V).
///
/// V = k × Z × e / r, where k = 1/(4πε₀) ≈ 8.9876e9 N·m²/C², e = 1.602e-19 C.
#[must_use]
#[inline]
pub fn nuclear_charge_potential(atomic_number: u32, distance_m: f64) -> f64 {
    if distance_m.is_nan() || distance_m <= 0.0 {
        return 0.0;
    }
    // CODATA 2018 (kept local: `bridge` compiles without the `field` feature).
    const K: f64 = 1.0 / (4.0 * std::f64::consts::PI * 8.854_187_812_8e-12);
    const E: f64 = 1.602_176_634e-19;
    K * atomic_number as f64 * E / distance_m
}

// ── Dravya bridges (material science) ──────────────────────────────────────

/// Piezoelectric blocking (actuation) stress (Pa) produced by an electric
/// field, for a single piezoelectric stress coefficient.
///
/// σ = e × E, with `piezo_coeff_c_m2` the stress constant e (C/m², e.g.
/// e₃₃) and `e_field_v_m` the field E (V/m) along the corresponding axis.
///
/// **Sign convention:** the result is the stress the material exerts on a
/// rigid clamp, positive when the field drives expansion (e > 0, E along
/// the poling direction). In the IEEE Std 176 constitutive form
/// T = c^E S − e E, the stress *inside* a clamped crystal (S = 0) is
/// T = −e E, i.e. the negative of this value; negate the result if you need
/// the IEEE T.
#[must_use]
#[inline]
pub fn e_field_to_piezo_stress(e_field_v_m: f64, piezo_coeff_c_m2: f64) -> f64 {
    e_field_v_m * piezo_coeff_c_m2
}

/// Magnetostrictive strain from magnetic flux density, quadratic-moment model.
///
/// ε = (3/2) λ_s (M/M_s)², with M/M_s approximated by B/B_sat and clamped to
/// |B/B_sat| ≤ 1 (the strain saturates once the material is saturated;
/// B > B_sat only adds μ₀H).
///
/// λ_s is the isotropic (polycrystalline) saturation magnetostriction
/// constant in λ(θ) = (3/2) λ_s (cos²θ − 1/3). The (3/2) prefactor measures
/// strain relative to a demagnetized state whose domains are perpendicular
/// to the measuring axis (e.g. under compressive prestress, as in
/// Jiles-type and Terfenol-D models); relative to a randomly oriented
/// demagnetized state the saturation strain is λ_s, so multiply by 2/3 for
/// that reference. This is not a single-crystal (cubic λ₁₀₀/λ₁₁₁) formula.
///
/// - `saturation_magnetostriction`: λ_s, dimensionless and signed (e.g.
///   ≈ −7e-6 for polycrystalline iron, ≈ −33e-6 for nickel; negative means
///   contraction along the field).
/// - `saturation_flux_t`: saturation flux density B_sat (e.g. ≈ 2.15 T for
///   iron).
///
/// Returns 0.0 if `saturation_flux_t` is not positive (including NaN).
#[must_use]
#[inline]
pub fn magnetic_flux_to_magnetostrictive_strain(
    flux_density_t: f64,
    saturation_magnetostriction: f64,
    saturation_flux_t: f64,
) -> f64 {
    if saturation_flux_t.is_nan() || saturation_flux_t <= 0.0 {
        return 0.0;
    }
    let ratio = (flux_density_t / saturation_flux_t).clamp(-1.0, 1.0);
    1.5 * saturation_magnetostriction * ratio * ratio
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Ushma ──────────────────────────────────────────────────────────

    #[test]
    fn joule_heating_basic() {
        // 1 A/m², 1e-7 Ω·m (copper) → 1e-7 W/m³
        let p = joule_heating_power_density(1.0, 1e-7);
        assert!((p - 1e-7).abs() < 1e-12);
    }

    #[test]
    fn heating_rate_basic() {
        let rate = em_energy_to_heating_rate(1000.0, 1000.0, 4186.0);
        assert!((rate - 1000.0 / (1000.0 * 4186.0)).abs() < 1e-10);
    }

    #[test]
    fn heating_rate_zero_density() {
        assert_eq!(em_energy_to_heating_rate(1000.0, 0.0, 4186.0), 0.0);
        assert_eq!(em_energy_to_heating_rate(1000.0, f64::NAN, 4186.0), 0.0);
    }

    #[test]
    fn joule_heating_chain_copper() {
        // J = 1e6 A/m² in copper (ρ = 1.68e-8 Ω·m) → p = 1.68e4 W/m³;
        // ρ_m = 8960 kg/m³, c_p = 385 J/(kg·K) → dT/dt = 4.87013e-3 K/s.
        let p = joule_heating_power_density(1e6, 1.68e-8);
        assert!((p - 1.68e4).abs() < 1e-9);
        let rate = em_energy_to_heating_rate(p, 8960.0, 385.0);
        assert!((rate - 0.004_870_129_870_129_87).abs() < 1e-15);
    }

    // ── Prakash ────────────────────────────────────────────────────────

    #[test]
    fn bk7_visible_range() {
        let n = wavelength_to_refractive_index_bk7(550.0);
        assert!(n > 1.5 && n < 1.55, "BK7 at 550nm: got {n}");
    }

    #[test]
    fn bk7_dispersion() {
        let n_blue = wavelength_to_refractive_index_bk7(450.0);
        let n_red = wavelength_to_refractive_index_bk7(650.0);
        assert!(n_blue > n_red, "blue should have higher n than red");
    }

    #[test]
    fn bk7_zero_wavelength() {
        assert_eq!(wavelength_to_refractive_index_bk7(0.0), 1.0);
        assert_eq!(wavelength_to_refractive_index_bk7(f64::NAN), 1.0);
    }

    #[test]
    fn bk7_d_line() {
        // Schott N-BK7 catalog n_d = 1.51680 at 587.5618 nm.
        let n = wavelength_to_refractive_index_bk7(587.5618);
        assert!((n - 1.5168).abs() < 1e-5, "got {n}");
    }

    #[test]
    fn e_field_intensity_roundtrip() {
        let e = 100.0;
        let i = e_field_to_optical_intensity(e);
        assert!(i > 0.0);
        // I = ε₀cE²/2 ≈ 8.854e-12 * 3e8 * 10000 / 2 ≈ 13.28 W/m²
        assert!((i - 13.28).abs() < 0.1);
    }

    // ── Tanmatra ───────────────────────────────────────────────────────

    #[test]
    fn hydrogen_lyman_alpha() {
        // Lyman-α: 10.2 eV → ~121.5 nm
        let nm = energy_level_to_wavelength_nm(10.2);
        assert!((nm - 121.6).abs() < 1.0);
    }

    #[test]
    fn energy_zero() {
        assert_eq!(energy_level_to_wavelength_nm(0.0), 0.0);
        assert_eq!(energy_level_to_wavelength_nm(f64::NAN), 0.0);
    }

    #[test]
    fn energy_one_ev_is_hc() {
        // hc/e = 6.62607015e-34 × 299792458 / 1.602176634e-19 m = 1239.8419843320025 nm
        let nm = energy_level_to_wavelength_nm(1.0);
        assert!((nm - 1_239.841_984_332_002_5).abs() < 1e-9);
    }

    #[test]
    fn nuclear_potential_hydrogen() {
        // Z=1 at 1 Bohr radius (5.29e-11 m) ≈ 27.2 V
        let v = nuclear_charge_potential(1, 5.29e-11);
        assert!((v - 27.2).abs() < 0.5);
    }

    #[test]
    fn nuclear_potential_zero_distance() {
        assert_eq!(nuclear_charge_potential(1, 0.0), 0.0);
    }

    // ── Dravya ─────────────────────────────────────────────────────────

    #[test]
    fn piezo_stress_basic() {
        // PZT: e ≈ 15 C/m², E = 1000 V/m → σ = 15000 Pa
        let s = e_field_to_piezo_stress(1000.0, 15.0);
        assert!((s - 15000.0).abs() < 0.1);
    }

    #[test]
    fn magnetostriction_iron() {
        // Polycrystalline iron: λ_s ≈ −7e-6, B_sat ≈ 2.15 T, B = 1.0 T →
        // ε = 1.5 × (−7e-6) × (1/2.15)² = −2.2714981e-6 (contraction).
        let e = magnetic_flux_to_magnetostrictive_strain(1.0, -7e-6, 2.15);
        assert!((e - (-2.271_498_107_084_910_6e-6)).abs() < 1e-18);
    }

    #[test]
    fn magnetostriction_saturates_beyond_bsat() {
        // Beyond saturation the strain stays at (3/2)λ_s (was (3/2)λ_s(B/B_sat)²).
        let lambda_s = -7e-6;
        let at_sat = 1.5 * lambda_s;
        for b in [2.15, 3.0, 10.0, -3.0] {
            let e = magnetic_flux_to_magnetostrictive_strain(b, lambda_s, 2.15);
            assert!((e - at_sat).abs() < 1e-18, "B = {b}: {e}");
        }
    }

    #[test]
    fn magnetostriction_zero_saturation() {
        assert_eq!(
            magnetic_flux_to_magnetostrictive_strain(1.0, 35e-6, 0.0),
            0.0
        );
    }
}
