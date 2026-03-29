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

/// Convert current density (A/m²), resistance per unit length (Ω/m), and
/// conductor cross-section (m²) to Joule heating power per unit volume (W/m³).
///
/// P/V = J² × ρ_resistivity = J² / σ
#[must_use]
#[inline]
pub fn joule_heating_power_density(current_density_a_m2: f64, resistivity_ohm_m: f64) -> f64 {
    current_density_a_m2 * current_density_a_m2 * resistivity_ohm_m
}

/// Convert EM field energy density (J/m³) to temperature rise rate (K/s)
/// given material properties.
///
/// dT/dt = u / (ρ × c_p) where u is volumetric energy deposition rate.
#[must_use]
#[inline]
pub fn em_energy_to_heating_rate(
    energy_density_j_m3: f64,
    material_density_kg_m3: f64,
    specific_heat_j_per_kg_k: f64,
) -> f64 {
    let denom = material_density_kg_m3 * specific_heat_j_per_kg_k;
    if denom <= 0.0 {
        return 0.0;
    }
    energy_density_j_m3 / denom
}

// ── Prakash bridges (optics) ───────────────────────────────────────────────

/// Convert wavelength (nm) to refractive index using the Sellmeier equation
/// for BK7 glass.
///
/// n² = 1 + B₁λ²/(λ²-C₁) + B₂λ²/(λ²-C₂) + B₃λ²/(λ²-C₃)
#[must_use]
pub fn wavelength_to_refractive_index_bk7(wavelength_nm: f64) -> f64 {
    if wavelength_nm <= 0.0 {
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
/// I = ε₀ × c × E² / 2
#[must_use]
#[inline]
pub fn e_field_to_optical_intensity(e_amplitude_v_m: f64) -> f64 {
    const EPSILON_0: f64 = 8.854_187_817e-12;
    const C: f64 = 299_792_458.0;
    0.5 * EPSILON_0 * C * e_amplitude_v_m * e_amplitude_v_m
}

// ── Tanmatra bridges (atomic/nuclear physics) ──────────────────────────────

/// Convert electron energy level transition (eV) to photon emission
/// wavelength (nm).
///
/// λ = hc / E, where h = 4.1357e-15 eV·s, c = 3e8 m/s.
#[must_use]
#[inline]
pub fn energy_level_to_wavelength_nm(energy_ev: f64) -> f64 {
    if energy_ev <= 0.0 {
        return 0.0;
    }
    // hc in eV·nm = 1239.8419843320025
    1239.842 / energy_ev
}

/// Convert nuclear charge number (Z) to Coulomb potential at distance r (V).
///
/// V = k × Z × e / r, where k = 8.9876e9, e = 1.602e-19 C.
#[must_use]
#[inline]
pub fn nuclear_charge_potential(atomic_number: u32, distance_m: f64) -> f64 {
    if distance_m <= 0.0 {
        return 0.0;
    }
    const K: f64 = 8.987_551_792e9;
    const E: f64 = 1.602_176_634e-19;
    K * atomic_number as f64 * E / distance_m
}

// ── Dravya bridges (material science) ──────────────────────────────────────

/// Convert electric field strength (V/m) to piezoelectric stress (Pa)
/// for a given piezoelectric coefficient.
///
/// σ = e × E, where e is the piezoelectric stress constant (C/m²).
#[must_use]
#[inline]
pub fn e_field_to_piezo_stress(e_field_v_m: f64, piezo_coeff_c_m2: f64) -> f64 {
    e_field_v_m * piezo_coeff_c_m2
}

/// Convert magnetic flux density (T) to magnetostrictive strain for
/// a given magnetostrictive coefficient.
///
/// ε = λ_s × (3/2) × (B/B_sat)² for cubic symmetry.
/// `saturation_magnetostriction`: λ_s (e.g. ~35e-6 for iron).
/// `saturation_flux_t`: saturation flux density (e.g. ~2.15 T for iron).
#[must_use]
#[inline]
pub fn magnetic_flux_to_magnetostrictive_strain(
    flux_density_t: f64,
    saturation_magnetostriction: f64,
    saturation_flux_t: f64,
) -> f64 {
    if saturation_flux_t <= 0.0 {
        return 0.0;
    }
    let ratio = flux_density_t / saturation_flux_t;
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
        // Iron: λ_s ≈ 35e-6, B_sat ≈ 2.15 T, B = 1.0 T
        let e = magnetic_flux_to_magnetostrictive_strain(1.0, 35e-6, 2.15);
        let expected = 1.5 * 35e-6 * (1.0 / 2.15_f64).powi(2);
        assert!((e - expected).abs() < 1e-10);
    }

    #[test]
    fn magnetostriction_zero_saturation() {
        assert_eq!(
            magnetic_flux_to_magnetostrictive_strain(1.0, 35e-6, 0.0),
            0.0
        );
    }
}
