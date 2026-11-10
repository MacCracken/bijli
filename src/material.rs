//! Dielectric and magnetic materials — polarization, bound charges, magnetization.
//!
//! All SI units: C/m², A/m, F/m, H/m.

use crate::error::{BijliError, Result};
use crate::field::{EPSILON_0, FieldVector, MU_0};

// ── Dielectrics ────────────────────────────────────────────────────

/// Electric susceptibility from relative permittivity: χ_e = ε_r − 1.
#[inline]
#[must_use]
pub fn electric_susceptibility(relative_permittivity: f64) -> f64 {
    relative_permittivity - 1.0
}

/// Relative permittivity from electric susceptibility: ε_r = 1 + χ_e.
#[inline]
#[must_use]
pub fn relative_permittivity(electric_susceptibility: f64) -> f64 {
    1.0 + electric_susceptibility
}

/// Absolute permittivity: ε = ε_r × ε₀.
#[inline]
#[must_use]
pub fn absolute_permittivity(relative_permittivity: f64) -> f64 {
    relative_permittivity * EPSILON_0
}

/// Polarization vector: **P** = ε₀ χ_e **E** = (ε_r − 1)ε₀ **E**.
#[inline]
#[must_use]
pub fn polarization(relative_permittivity: f64, electric_field: &FieldVector) -> FieldVector {
    electric_field.scale((relative_permittivity - 1.0) * EPSILON_0)
}

/// Displacement field: **D** = ε₀**E** + **P** = ε**E**.
#[inline]
#[must_use]
pub fn displacement_field(relative_permittivity: f64, electric_field: &FieldVector) -> FieldVector {
    electric_field.scale(relative_permittivity * EPSILON_0)
}

/// Bound surface charge density: σ_b = **P** ⋅ **n̂**.
///
/// `polarization` is the polarization vector, `normal` is the outward surface normal.
#[inline]
#[must_use]
pub fn bound_surface_charge(polarization: &FieldVector, normal: &FieldVector) -> f64 {
    polarization.dot(normal)
}

/// Bound volume charge density: ρ_b = −∇⋅**P**.
///
/// For uniform polarization, this is zero. This function computes it from
/// a finite-difference divergence given ∂Px/∂x, ∂Py/∂y, ∂Pz/∂z.
#[inline]
#[must_use]
pub fn bound_volume_charge(dp_dx: f64, dp_dy: f64, dp_dz: f64) -> f64 {
    -(dp_dx + dp_dy + dp_dz)
}

/// Energy density in a dielectric: u = ε E²/2.
#[inline]
#[must_use]
pub fn dielectric_energy_density(relative_permittivity: f64, e_magnitude: f64) -> f64 {
    0.5 * relative_permittivity * EPSILON_0 * e_magnitude * e_magnitude
}

/// Clausius-Mossotti relation: (ε_r − 1)/(ε_r + 2) = Nα/(3ε₀).
///
/// Given molecular polarizability α (SI, C·m²/V = F·m²) and number density
/// N (m⁻³), returns ε_r = (1 + 2X)/(1 − X) with X = Nα/(3ε₀).
///
/// # Errors
/// [`BijliError::InvalidParameter`] for non-finite inputs, N < 0, or
/// X ≥ 1 (the relation diverges at X = 1 — the "polarization catastrophe" —
/// and gives a negative, unphysical ε_r beyond it);
/// [`BijliError::InvalidPermittivity`] if X ≤ −1/2 (ε_r ≤ 0, only possible
/// for negative α).
#[inline]
pub fn clausius_mossotti(polarizability: f64, number_density: f64) -> Result<f64> {
    if !polarizability.is_finite() || !number_density.is_finite() || number_density < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "polarizability and number density must be finite with N ≥ 0, got α = {polarizability}, N = {number_density}"
            ),
        });
    }
    let chi = number_density * polarizability / (3.0 * EPSILON_0);
    if chi.is_nan() || chi >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "Nα/(3ε₀) = {chi} ≥ 1: Clausius-Mossotti diverges (polarization catastrophe), no physical ε_r"
            ),
        });
    }
    let eps_r = (1.0 + 2.0 * chi) / (1.0 - chi);
    if eps_r.is_nan() || eps_r <= 0.0 {
        return Err(BijliError::InvalidPermittivity { value: eps_r });
    }
    Ok(eps_r)
}

// ── Magnetic materials ─────────────────────────────────────────────

/// Magnetic susceptibility from relative permeability: χ_m = μ_r − 1.
#[inline]
#[must_use]
pub fn magnetic_susceptibility(relative_permeability: f64) -> f64 {
    relative_permeability - 1.0
}

/// Relative permeability from magnetic susceptibility: μ_r = 1 + χ_m.
#[inline]
#[must_use]
pub fn relative_permeability(magnetic_susceptibility: f64) -> f64 {
    1.0 + magnetic_susceptibility
}

/// Absolute permeability: μ = μ_r × μ₀.
#[inline]
#[must_use]
pub fn absolute_permeability(relative_permeability: f64) -> f64 {
    relative_permeability * MU_0
}

/// Magnetization: **M** = χ_m **H**.
///
/// `h_field` is the auxiliary magnetic field **H** (A/m).
#[inline]
#[must_use]
pub fn magnetization(relative_permeability: f64, h_field: &FieldVector) -> FieldVector {
    h_field.scale(relative_permeability - 1.0)
}

/// Auxiliary field **H** from **B**: **H** = **B**/μ₀ − **M** = **B**/(μ_r μ₀).
#[inline]
pub fn h_field_from_b(relative_permeability: f64, b_field: &FieldVector) -> Result<FieldVector> {
    if relative_permeability.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "relative permeability cannot be zero".into(),
        });
    }
    Ok(b_field.scale(1.0 / (relative_permeability * MU_0)))
}

/// **B** from **H**: **B** = μ_r μ₀ **H**.
#[inline]
#[must_use]
pub fn b_field_from_h(relative_permeability: f64, h_field: &FieldVector) -> FieldVector {
    h_field.scale(relative_permeability * MU_0)
}

/// Magnetic energy density in a material: u = B²/(2μ).
#[inline]
pub fn magnetic_energy_density_material(
    relative_permeability: f64,
    b_magnitude: f64,
) -> Result<f64> {
    let mu = relative_permeability * MU_0;
    if mu.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "permeability cannot be zero for energy density".into(),
        });
    }
    Ok(b_magnitude * b_magnitude / (2.0 * mu))
}

/// Bound surface current density: **K**_b = **M** × **n̂** (A/m).
#[inline]
#[must_use]
pub fn bound_surface_current(magnetization: &FieldVector, normal: &FieldVector) -> FieldVector {
    magnetization.cross(normal)
}

/// Curie's law for paramagnetic susceptibility: χ_m = C/T.
///
/// `curie_constant` is the material's Curie constant (K),
/// `temperature` is in kelvin.
#[inline]
pub fn curie_law(curie_constant: f64, temperature: f64) -> Result<f64> {
    if temperature <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("temperature must be positive, got {temperature} K"),
        });
    }
    Ok(curie_constant / temperature)
}

/// Curie-Weiss law for ferromagnetic susceptibility: χ_m = C/(T − T_c).
///
/// Valid only for T > T_c (above Curie temperature).
#[inline]
pub fn curie_weiss(curie_constant: f64, temperature: f64, curie_temperature: f64) -> Result<f64> {
    let denom = temperature - curie_temperature;
    if denom <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "temperature ({temperature} K) must be above Curie temperature ({curie_temperature} K)"
            ),
        });
    }
    Ok(curie_constant / denom)
}

/// Classify magnetic material type from susceptibility.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum MagneticType {
    Diamagnetic,
    Paramagnetic,
    Ferromagnetic,
}

/// Classify a material's magnetic behavior from its susceptibility.
#[inline]
#[must_use]
pub fn classify_magnetic(susceptibility: f64) -> MagneticType {
    if susceptibility < -1e-15 {
        MagneticType::Diamagnetic
    } else if susceptibility > 1.0 {
        MagneticType::Ferromagnetic
    } else {
        MagneticType::Paramagnetic
    }
}

// ── Unified Material struct ────────────────────────────────────────

/// Unified electromagnetic material properties.
///
/// Combines dielectric, magnetic, and conductive properties in one struct,
/// bridging the `material`, `fdtd`, and `wave` modules.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Material {
    /// Relative permittivity ε_r.
    pub eps_r: f64,
    /// Relative permeability μ_r.
    pub mu_r: f64,
    /// Conductivity σ (S/m). Zero for lossless.
    pub conductivity: f64,
}

impl Material {
    /// Create a material with given properties.
    #[inline]
    #[must_use]
    pub fn new(eps_r: f64, mu_r: f64, conductivity: f64) -> Self {
        Self {
            eps_r,
            mu_r,
            conductivity,
        }
    }

    /// Vacuum.
    #[inline]
    #[must_use]
    pub fn vacuum() -> Self {
        Self {
            eps_r: 1.0,
            mu_r: 1.0,
            conductivity: 0.0,
        }
    }

    /// Simple dielectric (non-magnetic, lossless).
    #[inline]
    #[must_use]
    pub fn dielectric(eps_r: f64) -> Self {
        Self {
            eps_r,
            mu_r: 1.0,
            conductivity: 0.0,
        }
    }

    /// Conductor with given conductivity.
    #[inline]
    #[must_use]
    pub fn conductor(conductivity: f64) -> Self {
        Self {
            eps_r: 1.0,
            mu_r: 1.0,
            conductivity,
        }
    }

    /// Refractive index n = √(ε_r μ_r) (lossless approximation).
    #[inline]
    #[must_use]
    pub fn refractive_index(&self) -> f64 {
        (self.eps_r * self.mu_r).sqrt()
    }

    /// Absolute permittivity ε = ε_r ε₀.
    #[inline]
    #[must_use]
    pub fn permittivity(&self) -> f64 {
        self.eps_r * EPSILON_0
    }

    /// Absolute permeability μ = μ_r μ₀.
    #[inline]
    #[must_use]
    pub fn permeability(&self) -> f64 {
        self.mu_r * MU_0
    }

    /// Lossless wave impedance η = √(μ/ε) (Ω). Conductivity is ignored.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] unless ε_r is positive and finite;
    /// [`BijliError::InvalidPermeability`] unless μ_r is positive and finite.
    #[inline]
    pub fn impedance(&self) -> Result<f64> {
        if !(self.eps_r.is_finite() && self.eps_r > 0.0) {
            return Err(BijliError::InvalidPermittivity { value: self.eps_r });
        }
        if !(self.mu_r.is_finite() && self.mu_r > 0.0) {
            return Err(BijliError::InvalidPermeability { value: self.mu_r });
        }
        Ok((self.mu_r * MU_0 / (self.eps_r * EPSILON_0)).sqrt())
    }

    /// Refractive index n = √(ε_r μ_r) and relative impedance η/η₀ = √(μ_r/ε_r)
    /// for **lossless** optics (Fresnel coefficients).
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] / [`BijliError::InvalidPermeability`]
    /// unless ε_r / μ_r are positive and finite;
    /// [`BijliError::InvalidParameter`] if σ ≠ 0 — a conductor needs the
    /// frequency-dependent complex index n² = ε_r μ_r − iσμ_r/(ωε₀), which a
    /// frequency-free `Material` cannot provide.
    #[inline]
    pub(crate) fn lossless_optical_constants(&self) -> Result<(f64, f64)> {
        if !(self.eps_r.is_finite() && self.eps_r > 0.0) {
            return Err(BijliError::InvalidPermittivity { value: self.eps_r });
        }
        if !(self.mu_r.is_finite() && self.mu_r > 0.0) {
            return Err(BijliError::InvalidPermeability { value: self.mu_r });
        }
        if self.conductivity != 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "Fresnel coefficients from Material are lossless-only; conductivity σ = {} S/m \
                     needs a frequency-dependent complex refractive index",
                    self.conductivity
                ),
            });
        }
        Ok((
            (self.eps_r * self.mu_r).sqrt(),
            (self.mu_r / self.eps_r).sqrt(),
        ))
    }

    /// Wave speed v = 1/√(εμ).
    #[inline]
    pub fn wave_speed(&self) -> Result<f64> {
        if !(self.eps_r > 0.0 && self.mu_r > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: "ε_r and μ_r must be positive".into(),
            });
        }
        Ok(1.0 / (self.eps_r * EPSILON_0 * self.mu_r * MU_0).sqrt())
    }

    /// Whether this material is lossy (has conductivity).
    #[inline]
    #[must_use]
    pub fn is_lossy(&self) -> bool {
        self.conductivity > 0.0
    }

    /// Normal-incidence reflectance at the interface with another lossless
    /// material.
    ///
    /// R = ((η₂ − η₁)/(η₂ + η₁))², η = √(μ_r/ε_r) (relative wave impedance).
    /// For non-magnetic media this reduces to ((n₁ − n₂)/(n₁ + n₂))²; for
    /// magnetic media it is the correct form (e.g. ε_r = μ_r = 4 is
    /// impedance-matched to vacuum: R = 0, although n = 4).
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] / [`BijliError::InvalidPermeability`]
    /// for non-positive or non-finite ε_r / μ_r; [`BijliError::InvalidParameter`]
    /// if either material has σ ≠ 0 (lossless only — a conductor's
    /// reflectance depends on frequency).
    #[inline]
    pub fn try_reflectance_at(&self, other: &Self) -> Result<f64> {
        let (_, eta1) = self.lossless_optical_constants()?;
        let (_, eta2) = other.lossless_optical_constants()?;
        let r = (eta2 - eta1) / (eta2 + eta1);
        Ok(r * r)
    }

    /// Normal-incidence reflectance at the interface with another material.
    ///
    /// Same as [`try_reflectance_at`](Self::try_reflectance_at), but returns
    /// `NaN` where that returns an error (conductive or invalid materials).
    #[deprecated(
        since = "1.2.0",
        note = "cannot report conductive (σ ≠ 0) or invalid materials; use `try_reflectance_at`, which returns Result"
    )]
    #[inline]
    #[must_use]
    pub fn reflectance_at(&self, other: &Self) -> f64 {
        self.try_reflectance_at(other).unwrap_or(f64::NAN)
    }
}

impl Default for Material {
    /// Default is vacuum.
    #[inline]
    fn default() -> Self {
        Self::vacuum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dielectric tests ───────────────────────────────────────────

    #[test]
    fn test_electric_susceptibility() {
        // Glass: ε_r ≈ 4.0 → χ_e = 3.0
        assert!((electric_susceptibility(4.0) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_relative_permittivity() {
        assert!((relative_permittivity(3.0) - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_absolute_permittivity() {
        let eps = absolute_permittivity(4.0);
        assert!((eps - 4.0 * EPSILON_0).abs() < 1e-25);
    }

    #[test]
    fn test_polarization_vacuum() {
        // Vacuum: ε_r = 1 → P = 0
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let p = polarization(1.0, &e);
        assert!(p.magnitude() < 1e-25);
    }

    #[test]
    fn test_polarization_dielectric() {
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let p = polarization(4.0, &e);
        // P = 3ε₀ × 1000
        assert!((p.x - 3.0 * EPSILON_0 * 1000.0).abs() < 1e-15);
    }

    #[test]
    fn test_displacement_field() {
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let d = displacement_field(4.0, &e);
        assert!((d.x - 4.0 * EPSILON_0 * 1000.0).abs() < 1e-15);
    }

    #[test]
    fn test_bound_surface_charge() {
        let p = FieldVector::new(1e-6, 0.0, 0.0);
        let n = FieldVector::new(1.0, 0.0, 0.0);
        assert!((bound_surface_charge(&p, &n) - 1e-6).abs() < 1e-20);
    }

    #[test]
    fn test_bound_volume_charge_uniform() {
        // Uniform polarization → zero bound volume charge
        assert!(bound_volume_charge(0.0, 0.0, 0.0).abs() < 1e-30);
    }

    #[test]
    fn test_dielectric_energy_density() {
        let u = dielectric_energy_density(4.0, 1000.0);
        assert!((u - 0.5 * 4.0 * EPSILON_0 * 1e6).abs() < 1e-10);
    }

    #[test]
    fn test_clausius_mossotti() {
        // Small polarizability → ε_r ≈ 1
        let eps_r = clausius_mossotti(1e-40, 1e28).unwrap();
        assert!(eps_r > 1.0);
        assert!(eps_r < 2.0);
    }

    // ── Magnetic material tests ────────────────────────────────────

    #[test]
    fn test_magnetic_susceptibility() {
        assert!((magnetic_susceptibility(1.001) - 0.001).abs() < 1e-15);
    }

    #[test]
    fn test_relative_permeability() {
        assert!((relative_permeability(0.001) - 1.001).abs() < 1e-15);
    }

    #[test]
    fn test_absolute_permeability() {
        let mu = absolute_permeability(1000.0);
        assert!((mu - 1000.0 * MU_0).abs() < 1e-15);
    }

    #[test]
    fn test_magnetization() {
        let h = FieldVector::new(1000.0, 0.0, 0.0);
        let m = magnetization(1000.0, &h);
        // M = (μ_r - 1)H = 999 × 1000
        assert!((m.x - 999_000.0).abs() < 1e-6);
    }

    #[test]
    fn test_h_b_roundtrip() {
        let mu_r = 1000.0;
        let h = FieldVector::new(100.0, 0.0, 0.0);
        let b = b_field_from_h(mu_r, &h);
        let h_back = h_field_from_b(mu_r, &b).unwrap();
        assert!((h_back.x - h.x).abs() < 1e-10);
    }

    #[test]
    fn test_magnetic_energy_density() {
        let u = magnetic_energy_density_material(1.0, 1.0).unwrap();
        // u = B²/(2μ₀) — same as free space
        let expected = 1.0 / (2.0 * MU_0);
        assert!((u - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_bound_surface_current() {
        let m = FieldVector::new(1000.0, 0.0, 0.0);
        let n = FieldVector::new(0.0, 1.0, 0.0);
        let k = bound_surface_current(&m, &n);
        // M × n̂ = (1000,0,0) × (0,1,0) = (0,0,1000)
        assert!((k.z - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_curie_law() {
        let chi = curie_law(1.0, 300.0).unwrap();
        assert!((chi - 1.0 / 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_curie_law_zero_temp() {
        assert!(curie_law(1.0, 0.0).is_err());
    }

    #[test]
    fn test_curie_weiss() {
        let chi = curie_weiss(1.0, 1000.0, 770.0).unwrap();
        assert!((chi - 1.0 / 230.0).abs() < 1e-10);
    }

    #[test]
    fn test_curie_weiss_below_tc() {
        assert!(curie_weiss(1.0, 700.0, 770.0).is_err());
    }

    #[test]
    fn test_classify_diamagnetic() {
        assert_eq!(classify_magnetic(-1e-5), MagneticType::Diamagnetic);
    }

    #[test]
    fn test_classify_paramagnetic() {
        assert_eq!(classify_magnetic(1e-3), MagneticType::Paramagnetic);
    }

    #[test]
    fn test_classify_ferromagnetic() {
        assert_eq!(classify_magnetic(1000.0), MagneticType::Ferromagnetic);
    }

    // ── Unified Material tests ────────────────────────────────────

    #[test]
    fn test_material_vacuum() {
        let m = Material::vacuum();
        assert!((m.eps_r - 1.0).abs() < 1e-15);
        assert!((m.mu_r - 1.0).abs() < 1e-15);
        assert!(m.conductivity.abs() < 1e-30);
        assert!(!m.is_lossy());
    }

    #[test]
    fn test_material_dielectric() {
        let m = Material::dielectric(4.0);
        assert!((m.refractive_index() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_material_conductor() {
        let m = Material::conductor(5.96e7); // copper
        assert!(m.is_lossy());
    }

    #[test]
    fn test_material_wave_speed_vacuum() {
        let v = Material::vacuum().wave_speed().unwrap();
        assert!((v - crate::field::SPEED_OF_LIGHT).abs() / crate::field::SPEED_OF_LIGHT < 1e-6);
    }

    #[test]
    fn test_material_impedance_vacuum() {
        let z = Material::vacuum().impedance().unwrap();
        // Z₀ ≈ 376.73 Ω
        assert!((z - 376.73).abs() < 0.1);
    }

    #[test]
    #[allow(deprecated)]
    fn test_material_reflectance_same() {
        let m = Material::dielectric(2.25);
        assert!(m.reflectance_at(&m) < 1e-10); // no reflection at same material
    }

    #[test]
    #[allow(deprecated)]
    fn test_material_reflectance_air_glass() {
        let air = Material::vacuum();
        let glass = Material::dielectric(2.25); // n = 1.5
        let r = air.reflectance_at(&glass);
        // R = ((1-1.5)/(1+1.5))² = 0.04
        assert!((r - 0.04).abs() < 1e-6);
    }

    #[test]
    fn test_material_default_is_vacuum() {
        assert_eq!(Material::default(), Material::vacuum());
    }

    #[test]
    fn test_material_permittivity() {
        let m = Material::dielectric(4.0);
        assert!((m.permittivity() - 4.0 * EPSILON_0).abs() < 1e-25);
    }

    #[test]
    fn test_material_permeability() {
        let m = Material::new(1.0, 1000.0, 0.0);
        assert!((m.permeability() - 1000.0 * MU_0).abs() < 1e-15);
    }

    // ── Regression tests ──────────────────────────────────────────

    #[test]
    fn test_try_reflectance_magnetic_impedance_matched() {
        // ε_r = μ_r = 4 → η = η₀: no reflection from vacuum at normal
        // incidence (old n-only formula gave ((1−4)/(1+4))² = 0.36).
        let m = Material::new(4.0, 4.0, 0.0);
        let r = Material::vacuum().try_reflectance_at(&m).unwrap();
        assert!(r.abs() < 1e-15, "R = {r}");
        // Pure magnetic medium μ_r = 2.25, ε_r = 1: η = 1.5 η₀ → R = (0.5/2.5)² = 0.04
        let mag = Material::new(1.0, 2.25, 0.0);
        let r = Material::vacuum().try_reflectance_at(&mag).unwrap();
        assert!((r - 0.04).abs() < 1e-15);
    }

    #[test]
    fn test_try_reflectance_dielectric_matches_n_formula() {
        let r = Material::vacuum()
            .try_reflectance_at(&Material::dielectric(2.25))
            .unwrap();
        assert!((r - 0.04).abs() < 1e-15);
    }

    #[test]
    fn test_try_reflectance_rejects_conductor_and_invalid() {
        // Copper: old code returned R = 0 (σ ignored).
        let cu = Material::conductor(5.96e7);
        assert!(Material::vacuum().try_reflectance_at(&cu).is_err());
        assert!(cu.try_reflectance_at(&Material::vacuum()).is_err());
        assert!(matches!(
            Material::vacuum().try_reflectance_at(&Material::dielectric(-2.0)),
            Err(BijliError::InvalidPermittivity { .. })
        ));
        assert!(matches!(
            Material::vacuum().try_reflectance_at(&Material::new(1.0, -1.0, 0.0)),
            Err(BijliError::InvalidPermeability { .. })
        ));
        assert!(
            Material::vacuum()
                .try_reflectance_at(&Material::dielectric(f64::NAN))
                .is_err()
        );
    }

    #[test]
    #[allow(deprecated)]
    fn test_reflectance_at_deprecated_never_wrong() {
        let m = Material::new(4.0, 4.0, 0.0);
        assert!(Material::vacuum().reflectance_at(&m).abs() < 1e-15);
        assert!(
            Material::vacuum()
                .reflectance_at(&Material::conductor(5.96e7))
                .is_nan()
        );
    }

    #[test]
    fn test_impedance_validates_permeability() {
        // Old code: μ_r = −1 → NaN, μ_r = 0 → 0 Ω.
        assert!(matches!(
            Material::new(1.0, -1.0, 0.0).impedance(),
            Err(BijliError::InvalidPermeability { .. })
        ));
        assert!(matches!(
            Material::new(1.0, 0.0, 0.0).impedance(),
            Err(BijliError::InvalidPermeability { .. })
        ));
        assert!(matches!(
            Material::dielectric(f64::NAN).impedance(),
            Err(BijliError::InvalidPermittivity { .. })
        ));
        // η = η₀ √(μ_r/ε_r): μ_r = 4, ε_r = 1 → 2η₀
        let z0 = Material::vacuum().impedance().unwrap();
        let z = Material::new(1.0, 4.0, 0.0).impedance().unwrap();
        assert!((z / z0 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_wave_speed_nan_rejected() {
        assert!(Material::dielectric(f64::NAN).wave_speed().is_err());
        assert!(Material::new(1.0, f64::NAN, 0.0).wave_speed().is_err());
    }

    #[test]
    fn test_clausius_mossotti_closed_form() {
        // X = Nα/(3ε₀) = 0.5 → ε_r = (1 + 1)/(1 − 0.5) = 4
        let n = 1e28;
        let alpha = 1.5 * EPSILON_0 / n;
        let eps = clausius_mossotti(alpha, n).unwrap();
        assert!((eps - 4.0).abs() < 1e-12, "{eps}");
    }

    #[test]
    fn test_clausius_mossotti_catastrophe_rejected() {
        // X = 1.2: old code returned ε_r = (1 + 2.4)/(1 − 1.2) = −17 as Ok.
        let n = 1e28;
        assert!(clausius_mossotti(3.6 * EPSILON_0 / n, n).is_err());
        assert!(clausius_mossotti(3.0 * EPSILON_0 * (1.0 + 1e-9) / n, n).is_err()); // X ≈ 1⁺
        assert!(clausius_mossotti(f64::NAN, n).is_err());
        assert!(clausius_mossotti(1e-40, -1.0).is_err());
        // X = −0.6 ≤ −1/2 → ε_r < 0
        assert!(matches!(
            clausius_mossotti(-1.8 * EPSILON_0 / n, n),
            Err(BijliError::InvalidPermittivity { .. })
        ));
    }
}
