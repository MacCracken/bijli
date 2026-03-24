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

/// Clausius-Mossotti relation: (ε_r − 1)/(ε_r + 2) = nα/(3ε₀).
///
/// Given molecular polarizability α and number density n, returns ε_r.
#[inline]
pub fn clausius_mossotti(polarizability: f64, number_density: f64) -> Result<f64> {
    let chi = number_density * polarizability / (3.0 * EPSILON_0);
    let denom = 1.0 - chi;
    if denom.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Clausius-Mossotti divergence (ferroelectric transition)".into(),
        });
    }
    Ok((1.0 + 2.0 * chi) / denom)
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
}
