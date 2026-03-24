//! Maxwell's equations — divergence, curl, wave equation.
//!
//! Provides finite-difference approximations for Maxwell's equations
//! on a uniform 3D grid.

use crate::error::{BijliError, Result};
use crate::field::{EPSILON_0, FieldVector, MU_0};

/// Gauss's law for electricity: ∇⋅E = ρ/ε₀.
///
/// Given charge density ρ (C/m³), returns the required divergence of E.
#[inline]
#[must_use]
pub fn gauss_electric_divergence(charge_density: f64) -> f64 {
    charge_density / EPSILON_0
}

/// Gauss's law for magnetism: ∇⋅B = 0.
///
/// Returns 0.0 (no magnetic monopoles).
#[inline]
#[must_use]
pub fn gauss_magnetic_divergence() -> f64 {
    0.0
}

/// Displacement current density: J_d = ε₀ ∂E/∂t.
///
/// Given the time derivative of E (V/m/s), returns the displacement current (A/m²).
#[inline]
#[must_use]
pub fn displacement_current(de_dt: &FieldVector) -> FieldVector {
    de_dt.scale(EPSILON_0)
}

/// EM wave speed in a medium: v = 1/√(εμ).
#[inline]
pub fn wave_speed(permittivity: f64, permeability: f64) -> Result<f64> {
    if permittivity <= 0.0 {
        return Err(BijliError::InvalidPermittivity {
            value: permittivity,
        });
    }
    if permeability <= 0.0 {
        return Err(BijliError::InvalidPermeability {
            value: permeability,
        });
    }
    Ok(1.0 / (permittivity * permeability).sqrt())
}

/// Refractive index: n = c/v = √(ε_r × μ_r).
///
/// # Errors
///
/// Returns [`BijliError::InvalidPermittivity`] or [`BijliError::InvalidPermeability`]
/// if the relative values are not positive.
#[inline]
pub fn refractive_index(relative_permittivity: f64, relative_permeability: f64) -> Result<f64> {
    if relative_permittivity <= 0.0 {
        return Err(BijliError::InvalidPermittivity {
            value: relative_permittivity,
        });
    }
    if relative_permeability <= 0.0 {
        return Err(BijliError::InvalidPermeability {
            value: relative_permeability,
        });
    }
    Ok((relative_permittivity * relative_permeability).sqrt())
}

/// Impedance of a medium: η = √(μ/ε) (ohms).
#[inline]
pub fn impedance(permittivity: f64, permeability: f64) -> Result<f64> {
    if permittivity <= 0.0 {
        return Err(BijliError::InvalidPermittivity {
            value: permittivity,
        });
    }
    if permeability <= 0.0 {
        return Err(BijliError::InvalidPermeability {
            value: permeability,
        });
    }
    Ok((permeability / permittivity).sqrt())
}

/// Impedance of free space: η₀ = √(μ₀/ε₀) ≈ 376.73 Ω.
#[inline]
#[must_use]
pub fn free_space_impedance() -> f64 {
    (MU_0 / EPSILON_0).sqrt()
}

/// Wavelength from frequency: λ = v/f.
#[inline]
pub fn wavelength(frequency: f64, velocity: f64) -> Result<f64> {
    if frequency <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("frequency must be positive, got {frequency}"),
        });
    }
    if velocity <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("velocity must be positive, got {velocity}"),
        });
    }
    Ok(velocity / frequency)
}

/// Frequency from wavelength: f = v/λ.
#[inline]
pub fn frequency(wavelength_m: f64, velocity: f64) -> Result<f64> {
    if wavelength_m <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive, got {wavelength_m}"),
        });
    }
    if velocity <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("velocity must be positive, got {velocity}"),
        });
    }
    Ok(velocity / wavelength_m)
}

/// Skin depth: δ = √(2/(ωμσ)) (meters).
///
/// How far EM waves penetrate into a conductor.
#[inline]
pub fn skin_depth(angular_freq: f64, permeability: f64, conductivity: f64) -> Result<f64> {
    let denom = angular_freq * permeability * conductivity;
    if denom <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: "all parameters must be positive for skin depth".into(),
        });
    }
    Ok((2.0 / denom).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::SPEED_OF_LIGHT;

    #[test]
    fn test_gauss_electric() {
        let div = gauss_electric_divergence(1e-6);
        // ρ/ε₀ should be large
        assert!(div > 1e4);
    }

    #[test]
    fn test_gauss_magnetic() {
        assert_eq!(gauss_magnetic_divergence(), 0.0);
    }

    #[test]
    fn test_wave_speed_vacuum() {
        let v = wave_speed(EPSILON_0, MU_0).unwrap();
        assert!((v - SPEED_OF_LIGHT).abs() / SPEED_OF_LIGHT < 1e-6);
    }

    #[test]
    fn test_refractive_index_vacuum() {
        let n = refractive_index(1.0, 1.0).unwrap();
        assert!((n - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_refractive_index_glass() {
        // Glass: ε_r ≈ 2.25, μ_r ≈ 1.0 → n ≈ 1.5
        let n = refractive_index(2.25, 1.0).unwrap();
        assert!((n - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_refractive_index_negative_permittivity() {
        assert!(refractive_index(-1.0, 1.0).is_err());
    }

    #[test]
    fn test_wavelength_zero_velocity() {
        assert!(wavelength(5e14, 0.0).is_err());
    }

    #[test]
    fn test_frequency_zero_velocity() {
        assert!(frequency(530e-9, 0.0).is_err());
    }

    #[test]
    fn test_free_space_impedance() {
        let z0 = free_space_impedance();
        assert!((z0 - 376.73).abs() < 0.1);
    }

    #[test]
    fn test_wavelength_visible_light() {
        // Green light: f ≈ 5.66e14 Hz → λ ≈ 530 nm
        let lambda = wavelength(5.66e14, SPEED_OF_LIGHT).unwrap();
        assert!((lambda * 1e9 - 530.0).abs() < 1.0);
    }

    #[test]
    fn test_frequency_from_wavelength() {
        let f = frequency(530e-9, SPEED_OF_LIGHT).unwrap();
        assert!((f - 5.66e14).abs() / 5.66e14 < 0.01);
    }

    #[test]
    fn test_skin_depth() {
        // Copper at 60 Hz: σ ≈ 5.96e7, μ ≈ μ₀
        let delta = skin_depth(2.0 * std::f64::consts::PI * 60.0, MU_0, 5.96e7).unwrap();
        // Should be around 8.5 mm
        assert!((delta * 1000.0 - 8.5).abs() < 1.0);
    }

    #[test]
    fn test_displacement_current() {
        let de_dt = FieldVector::new(1e6, 0.0, 0.0);
        let jd = displacement_current(&de_dt);
        assert!((jd.x - EPSILON_0 * 1e6).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_permittivity() {
        assert!(wave_speed(-1.0, MU_0).is_err());
    }

    #[test]
    fn test_invalid_permeability() {
        assert!(impedance(EPSILON_0, -1.0).is_err());
    }
}
