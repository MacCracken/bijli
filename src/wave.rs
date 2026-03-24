//! Electromagnetic wave propagation, polarization, Poynting vector.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};
use crate::field::{FieldVector, EPSILON_0, MU_0, SPEED_OF_LIGHT};

/// Poynting vector S = E × B / μ₀ (W/m²).
///
/// Represents energy flux of the electromagnetic field.
#[must_use]
pub fn poynting_vector(e: &FieldVector, b: &FieldVector) -> FieldVector {
    e.cross(b).scale(1.0 / MU_0)
}

/// Intensity (time-averaged power per unit area) of a plane wave.
///
/// I = c⋅ε₀⋅E₀²/2 (W/m²)
#[inline]
#[must_use]
pub fn plane_wave_intensity(e_amplitude: f64) -> f64 {
    0.5 * SPEED_OF_LIGHT * EPSILON_0 * e_amplitude * e_amplitude
}

/// Radiation pressure for a perfectly absorbed wave.
///
/// P = I/c (Pa)
#[inline]
#[must_use]
pub fn radiation_pressure_absorbed(intensity: f64) -> f64 {
    intensity / SPEED_OF_LIGHT
}

/// Radiation pressure for a perfectly reflected wave.
///
/// P = 2I/c (Pa)
#[inline]
#[must_use]
pub fn radiation_pressure_reflected(intensity: f64) -> f64 {
    2.0 * intensity / SPEED_OF_LIGHT
}

/// E-field amplitude from B-field amplitude: E₀ = cB₀.
#[inline]
#[must_use]
pub fn e_from_b(b_amplitude: f64) -> f64 {
    SPEED_OF_LIGHT * b_amplitude
}

/// B-field amplitude from E-field amplitude: B₀ = E₀/c.
#[inline]
#[must_use]
pub fn b_from_e(e_amplitude: f64) -> f64 {
    e_amplitude / SPEED_OF_LIGHT
}

/// Plane wave E-field at position and time.
///
/// E(x,t) = E₀ sin(kx - ωt + φ) (for wave propagating in +x)
#[inline]
#[must_use]
pub fn plane_wave_e(e_amplitude: f64, k: f64, x: f64, omega: f64, t: f64, phase: f64) -> f64 {
    e_amplitude * (k * x - omega * t + phase).sin()
}

/// Wave number from frequency: k = ω/v = 2πf/v.
pub fn wave_number(frequency: f64, velocity: f64) -> Result<f64> {
    if velocity.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "velocity cannot be zero for wave number".into(),
        });
    }
    Ok(2.0 * std::f64::consts::PI * frequency / velocity)
}

/// Angular frequency from frequency: ω = 2πf.
#[inline]
#[must_use]
pub fn angular_frequency(frequency: f64) -> f64 {
    2.0 * std::f64::consts::PI * frequency
}

/// Electromagnetic momentum density: g = S/c² (kg/(m²⋅s)).
#[must_use]
pub fn momentum_density(e: &FieldVector, b: &FieldVector) -> FieldVector {
    poynting_vector(e, b).scale(1.0 / (SPEED_OF_LIGHT * SPEED_OF_LIGHT))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poynting_vector() {
        let e = FieldVector::new(0.0, 1000.0, 0.0); // E in y
        let b = FieldVector::new(0.0, 0.0, 1e-6); // B in z
        let s = poynting_vector(&e, &b);
        // S = E × B / μ₀ should be in +x direction
        assert!(s.x > 0.0);
    }

    #[test]
    fn test_plane_wave_intensity() {
        let i = plane_wave_intensity(1000.0);
        assert!(i > 0.0);
        // I = 0.5 * c * ε₀ * E₀²
        let expected = 0.5 * SPEED_OF_LIGHT * EPSILON_0 * 1e6;
        assert!((i - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_radiation_pressure() {
        let i = 1000.0; // 1 kW/m²
        let p_abs = radiation_pressure_absorbed(i);
        let p_ref = radiation_pressure_reflected(i);
        assert!((p_ref - 2.0 * p_abs).abs() < 1e-15);
    }

    #[test]
    fn test_e_b_relationship() {
        let e0 = 1000.0;
        let b0 = b_from_e(e0);
        let e_back = e_from_b(b0);
        assert!((e_back - e0).abs() < 1e-6);
    }

    #[test]
    fn test_plane_wave() {
        // At t=0, x=0, phase=0: sin(0) = 0
        let e = plane_wave_e(100.0, 1.0, 0.0, 1.0, 0.0, 0.0);
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn test_wave_number() {
        // Visible light: f ≈ 5e14 Hz → k ≈ 1.05e7 rad/m
        let k = wave_number(5e14, SPEED_OF_LIGHT).unwrap();
        assert!((k - 1.047e7).abs() / 1.047e7 < 0.01);
    }

    #[test]
    fn test_angular_frequency() {
        let omega = angular_frequency(60.0);
        assert!((omega - 120.0 * std::f64::consts::PI).abs() < 1e-10);
    }
}
