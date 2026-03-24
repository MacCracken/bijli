//! Electromagnetic wave propagation, polarization, Poynting vector.

use crate::error::{BijliError, Result};
use crate::field::{EPSILON_0, FieldVector, MU_0, SPEED_OF_LIGHT};

/// Poynting vector S = E × B / μ₀ (W/m²).
///
/// Represents energy flux of the electromagnetic field.
#[inline]
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
#[inline]
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
#[inline]
#[must_use]
pub fn momentum_density(e: &FieldVector, b: &FieldVector) -> FieldVector {
    poynting_vector(e, b).scale(1.0 / (SPEED_OF_LIGHT * SPEED_OF_LIGHT))
}

// ── Refraction, Reflection, Transmission ───────────────────────────

/// Snell's law: n₁ sin θ₁ = n₂ sin θ₂.
///
/// Returns the refracted angle θ₂ in radians, or `None` for total internal reflection.
#[inline]
pub fn snell_refraction_angle(n1: f64, n2: f64, theta_i: f64) -> Result<Option<f64>> {
    if n2.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "refractive index n₂ cannot be zero".into(),
        });
    }
    let sin_t = n1 * theta_i.sin() / n2;
    if sin_t.abs() > 1.0 {
        Ok(None) // total internal reflection
    } else {
        Ok(Some(sin_t.asin()))
    }
}

/// Critical angle for total internal reflection: θ_c = arcsin(n₂/n₁).
///
/// Only exists when n₁ > n₂.
#[inline]
pub fn critical_angle(n1: f64, n2: f64) -> Result<Option<f64>> {
    if n1.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "refractive index n₁ cannot be zero".into(),
        });
    }
    let ratio = n2 / n1;
    if ratio >= 1.0 {
        Ok(None) // no total internal reflection when n₂ ≥ n₁
    } else {
        Ok(Some(ratio.asin()))
    }
}

/// Brewster's angle: θ_B = arctan(n₂/n₁).
///
/// At this angle, reflected light is completely s-polarized.
#[inline]
pub fn brewster_angle(n1: f64, n2: f64) -> Result<f64> {
    if n1.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "refractive index n₁ cannot be zero".into(),
        });
    }
    Ok((n2 / n1).atan())
}

/// Fresnel reflection coefficient for s-polarization (TE).
///
/// r_s = (n₁ cos θ_i − n₂ cos θ_t) / (n₁ cos θ_i + n₂ cos θ_t)
#[inline]
pub fn fresnel_rs(n1: f64, theta_i: f64, n2: f64, theta_t: f64) -> Result<f64> {
    let num = n1 * theta_i.cos() - n2 * theta_t.cos();
    let den = n1 * theta_i.cos() + n2 * theta_t.cos();
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Fresnel denominator is zero".into(),
        });
    }
    Ok(num / den)
}

/// Fresnel reflection coefficient for p-polarization (TM).
///
/// r_p = (n₂ cos θ_i − n₁ cos θ_t) / (n₂ cos θ_i + n₁ cos θ_t)
#[inline]
pub fn fresnel_rp(n1: f64, theta_i: f64, n2: f64, theta_t: f64) -> Result<f64> {
    let num = n2 * theta_i.cos() - n1 * theta_t.cos();
    let den = n2 * theta_i.cos() + n1 * theta_t.cos();
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Fresnel denominator is zero".into(),
        });
    }
    Ok(num / den)
}

/// Reflectance (power) for s-polarization: R_s = |r_s|².
#[inline]
pub fn reflectance_s(n1: f64, theta_i: f64, n2: f64, theta_t: f64) -> Result<f64> {
    let r = fresnel_rs(n1, theta_i, n2, theta_t)?;
    Ok(r * r)
}

/// Reflectance (power) for p-polarization: R_p = |r_p|².
#[inline]
pub fn reflectance_p(n1: f64, theta_i: f64, n2: f64, theta_t: f64) -> Result<f64> {
    let r = fresnel_rp(n1, theta_i, n2, theta_t)?;
    Ok(r * r)
}

/// Normal-incidence reflectance: R = ((n₁ − n₂)/(n₁ + n₂))².
#[inline]
pub fn reflectance_normal(n1: f64, n2: f64) -> Result<f64> {
    let den = n1 + n2;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "sum of refractive indices cannot be zero".into(),
        });
    }
    let r = (n1 - n2) / den;
    Ok(r * r)
}

/// Normal-incidence transmittance: T = 1 − R.
#[inline]
pub fn transmittance_normal(n1: f64, n2: f64) -> Result<f64> {
    Ok(1.0 - reflectance_normal(n1, n2)?)
}

// ── Waveguides ─────────────────────────────────────────────────────

/// Cutoff frequency of a rectangular waveguide TE_mn or TM_mn mode.
///
/// f_c = (c/2)√((m/a)² + (n/b)²)
///
/// `a` and `b` are the waveguide dimensions (a ≥ b by convention).
#[inline]
pub fn rectangular_waveguide_cutoff(a: f64, b: f64, m: u32, n: u32, velocity: f64) -> Result<f64> {
    if a <= 0.0 || b <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("waveguide dimensions must be positive, got a={a}, b={b}"),
        });
    }
    let ma = m as f64 / a;
    let nb = n as f64 / b;
    Ok(0.5 * velocity * (ma * ma + nb * nb).sqrt())
}

/// Cutoff wavelength of a rectangular waveguide mode.
///
/// λ_c = 2/√((m/a)² + (n/b)²)
#[inline]
pub fn rectangular_waveguide_cutoff_wavelength(a: f64, b: f64, m: u32, n: u32) -> Result<f64> {
    if a <= 0.0 || b <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("waveguide dimensions must be positive, got a={a}, b={b}"),
        });
    }
    let ma = m as f64 / a;
    let nb = n as f64 / b;
    let denom = (ma * ma + nb * nb).sqrt();
    if denom < 1e-30 {
        return Err(BijliError::InvalidParameter {
            reason: "both m and n cannot be zero".into(),
        });
    }
    Ok(2.0 / denom)
}

/// Guide wavelength: λ_g = λ/√(1 − (λ/λ_c)²).
///
/// Only valid when λ < λ_c (above cutoff).
#[inline]
pub fn guide_wavelength(free_space_wavelength: f64, cutoff_wavelength: f64) -> Result<f64> {
    if cutoff_wavelength.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "cutoff wavelength cannot be zero".into(),
        });
    }
    let ratio = free_space_wavelength / cutoff_wavelength;
    let radicand = 1.0 - ratio * ratio;
    if radicand <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("frequency below cutoff: λ/λ_c = {ratio:.4} (must be < 1)"),
        });
    }
    Ok(free_space_wavelength / radicand.sqrt())
}

/// Cutoff frequency of a cylindrical waveguide TE_mn mode.
///
/// f_c = x'_mn c / (2πa) where x'_mn is the n-th root of J'_m.
///
/// `bessel_root` is the appropriate root (e.g., 1.8412 for TE₁₁).
#[inline]
pub fn cylindrical_waveguide_cutoff(radius: f64, bessel_root: f64, velocity: f64) -> Result<f64> {
    if radius <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("waveguide radius must be positive, got {radius}"),
        });
    }
    Ok(bessel_root * velocity / (2.0 * std::f64::consts::PI * radius))
}

// ── Antenna patterns ───────────────────────────────────────────────

/// Radiation intensity of a Hertzian (short) dipole.
///
/// U(θ) ∝ sin²θ. Returns normalized pattern (peak = 1.0).
#[inline]
#[must_use]
pub fn hertzian_dipole_pattern(theta: f64) -> f64 {
    let s = theta.sin();
    s * s
}

/// Radiation intensity of a half-wave dipole antenna.
///
/// U(θ) ∝ [cos(π/2 cos θ)/sin θ]². Returns normalized pattern.
#[inline]
#[must_use]
pub fn half_wave_dipole_pattern(theta: f64) -> f64 {
    let sin_t = theta.sin();
    if sin_t.abs() < 1e-15 {
        return 0.0; // on-axis nulls
    }
    let num = (std::f64::consts::FRAC_PI_2 * theta.cos()).cos();
    let f = num / sin_t;
    f * f
}

/// Directivity of a Hertzian dipole: D = 1.5 (1.76 dBi).
#[inline]
#[must_use]
pub fn hertzian_dipole_directivity() -> f64 {
    1.5
}

/// Directivity of a half-wave dipole: D ≈ 1.64 (2.15 dBi).
#[inline]
#[must_use]
pub fn half_wave_dipole_directivity() -> f64 {
    1.6409 // exact value: 4/Cin(2π) where Cin is the cosine integral
}

/// Radiation resistance of a Hertzian dipole: R_rad = 80π²(dl/λ)².
#[inline]
pub fn hertzian_dipole_radiation_resistance(dipole_length: f64, wavelength: f64) -> Result<f64> {
    if wavelength.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "wavelength cannot be zero".into(),
        });
    }
    let ratio = dipole_length / wavelength;
    Ok(80.0 * std::f64::consts::PI * std::f64::consts::PI * ratio * ratio)
}

/// Radiation resistance of a half-wave dipole: R_rad ≈ 73.1 Ω.
#[inline]
#[must_use]
pub fn half_wave_dipole_radiation_resistance() -> f64 {
    73.1
}

/// Effective aperture of an antenna: A_e = λ²D/(4π).
#[inline]
pub fn effective_aperture(wavelength: f64, directivity: f64) -> Result<f64> {
    if wavelength < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be non-negative, got {wavelength}"),
        });
    }
    Ok(wavelength * wavelength * directivity / (4.0 * std::f64::consts::PI))
}

/// Friis transmission equation: P_r/P_t = G_t G_r (λ/(4πd))².
///
/// Returns the power ratio (dimensionless).
#[inline]
pub fn friis_transmission(
    gain_tx: f64,
    gain_rx: f64,
    wavelength: f64,
    distance: f64,
) -> Result<f64> {
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let ratio = wavelength / (4.0 * std::f64::consts::PI * distance);
    Ok(gain_tx * gain_rx * ratio * ratio)
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

    #[test]
    fn test_wave_number_zero_velocity_fails() {
        assert!(wave_number(5e14, 0.0).is_err());
    }

    #[test]
    fn test_momentum_density() {
        let e = FieldVector::new(0.0, 1000.0, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1e-6);
        let g = momentum_density(&e, &b);
        // g = S/c² — should be in +x direction (same as Poynting vector)
        assert!(g.x > 0.0);
        assert!(g.y.abs() < 1e-30);
        assert!(g.z.abs() < 1e-30);
    }

    #[test]
    fn test_radiation_pressure_absorbed_value() {
        // 1 kW/m² → P = I/c ≈ 3.34e-6 Pa
        let p = radiation_pressure_absorbed(1000.0);
        assert!((p - 1000.0 / SPEED_OF_LIGHT).abs() < 1e-15);
    }

    // ── Refraction / reflection tests ──────────────────────────────

    #[test]
    fn test_snell_normal_incidence() {
        let theta_t = snell_refraction_angle(1.0, 1.5, 0.0).unwrap().unwrap();
        assert!(theta_t.abs() < 1e-10);
    }

    #[test]
    fn test_snell_air_to_glass() {
        // 45° in air → ~28.1° in glass (n=1.5)
        let theta_i = std::f64::consts::FRAC_PI_4;
        let theta_t = snell_refraction_angle(1.0, 1.5, theta_i).unwrap().unwrap();
        let expected = (theta_i.sin() / 1.5).asin();
        assert!((theta_t - expected).abs() < 1e-10);
    }

    #[test]
    fn test_total_internal_reflection() {
        // Glass to air at steep angle
        let theta_i = 1.0; // ~57° > critical angle for n=1.5→1.0
        let result = snell_refraction_angle(1.5, 1.0, theta_i).unwrap();
        assert!(result.is_none()); // TIR
    }

    #[test]
    fn test_critical_angle() {
        // Glass (n=1.5) to air (n=1): θ_c = arcsin(1/1.5) ≈ 41.8°
        let theta_c = critical_angle(1.5, 1.0).unwrap().unwrap();
        assert!((theta_c - (1.0_f64 / 1.5).asin()).abs() < 1e-10);
    }

    #[test]
    fn test_critical_angle_none_for_dense_medium() {
        // Air to glass: no total internal reflection
        assert!(critical_angle(1.0, 1.5).unwrap().is_none());
    }

    #[test]
    fn test_brewster_angle() {
        // Air to glass: θ_B = arctan(1.5) ≈ 56.3°
        let theta_b = brewster_angle(1.0, 1.5).unwrap();
        assert!((theta_b - 1.5_f64.atan()).abs() < 1e-10);
    }

    #[test]
    fn test_fresnel_normal_incidence() {
        // Normal incidence: r_s = r_p = (n₁-n₂)/(n₁+n₂)
        let rs = fresnel_rs(1.0, 0.0, 1.5, 0.0).unwrap();
        let rp = fresnel_rp(1.0, 0.0, 1.5, 0.0).unwrap();
        let expected = (1.0 - 1.5) / (1.0 + 1.5);
        assert!((rs - expected).abs() < 1e-10);
        assert!((rp - (-expected)).abs() < 1e-10); // rp has opposite sign convention
    }

    #[test]
    fn test_reflectance_normal() {
        // Air to glass: R = ((1-1.5)/(1+1.5))² = 0.04
        let r = reflectance_normal(1.0, 1.5).unwrap();
        assert!((r - 0.04).abs() < 1e-10);
    }

    #[test]
    fn test_transmittance_complement() {
        let r = reflectance_normal(1.0, 1.5).unwrap();
        let t = transmittance_normal(1.0, 1.5).unwrap();
        assert!((r + t - 1.0).abs() < 1e-10);
    }

    // ── Waveguide tests ────────────────────────────────────────────

    #[test]
    fn test_rectangular_waveguide_te10() {
        // Standard WR-90: a = 22.86mm, b = 10.16mm
        // TE₁₀ cutoff: f_c = c/(2a) ≈ 6.56 GHz
        let a = 22.86e-3;
        let b = 10.16e-3;
        let fc = rectangular_waveguide_cutoff(a, b, 1, 0, SPEED_OF_LIGHT).unwrap();
        assert!((fc / 1e9 - 6.56).abs() < 0.01);
    }

    #[test]
    fn test_rectangular_waveguide_cutoff_wavelength() {
        let a = 22.86e-3;
        let b = 10.16e-3;
        let lambda_c = rectangular_waveguide_cutoff_wavelength(a, b, 1, 0).unwrap();
        assert!((lambda_c - 2.0 * a).abs() < 1e-10); // TE₁₀: λ_c = 2a
    }

    #[test]
    fn test_guide_wavelength() {
        let lambda_c = 2.0 * 22.86e-3; // TE₁₀ cutoff wavelength
        let lambda = 30e-3; // 10 GHz in free space
        let lambda_g = guide_wavelength(lambda, lambda_c).unwrap();
        assert!(lambda_g > lambda); // guide wavelength always longer
    }

    #[test]
    fn test_guide_wavelength_below_cutoff() {
        assert!(guide_wavelength(50e-3, 40e-3).is_err()); // λ > λ_c
    }

    #[test]
    fn test_cylindrical_waveguide_cutoff() {
        // TE₁₁ mode: x'₁₁ = 1.8412
        let r = 10e-3; // 10mm radius
        let fc = cylindrical_waveguide_cutoff(r, 1.8412, SPEED_OF_LIGHT).unwrap();
        // f_c ≈ 8.79 GHz
        assert!((fc / 1e9 - 8.79).abs() < 0.1);
    }

    // ── Antenna tests ──────────────────────────────────────────────

    #[test]
    fn test_hertzian_dipole_pattern_peak() {
        // Peak at θ = π/2 (equatorial plane)
        assert!((hertzian_dipole_pattern(std::f64::consts::FRAC_PI_2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hertzian_dipole_pattern_null() {
        // Null along axis (θ = 0)
        assert!(hertzian_dipole_pattern(0.0).abs() < 1e-10);
    }

    #[test]
    fn test_half_wave_dipole_pattern_peak() {
        let peak = half_wave_dipole_pattern(std::f64::consts::FRAC_PI_2);
        assert!((peak - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_half_wave_dipole_pattern_null() {
        assert!(half_wave_dipole_pattern(0.0).abs() < 1e-10);
    }

    #[test]
    fn test_hertzian_radiation_resistance() {
        // dl/λ = 0.01 → R_rad = 80π² × 0.0001 ≈ 0.0789 Ω
        let r = hertzian_dipole_radiation_resistance(0.01, 1.0).unwrap();
        assert!((r - 80.0 * std::f64::consts::PI * std::f64::consts::PI * 1e-4).abs() < 1e-6);
    }

    #[test]
    fn test_effective_aperture() {
        // Half-wave dipole at 1m wavelength
        let ae = effective_aperture(1.0, half_wave_dipole_directivity()).unwrap();
        // A_e = λ²D/(4π) ≈ 0.1305 m²
        assert!((ae - 0.1305).abs() < 0.001);
    }

    #[test]
    fn test_friis_inverse_square() {
        // Doubling distance → 1/4 power
        let p1 = friis_transmission(1.0, 1.0, 1.0, 1.0).unwrap();
        let p2 = friis_transmission(1.0, 1.0, 1.0, 2.0).unwrap();
        assert!((p1 / p2 - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_friis_singularity() {
        assert!(friis_transmission(1.0, 1.0, 1.0, 0.0).is_err());
    }
}
