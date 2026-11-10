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
    let cos_i = theta_i.cos();
    let cos_t = theta_t.cos();
    let num = n1 * cos_i - n2 * cos_t;
    let den = n1 * cos_i + n2 * cos_t;
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
    let cos_i = theta_i.cos();
    let cos_t = theta_t.cos();
    let num = n2 * cos_i - n1 * cos_t;
    let den = n2 * cos_i + n1 * cos_t;
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

// ── Trig-free Fresnel (direct cosine interface) ──────────────────
//
// Conventions shared by the `*_direct`, unpolarized, Schlick and
// material-based functions below:
// - n₁ (incident side) and n₂ must be positive and finite (lossless media).
// - cos θ_i is the cosine of the angle between the incident ray and the
//   surface normal on the incident side, in [0, 1]. Values within
//   `COS_THETA_TOLERANCE` (1e-6) outside that range — rounding from
//   normalized f32/f64 vectors — are clamped; anything else (including NaN or
//   a cosine computed against the wrong-side normal) is an error. No
//   automatic side flip is done: swap n₁ and n₂ yourself for rays leaving
//   the denser medium.

/// Tolerance for cos θ_i slightly outside [0, 1] (rounding); clamped.
const COS_THETA_TOLERANCE: f64 = 1e-6;

/// Validate refractive indices for the trig-free Fresnel family.
#[inline]
fn check_indices(n1: f64, n2: f64) -> Result<()> {
    if n2 == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "refractive index n₂ cannot be zero".into(),
        });
    }
    if !(n1.is_finite() && n1 > 0.0 && n2.is_finite() && n2 > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("refractive indices must be positive and finite, got n₁={n1}, n₂={n2}"),
        });
    }
    Ok(())
}

/// Validate cos θ_i ∈ [0, 1] (clamping within `COS_THETA_TOLERANCE`).
#[inline]
fn check_cos_theta_i(cos_theta_i: f64) -> Result<f64> {
    if (0.0..=1.0).contains(&cos_theta_i) {
        Ok(cos_theta_i)
    } else if (-COS_THETA_TOLERANCE..=1.0 + COS_THETA_TOLERANCE).contains(&cos_theta_i) {
        Ok(cos_theta_i.clamp(0.0, 1.0))
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!(
                "cos θ_i must lie in [0, 1] (measured against the incident-side normal), got {cos_theta_i}"
            ),
        })
    }
}

/// Validated (cos θ_i, cos θ_t); cos θ_t is `None` under total internal reflection.
#[inline]
fn fresnel_cosines(n1: f64, n2: f64, cos_theta_i: f64) -> Result<(f64, Option<f64>)> {
    check_indices(n1, n2)?;
    let cos_i = check_cos_theta_i(cos_theta_i)?;
    let ratio = n1 / n2;
    let sin2_t = ratio * ratio * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        Ok((cos_i, None)) // total internal reflection
    } else {
        Ok((cos_i, Some((1.0 - sin2_t).sqrt())))
    }
}

/// a/b for a Fresnel ratio whose denominator is a sum of non-negative terms.
#[inline]
fn fresnel_ratio(num: f64, den: f64) -> Result<f64> {
    if den == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "Fresnel denominator is zero".into(),
        });
    }
    Ok(num / den)
}

/// Compute cos θ_t from cos θ_i via Snell's law (no trig).
///
/// cos θ_t = √(1 − (n₁/n₂)²(1 − cos²θ_i))
///
/// Returns `None` for total internal reflection (negative radicand).
///
/// # Errors
/// [`BijliError::DivisionByZero`] for n₂ = 0; [`BijliError::InvalidParameter`]
/// for non-positive/non-finite indices or cos θ_i outside [0, 1] (beyond a
/// 1e-6 rounding tolerance, which is clamped).
#[inline]
pub fn snell_cos_theta_t(n1: f64, n2: f64, cos_theta_i: f64) -> Result<Option<f64>> {
    Ok(fresnel_cosines(n1, n2, cos_theta_i)?.1)
}

/// Trig-free Fresnel reflection coefficient for s-polarization (TE).
///
/// Takes cos θ_i directly; computes cos θ_t via Snell's law internally.
/// Returns `None` for total internal reflection. For n₁ = n₂ it is exactly
/// 0 at every angle (including grazing incidence, where the formula is 0/0).
///
/// r_s = (n₁ cos θ_i − n₂ cos θ_t) / (n₁ cos θ_i + n₂ cos θ_t)
///
/// # Errors
/// As [`snell_cos_theta_t`].
#[inline]
pub fn fresnel_rs_direct(n1: f64, n2: f64, cos_theta_i: f64) -> Result<Option<f64>> {
    let (cos_i, cos_t) = fresnel_cosines(n1, n2, cos_theta_i)?;
    let Some(cos_t) = cos_t else {
        return Ok(None); // TIR → |r| = 1
    };
    if n1 == n2 {
        return Ok(Some(0.0));
    }
    let num = n1 * cos_i - n2 * cos_t;
    let den = n1 * cos_i + n2 * cos_t;
    fresnel_ratio(num, den).map(Some)
}

/// Trig-free Fresnel reflection coefficient for p-polarization (TM).
///
/// Takes cos θ_i directly; computes cos θ_t via Snell's law internally.
/// Returns `None` for total internal reflection. For n₁ = n₂ it is exactly
/// 0 at every angle.
///
/// r_p = (n₂ cos θ_i − n₁ cos θ_t) / (n₂ cos θ_i + n₁ cos θ_t)
///
/// # Errors
/// As [`snell_cos_theta_t`].
#[inline]
pub fn fresnel_rp_direct(n1: f64, n2: f64, cos_theta_i: f64) -> Result<Option<f64>> {
    let (cos_i, cos_t) = fresnel_cosines(n1, n2, cos_theta_i)?;
    let Some(cos_t) = cos_t else {
        return Ok(None); // TIR → |r| = 1
    };
    if n1 == n2 {
        return Ok(Some(0.0));
    }
    let num = n2 * cos_i - n1 * cos_t;
    let den = n2 * cos_i + n1 * cos_t;
    fresnel_ratio(num, den).map(Some)
}

/// Trig-free reflectance for s-polarization: R_s = |r_s|².
///
/// Returns 1.0 for total internal reflection.
#[inline]
pub fn reflectance_s_direct(n1: f64, n2: f64, cos_theta_i: f64) -> Result<f64> {
    match fresnel_rs_direct(n1, n2, cos_theta_i)? {
        Some(r) => Ok(r * r),
        None => Ok(1.0), // TIR
    }
}

/// Trig-free reflectance for p-polarization: R_p = |r_p|².
///
/// Returns 1.0 for total internal reflection.
#[inline]
pub fn reflectance_p_direct(n1: f64, n2: f64, cos_theta_i: f64) -> Result<f64> {
    match fresnel_rp_direct(n1, n2, cos_theta_i)? {
        Some(r) => Ok(r * r),
        None => Ok(1.0), // TIR
    }
}

/// Trig-free unpolarized reflectance: R = (R_s + R_p) / 2, in [0, 1].
///
/// Returns 1.0 for total internal reflection. cos θ_t is computed once and
/// shared by both polarizations.
///
/// # Errors
/// As [`snell_cos_theta_t`].
#[inline]
pub fn reflectance_unpolarized(n1: f64, n2: f64, cos_theta_i: f64) -> Result<f64> {
    let (cos_i, cos_t) = fresnel_cosines(n1, n2, cos_theta_i)?;
    let Some(cos_t) = cos_t else {
        return Ok(1.0); // TIR
    };
    if n1 == n2 {
        return Ok(0.0);
    }
    let rs = fresnel_ratio(n1 * cos_i - n2 * cos_t, n1 * cos_i + n2 * cos_t)?;
    let rp = fresnel_ratio(n2 * cos_i - n1 * cos_t, n2 * cos_i + n1 * cos_t)?;
    Ok(0.5 * (rs * rs + rp * rp))
}

/// Trig-free unpolarized transmittance: T = 1 − R, in [0, 1].
#[inline]
pub fn transmittance_unpolarized(n1: f64, n2: f64, cos_theta_i: f64) -> Result<f64> {
    Ok(1.0 - reflectance_unpolarized(n1, n2, cos_theta_i)?)
}

/// Schlick's approximation for reflectance (common in rendering).
///
/// R ≈ R₀ + (1 − R₀)(1 − cos θ)⁵, R₀ = ((n₁−n₂)/(n₁+n₂))²,
///
/// where θ is the angle on the **optically thinner** side: θ = θ_i when
/// n₁ ≤ n₂ and θ = θ_t when n₁ > n₂. Under total internal reflection the
/// result is exactly 1.0 (de Greve, *Reflections and Refractions in Ray
/// Tracing*, 2006). Using cos θ_i for n₁ > n₂ underestimates R badly near
/// the critical angle.
///
/// # Errors
/// As [`snell_cos_theta_t`].
#[inline]
pub fn schlick_reflectance(n1: f64, n2: f64, cos_theta_i: f64) -> Result<f64> {
    let (cos_i, cos_t) = fresnel_cosines(n1, n2, cos_theta_i)?;
    let cos = if n1 > n2 {
        match cos_t {
            Some(c) => c,
            None => return Ok(1.0), // TIR
        }
    } else {
        cos_i
    };
    let r0 = reflectance_normal(n1, n2)?;
    let one_minus_cos = 1.0 - cos;
    let omc2 = one_minus_cos * one_minus_cos;
    let omc5 = omc2 * omc2 * one_minus_cos;
    Ok(r0 + (1.0 - r0) * omc5)
}

// ── Material-based Fresnel interface ───────────────────────────────

/// Trig-free unpolarized reflectance at the interface between two
/// **lossless** materials, valid for magnetic media.
///
/// Snell's law uses n = √(ε_r μ_r); the amplitudes use the relative wave
/// impedance η = √(μ_r/ε_r):
///
/// ```text
/// r_s = (η₂ cos θ_i − η₁ cos θ_t) / (η₂ cos θ_i + η₁ cos θ_t)
/// r_p = (η₁ cos θ_i − η₂ cos θ_t) / (η₁ cos θ_i + η₂ cos θ_t)
/// R   = (r_s² + r_p²)/2          (1.0 under total internal reflection)
/// ```
///
/// For non-magnetic media this equals [`reflectance_unpolarized`] with
/// n = √ε_r. `cos_theta_i` follows the same rules as there.
///
/// # Errors
/// [`BijliError::InvalidPermittivity`] / [`BijliError::InvalidPermeability`]
/// for non-positive or non-finite ε_r / μ_r; [`BijliError::InvalidParameter`]
/// if either material has conductivity σ ≠ 0 (a conductor's Fresnel
/// coefficients need its frequency-dependent complex index) or for an
/// out-of-range cos θ_i.
#[inline]
pub fn reflectance_at_interface(
    mat1: &crate::material::Material,
    mat2: &crate::material::Material,
    cos_theta_i: f64,
) -> Result<f64> {
    let (n1, eta1) = mat1.lossless_optical_constants()?;
    let (n2, eta2) = mat2.lossless_optical_constants()?;
    let (cos_i, cos_t) = fresnel_cosines(n1, n2, cos_theta_i)?;
    let Some(cos_t) = cos_t else {
        return Ok(1.0); // TIR
    };
    if n1 == n2 {
        // cos θ_t = cos θ_i at every angle: R = Γ², Γ = (η₂ − η₁)/(η₂ + η₁)
        // (also covers grazing incidence, where the general form is 0/0).
        let gamma = (eta2 - eta1) / (eta2 + eta1);
        return Ok(gamma * gamma);
    }
    let rs = fresnel_ratio(eta2 * cos_i - eta1 * cos_t, eta2 * cos_i + eta1 * cos_t)?;
    let rp = fresnel_ratio(eta1 * cos_i - eta2 * cos_t, eta1 * cos_i + eta2 * cos_t)?;
    Ok(0.5 * (rs * rs + rp * rp))
}

/// Normal-incidence reflectance between two lossless materials:
/// R = ((η₂ − η₁)/(η₂ + η₁))², η = √(μ_r/ε_r).
///
/// Equivalent to [`crate::material::Material::try_reflectance_at`].
///
/// # Errors
/// As [`reflectance_at_interface`].
#[inline]
pub fn reflectance_normal_materials(
    mat1: &crate::material::Material,
    mat2: &crate::material::Material,
) -> Result<f64> {
    mat1.try_reflectance_at(mat2)
}

/// Transmission line reflection coefficient from material impedances.
///
/// Bridges `circuit`/`rf` reflection coefficient with `wave`/`material` modules.
/// Γ = (η₂ − η₁)/(η₂ + η₁) where η = √(μ/ε) is the (lossless) wave impedance.
///
/// # Errors
/// [`BijliError::InvalidPermittivity`] / [`BijliError::InvalidPermeability`]
/// for non-positive or non-finite ε_r / μ_r; [`BijliError::InvalidParameter`]
/// if either material has σ ≠ 0 (a conductor's impedance is complex and
/// frequency-dependent).
#[inline]
pub fn material_reflection_coefficient(
    mat1: &crate::material::Material,
    mat2: &crate::material::Material,
) -> Result<f64> {
    // Validates ε_r, μ_r and rejects conductors.
    mat1.lossless_optical_constants()?;
    mat2.lossless_optical_constants()?;
    let eta1 = mat1.impedance()?;
    let eta2 = mat2.impedance()?;
    let den = eta2 + eta1;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "η₁ + η₂ = 0".into(),
        });
    }
    Ok((eta2 - eta1) / den)
}

// ── Waveguides ─────────────────────────────────────────────────────

/// Validate rectangular waveguide dimensions and mode indices.
#[inline]
fn check_rectangular_waveguide(a: f64, b: f64, m: u32, n: u32) -> Result<()> {
    if !(a.is_finite() && a > 0.0 && b.is_finite() && b > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("waveguide dimensions must be positive and finite, got a={a}, b={b}"),
        });
    }
    if m == 0 && n == 0 {
        return Err(BijliError::InvalidParameter {
            reason: "both m and n cannot be zero (no TE₀₀/TM₀₀ mode)".into(),
        });
    }
    Ok(())
}

/// Cutoff frequency of a rectangular waveguide TE_mn or TM_mn mode.
///
/// f_c = (v/2)√((m/a)² + (n/b)²)
///
/// `a` and `b` are the waveguide dimensions (a ≥ b by convention); `velocity`
/// is the phase velocity of the filling medium (c for air).
///
/// # Errors
/// [`BijliError::InvalidParameter`] for non-positive/non-finite `a`, `b` or
/// `velocity`, or m = n = 0 (consistent with
/// [`rectangular_waveguide_cutoff_wavelength`]).
#[inline]
pub fn rectangular_waveguide_cutoff(a: f64, b: f64, m: u32, n: u32, velocity: f64) -> Result<f64> {
    check_rectangular_waveguide(a, b, m, n)?;
    if !(velocity.is_finite() && velocity > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("wave velocity must be positive and finite, got {velocity}"),
        });
    }
    let ma = m as f64 / a;
    let nb = n as f64 / b;
    Ok(0.5 * velocity * (ma * ma + nb * nb).sqrt())
}

/// Cutoff wavelength of a rectangular waveguide mode.
///
/// λ_c = 2/√((m/a)² + (n/b)²)
///
/// # Errors
/// [`BijliError::InvalidParameter`] for non-positive/non-finite `a`, `b`, or
/// m = n = 0.
#[inline]
pub fn rectangular_waveguide_cutoff_wavelength(a: f64, b: f64, m: u32, n: u32) -> Result<f64> {
    check_rectangular_waveguide(a, b, m, n)?;
    let ma = m as f64 / a;
    let nb = n as f64 / b;
    let denom = (ma * ma + nb * nb).sqrt();
    if denom.is_nan() || denom <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("cutoff wavelength not representable for a={a}, b={b}"),
        });
    }
    Ok(2.0 / denom)
}

/// Guide wavelength: λ_g = λ/√(1 − (λ/λ_c)²).
///
/// Only valid when λ < λ_c (above cutoff). λ_c = +∞ (a TEM line) gives λ_g = λ.
///
/// # Errors
/// [`BijliError::DivisionByZero`] for λ_c = 0; [`BijliError::InvalidParameter`]
/// for non-positive/NaN wavelengths, a non-finite λ, or λ ≥ λ_c.
#[inline]
pub fn guide_wavelength(free_space_wavelength: f64, cutoff_wavelength: f64) -> Result<f64> {
    if cutoff_wavelength == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "cutoff wavelength cannot be zero".into(),
        });
    }
    if cutoff_wavelength.is_nan() || cutoff_wavelength < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("cutoff wavelength must be positive, got {cutoff_wavelength}"),
        });
    }
    if !(free_space_wavelength.is_finite() && free_space_wavelength > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive and finite, got {free_space_wavelength}"),
        });
    }
    let ratio = free_space_wavelength / cutoff_wavelength;
    let radicand = 1.0 - ratio * ratio;
    if radicand.is_nan() || radicand <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("frequency below cutoff: λ/λ_c = {ratio:.4} (must be < 1)"),
        });
    }
    Ok(free_space_wavelength / radicand.sqrt())
}

/// Cutoff frequency of a circular (cylindrical) waveguide mode.
///
/// f_c = x · v / (2πa)
///
/// where `bessel_root` x is the mode's characteristic root and `radius` is a:
/// - **TE_mn**: x = x′_mn, the n-th root of J′_m — e.g. TE₁₁ 1.8412 (the
///   dominant mode), TE₂₁ 3.0542, TE₀₁ 3.8317;
/// - **TM_mn**: x = x_mn, the n-th root of J_m — e.g. TM₀₁ 2.4048,
///   TM₁₁ 3.8317.
///
/// `velocity` is the phase velocity of the filling medium (c for air).
///
/// # Errors
/// [`BijliError::InvalidParameter`] unless `radius`, `bessel_root` and
/// `velocity` are positive and finite.
#[inline]
pub fn cylindrical_waveguide_cutoff(radius: f64, bessel_root: f64, velocity: f64) -> Result<f64> {
    if !(radius.is_finite() && radius > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("waveguide radius must be positive and finite, got {radius}"),
        });
    }
    if !(bessel_root.is_finite() && bessel_root > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Bessel root must be positive and finite, got {bessel_root}"),
        });
    }
    if !(velocity.is_finite() && velocity > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("wave velocity must be positive and finite, got {velocity}"),
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

    // ── Trig-free Fresnel tests ───────────────────────────────────

    #[test]
    fn test_snell_cos_theta_t_normal() {
        // Normal incidence: cos θ_i = 1 → cos θ_t = 1
        let cos_t = snell_cos_theta_t(1.0, 1.5, 1.0).unwrap().unwrap();
        assert!((cos_t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_snell_cos_theta_t_45_deg() {
        // Air to glass at 45°: cos(45°) = √2/2
        let cos_i = std::f64::consts::FRAC_PI_4.cos();
        let cos_t = snell_cos_theta_t(1.0, 1.5, cos_i).unwrap().unwrap();
        // Verify via trig: sin θ_t = sin(45°)/1.5, cos θ_t = √(1 - sin²θ_t)
        let sin_t = std::f64::consts::FRAC_PI_4.sin() / 1.5;
        let expected = (1.0 - sin_t * sin_t).sqrt();
        assert!((cos_t - expected).abs() < 1e-10);
    }

    #[test]
    fn test_snell_cos_theta_t_tir() {
        // Glass to air at steep angle → TIR
        let cos_i = 0.3; // steep angle (cos θ_i small → θ_i large)
        assert!(snell_cos_theta_t(1.5, 1.0, cos_i).unwrap().is_none());
    }

    #[test]
    fn test_fresnel_rs_direct_normal() {
        // Normal incidence: r_s = (n₁ - n₂)/(n₁ + n₂)
        let r = fresnel_rs_direct(1.0, 1.5, 1.0).unwrap().unwrap();
        let expected = (1.0 - 1.5) / (1.0 + 1.5);
        assert!((r - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fresnel_rp_direct_normal() {
        // Normal incidence: r_p = (n₂ - n₁)/(n₂ + n₁)
        let r = fresnel_rp_direct(1.0, 1.5, 1.0).unwrap().unwrap();
        let expected = (1.5 - 1.0) / (1.5 + 1.0);
        assert!((r - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fresnel_direct_matches_trig() {
        // Verify trig-free matches trig-based at 30°
        let theta_i = std::f64::consts::PI / 6.0;
        let cos_i = theta_i.cos();
        let theta_t = snell_refraction_angle(1.0, 1.5, theta_i).unwrap().unwrap();

        let rs_trig = fresnel_rs(1.0, theta_i, 1.5, theta_t).unwrap();
        let rs_direct = fresnel_rs_direct(1.0, 1.5, cos_i).unwrap().unwrap();
        assert!((rs_trig - rs_direct).abs() < 1e-10);

        let rp_trig = fresnel_rp(1.0, theta_i, 1.5, theta_t).unwrap();
        let rp_direct = fresnel_rp_direct(1.0, 1.5, cos_i).unwrap().unwrap();
        assert!((rp_trig - rp_direct).abs() < 1e-10);
    }

    #[test]
    fn test_fresnel_direct_tir() {
        // Glass to air, steep angle → TIR, returns None
        assert!(fresnel_rs_direct(1.5, 1.0, 0.3).unwrap().is_none());
        assert!(fresnel_rp_direct(1.5, 1.0, 0.3).unwrap().is_none());
    }

    #[test]
    fn test_reflectance_s_direct_tir() {
        // TIR → reflectance = 1.0
        let r = reflectance_s_direct(1.5, 1.0, 0.3).unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflectance_p_direct_tir() {
        let r = reflectance_p_direct(1.5, 1.0, 0.3).unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflectance_unpolarized_normal() {
        // At normal incidence, R_s = R_p, so unpolarized = same as normal
        let r = reflectance_unpolarized(1.0, 1.5, 1.0).unwrap();
        let r_normal = reflectance_normal(1.0, 1.5).unwrap();
        assert!((r - r_normal).abs() < 1e-10);
    }

    #[test]
    fn test_transmittance_unpolarized_complement() {
        let cos_i = std::f64::consts::FRAC_PI_4.cos();
        let r = reflectance_unpolarized(1.0, 1.5, cos_i).unwrap();
        let t = transmittance_unpolarized(1.0, 1.5, cos_i).unwrap();
        assert!((r + t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_schlick_normal_matches_exact() {
        // At normal incidence, Schlick = exact
        let r_schlick = schlick_reflectance(1.0, 1.5, 1.0).unwrap();
        let r_exact = reflectance_normal(1.0, 1.5).unwrap();
        assert!((r_schlick - r_exact).abs() < 1e-10);
    }

    #[test]
    fn test_schlick_grazing() {
        // At grazing incidence (cos θ → 0), Schlick → 1.0
        let r = schlick_reflectance(1.0, 1.5, 0.0).unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_schlick_approximates_fresnel() {
        // Schlick should be close to exact unpolarized Fresnel
        let cos_i = std::f64::consts::FRAC_PI_4.cos();
        let r_schlick = schlick_reflectance(1.0, 1.5, cos_i).unwrap();
        let r_exact = reflectance_unpolarized(1.0, 1.5, cos_i).unwrap();
        // Schlick is an approximation — within a few percent
        assert!((r_schlick - r_exact).abs() < 0.02);
    }

    #[test]
    fn test_snell_cos_theta_t_zero_n2() {
        assert!(snell_cos_theta_t(1.0, 0.0, 0.5).is_err());
    }

    // ── Material-based interface tests ────────────────────────────

    #[test]
    fn test_reflectance_at_interface_normal() {
        use crate::material::Material;
        let air = Material::vacuum();
        let glass = Material::dielectric(2.25);
        let r = reflectance_at_interface(&air, &glass, 1.0).unwrap();
        let expected = reflectance_normal(1.0, 1.5).unwrap();
        assert!((r - expected).abs() < 1e-6);
    }

    #[test]
    fn test_reflectance_normal_materials_air_glass() {
        use crate::material::Material;
        let air = Material::vacuum();
        let glass = Material::dielectric(2.25);
        let r = reflectance_normal_materials(&air, &glass).unwrap();
        assert!((r - 0.04).abs() < 1e-6);
    }

    #[test]
    fn test_material_reflection_coefficient() {
        use crate::material::Material;
        let air = Material::vacuum();
        let glass = Material::dielectric(2.25);
        // Γ = (η₂ - η₁)/(η₂ + η₁), and η ∝ 1/n for non-magnetic
        // η_glass = η₀/1.5, Γ = (η₀/1.5 - η₀)/(η₀/1.5 + η₀) = (1/1.5-1)/(1/1.5+1) = (-1/3)/(5/3) = -0.2
        let g = material_reflection_coefficient(&air, &glass).unwrap();
        assert!((g - (-0.2)).abs() < 1e-6);
    }

    #[test]
    fn test_material_reflection_coefficient_same_material() {
        use crate::material::Material;
        let m = Material::dielectric(4.0);
        let g = material_reflection_coefficient(&m, &m).unwrap();
        assert!(g.abs() < 1e-10);
    }

    // ── Fresnel regression tests ──────────────────────────────────
    //
    // Reference values: exact Fresnel equations evaluated independently in
    // Python from Jackson (3rd ed.) eqs. (7.39) and (7.41), which are written
    // in terms of n and μ (not impedances).

    #[test]
    fn test_schlick_dense_to_rare_uses_cos_theta_t() {
        // Glass → air near the critical angle (41.8°). Exact R(cos θ_i = 0.76)
        // = 0.30233831088147745; old Schlick (cos θ_i) gave 0.0408.
        let r = schlick_reflectance(1.5, 1.0, 0.76).unwrap();
        assert!((r - 0.302_338_310_881_477_45).abs() < 0.02, "{r}");
        let r = schlick_reflectance(1.5, 1.0, 0.9).unwrap();
        assert!((r - 0.046_332_647_954_037_66).abs() < 0.01, "{r}");
        // Under TIR R = 1 exactly (old: 0.07).
        assert_eq!(schlick_reflectance(1.5, 1.0, 0.5).unwrap(), 1.0);
        // Rare → dense is unchanged.
        let r0 = reflectance_normal(1.0, 1.5).unwrap();
        let expected = r0 + (1.0 - r0) * 0.5_f64.powi(5);
        assert!((schlick_reflectance(1.0, 1.5, 0.5).unwrap() - expected).abs() < 1e-15);
    }

    #[test]
    fn test_cos_theta_i_range_checked() {
        // Old code: reflectance_unpolarized(1, 1.5, −1) = 25.
        assert!(reflectance_unpolarized(1.0, 1.5, -1.0).is_err());
        assert!(reflectance_unpolarized(1.0, 1.5, 1.5).is_err());
        assert!(reflectance_unpolarized(1.0, 1.5, f64::NAN).is_err());
        assert!(transmittance_unpolarized(1.0, 1.5, -0.5).is_err());
        assert!(schlick_reflectance(1.0, 1.5, -0.2).is_err());
        assert!(fresnel_rs_direct(1.0, 1.5, 2.0).is_err());
        assert!(snell_cos_theta_t(1.0, 1.5, -0.1).is_err());
        // Rounding-level excursions are clamped.
        let normal = reflectance_unpolarized(1.0, 1.5, 1.0).unwrap();
        assert_eq!(
            reflectance_unpolarized(1.0, 1.5, 1.0 + 1e-9).unwrap(),
            normal
        );
        let grazing = reflectance_unpolarized(1.0, 1.5, 0.0).unwrap();
        assert_eq!(reflectance_unpolarized(1.0, 1.5, -1e-9).unwrap(), grazing);
        assert!((grazing - 1.0).abs() < 1e-15); // R → 1 at grazing incidence
    }

    #[test]
    fn test_invalid_indices_rejected() {
        assert!(reflectance_unpolarized(-1.0, 1.5, 0.5).is_err());
        assert!(reflectance_unpolarized(1.0, f64::NAN, 0.5).is_err());
        assert!(snell_cos_theta_t(f64::INFINITY, 1.5, 0.5).is_err());
    }

    #[test]
    fn test_reflectance_and_transmittance_in_unit_interval() {
        for (n1, n2) in [(1.0, 1.5), (1.5, 1.0), (1.0, 3.5), (2.4, 1.0), (1.0, 1.0)] {
            for i in 0..=1000 {
                let c = f64::from(i) / 1000.0;
                let r = reflectance_unpolarized(n1, n2, c).unwrap();
                let t = transmittance_unpolarized(n1, n2, c).unwrap();
                assert!((0.0..=1.0).contains(&r), "R={r} n1={n1} n2={n2} c={c}");
                assert!((0.0..=1.0).contains(&t), "T={t} n1={n1} n2={n2} c={c}");
                let s = schlick_reflectance(n1, n2, c).unwrap();
                assert!((0.0..=1.0).contains(&s), "Schlick={s}");
            }
        }
    }

    #[test]
    fn test_index_matched_grazing_is_zero() {
        // Old code: Err(DivisionByZero) from the 0/0 at cos θ_i = 0.
        assert_eq!(fresnel_rs_direct(1.0, 1.0, 0.0).unwrap(), Some(0.0));
        assert_eq!(fresnel_rp_direct(1.0, 1.0, 0.0).unwrap(), Some(0.0));
        assert_eq!(reflectance_unpolarized(1.33, 1.33, 0.0).unwrap(), 0.0);
        assert_eq!(transmittance_unpolarized(1.33, 1.33, 0.0).unwrap(), 1.0);
    }

    #[test]
    fn test_reflectance_unpolarized_single_pass_matches_per_polarization() {
        // The shared-cos θ_t fast path must reproduce (R_s + R_p)/2 exactly.
        for (n1, n2) in [(1.0, 1.5), (1.5, 1.0), (1.0, 2.42), (1.33, 1.0)] {
            for i in 0..=200 {
                let c = f64::from(i) / 200.0;
                let r = reflectance_unpolarized(n1, n2, c).unwrap();
                let rs = reflectance_s_direct(n1, n2, c).unwrap();
                let rp = reflectance_p_direct(n1, n2, c).unwrap();
                assert_eq!(r, 0.5 * (rs + rp), "n1={n1} n2={n2} c={c}");
            }
        }
    }

    #[test]
    fn test_reflectance_at_interface_magnetic_media() {
        use crate::material::Material;
        let vac = Material::vacuum();
        let matched = Material::new(4.0, 4.0, 0.0); // n = 4, η = η₀
        // Old (n-only) values: 0.36, 0.3595, 0.3725.
        assert!(reflectance_at_interface(&vac, &matched, 1.0).unwrap().abs() < 1e-15);
        let r = reflectance_at_interface(&vac, &matched, 0.8).unwrap();
        assert!((r - 0.011_127_846_503_534_383).abs() < 1e-12, "{r}");
        let r = reflectance_at_interface(&vac, &matched, 0.3).unwrap();
        assert!((r - 0.278_768_233_115_061_8).abs() < 1e-12, "{r}");
        let r = reflectance_at_interface(&matched, &vac, 0.99).unwrap();
        assert!((r - 0.008_200_057_075_464_767).abs() < 1e-12, "{r}");
        // ε_r = 1, μ_r = 2.25 at 45°
        let mag = Material::new(1.0, 2.25, 0.0);
        let r = reflectance_at_interface(&vac, &mag, std::f64::consts::FRAC_PI_4.cos()).unwrap();
        assert!((r - 0.050_239_911_012_235_95).abs() < 1e-12, "{r}");
        // Same index, different impedance, grazing: R = Γ² (no 0/0 error).
        let a = Material::new(4.0, 4.0, 0.0);
        let b = Material::new(16.0, 1.0, 0.0); // n = 4, η = η₀/4
        let r = reflectance_at_interface(&a, &b, 0.0).unwrap();
        assert!((r - 0.36).abs() < 1e-15, "{r}");
        // Normal incidence wrapper agrees.
        assert!(reflectance_normal_materials(&vac, &matched).unwrap().abs() < 1e-15);
    }

    #[test]
    fn test_reflectance_at_interface_dielectric_matches_index_form() {
        use crate::material::Material;
        let air = Material::vacuum();
        let glass = Material::dielectric(2.25);
        for i in 0..=50 {
            let c = f64::from(i) / 50.0;
            let a = reflectance_at_interface(&air, &glass, c).unwrap();
            let b = reflectance_unpolarized(1.0, 1.5, c).unwrap();
            assert!((a - b).abs() < 1e-14, "c={c}: {a} vs {b}");
            let a = reflectance_at_interface(&glass, &air, c).unwrap();
            let b = reflectance_unpolarized(1.5, 1.0, c).unwrap();
            assert!((a - b).abs() < 1e-14, "c={c}: {a} vs {b}");
        }
    }

    #[test]
    fn test_material_fresnel_rejects_conductors_and_invalid() {
        use crate::error::BijliError;
        use crate::material::Material;
        let vac = Material::vacuum();
        let cu = Material::conductor(5.96e7); // old code: R = 0, Γ = 0
        assert!(reflectance_at_interface(&vac, &cu, 1.0).is_err());
        assert!(reflectance_normal_materials(&vac, &cu).is_err());
        assert!(material_reflection_coefficient(&vac, &cu).is_err());
        let neg = Material::dielectric(-2.0);
        assert!(matches!(
            reflectance_at_interface(&vac, &neg, 1.0),
            Err(BijliError::InvalidPermittivity { .. })
        ));
        assert!(matches!(
            reflectance_normal_materials(&vac, &neg),
            Err(BijliError::InvalidPermittivity { .. })
        ));
        assert!(matches!(
            material_reflection_coefficient(&vac, &neg),
            Err(BijliError::InvalidPermittivity { .. })
        ));
        let bad_mu = Material::new(1.0, -1.0, 0.0);
        assert!(matches!(
            reflectance_at_interface(&vac, &bad_mu, 1.0),
            Err(BijliError::InvalidPermeability { .. })
        ));
        assert!(reflectance_at_interface(&vac, &Material::dielectric(2.25), -0.5).is_err());
    }

    // ── Waveguide regression tests ────────────────────────────────

    #[test]
    fn test_rectangular_waveguide_cutoff_validation() {
        let (a, b) = (22.86e-3, 10.16e-3);
        // Old code: Ok(0.0) for m = n = 0 (the wavelength sibling errs).
        assert!(rectangular_waveguide_cutoff(a, b, 0, 0, SPEED_OF_LIGHT).is_err());
        assert!(rectangular_waveguide_cutoff(a, b, 1, 0, 0.0).is_err());
        assert!(rectangular_waveguide_cutoff(a, b, 1, 0, -SPEED_OF_LIGHT).is_err());
        assert!(rectangular_waveguide_cutoff(a, b, 1, 0, f64::NAN).is_err());
        assert!(rectangular_waveguide_cutoff(f64::NAN, b, 1, 0, SPEED_OF_LIGHT).is_err());
        assert!(rectangular_waveguide_cutoff_wavelength(a, f64::NAN, 1, 0).is_err());
        // TE₁₀ closed form f_c = c/(2a).
        let fc = rectangular_waveguide_cutoff(a, b, 1, 0, SPEED_OF_LIGHT).unwrap();
        assert!((fc - SPEED_OF_LIGHT / (2.0 * a)).abs() / fc < 1e-15);
    }

    #[test]
    fn test_guide_wavelength_validation() {
        assert!(guide_wavelength(f64::NAN, 0.05).is_err());
        assert!(guide_wavelength(0.03, f64::NAN).is_err());
        assert!(guide_wavelength(-0.03, 0.05).is_err());
        assert!(guide_wavelength(0.03, 0.0).is_err());
        // TEM line (no cutoff): λ_g = λ.
        assert_eq!(guide_wavelength(0.03, f64::INFINITY).unwrap(), 0.03);
    }

    #[test]
    fn test_cylindrical_waveguide_tm01_and_validation() {
        // TM₀₁ uses the first root of J₀, 2.4048: f_c = 2.4048 c/(2πa).
        let r = 10e-3;
        let fc = cylindrical_waveguide_cutoff(r, 2.4048, SPEED_OF_LIGHT).unwrap();
        let expected = 2.4048 * SPEED_OF_LIGHT / (2.0 * std::f64::consts::PI * r);
        assert!((fc - expected).abs() / expected < 1e-15);
        assert!((fc / 1e9 - 11.474).abs() < 0.01);
        assert!(cylindrical_waveguide_cutoff(f64::NAN, 1.8412, SPEED_OF_LIGHT).is_err());
        assert!(cylindrical_waveguide_cutoff(r, -1.8412, SPEED_OF_LIGHT).is_err());
        assert!(cylindrical_waveguide_cutoff(r, 1.8412, f64::NAN).is_err());
    }
}
