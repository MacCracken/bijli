//! Mie and Rayleigh scattering — cross-sections, efficiencies, phase functions.
//!
//! Mie theory gives the exact solution for scattering of a plane wave by a
//! homogeneous sphere. Rayleigh scattering is the small-particle limit (x ≪ 1).
//!
//! Size parameter: x = 2πr/λ where r is the sphere radius and λ is the wavelength.

use crate::error::{BijliError, Result};
use crate::polarization::{Complex, ComplexExt, complex_real, complex_zero};

// ── Mie scattering ────────────────────────────────────────────────

/// Mie scattering result for a homogeneous sphere.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MieResult {
    /// Scattering efficiency Q_sca.
    pub q_sca: f64,
    /// Extinction efficiency Q_ext.
    pub q_ext: f64,
    /// Absorption efficiency Q_abs = Q_ext − Q_sca.
    pub q_abs: f64,
    /// Asymmetry parameter g = ⟨cos θ⟩.
    pub g: f64,
    /// Mie coefficients a_n.
    pub a: Vec<Complex>,
    /// Mie coefficients b_n.
    pub b: Vec<Complex>,
}

/// Compute Mie scattering for a sphere.
///
/// # Parameters
/// - `x`: size parameter 2πr/λ (must be > 0)
/// - `m`: complex refractive index of the sphere relative to the medium
///
/// # Returns
/// [`MieResult`] with efficiencies, asymmetry parameter, and Mie coefficients.
pub fn mie(x: f64, m: Complex) -> Result<MieResult> {
    if x <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("size parameter x must be positive, got {x}"),
        });
    }

    tracing::debug!(x, m_re = m.re, m_im = m.im, "computing Mie scattering");

    let mx = Complex::new(m.re * x, m.im * x);

    // Number of terms: Wiscombe's criterion
    let n_max = (x + 4.0 * x.cbrt() + 2.0).ceil() as usize;
    let n_max = n_max.max(4);

    // Compute ψ_n(x) and ξ_n(x) via upward recurrence for Riccati-Bessel functions
    // ψ_n(z) = z j_n(z), ξ_n(z) = z h_n^(1)(z)
    // Use logarithmic derivative D_n(z) = d/dz [ln(z j_n(z))] = ψ_n'/ψ_n

    // Logarithmic derivatives D_n(mx) via downward recurrence (stable)
    let n_stop = n_max + 15;
    let mut d_mx = vec![complex_zero(); n_stop + 1];
    // Start from large n and recur downward: D_{n-1}(z) = n/z - 1/(D_n(z) + n/z)
    for n in (1..=n_stop).rev() {
        let nf = complex_real(n as f64);
        let ratio = nf * inv_complex(mx)?;
        d_mx[n - 1] = ratio - inv_complex(d_mx[n] + ratio)?;
    }

    // Riccati-Bessel functions ψ_n(x), ξ_n(x) via upward recurrence
    let mut psi_prev = x.sin(); // ψ_0(x) = sin(x)
    let mut psi_curr = x.sin() / x - x.cos(); // ψ_1(x)
    let mut xi_prev = Complex::new(x.sin(), x.cos()); // ξ_0(x) = sin(x) + i cos(x)
    let mut xi_curr = Complex::new(x.sin() / x - x.cos(), x.cos() / x + x.sin());

    let mut a_coeffs = Vec::with_capacity(n_max);
    let mut b_coeffs = Vec::with_capacity(n_max);
    let mut q_sca = 0.0;
    let mut q_ext = 0.0;
    let mut g_num = 0.0; // numerator for asymmetry parameter

    for n in 1..=n_max {
        if n > 1 {
            let nf = (2 * n - 1) as f64;
            let psi_next = nf / x * psi_curr - psi_prev;
            psi_prev = psi_curr;
            psi_curr = psi_next;

            let xi_next = complex_real(nf / x) * xi_curr - xi_prev;
            xi_prev = xi_curr;
            xi_curr = xi_next;
        }

        // Bohren & Huffman (1983) Mie coefficients:
        // a_n = (m D_n(mx) + n/x) ψ_n - ψ_{n-1} / ((m D_n(mx) + n/x) ξ_n - ξ_{n-1})
        // b_n = (D_n(mx)/m + n/x) ψ_n - ψ_{n-1} / ((D_n(mx)/m + n/x) ξ_n - ξ_{n-1})

        let d_n = d_mx[n]; // D_n(mx)
        let n_over_x = complex_real(n as f64 / x);

        let a_factor = complex_div(d_n, m)? + n_over_x;
        let a_num = a_factor * complex_real(psi_curr) - complex_real(psi_prev);
        let a_den = a_factor * xi_curr - xi_prev;
        let a_n = complex_div(a_num, a_den)?;

        let b_factor = m * d_n + n_over_x;
        let b_num = b_factor * complex_real(psi_curr) - complex_real(psi_prev);
        let b_den = b_factor * xi_curr - xi_prev;
        let b_n = complex_div(b_num, b_den)?;

        a_coeffs.push(a_n);
        b_coeffs.push(b_n);

        let nf = n as f64;
        let weight = 2.0 * nf + 1.0;

        q_sca += weight * (a_n.norm_sq() + b_n.norm_sq());
        q_ext += weight * (a_n.re + b_n.re);

        // Asymmetry parameter contributions
        if n > 1 {
            let a_prev = a_coeffs[n - 2];
            let b_prev = b_coeffs[n - 2];
            let n1 = (n - 1) as f64;
            g_num += n1 * (n1 + 2.0) / (n1 + 1.0) * (a_prev.conj() * a_n + b_prev.conj() * b_n).re;
            g_num += (2.0 * n1 + 1.0) / (n1 * (n1 + 1.0)) * (a_prev.conj() * b_prev).re;
        }
    }

    let x2 = x * x;
    q_sca *= 2.0 / x2;
    q_ext *= 2.0 / x2;
    let q_abs = q_ext - q_sca;

    let g = if q_sca > 1e-30 {
        4.0 / (x2 * q_sca) * g_num
    } else {
        0.0
    };

    Ok(MieResult {
        q_sca,
        q_ext,
        q_abs,
        g,
        a: a_coeffs,
        b: b_coeffs,
    })
}

/// Scattering cross-section from efficiency: σ_sca = Q_sca × πr².
#[inline]
#[must_use]
pub fn cross_section_from_efficiency(q: f64, radius: f64) -> f64 {
    q * std::f64::consts::PI * radius * radius
}

// ── Rayleigh scattering ───────────────────────────────────────────

/// Rayleigh scattering cross-section for a small sphere.
///
/// σ = (8π/3)(2πr/λ)⁴ r² |K|² where K = (m²−1)/(m²+2)
///
/// Valid when x = 2πr/λ ≪ 1 (particle much smaller than wavelength).
#[inline]
pub fn rayleigh_cross_section(radius: f64, wavelength: f64, m: Complex) -> Result<f64> {
    if radius <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("radius must be positive, got {radius}"),
        });
    }
    if wavelength <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive, got {wavelength}"),
        });
    }

    let x = 2.0 * std::f64::consts::PI * radius / wavelength;
    let m2 = m * m;
    let k = complex_div(m2 - complex_real(1.0), m2 + complex_real(2.0))?;
    let x4 = x * x * x * x;

    Ok(8.0 / 3.0 * std::f64::consts::PI * radius * radius * x4 * k.norm_sq())
}

/// Rayleigh scattering efficiency Q_sca = σ/(πr²).
#[inline]
pub fn rayleigh_efficiency(x: f64, m: Complex) -> Result<f64> {
    if x <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("size parameter must be positive, got {x}"),
        });
    }
    let m2 = m * m;
    let k = complex_div(m2 - complex_real(1.0), m2 + complex_real(2.0))?;
    Ok(8.0 / 3.0 * x * x * x * x * k.norm_sq())
}

/// Rayleigh phase function: P(θ) = (3/4)(1 + cos²θ).
///
/// Normalized so that integral over 4π solid angle = 4π.
#[inline]
#[must_use]
pub fn rayleigh_phase_function(cos_theta: f64) -> f64 {
    0.75 * (1.0 + cos_theta * cos_theta)
}

/// Rayleigh scattering intensity ratio I/I₀ at distance r and angle θ.
///
/// I/I₀ = (8π⁴α² / (r²λ⁴)) × (1 + cos²θ)/2
///
/// For a dielectric sphere, α = 3ε₀V(ε_r − 1)/(ε_r + 2).
/// This simplified version uses the cross-section directly.
#[inline]
pub fn rayleigh_angular_intensity(
    cross_section: f64,
    distance: f64,
    cos_theta: f64,
) -> Result<f64> {
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(
        cross_section / (distance * distance) * rayleigh_phase_function(cos_theta)
            / (4.0 * std::f64::consts::PI),
    )
}

// ── Complex arithmetic helpers ────────────────────────────────────

/// 1/z for a complex number.
#[inline]
fn inv_complex(z: Complex) -> Result<Complex> {
    let d = z.norm_sq();
    if d < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "complex division by zero in Mie computation".into(),
        });
    }
    Ok(Complex::new(z.re / d, -z.im / d))
}

/// a/b for complex numbers.
#[inline]
fn complex_div(a: Complex, b: Complex) -> Result<Complex> {
    let d = b.norm_sq();
    if d < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "complex division by zero in Mie computation".into(),
        });
    }
    Ok(Complex::new(
        (a.re * b.re + a.im * b.im) / d,
        (a.im * b.re - a.re * b.im) / d,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Mie tests ─────────────────────────────────────────────────

    #[test]
    fn test_mie_small_sphere_approaches_rayleigh() {
        // For very small x, Mie should approximate Rayleigh
        let x = 0.01;
        let m = Complex::new(1.5, 0.0);
        let mie_result = mie(x, m).unwrap();
        let ray_q = rayleigh_efficiency(x, m).unwrap();
        assert!(
            (mie_result.q_sca - ray_q).abs() / ray_q < 0.01,
            "Mie Q_sca={} vs Rayleigh Q={}",
            mie_result.q_sca,
            ray_q
        );
    }

    #[test]
    fn test_mie_nonabsorbing_q_ext_equals_q_sca() {
        // For real refractive index (no absorption), Q_abs = 0
        let x = 1.0;
        let m = complex_real(1.5);
        let result = mie(x, m).unwrap();
        assert!(result.q_abs.abs() < 1e-10);
    }

    #[test]
    fn test_mie_absorbing_sphere() {
        // For an absorbing sphere (im(m) < 0 in e^{+ikr} convention),
        // Q_ext should exceed Q_sca and Q_abs > 0.
        let x = 1.0;
        let m = Complex::new(1.5, -0.1); // negative im = absorption
        let result = mie(x, m).unwrap();
        assert!(result.q_ext > 0.0, "Q_ext = {}", result.q_ext);
        assert!(result.q_sca > 0.0, "Q_sca = {}", result.q_sca);
        assert!(result.q_abs > 0.0, "Q_abs = {}", result.q_abs);
    }

    #[test]
    fn test_mie_extinction_positive() {
        let x = 2.0;
        let m = Complex::new(1.33, 0.0);
        let result = mie(x, m).unwrap();
        assert!(result.q_ext > 0.0);
        assert!(result.q_sca > 0.0);
    }

    #[test]
    fn test_mie_asymmetry_small_sphere() {
        // Small sphere: g ≈ 0 (isotropic)
        let x = 0.01;
        let m = complex_real(1.5);
        let result = mie(x, m).unwrap();
        assert!(result.g.abs() < 0.01);
    }

    #[test]
    fn test_mie_asymmetry_large_sphere() {
        // Large sphere: g > 0 (forward scattering dominant)
        let x = 10.0;
        let m = complex_real(1.5);
        let result = mie(x, m).unwrap();
        assert!(result.g > 0.0);
    }

    #[test]
    fn test_mie_invalid_x() {
        assert!(mie(0.0, complex_real(1.5)).is_err());
        assert!(mie(-1.0, complex_real(1.5)).is_err());
    }

    #[test]
    fn test_mie_water_droplet() {
        // Water droplet (n = 1.33) at visible wavelength
        // x = 2πr/λ, r = 1μm, λ = 0.5μm → x ≈ 12.57
        let x = 2.0 * std::f64::consts::PI * 1.0 / 0.5;
        let m = Complex::new(1.33, 0.0);
        let result = mie(x, m).unwrap();
        // Q_ext should be around 2 for large dielectric spheres (extinction paradox)
        assert!(result.q_ext > 1.0);
        assert!(result.q_ext < 5.0);
    }

    #[test]
    fn test_mie_coefficients_count() {
        let x = 5.0;
        let m = complex_real(1.5);
        let result = mie(x, m).unwrap();
        assert!(result.a.len() >= 4);
        assert_eq!(result.a.len(), result.b.len());
    }

    // ── Rayleigh tests ────────────────────────────────────────────

    #[test]
    fn test_rayleigh_cross_section_lambda4() {
        // σ ∝ λ⁻⁴: halving wavelength → 16× cross-section
        let r = 1e-8;
        let m = complex_real(1.5);
        let s1 = rayleigh_cross_section(r, 500e-9, m).unwrap();
        let s2 = rayleigh_cross_section(r, 250e-9, m).unwrap();
        assert!((s2 / s1 - 16.0).abs() < 0.1);
    }

    #[test]
    fn test_rayleigh_cross_section_r6() {
        // σ ∝ r⁶: doubling radius → 64× cross-section
        let m = complex_real(1.5);
        let s1 = rayleigh_cross_section(1e-8, 500e-9, m).unwrap();
        let s2 = rayleigh_cross_section(2e-8, 500e-9, m).unwrap();
        assert!((s2 / s1 - 64.0).abs() < 0.1);
    }

    #[test]
    fn test_rayleigh_phase_function_forward_backward_equal() {
        // Rayleigh scattering is symmetric: P(θ) = P(π−θ)
        let p_fwd = rayleigh_phase_function(0.5);
        let p_bwd = rayleigh_phase_function(-0.5);
        assert!((p_fwd - p_bwd).abs() < TOL);
    }

    #[test]
    fn test_rayleigh_phase_function_90_deg() {
        // P(90°) = 3/4 (minimum, cos θ = 0)
        assert!((rayleigh_phase_function(0.0) - 0.75).abs() < TOL);
    }

    #[test]
    fn test_rayleigh_phase_function_forward() {
        // P(0°) = 3/2 (maximum, cos θ = 1)
        assert!((rayleigh_phase_function(1.0) - 1.5).abs() < TOL);
    }

    #[test]
    fn test_rayleigh_efficiency_positive() {
        let q = rayleigh_efficiency(0.1, complex_real(1.5)).unwrap();
        assert!(q > 0.0);
    }

    #[test]
    fn test_rayleigh_invalid_params() {
        assert!(rayleigh_cross_section(0.0, 500e-9, complex_real(1.5)).is_err());
        assert!(rayleigh_cross_section(1e-8, 0.0, complex_real(1.5)).is_err());
        assert!(rayleigh_efficiency(0.0, complex_real(1.5)).is_err());
    }

    #[test]
    fn test_rayleigh_angular_intensity_inverse_square() {
        let sigma = 1e-20;
        let cos_t = 0.5;
        let i1 = rayleigh_angular_intensity(sigma, 1.0, cos_t).unwrap();
        let i2 = rayleigh_angular_intensity(sigma, 2.0, cos_t).unwrap();
        assert!((i1 / i2 - 4.0).abs() < TOL);
    }

    #[test]
    fn test_rayleigh_angular_intensity_singularity() {
        assert!(rayleigh_angular_intensity(1e-20, 0.0, 0.5).is_err());
    }

    #[test]
    fn test_cross_section_from_efficiency() {
        let r = 1e-6;
        let q = 2.0;
        let sigma = cross_section_from_efficiency(q, r);
        let expected = q * std::f64::consts::PI * r * r;
        assert!((sigma - expected).abs() < 1e-30);
    }

    // ── Mie-Rayleigh consistency ──────────────────────────────────

    #[test]
    fn test_mie_rayleigh_cross_section_consistency() {
        // Cross-section from Mie Q_sca should match Rayleigh formula for small x
        let r = 10e-9;
        let wavelength = 500e-9;
        let x = 2.0 * std::f64::consts::PI * r / wavelength;
        let m = complex_real(1.5);

        let mie_sigma = cross_section_from_efficiency(mie(x, m).unwrap().q_sca, r);
        let ray_sigma = rayleigh_cross_section(r, wavelength, m).unwrap();

        assert!(
            (mie_sigma - ray_sigma).abs() / ray_sigma < 0.01,
            "Mie σ={mie_sigma} vs Rayleigh σ={ray_sigma}"
        );
    }
}
