//! Mie and Rayleigh scattering — cross-sections, efficiencies, phase functions.
//!
//! Mie theory gives the exact solution for scattering of a plane wave by a
//! homogeneous sphere. Rayleigh scattering is the small-particle limit (x ≪ 1).
//!
//! Size parameter: x = 2πr/λ where r is the sphere radius and λ is the
//! wavelength *in the surrounding medium*; the refractive index `m` is the
//! sphere's index relative to that medium.
//!
//! # Sign convention
//!
//! [`mie`] uses the **e^{+iωt}** (engineering) time convention: the outgoing
//! Riccati–Hankel function is ξₙ(x) = ψₙ(x) + iχₙ(x) = x hₙ⁽²⁾(x), so
//! ξ₀(x) = sin x + i cos x. In this convention an **absorbing** sphere has
//! **Im(m) < 0** (m = n − iκ, κ ≥ 0); `mie` rejects Im(m) > 0 (a gain medium,
//! or — far more often — an index written in the other convention).
//!
//! Bohren & Huffman (1983) and most optics texts use e^{−iωt} with
//! m = n + iκ. To use such data pass `m.conj()`. The returned efficiencies and
//! asymmetry parameter are convention-independent; the Mie coefficients
//! `a`, `b` are the complex conjugates of Bohren & Huffman's aₙ, bₙ.
//!
//! The Jones calculus in [`crate::polarization`] uses the opposite
//! (e^{−iωt}) convention — see that module's docs before mixing phases
//! between the two modules.

use crate::error::{BijliError, Result};
use crate::polarization::{Complex, complex_real, complex_zero};

// ── Mie scattering ────────────────────────────────────────────────

/// Largest size parameter x accepted by [`mie`].
///
/// The series needs ≈ x + 4x^{1/3} + 2 terms, so this bounds both run time
/// and the size of the returned coefficient vectors (≈ 10⁵ entries).
pub const MIE_MAX_SIZE_PARAMETER: f64 = 1.0e5;

/// Largest |m|·x accepted by [`mie`].
///
/// The logarithmic derivative Dₙ(mx) is computed by downward recurrence
/// from n ≈ max(n_max, |mx| + 4|mx|^{1/3}) + 15, so |mx| bounds the run time.
pub const MIE_MAX_MX: f64 = 1.0e7;

/// Below this size parameter ψ₁(x) is evaluated by its power series to avoid
/// the cancellation in sin x / x − cos x.
const PSI1_SERIES_THRESHOLD: f64 = 0.1;

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
    /// Mie coefficients a_n (n = 1…n_max), in the e^{+iωt} convention
    /// (complex conjugates of Bohren & Huffman's aₙ).
    pub a: Vec<Complex>,
    /// Mie coefficients b_n (n = 1…n_max), in the e^{+iωt} convention
    /// (complex conjugates of Bohren & Huffman's bₙ).
    pub b: Vec<Complex>,
}

/// Compute Mie scattering for a homogeneous sphere.
///
/// Uses the Bohren & Huffman (1983, App. A) formulation:
///
/// ```text
/// aₙ = [(Dₙ(mx)/m + n/x) ψₙ(x) − ψₙ₋₁(x)] / [(Dₙ(mx)/m + n/x) ξₙ(x) − ξₙ₋₁(x)]
/// bₙ = [(m Dₙ(mx) + n/x) ψₙ(x) − ψₙ₋₁(x)] / [(m Dₙ(mx) + n/x) ξₙ(x) − ξₙ₋₁(x)]
///
/// Q_sca = (2/x²) Σ (2n+1)(|aₙ|² + |bₙ|²)
/// Q_ext = (2/x²) Σ (2n+1) Re(aₙ + bₙ)
/// g     = (4/(x² Q_sca)) Σ [ n(n+2)/(n+1) Re(aₙa*ₙ₊₁ + bₙb*ₙ₊₁)
///                            + (2n+1)/(n(n+1)) Re(aₙb*ₙ) ]
/// ```
///
/// with the series truncated at Wiscombe's n_max = ⌈x + 4x^{1/3} + 2⌉ (at
/// least 4). Dₙ(mx) = ψₙ′(mx)/ψₙ(mx) is obtained by downward recurrence
/// started at max(n_max, |mx| + 4|mx|^{1/3}) + 15, which keeps it accurate
/// when |mx| ≫ x (high-index spheres). ψₙ and χₙ use upward recurrence; ψ₁ is
/// taken from its power series for x < 0.1.
///
/// # Conventions
///
/// e^{+iωt} time dependence: **absorption requires Im(m) < 0**
/// (m = n − iκ). For an index quoted as n + iκ (Bohren & Huffman, e^{−iωt})
/// pass `m.conj()`. See the [module docs](self).
///
/// # Accuracy
///
/// Q_ext, Q_sca and g agree with high-precision references to ~1e-12
/// relative for 1e-8 ≤ x ≤ 200. For x ≲ 1e-4 the asymmetry parameter (which
/// scales as x²) is accurate only in absolute terms (|Δg| ≲ 1e-15), as are the
/// tiny higher-order coefficients (b₁, a₂, …).
///
/// # Parameters
/// - `x`: size parameter 2πr/λ (finite, 0 < x ≤ [`MIE_MAX_SIZE_PARAMETER`])
/// - `m`: complex refractive index of the sphere relative to the medium
///   (finite, non-zero, Im(m) ≤ 0, |m|·x ≤ [`MIE_MAX_MX`])
///
/// # Errors
/// [`BijliError::InvalidParameter`] for out-of-range or non-finite `x`/`m`,
/// for Im(m) > 0, or if the series over/underflows (only for vanishingly
/// small x).
///
/// # Returns
/// [`MieResult`] with efficiencies, asymmetry parameter, and Mie coefficients.
pub fn mie(x: f64, m: Complex) -> Result<MieResult> {
    if !(x.is_finite() && x > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("size parameter x must be positive and finite, got {x}"),
        });
    }
    if x > MIE_MAX_SIZE_PARAMETER {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "size parameter x = {x} exceeds the supported maximum {MIE_MAX_SIZE_PARAMETER}"
            ),
        });
    }
    if !m.re.is_finite() || !m.im.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("refractive index must be finite, got {m}"),
        });
    }
    if m.re == 0.0 && m.im == 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: "refractive index cannot be zero".into(),
        });
    }
    if m.im > 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "Im(m) = {} > 0: bijli's Mie solver uses the e^{{+iωt}} convention, where absorbing \
                 media have Im(m) < 0 (m = n − iκ); for an index written as n + iκ pass m.conj()",
                m.im
            ),
        });
    }
    let mx_abs = m.re.hypot(m.im) * x;
    if mx_abs.is_nan() || mx_abs > MIE_MAX_MX {
        return Err(BijliError::InvalidParameter {
            reason: format!("|m|·x = {mx_abs} exceeds the supported maximum {MIE_MAX_MX}"),
        });
    }

    let mx = m * x;

    // Number of terms: Wiscombe's criterion. x ≤ 1e5, so the cast is exact.
    let n_max = ((x + 4.0 * x.cbrt() + 2.0).ceil() as usize).max(4);
    // Start of the Dₙ(mx) downward recurrence. It must sit well above both
    // n_max and |mx|: the start-value error only decays once n < |mx|, so the
    // |mx| branch needs an extra 4|mx|^{1/3} (Bohren & Huffman's
    // max(n_max, |mx|) + 15 is not enough for |mx| ≳ 50), while the n_max
    // branch keeps the classic +15. |mx| ≤ 1e7, so the cast is exact.
    let n_start = ((n_max as f64).max(mx_abs + 4.0 * mx_abs.cbrt()) + 15.0).ceil() as usize;

    tracing::debug!(
        x,
        m_re = m.re,
        m_im = m.im,
        n_max,
        n_start,
        "computing Mie scattering"
    );

    // Logarithmic derivative Dₙ(z) = ψₙ′(z)/ψₙ(z), z = mx, by downward
    // recurrence (stable): D_{n−1}(z) = n/z − 1/(Dₙ(z) + n/z), D_{n_start} = 0.
    // Only n ≤ n_max is stored.
    let inv_mx = complex_recip(mx)?;
    let mut d_mx = vec![complex_zero(); n_max + 1];
    let mut d_n = complex_zero();
    for n in (1..=n_start).rev() {
        let ratio = inv_mx * n as f64;
        d_n = ratio - complex_recip(d_n + ratio)?;
        if let Some(slot) = d_mx.get_mut(n - 1) {
            *slot = d_n;
        }
    }

    // Riccati–Bessel functions by upward recurrence:
    // ψₙ(x) = x jₙ(x), χₙ(x) = −x yₙ(x), ξₙ(x) = ψₙ(x) + iχₙ(x) = x hₙ⁽²⁾(x).
    let (sin_x, cos_x) = x.sin_cos();
    let mut psi_prev = sin_x; // ψ₀
    let mut psi_curr = psi1(x, sin_x, cos_x); // ψ₁
    let mut chi_prev = cos_x; // χ₀
    let mut chi_curr = cos_x / x + sin_x; // χ₁

    let inv_m = complex_recip(m)?;
    let mut a_coeffs = Vec::with_capacity(n_max);
    let mut b_coeffs = Vec::with_capacity(n_max);
    let mut sca_sum = 0.0;
    let mut ext_sum = 0.0;
    let mut g_sum = 0.0;
    let mut a_prev = complex_zero();
    let mut b_prev = complex_zero();

    for (n, &d_n) in d_mx.iter().enumerate().skip(1) {
        if n > 1 {
            let f = (2 * n - 1) as f64 / x;
            let psi_next = f * psi_curr - psi_prev;
            psi_prev = psi_curr;
            psi_curr = psi_next;
            let chi_next = f * chi_curr - chi_prev;
            chi_prev = chi_curr;
            chi_curr = chi_next;
        }
        let xi_curr = Complex::new(psi_curr, chi_curr);
        let xi_prev = Complex::new(psi_prev, chi_prev);
        let n_over_x = complex_real(n as f64 / x);

        // aₙ uses Dₙ/m + n/x, bₙ uses m Dₙ + n/x (Bohren & Huffman eq. 4.88).
        let a_factor = d_n * inv_m + n_over_x;
        let a_n = complex_div(
            a_factor * psi_curr - complex_real(psi_prev),
            a_factor * xi_curr - xi_prev,
        )?;

        let b_factor = m * d_n + n_over_x;
        let b_n = complex_div(
            b_factor * psi_curr - complex_real(psi_prev),
            b_factor * xi_curr - xi_prev,
        )?;

        let nf = n as f64;
        let weight = 2.0 * nf + 1.0;
        sca_sum += weight * (a_n.norm_sq() + b_n.norm_sq());
        ext_sum += weight * (a_n.re + b_n.re);

        // Asymmetry: (2n+1)/(n(n+1)) Re(aₙb*ₙ) for every n, plus the
        // (n−1)(n+1)/n Re(aₙ₋₁a*ₙ + bₙ₋₁b*ₙ) cross term for n > 1.
        g_sum += weight / (nf * (nf + 1.0)) * (a_n.re * b_n.re + a_n.im * b_n.im);
        if n > 1 {
            g_sum += (nf - 1.0) * (nf + 1.0) / nf
                * (a_prev.re * a_n.re
                    + a_prev.im * a_n.im
                    + b_prev.re * b_n.re
                    + b_prev.im * b_n.im);
        }

        a_prev = a_n;
        b_prev = b_n;
        a_coeffs.push(a_n);
        b_coeffs.push(b_n);
    }

    let scale = 2.0 / (x * x);
    let q_sca = scale * sca_sum;
    let q_ext = scale * ext_sum;
    let q_abs = q_ext - q_sca;
    // g = (4/(x² Q_sca)) Σ… = 2 Σ… / Σ(2n+1)(|a|²+|b|²); x² cancels.
    let g = if sca_sum > 0.0 {
        2.0 * g_sum / sca_sum
    } else {
        0.0
    };

    if !(q_sca.is_finite() && q_ext.is_finite() && g.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Mie series over/underflowed for x = {x}, m = {m}"),
        });
    }

    tracing::debug!(q_ext, q_sca, g, "Mie scattering computed");

    Ok(MieResult {
        q_sca,
        q_ext,
        q_abs,
        g,
        a: a_coeffs,
        b: b_coeffs,
    })
}

/// ψ₁(x) = sin x / x − cos x, via its Maclaurin series for small x where the
/// direct form cancels catastrophically (relative error ≈ 3ε/x²).
///
/// ψ₁(x) = x²/3 − x⁴/30 + x⁶/840 − x⁸/45360 + x¹⁰/3991680 − …
#[inline]
fn psi1(x: f64, sin_x: f64, cos_x: f64) -> f64 {
    if x < PSI1_SERIES_THRESHOLD {
        let x2 = x * x;
        x2 / 3.0 * (1.0 - x2 / 10.0 * (1.0 - x2 / 28.0 * (1.0 - x2 / 54.0 * (1.0 - x2 / 88.0))))
    } else {
        sin_x / x - cos_x
    }
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
/// Depends on m only through |K|², so it is independent of the sign
/// convention for Im(m).
#[inline]
pub fn rayleigh_cross_section(radius: f64, wavelength: f64, m: Complex) -> Result<f64> {
    if !(radius.is_finite() && radius > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("radius must be positive and finite, got {radius}"),
        });
    }
    if !(wavelength.is_finite() && wavelength > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive and finite, got {wavelength}"),
        });
    }

    let x = 2.0 * std::f64::consts::PI * radius / wavelength;
    let m2 = m * m;
    let k = complex_div(m2 - complex_real(1.0), m2 + complex_real(2.0))?;
    let x4 = x * x * x * x;

    Ok(8.0 / 3.0 * std::f64::consts::PI * radius * radius * x4 * k.norm_sq())
}

/// Rayleigh scattering efficiency Q_sca = σ/(πr²) = (8/3) x⁴ |K|².
#[inline]
pub fn rayleigh_efficiency(x: f64, m: Complex) -> Result<f64> {
    if !(x.is_finite() && x > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("size parameter must be positive and finite, got {x}"),
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

/// Rayleigh-scattered intensity ratio I/I₀ at distance r and angle θ for
/// unpolarized incident light.
///
/// ```text
/// I/I₀ = (σ / r²) · P(θ)/(4π) = 3σ (1 + cos²θ) / (16π r²)
///      = 8π⁴ α_v² (1 + cos²θ) / (λ⁴ r²)
/// ```
///
/// where σ is the total scattering cross-section (m², e.g. from
/// [`rayleigh_cross_section`]) and α_v = α/(4πε₀) is the polarizability
/// *volume* (m³); for a sphere of radius a, α_v = a³ (m²−1)/(m²+2).
#[inline]
pub fn rayleigh_angular_intensity(
    cross_section: f64,
    distance: f64,
    cos_theta: f64,
) -> Result<f64> {
    if distance.is_nan() {
        return Err(BijliError::InvalidParameter {
            reason: "distance must not be NaN".into(),
        });
    }
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(
        cross_section / (distance * distance) * rayleigh_phase_function(cos_theta)
            / (4.0 * std::f64::consts::PI),
    )
}

// ── Complex arithmetic helpers ────────────────────────────────────

/// Reciprocal 1/z: a single division when |z|² is a normal float, Smith's
/// algorithm (via [`complex_div`]) when |z|² would overflow or underflow.
#[inline]
fn complex_recip(z: Complex) -> Result<Complex> {
    let d = z.re * z.re + z.im * z.im;
    if d.is_normal() {
        let inv = 1.0 / d;
        Ok(Complex::new(z.re * inv, -z.im * inv))
    } else {
        complex_div(complex_real(1.0), z)
    }
}

/// a/b for complex numbers (Smith's algorithm — no intermediate |b|², so no
/// spurious over/underflow).
#[inline]
fn complex_div(a: Complex, b: Complex) -> Result<Complex> {
    if b.re == 0.0 && b.im == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "complex division by zero in scattering computation".into(),
        });
    }
    if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let den = b.re + b.im * r;
        Ok(Complex::new(
            (a.re + a.im * r) / den,
            (a.im - a.re * r) / den,
        ))
    } else {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        Ok(Complex::new(
            (a.re * r + a.im) / den,
            (a.im * r - a.re) / den,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Mie tests ─────────────────────────────────────────────────

    #[test]
    fn test_mie_nmax_dominated_large_x_matches_reference() {
        // Low-contrast large spheres (|mx| < n_max): the D_n start is set by
        // the n_max branch. References: 40-digit mpmath evaluation from
        // Bessel functions (independent of the recurrences used here).
        for (x, m, q_ref, g_ref) in [
            (100.0, 1.05, 2.379_520_657_428_609, 0.986_193_953_491_216_8),
            (60.0, 1.02, 2.120_134_468_004_568, 0.998_148_563_737_599),
        ] {
            let r = mie(x, complex_real(m)).unwrap();
            assert!(
                ((r.q_ext - q_ref) / q_ref).abs() < 1e-12,
                "x={x} m={m}: Q_ext {} vs {q_ref}",
                r.q_ext
            );
            assert!(
                ((r.g - g_ref) / g_ref).abs() < 1e-12,
                "x={x} m={m}: g {} vs {g_ref}",
                r.g
            );
        }
    }

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

    // ── Reference-value regression tests ──────────────────────────
    //
    // Reference values come from an independent 40-digit mpmath evaluation of
    // the Mie series using Bessel functions directly (no recurrences), in
    // Bohren & Huffman's e^{−iωt} convention (bijli input = conj of BH m).

    fn rel(a: f64, b: f64) -> f64 {
        ((a - b) / b).abs()
    }

    #[test]
    fn test_mie_bohren_huffman_sample() {
        // Bohren & Huffman App. A sample: r = 0.525 µm, λ = 0.6328 µm, m = 1.55.
        // BH tabulate Q_ext = Q_sca = 3.10543, g = 0.63314.
        let x = 2.0 * std::f64::consts::PI * 0.525 / 0.6328;
        let r = mie(x, complex_real(1.55)).unwrap();
        assert!((r.q_ext - 3.10543).abs() < 5e-6, "Q_ext = {}", r.q_ext);
        assert!((r.g - 0.63314).abs() < 5e-6, "g = {}", r.g);
        // 40-digit reference: Q_ext = 3.10542553146588, g = 0.633136758040895
        assert!(rel(r.q_ext, 3.105_425_531_465_88) < 1e-11);
        assert!(rel(r.q_sca, 3.105_425_531_465_88) < 1e-11);
        assert!(rel(r.g, 0.633_136_758_040_895) < 1e-11);
    }

    #[test]
    fn test_mie_high_index_large_mx() {
        // |mx| = 80 ≫ n_max ≈ 33: the old D_n start (n_max + 15) gave
        // Q_ext = 2.0758, g = 0.4230.
        let r = mie(20.0, complex_real(4.0)).unwrap();
        assert!(
            rel(r.q_ext, 2.089_041_014_033_205) < 1e-10,
            "Q_ext = {}",
            r.q_ext
        );
        assert!(rel(r.g, 0.500_862_780_527_793) < 1e-10, "g = {}", r.g);

        // |mx| = 175: old code gave Q_ext = 2.2122.
        let r = mie(50.0, complex_real(3.5)).unwrap();
        assert!(
            rel(r.q_ext, 2.080_685_067_202_63) < 1e-10,
            "Q_ext = {}",
            r.q_ext
        );
        assert!(rel(r.g, 0.514_978_310_390_713_1) < 1e-10, "g = {}", r.g);
    }

    #[test]
    fn test_mie_large_sphere() {
        // x = 200, |mx| = 300: old code gave Q_ext = 2.0150.
        let r = mie(200.0, complex_real(1.5)).unwrap();
        assert!(
            rel(r.q_ext, 2.092_092_687_530_17) < 1e-10,
            "Q_ext = {}",
            r.q_ext
        );
        assert!(rel(r.g, 0.821_956_642_316_233_5) < 1e-10, "g = {}", r.g);
    }

    #[test]
    fn test_mie_absorbing_reference() {
        // BH m = 1.33 + 0.1i  ↔  bijli m = 1.33 − 0.1i
        let r = mie(10.0, Complex::new(1.33, -0.1)).unwrap();
        assert!(rel(r.q_ext, 2.398_526_237_960_563) < 1e-10);
        assert!(rel(r.q_sca, 1.187_346_283_695_99) < 1e-10);
        assert!(rel(r.q_abs, 1.211_179_954_264_573) < 1e-10);
        assert!(rel(r.g, 0.919_321_731_299_517) < 1e-10);

        // Metal-like, strongly absorbing: BH m = 0.2 + 3i
        let r = mie(5.0, Complex::new(0.2, -3.0)).unwrap();
        assert!(rel(r.q_ext, 3.191_339_268_251_32) < 1e-10);
        assert!(rel(r.q_sca, 2.982_937_412_023_67) < 1e-10);
        assert!(rel(r.g, 0.551_568_886_029_315_8) < 1e-10);
    }

    #[test]
    fn test_mie_coefficients_are_bh_conjugates() {
        // BH (m = 1.5 + 0.1i, x = 1): a₁ = 0.0682287821494085 − 0.17068948273117i
        //                              b₁ = 0.00864512703725183 − 0.0272424021476631i
        let r = mie(1.0, Complex::new(1.5, -0.1)).unwrap();
        assert!((r.a[0].re - 0.068_228_782_149_408_5).abs() < 1e-13);
        assert!((r.a[0].im - 0.170_689_482_731_17).abs() < 1e-13);
        assert!((r.b[0].re - 0.008_645_127_037_251_83).abs() < 1e-13);
        assert!((r.b[0].im - 0.027_242_402_147_663_1).abs() < 1e-13);
        assert!(rel(r.q_abs, 0.273_630_438_032_15) < 1e-10);
    }

    #[test]
    fn test_mie_tiny_sphere_matches_series() {
        // x = 1e-8: ψ₁ = sin x/x − cos x cancels completely; the old code
        // returned Q_sca ≈ 1.17e-31 (50× the true value).
        // Reference 2.30680507497117e-33 (= (8/3) x⁴ |K|² to O(x²)).
        let r = mie(1e-8, complex_real(1.5)).unwrap();
        assert!(
            rel(r.q_sca, 2.306_805_074_971_17e-33) < 1e-12,
            "Q_sca = {}",
            r.q_sca
        );
        assert!(rel(r.q_ext, 2.306_805_074_971_17e-33) < 1e-12);
        assert!(r.g.abs() < 1e-15);
        let ray = rayleigh_efficiency(1e-8, complex_real(1.5)).unwrap();
        assert!(rel(r.q_sca, ray) < 1e-12);

        // Small absorbing sphere: BH m = 1.5 + 0.01i, x = 1e-3
        let r = mie(1e-3, Complex::new(1.5, -0.01)).unwrap();
        assert!(rel(r.q_ext, 1.993_075_206_713_3e-5) < 1e-10);
        assert!(rel(r.q_sca, 2.307_758_494_063_458e-13) < 1e-10);
    }

    #[test]
    fn test_mie_psi1_series_continuity() {
        // Results must be smooth across the ψ₁ series threshold (x = 0.1).
        // References: Q_sca(0.0999) = 2.299186415536038e-5,
        //             Q_sca(0.1001) = 2.317660029332435e-5.
        let lo = mie(0.0999, complex_real(1.5)).unwrap();
        let hi = mie(0.1001, complex_real(1.5)).unwrap();
        assert!(rel(lo.q_sca, 2.299_186_415_536_038e-5) < 1e-12);
        assert!(rel(hi.q_sca, 2.317_660_029_332_435e-5) < 1e-12);
    }

    #[test]
    fn test_mie_rejects_gain_convention() {
        // Im(m) > 0 is gain in the e^{+iωt} convention — almost always a
        // convention mix-up. Old code returned Q_abs < 0 silently.
        assert!(mie(1.0, Complex::new(1.5, 0.1)).is_err());
        assert!(mie(1.0, Complex::new(1.5, 0.1).conj()).is_ok());
    }

    #[test]
    fn test_mie_rejects_non_finite_and_huge() {
        let m = complex_real(1.5);
        assert!(mie(f64::NAN, m).is_err());
        assert!(mie(f64::INFINITY, m).is_err());
        assert!(mie(1e10, m).is_err());
        assert!(mie(MIE_MAX_SIZE_PARAMETER * 1.0001, m).is_err());
        assert!(mie(1.0, Complex::new(f64::NAN, 0.0)).is_err());
        assert!(mie(1.0, Complex::new(1.5, f64::NAN)).is_err());
        assert!(mie(1.0, Complex::new(f64::INFINITY, 0.0)).is_err());
        assert!(mie(1.0, complex_zero()).is_err());
        assert!(mie(10.0, complex_real(2e6)).is_err()); // |mx| > MIE_MAX_MX
        assert!(mie(1e-320, m).is_err()); // series overflows → Err, not NaN
    }

    #[test]
    fn test_rayleigh_angular_intensity_closed_form() {
        // I/I₀ = 8π⁴ α_v² (1 + cos²θ)/(λ⁴ r²), α_v = a³ (m²−1)/(m²+2)
        let a = 10e-9;
        let lambda: f64 = 500e-9;
        let dist = 2.0;
        let cos_t = 0.3;
        let k_fac = (1.5_f64 * 1.5 - 1.0) / (1.5 * 1.5 + 2.0);
        let alpha_v = a * a * a * k_fac;
        let pi = std::f64::consts::PI;
        let expected = 8.0 * pi.powi(4) * alpha_v * alpha_v * (1.0 + cos_t * cos_t)
            / (lambda.powi(4) * dist * dist);
        let sigma = rayleigh_cross_section(a, lambda, complex_real(1.5)).unwrap();
        let got = rayleigh_angular_intensity(sigma, dist, cos_t).unwrap();
        assert!(rel(got, expected) < 1e-12, "{got} vs {expected}");
    }

    #[test]
    fn test_rayleigh_nan_inputs_rejected() {
        let m = complex_real(1.5);
        assert!(rayleigh_cross_section(f64::NAN, 500e-9, m).is_err());
        assert!(rayleigh_cross_section(1e-8, f64::NAN, m).is_err());
        assert!(rayleigh_efficiency(f64::NAN, m).is_err());
        assert!(rayleigh_angular_intensity(1e-20, f64::NAN, 0.5).is_err());
    }
}
