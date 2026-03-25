//! Gaussian beams, ABCD ray transfer matrices, resonator stability,
//! and higher-order Hermite-Gaussian / Laguerre-Gaussian modes.
//!
//! All SI units. Wavelength in meters, distances in meters.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

// ── Gaussian beam parameter ───────────────────────────────────────

/// Gaussian beam parameters at a reference plane.
///
/// A fundamental Gaussian beam is fully described by its wavelength λ and
/// complex beam parameter q(z) = z - z_waist + i z_R, where z_R = πw₀²/λ.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GaussianBeam {
    /// Wavelength λ (m).
    pub wavelength: f64,
    /// Beam waist radius w₀ (m) — the 1/e² intensity radius at the focus.
    pub waist: f64,
    /// Rayleigh range z_R = πw₀²/λ (m).
    pub rayleigh_range: f64,
}

impl GaussianBeam {
    /// Create a Gaussian beam from wavelength and waist radius.
    #[inline]
    pub fn new(wavelength: f64, waist: f64) -> Result<Self> {
        if wavelength <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("wavelength must be positive, got {wavelength}"),
            });
        }
        if waist <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("beam waist must be positive, got {waist}"),
            });
        }
        let rayleigh_range = std::f64::consts::PI * waist * waist / wavelength;
        Ok(Self {
            wavelength,
            waist,
            rayleigh_range,
        })
    }

    /// Complex beam parameter q(z) = (z - z_waist) + i z_R.
    ///
    /// `z` is the distance from the waist along the propagation axis.
    #[inline]
    #[must_use]
    pub fn q(&self, z: f64) -> [f64; 2] {
        [z, self.rayleigh_range]
    }

    /// Spot size (1/e² intensity radius) at distance z from waist.
    ///
    /// w(z) = w₀ √(1 + (z/z_R)²)
    #[inline]
    #[must_use]
    pub fn spot_size(&self, z: f64) -> f64 {
        let ratio = z / self.rayleigh_range;
        self.waist * (1.0 + ratio * ratio).sqrt()
    }

    /// Radius of curvature of the wavefront at distance z from waist.
    ///
    /// R(z) = z(1 + (z_R/z)²). At z=0 (waist), R → ∞ (flat wavefront).
    #[inline]
    pub fn radius_of_curvature(&self, z: f64) -> Result<f64> {
        if z.abs() < 1e-30 {
            return Err(BijliError::InvalidParameter {
                reason: "radius of curvature is infinite at the beam waist".into(),
            });
        }
        let ratio = self.rayleigh_range / z;
        Ok(z * (1.0 + ratio * ratio))
    }

    /// Gouy phase at distance z from waist: ζ(z) = arctan(z/z_R).
    #[inline]
    #[must_use]
    pub fn gouy_phase(&self, z: f64) -> f64 {
        (z / self.rayleigh_range).atan()
    }

    /// Far-field half-angle divergence: θ = λ/(πw₀).
    #[inline]
    #[must_use]
    pub fn divergence(&self) -> f64 {
        self.wavelength / (std::f64::consts::PI * self.waist)
    }

    /// Peak intensity relative to total power: I₀ = 2P/(πw²).
    ///
    /// Returns the factor 2/(πw(z)²) that multiplies total power P.
    #[inline]
    #[must_use]
    pub fn peak_intensity_factor(&self, z: f64) -> f64 {
        let w = self.spot_size(z);
        2.0 / (std::f64::consts::PI * w * w)
    }

    /// Confocal parameter (twice the Rayleigh range): b = 2z_R.
    #[inline]
    #[must_use]
    pub fn confocal_parameter(&self) -> f64 {
        2.0 * self.rayleigh_range
    }

    /// Depth of focus (distance over which spot size grows by √2): DOF = 2z_R.
    #[inline]
    #[must_use]
    pub fn depth_of_focus(&self) -> f64 {
        2.0 * self.rayleigh_range
    }
}

// ── ABCD ray transfer matrices ────────────────────────────────────

/// ABCD ray transfer matrix — 2×2 real matrix for paraxial ray tracing.
///
/// A ray is described by [y, θ] where y is the height from the axis
/// and θ is the angle (paraxial, so θ ≈ tan θ ≈ sin θ).
///
/// ```text
/// [y']   [A  B] [y]
/// [θ'] = [C  D] [θ]
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AbcdMatrix {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl AbcdMatrix {
    /// Create an ABCD matrix from components.
    #[inline]
    #[must_use]
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }

    /// Identity matrix (no optical element).
    #[inline]
    #[must_use]
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Free-space propagation over distance d.
    #[inline]
    #[must_use]
    pub fn free_space(distance: f64) -> Self {
        Self::new(1.0, distance, 0.0, 1.0)
    }

    /// Thin lens with focal length f.
    ///
    /// Positive f = converging, negative f = diverging.
    #[inline]
    pub fn thin_lens(focal_length: f64) -> Result<Self> {
        if focal_length.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "focal length cannot be zero".into(),
            });
        }
        Ok(Self::new(1.0, 0.0, -1.0 / focal_length, 1.0))
    }

    /// Curved mirror with radius of curvature R.
    ///
    /// Concave (R > 0) focuses, convex (R < 0) diverges.
    /// Equivalent to a thin lens with f = R/2.
    #[inline]
    pub fn curved_mirror(radius: f64) -> Result<Self> {
        if radius.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "mirror radius cannot be zero".into(),
            });
        }
        Ok(Self::new(1.0, 0.0, -2.0 / radius, 1.0))
    }

    /// Flat interface between media with refractive indices n1 → n2.
    #[inline]
    pub fn flat_interface(n1: f64, n2: f64) -> Result<Self> {
        if n2.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "refractive index n₂ cannot be zero".into(),
            });
        }
        Ok(Self::new(1.0, 0.0, 0.0, n1 / n2))
    }

    /// Curved interface (radius R) between media n1 → n2.
    ///
    /// R > 0 when center of curvature is on the transmission side.
    #[inline]
    pub fn curved_interface(n1: f64, n2: f64, radius: f64) -> Result<Self> {
        if n2.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "refractive index n₂ cannot be zero".into(),
            });
        }
        if radius.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "interface radius cannot be zero".into(),
            });
        }
        Ok(Self::new(1.0, 0.0, (n1 - n2) / (n2 * radius), n1 / n2))
    }

    /// Apply this matrix to a ray [y, θ].
    #[inline]
    #[must_use]
    pub fn apply_ray(&self, y: f64, theta: f64) -> (f64, f64) {
        (self.a * y + self.b * theta, self.c * y + self.d * theta)
    }

    /// Compose two ABCD matrices: self × other (self applied after other).
    #[inline]
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        Self {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }

    /// Determinant AD - BC. For a system in the same medium, det = 1.
    /// For a system spanning media n1 → n2, det = n1/n2.
    #[inline]
    #[must_use]
    pub fn determinant(&self) -> f64 {
        self.a * self.d - self.b * self.c
    }

    /// Propagate a complex beam parameter through this ABCD matrix.
    ///
    /// q_out = (A q_in + B) / (C q_in + D)
    ///
    /// `q` is [Re(q), Im(q)] = [z_from_waist, z_R].
    /// Returns the transformed [Re(q'), Im(q')].
    #[inline]
    pub fn propagate_beam(&self, q: [f64; 2]) -> Result<[f64; 2]> {
        // q = q[0] + i*q[1]
        // numerator = A*q + B = (A*q[0] + B) + i*(A*q[1])
        let num_re = self.a * q[0] + self.b;
        let num_im = self.a * q[1];
        // denominator = C*q + D = (C*q[0] + D) + i*(C*q[1])
        let den_re = self.c * q[0] + self.d;
        let den_im = self.c * q[1];
        let den_sq = den_re * den_re + den_im * den_im;
        if den_sq < 1e-60 {
            return Err(BijliError::DivisionByZero {
                context: "ABCD beam propagation denominator is zero".into(),
            });
        }
        Ok([
            (num_re * den_re + num_im * den_im) / den_sq,
            (num_im * den_re - num_re * den_im) / den_sq,
        ])
    }

    /// Extract beam parameters after propagation through this matrix.
    ///
    /// Given a `GaussianBeam` at its waist, propagate through this system
    /// and return the new waist size and waist location (relative to output plane).
    pub fn transform_beam(&self, beam: &GaussianBeam) -> Result<GaussianBeam> {
        let q_in = beam.q(0.0); // at waist: q = i*z_R
        let q_out = self.propagate_beam(q_in)?;
        // Im(q_out) = z_R_new = π w₀_new² / λ
        let z_r_new = q_out[1];
        if z_r_new <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: "transformed beam has non-physical Rayleigh range".into(),
            });
        }
        let w0_new = (z_r_new * beam.wavelength / std::f64::consts::PI).sqrt();
        GaussianBeam::new(beam.wavelength, w0_new)
    }
}

impl Default for AbcdMatrix {
    #[inline]
    fn default() -> Self {
        Self::identity()
    }
}

// ── Resonator stability ───────────────────────────────────────────

/// Resonator stability analysis for a two-mirror cavity.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ResonatorStability {
    /// g-parameter for mirror 1: g₁ = 1 − L/R₁.
    pub g1: f64,
    /// g-parameter for mirror 2: g₂ = 1 − L/R₂.
    pub g2: f64,
    /// Product g₁g₂ (stable if 0 ≤ g₁g₂ ≤ 1).
    pub g_product: f64,
    /// Whether the cavity is stable.
    pub is_stable: bool,
}

/// Analyze stability of a two-mirror resonator.
///
/// # Parameters
/// - `length`: cavity length L (m)
/// - `r1`: radius of curvature of mirror 1 (m, positive = concave)
/// - `r2`: radius of curvature of mirror 2 (m, positive = concave)
#[inline]
pub fn resonator_stability(length: f64, r1: f64, r2: f64) -> Result<ResonatorStability> {
    if length <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("cavity length must be positive, got {length}"),
        });
    }
    if r1.abs() < 1e-30 || r2.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "mirror radius of curvature cannot be zero".into(),
        });
    }
    let g1 = 1.0 - length / r1;
    let g2 = 1.0 - length / r2;
    let g_product = g1 * g2;
    let is_stable = (0.0..=1.0).contains(&g_product);
    Ok(ResonatorStability {
        g1,
        g2,
        g_product,
        is_stable,
    })
}

/// Beam waist of the fundamental mode inside a stable resonator.
///
/// w₀² = (Lλ/π) √(g₁g₂(1−g₁g₂)) / (g₁+g₂−2g₁g₂)²)
///
/// Returns the waist radius in meters, or error if cavity is unstable.
pub fn resonator_beam_waist(length: f64, r1: f64, r2: f64, wavelength: f64) -> Result<f64> {
    let stab = resonator_stability(length, r1, r2)?;
    if !stab.is_stable {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "cavity is unstable (g₁g₂ = {:.4}), no confined mode exists",
                stab.g_product
            ),
        });
    }
    if wavelength <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive, got {wavelength}"),
        });
    }

    let g1 = stab.g1;
    let g2 = stab.g2;
    let gp = stab.g_product;

    let denom = (g1 + g2 - 2.0 * gp).powi(2);
    if denom < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "degenerate cavity geometry".into(),
        });
    }

    let w0_sq = length * wavelength / std::f64::consts::PI * (gp * (1.0 - gp) / denom).abs().sqrt();
    Ok(w0_sq.sqrt())
}

/// Transverse mode spacing (frequency difference between adjacent modes).
///
/// Δν_transverse = (Δν_FSR / π) × arccos(±√(g₁g₂))
///
/// `fsr` is the free spectral range c/(2L).
#[inline]
pub fn transverse_mode_spacing(fsr: f64, g1: f64, g2: f64) -> Result<f64> {
    let gp = g1 * g2;
    if !(0.0..=1.0).contains(&gp) {
        return Err(BijliError::InvalidParameter {
            reason: format!("cavity is unstable (g₁g₂ = {gp:.4})"),
        });
    }
    Ok(fsr / std::f64::consts::PI * gp.sqrt().acos())
}

// ── Hermite-Gaussian modes ────────────────────────────────────────

/// Hermite-Gaussian beam amplitude HG_mn at position (x, y, z).
///
/// Returns the normalized field amplitude (dimensionless, relative to peak of TEM₀₀).
///
/// # Parameters
/// - `beam`: the fundamental Gaussian beam parameters
/// - `m`, `n`: mode indices (m for x, n for y)
/// - `x`, `y`: transverse coordinates (m)
/// - `z`: distance from waist (m)
pub fn hermite_gaussian(beam: &GaussianBeam, m: u32, n: u32, x: f64, y: f64, z: f64) -> f64 {
    let w = beam.spot_size(z);
    let w0 = beam.waist;

    // Normalized coordinates
    let sqrt2_x = std::f64::consts::SQRT_2 * x / w;
    let sqrt2_y = std::f64::consts::SQRT_2 * y / w;

    // Hermite polynomials
    let hm = hermite_polynomial(m, sqrt2_x);
    let hn = hermite_polynomial(n, sqrt2_y);

    // Gaussian envelope
    let r_sq = x * x + y * y;
    let envelope = (-r_sq / (w * w)).exp();

    // Amplitude scaling: w₀/w
    let scale = w0 / w;

    // Gouy phase contribution: (m + n + 1) × ζ(z)
    // (phase is real part of complex amplitude — for intensity, this cancels)
    // For amplitude we include it
    let gouy = (m + n + 1) as f64 * beam.gouy_phase(z);

    scale * hm * hn * envelope * gouy.cos()
}

/// Laguerre-Gaussian beam intensity profile LG_pl at position (r, z).
///
/// Returns the normalized intensity (dimensionless, relative to peak of TEM₀₀).
///
/// # Parameters
/// - `beam`: the fundamental Gaussian beam parameters
/// - `p`: radial mode index (number of radial nodes)
/// - `l`: azimuthal mode index (orbital angular momentum)
/// - `r`: radial distance from axis (m)
/// - `z`: distance from waist (m)
#[must_use]
pub fn laguerre_gaussian_intensity(beam: &GaussianBeam, p: u32, l: i32, r: f64, z: f64) -> f64 {
    let w = beam.spot_size(z);
    let w0 = beam.waist;

    let rho = std::f64::consts::SQRT_2 * r / w;
    let rho_sq = rho * rho;

    // Generalized Laguerre polynomial L_p^|l|(2r²/w²)
    let lp = laguerre_polynomial(p, l.unsigned_abs(), rho_sq);

    // Intensity profile (no azimuthal dependence when computing |E|²)
    let envelope = (-2.0 * r * r / (w * w)).exp();
    let scale = (w0 * w0) / (w * w);

    scale * rho_sq.powi(l.unsigned_abs() as i32) * lp * lp * envelope
}

// ── Polynomial helpers ────────────────────────────────────────────

/// Hermite polynomial H_n(x) via recurrence: H₀=1, H₁=2x, H_n = 2x H_{n-1} - 2(n-1) H_{n-2}.
#[must_use]
fn hermite_polynomial(n: u32, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut h_prev = 1.0;
            let mut h_curr = 2.0 * x;
            for k in 2..=n {
                let h_next = 2.0 * x * h_curr - 2.0 * (k - 1) as f64 * h_prev;
                h_prev = h_curr;
                h_curr = h_next;
            }
            h_curr
        }
    }
}

/// Generalized Laguerre polynomial L_p^l(x) via recurrence.
///
/// L₀ˡ = 1, L₁ˡ = 1+l-x, L_pˡ = ((2p+l-1-x) L_{p-1} - (p+l-1) L_{p-2}) / p
#[must_use]
fn laguerre_polynomial(p: u32, l: u32, x: f64) -> f64 {
    match p {
        0 => 1.0,
        1 => 1.0 + l as f64 - x,
        _ => {
            let mut l_prev = 1.0;
            let mut l_curr = 1.0 + l as f64 - x;
            for k in 2..=p {
                let kf = k as f64;
                let lf = l as f64;
                let l_next = ((2.0 * kf + lf - 1.0 - x) * l_curr - (kf + lf - 1.0) * l_prev) / kf;
                l_prev = l_curr;
                l_curr = l_next;
            }
            l_curr
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Gaussian beam tests ───────────────────────────────────────

    #[test]
    fn test_beam_creation() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        // z_R = πw₀²/λ ≈ 0.0295 m
        let expected_zr = std::f64::consts::PI * (100e-6_f64).powi(2) / 1064e-9;
        assert!((beam.rayleigh_range - expected_zr).abs() / expected_zr < TOL);
    }

    #[test]
    fn test_beam_invalid() {
        assert!(GaussianBeam::new(0.0, 100e-6).is_err());
        assert!(GaussianBeam::new(1064e-9, -1.0).is_err());
    }

    #[test]
    fn test_spot_size_at_waist() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        assert!((beam.spot_size(0.0) - beam.waist).abs() < 1e-20);
    }

    #[test]
    fn test_spot_size_at_rayleigh() {
        // w(z_R) = w₀ √2
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let w_zr = beam.spot_size(beam.rayleigh_range);
        assert!((w_zr / beam.waist - std::f64::consts::SQRT_2).abs() < TOL);
    }

    #[test]
    fn test_spot_size_symmetric() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        assert!((beam.spot_size(0.01) - beam.spot_size(-0.01)).abs() < 1e-20);
    }

    #[test]
    fn test_radius_of_curvature_far_field() {
        // For z >> z_R, R(z) ≈ z
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let z = 100.0 * beam.rayleigh_range;
        let r = beam.radius_of_curvature(z).unwrap();
        assert!((r / z - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_radius_of_curvature_at_waist_fails() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        assert!(beam.radius_of_curvature(0.0).is_err());
    }

    #[test]
    fn test_gouy_phase_at_rayleigh() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let gouy = beam.gouy_phase(beam.rayleigh_range);
        assert!((gouy - std::f64::consts::FRAC_PI_4).abs() < TOL);
    }

    #[test]
    fn test_divergence() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let theta = beam.divergence();
        let expected = 1064e-9 / (std::f64::consts::PI * 100e-6);
        assert!((theta - expected).abs() / expected < TOL);
    }

    #[test]
    fn test_confocal_equals_dof() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        assert!((beam.confocal_parameter() - beam.depth_of_focus()).abs() < 1e-30);
    }

    // ── ABCD matrix tests ─────────────────────────────────────────

    #[test]
    fn test_identity_ray() {
        let (y, theta) = AbcdMatrix::identity().apply_ray(1.0, 0.1);
        assert!((y - 1.0).abs() < TOL);
        assert!((theta - 0.1).abs() < TOL);
    }

    #[test]
    fn test_free_space_propagation() {
        let m = AbcdMatrix::free_space(1.0);
        let (y, theta) = m.apply_ray(0.0, 0.1);
        // y = 0 + 1.0 * 0.1 = 0.1, theta unchanged
        assert!((y - 0.1).abs() < TOL);
        assert!((theta - 0.1).abs() < TOL);
    }

    #[test]
    fn test_thin_lens_focuses() {
        let lens = AbcdMatrix::thin_lens(0.1).unwrap();
        // Parallel ray at height 1: (1, 0) → (1, -10)
        let (y, theta) = lens.apply_ray(1.0, 0.0);
        assert!((y - 1.0).abs() < TOL);
        assert!((theta - (-10.0)).abs() < TOL);
    }

    #[test]
    fn test_thin_lens_zero_focal() {
        assert!(AbcdMatrix::thin_lens(0.0).is_err());
    }

    #[test]
    fn test_curved_mirror_equiv_lens() {
        // Curved mirror with R = 2f has same ABCD as thin lens with f
        let f = 0.5;
        let lens = AbcdMatrix::thin_lens(f).unwrap();
        let mirror = AbcdMatrix::curved_mirror(2.0 * f).unwrap();
        assert!((lens.c - mirror.c).abs() < TOL);
    }

    #[test]
    fn test_flat_interface() {
        let m = AbcdMatrix::flat_interface(1.0, 1.5).unwrap();
        // Snell's law in paraxial: θ₂ = (n₁/n₂)θ₁
        let (y, theta) = m.apply_ray(0.0, 0.3);
        assert!(y.abs() < TOL);
        assert!((theta - 0.3 / 1.5).abs() < TOL);
    }

    #[test]
    fn test_compose_free_space() {
        // Two free-space propagations = one longer propagation
        let m1 = AbcdMatrix::free_space(1.0);
        let m2 = AbcdMatrix::free_space(2.0);
        let composed = m2.compose(&m1);
        let direct = AbcdMatrix::free_space(3.0);
        assert!((composed.b - direct.b).abs() < TOL);
    }

    #[test]
    fn test_determinant_free_space() {
        assert!((AbcdMatrix::free_space(1.0).determinant() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_determinant_interface() {
        let m = AbcdMatrix::flat_interface(1.0, 1.5).unwrap();
        assert!((m.determinant() - 1.0 / 1.5).abs() < TOL);
    }

    #[test]
    fn test_beam_propagation_free_space() {
        // Free-space propagation: waist doesn't change, only position
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let m = AbcdMatrix::free_space(0.01);
        let q_out = m.propagate_beam(beam.q(0.0)).unwrap();
        // Im(q) = z_R should be unchanged
        assert!((q_out[1] - beam.rayleigh_range).abs() / beam.rayleigh_range < TOL);
        // Re(q) should be 0.01 (we moved 1cm from waist)
        assert!((q_out[0] - 0.01).abs() < TOL);
    }

    #[test]
    fn test_beam_transform_through_lens() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let lens = AbcdMatrix::thin_lens(0.1).unwrap();
        let new_beam = lens.transform_beam(&beam).unwrap();
        // Lens at waist doesn't change waist size (R=∞ at waist → still ∞ after lens for w₀ << f)
        // but for very tight focus this holds approximately
        assert!(new_beam.waist > 0.0);
        assert!(new_beam.wavelength == beam.wavelength);
    }

    // ── Resonator stability tests ─────────────────────────────────

    #[test]
    fn test_symmetric_confocal() {
        // Symmetric confocal: R₁ = R₂ = L → g₁ = g₂ = 0 → g₁g₂ = 0 (marginally stable)
        let stab = resonator_stability(0.5, 0.5, 0.5).unwrap();
        assert!(stab.g1.abs() < TOL);
        assert!(stab.is_stable);
    }

    #[test]
    fn test_planar_cavity() {
        // Planar cavity: R₁ = R₂ = ∞ → g₁ = g₂ = 1 → g₁g₂ = 1 (marginally stable)
        let stab = resonator_stability(0.1, 1e10, 1e10).unwrap();
        assert!((stab.g_product - 1.0).abs() < TOL);
        assert!(stab.is_stable);
    }

    #[test]
    fn test_concentric_unstable() {
        // Concentric: R₁ = R₂ = L/2 → g₁ = g₂ = -1 → g₁g₂ = 1 (marginally stable, edge)
        // Just past concentric → unstable
        let stab = resonator_stability(1.0, 0.4, 0.4).unwrap();
        assert!(!stab.is_stable);
    }

    #[test]
    fn test_stable_cavity() {
        // R₁ = R₂ = 2L → g₁ = g₂ = 0.5 → g₁g₂ = 0.25
        let stab = resonator_stability(0.1, 0.2, 0.2).unwrap();
        assert!((stab.g_product - 0.25).abs() < TOL);
        assert!(stab.is_stable);
    }

    #[test]
    fn test_resonator_beam_waist() {
        let w = resonator_beam_waist(0.1, 0.2, 0.2, 1064e-9).unwrap();
        assert!(w > 0.0);
        assert!(w < 1e-3); // should be sub-mm for typical laser cavities
    }

    #[test]
    fn test_resonator_beam_waist_unstable_fails() {
        assert!(resonator_beam_waist(1.0, 0.4, 0.4, 1064e-9).is_err());
    }

    #[test]
    fn test_transverse_mode_spacing() {
        let fsr = 1.5e9; // 1.5 GHz
        let spacing = transverse_mode_spacing(fsr, 0.5, 0.5).unwrap();
        // For g₁g₂ = 0.25, arccos(0.5) = π/3, spacing = FSR/3
        assert!((spacing - fsr / 3.0).abs() / fsr < TOL);
    }

    #[test]
    fn test_transverse_mode_spacing_unstable_fails() {
        assert!(transverse_mode_spacing(1.5e9, -2.0, -2.0).is_err());
    }

    #[test]
    fn test_resonator_invalid_length() {
        assert!(resonator_stability(0.0, 1.0, 1.0).is_err());
        assert!(resonator_stability(-1.0, 1.0, 1.0).is_err());
    }

    // ── HG / LG mode tests ───────────────────────────────────────

    #[test]
    fn test_hermite_polynomial_values() {
        assert!((hermite_polynomial(0, 1.0) - 1.0).abs() < TOL);
        assert!((hermite_polynomial(1, 1.0) - 2.0).abs() < TOL);
        // H₂(x) = 4x² - 2
        assert!((hermite_polynomial(2, 1.0) - 2.0).abs() < TOL);
        // H₃(x) = 8x³ - 12x, H₃(1) = -4
        assert!((hermite_polynomial(3, 1.0) - (-4.0)).abs() < TOL);
    }

    #[test]
    fn test_laguerre_polynomial_values() {
        // L₀⁰ = 1
        assert!((laguerre_polynomial(0, 0, 1.0) - 1.0).abs() < TOL);
        // L₁⁰ = 1 - x
        assert!((laguerre_polynomial(1, 0, 0.5) - 0.5).abs() < TOL);
        // L₂⁰ = 1 - 2x + x²/2, L₂⁰(1) = 1 - 2 + 0.5 = -0.5
        assert!((laguerre_polynomial(2, 0, 1.0) - (-0.5)).abs() < TOL);
        // L₁¹(x) = 2 - x, L₁¹(1) = 1
        assert!((laguerre_polynomial(1, 1, 1.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_hg00_is_gaussian() {
        // HG₀₀ at waist should be a simple Gaussian peak at center
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let peak = hermite_gaussian(&beam, 0, 0, 0.0, 0.0, 0.0);
        assert!((peak - 1.0).abs() < TOL); // normalized to 1 at center of waist

        // Off-axis should be less
        let off = hermite_gaussian(&beam, 0, 0, 50e-6, 0.0, 0.0);
        assert!(off < peak);
        assert!(off > 0.0);
    }

    #[test]
    fn test_hg10_has_node_at_center() {
        // HG₁₀ has a node at x=0
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let at_center = hermite_gaussian(&beam, 1, 0, 0.0, 0.0, 0.0);
        assert!(at_center.abs() < TOL);
    }

    #[test]
    fn test_lg00_is_gaussian() {
        // LG₀₀ = fundamental Gaussian
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let peak = laguerre_gaussian_intensity(&beam, 0, 0, 0.0, 0.0);
        assert!((peak - 1.0).abs() < TOL);
    }

    #[test]
    fn test_lg_donut_mode() {
        // LG₀₁ (vortex/donut beam): zero on axis, ring pattern
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let on_axis = laguerre_gaussian_intensity(&beam, 0, 1, 0.0, 0.0);
        assert!(on_axis.abs() < TOL); // zero at center

        let off_axis = laguerre_gaussian_intensity(&beam, 0, 1, 70e-6, 0.0);
        assert!(off_axis > 0.0); // nonzero on ring
    }

    #[test]
    fn test_abcd_default_is_identity() {
        let m = AbcdMatrix::default();
        assert!((m.a - 1.0).abs() < TOL);
        assert!((m.d - 1.0).abs() < TOL);
        assert!(m.b.abs() < TOL);
        assert!(m.c.abs() < TOL);
    }
}
