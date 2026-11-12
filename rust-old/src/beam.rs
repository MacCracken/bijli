//! Gaussian beams, ABCD ray transfer matrices, resonator stability,
//! and higher-order Hermite-Gaussian / Laguerre-Gaussian modes.
//!
//! All SI units. Wavelength in meters, distances in meters.
//!
//! # Conventions
//!
//! - [`GaussianBeam::wavelength`] is the wavelength **in the medium** the beam
//!   currently propagates in (λ₀/n); the complex beam parameter is
//!   q = (z − z_waist) + i z_R with z_R = πw₀²/λ.
//! - ABCD matrices act on (height, actual ray angle), so a system that ends
//!   in a different medium has det = n_in/n_out, and the output wavelength is
//!   λ_out = λ_in · det.
//! - Complex mode fields ([`hermite_gaussian_field`]) use the physics
//!   **e^{−iωt}** convention (field ∝ e^{+ikz}), matching
//!   [`crate::polarization`]. Siegman's *Lasers* uses e^{+jωt}; conjugate to
//!   compare.

use serde::{Deserialize, Serialize};

use hisab::num::Complex;

use crate::error::{BijliError, Result};

// ── Gaussian beam parameter ───────────────────────────────────────

/// Gaussian beam parameters at a reference plane.
///
/// A fundamental Gaussian beam is fully described by its wavelength λ and
/// complex beam parameter q(z) = z - z_waist + i z_R, where z_R = πw₀²/λ.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GaussianBeam {
    /// Wavelength λ in the propagation medium (m), i.e. λ₀/n.
    pub wavelength: f64,
    /// Beam waist radius w₀ (m) — the 1/e² intensity radius at the focus.
    pub waist: f64,
    /// Rayleigh range z_R = πw₀²/λ (m).
    pub rayleigh_range: f64,
}

impl GaussianBeam {
    /// Create a Gaussian beam from wavelength (in the medium) and waist radius.
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] unless both are positive and finite
    /// and the resulting Rayleigh range is positive and finite.
    #[inline]
    pub fn new(wavelength: f64, waist: f64) -> Result<Self> {
        if !(wavelength.is_finite() && wavelength > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!("wavelength must be positive and finite, got {wavelength}"),
            });
        }
        if !(waist.is_finite() && waist > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!("beam waist must be positive and finite, got {waist}"),
            });
        }
        let rayleigh_range = std::f64::consts::PI * waist * waist / wavelength;
        if !(rayleigh_range.is_finite() && rayleigh_range > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "Rayleigh range πw₀²/λ = {rayleigh_range} is not representable \
                     (w₀ = {waist}, λ = {wavelength})"
                ),
            });
        }
        tracing::debug!(wavelength, waist, rayleigh_range, "Gaussian beam created");
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
        if focal_length.is_nan() {
            return Err(BijliError::InvalidParameter {
                reason: "focal length must not be NaN".into(),
            });
        }
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
        if radius.is_nan() {
            return Err(BijliError::InvalidParameter {
                reason: "mirror radius must not be NaN".into(),
            });
        }
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
        if n1.is_nan() || n2.is_nan() {
            return Err(BijliError::InvalidParameter {
                reason: "refractive indices must not be NaN".into(),
            });
        }
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
        if n1.is_nan() || n2.is_nan() || radius.is_nan() {
            return Err(BijliError::InvalidParameter {
                reason: "refractive indices and radius must not be NaN".into(),
            });
        }
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

    /// Output beam of this system for an input beam whose **waist lies at the
    /// input plane**.
    ///
    /// q_out = (A·i z_R + B)/(C·i z_R + D); the output beam has Rayleigh range
    /// z_R′ = Im(q_out), wavelength λ′ = λ·det(M) (det = n_in/n_out, so a
    /// system ending in a denser medium shortens λ) and waist
    /// w₀′ = √(z_R′ λ′/π).
    ///
    /// The waist *location* is not part of [`GaussianBeam`]; use
    /// [`transform_beam_with_waist_position`](Self::transform_beam_with_waist_position)
    /// to obtain it.
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] if det(M) ≤ 0 or the output Rayleigh
    /// range is not positive; [`BijliError::DivisionByZero`] if
    /// C·q + D = 0.
    #[inline]
    pub fn transform_beam(&self, beam: &GaussianBeam) -> Result<GaussianBeam> {
        self.transform_beam_with_waist_position(beam)
            .map(|(out, _)| out)
    }

    /// Like [`transform_beam`](Self::transform_beam), but also returns the
    /// signed distance from the **output plane to the new waist**:
    /// d = −Re(q_out). Positive d means the waist lies after (downstream of)
    /// the output plane (converging beam); negative means it lies before it
    /// (diverging beam, virtual waist).
    ///
    /// Example: a thin lens of focal length f placed at the input waist gives
    /// d = f z_R²/(z_R² + f²) and w₀′ = w₀ f/√(f² + z_R²) (Self, 1983).
    ///
    /// # Errors
    /// Same as [`transform_beam`](Self::transform_beam).
    pub fn transform_beam_with_waist_position(
        &self,
        beam: &GaussianBeam,
    ) -> Result<(GaussianBeam, f64)> {
        let det = self.determinant();
        if !(det.is_finite() && det > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "ABCD determinant must be positive and finite (det = n_in/n_out), got {det}"
                ),
            });
        }
        let q_in = beam.q(0.0); // at waist: q = i*z_R
        let q_out = self.propagate_beam(q_in)?;
        // Im(q_out) = z_R′ = π w₀′² / λ′
        let z_r_new = q_out[1];
        if !(z_r_new.is_finite() && z_r_new > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: "transformed beam has non-physical Rayleigh range".into(),
            });
        }
        // Same-medium systems have det = 1 analytically; treat rounding noise
        // as exactly 1 so the wavelength is preserved bit-for-bit.
        let wavelength = if (det - 1.0).abs() <= 1e-12 {
            beam.wavelength
        } else {
            beam.wavelength * det
        };
        let w0_new = (z_r_new * wavelength / std::f64::consts::PI).sqrt();
        let out = GaussianBeam::new(wavelength, w0_new)?;
        Ok((out, -q_out[0]))
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
/// - `length`: cavity length L (m, positive and finite)
/// - `r1`: radius of curvature of mirror 1 (m, positive = concave; ±∞ = flat)
/// - `r2`: radius of curvature of mirror 2 (m, positive = concave; ±∞ = flat)
#[inline]
pub fn resonator_stability(length: f64, r1: f64, r2: f64) -> Result<ResonatorStability> {
    if !(length.is_finite() && length > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("cavity length must be positive and finite, got {length}"),
        });
    }
    if r1.is_nan() || r2.is_nan() {
        return Err(BijliError::InvalidParameter {
            reason: "mirror radii must not be NaN".into(),
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
/// ```text
/// w₀² = (Lλ/π) · √[ g₁g₂(1 − g₁g₂) / (g₁ + g₂ − 2g₁g₂)² ]
/// ```
///
/// (Siegman, *Lasers*, ch. 19). The symmetric confocal cavity
/// (g₁ = g₂ = 0), where the expression is 0/0, returns its limit
/// w₀ = √(Lλ/2π). On the other stability edges the formula's limit is
/// returned (w₀ = 0 when exactly one gᵢ = 0 or g₁g₂ = 1 with gᵢ ≠ 1).
///
/// Returns the waist radius in meters.
///
/// # Errors
/// [`BijliError::InvalidParameter`] if the cavity is unstable or the
/// wavelength is not positive and finite; [`BijliError::DivisionByZero`]
/// for the planar cavity (g₁ = g₂ = 1), which has no confined mode.
#[inline]
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
    if !(wavelength.is_finite() && wavelength > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("wavelength must be positive and finite, got {wavelength}"),
        });
    }

    let g1 = stab.g1;
    let g2 = stab.g2;
    let gp = stab.g_product;

    // Inside the stability region s = g₁ + g₂ − 2g₁g₂ vanishes only at
    // (0, 0) (confocal) and (1, 1) (planar). Dividing by |s| rather than
    // comparing s² to a threshold keeps near-confocal cavities accurate.
    let s = g1 + g2 - 2.0 * gp;
    if s == 0.0 {
        if g1 == 0.0 && g2 == 0.0 {
            return Ok((length * wavelength / (2.0 * std::f64::consts::PI)).sqrt());
        }
        return Err(BijliError::DivisionByZero {
            context: "planar cavity (g₁ = g₂ = 1) has no confined mode".into(),
        });
    }

    // For a stable cavity (0 ≤ gp ≤ 1), gp*(1-gp) is always non-negative.
    let w0_sq = length * wavelength / std::f64::consts::PI * (gp * (1.0 - gp)).sqrt() / s.abs();
    Ok(w0_sq.sqrt())
}

/// Transverse mode spacing (frequency step per unit of m + n).
///
/// ```text
/// Δν_T = (Δν_FSR / π) · arccos(±√(g₁g₂))
/// ```
///
/// with the **+** sign when g₁, g₂ ≥ 0 and the **−** sign when g₁, g₂ < 0
/// (Siegman, *Lasers*, ch. 19), so that the TEM_mnq resonances lie at
/// ν = q·Δν_FSR + (m + n + 1)·Δν_T, with Δν_T ∈ [0, Δν_FSR].
///
/// `fsr` is the free spectral range c/(2L) (positive, finite).
///
/// # Errors
/// [`BijliError::InvalidParameter`] for an unstable cavity (g₁g₂ ∉ [0, 1]),
/// NaN g-parameters, or a non-positive/non-finite `fsr`.
#[inline]
pub fn transverse_mode_spacing(fsr: f64, g1: f64, g2: f64) -> Result<f64> {
    if !(fsr.is_finite() && fsr > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("free spectral range must be positive and finite, got {fsr}"),
        });
    }
    let gp = g1 * g2;
    if !(0.0..=1.0).contains(&gp) {
        return Err(BijliError::InvalidParameter {
            reason: format!("cavity is unstable (g₁g₂ = {gp:.4})"),
        });
    }
    let root = gp.sqrt();
    // Both g < 0 (e.g. near-concentric): lower sign. If either g is 0 the
    // root is 0 and the sign is irrelevant.
    let arg = if g1 < 0.0 || g2 < 0.0 { -root } else { root };
    Ok(fsr / std::f64::consts::PI * arg.acos())
}

// ── Hermite-Gaussian modes ────────────────────────────────────────

/// Hermite-Gaussian real amplitude envelope HG_mn at position (x, y, z).
///
/// ```text
/// u_mn(x, y, z) = (w₀/w) · H_m(√2 x/w) · H_n(√2 y/w) · exp(−(x² + y²)/w²)
/// ```
///
/// This is the field with all propagation phases (kz, k r²/2R, Gouy phase)
/// removed: |u_mn| = |E| and u_mn² is the intensity. The sign of the Hermite
/// lobes is kept. Use [`hermite_gaussian_field`] for the complex field.
///
/// Normalization: the TEM₀₀ prefactor is used for every mode, so
/// HG₀₀(0, 0, 0) = 1. Modes are **not** power-normalized: HG_mn carries
/// 2^{m+n} m! n! times the power of TEM₀₀; multiply by 1/√(2^{m+n} m! n!)
/// to compare modes at equal power.
///
/// # Parameters
/// - `beam`: the fundamental Gaussian beam parameters
/// - `m`, `n`: mode indices (m for x, n for y)
/// - `x`, `y`: transverse coordinates (m)
/// - `z`: distance from waist (m)
#[inline]
#[must_use]
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

    scale * hm * hn * envelope
}

/// Complex Hermite-Gaussian field HG_mn at (x, y, z), e^{−iωt} convention.
///
/// ```text
/// E_mn = u_mn(x, y, z) · exp(i [k z + k (x² + y²)/(2R(z)) − (m + n + 1) ζ(z)])
/// ```
///
/// with u_mn = [`hermite_gaussian`], k = 2π/λ (λ in the medium),
/// 1/R(z) = z/(z² + z_R²) (finite at the waist) and ζ(z) = arctan(z/z_R).
/// For Siegman's e^{+jωt} convention take the complex conjugate.
/// Normalization as in [`hermite_gaussian`].
#[inline]
#[must_use]
pub fn hermite_gaussian_field(
    beam: &GaussianBeam,
    m: u32,
    n: u32,
    x: f64,
    y: f64,
    z: f64,
) -> Complex {
    let amplitude = hermite_gaussian(beam, m, n, x, y, z);
    let k = 2.0 * std::f64::consts::PI / beam.wavelength;
    let z_r = beam.rayleigh_range;
    let inv_r = z / (z * z + z_r * z_r);
    let order = f64::from(m) + f64::from(n) + 1.0;
    let phase = k * z + 0.5 * k * (x * x + y * y) * inv_r - order * beam.gouy_phase(z);
    Complex::from_polar(amplitude, phase)
}

/// Laguerre-Gaussian beam intensity profile LG_pl at position (r, z).
///
/// ```text
/// I_pl(r, z) = (w₀/w)² · (2r²/w²)^|l| · [L_p^|l|(2r²/w²)]² · exp(−2r²/w²)
/// ```
///
/// This is |E|² (phases, including the e^{ilφ} vortex and Gouy phase
/// (2p + |l| + 1)ζ, drop out). Normalization: the TEM₀₀ prefactor is used for
/// every mode, so LG₀₀(0, 0) = 1. Modes are **not** power-normalized: LG_pl
/// carries (p + |l|)!/p! times the power of TEM₀₀; multiply by
/// p!/(p + |l|)! to compare modes at equal power.
///
/// # Parameters
/// - `beam`: the fundamental Gaussian beam parameters
/// - `p`: radial mode index (number of radial nodes)
/// - `l`: azimuthal mode index (orbital angular momentum)
/// - `r`: radial distance from axis (m)
/// - `z`: distance from waist (m)
#[inline]
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

    // ── Audit-driven edge case tests ──────────────────────────────

    #[test]
    fn test_curved_interface() {
        // Convex surface n1=1 → n2=1.5, R=0.1m
        let m = AbcdMatrix::curved_interface(1.0, 1.5, 0.1).unwrap();
        // det = n1/n2
        assert!((m.determinant() - 1.0 / 1.5).abs() < TOL);
        // C = (n1 - n2)/(n2 * R) = -0.5/0.15 = -10/3
        assert!((m.c - (-10.0 / 3.0)).abs() < TOL);
    }

    #[test]
    fn test_curved_interface_zero_n2() {
        assert!(AbcdMatrix::curved_interface(1.0, 0.0, 0.1).is_err());
    }

    #[test]
    fn test_curved_interface_zero_radius() {
        assert!(AbcdMatrix::curved_interface(1.0, 1.5, 0.0).is_err());
    }

    #[test]
    fn test_curved_mirror_zero_radius() {
        assert!(AbcdMatrix::curved_mirror(0.0).is_err());
    }

    #[test]
    fn test_flat_interface_zero_n2() {
        assert!(AbcdMatrix::flat_interface(1.0, 0.0).is_err());
    }

    #[test]
    fn test_peak_intensity_at_waist_vs_far() {
        // Peak intensity decreases away from waist
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let i_waist = beam.peak_intensity_factor(0.0);
        let i_far = beam.peak_intensity_factor(10.0 * beam.rayleigh_range);
        assert!(i_waist > i_far);
    }

    #[test]
    fn test_gouy_phase_limits() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        // At z → +∞, Gouy → π/2
        let far = beam.gouy_phase(1e10 * beam.rayleigh_range);
        assert!((far - std::f64::consts::FRAC_PI_2).abs() < 1e-4);
        // At z → -∞, Gouy → -π/2
        let neg_far = beam.gouy_phase(-1e10 * beam.rayleigh_range);
        assert!((neg_far - (-std::f64::consts::FRAC_PI_2)).abs() < 1e-4);
    }

    #[test]
    fn test_radius_of_curvature_at_rayleigh() {
        // R(z_R) = 2 z_R (minimum radius of curvature)
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let r = beam.radius_of_curvature(beam.rayleigh_range).unwrap();
        assert!((r - 2.0 * beam.rayleigh_range).abs() / r < TOL);
    }

    #[test]
    fn test_lens_imaging_condition() {
        // Object at 2f, image at 2f: compose free_space(2f) → lens(f) → free_space(2f)
        // Should give magnification = -1 (A = -1)
        let f = 0.1;
        let system = AbcdMatrix::free_space(2.0 * f)
            .compose(&AbcdMatrix::thin_lens(f).unwrap())
            .compose(&AbcdMatrix::free_space(2.0 * f));
        assert!((system.a - (-1.0)).abs() < TOL);
        assert!(system.b.abs() < TOL); // imaging condition: B = 0
    }

    #[test]
    fn test_resonator_beam_waist_invalid_wavelength() {
        assert!(resonator_beam_waist(0.1, 0.2, 0.2, 0.0).is_err());
        assert!(resonator_beam_waist(0.1, 0.2, 0.2, -1.0).is_err());
    }

    #[test]
    fn test_hg_higher_order_symmetry() {
        // HG₂₀ should be symmetric in x: HG(x) = HG(-x)
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let pos = hermite_gaussian(&beam, 2, 0, 50e-6, 0.0, 0.0);
        let neg = hermite_gaussian(&beam, 2, 0, -50e-6, 0.0, 0.0);
        assert!((pos - neg).abs() < 1e-15);
    }

    #[test]
    fn test_lg_radial_symmetry() {
        // LG intensity depends only on r, not angle
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let i1 = laguerre_gaussian_intensity(&beam, 1, 0, 50e-6, 0.0);
        let i2 = laguerre_gaussian_intensity(&beam, 1, 0, 50e-6, 0.01);
        // Different z but same r — should differ (z affects w)
        // But at same z, same r, intensity is same regardless of azimuthal angle
        // (which isn't a parameter — the function is already radially symmetric)
        assert!(i1 > 0.0);
        assert!(i2 > 0.0);
    }

    // ── Regression tests (closed-form references) ─────────────────

    #[test]
    fn test_resonator_beam_waist_symmetric_confocal() {
        // g₁ = g₂ = 0 is stable; Siegman: w₀ = √(Lλ/2π). Old code: Err.
        let (l, lambda) = (0.5, 1064e-9);
        let w = resonator_beam_waist(l, 0.5, 0.5, lambda).unwrap();
        let expected = (l * lambda / (2.0 * std::f64::consts::PI)).sqrt();
        assert!((w - expected).abs() / expected < 1e-12, "{w} vs {expected}");
    }

    #[test]
    fn test_resonator_beam_waist_symmetric_closed_form() {
        // Symmetric cavity: w₀² = (Lλ/2π)·√((1+g)/(1−g)); g = 0.5 → √3.
        let (l, lambda) = (0.1, 1064e-9);
        let w = resonator_beam_waist(l, 0.2, 0.2, lambda).unwrap();
        let expected = (l * lambda / (2.0 * std::f64::consts::PI) * 3.0_f64.sqrt()).sqrt();
        assert!((w - expected).abs() / expected < 1e-12);
        // Near-confocal (g ≈ ±1e-9) stays continuous with the confocal limit.
        let near =
            resonator_beam_waist(0.5, 0.5 * (1.0 + 1e-9), 0.5 * (1.0 + 1e-9), lambda).unwrap();
        let conf = resonator_beam_waist(0.5, 0.5, 0.5, lambda).unwrap();
        assert!((near - conf).abs() / conf < 1e-6);
    }

    #[test]
    fn test_resonator_beam_waist_planar_errs() {
        // Flat–flat cavity (R = ∞ → g = 1): no confined mode.
        assert!(resonator_beam_waist(0.1, f64::INFINITY, f64::INFINITY, 1064e-9).is_err());
        let stab = resonator_stability(0.1, f64::INFINITY, f64::INFINITY).unwrap();
        assert!((stab.g_product - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_transverse_mode_spacing_negative_g() {
        // Siegman: lower sign for g₁, g₂ < 0. g = −0.5: arccos(−0.5) = 2π/3
        // → Δν_T = 2·FSR/3. Old code returned FSR/3.
        let fsr = 1.5e9;
        let s = transverse_mode_spacing(fsr, -0.5, -0.5).unwrap();
        assert!((s - 2.0 * fsr / 3.0).abs() / fsr < 1e-12, "{s}");
        // Near-concentric (g → −1): Δν_T → FSR; near-planar (g → +1): → 0.
        assert!((transverse_mode_spacing(fsr, -1.0, -1.0).unwrap() - fsr).abs() / fsr < 1e-12);
        assert!(transverse_mode_spacing(fsr, 1.0, 1.0).unwrap().abs() < 1e-3);
        // Confocal: Δν_T = FSR/2 regardless of sign.
        assert!((transverse_mode_spacing(fsr, 0.0, 0.0).unwrap() - fsr / 2.0).abs() / fsr < 1e-12);
    }

    #[test]
    fn test_transverse_mode_spacing_invalid() {
        assert!(transverse_mode_spacing(f64::NAN, 0.5, 0.5).is_err());
        assert!(transverse_mode_spacing(0.0, 0.5, 0.5).is_err());
        assert!(transverse_mode_spacing(-1.0, 0.5, 0.5).is_err());
        assert!(transverse_mode_spacing(1e9, f64::NAN, 0.5).is_err());
    }

    #[test]
    fn test_hg_envelope_at_rayleigh_range() {
        // At z = z_R, w = √2 w₀. The old cos((m+n+1)ζ) factor made HG₁₀
        // cos(π/2) ≈ 0 everywhere and halved HG₀₀ on axis.
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let z_r = beam.rayleigh_range;
        // HG₁₀ at x = w₀: (1/√2)·H₁(1)·e^{−1/2} = √2 e^{−1/2}
        let v = hermite_gaussian(&beam, 1, 0, beam.waist, 0.0, z_r);
        let expected = std::f64::consts::SQRT_2 * (-0.5_f64).exp();
        assert!((v - expected).abs() < 1e-12, "{v} vs {expected}");
        // HG₀₀ on axis: |E| ∝ w₀/w = 1/√2 (intensity halves at z_R)
        let v = hermite_gaussian(&beam, 0, 0, 0.0, 0.0, z_r);
        assert!((v - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_gaussian_field_closed_form() {
        // λ = 1, w₀ = 1 → z_R = π, k = 2π. At z = z_R: w² = 2, R = 2z_R, ζ = π/4.
        let beam = GaussianBeam::new(1.0, 1.0).unwrap();
        let pi = std::f64::consts::PI;
        let z = pi;
        // HG₀₀ at x = 0.5: amplitude (1/√2)e^{−0.125},
        // phase k z + k x²/(2R) − ζ = 2π² + 2π·0.25/(4π) − π/4
        let e = hermite_gaussian_field(&beam, 0, 0, 0.5, 0.0, z);
        let amp = std::f64::consts::FRAC_1_SQRT_2 * (-0.125_f64).exp();
        let ph = 2.0 * pi * pi + 0.125 - pi / 4.0;
        assert!((e.re - amp * ph.cos()).abs() < 1e-12);
        assert!((e.im - amp * ph.sin()).abs() < 1e-12);
        // HG₂₀ on axis: H₂(0) = −2 → amplitude −√2, Gouy factor 3ζ.
        let e = hermite_gaussian_field(&beam, 2, 0, 0.0, 0.0, z);
        let ph = 2.0 * pi * pi - 3.0 * pi / 4.0;
        assert!((e.re + std::f64::consts::SQRT_2 * ph.cos()).abs() < 1e-12);
        assert!((e.im + std::f64::consts::SQRT_2 * ph.sin()).abs() < 1e-12);
        // At the waist the field is real and equals the envelope.
        let e = hermite_gaussian_field(&beam, 1, 1, 0.3, -0.2, 0.0);
        let u = hermite_gaussian(&beam, 1, 1, 0.3, -0.2, 0.0);
        assert!((e.re - u).abs() < 1e-15 && e.im.abs() < 1e-15);
        // |E| equals the envelope magnitude everywhere.
        let e = hermite_gaussian_field(&beam, 3, 1, 0.4, 0.7, 2.3);
        let u = hermite_gaussian(&beam, 3, 1, 0.4, 0.7, 2.3);
        assert!((e.abs() - u.abs()).abs() < 1e-12);
    }

    #[test]
    fn test_transform_beam_flat_interface() {
        // Air → glass (n = 1.5) at the waist: w₀ is continuous across the
        // interface, λ′ = λ/1.5, z_R′ = 1.5 z_R. Old code returned √1.5·w₀.
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let m = AbcdMatrix::flat_interface(1.0, 1.5).unwrap();
        let out = m.transform_beam(&beam).unwrap();
        assert!(
            (out.waist - beam.waist).abs() / beam.waist < 1e-12,
            "{}",
            out.waist
        );
        assert!((out.wavelength - 1064e-9 / 1.5).abs() / out.wavelength < 1e-12);
        assert!(
            (out.rayleigh_range - 1.5 * beam.rayleigh_range).abs() / out.rayleigh_range < 1e-12
        );
    }

    #[test]
    fn test_transform_beam_lens_waist_position() {
        // Self (1983), input waist at a thin lens:
        //   d = f z_R²/(z_R² + f²),  w₀′ = w₀ f/√(f² + z_R²)
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        let f = 0.1;
        let z_r = beam.rayleigh_range;
        let lens = AbcdMatrix::thin_lens(f).unwrap();
        let (out, d) = lens.transform_beam_with_waist_position(&beam).unwrap();
        let d_expected = f * z_r * z_r / (z_r * z_r + f * f);
        let w_expected = beam.waist * f / (f * f + z_r * z_r).sqrt();
        assert!(
            (d - d_expected).abs() / d_expected < 1e-12,
            "{d} vs {d_expected}"
        );
        assert!((out.waist - w_expected).abs() / w_expected < 1e-12);
        assert_eq!(out.wavelength, beam.wavelength);
        // Free space: the waist is a distance L behind the output plane.
        let (out, d) = AbcdMatrix::free_space(0.25)
            .transform_beam_with_waist_position(&beam)
            .unwrap();
        assert!((d + 0.25).abs() < 1e-15);
        assert!((out.waist - beam.waist).abs() / beam.waist < 1e-12);
    }

    #[test]
    fn test_transform_beam_invalid_determinant() {
        let beam = GaussianBeam::new(1064e-9, 100e-6).unwrap();
        assert!(
            AbcdMatrix::new(1.0, 0.0, 0.0, -1.0)
                .transform_beam(&beam)
                .is_err()
        );
        assert!(
            AbcdMatrix::new(1.0, 0.0, 0.0, f64::NAN)
                .transform_beam(&beam)
                .is_err()
        );
    }

    #[test]
    fn test_nan_inputs_rejected() {
        assert!(GaussianBeam::new(f64::NAN, 100e-6).is_err());
        assert!(GaussianBeam::new(1064e-9, f64::NAN).is_err());
        assert!(GaussianBeam::new(f64::INFINITY, 100e-6).is_err());
        assert!(GaussianBeam::new(1064e-9, 1e-200).is_err()); // z_R underflows to 0
        assert!(resonator_stability(f64::NAN, 1.0, 1.0).is_err());
        assert!(resonator_stability(0.1, f64::NAN, 1.0).is_err());
        assert!(resonator_stability(f64::INFINITY, 1.0, 1.0).is_err());
        assert!(resonator_beam_waist(0.1, 0.2, 0.2, f64::NAN).is_err());
        assert!(AbcdMatrix::thin_lens(f64::NAN).is_err());
        assert!(AbcdMatrix::curved_mirror(f64::NAN).is_err());
        assert!(AbcdMatrix::flat_interface(f64::NAN, 1.0).is_err());
        assert!(AbcdMatrix::curved_interface(1.0, 1.5, f64::NAN).is_err());
    }
}
