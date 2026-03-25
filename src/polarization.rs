//! Polarization formalism — Jones vectors/matrices, Stokes parameters, Mueller matrices.
//!
//! Jones formalism handles fully polarized coherent light. Stokes/Mueller
//! handles partially polarized and incoherent light.
//!
//! All angles in radians. Jones vectors use complex amplitudes (Re, Im pairs).

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

// ── Complex helper ────────────────────────────────────────────────

/// Minimal complex number for polarization math (avoids external dep).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    /// Create a complex number.
    #[inline]
    #[must_use]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Real number (zero imaginary part).
    #[inline]
    #[must_use]
    pub fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }

    /// Pure imaginary number.
    #[inline]
    #[must_use]
    pub fn imag(im: f64) -> Self {
        Self { re: 0.0, im }
    }

    /// Zero.
    #[inline]
    #[must_use]
    pub fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    /// Complex conjugate.
    #[inline]
    #[must_use]
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Squared magnitude |z|².
    #[inline]
    #[must_use]
    pub fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Magnitude |z|.
    #[inline]
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// e^(iθ) = cos θ + i sin θ.
    #[inline]
    #[must_use]
    pub fn from_polar(magnitude: f64, phase: f64) -> Self {
        Self {
            re: magnitude * phase.cos(),
            im: magnitude * phase.sin(),
        }
    }

    /// Phase angle arg(z).
    #[inline]
    #[must_use]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

impl std::ops::Add for Complex {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl std::ops::Mul<Complex> for f64 {
    type Output = Complex;
    #[inline]
    fn mul(self, rhs: Complex) -> Complex {
        Complex {
            re: self * rhs.re,
            im: self * rhs.im,
        }
    }
}

impl std::ops::Neg for Complex {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{}+{}i", self.re, self.im)
        } else {
            write!(f, "{}{}i", self.re, self.im)
        }
    }
}

// ── Jones vector ──────────────────────────────────────────────────

/// Jones vector — represents the polarization state of fully polarized light.
///
/// Components are the complex amplitudes of the x and y electric field:
/// **J** = [E_x, E_y]^T
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JonesVector {
    pub x: Complex,
    pub y: Complex,
}

impl JonesVector {
    /// Create a Jones vector from complex components.
    #[inline]
    #[must_use]
    pub fn new(x: Complex, y: Complex) -> Self {
        Self { x, y }
    }

    /// Horizontal linear polarization: [1, 0].
    #[inline]
    #[must_use]
    pub fn horizontal() -> Self {
        Self {
            x: Complex::real(1.0),
            y: Complex::zero(),
        }
    }

    /// Vertical linear polarization: [0, 1].
    #[inline]
    #[must_use]
    pub fn vertical() -> Self {
        Self {
            x: Complex::zero(),
            y: Complex::real(1.0),
        }
    }

    /// Diagonal (+45°) linear polarization: [1, 1]/√2.
    #[inline]
    #[must_use]
    pub fn diagonal() -> Self {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        Self {
            x: Complex::real(s),
            y: Complex::real(s),
        }
    }

    /// Anti-diagonal (−45°) linear polarization: [1, −1]/√2.
    #[inline]
    #[must_use]
    pub fn anti_diagonal() -> Self {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        Self {
            x: Complex::real(s),
            y: Complex::real(-s),
        }
    }

    /// Right circular polarization: [1, −i]/√2.
    #[inline]
    #[must_use]
    pub fn right_circular() -> Self {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        Self {
            x: Complex::real(s),
            y: Complex::new(0.0, -s),
        }
    }

    /// Left circular polarization: [1, +i]/√2.
    #[inline]
    #[must_use]
    pub fn left_circular() -> Self {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        Self {
            x: Complex::real(s),
            y: Complex::new(0.0, s),
        }
    }

    /// Linear polarization at angle θ from horizontal: [cos θ, sin θ].
    #[inline]
    #[must_use]
    pub fn linear(angle: f64) -> Self {
        Self {
            x: Complex::real(angle.cos()),
            y: Complex::real(angle.sin()),
        }
    }

    /// Intensity (squared magnitude): |E_x|² + |E_y|².
    #[inline]
    #[must_use]
    pub fn intensity(&self) -> f64 {
        self.x.norm_sq() + self.y.norm_sq()
    }

    /// Normalize to unit intensity.
    #[inline]
    pub fn normalized(&self) -> Result<Self> {
        let i = self.intensity();
        if i < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "cannot normalize zero Jones vector".into(),
            });
        }
        let s = 1.0 / i.sqrt();
        Ok(Self {
            x: self.x * s,
            y: self.y * s,
        })
    }

    /// Inner product ⟨J₁|J₂⟩ = x₁*x₂ + y₁*y₂ (conjugate-linear in first arg).
    #[inline]
    #[must_use]
    pub fn inner(&self, other: &Self) -> Complex {
        self.x.conj() * other.x + self.y.conj() * other.y
    }

    /// Degree of orthogonality: |⟨J₁|J₂⟩|² / (I₁ × I₂).
    ///
    /// Returns 0 for orthogonal states, 1 for parallel states.
    #[inline]
    pub fn overlap(&self, other: &Self) -> Result<f64> {
        let denom = self.intensity() * other.intensity();
        if denom < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "cannot compute overlap with zero-intensity state".into(),
            });
        }
        Ok(self.inner(other).norm_sq() / denom)
    }
}

// ── Jones matrix ──────────────────────────────────────────────────

/// Jones matrix — 2×2 complex matrix representing an optical element.
///
/// Stored as [[a, b], [c, d]] where the matrix acts as:
/// ```text
/// [E_x']   [a  b] [E_x]
/// [E_y'] = [c  d] [E_y]
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JonesMatrix {
    pub m: [[Complex; 2]; 2],
}

impl JonesMatrix {
    /// Create a Jones matrix from components.
    #[inline]
    #[must_use]
    pub fn new(a: Complex, b: Complex, c: Complex, d: Complex) -> Self {
        Self {
            m: [[a, b], [c, d]],
        }
    }

    /// Identity matrix (no effect on polarization).
    #[inline]
    #[must_use]
    pub fn identity() -> Self {
        Self::new(
            Complex::real(1.0),
            Complex::zero(),
            Complex::zero(),
            Complex::real(1.0),
        )
    }

    /// Apply this matrix to a Jones vector.
    #[inline]
    #[must_use]
    pub fn apply(&self, v: &JonesVector) -> JonesVector {
        JonesVector {
            x: self.m[0][0] * v.x + self.m[0][1] * v.y,
            y: self.m[1][0] * v.x + self.m[1][1] * v.y,
        }
    }

    /// Compose two Jones matrices: self × other (self applied after other).
    #[inline]
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        let a = self.m[0][0] * other.m[0][0] + self.m[0][1] * other.m[1][0];
        let b = self.m[0][0] * other.m[0][1] + self.m[0][1] * other.m[1][1];
        let c = self.m[1][0] * other.m[0][0] + self.m[1][1] * other.m[1][0];
        let d = self.m[1][0] * other.m[0][1] + self.m[1][1] * other.m[1][1];
        Self::new(a, b, c, d)
    }

    /// Horizontal linear polarizer.
    #[inline]
    #[must_use]
    pub fn horizontal_polarizer() -> Self {
        Self::new(
            Complex::real(1.0),
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
        )
    }

    /// Vertical linear polarizer.
    #[inline]
    #[must_use]
    pub fn vertical_polarizer() -> Self {
        Self::new(
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
            Complex::real(1.0),
        )
    }

    /// Linear polarizer at angle θ from horizontal.
    #[inline]
    #[must_use]
    pub fn linear_polarizer(angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(
            Complex::real(c * c),
            Complex::real(c * s),
            Complex::real(c * s),
            Complex::real(s * s),
        )
    }

    /// Half-wave plate (retarder with δ = π) with fast axis at angle θ.
    ///
    /// Flips polarization about the fast axis.
    #[inline]
    #[must_use]
    pub fn half_wave_plate(fast_axis_angle: f64) -> Self {
        Self::waveplate(fast_axis_angle, std::f64::consts::PI)
    }

    /// Quarter-wave plate (retarder with δ = π/2) with fast axis at angle θ.
    ///
    /// Converts linear to circular polarization (and vice versa) at 45°.
    #[inline]
    #[must_use]
    pub fn quarter_wave_plate(fast_axis_angle: f64) -> Self {
        Self::waveplate(fast_axis_angle, std::f64::consts::FRAC_PI_2)
    }

    /// General waveplate (linear retarder) with fast axis at angle θ
    /// and retardance δ (phase difference between fast and slow axes).
    #[inline]
    #[must_use]
    pub fn waveplate(fast_axis_angle: f64, retardance: f64) -> Self {
        let c = fast_axis_angle.cos();
        let s = fast_axis_angle.sin();
        let half_d = retardance / 2.0;
        let cos_d = half_d.cos();
        let sin_d = half_d.sin();

        // Rotation into fast-axis frame, apply phase, rotate back
        // M = R(-θ) × diag(e^{-iδ/2}, e^{iδ/2}) × R(θ)
        let a = Complex::new(cos_d, -sin_d * (c * c - s * s));
        let b = Complex::new(0.0, -sin_d * 2.0 * c * s);
        let d = Complex::new(cos_d, sin_d * (c * c - s * s));

        Self::new(a, b, b, d)
    }

    /// Faraday rotator — rotates polarization by angle θ.
    #[inline]
    #[must_use]
    pub fn rotator(angle: f64) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::new(
            Complex::real(c),
            Complex::real(-s),
            Complex::real(s),
            Complex::real(c),
        )
    }

    /// Determinant of the Jones matrix.
    #[inline]
    #[must_use]
    pub fn determinant(&self) -> Complex {
        self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]
    }
}

// ── Stokes vector ─────────────────────────────────────────────────

/// Stokes vector — describes polarization state including partially polarized light.
///
/// S = [S₀, S₁, S₂, S₃] where:
/// - S₀ = total intensity
/// - S₁ = horizontal vs vertical preference
/// - S₂ = diagonal vs anti-diagonal preference
/// - S₃ = right vs left circular preference
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StokesVector {
    pub s: [f64; 4],
}

impl StokesVector {
    /// Create a Stokes vector from components.
    #[inline]
    #[must_use]
    pub fn new(s0: f64, s1: f64, s2: f64, s3: f64) -> Self {
        Self {
            s: [s0, s1, s2, s3],
        }
    }

    /// Convert from a Jones vector (fully polarized).
    ///
    /// S₀ = |E_x|² + |E_y|²
    /// S₁ = |E_x|² − |E_y|²
    /// S₂ = 2 Re(E_x E_y*)
    /// S₃ = 2 Im(E_x E_y*)
    #[inline]
    #[must_use]
    pub fn from_jones(j: &JonesVector) -> Self {
        let ix = j.x.norm_sq();
        let iy = j.y.norm_sq();
        let cross = j.x * j.y.conj();
        Self {
            s: [ix + iy, ix - iy, 2.0 * cross.re, 2.0 * cross.im],
        }
    }

    /// Unpolarized light with given intensity.
    #[inline]
    #[must_use]
    pub fn unpolarized(intensity: f64) -> Self {
        Self {
            s: [intensity, 0.0, 0.0, 0.0],
        }
    }

    /// Horizontal linear polarization.
    #[inline]
    #[must_use]
    pub fn horizontal(intensity: f64) -> Self {
        Self {
            s: [intensity, intensity, 0.0, 0.0],
        }
    }

    /// Vertical linear polarization.
    #[inline]
    #[must_use]
    pub fn vertical(intensity: f64) -> Self {
        Self {
            s: [intensity, -intensity, 0.0, 0.0],
        }
    }

    /// Right circular polarization.
    #[inline]
    #[must_use]
    pub fn right_circular(intensity: f64) -> Self {
        Self {
            s: [intensity, 0.0, 0.0, intensity],
        }
    }

    /// Left circular polarization.
    #[inline]
    #[must_use]
    pub fn left_circular(intensity: f64) -> Self {
        Self {
            s: [intensity, 0.0, 0.0, -intensity],
        }
    }

    /// Total intensity S₀.
    #[inline]
    #[must_use]
    pub fn intensity(&self) -> f64 {
        self.s[0]
    }

    /// Degree of polarization: p = √(S₁² + S₂² + S₃²) / S₀.
    ///
    /// 0 = unpolarized, 1 = fully polarized.
    #[inline]
    pub fn degree_of_polarization(&self) -> Result<f64> {
        if self.s[0].abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "cannot compute DOP for zero-intensity light".into(),
            });
        }
        let pol = (self.s[1] * self.s[1] + self.s[2] * self.s[2] + self.s[3] * self.s[3]).sqrt();
        Ok(pol / self.s[0])
    }

    /// Ellipticity angle χ: tan(2χ) = S₃ / √(S₁² + S₂²).
    ///
    /// χ = 0 for linear, ±π/4 for circular.
    #[inline]
    #[must_use]
    pub fn ellipticity_angle(&self) -> f64 {
        let lin = (self.s[1] * self.s[1] + self.s[2] * self.s[2]).sqrt();
        0.5 * self.s[3].atan2(lin)
    }

    /// Orientation angle ψ of the polarization ellipse: tan(2ψ) = S₂/S₁.
    #[inline]
    #[must_use]
    pub fn orientation_angle(&self) -> f64 {
        0.5 * self.s[2].atan2(self.s[1])
    }

    /// Add two Stokes vectors (incoherent superposition).
    #[inline]
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            s: [
                self.s[0] + other.s[0],
                self.s[1] + other.s[1],
                self.s[2] + other.s[2],
                self.s[3] + other.s[3],
            ],
        }
    }
}

impl std::ops::Add for StokesVector {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            s: [
                self.s[0] + rhs.s[0],
                self.s[1] + rhs.s[1],
                self.s[2] + rhs.s[2],
                self.s[3] + rhs.s[3],
            ],
        }
    }
}

// ── Mueller matrix ────────────────────────────────────────────────

/// Mueller matrix — 4×4 real matrix acting on Stokes vectors.
///
/// Handles partially polarized and incoherent light, unlike Jones.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MuellerMatrix {
    pub m: [[f64; 4]; 4],
}

impl MuellerMatrix {
    /// Create from a 4×4 array.
    #[inline]
    #[must_use]
    pub fn new(m: [[f64; 4]; 4]) -> Self {
        Self { m }
    }

    /// Identity Mueller matrix.
    #[inline]
    #[must_use]
    pub fn identity() -> Self {
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Apply this Mueller matrix to a Stokes vector.
    #[inline]
    #[must_use]
    pub fn apply(&self, s: &StokesVector) -> StokesVector {
        let mut out = [0.0; 4];
        for (out_i, row) in out.iter_mut().zip(self.m.iter()) {
            *out_i = row.iter().zip(s.s.iter()).map(|(a, b)| a * b).sum();
        }
        StokesVector { s: out }
    }

    /// Compose two Mueller matrices: self × other.
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        let mut m = [[0.0; 4]; 4];
        for (m_row, self_row) in m.iter_mut().zip(self.m.iter()) {
            for (m_ij, j) in m_row.iter_mut().zip(0..) {
                *m_ij = self_row
                    .iter()
                    .zip(other.m.iter())
                    .map(|(a, other_row)| a * other_row[j])
                    .sum();
            }
        }
        Self { m }
    }

    /// Horizontal linear polarizer.
    #[inline]
    #[must_use]
    pub fn horizontal_polarizer() -> Self {
        Self {
            m: [
                [0.5, 0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ],
        }
    }

    /// Vertical linear polarizer.
    #[inline]
    #[must_use]
    pub fn vertical_polarizer() -> Self {
        Self {
            m: [
                [0.5, -0.5, 0.0, 0.0],
                [-0.5, 0.5, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ],
        }
    }

    /// Linear polarizer at angle θ.
    #[inline]
    #[must_use]
    pub fn linear_polarizer(angle: f64) -> Self {
        let c2 = (2.0 * angle).cos();
        let s2 = (2.0 * angle).sin();
        Self {
            m: [
                [0.5, 0.5 * c2, 0.5 * s2, 0.0],
                [0.5 * c2, 0.5 * c2 * c2, 0.5 * c2 * s2, 0.0],
                [0.5 * s2, 0.5 * c2 * s2, 0.5 * s2 * s2, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ],
        }
    }

    /// General waveplate (retarder) with fast axis at angle θ and retardance δ.
    #[inline]
    #[must_use]
    pub fn waveplate(fast_axis_angle: f64, retardance: f64) -> Self {
        let c2 = (2.0 * fast_axis_angle).cos();
        let s2 = (2.0 * fast_axis_angle).sin();
        let cd = retardance.cos();
        let sd = retardance.sin();
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, c2 * c2 + s2 * s2 * cd, c2 * s2 * (1.0 - cd), -s2 * sd],
                [0.0, c2 * s2 * (1.0 - cd), s2 * s2 + c2 * c2 * cd, c2 * sd],
                [0.0, s2 * sd, -c2 * sd, cd],
            ],
        }
    }

    /// Quarter-wave plate with fast axis at angle θ.
    #[inline]
    #[must_use]
    pub fn quarter_wave_plate(fast_axis_angle: f64) -> Self {
        Self::waveplate(fast_axis_angle, std::f64::consts::FRAC_PI_2)
    }

    /// Half-wave plate with fast axis at angle θ.
    #[inline]
    #[must_use]
    pub fn half_wave_plate(fast_axis_angle: f64) -> Self {
        Self::waveplate(fast_axis_angle, std::f64::consts::PI)
    }

    /// Rotation matrix (coordinate rotation by angle θ).
    #[inline]
    #[must_use]
    pub fn rotation(angle: f64) -> Self {
        let c2 = (2.0 * angle).cos();
        let s2 = (2.0 * angle).sin();
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, c2, s2, 0.0],
                [0.0, -s2, c2, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Ideal depolarizer: output is unpolarized regardless of input.
    #[inline]
    #[must_use]
    pub fn depolarizer() -> Self {
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ],
        }
    }

    /// Attenuator (neutral density filter): scales intensity by factor.
    #[inline]
    #[must_use]
    pub fn attenuator(transmission: f64) -> Self {
        Self {
            m: [
                [transmission, 0.0, 0.0, 0.0],
                [0.0, transmission, 0.0, 0.0],
                [0.0, 0.0, transmission, 0.0],
                [0.0, 0.0, 0.0, transmission],
            ],
        }
    }

    /// Right circular polarizer.
    #[inline]
    #[must_use]
    pub fn right_circular_polarizer() -> Self {
        Self {
            m: [
                [0.5, 0.0, 0.0, 0.5],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0, 0.5],
            ],
        }
    }

    /// Left circular polarizer.
    #[inline]
    #[must_use]
    pub fn left_circular_polarizer() -> Self {
        Self {
            m: [
                [0.5, 0.0, 0.0, -0.5],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [-0.5, 0.0, 0.0, 0.5],
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    // ── Complex tests ─────────────────────────────────────────────

    #[test]
    fn test_complex_mul() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a * b;
        // (1+2i)(3+4i) = 3+4i+6i+8i² = -5+10i
        assert!((c.re - (-5.0)).abs() < TOL);
        assert!((c.im - 10.0).abs() < TOL);
    }

    #[test]
    fn test_complex_conj() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.conj().im - (-4.0)).abs() < TOL);
    }

    #[test]
    fn test_complex_norm() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.norm() - 5.0).abs() < TOL);
    }

    #[test]
    fn test_complex_from_polar() {
        let z = Complex::from_polar(1.0, std::f64::consts::FRAC_PI_2);
        assert!(z.re.abs() < TOL);
        assert!((z.im - 1.0).abs() < TOL);
    }

    // ── Jones vector tests ────────────────────────────────────────

    #[test]
    fn test_standard_states_normalized() {
        let states = [
            JonesVector::horizontal(),
            JonesVector::vertical(),
            JonesVector::diagonal(),
            JonesVector::anti_diagonal(),
            JonesVector::right_circular(),
            JonesVector::left_circular(),
        ];
        for s in &states {
            assert!((s.intensity() - 1.0).abs() < TOL, "state not normalized");
        }
    }

    #[test]
    fn test_hv_orthogonal() {
        let h = JonesVector::horizontal();
        let v = JonesVector::vertical();
        assert!(h.inner(&v).norm() < TOL);
        assert!((h.overlap(&v).unwrap()).abs() < TOL);
    }

    #[test]
    fn test_da_orthogonal() {
        let d = JonesVector::diagonal();
        let a = JonesVector::anti_diagonal();
        assert!(d.inner(&a).norm() < TOL);
    }

    #[test]
    fn test_rl_orthogonal() {
        let r = JonesVector::right_circular();
        let l = JonesVector::left_circular();
        assert!(r.inner(&l).norm() < TOL);
    }

    #[test]
    fn test_linear_at_zero_is_horizontal() {
        let h = JonesVector::horizontal();
        let l = JonesVector::linear(0.0);
        assert!((l.x.re - h.x.re).abs() < TOL);
        assert!((l.y.re - h.y.re).abs() < TOL);
    }

    #[test]
    fn test_linear_at_90_is_vertical() {
        let v = JonesVector::vertical();
        let l = JonesVector::linear(std::f64::consts::FRAC_PI_2);
        assert!((l.x.re - v.x.re).abs() < TOL);
        assert!((l.y.re - v.y.re).abs() < TOL);
    }

    #[test]
    fn test_normalize_zero_fails() {
        let z = JonesVector::new(Complex::zero(), Complex::zero());
        assert!(z.normalized().is_err());
    }

    // ── Jones matrix tests ────────────────────────────────────────

    #[test]
    fn test_identity_preserves() {
        let h = JonesVector::horizontal();
        let result = JonesMatrix::identity().apply(&h);
        assert!((result.x.re - h.x.re).abs() < TOL);
        assert!((result.y.re - h.y.re).abs() < TOL);
    }

    #[test]
    fn test_h_polarizer_passes_h() {
        let h = JonesVector::horizontal();
        let result = JonesMatrix::horizontal_polarizer().apply(&h);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_h_polarizer_blocks_v() {
        let v = JonesVector::vertical();
        let result = JonesMatrix::horizontal_polarizer().apply(&v);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_v_polarizer_passes_v() {
        let v = JonesVector::vertical();
        let result = JonesMatrix::vertical_polarizer().apply(&v);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_v_polarizer_blocks_h() {
        let h = JonesVector::horizontal();
        let result = JonesMatrix::vertical_polarizer().apply(&h);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_linear_polarizer_malus_law() {
        // Malus's law: I = I₀ cos²θ
        let h = JonesVector::horizontal();
        let angle = std::f64::consts::PI / 6.0; // 30°
        let pol = JonesMatrix::linear_polarizer(angle);
        let result = pol.apply(&h);
        let expected = angle.cos() * angle.cos();
        assert!((result.intensity() - expected).abs() < TOL);
    }

    #[test]
    fn test_crossed_polarizers_block() {
        // H polarizer then V polarizer → zero
        let h = JonesVector::horizontal();
        let system =
            JonesMatrix::vertical_polarizer().compose(&JonesMatrix::horizontal_polarizer());
        let result = system.apply(&h);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_qwp_h_to_circular() {
        // QWP at 45° converts H to right circular
        let h = JonesVector::horizontal();
        let qwp = JonesMatrix::quarter_wave_plate(std::f64::consts::FRAC_PI_4);
        let result = qwp.apply(&h);
        // Should be circularly polarized: |E_x| = |E_y|
        assert!((result.x.norm() - result.y.norm()).abs() < TOL);
        // Phase difference should be π/2
        let phase_diff = result.y.arg() - result.x.arg();
        assert!((phase_diff.abs() - std::f64::consts::FRAC_PI_2).abs() < TOL);
    }

    #[test]
    fn test_hwp_flips_polarization() {
        // HWP at 22.5° converts H to D
        let h = JonesVector::horizontal();
        let hwp = JonesMatrix::half_wave_plate(std::f64::consts::PI / 8.0);
        let result = hwp.apply(&h);
        let d = JonesVector::diagonal();
        // Should be diagonal (up to global phase)
        let overlap = result.normalized().unwrap().overlap(&d).unwrap();
        assert!((overlap - 1.0).abs() < TOL);
    }

    #[test]
    fn test_hwp_preserves_intensity() {
        let d = JonesVector::diagonal();
        let hwp = JonesMatrix::half_wave_plate(0.0);
        let result = hwp.apply(&d);
        assert!((result.intensity() - d.intensity()).abs() < TOL);
    }

    #[test]
    fn test_qwp_preserves_intensity() {
        let d = JonesVector::diagonal();
        let qwp = JonesMatrix::quarter_wave_plate(0.0);
        let result = qwp.apply(&d);
        assert!((result.intensity() - d.intensity()).abs() < TOL);
    }

    #[test]
    fn test_rotator() {
        // 90° rotation: H → V
        let h = JonesVector::horizontal();
        let rot = JonesMatrix::rotator(std::f64::consts::FRAC_PI_2);
        let result = rot.apply(&h);
        assert!(result.x.norm() < TOL);
        assert!((result.y.re - 1.0).abs() < TOL);
    }

    #[test]
    fn test_compose_associative() {
        let a = JonesMatrix::linear_polarizer(0.3);
        let b = JonesMatrix::quarter_wave_plate(0.5);
        let c = JonesMatrix::rotator(0.7);

        let ab_c = a.compose(&b).compose(&c);
        let a_bc = a.compose(&b.compose(&c));

        // Matrix multiplication is not commutative but is associative
        // We check (AB)C = A(BC) — wait, compose is self × other,
        // so a.compose(&b) = a × b. Then (a×b).compose(&c) = (a×b) × c
        // and a.compose(&(b.compose(&c))) = a × (b×c). Same thing.
        let v = JonesVector::diagonal();
        let r1 = ab_c.apply(&v);
        let r2 = a_bc.apply(&v);
        assert!((r1.x.re - r2.x.re).abs() < TOL);
        assert!((r1.x.im - r2.x.im).abs() < TOL);
        assert!((r1.y.re - r2.y.re).abs() < TOL);
        assert!((r1.y.im - r2.y.im).abs() < TOL);
    }

    #[test]
    fn test_determinant_polarizer() {
        // Polarizer is singular (det = 0)
        let det = JonesMatrix::horizontal_polarizer().determinant();
        assert!(det.norm() < TOL);
    }

    #[test]
    fn test_determinant_waveplate() {
        // Waveplate is unitary (|det| = 1)
        let det = JonesMatrix::quarter_wave_plate(0.3).determinant();
        assert!((det.norm() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_determinant_rotator() {
        // Rotator is unitary (|det| = 1)
        let det = JonesMatrix::rotator(0.7).determinant();
        assert!((det.norm() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_waveplate_zero_retardance_is_identity() {
        let wp = JonesMatrix::waveplate(0.5, 0.0);
        let h = JonesVector::horizontal();
        let result = wp.apply(&h);
        assert!((result.x.re - 1.0).abs() < TOL);
        assert!(result.y.norm() < TOL);
    }

    // ── Stokes vector tests ───────────────────────────────────────

    #[test]
    fn test_stokes_from_jones_horizontal() {
        let s = StokesVector::from_jones(&JonesVector::horizontal());
        assert!((s.s[0] - 1.0).abs() < TOL); // S₀ = 1
        assert!((s.s[1] - 1.0).abs() < TOL); // S₁ = 1 (horizontal)
        assert!(s.s[2].abs() < TOL);
        assert!(s.s[3].abs() < TOL);
    }

    #[test]
    fn test_stokes_from_jones_vertical() {
        let s = StokesVector::from_jones(&JonesVector::vertical());
        assert!((s.s[0] - 1.0).abs() < TOL);
        assert!((s.s[1] - (-1.0)).abs() < TOL); // S₁ = -1 (vertical)
        assert!(s.s[2].abs() < TOL);
        assert!(s.s[3].abs() < TOL);
    }

    #[test]
    fn test_stokes_from_jones_diagonal() {
        let s = StokesVector::from_jones(&JonesVector::diagonal());
        assert!((s.s[0] - 1.0).abs() < TOL);
        assert!(s.s[1].abs() < TOL);
        assert!((s.s[2] - 1.0).abs() < TOL); // S₂ = 1 (diagonal)
        assert!(s.s[3].abs() < TOL);
    }

    #[test]
    fn test_stokes_from_jones_right_circular() {
        let s = StokesVector::from_jones(&JonesVector::right_circular());
        assert!((s.s[0] - 1.0).abs() < TOL);
        assert!(s.s[1].abs() < TOL);
        assert!(s.s[2].abs() < TOL);
        assert!((s.s[3] - 1.0).abs() < TOL); // S₃ = +1 (RCP: E_y lags E_x by π/2)
    }

    #[test]
    fn test_stokes_fully_polarized_dop() {
        let s = StokesVector::from_jones(&JonesVector::horizontal());
        assert!((s.degree_of_polarization().unwrap() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_stokes_unpolarized_dop() {
        let s = StokesVector::unpolarized(1.0);
        assert!(s.degree_of_polarization().unwrap().abs() < TOL);
    }

    #[test]
    fn test_stokes_partially_polarized() {
        // Mix H and unpolarized: S = [2, 1, 0, 0] → DOP = 0.5
        let s = StokesVector::horizontal(1.0).add(&StokesVector::unpolarized(1.0));
        assert!((s.degree_of_polarization().unwrap() - 0.5).abs() < TOL);
    }

    #[test]
    fn test_stokes_orientation_horizontal() {
        let s = StokesVector::horizontal(1.0);
        assert!(s.orientation_angle().abs() < TOL); // ψ = 0 for horizontal
    }

    #[test]
    fn test_stokes_ellipticity_linear() {
        let s = StokesVector::horizontal(1.0);
        assert!(s.ellipticity_angle().abs() < TOL); // χ = 0 for linear
    }

    #[test]
    fn test_stokes_ellipticity_circular() {
        let s = StokesVector::right_circular(1.0);
        assert!((s.ellipticity_angle() - std::f64::consts::FRAC_PI_4).abs() < TOL);
    }

    #[test]
    fn test_stokes_zero_intensity_dop_fails() {
        let s = StokesVector::new(0.0, 0.0, 0.0, 0.0);
        assert!(s.degree_of_polarization().is_err());
    }

    // ── Mueller matrix tests ──────────────────────────────────────

    #[test]
    fn test_mueller_identity() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::identity().apply(&s);
        for i in 0..4 {
            assert!((result.s[i] - s.s[i]).abs() < TOL);
        }
    }

    #[test]
    fn test_mueller_h_polarizer_passes_h() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::horizontal_polarizer().apply(&s);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_mueller_h_polarizer_blocks_v() {
        let s = StokesVector::vertical(1.0);
        let result = MuellerMatrix::horizontal_polarizer().apply(&s);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_mueller_v_polarizer_blocks_h() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::vertical_polarizer().apply(&s);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_mueller_polarizer_halves_unpolarized() {
        let s = StokesVector::unpolarized(1.0);
        let result = MuellerMatrix::horizontal_polarizer().apply(&s);
        assert!((result.intensity() - 0.5).abs() < TOL);
    }

    #[test]
    fn test_mueller_waveplate_preserves_intensity() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::quarter_wave_plate(0.3).apply(&s);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_mueller_depolarizer() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::depolarizer().apply(&s);
        assert!((result.intensity() - 1.0).abs() < TOL);
        assert!(result.s[1].abs() < TOL);
        assert!(result.s[2].abs() < TOL);
        assert!(result.s[3].abs() < TOL);
    }

    #[test]
    fn test_mueller_attenuator() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::attenuator(0.5).apply(&s);
        assert!((result.intensity() - 0.5).abs() < TOL);
    }

    #[test]
    fn test_mueller_rcp_polarizer() {
        let s = StokesVector::right_circular(1.0);
        let result = MuellerMatrix::right_circular_polarizer().apply(&s);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_mueller_rcp_blocks_lcp() {
        let s = StokesVector::left_circular(1.0);
        let result = MuellerMatrix::right_circular_polarizer().apply(&s);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_mueller_crossed_polarizers() {
        let system =
            MuellerMatrix::vertical_polarizer().compose(&MuellerMatrix::horizontal_polarizer());
        let s = StokesVector::horizontal(1.0);
        let result = system.apply(&s);
        assert!(result.intensity() < TOL);
    }

    #[test]
    fn test_mueller_rotation_preserves_intensity() {
        let s = StokesVector::horizontal(1.0);
        let result = MuellerMatrix::rotation(0.5).apply(&s);
        assert!((result.intensity() - 1.0).abs() < TOL);
    }

    #[test]
    fn test_jones_mueller_agreement_polarizer() {
        // Jones and Mueller should agree for fully polarized light
        let j = JonesVector::diagonal();
        let jones_result = JonesMatrix::horizontal_polarizer().apply(&j);
        let jones_stokes = StokesVector::from_jones(&jones_result);

        let stokes_in = StokesVector::from_jones(&j);
        let mueller_result = MuellerMatrix::horizontal_polarizer().apply(&stokes_in);

        for i in 0..4 {
            assert!(
                (jones_stokes.s[i] - mueller_result.s[i]).abs() < TOL,
                "mismatch at S[{i}]"
            );
        }
    }

    // ── Audit-driven edge case tests ──────────────────────────────

    #[test]
    fn test_stokes_add_operator() {
        let a = StokesVector::horizontal(1.0);
        let b = StokesVector::vertical(1.0);
        let c = a + b;
        // H + V = unpolarized with intensity 2
        assert!((c.intensity() - 2.0).abs() < TOL);
        assert!(c.s[1].abs() < TOL); // S₁ cancels
    }

    #[test]
    fn test_complex_display() {
        let pos = Complex::new(1.0, 2.0);
        assert_eq!(format!("{pos}"), "1+2i");
        let neg = Complex::new(1.0, -2.0);
        assert_eq!(format!("{neg}"), "1-2i");
    }

    #[test]
    fn test_complex_sub() {
        let a = Complex::new(5.0, 3.0);
        let b = Complex::new(2.0, 1.0);
        let c = a - b;
        assert!((c.re - 3.0).abs() < TOL);
        assert!((c.im - 2.0).abs() < TOL);
    }

    #[test]
    fn test_complex_neg() {
        let z = Complex::new(3.0, -4.0);
        let neg = -z;
        assert!((neg.re - (-3.0)).abs() < TOL);
        assert!((neg.im - 4.0).abs() < TOL);
    }

    #[test]
    fn test_complex_arg() {
        let z = Complex::new(0.0, 1.0);
        assert!((z.arg() - std::f64::consts::FRAC_PI_2).abs() < TOL);
    }

    #[test]
    fn test_jones_mueller_agreement_waveplate() {
        // Jones and Mueller QWP should agree for fully polarized light
        let j = JonesVector::horizontal();
        let angle = 0.3;

        let jones_out = JonesMatrix::quarter_wave_plate(angle).apply(&j);
        let jones_stokes = StokesVector::from_jones(&jones_out);

        let stokes_in = StokesVector::from_jones(&j);
        let mueller_out = MuellerMatrix::quarter_wave_plate(angle).apply(&stokes_in);

        for i in 0..4 {
            assert!(
                (jones_stokes.s[i] - mueller_out.s[i]).abs() < 1e-8,
                "QWP mismatch at S[{i}]: jones={} mueller={}",
                jones_stokes.s[i],
                mueller_out.s[i]
            );
        }
    }

    #[test]
    fn test_mueller_compose_identity() {
        let m = MuellerMatrix::horizontal_polarizer();
        let composed = MuellerMatrix::identity().compose(&m);
        let s = StokesVector::horizontal(1.0);
        let r1 = m.apply(&s);
        let r2 = composed.apply(&s);
        for i in 0..4 {
            assert!((r1.s[i] - r2.s[i]).abs() < TOL);
        }
    }

    #[test]
    fn test_mie_result_serde() {
        // MieResult should roundtrip through serde
        let result = crate::scattering::mie(1.0, Complex::real(1.5)).unwrap();
        let json = serde_json::to_string(&result).unwrap();
        let back: crate::scattering::MieResult = serde_json::from_str(&json).unwrap();
        assert!((back.q_ext - result.q_ext).abs() < TOL);
    }
}
