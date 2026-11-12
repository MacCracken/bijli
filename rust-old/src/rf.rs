//! RF & Microwave — transmission lines, S-parameters, Smith chart,
//! matching networks, Touchstone I/O, antenna arrays.
//!
//! All SI units unless noted: impedances in ohms, admittances in siemens,
//! frequencies in Hz, lengths in metres, angles in radians. Phasors use the
//! e^{+jωt} convention (inductive reactance is +jωL, a lossless line has
//! γ = jβ).
//!
//! # Contents
//!
//! - **Transmission lines** — lossless/lossy characteristic impedance,
//!   propagation constant, phase velocity, input impedance of a terminated
//!   lossless line ([`input_impedance_lossless_complex`]).
//! - **Reflection** — Γ (real and complex), VSWR, return / insertion /
//!   mismatch loss.
//! - **Smith chart** — impedance normalisation and Γ ↔ z mapping.
//! - **Matching** — quarter-wave transformer, single shunt short-circuited
//!   stub ([`single_stub_match`], [`single_stub_match_complex`]).
//! - **S-parameters** — N-port [`SMatrix`] container; for 2-ports, cascading
//!   ([`SMatrix::cascade`], computed internally via transfer (T) parameters)
//!   and S → Z conversion ([`s_to_z`]). T-parameters are not exposed and no
//!   Z → S, Y or ABCD conversions are provided. Conversions assume all ports
//!   share one real reference impedance.
//! - **Touchstone v1.1** — read/write of `.sNp` S-parameter data
//!   ([`parse_touchstone`], [`try_write_touchstone`]).
//! - **Antenna arrays** — linear, planar and circular array factors.

use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};
use crate::polarization::{Complex, ComplexExt, complex_from_polar, complex_real, complex_zero};

/// Maximum port count accepted by the Touchstone reader/writer.
///
/// Bounds the per-record buffer (1 + 2n² values) so that a caller-supplied
/// port count derived from untrusted input cannot request a huge allocation.
const MAX_TOUCHSTONE_PORTS: usize = 256;

// ── Internal helpers ──────────────────────────────────────────────

#[inline]
fn is_finite_complex(c: Complex) -> bool {
    c.re.is_finite() && c.im.is_finite()
}

/// Relative equality used for reference-impedance consistency checks.
#[inline]
fn z0_equal(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}

#[inline]
fn require_positive_finite(name: &str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!("{name} must be finite and positive, got {value}"),
        })
    }
}

/// Complex division `a / b` using Smith's algorithm, which avoids the
/// premature overflow/underflow of the textbook `a·b̄ / |b|²` form.
///
/// The caller must ensure `b ≠ 0` (the result is NaN otherwise).
#[inline]
fn cdiv(a: Complex, b: Complex) -> Complex {
    if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let d = b.re + b.im * r;
        Complex::new((a.re + a.im * r) / d, (a.im - a.re * r) / d)
    } else {
        let r = b.re / b.im;
        let d = b.re * r + b.im;
        Complex::new((a.re * r + a.im) / d, (a.im * r - a.re) / d)
    }
}

/// Reciprocal 1/z: a single division when |z|² is a normal float, falling
/// back to Smith's algorithm ([`cdiv`]) when |z|² would over/underflow.
#[inline]
fn crecip(z: Complex) -> Complex {
    let d = z.re * z.re + z.im * z.im;
    if d.is_normal() {
        let inv = 1.0 / d;
        Complex::new(z.re * inv, -z.im * inv)
    } else {
        cdiv(complex_real(1.0), z)
    }
}

/// Complex infinity (the point at infinity of the Smith chart / Riemann
/// sphere), represented as `+∞ + 0j`.
const COMPLEX_INFINITY: Complex = Complex::new(f64::INFINITY, 0.0);

// ── Transmission lines ────────────────────────────────────────────

/// Characteristic impedance of a lossless transmission line: Z₀ = √(L'/C').
///
/// `inductance_per_m` (H/m) and `capacitance_per_m` (F/m) must be finite and
/// positive.
#[inline]
pub fn characteristic_impedance_lossless(
    inductance_per_m: f64,
    capacitance_per_m: f64,
) -> Result<f64> {
    if !(inductance_per_m.is_finite()
        && inductance_per_m > 0.0
        && capacitance_per_m.is_finite()
        && capacitance_per_m > 0.0)
    {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "L' and C' must be finite and positive, got L'={inductance_per_m}, C'={capacitance_per_m}"
            ),
        });
    }
    Ok((inductance_per_m / capacitance_per_m).sqrt())
}

/// Complex characteristic impedance of a lossy line:
/// Z₀ = √((R' + jωL') / (G' + jωC')), principal branch (Re Z₀ ≥ 0).
///
/// Per-unit-length parameters: `r` (Ω/m), `l` (H/m), `g` (S/m), `c` (F/m);
/// `omega` is the angular frequency (rad/s).
///
/// # Singular cases
///
/// The expression requires G' + jωC' ≠ 0. At DC (`omega == 0`) with
/// `g == 0` and `c > 0` the ω → 0⁺ limit is returned instead:
/// √(L'/C') when `r == 0` (lossless line), or complex infinity
/// `+∞ − j∞` (|Z₀| → ∞ at phase −45°) when `r > 0`. Any other case with
/// G' + jωC' = 0 is undefined and yields NaN components.
#[inline]
#[must_use]
pub fn characteristic_impedance_lossy(r: f64, l: f64, g: f64, c: f64, omega: f64) -> Complex {
    let num = Complex::new(r, omega * l);
    let den = Complex::new(g, omega * c);
    if den.re == 0.0 && den.im == 0.0 {
        if omega == 0.0 && g == 0.0 && c > 0.0 {
            if r == 0.0 {
                return complex_real(l / c).sqrt();
            }
            if r > 0.0 {
                return Complex::new(f64::INFINITY, f64::NEG_INFINITY);
            }
        }
        return Complex::new(f64::NAN, f64::NAN);
    }
    let ratio = cdiv(num, den);
    let mag = ratio.norm().sqrt();
    let phase = ratio.arg() / 2.0;
    complex_from_polar(mag, phase)
}

/// Complex propagation constant: γ = α + jβ = √((R'+jωL')(G'+jωC')),
/// principal branch (α ≥ 0, β ≥ 0 for a passive line at ω > 0).
///
/// α in Np/m, β in rad/m.
#[inline]
#[must_use]
pub fn propagation_constant(r: f64, l: f64, g: f64, c: f64, omega: f64) -> Complex {
    let a = Complex::new(r, omega * l);
    let b = Complex::new(g, omega * c);
    let prod = a * b;
    let mag = prod.norm().sqrt();
    let phase = prod.arg() / 2.0;
    complex_from_polar(mag, phase)
}

/// Phase velocity on a lossless line: v_p = 1/√(L'C') (m/s).
#[inline]
pub fn phase_velocity(inductance_per_m: f64, capacitance_per_m: f64) -> Result<f64> {
    if !(inductance_per_m.is_finite()
        && inductance_per_m > 0.0
        && capacitance_per_m.is_finite()
        && capacitance_per_m > 0.0)
    {
        return Err(BijliError::InvalidParameter {
            reason: "L' and C' must be finite and positive".into(),
        });
    }
    Ok(1.0 / (inductance_per_m * capacitance_per_m).sqrt())
}

/// Input impedance of a lossless transmission line of characteristic
/// impedance `z0` terminated in a complex load `z_load`:
///
/// Z_in = Z₀ (Z_L + jZ₀ tan βl) / (Z₀ + jZ_L tan βl)
///
/// evaluated in the equivalent form
/// Z₀ (Z_L cos βl + jZ₀ sin βl) / (Z₀ cos βl + jZ_L sin βl), which stays
/// finite at βl = π/2. `beta` is the phase constant (rad/m) and `length`
/// the line length (m); only the product βl matters, so `beta = 2π` with
/// `length` in wavelengths also works.
///
/// # Errors
///
/// - [`BijliError::InvalidParameter`] if `z0` is not finite and positive, or
///   `z_load`, `beta·length` are not finite.
/// - [`BijliError::DivisionByZero`] if the denominator vanishes (relative to
///   |Z₀| + |Z_L|), i.e. Z_in is unbounded — e.g. a short-circuited
///   quarter-wave line.
#[inline]
pub fn input_impedance_lossless_complex(
    z0: f64,
    z_load: Complex,
    beta: f64,
    length: f64,
) -> Result<Complex> {
    require_positive_finite("Z₀", z0)?;
    if !is_finite_complex(z_load) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Z_L must be finite, got {z_load}"),
        });
    }
    let bl = beta * length;
    if !bl.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("electrical length βl must be finite, got β={beta}, l={length}"),
        });
    }
    let (s, c) = bl.sin_cos();
    // num = Z_L cos βl + j Z₀ sin βl
    let num = Complex::new(z_load.re * c, z_load.im * c + z0 * s);
    // den = Z₀ cos βl + j Z_L sin βl
    let den = Complex::new(z0 * c - z_load.im * s, z_load.re * s);
    let scale = z0 + z_load.norm();
    let den_mag = den.norm();
    if den_mag.is_nan() || den_mag <= 1e-15 * scale {
        return Err(BijliError::DivisionByZero {
            context: "transmission line input impedance is unbounded (Z₀cos βl + jZ_L sin βl = 0)"
                .into(),
        });
    }
    Ok(cdiv(num, den) * z0)
}

/// Input impedance of a lossless line with a *real* load, as a real number.
///
/// Z_in is complex in general (Z₀ (Z_L + jZ₀ tan βl)/(Z₀ + jZ_L tan βl)), so
/// this only succeeds when the imaginary part is negligible
/// (|Im Z_in| ≤ 10⁻⁹ |Z_in|) — e.g. a matched load or a line that is a
/// multiple of λ/4 long — and returns `Re Z_in`.
///
/// # Errors
///
/// Everything [`input_impedance_lossless_complex`] returns, plus
/// [`BijliError::InvalidParameter`] when Z_in has a non-negligible
/// imaginary part.
#[deprecated(
    since = "1.2.0",
    note = "Z_in is complex in general; use input_impedance_lossless_complex"
)]
#[inline]
pub fn input_impedance_lossless(z0: f64, z_load: f64, beta: f64, length: f64) -> Result<f64> {
    let z = input_impedance_lossless_complex(z0, complex_real(z_load), beta, length)?;
    if z.im.abs() <= 1e-9 * z.norm() {
        Ok(z.re)
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!(
                "input impedance is complex ({} {:+} j Ω); use input_impedance_lossless_complex",
                z.re, z.im
            ),
        })
    }
}

// ── Reflection coefficient & VSWR ─────────────────────────────────

/// Voltage reflection coefficient: Γ = (Z_L − Z₀) / (Z_L + Z₀).
#[inline]
pub fn reflection_coefficient(z_load: f64, z0: f64) -> Result<f64> {
    if !(z_load.is_finite() && z0.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Z_L and Z₀ must be finite, got Z_L={z_load}, Z₀={z0}"),
        });
    }
    let den = z_load + z0;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Z_L + Z₀ = 0".into(),
        });
    }
    Ok((z_load - z0) / den)
}

/// Complex reflection coefficient for complex impedances:
/// Γ = (Z_L − Z₀) / (Z_L + Z₀).
#[inline]
pub fn reflection_coefficient_complex(z_load: Complex, z0: Complex) -> Result<Complex> {
    if !(is_finite_complex(z_load) && is_finite_complex(z0)) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Z_L and Z₀ must be finite, got Z_L={z_load}, Z₀={z0}"),
        });
    }
    let den = z_load + z0;
    if den.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "Z_L + Z₀ = 0".into(),
        });
    }
    Ok(cdiv(z_load - z0, den))
}

/// VSWR from reflection coefficient magnitude: VSWR = (1 + |Γ|) / (1 − |Γ|).
///
/// Requires 0 ≤ |Γ| < 1 (|Γ| = 1 gives infinite VSWR).
#[inline]
pub fn vswr(gamma_mag: f64) -> Result<f64> {
    if gamma_mag >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be < 1, got {gamma_mag}"),
        });
    }
    if gamma_mag.is_nan() || gamma_mag < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be ≥ 0, got {gamma_mag}"),
        });
    }
    Ok((1.0 + gamma_mag) / (1.0 - gamma_mag))
}

/// Return loss in dB: RL = −20 log₁₀|Γ|.
///
/// Requires 0 < |Γ| ≤ 1 (a perfect match has infinite return loss; |Γ| > 1
/// is not a passive load and is rejected, consistent with [`vswr`]).
#[inline]
pub fn return_loss(gamma_mag: f64) -> Result<f64> {
    if !(gamma_mag > 0.0 && gamma_mag <= 1.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be in (0, 1] for return loss, got {gamma_mag}"),
        });
    }
    Ok(-20.0 * gamma_mag.log10())
}

/// Insertion loss in dB due to mismatch alone: IL = −10 log₁₀(1 − |Γ|²).
///
/// Requires 0 ≤ |Γ| < 1 (|Γ| = 1 is total reflection, infinite loss).
#[inline]
pub fn insertion_loss(gamma_mag: f64) -> Result<f64> {
    if gamma_mag >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: "total reflection — infinite insertion loss".into(),
        });
    }
    if gamma_mag.is_nan() || gamma_mag < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be in [0, 1), got {gamma_mag}"),
        });
    }
    Ok(-10.0 * (1.0 - gamma_mag * gamma_mag).log10())
}

/// Mismatch loss in dB: ML = −10 log₁₀(1 − |Γ|²).
///
/// Same as insertion loss for a lossless network.
#[inline]
pub fn mismatch_loss(gamma_mag: f64) -> Result<f64> {
    insertion_loss(gamma_mag)
}

// ── Smith chart ───────────────────────────────────────────────────

/// Normalize impedance: z = Z / Z₀.
#[inline]
pub fn normalize_impedance(z: Complex, z0: f64) -> Result<Complex> {
    if !z0.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("reference impedance Z₀ must be finite, got {z0}"),
        });
    }
    if z0.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "reference impedance Z₀ cannot be zero".into(),
        });
    }
    Ok(z * (1.0 / z0))
}

/// Denormalize impedance: Z = z × Z₀.
#[inline]
#[must_use]
pub fn denormalize_impedance(z_norm: Complex, z0: f64) -> Complex {
    z_norm * z0
}

/// Reflection coefficient from normalized impedance: Γ = (z − 1) / (z + 1).
///
/// Poles and infinities: z = −1 maps to complex infinity, returned as
/// `+∞ + 0j`; an infinite z (open circuit, either component ±∞) maps to
/// Γ = 1. NaN input yields NaN components. This makes `gamma_from_z` and
/// [`z_from_gamma`] mutually inverse at their poles.
#[inline]
#[must_use]
pub fn gamma_from_z(z_norm: Complex) -> Complex {
    if z_norm.re.is_nan() || z_norm.im.is_nan() {
        return Complex::new(f64::NAN, f64::NAN);
    }
    if z_norm.re.is_infinite() || z_norm.im.is_infinite() {
        return complex_real(1.0);
    }
    let one = complex_real(1.0);
    let den = z_norm + one;
    if den.re == 0.0 && den.im == 0.0 {
        return COMPLEX_INFINITY;
    }
    cdiv(z_norm - one, den)
}

/// Normalized impedance from reflection coefficient: z = (1 + Γ) / (1 − Γ).
///
/// Poles and infinities: Γ = 1 (open circuit) maps to complex infinity,
/// returned as `+∞ + 0j`; an infinite Γ maps to z = −1. NaN input yields
/// NaN components.
#[inline]
#[must_use]
pub fn z_from_gamma(gamma: Complex) -> Complex {
    if gamma.re.is_nan() || gamma.im.is_nan() {
        return Complex::new(f64::NAN, f64::NAN);
    }
    if gamma.re.is_infinite() || gamma.im.is_infinite() {
        return complex_real(-1.0);
    }
    let one = complex_real(1.0);
    let den = one - gamma;
    if den.re == 0.0 && den.im == 0.0 {
        return COMPLEX_INFINITY;
    }
    cdiv(one + gamma, den)
}

// ── Matching networks ─────────────────────────────────────────────

/// Quarter-wave transformer impedance: Z_t = √(Z₀ × Z_L).
///
/// Matches a real load `z_load` to a line `z0`; both must be finite and
/// positive.
#[inline]
pub fn quarter_wave_transformer(z0: f64, z_load: f64) -> Result<f64> {
    if !(z0.is_finite() && z0 > 0.0 && z_load.is_finite() && z_load > 0.0) {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "Z₀ and Z_L must be finite and positive for a quarter-wave transformer, got Z₀={z0}, Z_L={z_load}"
            ),
        });
    }
    Ok((z0 * z_load).sqrt())
}

/// Single shunt short-circuited stub match for a real load (Pozar §5.2).
///
/// Returns `(d, ℓ)` in wavelengths: `d` is the distance from the load
/// (toward the generator) at which a short-circuited shunt stub of the same
/// Z₀ and length `ℓ` is connected so that the line is matched. Of the two
/// solutions, the one with the smaller `d` is returned; see
/// [`single_stub_match_complex`] for both.
///
/// A matched load (`z_load == z0`) returns `(0, 0.25)`: a λ/4
/// short-circuited stub presents zero admittance, i.e. it is equivalent to
/// no stub.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `z_load` or `z0` is not finite and
/// positive.
#[inline]
pub fn single_stub_match(z_load: f64, z0: f64) -> Result<(f64, f64)> {
    require_positive_finite("Z_L", z_load)?;
    let [first, _] = single_stub_match_complex(complex_real(z_load), z0)?;
    Ok(first)
}

/// Single shunt short-circuited stub match for a complex load (Pozar §5.2,
/// eqs. 5.9–5.11).
///
/// With Z_L = R_L + jX_L and t = tan βd:
///
/// - R_L ≠ Z₀: t = [X_L ± √(R_L((Z₀ − R_L)² + X_L²)/Z₀)] / (R_L − Z₀)
/// - R_L = Z₀: t = −X_L/(2Z₀), plus the t → ∞ solution d = λ/4
///
/// d/λ = atan(t)/2π wrapped into [0, ½). At that point the line admittance
/// is Y = Y₀ + jB, and a short-circuited stub with input admittance
/// −jY₀ cot βℓ = −jB cancels B: ℓ/λ = atan2(Y₀, B)/2π ∈ (0, ½).
///
/// Returns both `(d, ℓ)` solutions in wavelengths, sorted by increasing `d`.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `z0` is not finite and positive, or
/// `z_load` is not finite or has R_L ≤ 0 (a lossless stub cannot match a
/// purely reactive or active load).
pub fn single_stub_match_complex(z_load: Complex, z0: f64) -> Result<[(f64, f64); 2]> {
    require_positive_finite("Z₀", z0)?;
    if !is_finite_complex(z_load) {
        return Err(BijliError::InvalidParameter {
            reason: format!("Z_L must be finite, got {z_load}"),
        });
    }
    let (rl, xl) = (z_load.re, z_load.im);
    if rl <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("load resistance must be positive for a stub match, got {rl}"),
        });
    }

    // d/λ from t = tan βd, wrapped into [0, 0.5); `+ 0.0` normalises −0.0.
    let d_from_t = |t: f64| {
        let d = t.atan() / TAU;
        if d < 0.0 { d + 0.5 } else { d + 0.0 }
    };
    let (d1, d2) = if (rl - z0).abs() <= 1e-12 * z0 {
        (d_from_t(-xl / (2.0 * z0)), 0.25)
    } else {
        let disc = (rl * ((z0 - rl) * (z0 - rl) + xl * xl) / z0).sqrt();
        (
            d_from_t((xl + disc) / (rl - z0)),
            d_from_t((xl - disc) / (rl - z0)),
        )
    };

    let y0 = 1.0 / z0;
    let stub_length = |d: f64| -> Result<f64> {
        // β = 2π per wavelength, length in wavelengths.
        let z_in = input_impedance_lossless_complex(z0, z_load, TAU, d)?;
        let y = crecip(z_in);
        Ok(y0.atan2(y.im) / TAU)
    };
    let mut solutions = [(d1, stub_length(d1)?), (d2, stub_length(d2)?)];
    if solutions[1].0 < solutions[0].0 {
        solutions.swap(0, 1);
    }
    tracing::debug!(
        z_load = %z_load,
        z0,
        d1 = solutions[0].0,
        l1 = solutions[0].1,
        d2 = solutions[1].0,
        l2 = solutions[1].1,
        "single shunt stub match"
    );
    Ok(solutions)
}

// ── S-parameters ──────────────────────────────────────────────────

/// S-parameter matrix for an N-port network.
///
/// Invariant expected by the methods (the fields are public, so it is not
/// enforced by construction): `s` is `ports × ports` and `z0` has `ports`
/// entries. Methods return [`BijliError::InvalidParameter`] instead of
/// panicking when it does not hold.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SMatrix {
    /// Number of ports.
    pub ports: usize,
    /// S-parameter data: `s[i][j]` = S_{i+1,j+1}.
    pub s: Vec<Vec<Complex>>,
    /// Reference impedance per port (ohms).
    pub z0: Vec<f64>,
}

impl SMatrix {
    /// Create an N-port S-matrix initialized to zero.
    #[must_use]
    pub fn new(ports: usize, z0: f64) -> Self {
        Self {
            ports,
            s: vec![vec![complex_zero(); ports]; ports],
            z0: vec![z0; ports],
        }
    }

    /// Create a 2-port S-matrix from individual parameters.
    #[must_use]
    pub fn two_port(s11: Complex, s12: Complex, s21: Complex, s22: Complex, z0: f64) -> Self {
        Self {
            ports: 2,
            s: vec![vec![s11, s12], vec![s21, s22]],
            z0: vec![z0; 2],
        }
    }

    /// Check that `s` is `ports × ports` and `z0` has `ports` entries.
    fn check_shape(&self) -> Result<()> {
        let n = self.ports;
        if self.s.len() != n || self.s.iter().any(|row| row.len() != n) || self.z0.len() != n {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "inconsistent SMatrix: ports = {n}, but s has {} rows and z0 has {} entries \
                     (expected {n}×{n} and {n})",
                    self.s.len(),
                    self.z0.len()
                ),
            });
        }
        Ok(())
    }

    /// The common reference impedance of all ports (finite, positive, equal).
    fn uniform_z0(&self) -> Result<f64> {
        let first = *self
            .z0
            .first()
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: "SMatrix has no reference impedances".into(),
            })?;
        require_positive_finite("port reference impedance", first)?;
        if self.z0.iter().any(|&z| !z0_equal(z, first)) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "all port reference impedances must be equal, got {:?}",
                    self.z0
                ),
            });
        }
        Ok(first)
    }

    /// `[S11, S12, S21, S22]` of a well-formed 2-port.
    fn two_port_params(&self) -> Result<[Complex; 4]> {
        if self.ports != 2 {
            return Err(BijliError::InvalidParameter {
                reason: format!("expected a 2-port network, got {} ports", self.ports),
            });
        }
        self.check_shape()?;
        let at = |i: usize, j: usize| self.s.get(i).and_then(|row| row.get(j)).copied();
        match (at(0, 0), at(0, 1), at(1, 0), at(1, 1)) {
            (Some(s11), Some(s12), Some(s21), Some(s22)) => Ok([s11, s12, s21, s22]),
            _ => Err(BijliError::InvalidParameter {
                reason: "inconsistent 2-port SMatrix data".into(),
            }),
        }
    }

    /// Get S_{ij} (1-indexed ports).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Result<Complex> {
        if i == 0 || j == 0 || i > self.ports || j > self.ports {
            return Err(BijliError::InvalidParameter {
                reason: format!("port indices must be 1..={}, got ({i}, {j})", self.ports),
            });
        }
        self.s
            .get(i - 1)
            .and_then(|row| row.get(j - 1))
            .copied()
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("inconsistent SMatrix: no entry for S{i}{j}"),
            })
    }

    /// Set S_{ij} (1-indexed ports).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: Complex) -> Result<()> {
        if i == 0 || j == 0 || i > self.ports || j > self.ports {
            return Err(BijliError::InvalidParameter {
                reason: format!("port indices must be 1..={}, got ({i}, {j})", self.ports),
            });
        }
        let slot = self
            .s
            .get_mut(i - 1)
            .and_then(|row| row.get_mut(j - 1))
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("inconsistent SMatrix: no entry for S{i}{j}"),
            })?;
        *slot = value;
        Ok(())
    }

    /// Cascade two 2-port networks (self followed by `other`, port 2 of
    /// `self` connected to port 1 of `other`).
    ///
    /// Computed by multiplying transfer matrices T with
    /// `[b₁, a₁]ᵀ = T [a₂, b₂]ᵀ`, T = (1/S21) [[−ΔS, S11], [−S22, 1]].
    ///
    /// # Errors
    ///
    /// - [`BijliError::InvalidParameter`] if either network is not a
    ///   well-formed 2-port, or the four port reference impedances are not
    ///   all equal (finite, positive).
    /// - [`BijliError::DivisionByZero`] if S21 of either network is zero (no
    ///   T-parameter representation) or the cascade is singular.
    pub fn cascade(&self, other: &Self) -> Result<Self> {
        if self.ports != 2 || other.ports != 2 {
            return Err(BijliError::InvalidParameter {
                reason: "cascade requires two 2-port networks".into(),
            });
        }
        let a = self.two_port_params()?;
        let b = other.two_port_params()?;
        let z0 = self.uniform_z0()?;
        let z0_other = other.uniform_z0()?;
        if !z0_equal(z0, z0_other) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "cascade requires equal reference impedances, got {z0} Ω and {z0_other} Ω"
                ),
            });
        }
        tracing::trace!(z0, "cascading 2-port S-matrices");

        let t1 = s_to_t(a)?;
        let t2 = s_to_t(b)?;

        // T = T1 × T2 (2×2 complex matrix multiply)
        let t = [
            [
                t1[0][0] * t2[0][0] + t1[0][1] * t2[1][0],
                t1[0][0] * t2[0][1] + t1[0][1] * t2[1][1],
            ],
            [
                t1[1][0] * t2[0][0] + t1[1][1] * t2[1][0],
                t1[1][0] * t2[0][1] + t1[1][1] * t2[1][1],
            ],
        ];

        t_to_s(&t, z0)
    }
}

/// Convert 2-port S-parameters `[S11, S12, S21, S22]` to T (transfer)
/// parameters.
fn s_to_t([s11, s12, s21, s22]: [Complex; 4]) -> Result<[[Complex; 2]; 2]> {
    if s21.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "S21 = 0, cannot convert to T-parameters".into(),
        });
    }
    let inv_s21 = crecip(s21);
    let det = s11 * s22 - s12 * s21;
    Ok([
        [-(det * inv_s21), s11 * inv_s21],
        [-(s22 * inv_s21), inv_s21],
    ])
}

/// Convert T (transfer) parameters to 2-port S-parameters.
fn t_to_s(t: &[[Complex; 2]; 2], z0: f64) -> Result<SMatrix> {
    let t22 = t[1][1];
    if t22.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "T22 = 0, cannot convert to S-parameters".into(),
        });
    }
    let inv_t22 = crecip(t22);
    let s11 = t[0][1] * inv_t22;
    let s12 = (t[0][0] * t[1][1] - t[0][1] * t[1][0]) * inv_t22;
    let s21 = inv_t22;
    let s22 = -(t[1][0] * inv_t22);
    Ok(SMatrix::two_port(s11, s12, s21, s22, z0))
}

/// Convert 2-port S-parameters to Z-parameters (ohms).
///
/// Z = Z₀ (I + S)(I − S)⁻¹ for a common real reference impedance Z₀
/// (Pozar Table 4.2). Returns `[[Z11, Z12], [Z21, Z22]]`.
///
/// # Errors
///
/// - [`BijliError::InvalidParameter`] if `s` is not a well-formed 2-port or
///   its port reference impedances are not all equal (finite, positive).
/// - [`BijliError::DivisionByZero`] if (I − S) is singular (e.g. an ideal
///   thru, which has no Z-parameters).
pub fn s_to_z(s: &SMatrix) -> Result<[[Complex; 2]; 2]> {
    if s.ports != 2 {
        return Err(BijliError::InvalidParameter {
            reason: "S→Z conversion requires 2-port network".into(),
        });
    }
    let [s11, s12, s21, s22] = s.two_port_params()?;
    let z0_value = s.uniform_z0()?;
    tracing::trace!(z0 = z0_value, "converting 2-port S-parameters to Z");
    let one = complex_real(1.0);
    let z0 = complex_real(z0_value);

    let det = (one - s11) * (one - s22) - s12 * s21;
    if det.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "singular S→Z conversion".into(),
        });
    }
    let inv_det = crecip(det);

    Ok([
        [
            z0 * ((one + s11) * (one - s22) + s12 * s21) * inv_det,
            z0 * complex_real(2.0) * s12 * inv_det,
        ],
        [
            z0 * complex_real(2.0) * s21 * inv_det,
            z0 * ((one - s11) * (one + s22) + s12 * s21) * inv_det,
        ],
    ])
}

// ── Touchstone I/O ────────────────────────────────────────────────

/// Touchstone data format options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TouchstoneFormat {
    /// Magnitude and angle (degrees).
    MagnitudeAngle,
    /// Real and imaginary parts.
    RealImag,
    /// dB magnitude (20 log₁₀|S|) and angle (degrees).
    DbAngle,
}

/// Parsed Touchstone data for an N-port network.
///
/// Consistency expected by [`try_write_touchstone`]: `frequencies` and
/// `s_matrices` have equal length, and every matrix is a well-formed
/// `ports`-port [`SMatrix`] whose reference impedances all equal `z0`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TouchstoneData {
    /// Number of ports.
    pub ports: usize,
    /// Reference impedance (ohms) — the option line's `R` value.
    pub z0: f64,
    /// Frequency points (Hz).
    pub frequencies: Vec<f64>,
    /// S-parameter matrices at each frequency.
    pub s_matrices: Vec<SMatrix>,
}

fn invalid_port_count(ports: usize) -> BijliError {
    BijliError::InvalidParameter {
        reason: format!("Touchstone port count must be in 1..={MAX_TOUCHSTONE_PORTS}, got {ports}"),
    }
}

/// Number of values in one network-data record: 1 + 2n².
fn touchstone_record_len(ports: usize) -> Result<usize> {
    if ports == 0 || ports > MAX_TOUCHSTONE_PORTS {
        return Err(invalid_port_count(ports));
    }
    ports
        .checked_mul(ports)
        .and_then(|nn| nn.checked_mul(2))
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| invalid_port_count(ports))
}

/// Truncate untrusted text for inclusion in an error message.
fn snippet(s: &str) -> String {
    const MAX_CHARS: usize = 32;
    if s.chars().count() > MAX_CHARS {
        let mut out: String = s.chars().take(MAX_CHARS).collect();
        out.push('…');
        out
    } else {
        s.to_string()
    }
}

/// Apply a Touchstone option line (text after `#`).
fn parse_option_line(
    opts: &str,
    line_no: usize,
    freq_mult: &mut f64,
    format: &mut TouchstoneFormat,
    z0: &mut f64,
) -> Result<()> {
    let mut tokens = opts.split_whitespace();
    while let Some(tok) = tokens.next() {
        match tok.to_ascii_uppercase().as_str() {
            "HZ" => *freq_mult = 1.0,
            "KHZ" => *freq_mult = 1e3,
            "MHZ" => *freq_mult = 1e6,
            "GHZ" => *freq_mult = 1e9,
            "S" => {}
            "Y" | "Z" | "H" | "G" => {
                return Err(BijliError::InvalidParameter {
                    reason: format!(
                        "Touchstone line {line_no}: only S-parameter data is supported, got parameter type {tok}"
                    ),
                });
            }
            "MA" => *format = TouchstoneFormat::MagnitudeAngle,
            "RI" => *format = TouchstoneFormat::RealImag,
            "DB" => *format = TouchstoneFormat::DbAngle,
            "R" => {
                let value = tokens.next().ok_or_else(|| BijliError::InvalidParameter {
                    reason: format!(
                        "Touchstone line {line_no}: option R must be followed by the reference impedance"
                    ),
                })?;
                let r = value
                    .parse::<f64>()
                    .ok()
                    .filter(|r| r.is_finite() && *r > 0.0)
                    .ok_or_else(|| BijliError::InvalidParameter {
                        reason: format!(
                            "Touchstone line {line_no}: reference impedance must be a finite positive number, got {:?}",
                            snippet(value)
                        ),
                    })?;
                *z0 = r;
            }
            _ => {
                return Err(BijliError::InvalidParameter {
                    reason: format!(
                        "Touchstone line {line_no}: unrecognized option {:?}",
                        snippet(tok)
                    ),
                });
            }
        }
    }
    Ok(())
}

/// Build one S-matrix from a complete record `[f, v₁, v₂, …]`.
fn touchstone_record(
    record: &[f64],
    ports: usize,
    format: TouchstoneFormat,
    z0: f64,
    line_no: usize,
) -> Result<SMatrix> {
    let values = record.get(1..).unwrap_or(&[]);
    let mut sm = SMatrix::new(ports, z0);
    let (pairs, _) = values.as_chunks::<2>();
    for (k, &[v1, v2]) in pairs.iter().enumerate() {
        // v1.1: 2-port data is ordered N11 N21 N12 N22 (column-major);
        // every other port count is row-major.
        let (i, j) = if ports == 2 {
            (k % 2, k / 2)
        } else {
            (k / ports, k % ports)
        };
        let c = parse_complex_pair(v1, v2, format);
        if !is_finite_complex(c) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "Touchstone record ending at line {line_no}: S{}{} = ({v1}, {v2}) is not representable",
                    i + 1,
                    j + 1
                ),
            });
        }
        let slot =
            sm.s.get_mut(i)
                .and_then(|row| row.get_mut(j))
                .ok_or_else(|| BijliError::InvalidParameter {
                    reason: "Touchstone record has more values than the port count allows".into(),
                })?;
        *slot = c;
    }
    Ok(sm)
}

/// Parse Touchstone v1.1 (`.sNp`) S-parameter content for an N-port network.
///
/// - The option line `# <unit> <param> <format> R <z0>` is case-insensitive;
///   defaults are `GHz S MA R 50`. Only the first option line is used (later
///   ones are ignored, per the spec) and it must precede the data. Units:
///   Hz/kHz/MHz/GHz; formats: RI, MA, DB (angles in degrees). Parameter types
///   other than S are rejected.
/// - `!` starts a comment anywhere on a line.
/// - Values are read as a token stream, so records may span several lines
///   (as v1.1 requires for 3+ ports: one matrix row per line, at most four
///   pairs per line). Each record is 1 + 2n² numbers.
/// - Data order: 2-port records are `f N11 N21 N12 N22`; all other port
///   counts are row-major (`N11 N12 … N1n N21 …`).
/// - For 2-ports, a record whose frequency does not increase marks the start
///   of the noise-parameter block, which is skipped.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `ports` is not in 1..=256, an option
/// is unrecognized or malformed (including a missing/non-positive `R`), the
/// parameter type is not S, a Touchstone 2.0 keyword is present, any data
/// token is not a finite number, a frequency is negative, the data ends with
/// a partial record, or no records are found.
pub fn parse_touchstone(content: &str, ports: usize) -> Result<TouchstoneData> {
    let record_len = touchstone_record_len(ports)?;
    tracing::debug!(ports, bytes = content.len(), "parsing Touchstone data");

    let mut frequencies: Vec<f64> = Vec::new();
    let mut s_matrices = Vec::new();
    let mut z0 = 50.0;
    let mut format = TouchstoneFormat::MagnitudeAngle;
    let mut freq_mult = 1e9; // default GHz
    let mut seen_option = false;
    let mut in_noise_block = false;
    let mut record: Vec<f64> = Vec::with_capacity(record_len);

    for (idx, raw) in content.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw.split('!').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if let Some(opts) = line.strip_prefix('#') {
            if seen_option {
                tracing::debug!(line = line_no, "ignoring additional Touchstone option line");
                continue;
            }
            if !frequencies.is_empty() || !record.is_empty() || in_noise_block {
                return Err(BijliError::InvalidParameter {
                    reason: format!(
                        "Touchstone line {line_no}: option line must precede the network data"
                    ),
                });
            }
            seen_option = true;
            parse_option_line(opts, line_no, &mut freq_mult, &mut format, &mut z0)?;
            continue;
        }
        if line.starts_with('[') {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "Touchstone line {line_no}: Touchstone 2.0 keywords are not supported (v1.1 only)"
                ),
            });
        }

        for tok in line.split_whitespace() {
            let value = tok
                .parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| BijliError::InvalidParameter {
                    reason: format!(
                        "Touchstone line {line_no}: expected a finite number, got {:?}",
                        snippet(tok)
                    ),
                })?;
            if in_noise_block {
                continue;
            }
            if record.is_empty() {
                // First value of a record: its frequency.
                let freq = value * freq_mult;
                if !(freq.is_finite() && freq >= 0.0) {
                    return Err(BijliError::InvalidParameter {
                        reason: format!(
                            "Touchstone line {line_no}: frequency must be finite and non-negative, got {value}"
                        ),
                    });
                }
                if ports == 2 && frequencies.last().is_some_and(|&last| freq <= last) {
                    tracing::debug!(line = line_no, "skipping 2-port noise parameter block");
                    in_noise_block = true;
                    continue;
                }
            }
            record.push(value);
            if record.len() == record_len {
                let sm = touchstone_record(&record, ports, format, z0, line_no)?;
                let freq = record.first().copied().unwrap_or(0.0) * freq_mult;
                frequencies.push(freq);
                s_matrices.push(sm);
                record.clear();
            }
        }
    }

    if !record.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "Touchstone data ends with an incomplete record: {} of {record_len} values",
                record.len()
            ),
        });
    }
    if frequencies.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "no valid data found in Touchstone content".into(),
        });
    }

    tracing::debug!(
        ports,
        points = frequencies.len(),
        z0,
        ?format,
        "parsed Touchstone data"
    );
    Ok(TouchstoneData {
        ports,
        z0,
        frequencies,
        s_matrices,
    })
}

/// Check the consistency that [`try_write_touchstone`] requires.
fn validate_touchstone(data: &TouchstoneData) -> Result<()> {
    touchstone_record_len(data.ports)?;
    require_positive_finite("Touchstone reference impedance", data.z0)?;
    if data.frequencies.len() != data.s_matrices.len() {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "Touchstone data has {} frequencies but {} S-matrices",
                data.frequencies.len(),
                data.s_matrices.len()
            ),
        });
    }
    let mut previous: Option<f64> = None;
    for (k, (&freq, sm)) in data.frequencies.iter().zip(&data.s_matrices).enumerate() {
        if !(freq.is_finite() && freq >= 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!("frequency #{k} must be finite and non-negative, got {freq}"),
            });
        }
        if previous.is_some_and(|p| freq <= p) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "Touchstone frequencies must be strictly increasing (#{k}: {freq})"
                ),
            });
        }
        previous = Some(freq);
        if sm.ports != data.ports {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "S-matrix #{k} has {} ports, expected {}",
                    sm.ports, data.ports
                ),
            });
        }
        sm.check_shape()?;
        if sm.z0.iter().any(|&z| !z0_equal(z, data.z0)) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "S-matrix #{k} reference impedances {:?} differ from the file's R {}",
                    sm.z0, data.z0
                ),
            });
        }
        if !sm.s.iter().flatten().all(|&c| is_finite_complex(c)) {
            return Err(BijliError::InvalidParameter {
                reason: format!("S-matrix #{k} contains non-finite values"),
            });
        }
    }
    Ok(())
}

fn push_touchstone_pair(out: &mut String, c: Complex, format: TouchstoneFormat) -> Result<()> {
    let (v1, v2) = format_complex_pair(c, format);
    if !(v1.is_finite() && v2.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("S-parameter {c} cannot be written in {format:?} format"),
        });
    }
    out.push_str(&format!("  {v1:e}  {v2:e}"));
    Ok(())
}

/// Format S-parameter data as Touchstone v1.1 text.
///
/// Writes `# GHz S <fmt> R <z0>` followed by one record per frequency.
/// Numbers use the shortest round-trip scientific notation, so the output
/// re-parses (with [`parse_touchstone`]) to the same values up to the
/// GHz ↔ Hz scaling. 1- and 2-port records occupy one line (2-port order
/// `f N11 N21 N12 N22`); for 3+ ports each matrix row starts a new line
/// with at most four pairs per line. In DB format |S| is clamped to
/// [`f64::MIN_POSITIVE`] so a zero magnitude writes as a finite dB value.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if the data is inconsistent (see
/// [`TouchstoneData`]), `ports` is not in 1..=256, `z0` is not finite and
/// positive, frequencies are not finite, non-negative and strictly
/// increasing, or any S value is non-finite.
pub fn try_write_touchstone(data: &TouchstoneData, format: TouchstoneFormat) -> Result<String> {
    validate_touchstone(data)?;
    tracing::debug!(
        ports = data.ports,
        points = data.frequencies.len(),
        ?format,
        "writing Touchstone data"
    );

    let mut out = String::new();
    out.push_str("! Touchstone file generated by bijli\n");

    let fmt_str = match format {
        TouchstoneFormat::MagnitudeAngle => "MA",
        TouchstoneFormat::RealImag => "RI",
        TouchstoneFormat::DbAngle => "DB",
    };
    out.push_str(&format!("# GHz S {fmt_str} R {}\n", data.z0));

    for (freq, sm) in data.frequencies.iter().zip(&data.s_matrices) {
        out.push_str(&format!("{:e}", freq / 1e9));
        if data.ports == 2 {
            let [s11, s12, s21, s22] = sm.two_port_params()?;
            for c in [s11, s21, s12, s22] {
                push_touchstone_pair(&mut out, c, format)?;
            }
        } else {
            let mut first_line = true;
            for row in &sm.s {
                for chunk in row.chunks(4) {
                    if !first_line {
                        out.push_str("\n ");
                    }
                    first_line = false;
                    for &c in chunk {
                        push_touchstone_pair(&mut out, c, format)?;
                    }
                }
            }
        }
        out.push('\n');
    }

    Ok(out)
}

/// Format Touchstone output string.
///
/// Returns the same text as [`try_write_touchstone`] for consistent data.
/// Because this signature cannot report errors, inconsistent data produces
/// an empty string (and a `tracing` warning) instead of panicking.
#[deprecated(
    since = "1.2.0",
    note = "cannot report inconsistent data; use try_write_touchstone"
)]
#[must_use]
pub fn write_touchstone(data: &TouchstoneData, format: TouchstoneFormat) -> String {
    try_write_touchstone(data, format).unwrap_or_else(|e| {
        tracing::warn!(error = %e, "write_touchstone: inconsistent data, nothing written");
        String::new()
    })
}

fn parse_complex_pair(v1: f64, v2: f64, format: TouchstoneFormat) -> Complex {
    match format {
        TouchstoneFormat::RealImag => Complex::new(v1, v2),
        TouchstoneFormat::MagnitudeAngle => complex_from_polar(v1, v2.to_radians()),
        TouchstoneFormat::DbAngle => {
            let mag = 10.0_f64.powf(v1 / 20.0);
            complex_from_polar(mag, v2.to_radians())
        }
    }
}

fn format_complex_pair(c: Complex, format: TouchstoneFormat) -> (f64, f64) {
    match format {
        TouchstoneFormat::RealImag => (c.re, c.im),
        TouchstoneFormat::MagnitudeAngle => (c.norm(), c.arg().to_degrees()),
        TouchstoneFormat::DbAngle => (
            20.0 * c.norm().max(f64::MIN_POSITIVE).log10(),
            c.arg().to_degrees(),
        ),
    }
}

// ── Antenna arrays ────────────────────────────────────────────────

/// Array factor for a uniform linear array.
///
/// AF(θ) = sin(Nψ/2) / (N sin(ψ/2)) where ψ = kd cos θ + β
///
/// - `n_elements`: number of array elements
/// - `spacing`: element spacing d (meters)
/// - `wavenumber`: k = 2π/λ
/// - `progressive_phase`: β (radians), phase shift between elements
/// - `theta`: observation angle from array axis (radians)
#[inline]
#[must_use]
pub fn linear_array_factor(
    n_elements: usize,
    spacing: f64,
    wavenumber: f64,
    progressive_phase: f64,
    theta: f64,
) -> f64 {
    if n_elements == 0 {
        return 0.0;
    }
    if n_elements == 1 {
        return 1.0;
    }

    let psi = wavenumber * spacing * theta.cos() + progressive_phase;
    let half_psi = psi / 2.0;
    let n = n_elements as f64;

    // AF = sin(N ψ/2) / (N sin(ψ/2))
    let sin_n_half_psi = (n * half_psi).sin();
    let sin_half_psi = half_psi.sin();

    if sin_half_psi.abs() < 1e-15 {
        return 1.0; // main beam peak (L'Hôpital)
    }

    (sin_n_half_psi / (n * sin_half_psi)).abs()
}

/// Array factor for a planar (rectangular) array.
///
/// AF = AF_x(θ,φ) × AF_y(θ,φ)
///
/// Returns the product of x-axis and y-axis linear array factors.
#[inline]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn planar_array_factor(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    wavenumber: f64,
    beta_x: f64,
    beta_y: f64,
    theta: f64,
    phi: f64,
) -> f64 {
    let af_x = linear_array_factor(
        nx,
        dx,
        wavenumber,
        beta_x,
        // Effective angle for x-axis: project onto x
        (theta.sin() * phi.cos()).acos(),
    );
    let af_y = linear_array_factor(
        ny,
        dy,
        wavenumber,
        beta_y,
        // Effective angle for y-axis: project onto y
        (theta.sin() * phi.sin()).acos(),
    );
    af_x * af_y
}

/// Array factor for a uniform circular array.
///
/// N elements equally spaced on a circle of radius a.
/// AF(θ,φ) = Σ exp(j(ka sin θ cos(φ − φₙ) + βₙ))
///
/// Returns the magnitude of the array factor, normalized to N.
#[must_use]
pub fn circular_array_factor(
    n_elements: usize,
    radius: f64,
    wavenumber: f64,
    theta: f64,
    phi: f64,
    progressive_phase: f64,
) -> f64 {
    if n_elements == 0 {
        return 0.0;
    }

    let n = n_elements as f64;
    let mut sum_re = 0.0;
    let mut sum_im = 0.0;

    for i in 0..n_elements {
        let phi_n = 2.0 * std::f64::consts::PI * i as f64 / n;
        let psi =
            wavenumber * radius * theta.sin() * (phi - phi_n).cos() + progressive_phase * i as f64;
        sum_re += psi.cos();
        sum_im += psi.sin();
    }

    (sum_re * sum_re + sum_im * sum_im).sqrt() / n
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    fn close(a: Complex, b: Complex, tol: f64) -> bool {
        (a.re - b.re).abs() <= tol && (a.im - b.im).abs() <= tol
    }

    // ── Transmission line tests ───────────────────────────────────

    #[test]
    fn test_z0_lossless() {
        // Coax: L'=250nH/m, C'=100pF/m → Z₀ = 50Ω
        let z0 = characteristic_impedance_lossless(250e-9, 100e-12).unwrap();
        assert!((z0 - 50.0).abs() < 0.1);
    }

    #[test]
    fn test_z0_invalid() {
        assert!(characteristic_impedance_lossless(-1.0, 100e-12).is_err());
        assert!(characteristic_impedance_lossless(250e-9, 0.0).is_err());
        assert!(characteristic_impedance_lossless(f64::NAN, 100e-12).is_err());
        assert!(phase_velocity(f64::NAN, 100e-12).is_err());
    }

    #[test]
    fn test_phase_velocity() {
        let vp = phase_velocity(250e-9, 100e-12).unwrap();
        // v = 1/√(LC) ≈ 2e8 m/s for coax with polyethylene dielectric
        assert!(vp > 1e8);
        assert!(vp < 3e8);
    }

    #[test]
    #[allow(deprecated)]
    fn test_input_impedance_matched() {
        // Matched load: Z_in = Z₀ regardless of length
        let z_in = input_impedance_lossless(50.0, 50.0, 10.0, 0.1).unwrap();
        assert!((z_in - 50.0).abs() < 0.01);
    }

    #[test]
    #[allow(deprecated)]
    fn test_input_impedance_quarter_wave() {
        // λ/4 transformer: Z_in = Z₀²/Z_L
        let beta = 2.0 * std::f64::consts::PI; // β = 2π/λ, so βl = π/2 at l = λ/4
        let z_in = input_impedance_lossless(50.0, 100.0, beta, 0.25).unwrap();
        // Z_in = 50² / 100 = 25
        assert!((z_in - 25.0).abs() < 0.1);
    }

    #[test]
    fn test_input_impedance_complex_eighth_wave() {
        // Z₀=50, Z_L=100, βl=π/4: Z_in = 50(100 + j50)/(50 + j100) = 40 − j30
        let z = input_impedance_lossless_complex(
            50.0,
            complex_real(100.0),
            1.0,
            std::f64::consts::FRAC_PI_4,
        )
        .unwrap();
        assert!(close(z, Complex::new(40.0, -30.0), 1e-9), "got {z}");
    }

    #[test]
    fn test_input_impedance_complex_load() {
        // Reference (python, tan form): Z₀=50, Z_L=30+j40, βl=0.6π
        let z = input_impedance_lossless_complex(50.0, Complex::new(30.0, 40.0), TAU, 0.3).unwrap();
        let expected = Complex::new(20.404995607197357, -22.010650955158418);
        assert!(close(z, expected, 1e-9), "got {z}");
    }

    #[test]
    fn test_input_impedance_short_quarter_wave_is_unbounded() {
        let r = input_impedance_lossless_complex(50.0, complex_zero(), TAU, 0.25);
        assert!(matches!(r, Err(BijliError::DivisionByZero { .. })));
    }

    #[test]
    fn test_input_impedance_complex_invalid() {
        assert!(input_impedance_lossless_complex(0.0, complex_real(50.0), 1.0, 1.0).is_err());
        assert!(input_impedance_lossless_complex(f64::NAN, complex_real(50.0), 1.0, 1.0).is_err());
        assert!(
            input_impedance_lossless_complex(50.0, Complex::new(f64::NAN, 0.0), 1.0, 1.0).is_err()
        );
        assert!(
            input_impedance_lossless_complex(50.0, complex_real(50.0), f64::INFINITY, 1.0).is_err()
        );
    }

    #[test]
    #[allow(deprecated)]
    fn test_input_impedance_real_rejects_complex_result() {
        // The true Z_in is 40 − j30; the old real-valued version returned 50.
        let r = input_impedance_lossless(50.0, 100.0, 1.0, std::f64::consts::FRAC_PI_4);
        assert!(matches!(r, Err(BijliError::InvalidParameter { .. })));
    }

    // ── Reflection / VSWR tests ───────────────────────────────────

    #[test]
    fn test_reflection_matched() {
        let gamma = reflection_coefficient(50.0, 50.0).unwrap();
        assert!(gamma.abs() < TOL);
    }

    #[test]
    fn test_reflection_open() {
        // Open circuit: Z_L = ∞ → Γ = 1
        let gamma = reflection_coefficient(1e15, 50.0).unwrap();
        assert!((gamma - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_reflection_short() {
        // Short circuit: Z_L = 0 → Γ = -1
        let gamma = reflection_coefficient(0.0, 50.0).unwrap();
        assert!((gamma - (-1.0)).abs() < TOL);
    }

    #[test]
    fn test_reflection_nan_rejected() {
        assert!(reflection_coefficient(f64::NAN, 50.0).is_err());
        assert!(
            reflection_coefficient_complex(Complex::new(f64::NAN, 0.0), complex_real(50.0))
                .is_err()
        );
    }

    #[test]
    fn test_vswr_matched() {
        let v = vswr(0.0).unwrap();
        assert!((v - 1.0).abs() < TOL);
    }

    #[test]
    fn test_vswr_2_to_1() {
        // |Γ| = 1/3 → VSWR = 2:1
        let v = vswr(1.0 / 3.0).unwrap();
        assert!((v - 2.0).abs() < TOL);
    }

    #[test]
    fn test_vswr_invalid() {
        assert!(vswr(1.0).is_err());
        assert!(vswr(-0.1).is_err());
        assert!(vswr(f64::NAN).is_err());
    }

    #[test]
    fn test_return_loss() {
        // |Γ| = 0.1 → RL = 20 dB
        let rl = return_loss(0.1).unwrap();
        assert!((rl - 20.0).abs() < 0.01);
        // Total reflection: 0 dB
        assert!(return_loss(1.0).unwrap().abs() < TOL);
    }

    #[test]
    fn test_return_loss_invalid() {
        assert!(return_loss(f64::NAN).is_err());
        assert!(return_loss(1.5).is_err());
    }

    #[test]
    fn test_insertion_loss_matched() {
        let il = insertion_loss(0.0).unwrap();
        assert!(il.abs() < TOL); // 0 dB for perfect match
    }

    #[test]
    fn test_insertion_loss_invalid() {
        assert!(insertion_loss(f64::NAN).is_err());
        assert!(insertion_loss(-0.5).is_err());
        assert!(mismatch_loss(f64::NAN).is_err());
    }

    // ── Smith chart tests ─────────────────────────────────────────

    #[test]
    fn test_normalize_denormalize() {
        let z = Complex::new(100.0, 50.0);
        let z_n = normalize_impedance(z, 50.0).unwrap();
        let z_back = denormalize_impedance(z_n, 50.0);
        assert!((z_back.re - z.re).abs() < TOL);
        assert!((z_back.im - z.im).abs() < TOL);
    }

    #[test]
    fn test_gamma_z_roundtrip() {
        let z = Complex::new(2.0, 1.0); // normalized
        let g = gamma_from_z(z);
        let z_back = z_from_gamma(g);
        assert!((z_back.re - z.re).abs() < TOL);
        assert!((z_back.im - z.im).abs() < TOL);
    }

    #[test]
    fn test_gamma_matched() {
        let g = gamma_from_z(complex_real(1.0));
        assert!(g.norm() < TOL);
    }

    #[test]
    fn test_smith_poles_are_infinite_not_nan() {
        let g = gamma_from_z(complex_real(-1.0));
        assert!(g.re.is_infinite() && !g.im.is_nan(), "got {g}");
        let z = z_from_gamma(complex_real(1.0));
        assert!(z.re.is_infinite() && !z.im.is_nan(), "got {z}");
        // Open circuit z = ∞ ↔ Γ = 1, and the inverse pole Γ = ∞ ↔ z = −1.
        assert!(close(gamma_from_z(z), complex_real(1.0), 0.0));
        assert!(close(z_from_gamma(g), complex_real(-1.0), 0.0));
    }

    #[test]
    fn test_lossy_z0_dc_limits() {
        // R = G = ω = 0: lossless DC limit √(L/C) = 50 Ω (was 0/0 = NaN).
        let z = characteristic_impedance_lossy(0.0, 250e-9, 0.0, 100e-12, 0.0);
        assert!(close(z, complex_real(50.0), 1e-9), "got {z}");
        // R > 0, G = ω = 0: |Z₀| → ∞.
        let z = characteristic_impedance_lossy(1.0, 250e-9, 0.0, 100e-12, 0.0);
        assert!(z.re.is_infinite(), "got {z}");
        // G > 0 at DC: √(R/G).
        let z = characteristic_impedance_lossy(4.0, 250e-9, 1.0, 100e-12, 0.0);
        assert!(close(z, complex_real(2.0), 1e-12), "got {z}");
    }

    // ── Matching tests ────────────────────────────────────────────

    #[test]
    fn test_quarter_wave_transformer() {
        // Match 50Ω to 100Ω: Z_t = √(5000) ≈ 70.71Ω
        let zt = quarter_wave_transformer(50.0, 100.0).unwrap();
        assert!((zt - 70.71).abs() < 0.01);
    }

    #[test]
    fn test_quarter_wave_both_negative() {
        assert!(quarter_wave_transformer(-50.0, -100.0).is_err());
        assert!(quarter_wave_transformer(50.0, 0.0).is_err());
        assert!(quarter_wave_transformer(f64::NAN, 100.0).is_err());
    }

    #[test]
    fn test_single_stub_already_matched() {
        // Matched: no stub needed, i.e. a λ/4 short stub (zero admittance).
        let (d, l) = single_stub_match(50.0, 50.0).unwrap();
        assert!(d.abs() < TOL);
        assert!((l - 0.25).abs() < TOL);
    }

    #[test]
    fn test_single_stub_real_load_reference() {
        // Reference values from an independent brute-force root search of
        // Re Y(d) = Y₀ (python): 100 Ω on 50 Ω → (0.15204, 0.15204).
        let (d, l) = single_stub_match(100.0, 50.0).unwrap();
        assert!((d - 0.1520433619923482).abs() < 1e-9, "d = {d}");
        assert!((l - 0.1520433619923482).abs() < 1e-9, "l = {l}");
        // 25 Ω on 50 Ω → (0.09796, 0.34796).
        let (d, l) = single_stub_match(25.0, 50.0).unwrap();
        assert!((d - 0.09795663800765184).abs() < 1e-9, "d = {d}");
        assert!((l - 0.3479566380076518).abs() < 1e-9, "l = {l}");
    }

    #[test]
    fn test_single_stub_pozar_example() {
        // Pozar Example 5.2: Z_L = 60 − j80 Ω, Z₀ = 50 Ω, short-circuited
        // shunt stubs: d₁ = 0.110λ, ℓ₁ = 0.095λ; d₂ = 0.260λ, ℓ₂ = 0.405λ.
        let [(d1, l1), (d2, l2)] =
            single_stub_match_complex(Complex::new(60.0, -80.0), 50.0).unwrap();
        assert!((d1 - 0.11042321863830029).abs() < 1e-9);
        assert!((l1 - 0.09497462163589147).abs() < 1e-9);
        assert!((d2 - 0.2594445306228258).abs() < 1e-9);
        assert!((l2 - 0.4050253783641085).abs() < 1e-9);
    }

    #[test]
    fn test_single_stub_admittance_is_matched() {
        // Independent check with the tan form: Y(d) + Y_stub(ℓ) = Y₀.
        let z0 = 50.0;
        let y0 = 1.0 / z0;
        for zl in [
            Complex::new(100.0, 0.0),
            Complex::new(25.0, 0.0),
            Complex::new(60.0, -80.0),
            Complex::new(50.0, 50.0),
            Complex::new(10.0, 120.0),
        ] {
            for (d, l) in single_stub_match_complex(zl, z0).unwrap() {
                assert!((0.0..0.5).contains(&d) && l > 0.0 && l < 0.5);
                let t = (TAU * d).tan();
                let j = Complex::new(0.0, 1.0);
                let z_in = (zl + j * (z0 * t)) / (complex_real(z0) + j * zl * t) * z0;
                let y_line = complex_real(1.0) / z_in;
                let y_stub = Complex::new(0.0, -y0 / (TAU * l).tan());
                let y = y_line + y_stub;
                assert!(close(y, complex_real(y0), 1e-9), "Z_L={zl}: Y={y}");
            }
        }
    }

    #[test]
    fn test_single_stub_invalid() {
        assert!(single_stub_match(-100.0, 50.0).is_err());
        assert!(single_stub_match(100.0, 0.0).is_err());
        assert!(single_stub_match(f64::NAN, 50.0).is_err());
        assert!(single_stub_match_complex(Complex::new(0.0, 30.0), 50.0).is_err());
    }

    // ── S-parameter tests ─────────────────────────────────────────

    #[test]
    fn test_smatrix_creation() {
        let s = SMatrix::new(2, 50.0);
        assert_eq!(s.ports, 2);
        assert!(s.get(1, 1).unwrap().norm() < TOL);
    }

    #[test]
    fn test_smatrix_get_set() {
        let mut s = SMatrix::new(2, 50.0);
        s.set(1, 2, Complex::new(0.5, 0.1)).unwrap();
        let v = s.get(1, 2).unwrap();
        assert!((v.re - 0.5).abs() < TOL);
    }

    #[test]
    fn test_smatrix_invalid_index() {
        let s = SMatrix::new(2, 50.0);
        assert!(s.get(0, 1).is_err());
        assert!(s.get(3, 1).is_err());
    }

    #[test]
    fn test_smatrix_inconsistent_does_not_panic() {
        let mut bad = SMatrix {
            ports: 2,
            s: vec![],
            z0: vec![50.0; 2],
        };
        assert!(bad.get(1, 1).is_err());
        assert!(bad.set(2, 2, complex_real(1.0)).is_err());
        assert!(bad.cascade(&SMatrix::new(2, 50.0)).is_err());
        assert!(SMatrix::new(2, 50.0).cascade(&bad).is_err());
        assert!(s_to_z(&bad).is_err());
        let ragged = SMatrix {
            ports: 2,
            s: vec![vec![complex_real(0.1)], vec![complex_real(0.2)]],
            z0: vec![],
        };
        assert!(s_to_z(&ragged).is_err());
        assert!(ragged.cascade(&ragged).is_err());
    }

    #[test]
    fn test_cascade_thru() {
        // Two ideal thru connections: S21 = S12 = 1, S11 = S22 = 0
        let thru = SMatrix::two_port(
            complex_zero(),
            complex_real(1.0),
            complex_real(1.0),
            complex_zero(),
            50.0,
        );
        let cascaded = thru.cascade(&thru).unwrap();
        // Thru × Thru = Thru
        assert!((cascaded.s[1][0].re - 1.0).abs() < TOL); // S21 = 1
        assert!(cascaded.s[0][0].norm() < TOL); // S11 = 0
    }

    #[test]
    fn test_cascade_matches_star_product() {
        // Independent reference: Redheffer star product (python),
        // S11 = A11 + A12 A21 B11 / (1 − A22 B11), S21 = A21 B21 / (1 − A22 B11), …
        let a = SMatrix::two_port(
            Complex::new(0.1, 0.2),
            Complex::new(0.3, -0.1),
            Complex::new(0.8, 0.05),
            Complex::new(-0.2, 0.1),
            50.0,
        );
        let b = SMatrix::two_port(
            Complex::new(-0.15, 0.05),
            Complex::new(0.6, 0.2),
            Complex::new(0.7, -0.3),
            Complex::new(0.25, -0.05),
            50.0,
        );
        let c = a.cascade(&b).unwrap();
        let tol = 1e-12;
        assert!(close(
            c.s[0][0],
            Complex::new(0.06624178712220763, 0.22342969776609725),
            tol
        ));
        assert!(close(
            c.s[0][1],
            Complex::new(0.20499342969776613, -0.005256241787122209),
            tol
        ));
        assert!(close(
            c.s[1][0],
            Complex::new(0.5839684625492773, -0.2252299605781866),
            tol
        ));
        assert!(close(
            c.s[1][1],
            Complex::new(0.15717477003942182, 0.009816031537450719),
            tol
        ));
        assert_eq!(c.z0, vec![50.0, 50.0]);
    }

    #[test]
    fn test_cascade_requires_equal_z0() {
        let a = SMatrix::two_port(
            complex_zero(),
            complex_real(1.0),
            complex_real(1.0),
            complex_zero(),
            50.0,
        );
        let mut mixed = a.clone();
        mixed.z0 = vec![50.0, 75.0];
        assert!(mixed.cascade(&a).is_err());
        assert!(a.cascade(&mixed).is_err());
        let mut other = a.clone();
        other.z0 = vec![75.0, 75.0];
        assert!(a.cascade(&other).is_err());
    }

    #[test]
    fn test_s_to_z_conversion() {
        // 3dB attenuator: S11=S22=0, S12=S21=1/√2 in 50Ω
        let s = SMatrix::two_port(
            complex_zero(),
            complex_real(std::f64::consts::FRAC_1_SQRT_2),
            complex_real(std::f64::consts::FRAC_1_SQRT_2),
            complex_zero(),
            50.0,
        );
        let z = s_to_z(&s).unwrap();
        // Just verify it produces finite values and Z12 = Z21 (reciprocal)
        assert!(z[0][0].norm() > 0.0);
        assert!((z[0][1].re - z[1][0].re).abs() < TOL);
        assert!((z[0][1].im - z[1][0].im).abs() < TOL);
    }

    #[test]
    fn test_s_to_z_shunt_resistor() {
        // Shunt R = 100 Ω in 50 Ω: S11 = S22 = −Z₀/(2R+Z₀) = −0.2,
        // S21 = S12 = 2R/(2R+Z₀) = 0.8; every Z-parameter equals R.
        let s = SMatrix::two_port(
            complex_real(-0.2),
            complex_real(0.8),
            complex_real(0.8),
            complex_real(-0.2),
            50.0,
        );
        let z = s_to_z(&s).unwrap();
        for row in z {
            for zij in row {
                assert!(close(zij, complex_real(100.0), 1e-9), "got {zij}");
            }
        }
    }

    #[test]
    fn test_s_to_z_requires_equal_z0() {
        let mut s = SMatrix::two_port(
            complex_real(-0.2),
            complex_real(0.8),
            complex_real(0.8),
            complex_real(-0.2),
            50.0,
        );
        s.z0 = vec![50.0, 75.0];
        assert!(s_to_z(&s).is_err());
    }

    #[test]
    fn test_s_to_z_thru_is_singular() {
        // Ideal thru has no Z-parameter representation
        let thru = SMatrix::two_port(
            complex_zero(),
            complex_real(1.0),
            complex_real(1.0),
            complex_zero(),
            50.0,
        );
        assert!(s_to_z(&thru).is_err());
    }

    // ── Touchstone tests ──────────────────────────────────────────

    fn one_port_data() -> TouchstoneData {
        let mut a = SMatrix::new(1, 50.0);
        a.s[0][0] = Complex::new(0.5, 0.1);
        let mut b = SMatrix::new(1, 50.0);
        b.s[0][0] = Complex::new(0.3, -0.2);
        TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1e9, 2e9],
            s_matrices: vec![a, b],
        }
    }

    #[test]
    fn test_touchstone_roundtrip() {
        let data = one_port_data();
        for format in [
            TouchstoneFormat::RealImag,
            TouchstoneFormat::MagnitudeAngle,
            TouchstoneFormat::DbAngle,
        ] {
            let content = try_write_touchstone(&data, format).unwrap();
            assert!(content.contains("# GHz S"));
            assert!(content.contains("R 50\n"));
            let back = parse_touchstone(&content, 1).unwrap();
            assert_eq!(back.frequencies, data.frequencies);
            for (x, y) in back.s_matrices.iter().zip(&data.s_matrices) {
                assert!(close(x.s[0][0], y.s[0][0], 1e-12), "{format:?}");
            }
        }
    }

    #[test]
    #[allow(deprecated)]
    fn test_write_touchstone_deprecated_matches_try() {
        let data = one_port_data();
        let a = write_touchstone(&data, TouchstoneFormat::RealImag);
        let b = try_write_touchstone(&data, TouchstoneFormat::RealImag).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    #[allow(deprecated)]
    fn test_write_touchstone_inconsistent_does_not_panic() {
        let data = TouchstoneData {
            ports: 2,
            z0: 50.0,
            frequencies: vec![1e9],
            s_matrices: vec![],
        };
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());
        assert!(write_touchstone(&data, TouchstoneFormat::RealImag).is_empty());

        // Matrix size disagrees with `ports`.
        let data = TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1e9],
            s_matrices: vec![SMatrix::new(2, 50.0)],
        };
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());

        // Matrix reference impedance disagrees with the file's R.
        let data = TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1e9],
            s_matrices: vec![SMatrix::new(1, 75.0)],
        };
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());

        // ports = 0 and non-finite values.
        let mut data = one_port_data();
        data.ports = 0;
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());
        let mut data = one_port_data();
        data.s_matrices[0].s[0][0] = Complex::new(f64::NAN, 0.0);
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());
        let mut data = one_port_data();
        data.frequencies = vec![2e9, 1e9];
        assert!(try_write_touchstone(&data, TouchstoneFormat::RealImag).is_err());
    }

    #[test]
    fn test_parse_touchstone_ri() {
        let content = "# GHz S RI R 50\n1.0  0.5  0.1\n2.0  0.3  -0.2\n";
        let data = parse_touchstone(content, 1).unwrap();
        assert_eq!(data.frequencies.len(), 2);
        assert!((data.frequencies[0] - 1e9).abs() < 1.0);
        assert!((data.s_matrices[0].s[0][0].re - 0.5).abs() < TOL);
    }

    #[test]
    fn test_parse_touchstone_empty() {
        assert!(parse_touchstone("! comment only\n", 1).is_err());
    }

    #[test]
    fn test_parse_touchstone_two_port_order() {
        // v1.1 2-port order is f N11 N21 N12 N22 (hand-written, S21 ≠ S12).
        let content = "! 2-port\n# GHZ S RI R 50\n1.0  0.1 0.01  0.9 0.09  0.2 0.02  0.3 0.03\n";
        let data = parse_touchstone(content, 2).unwrap();
        let s = &data.s_matrices[0];
        assert!(close(s.get(1, 1).unwrap(), Complex::new(0.1, 0.01), 1e-15));
        assert!(close(s.get(2, 1).unwrap(), Complex::new(0.9, 0.09), 1e-15));
        assert!(close(s.get(1, 2).unwrap(), Complex::new(0.2, 0.02), 1e-15));
        assert!(close(s.get(2, 2).unwrap(), Complex::new(0.3, 0.03), 1e-15));
    }

    #[test]
    fn test_write_touchstone_two_port_order() {
        let sm = SMatrix::two_port(
            complex_real(0.1),
            complex_real(0.2),
            complex_real(0.9),
            complex_real(0.3),
            50.0,
        );
        let data = TouchstoneData {
            ports: 2,
            z0: 50.0,
            frequencies: vec![1e9],
            s_matrices: vec![sm],
        };
        let text = try_write_touchstone(&data, TouchstoneFormat::RealImag).unwrap();
        let line = text.lines().last().unwrap();
        let v: Vec<f64> = line
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        // f, N11, N21, N12, N22
        assert_eq!(v, vec![1.0, 0.1, 0.0, 0.9, 0.0, 0.2, 0.0, 0.3, 0.0]);
    }

    #[test]
    fn test_parse_touchstone_three_port_multiline() {
        // v1.1: one matrix row per line for 3- and 4-ports.
        let content = "# GHz S RI R 50\n\
            1.0  0.11 -0.11  0.12 -0.12  0.13 -0.13\n\
                 0.21 -0.21  0.22 -0.22  0.23 -0.23\n\
                 0.31 -0.31  0.32 -0.32  0.33 -0.33 ! row 3\n\
            2.0  0.11 0 0.12 0 0.13 0\n 0.21 0 0.22 0 0.23 0\n 0.31 0 0.32 0 0.33 0\n";
        let data = parse_touchstone(content, 3).unwrap();
        assert_eq!(data.frequencies, vec![1e9, 2e9]);
        let s = &data.s_matrices[0];
        for i in 0..3 {
            for j in 0..3 {
                let v = 0.1 * (i + 1) as f64 + 0.01 * (j + 1) as f64;
                assert!(
                    close(s.s[i][j], Complex::new(v, -v), 1e-12),
                    "S{}{}",
                    i + 1,
                    j + 1
                );
            }
        }
    }

    #[test]
    fn test_parse_touchstone_five_port_wrapped() {
        // ≥ 5 ports: each row wraps after four pairs.
        let n = 5;
        let mut content = String::from("# MHz S MA R 75\n");
        for (fi, f) in [100.0, 200.0].iter().enumerate() {
            for i in 0..n {
                let pairs: Vec<String> = (0..n)
                    .map(|j| format!("{} {}", (10 * (i + 1) + j + 1) as f64 / 100.0, fi as f64))
                    .collect();
                let head = if i == 0 {
                    format!("{f}")
                } else {
                    String::new()
                };
                content.push_str(&format!("{head} {}\n", pairs[..4].join(" ")));
                content.push_str(&format!("   {}\n", pairs[4]));
            }
        }
        let data = parse_touchstone(&content, n).unwrap();
        assert_eq!(data.frequencies, vec![100e6, 200e6]);
        assert!((data.z0 - 75.0).abs() < TOL);
        let s = &data.s_matrices[0];
        assert!(close(s.s[2][4], complex_real(0.35), 1e-12));
        assert!(close(s.s[4][0], complex_real(0.51), 1e-12));
    }

    #[test]
    fn test_touchstone_multiport_write_layout_and_roundtrip() {
        for n in [3usize, 4, 5, 9] {
            let mut sm = SMatrix::new(n, 50.0);
            for i in 0..n {
                for j in 0..n {
                    sm.s[i][j] = Complex::new(0.01 * (i * n + j) as f64, -0.003 * j as f64);
                }
            }
            let data = TouchstoneData {
                ports: n,
                z0: 50.0,
                frequencies: vec![1.5e9],
                s_matrices: vec![sm.clone()],
            };
            let text = try_write_touchstone(&data, TouchstoneFormat::RealImag).unwrap();
            let lines_per_row = n.div_ceil(4);
            let data_lines: Vec<&str> = text.lines().skip(2).collect();
            assert_eq!(data_lines.len(), n * lines_per_row, "n = {n}");
            for line in &data_lines {
                let tokens = line.split_whitespace().count();
                assert!(tokens <= 1 + 8, "n = {n}: {line}");
            }
            let back = parse_touchstone(&text, n).unwrap();
            assert_eq!(back.s_matrices[0], sm, "n = {n}");
        }
    }

    #[test]
    fn test_touchstone_write_precision() {
        // Low frequency and tiny S values must survive (were "0.000001" / "0.000000").
        let mut sm = SMatrix::new(1, 50.0);
        sm.s[0][0] = Complex::new(1e-9, -2e-12);
        let data = TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1234.5, 1e6 + 0.25],
            s_matrices: vec![sm.clone(), sm],
        };
        for format in [
            TouchstoneFormat::RealImag,
            TouchstoneFormat::MagnitudeAngle,
            TouchstoneFormat::DbAngle,
        ] {
            let text = try_write_touchstone(&data, format).unwrap();
            let back = parse_touchstone(&text, 1).unwrap();
            for (f, g) in back.frequencies.iter().zip(&data.frequencies) {
                assert!((f - g).abs() <= 1e-12 * g, "{format:?}: {f} vs {g}");
            }
            let s = back.s_matrices[0].s[0][0];
            assert!(
                (s.re - 1e-9).abs() <= 1e-21 && (s.im + 2e-12).abs() <= 1e-21,
                "{format:?}: {s}"
            );
        }
    }

    #[test]
    fn test_touchstone_db_zero_magnitude_is_finite() {
        let data = TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1e9],
            s_matrices: vec![SMatrix::new(1, 50.0)],
        };
        let text = try_write_touchstone(&data, TouchstoneFormat::DbAngle).unwrap();
        assert!(!text.to_lowercase().contains("inf"), "{text}");
        let back = parse_touchstone(&text, 1).unwrap();
        assert!(back.s_matrices[0].s[0][0].norm() < 1e-300);
    }

    #[test]
    fn test_parse_touchstone_inline_comment() {
        let data = parse_touchstone("# GHz S RI R 50 ! options\n1.0 0.5 0.1 ! first\n", 1).unwrap();
        assert!(close(
            data.s_matrices[0].s[0][0],
            Complex::new(0.5, 0.1),
            1e-15
        ));
    }

    #[test]
    fn test_parse_touchstone_rejects_bad_tokens() {
        // A non-numeric token used to be dropped, silently shifting values.
        assert!(parse_touchstone("# GHz S RI R 50\n1.0 0.5 abc 0.1\n", 1).is_err());
        assert!(parse_touchstone("# GHz S RI R 50\n1.0 nan 0.1\n", 1).is_err());
        assert!(parse_touchstone("# GHz S RI R 50\n1.0 inf 0.1\n", 1).is_err());
        assert!(parse_touchstone("# GHz S RI R 50\n-1.0 0.5 0.1\n", 1).is_err());
        assert!(parse_touchstone("# GHz S DB R 50\n1.0 7000 0\n", 1).is_err());
    }

    #[test]
    fn test_parse_touchstone_partial_record() {
        assert!(parse_touchstone("# GHz S RI R 50\n1.0 0.5 0.1\n2.0 0.3\n", 1).is_err());
        assert!(parse_touchstone("# GHz S RI R 50\n1.0 0.1 0 0.9 0 0.2 0\n", 2).is_err());
    }

    #[test]
    fn test_parse_touchstone_option_line_rules() {
        // Only the first option line counts.
        let d = parse_touchstone("# GHz S RI R 50\n# MHz S MA R 75\n1.0 0.5 0.1\n", 1).unwrap();
        assert!((d.z0 - 50.0).abs() < TOL);
        assert!((d.frequencies[0] - 1e9).abs() < 1.0);
        assert!(close(
            d.s_matrices[0].s[0][0],
            Complex::new(0.5, 0.1),
            1e-15
        ));
        // Non-S parameter types.
        for p in ["Y", "Z", "H", "G"] {
            let content = format!("# GHz {p} RI R 50\n1.0 0.5 0.1\n");
            assert!(parse_touchstone(&content, 1).is_err(), "{p}");
        }
        // R must be followed by a finite positive number.
        for r in ["", "abc", "-50", "0", "nan", "inf"] {
            let content = format!("# GHz S RI R {r}\n1.0 0.5 0.1\n");
            assert!(parse_touchstone(&content, 1).is_err(), "R {r:?}");
        }
        // Unknown option token, option line after data, v2 keywords.
        assert!(parse_touchstone("# GHz S RI XYZ R 50\n1.0 0.5 0.1\n", 1).is_err());
        assert!(parse_touchstone("1.0 0.5 0.1\n# GHz S RI R 50\n2.0 0.5 0.1\n", 1).is_err());
        assert!(parse_touchstone("[Version] 2.0\n# GHz S RI R 50\n1.0 0.5 0.1\n", 1).is_err());
        // Defaults: GHz S MA R 50.
        let d = parse_touchstone("1.0 0.5 90\n", 1).unwrap();
        assert!((d.frequencies[0] - 1e9).abs() < 1.0);
        assert!(close(
            d.s_matrices[0].s[0][0],
            Complex::new(0.0, 0.5),
            1e-12
        ));
    }

    #[test]
    fn test_parse_touchstone_port_limits() {
        assert!(parse_touchstone("1.0 0.5 0.1\n", 0).is_err());
        assert!(parse_touchstone("1.0 0.5 0.1\n", 100_000).is_err());
        assert!(parse_touchstone("1.0 0.5 0.1\n", usize::MAX).is_err());
    }

    #[test]
    fn test_parse_touchstone_two_port_noise_block_skipped() {
        let content = "# GHZ S MA R 50\n\
            1.0 0.5 10 0.9 -20 0.1 30 0.4 40\n\
            2.0 0.6 11 0.8 -21 0.1 31 0.3 41\n\
            ! noise parameters\n\
            1.0 1.5 0.3 45 0.2\n\
            2.0 1.8 0.35 50 0.25\n";
        let data = parse_touchstone(content, 2).unwrap();
        assert_eq!(data.frequencies, vec![1e9, 2e9]);
        // MA: S21 at 1 GHz is 0.9∠−20°.
        let s21 = data.s_matrices[0].get(2, 1).unwrap();
        assert!((s21.norm() - 0.9).abs() < 1e-12);
        assert!((s21.arg().to_degrees() + 20.0).abs() < 1e-9);
    }

    // ── Array factor tests ────────────────────────────────────────

    #[test]
    fn test_linear_array_broadside() {
        // Broadside array (β=0): peak at θ = 90°
        let k = 2.0 * std::f64::consts::PI;
        let af = linear_array_factor(4, 0.5, k, 0.0, std::f64::consts::FRAC_PI_2);
        assert!((af - 1.0).abs() < TOL);
    }

    #[test]
    fn test_linear_array_single_element() {
        assert!((linear_array_factor(1, 0.5, 10.0, 0.0, 0.5) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_linear_array_zero_elements() {
        assert!(linear_array_factor(0, 0.5, 10.0, 0.0, 0.5).abs() < TOL);
    }

    #[test]
    fn test_circular_array_single() {
        assert!((circular_array_factor(1, 0.5, 10.0, 0.5, 0.0, 0.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_circular_array_broadside_peak() {
        // Uniform circular array with β=0: peak at θ=0 (endfire)
        let af = circular_array_factor(8, 0.25, 2.0 * std::f64::consts::PI, 0.0, 0.0, 0.0);
        assert!((af - 1.0).abs() < TOL); // all elements in phase at θ=0
    }

    #[test]
    fn test_planar_array_nonzero() {
        let k = 2.0 * std::f64::consts::PI;
        let af = planar_array_factor(4, 4, 0.5, 0.5, k, 0.0, 0.0, 0.3, 0.5);
        // Should produce a finite positive value
        assert!(af > 0.0);
        assert!(af <= 1.0);
    }

    // ── Audit-driven edge case tests ──────────────────────────────

    #[test]
    fn test_normalize_zero_z0() {
        assert!(normalize_impedance(Complex::new(50.0, 0.0), 0.0).is_err());
        assert!(normalize_impedance(Complex::new(50.0, 0.0), f64::NAN).is_err());
    }

    #[test]
    fn test_reflection_coefficient_complex_matched() {
        let z = Complex::new(50.0, 0.0);
        let g = reflection_coefficient_complex(z, z).unwrap();
        assert!(g.norm() < TOL);
    }

    #[test]
    fn test_return_loss_zero_gamma() {
        // |Γ| = 0 → RL = ∞, but log(0) = -inf, so we error
        assert!(return_loss(0.0).is_err());
    }

    #[test]
    fn test_insertion_loss_total_reflection() {
        assert!(insertion_loss(1.0).is_err());
    }

    #[test]
    fn test_cascade_non_2port() {
        let s3 = SMatrix::new(3, 50.0);
        let s2 = SMatrix::new(2, 50.0);
        assert!(s3.cascade(&s2).is_err());
    }

    #[test]
    fn test_s_to_z_non_2port() {
        let s3 = SMatrix::new(3, 50.0);
        assert!(s_to_z(&s3).is_err());
    }

    #[test]
    fn test_touchstone_parse_mhz() {
        let content = "# MHZ S RI R 75\n100.0  0.1  0.2\n";
        let data = parse_touchstone(content, 1).unwrap();
        assert!((data.z0 - 75.0).abs() < TOL);
        assert!((data.frequencies[0] - 100e6).abs() < 1.0);
    }

    #[test]
    fn test_quarter_wave_negative() {
        assert!(quarter_wave_transformer(-50.0, 100.0).is_err());
    }

    #[test]
    fn test_linear_array_factor_bounded() {
        // AF should always be in [0, 1] for normalized array
        let k = 2.0 * std::f64::consts::PI;
        for i in 0..100 {
            let theta = i as f64 * std::f64::consts::PI / 100.0;
            let af = linear_array_factor(8, 0.5, k, 0.0, theta);
            assert!(
                (-TOL..=1.0 + TOL).contains(&af),
                "AF out of range at θ={theta}: {af}"
            );
        }
    }

    #[test]
    fn test_propagation_constant_lossless() {
        // Lossless (R=0, G=0): γ = jβ = jω√(LC)
        let gamma = propagation_constant(0.0, 250e-9, 0.0, 100e-12, 2e9 * std::f64::consts::PI);
        assert!(gamma.re.abs() < 1e-3); // α ≈ 0 (lossless)
        assert!(gamma.im > 0.0); // β > 0
    }

    #[test]
    fn test_characteristic_impedance_lossy_approaches_lossless() {
        // With R, G very small, lossy Z₀ ≈ lossless Z₀
        let omega = 2e9 * std::f64::consts::PI;
        let z_lossy = characteristic_impedance_lossy(0.001, 250e-9, 1e-6, 100e-12, omega);
        let z_lossless = characteristic_impedance_lossless(250e-9, 100e-12).unwrap();
        assert!((z_lossy.re - z_lossless).abs() / z_lossless < 0.01);
        assert!(z_lossy.im.abs() < 1.0); // small imaginary part
    }

    #[test]
    fn test_feature_gate_rf_compiles_alone() {
        // This test existing proves rf feature compiles
        let _ = reflection_coefficient(75.0, 50.0);
    }
}
