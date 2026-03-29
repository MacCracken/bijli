//! RF & Microwave — transmission lines, S-parameters, Smith chart,
//! matching networks, Touchstone I/O, antenna arrays.
//!
//! All SI units unless noted. Impedances in ohms, frequencies in Hz.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};
use crate::polarization::{Complex, ComplexExt, complex_from_polar, complex_real, complex_zero};

// ── Transmission lines ────────────────────────────────────────────

/// Characteristic impedance of a lossless transmission line: Z₀ = √(L/C).
#[inline]
pub fn characteristic_impedance_lossless(
    inductance_per_m: f64,
    capacitance_per_m: f64,
) -> Result<f64> {
    if inductance_per_m <= 0.0 || capacitance_per_m <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "L' and C' must be positive, got L'={inductance_per_m}, C'={capacitance_per_m}"
            ),
        });
    }
    Ok((inductance_per_m / capacitance_per_m).sqrt())
}

/// Complex characteristic impedance of a lossy line: Z₀ = √((R'+jωL')/(G'+jωC')).
#[inline]
#[must_use]
pub fn characteristic_impedance_lossy(r: f64, l: f64, g: f64, c: f64, omega: f64) -> Complex {
    let num = Complex::new(r, omega * l);
    let den = Complex::new(g, omega * c);
    // sqrt(num/den)
    let ratio = num / den;
    let mag = ratio.norm().sqrt();
    let phase = ratio.arg() / 2.0;
    complex_from_polar(mag, phase)
}

/// Complex propagation constant: γ = α + jβ = √((R'+jωL')(G'+jωC')).
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

/// Phase velocity on a lossless line: v_p = 1/√(LC).
#[inline]
pub fn phase_velocity(inductance_per_m: f64, capacitance_per_m: f64) -> Result<f64> {
    if inductance_per_m <= 0.0 || capacitance_per_m <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: "L' and C' must be positive".into(),
        });
    }
    Ok(1.0 / (inductance_per_m * capacitance_per_m).sqrt())
}

/// Input impedance of a transmission line terminated with Z_L.
///
/// Z_in = Z₀ (Z_L + Z₀ tanh(γl)) / (Z₀ + Z_L tanh(γl))
///
/// For lossless line: γ = jβ, tanh(jβl) = j tan(βl).
#[inline]
pub fn input_impedance_lossless(z0: f64, z_load: f64, beta: f64, length: f64) -> Result<f64> {
    let tan_bl = (beta * length).tan();
    let den = z0 + z_load * tan_bl;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "transmission line input impedance denominator is zero".into(),
        });
    }
    Ok(z0 * (z_load + z0 * tan_bl) / den)
}

// ── Reflection coefficient & VSWR ─────────────────────────────────

/// Voltage reflection coefficient: Γ = (Z_L − Z₀) / (Z_L + Z₀).
#[inline]
pub fn reflection_coefficient(z_load: f64, z0: f64) -> Result<f64> {
    let den = z_load + z0;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Z_L + Z₀ = 0".into(),
        });
    }
    Ok((z_load - z0) / den)
}

/// Complex reflection coefficient for complex impedances.
#[inline]
pub fn reflection_coefficient_complex(z_load: Complex, z0: Complex) -> Result<Complex> {
    let den = z_load + z0;
    if den.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "Z_L + Z₀ = 0".into(),
        });
    }
    Ok((z_load - z0) / den)
}

/// VSWR from reflection coefficient magnitude: VSWR = (1 + |Γ|) / (1 − |Γ|).
#[inline]
pub fn vswr(gamma_mag: f64) -> Result<f64> {
    if gamma_mag >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be < 1, got {gamma_mag}"),
        });
    }
    if gamma_mag < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be ≥ 0, got {gamma_mag}"),
        });
    }
    Ok((1.0 + gamma_mag) / (1.0 - gamma_mag))
}

/// Return loss in dB: RL = −20 log₁₀|Γ|.
#[inline]
pub fn return_loss(gamma_mag: f64) -> Result<f64> {
    if gamma_mag <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|Γ| must be positive for return loss, got {gamma_mag}"),
        });
    }
    Ok(-20.0 * gamma_mag.log10())
}

/// Insertion loss in dB: IL = −10 log₁₀(1 − |Γ|²).
#[inline]
pub fn insertion_loss(gamma_mag: f64) -> Result<f64> {
    let g2 = gamma_mag * gamma_mag;
    if g2 >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: "total reflection — infinite insertion loss".into(),
        });
    }
    Ok(-10.0 * (1.0 - g2).log10())
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
#[inline]
#[must_use]
pub fn gamma_from_z(z_norm: Complex) -> Complex {
    let one = complex_real(1.0);
    (z_norm - one) / (z_norm + one)
}

/// Normalized impedance from reflection coefficient: z = (1 + Γ) / (1 − Γ).
#[inline]
#[must_use]
pub fn z_from_gamma(gamma: Complex) -> Complex {
    let one = complex_real(1.0);
    (one + gamma) / (one - gamma)
}

// ── Matching networks ─────────────────────────────────────────────

/// Quarter-wave transformer impedance: Z_t = √(Z₀ × Z_L).
#[inline]
pub fn quarter_wave_transformer(z0: f64, z_load: f64) -> Result<f64> {
    let product = z0 * z_load;
    if product < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: "Z₀ × Z_L must be non-negative for quarter-wave transformer".into(),
        });
    }
    Ok(product.sqrt())
}

/// Single-stub matching — required stub length and position.
///
/// Returns (distance_from_load, stub_length) in units of wavelength.
/// Uses shunt short-circuited stub.
#[inline]
pub fn single_stub_match(z_load: f64, z0: f64) -> Result<(f64, f64)> {
    if z0.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Z₀ cannot be zero".into(),
        });
    }
    let r = z_load / z0;
    if (r - 1.0).abs() < 1e-10 {
        // Already matched
        return Ok((0.0, 0.0));
    }

    // Distance from load (in wavelengths)
    let d = if r > 1.0 {
        (1.0 / (r.sqrt())).atan() / (2.0 * std::f64::consts::PI)
    } else {
        (std::f64::consts::PI - (1.0 / (r.sqrt())).atan()) / (2.0 * std::f64::consts::PI)
    };

    // Stub length (short-circuit stub, in wavelengths)
    let b = ((r - 1.0).abs()).sqrt() / z0;
    let stub_len = (1.0 / b).atan() / (2.0 * std::f64::consts::PI);

    Ok((d, stub_len.abs()))
}

// ── S-parameters ──────────────────────────────────────────────────

/// S-parameter matrix for an N-port network.
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

    /// Get S_{ij} (1-indexed ports).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Result<Complex> {
        if i == 0 || j == 0 || i > self.ports || j > self.ports {
            return Err(BijliError::InvalidParameter {
                reason: format!("port indices must be 1..={}, got ({i}, {j})", self.ports),
            });
        }
        Ok(self.s[i - 1][j - 1])
    }

    /// Set S_{ij} (1-indexed ports).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: Complex) -> Result<()> {
        if i == 0 || j == 0 || i > self.ports || j > self.ports {
            return Err(BijliError::InvalidParameter {
                reason: format!("port indices must be 1..={}, got ({i}, {j})", self.ports),
            });
        }
        self.s[i - 1][j - 1] = value;
        Ok(())
    }

    /// Cascade two 2-port networks using T-parameter conversion.
    ///
    /// Returns the combined S-matrix of self followed by other.
    pub fn cascade(&self, other: &Self) -> Result<Self> {
        if self.ports != 2 || other.ports != 2 {
            return Err(BijliError::InvalidParameter {
                reason: "cascade requires two 2-port networks".into(),
            });
        }

        // Convert to T-parameters, multiply, convert back
        let t1 = s_to_t(self)?;
        let t2 = s_to_t(other)?;

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

        t_to_s(&t, self.z0[0])
    }
}

/// Convert 2-port S-parameters to T (transfer) parameters.
fn s_to_t(s: &SMatrix) -> Result<[[Complex; 2]; 2]> {
    let s21 = s.s[1][0];
    if s21.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "S21 = 0, cannot convert to T-parameters".into(),
        });
    }
    let s11 = s.s[0][0];
    let s12 = s.s[0][1];
    let s22 = s.s[1][1];
    let inv_s21 = complex_real(1.0) / s21;

    let det = s11 * s22 - s12 * s21;
    Ok([
        [inv_s21 * complex_real(-1.0) * det, s11 * inv_s21],
        [complex_real(-1.0) * s22 * inv_s21, inv_s21],
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
    let inv_t22 = complex_real(1.0) / t22;
    let s11 = t[0][1] * inv_t22;
    let s12 = (t[0][0] * t[1][1] - t[0][1] * t[1][0]) * inv_t22;
    let s21 = inv_t22;
    let s22 = complex_real(-1.0) * t[1][0] * inv_t22;
    Ok(SMatrix::two_port(s11, s12, s21, s22, z0))
}

/// Convert 2-port S-parameters to Z-parameters.
pub fn s_to_z(s: &SMatrix) -> Result<[[Complex; 2]; 2]> {
    if s.ports != 2 {
        return Err(BijliError::InvalidParameter {
            reason: "S→Z conversion requires 2-port network".into(),
        });
    }
    let one = complex_real(1.0);
    let s11 = s.s[0][0];
    let s12 = s.s[0][1];
    let s21 = s.s[1][0];
    let s22 = s.s[1][1];
    let z0 = complex_real(s.z0[0]);

    let det = (one - s11) * (one - s22) - s12 * s21;
    if det.norm_sq() < 1e-60 {
        return Err(BijliError::DivisionByZero {
            context: "singular S→Z conversion".into(),
        });
    }
    let inv_det = one / det;

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
    /// dB magnitude and angle (degrees).
    DbAngle,
}

/// Parsed Touchstone data for an N-port network.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TouchstoneData {
    /// Number of ports.
    pub ports: usize,
    /// Reference impedance (ohms).
    pub z0: f64,
    /// Frequency points (Hz).
    pub frequencies: Vec<f64>,
    /// S-parameter matrices at each frequency.
    pub s_matrices: Vec<SMatrix>,
}

/// Parse Touchstone (.snp) file content.
///
/// Supports S-parameter data in RI (real/imaginary), MA (magnitude/angle),
/// and DB (dB/angle) formats.
pub fn parse_touchstone(content: &str, ports: usize) -> Result<TouchstoneData> {
    let mut frequencies = Vec::new();
    let mut s_matrices = Vec::new();
    let mut z0 = 50.0;
    let mut format = TouchstoneFormat::MagnitudeAngle;
    let mut freq_mult = 1e9; // default GHz

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('!') {
            continue;
        }
        if let Some(opts) = line.strip_prefix('#') {
            // Option line: # GHz S MA R 50
            let parts: Vec<&str> = opts.split_whitespace().collect();
            for (i, &p) in parts.iter().enumerate() {
                match p.to_uppercase().as_str() {
                    "HZ" => freq_mult = 1.0,
                    "KHZ" => freq_mult = 1e3,
                    "MHZ" => freq_mult = 1e6,
                    "GHZ" => freq_mult = 1e9,
                    "MA" => format = TouchstoneFormat::MagnitudeAngle,
                    "RI" => format = TouchstoneFormat::RealImag,
                    "DB" => format = TouchstoneFormat::DbAngle,
                    "R" => {
                        if let Some(Ok(v)) = parts.get(i + 1).map(|s| s.parse::<f64>()) {
                            z0 = v;
                        }
                    }
                    _ => {}
                }
            }
            continue;
        }

        // Data line
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        if values.is_empty() {
            continue;
        }

        let freq = values[0] * freq_mult;
        let n = ports;
        let expected = 1 + 2 * n * n; // freq + n*n complex values
        if values.len() < expected {
            continue; // skip incomplete lines
        }

        let mut sm = SMatrix::new(n, z0);
        let mut idx = 1;
        for i in 0..n {
            for j in 0..n {
                let c = parse_complex_pair(values[idx], values[idx + 1], format);
                sm.s[i][j] = c;
                idx += 2;
            }
        }

        frequencies.push(freq);
        s_matrices.push(sm);
    }

    if frequencies.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "no valid data found in Touchstone content".into(),
        });
    }

    Ok(TouchstoneData {
        ports,
        z0,
        frequencies,
        s_matrices,
    })
}

/// Format Touchstone output string.
#[must_use]
pub fn write_touchstone(data: &TouchstoneData, format: TouchstoneFormat) -> String {
    let mut out = String::new();
    out.push_str("! Touchstone file generated by bijli\n");

    let fmt_str = match format {
        TouchstoneFormat::MagnitudeAngle => "MA",
        TouchstoneFormat::RealImag => "RI",
        TouchstoneFormat::DbAngle => "DB",
    };
    out.push_str(&format!("# GHz S {fmt_str} R {}\n", data.z0));

    for (fi, freq) in data.frequencies.iter().enumerate() {
        let freq_ghz = freq / 1e9;
        out.push_str(&format!("{freq_ghz:.6}"));

        let sm = &data.s_matrices[fi];
        for i in 0..data.ports {
            for j in 0..data.ports {
                let c = sm.s[i][j];
                let (v1, v2) = format_complex_pair(c, format);
                out.push_str(&format!("  {v1:.6}  {v2:.6}"));
            }
        }
        out.push('\n');
    }

    out
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
        TouchstoneFormat::DbAngle => (20.0 * c.norm().log10(), c.arg().to_degrees()),
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
    }

    #[test]
    fn test_phase_velocity() {
        let vp = phase_velocity(250e-9, 100e-12).unwrap();
        // v = 1/√(LC) ≈ 2e8 m/s for coax with polyethylene dielectric
        assert!(vp > 1e8);
        assert!(vp < 3e8);
    }

    #[test]
    fn test_input_impedance_matched() {
        // Matched load: Z_in = Z₀ regardless of length
        let z_in = input_impedance_lossless(50.0, 50.0, 10.0, 0.1).unwrap();
        assert!((z_in - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_input_impedance_quarter_wave() {
        // λ/4 transformer: Z_in = Z₀²/Z_L
        let beta = 2.0 * std::f64::consts::PI; // β = 2π/λ, so βl = π/2 at l = λ/4
        let z_in = input_impedance_lossless(50.0, 100.0, beta, 0.25).unwrap();
        // Z_in = 50² / 100 = 25
        assert!((z_in - 25.0).abs() < 0.1);
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
    }

    #[test]
    fn test_return_loss() {
        // |Γ| = 0.1 → RL = 20 dB
        let rl = return_loss(0.1).unwrap();
        assert!((rl - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_insertion_loss_matched() {
        let il = insertion_loss(0.0).unwrap();
        assert!(il.abs() < TOL); // 0 dB for perfect match
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

    // ── Matching tests ────────────────────────────────────────────

    #[test]
    fn test_quarter_wave_transformer() {
        // Match 50Ω to 100Ω: Z_t = √(5000) ≈ 70.71Ω
        let zt = quarter_wave_transformer(50.0, 100.0).unwrap();
        assert!((zt - 70.71).abs() < 0.01);
    }

    #[test]
    fn test_single_stub_already_matched() {
        let (d, l) = single_stub_match(50.0, 50.0).unwrap();
        assert!(d.abs() < TOL);
        assert!(l.abs() < TOL);
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
        // For S11=S22=0: Z11 = Z22 = Z₀(1+S12*S21)/(1-S12*S21) ... complex
        // Just verify it produces finite values and Z12 = Z21 (reciprocal)
        assert!(z[0][0].norm() > 0.0);
        assert!((z[0][1].re - z[1][0].re).abs() < TOL);
        assert!((z[0][1].im - z[1][0].im).abs() < TOL);
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

    #[test]
    fn test_touchstone_roundtrip() {
        let mut data = TouchstoneData {
            ports: 1,
            z0: 50.0,
            frequencies: vec![1e9, 2e9],
            s_matrices: vec![
                SMatrix::two_port(
                    Complex::new(0.5, 0.1),
                    complex_zero(),
                    complex_zero(),
                    complex_zero(),
                    50.0,
                ),
                SMatrix::two_port(
                    Complex::new(0.3, -0.2),
                    complex_zero(),
                    complex_zero(),
                    complex_zero(),
                    50.0,
                ),
            ],
        };
        data.ports = 1;
        // Write as RI format
        let content = write_touchstone(&data, TouchstoneFormat::RealImag);
        assert!(content.contains("# GHz S RI R 50"));
        assert!(content.contains("1.000000"));
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
