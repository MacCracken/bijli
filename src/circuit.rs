//! Ohm's law, resistance, capacitance, inductance, and basic circuits.
//!
//! All SI units: ohms, farads, henries, volts, amperes, seconds.

use crate::error::{BijliError, Result};

// ── Ohm's law ──────────────────────────────────────────────────────

/// Voltage from current and resistance: V = IR.
#[inline]
#[must_use]
pub fn voltage(current: f64, resistance: f64) -> f64 {
    current * resistance
}

/// Current from voltage and resistance: I = V/R.
#[inline]
pub fn current(voltage: f64, resistance: f64) -> Result<f64> {
    if resistance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "resistance cannot be zero for Ohm's law".into(),
        });
    }
    Ok(voltage / resistance)
}

/// Resistance from voltage and current: R = V/I.
#[inline]
pub fn resistance(voltage: f64, current: f64) -> Result<f64> {
    if current.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "current cannot be zero for resistance calculation".into(),
        });
    }
    Ok(voltage / current)
}

/// Power dissipated: P = IV = I²R = V²/R.
#[inline]
#[must_use]
pub fn power(current: f64, voltage: f64) -> f64 {
    current * voltage
}

/// Resistivity to resistance: R = ρL/A.
#[inline]
pub fn resistance_from_geometry(
    resistivity: f64,
    length: f64,
    cross_section_area: f64,
) -> Result<f64> {
    if cross_section_area <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("cross-section area must be positive, got {cross_section_area}"),
        });
    }
    Ok(resistivity * length / cross_section_area)
}

// ── Series / parallel ──────────────────────────────────────────────

/// Total resistance of resistors in series: R = R₁ + R₂ + ...
///
/// # Errors
///
/// Returns [`BijliError::InvalidParameter`] if the slice is empty.
#[inline]
pub fn resistance_series(resistances: &[f64]) -> Result<f64> {
    if resistances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one resistance value is required".into(),
        });
    }
    Ok(resistances.iter().sum())
}

/// Total resistance of resistors in parallel: 1/R = 1/R₁ + 1/R₂ + ...
///
/// # Errors
///
/// Returns an error if the slice is empty or any resistance is zero.
#[inline]
pub fn resistance_parallel(resistances: &[f64]) -> Result<f64> {
    if resistances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one resistance value is required".into(),
        });
    }
    let mut sum_inv = 0.0;
    for &r in resistances {
        if r.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "resistance cannot be zero in parallel combination".into(),
            });
        }
        sum_inv += 1.0 / r;
    }
    if sum_inv.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "sum of inverse resistances is zero".into(),
        });
    }
    Ok(1.0 / sum_inv)
}

/// Total capacitance of capacitors in series: 1/C = 1/C₁ + 1/C₂ + ...
///
/// # Errors
///
/// Returns an error if the slice is empty or any capacitance is zero.
#[inline]
pub fn capacitance_series(capacitances: &[f64]) -> Result<f64> {
    if capacitances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one capacitance value is required".into(),
        });
    }
    let mut sum_inv = 0.0;
    for &c in capacitances {
        if c.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "capacitance cannot be zero in series combination".into(),
            });
        }
        sum_inv += 1.0 / c;
    }
    if sum_inv.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "sum of inverse capacitances is zero".into(),
        });
    }
    Ok(1.0 / sum_inv)
}

/// Total capacitance of capacitors in parallel: C = C₁ + C₂ + ...
///
/// # Errors
///
/// Returns [`BijliError::InvalidParameter`] if the slice is empty.
#[inline]
pub fn capacitance_parallel(capacitances: &[f64]) -> Result<f64> {
    if capacitances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one capacitance value is required".into(),
        });
    }
    Ok(capacitances.iter().sum())
}

/// Total inductance of inductors in series: L = L₁ + L₂ + ...
///
/// # Errors
///
/// Returns [`BijliError::InvalidParameter`] if the slice is empty.
#[inline]
pub fn inductance_series(inductances: &[f64]) -> Result<f64> {
    if inductances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one inductance value is required".into(),
        });
    }
    Ok(inductances.iter().sum())
}

/// Total inductance of inductors in parallel: 1/L = 1/L₁ + 1/L₂ + ...
///
/// # Errors
///
/// Returns an error if the slice is empty or any inductance is zero.
#[inline]
pub fn inductance_parallel(inductances: &[f64]) -> Result<f64> {
    if inductances.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "at least one inductance value is required".into(),
        });
    }
    let mut sum_inv = 0.0;
    for &l in inductances {
        if l.abs() < 1e-30 {
            return Err(BijliError::DivisionByZero {
                context: "inductance cannot be zero in parallel combination".into(),
            });
        }
        sum_inv += 1.0 / l;
    }
    if sum_inv.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "sum of inverse inductances is zero".into(),
        });
    }
    Ok(1.0 / sum_inv)
}

// ── Capacitance & Inductance ───────────────────────────────────────

/// Parallel plate capacitance: C = εA/d.
#[inline]
pub fn parallel_plate_capacitance(permittivity: f64, area: f64, separation: f64) -> Result<f64> {
    if separation <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("plate separation must be positive, got {separation}"),
        });
    }
    Ok(permittivity * area / separation)
}

/// Energy stored in a capacitor: U = CV²/2.
#[inline]
#[must_use]
pub fn capacitor_energy(capacitance: f64, voltage: f64) -> f64 {
    0.5 * capacitance * voltage * voltage
}

/// Charge on a capacitor: Q = CV.
#[inline]
#[must_use]
pub fn capacitor_charge(capacitance: f64, voltage: f64) -> f64 {
    capacitance * voltage
}

/// Solenoid inductance: L = μ₀n²Al.
///
/// `n` is turns per unit length, `area` is cross-section, `length` is solenoid length.
#[inline]
#[must_use]
pub fn solenoid_inductance(
    permeability: f64,
    turns_per_length: f64,
    area: f64,
    length: f64,
) -> f64 {
    permeability * turns_per_length * turns_per_length * area * length
}

/// Energy stored in an inductor: U = LI²/2.
#[inline]
#[must_use]
pub fn inductor_energy(inductance: f64, current: f64) -> f64 {
    0.5 * inductance * current * current
}

// ── RC circuit ─────────────────────────────────────────────────────

/// RC time constant: τ = RC.
#[inline]
#[must_use]
pub fn rc_time_constant(resistance: f64, capacitance: f64) -> f64 {
    resistance * capacitance
}

/// RC charging voltage: V(t) = V₀(1 − e^(−t/RC)).
#[inline]
#[must_use]
pub fn rc_charging_voltage(v0: f64, resistance: f64, capacitance: f64, t: f64) -> f64 {
    let tau = resistance * capacitance;
    if tau.abs() < 1e-30 {
        return v0; // instant charge
    }
    v0 * (1.0 - (-t / tau).exp())
}

/// RC discharging voltage: V(t) = V₀ e^(−t/RC).
#[inline]
#[must_use]
pub fn rc_discharging_voltage(v0: f64, resistance: f64, capacitance: f64, t: f64) -> f64 {
    let tau = resistance * capacitance;
    if tau.abs() < 1e-30 {
        return 0.0; // instant discharge
    }
    v0 * (-t / tau).exp()
}

// ── RL circuit ─────────────────────────────────────────────────────

/// RL time constant: τ = L/R.
#[inline]
pub fn rl_time_constant(inductance: f64, resistance: f64) -> Result<f64> {
    if resistance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "resistance cannot be zero for RL time constant".into(),
        });
    }
    Ok(inductance / resistance)
}

/// RL circuit current rise: I(t) = (V/R)(1 − e^(−tR/L)).
#[inline]
pub fn rl_current_rise(voltage: f64, resistance: f64, inductance: f64, t: f64) -> Result<f64> {
    if resistance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "resistance cannot be zero for RL current".into(),
        });
    }
    let tau = inductance / resistance;
    if tau.abs() < 1e-30 {
        return Ok(voltage / resistance); // instant rise
    }
    Ok((voltage / resistance) * (1.0 - (-t / tau).exp()))
}

/// RL circuit current decay: I(t) = I₀ e^(−tR/L).
#[inline]
pub fn rl_current_decay(i0: f64, resistance: f64, inductance: f64, t: f64) -> Result<f64> {
    if resistance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "resistance cannot be zero for RL decay".into(),
        });
    }
    let tau = inductance / resistance;
    if tau.abs() < 1e-30 {
        return Ok(0.0);
    }
    Ok(i0 * (-t / tau).exp())
}

// ── RLC circuit ────────────────────────────────────────────────────

/// Resonant frequency of an LC/RLC circuit: f₀ = 1/(2π√(LC)).
#[inline]
pub fn resonant_frequency(inductance: f64, capacitance: f64) -> Result<f64> {
    let product = inductance * capacitance;
    if product <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("LC product must be positive, got L={inductance}, C={capacitance}"),
        });
    }
    Ok(1.0 / (2.0 * std::f64::consts::PI * product.sqrt()))
}

/// Angular resonant frequency: ω₀ = 1/√(LC).
#[inline]
pub fn resonant_angular_frequency(inductance: f64, capacitance: f64) -> Result<f64> {
    let product = inductance * capacitance;
    if product <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("LC product must be positive, got L={inductance}, C={capacitance}"),
        });
    }
    Ok(1.0 / product.sqrt())
}

/// Quality factor of an RLC circuit: Q = (1/R)√(L/C).
#[inline]
pub fn quality_factor(resistance: f64, inductance: f64, capacitance: f64) -> Result<f64> {
    if resistance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "resistance cannot be zero for quality factor".into(),
        });
    }
    if capacitance <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("capacitance must be positive, got {capacitance}"),
        });
    }
    Ok((inductance / capacitance).sqrt() / resistance)
}

/// Damping ratio of a series RLC circuit: ζ = R/(2)√(C/L).
#[inline]
pub fn damping_ratio(resistance: f64, inductance: f64, capacitance: f64) -> Result<f64> {
    if inductance <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("inductance must be positive, got {inductance}"),
        });
    }
    Ok(resistance * 0.5 * (capacitance / inductance).sqrt())
}

/// Impedance of an RLC series circuit at angular frequency ω.
///
/// Z = √(R² + (ωL − 1/(ωC))²)
#[inline]
pub fn rlc_impedance(
    resistance: f64,
    inductance: f64,
    capacitance: f64,
    angular_freq: f64,
) -> Result<f64> {
    if angular_freq.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "angular frequency cannot be zero for impedance".into(),
        });
    }
    if capacitance.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "capacitance cannot be zero for impedance".into(),
        });
    }
    let x = angular_freq * inductance - 1.0 / (angular_freq * capacitance);
    Ok((resistance * resistance + x * x).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Ohm's law ──────────────────────────────────────────────────

    #[test]
    fn test_voltage() {
        assert!((voltage(2.0, 3.0) - 6.0).abs() < 1e-15);
    }

    #[test]
    fn test_current() {
        assert!((current(6.0, 3.0).unwrap() - 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_current_zero_resistance() {
        assert!(current(6.0, 0.0).is_err());
    }

    #[test]
    fn test_resistance() {
        assert!((resistance(6.0, 2.0).unwrap() - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_power() {
        assert!((power(2.0, 6.0) - 12.0).abs() < 1e-15);
    }

    #[test]
    fn test_resistance_from_geometry() {
        // Copper wire: ρ ≈ 1.68e-8, L = 1m, A = 1mm² = 1e-6 m²
        let r = resistance_from_geometry(1.68e-8, 1.0, 1e-6).unwrap();
        assert!((r - 0.0168).abs() < 1e-4);
    }

    // ── Series / parallel ──────────────────────────────────────────

    #[test]
    fn test_resistance_series() {
        assert!((resistance_series(&[10.0, 20.0, 30.0]).unwrap() - 60.0).abs() < 1e-15);
    }

    #[test]
    fn test_resistance_parallel() {
        // Two 10Ω in parallel = 5Ω
        let r = resistance_parallel(&[10.0, 10.0]).unwrap();
        assert!((r - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_resistance_parallel_zero() {
        assert!(resistance_parallel(&[10.0, 0.0]).is_err());
    }

    #[test]
    fn test_capacitance_series() {
        // Two 10μF in series = 5μF
        let c = capacitance_series(&[10e-6, 10e-6]).unwrap();
        assert!((c - 5e-6).abs() / 5e-6 < 1e-10);
    }

    #[test]
    fn test_capacitance_parallel() {
        assert!((capacitance_parallel(&[10e-6, 20e-6]).unwrap() - 30e-6).abs() < 1e-20);
    }

    #[test]
    fn test_inductance_series() {
        assert!((inductance_series(&[1e-3, 2e-3]).unwrap() - 3e-3).abs() < 1e-18);
    }

    #[test]
    fn test_inductance_parallel() {
        let l = inductance_parallel(&[2e-3, 2e-3]).unwrap();
        assert!((l - 1e-3).abs() / 1e-3 < 1e-10);
    }

    // ── Capacitance & Inductance ───────────────────────────────────

    #[test]
    fn test_parallel_plate_capacitance() {
        use crate::field::EPSILON_0;
        // 1m² plates, 1mm apart, vacuum
        let c = parallel_plate_capacitance(EPSILON_0, 1.0, 1e-3).unwrap();
        assert!((c - EPSILON_0 * 1000.0).abs() / c < 1e-6);
    }

    #[test]
    fn test_capacitor_energy() {
        // 1F at 1V = 0.5J
        assert!((capacitor_energy(1.0, 1.0) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn test_capacitor_charge() {
        assert!((capacitor_charge(1e-6, 5.0) - 5e-6).abs() < 1e-20);
    }

    #[test]
    fn test_solenoid_inductance() {
        use crate::field::MU_0;
        // 1000 turns/m, 1cm² cross-section, 10cm long
        let l = solenoid_inductance(MU_0, 1000.0, 1e-4, 0.1);
        // L = μ₀ × 10⁶ × 10⁻⁴ × 0.1 ≈ 1.257e-5 H
        assert!((l - 1.257e-5).abs() / 1.257e-5 < 0.01);
    }

    #[test]
    fn test_inductor_energy() {
        // 1H at 1A = 0.5J
        assert!((inductor_energy(1.0, 1.0) - 0.5).abs() < 1e-15);
    }

    // ── RC circuit ─────────────────────────────────────────────────

    #[test]
    fn test_rc_time_constant() {
        assert!((rc_time_constant(1e3, 1e-6) - 1e-3).abs() < 1e-18);
    }

    #[test]
    fn test_rc_charging_at_zero() {
        let v = rc_charging_voltage(10.0, 1e3, 1e-6, 0.0);
        assert!(v.abs() < 1e-10);
    }

    #[test]
    fn test_rc_charging_at_one_tau() {
        let tau = 1e-3;
        let v = rc_charging_voltage(10.0, 1e3, 1e-6, tau);
        // V(τ) = V₀(1 - 1/e) ≈ 0.6321 × V₀
        assert!((v / 10.0 - 0.6321).abs() < 0.001);
    }

    #[test]
    fn test_rc_discharging_at_one_tau() {
        let tau = 1e-3;
        let v = rc_discharging_voltage(10.0, 1e3, 1e-6, tau);
        // V(τ) = V₀/e ≈ 3.679
        assert!((v - 10.0 / std::f64::consts::E).abs() < 0.01);
    }

    // ── RL circuit ─────────────────────────────────────────────────

    #[test]
    fn test_rl_time_constant() {
        let tau = rl_time_constant(1e-3, 1.0).unwrap();
        assert!((tau - 1e-3).abs() < 1e-18);
    }

    #[test]
    fn test_rl_current_rise_at_one_tau() {
        let v = 10.0;
        let r = 1.0;
        let l = 1e-3;
        let tau = l / r;
        let i = rl_current_rise(v, r, l, tau).unwrap();
        // I(τ) = (V/R)(1 - 1/e) ≈ 6.321 A
        assert!((i / (v / r) - 0.6321).abs() < 0.001);
    }

    #[test]
    fn test_rl_current_decay() {
        let i = rl_current_decay(10.0, 1.0, 1e-3, 1e-3).unwrap();
        assert!((i - 10.0 / std::f64::consts::E).abs() < 0.01);
    }

    // ── RLC circuit ────────────────────────────────────────────────

    #[test]
    fn test_resonant_frequency() {
        // L = 1mH, C = 1μF → f₀ ≈ 5033 Hz
        let f = resonant_frequency(1e-3, 1e-6).unwrap();
        assert!((f - 5032.9).abs() < 1.0);
    }

    #[test]
    fn test_resonant_angular_frequency() {
        let omega = resonant_angular_frequency(1e-3, 1e-6).unwrap();
        let f = resonant_frequency(1e-3, 1e-6).unwrap();
        assert!((omega - 2.0 * std::f64::consts::PI * f).abs() < 1.0);
    }

    #[test]
    fn test_quality_factor() {
        // R=10Ω, L=1mH, C=1μF → Q = √(L/C)/R ≈ 3.162
        let q = quality_factor(10.0, 1e-3, 1e-6).unwrap();
        assert!((q - 3.162).abs() < 0.01);
    }

    #[test]
    fn test_damping_ratio() {
        let zeta = damping_ratio(10.0, 1e-3, 1e-6).unwrap();
        // ζ = R/2 × √(C/L) = 5 × √(1e-3) ≈ 0.158
        assert!((zeta - 0.158).abs() < 0.01);
    }

    #[test]
    fn test_rlc_impedance_at_resonance() {
        // At resonance, ωL = 1/(ωC), so Z = R
        let omega = resonant_angular_frequency(1e-3, 1e-6).unwrap();
        let z = rlc_impedance(10.0, 1e-3, 1e-6, omega).unwrap();
        assert!((z - 10.0).abs() < 0.01);
    }

    #[test]
    fn test_rlc_impedance_zero_freq() {
        assert!(rlc_impedance(10.0, 1e-3, 1e-6, 0.0).is_err());
    }

    // ── Empty slice edge cases ────────────────────────────────────

    #[test]
    fn test_resistance_series_empty() {
        assert!(resistance_series(&[]).is_err());
    }

    #[test]
    fn test_resistance_parallel_empty() {
        assert!(resistance_parallel(&[]).is_err());
    }

    #[test]
    fn test_capacitance_series_empty() {
        assert!(capacitance_series(&[]).is_err());
    }

    #[test]
    fn test_capacitance_parallel_empty() {
        assert!(capacitance_parallel(&[]).is_err());
    }

    #[test]
    fn test_inductance_series_empty() {
        assert!(inductance_series(&[]).is_err());
    }

    #[test]
    fn test_inductance_parallel_empty() {
        assert!(inductance_parallel(&[]).is_err());
    }

    // ── Single-element combination ────────────────────────────────

    #[test]
    fn test_resistance_series_single() {
        assert!((resistance_series(&[47.0]).unwrap() - 47.0).abs() < 1e-15);
    }

    #[test]
    fn test_resistance_parallel_single() {
        assert!((resistance_parallel(&[47.0]).unwrap() - 47.0).abs() < 1e-10);
    }

    #[test]
    fn test_capacitance_parallel_single() {
        assert!((capacitance_parallel(&[10e-6]).unwrap() - 10e-6).abs() < 1e-20);
    }

    // ── RC/RL zero-tau edge cases ─────────────────────────────────

    #[test]
    fn test_rc_charging_zero_tau() {
        // tau ≈ 0 → instant charge
        let v = rc_charging_voltage(10.0, 0.0, 1e-6, 1.0);
        assert!((v - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_rc_discharging_zero_tau() {
        // tau ≈ 0 → instant discharge
        let v = rc_discharging_voltage(10.0, 0.0, 1e-6, 1.0);
        assert!(v.abs() < 1e-10);
    }

    // ── Geometry edge cases ───────────────────────────────────────

    #[test]
    fn test_resistance_from_geometry_zero_area() {
        assert!(resistance_from_geometry(1.68e-8, 1.0, 0.0).is_err());
    }

    #[test]
    fn test_parallel_plate_capacitance_zero_separation() {
        use crate::field::EPSILON_0;
        assert!(parallel_plate_capacitance(EPSILON_0, 1.0, 0.0).is_err());
    }

    #[test]
    fn test_quality_factor_zero_resistance() {
        assert!(quality_factor(0.0, 1e-3, 1e-6).is_err());
    }

    #[test]
    fn test_damping_ratio_zero_inductance() {
        assert!(damping_ratio(10.0, 0.0, 1e-6).is_err());
    }

    #[test]
    fn test_resonant_frequency_negative_lc() {
        assert!(resonant_frequency(-1e-3, 1e-6).is_err());
    }
}
