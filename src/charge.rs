//! Point charges, Coulomb's law, Lorentz force, charge distributions.
//!
//! All SI units: coulombs, newtons, meters, seconds.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};
use crate::field::{COULOMB_K, FieldVector};

/// Elementary charge e (coulombs).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Electron mass (kg).
pub const ELECTRON_MASS: f64 = 9.109_383_702e-31;

/// Proton mass (kg).
pub const PROTON_MASS: f64 = 1.672_621_924e-27;

/// A point charge with position and velocity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointCharge {
    /// Charge in coulombs.
    pub charge: f64,
    /// Mass in kilograms.
    pub mass: f64,
    /// Position [x, y, z] in meters.
    pub position: [f64; 3],
    /// Velocity [vx, vy, vz] in m/s.
    pub velocity: [f64; 3],
}

impl PointCharge {
    /// Create a new point charge.
    ///
    /// # Errors
    ///
    /// Returns [`BijliError::InvalidParameter`] if `mass` is not positive (or NaN).
    #[inline]
    pub fn new(charge: f64, mass: f64, position: [f64; 3], velocity: [f64; 3]) -> Result<Self> {
        if mass.is_nan() || mass <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("mass must be positive, got {mass}"),
            });
        }
        Ok(Self {
            charge,
            mass,
            position,
            velocity,
        })
    }

    /// Create an electron at rest at a position.
    #[inline]
    #[must_use]
    pub fn electron(position: [f64; 3]) -> Self {
        // SAFETY: ELECTRON_MASS is a known positive constant.
        Self {
            charge: -ELEMENTARY_CHARGE,
            mass: ELECTRON_MASS,
            position,
            velocity: [0.0; 3],
        }
    }

    /// Create a proton at rest at a position.
    #[inline]
    #[must_use]
    pub fn proton(position: [f64; 3]) -> Self {
        // SAFETY: PROTON_MASS is a known positive constant.
        Self {
            charge: ELEMENTARY_CHARGE,
            mass: PROTON_MASS,
            position,
            velocity: [0.0; 3],
        }
    }

    /// Electric field at a point due to this charge.
    ///
    /// Delegates to [`crate::field::electric_field_point_charge`].
    #[inline]
    pub fn electric_field_at(&self, field_pos: [f64; 3]) -> Result<FieldVector> {
        crate::field::electric_field_point_charge(self.charge, self.position, field_pos)
    }

    /// Electric potential at a point due to this charge.
    ///
    /// Delegates to [`crate::field::electric_potential_point_charge`].
    #[inline]
    pub fn electric_potential_at(&self, field_pos: [f64; 3]) -> Result<f64> {
        crate::field::electric_potential_point_charge(self.charge, self.position, field_pos)
    }

    /// Magnetic field at a point due to this moving charge (Biot-Savart).
    ///
    /// Delegates to [`crate::field::magnetic_field_moving_charge`].
    #[inline]
    pub fn magnetic_field_at(&self, field_pos: [f64; 3]) -> Result<FieldVector> {
        crate::field::magnetic_field_moving_charge(
            self.charge,
            self.position,
            self.velocity,
            field_pos,
        )
    }

    /// Lorentz force on this charge in given E and B fields.
    #[inline]
    #[must_use]
    pub fn lorentz_force_in(
        &self,
        electric_field: &FieldVector,
        magnetic_field: &FieldVector,
    ) -> FieldVector {
        let v = FieldVector::new(self.velocity[0], self.velocity[1], self.velocity[2]);
        lorentz_force(self.charge, &v, electric_field, magnetic_field)
    }
}

/// Coulomb force between two point charges.
///
/// F = kq₁q₂/r² r̂ (newtons, directed from charge2 toward charge1 if repulsive)
#[inline]
pub fn coulomb_force(q1: &PointCharge, q2: &PointCharge) -> Result<FieldVector> {
    let dx = q1.position[0] - q2.position[0];
    let dy = q1.position[1] - q2.position[1];
    let dz = q1.position[2] - q2.position[2];
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    // kq₁q₂ / r³ — avoids separate sqrt + division
    let inv_r = 1.0 / r_sq.sqrt();
    let factor = COULOMB_K * q1.charge * q2.charge * inv_r * inv_r * inv_r;

    Ok(FieldVector {
        x: factor * dx,
        y: factor * dy,
        z: factor * dz,
    })
}

/// Lorentz force on a charged particle: F = q(E + v × B).
#[inline]
#[must_use]
pub fn lorentz_force(
    charge: f64,
    velocity: &FieldVector,
    electric_field: &FieldVector,
    magnetic_field: &FieldVector,
) -> FieldVector {
    (*electric_field + velocity.cross(magnetic_field)) * charge
}

/// Electric dipole moment: p = qd (C⋅m).
#[inline]
#[must_use]
pub fn dipole_moment(charge: f64, separation: f64) -> f64 {
    charge * separation
}

/// Potential energy of a charge in an electric potential.
///
/// U = qV (joules)
#[inline]
#[must_use]
pub fn potential_energy(charge: f64, potential: f64) -> f64 {
    charge * potential
}

/// Coulomb potential energy between two point charges.
///
/// U = kq₁q₂/r (joules), `distance` r ≥ 0 in metres.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `distance` is negative or NaN,
/// [`BijliError::Singularity`] if it is (numerically) zero.
#[inline]
pub fn coulomb_potential_energy(q1: f64, q2: f64, distance: f64) -> Result<f64> {
    crate::field::check_radial_distance(distance)?;
    Ok(COULOMB_K * q1 * q2 / distance)
}

/// Validate a particle mass: zero → `DivisionByZero`, negative/NaN → `InvalidParameter`.
#[inline]
fn check_mass(mass: f64, what: &str) -> Result<()> {
    if mass == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: format!("mass cannot be zero for {what}"),
        });
    }
    if mass.is_nan() || mass <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("mass must be positive for {what}, got {mass}"),
        });
    }
    Ok(())
}

/// Cyclotron frequency for a charged particle in a magnetic field.
///
/// ω_c = |q|B/m (rad/s), with `b_magnitude` = |B| ≥ 0 in tesla.
///
/// # Errors
///
/// - [`BijliError::DivisionByZero`] if `mass` is zero.
/// - [`BijliError::InvalidParameter`] if `mass` is negative, `b_magnitude`
///   is negative, or any argument is NaN.
#[inline]
pub fn cyclotron_frequency(charge: f64, mass: f64, b_magnitude: f64) -> Result<f64> {
    check_mass(mass, "cyclotron frequency")?;
    if charge.is_nan() {
        return Err(BijliError::InvalidParameter {
            reason: "charge must not be NaN".into(),
        });
    }
    if b_magnitude.is_nan() || b_magnitude < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|B| must be non-negative, got {b_magnitude}"),
        });
    }
    Ok(charge.abs() * b_magnitude / mass)
}

/// Larmor radius (gyroradius) of a charged particle.
///
/// r_L = m|v_⊥|/(|q|B) (meters), with `b_magnitude` = |B| > 0 in tesla. Any
/// non-zero charge and field strength is accepted (e.g. an electron in a
/// 10⁻¹² T interplanetary field).
///
/// # Errors
///
/// - [`BijliError::DivisionByZero`] if `mass`, `charge` or `b_magnitude` is zero.
/// - [`BijliError::InvalidParameter`] if `mass` or `b_magnitude` is negative,
///   or any argument is NaN.
#[inline]
pub fn larmor_radius(mass: f64, v_perp: f64, charge: f64, b_magnitude: f64) -> Result<f64> {
    check_mass(mass, "Larmor radius")?;
    if v_perp.is_nan() || charge.is_nan() {
        return Err(BijliError::InvalidParameter {
            reason: "v_perp and charge must not be NaN".into(),
        });
    }
    if charge == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "charge cannot be zero for Larmor radius".into(),
        });
    }
    if b_magnitude == 0.0 {
        return Err(BijliError::DivisionByZero {
            context: "|B| cannot be zero for Larmor radius".into(),
        });
    }
    if b_magnitude.is_nan() || b_magnitude <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("|B| must be positive, got {b_magnitude}"),
        });
    }
    // Divide in two steps so a tiny |q|·B product cannot underflow to zero.
    Ok(mass * v_perp.abs() / charge.abs() / b_magnitude)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coulomb_force_repulsion() {
        let q1 = PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]).unwrap();
        let q2 = PointCharge::new(1e-6, 1.0, [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
        let f = coulomb_force(&q1, &q2).unwrap();
        // Like charges repel: force on q1 is in -x direction (away from q2)
        assert!(f.x < 0.0);
    }

    #[test]
    fn test_coulomb_force_attraction() {
        let q1 = PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]).unwrap();
        let q2 = PointCharge::new(-1e-6, 1.0, [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
        let f = coulomb_force(&q1, &q2).unwrap();
        // Opposite charges attract: force on q1 is in +x direction (toward q2)
        assert!(f.x > 0.0);
    }

    #[test]
    fn test_lorentz_force_electric_only() {
        let v = FieldVector::zero();
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let f = lorentz_force(ELEMENTARY_CHARGE, &v, &e, &b);
        assert!((f.x - ELEMENTARY_CHARGE * 1000.0).abs() < 1e-25);
    }

    #[test]
    fn test_lorentz_force_magnetic_only() {
        let v = FieldVector::new(1e6, 0.0, 0.0);
        let e = FieldVector::zero();
        let b = FieldVector::new(0.0, 0.0, 1.0);
        let f = lorentz_force(ELEMENTARY_CHARGE, &v, &e, &b);
        // v × B = (1e6, 0, 0) × (0, 0, 1) = (0, -1e6, 0)
        assert!(f.y < 0.0);
    }

    #[test]
    fn test_electron_proton() {
        let e = PointCharge::electron([0.0; 3]);
        let p = PointCharge::proton([1e-10, 0.0, 0.0]);
        let f = coulomb_force(&e, &p).unwrap();
        // Electron attracted to proton
        assert!(f.x > 0.0);
    }

    #[test]
    fn test_cyclotron_frequency() {
        let omega = cyclotron_frequency(ELEMENTARY_CHARGE, ELECTRON_MASS, 1.0).unwrap();
        // ω_c = eB/m ≈ 1.76e11 rad/s for electron in 1T
        assert!((omega - 1.76e11).abs() / 1.76e11 < 0.01);
    }

    #[test]
    fn test_coulomb_potential_energy() {
        let u = coulomb_potential_energy(ELEMENTARY_CHARGE, -ELEMENTARY_CHARGE, 5.3e-11).unwrap();
        // Hydrogen ground state: ≈ -4.36e-18 J ≈ -27.2 eV
        assert!(u < 0.0);
    }

    #[test]
    fn test_dipole_moment() {
        let p = dipole_moment(ELEMENTARY_CHARGE, 1e-10);
        assert!((p - ELEMENTARY_CHARGE * 1e-10).abs() < 1e-40);
    }

    #[test]
    fn test_new_rejects_zero_mass() {
        let result = PointCharge::new(1e-6, 0.0, [0.0; 3], [0.0; 3]);
        assert!(result.is_err());
    }

    #[test]
    fn test_new_rejects_negative_mass() {
        let result = PointCharge::new(1e-6, -1.0, [0.0; 3], [0.0; 3]);
        assert!(result.is_err());
    }

    #[test]
    fn test_potential_energy() {
        // 1C at 1V → 1J
        let u = potential_energy(1.0, 1.0);
        assert!((u - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_larmor_radius() {
        // Electron at 1e6 m/s in 1T field
        let r = larmor_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, 1.0).unwrap();
        // r_L = mv/(eB) ≈ 5.69e-6 m
        assert!((r - 5.69e-6).abs() / 5.69e-6 < 0.01);
    }

    #[test]
    fn test_larmor_radius_zero_b_fails() {
        let result = larmor_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_cyclotron_frequency_zero_mass_fails() {
        let result = cyclotron_frequency(ELEMENTARY_CHARGE, 0.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_coulomb_potential_energy_singularity() {
        let result = coulomb_potential_energy(1.0, 1.0, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_coulomb_force_singularity() {
        let q1 = PointCharge::new(1e-6, 1.0, [0.0; 3], [0.0; 3]).unwrap();
        let q2 = PointCharge::new(1e-6, 1.0, [0.0; 3], [0.0; 3]).unwrap();
        assert!(coulomb_force(&q1, &q2).is_err());
    }

    // ── PointCharge convenience method tests ──────────────────────

    #[test]
    fn test_point_charge_electric_field_at() {
        let q = PointCharge::new(1e-6, 1.0, [0.0; 3], [0.0; 3]).unwrap();
        let e = q.electric_field_at([1.0, 0.0, 0.0]).unwrap();
        assert!(e.x > 0.0); // positive charge → field points away
    }

    #[test]
    fn test_point_charge_electric_potential_at() {
        let q = PointCharge::new(1e-6, 1.0, [0.0; 3], [0.0; 3]).unwrap();
        let v = q.electric_potential_at([1.0, 0.0, 0.0]).unwrap();
        assert!(v > 0.0); // positive charge → positive potential
    }

    #[test]
    fn test_point_charge_magnetic_field_at() {
        let q = PointCharge::new(1e-6, 1.0, [0.0; 3], [1e6, 0.0, 0.0]).unwrap();
        let b = q.magnetic_field_at([0.0, 1.0, 0.0]).unwrap();
        assert!(b.z > 0.0); // v×r̂ = x̂×ŷ = ẑ
    }

    #[test]
    fn test_point_charge_lorentz_force() {
        let q =
            PointCharge::new(ELEMENTARY_CHARGE, PROTON_MASS, [0.0; 3], [1e6, 0.0, 0.0]).unwrap();
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let f = q.lorentz_force_in(&e, &b);
        assert!(f.x > 0.0); // qE in +x
    }

    #[test]
    fn test_point_charge_singularity() {
        let q = PointCharge::new(1e-6, 1.0, [0.0; 3], [0.0; 3]).unwrap();
        assert!(q.electric_field_at([0.0; 3]).is_err());
    }

    // ── Audit regression tests (1.2.0) ────────────────────────────

    #[test]
    fn test_larmor_radius_weak_field() {
        // Electron at 1e3 m/s in 1e-12 T (|q|B ≈ 1.6e-31 < old 1e-30 floor).
        let b = 1e-12;
        let v = 1e3;
        let r = larmor_radius(ELECTRON_MASS, v, -ELEMENTARY_CHARGE, b).unwrap();
        // r_L = m v / (|q| B) = 9.109383702e-31·1e3 / (1.602176634e-19·1e-12)
        let expected = 5.685_630_103_877_797e3;
        assert!((r - expected).abs() / expected < 1e-9, "{r}");
    }

    #[test]
    fn test_larmor_radius_rejects_invalid() {
        assert!(larmor_radius(ELECTRON_MASS, 1e6, 0.0, 1.0).is_err());
        assert!(larmor_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, -1.0).is_err());
        assert!(larmor_radius(-ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, 1.0).is_err());
        assert!(larmor_radius(ELECTRON_MASS, f64::NAN, ELEMENTARY_CHARGE, 1.0).is_err());
        assert!(larmor_radius(ELECTRON_MASS, 1e6, f64::NAN, 1.0).is_err());
        assert!(larmor_radius(ELECTRON_MASS, 1e6, ELEMENTARY_CHARGE, f64::NAN).is_err());
    }

    #[test]
    fn test_cyclotron_frequency_rejects_invalid() {
        assert!(cyclotron_frequency(ELEMENTARY_CHARGE, -ELECTRON_MASS, 1.0).is_err());
        assert!(cyclotron_frequency(ELEMENTARY_CHARGE, f64::NAN, 1.0).is_err());
        assert!(cyclotron_frequency(ELEMENTARY_CHARGE, ELECTRON_MASS, -1.0).is_err());
        assert!(cyclotron_frequency(f64::NAN, ELECTRON_MASS, 1.0).is_err());
        // B = 0 is a valid (zero) frequency.
        assert_eq!(
            cyclotron_frequency(ELEMENTARY_CHARGE, ELECTRON_MASS, 0.0).unwrap(),
            0.0
        );
    }

    #[test]
    fn test_coulomb_potential_energy_rejects_negative_distance() {
        assert!(coulomb_potential_energy(1.0, 1.0, -1.0).is_err());
        assert!(coulomb_potential_energy(1.0, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn test_new_rejects_nan_mass() {
        assert!(PointCharge::new(1e-6, f64::NAN, [0.0; 3], [0.0; 3]).is_err());
    }
}
