//! Point charges, Coulomb's law, Lorentz force, charge distributions.
//!
//! All SI units: coulombs, newtons, meters, seconds.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};
use crate::field::{FieldVector, COULOMB_K};

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
    #[must_use]
    pub fn new(charge: f64, mass: f64, position: [f64; 3], velocity: [f64; 3]) -> Self {
        Self {
            charge,
            mass,
            position,
            velocity,
        }
    }

    /// Create an electron at rest at a position.
    #[must_use]
    pub fn electron(position: [f64; 3]) -> Self {
        Self::new(-ELEMENTARY_CHARGE, ELECTRON_MASS, position, [0.0; 3])
    }

    /// Create a proton at rest at a position.
    #[must_use]
    pub fn proton(position: [f64; 3]) -> Self {
        Self::new(ELEMENTARY_CHARGE, PROTON_MASS, position, [0.0; 3])
    }
}

/// Coulomb force between two point charges.
///
/// F = kq₁q₂/r² r̂ (newtons, directed from charge2 toward charge1 if repulsive)
pub fn coulomb_force(q1: &PointCharge, q2: &PointCharge) -> Result<FieldVector> {
    let dx = q1.position[0] - q2.position[0];
    let dy = q1.position[1] - q2.position[1];
    let dz = q1.position[2] - q2.position[2];
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let r = r_sq.sqrt();
    let f_mag = COULOMB_K * q1.charge * q2.charge / r_sq;

    Ok(FieldVector {
        x: f_mag * dx / r,
        y: f_mag * dy / r,
        z: f_mag * dz / r,
    })
}

/// Lorentz force on a charged particle: F = q(E + v × B).
#[must_use]
pub fn lorentz_force(
    charge: f64,
    velocity: &FieldVector,
    electric_field: &FieldVector,
    magnetic_field: &FieldVector,
) -> FieldVector {
    let v_cross_b = velocity.cross(magnetic_field);
    let total = electric_field.add(&v_cross_b);
    total.scale(charge)
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
/// U = kq₁q₂/r (joules)
pub fn coulomb_potential_energy(q1: f64, q2: f64, distance: f64) -> Result<f64> {
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(COULOMB_K * q1 * q2 / distance)
}

/// Cyclotron frequency for a charged particle in a magnetic field.
///
/// ω_c = |q|B/m (rad/s)
pub fn cyclotron_frequency(charge: f64, mass: f64, b_magnitude: f64) -> Result<f64> {
    if mass.abs() < 1e-50 {
        return Err(BijliError::DivisionByZero {
            context: "mass cannot be zero for cyclotron frequency".into(),
        });
    }
    Ok(charge.abs() * b_magnitude / mass)
}

/// Larmor radius (gyroradius) of a charged particle.
///
/// r_L = mv_⊥/(|q|B) (meters)
pub fn larmor_radius(mass: f64, v_perp: f64, charge: f64, b_magnitude: f64) -> Result<f64> {
    let denom = charge.abs() * b_magnitude;
    if denom < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "charge×B cannot be zero for Larmor radius".into(),
        });
    }
    Ok(mass * v_perp.abs() / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coulomb_force_repulsion() {
        let q1 = PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]);
        let q2 = PointCharge::new(1e-6, 1.0, [1.0, 0.0, 0.0], [0.0; 3]);
        let f = coulomb_force(&q1, &q2).unwrap();
        // Like charges repel: force on q1 is in -x direction (away from q2)
        assert!(f.x < 0.0);
    }

    #[test]
    fn test_coulomb_force_attraction() {
        let q1 = PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]);
        let q2 = PointCharge::new(-1e-6, 1.0, [1.0, 0.0, 0.0], [0.0; 3]);
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
}
