//! Electric and magnetic vector fields, field lines, potentials.
//!
//! All SI units: meters, coulombs, volts/meter, tesla.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

/// Permittivity of free space ε₀ (F/m).
pub const EPSILON_0: f64 = 8.854_187_817e-12;

/// Permeability of free space μ₀ (H/m).
pub const MU_0: f64 = 1.256_637_062e-6;

/// Speed of light c = 1/√(ε₀μ₀) (m/s).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Coulomb's constant k = 1/(4πε₀) (N⋅m²/C²).
pub const COULOMB_K: f64 = 8.987_551_792e9;

/// A 3D vector field value at a point.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FieldVector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl FieldVector {
    /// Create a new field vector.
    #[inline]
    #[must_use]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    #[inline]
    #[must_use]
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Magnitude of the vector.
    #[inline]
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Squared magnitude (avoids sqrt).
    #[inline]
    #[must_use]
    pub fn magnitude_sq(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Unit vector in this direction.
    pub fn normalized(&self) -> Result<Self> {
        let mag = self.magnitude();
        if mag < 1e-15 {
            return Err(BijliError::DivisionByZero {
                context: "zero vector cannot be normalized".into(),
            });
        }
        Ok(Self {
            x: self.x / mag,
            y: self.y / mag,
            z: self.z / mag,
        })
    }

    /// Dot product.
    #[inline]
    #[must_use]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    #[inline]
    #[must_use]
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Scale by a constant.
    #[inline]
    #[must_use]
    pub fn scale(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    /// Add two field vectors.
    #[inline]
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

/// Electric field E at a point due to a point charge.
///
/// E = kq/r² r̂ (V/m)
pub fn electric_field_point_charge(
    charge_q: f64,
    charge_pos: [f64; 3],
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let dx = field_pos[0] - charge_pos[0];
    let dy = field_pos[1] - charge_pos[1];
    let dz = field_pos[2] - charge_pos[2];
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let r = r_sq.sqrt();
    let e_mag = COULOMB_K * charge_q / r_sq;

    Ok(FieldVector {
        x: e_mag * dx / r,
        y: e_mag * dy / r,
        z: e_mag * dz / r,
    })
}

/// Electric potential V at a point due to a point charge.
///
/// V = kq/r (volts)
pub fn electric_potential_point_charge(
    charge_q: f64,
    charge_pos: [f64; 3],
    field_pos: [f64; 3],
) -> Result<f64> {
    let dx = field_pos[0] - charge_pos[0];
    let dy = field_pos[1] - charge_pos[1];
    let dz = field_pos[2] - charge_pos[2];
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    Ok(COULOMB_K * charge_q / r_sq.sqrt())
}

/// Magnetic field B at a point due to a moving point charge (Biot-Savart).
///
/// B = (μ₀/4π) × q(v × r̂)/r² (tesla)
pub fn magnetic_field_moving_charge(
    charge_q: f64,
    charge_pos: [f64; 3],
    velocity: [f64; 3],
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let dx = field_pos[0] - charge_pos[0];
    let dy = field_pos[1] - charge_pos[1];
    let dz = field_pos[2] - charge_pos[2];
    let r_sq = dx * dx + dy * dy + dz * dz;

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let r = r_sq.sqrt();
    let r_hat = [dx / r, dy / r, dz / r];

    // v × r̂
    let cross_x = velocity[1] * r_hat[2] - velocity[2] * r_hat[1];
    let cross_y = velocity[2] * r_hat[0] - velocity[0] * r_hat[2];
    let cross_z = velocity[0] * r_hat[1] - velocity[1] * r_hat[0];

    let factor = (MU_0 / (4.0 * std::f64::consts::PI)) * charge_q / r_sq;

    Ok(FieldVector {
        x: factor * cross_x,
        y: factor * cross_y,
        z: factor * cross_z,
    })
}

/// Superposition of electric fields from multiple point charges.
pub fn electric_field_superposition(
    charges: &[(f64, [f64; 3])],
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let mut total = FieldVector::zero();
    for &(q, pos) in charges {
        let e = electric_field_point_charge(q, pos, field_pos)?;
        total = total.add(&e);
    }
    Ok(total)
}

/// Energy stored in an electric field (energy density u = ε₀E²/2).
#[inline]
#[must_use]
pub fn electric_energy_density(e_magnitude: f64) -> f64 {
    0.5 * EPSILON_0 * e_magnitude * e_magnitude
}

/// Energy stored in a magnetic field (energy density u = B²/(2μ₀)).
#[inline]
#[must_use]
pub fn magnetic_energy_density(b_magnitude: f64) -> f64 {
    b_magnitude * b_magnitude / (2.0 * MU_0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_vector_magnitude() {
        let v = FieldVector::new(3.0, 4.0, 0.0);
        assert!((v.magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_field_vector_dot() {
        let a = FieldVector::new(1.0, 0.0, 0.0);
        let b = FieldVector::new(0.0, 1.0, 0.0);
        assert!(a.dot(&b).abs() < 1e-10);
    }

    #[test]
    fn test_field_vector_cross() {
        let x = FieldVector::new(1.0, 0.0, 0.0);
        let y = FieldVector::new(0.0, 1.0, 0.0);
        let z = x.cross(&y);
        assert!((z.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_electric_field_positive_charge() {
        let e = electric_field_point_charge(1e-6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        // E should point away from positive charge (positive x)
        assert!(e.x > 0.0);
        assert!(e.y.abs() < 1e-10);
    }

    #[test]
    fn test_electric_field_negative_charge() {
        let e = electric_field_point_charge(-1e-6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        // E should point toward negative charge (negative x)
        assert!(e.x < 0.0);
    }

    #[test]
    fn test_singularity_at_source() {
        let result = electric_field_point_charge(1e-6, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(result.is_err());
    }

    #[test]
    fn test_electric_potential() {
        let v = electric_potential_point_charge(1e-6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        // V = kq/r ≈ 8987.55 V
        assert!((v - 8987.551_792).abs() < 1.0);
    }

    #[test]
    fn test_inverse_square_law() {
        let e1 = electric_field_point_charge(1e-6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        let e2 = electric_field_point_charge(1e-6, [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]).unwrap();
        // E at 2r should be E at r / 4
        let ratio = e1.magnitude() / e2.magnitude();
        assert!((ratio - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_superposition_cancellation() {
        let charges = vec![
            (1e-6, [-1.0, 0.0, 0.0]),
            (1e-6, [1.0, 0.0, 0.0]),
        ];
        let e = electric_field_superposition(&charges, [0.0, 0.0, 0.0]).unwrap();
        // Symmetric charges: E_x cancels at midpoint
        assert!(e.x.abs() < 1e-6);
    }

    #[test]
    fn test_magnetic_field_moving_charge() {
        let b = magnetic_field_moving_charge(
            1e-6,
            [0.0, 0.0, 0.0],
            [1e6, 0.0, 0.0], // moving in +x
            [0.0, 1.0, 0.0], // field point at +y
        ).unwrap();
        // v × r̂ = x̂ × ŷ = ẑ → B should be in +z direction
        assert!(b.z > 0.0);
    }

    #[test]
    fn test_energy_density_electric() {
        let u = electric_energy_density(1000.0); // 1 kV/m
        assert!(u > 0.0);
        assert!((u - 0.5 * EPSILON_0 * 1e6).abs() < 1e-10);
    }

    #[test]
    fn test_energy_density_magnetic() {
        let u = magnetic_energy_density(1.0); // 1 tesla
        assert!(u > 0.0);
    }

    #[test]
    fn test_normalized() {
        let v = FieldVector::new(3.0, 4.0, 0.0);
        let n = v.normalized().unwrap();
        assert!((n.magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_vector_normalize_fails() {
        let v = FieldVector::zero();
        assert!(v.normalized().is_err());
    }
}
