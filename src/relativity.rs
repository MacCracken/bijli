//! Relativistic electrodynamics — EM tensor, Lorentz field transforms,
//! retarded potentials, Liénard-Wiechert fields.
//!
//! All SI units. Uses the (−,+,+,+) metric signature.
//! Indices: 0=t, 1=x, 2=y, 3=z.

use crate::error::{BijliError, Result};
use crate::field::{FieldVector, MU_0, SPEED_OF_LIGHT, SPEED_OF_LIGHT_SQ};

/// Speed of light squared, alias for readability.
const C2: f64 = SPEED_OF_LIGHT_SQ;

/// Speed of light cubed.
const C3: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT_SQ;

// ── Lorentz factor ─────────────────────────────────────────────────

/// Lorentz factor γ = 1/√(1 − v²/c²).
#[inline]
pub fn lorentz_factor(speed: f64) -> Result<f64> {
    let beta_sq = speed * speed / C2;
    if beta_sq >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("speed ({speed} m/s) must be less than c"),
        });
    }
    Ok(1.0 / (1.0 - beta_sq).sqrt())
}

/// Velocity ratio β = v/c.
#[inline]
#[must_use]
pub fn beta(speed: f64) -> f64 {
    speed / SPEED_OF_LIGHT
}

// ── Electromagnetic tensor ─────────────────────────────────────────

/// The electromagnetic field tensor F^μν (contravariant).
///
/// A 4×4 antisymmetric matrix encoding E and B fields:
/// ```text
///        ⎡  0    -Ex/c  -Ey/c  -Ez/c ⎤
/// F^μν = ⎢ Ex/c    0    -Bz     By   ⎥
///        ⎢ Ey/c   Bz      0    -Bx   ⎥
///        ⎣ Ez/c  -By     Bx      0    ⎦
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmTensor {
    /// The 4×4 matrix F^μν stored in row-major order.
    pub components: [[f64; 4]; 4],
}

impl EmTensor {
    /// Construct the EM tensor from E and B fields.
    #[inline]
    #[must_use]
    pub fn from_fields(e: &FieldVector, b: &FieldVector) -> Self {
        let ex_c = e.x / SPEED_OF_LIGHT;
        let ey_c = e.y / SPEED_OF_LIGHT;
        let ez_c = e.z / SPEED_OF_LIGHT;

        Self {
            components: [
                [0.0, -ex_c, -ey_c, -ez_c],
                [ex_c, 0.0, -b.z, b.y],
                [ey_c, b.z, 0.0, -b.x],
                [ez_c, -b.y, b.x, 0.0],
            ],
        }
    }

    /// Extract the electric field from the tensor.
    #[inline]
    #[must_use]
    pub fn electric_field(&self) -> FieldVector {
        let f = &self.components;
        FieldVector::new(
            f[1][0] * SPEED_OF_LIGHT,
            f[2][0] * SPEED_OF_LIGHT,
            f[3][0] * SPEED_OF_LIGHT,
        )
    }

    /// Extract the magnetic field from the tensor.
    #[inline]
    #[must_use]
    pub fn magnetic_field(&self) -> FieldVector {
        let f = &self.components;
        FieldVector::new(f[3][2], f[1][3], f[2][1])
    }

    /// First Lorentz invariant: F_μν F^μν = 2(B² − E²/c²).
    #[inline]
    #[must_use]
    pub fn first_invariant(&self) -> f64 {
        let e = self.electric_field();
        let b = self.magnetic_field();
        2.0 * (b.magnitude_sq() - e.magnitude_sq() / C2)
    }

    /// Second Lorentz invariant: ε_μνρσ F^μν F^ρσ / 8 = −E⋅B/c.
    ///
    /// Proportional to E⋅B — vanishes for radiation fields.
    #[inline]
    #[must_use]
    pub fn second_invariant(&self) -> f64 {
        let e = self.electric_field();
        let b = self.magnetic_field();
        -e.dot(&b) / SPEED_OF_LIGHT
    }

    /// Dual tensor *F^μν (Hodge dual), swaps E↔B roles.
    #[inline]
    #[must_use]
    pub fn dual(&self) -> Self {
        let e = self.electric_field();
        let b = self.magnetic_field();
        // Dual: E → cB, B → −E/c
        let e_dual = b.scale(SPEED_OF_LIGHT);
        let b_dual = e.scale(-1.0 / SPEED_OF_LIGHT);
        Self::from_fields(&e_dual, &b_dual)
    }
}

// ── Lorentz transformation of fields ───────────────────────────────

/// Transform E and B fields to a frame moving with velocity v along an
/// arbitrary direction.
///
/// Returns (E', B') in the boosted frame.
///
/// Uses the general field transformation:
///   E'∥ = E∥
///   E'⊥ = γ(E⊥ + v × B)
///   B'∥ = B∥
///   B'⊥ = γ(B⊥ − v × E/c²)
#[inline]
pub fn lorentz_transform_fields(
    e: &FieldVector,
    b: &FieldVector,
    velocity: &FieldVector,
) -> Result<(FieldVector, FieldVector)> {
    let v_sq = velocity.magnitude_sq();
    if v_sq >= C2 {
        return Err(BijliError::InvalidParameter {
            reason: "boost velocity must be less than c".into(),
        });
    }

    if v_sq < 1e-30 {
        // No boost — return fields unchanged
        return Ok((*e, *b));
    }

    let gamma = 1.0 / (1.0 - v_sq / C2).sqrt();
    let v_hat = FieldVector::new(
        velocity.x / v_sq.sqrt(),
        velocity.y / v_sq.sqrt(),
        velocity.z / v_sq.sqrt(),
    );

    // Parallel components (along v)
    let e_par_mag = e.dot(&v_hat);
    let b_par_mag = b.dot(&v_hat);
    let e_par = v_hat.scale(e_par_mag);
    let b_par = v_hat.scale(b_par_mag);

    // Perpendicular components
    let e_perp = *e - e_par;
    let b_perp = *b - b_par;

    // v × B and v × E
    let v_cross_b = velocity.cross(b);
    let v_cross_e = velocity.cross(e);

    // Transform
    let e_prime = e_par + (e_perp + v_cross_b) * gamma;
    let b_prime = b_par + (b_perp - v_cross_e.scale(1.0 / C2)) * gamma;

    Ok((e_prime, b_prime))
}

/// Special case: boost along x-axis with speed v.
///
/// E'_x = E_x           B'_x = B_x
/// E'_y = γ(E_y − vB_z) B'_y = γ(B_y + vE_z/c²)
/// E'_z = γ(E_z + vB_y) B'_z = γ(B_z − vE_y/c²)
#[inline]
pub fn lorentz_transform_x(
    e: &FieldVector,
    b: &FieldVector,
    speed: f64,
) -> Result<(FieldVector, FieldVector)> {
    let gamma = lorentz_factor(speed)?;

    let e_prime = FieldVector::new(
        e.x,
        gamma * (e.y - speed * b.z),
        gamma * (e.z + speed * b.y),
    );
    let b_prime = FieldVector::new(
        b.x,
        gamma * (b.y + speed * e.z / C2),
        gamma * (b.z - speed * e.y / C2),
    );

    Ok((e_prime, b_prime))
}

// ── Four-vectors ───────────────────────────────────────────────────

/// A four-vector (ct, x, y, z).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourVector {
    pub t: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl FourVector {
    /// Create a new four-vector from (ct, x, y, z).
    #[inline]
    #[must_use]
    pub fn new(ct: f64, x: f64, y: f64, z: f64) -> Self {
        Self { t: ct, x, y, z }
    }

    /// Create a four-position from time and spatial coordinates.
    #[inline]
    #[must_use]
    pub fn from_time_position(time: f64, pos: [f64; 3]) -> Self {
        Self {
            t: SPEED_OF_LIGHT * time,
            x: pos[0],
            y: pos[1],
            z: pos[2],
        }
    }

    /// Minkowski inner product: η_μν A^μ B^ν = −A⁰B⁰ + A¹B¹ + A²B² + A³B³.
    #[inline]
    #[must_use]
    pub fn dot(&self, other: &Self) -> f64 {
        -self.t * other.t + self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Squared interval: s² = −(ct)² + x² + y² + z².
    #[inline]
    #[must_use]
    pub fn interval_sq(&self) -> f64 {
        -self.t * self.t + self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Spatial magnitude.
    #[inline]
    #[must_use]
    pub fn spatial_magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Boost along x-axis with speed v.
    #[inline]
    pub fn boost_x(&self, speed: f64) -> Result<Self> {
        let gamma = lorentz_factor(speed)?;
        let beta_val = speed / SPEED_OF_LIGHT;
        Ok(Self {
            t: gamma * (self.t - beta_val * self.x),
            x: gamma * (self.x - beta_val * self.t),
            y: self.y,
            z: self.z,
        })
    }
}

// ── Four-potential and retarded potentials ──────────────────────────

/// Four-potential A^μ = (φ/c, A_x, A_y, A_z).
#[derive(Debug, Clone, Copy)]
pub struct FourPotential {
    /// Scalar potential φ/c.
    pub phi_over_c: f64,
    /// Vector potential components.
    pub a: FieldVector,
}

/// Retarded scalar potential from a point charge.
///
/// φ(r, t) = q / (4πε₀ |r − r'(t_r)|)
///
/// where t_r is the retarded time: t_r = t − |r − r'(t_r)|/c.
///
/// For a charge at rest at `charge_pos`, this reduces to the Coulomb potential.
/// For a moving charge, `retarded_distance` is |r − r'(t_r)|.
#[inline]
pub fn retarded_scalar_potential(charge: f64, retarded_distance: f64) -> Result<f64> {
    if retarded_distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(charge / (4.0 * std::f64::consts::PI * crate::field::EPSILON_0 * retarded_distance))
}

/// Retarded vector potential from a moving point charge.
///
/// A(r, t) = μ₀ q v / (4π |r − r'(t_r)|)
#[inline]
pub fn retarded_vector_potential(
    charge: f64,
    velocity: &FieldVector,
    retarded_distance: f64,
) -> Result<FieldVector> {
    if retarded_distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let factor = MU_0 * charge / (4.0 * std::f64::consts::PI * retarded_distance);
    Ok(velocity.scale(factor))
}

// ── Liénard-Wiechert fields ────────────────────────────────────────

/// Liénard-Wiechert electric field from a moving point charge.
///
/// E = (q/4πε₀) × (n̂ − β) / (γ²(1 − n̂⋅β)³ R²)
///   + (q/4πε₀c) × n̂ × [(n̂ − β) × β̇] / ((1 − n̂⋅β)³ R)
///
/// The first term is the "velocity field" (Coulomb-like, ∝ 1/R²).
/// The second term is the "acceleration field" (radiation, ∝ 1/R).
///
/// # Parameters
/// - `charge`: charge in coulombs
/// - `r_vec`: vector from retarded charge position to field point (R = |r_vec|)
/// - `beta_vec`: velocity/c at retarded time
/// - `beta_dot`: d(β)/dt at retarded time (acceleration/c)
///
/// # Returns
/// The electric field at the field point.
pub fn lienard_wiechert_e(
    charge: f64,
    r_vec: &FieldVector,
    beta_vec: &FieldVector,
    beta_dot: &FieldVector,
) -> Result<FieldVector> {
    let r_mag = r_vec.magnitude();
    if r_mag < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let beta_sq = beta_vec.magnitude_sq();
    if beta_sq >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: "β must be < 1".into(),
        });
    }

    let n_hat = r_vec.scale(1.0 / r_mag);
    let n_dot_beta = n_hat.dot(beta_vec);
    let kappa = 1.0 - n_dot_beta; // (1 − n̂⋅β)

    if kappa.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "Liénard-Wiechert denominator (1 − n̂⋅β) is zero".into(),
        });
    }

    let gamma_sq = 1.0 / (1.0 - beta_sq);
    let kappa3 = kappa * kappa * kappa;
    let prefactor = charge / (4.0 * std::f64::consts::PI * crate::field::EPSILON_0);

    // Velocity (Coulomb) term: (n̂ − β) / (γ²κ³R²)
    let n_minus_beta = n_hat - *beta_vec;
    let velocity_term = n_minus_beta.scale(prefactor / (gamma_sq * kappa3 * r_mag * r_mag));

    // Acceleration (radiation) term: n̂ × [(n̂ − β) × β̇] / (κ³Rc)
    let cross_inner = n_minus_beta.cross(beta_dot);
    let cross_outer = n_hat.cross(&cross_inner);
    let accel_term = cross_outer.scale(prefactor / (kappa3 * r_mag * SPEED_OF_LIGHT));

    Ok(velocity_term + accel_term)
}

/// Liénard-Wiechert magnetic field: B = n̂ × E / c.
///
/// The B field is always perpendicular to both n̂ and E.
#[inline]
#[must_use]
pub fn lienard_wiechert_b(n_hat: &FieldVector, e_field: &FieldVector) -> FieldVector {
    n_hat.cross(e_field).scale(1.0 / SPEED_OF_LIGHT)
}

/// Larmor formula: power radiated by an accelerating charge.
///
/// P = q²a² / (6πε₀c³) (non-relativistic)
#[inline]
#[must_use]
pub fn larmor_power(charge: f64, acceleration: f64) -> f64 {
    let c3 = C3;
    charge * charge * acceleration * acceleration
        / (6.0 * std::f64::consts::PI * crate::field::EPSILON_0 * c3)
}

/// Relativistic Larmor formula (Liénard extension).
///
/// P = q²γ⁶ / (6πε₀c) × [(a)² − (v × a)²/c²]
///
/// In the rest frame this reduces to the non-relativistic Larmor formula.
#[inline]
pub fn relativistic_larmor_power(
    charge: f64,
    velocity: &FieldVector,
    acceleration: &FieldVector,
) -> Result<f64> {
    let v_sq = velocity.magnitude_sq();
    if v_sq >= C2 {
        return Err(BijliError::InvalidParameter {
            reason: "velocity must be less than c".into(),
        });
    }

    let gamma = 1.0 / (1.0 - v_sq / C2).sqrt();
    let gamma2 = gamma * gamma;
    let gamma6 = gamma2 * gamma2 * gamma2;
    let a_sq = acceleration.magnitude_sq();
    let v_cross_a = velocity.cross(acceleration);
    let v_cross_a_sq = v_cross_a.magnitude_sq();

    let c3 = C3;
    let prefactor = charge * charge / (6.0 * std::f64::consts::PI * crate::field::EPSILON_0 * c3);

    Ok(prefactor * gamma6 * (a_sq - v_cross_a_sq / C2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::EPSILON_0;

    // ── Lorentz factor ─────────────────────────────────────────────

    #[test]
    fn test_lorentz_factor_zero() {
        let gamma = lorentz_factor(0.0).unwrap();
        assert!((gamma - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_lorentz_factor_half_c() {
        let gamma = lorentz_factor(0.5 * SPEED_OF_LIGHT).unwrap();
        let expected = 1.0 / (1.0 - 0.25_f64).sqrt();
        assert!((gamma - expected).abs() < 1e-10);
    }

    #[test]
    fn test_lorentz_factor_superluminal() {
        assert!(lorentz_factor(SPEED_OF_LIGHT).is_err());
        assert!(lorentz_factor(1.1 * SPEED_OF_LIGHT).is_err());
    }

    #[test]
    fn test_beta() {
        assert!((beta(SPEED_OF_LIGHT) - 1.0).abs() < 1e-15);
        assert!((beta(0.5 * SPEED_OF_LIGHT) - 0.5).abs() < 1e-15);
    }

    // ── EM tensor ──────────────────────────────────────────────────

    #[test]
    fn test_em_tensor_roundtrip() {
        let e = FieldVector::new(1e3, 2e3, 3e3);
        let b = FieldVector::new(1e-3, 2e-3, 3e-3);
        let f = EmTensor::from_fields(&e, &b);
        let e_back = f.electric_field();
        let b_back = f.magnetic_field();
        assert!((e_back.x - e.x).abs() < 1e-6);
        assert!((e_back.y - e.y).abs() < 1e-6);
        assert!((e_back.z - e.z).abs() < 1e-6);
        assert!((b_back.x - b.x).abs() < 1e-10);
        assert!((b_back.y - b.y).abs() < 1e-10);
        assert!((b_back.z - b.z).abs() < 1e-10);
    }

    #[test]
    fn test_em_tensor_antisymmetric() {
        let e = FieldVector::new(1e3, 0.0, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1e-3);
        let f = EmTensor::from_fields(&e, &b);
        for mu in 0..4 {
            for nu in 0..4 {
                assert!(
                    (f.components[mu][nu] + f.components[nu][mu]).abs() < 1e-20,
                    "F^{mu}{nu} + F^{nu}{mu} != 0"
                );
            }
        }
    }

    #[test]
    fn test_first_invariant_radiation() {
        // For radiation, |E| = c|B|, so B² - E²/c² = 0
        let e = FieldVector::new(0.0, SPEED_OF_LIGHT, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1.0);
        let f = EmTensor::from_fields(&e, &b);
        assert!(f.first_invariant().abs() < 1e-10);
    }

    #[test]
    fn test_second_invariant_radiation() {
        // For radiation, E ⊥ B → E⋅B = 0
        let e = FieldVector::new(0.0, SPEED_OF_LIGHT, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1.0);
        let f = EmTensor::from_fields(&e, &b);
        assert!(f.second_invariant().abs() < 1e-10);
    }

    #[test]
    fn test_dual_tensor() {
        let e = FieldVector::new(1e3, 0.0, 0.0);
        let b = FieldVector::new(0.0, 1e-3, 0.0);
        let f = EmTensor::from_fields(&e, &b);
        let f_dual = f.dual();
        // Dual of dual should give back -F
        let f_dd = f_dual.dual();
        for mu in 0..4 {
            for nu in 0..4 {
                assert!(
                    (f_dd.components[mu][nu] + f.components[mu][nu]).abs() < 1e-10,
                    "dual(dual(F)) != -F"
                );
            }
        }
    }

    // ── Lorentz transform ──────────────────────────────────────────

    #[test]
    fn test_lorentz_transform_zero_velocity() {
        let e = FieldVector::new(1e3, 0.0, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1e-3);
        let (ep, bp) = lorentz_transform_x(&e, &b, 0.0).unwrap();
        assert!((ep.x - e.x).abs() < 1e-10);
        assert!((bp.z - b.z).abs() < 1e-15);
    }

    #[test]
    fn test_lorentz_transform_pure_e_creates_b() {
        // Boost along x with pure E in y → should create B in z
        let e = FieldVector::new(0.0, 1e3, 0.0);
        let b = FieldVector::zero();
        let v = 0.5 * SPEED_OF_LIGHT;
        let (_, bp) = lorentz_transform_x(&e, &b, v).unwrap();
        assert!(bp.z.abs() > 0.0);
    }

    #[test]
    fn test_lorentz_transform_preserves_parallel() {
        // E along x (parallel to boost) should be unchanged
        let e = FieldVector::new(1e3, 0.0, 0.0);
        let b = FieldVector::zero();
        let (ep, _) = lorentz_transform_x(&e, &b, 0.5 * SPEED_OF_LIGHT).unwrap();
        assert!((ep.x - e.x).abs() < 1e-10);
    }

    #[test]
    fn test_lorentz_transform_invariant_preserved() {
        let e = FieldVector::new(1e3, 2e3, 3e3);
        let b = FieldVector::new(1e-6, 2e-6, 3e-6);
        let f1 = EmTensor::from_fields(&e, &b);
        let inv1 = f1.first_invariant();

        let (ep, bp) = lorentz_transform_x(&e, &b, 0.8 * SPEED_OF_LIGHT).unwrap();
        let f2 = EmTensor::from_fields(&ep, &bp);
        let inv2 = f2.first_invariant();

        assert!((inv1 - inv2).abs() / (inv1.abs() + 1e-30) < 1e-6);
    }

    #[test]
    fn test_general_lorentz_transform_x_matches_special() {
        let e = FieldVector::new(0.0, 1e3, 0.0);
        let b = FieldVector::new(0.0, 0.0, 1e-3);
        let v = 0.3 * SPEED_OF_LIGHT;

        let (ep1, bp1) = lorentz_transform_x(&e, &b, v).unwrap();
        let vel = FieldVector::new(v, 0.0, 0.0);
        let (ep2, bp2) = lorentz_transform_fields(&e, &b, &vel).unwrap();

        assert!((ep1.x - ep2.x).abs() < 1e-6);
        assert!((ep1.y - ep2.y).abs() < 1e-6);
        assert!((ep1.z - ep2.z).abs() < 1e-6);
        assert!((bp1.x - bp2.x).abs() < 1e-12);
        assert!((bp1.y - bp2.y).abs() < 1e-12);
        assert!((bp1.z - bp2.z).abs() < 1e-12);
    }

    #[test]
    fn test_lorentz_transform_superluminal() {
        let e = FieldVector::new(1e3, 0.0, 0.0);
        let b = FieldVector::zero();
        assert!(lorentz_transform_x(&e, &b, SPEED_OF_LIGHT).is_err());
    }

    // ── Four-vector ────────────────────────────────────────────────

    #[test]
    fn test_four_vector_lightlike() {
        // Photon: s² = 0
        let fv = FourVector::new(SPEED_OF_LIGHT, SPEED_OF_LIGHT, 0.0, 0.0);
        assert!(fv.interval_sq().abs() < 1e-10);
    }

    #[test]
    fn test_four_vector_timelike() {
        // Massive particle at rest: s² < 0
        let fv = FourVector::new(SPEED_OF_LIGHT, 0.0, 0.0, 0.0);
        assert!(fv.interval_sq() < 0.0);
    }

    #[test]
    fn test_four_vector_boost_roundtrip() {
        let fv = FourVector::new(SPEED_OF_LIGHT * 1.0, 1.0, 2.0, 3.0);
        let v = 0.5 * SPEED_OF_LIGHT;
        let boosted = fv.boost_x(v).unwrap();
        let back = boosted.boost_x(-v).unwrap();
        assert!((back.t - fv.t).abs() < 1e-6);
        assert!((back.x - fv.x).abs() < 1e-6);
        assert!((back.y - fv.y).abs() < 1e-10);
    }

    #[test]
    fn test_four_vector_interval_invariant() {
        let fv = FourVector::new(SPEED_OF_LIGHT * 2.0, 1.0, 1.0, 1.0);
        let s2_orig = fv.interval_sq();
        let boosted = fv.boost_x(0.6 * SPEED_OF_LIGHT).unwrap();
        let s2_boost = boosted.interval_sq();
        assert!((s2_orig - s2_boost).abs() < 1e-6);
    }

    // ── Retarded potentials ────────────────────────────────────────

    #[test]
    fn test_retarded_scalar_coulomb() {
        // Static charge: retarded potential = Coulomb
        let q = 1e-6;
        let r = 1.0;
        let phi = retarded_scalar_potential(q, r).unwrap();
        let expected = q / (4.0 * std::f64::consts::PI * EPSILON_0 * r);
        assert!((phi - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_retarded_scalar_singularity() {
        assert!(retarded_scalar_potential(1e-6, 0.0).is_err());
    }

    #[test]
    fn test_retarded_vector_potential() {
        let q = 1e-6;
        let v = FieldVector::new(1e6, 0.0, 0.0);
        let r = 1.0;
        let a = retarded_vector_potential(q, &v, r).unwrap();
        let expected = MU_0 * q * 1e6 / (4.0 * std::f64::consts::PI * r);
        assert!((a.x - expected).abs() / expected < 1e-10);
    }

    // ── Liénard-Wiechert ───────────────────────────────────────────

    #[test]
    fn test_lw_static_charge_is_coulomb() {
        // Static charge: LW should reduce to Coulomb field
        let q = 1e-6;
        let r = FieldVector::new(1.0, 0.0, 0.0);
        let beta_vec = FieldVector::zero();
        let beta_dot = FieldVector::zero();
        let e = lienard_wiechert_e(q, &r, &beta_vec, &beta_dot).unwrap();
        let expected = q / (4.0 * std::f64::consts::PI * EPSILON_0);
        assert!((e.x - expected).abs() / expected < 1e-6);
        assert!(e.y.abs() < 1e-6);
        assert!(e.z.abs() < 1e-6);
    }

    #[test]
    fn test_lw_b_perpendicular_to_e() {
        let q = 1e-6;
        let r = FieldVector::new(1.0, 0.0, 0.0);
        let beta_vec = FieldVector::new(0.1, 0.0, 0.0);
        let beta_dot = FieldVector::zero();
        let e = lienard_wiechert_e(q, &r, &beta_vec, &beta_dot).unwrap();
        let n_hat = r.scale(1.0 / r.magnitude());
        let b = lienard_wiechert_b(&n_hat, &e);
        // B should be perpendicular to n̂
        assert!(b.dot(&n_hat).abs() < 1e-10);
    }

    #[test]
    fn test_lw_singularity() {
        let r = FieldVector::zero();
        assert!(lienard_wiechert_e(1e-6, &r, &FieldVector::zero(), &FieldVector::zero()).is_err());
    }

    // ── Larmor radiation ───────────────────────────────────────────

    #[test]
    fn test_larmor_power_positive() {
        let p = larmor_power(1e-6, 1e10);
        assert!(p > 0.0);
    }

    #[test]
    fn test_larmor_zero_acceleration() {
        assert!(larmor_power(1e-6, 0.0).abs() < 1e-50);
    }

    #[test]
    fn test_relativistic_larmor_rest_frame() {
        // At v=0, relativistic formula should equal non-relativistic
        let q = 1e-6;
        let a = 1e10;
        let v = FieldVector::zero();
        let acc = FieldVector::new(a, 0.0, 0.0);
        let p_rel = relativistic_larmor_power(q, &v, &acc).unwrap();
        let p_nr = larmor_power(q, a);
        assert!((p_rel - p_nr).abs() / p_nr < 1e-10);
    }

    #[test]
    fn test_relativistic_larmor_exceeds_nonrelativistic() {
        // At high speeds, relativistic power >> non-relativistic
        let q = 1e-6;
        let a = 1e10;
        let v = FieldVector::new(0.9 * SPEED_OF_LIGHT, 0.0, 0.0);
        let acc = FieldVector::new(0.0, a, 0.0); // perpendicular acceleration
        let p_rel = relativistic_larmor_power(q, &v, &acc).unwrap();
        let p_nr = larmor_power(q, a);
        assert!(p_rel > p_nr);
    }

    #[test]
    fn test_relativistic_larmor_superluminal() {
        let v = FieldVector::new(SPEED_OF_LIGHT, 0.0, 0.0);
        let acc = FieldVector::new(0.0, 1e10, 0.0);
        assert!(relativistic_larmor_power(1e-6, &v, &acc).is_err());
    }
}
