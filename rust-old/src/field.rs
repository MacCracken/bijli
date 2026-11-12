//! Electric and magnetic vector fields, field lines, potentials.
//!
//! All SI units: meters, coulombs, volts/meter, tesla.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

/// Permittivity of free space ε₀ (F/m), CODATA 2018: 8.854 187 8128(13) × 10⁻¹² F/m.
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// Permeability of free space μ₀ (H/m), CODATA 2018: 1.256 637 062 12(19) × 10⁻⁶ H/m.
pub const MU_0: f64 = 1.256_637_062_12e-6;

/// Speed of light in vacuum c (m/s), exact by definition of the metre.
///
/// c = 1/√(ε₀μ₀); with the CODATA 2018 [`EPSILON_0`] and [`MU_0`] above,
/// ε₀μ₀c² − 1 ≈ −4 × 10⁻¹⁴ (well inside their 1.5 × 10⁻¹⁰ uncertainty).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Coulomb's constant k = 1/(4πε₀) ≈ 8.987 551 7923 × 10⁹ N⋅m²/C² (CODATA 2018).
///
/// Computed from [`EPSILON_0`] at compile time, so `COULOMB_K` and
/// `1/(4π EPSILON_0)` agree to the last bit.
pub const COULOMB_K: f64 = 1.0 / (4.0 * std::f64::consts::PI * EPSILON_0);

/// Speed of light squared c² (m²/s²).
pub const SPEED_OF_LIGHT_SQ: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

/// A 3D vector field value at a point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
    ///
    /// Any finite, non-zero vector can be normalized, whatever its magnitude:
    /// SI-scale values such as 1.6 × 10⁻¹⁹ N or 6 × 10⁻³⁰ C⋅m work directly, and
    /// vectors whose squared magnitude would underflow or overflow are first
    /// rescaled by their largest component.
    ///
    /// # Errors
    ///
    /// - [`BijliError::DivisionByZero`] for the zero vector.
    /// - [`BijliError::InvalidParameter`] if any component is NaN or infinite.
    #[inline]
    pub fn normalized(&self) -> Result<Self> {
        let mag_sq = self.magnitude_sq();
        // Fast path: |v|² is a normal finite number, so |v| is accurate.
        if mag_sq.is_finite() && mag_sq >= f64::MIN_POSITIVE {
            let mag = mag_sq.sqrt();
            return Ok(Self {
                x: self.x / mag,
                y: self.y / mag,
                z: self.z / mag,
            });
        }
        self.normalized_rescaled()
    }

    /// Slow path of [`normalized`](Self::normalized): zero, non-finite, or
    /// components whose squares under/overflow.
    #[cold]
    fn normalized_rescaled(&self) -> Result<Self> {
        if !(self.x.is_finite() && self.y.is_finite() && self.z.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("cannot normalize non-finite vector {self}"),
            });
        }
        let s = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if s == 0.0 {
            return Err(BijliError::DivisionByZero {
                context: "zero vector cannot be normalized".into(),
            });
        }
        let v = Self::new(self.x / s, self.y / s, self.z / s);
        let mag = v.magnitude(); // in [1, √3]
        Ok(Self::new(v.x / mag, v.y / mag, v.z / mag))
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
}

impl Default for FieldVector {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl From<[f64; 3]> for FieldVector {
    #[inline]
    fn from(arr: [f64; 3]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }
}

impl From<FieldVector> for [f64; 3] {
    #[inline]
    fn from(v: FieldVector) -> Self {
        [v.x, v.y, v.z]
    }
}

impl std::ops::Add for FieldVector {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Sub for FieldVector {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for FieldVector {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Neg for FieldVector {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl std::ops::AddAssign for FieldVector {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl std::ops::SubAssign for FieldVector {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl std::ops::MulAssign<f64> for FieldVector {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl std::ops::Mul<FieldVector> for f64 {
    type Output = FieldVector;
    #[inline]
    fn mul(self, rhs: FieldVector) -> FieldVector {
        FieldVector {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl std::fmt::Display for FieldVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

/// Electric field E at a point due to a point charge.
///
/// E = kq/r² r̂ (V/m)
#[inline]
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

    // kq / r³ — avoids separate sqrt + division
    let inv_r = 1.0 / r_sq.sqrt();
    let factor = COULOMB_K * charge_q * inv_r * inv_r * inv_r;

    Ok(FieldVector {
        x: factor * dx,
        y: factor * dy,
        z: factor * dz,
    })
}

/// Electric potential V at a point due to a point charge.
///
/// V = kq/r (volts)
#[inline]
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
#[inline]
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

    let inv_r = 1.0 / r_sq.sqrt();
    let r_hat = [dx * inv_r, dy * inv_r, dz * inv_r];

    // v × r̂
    let cross_x = velocity[1] * r_hat[2] - velocity[2] * r_hat[1];
    let cross_y = velocity[2] * r_hat[0] - velocity[0] * r_hat[2];
    let cross_z = velocity[0] * r_hat[1] - velocity[1] * r_hat[0];

    let factor = (MU_0 / (4.0 * std::f64::consts::PI)) * charge_q * inv_r * inv_r;

    Ok(FieldVector {
        x: factor * cross_x,
        y: factor * cross_y,
        z: factor * cross_z,
    })
}

/// Superposition of electric fields from multiple point charges.
#[inline]
pub fn electric_field_superposition(
    charges: &[(f64, [f64; 3])],
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let mut total = FieldVector::zero();
    for &(q, pos) in charges {
        total += electric_field_point_charge(q, pos, field_pos)?;
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

// ── Dipole fields ──────────────────────────────────────────────────

/// Electric dipole field at position r from a dipole at the origin.
///
/// For a dipole with moment **p** = p ẑ:
///   E_r = 2p cos θ / (4πε₀r³)
///   E_θ = p sin θ / (4πε₀r³)
///
/// General form: **E** = (1/4πε₀r³)[3(**p**⋅**r̂**)**r̂** − **p**]
#[inline]
pub fn electric_dipole_field(
    dipole_moment: &FieldVector,
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let r_sq =
        field_pos[0] * field_pos[0] + field_pos[1] * field_pos[1] + field_pos[2] * field_pos[2];

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let inv_r = 1.0 / r_sq.sqrt();
    let inv_r3 = inv_r * inv_r * inv_r;
    let r_hat = FieldVector::new(
        field_pos[0] * inv_r,
        field_pos[1] * inv_r,
        field_pos[2] * inv_r,
    );

    // E = (1/4πε₀r³)[3(p⋅r̂)r̂ − p]
    let p_dot_rhat = dipole_moment.dot(&r_hat);
    let factor = inv_r3 / (4.0 * std::f64::consts::PI * EPSILON_0);

    Ok(r_hat.scale(3.0 * p_dot_rhat * factor) - dipole_moment.scale(factor))
}

/// Electric dipole potential at position r.
///
/// V = p⋅r̂ / (4πε₀r²)
#[inline]
pub fn electric_dipole_potential(dipole_moment: &FieldVector, field_pos: [f64; 3]) -> Result<f64> {
    let r_sq =
        field_pos[0] * field_pos[0] + field_pos[1] * field_pos[1] + field_pos[2] * field_pos[2];

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let inv_r = 1.0 / r_sq.sqrt();
    let r_hat = FieldVector::new(
        field_pos[0] * inv_r,
        field_pos[1] * inv_r,
        field_pos[2] * inv_r,
    );

    Ok(dipole_moment.dot(&r_hat) * inv_r * inv_r / (4.0 * std::f64::consts::PI * EPSILON_0))
}

/// Magnetic dipole field at position r from a dipole at the origin.
///
/// General form: **B** = (μ₀/4πr³)[3(**m**⋅**r̂**)**r̂** − **m**]
///
/// Same structure as electric dipole but with μ₀/4π instead of 1/(4πε₀).
#[inline]
pub fn magnetic_dipole_field(
    magnetic_moment: &FieldVector,
    field_pos: [f64; 3],
) -> Result<FieldVector> {
    let r_sq =
        field_pos[0] * field_pos[0] + field_pos[1] * field_pos[1] + field_pos[2] * field_pos[2];

    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let inv_r = 1.0 / r_sq.sqrt();
    let inv_r3 = inv_r * inv_r * inv_r;
    let r_hat = FieldVector::new(
        field_pos[0] * inv_r,
        field_pos[1] * inv_r,
        field_pos[2] * inv_r,
    );

    let m_dot_rhat = magnetic_moment.dot(&r_hat);
    let factor = MU_0 * inv_r3 / (4.0 * std::f64::consts::PI);

    Ok(r_hat.scale(3.0 * m_dot_rhat * factor) - magnetic_moment.scale(factor))
}

// ── Gauss's law applications ───────────────────────────────────────

/// Electric field from an infinite plane of uniform surface charge density.
///
/// |E| = σ/(2ε₀), directed away from the plane.
#[inline]
#[must_use]
pub fn electric_field_infinite_plane(surface_charge_density: f64) -> f64 {
    surface_charge_density / (2.0 * EPSILON_0)
}

/// Validate a radial distance argument (m).
///
/// NaN or negative → [`BijliError::InvalidParameter`]; below 10⁻³⁰ m →
/// [`BijliError::Singularity`].
#[inline]
pub(crate) fn check_radial_distance(distance: f64) -> Result<()> {
    if distance.is_nan() || distance < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("distance must be non-negative, got {distance}"),
        });
    }
    if distance < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(())
}

/// Electric field outside a uniformly charged sphere (r > R).
///
/// E_r = Q/(4πε₀r²) (V/m), same as a point charge; `distance` is r ≥ 0 (m).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `distance` is negative or NaN,
/// [`BijliError::Singularity`] if it is (numerically) zero.
#[inline]
pub fn electric_field_charged_sphere(total_charge: f64, distance: f64) -> Result<f64> {
    check_radial_distance(distance)?;
    Ok(COULOMB_K * total_charge / (distance * distance))
}

/// Electric field inside a uniformly charged sphere (0 ≤ r ≤ R).
///
/// E_r = Qr/(4πε₀R³) (V/m) — grows linearly with r and equals the outside
/// value Q/(4πε₀R²) at r = R. For r > R use [`electric_field_charged_sphere`].
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `sphere_radius` is not positive and
/// finite, or `distance` is outside [0, R] (or NaN).
#[inline]
pub fn electric_field_inside_sphere(
    total_charge: f64,
    sphere_radius: f64,
    distance: f64,
) -> Result<f64> {
    if !(sphere_radius > 0.0 && sphere_radius.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("sphere radius must be positive and finite, got {sphere_radius}"),
        });
    }
    if !(distance >= 0.0 && distance <= sphere_radius) {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "distance must lie in [0, R] = [0, {sphere_radius}] for the interior field, got {distance}"
            ),
        });
    }
    Ok(COULOMB_K * total_charge * (distance / sphere_radius) / (sphere_radius * sphere_radius))
}

/// Electric field from an infinite line of uniform linear charge density.
///
/// E_r = λ/(2πε₀r) (V/m), directed radially outward for λ > 0; `distance` is
/// the perpendicular distance r ≥ 0 (m).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `distance` is negative or NaN,
/// [`BijliError::Singularity`] if it is (numerically) zero.
#[inline]
pub fn electric_field_infinite_line(linear_charge_density: f64, distance: f64) -> Result<f64> {
    check_radial_distance(distance)?;
    Ok(linear_charge_density / (2.0 * std::f64::consts::PI * EPSILON_0 * distance))
}

/// Electric field outside an infinite uniformly charged cylinder (r > R).
///
/// E_r = λ/(2πε₀r) where λ is the charge per unit length.
/// Same form (and same errors) as [`electric_field_infinite_line`].
#[inline]
pub fn electric_field_outside_cylinder(linear_charge_density: f64, distance: f64) -> Result<f64> {
    electric_field_infinite_line(linear_charge_density, distance)
}

/// Electric field inside an infinite uniformly charged cylinder (r < R).
///
/// E_r = ρr/(2ε₀) where ρ is the volume charge density.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `distance` is negative or NaN.
#[inline]
pub fn electric_field_inside_cylinder(volume_charge_density: f64, distance: f64) -> Result<f64> {
    if distance.is_nan() || distance < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("distance must be non-negative, got {distance}"),
        });
    }
    Ok(volume_charge_density * distance / (2.0 * EPSILON_0))
}

// ── Charge distributions ───────────────────────────────────────────

/// Electric field on the axis of a uniformly charged ring.
///
/// E_z = Qz / (4πε₀(z² + R²)^(3/2)) (V/m)
///
/// The ring (radius R, total charge Q) lies in the plane z = 0, centred on the
/// axis; `axial_distance` is the signed coordinate z (m), so E_z is odd in z
/// and vanishes at the centre.
///
/// # Errors
///
/// - [`BijliError::InvalidParameter`] if `ring_radius` is negative, NaN or
///   infinite, or `axial_distance` is not finite.
/// - [`BijliError::Singularity`] only for a degenerate ring (R = 0, a point
///   charge) evaluated at its own position (z = 0).
#[inline]
pub fn electric_field_ring_axis(
    total_charge: f64,
    ring_radius: f64,
    axial_distance: f64,
) -> Result<f64> {
    if !(ring_radius >= 0.0 && ring_radius.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("ring radius must be non-negative and finite, got {ring_radius}"),
        });
    }
    if !axial_distance.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("axial distance must be finite, got {axial_distance}"),
        });
    }
    if axial_distance == 0.0 {
        // Centre of the ring: E_z = 0 by symmetry (R > 0); R = 0 is a point charge.
        return if ring_radius == 0.0 {
            Err(BijliError::Singularity)
        } else {
            Ok(0.0)
        };
    }
    let d_sq = axial_distance * axial_distance + ring_radius * ring_radius;
    let denom = d_sq * d_sq.sqrt(); // (z² + R²)^(3/2)
    Ok(COULOMB_K * total_charge * axial_distance / denom)
}

/// Electric field on the axis of a uniformly charged disk.
///
/// ```text
/// E_z = sgn(z)·(σ/2ε₀)·[1 − |z|/√(z² + R²)]
///     = sgn(z)·(σ/2ε₀)·R² / (r(r + |z|)),      r = √(z² + R²)
/// ```
///
/// The disk (radius R, surface charge density σ in C/m²) lies in the plane
/// z = 0; `axial_distance` is the signed coordinate z (m) and E_z (V/m) is odd
/// in z. The second, cancellation-free form is evaluated, so the far field
/// correctly approaches kσπR²/z² instead of rounding to zero.
///
/// The field jumps by σ/ε₀ across the disk. At z = 0 the one-sided limit
/// selected by the sign bit of `axial_distance` is returned: `+0.0` gives
/// +σ/(2ε₀) (limit z → 0⁺), `-0.0` gives −σ/(2ε₀) (limit z → 0⁻).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `disk_radius` is not positive and
/// finite, or `axial_distance` is not finite.
#[inline]
pub fn electric_field_disk_axis(
    surface_charge_density: f64,
    disk_radius: f64,
    axial_distance: f64,
) -> Result<f64> {
    if !(disk_radius > 0.0 && disk_radius.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("disk radius must be positive and finite, got {disk_radius}"),
        });
    }
    if !axial_distance.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("axial distance must be finite, got {axial_distance}"),
        });
    }
    let abs_z = axial_distance.abs();
    let r = (axial_distance * axial_distance + disk_radius * disk_radius).sqrt();
    // R²/(r(r + |z|)) written as (R/r)(R/(r + |z|)): no cancellation, no R² overflow.
    let geometric = (disk_radius / r) * (disk_radius / (r + abs_z));
    let sign = 1.0_f64.copysign(axial_distance);
    Ok(sign * surface_charge_density / (2.0 * EPSILON_0) * geometric)
}

// ── Field line tracing ─────────────────────────────────────────────

/// Upper bound on the initial allocation of [`trace_field_line`]; the point
/// buffer still grows on demand beyond it.
const TRACE_INITIAL_CAPACITY: usize = 4096;

/// Trace a single field line from a starting point using forward-Euler steps.
///
/// Each step moves a distance `step_size` (m) along the local unit field
/// direction: x_{n+1} = x_n + h·F(x_n)/|F(x_n)|. Returns the visited
/// positions, starting with `start` (at most `max_steps + 1` points).
///
/// `field_fn` computes the field at a given position. Tracing stops early and
/// returns the points collected so far when
/// - the field is exactly zero (a null point), or
/// - the line reaches a source or sink: `field_fn` returns
///   [`BijliError::Singularity`] (as the point-charge functions in this module
///   do) or a field with an infinite component.
///
/// # Errors
///
/// - [`BijliError::InvalidParameter`] if `step_size` is not positive and
///   finite, `start` is not finite, or `field_fn` returns a NaN component.
/// - Any error from `field_fn` other than [`BijliError::Singularity`] is
///   propagated.
pub fn trace_field_line<F>(
    start: [f64; 3],
    step_size: f64,
    max_steps: usize,
    field_fn: F,
) -> Result<Vec<[f64; 3]>>
where
    F: Fn([f64; 3]) -> Result<FieldVector>,
{
    if !(step_size > 0.0 && step_size.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("step size must be positive and finite, got {step_size}"),
        });
    }
    if !start.iter().all(|c| c.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("start position must be finite, got {start:?}"),
        });
    }

    tracing::debug!(?start, step_size, max_steps, "tracing field line");

    let mut points = Vec::with_capacity(max_steps.saturating_add(1).min(TRACE_INITIAL_CAPACITY));
    let mut pos = start;
    points.push(pos);
    let mut stop_reason = "max_steps reached";

    for _ in 0..max_steps {
        let field = match field_fn(pos) {
            Ok(field) => field,
            Err(BijliError::Singularity) => {
                stop_reason = "singularity (source/sink)";
                break;
            }
            Err(e) => return Err(e),
        };
        if field.x.is_nan() || field.y.is_nan() || field.z.is_nan() {
            return Err(BijliError::InvalidParameter {
                reason: format!("field function returned NaN at {pos:?}"),
            });
        }
        if field.x.is_infinite() || field.y.is_infinite() || field.z.is_infinite() {
            stop_reason = "infinite field (source/sink)";
            break;
        }
        // Finite field: normalization fails only for the zero vector.
        let Ok(dir) = field.normalized() else {
            stop_reason = "null point";
            break;
        };
        pos = [
            pos[0] + step_size * dir.x,
            pos[1] + step_size * dir.y,
            pos[2] + step_size * dir.z,
        ];
        points.push(pos);
    }

    tracing::debug!(points = points.len(), stop_reason, "field line traced");
    Ok(points)
}

// ── Multipole expansion ───────────────────────────────────────────

/// Electric quadrupole potential at position r from a quadrupole at the origin.
///
/// ```text
/// V = (1/4πε₀) · ½ Σᵢⱼ Qᵢⱼ xᵢ xⱼ / r⁵        (volts; Jackson eq. 4.10)
/// ```
///
/// `q_tensor` is the traceless quadrupole moment tensor (C⋅m²) in Jackson's
/// convention, Qᵢⱼ = Σₖ qₖ (3xₖᵢ xₖⱼ − rₖ² δᵢⱼ) (or ∫ρ(3xᵢxⱼ − r²δᵢⱼ) d³x), as a
/// 3×3 array. Some texts define Qᵢⱼ with an extra factor ½; double such a
/// tensor before passing it. Example: +q at ±d ẑ and −2q at the origin give
/// Q_zz = 4qd², Q_xx = Q_yy = −2qd², and V = 2qd²P₂(cos θ)/(4πε₀r³).
///
/// # Errors
///
/// [`BijliError::Singularity`] at (numerically) the origin.
#[inline]
pub fn electric_quadrupole_potential(q_tensor: &[[f64; 3]; 3], field_pos: [f64; 3]) -> Result<f64> {
    let r_sq =
        field_pos[0] * field_pos[0] + field_pos[1] * field_pos[1] + field_pos[2] * field_pos[2];
    if r_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let r = r_sq.sqrt();
    let r5 = r_sq * r_sq * r;

    let mut sum = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            sum += q_tensor[i][j] * field_pos[i] * field_pos[j];
        }
    }

    Ok(sum / (8.0 * std::f64::consts::PI * EPSILON_0 * r5))
}

/// Unnormalized real (tesseral) multipole moment of a set of point charges.
///
/// Returns the real solid-harmonic moment (C⋅mˡ)
///
/// ```text
/// M_lm = Σᵢ qᵢ rᵢˡ P_l^|m|(cos θᵢ) × { cos(m φᵢ)     m ≥ 0
///                                   { sin(|m| φᵢ)   m < 0
/// ```
///
/// where P_l^m is the associated Legendre function *including* the
/// Condon–Shortley phase (−1)^m (so P₁¹(x) = −√(1 − x²)), and θᵢ, φᵢ are the
/// polar and azimuthal angles of each charge about the origin. No
/// spherical-harmonic normalization is applied, so this is **not** the
/// complex Jackson moment q_lm = Σ q rˡ Y*_lm; for m ≥ 0 the two are related by
///
/// ```text
/// q_lm = N_lm (M_lm − i M_{l,−m}),   N_lm = √((2l+1)/(4π) · (l−m)!/(l+m)!),
/// q_{l,−m} = (−1)^m q*_lm,           (M_{l,−0} ≡ 0)
/// ```
///
/// Use [`spherical_multipole_moment`] for the normalized complex q_lm. For one
/// charge q at (x, y, z): M₁₀ = qz, M₁₁ = −qx, M₁,₋₁ = −qy, M₂₂ = 3q(x² − y²).
///
/// `l`: multipole order (0 = monopole, 1 = dipole, 2 = quadrupole, ...);
/// `m`: azimuthal index (−l..=l); `charges`: slice of (charge, [x, y, z]).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if |m| > l.
pub fn multipole_moment(l: u32, m: i32, charges: &[(f64, [f64; 3])]) -> Result<f64> {
    if m.unsigned_abs() > l {
        return Err(BijliError::InvalidParameter {
            reason: format!("|m| must be ≤ l, got l={l}, m={m}"),
        });
    }
    let m_abs = m.unsigned_abs();
    let m_f = f64::from(m_abs);

    let mut moment = 0.0;
    for &(q, pos) in charges {
        let Some((r, cos_t, sin_t, phi)) = spherical_coords(pos) else {
            // Charge at the origin: only the monopole (rˡ = 1 for l = 0) survives.
            if l == 0 {
                moment += q;
            }
            continue;
        };
        let azimuthal = if m >= 0 {
            (m_f * phi).cos()
        } else {
            (m_f * phi).sin()
        };
        moment += q * radial_power(r, l) * associated_legendre(l, m_abs, cos_t, sin_t) * azimuthal;
    }

    Ok(moment)
}

/// Normalized complex spherical multipole moment q_lm (Jackson eq. 4.3).
///
/// ```text
/// q_lm = Σᵢ qᵢ rᵢˡ Y*_lm(θᵢ, φᵢ)                                   (C⋅mˡ)
/// Y_lm(θ, φ) = √((2l+1)/(4π) · (l−m)!/(l+m)!) P_l^m(cos θ) e^{imφ}
/// ```
///
/// with orthonormal spherical harmonics including the Condon–Shortley phase
/// (Y_{l,−m} = (−1)^m Y*_lm, hence q_{l,−m} = (−1)^m q*_lm). The potential
/// outside the charges is then
///
/// ```text
/// Φ(r) = (1/ε₀) Σ_lm q_lm Y_lm(θ, φ) / ((2l + 1) r^{l+1})        (Jackson eq. 4.1)
/// ```
///
/// Low orders (Jackson eqs. 4.4–4.6): q₀₀ = Q/√(4π), q₁₀ = √(3/4π) p_z,
/// q₁₁ = −√(3/8π)(p_x − i p_y), q₂₀ = ½√(5/4π) Q₃₃.
///
/// Returns `(Re q_lm, Im q_lm)`. The normalized associated Legendre functions
/// are generated by a stable recurrence, so high orders do not overflow.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if |m| > l.
pub fn spherical_multipole_moment(
    l: u32,
    m: i32,
    charges: &[(f64, [f64; 3])],
) -> Result<(f64, f64)> {
    if m.unsigned_abs() > l {
        return Err(BijliError::InvalidParameter {
            reason: format!("|m| must be ≤ l, got l={l}, m={m}"),
        });
    }
    let m_abs = m.unsigned_abs();
    let m_f = f64::from(m_abs);
    // Y*_{l,−|m|} = (−1)^|m| Y_{l,|m|}
    let neg_sign = if m_abs.is_multiple_of(2) { 1.0 } else { -1.0 };

    let (mut re, mut im) = (0.0, 0.0);
    for &(q, pos) in charges {
        let Some((r, cos_t, sin_t, phi)) = spherical_coords(pos) else {
            if l == 0 {
                re += q * 0.5 / std::f64::consts::PI.sqrt(); // Y₀₀ = 1/√(4π)
            }
            continue;
        };
        let w = q * radial_power(r, l) * normalized_associated_legendre(l, m_abs, cos_t, sin_t);
        let (s, c) = (m_f * phi).sin_cos();
        if m >= 0 {
            // Y*_lm = N P e^{−imφ}
            re += w * c;
            im -= w * s;
        } else {
            re += neg_sign * w * c;
            im += neg_sign * w * s;
        }
    }

    Ok((re, im))
}

/// Spherical coordinates (r, cos θ, sin θ, φ) of a point, or `None` at the
/// origin. cos θ and sin θ are computed directly from the Cartesian
/// components (no acos/cos round trip), so both stay accurate near the poles.
#[inline]
#[must_use]
fn spherical_coords(pos: [f64; 3]) -> Option<(f64, f64, f64, f64)> {
    let rho_sq = pos[0] * pos[0] + pos[1] * pos[1];
    let r = (rho_sq + pos[2] * pos[2]).sqrt();
    if r == 0.0 {
        return None;
    }
    let cos_t = (pos[2] / r).clamp(-1.0, 1.0);
    let sin_t = (rho_sq.sqrt() / r).clamp(0.0, 1.0);
    Some((r, cos_t, sin_t, pos[1].atan2(pos[0])))
}

/// rˡ for any `u32` order.
#[inline]
#[must_use]
fn radial_power(r: f64, l: u32) -> f64 {
    match i32::try_from(l) {
        Ok(li) => r.powi(li),
        Err(_) => r.powf(f64::from(l)),
    }
}

/// Associated Legendre function P_l^m(x) via the standard upward recurrence.
///
/// Physics convention, including the Condon–Shortley phase (−1)^m.
/// `s` = √(1 − x²) = sin θ is passed separately so it can be computed
/// without cancellation near |x| = 1.
#[must_use]
fn associated_legendre(l: u32, m: u32, x: f64, s: f64) -> f64 {
    if m > l {
        return 0.0;
    }

    // P_m^m(x) = (-1)^m (2m-1)!! (1-x²)^{m/2}
    let mut pmm = 1.0;
    let mut fact = 1.0;
    for _ in 0..m {
        pmm *= -fact * s;
        fact += 2.0;
    }

    if l == m {
        return pmm;
    }

    let m_f = f64::from(m);
    // P_{m+1}^m(x) = x(2m+1) P_m^m
    let mut pmmp1 = x * (2.0 * m_f + 1.0) * pmm;
    if l == m + 1 {
        return pmmp1;
    }

    // Recurrence: (l-m) P_l^m = x(2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    let mut pll = 0.0;
    for ll in (m + 2)..=l {
        let ll_f = f64::from(ll);
        pll = (x * (2.0 * ll_f - 1.0) * pmmp1 - (ll_f + m_f - 1.0) * pmm) / (ll_f - m_f);
        pmm = pmmp1;
        pmmp1 = pll;
    }

    pll
}

/// Normalized associated Legendre function
/// P̄_l^m(x) = √((2l+1)/(4π) · (l−m)!/(l+m)!) P_l^m(x), CS phase included,
/// so that Y_lm = P̄_l^m(cos θ) e^{imφ}.
///
/// Stable recurrence (Numerical Recipes, 3rd ed., §6.7): no factorials, no
/// overflow at high order. `s` = sin θ as in [`associated_legendre`].
#[must_use]
fn normalized_associated_legendre(l: u32, m: u32, x: f64, s: f64) -> f64 {
    if m > l {
        return 0.0;
    }
    // P̄_m^m = (−1)^m √((2m+1)/(4π)) √((2m−1)!!/(2m)!!) sᵐ
    let mut pmm = 1.0;
    for k in 1..=m {
        let k_f = f64::from(k);
        pmm *= -s * ((2.0 * k_f - 1.0) / (2.0 * k_f)).sqrt();
    }
    let m_f = f64::from(m);
    pmm *= ((2.0 * m_f + 1.0) / (4.0 * std::f64::consts::PI)).sqrt();
    if l == m {
        return pmm;
    }

    // P̄_{m+1}^m = x √(2m+3) P̄_m^m
    let mut pmmp1 = x * (2.0 * m_f + 3.0).sqrt() * pmm;
    if l == m + 1 {
        return pmmp1;
    }

    // P̄_l^m = a_l [x P̄_{l−1}^m − P̄_{l−2}^m / a_{l−1}],  a_l = √((4l²−1)/(l²−m²))
    let m_sq = m_f * m_f;
    let mut a_prev =
        ((4.0 * (m_f + 1.0) * (m_f + 1.0) - 1.0) / ((m_f + 1.0) * (m_f + 1.0) - m_sq)).sqrt();
    let mut pll = 0.0;
    for ll in (m + 2)..=l {
        let ll_f = f64::from(ll);
        let a = ((4.0 * ll_f * ll_f - 1.0) / (ll_f * ll_f - m_sq)).sqrt();
        pll = a * (x * pmmp1 - pmm / a_prev);
        pmm = pmmp1;
        pmmp1 = pll;
        a_prev = a;
    }
    pll
}

// ── Method of images ──────────────────────────────────────────────

/// Image charge for a point charge near an infinite conducting plane.
///
/// The plane is at z = `plane_z`. Returns (image_charge, image_position).
#[inline]
#[must_use]
pub fn image_charge_plane(charge: f64, position: [f64; 3], plane_z: f64) -> (f64, [f64; 3]) {
    let image_z = 2.0 * plane_z - position[2];
    (-charge, [position[0], position[1], image_z])
}

/// Image charge for a point charge near a grounded conducting sphere.
///
/// Sphere centered at origin with radius R. Returns (image_charge, image_position).
///
/// q' = −(R/d)q at position (R²/d) r̂, where d = |r| > R (Griffiths Ex. 3.2).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `sphere_radius` is not positive and
/// finite, or the charge is not strictly outside the sphere (including a NaN
/// position).
#[inline]
pub fn image_charge_sphere(
    charge: f64,
    position: [f64; 3],
    sphere_radius: f64,
) -> Result<(f64, [f64; 3])> {
    if !(sphere_radius > 0.0 && sphere_radius.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("sphere radius must be positive and finite, got {sphere_radius}"),
        });
    }
    let d_sq = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
    let d = d_sq.sqrt();
    if d.is_nan() || d <= sphere_radius {
        return Err(BijliError::InvalidParameter {
            reason: format!("charge must be outside the sphere (|r| = {d} > R = {sphere_radius})"),
        });
    }

    let ratio = sphere_radius / d;
    let q_image = -charge * ratio;
    let scale = ratio * ratio; // R²/d²
    let pos_image = [
        position[0] * scale,
        position[1] * scale,
        position[2] * scale,
    ];

    Ok((q_image, pos_image))
}

/// Image charge for a point charge near a dielectric half-space (Jackson §4.4).
///
/// The charge q sits at height z > 0 in medium 1 (relative permittivity ε₁,
/// region z > 0); medium 2 (ε₂) fills z < 0; the interface is z = 0.
///
/// Returns (image_charge_for_field_in_medium1, image_position): the image
/// q' = −q(ε₂ − ε₁)/(ε₂ + ε₁) at the mirror point (x, y, −z). The potential in
/// medium 1 is Φ₁ = (1/4πε₀ε₁)(q/R₁ + q'/R₂) — note the ε₁: divide results of
/// the vacuum point-charge functions (e.g. [`electric_potential_point_charge`])
/// by ε₁. (In medium 2 the field is that of q'' = 2ε₂q/(ε₂ + ε₁) at the real
/// charge position, in a uniform medium ε₂.)
///
/// # Errors
///
/// - [`BijliError::InvalidParameter`] if the charge is not strictly above the
///   interface (z ≤ 0, or non-finite) — at z = 0 the image would coincide with
///   the charge.
/// - [`BijliError::InvalidPermittivity`] if ε₁ or ε₂ is not positive and finite.
#[inline]
pub fn image_charge_dielectric(
    charge: f64,
    position: [f64; 3],
    eps_r1: f64,
    eps_r2: f64,
) -> Result<(f64, [f64; 3])> {
    if !(position[2] > 0.0 && position[2].is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "charge must be in medium 1 (z > 0) for dielectric image method, got z = {}",
                position[2]
            ),
        });
    }
    for eps in [eps_r1, eps_r2] {
        if !(eps > 0.0 && eps.is_finite()) {
            return Err(BijliError::InvalidPermittivity { value: eps });
        }
    }
    let q_image = -charge * (eps_r2 - eps_r1) / (eps_r1 + eps_r2);
    let pos_image = [position[0], position[1], -position[2]];
    Ok((q_image, pos_image))
}

// ── Green's functions ─────────────────────────────────────────────

/// Free-space scalar Green's function: G(r,r') = e^{ikR}/(4πR).
///
/// Returns (Re(G), Im(G)) since the result is complex.
/// `k` is the wavenumber (2π/λ).
#[inline]
pub fn greens_function_scalar(
    field_pos: [f64; 3],
    source_pos: [f64; 3],
    k: f64,
) -> Result<(f64, f64)> {
    let dx = field_pos[0] - source_pos[0];
    let dy = field_pos[1] - source_pos[1];
    let dz = field_pos[2] - source_pos[2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    if r < 1e-30 {
        return Err(BijliError::Singularity);
    }

    let factor = 1.0 / (4.0 * std::f64::consts::PI * r);
    let kr = k * r;
    Ok((factor * kr.cos(), factor * kr.sin()))
}

/// Static Green's function (k=0): G(r,r') = 1/(4πR).
#[inline]
pub fn greens_function_static(field_pos: [f64; 3], source_pos: [f64; 3]) -> Result<f64> {
    let dx = field_pos[0] - source_pos[0];
    let dy = field_pos[1] - source_pos[1];
    let dz = field_pos[2] - source_pos[2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    if r < 1e-30 {
        return Err(BijliError::Singularity);
    }

    Ok(1.0 / (4.0 * std::f64::consts::PI * r))
}

// ── Maxwell stress tensor ─────────────────────────────────────────

/// Maxwell stress tensor component T_{ij}.
///
/// T_{ij} = ε₀(E_iE_j − ½δ_{ij}E²) + (1/μ₀)(B_iB_j − ½δ_{ij}B²)
///
/// Returns the 3×3 stress tensor (N/m² = Pa).
#[inline]
#[must_use]
pub fn maxwell_stress_tensor(e: &FieldVector, b: &FieldVector) -> [[f64; 3]; 3] {
    let e_arr = [e.x, e.y, e.z];
    let b_arr = [b.x, b.y, b.z];
    let e_sq = e.magnitude_sq();
    let b_sq = b.magnitude_sq();
    let inv_mu0 = 1.0 / MU_0;

    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            t[i][j] = EPSILON_0 * (e_arr[i] * e_arr[j] - 0.5 * delta * e_sq)
                + inv_mu0 * (b_arr[i] * b_arr[j] - 0.5 * delta * b_sq);
        }
    }
    t
}

/// Electromagnetic force per unit area (radiation pressure) from stress tensor.
///
/// F_i = Σ_j T_{ij} n_j where n is the outward surface normal.
#[inline]
#[must_use]
pub fn stress_tensor_force(tensor: &[[f64; 3]; 3], normal: &FieldVector) -> FieldVector {
    let n = [normal.x, normal.y, normal.z];
    let mut f = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            f[i] += tensor[i][j] * n[j];
        }
    }
    FieldVector::new(f[0], f[1], f[2])
}

/// Total Poynting flux through a closed surface (discrete integration).
///
/// P = ∮ S⃗ · n̂ dA ≈ Σ (E × B / μ₀) · n̂ ΔA
///
/// `samples`: slice of (E, B, outward_normal, area_element) at surface points.
#[must_use]
pub fn poynting_flux_surface(samples: &[(FieldVector, FieldVector, FieldVector, f64)]) -> f64 {
    let inv_mu0 = 1.0 / MU_0;
    let mut total = 0.0;
    for &(ref e, ref b, ref normal, da) in samples {
        let s = e.cross(b).scale(inv_mu0);
        total += s.dot(normal) * da;
    }
    total
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
        assert!((v - 8_987.551_792).abs() < 1.0);
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
        let charges = vec![(1e-6, [-1.0, 0.0, 0.0]), (1e-6, [1.0, 0.0, 0.0])];
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
        )
        .unwrap();
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

    #[test]
    fn test_ops_add() {
        let a = FieldVector::new(1.0, 2.0, 3.0);
        let b = FieldVector::new(4.0, 5.0, 6.0);
        let c = a + b;
        assert_eq!(c, FieldVector::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn test_ops_sub() {
        let a = FieldVector::new(4.0, 5.0, 6.0);
        let b = FieldVector::new(1.0, 2.0, 3.0);
        let c = a - b;
        assert_eq!(c, FieldVector::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_ops_mul() {
        let a = FieldVector::new(1.0, 2.0, 3.0);
        let c = a * 2.0;
        assert_eq!(c, FieldVector::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn test_ops_neg() {
        let a = FieldVector::new(1.0, -2.0, 3.0);
        let c = -a;
        assert_eq!(c, FieldVector::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn test_partial_eq() {
        let a = FieldVector::new(1.0, 2.0, 3.0);
        let b = FieldVector::new(1.0, 2.0, 3.0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_magnitude_sq() {
        let v = FieldVector::new(1.0, 2.0, 3.0);
        assert!((v.magnitude_sq() - 14.0).abs() < 1e-10);
    }

    #[test]
    fn test_ops_add_assign() {
        let mut a = FieldVector::new(1.0, 2.0, 3.0);
        a += FieldVector::new(4.0, 5.0, 6.0);
        assert_eq!(a, FieldVector::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn test_ops_sub_assign() {
        let mut a = FieldVector::new(4.0, 5.0, 6.0);
        a -= FieldVector::new(1.0, 2.0, 3.0);
        assert_eq!(a, FieldVector::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_ops_mul_assign() {
        let mut a = FieldVector::new(1.0, 2.0, 3.0);
        a *= 2.0;
        assert_eq!(a, FieldVector::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn test_ops_scalar_mul_lhs() {
        let a = FieldVector::new(1.0, 2.0, 3.0);
        let c = 2.0 * a;
        assert_eq!(c, FieldVector::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn test_display() {
        let v = FieldVector::new(1.0, 2.5, -3.0);
        let s = format!("{v}");
        assert_eq!(s, "(1, 2.5, -3)");
    }

    #[test]
    fn test_electric_field_singularity() {
        assert!(electric_field_point_charge(1e-6, [0.0; 3], [0.0; 3]).is_err());
    }

    #[test]
    fn test_electric_potential_singularity() {
        assert!(electric_potential_point_charge(1e-6, [0.0; 3], [0.0; 3]).is_err());
    }

    #[test]
    fn test_magnetic_field_singularity() {
        assert!(magnetic_field_moving_charge(1e-6, [0.0; 3], [1e6, 0.0, 0.0], [0.0; 3]).is_err());
    }

    // ── Dipole field tests ─────────────────────────────────────────

    #[test]
    fn test_electric_dipole_on_axis() {
        // Dipole p = p_z ẑ, field point on +z axis at distance d
        // E should be in +z direction (along dipole axis), E_z = 2p/(4πε₀d³)
        let p = FieldVector::new(0.0, 0.0, 1e-30); // tiny dipole moment
        let d = 1.0;
        let e = electric_dipole_field(&p, [0.0, 0.0, d]).unwrap();
        let expected = 2.0 * 1e-30 / (4.0 * std::f64::consts::PI * EPSILON_0 * d * d * d);
        assert!((e.z - expected).abs() / expected < 1e-6);
        assert!(e.x.abs() < 1e-40);
    }

    #[test]
    fn test_electric_dipole_equatorial() {
        // On equatorial plane (perpendicular to dipole axis), field is antiparallel to p
        // E = -p/(4πε₀r³)
        let p = FieldVector::new(0.0, 0.0, 1e-30);
        let d = 1.0;
        let e = electric_dipole_field(&p, [d, 0.0, 0.0]).unwrap();
        let expected = -1e-30 / (4.0 * std::f64::consts::PI * EPSILON_0 * d * d * d);
        assert!((e.z - expected).abs() / expected.abs() < 1e-6);
    }

    #[test]
    fn test_electric_dipole_singularity() {
        let p = FieldVector::new(0.0, 0.0, 1e-30);
        assert!(electric_dipole_field(&p, [0.0; 3]).is_err());
    }

    #[test]
    fn test_electric_dipole_potential_on_axis() {
        let p = FieldVector::new(0.0, 0.0, 1e-30);
        let v = electric_dipole_potential(&p, [0.0, 0.0, 1.0]).unwrap();
        // V = p/(4πε₀r²) on axis
        let expected = 1e-30 / (4.0 * std::f64::consts::PI * EPSILON_0);
        assert!((v - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_electric_dipole_potential_equatorial_zero() {
        // On equatorial plane, p⋅r̂ = 0 → V = 0
        let p = FieldVector::new(0.0, 0.0, 1e-30);
        let v = electric_dipole_potential(&p, [1.0, 0.0, 0.0]).unwrap();
        assert!(v.abs() < 1e-50);
    }

    #[test]
    fn test_magnetic_dipole_on_axis() {
        // Same structure as electric dipole: B_z = 2μ₀m/(4πd³) on axis
        let m = FieldVector::new(0.0, 0.0, 1.0); // 1 A⋅m²
        let d = 1.0;
        let b = magnetic_dipole_field(&m, [0.0, 0.0, d]).unwrap();
        let expected = 2.0 * MU_0 / (4.0 * std::f64::consts::PI * d * d * d);
        assert!((b.z - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_magnetic_dipole_singularity() {
        let m = FieldVector::new(0.0, 0.0, 1.0);
        assert!(magnetic_dipole_field(&m, [0.0; 3]).is_err());
    }

    // ── Gauss's law tests ──────────────────────────────────────────

    #[test]
    fn test_infinite_plane() {
        let sigma = 1e-6; // 1 μC/m²
        let e = electric_field_infinite_plane(sigma);
        assert!((e - sigma / (2.0 * EPSILON_0)).abs() < 1.0);
    }

    #[test]
    fn test_charged_sphere_matches_point_charge() {
        // Outside a sphere, field = point charge field
        let q = 1e-6;
        let r = 2.0;
        let e_sphere = electric_field_charged_sphere(q, r).unwrap();
        let e_point = electric_field_point_charge(q, [0.0; 3], [r, 0.0, 0.0])
            .unwrap()
            .magnitude();
        assert!((e_sphere - e_point).abs() / e_point < 1e-6);
    }

    #[test]
    fn test_charged_sphere_singularity() {
        assert!(electric_field_charged_sphere(1e-6, 0.0).is_err());
    }

    #[test]
    fn test_inside_sphere_zero_at_center() {
        let e = electric_field_inside_sphere(1e-6, 1.0, 0.0).unwrap();
        assert!(e.abs() < 1e-30);
    }

    #[test]
    fn test_inside_sphere_linear_growth() {
        let q = 1e-6;
        let big_r = 1.0;
        let e1 = electric_field_inside_sphere(q, big_r, 0.25).unwrap();
        let e2 = electric_field_inside_sphere(q, big_r, 0.50).unwrap();
        // E grows linearly with r → ratio should be 2.0
        assert!((e2 / e1 - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_infinite_line() {
        let lambda = 1e-6; // 1 μC/m
        let r = 1.0;
        let e = electric_field_infinite_line(lambda, r).unwrap();
        let expected = lambda / (2.0 * std::f64::consts::PI * EPSILON_0 * r);
        assert!((e - expected).abs() / expected < 1e-6);
    }

    #[test]
    fn test_infinite_line_singularity() {
        assert!(electric_field_infinite_line(1e-6, 0.0).is_err());
    }

    #[test]
    fn test_outside_cylinder_equals_line() {
        let lambda = 1e-6;
        let r = 2.0;
        let e_cyl = electric_field_outside_cylinder(lambda, r).unwrap();
        let e_line = electric_field_infinite_line(lambda, r).unwrap();
        assert!((e_cyl - e_line).abs() < 1e-30);
    }

    #[test]
    fn test_inside_cylinder() {
        let rho = 1e-6; // 1 μC/m³
        let r = 0.5;
        let e = electric_field_inside_cylinder(rho, r).unwrap();
        let expected = rho * r / (2.0 * EPSILON_0);
        assert!((e - expected).abs() / expected < 1e-6);
    }

    // ── Charge distribution tests ──────────────────────────────────

    #[test]
    fn test_ring_axis_center() {
        // At center of ring (z=0), field is zero by symmetry
        let e = electric_field_ring_axis(1e-6, 1.0, 0.0).unwrap();
        assert!(e.abs() < 1e-30);
    }

    #[test]
    fn test_ring_axis_far_field() {
        // Far from ring, should approximate point charge
        let q = 1e-6;
        let big_r = 0.01; // tiny ring
        let z = 10.0; // far away
        let e_ring = electric_field_ring_axis(q, big_r, z).unwrap();
        let e_point = COULOMB_K * q / (z * z);
        assert!((e_ring - e_point).abs() / e_point < 1e-4);
    }

    #[test]
    fn test_disk_axis_far_field_approaches_point() {
        // Very far from disk → approaches point charge: E ≈ Q/(4πε₀z²)
        let sigma = 1e-6;
        let big_r = 0.01;
        let z = 100.0;
        let e_disk = electric_field_disk_axis(sigma, big_r, z).unwrap();
        // Total charge Q = σπR²
        let q = sigma * std::f64::consts::PI * big_r * big_r;
        let e_point = COULOMB_K * q / (z * z);
        assert!((e_disk - e_point).abs() / e_point < 0.01);
    }

    #[test]
    fn test_disk_axis_close_approaches_infinite_plane() {
        // Very close to large disk → approaches σ/(2ε₀)
        let sigma = 1e-6;
        let big_r = 1000.0; // huge disk
        let z = 0.001; // very close
        let e_disk = electric_field_disk_axis(sigma, big_r, z).unwrap();
        let e_plane = electric_field_infinite_plane(sigma);
        assert!((e_disk - e_plane).abs() / e_plane < 0.01);
    }

    // ── Field line tracing tests ───────────────────────────────────

    #[test]
    fn test_trace_field_line_uniform() {
        // Uniform field in +x: line should be straight
        let points = trace_field_line([0.0, 0.0, 0.0], 0.1, 10, |_pos| {
            Ok(FieldVector::new(1.0, 0.0, 0.0))
        })
        .unwrap();
        assert_eq!(points.len(), 11);
        assert!((points[10][0] - 1.0).abs() < 1e-10);
        assert!(points[10][1].abs() < 1e-10);
    }

    #[test]
    fn test_trace_field_line_from_charge() {
        // Trace radially outward from a positive charge
        let q = 1e-6;
        let points = trace_field_line([0.1, 0.0, 0.0], 0.01, 5, |pos| {
            electric_field_point_charge(q, [0.0; 3], pos)
        })
        .unwrap();
        // Should move in +x direction
        assert!(points.last().unwrap()[0] > points[0][0]);
    }

    #[test]
    fn test_trace_field_line_invalid_step() {
        assert!(
            trace_field_line([0.0; 3], -1.0, 10, |_| Ok(FieldVector::new(1.0, 0.0, 0.0)),).is_err()
        );
    }

    // ── Trait impl tests ──────────────────────────────────────────

    #[test]
    fn test_default_is_zero() {
        let v = FieldVector::default();
        assert_eq!(v, FieldVector::zero());
    }

    #[test]
    fn test_from_array() {
        let v = FieldVector::from([1.0, 2.0, 3.0]);
        assert_eq!(v, FieldVector::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_into_array() {
        let v = FieldVector::new(1.0, 2.0, 3.0);
        let arr: [f64; 3] = v.into();
        assert_eq!(arr, [1.0, 2.0, 3.0]);
    }

    // ── Multipole expansion tests ─────────────────────────────────

    #[test]
    fn test_quadrupole_potential_singularity() {
        let q = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        assert!(electric_quadrupole_potential(&q, [0.0; 3]).is_err());
    }

    #[test]
    fn test_quadrupole_potential_decays_faster_than_dipole() {
        // Quadrupole ∝ 1/r³, dipole ∝ 1/r²
        let q = [[1e-30, 0.0, 0.0], [0.0, -1e-30, 0.0], [0.0, 0.0, 0.0]];
        let v1 = electric_quadrupole_potential(&q, [1.0, 0.0, 0.0])
            .unwrap()
            .abs();
        let v2 = electric_quadrupole_potential(&q, [2.0, 0.0, 0.0])
            .unwrap()
            .abs();
        // V ∝ 1/r³ → ratio ≈ 8
        let ratio = v1 / v2;
        assert!((ratio - 8.0).abs() < 0.1);
    }

    #[test]
    fn test_multipole_moment_monopole() {
        // l=0: monopole moment = total charge
        let charges = vec![(1e-6, [0.0, 0.0, 0.0]), (2e-6, [1.0, 0.0, 0.0])];
        let q00 = multipole_moment(0, 0, &charges).unwrap();
        // P_0^0(x) = 1, Y_00 ∝ 1, so monopole ∝ Σq
        assert!(q00 > 0.0);
    }

    #[test]
    fn test_multipole_moment_invalid_m() {
        assert!(multipole_moment(1, 2, &[]).is_err());
    }

    #[test]
    fn test_associated_legendre_basic() {
        let s = (1.0 - 0.25_f64).sqrt();
        // P_0^0 = 1
        assert!((associated_legendre(0, 0, 0.5, s) - 1.0).abs() < 1e-10);
        // P_1^0 = x
        assert!((associated_legendre(1, 0, 0.5, s) - 0.5).abs() < 1e-10);
        // P_1^1 = -√(1-x²)
        let p11 = associated_legendre(1, 1, 0.5, s);
        assert!((p11 - (-(1.0 - 0.25_f64).sqrt())).abs() < 1e-10);
    }

    // ── Method of images tests ────────────────────────────────────

    #[test]
    fn test_image_charge_plane() {
        let (q_img, pos_img) = image_charge_plane(1e-6, [0.0, 0.0, 1.0], 0.0);
        assert!((q_img - (-1e-6)).abs() < 1e-20);
        assert!((pos_img[2] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_image_charge_plane_tangential_zero() {
        // On the conducting plane (z=0), tangential E (x,y) should be zero
        let q = 1e-6;
        let pos = [0.0, 0.0, 1.0];
        let (q_img, pos_img) = image_charge_plane(q, pos, 0.0);
        let point_on_plane = [1.0, 0.0, 0.0];
        let e_real = electric_field_point_charge(q, pos, point_on_plane).unwrap();
        let e_img = electric_field_point_charge(q_img, pos_img, point_on_plane).unwrap();
        // Tangential components (x,y) should cancel by symmetry
        assert!((e_real.x + e_img.x).abs() < 1e-10);
        assert!((e_real.y + e_img.y).abs() < 1e-10);
    }

    #[test]
    fn test_image_charge_sphere() {
        let (q_img, pos_img) = image_charge_sphere(1e-6, [0.0, 0.0, 2.0], 1.0).unwrap();
        // q' = -R/d * q = -0.5e-6
        assert!((q_img - (-0.5e-6)).abs() < 1e-20);
        // pos = R²/d² * pos = 0.25 * [0,0,2] = [0,0,0.5]
        assert!((pos_img[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_image_charge_sphere_inside_fails() {
        assert!(image_charge_sphere(1e-6, [0.0, 0.0, 0.5], 1.0).is_err());
    }

    #[test]
    fn test_image_charge_dielectric() {
        let (q_img, pos_img) = image_charge_dielectric(1e-6, [0.0, 0.0, 1.0], 1.0, 4.0).unwrap();
        // q' = -q(ε₂-ε₁)/(ε₂+ε₁) = -q(3/5) = -0.6e-6
        assert!((q_img - (-0.6e-6)).abs() < 1e-12);
        assert!((pos_img[2] - (-1.0)).abs() < 1e-10);
    }

    // ── Green's function tests ────────────────────────────────────

    #[test]
    fn test_greens_static_coulomb() {
        // Static Green's function at r=1 should be 1/(4π)
        let g = greens_function_static([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
        assert!((g - 1.0 / (4.0 * std::f64::consts::PI)).abs() < 1e-10);
    }

    #[test]
    fn test_greens_static_singularity() {
        assert!(greens_function_static([0.0; 3], [0.0; 3]).is_err());
    }

    #[test]
    fn test_greens_scalar_k0_matches_static() {
        let (re, im) = greens_function_scalar([2.0, 0.0, 0.0], [0.0; 3], 0.0).unwrap();
        let g_static = greens_function_static([2.0, 0.0, 0.0], [0.0; 3]).unwrap();
        assert!((re - g_static).abs() < 1e-10);
        assert!(im.abs() < 1e-10);
    }

    // ── Stress tensor tests ───────────────────────────────────────

    #[test]
    fn test_stress_tensor_pure_e_field() {
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let t = maxwell_stress_tensor(&e, &b);
        // T_xx = ε₀(E_x² - E²/2) = ε₀ E²/2
        let expected = 0.5 * EPSILON_0 * 1e6;
        assert!((t[0][0] - expected).abs() / expected < 1e-6);
        // T_yy = T_zz = -ε₀ E²/2
        assert!((t[1][1] - (-expected)).abs() / expected < 1e-6);
    }

    #[test]
    fn test_stress_tensor_force_on_surface() {
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let t = maxwell_stress_tensor(&e, &b);
        let normal = FieldVector::new(1.0, 0.0, 0.0); // surface facing +x
        let force = stress_tensor_force(&t, &normal);
        // Force should be in +x (tension along field lines)
        assert!(force.x > 0.0);
    }

    #[test]
    fn test_poynting_flux_plane_wave() {
        // Plane wave: S = E × B / μ₀, all pointing in same direction
        let e = FieldVector::new(0.0, 100.0, 0.0); // E in y
        let b = FieldVector::new(0.0, 0.0, 100.0 / SPEED_OF_LIGHT); // B in z
        let normal = FieldVector::new(1.0, 0.0, 0.0); // surface normal +x
        let da = 0.01; // 1 cm²
        let flux = poynting_flux_surface(&[(e, b, normal, da)]);
        assert!(flux > 0.0);
    }

    // ── Audit-driven V1.5 edge case tests ─────────────────────────

    #[test]
    fn test_stress_tensor_symmetric() {
        // Maxwell stress tensor should be symmetric: T_{ij} = T_{ji}
        let e = FieldVector::new(100.0, 200.0, 300.0);
        let b = FieldVector::new(1e-4, 2e-4, 3e-4);
        let t = maxwell_stress_tensor(&e, &b);
        // Check off-diagonal symmetry
        assert!((t[0][1] - t[1][0]).abs() < 1e-15);
        assert!((t[0][2] - t[2][0]).abs() < 1e-15);
        assert!((t[1][2] - t[2][1]).abs() < 1e-15);
    }

    #[test]
    fn test_stress_tensor_trace() {
        // Trace = ε₀E²/2 + B²/(2μ₀) - ε₀E²/2 - B²/(2μ₀) ... no,
        // Tr(T) = ε₀(E²-3E²/2) + (B²-3B²/2)/μ₀ = -ε₀E²/2 - B²/(2μ₀)
        let e = FieldVector::new(1000.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let t = maxwell_stress_tensor(&e, &b);
        let trace = t[0][0] + t[1][1] + t[2][2];
        let expected = -0.5 * EPSILON_0 * e.magnitude_sq();
        assert!((trace - expected).abs() / expected.abs() < 1e-6);
    }

    #[test]
    fn test_poynting_flux_empty_surface() {
        let flux = poynting_flux_surface(&[]);
        assert!(flux.abs() < 1e-30);
    }

    #[test]
    fn test_greens_function_inverse_r() {
        // G should fall off as 1/r
        let g1 = greens_function_static([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
        let g2 = greens_function_static([2.0, 0.0, 0.0], [0.0; 3]).unwrap();
        assert!((g1 / g2 - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_greens_scalar_oscillates() {
        // For k > 0, the imaginary part should oscillate
        let k = 100.0;
        let (_, im1) = greens_function_scalar([0.01, 0.0, 0.0], [0.0; 3], k).unwrap();
        let (_, im2) = greens_function_scalar([0.05, 0.0, 0.0], [0.0; 3], k).unwrap();
        // Different signs possible (oscillation)
        assert!(im1.abs() > 0.0);
        assert!(im2.abs() > 0.0);
    }

    #[test]
    fn test_image_charge_dielectric_wrong_halfspace() {
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, -1.0], 1.0, 4.0).is_err());
    }

    #[test]
    fn test_image_charge_dielectric_equal_media() {
        // ε₁ = ε₂ → no image charge (q' = 0)
        let (q_img, _) = image_charge_dielectric(1e-6, [0.0, 0.0, 1.0], 2.0, 2.0).unwrap();
        assert!(q_img.abs() < 1e-20);
    }

    #[test]
    fn test_image_charge_sphere_on_surface() {
        // Charge exactly on sphere surface → error (not outside)
        assert!(image_charge_sphere(1e-6, [0.0, 0.0, 1.0], 1.0).is_err());
    }

    #[test]
    fn test_associated_legendre_p20() {
        // P_2^0(x) = (3x²-1)/2
        let x = 0.5;
        let p = associated_legendre(2, 0, x, (1.0 - x * x).sqrt());
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((p - expected).abs() < 1e-10);
    }

    #[test]
    fn test_multipole_moment_empty_charges() {
        let m = multipole_moment(1, 0, &[]).unwrap();
        assert!(m.abs() < 1e-30);
    }

    // ── Audit regression tests (1.2.0) ────────────────────────────

    /// Direct Coulomb sum used as an independent reference.
    fn direct_potential(charges: &[(f64, [f64; 3])], p: [f64; 3]) -> f64 {
        charges
            .iter()
            .map(|&(q, c)| {
                let d =
                    ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt();
                COULOMB_K * q / d
            })
            .sum()
    }

    #[test]
    fn test_quadrupole_potential_linear_quadrupole_far_field() {
        // +q at ±d ẑ, −2q at origin. Jackson 4.10: Q_ij = Σ q(3xᵢxⱼ − r²δᵢⱼ)
        // → Q_zz = 4qd², Q_xx = Q_yy = −2qd²; V → 2kqd² P₂(cosθ)/r³.
        let q = 1e-9;
        let d = 1e-3;
        let charges = [
            (q, [0.0, 0.0, d]),
            (q, [0.0, 0.0, -d]),
            (-2.0 * q, [0.0; 3]),
        ];
        let qd2 = q * d * d;
        let q_tensor = [
            [-2.0 * qd2, 0.0, 0.0],
            [0.0, -2.0 * qd2, 0.0],
            [0.0, 0.0, 4.0 * qd2],
        ];
        for p in [
            [0.0, 0.0, 1.0],
            [0.3, -0.4, 0.5],
            [1.0, 0.0, 0.0],
            [0.2, 0.7, -0.9],
        ] {
            let v = electric_quadrupole_potential(&q_tensor, p).unwrap();
            let v_ref = direct_potential(&charges, p);
            // Next multipole (l = 4) is O(d²/r²) ≈ 1e-6 relative.
            assert!(
                (v - v_ref).abs() / v_ref.abs() < 1e-5,
                "at {p:?}: quadrupole {v} vs direct {v_ref}"
            );
        }
    }

    #[test]
    fn test_disk_axis_below_disk_antisymmetric() {
        let sigma = 1e-6;
        let big_r = 1.0;
        let z = 0.5;
        let above = electric_field_disk_axis(sigma, big_r, z).unwrap();
        let below = electric_field_disk_axis(sigma, big_r, -z).unwrap();
        let expected = sigma / (2.0 * EPSILON_0) * (1.0 - z / (z * z + big_r * big_r).sqrt());
        assert!((above - expected).abs() / expected < 1e-12);
        assert!(
            (below + expected).abs() / expected < 1e-12,
            "below = {below}"
        );
    }

    #[test]
    fn test_disk_axis_far_field_no_cancellation() {
        // z/R = 1e8: E_z = (σ/2ε₀)·R²/(r(r+z)) ≈ kσπR²/z² (1 − 3R²/(4z²)…)
        let sigma = 1e-6;
        let big_r = 1.0;
        let z = 1e8;
        let e = electric_field_disk_axis(sigma, big_r, z).unwrap();
        let e_point = COULOMB_K * sigma * std::f64::consts::PI * big_r * big_r / (z * z);
        assert!((e - e_point).abs() / e_point < 1e-12, "{e} vs {e_point}");
    }

    #[test]
    fn test_disk_axis_on_surface_one_sided_limits() {
        let sigma = 1e-6;
        let plane = sigma / (2.0 * EPSILON_0);
        let up = electric_field_disk_axis(sigma, 1.0, 0.0).unwrap();
        let down = electric_field_disk_axis(sigma, 1.0, -0.0).unwrap();
        assert!((up - plane).abs() / plane < 1e-15);
        assert!((down + plane).abs() / plane < 1e-15);
    }

    #[test]
    fn test_ring_axis_atomic_scale() {
        // Bohr-radius ring: (z² + R²)^{3/2} ≈ 1.6e-31 must not be a singularity.
        let q = -1.602_176_634e-19;
        let big_r = 5.29e-11;
        let z = 1e-11;
        let e = electric_field_ring_axis(q, big_r, z).unwrap();
        let expected = COULOMB_K * q * z / (z * z + big_r * big_r).powf(1.5);
        assert!((e - expected).abs() / expected.abs() < 1e-12);
    }

    #[test]
    fn test_ring_axis_true_singularity() {
        assert!(matches!(
            electric_field_ring_axis(1e-9, 0.0, 0.0),
            Err(BijliError::Singularity)
        ));
        assert!(electric_field_ring_axis(1e-9, f64::NAN, 1.0).is_err());
    }

    #[test]
    fn test_inside_sphere_rejects_outside_domain() {
        assert!(electric_field_inside_sphere(1e-6, 1.0, 2.0).is_err());
        assert!(electric_field_inside_sphere(1e-6, 1.0, -0.5).is_err());
        assert!(electric_field_inside_sphere(1e-6, 1.0, f64::NAN).is_err());
        assert!(electric_field_inside_sphere(1e-6, f64::NAN, 0.5).is_err());
        // Surface value matches the outside solution (continuity).
        let at_surface = electric_field_inside_sphere(1e-6, 1.0, 1.0).unwrap();
        let outside = electric_field_charged_sphere(1e-6, 1.0).unwrap();
        assert!((at_surface - outside).abs() / outside < 1e-15);
    }

    #[test]
    fn test_image_charge_sphere_rejects_bad_radius() {
        assert!(image_charge_sphere(1e-6, [0.0, 0.0, 2.0], -1.0).is_err());
        assert!(image_charge_sphere(1e-6, [0.0, 0.0, 2.0], 0.0).is_err());
        assert!(image_charge_sphere(1e-6, [0.0, 0.0, 2.0], f64::NAN).is_err());
        assert!(image_charge_sphere(1e-6, [f64::NAN, 0.0, 2.0], 1.0).is_err());
    }

    #[test]
    fn test_image_charge_sphere_grounded_surface_potential_zero() {
        // Grounded sphere: V(real) + V(image) = 0 everywhere on |r| = R.
        let (q, pos, big_r) = (1e-9, [0.3, -1.1, 2.0], 0.7);
        let (qi, pi) = image_charge_sphere(q, pos, big_r).unwrap();
        for p in [
            [big_r, 0.0, 0.0],
            [0.0, 0.0, -big_r],
            [0.4, 0.4, big_r * 0.5],
        ] {
            let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let s = [p[0] * big_r / n, p[1] * big_r / n, p[2] * big_r / n];
            let v = direct_potential(&[(q, pos), (qi, pi)], s);
            let scale = direct_potential(&[(q, pos)], s);
            assert!(v.abs() / scale < 1e-12);
        }
    }

    #[test]
    fn test_image_charge_dielectric_rejects_invalid() {
        // Charge on the interface: image would coincide with the charge.
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, 0.0], 1.0, 4.0).is_err());
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, f64::NAN], 1.0, 4.0).is_err());
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, 1.0], 2.0, -1.0).is_err());
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, 1.0], 0.0, 4.0).is_err());
        assert!(image_charge_dielectric(1e-6, [0.0, 0.0, 1.0], 1.0, f64::NAN).is_err());
    }

    #[test]
    fn test_normalized_small_si_magnitudes() {
        // 1.6e-19 N (force on e in 1 V/m) and 6e-30 C·m (water dipole).
        for s in [1.6e-19, 6.2e-30, 1e-200, 1e-300] {
            let v = FieldVector::new(3.0 * s, 0.0, 4.0 * s);
            let n = v.normalized().unwrap();
            assert!((n.x - 0.6).abs() < 1e-15, "scale {s}: {n}");
            assert!((n.z - 0.8).abs() < 1e-15, "scale {s}: {n}");
        }
    }

    #[test]
    fn test_normalized_huge_components() {
        // |v|² overflows to ∞; the direction is still well defined.
        let v = FieldVector::new(1e200, 1e200, 0.0);
        let n = v.normalized().unwrap();
        let h = std::f64::consts::FRAC_1_SQRT_2;
        assert!((n.x - h).abs() < 1e-15 && (n.y - h).abs() < 1e-15, "{n}");
    }

    #[test]
    fn test_normalized_non_finite_fails() {
        assert!(FieldVector::new(f64::NAN, 1.0, 0.0).normalized().is_err());
        assert!(
            FieldVector::new(f64::INFINITY, 1.0, 0.0)
                .normalized()
                .is_err()
        );
    }

    #[test]
    fn test_trace_field_line_huge_max_steps() {
        // Must not try to reserve usize::MAX + 1 elements.
        let points =
            trace_field_line([0.0; 3], 0.1, usize::MAX, |_| Ok(FieldVector::zero())).unwrap();
        assert_eq!(points.len(), 1);
    }

    #[test]
    fn test_trace_field_line_nan_step_rejected() {
        let f = |_: [f64; 3]| Ok(FieldVector::new(1.0, 0.0, 0.0));
        assert!(trace_field_line([0.0; 3], f64::NAN, 10, f).is_err());
        assert!(trace_field_line([0.0; 3], f64::INFINITY, 10, f).is_err());
        assert!(trace_field_line([f64::NAN, 0.0, 0.0], 0.1, 10, f).is_err());
    }

    #[test]
    fn test_trace_field_line_stops_at_sink() {
        // Field function reports a singularity once x ≥ 0.25: the line ends there.
        let points = trace_field_line([0.0; 3], 0.1, 100, |p| {
            if p[0] >= 0.25 {
                Err(BijliError::Singularity)
            } else {
                Ok(FieldVector::new(1.0, 0.0, 0.0))
            }
        })
        .unwrap();
        assert_eq!(points.len(), 4); // x = 0, 0.1, 0.2, 0.3
        assert!((points[3][0] - 0.3).abs() < 1e-12);
        // Other errors still propagate.
        assert!(trace_field_line([0.0; 3], 0.1, 10, |_| Err(BijliError::ZeroCharge)).is_err());
    }

    #[test]
    fn test_trace_field_line_tiny_field_still_traced() {
        // Physical fields below 1e-30 (in SI) still have a direction.
        let points =
            trace_field_line([0.0; 3], 0.5, 2, |_| Ok(FieldVector::new(0.0, 1e-35, 0.0))).unwrap();
        assert_eq!(points.len(), 3);
        assert!((points[2][1] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_gauss_law_rejects_negative_distance() {
        assert!(electric_field_charged_sphere(1e-6, -2.0).is_err());
        assert!(electric_field_infinite_line(1e-6, -1.0).is_err());
        assert!(electric_field_outside_cylinder(1e-6, -1.0).is_err());
        assert!(electric_field_charged_sphere(1e-6, f64::NAN).is_err());
        assert!(electric_field_infinite_line(1e-6, f64::NAN).is_err());
        assert!(electric_field_inside_cylinder(1e-6, f64::NAN).is_err());
    }

    #[test]
    fn test_constants_codata_2018_consistent() {
        // ε₀μ₀c² = 1 and k = 1/(4πε₀) to CODATA 2018 precision.
        let one = EPSILON_0 * MU_0 * SPEED_OF_LIGHT_SQ;
        assert!((one - 1.0).abs() < 1e-12, "ε₀μ₀c² − 1 = {}", one - 1.0);
        let k = 1.0 / (4.0 * std::f64::consts::PI * EPSILON_0);
        assert!((COULOMB_K - k).abs() / k < 1e-15);
        assert!((EPSILON_0 - 8.854_187_812_8e-12).abs() < 1e-22);
        assert!((COULOMB_K - 8.987_551_792_3e9).abs() / COULOMB_K < 1e-10);
    }

    #[test]
    fn test_multipole_moment_near_pole_precision() {
        // Charge at (1e-9, 0, 1): cosθ rounds to 1, so 1 − cos²θ = 0 loses sinθ.
        // M_11 = q r P₁¹(cosθ) cos φ = −q x exactly.
        let q = 2e-9;
        let x = 1e-9;
        let m11 = multipole_moment(1, 1, &[(q, [x, 0.0, 1.0])]).unwrap();
        assert!((m11 - (-q * x)).abs() / (q * x) < 1e-12, "{m11}");
    }

    #[test]
    fn test_multipole_moment_documented_definition() {
        // Unnormalized real solid harmonics with Condon–Shortley phase:
        //   l=1: m=0 → q z, m=1 → −q x, m=−1 → −q y
        //   l=2: m=0 → q(3z² − r²)/2, m=2 → 3q(x² − y²), m=−2 → 6q x y, m=1 → −3q x z
        let q = 1.5e-9;
        let (x, y, z) = (0.3, -0.7, 0.4);
        let r2 = x * x + y * y + z * z;
        let c = [(q, [x, y, z])];
        let cases = [
            (1, 0, q * z),
            (1, 1, -q * x),
            (1, -1, -q * y),
            (2, 0, q * (3.0 * z * z - r2) / 2.0),
            (2, 1, -3.0 * q * x * z),
            (2, -1, -3.0 * q * y * z),
            (2, 2, 3.0 * q * (x * x - y * y)),
            (2, -2, 6.0 * q * x * y),
        ];
        for (l, m, expected) in cases {
            let got = multipole_moment(l, m, &c).unwrap();
            assert!(
                (got - expected).abs() < 1e-12 * q,
                "l={l} m={m}: {got} vs {expected}"
            );
        }
    }

    #[test]
    fn test_spherical_multipole_moment_jackson_low_orders() {
        // Jackson 4.4–4.6 for one charge q at (x, y, z):
        //   q₀₀ = q/√(4π), q₁₀ = √(3/4π) qz, q₁₁ = −√(3/8π) q(x − iy),
        //   q₂₀ = ½√(5/4π) q(3z² − r²), q₂₁ = −√(15/8π) qz(x − iy),
        //   q₂₂ = ¼√(15/2π) q(x − iy)²,  q_{l,−m} = (−1)^m q*_lm.
        use std::f64::consts::PI;
        let q = 2.5e-9;
        let (x, y, z) = (0.3, -0.7, 0.4);
        let r2 = x * x + y * y + z * z;
        let c = [(q, [x, y, z])];
        let (xm_re, xm_im) = (x, -y); // x − iy
        let sq = (xm_re * xm_re - xm_im * xm_im, 2.0 * xm_re * xm_im); // (x − iy)²
        let cases = [
            (0, 0, (q / (4.0 * PI).sqrt(), 0.0)),
            (1, 0, ((3.0 / (4.0 * PI)).sqrt() * q * z, 0.0)),
            (1, 1, {
                let k = -(3.0 / (8.0 * PI)).sqrt() * q;
                (k * xm_re, k * xm_im)
            }),
            (1, -1, {
                let k = (3.0 / (8.0 * PI)).sqrt() * q;
                (k * xm_re, -k * xm_im)
            }),
            (
                2,
                0,
                (
                    0.5 * (5.0 / (4.0 * PI)).sqrt() * q * (3.0 * z * z - r2),
                    0.0,
                ),
            ),
            (2, 1, {
                let k = -(15.0 / (8.0 * PI)).sqrt() * q * z;
                (k * xm_re, k * xm_im)
            }),
            (2, -2, {
                let k = 0.25 * (15.0 / (2.0 * PI)).sqrt() * q;
                (k * sq.0, -k * sq.1)
            }),
            (2, 2, {
                let k = 0.25 * (15.0 / (2.0 * PI)).sqrt() * q;
                (k * sq.0, k * sq.1)
            }),
        ];
        for (l, m, (re_ref, im_ref)) in cases {
            let (re, im) = spherical_multipole_moment(l, m, &c).unwrap();
            assert!(
                (re - re_ref).abs() < 1e-14 * q && (im - im_ref).abs() < 1e-14 * q,
                "l={l} m={m}: ({re}, {im}) vs ({re_ref}, {im_ref})"
            );
        }
    }

    #[test]
    fn test_spherical_multipole_moment_relation_to_real_moment() {
        // q_lm = N_lm (M_lm − i M_{l,−m}) for m > 0.
        let charges = [
            (1e-9, [0.2, 0.5, -0.3]),
            (-3e-9, [-0.4, 0.1, 0.6]),
            (2e-9, [0.0, -0.8, 0.1]),
        ];
        let (l, m) = (5_u32, 3_i32);
        // N_53 = √(11/(4π) · 2!/8!)
        let n = (11.0 / (4.0 * std::f64::consts::PI) * 2.0 / 40_320.0).sqrt();
        let m_c = multipole_moment(l, m, &charges).unwrap();
        let m_s = multipole_moment(l, -m, &charges).unwrap();
        let (re, im) = spherical_multipole_moment(l, m, &charges).unwrap();
        let scale = m_c.abs().max(m_s.abs()) * n;
        assert!((re - n * m_c).abs() < 1e-13 * scale);
        assert!((im + n * m_s).abs() < 1e-13 * scale);
    }

    #[test]
    fn test_spherical_multipole_moment_high_order_no_overflow() {
        // Unnormalized P_200^200 = 399!! overflows f64; the normalized
        // recurrence must not. Y_ll(π/2, 0) = √((2l+1)!/(4π))/(2ˡ l!) for even l
        // = 1.127_941_485_749_480_3 (log-gamma reference) at l = 200.
        let (re, im) = spherical_multipole_moment(200, 200, &[(1.0, [1.0, 0.0, 0.0])]).unwrap();
        assert!((re - 1.127_941_485_749_480_3).abs() < 1e-12, "{re}");
        assert!(im.abs() < 1e-15);
        // m = 0 on the unit sphere at cos θ = 0.3: √(401/4π) P₂₀₀(0.3)
        // = −0.055_131_838_971_581_58 (Bonnet recurrence reference).
        let s = (1.0 - 0.09_f64).sqrt();
        let (re, _) = spherical_multipole_moment(200, 0, &[(1.0, [s, 0.0, 0.3])]).unwrap();
        assert!((re - (-0.055_131_838_971_581_58)).abs() < 1e-12, "{re}");
    }

    #[test]
    fn test_spherical_multipole_moment_invalid_and_origin() {
        assert!(spherical_multipole_moment(1, -2, &[]).is_err());
        // A charge at the origin contributes only to the monopole.
        let c = [(4.0, [0.0; 3])];
        let (re, im) = spherical_multipole_moment(0, 0, &c).unwrap();
        assert!((re - 4.0 / (4.0 * std::f64::consts::PI).sqrt()).abs() < 1e-15 && im == 0.0);
        assert_eq!(spherical_multipole_moment(2, 1, &c).unwrap(), (0.0, 0.0));
        assert!((multipole_moment(0, 0, &c).unwrap() - 4.0).abs() < 1e-15);
    }
}
