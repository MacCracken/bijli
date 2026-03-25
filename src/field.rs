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
    #[inline]
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

/// Electric field outside a uniformly charged sphere (r > R).
///
/// E = Q/(4πε₀r²), same as a point charge.
#[inline]
pub fn electric_field_charged_sphere(total_charge: f64, distance: f64) -> Result<f64> {
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(COULOMB_K * total_charge / (distance * distance))
}

/// Electric field inside a uniformly charged sphere (r < R).
///
/// E = Qr/(4πε₀R³) — grows linearly with r.
#[inline]
pub fn electric_field_inside_sphere(
    total_charge: f64,
    sphere_radius: f64,
    distance: f64,
) -> Result<f64> {
    if sphere_radius <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("sphere radius must be positive, got {sphere_radius}"),
        });
    }
    Ok(COULOMB_K * total_charge * distance / (sphere_radius * sphere_radius * sphere_radius))
}

/// Electric field from an infinite line of uniform linear charge density.
///
/// E = λ/(2πε₀r), directed radially outward.
#[inline]
pub fn electric_field_infinite_line(linear_charge_density: f64, distance: f64) -> Result<f64> {
    if distance.abs() < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(linear_charge_density / (2.0 * std::f64::consts::PI * EPSILON_0 * distance))
}

/// Electric field outside an infinite uniformly charged cylinder (r > R).
///
/// E = λ/(2πε₀r) where λ is the charge per unit length.
/// Same form as infinite line charge.
#[inline]
pub fn electric_field_outside_cylinder(linear_charge_density: f64, distance: f64) -> Result<f64> {
    electric_field_infinite_line(linear_charge_density, distance)
}

/// Electric field inside an infinite uniformly charged cylinder (r < R).
///
/// E = ρr/(2ε₀) where ρ is the volume charge density.
#[inline]
pub fn electric_field_inside_cylinder(volume_charge_density: f64, distance: f64) -> Result<f64> {
    if distance < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("distance must be non-negative, got {distance}"),
        });
    }
    Ok(volume_charge_density * distance / (2.0 * EPSILON_0))
}

// ── Charge distributions ───────────────────────────────────────────

/// Electric field on the axis of a uniformly charged ring.
///
/// E_z = Qz / (4πε₀(z² + R²)^(3/2))
#[inline]
pub fn electric_field_ring_axis(
    total_charge: f64,
    ring_radius: f64,
    axial_distance: f64,
) -> Result<f64> {
    if ring_radius < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("ring radius must be non-negative, got {ring_radius}"),
        });
    }
    let denom_sq = axial_distance * axial_distance + ring_radius * ring_radius;
    let denom = denom_sq * denom_sq.sqrt(); // (z² + R²)^(3/2)
    if denom < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(COULOMB_K * total_charge * axial_distance / denom)
}

/// Electric field on the axis of a uniformly charged disk.
///
/// E_z = (σ/2ε₀)[1 − z/√(z² + R²)]
#[inline]
pub fn electric_field_disk_axis(
    surface_charge_density: f64,
    disk_radius: f64,
    axial_distance: f64,
) -> Result<f64> {
    if disk_radius <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("disk radius must be positive, got {disk_radius}"),
        });
    }
    let r = (axial_distance * axial_distance + disk_radius * disk_radius).sqrt();
    if r < 1e-30 {
        return Err(BijliError::Singularity);
    }
    Ok(surface_charge_density / (2.0 * EPSILON_0) * (1.0 - axial_distance / r))
}

// ── Field line tracing ─────────────────────────────────────────────

/// Trace a single field line from a starting point using Euler integration.
///
/// Steps along the local field direction with the given step size.
/// Returns a vector of positions along the field line.
///
/// `field_fn` computes the field at a given position.
pub fn trace_field_line<F>(
    start: [f64; 3],
    step_size: f64,
    max_steps: usize,
    field_fn: F,
) -> Result<Vec<[f64; 3]>>
where
    F: Fn([f64; 3]) -> Result<FieldVector>,
{
    if step_size <= 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: format!("step size must be positive, got {step_size}"),
        });
    }

    tracing::debug!(?start, step_size, max_steps, "tracing field line");

    let mut points = Vec::with_capacity(max_steps + 1);
    let mut pos = start;
    points.push(pos);

    for _ in 0..max_steps {
        let field = field_fn(pos)?;
        let mag = field.magnitude();
        if mag < 1e-30 {
            break; // field vanishes — stop tracing
        }
        // Step along the unit field direction
        let inv_mag = step_size / mag;
        pos = [
            pos[0] + field.x * inv_mag,
            pos[1] + field.y * inv_mag,
            pos[2] + field.z * inv_mag,
        ];
        points.push(pos);
    }

    Ok(points)
}

// ── Multipole expansion ───────────────────────────────────────────

/// Electric quadrupole potential at position r from a quadrupole at the origin.
///
/// V = (1/4πε₀) Σᵢⱼ Qᵢⱼ rᵢ rⱼ / (4r⁵)
///
/// `q_tensor` is the traceless quadrupole moment tensor Q_{ij} as a 3×3 array.
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

    Ok(sum / (4.0 * std::f64::consts::PI * EPSILON_0 * 4.0 * r5))
}

/// Spherical harmonic expansion coefficient for a charge distribution.
///
/// Computes the multipole moment q_{lm} = Σ qᵢ rᵢˡ Y*_{lm}(θᵢ,φᵢ)
/// for point charges at given positions.
///
/// Returns the real part of the expansion (sufficient for real charge distributions
/// with the convention q_{l,-m} = (-1)^m q*_{lm}).
///
/// `l`: multipole order (0=monopole, 1=dipole, 2=quadrupole, ...)
/// `m`: azimuthal index (-l..=l)
/// `charges`: slice of (charge, [x, y, z]) tuples
pub fn multipole_moment(l: u32, m: i32, charges: &[(f64, [f64; 3])]) -> Result<f64> {
    if m.unsigned_abs() > l {
        return Err(BijliError::InvalidParameter {
            reason: format!("|m| must be ≤ l, got l={l}, m={m}"),
        });
    }

    let mut moment = 0.0;
    for &(q, pos) in charges {
        let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        if r < 1e-30 && l > 0 {
            continue; // r^l = 0 for l > 0 at origin
        }
        let theta = if r > 1e-30 { (pos[2] / r).acos() } else { 0.0 };
        let phi = pos[1].atan2(pos[0]);

        let r_l = r.powi(l as i32);
        // Real spherical harmonic (simplified: cos(mφ) for m≥0, sin(|m|φ) for m<0)
        let angular = if m >= 0 {
            associated_legendre(l, m as u32, theta.cos()) * (m as f64 * phi).cos()
        } else {
            associated_legendre(l, m.unsigned_abs(), theta.cos())
                * (m.unsigned_abs() as f64 * phi).sin()
        };

        moment += q * r_l * angular;
    }

    Ok(moment)
}

/// Associated Legendre polynomial P_l^m(x) via recurrence.
///
/// Uses the standard physics convention (includes Condon-Shortley phase).
#[must_use]
fn associated_legendre(l: u32, m: u32, x: f64) -> f64 {
    if m > l {
        return 0.0;
    }

    // P_m^m(x) = (-1)^m (2m-1)!! (1-x²)^{m/2}
    let mut pmm = 1.0;
    if m > 0 {
        let somx2 = (1.0 - x * x).sqrt();
        let mut fact = 1.0;
        for _ in 1..=m {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    if l == m {
        return pmm;
    }

    // P_{m+1}^m(x) = x(2m+1) P_m^m
    let mut pmmp1 = x * (2 * m + 1) as f64 * pmm;
    if l == m + 1 {
        return pmmp1;
    }

    // Recurrence: (l-m) P_l^m = x(2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    let mut pll = 0.0;
    for ll in (m + 2)..=l {
        pll = (x * (2 * ll - 1) as f64 * pmmp1 - (ll + m - 1) as f64 * pmm) / (ll - m) as f64;
        pmm = pmmp1;
        pmmp1 = pll;
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
/// q' = -(R/d)q at position (R²/d) r̂
#[inline]
pub fn image_charge_sphere(
    charge: f64,
    position: [f64; 3],
    sphere_radius: f64,
) -> Result<(f64, [f64; 3])> {
    let d_sq = position[0] * position[0] + position[1] * position[1] + position[2] * position[2];
    if d_sq < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let d = d_sq.sqrt();
    if d <= sphere_radius {
        return Err(BijliError::InvalidParameter {
            reason: "charge must be outside the sphere".into(),
        });
    }

    let q_image = -charge * sphere_radius / d;
    let scale = sphere_radius * sphere_radius / d_sq;
    let pos_image = [
        position[0] * scale,
        position[1] * scale,
        position[2] * scale,
    ];

    Ok((q_image, pos_image))
}

/// Image charge for a point charge near a dielectric half-space.
///
/// Charge in medium with ε₁ (z > 0), dielectric with ε₂ (z < 0).
/// Interface at z = 0.
///
/// Returns (image_charge_for_field_in_medium1, image_position).
/// The image charge is q' = -q(ε₂ - ε₁)/(ε₂ + ε₁) at the mirror position.
#[inline]
pub fn image_charge_dielectric(
    charge: f64,
    position: [f64; 3],
    eps_r1: f64,
    eps_r2: f64,
) -> Result<(f64, [f64; 3])> {
    if position[2] < 0.0 {
        return Err(BijliError::InvalidParameter {
            reason: "charge must be in medium 1 (z > 0) for dielectric image method".into(),
        });
    }
    let den = eps_r1 + eps_r2;
    if den.abs() < 1e-30 {
        return Err(BijliError::DivisionByZero {
            context: "ε₁ + ε₂ = 0".into(),
        });
    }
    let q_image = -charge * (eps_r2 - eps_r1) / den;
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
/// Returns the 3×3 stress tensor.
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
        // P_0^0 = 1
        assert!((associated_legendre(0, 0, 0.5) - 1.0).abs() < 1e-10);
        // P_1^0 = x
        assert!((associated_legendre(1, 0, 0.5) - 0.5).abs() < 1e-10);
        // P_1^1 = -√(1-x²)
        let p11 = associated_legendre(1, 1, 0.5);
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
        let p = associated_legendre(2, 0, x);
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((p - expected).abs() < 1e-10);
    }

    #[test]
    fn test_multipole_moment_empty_charges() {
        let m = multipole_moment(1, 0, &[]).unwrap();
        assert!(m.abs() < 1e-30);
    }
}
