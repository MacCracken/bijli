//! Relativistic electrodynamics — EM tensor, Lorentz field transforms,
//! retarded potentials, Liénard-Wiechert fields.
//!
//! All SI units. Indices: 0 = ct, 1 = x, 2 = y, 3 = z.
//!
//! # Sign conventions
//!
//! - [`FourVector::dot`] and [`FourVector::interval_sq`] use the (−,+,+,+)
//!   metric η = diag(−1, 1, 1, 1).
//! - [`EmTensor::components`] stores F^μν in the Jackson/Griffiths layout
//!   (F^{0i} = −Eᵢ/c, F^{12} = −B_z), which is the contravariant field tensor of
//!   the **(+,−,−,−)** signature — the negative of F^μν under η. Contract it
//!   with g = diag(1, −1, −1, −1), or use [`EmTensor::four_force`]; see the
//!   [`EmTensor`] docs.

use crate::error::{BijliError, Result};
use crate::field::{FieldVector, MU_0, SPEED_OF_LIGHT, SPEED_OF_LIGHT_SQ};

/// Speed of light squared, alias for readability.
const C2: f64 = SPEED_OF_LIGHT_SQ;

/// Speed of light cubed.
const C3: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT_SQ;

// ── Lorentz factor ─────────────────────────────────────────────────

/// Lorentz factor γ = 1/√(1 − v²/c²).
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if |v| ≥ c or `speed` is NaN.
#[inline]
pub fn lorentz_factor(speed: f64) -> Result<f64> {
    let beta_sq = speed * speed / C2;
    if beta_sq.is_nan() || beta_sq >= 1.0 {
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
/// A 4×4 antisymmetric matrix encoding E and B fields (indices 0 = ct, 1..3 = x, y, z):
/// ```text
///        ⎡  0    -Ex/c  -Ey/c  -Ez/c ⎤
/// F^μν = ⎢ Ex/c    0    -Bz     By   ⎥
///        ⎢ Ey/c   Bz      0    -Bx   ⎥
///        ⎣ Ez/c  -By     Bx      0    ⎦
/// ```
///
/// # Sign convention
///
/// This is the layout of Jackson (eq. 11.137, in SI form) and Griffiths
/// (eq. 12.119): F^{0i} = −Eᵢ/c, F^{ij} = −ε_{ijk}B_k. It is the contravariant
/// tensor F^μν = ∂^μA^ν − ∂^νA^μ of the **(+,−,−,−)** signature, *not* of the
/// (−,+,+,+) metric used by [`FourVector::dot`] — under that metric
/// F^{0i} = +Eᵢ/c and F^{ij} = +ε_{ijk}B_k, i.e. the negative of this matrix.
///
/// Contract these components with g = diag(1, −1, −1, −1). For example the
/// Lorentz four-force on a charge q with four-velocity u = γ(c, v) is
///
/// ```text
/// K^μ = q F^μν g_νλ u^λ = q (F^{μ0}u⁰ − F^{μ1}u¹ − F^{μ2}u² − F^{μ3}u³)
///     = (γq E⋅v/c,  γq (E + v × B)),
/// ```
///
/// which [`EmTensor::four_force`] implements. Contracting with η instead would
/// give −K. The invariants below are independent of the signature choice.
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

    /// Second Lorentz invariant: ε_μνρσ F^μν F^ρσ / 8 = −E⋅B/c,
    /// with Jackson's Levi-Civita convention ε^{0123} = +1 (so ε_{0123} = −1).
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

    /// Lorentz four-force K^μ = dp^μ/dτ = q F^μ_ν u^ν on a charge (N).
    ///
    /// `four_velocity` is u^μ = γ(c, v) (m/s), stored in the (ct, x, y, z) slots
    /// of [`FourVector`]. Returns K = (γq E⋅v/c, γq(E + v × B)); K⁰ = (γ/c) dW/dt
    /// is the power term, and K⋅u = 0 under [`FourVector::dot`].
    ///
    /// Uses the stored components with the (+,−,−,−) metric (see the type docs):
    /// K^μ = q (F^{μ0}u⁰ − F^{μ1}u¹ − F^{μ2}u² − F^{μ3}u³).
    #[inline]
    #[must_use]
    pub fn four_force(&self, charge: f64, four_velocity: &FourVector) -> FourVector {
        let f = &self.components;
        // Covariant u_ν = g_νλ u^λ with g = diag(1, −1, −1, −1).
        let u_low = [
            four_velocity.t,
            -four_velocity.x,
            -four_velocity.y,
            -four_velocity.z,
        ];
        let mut k = [0.0; 4];
        for (k_mu, row) in k.iter_mut().zip(f.iter()) {
            *k_mu = charge
                * (row[0] * u_low[0] + row[1] * u_low[1] + row[2] * u_low[2] + row[3] * u_low[3]);
        }
        FourVector::new(k[0], k[1], k[2], k[3])
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
///
/// evaluated in the equivalent normalization-free form (Jackson eq. 11.149)
///   E' = γ(E + v × B) − γ²/(γ+1) (v⋅E) v/c²
///   B' = γ(B − v × E/c²) − γ²/(γ+1) (v⋅B) v/c²,
/// which stays exact for arbitrarily small (non-zero) |v|.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if |v| ≥ c or `velocity` has a NaN component.
#[inline]
pub fn lorentz_transform_fields(
    e: &FieldVector,
    b: &FieldVector,
    velocity: &FieldVector,
) -> Result<(FieldVector, FieldVector)> {
    let v_sq = velocity.magnitude_sq();
    if v_sq.is_nan() || v_sq >= C2 {
        return Err(BijliError::InvalidParameter {
            reason: "boost velocity must be less than c".into(),
        });
    }

    let gamma = 1.0 / (1.0 - v_sq / C2).sqrt();
    let k = gamma * gamma / ((gamma + 1.0) * C2);

    let e_prime = (*e + velocity.cross(b)) * gamma - velocity.scale(k * velocity.dot(e));
    let b_prime =
        (*b - velocity.cross(e).scale(1.0 / C2)) * gamma - velocity.scale(k * velocity.dot(b));

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

/// Retarded scalar potential from a point charge, given the effective distance.
///
/// φ(r, t) = q / (4πε₀ R_eff)   (V)
///
/// `retarded_distance` must be the Liénard–Wiechert *effective* distance
///
/// ```text
/// R_eff = |R| − R⋅β = |R| (1 − n̂⋅β),   R = r − r'(t_r),  β = v(t_r)/c,
/// ```
///
/// evaluated at the retarded time t_r = t − |R|/c (see [`solve_retarded_time`]).
/// With that argument the result is the exact Liénard–Wiechert potential
/// (Griffiths eq. 10.46, Jackson eq. 14.8). For a charge at rest R_eff = |R|
/// and this is the Coulomb potential; passing the plain distance |R| for a
/// moving charge omits the 1/(1 − n̂⋅β) factor, an O(β) error. To compute R_eff
/// from vectors, use [`lienard_wiechert_scalar_potential`].
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `retarded_distance` is negative or NaN,
/// [`BijliError::Singularity`] if it is (numerically) zero.
#[inline]
pub fn retarded_scalar_potential(charge: f64, retarded_distance: f64) -> Result<f64> {
    crate::field::check_radial_distance(retarded_distance)?;
    Ok(charge / (4.0 * std::f64::consts::PI * crate::field::EPSILON_0 * retarded_distance))
}

/// Retarded vector potential from a moving point charge, given the effective distance.
///
/// A(r, t) = μ₀ q v / (4π R_eff)   (T⋅m)
///
/// `velocity` is v(t_r) and `retarded_distance` the Liénard–Wiechert
/// effective distance R_eff = |R| − R⋅v/c, both at the retarded time, exactly
/// as for [`retarded_scalar_potential`] (then A = (v/c²) φ). To compute R_eff
/// from vectors, use [`lienard_wiechert_vector_potential`].
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if `retarded_distance` is negative or NaN,
/// [`BijliError::Singularity`] if it is (numerically) zero.
#[inline]
pub fn retarded_vector_potential(
    charge: f64,
    velocity: &FieldVector,
    retarded_distance: f64,
) -> Result<FieldVector> {
    crate::field::check_radial_distance(retarded_distance)?;
    let factor = MU_0 * charge / (4.0 * std::f64::consts::PI * retarded_distance);
    Ok(velocity.scale(factor))
}

/// Liénard–Wiechert effective distance R_eff = |R| − R⋅β (> 0 for |β| < 1).
#[inline]
fn lienard_wiechert_effective_distance(r_vec: &FieldVector, beta_vec: &FieldVector) -> Result<f64> {
    let r_mag = r_vec.magnitude();
    if r_mag.is_nan() {
        return Err(BijliError::InvalidParameter {
            reason: "R vector must not contain NaN".into(),
        });
    }
    if r_mag < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let beta_sq = beta_vec.magnitude_sq();
    if beta_sq.is_nan() || beta_sq >= 1.0 {
        return Err(BijliError::InvalidParameter {
            reason: "β must be < 1".into(),
        });
    }
    Ok(r_mag - r_vec.dot(beta_vec))
}

/// Liénard–Wiechert scalar potential of a moving point charge.
///
/// φ(r, t) = q / (4πε₀ (|R| − R⋅β)) = q / (4πε₀ |R| (1 − n̂⋅β))   (V)
///
/// (Griffiths eq. 10.46, Jackson eq. 14.8.)
///
/// # Parameters
/// - `charge`: charge in coulombs
/// - `r_vec`: R = r − r'(t_r), from the retarded charge position to the field point (m)
/// - `beta_vec`: β = v(t_r)/c at the retarded time
///
/// # Errors
///
/// [`BijliError::Singularity`] if R = 0, [`BijliError::InvalidParameter`] if
/// |β| ≥ 1 or an input is NaN.
#[inline]
pub fn lienard_wiechert_scalar_potential(
    charge: f64,
    r_vec: &FieldVector,
    beta_vec: &FieldVector,
) -> Result<f64> {
    let r_eff = lienard_wiechert_effective_distance(r_vec, beta_vec)?;
    Ok(charge / (4.0 * std::f64::consts::PI * crate::field::EPSILON_0 * r_eff))
}

/// Liénard–Wiechert vector potential of a moving point charge.
///
/// A(r, t) = μ₀ q c β / (4π (|R| − R⋅β)) = (β/c) φ   (T⋅m)
///
/// (Griffiths eq. 10.47.) Parameters and errors as for
/// [`lienard_wiechert_scalar_potential`].
#[inline]
pub fn lienard_wiechert_vector_potential(
    charge: f64,
    r_vec: &FieldVector,
    beta_vec: &FieldVector,
) -> Result<FieldVector> {
    let r_eff = lienard_wiechert_effective_distance(r_vec, beta_vec)?;
    let factor = MU_0 * charge * SPEED_OF_LIGHT / (4.0 * std::f64::consts::PI * r_eff);
    Ok(beta_vec.scale(factor))
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
    if beta_sq.is_nan() || beta_sq >= 1.0 {
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
/// P = q²γ⁶ / (6πε₀c³) × [a² − |v × a|²/c²]   (W)
///
/// In the rest frame this reduces to the non-relativistic Larmor formula.
#[inline]
pub fn relativistic_larmor_power(
    charge: f64,
    velocity: &FieldVector,
    acceleration: &FieldVector,
) -> Result<f64> {
    let v_sq = velocity.magnitude_sq();
    if v_sq.is_nan() || v_sq >= C2 {
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

// ── Retarded time solver ──────────────────────────────────────────

/// Maximum number of bracket doublings in [`solve_retarded_time`]; 2²¹⁰⁰
/// spans the whole `f64` range, so this only stops unbounded searches.
const RETARDED_BRACKET_MAX_EXPANSIONS: usize = 2100;

/// Refinement budget used by [`lienard_wiechert_fields`].
const LW_MAX_ITERATIONS: usize = 200;

/// Solve for the retarded time t_r given observer position, observation time,
/// and a charge trajectory.
///
/// Finds the unique t_r ≤ t with |r_obs − r_charge(t_r)| = c(t − t_r) for a
/// charge moving slower than light.
///
/// # Method
///
/// Works on the delay τ = t − t_r ≥ 0 with g(τ) = |r_obs − r(t − τ)| − cτ.
/// g(0) ≥ 0 and, for |v| < c, g'(τ) = R̂⋅v − c < 0, so the root is unique.
/// Starting from the static estimate τ₀ = |r_obs − r(t)|/c, a bracket
/// [τ_lo, τ_hi] with g(τ_lo) > 0 ≥ g(τ_hi) is found by doubling τ_hi, then
/// refined by Newton steps safeguarded by bisection (a Newton step that leaves
/// the bracket, or does not at least halve the previous step, is replaced by
/// bisection). The iteration stops once the step (or half the bracket) is
/// ≤ 4ε·max(|t|, τ) — the resolution with which t_r = t − τ is representable —
/// or once |g| has dropped to the rounding level of its own inputs,
/// 2ε·max(|R|, cτ, ‖r_obs‖∞, ‖r(t_r)‖∞, ‖v‖∞|t_r|), after which one last
/// Newton correction is applied. (The second test matters when a charge
/// approaches the observer at nearly c, where |g'| = c − R̂⋅v is small.) No
/// fixed absolute tolerance on g is used: c(t − t_r) alone carries rounding
/// noise of ≈ c·ulp(t) (≈ 7 × 10⁻⁸ m at t = 1 s), which made the previous
/// `|g| < 10⁻¹² m` criterion unreachable for ordinary t.
///
/// # Parameters
///
/// - `trajectory_fn(t)` returns the charge's (position in m, velocity in m/s)
///   at time t; the motion must be subluminal.
/// - `max_iterations`: maximum number of Newton/bisection refinement steps
///   (one trajectory evaluation each; quadratic convergence typically needs
///   2–5, pure bisection at most ≈ 55). Bracketing costs an additional
///   ≈ 2 + log₂(1/(1 − β)) evaluations for speed βc that are not counted.
///
/// # Errors
///
/// [`BijliError::InvalidParameter`] if the observer position or time is not
/// finite, the trajectory returns a non-finite state, no bracket exists
/// (e.g. superluminal motion), or the refinement has not converged after
/// `max_iterations` steps.
pub fn solve_retarded_time<F>(
    observer_pos: [f64; 3],
    observation_time: f64,
    trajectory_fn: F,
    max_iterations: usize,
) -> Result<f64>
where
    F: Fn(f64) -> ([f64; 3], [f64; 3]),
{
    if !(observation_time.is_finite() && observer_pos.iter().all(|c| c.is_finite())) {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "observer position and time must be finite, got {observer_pos:?} at t = {observation_time}"
            ),
        });
    }
    tracing::debug!(
        ?observer_pos,
        observation_time,
        max_iterations,
        "solving retarded time"
    );
    let t = observation_time;

    // g, dg/dτ and the rounding scale of g at t_r = t − τ. The light-travel
    // term uses c·(t − t_r) with the rounded t_r, so the root satisfies the
    // light-cone condition for the value actually returned.
    let eval = |tau: f64| -> Result<LightConeSample> {
        let t_r = t - tau;
        let (pos, vel) = trajectory_fn(t_r);
        if !pos.iter().chain(vel.iter()).all(|c| c.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("trajectory returned a non-finite state at t = {t_r}"),
            });
        }
        let dx = observer_pos[0] - pos[0];
        let dy = observer_pos[1] - pos[1];
        let dz = observer_pos[2] - pos[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let light = SPEED_OF_LIGHT * (t - t_r);
        let dg = if dist > 0.0 {
            (dx * vel[0] + dy * vel[1] + dz * vel[2]) / dist - SPEED_OF_LIGHT
        } else {
            -SPEED_OF_LIGHT
        };
        let max_abs = |v: &[f64; 3]| v[0].abs().max(v[1].abs()).max(v[2].abs());
        Ok(LightConeSample {
            tau,
            g: dist - light,
            dg,
            // |g| below 2ε·max(...) is indistinguishable from zero: the inputs to g
            // (positions, distance, light-travel term) carry that much rounding.
            // (|v|·|t_r| covers time rounding inside `trajectory_fn` itself.)
            noise: 2.0
                * f64::EPSILON
                * dist
                    .max(light)
                    .max(max_abs(&observer_pos))
                    .max(max_abs(&pos))
                    .max(max_abs(&vel) * t_r.abs()),
        })
    };
    let tol = |tau: f64| 4.0 * f64::EPSILON * t.abs().max(tau);

    let s0 = eval(0.0)?;
    if s0.g == 0.0 {
        // Observer coincides with the charge at time t.
        tracing::debug!(t_r = t, "retarded time: observer on the trajectory");
        return Ok(t);
    }

    // ── Bracket: g(lo) > 0 ≥ g(hi) ──
    let mut lo = s0;
    let mut hi = eval((s0.g / SPEED_OF_LIGHT).max(f64::MIN_POSITIVE))?;
    let mut expansions = 0_usize;
    while hi.g > 0.0 {
        let next_tau = 2.0 * hi.tau;
        if expansions >= RETARDED_BRACKET_MAX_EXPANSIONS || !next_tau.is_finite() {
            return Err(BijliError::InvalidParameter {
                reason: "retarded time solver could not bracket the light-cone root \
                         (is the trajectory faster than light?)"
                    .into(),
            });
        }
        lo = hi;
        hi = eval(next_tau)?;
        expansions += 1;
    }

    // ── Safeguarded Newton / bisection on [lo.tau, hi.tau] ──
    let mut cur = if lo.g.abs() < hi.g.abs() { lo } else { hi };
    let mut last_step = hi.tau - lo.tau;
    let mut iterations = 0_usize;
    let converged_tau = loop {
        if cur.g == 0.0 {
            break Some(cur.tau);
        }
        if cur.g.abs() <= cur.noise {
            // Residual at the rounding floor: apply one final (unevaluated)
            // Newton correction, kept inside the bracket.
            let polished = if cur.dg < 0.0 {
                (cur.tau - cur.g / cur.dg).max(lo.tau).min(hi.tau)
            } else {
                cur.tau
            };
            break Some(polished);
        }
        if iterations == max_iterations {
            break None;
        }
        iterations += 1;

        let newton = cur.tau - cur.g / cur.dg;
        let use_newton = cur.dg < 0.0
            && newton >= lo.tau
            && newton <= hi.tau
            && (2.0 * cur.g).abs() <= (last_step * cur.dg).abs();
        let (next, step) = if use_newton {
            (newton, (cur.g / cur.dg).abs())
        } else {
            let half = 0.5 * (hi.tau - lo.tau);
            (lo.tau + half, half)
        };
        if step <= tol(next) || next == cur.tau {
            break Some(next);
        }
        last_step = step;
        cur = eval(next)?;
        if cur.g > 0.0 {
            lo = cur;
        } else {
            hi = cur;
        }
    };

    match converged_tau {
        Some(tau) => {
            tracing::debug!(
                t_r = t - tau,
                tau,
                iterations,
                expansions,
                "retarded time converged"
            );
            Ok(t - tau)
        }
        None => {
            tracing::debug!(
                max_iterations,
                lo = lo.tau,
                hi = hi.tau,
                "retarded time solver did not converge"
            );
            Err(BijliError::InvalidParameter {
                reason: format!(
                    "retarded time solver did not converge in {max_iterations} iterations"
                ),
            })
        }
    }
}

/// One evaluation of the light-cone function g(τ) = |r_obs − r(t − τ)| − cτ.
#[derive(Debug, Clone, Copy)]
struct LightConeSample {
    tau: f64,
    g: f64,
    dg: f64,
    noise: f64,
}

/// Complete Liénard-Wiechert fields (E and B) at an observation point.
///
/// Combines the retarded time solver ([`solve_retarded_time`], with an
/// internal budget of 200 refinement steps) with [`lienard_wiechert_e`] and
/// [`lienard_wiechert_b`]. Returns (E, B) at the observer position (V/m, T).
///
/// `trajectory_fn` returns (position, velocity, acceleration) at time t (SI);
/// the motion must be subluminal.
///
/// # Errors
///
/// Propagates the errors of [`solve_retarded_time`] and
/// [`lienard_wiechert_e`]; [`BijliError::Singularity`] if the observer sits
/// on the charge's retarded position; [`BijliError::InvalidParameter`] if the
/// acceleration is not finite.
pub fn lienard_wiechert_fields<F>(
    charge: f64,
    observer_pos: [f64; 3],
    observation_time: f64,
    trajectory_fn: F,
) -> Result<(FieldVector, FieldVector)>
where
    F: Fn(f64) -> ([f64; 3], [f64; 3], [f64; 3]),
{
    tracing::trace!(
        charge,
        ?observer_pos,
        observation_time,
        "Liénard-Wiechert fields"
    );
    // Solve for retarded time using position+velocity only
    let t_r = solve_retarded_time(
        observer_pos,
        observation_time,
        |t| {
            let (pos, vel, _) = trajectory_fn(t);
            (pos, vel)
        },
        LW_MAX_ITERATIONS,
    )?;

    let (pos_r, vel_r, acc_r) = trajectory_fn(t_r);
    if !acc_r.iter().all(|a| a.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!("trajectory returned a non-finite acceleration at t = {t_r}"),
        });
    }

    // r⃗ = observer - retarded position
    let r_vec = FieldVector::new(
        observer_pos[0] - pos_r[0],
        observer_pos[1] - pos_r[1],
        observer_pos[2] - pos_r[2],
    );
    let r_mag = r_vec.magnitude();
    if r_mag < 1e-30 {
        return Err(BijliError::Singularity);
    }
    let n_hat = r_vec.scale(1.0 / r_mag);

    // β⃗ = v⃗/c, β̇ = a⃗/c
    let beta_vec = FieldVector::new(
        vel_r[0] / SPEED_OF_LIGHT,
        vel_r[1] / SPEED_OF_LIGHT,
        vel_r[2] / SPEED_OF_LIGHT,
    );
    let beta_dot = FieldVector::new(
        acc_r[0] / SPEED_OF_LIGHT,
        acc_r[1] / SPEED_OF_LIGHT,
        acc_r[2] / SPEED_OF_LIGHT,
    );

    let e = lienard_wiechert_e(charge, &r_vec, &beta_vec, &beta_dot)?;
    let b = lienard_wiechert_b(&n_hat, &e);

    Ok((e, b))
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

    // ── Retarded time solver tests ────────────────────────────────

    #[test]
    fn test_retarded_time_static_charge() {
        // Static charge at origin: t_r = t - r/c
        let obs = [1.0, 0.0, 0.0];
        let t_obs = 1e-8;
        let t_r = solve_retarded_time(obs, t_obs, |_t| ([0.0; 3], [0.0; 3]), 50).unwrap();
        let expected = t_obs - 1.0 / SPEED_OF_LIGHT;
        assert!((t_r - expected).abs() < 1e-18);
    }

    #[test]
    fn test_retarded_time_moving_charge() {
        // Charge moving in +x at constant velocity
        let v = 0.1 * SPEED_OF_LIGHT;
        let obs = [2.0, 0.0, 0.0];
        let t_obs = 1e-8;
        let t_r =
            solve_retarded_time(obs, t_obs, |t| ([v * t, 0.0, 0.0], [v, 0.0, 0.0]), 50).unwrap();
        // Verify light-cone condition: |obs - pos(t_r)| = c(t_obs - t_r)
        let pos_r = v * t_r;
        let r = (obs[0] - pos_r).abs();
        let c_dt = SPEED_OF_LIGHT * (t_obs - t_r);
        assert!((r - c_dt).abs() < 1e-10);
    }

    #[test]
    fn test_lw_fields_static() {
        // Static charge: convenience function should match direct computation
        let q = 1e-6;
        let obs = [0.0, 1.0, 0.0];
        let t = 1e-8;
        let (e, b) =
            lienard_wiechert_fields(q, obs, t, |_t| ([0.0; 3], [0.0; 3], [0.0; 3])).unwrap();
        // E should be Coulomb field pointing in +y
        assert!(e.y > 0.0);
        assert!(e.x.abs() < 1e-3);
        // B should be ~zero for static charge
        assert!(b.magnitude() < 1e-10);
    }

    // ── Audit regression tests (1.2.0) ────────────────────────────

    const Q_E: f64 = 1.602_176_634e-19;

    /// Closed-form delay τ = t − t_r for uniform motion, where `w` is the
    /// vector from the charge's *present* position to the observer:
    /// (c² − v²)τ² − 2(w⋅v)τ − w² = 0, positive root (cancellation-free form).
    fn uniform_motion_delay(w: [f64; 3], v: [f64; 3]) -> f64 {
        let wv = w[0] * v[0] + w[1] * v[1] + w[2] * v[2];
        let w2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let s = (wv * wv + (C2 - v2) * w2).sqrt();
        if wv >= 0.0 {
            (wv + s) / (C2 - v2)
        } else {
            w2 / (s - wv)
        }
    }

    #[test]
    fn test_retarded_time_static_charge_large_t() {
        // Rounding noise in c(t − t_r) is ≈ c·ulp(t) ≫ 1e-12 m for t ≳ 1e-4 s;
        // the solver must still converge.
        let obs = [1.0, 0.0, 0.0];
        for t_obs in [1e-3, 1.0, 10.0, 1e4] {
            let t_r = solve_retarded_time(obs, t_obs, |_t| ([0.0; 3], [0.0; 3]), 50).unwrap();
            let expected = t_obs - 1.0 / SPEED_OF_LIGHT;
            assert!(t_r <= t_obs);
            assert!(
                (t_r - expected).abs() <= 8.0 * f64::EPSILON * t_obs,
                "t = {t_obs}: t_r = {t_r}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_retarded_time_uniform_motion_closed_form() {
        let v = [
            0.8 * SPEED_OF_LIGHT,
            -0.3 * SPEED_OF_LIGHT,
            0.1 * SPEED_OF_LIGHT,
        ];
        let r0 = [0.1, -0.2, 0.05];
        let obs = [1.0, 0.5, 0.2];
        // Charge is at r0 at time t0; test small and large observation times.
        for (t0, t_obs) in [(0.0, 2e-9), (1.0, 1.0 + 3e-9), (-5.0, -5.0 + 1e-8)] {
            let traj = |t: f64| {
                let dt = t - t0;
                ([r0[0] + v[0] * dt, r0[1] + v[1] * dt, r0[2] + v[2] * dt], v)
            };
            let t_r = solve_retarded_time(obs, t_obs, traj, 50).unwrap();
            let (p_now, _) = traj(t_obs);
            let w = [obs[0] - p_now[0], obs[1] - p_now[1], obs[2] - p_now[2]];
            let tau = uniform_motion_delay(w, v);
            let expected = t_obs - tau;
            assert!(t_r <= t_obs);
            assert!(
                (t_r - expected).abs() <= 1e-13 * tau + 8.0 * f64::EPSILON * t_obs.abs(),
                "t = {t_obs}: t_r = {t_r}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_retarded_time_max_iterations_respected() {
        let traj = |t: f64| {
            (
                [0.5 * SPEED_OF_LIGHT * t, 0.0, 0.0],
                [0.5 * SPEED_OF_LIGHT, 0.0, 0.0],
            )
        };
        assert!(solve_retarded_time([1.0, 0.3, 0.0], 1e-8, traj, 0).is_err());
        assert!(solve_retarded_time([1.0, 0.3, 0.0], 1e-8, traj, 50).is_ok());
    }

    #[test]
    fn test_retarded_time_rejects_superluminal_and_nan() {
        // Trajectory moving away faster than light: no retarded solution.
        let fast = |t: f64| {
            (
                [-2.0 * SPEED_OF_LIGHT * t, 0.0, 0.0],
                [-2.0 * SPEED_OF_LIGHT, 0.0, 0.0],
            )
        };
        assert!(solve_retarded_time([1.0, 0.0, 0.0], 1.0, fast, 100).is_err());
        let bad = |_t: f64| ([f64::NAN, 0.0, 0.0], [0.0; 3]);
        assert!(solve_retarded_time([1.0, 0.0, 0.0], 1.0, bad, 100).is_err());
        let still = |_t: f64| ([0.0; 3], [0.0; 3]);
        assert!(solve_retarded_time([1.0, 0.0, 0.0], f64::NAN, still, 100).is_err());
    }

    #[test]
    fn test_lw_fields_uniform_motion_heaviside() {
        // Griffiths 10.75/10.76: field of a uniformly moving charge, in terms of
        // the present position: E = q(1−β²)/(4πε₀(1 − β²sin²θ)^{3/2}) ŵ/w²,
        // B = v × E / c².
        let q = 1e-9;
        let v = [0.9 * SPEED_OF_LIGHT, 0.0, 0.0];
        let obs = [0.3, 0.8, -0.4];
        for t_obs in [3e-9, 2.0] {
            let traj = |t: f64| ([v[0] * (t - t_obs), 0.0, 0.0], v, [0.0; 3]);
            let (e, b) = lienard_wiechert_fields(q, obs, t_obs, traj).unwrap();
            let (p_now, _, _) = traj(t_obs);
            let w = FieldVector::new(obs[0] - p_now[0], obs[1] - p_now[1], obs[2] - p_now[2]);
            let wm = w.magnitude();
            let beta2 = (v[0] / SPEED_OF_LIGHT).powi(2);
            let cos_t = w.x / wm;
            let sin2 = 1.0 - cos_t * cos_t;
            let mag = q / (4.0 * std::f64::consts::PI * EPSILON_0) * (1.0 - beta2)
                / ((1.0 - beta2 * sin2).powf(1.5) * wm * wm);
            let e_ref = w.scale(mag / wm);
            let b_ref = FieldVector::new(v[0], v[1], v[2])
                .cross(&e_ref)
                .scale(1.0 / C2);
            assert!(
                (e - e_ref).magnitude() / e_ref.magnitude() < 1e-9,
                "t = {t_obs}: E = {e}, ref = {e_ref}"
            );
            assert!((b - b_ref).magnitude() / b_ref.magnitude() < 1e-9);
        }
    }

    #[test]
    fn test_lw_fields_circular_orbit_beta_099() {
        // R = 1 m, β = 0.99 circular orbit in the xy-plane. References for t_r,
        // E and B from a 60-digit bisection + Griffiths 10.72 (python decimal).
        // The first two cases made the old Newton solver exhaust 50 iterations.
        let big_r = 1.0;
        let omega = 0.99 * SPEED_OF_LIGHT / big_r;
        let traj_pv = |t: f64| {
            let (s, c) = (omega * t).sin_cos();
            (
                [big_r * c, big_r * s, 0.0],
                [-big_r * omega * s, big_r * omega * c, 0.0],
            )
        };
        let traj_pva = |t: f64| {
            let (p, v) = traj_pv(t);
            let (s, c) = (omega * t).sin_cos();
            let a2 = big_r * omega * omega;
            (p, v, [-a2 * c, -a2 * s, 0.0])
        };
        type Case = ([f64; 3], f64, f64, [f64; 3], [f64; 3]);
        let cases: [Case; 3] = [
            (
                [0.48, -1.02, 0.0],
                8.37e-8,
                8.222_399_231_049_45e-8,
                [6.227_148_912_038_439e-10, -5.141_662_384_120_919e-10, 0.0],
                [0.0, 0.0, 2.690_862_240_685_010_3e-18],
            ),
            (
                [1.23, -0.59, 0.01],
                2e-10,
                -1.102_167_177_929_938_3e-9,
                [
                    1.266_102_956_952_649_4e-9,
                    1.196_125_255_494_040_6e-9,
                    -2.799_254_661_512_416e-11,
                ],
                [
                    -3.793_853_944_382_319e-20,
                    1.758_794_283_294_848e-19,
                    5.799_391_223_903_531e-18,
                ],
            ),
            (
                [0.89, 0.805, -0.081],
                6.93e-8,
                6.751_098_042_729_196e-8,
                [
                    1.826_043_746_748_316_8e-10,
                    6.531_739_882_140_894e-10,
                    5.759_862_144_199_076e-12,
                ],
                [
                    3.246_628_266_085_001e-19,
                    -1.104_694_176_817_058_7e-19,
                    2.234_584_329_904_851e-18,
                ],
            ),
        ];
        for (obs, t_obs, t_r_ref, e_ref, b_ref) in cases {
            let t_r = solve_retarded_time(obs, t_obs, traj_pv, 50).unwrap();
            assert!(t_r <= t_obs);
            assert!((t_r - t_r_ref).abs() < 1e-21, "obs {obs:?}: t_r = {t_r}");

            let (e, b) = lienard_wiechert_fields(Q_E, obs, t_obs, traj_pva).unwrap();
            let e_ref = FieldVector::from(e_ref);
            let b_ref = FieldVector::from(b_ref);
            assert!(
                (e - e_ref).magnitude() / e_ref.magnitude() < 1e-8,
                "E = {e}"
            );
            assert!(
                (b - b_ref).magnitude() / b_ref.magnitude() < 1e-8,
                "B = {b}"
            );
        }
    }

    #[test]
    fn test_lorentz_factor_nan_rejected() {
        assert!(lorentz_factor(f64::NAN).is_err());
        assert!(lorentz_factor(f64::INFINITY).is_err());
    }

    #[test]
    fn test_lorentz_transform_fields_nan_rejected() {
        let e = FieldVector::new(1.0, 0.0, 0.0);
        let b = FieldVector::zero();
        let v = FieldVector::new(f64::NAN, 0.0, 0.0);
        assert!(lorentz_transform_fields(&e, &b, &v).is_err());
        let v = FieldVector::new(1.1 * SPEED_OF_LIGHT, 0.0, 0.0);
        assert!(lorentz_transform_fields(&e, &b, &v).is_err());
    }

    #[test]
    fn test_lorentz_transform_fields_tiny_velocity() {
        // Pure B = 1 T, v = 1e-16 m/s ŷ: E' = v × B = 1e-16 V/m x̂ (γ = 1).
        let e = FieldVector::zero();
        let b = FieldVector::new(0.0, 0.0, 1.0);
        let v = FieldVector::new(0.0, 1e-16, 0.0);
        let (ep, bp) = lorentz_transform_fields(&e, &b, &v).unwrap();
        assert!((ep.x - 1e-16).abs() < 1e-30, "E' = {ep}");
        assert!((bp.z - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_lorentz_transform_fields_oblique_matches_rotated_x_boost() {
        // Boost along (1,1,0)/√2: compare with lorentz_transform_x in a frame
        // rotated by 45° about z.
        let e = FieldVector::new(1e3, -2e3, 5e2);
        let b = FieldVector::new(2e-6, 1e-6, -3e-6);
        let speed = 0.75 * SPEED_OF_LIGHT;
        let h = std::f64::consts::FRAC_1_SQRT_2;
        let rot = |v: &FieldVector| FieldVector::new(h * (v.x + v.y), h * (v.y - v.x), v.z);
        let unrot = |v: &FieldVector| FieldVector::new(h * (v.x - v.y), h * (v.x + v.y), v.z);
        let (ex, bx) = lorentz_transform_x(&rot(&e), &rot(&b), speed).unwrap();
        let (e_ref, b_ref) = (unrot(&ex), unrot(&bx));
        let vel = FieldVector::new(h * speed, h * speed, 0.0);
        let (ep, bp) = lorentz_transform_fields(&e, &b, &vel).unwrap();
        assert!((ep - e_ref).magnitude() / e_ref.magnitude() < 1e-12);
        assert!((bp - b_ref).magnitude() / b_ref.magnitude() < 1e-12);
    }

    #[test]
    fn test_relativistic_larmor_reference_value() {
        // P = q²γ⁶(a² − |v×a|²/c²)/(6πε₀c³), v = 0.6c x̂, a = 1e18 ŷ, q = e
        // → 1.393_634_464_118_5e-17 W (python decimal reference).
        let v = FieldVector::new(0.6 * SPEED_OF_LIGHT, 0.0, 0.0);
        let a = FieldVector::new(0.0, 1e18, 0.0);
        let p = relativistic_larmor_power(Q_E, &v, &a).unwrap();
        assert!((p - 1.393_634_464_118_5e-17).abs() / p < 1e-12, "{p}");
    }

    #[test]
    fn test_four_force_matches_lorentz_force() {
        // K = (γq E⋅v/c, γq(E + v × B)) with u = γ(c, v); also K⋅u = 0 under η.
        let q = -Q_E;
        let e = FieldVector::new(3e4, -1e4, 2e4);
        let b = FieldVector::new(0.2, 0.5, -0.1);
        let v = FieldVector::new(
            0.3 * SPEED_OF_LIGHT,
            -0.5 * SPEED_OF_LIGHT,
            0.2 * SPEED_OF_LIGHT,
        );
        let gamma = 1.0 / (1.0 - v.magnitude_sq() / C2).sqrt();
        let u = FourVector::new(
            gamma * SPEED_OF_LIGHT,
            gamma * v.x,
            gamma * v.y,
            gamma * v.z,
        );
        let f = EmTensor::from_fields(&e, &b);
        // Layout pinned by the docs: F^{01} = −E_x/c, F^{12} = −B_z.
        assert_eq!(f.components[0][1], -e.x / SPEED_OF_LIGHT);
        assert_eq!(f.components[1][2], -b.z);

        let k = f.four_force(q, &u);
        // Independent 3-vector Lorentz force q(E + v × B)
        let v_cross_b = FieldVector::new(
            v.y * b.z - v.z * b.y,
            v.z * b.x - v.x * b.z,
            v.x * b.y - v.y * b.x,
        );
        let f3 = (e + v_cross_b).scale(q);
        let scale = gamma * f3.magnitude();
        assert!((k.x - gamma * f3.x).abs() < 1e-12 * scale);
        assert!((k.y - gamma * f3.y).abs() < 1e-12 * scale);
        assert!((k.z - gamma * f3.z).abs() < 1e-12 * scale);
        let k0 = gamma * q * e.dot(&v) / SPEED_OF_LIGHT;
        assert!((k.t - k0).abs() < 1e-12 * scale);
        // Four-force is orthogonal to the four-velocity: η(K, u) = 0.
        assert!(k.dot(&u).abs() < 1e-9 * scale * u.t);
    }

    #[test]
    fn test_lienard_wiechert_potentials_uniform_motion_closed_form() {
        // Griffiths Ex. 10.4 (charge through the origin at t = 0, velocity v):
        //   V = (1/4πε₀) qc / √((c²t − r⋅v)² + (c² − v²)(r² − c²t²)),  A = (v/c²) V
        let q = 1e-9;
        let v = [
            0.6 * SPEED_OF_LIGHT,
            0.2 * SPEED_OF_LIGHT,
            -0.1 * SPEED_OF_LIGHT,
        ];
        let vv = FieldVector::from(v);
        let r = [0.5, -0.3, 0.8];
        let t = 4e-9;
        let rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let v2 = vv.magnitude_sq();
        let v_ref = q * SPEED_OF_LIGHT
            / (4.0 * std::f64::consts::PI * EPSILON_0)
            / ((C2 * t - rv).powi(2) + (C2 - v2) * (r2 - C2 * t * t)).sqrt();
        let a_ref = vv.scale(v_ref / C2);

        // Retarded geometry from the closed-form delay.
        let w = [r[0] - v[0] * t, r[1] - v[1] * t, r[2] - v[2] * t];
        let tau = uniform_motion_delay(w, v);
        let t_r = t - tau;
        let r_vec = FieldVector::new(r[0] - v[0] * t_r, r[1] - v[1] * t_r, r[2] - v[2] * t_r);
        let beta = vv.scale(1.0 / SPEED_OF_LIGHT);

        let phi = lienard_wiechert_scalar_potential(q, &r_vec, &beta).unwrap();
        assert!(
            (phi - v_ref).abs() / v_ref < 1e-12,
            "φ = {phi}, ref {v_ref}"
        );
        let a = lienard_wiechert_vector_potential(q, &r_vec, &beta).unwrap();
        assert!((a - a_ref).magnitude() / a_ref.magnitude() < 1e-12);

        // The effective-distance functions agree when given R_eff = |R| − R⋅β.
        let r_eff = r_vec.magnitude() - r_vec.dot(&beta);
        let phi_eff = retarded_scalar_potential(q, r_eff).unwrap();
        assert!((phi_eff - v_ref).abs() / v_ref < 1e-12);
        let a_eff = retarded_vector_potential(q, &vv, r_eff).unwrap();
        assert!((a_eff - a_ref).magnitude() / a_ref.magnitude() < 1e-12);
        // …while the plain distance |R| is off by the factor 1 − n̂⋅β.
        let phi_plain = retarded_scalar_potential(q, r_vec.magnitude()).unwrap();
        assert!((phi_plain - v_ref).abs() / v_ref > 1e-2);
    }

    #[test]
    fn test_lienard_wiechert_potentials_errors() {
        let b = FieldVector::new(0.5, 0.0, 0.0);
        assert!(matches!(
            lienard_wiechert_scalar_potential(1.0, &FieldVector::zero(), &b),
            Err(BijliError::Singularity)
        ));
        let r = FieldVector::new(1.0, 0.0, 0.0);
        assert!(
            lienard_wiechert_scalar_potential(1.0, &r, &FieldVector::new(1.0, 0.0, 0.0)).is_err()
        );
        assert!(
            lienard_wiechert_vector_potential(1.0, &r, &FieldVector::new(f64::NAN, 0.0, 0.0))
                .is_err()
        );
        assert!(
            lienard_wiechert_vector_potential(1.0, &FieldVector::new(f64::NAN, 0.0, 0.0), &b)
                .is_err()
        );
        // Static limit: Coulomb potential, zero vector potential.
        let phi = lienard_wiechert_scalar_potential(1e-9, &r, &FieldVector::zero()).unwrap();
        assert!((phi - 1e-9 / (4.0 * std::f64::consts::PI * EPSILON_0)).abs() / phi < 1e-15);
    }

    #[test]
    fn test_retarded_potentials_reject_negative_distance() {
        assert!(retarded_scalar_potential(1e-6, -1.0).is_err());
        assert!(retarded_scalar_potential(1e-6, f64::NAN).is_err());
        let v = FieldVector::new(1e6, 0.0, 0.0);
        assert!(retarded_vector_potential(1e-6, &v, -1.0).is_err());
    }
}
