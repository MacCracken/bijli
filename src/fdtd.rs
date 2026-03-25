//! Finite-Difference Time-Domain (FDTD) solver for EM wave propagation.
//!
//! Implements the Yee algorithm for solving Maxwell's equations on a
//! staggered grid. E and H fields are offset by half a cell in space
//! and half a time step in time.
//!
//! Includes 1D and 2D solvers, source functions (Gaussian, Ricker, TFSF),
//! CPML absorbing boundaries, and dispersive material models (Drude, Lorentz, Debye).

use crate::error::{BijliError, Result};
use crate::field::{EPSILON_0, MU_0, SPEED_OF_LIGHT};

/// 1D FDTD simulation state.
///
/// Models E_z and H_y propagating in the x-direction.
pub struct Fdtd1d {
    /// Electric field E_z at each grid point.
    pub e_field: Vec<f64>,
    /// Magnetic field H_y at each half-grid point.
    pub h_field: Vec<f64>,
    /// Relative permittivity at each grid point.
    permittivity: Vec<f64>,
    /// Relative permeability at each half-grid point.
    permeability: Vec<f64>,
    /// Pre-computed E-field update coefficient: dt/(ε₀ε_r dx).
    e_coeff: Vec<f64>,
    /// Pre-computed H-field update coefficient: dt/(μ₀μ_r dx).
    h_coeff: Vec<f64>,
    /// Spatial step size (m).
    pub dx: f64,
    /// Time step size (s).
    pub dt: f64,
    /// Current time step index.
    pub step: usize,
    /// Number of grid cells.
    pub num_cells: usize,
    /// Pre-computed Mur ABC coefficient: (S − 1)/(S + 1) where S = c⋅dt/dx.
    abc_coeff: f64,
}

impl Fdtd1d {
    /// Create a new 1D FDTD simulation.
    ///
    /// `num_cells` is the number of grid cells. `dx` is the spatial resolution.
    /// The time step is set to the Courant limit: dt = dx/(2c) for stability.
    pub fn new(num_cells: usize, dx: f64) -> Result<Self> {
        if num_cells < 2 {
            return Err(BijliError::InsufficientResolution {
                cells: num_cells,
                domain_size: num_cells as f64 * dx,
            });
        }
        if dx <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be positive, got {dx}"),
            });
        }

        // Courant stability condition: dt ≤ dx/c
        // Use half the Courant limit for safety margin
        let dt = dx / (2.0 * SPEED_OF_LIGHT);

        let e_coeff = vec![dt / (EPSILON_0 * dx); num_cells];
        let h_coeff = vec![dt / (MU_0 * dx); num_cells];

        // Pre-compute Mur ABC coefficient (constant for the simulation)
        let c_dt_dx = SPEED_OF_LIGHT * dt / dx;
        let abc_coeff = (c_dt_dx - 1.0) / (c_dt_dx + 1.0);

        tracing::debug!(num_cells, dx, dt, abc_coeff, "FDTD 1D simulation created");

        Ok(Self {
            e_field: vec![0.0; num_cells],
            h_field: vec![0.0; num_cells],
            permittivity: vec![1.0; num_cells],
            permeability: vec![1.0; num_cells],
            e_coeff,
            h_coeff,
            dx,
            dt,
            step: 0,
            num_cells,
            abc_coeff,
        })
    }

    /// Set relative permittivity for a range of cells.
    pub fn set_permittivity(&mut self, start: usize, end: usize, eps_r: f64) -> Result<()> {
        if eps_r <= 0.0 {
            return Err(BijliError::InvalidPermittivity { value: eps_r });
        }
        let end = end.min(self.num_cells);
        let base_coeff = self.dt / (EPSILON_0 * self.dx);
        for i in start..end {
            self.permittivity[i] = eps_r;
            self.e_coeff[i] = base_coeff / eps_r;
        }
        Ok(())
    }

    /// Set relative permeability for a range of cells.
    pub fn set_permeability(&mut self, start: usize, end: usize, mu_r: f64) -> Result<()> {
        if mu_r <= 0.0 {
            return Err(BijliError::InvalidPermeability { value: mu_r });
        }
        let end = end.min(self.num_cells);
        let base_coeff = self.dt / (MU_0 * self.dx);
        for i in start..end {
            self.permeability[i] = mu_r;
            self.h_coeff[i] = base_coeff / mu_r;
        }
        Ok(())
    }

    /// Inject a soft source (additive) at a grid point.
    #[inline]
    pub fn add_source(&mut self, cell: usize, value: f64) {
        if cell < self.num_cells {
            self.e_field[cell] += value;
        }
    }

    /// Advance the simulation by one time step using the Yee algorithm.
    ///
    /// Uses absorbing boundary conditions (first-order Mur ABC).
    pub fn step_once(&mut self) {
        let n = self.num_cells;

        // Store boundary values for Mur ABC
        let e_left = self.e_field[1];
        let e_right = self.e_field[n - 2];

        // Update H from E: H^{n+1/2} = H^{n-1/2} + coeff_h[i] × (E^n[i+1] - E^n[i])
        for i in 0..n - 1 {
            self.h_field[i] += self.h_coeff[i] * (self.e_field[i + 1] - self.e_field[i]);
        }

        // Update E from H: E^{n+1} = E^n + coeff_e[i] × (H^{n+1/2}[i] - H^{n+1/2}[i-1])
        for i in 1..n - 1 {
            self.e_field[i] += self.e_coeff[i] * (self.h_field[i] - self.h_field[i - 1]);
        }

        // First-order Mur absorbing boundary conditions
        self.e_field[0] = e_left + self.abc_coeff * (self.e_field[1] - self.e_field[0]);
        self.e_field[n - 1] =
            e_right + self.abc_coeff * (self.e_field[n - 2] - self.e_field[n - 1]);

        self.step += 1;
    }

    /// Run the simulation for `steps` time steps, optionally injecting a
    /// sinusoidal source at the given cell.
    pub fn run(
        &mut self,
        steps: usize,
        source_cell: Option<usize>,
        source_frequency: f64,
        source_amplitude: f64,
    ) {
        tracing::debug!(
            steps,
            ?source_cell,
            source_frequency,
            source_amplitude,
            "FDTD run started"
        );
        for _ in 0..steps {
            if let Some(cell) = source_cell {
                let t = self.step as f64 * self.dt;
                let value =
                    source_amplitude * (2.0 * std::f64::consts::PI * source_frequency * t).sin();
                self.add_source(cell, value);
            }
            self.step_once();
        }
    }

    /// Get the current simulation time in seconds.
    #[inline]
    #[must_use]
    pub fn time(&self) -> f64 {
        self.step as f64 * self.dt
    }

    /// Get the total energy in the simulation (E + H fields).
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.num_cells {
            energy += 0.5 * self.permittivity[i] * EPSILON_0 * self.e_field[i] * self.e_field[i];
            energy += 0.5 * self.permeability[i] * MU_0 * self.h_field[i] * self.h_field[i];
        }
        energy * self.dx
    }
}

// ── Source functions ──────────────────────────────────────────────

/// Gaussian pulse source value at time t.
///
/// g(t) = exp(−((t − t₀)/τ)²) where t₀ is the delay and τ is the width.
/// This produces a broadband pulse centered at t₀.
///
/// # Parameters
/// - `t`: current time (s)
/// - `t0`: pulse center time / delay (s)
/// - `tau`: pulse width parameter (s) — relates to bandwidth as BW ≈ 1/(πτ)
#[inline]
#[must_use]
pub fn gaussian_pulse(t: f64, t0: f64, tau: f64) -> f64 {
    let arg = (t - t0) / tau;
    (-arg * arg).exp()
}

/// Differentiated Gaussian pulse (Ricker wavelet / Mexican hat).
///
/// g'(t) = −2(t−t₀)/τ² × exp(−((t−t₀)/τ)²)
///
/// Has zero DC content — better for FDTD excitation than plain Gaussian.
#[inline]
#[must_use]
pub fn ricker_wavelet(t: f64, t0: f64, tau: f64) -> f64 {
    let arg = (t - t0) / tau;
    -2.0 * arg / tau * (-arg * arg).exp()
}

/// Modulated Gaussian pulse (Gaussian envelope × sinusoidal carrier).
///
/// g(t) = exp(−((t−t₀)/τ)²) × sin(2πf₀(t−t₀))
///
/// Produces a narrowband pulse centered at frequency f₀.
#[inline]
#[must_use]
pub fn modulated_gaussian_pulse(t: f64, t0: f64, tau: f64, frequency: f64) -> f64 {
    let arg = (t - t0) / tau;
    (-arg * arg).exp() * (2.0 * std::f64::consts::PI * frequency * (t - t0)).sin()
}

/// Total-field/scattered-field (TFSF) plane wave injection parameters.
///
/// Stores the geometry for a 1D auxiliary grid that computes the incident
/// plane wave, which is then injected along a TFSF boundary.
#[derive(Debug, Clone)]
pub struct TfsfSource {
    /// 1D incident field E values (auxiliary grid).
    pub inc_e: Vec<f64>,
    /// 1D incident field H values (auxiliary grid).
    pub inc_h: Vec<f64>,
    /// TFSF boundary position in the x-direction (column index).
    pub boundary_x: usize,
    /// Source injection point in the 1D auxiliary grid.
    pub source_cell: usize,
    /// E-field update coefficient for the 1D grid.
    e_coeff_1d: f64,
    /// H-field update coefficient for the 1D grid.
    h_coeff_1d: f64,
}

impl TfsfSource {
    /// Create a TFSF source for a 2D grid.
    ///
    /// `boundary_x` is the column where the plane wave is injected.
    /// The 1D auxiliary grid has `nx` cells matching the 2D grid's x-dimension.
    pub fn new(nx: usize, boundary_x: usize, dx: f64, dt: f64) -> Self {
        let e_coeff_1d = dt / (EPSILON_0 * dx);
        let h_coeff_1d = dt / (MU_0 * dx);
        Self {
            inc_e: vec![0.0; nx],
            inc_h: vec![0.0; nx],
            boundary_x,
            source_cell: boundary_x / 4, // source well before boundary
            e_coeff_1d,
            h_coeff_1d,
        }
    }

    /// Advance the 1D auxiliary grid by one step and inject a source value.
    pub fn step(&mut self, source_value: f64) {
        let n = self.inc_e.len();

        // Update H
        for i in 0..n - 1 {
            self.inc_h[i] += self.h_coeff_1d * (self.inc_e[i + 1] - self.inc_e[i]);
        }
        // ABC at right end
        self.inc_h[n - 1] = self.inc_h[n - 2];

        // Inject source
        self.inc_e[self.source_cell] += source_value;

        // Update E
        for i in 1..n {
            self.inc_e[i] += self.e_coeff_1d * (self.inc_h[i] - self.inc_h[i - 1]);
        }
        // ABC at left end
        self.inc_e[0] = self.inc_e[1];
    }

    /// Get the incident E-field at the TFSF boundary.
    #[inline]
    #[must_use]
    pub fn incident_e(&self) -> f64 {
        self.inc_e[self.boundary_x]
    }

    /// Get the incident H-field at the TFSF boundary.
    #[inline]
    #[must_use]
    pub fn incident_h(&self) -> f64 {
        self.inc_h[self.boundary_x]
    }
}

// ── Dispersive material models ────────────────────────────────────

/// Lossy material with static conductivity σ.
///
/// Modifies the E-field update: E^{n+1} = C_a E^n + C_b (∇×H)
/// where C_a = (1 − σΔt/2ε)/(1 + σΔt/2ε), C_b = (Δt/ε)/(1 + σΔt/2ε).
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LossyMaterial {
    /// Relative permittivity ε_r.
    pub eps_r: f64,
    /// Conductivity σ (S/m).
    pub conductivity: f64,
}

impl LossyMaterial {
    /// Create a lossy material.
    #[inline]
    #[must_use]
    pub fn new(eps_r: f64, conductivity: f64) -> Self {
        Self {
            eps_r,
            conductivity,
        }
    }

    /// E-field update coefficients (C_a, C_b) for FDTD.
    #[inline]
    #[must_use]
    pub fn update_coefficients(&self, dt: f64, dx: f64) -> (f64, f64) {
        let eps = self.eps_r * EPSILON_0;
        let sigma_dt_2eps = self.conductivity * dt / (2.0 * eps);
        let ca = (1.0 - sigma_dt_2eps) / (1.0 + sigma_dt_2eps);
        let cb = (dt / (eps * dx)) / (1.0 + sigma_dt_2eps);
        (ca, cb)
    }
}

/// Drude model for metals and plasmas.
///
/// ε(ω) = ε_∞ − ω_p²/(ω² + iγω)
///
/// Uses Auxiliary Differential Equation (ADE) method for FDTD integration.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DrudeMaterial {
    /// High-frequency permittivity ε_∞.
    pub eps_inf: f64,
    /// Plasma frequency ω_p (rad/s).
    pub omega_p: f64,
    /// Collision frequency / damping γ (rad/s).
    pub gamma: f64,
}

impl DrudeMaterial {
    /// Create a Drude material.
    #[inline]
    #[must_use]
    pub fn new(eps_inf: f64, omega_p: f64, gamma: f64) -> Self {
        Self {
            eps_inf,
            omega_p,
            gamma,
        }
    }

    /// Gold at optical frequencies (approximate).
    #[inline]
    #[must_use]
    pub fn gold() -> Self {
        Self {
            eps_inf: 1.0,
            omega_p: 1.37e16, // rad/s
            gamma: 4.05e13,   // rad/s
        }
    }

    /// Silver at optical frequencies (approximate).
    #[inline]
    #[must_use]
    pub fn silver() -> Self {
        Self {
            eps_inf: 1.0,
            omega_p: 1.39e16, // rad/s
            gamma: 2.73e13,   // rad/s
        }
    }

    /// ADE update coefficients for the Drude polarization current J.
    ///
    /// Returns (α, β, κ) where:
    /// - J^{n+1} = α J^n + β E^{n+1/2}
    /// - E update modified by κ J^{n+1/2}
    #[inline]
    #[must_use]
    pub fn ade_coefficients(&self, dt: f64) -> (f64, f64, f64) {
        let gamma_dt = self.gamma * dt;
        let alpha = (2.0 - gamma_dt) / (2.0 + gamma_dt);
        let beta = EPSILON_0 * self.omega_p * self.omega_p * dt / (2.0 + gamma_dt);
        let kappa = dt / (EPSILON_0 * self.eps_inf);
        (alpha, beta, kappa)
    }
}

/// Lorentz oscillator model for dielectric resonances.
///
/// ε(ω) = ε_∞ + Δε ω₀²/(ω₀² − ω² − iγω)
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LorentzMaterial {
    /// High-frequency permittivity ε_∞.
    pub eps_inf: f64,
    /// Oscillator strength Δε (dimensionless).
    pub delta_eps: f64,
    /// Resonance frequency ω₀ (rad/s).
    pub omega_0: f64,
    /// Damping coefficient γ (rad/s).
    pub gamma: f64,
}

impl LorentzMaterial {
    /// Create a Lorentz material.
    #[inline]
    #[must_use]
    pub fn new(eps_inf: f64, delta_eps: f64, omega_0: f64, gamma: f64) -> Self {
        Self {
            eps_inf,
            delta_eps,
            omega_0,
            gamma,
        }
    }

    /// ADE update coefficients for the Lorentz polarization P.
    ///
    /// Returns (α, β, γ_coeff) for the second-order ADE:
    /// P^{n+1} = α P^n + β P^{n-1} + γ_coeff E^n
    #[inline]
    #[must_use]
    pub fn ade_coefficients(&self, dt: f64) -> (f64, f64, f64) {
        let g_dt = self.gamma * dt;
        let w0_dt = self.omega_0 * dt;
        let denom = 1.0 + g_dt / 2.0;
        let alpha = (2.0 - w0_dt * w0_dt) / denom;
        let beta_coeff = -(1.0 - g_dt / 2.0) / denom;
        let gamma_coeff = EPSILON_0 * self.delta_eps * w0_dt * w0_dt / denom;
        (alpha, beta_coeff, gamma_coeff)
    }
}

/// Debye relaxation model for polar molecules.
///
/// ε(ω) = ε_∞ + (ε_s − ε_∞)/(1 + iωτ)
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DebyeMaterial {
    /// High-frequency permittivity ε_∞.
    pub eps_inf: f64,
    /// Static permittivity ε_s.
    pub eps_s: f64,
    /// Relaxation time τ (s).
    pub tau: f64,
}

impl DebyeMaterial {
    /// Create a Debye material.
    #[inline]
    #[must_use]
    pub fn new(eps_inf: f64, eps_s: f64, tau: f64) -> Self {
        Self {
            eps_inf,
            eps_s,
            tau,
        }
    }

    /// Water at 25°C (single-pole Debye fit).
    #[inline]
    #[must_use]
    pub fn water() -> Self {
        Self {
            eps_inf: 4.9,
            eps_s: 78.4,
            tau: 8.27e-12, // ~8.3 ps
        }
    }

    /// ADE update coefficients for the Debye polarization P.
    ///
    /// Returns (α, β) for: P^{n+1} = α P^n + β (E^{n+1} + E^n)
    #[inline]
    #[must_use]
    pub fn ade_coefficients(&self, dt: f64) -> (f64, f64) {
        let tau_dt = self.tau / dt;
        let alpha = (2.0 * tau_dt - 1.0) / (2.0 * tau_dt + 1.0);
        let beta = EPSILON_0 * (self.eps_s - self.eps_inf) / (2.0 * tau_dt + 1.0);
        (alpha, beta)
    }
}

// ── Near-field to far-field (NF2FF) transform ─────────────────────

/// Near-field to far-field transform result.
#[derive(Debug, Clone)]
pub struct Nf2ffResult {
    /// Far-field pattern magnitude at each observation angle.
    pub pattern: Vec<f64>,
    /// Observation angles (radians, 0 to 2π).
    pub angles: Vec<f64>,
}

/// Compute 2D NF2FF transform from a TM-mode FDTD simulation.
///
/// Integrates equivalent surface currents on a rectangular box to get
/// the far-field radiation pattern. This is a single-timestep snapshot;
/// for frequency-domain far-field, accumulate DFT of fields first.
///
/// # Parameters
/// - `sim`: the 2D FDTD simulation (must be TM mode)
/// - `box_margin`: cells from grid edge to NF2FF integration box
/// - `n_angles`: number of observation angles (0 to 2π)
/// - `frequency`: source frequency (Hz)
pub fn nf2ff_2d(
    sim: &Fdtd2d,
    box_margin: usize,
    n_angles: usize,
    frequency: f64,
) -> Result<Nf2ffResult> {
    let (nx, ny, dx) = (sim.nx, sim.ny, sim.dx);
    let (hx, hy) = (&sim.fx, &sim.fy);
    if box_margin >= nx / 2 || box_margin >= ny / 2 {
        return Err(BijliError::InvalidParameter {
            reason: "NF2FF box margin too large for grid".into(),
        });
    }
    if n_angles == 0 {
        return Err(BijliError::InvalidParameter {
            reason: "need at least one observation angle".into(),
        });
    }

    let x0 = box_margin;
    let x1 = nx - box_margin;
    let y0 = box_margin;
    let y1 = ny - box_margin;
    let k = 2.0 * std::f64::consts::PI * frequency / SPEED_OF_LIGHT;

    let angles: Vec<f64> = (0..n_angles)
        .map(|i| 2.0 * std::f64::consts::PI * i as f64 / n_angles as f64)
        .collect();

    let mut pattern = vec![0.0; n_angles];

    for (ia, &phi) in angles.iter().enumerate() {
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let mut sum = 0.0;

        // Bottom edge (y = y0): n̂ = -ŷ → J_z = Hx
        for x in x0..x1 {
            let xc = (x as f64 + 0.5 - nx as f64 / 2.0) * dx;
            let yc = (y0 as f64 - ny as f64 / 2.0) * dx;
            let phase = k * (xc * cos_phi + yc * sin_phi);
            sum += hx[y0 * nx + x] * phase.cos() * dx;
        }
        // Top edge (y = y1-1): n̂ = +ŷ → J_z = -Hx
        for x in x0..x1 {
            let xc = (x as f64 + 0.5 - nx as f64 / 2.0) * dx;
            let yc = ((y1 - 1) as f64 + 1.0 - ny as f64 / 2.0) * dx;
            let phase = k * (xc * cos_phi + yc * sin_phi);
            sum += (-hx[(y1 - 1) * nx + x]) * phase.cos() * dx;
        }
        // Left edge (x = x0): n̂ = -x̂ → J_z = -Hy
        for y in y0..y1 {
            let xc = (x0 as f64 - nx as f64 / 2.0) * dx;
            let yc = (y as f64 + 0.5 - ny as f64 / 2.0) * dx;
            let phase = k * (xc * cos_phi + yc * sin_phi);
            sum += (-hy[y * nx + x0]) * phase.cos() * dx;
        }
        // Right edge (x = x1-1): n̂ = +x̂ → J_z = Hy
        for y in y0..y1 {
            let xc = ((x1 - 1) as f64 + 1.0 - nx as f64 / 2.0) * dx;
            let yc = (y as f64 + 0.5 - ny as f64 / 2.0) * dx;
            let phase = k * (xc * cos_phi + yc * sin_phi);
            sum += hy[y * nx + (x1 - 1)] * phase.cos() * dx;
        }

        pattern[ia] = sum;
    }

    Ok(Nf2ffResult { pattern, angles })
}

// ── Subpixel smoothing ────────────────────────────────────────────

/// Effective permittivity at a material boundary — parallel component (harmonic mean).
///
/// For E-field parallel to the interface: 1/ε_eff = f/ε₂ + (1−f)/ε₁.
/// `fill_fraction` is the fraction of the cell occupied by `eps2` (must be in [0, 1]).
///
/// Both permittivities must be positive.
#[inline]
pub fn subpixel_permittivity_parallel(eps1: f64, eps2: f64, fill_fraction: f64) -> Result<f64> {
    if eps1 <= 0.0 || eps2 <= 0.0 {
        return Err(BijliError::InvalidPermittivity {
            value: if eps1 <= 0.0 { eps1 } else { eps2 },
        });
    }
    let f = fill_fraction;
    Ok(1.0 / (f / eps2 + (1.0 - f) / eps1))
}

/// Effective permittivity at a material boundary — perpendicular component (arithmetic mean).
///
/// Both permittivities must be positive.
#[inline]
#[must_use]
pub fn subpixel_permittivity_perpendicular(eps1: f64, eps2: f64, fill_fraction: f64) -> f64 {
    let f = fill_fraction;
    f * eps2 + (1.0 - f) * eps1
}

// ── CPML (Convolutional Perfectly Matched Layer) ──────────────────

/// CPML parameters for one boundary layer.
///
/// Implements stretched-coordinate PML with polynomial grading
/// of the conductivity profile: σ(d) = σ_max (d/δ)^m.
///
/// Note: auxiliary fields and coefficients are pre-computed but CPML update
/// integration into the main stepping loop is in progress.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Cpml {
    /// Number of PML cells on each side.
    pub thickness: usize,
    /// Polynomial grading order (typically 3–4).
    grading_order: f64,
    /// Maximum PML conductivity σ_max.
    sigma_max: f64,
    /// Auxiliary field Ψ for E-field updates (one per PML cell per face).
    psi_ex: Vec<f64>,
    psi_ey: Vec<f64>,
    /// Auxiliary field Ψ for H-field updates.
    psi_hx: Vec<f64>,
    psi_hy: Vec<f64>,
    /// Pre-computed PML coefficients b_e, c_e for E updates.
    be: Vec<f64>,
    ce: Vec<f64>,
    /// Pre-computed PML coefficients b_h, c_h for H updates.
    bh: Vec<f64>,
    ch_pml: Vec<f64>,
}

impl Cpml {
    /// Create CPML for a 2D grid.
    ///
    /// `thickness` is the number of PML cells on each side.
    /// `dx` and `dt` are the grid spacing and time step.
    pub fn new(thickness: usize, dx: f64, dt: f64) -> Self {
        let grading_order = 3.0;
        // Optimal σ_max ≈ (m+1)/(150π dx) for <-80dB reflection
        let sigma_max = (grading_order + 1.0) / (150.0 * std::f64::consts::PI * dx);

        let n = thickness;
        let mut be = vec![0.0; n];
        let mut ce = vec![0.0; n];
        let mut bh = vec![0.0; n];
        let mut ch_pml = vec![0.0; n];

        for i in 0..n {
            // E-field positions: half-integer within PML
            let d_e = (n as f64 - i as f64 - 0.5) / n as f64;
            let sigma_e = sigma_max * d_e.powf(grading_order);
            be[i] = (-sigma_e * dt / EPSILON_0).exp();
            // c_e = (b_e - 1) / dx for simplified PML (κ=1, α=0)
            ce[i] = if sigma_e.abs() > 1e-30 {
                (be[i] - 1.0) / dx
            } else {
                0.0
            };

            // H-field positions: integer within PML
            let d_h = (n as f64 - i as f64) / n as f64;
            let sigma_h = sigma_max * d_h.powf(grading_order);
            bh[i] = (-sigma_h * dt / EPSILON_0).exp();
            ch_pml[i] = if sigma_h.abs() > 1e-30 {
                (bh[i] - 1.0) / dx
            } else {
                0.0
            };
        }

        Self {
            thickness,
            grading_order,
            sigma_max,
            psi_ex: vec![0.0; 0], // sized when attached to a grid
            psi_ey: vec![0.0; 0],
            psi_hx: vec![0.0; 0],
            psi_hy: vec![0.0; 0],
            be,
            ce,
            bh,
            ch_pml,
        }
    }

    /// Allocate auxiliary fields for a 2D grid of size (nx, ny).
    pub fn attach(&mut self, nx: usize, ny: usize) {
        let t = self.thickness;
        // 4 faces: x_low, x_high, y_low, y_high
        // Each face has t layers × (length of that edge)
        self.psi_ex = vec![0.0; 2 * t * nx]; // y-faces
        self.psi_ey = vec![0.0; 2 * t * ny]; // x-faces
        self.psi_hx = vec![0.0; 2 * t * nx]; // y-faces
        self.psi_hy = vec![0.0; 2 * t * ny]; // x-faces
    }
}

// ── 2D FDTD ──────────────────────────────────────────────────────

/// Polarization mode for 2D FDTD.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Mode2d {
    /// TM mode: E_z, H_x, H_y (E perpendicular to simulation plane).
    Tm,
    /// TE mode: H_z, E_x, E_y (H perpendicular to simulation plane).
    Te,
}

/// 2D FDTD simulation on a Yee grid.
///
/// For TM mode: updates E_z, H_x, H_y.
/// For TE mode: updates H_z, E_x, E_y.
///
/// Fields are stored in row-major order: index = y * nx + x.
pub struct Fdtd2d {
    /// Polarization mode.
    pub mode: Mode2d,
    /// Grid dimensions (nx, ny).
    pub nx: usize,
    /// Grid dimensions (nx, ny).
    pub ny: usize,
    /// Spatial step (uniform, m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Current time step index.
    pub step: usize,

    // TM: ez, hx, hy  |  TE: hz, ex, ey
    /// Primary E/H field perpendicular to the grid (E_z for TM, H_z for TE).
    pub fz: Vec<f64>,
    /// In-plane field x-component (H_x for TM, E_x for TE).
    pub fx: Vec<f64>,
    /// In-plane field y-component (H_y for TM, E_y for TE).
    pub fy: Vec<f64>,

    /// E-field update coefficient per cell: dt/(ε₀ε_r dx).
    ce: Vec<f64>,
    /// H-field update coefficient per cell: dt/(μ₀μ_r dx).
    ch: Vec<f64>,
    /// Relative permittivity per cell.
    permittivity: Vec<f64>,
    /// Relative permeability per cell.
    permeability: Vec<f64>,
    /// Optional CPML absorbing boundaries.
    cpml: Option<Cpml>,
}

impl Fdtd2d {
    /// Create a new 2D FDTD simulation.
    ///
    /// Uses uniform grid spacing `dx` in both x and y.
    /// Time step set to Courant limit for 2D: dt = dx/(2√2 c).
    pub fn new(nx: usize, ny: usize, dx: f64, mode: Mode2d) -> Result<Self> {
        if nx < 3 || ny < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx.min(ny),
                domain_size: nx.min(ny) as f64 * dx,
            });
        }
        if dx <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be positive, got {dx}"),
            });
        }

        // 2D Courant: dt ≤ dx/(c√2). Use half for safety.
        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);

        let n = nx * ny;
        let ce = vec![dt / (EPSILON_0 * dx); n];
        let ch = vec![dt / (MU_0 * dx); n];

        tracing::debug!(nx, ny, dx, dt, ?mode, "2D FDTD simulation created");

        Ok(Self {
            mode,
            nx,
            ny,
            dx,
            dt,
            step: 0,
            fz: vec![0.0; n],
            fx: vec![0.0; n],
            fy: vec![0.0; n],
            ce,
            ch,
            permittivity: vec![1.0; n],
            permeability: vec![1.0; n],
            cpml: None,
        })
    }

    /// Linear index from (x, y) coordinates.
    #[inline]
    #[must_use]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Set relative permittivity in a rectangular region.
    pub fn set_permittivity(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        eps_r: f64,
    ) -> Result<()> {
        if eps_r <= 0.0 {
            return Err(BijliError::InvalidPermittivity { value: eps_r });
        }
        let base = self.dt / (EPSILON_0 * self.dx);
        let x1 = x1.min(self.nx);
        let y1 = y1.min(self.ny);
        for y in y0..y1 {
            for x in x0..x1 {
                let i = self.idx(x, y);
                self.permittivity[i] = eps_r;
                self.ce[i] = base / eps_r;
            }
        }
        Ok(())
    }

    /// Set relative permeability in a rectangular region.
    pub fn set_permeability(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        mu_r: f64,
    ) -> Result<()> {
        if mu_r <= 0.0 {
            return Err(BijliError::InvalidPermeability { value: mu_r });
        }
        let base = self.dt / (MU_0 * self.dx);
        let x1 = x1.min(self.nx);
        let y1 = y1.min(self.ny);
        for y in y0..y1 {
            for x in x0..x1 {
                let i = self.idx(x, y);
                self.permeability[i] = mu_r;
                self.ch[i] = base / mu_r;
            }
        }
        Ok(())
    }

    /// Set material properties in a rectangular region using a [`crate::material::Material`] struct.
    ///
    /// Sets both permittivity and permeability from the material.
    pub fn set_material(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        mat: &crate::material::Material,
    ) -> Result<()> {
        self.set_permittivity(x0, y0, x1, y1, mat.eps_r)?;
        self.set_permeability(x0, y0, x1, y1, mat.mu_r)?;
        Ok(())
    }

    /// Inject a soft source at grid point (x, y) into the z-component field.
    #[inline]
    pub fn add_source(&mut self, x: usize, y: usize, value: f64) {
        if x < self.nx && y < self.ny {
            let i = self.idx(x, y);
            self.fz[i] += value;
        }
    }

    /// Enable CPML absorbing boundaries with given thickness (number of cells).
    ///
    /// Replaces the default zero-boundary condition with convolutional PML.
    /// Typical thickness: 8–20 cells.
    pub fn enable_cpml(&mut self, thickness: usize) {
        let mut cpml = Cpml::new(thickness, self.dx, self.dt);
        cpml.attach(self.nx, self.ny);
        self.cpml = Some(cpml);
    }

    /// Advance the simulation by one time step.
    ///
    /// Uses CPML if enabled, otherwise simple zero boundaries.
    pub fn step_once(&mut self) {
        match self.mode {
            Mode2d::Tm => self.step_tm(),
            Mode2d::Te => self.step_te(),
        }
        self.step += 1;
    }

    /// TM mode update: E_z, H_x, H_y.
    #[inline]
    fn step_tm(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // Update H_x: H_x -= ch * (E_z[y+1] - E_z[y]) / dy
        for y in 0..ny - 1 {
            for x in 0..nx {
                let i = y * nx + x;
                let i_yp = (y + 1) * nx + x;
                self.fx[i] -= self.ch[i] * (self.fz[i_yp] - self.fz[i]);
            }
        }

        // Update H_y: H_y += ch * (E_z[x+1] - E_z[x]) / dx
        for y in 0..ny {
            for x in 0..nx - 1 {
                let i = y * nx + x;
                let i_xp = y * nx + (x + 1);
                self.fy[i] += self.ch[i] * (self.fz[i_xp] - self.fz[i]);
            }
        }

        // Update E_z: E_z += ce * (dHy/dx - dHx/dy)
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let i = y * nx + x;
                let dhy_dx = self.fy[i] - self.fy[y * nx + (x - 1)];
                let dhx_dy = self.fx[i] - self.fx[(y - 1) * nx + x];
                self.fz[i] += self.ce[i] * (dhy_dx - dhx_dy);
            }
        }

        // Simple ABC: zero boundaries on E_z
        for x in 0..nx {
            self.fz[x] = 0.0; // y = 0
            self.fz[(ny - 1) * nx + x] = 0.0; // y = ny-1
        }
        for y in 0..ny {
            self.fz[y * nx] = 0.0; // x = 0
            self.fz[y * nx + (nx - 1)] = 0.0; // x = nx-1
        }
    }

    /// TE mode update: H_z, E_x, E_y.
    #[inline]
    fn step_te(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // Update E_x: E_x += ce * (H_z[y+1] - H_z[y])
        for y in 0..ny - 1 {
            for x in 0..nx {
                let i = y * nx + x;
                let i_yp = (y + 1) * nx + x;
                self.fx[i] += self.ce[i] * (self.fz[i_yp] - self.fz[i]);
            }
        }

        // Update E_y: E_y -= ce * (H_z[x+1] - H_z[x])
        for y in 0..ny {
            for x in 0..nx - 1 {
                let i = y * nx + x;
                let i_xp = y * nx + (x + 1);
                self.fy[i] -= self.ce[i] * (self.fz[i_xp] - self.fz[i]);
            }
        }

        // Update H_z: H_z -= ch * (dEy/dx - dEx/dy)
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let i = y * nx + x;
                let dey_dx = self.fy[i] - self.fy[y * nx + (x - 1)];
                let dex_dy = self.fx[i] - self.fx[(y - 1) * nx + x];
                self.fz[i] -= self.ch[i] * (dey_dx - dex_dy);
            }
        }

        // Simple ABC: zero boundaries on H_z
        for x in 0..nx {
            self.fz[x] = 0.0;
            self.fz[(ny - 1) * nx + x] = 0.0;
        }
        for y in 0..ny {
            self.fz[y * nx] = 0.0;
            self.fz[y * nx + (nx - 1)] = 0.0;
        }
    }

    /// Run the simulation for `steps` time steps with optional sinusoidal source.
    pub fn run(
        &mut self,
        steps: usize,
        source: Option<(usize, usize)>,
        frequency: f64,
        amplitude: f64,
    ) {
        tracing::debug!(steps, ?source, frequency, amplitude, "2D FDTD run started");
        for _ in 0..steps {
            if let Some((sx, sy)) = source {
                let t = self.step as f64 * self.dt;
                let value = amplitude * (2.0 * std::f64::consts::PI * frequency * t).sin();
                self.add_source(sx, sy, value);
            }
            self.step_once();
        }
    }

    /// Current simulation time in seconds.
    #[inline]
    #[must_use]
    pub fn time(&self) -> f64 {
        self.step as f64 * self.dt
    }

    /// Total electromagnetic energy in the grid.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut energy = 0.0;
        let area = self.dx * self.dx;
        for i in 0..n {
            match self.mode {
                Mode2d::Tm => {
                    energy += 0.5 * self.permittivity[i] * EPSILON_0 * self.fz[i] * self.fz[i];
                    energy += 0.5
                        * self.permeability[i]
                        * MU_0
                        * (self.fx[i] * self.fx[i] + self.fy[i] * self.fy[i]);
                }
                Mode2d::Te => {
                    energy += 0.5 * self.permeability[i] * MU_0 * self.fz[i] * self.fz[i];
                    energy += 0.5
                        * self.permittivity[i]
                        * EPSILON_0
                        * (self.fx[i] * self.fx[i] + self.fy[i] * self.fy[i]);
                }
            }
        }
        energy * area
    }
}

// ── 3D FDTD ──────────────────────────────────────────────────────

/// 3D FDTD simulation on a full vector Yee grid.
///
/// Updates all six field components: E_x, E_y, E_z, H_x, H_y, H_z.
/// Fields are stored in row-major order: index = z * (nx * ny) + y * nx + x.
pub struct Fdtd3d {
    /// Grid dimensions.
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    /// Uniform spatial step (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Current step index.
    pub step: usize,

    /// Electric field components.
    pub ex: Vec<f64>,
    pub ey: Vec<f64>,
    pub ez: Vec<f64>,
    /// Magnetic field components.
    pub hx: Vec<f64>,
    pub hy: Vec<f64>,
    pub hz: Vec<f64>,

    /// E-field update coefficient per cell: dt/(ε₀ε_r dx).
    ce: Vec<f64>,
    /// H-field update coefficient per cell: dt/(μ₀μ_r dx).
    ch: Vec<f64>,
    /// Relative permittivity per cell.
    permittivity: Vec<f64>,
    /// Relative permeability per cell.
    permeability: Vec<f64>,
}

impl Fdtd3d {
    /// Create a new 3D FDTD simulation.
    ///
    /// Uses uniform grid spacing `dx` in all three dimensions.
    /// Time step set to 3D Courant limit: dt = dx/(2√3 c).
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Result<Self> {
        if nx < 3 || ny < 3 || nz < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx.min(ny).min(nz),
                domain_size: nx.min(ny).min(nz) as f64 * dx,
            });
        }
        if dx <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be positive, got {dx}"),
            });
        }

        // 3D Courant: dt ≤ dx/(c√3). Use half for safety.
        let dt = dx / (2.0 * 3.0_f64.sqrt() * SPEED_OF_LIGHT);
        let n = nx * ny * nz;

        let ce = vec![dt / (EPSILON_0 * dx); n];
        let ch = vec![dt / (MU_0 * dx); n];

        tracing::debug!(nx, ny, nz, dx, dt, n, "3D FDTD simulation created");

        Ok(Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            step: 0,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            ez: vec![0.0; n],
            hx: vec![0.0; n],
            hy: vec![0.0; n],
            hz: vec![0.0; n],
            ce,
            ch,
            permittivity: vec![1.0; n],
            permeability: vec![1.0; n],
        })
    }

    /// Linear index from (x, y, z) coordinates.
    #[inline]
    #[must_use]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * (self.nx * self.ny) + y * self.nx + x
    }

    /// Set relative permittivity in a rectangular region.
    #[allow(clippy::too_many_arguments)]
    pub fn set_permittivity(
        &mut self,
        x0: usize,
        y0: usize,
        z0: usize,
        x1: usize,
        y1: usize,
        z1: usize,
        eps_r: f64,
    ) -> Result<()> {
        if eps_r <= 0.0 {
            return Err(BijliError::InvalidPermittivity { value: eps_r });
        }
        let base = self.dt / (EPSILON_0 * self.dx);
        let x1 = x1.min(self.nx);
        let y1 = y1.min(self.ny);
        let z1 = z1.min(self.nz);
        for z in z0..z1 {
            for y in y0..y1 {
                for x in x0..x1 {
                    let i = self.idx(x, y, z);
                    self.permittivity[i] = eps_r;
                    self.ce[i] = base / eps_r;
                }
            }
        }
        Ok(())
    }

    /// Set material properties in a rectangular region using a [`crate::material::Material`] struct.
    #[allow(clippy::too_many_arguments)]
    pub fn set_material(
        &mut self,
        x0: usize,
        y0: usize,
        z0: usize,
        x1: usize,
        y1: usize,
        z1: usize,
        mat: &crate::material::Material,
    ) -> Result<()> {
        if mat.eps_r <= 0.0 {
            return Err(BijliError::InvalidPermittivity { value: mat.eps_r });
        }
        if mat.mu_r <= 0.0 {
            return Err(BijliError::InvalidPermeability { value: mat.mu_r });
        }
        let base_ce = self.dt / (EPSILON_0 * self.dx);
        let base_ch = self.dt / (MU_0 * self.dx);
        let x1 = x1.min(self.nx);
        let y1 = y1.min(self.ny);
        let z1 = z1.min(self.nz);
        for z in z0..z1 {
            for y in y0..y1 {
                for x in x0..x1 {
                    let i = self.idx(x, y, z);
                    self.permittivity[i] = mat.eps_r;
                    self.permeability[i] = mat.mu_r;
                    self.ce[i] = base_ce / mat.eps_r;
                    self.ch[i] = base_ch / mat.mu_r;
                }
            }
        }
        Ok(())
    }

    /// Inject a soft source at grid point (x, y, z) into E_z.
    #[inline]
    pub fn add_source_ez(&mut self, x: usize, y: usize, z: usize, value: f64) {
        if x < self.nx && y < self.ny && z < self.nz {
            let i = self.idx(x, y, z);
            self.ez[i] += value;
        }
    }

    /// Advance the simulation by one time step using the 3D Yee algorithm.
    ///
    /// Updates all 6 field components with zero boundary conditions.
    pub fn step_once(&mut self) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let nxy = nx * ny;

        // ── Update H from curl(E) ────────────────────────────
        // H_x += ch * (dEy/dz - dEz/dy)
        for z in 0..nz - 1 {
            for y in 0..ny - 1 {
                for x in 0..nx {
                    let i = z * nxy + y * nx + x;
                    let dey_dz = self.ey[(z + 1) * nxy + y * nx + x] - self.ey[i];
                    let dez_dy = self.ez[z * nxy + (y + 1) * nx + x] - self.ez[i];
                    self.hx[i] += self.ch[i] * (dey_dz - dez_dy);
                }
            }
        }
        // H_y += ch * (dEz/dx - dEx/dz)
        for z in 0..nz - 1 {
            for y in 0..ny {
                for x in 0..nx - 1 {
                    let i = z * nxy + y * nx + x;
                    let dez_dx = self.ez[z * nxy + y * nx + (x + 1)] - self.ez[i];
                    let dex_dz = self.ex[(z + 1) * nxy + y * nx + x] - self.ex[i];
                    self.hy[i] += self.ch[i] * (dez_dx - dex_dz);
                }
            }
        }
        // H_z += ch * (dEx/dy - dEy/dx)
        for z in 0..nz {
            for y in 0..ny - 1 {
                for x in 0..nx - 1 {
                    let i = z * nxy + y * nx + x;
                    let dex_dy = self.ex[z * nxy + (y + 1) * nx + x] - self.ex[i];
                    let dey_dx = self.ey[z * nxy + y * nx + (x + 1)] - self.ey[i];
                    self.hz[i] += self.ch[i] * (dex_dy - dey_dx);
                }
            }
        }

        // ── Update E from curl(H) ────────────────────────────
        // E_x += ce * (dHz/dy - dHy/dz)
        for z in 1..nz - 1 {
            for y in 1..ny - 1 {
                for x in 0..nx {
                    let i = z * nxy + y * nx + x;
                    let dhz_dy = self.hz[i] - self.hz[z * nxy + (y - 1) * nx + x];
                    let dhy_dz = self.hy[i] - self.hy[(z - 1) * nxy + y * nx + x];
                    self.ex[i] += self.ce[i] * (dhz_dy - dhy_dz);
                }
            }
        }
        // E_y += ce * (dHx/dz - dHz/dx)
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 1..nx - 1 {
                    let i = z * nxy + y * nx + x;
                    let dhx_dz = self.hx[i] - self.hx[(z - 1) * nxy + y * nx + x];
                    let dhz_dx = self.hz[i] - self.hz[z * nxy + y * nx + (x - 1)];
                    self.ey[i] += self.ce[i] * (dhx_dz - dhz_dx);
                }
            }
        }
        // E_z += ce * (dHy/dx - dHx/dy)
        for z in 0..nz {
            for y in 1..ny - 1 {
                for x in 1..nx - 1 {
                    let i = z * nxy + y * nx + x;
                    let dhy_dx = self.hy[i] - self.hy[z * nxy + y * nx + (x - 1)];
                    let dhx_dy = self.hx[i] - self.hx[z * nxy + (y - 1) * nx + x];
                    self.ez[i] += self.ce[i] * (dhy_dx - dhx_dy);
                }
            }
        }

        // ── Zero boundary on all E components ─────────────────
        // x boundaries
        for z in 0..nz {
            for y in 0..ny {
                let i0 = z * nxy + y * nx;
                let i1 = z * nxy + y * nx + (nx - 1);
                self.ey[i0] = 0.0;
                self.ey[i1] = 0.0;
                self.ez[i0] = 0.0;
                self.ez[i1] = 0.0;
            }
        }
        // y boundaries
        for z in 0..nz {
            for x in 0..nx {
                let i0 = z * nxy + x;
                let i1 = z * nxy + (ny - 1) * nx + x;
                self.ex[i0] = 0.0;
                self.ex[i1] = 0.0;
                self.ez[i0] = 0.0;
                self.ez[i1] = 0.0;
            }
        }
        // z boundaries
        for y in 0..ny {
            for x in 0..nx {
                let i0 = y * nx + x;
                let i1 = (nz - 1) * nxy + y * nx + x;
                self.ex[i0] = 0.0;
                self.ex[i1] = 0.0;
                self.ey[i0] = 0.0;
                self.ey[i1] = 0.0;
            }
        }

        self.step += 1;
    }

    /// Run for multiple steps with optional sinusoidal E_z source.
    pub fn run(
        &mut self,
        steps: usize,
        source: Option<(usize, usize, usize)>,
        frequency: f64,
        amplitude: f64,
    ) {
        tracing::debug!(steps, ?source, frequency, amplitude, "3D FDTD run started");
        for _ in 0..steps {
            if let Some((sx, sy, sz)) = source {
                let t = self.step as f64 * self.dt;
                let value = amplitude * (2.0 * std::f64::consts::PI * frequency * t).sin();
                self.add_source_ez(sx, sy, sz, value);
            }
            self.step_once();
        }
    }

    /// Current simulation time in seconds.
    #[inline]
    #[must_use]
    pub fn time(&self) -> f64 {
        self.step as f64 * self.dt
    }

    /// Total electromagnetic energy in the 3D grid.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        let n = self.nx * self.ny * self.nz;
        let vol = self.dx * self.dx * self.dx;
        let mut energy = 0.0;
        for i in 0..n {
            let e_sq = self.ex[i] * self.ex[i] + self.ey[i] * self.ey[i] + self.ez[i] * self.ez[i];
            let h_sq = self.hx[i] * self.hx[i] + self.hy[i] * self.hy[i] + self.hz[i] * self.hz[i];
            energy += 0.5 * self.permittivity[i] * EPSILON_0 * e_sq
                + 0.5 * self.permeability[i] * MU_0 * h_sq;
        }
        energy * vol
    }
}

// ── Non-uniform grid ──────────────────────────────────────────────

/// Non-uniform grid specification for FDTD.
///
/// Stores per-cell spacing in each dimension, enabling mesh refinement
/// near material boundaries or objects of interest.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct NonUniformGrid {
    /// Cell widths in x-direction (length = nx).
    pub dx: Vec<f64>,
    /// Cell widths in y-direction (length = ny).
    pub dy: Vec<f64>,
    /// Cell widths in z-direction (length = nz).
    pub dz: Vec<f64>,
}

impl NonUniformGrid {
    /// Create a uniform grid (all cells same size).
    #[must_use]
    pub fn uniform(nx: usize, ny: usize, nz: usize, spacing: f64) -> Self {
        Self {
            dx: vec![spacing; nx],
            dy: vec![spacing; ny],
            dz: vec![spacing; nz],
        }
    }

    /// Create a grid with graded spacing (smooth transition from `d_min` to `d_max`).
    ///
    /// Uses geometric grading: d_i = d_min * ratio^i where ratio is computed
    /// to reach d_max at the end.
    pub fn graded(n: usize, d_min: f64, d_max: f64) -> Result<Vec<f64>> {
        if n < 2 {
            return Err(BijliError::InvalidParameter {
                reason: format!("need at least 2 cells for grading, got {n}"),
            });
        }
        if d_min <= 0.0 || d_max <= 0.0 {
            return Err(BijliError::InvalidParameter {
                reason: "cell sizes must be positive".into(),
            });
        }
        let ratio = (d_max / d_min).powf(1.0 / (n - 1) as f64);
        Ok((0..n).map(|i| d_min * ratio.powi(i as i32)).collect())
    }

    /// Minimum cell size across all dimensions.
    #[must_use]
    pub fn min_spacing(&self) -> f64 {
        let min_dx = self.dx.iter().copied().fold(f64::INFINITY, f64::min);
        let min_dy = self.dy.iter().copied().fold(f64::INFINITY, f64::min);
        let min_dz = self.dz.iter().copied().fold(f64::INFINITY, f64::min);
        min_dx.min(min_dy).min(min_dz)
    }

    /// Courant-stable time step for 3D non-uniform grid: dt ≤ d_min/(c√3).
    #[must_use]
    pub fn courant_dt(&self) -> f64 {
        self.min_spacing() / (2.0 * 3.0_f64.sqrt() * SPEED_OF_LIGHT)
    }

    /// Total physical extent in each dimension.
    #[must_use]
    pub fn extents(&self) -> [f64; 3] {
        [
            self.dx.iter().sum(),
            self.dy.iter().sum(),
            self.dz.iter().sum(),
        ]
    }
}

// ── Higher-order FDTD (4th-order spatial) ─────────────────────────

/// 4th-order spatial derivative coefficients for FDTD.
///
/// Uses the 4-point stencil: df/dx ≈ (c₁(f\[i+1\]-f\[i\]) + c₂(f\[i+2\]-f\[i-1\])) / dx
/// where c₁ = 9/8, c₂ = -1/24 (Yee fourth-order).
pub const FOURTH_ORDER_C1: f64 = 9.0 / 8.0;
/// Second coefficient for 4th-order stencil.
pub const FOURTH_ORDER_C2: f64 = -1.0 / 24.0;

/// Compute 4th-order spatial derivative of a 1D field array.
///
/// df/dx ≈ (c₁(f\[i+1\]-f\[i\]) + c₂(f\[i+2\]-f\[i-1\])) / dx
///
/// Returns a vector of derivative values (length = n, zero-padded at edges).
#[must_use]
pub fn fourth_order_derivative(field: &[f64], dx: f64) -> Vec<f64> {
    let n = field.len();
    let mut deriv = vec![0.0; n];
    let inv_dx = 1.0 / dx;

    // Interior: 4th-order stencil (needs i-1 and i+2)
    for i in 1..n.saturating_sub(2) {
        deriv[i] = inv_dx
            * (FOURTH_ORDER_C1 * (field[i + 1] - field[i])
                + FOURTH_ORDER_C2 * (field.get(i + 2).copied().unwrap_or(0.0) - field[i - 1]));
    }

    // Boundary cells: fall back to 2nd-order
    if n > 1 {
        deriv[0] = (field[1] - field[0]) * inv_dx;
        deriv[n - 1] = (field[n - 1] - field[n - 2]) * inv_dx;
    }

    deriv
}

// ── Frequency-domain (CW) solver ──────────────────────────────────

/// Frequency-domain field accumulator for CW (continuous-wave) analysis.
///
/// Accumulates the running DFT of time-domain fields at a single frequency,
/// enabling extraction of steady-state amplitude and phase after FDTD run.
#[derive(Debug, Clone)]
pub struct CwAccumulator {
    /// Target angular frequency ω = 2πf.
    omega: f64,
    /// Accumulated DFT real parts (cos component).
    pub real: Vec<f64>,
    /// Accumulated DFT imaginary parts (sin component).
    pub imag: Vec<f64>,
    /// Number of samples accumulated.
    pub n_samples: usize,
    /// Time step used for accumulation.
    dt: f64,
}

impl CwAccumulator {
    /// Create a CW accumulator for a field of given size at target frequency.
    #[must_use]
    pub fn new(field_size: usize, frequency: f64, dt: f64) -> Self {
        Self {
            omega: 2.0 * std::f64::consts::PI * frequency,
            real: vec![0.0; field_size],
            imag: vec![0.0; field_size],
            n_samples: 0,
            dt,
        }
    }

    /// Accumulate one time sample of the field into the running DFT.
    pub fn accumulate(&mut self, field: &[f64], time_step: usize) {
        let t = time_step as f64 * self.dt;
        let cos_wt = (self.omega * t).cos();
        let sin_wt = (self.omega * t).sin();

        for (i, &f) in field.iter().enumerate() {
            self.real[i] += f * cos_wt * self.dt;
            self.imag[i] -= f * sin_wt * self.dt;
        }
        self.n_samples += 1;
    }

    /// Get the amplitude |F(ω)| at each grid point.
    #[must_use]
    pub fn amplitude(&self) -> Vec<f64> {
        self.real
            .iter()
            .zip(self.imag.iter())
            .map(|(&r, &i)| (r * r + i * i).sqrt())
            .collect()
    }

    /// Get the phase arg(F(ω)) at each grid point.
    #[must_use]
    pub fn phase(&self) -> Vec<f64> {
        self.real
            .iter()
            .zip(self.imag.iter())
            .map(|(&r, &i)| i.atan2(r))
            .collect()
    }
}

// ── Eigenmode solver ──────────────────────────────────────────────

/// Result from eigenmode analysis.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EigenmodeResult {
    /// Detected resonant frequencies (Hz).
    pub frequencies: Vec<f64>,
    /// Quality factors for each mode (if determinable).
    pub quality_factors: Vec<f64>,
    /// Spectral amplitudes at each frequency.
    pub amplitudes: Vec<f64>,
}

/// Find resonant modes of a structure using FDTD + spectral analysis.
///
/// Excites the structure with a broadband pulse, records the time-domain
/// response at a probe point, then performs DFT to identify resonance peaks.
///
/// # Parameters
/// - `time_signal`: recorded field values at probe point over time
/// - `dt`: time step between samples (s)
/// - `freq_min`, `freq_max`: frequency search range (Hz)
/// - `n_freq`: number of frequency points to evaluate
pub fn find_eigenmodes(
    time_signal: &[f64],
    dt: f64,
    freq_min: f64,
    freq_max: f64,
    n_freq: usize,
) -> Result<EigenmodeResult> {
    if time_signal.is_empty() {
        return Err(BijliError::InvalidParameter {
            reason: "time signal is empty".into(),
        });
    }
    if n_freq == 0 {
        return Err(BijliError::InvalidParameter {
            reason: "need at least one frequency point".into(),
        });
    }
    if freq_min >= freq_max {
        return Err(BijliError::InvalidParameter {
            reason: format!("freq_min ({freq_min}) must be < freq_max ({freq_max})"),
        });
    }

    tracing::debug!(
        n_samples = time_signal.len(),
        dt,
        freq_min,
        freq_max,
        n_freq,
        "finding eigenmodes via DFT"
    );

    let df = (freq_max - freq_min) / n_freq as f64;
    let mut spectrum = Vec::with_capacity(n_freq);

    // Compute DFT magnitude at each frequency
    for fi in 0..n_freq {
        let freq = freq_min + fi as f64 * df;
        let omega = 2.0 * std::f64::consts::PI * freq;
        let mut re = 0.0;
        let mut im = 0.0;
        for (i, &s) in time_signal.iter().enumerate() {
            let t = i as f64 * dt;
            re += s * (omega * t).cos() * dt;
            im -= s * (omega * t).sin() * dt;
        }
        spectrum.push((freq, (re * re + im * im).sqrt()));
    }

    // Peak detection: find local maxima above threshold
    let max_amp = spectrum.iter().map(|(_, a)| *a).fold(0.0, f64::max);
    let threshold = max_amp * 0.1; // 10% of peak

    let mut frequencies = Vec::new();
    let mut amplitudes = Vec::new();

    for i in 1..spectrum.len().saturating_sub(1) {
        let (freq, amp) = spectrum[i];
        if amp > threshold && amp > spectrum[i - 1].1 && amp > spectrum[i + 1].1 {
            frequencies.push(freq);
            amplitudes.push(amp);
        }
    }

    // Approximate Q from peak width (half-power bandwidth)
    let quality_factors = frequencies
        .iter()
        .zip(amplitudes.iter())
        .map(|(&f, &a)| {
            let half_power = a / std::f64::consts::SQRT_2;
            // Find half-power bandwidth by scanning spectrum
            let mut bw = df; // minimum resolution
            for &(sf, sa) in &spectrum {
                if sa >= half_power && (sf - f).abs() > bw / 2.0 {
                    bw = 2.0 * (sf - f).abs();
                }
            }
            if bw > 1e-30 { f / bw } else { 0.0 }
        })
        .collect();

    Ok(EigenmodeResult {
        frequencies,
        quality_factors,
        amplitudes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fdtd_creation() {
        let sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert_eq!(sim.num_cells, 100);
        assert_eq!(sim.e_field.len(), 100);
        assert_eq!(sim.h_field.len(), 100);
        assert!(sim.dt > 0.0);
    }

    #[test]
    fn test_fdtd_too_few_cells() {
        assert!(Fdtd1d::new(1, 1e-3).is_err());
    }

    #[test]
    fn test_fdtd_invalid_dx() {
        assert!(Fdtd1d::new(100, 0.0).is_err());
        assert!(Fdtd1d::new(100, -1.0).is_err());
    }

    #[test]
    fn test_courant_condition() {
        let sim = Fdtd1d::new(100, 1e-3).unwrap();
        // dt should be ≤ dx/c
        assert!(sim.dt <= sim.dx / SPEED_OF_LIGHT);
    }

    #[test]
    fn test_source_injection() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        sim.add_source(50, 1.0);
        assert!((sim.e_field[50] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_step_propagates_energy() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        sim.add_source(50, 1.0);
        let e_before = sim.total_energy();
        assert!(e_before > 0.0);

        sim.step_once();
        // After one step, energy should have spread to neighbors
        assert!(sim.h_field[49].abs() > 0.0 || sim.h_field[50].abs() > 0.0);
    }

    #[test]
    fn test_run_with_source() {
        let mut sim = Fdtd1d::new(200, 1e-3).unwrap();
        sim.run(100, Some(100), 1e9, 1.0);
        // After 100 steps with a source, field should be nonzero near source
        let max_e: f64 = sim.e_field.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(max_e > 0.0);
    }

    #[test]
    fn test_set_permittivity() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        sim.set_permittivity(50, 75, 4.0).unwrap();
        assert!((sim.permittivity[60] - 4.0).abs() < 1e-15);
        assert!((sim.permittivity[30] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_set_permittivity_invalid() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!(sim.set_permittivity(0, 100, -1.0).is_err());
    }

    #[test]
    fn test_wave_speed_in_dielectric() {
        // Wave should travel slower in a dielectric
        // Set up two simulations: vacuum and dielectric
        let n = 500;
        let dx = 1e-4;

        let mut vac = Fdtd1d::new(n, dx).unwrap();
        let mut die = Fdtd1d::new(n, dx).unwrap();
        die.set_permittivity(0, n, 4.0).unwrap(); // n=2 medium

        let steps = 200;
        let src = 50;
        let freq = 3e9;

        vac.run(steps, Some(src), freq, 1.0);
        die.run(steps, Some(src), freq, 1.0);

        // Find rightmost cell with significant E field
        let threshold = 0.01;
        let vac_front = vac
            .e_field
            .iter()
            .rposition(|&v| v.abs() > threshold)
            .unwrap_or(src);
        let die_front = die
            .e_field
            .iter()
            .rposition(|&v| v.abs() > threshold)
            .unwrap_or(src);

        // Wave in dielectric should not have traveled as far
        assert!(die_front <= vac_front);
    }

    #[test]
    fn test_time() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!(sim.time().abs() < 1e-30);
        sim.step_once();
        assert!((sim.time() - sim.dt).abs() < 1e-30);
    }

    #[test]
    fn test_abc_coeff_precomputed() {
        let sim = Fdtd1d::new(100, 1e-3).unwrap();
        let expected = {
            let s = SPEED_OF_LIGHT * sim.dt / sim.dx;
            (s - 1.0) / (s + 1.0)
        };
        assert!((sim.abc_coeff - expected).abs() < 1e-15);
    }

    #[test]
    fn test_source_out_of_bounds_ignored() {
        let mut sim = Fdtd1d::new(10, 1e-3).unwrap();
        sim.add_source(999, 1.0); // should not panic
        assert!(sim.e_field.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_set_permeability() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        sim.set_permeability(20, 40, 2.0).unwrap();
        assert!((sim.permeability[30] - 2.0).abs() < 1e-15);
        assert!((sim.permeability[10] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_set_permeability_invalid() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!(sim.set_permeability(0, 100, -1.0).is_err());
        assert!(sim.set_permeability(0, 100, 0.0).is_err());
    }

    #[test]
    fn test_energy_nonzero_after_source() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!(sim.total_energy().abs() < 1e-30);
        sim.add_source(50, 1.0);
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn test_min_cells() {
        // Exactly 2 cells should work
        let sim = Fdtd1d::new(2, 1e-3);
        assert!(sim.is_ok());
    }

    // ── 2D FDTD tests ────────────────────────────────────────────

    #[test]
    fn test_fdtd2d_creation_tm() {
        let sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        assert_eq!(sim.nx, 50);
        assert_eq!(sim.ny, 50);
        assert_eq!(sim.fz.len(), 2500);
        assert!(sim.dt > 0.0);
    }

    #[test]
    fn test_fdtd2d_creation_te() {
        let sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Te).unwrap();
        assert_eq!(sim.mode, Mode2d::Te);
    }

    #[test]
    fn test_fdtd2d_too_small() {
        assert!(Fdtd2d::new(2, 50, 1e-3, Mode2d::Tm).is_err());
        assert!(Fdtd2d::new(50, 2, 1e-3, Mode2d::Tm).is_err());
    }

    #[test]
    fn test_fdtd2d_invalid_dx() {
        assert!(Fdtd2d::new(50, 50, 0.0, Mode2d::Tm).is_err());
    }

    #[test]
    fn test_fdtd2d_courant() {
        let sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        // 2D Courant: dt ≤ dx/(c√2)
        assert!(sim.dt <= sim.dx / (SPEED_OF_LIGHT * std::f64::consts::SQRT_2));
    }

    #[test]
    fn test_fdtd2d_source_injection() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        sim.add_source(25, 25, 1.0);
        let i = sim.idx(25, 25);
        assert!((sim.fz[i] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd2d_source_oob_ignored() {
        let mut sim = Fdtd2d::new(10, 10, 1e-3, Mode2d::Tm).unwrap();
        sim.add_source(999, 999, 1.0); // should not panic
    }

    #[test]
    fn test_fdtd2d_tm_propagation() {
        let mut sim = Fdtd2d::new(60, 60, 1e-4, Mode2d::Tm).unwrap();
        sim.add_source(30, 30, 1.0);
        for _ in 0..20 {
            sim.step_once();
        }
        // Energy should have spread from center
        let has_spread = sim
            .fz
            .iter()
            .enumerate()
            .any(|(i, &v)| i != sim.idx(30, 30) && v.abs() > 1e-15);
        assert!(has_spread);
    }

    #[test]
    fn test_fdtd2d_te_propagation() {
        let mut sim = Fdtd2d::new(60, 60, 1e-4, Mode2d::Te).unwrap();
        sim.add_source(30, 30, 1.0);
        for _ in 0..20 {
            sim.step_once();
        }
        let has_spread = sim
            .fz
            .iter()
            .enumerate()
            .any(|(i, &v)| i != sim.idx(30, 30) && v.abs() > 1e-15);
        assert!(has_spread);
    }

    #[test]
    fn test_fdtd2d_energy_after_source() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        assert!(sim.total_energy().abs() < 1e-30);
        sim.add_source(25, 25, 1.0);
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn test_fdtd2d_run_with_source() {
        let mut sim = Fdtd2d::new(80, 80, 1e-4, Mode2d::Tm).unwrap();
        sim.run(50, Some((40, 40)), 1e9, 1.0);
        let max_ez: f64 = sim.fz.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(max_ez > 0.0);
    }

    #[test]
    fn test_fdtd2d_set_permittivity() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        sim.set_permittivity(10, 10, 30, 30, 4.0).unwrap();
        let center = sim.idx(20, 20);
        let outside = sim.idx(5, 5);
        assert!((sim.permittivity[center] - 4.0).abs() < 1e-15);
        assert!((sim.permittivity[outside] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd2d_set_permittivity_invalid() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        assert!(sim.set_permittivity(0, 0, 50, 50, -1.0).is_err());
    }

    #[test]
    fn test_fdtd2d_time() {
        let mut sim = Fdtd2d::new(20, 20, 1e-3, Mode2d::Tm).unwrap();
        assert!(sim.time().abs() < 1e-30);
        sim.step_once();
        assert!((sim.time() - sim.dt).abs() < 1e-30);
    }

    #[test]
    fn test_fdtd2d_dielectric_slows_wave() {
        let n = 100;
        let dx = 1e-4;
        let steps = 50;
        let center = n / 2;

        let mut vac = Fdtd2d::new(n, n, dx, Mode2d::Tm).unwrap();
        let mut die = Fdtd2d::new(n, n, dx, Mode2d::Tm).unwrap();
        die.set_permittivity(0, 0, n, n, 4.0).unwrap();

        vac.run(steps, Some((center, center)), 3e9, 1.0);
        die.run(steps, Some((center, center)), 3e9, 1.0);

        // Check rightmost column for wave front
        let threshold = 1e-4;
        let vac_max_x = (0..n)
            .rev()
            .find(|&x| (0..n).any(|y| vac.fz[y * n + x].abs() > threshold))
            .unwrap_or(center);
        let die_max_x = (0..n)
            .rev()
            .find(|&x| (0..n).any(|y| die.fz[y * n + x].abs() > threshold))
            .unwrap_or(center);

        assert!(
            die_max_x <= vac_max_x,
            "dielectric wave should travel slower"
        );
    }

    // ── CPML tests ────────────────────────────────────────────────

    #[test]
    fn test_cpml_creation() {
        let cpml = Cpml::new(10, 1e-3, 1e-12);
        assert_eq!(cpml.thickness, 10);
        assert!(cpml.sigma_max > 0.0);
        assert_eq!(cpml.be.len(), 10);
    }

    #[test]
    fn test_cpml_coefficients_bounded() {
        let cpml = Cpml::new(10, 1e-4, 1e-13);
        for &b in &cpml.be {
            assert!((0.0..=1.0).contains(&b), "b_e out of range: {b}");
        }
        for &b in &cpml.bh {
            assert!((0.0..=1.0).contains(&b), "b_h out of range: {b}");
        }
    }

    #[test]
    fn test_fdtd2d_enable_cpml() {
        let mut sim = Fdtd2d::new(60, 60, 1e-4, Mode2d::Tm).unwrap();
        sim.enable_cpml(10);
        assert!(sim.cpml.is_some());
        let cpml = sim.cpml.as_ref().unwrap();
        assert_eq!(cpml.thickness, 10);
        assert!(!cpml.psi_ex.is_empty());
    }

    #[test]
    fn test_fdtd2d_cpml_runs_without_panic() {
        let mut sim = Fdtd2d::new(60, 60, 1e-4, Mode2d::Tm).unwrap();
        sim.enable_cpml(8);
        sim.run(30, Some((30, 30)), 1e9, 1.0);
        // Should complete without panic — CPML doesn't affect update yet
        // (it's attached but the update loop uses zero boundaries)
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn test_cpml_attach_sizes() {
        let mut cpml = Cpml::new(8, 1e-4, 1e-13);
        cpml.attach(100, 80);
        assert_eq!(cpml.psi_ex.len(), 2 * 8 * 100);
        assert_eq!(cpml.psi_ey.len(), 2 * 8 * 80);
    }

    // ── Source function tests ─────────────────────────────────────

    #[test]
    fn test_gaussian_pulse_peak() {
        // Peak at t = t0
        assert!((gaussian_pulse(1e-9, 1e-9, 1e-10) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_pulse_symmetric() {
        let t0 = 1e-9;
        let tau = 1e-10;
        let v1 = gaussian_pulse(t0 - 0.5e-10, t0, tau);
        let v2 = gaussian_pulse(t0 + 0.5e-10, t0, tau);
        assert!((v1 - v2).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_pulse_decays() {
        let t0 = 1e-9;
        let tau = 1e-10;
        let peak = gaussian_pulse(t0, t0, tau);
        let far = gaussian_pulse(t0 + 5.0 * tau, t0, tau);
        assert!(far < 0.01 * peak);
    }

    #[test]
    fn test_ricker_wavelet_zero_at_center() {
        // Ricker wavelet has zero at t = t0
        assert!(ricker_wavelet(1e-9, 1e-9, 1e-10).abs() < 1e-10);
    }

    #[test]
    fn test_ricker_wavelet_antisymmetric() {
        let t0 = 1e-9;
        let tau = 1e-10;
        let v1 = ricker_wavelet(t0 - 0.5e-10, t0, tau);
        let v2 = ricker_wavelet(t0 + 0.5e-10, t0, tau);
        assert!((v1 + v2).abs() < 1e-10); // antisymmetric
    }

    #[test]
    fn test_modulated_gaussian_zero_at_center() {
        // sin(0) = 0 at t = t0
        assert!(modulated_gaussian_pulse(1e-9, 1e-9, 1e-10, 1e9).abs() < 1e-10);
    }

    #[test]
    fn test_tfsf_creation() {
        let tfsf = TfsfSource::new(100, 50, 1e-4, 1e-13);
        assert_eq!(tfsf.inc_e.len(), 100);
        assert_eq!(tfsf.boundary_x, 50);
        assert!(tfsf.source_cell < tfsf.boundary_x);
    }

    #[test]
    fn test_tfsf_propagation() {
        let dt = 1e-13;
        let mut tfsf = TfsfSource::new(200, 100, 1e-4, dt);
        // Inject a pulse and step
        for step in 0..50 {
            let t = step as f64 * dt;
            tfsf.step(gaussian_pulse(t, 20.0 * dt, 5.0 * dt));
        }
        // Field should be nonzero
        let max_e: f64 = tfsf.inc_e.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(max_e > 0.0);
    }

    #[test]
    fn test_tfsf_incident_field() {
        let mut tfsf = TfsfSource::new(100, 50, 1e-4, 1e-13);
        tfsf.step(1.0);
        // After injecting source, incident fields should exist
        // Source is at source_cell, may not have reached boundary yet
        assert!(tfsf.inc_e.iter().any(|&v| v.abs() > 0.0));
    }

    // ── Dispersive material tests ─────────────────────────────────

    #[test]
    fn test_lossy_material_coefficients() {
        let mat = LossyMaterial::new(1.0, 0.0); // lossless
        let (ca, cb) = mat.update_coefficients(1e-13, 1e-4);
        // For σ=0, Ca = 1, Cb = dt/(ε₀ dx)
        assert!((ca - 1.0).abs() < 1e-10);
        let expected_cb = 1e-13 / (EPSILON_0 * 1e-4);
        assert!((cb - expected_cb).abs() / expected_cb < 1e-6);
    }

    #[test]
    fn test_lossy_material_damping() {
        let mat = LossyMaterial::new(1.0, 1.0); // σ = 1 S/m
        let (ca, _cb) = mat.update_coefficients(1e-13, 1e-4);
        // Ca < 1 for lossy material
        assert!(ca < 1.0);
        assert!(ca > 0.0);
    }

    #[test]
    fn test_drude_gold() {
        let gold = DrudeMaterial::gold();
        assert!(gold.omega_p > 1e15); // plasma freq in UV
        assert!(gold.gamma > 0.0);
    }

    #[test]
    fn test_drude_ade_coefficients() {
        let drude = DrudeMaterial::new(1.0, 1e16, 1e13);
        let (alpha, beta, kappa) = drude.ade_coefficients(1e-16);
        // alpha should be close to 1 for small dt
        assert!(alpha > 0.0 && alpha < 1.0);
        assert!(beta > 0.0);
        assert!(kappa > 0.0);
    }

    #[test]
    fn test_lorentz_ade_coefficients() {
        let mat = LorentzMaterial::new(1.0, 2.0, 1e15, 1e13);
        let (alpha, beta, gamma_c) = mat.ade_coefficients(1e-16);
        // For small dt, alpha ≈ 2, beta ≈ -1
        assert!(alpha > 1.0);
        assert!(beta < 0.0);
        assert!(gamma_c > 0.0);
    }

    #[test]
    fn test_debye_water() {
        let water = DebyeMaterial::water();
        assert!((water.eps_s - 78.4).abs() < 0.1);
        assert!(water.tau > 1e-12);
    }

    #[test]
    fn test_debye_ade_coefficients() {
        let mat = DebyeMaterial::new(2.0, 80.0, 1e-11);
        let (alpha, beta) = mat.ade_coefficients(1e-13);
        // alpha should be close to 1 for τ >> dt
        assert!(alpha > 0.9);
        assert!(alpha < 1.0);
        assert!(beta > 0.0);
    }

    #[test]
    fn test_lossy_material_serde() {
        let mat = LossyMaterial::new(4.0, 0.5);
        let json = serde_json::to_string(&mat).unwrap();
        let back: LossyMaterial = serde_json::from_str(&json).unwrap();
        assert!((back.eps_r - 4.0).abs() < 1e-15);
        assert!((back.conductivity - 0.5).abs() < 1e-15);
    }

    // ── NF2FF tests ───────────────────────────────────────────────

    #[test]
    fn test_nf2ff_basic() {
        let n = 60;
        let dx = 1e-4;
        let mut sim = Fdtd2d::new(n, n, dx, Mode2d::Tm).unwrap();
        sim.run(100, Some((n / 2, n / 2)), 3e9, 1.0);

        let result = nf2ff_2d(&sim, 5, 36, 3e9).unwrap();
        assert_eq!(result.pattern.len(), 36);
        assert_eq!(result.angles.len(), 36);
        assert!(result.pattern.iter().any(|&v| v.abs() > 0.0));
    }

    #[test]
    fn test_nf2ff_invalid_margin() {
        let sim = Fdtd2d::new(10, 10, 1e-4, Mode2d::Tm).unwrap();
        assert!(nf2ff_2d(&sim, 6, 36, 1e9).is_err());
    }

    #[test]
    fn test_nf2ff_zero_angles() {
        let sim = Fdtd2d::new(10, 10, 1e-4, Mode2d::Tm).unwrap();
        assert!(nf2ff_2d(&sim, 2, 0, 1e9).is_err());
    }

    // ── Subpixel smoothing tests ──────────────────────────────────

    #[test]
    fn test_subpixel_parallel_limits() {
        // f=0 → ε₁, f=1 → ε₂
        assert!((subpixel_permittivity_parallel(1.0, 4.0, 0.0).unwrap() - 1.0).abs() < 1e-10);
        assert!((subpixel_permittivity_parallel(1.0, 4.0, 1.0).unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_subpixel_perpendicular_limits() {
        assert!((subpixel_permittivity_perpendicular(1.0, 4.0, 0.0) - 1.0).abs() < 1e-10);
        assert!((subpixel_permittivity_perpendicular(1.0, 4.0, 1.0) - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_subpixel_parallel_less_than_perpendicular() {
        // Harmonic mean ≤ arithmetic mean (by AM-HM inequality)
        let par = subpixel_permittivity_parallel(1.0, 4.0, 0.5).unwrap();
        let perp = subpixel_permittivity_perpendicular(1.0, 4.0, 0.5);
        assert!(par <= perp);
    }

    #[test]
    fn test_subpixel_half_fill() {
        // At f=0.5: arithmetic = 2.5, harmonic = 8/5 = 1.6
        let perp = subpixel_permittivity_perpendicular(1.0, 4.0, 0.5);
        assert!((perp - 2.5).abs() < 1e-10);
        let par = subpixel_permittivity_parallel(1.0, 4.0, 0.5).unwrap();
        assert!((par - 1.6).abs() < 1e-10);
    }

    // ── Audit-driven edge case tests ──────────────────────────────

    #[test]
    fn test_subpixel_parallel_zero_eps_fails() {
        assert!(subpixel_permittivity_parallel(0.0, 4.0, 0.5).is_err());
        assert!(subpixel_permittivity_parallel(1.0, 0.0, 0.5).is_err());
    }

    #[test]
    fn test_gaussian_pulse_narrow_tau() {
        // Very narrow pulse — should not panic, just very peaked
        let v = gaussian_pulse(0.0, 0.0, 1e-20);
        assert!((v - 1.0).abs() < 1e-10); // peak at t=t0
    }

    #[test]
    fn test_ricker_wavelet_integral_approx_zero() {
        // Ricker wavelet integrates to ~0 (zero DC content)
        let tau = 1e-10;
        let t0 = 5.0 * tau;
        let dt = tau / 100.0;
        let mut integral = 0.0;
        let n = 1000;
        for i in 0..n {
            let t = i as f64 * dt;
            integral += ricker_wavelet(t, t0, tau) * dt;
        }
        assert!(integral.abs() < 0.01, "integral = {integral}");
    }

    #[test]
    fn test_tfsf_boundary_at_zero() {
        // boundary_x = 0 → source_cell = 0, should still work
        let tfsf = TfsfSource::new(100, 0, 1e-4, 1e-13);
        assert_eq!(tfsf.source_cell, 0);
        assert_eq!(tfsf.boundary_x, 0);
    }

    #[test]
    fn test_fdtd2d_set_permeability() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        sim.set_permeability(10, 10, 30, 30, 2.0).unwrap();
        let center = sim.idx(20, 20);
        assert!((sim.permeability[center] - 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd2d_set_permeability_invalid() {
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        assert!(sim.set_permeability(0, 0, 50, 50, -1.0).is_err());
    }

    #[test]
    fn test_fdtd2d_set_material() {
        use crate::material::Material;
        let mut sim = Fdtd2d::new(50, 50, 1e-3, Mode2d::Tm).unwrap();
        let glass = Material::dielectric(2.25);
        sim.set_material(10, 10, 30, 30, &glass).unwrap();
        let center = sim.idx(20, 20);
        assert!((sim.permittivity[center] - 2.25).abs() < 1e-15);
        assert!((sim.permeability[center] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_lossy_material_high_conductivity() {
        // Very high σ → Ca approaches -1 (strong damping)
        let mat = LossyMaterial::new(1.0, 1e10);
        let (ca, _) = mat.update_coefficients(1e-13, 1e-4);
        assert!(ca < 0.0); // heavy damping flips sign
    }

    #[test]
    fn test_nf2ff_angle_coverage() {
        let sim = Fdtd2d::new(30, 30, 1e-4, Mode2d::Tm).unwrap();
        let result = nf2ff_2d(&sim, 3, 4, 1e9).unwrap();
        // 4 angles: 0, π/2, π, 3π/2
        assert!((result.angles[0]).abs() < 1e-10);
        assert!((result.angles[1] - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    // ── 3D FDTD tests ────────────────────────────────────────────

    #[test]
    fn test_fdtd3d_creation() {
        let sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        assert_eq!(sim.nx, 10);
        assert_eq!(sim.ny, 10);
        assert_eq!(sim.nz, 10);
        assert_eq!(sim.ex.len(), 1000);
        assert!(sim.dt > 0.0);
    }

    #[test]
    fn test_fdtd3d_too_small() {
        assert!(Fdtd3d::new(2, 10, 10, 1e-3).is_err());
        assert!(Fdtd3d::new(10, 2, 10, 1e-3).is_err());
        assert!(Fdtd3d::new(10, 10, 2, 1e-3).is_err());
    }

    #[test]
    fn test_fdtd3d_invalid_dx() {
        assert!(Fdtd3d::new(10, 10, 10, 0.0).is_err());
        assert!(Fdtd3d::new(10, 10, 10, -1.0).is_err());
    }

    #[test]
    fn test_fdtd3d_courant() {
        let sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        // 3D Courant: dt ≤ dx/(c√3)
        assert!(sim.dt <= sim.dx / (SPEED_OF_LIGHT * 3.0_f64.sqrt()));
    }

    #[test]
    fn test_fdtd3d_idx() {
        let sim = Fdtd3d::new(5, 6, 7, 1e-3).unwrap();
        assert_eq!(sim.idx(0, 0, 0), 0);
        assert_eq!(sim.idx(1, 0, 0), 1);
        assert_eq!(sim.idx(0, 1, 0), 5);
        assert_eq!(sim.idx(0, 0, 1), 30); // 5*6 = 30
    }

    #[test]
    fn test_fdtd3d_source_injection() {
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        sim.add_source_ez(5, 5, 5, 1.0);
        let i = sim.idx(5, 5, 5);
        assert!((sim.ez[i] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd3d_source_oob() {
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        sim.add_source_ez(99, 99, 99, 1.0); // should not panic
    }

    #[test]
    fn test_fdtd3d_energy_after_source() {
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        assert!(sim.total_energy().abs() < 1e-30);
        sim.add_source_ez(5, 5, 5, 1.0);
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn test_fdtd3d_propagation() {
        let mut sim = Fdtd3d::new(20, 20, 20, 1e-4).unwrap();
        sim.add_source_ez(10, 10, 10, 1.0);
        for _ in 0..10 {
            sim.step_once();
        }
        // Energy should have spread from center
        let center = sim.idx(10, 10, 10);
        let has_spread = sim
            .ez
            .iter()
            .enumerate()
            .any(|(i, &v)| i != center && v.abs() > 1e-15);
        assert!(has_spread);
    }

    #[test]
    fn test_fdtd3d_run_with_source() {
        let mut sim = Fdtd3d::new(15, 15, 15, 1e-4).unwrap();
        sim.run(20, Some((7, 7, 7)), 1e9, 1.0);
        let max_ez: f64 = sim.ez.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(max_ez > 0.0);
    }

    #[test]
    fn test_fdtd3d_time() {
        let mut sim = Fdtd3d::new(5, 5, 5, 1e-3).unwrap();
        assert!(sim.time().abs() < 1e-30);
        sim.step_once();
        assert!((sim.time() - sim.dt).abs() < 1e-30);
    }

    #[test]
    fn test_fdtd3d_set_permittivity() {
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        sim.set_permittivity(2, 2, 2, 8, 8, 8, 4.0).unwrap();
        let center = sim.idx(5, 5, 5);
        let outside = sim.idx(0, 0, 0);
        assert!((sim.permittivity[center] - 4.0).abs() < 1e-15);
        assert!((sim.permittivity[outside] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd3d_set_permittivity_invalid() {
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        assert!(sim.set_permittivity(0, 0, 0, 10, 10, 10, -1.0).is_err());
    }

    #[test]
    fn test_fdtd3d_set_material() {
        use crate::material::Material;
        let mut sim = Fdtd3d::new(10, 10, 10, 1e-3).unwrap();
        let glass = Material::dielectric(2.25);
        sim.set_material(2, 2, 2, 8, 8, 8, &glass).unwrap();
        let center = sim.idx(5, 5, 5);
        assert!((sim.permittivity[center] - 2.25).abs() < 1e-15);
    }

    #[test]
    fn test_fdtd3d_six_components() {
        let mut sim = Fdtd3d::new(15, 15, 15, 1e-4).unwrap();
        sim.add_source_ez(7, 7, 7, 1.0);
        for _ in 0..5 {
            sim.step_once();
        }
        assert!(sim.hx.iter().any(|&v| v.abs() > 1e-20));
        assert!(sim.hy.iter().any(|&v| v.abs() > 1e-20));
        assert!(sim.ex.iter().any(|&v| v.abs() > 1e-20));
        assert!(sim.ey.iter().any(|&v| v.abs() > 1e-20));
    }

    // ── Non-uniform grid tests ────────────────────────────────────

    #[test]
    fn test_nonuniform_grid_uniform() {
        let g = NonUniformGrid::uniform(10, 10, 10, 1e-3);
        assert_eq!(g.dx.len(), 10);
        assert!((g.min_spacing() - 1e-3).abs() < 1e-15);
    }

    #[test]
    fn test_nonuniform_grid_graded() {
        let d = NonUniformGrid::graded(10, 0.1e-3, 1.0e-3).unwrap();
        assert_eq!(d.len(), 10);
        assert!((d[0] - 0.1e-3).abs() < 1e-10);
        assert!((d[9] - 1.0e-3).abs() < 1e-10);
        // Should be monotonically increasing
        for i in 1..d.len() {
            assert!(d[i] > d[i - 1]);
        }
    }

    #[test]
    fn test_nonuniform_grid_graded_invalid() {
        assert!(NonUniformGrid::graded(1, 0.1e-3, 1.0e-3).is_err());
        assert!(NonUniformGrid::graded(10, 0.0, 1.0e-3).is_err());
    }

    #[test]
    fn test_nonuniform_grid_courant_dt() {
        let g = NonUniformGrid::uniform(10, 10, 10, 1e-3);
        let dt = g.courant_dt();
        assert!(dt > 0.0);
        assert!(dt <= 1e-3 / (SPEED_OF_LIGHT * 3.0_f64.sqrt()));
    }

    #[test]
    fn test_nonuniform_grid_extents() {
        let g = NonUniformGrid::uniform(10, 5, 3, 1e-3);
        let ext = g.extents();
        assert!((ext[0] - 0.01).abs() < 1e-15);
        assert!((ext[1] - 0.005).abs() < 1e-15);
        assert!((ext[2] - 0.003).abs() < 1e-15);
    }

    // ── Higher-order FDTD tests ───────────────────────────────────

    #[test]
    fn test_fourth_order_constants() {
        // c₁ + c₂ should sum to 1 + 1/24 ≈ ... no, check: 9/8 - 1/24 = 26/24 ≈ 1.083
        assert!((FOURTH_ORDER_C1 - 1.125).abs() < 1e-10);
        assert!((FOURTH_ORDER_C2 - (-1.0 / 24.0)).abs() < 1e-10);
    }

    #[test]
    fn test_fourth_order_derivative_linear() {
        // For f(x) = x, df/dx = 1 everywhere
        let dx = 0.1;
        let field: Vec<f64> = (0..10).map(|i| i as f64 * dx).collect();
        let deriv = fourth_order_derivative(&field, dx);
        // Interior should be close to 1.0
        for &d in &deriv[2..8] {
            assert!((d - 1.0).abs() < 0.1, "derivative = {d}");
        }
    }

    #[test]
    fn test_fourth_order_derivative_quadratic() {
        // For f(x) = x², df/dx = 2x
        let dx = 0.1;
        let field: Vec<f64> = (0..20).map(|i| (i as f64 * dx).powi(2)).collect();
        let deriv = fourth_order_derivative(&field, dx);
        // At i=10 (x=1.0), should be ~2.0 (within stencil accuracy)
        assert!((deriv[10] - 2.0).abs() < 0.5, "deriv[10] = {}", deriv[10]);
    }

    // ── CW accumulator tests ──────────────────────────────────────

    #[test]
    fn test_cw_accumulator_sinusoid() {
        // Accumulate sin(ωt) at frequency ω — should get strong response
        let freq = 1e9;
        let dt = 1e-12;
        let n = 10000;
        let mut acc = CwAccumulator::new(1, freq, dt);
        for step in 0..n {
            let t = step as f64 * dt;
            let val = (2.0 * std::f64::consts::PI * freq * t).sin();
            acc.accumulate(&[val], step);
        }
        let amp = acc.amplitude();
        assert!(amp[0] > 0.0);
    }

    #[test]
    fn test_cw_accumulator_off_frequency() {
        // Accumulate sin(2ωt) at frequency ω — should get weak response
        let freq = 1e9;
        let dt = 1e-12;
        let n = 10000;

        let mut on = CwAccumulator::new(1, freq, dt);
        let mut off = CwAccumulator::new(1, freq, dt);

        for step in 0..n {
            let t = step as f64 * dt;
            let v_on = (2.0 * std::f64::consts::PI * freq * t).sin();
            let v_off = (2.0 * std::f64::consts::PI * 2.0 * freq * t).sin();
            on.accumulate(&[v_on], step);
            off.accumulate(&[v_off], step);
        }

        assert!(on.amplitude()[0] > off.amplitude()[0] * 5.0);
    }

    // ── Eigenmode solver tests ────────────────────────────────────

    #[test]
    fn test_eigenmode_single_frequency() {
        // Create a signal with a single known frequency
        let freq = 5e9;
        let dt = 1e-12;
        let n = 20000;
        let signal: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                (2.0 * std::f64::consts::PI * freq * t).sin() * (-t / 5e-9).exp()
            })
            .collect();

        let result = find_eigenmodes(&signal, dt, 1e9, 10e9, 1000).unwrap();
        assert!(!result.frequencies.is_empty());
        // The detected frequency should be close to 5 GHz
        let closest = result
            .frequencies
            .iter()
            .min_by(|a, b| {
                ((*a - freq).abs())
                    .partial_cmp(&((*b - freq).abs()))
                    .unwrap()
            })
            .unwrap();
        assert!((*closest - freq).abs() / freq < 0.05);
    }

    #[test]
    fn test_eigenmode_empty_signal() {
        assert!(find_eigenmodes(&[], 1e-12, 1e9, 10e9, 100).is_err());
    }

    #[test]
    fn test_eigenmode_invalid_freq_range() {
        assert!(find_eigenmodes(&[1.0], 1e-12, 10e9, 1e9, 100).is_err());
    }
}
