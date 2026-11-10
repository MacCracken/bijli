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
/// Models E_z and H_y propagating in the x-direction. `e_field[i]` lives at
/// x = i·dx and `h_field[i]` at x = (i + ½)·dx; `h_field[num_cells − 1]` lies
/// outside the E domain and is never updated.
#[derive(Debug, Clone)]
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
    /// Mur ABC coefficient at the left end: (v⋅dt − dx)/(v⋅dt + dx), with v the
    /// local phase velocity next to the boundary.
    abc_left: f64,
    /// Mur ABC coefficient at the right end (same form as `abc_left`).
    abc_right: f64,
}

/// First-order Mur coefficient (v·dt − dx)/(v·dt + dx) for a medium with
/// relative permittivity `eps_r` and permeability `mu_r`.
#[inline]
#[must_use]
fn mur_coefficient(dt: f64, dx: f64, eps_r: f64, mu_r: f64) -> f64 {
    let v_dt = SPEED_OF_LIGHT * dt / (eps_r * mu_r).sqrt();
    (v_dt - dx) / (v_dt + dx)
}

/// `true` when `x` is a finite, strictly positive number (NaN-safe).
#[inline]
#[must_use]
fn is_positive_finite(x: f64) -> bool {
    x > 0.0 && x.is_finite()
}

impl Fdtd1d {
    /// Create a new 1D FDTD simulation.
    ///
    /// `num_cells` is the number of grid cells. `dx` is the spatial resolution.
    /// The time step is set to **half** the 1D Courant limit, dt = dx/(2c)
    /// (Courant number S = c·dt/dx = 0.5), for a safety margin.
    ///
    /// Both ends use first-order Mur absorbing boundaries.
    ///
    /// # Errors
    /// - [`BijliError::InsufficientResolution`] if `num_cells < 2`.
    /// - [`BijliError::InvalidParameter`] if `dx` is not finite and positive.
    pub fn new(num_cells: usize, dx: f64) -> Result<Self> {
        if num_cells < 2 {
            return Err(BijliError::InsufficientResolution {
                cells: num_cells,
                domain_size: num_cells as f64 * dx,
            });
        }
        if !is_positive_finite(dx) {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be finite and positive, got {dx}"),
            });
        }

        // Courant stability condition: dt ≤ dx/c
        // Use half the Courant limit for safety margin
        let dt = dx / (2.0 * SPEED_OF_LIGHT);

        let e_coeff = vec![dt / (EPSILON_0 * dx); num_cells];
        let h_coeff = vec![dt / (MU_0 * dx); num_cells];

        // Vacuum Mur coefficient; recomputed whenever materials change.
        let abc_coeff = mur_coefficient(dt, dx, 1.0, 1.0);

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
            abc_left: abc_coeff,
            abc_right: abc_coeff,
        })
    }

    /// Recompute the per-end Mur coefficients from the material next to each
    /// boundary.
    ///
    /// The left end uses ε_r of the first interior E node (index 1) and μ_r of
    /// the first H node (index 0); the right end uses index `n − 2` for both —
    /// i.e. the medium the boundary-adjacent updates actually propagate in.
    fn update_abc_coefficients(&mut self) {
        let n = self
            .num_cells
            .min(self.permittivity.len())
            .min(self.permeability.len());
        if n < 2 {
            return;
        }
        self.abc_left =
            mur_coefficient(self.dt, self.dx, self.permittivity[1], self.permeability[0]);
        self.abc_right = mur_coefficient(
            self.dt,
            self.dx,
            self.permittivity[n - 2],
            self.permeability[n - 2],
        );
    }

    /// Set relative permittivity for cells `start..end` (`end` is clamped to the grid).
    ///
    /// The Mur boundary coefficients are recomputed so the ABCs stay matched to
    /// the medium next to each boundary.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] if `eps_r` is not finite and positive.
    pub fn set_permittivity(&mut self, start: usize, end: usize, eps_r: f64) -> Result<()> {
        if !is_positive_finite(eps_r) {
            return Err(BijliError::InvalidPermittivity { value: eps_r });
        }
        let end = end
            .min(self.num_cells)
            .min(self.permittivity.len())
            .min(self.e_coeff.len());
        let base_coeff = self.dt / (EPSILON_0 * self.dx);
        for i in start..end {
            self.permittivity[i] = eps_r;
            self.e_coeff[i] = base_coeff / eps_r;
        }
        self.update_abc_coefficients();
        Ok(())
    }

    /// Set relative permeability for cells `start..end` (`end` is clamped to the grid).
    ///
    /// The Mur boundary coefficients are recomputed so the ABCs stay matched to
    /// the medium next to each boundary.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermeability`] if `mu_r` is not finite and positive.
    pub fn set_permeability(&mut self, start: usize, end: usize, mu_r: f64) -> Result<()> {
        if !is_positive_finite(mu_r) {
            return Err(BijliError::InvalidPermeability { value: mu_r });
        }
        let end = end
            .min(self.num_cells)
            .min(self.permeability.len())
            .min(self.h_coeff.len());
        let base_coeff = self.dt / (MU_0 * self.dx);
        for i in start..end {
            self.permeability[i] = mu_r;
            self.h_coeff[i] = base_coeff / mu_r;
        }
        self.update_abc_coefficients();
        Ok(())
    }

    /// Inject a soft source (additive) at a grid point. Out-of-grid cells are ignored.
    #[inline]
    pub fn add_source(&mut self, cell: usize, value: f64) {
        if cell < self.num_cells
            && let Some(e) = self.e_field.get_mut(cell)
        {
            *e += value;
        }
    }

    /// Number of cells the update loops can safely touch (guards against the
    /// public `num_cells`/field vectors having been changed inconsistently).
    #[inline]
    fn active_cells(&self) -> usize {
        self.num_cells
            .min(self.e_field.len())
            .min(self.h_field.len())
            .min(self.e_coeff.len())
            .min(self.h_coeff.len())
    }

    /// Advance the simulation by one time step using the Yee algorithm.
    ///
    /// Uses first-order Mur absorbing boundary conditions at both ends
    /// (Mur 1981; Taflove & Hagness, ch. 7):
    /// E₀ⁿ⁺¹ = E₁ⁿ + k (E₁ⁿ⁺¹ − E₀ⁿ), k = (v·dt − dx)/(v·dt + dx).
    ///
    /// If the public fields have been resized inconsistently the update is
    /// restricted to the common length (and is a no-op below 2 cells).
    pub fn step_once(&mut self) {
        let n = self.active_cells();
        if n < 2 {
            self.step += 1;
            return;
        }
        let e = &mut self.e_field[..n];
        let h = &mut self.h_field[..n];
        let ce = &self.e_coeff[..n];
        let ch = &self.h_coeff[..n];

        // Store boundary-neighbour values (time n) for the Mur ABC
        let e_left = e[1];
        let e_right = e[n - 2];

        // Update H from E: H^{n+1/2} = H^{n-1/2} + coeff_h[i] × (E^n[i+1] - E^n[i])
        // (zipped iterators: no bounds checks, vectorizable)
        for ((hi, &chi), w) in h[..n - 1].iter_mut().zip(&ch[..n - 1]).zip(e.windows(2)) {
            *hi += chi * (w[1] - w[0]);
        }

        // Update E from H: E^{n+1} = E^n + coeff_e[i] × (H^{n+1/2}[i] - H^{n+1/2}[i-1])
        for ((ei, &cei), w) in e[1..n - 1]
            .iter_mut()
            .zip(&ce[1..n - 1])
            .zip(h[..n - 1].windows(2))
        {
            *ei += cei * (w[1] - w[0]);
        }

        // First-order Mur absorbing boundary conditions
        e[0] = e_left + self.abc_left * (e[1] - e[0]);
        e[n - 1] = e_right + self.abc_right * (e[n - 2] - e[n - 1]);

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
        let energy: f64 = self
            .e_field
            .iter()
            .zip(&self.h_field)
            .zip(self.permittivity.iter().zip(&self.permeability))
            .take(self.num_cells)
            .map(|((&e, &h), (&eps_r, &mu_r))| {
                0.5 * eps_r * EPSILON_0 * e * e + 0.5 * mu_r * MU_0 * h * h
            })
            .sum();
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
/// - `tau`: pulse width parameter (s) — relates to bandwidth as BW ≈ 1/(πτ).
///   Must be non-zero; `tau == 0` or NaN returns `0.0` instead of NaN.
#[inline]
#[must_use]
pub fn gaussian_pulse(t: f64, t0: f64, tau: f64) -> f64 {
    if tau == 0.0 || tau.is_nan() {
        return 0.0;
    }
    let arg = (t - t0) / tau;
    (-arg * arg).exp()
}

/// Ricker wavelet ("Mexican hat"): the normalized negative second derivative
/// of the Gaussian [`gaussian_pulse`].
///
/// r(t) = (1 − 2u²) exp(−u²), u = (t − t₀)/τ
///
/// Equivalently r(t) = (1 − 2π²f_p²(t−t₀)²) exp(−π²f_p²(t−t₀)²) with peak
/// spectral frequency f_p = 1/(πτ). The waveform is even about t₀, has peak
/// value 1 at t = t₀ and zero DC content.
///
/// `tau == 0` or NaN returns `0.0`.
///
/// **Changed in 1.2.0:** earlier versions returned the *first* derivative of
/// the Gaussian (odd, peak ≈ 0.858/τ). That waveform is now
/// [`differentiated_gaussian`].
#[inline]
#[must_use]
pub fn ricker_wavelet(t: f64, t0: f64, tau: f64) -> f64 {
    if tau == 0.0 || tau.is_nan() {
        return 0.0;
    }
    let arg = (t - t0) / tau;
    let arg2 = arg * arg;
    (1.0 - 2.0 * arg2) * (-arg2).exp()
}

/// Differentiated Gaussian pulse: the time derivative of [`gaussian_pulse`].
///
/// g'(t) = −2(t−t₀)/τ² × exp(−((t−t₀)/τ)²)
///
/// Odd about t₀, zero DC content, peak magnitude √2·e^{−½}/τ ≈ 0.858/τ.
/// This is the waveform `ricker_wavelet` returned before 1.2.0.
///
/// `tau == 0` or NaN returns `0.0`.
#[inline]
#[must_use]
pub fn differentiated_gaussian(t: f64, t0: f64, tau: f64) -> f64 {
    if tau == 0.0 || tau.is_nan() {
        return 0.0;
    }
    let arg = (t - t0) / tau;
    -2.0 * arg / tau * (-arg * arg).exp()
}

/// Modulated Gaussian pulse (Gaussian envelope × sinusoidal carrier).
///
/// g(t) = exp(−((t−t₀)/τ)²) × sin(2πf₀(t−t₀))
///
/// Produces a narrowband pulse centered at frequency f₀.
/// `tau == 0` or NaN returns `0.0`.
#[inline]
#[must_use]
pub fn modulated_gaussian_pulse(t: f64, t0: f64, tau: f64, frequency: f64) -> f64 {
    if tau == 0.0 || tau.is_nan() {
        return 0.0;
    }
    let arg = (t - t0) / tau;
    (-arg * arg).exp() * (2.0 * std::f64::consts::PI * frequency * (t - t0)).sin()
}

/// Total-field/scattered-field (TFSF) plane wave injection parameters.
///
/// Holds a 1D auxiliary grid that computes the incident plane wave (E_z, H_y
/// travelling in +x), whose values are used to correct the main-grid updates
/// next to a TFSF boundary (Schneider, *Understanding the FDTD Method*,
/// §3.10–3.11; Taflove & Hagness, ch. 5).
///
/// # Grid layout and timing
/// - `inc_e[i]` is E^inc at x = i·dx; `inc_h[i]` is H^inc at x = (i + ½)·dx
///   (`inc_h[nx − 1]` lies outside the grid and stays 0).
/// - After k calls to [`step`](Self::step), `inc_e` holds E^inc at time k·dt and
///   `inc_h` holds H^inc at time (k − ½)·dt.
/// - Both ends of the auxiliary grid use first-order Mur absorbing boundaries,
///   so the incident pulse leaves the grid instead of bouncing back.
///
/// # Left-edge TFSF recipe
/// For a boundary at column `bx` (cells `x ≥ bx` are total-field) and update
/// conventions `H[i] += c_h·(E[i+1] − E[i])`, `E[i] += c_e·(H[i] − H[i−1])` (as in
/// [`Fdtd1d`] and the H_y/E_z updates of [`Fdtd2d`] TM mode), each main step is:
/// 1. update main H; then `H[bx − 1] −= c_h · incident_e()` (E^inc at time n),
/// 2. call `step(source)` (advances H^inc to n + ½ and E^inc to n + 1),
/// 3. update main E; then `E[bx] −= c_e · incident_h()` (H^inc at bx − ½, n + ½).
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
    /// Mur ABC coefficient (S − 1)/(S + 1), S = c·dt/dx.
    abc_coeff: f64,
}

impl TfsfSource {
    /// Create a TFSF source for a 2D grid, sanitizing invalid parameters.
    ///
    /// `boundary_x` is the first total-field column. The 1D auxiliary grid has
    /// `nx` cells matching the 2D grid's x-dimension, and the soft source sits at
    /// `max(boundary_x / 4, 1)`.
    ///
    /// This constructor never panics: out-of-range inputs are clamped (with a
    /// `tracing` warning) — `nx` to ≥ 3, `boundary_x` to `1..nx`, an invalid
    /// `dx`/`dt` is replaced so the Courant number S = c·dt/dx is 1, and
    /// `dt > dx/c` is clamped to `dx/c`. Prefer [`try_new`](Self::try_new),
    /// which reports invalid input as an error instead.
    #[must_use]
    pub fn new(nx: usize, boundary_x: usize, dx: f64, dt: f64) -> Self {
        let nx_s = nx.max(3);
        if nx_s != nx {
            tracing::warn!(nx, "TfsfSource: nx < 3, using 3 cells");
        }
        let bx = boundary_x.clamp(1, nx_s - 1);
        if bx != boundary_x {
            tracing::warn!(
                boundary_x,
                bx,
                "TfsfSource: boundary_x out of 1..nx, clamped"
            );
        }
        let (dx_s, dt_s) = match (is_positive_finite(dx), is_positive_finite(dt)) {
            // Same small tolerance as `try_new`, so dt = dx/c passes unchanged.
            (true, true) if SPEED_OF_LIGHT * dt <= dx * (1.0 + 1e-12) => (dx, dt),
            (true, true) => {
                tracing::warn!(
                    dx,
                    dt,
                    "TfsfSource: dt exceeds 1D Courant limit, clamped to dx/c"
                );
                (dx, dx / SPEED_OF_LIGHT)
            }
            (true, false) => {
                tracing::warn!(dt, "TfsfSource: invalid dt, using dx/c");
                (dx, dx / SPEED_OF_LIGHT)
            }
            (false, true) => {
                tracing::warn!(dx, "TfsfSource: invalid dx, using c·dt");
                (SPEED_OF_LIGHT * dt, dt)
            }
            (false, false) => {
                tracing::warn!(
                    dx,
                    dt,
                    "TfsfSource: invalid dx and dt, using dx = 1, dt = 1/c"
                );
                (1.0, 1.0 / SPEED_OF_LIGHT)
            }
        };
        Self::build(nx_s, bx, dx_s, dt_s)
    }

    /// Create a TFSF source, validating all parameters.
    ///
    /// # Errors
    /// - [`BijliError::InsufficientResolution`] if `nx < 3`.
    /// - [`BijliError::InvalidParameter`] if `boundary_x` is not in `1..nx`,
    ///   `dx` is not finite and positive, or `dt` is not in `(0, dx/c]`
    ///   (the 1D Courant limit).
    pub fn try_new(nx: usize, boundary_x: usize, dx: f64, dt: f64) -> Result<Self> {
        if nx < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx,
                domain_size: nx as f64 * dx,
            });
        }
        if boundary_x == 0 || boundary_x >= nx {
            return Err(BijliError::InvalidParameter {
                reason: format!("TFSF boundary_x must be in 1..{nx}, got {boundary_x}"),
            });
        }
        if !is_positive_finite(dx) {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be finite and positive, got {dx}"),
            });
        }
        // Small relative tolerance so dt = dx/c computed in floating point passes.
        if !is_positive_finite(dt) || SPEED_OF_LIGHT * dt > dx * (1.0 + 1e-12) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "dt must be in (0, dx/c] = (0, {}], got {dt}",
                    dx / SPEED_OF_LIGHT
                ),
            });
        }
        Ok(Self::build(nx, boundary_x, dx, dt))
    }

    /// Construct from already-validated parameters.
    fn build(nx: usize, boundary_x: usize, dx: f64, dt: f64) -> Self {
        let e_coeff_1d = dt / (EPSILON_0 * dx);
        let h_coeff_1d = dt / (MU_0 * dx);
        let s = SPEED_OF_LIGHT * dt / dx;
        let abc_coeff = (s - 1.0) / (s + 1.0);
        // Never inject at node 0 (overwritten by the ABC); the source must be
        // an interior node of the auxiliary grid.
        let source_cell = (boundary_x / 4).clamp(1, nx - 2);
        tracing::debug!(
            nx,
            boundary_x,
            source_cell,
            dx,
            dt,
            courant = s,
            "TFSF source created"
        );
        Self {
            inc_e: vec![0.0; nx],
            inc_h: vec![0.0; nx],
            boundary_x,
            source_cell,
            e_coeff_1d,
            h_coeff_1d,
            abc_coeff,
        }
    }

    /// Advance the 1D auxiliary grid by one step and inject a soft (additive)
    /// source value at `source_cell`.
    ///
    /// Never panics: if the public fields were modified inconsistently, the
    /// update runs over the common length of `inc_e`/`inc_h` (no-op below 3
    /// cells) and the source is clamped to an interior node.
    pub fn step(&mut self, source_value: f64) {
        let n = self.inc_e.len().min(self.inc_h.len());
        if n < 3 {
            return;
        }
        let e = &mut self.inc_e[..n];
        let h = &mut self.inc_h[..n];

        // Boundary-neighbour values at time n for the Mur ABCs
        let e_left = e[1];
        let e_right = e[n - 2];

        // Update H (H[n-1] at n−½ lies outside the E domain and is not updated)
        for i in 0..n - 1 {
            h[i] += self.h_coeff_1d * (e[i + 1] - e[i]);
        }

        // Inject source at an interior node
        e[self.source_cell.clamp(1, n - 2)] += source_value;

        // Update interior E
        for i in 1..n - 1 {
            e[i] += self.e_coeff_1d * (h[i] - h[i - 1]);
        }

        // First-order Mur ABCs at both ends
        e[0] = e_left + self.abc_coeff * (e[1] - e[0]);
        e[n - 1] = e_right + self.abc_coeff * (e[n - 2] - e[n - 1]);
    }

    /// Incident E-field at the TFSF boundary node x = `boundary_x`·dx.
    ///
    /// Returns 0.0 if `boundary_x` is outside the auxiliary grid.
    #[inline]
    #[must_use]
    pub fn incident_e(&self) -> f64 {
        self.inc_e.get(self.boundary_x).copied().unwrap_or(0.0)
    }

    /// Incident H-field at x = (`boundary_x` − ½)·dx, i.e. `inc_h[boundary_x − 1]`.
    ///
    /// This is the scattered-field H node adjacent to the boundary, whose
    /// incident value corrects the update of E at `boundary_x` (Schneider §3.10).
    /// Returns 0.0 if that node does not exist (`boundary_x == 0` or out of range).
    ///
    /// **Changed in 1.2.0:** earlier versions returned `inc_h[boundary_x]`
    /// (x = (`boundary_x` + ½)·dx, inside the total-field region), which is not
    /// the node required by the left-edge TFSF correction.
    #[inline]
    #[must_use]
    pub fn incident_h(&self) -> f64 {
        self.boundary_x
            .checked_sub(1)
            .and_then(|i| self.inc_h.get(i))
            .copied()
            .unwrap_or(0.0)
    }
}

// ── Dispersive material models ────────────────────────────────────

/// Validate that `value` is finite and non-negative (NaN-safe).
#[inline]
fn check_non_negative(name: &str, value: f64) -> Result<()> {
    if value >= 0.0 && value.is_finite() {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!("{name} must be finite and non-negative, got {value}"),
        })
    }
}

/// Validate that `value` is finite and strictly positive (NaN-safe).
#[inline]
fn check_positive(name: &str, value: f64) -> Result<()> {
    if is_positive_finite(value) {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!("{name} must be finite and positive, got {value}"),
        })
    }
}

/// Validate a relative permittivity (finite, > 0).
#[inline]
fn check_permittivity(value: f64) -> Result<()> {
    if is_positive_finite(value) {
        Ok(())
    } else {
        Err(BijliError::InvalidPermittivity { value })
    }
}

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
    ///
    /// Does not validate its inputs (ε_r ≤ 0 makes the coefficients
    /// meaningless); use [`try_new`](Self::try_new) for checked construction.
    #[inline]
    #[must_use]
    pub fn new(eps_r: f64, conductivity: f64) -> Self {
        Self {
            eps_r,
            conductivity,
        }
    }

    /// Create a lossy material, validating its parameters.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] if `eps_r` is not finite and positive;
    /// [`BijliError::InvalidParameter`] if `conductivity` is negative or not finite.
    pub fn try_new(eps_r: f64, conductivity: f64) -> Result<Self> {
        check_permittivity(eps_r)?;
        check_non_negative("conductivity", conductivity)?;
        Ok(Self::new(eps_r, conductivity))
    }

    /// E-field update coefficients (C_a, C_b) for FDTD.
    ///
    /// Requires `dt > 0`, `dx > 0` and a valid material (see [`try_new`](Self::try_new)).
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
/// ε(ω) = ε_∞ − ω_p²/(ω² + iγω)   (e^{−iωt} convention)
///
/// Uses Auxiliary Differential Equation (ADE) method for FDTD integration:
/// the polarization current obeys dJ/dt + γJ = ε₀ω_p² E.
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
    ///
    /// Does not validate its inputs (e.g. `eps_inf = 0` makes κ infinite); use
    /// [`try_new`](Self::try_new) for checked construction.
    #[inline]
    #[must_use]
    pub fn new(eps_inf: f64, omega_p: f64, gamma: f64) -> Self {
        Self {
            eps_inf,
            omega_p,
            gamma,
        }
    }

    /// Create a Drude material, validating its parameters.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] if `eps_inf` is not finite and positive;
    /// [`BijliError::InvalidParameter`] if `omega_p` or `gamma` is negative or not finite.
    pub fn try_new(eps_inf: f64, omega_p: f64, gamma: f64) -> Result<Self> {
        check_permittivity(eps_inf)?;
        check_non_negative("plasma frequency omega_p", omega_p)?;
        check_non_negative("damping gamma", gamma)?;
        Ok(Self::new(eps_inf, omega_p, gamma))
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
    /// Returns (α, β, κ) with α = (2 − γΔt)/(2 + γΔt), β = ε₀ω_p²Δt/(2 + γΔt),
    /// κ = Δt/(ε₀ε_∞), for the centred (trapezoidal) discretization of
    /// dJ/dt + γJ = ε₀ω_p²E (Taflove & Hagness, ch. 9):
    ///
    /// - J^{n+1} = α J^n + β (E^{n+1} + E^n)
    /// - E^{n+1} = E^n + κ [(∇×H)^{n+½} − ½(J^{n+1} + J^n)]
    ///
    /// Substituting the first into the second gives the explicit E update
    /// E^{n+1} (1 + κβ/2) = E^n (1 − κβ/2) + κ [(∇×H)^{n+½} − ½(1 + α) J^n].
    ///
    /// Requires `dt > 0` and a valid material (see [`try_new`](Self::try_new)).
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
/// ε(ω) = ε_∞ + Δε ω₀²/(ω₀² − ω² − iγω)   (e^{−iωt} convention)
///
/// The polarization obeys d²P/dt² + γ dP/dt + ω₀² P = ε₀ Δε ω₀² E.
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
    ///
    /// Does not validate its inputs; use [`try_new`](Self::try_new) for
    /// checked construction.
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

    /// Create a Lorentz material, validating its parameters.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] if `eps_inf` is not finite and positive;
    /// [`BijliError::InvalidParameter`] if `delta_eps` or `gamma` is negative or
    /// not finite, or `omega_0` is not finite and positive.
    pub fn try_new(eps_inf: f64, delta_eps: f64, omega_0: f64, gamma: f64) -> Result<Self> {
        check_permittivity(eps_inf)?;
        check_non_negative("oscillator strength delta_eps", delta_eps)?;
        check_positive("resonance frequency omega_0", omega_0)?;
        check_non_negative("damping gamma", gamma)?;
        Ok(Self::new(eps_inf, delta_eps, omega_0, gamma))
    }

    /// ADE update coefficients for the Lorentz polarization P.
    ///
    /// Returns (α, β, γ_coeff) for the central-difference discretization of the
    /// oscillator equation:
    /// P^{n+1} = α P^n + β P^{n-1} + γ_coeff E^n, with
    /// α = (2 − ω₀²Δt²)/(1 + γΔt/2), β = −(1 − γΔt/2)/(1 + γΔt/2),
    /// γ_coeff = ε₀Δεω₀²Δt²/(1 + γΔt/2).
    ///
    /// Requires `dt > 0`; the recursion is stable only for ω₀Δt < 2.
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
/// ε(ω) = ε_∞ + (ε_s − ε_∞)/(1 − iωτ)   (e^{−iωt} convention as for Drude/Lorentz;
/// equivalently 1 + jωτ with the engineering e^{jωt} convention)
///
/// The polarization obeys τ dP/dt + P = ε₀(ε_s − ε_∞) E.
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
    ///
    /// Does not validate its inputs; use [`try_new`](Self::try_new) for
    /// checked construction.
    #[inline]
    #[must_use]
    pub fn new(eps_inf: f64, eps_s: f64, tau: f64) -> Self {
        Self {
            eps_inf,
            eps_s,
            tau,
        }
    }

    /// Create a Debye material, validating its parameters.
    ///
    /// # Errors
    /// [`BijliError::InvalidPermittivity`] if `eps_inf` or `eps_s` is not finite
    /// and positive; [`BijliError::InvalidParameter`] if `eps_s < eps_inf` or
    /// `tau` is not finite and positive.
    pub fn try_new(eps_inf: f64, eps_s: f64, tau: f64) -> Result<Self> {
        check_permittivity(eps_inf)?;
        check_permittivity(eps_s)?;
        if eps_s < eps_inf {
            return Err(BijliError::InvalidParameter {
                reason: format!("Debye eps_s ({eps_s}) must be >= eps_inf ({eps_inf})"),
            });
        }
        check_positive("relaxation time tau", tau)?;
        Ok(Self::new(eps_inf, eps_s, tau))
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
    /// Returns (α, β) for the trapezoidal discretization of
    /// τ dP/dt + P = ε₀(ε_s − ε_∞)E:
    /// P^{n+1} = α P^n + β (E^{n+1} + E^n), with
    /// α = (2τ/Δt − 1)/(2τ/Δt + 1), β = ε₀(ε_s − ε_∞)/(2τ/Δt + 1).
    ///
    /// Requires `dt > 0` and a valid material (see [`try_new`](Self::try_new)).
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
    /// Far-field pattern magnitude at each observation angle (always ≥ 0):
    /// lim_{ρ→∞} √ρ·|E_z(ρ, φ)|. For phasor inputs in V/m this is in V·m^{−½};
    /// it is linear in the input field amplitudes (see [`nf2ff_2d_phasor`]).
    pub pattern: Vec<f64>,
    /// Observation angles φ (radians, 0 to 2π, measured from +x towards +y).
    pub angles: Vec<f64>,
}

/// Minimal complex number for the NF2FF accumulation (the `fdtd` feature does
/// not depend on the `polarization` module's complex type).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct Cplx {
    re: f64,
    im: f64,
}

impl Cplx {
    #[inline]
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// e^{iθ}
    #[inline]
    fn cis(theta: f64) -> Self {
        let (s, c) = theta.sin_cos();
        Self { re: c, im: s }
    }

    #[inline]
    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }

    #[inline]
    fn add(self, o: Self) -> Self {
        Self {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }

    #[inline]
    fn scale(self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }

    #[inline]
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// One straight side of the NF2FF integration contour.
struct Nf2ffSide {
    /// Position of the first sample relative to the contour centre (m).
    start: (f64, f64),
    /// Unit step direction along the side.
    dir: (f64, f64),
    /// Quadrature-weighted equivalent currents (J_z, M_x, M_y) per sample.
    samples: Vec<(Cplx, Cplx, Cplx)>,
}

/// Validate the shared NF2FF parameters; returns the free-space wavenumber.
fn nf2ff_validate(
    nx: usize,
    ny: usize,
    dx: f64,
    frequency: f64,
    box_margin: usize,
    n_angles: usize,
) -> Result<f64> {
    if n_angles == 0 {
        return Err(BijliError::InvalidParameter {
            reason: "need at least one observation angle".into(),
        });
    }
    if box_margin == 0 || box_margin >= nx / 2 || box_margin >= ny / 2 {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "NF2FF box margin must be in 1..min(nx, ny)/2 for a {nx}x{ny} grid, got {box_margin}"
            ),
        });
    }
    if !is_positive_finite(dx) {
        return Err(BijliError::InvalidParameter {
            reason: format!("dx must be finite and positive, got {dx}"),
        });
    }
    if !is_positive_finite(frequency) {
        return Err(BijliError::InvalidParameter {
            reason: format!("NF2FF frequency must be finite and positive, got {frequency}"),
        });
    }
    Ok(2.0 * std::f64::consts::PI * frequency / SPEED_OF_LIGHT)
}

/// Check that a field slice holds at least `nx·ny` samples.
fn nf2ff_check_len(name: &str, len: usize, nx: usize, ny: usize) -> Result<()> {
    match nx.checked_mul(ny) {
        Some(n) if len >= n => Ok(()),
        _ => Err(BijliError::InvalidParameter {
            reason: format!("NF2FF field `{name}` has {len} samples, need {nx}x{ny}"),
        }),
    }
}

/// TM-mode NF2FF core (2D Love equivalence, e^{jωt} phasors).
///
/// The contour runs through the E_z nodes x = x₀, x₁ and y = y₀, y₁ with
/// x₀ = y₀ = `box_margin`, x₁ = nx − 1 − margin, y₁ = ny − 1 − margin.
/// Tangential H on the contour is the average of the two Yee H samples half a
/// cell on either side. Equivalent currents: J = n̂×H (J_z = n_x H_y − n_y H_x),
/// M = −n̂×E (M_x = −n_y E_z, M_y = n_x E_z). The far field is
/// E_z ≈ −jk G₀(ρ) [η N_z − L_φ], |G₀| = ¼√(2/(πkρ)), with
/// N_z = ∮J_z e^{jk ρ̂·ρ'} dl', L_φ = −L_x sin φ + L_y cos φ; the returned
/// pattern is √(k/(8π))·|η N_z − L_φ| = √ρ |E_z|.
///
/// Callers must have validated the parameters and slice lengths.
#[allow(clippy::too_many_arguments)]
fn nf2ff_tm_core(
    ez: impl Fn(usize) -> Cplx,
    hx: impl Fn(usize) -> Cplx,
    hy: impl Fn(usize) -> Cplx,
    h_rotation: Cplx,
    nx: usize,
    ny: usize,
    dx: f64,
    k: f64,
    box_margin: usize,
    n_angles: usize,
) -> Nf2ffResult {
    let (x0, x1) = (box_margin, nx - 1 - box_margin);
    let (y0, y1) = (box_margin, ny - 1 - box_margin);
    let xc = 0.5 * (x0 + x1) as f64;
    let yc = 0.5 * (y0 + y1) as f64;
    let idx = |x: usize, y: usize| y * nx + x;
    let half = 0.5;
    // Trapezoid weights along each side (corner samples get half weight per side).
    let weight = |i: usize, last: usize| if i == 0 || i == last { half * dx } else { dx };

    let mut sides = Vec::with_capacity(4);
    // Bottom, y = y₀, n̂ = −ŷ: J_z = H_x, M_x = E_z, M_y = 0.
    // Top,    y = y₁, n̂ = +ŷ: J_z = −H_x, M_x = −E_z, M_y = 0.
    for (y, sign, y_below) in [(y0, 1.0, y0 - 1), (y1, -1.0, y1 - 1)] {
        let last = x1 - x0;
        let samples = (x0..=x1)
            .enumerate()
            .map(|(i, x)| {
                let w = weight(i, last);
                let h_t = hx(idx(x, y_below)).add(hx(idx(x, y))).scale(half);
                let e = ez(idx(x, y));
                (
                    h_t.mul(h_rotation).scale(sign * w),
                    e.scale(sign * w),
                    Cplx::default(),
                )
            })
            .collect();
        sides.push(Nf2ffSide {
            start: ((x0 as f64 - xc) * dx, (y as f64 - yc) * dx),
            dir: (1.0, 0.0),
            samples,
        });
    }
    // Left,  x = x₀, n̂ = −x̂: J_z = −H_y, M_x = 0, M_y = −E_z.
    // Right, x = x₁, n̂ = +x̂: J_z = H_y,  M_x = 0, M_y = E_z.
    for (x, sign) in [(x0, -1.0), (x1, 1.0)] {
        let last = y1 - y0;
        let samples = (y0..=y1)
            .enumerate()
            .map(|(i, y)| {
                let w = weight(i, last);
                let h_t = hy(idx(x - 1, y)).add(hy(idx(x, y))).scale(half);
                let e = ez(idx(x, y));
                (
                    h_t.mul(h_rotation).scale(sign * w),
                    Cplx::default(),
                    e.scale(sign * w),
                )
            })
            .collect();
        sides.push(Nf2ffSide {
            start: ((x as f64 - xc) * dx, (y0 as f64 - yc) * dx),
            dir: (0.0, 1.0),
            samples,
        });
    }

    let eta0 = (MU_0 / EPSILON_0).sqrt();
    let norm = (k / (8.0 * std::f64::consts::PI)).sqrt();
    let angles: Vec<f64> = (0..n_angles)
        .map(|i| 2.0 * std::f64::consts::PI * i as f64 / n_angles as f64)
        .collect();
    let pattern = angles
        .iter()
        .map(|&phi| {
            let (sin_phi, cos_phi) = phi.sin_cos();
            let mut n_z = Cplx::default();
            let mut l_x = Cplx::default();
            let mut l_y = Cplx::default();
            for side in &sides {
                // e^{jk ρ̂·ρ'} advanced by a constant phasor along the side.
                let mut rot = Cplx::cis(k * (side.start.0 * cos_phi + side.start.1 * sin_phi));
                let step = Cplx::cis(k * dx * (side.dir.0 * cos_phi + side.dir.1 * sin_phi));
                for &(jz, mx, my) in &side.samples {
                    n_z = n_z.add(jz.mul(rot));
                    l_x = l_x.add(mx.mul(rot));
                    l_y = l_y.add(my.mul(rot));
                    rot = rot.mul(step);
                }
            }
            let l_phi = l_y.scale(cos_phi).add(l_x.scale(-sin_phi));
            norm * n_z.scale(eta0).add(l_phi.scale(-1.0)).abs()
        })
        .collect();

    Nf2ffResult { pattern, angles }
}

/// Compute a 2D NF2FF transform from a single TM-mode FDTD snapshot.
///
/// Integrates the equivalent surface currents J = n̂×H and M = −n̂×E (Love's
/// equivalence principle) over a rectangular contour through the E_z nodes
/// `box_margin` cells inside the grid edge, and returns |F(φ)| normalized as
/// described on [`Nf2ffResult::pattern`].
///
/// **Monochromatic approximation:** a time snapshot is real-valued, so the
/// instantaneous fields are treated as the real part of the phasor (imaginary
/// part zero), and H (stored half a time step earlier than E on the Yee grid)
/// is used as is. The result therefore depends on the snapshot instant and is
/// only a rough estimate of the true far-field pattern. For accurate patterns
/// accumulate E_z, H_x, H_y with [`CwAccumulator`] and use [`nf2ff_2d_phasor`].
///
/// # Parameters
/// - `sim`: the 2D FDTD simulation (must be TM mode)
/// - `box_margin`: cells from grid edge to NF2FF integration contour (≥ 1)
/// - `n_angles`: number of observation angles (0 to 2π)
/// - `frequency`: source frequency (Hz)
///
/// # Errors
/// [`BijliError::InvalidParameter`] if the simulation is not TM mode, the
/// margin is 0 or too large for the grid, `n_angles == 0`, `frequency` is not
/// finite and positive, or the field vectors are shorter than nx·ny.
pub fn nf2ff_2d(
    sim: &Fdtd2d,
    box_margin: usize,
    n_angles: usize,
    frequency: f64,
) -> Result<Nf2ffResult> {
    if sim.mode != Mode2d::Tm {
        return Err(BijliError::InvalidParameter {
            reason: format!("NF2FF requires a TM-mode simulation, got {:?}", sim.mode),
        });
    }
    let (nx, ny, dx) = (sim.nx, sim.ny, sim.dx);
    let k = nf2ff_validate(nx, ny, dx, frequency, box_margin, n_angles)?;
    nf2ff_check_len("ez", sim.fz.len(), nx, ny)?;
    nf2ff_check_len("hx", sim.fx.len(), nx, ny)?;
    nf2ff_check_len("hy", sim.fy.len(), nx, ny)?;
    tracing::debug!(
        nx,
        ny,
        box_margin,
        n_angles,
        frequency,
        "NF2FF (snapshot) transform"
    );
    let real = |v: &[f64], i: usize| Cplx::new(v.get(i).copied().unwrap_or(0.0), 0.0);
    Ok(nf2ff_tm_core(
        |i| real(&sim.fz, i),
        |i| real(&sim.fx, i),
        |i| real(&sim.fy, i),
        Cplx::new(1.0, 0.0),
        nx,
        ny,
        dx,
        k,
        box_margin,
        n_angles,
    ))
}

/// Compute a 2D TM-mode NF2FF transform from frequency-domain phasors.
///
/// `ez`, `hx`, `hy` hold the running DFTs (complex phasors, e^{jωt} convention
/// as produced by [`CwAccumulator`]) of E_z, H_x, H_y on an `nx`×`ny` Yee grid
/// laid out like [`Fdtd2d`] TM mode (index = y·nx + x; E_z at (x, y), H_x at
/// (x, y + ½), H_y at (x + ½, y)). All three accumulators must share the same
/// frequency and time step; the frequency is taken from them.
///
/// The accumulators are assumed to have been fed after each FDTD step with the
/// same `time_step` index, so H (which on the Yee grid lags E by dt/2) is
/// rotated by e^{+jω·dt/2} to align it with E.
///
/// The transform integrates J = n̂×H and M = −n̂×E over a rectangular contour
/// through the E_z nodes `box_margin` cells inside the grid edge and returns
/// √ρ|E_z(ρ, φ)| as ρ → ∞ (see [`Nf2ffResult::pattern`]). The output is linear
/// in the inputs: for a CW steady state a running DFT is ≈ phasor ×
/// `n_samples`·dt/2, so divide by that factor for absolute amplitudes.
///
/// # Errors
/// [`BijliError::InvalidParameter`] if the accumulators disagree on frequency
/// or time step, the frequency is not positive, any accumulator holds fewer
/// than nx·ny samples, `dx` is not finite and positive, the margin is 0 or too
/// large, or `n_angles == 0`.
#[allow(clippy::too_many_arguments)]
pub fn nf2ff_2d_phasor(
    ez: &CwAccumulator,
    hx: &CwAccumulator,
    hy: &CwAccumulator,
    nx: usize,
    ny: usize,
    dx: f64,
    box_margin: usize,
    n_angles: usize,
) -> Result<Nf2ffResult> {
    let same = |a: f64, b: f64| (a - b).abs() <= 1e-12 * a.abs().max(b.abs());
    if !(same(ez.omega, hx.omega)
        && same(ez.omega, hy.omega)
        && same(ez.dt, hx.dt)
        && same(ez.dt, hy.dt))
    {
        return Err(BijliError::InvalidParameter {
            reason: "NF2FF phasor accumulators must share frequency and dt".into(),
        });
    }
    let frequency = ez.omega / (2.0 * std::f64::consts::PI);
    let k = nf2ff_validate(nx, ny, dx, frequency, box_margin, n_angles)?;
    for (name, acc) in [("ez", ez), ("hx", hx), ("hy", hy)] {
        nf2ff_check_len(name, acc.real.len().min(acc.imag.len()), nx, ny)?;
    }
    if !ez.dt.is_finite() {
        return Err(BijliError::InvalidParameter {
            reason: format!("accumulator dt must be finite, got {}", ez.dt),
        });
    }
    tracing::debug!(
        nx,
        ny,
        box_margin,
        n_angles,
        frequency,
        n_samples = ez.n_samples,
        "NF2FF (phasor) transform"
    );
    let phasor = |acc: &CwAccumulator, i: usize| {
        Cplx::new(
            acc.real.get(i).copied().unwrap_or(0.0),
            acc.imag.get(i).copied().unwrap_or(0.0),
        )
    };
    Ok(nf2ff_tm_core(
        |i| phasor(ez, i),
        |i| phasor(hx, i),
        |i| phasor(hy, i),
        Cplx::cis(0.5 * ez.omega * ez.dt),
        nx,
        ny,
        dx,
        k,
        box_margin,
        n_angles,
    ))
}

// ── Subpixel smoothing ────────────────────────────────────────────

/// Validate subpixel-smoothing inputs (NaN-safe).
#[inline]
fn check_subpixel(eps1: f64, eps2: f64, fill_fraction: f64) -> Result<()> {
    check_permittivity(eps1)?;
    check_permittivity(eps2)?;
    if !(0.0..=1.0).contains(&fill_fraction) {
        return Err(BijliError::InvalidParameter {
            reason: format!("fill_fraction must be in [0, 1], got {fill_fraction}"),
        });
    }
    Ok(())
}

/// Effective permittivity for the E-field component **parallel** to a material
/// interface: the arithmetic mean ⟨ε⟩ = f·ε₂ + (1−f)·ε₁.
///
/// Tangential E is continuous across the interface, so the layers act like
/// capacitors in parallel (Farjadpour et al., Opt. Lett. 31, 2972 (2006)).
/// `fill_fraction` is the fraction of the cell occupied by `eps2`.
///
/// **Fixed in 1.2.0:** earlier versions returned the harmonic mean here (the
/// two subpixel functions were swapped).
///
/// # Errors
/// [`BijliError::InvalidPermittivity`] if either permittivity is not finite and
/// positive; [`BijliError::InvalidParameter`] if `fill_fraction` is not in [0, 1].
#[inline]
pub fn subpixel_permittivity_parallel(eps1: f64, eps2: f64, fill_fraction: f64) -> Result<f64> {
    check_subpixel(eps1, eps2, fill_fraction)?;
    let f = fill_fraction;
    Ok(f * eps2 + (1.0 - f) * eps1)
}

/// Effective permittivity for the E-field component **perpendicular** to a
/// material interface: the harmonic mean ⟨ε⁻¹⟩⁻¹ = 1/(f/ε₂ + (1−f)/ε₁).
///
/// Normal D is continuous across the interface, so the layers act like
/// capacitors in series (Farjadpour et al., Opt. Lett. 31, 2972 (2006)).
/// `fill_fraction` is the fraction of the cell occupied by `eps2`.
///
/// Preconditions: both permittivities finite and positive. `fill_fraction` is
/// clamped to [0, 1] (NaN is treated as 0); invalid permittivities yield NaN.
/// Use [`subpixel_permittivity_perpendicular_checked`] to get an error instead.
///
/// **Fixed in 1.2.0:** earlier versions returned the arithmetic mean here.
#[inline]
#[must_use]
pub fn subpixel_permittivity_perpendicular(eps1: f64, eps2: f64, fill_fraction: f64) -> f64 {
    if !(is_positive_finite(eps1) && is_positive_finite(eps2)) {
        return f64::NAN;
    }
    // f64::clamp propagates NaN; map NaN to 0 first.
    let f = if fill_fraction.is_nan() {
        0.0
    } else {
        fill_fraction.clamp(0.0, 1.0)
    };
    1.0 / (f / eps2 + (1.0 - f) / eps1)
}

/// Checked variant of [`subpixel_permittivity_perpendicular`] (harmonic mean).
///
/// # Errors
/// [`BijliError::InvalidPermittivity`] if either permittivity is not finite and
/// positive; [`BijliError::InvalidParameter`] if `fill_fraction` is not in [0, 1].
#[inline]
pub fn subpixel_permittivity_perpendicular_checked(
    eps1: f64,
    eps2: f64,
    fill_fraction: f64,
) -> Result<f64> {
    check_subpixel(eps1, eps2, fill_fraction)?;
    Ok(subpixel_permittivity_perpendicular(
        eps1,
        eps2,
        fill_fraction,
    ))
}

// ── CPML (Convolutional Perfectly Matched Layer) ──────────────────

/// Free-space wave impedance η₀ = μ₀c ≈ 376.73 Ω.
const ETA_0: f64 = MU_0 * SPEED_OF_LIGHT;

/// Minimum number of grid nodes required between opposite CPML layers:
/// a grid needs `2·thickness + CPML_MIN_INTERIOR ≤ min(nx, ny)` in 2D and
/// `≤ min(nx, ny, nz)` in 3D.
const CPML_MIN_INTERIOR: usize = 3;

/// Whether a CPML of `thickness` cells per side fits a grid whose smallest
/// dimension is `min_dim` nodes (`thickness ≥ 1`, `2·thickness + 3 ≤ min_dim`).
#[inline]
fn cpml_fits(thickness: usize, min_dim: usize) -> bool {
    thickness > 0
        && thickness
            .checked_mul(2)
            .and_then(|v| v.checked_add(CPML_MIN_INTERIOR))
            .is_some_and(|need| need <= min_dim)
}

/// Grading parameters for a [`Cpml`] absorbing layer.
///
/// With ρ ∈ (0, 1] the normalised depth into the layer (ρ = 1 at the outer
/// wall), the profiles follow Roden & Gedney (2000) and Taflove & Hagness
/// (3rd ed., §7.9–7.11):
///
/// - σ(ρ) = σ_max ρ^m, with σ_max = `sigma_ratio` · σ_opt and
///   σ_opt = 0.8 (m + 1) / (η₀ Δx)
/// - κ(ρ) = 1 + (κ_max − 1) ρ^m
/// - α(ρ) = α_max (1 − ρ) (complex-frequency-shifted PML, CFS)
///
/// The defaults of [`CpmlParams::new`] (m = 3, σ_max = σ_opt, κ_max = 1,
/// α_max = 0) are the classic CPML, which absorbs propagating waves over a
/// wide range of incidence angles (≈ −88 dB reflection measured for a
/// 10-cell layer around a 2D point source, ≈ −80 dB in 3D). Raising
/// `kappa_max` (5–15) and `alpha_max` improves absorption of
/// evanescent / grazing fields; α_max = 2π ε₀ f_c only affects frequencies
/// near or below f_c, so keep f_c well below the band of interest.
///
/// κ also shortens the wavelength inside the layer by up to κ_max, so it
/// needs well-resolved sources: with κ_max = 5 a 10-cell layer reflects
/// ≈ −60 dB for a Gaussian with c·τ ≈ 4 cells but ≈ −85 dB at c·τ ≳ 7 cells
/// (2D and 3D alike).
///
/// σ_opt assumes a vacuum background; scale `sigma_ratio` by 1/√(ε_r μ_r)
/// when the layer is filled with another material.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct CpmlParams {
    /// Number of PML cells on each side (≥ 1, typically 8–20).
    pub thickness: usize,
    /// Polynomial grading order m (finite, ≥ 1; typically 3–4).
    pub grading_order: f64,
    /// σ_max / σ_opt (finite, > 0; typically 0.7–1.5).
    pub sigma_ratio: f64,
    /// Maximum coordinate-stretching factor κ_max (finite, ≥ 1).
    pub kappa_max: f64,
    /// Maximum CFS shift α_max in S/m (finite, ≥ 0).
    pub alpha_max: f64,
}

impl CpmlParams {
    /// Classic CPML parameters (m = 3, σ_max = σ_opt, κ_max = 1, α_max = 0)
    /// for a layer `thickness` cells thick.
    #[inline]
    #[must_use]
    pub fn new(thickness: usize) -> Self {
        Self {
            thickness,
            grading_order: 3.0,
            sigma_ratio: 1.0,
            kappa_max: 1.0,
            alpha_max: 0.0,
        }
    }

    /// Check that every parameter is finite and in range.
    fn validate(&self) -> Result<()> {
        let invalid = |reason: String| Err(BijliError::InvalidParameter { reason });
        if self.thickness == 0 {
            return invalid("CPML thickness must be at least 1 cell".into());
        }
        if !(self.grading_order >= 1.0 && self.grading_order.is_finite()) {
            return invalid(format!(
                "CPML grading order must be finite and ≥ 1, got {}",
                self.grading_order
            ));
        }
        if !(self.sigma_ratio > 0.0 && self.sigma_ratio.is_finite()) {
            return invalid(format!(
                "CPML sigma_ratio must be finite and > 0, got {}",
                self.sigma_ratio
            ));
        }
        if !(self.kappa_max >= 1.0 && self.kappa_max.is_finite()) {
            return invalid(format!(
                "CPML kappa_max must be finite and ≥ 1, got {}",
                self.kappa_max
            ));
        }
        if !(self.alpha_max >= 0.0 && self.alpha_max.is_finite()) {
            return invalid(format!(
                "CPML alpha_max must be finite and ≥ 0, got {}",
                self.alpha_max
            ));
        }
        Ok(())
    }
}

/// Convolutional PML (CPML) absorbing layer for [`Fdtd2d`] and [`Fdtd3d`].
///
/// Implements the stretched-coordinate PML of Roden & Gedney (2000) with
/// recursive convolution. Inside the layer every spatial derivative ∂/∂u
/// normal to the layer becomes (1/κ_u) ∂/∂u + Ψ_u, with
///
/// Ψ_u ← b_u Ψ_u + c_u ∂/∂u,  b = exp(−(σ/κ + α) Δt/ε₀),
/// c = σ (b − 1) / (κ (σ + κα)),
///
/// where σ, κ, α are graded as described in [`CpmlParams`]. The profiles are
/// sampled at the true Yee positions of each corrected component along the
/// derivative axis — integer nodes for the E update in 3D and for the
/// out-of-plane field in 2D (E_z in TM, H_z in TE), half-cell nodes for the
/// H update in 3D and the in-plane components in 2D. Ψ is stored only inside
/// the PML: the four edge strips in 2D, the six face slabs in 3D (12 Ψ
/// fields, one per field component and normal derivative).
///
/// Normally created and attached through [`Fdtd2d::try_enable_cpml`] or
/// [`Fdtd3d::try_enable_cpml`]; [`Fdtd3d::cpml`] returns the layer with its
/// grading (its own 3D Ψ storage is internal to the simulation).
#[derive(Debug, Clone)]
pub struct Cpml {
    /// Number of PML cells on each side.
    pub thickness: usize,
    /// Polynomial grading order m.
    grading_order: f64,
    /// Maximum PML conductivity σ_max (S/m).
    sigma_max: f64,
    /// Maximum coordinate stretching κ_max.
    kappa_max: f64,
    /// Maximum CFS shift α_max (S/m).
    alpha_max: f64,
    /// b, c and 1/κ − 1 at integer node positions. Index `l` is the distance
    /// in cells from the outer wall (`l = 0` is the wall node itself).
    b_node: Vec<f64>,
    c_node: Vec<f64>,
    km1_node: Vec<f64>,
    /// b, c and 1/κ − 1 at half-cell positions, `l + ½` cells from the wall.
    b_half: Vec<f64>,
    c_half: Vec<f64>,
    km1_half: Vec<f64>,
    /// Grid size the auxiliary fields are allocated for (0 × 0 if unattached).
    nx: usize,
    ny: usize,
    /// Ψ of ∂f_z/∂y in the f_x update (y-strips): index `(side·t + l)·nx + x`.
    psi_fx_y: Vec<f64>,
    /// Ψ of ∂f_x/∂y in the f_z update (y-strips): index `(side·t + l)·nx + x`.
    psi_fz_y: Vec<f64>,
    /// Ψ of ∂f_z/∂x in the f_y update (x-strips): index `y·2t + side·t + l`.
    psi_fy_x: Vec<f64>,
    /// Ψ of ∂f_y/∂x in the f_z update (x-strips): index `y·2t + side·t + l`.
    psi_fz_x: Vec<f64>,
}

impl Cpml {
    /// Create a CPML layer with the classic grading of [`CpmlParams::new`].
    ///
    /// `thickness` is the number of PML cells on each side; `dx` and `dt` are
    /// the grid spacing and time step of the simulation it is meant for.
    ///
    /// Invalid input (zero thickness, non-finite or non-positive `dx`/`dt`)
    /// logs a warning and yields an inert layer that can never be attached;
    /// use [`Cpml::with_params`] to get an error instead.
    #[must_use]
    pub fn new(thickness: usize, dx: f64, dt: f64) -> Self {
        match Self::with_params(&CpmlParams::new(thickness), dx, dt) {
            Ok(cpml) => cpml,
            Err(e) => {
                tracing::warn!(thickness, dx, dt, error = %e, "invalid CPML parameters; layer is inert");
                Self::from_profiles(
                    &CpmlParams::new(thickness),
                    0.0,
                    Profile::default(),
                    Profile::default(),
                )
            }
        }
    }

    /// Create a CPML layer with explicit grading parameters.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if a parameter of `params` is out of
    /// range (see [`CpmlParams`]), if `dx` or `dt` is not finite and positive,
    /// or if the profiles cannot be allocated.
    pub fn with_params(params: &CpmlParams, dx: f64, dt: f64) -> Result<Self> {
        params.validate()?;
        if !(dx > 0.0 && dx.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("CPML dx must be finite and positive, got {dx}"),
            });
        }
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("CPML dt must be finite and positive, got {dt}"),
            });
        }
        let t = params.thickness;
        let m = params.grading_order;
        let sigma_max = params.sigma_ratio * 0.8 * (m + 1.0) / (ETA_0 * dx);

        // Integer nodes sit l cells from the wall, half nodes l + ½ cells.
        let node = Profile::sample(t, 0.0, params, sigma_max, dt)?;
        let half = Profile::sample(t, 0.5, params, sigma_max, dt)?;
        if !(node.is_finite() && half.is_finite() && sigma_max.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("CPML coefficients are not finite for dx = {dx}, dt = {dt}"),
            });
        }

        tracing::debug!(
            thickness = t,
            grading_order = m,
            sigma_max,
            kappa_max = params.kappa_max,
            alpha_max = params.alpha_max,
            "CPML profile computed"
        );
        Ok(Self::from_profiles(params, sigma_max, node, half))
    }

    /// Assemble an unattached layer from sampled profiles (empty profiles
    /// give an inert layer that `is_attached_to` always rejects).
    fn from_profiles(params: &CpmlParams, sigma_max: f64, node: Profile, half: Profile) -> Self {
        Self {
            thickness: params.thickness,
            grading_order: params.grading_order,
            sigma_max,
            kappa_max: params.kappa_max,
            alpha_max: params.alpha_max,
            b_node: node.b,
            c_node: node.c,
            km1_node: node.km1,
            b_half: half.b,
            c_half: half.c,
            km1_half: half.km1,
            nx: 0,
            ny: 0,
            psi_fx_y: Vec::new(),
            psi_fz_y: Vec::new(),
            psi_fy_x: Vec::new(),
            psi_fz_x: Vec::new(),
        }
    }

    /// Allocate (and zero) the auxiliary Ψ fields for an `nx` × `ny` 2D grid.
    ///
    /// Ψ is stored only inside the PML strips: 2·thickness rows of `nx` cells
    /// and 2·thickness columns of `ny` cells per field. If the sizes overflow
    /// or cannot be allocated, a warning is logged and the layer is left
    /// unattached. ([`Fdtd3d`] keeps its own face-slab storage and does not
    /// use this.)
    pub fn attach(&mut self, nx: usize, ny: usize) {
        let sizes = self
            .thickness
            .checked_mul(2)
            .and_then(|t2| Some((t2.checked_mul(nx)?, t2.checked_mul(ny)?)));
        let allocated = sizes.and_then(|(sx, sy)| {
            Some((
                try_filled(sx, 0.0).ok()?,
                try_filled(sx, 0.0).ok()?,
                try_filled(sy, 0.0).ok()?,
                try_filled(sy, 0.0).ok()?,
            ))
        });
        match allocated {
            Some((fx_y, fz_y, fy_x, fz_x)) => {
                self.psi_fx_y = fx_y;
                self.psi_fz_y = fz_y;
                self.psi_fy_x = fy_x;
                self.psi_fz_x = fz_x;
                self.nx = nx;
                self.ny = ny;
            }
            None => {
                tracing::warn!(
                    thickness = self.thickness,
                    nx,
                    ny,
                    "CPML auxiliary fields could not be allocated; layer left unattached"
                );
                self.psi_fx_y = Vec::new();
                self.psi_fz_y = Vec::new();
                self.psi_fy_x = Vec::new();
                self.psi_fz_x = Vec::new();
                self.nx = 0;
                self.ny = 0;
            }
        }
    }

    /// Maximum PML conductivity σ_max (S/m).
    #[inline]
    #[must_use]
    pub fn sigma_max(&self) -> f64 {
        self.sigma_max
    }

    /// Maximum coordinate stretching κ_max.
    #[inline]
    #[must_use]
    pub fn kappa_max(&self) -> f64 {
        self.kappa_max
    }

    /// Maximum CFS shift α_max (S/m).
    #[inline]
    #[must_use]
    pub fn alpha_max(&self) -> f64 {
        self.alpha_max
    }

    /// Polynomial grading order m.
    #[inline]
    #[must_use]
    pub fn grading_order(&self) -> f64 {
        self.grading_order
    }

    /// Whether the profiles and Ψ storage are consistent with an `nx` × `ny`
    /// grid; the strip updates rely on this (it bounds every index they use).
    fn is_attached_to(&self, nx: usize, ny: usize) -> bool {
        let t = self.thickness;
        if !cpml_fits(t, nx.min(ny)) || self.nx != nx || self.ny != ny || !self.has_profiles() {
            return false;
        }
        let (Some(sx), Some(sy)) = ((2 * t).checked_mul(nx), (2 * t).checked_mul(ny)) else {
            return false;
        };
        self.psi_fx_y.len() == sx
            && self.psi_fz_y.len() == sx
            && self.psi_fy_x.len() == sy
            && self.psi_fz_x.len() == sy
    }

    /// Whether every coefficient profile has one entry per layer (false for
    /// the inert layer built from invalid input).
    fn has_profiles(&self) -> bool {
        [
            &self.b_node,
            &self.c_node,
            &self.km1_node,
            &self.b_half,
            &self.c_half,
            &self.km1_half,
        ]
        .iter()
        .all(|v| v.len() == self.thickness)
    }

    /// One 3D CPML derivative term over both face slabs normal to `AXIS`
    /// (0 = x, 1 = y, 2 = z):
    ///
    /// Ψ ← b·Ψ + c·Δ,  f += sign·k·((1/κ − 1)·Δ + Ψ)
    ///
    /// With `HALF` the corrected field sits half a cell off the nodes along
    /// `AXIS` (the H update): Δ = src[u + 1] − src[u] and layer l = 0..t lies
    /// l + ½ cells from the wall. Otherwise (the E update) Δ = src[u] − src[u − 1]
    /// and layer l = 1..t lies on the integer node l cells from the wall
    /// (l = 0 is the PEC wall itself, reset after the update). `psi` has the
    /// grid's shape with the `AXIS` extent replaced by 2·t layers, indexed
    /// `side·t + l`, in the grid's x-fastest order.
    ///
    /// Caller guarantees `Cpml3d::is_attached_to(pass.dims)`, slice lengths
    /// `nx·ny·nz` and `pass.ranges` within the main-loop update ranges.
    #[inline]
    fn correct_slabs<const AXIS: usize, const HALF: bool>(
        &self,
        psi: &mut [f64],
        f: &mut [f64],
        src: &[f64],
        k: &[f64],
        pass: &SlabPass,
    ) {
        let [nx, ny, _] = pass.dims;
        let nxy = nx * ny;
        let n = pass.dims[AXIS];
        let t = self.thickness;
        let stride = [1, nx, nxy][AXIS];
        let (b, c, km1) = if HALF {
            (&self.b_half, &self.c_half, &self.km1_half)
        } else {
            (&self.b_node, &self.c_node, &self.km1_node)
        };
        // First layer, and the high-side position offset: layer l of the high
        // slab sits at index n − off − l along AXIS.
        let (l0, off) = if HALF { (0, 2) } else { (1, 1) };
        let mut pdims = pass.dims;
        pdims[AXIS] = 2 * t;
        let sign = pass.sign;

        for side in 0..2 {
            let layer = |u: usize| if side == 0 { u } else { n - off - u };
            let mut r = pass.ranges.clone();
            r[AXIS] = if side == 0 {
                l0..t
            } else {
                n + 1 - off - t..n + 1 - off - l0
            };
            for z in r[2].clone() {
                for y in r[1].clone() {
                    let row = z * nxy + y * nx;
                    if AXIS == 0 {
                        // x-slab: the layer changes along the row — ascending
                        // on the low side, descending on the high side.
                        let xs = row + r[0].start..row + r[0].end;
                        let (up, down) = if HALF {
                            (xs.start + 1..xs.end + 1, xs.clone())
                        } else {
                            (xs.clone(), xs.start - 1..xs.end - 1)
                        };
                        let prow = (z * pdims[1] + y) * pdims[0] + side * t;
                        let cells = f[xs.clone()]
                            .iter_mut()
                            .zip(&k[xs])
                            .zip(src[up].iter().zip(&src[down]));
                        let layers = b[l0..t]
                            .iter()
                            .zip(&c[l0..t])
                            .zip(&km1[l0..t])
                            .zip(psi[prow + l0..prow + t].iter_mut());
                        if side == 0 {
                            cpml_cells(cells, layers, sign);
                        } else {
                            cpml_cells(cells, layers.rev(), sign);
                        }
                    } else {
                        // y- or z-slab: one layer per row, contiguous in x.
                        let l = layer(if AXIS == 1 { y } else { z });
                        let s = side * t + l;
                        let (py, pz) = if AXIS == 1 { (s, z) } else { (y, s) };
                        let xs = r[0].clone();
                        let prow = (pz * pdims[1] + py) * pdims[0];
                        let (hi, lo) = if HALF {
                            (row + stride, row)
                        } else {
                            (row, row - stride)
                        };
                        let (bl, cl, kl) = (b[l], c[l], km1[l]);
                        let cells = f[row + xs.start..row + xs.end]
                            .iter_mut()
                            .zip(&k[row + xs.start..row + xs.end])
                            .zip(
                                src[hi + xs.start..hi + xs.end]
                                    .iter()
                                    .zip(&src[lo + xs.start..lo + xs.end]),
                            )
                            .zip(psi[prow + xs.start..prow + xs.end].iter_mut());
                        for (((fv, &kk), (&up, &down)), p) in cells {
                            let d = up - down;
                            *p = bl * *p + cl * d;
                            *fv += sign * kk * (kl * d + *p);
                        }
                    }
                }
            }
        }
    }

    /// CPML term for f_x at (x, y + ½) in the y-strips:
    /// f_x += sign·k·((1/κ_y − 1)·Δ_y f_z + Ψ).
    ///
    /// Caller guarantees `is_attached_to(nx, ny)` and slice lengths `nx·ny`.
    fn correct_fx(
        &mut self,
        fx: &mut [f64],
        fz: &[f64],
        k: &[f64],
        dims: (usize, usize),
        sign: f64,
    ) {
        let (nx, ny) = dims;
        let t = self.thickness;
        for side in 0..2 {
            for l in 0..t {
                let y = if side == 0 { l } else { ny - 2 - l };
                let (b, c, km1) = (self.b_half[l], self.c_half[l], self.km1_half[l]);
                let r0 = y * nx;
                let r1 = r0 + nx;
                let p0 = (side * t + l) * nx;
                let cells = fx[r0..r1]
                    .iter_mut()
                    .zip(&k[r0..r1])
                    .zip(fz[r1..r1 + nx].iter().zip(&fz[r0..r1]))
                    .zip(self.psi_fx_y[p0..p0 + nx].iter_mut());
                for (((f, &kk), (&up, &here)), psi) in cells {
                    let d = up - here;
                    *psi = b * *psi + c * d;
                    *f += sign * kk * (km1 * d + *psi);
                }
            }
        }
    }

    /// CPML term for f_y at (x + ½, y) in the x-strips:
    /// f_y += sign·k·((1/κ_x − 1)·Δ_x f_z + Ψ).
    ///
    /// Caller guarantees `is_attached_to(nx, ny)` and slice lengths `nx·ny`.
    fn correct_fy(&mut self, fy: &mut [f64], fz: &[f64], k: &[f64], nx: usize, sign: f64) {
        let t = self.thickness;
        let (b, c, km1) = (&self.b_half, &self.c_half, &self.km1_half);
        let rows = fy
            .chunks_exact_mut(nx)
            .zip(fz.chunks_exact(nx))
            .zip(k.chunks_exact(nx))
            .zip(self.psi_fy_x.chunks_exact_mut(2 * t));
        for (((fy_r, fz_r), k_r), psi_r) in rows {
            for l in 0..t {
                // Low side: node x = l (position l + ½); high side: x = nx − 2 − l.
                for (s, x) in [(l, l), (t + l, nx - 2 - l)] {
                    let d = fz_r[x + 1] - fz_r[x];
                    let psi = &mut psi_r[s];
                    *psi = b[l] * *psi + c[l] * d;
                    fy_r[x] += sign * k_r[x] * (km1[l] * d + *psi);
                }
            }
        }
    }

    /// CPML terms for f_z at integer nodes:
    /// f_z += sign·k·[((1/κ_x − 1)·Δ_x f_y + Ψ_x) − ((1/κ_y − 1)·Δ_y f_x + Ψ_y)].
    ///
    /// The wall nodes (l = 0) are skipped: they are reset by the boundary.
    /// Caller guarantees `is_attached_to(nx, ny)` and slice lengths `nx·ny`.
    fn correct_fz(
        &mut self,
        fz: &mut [f64],
        in_plane: (&[f64], &[f64]),
        k: &[f64],
        dims: (usize, usize),
        sign: f64,
    ) {
        let (fx, fy) = in_plane;
        let (nx, ny) = dims;
        let t = self.thickness;

        // x-strips: Δ_x f_y at nodes x = l and x = nx − 1 − l.
        {
            let (b, c, km1) = (&self.b_node, &self.c_node, &self.km1_node);
            let rows = fz
                .chunks_exact_mut(nx)
                .zip(fy.chunks_exact(nx))
                .zip(k.chunks_exact(nx))
                .zip(self.psi_fz_x.chunks_exact_mut(2 * t))
                .skip(1)
                .take(ny - 2);
            for (((fz_r, fy_r), k_r), psi_r) in rows {
                for l in 1..t {
                    for (s, x) in [(l, l), (t + l, nx - 1 - l)] {
                        let d = fy_r[x] - fy_r[x - 1];
                        let psi = &mut psi_r[s];
                        *psi = b[l] * *psi + c[l] * d;
                        fz_r[x] += sign * k_r[x] * (km1[l] * d + *psi);
                    }
                }
            }
        }

        // y-strips: Δ_y f_x at rows y = l and y = ny − 1 − l, interior columns.
        for side in 0..2 {
            for l in 1..t {
                let y = if side == 0 { l } else { ny - 1 - l };
                let (b, c, km1) = (self.b_node[l], self.c_node[l], self.km1_node[l]);
                let r0 = y * nx;
                let p0 = (side * t + l) * nx;
                let cells = fz[r0 + 1..r0 + nx - 1]
                    .iter_mut()
                    .zip(&k[r0 + 1..r0 + nx - 1])
                    .zip(fx[r0 + 1..r0 + nx - 1].iter().zip(&fx[r0 + 1 - nx..r0 - 1]))
                    .zip(self.psi_fz_y[p0 + 1..p0 + nx - 1].iter_mut());
                for (((f, &kk), (&here, &below)), psi) in cells {
                    let d = here - below;
                    *psi = b * *psi + c * d;
                    *f -= sign * kk * (km1 * d + *psi);
                }
            }
        }
    }
}

/// CPML update of consecutive cells paired with their layers:
/// Ψ ← b·Ψ + c·Δ, f += sign·k·((1/κ − 1)·Δ + Ψ), Δ = up − down.
#[inline]
fn cpml_cells<'a>(
    cells: impl Iterator<Item = ((&'a mut f64, &'a f64), (&'a f64, &'a f64))>,
    layers: impl Iterator<Item = (((&'a f64, &'a f64), &'a f64), &'a mut f64)>,
    sign: f64,
) {
    for (((f, &k), (&up, &down)), (((&b, &c), &km1), psi)) in cells.zip(layers) {
        let d = up - down;
        *psi = b * *psi + c * d;
        *f += sign * k * (km1 * d + *psi);
    }
}

/// CPML coefficient profile `(b, c, 1/κ − 1)` sampled per layer.
#[derive(Debug, Default)]
struct Profile {
    b: Vec<f64>,
    c: Vec<f64>,
    km1: Vec<f64>,
}

impl Profile {
    /// Sample layers `l = 0..t` located `l + offset` cells from the outer
    /// wall, i.e. at normalised depth ρ = (t − l − offset)/t.
    fn sample(t: usize, offset: f64, params: &CpmlParams, sigma_max: f64, dt: f64) -> Result<Self> {
        let mut p = Self {
            b: try_filled(t, 0.0)?,
            c: try_filled(t, 0.0)?,
            km1: try_filled(t, 0.0)?,
        };
        let tf = t as f64;
        let layers = p.b.iter_mut().zip(p.c.iter_mut()).zip(p.km1.iter_mut());
        for (l, ((b, c), km1)) in layers.enumerate() {
            let rho = (tf - l as f64 - offset) / tf;
            (*b, *c, *km1) = cpml_coefficients(rho, params, sigma_max, dt);
        }
        Ok(p)
    }

    /// Whether every coefficient is finite.
    fn is_finite(&self) -> bool {
        self.b
            .iter()
            .chain(&self.c)
            .chain(&self.km1)
            .all(|v| v.is_finite())
    }
}

/// CPML recursive-convolution coefficients `(b, c, 1/κ − 1)` at normalised
/// depth `rho` ∈ [0, 1].
#[inline]
fn cpml_coefficients(rho: f64, params: &CpmlParams, sigma_max: f64, dt: f64) -> (f64, f64, f64) {
    let grade = rho.powf(params.grading_order);
    let sigma = sigma_max * grade;
    let kappa = 1.0 + (params.kappa_max - 1.0) * grade;
    let alpha = params.alpha_max * (1.0 - rho);
    let b = (-(sigma / kappa + alpha) * dt / EPSILON_0).exp();
    let denom = kappa * (sigma + kappa * alpha);
    let c = if sigma > 0.0 && denom > 0.0 {
        sigma * (b - 1.0) / denom
    } else {
        0.0
    };
    (b, c, 1.0 / kappa - 1.0)
}

// ── Grid helpers (shared by Fdtd2d / Fdtd3d) ─────────────────────

/// Allocate `n` copies of `value`, returning an error instead of aborting or
/// panicking when the allocation is impossible.
fn try_filled(n: usize, value: f64) -> Result<Vec<f64>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| BijliError::InvalidParameter {
            reason: format!("cannot allocate {n} grid values"),
        })?;
    v.resize(n, value);
    Ok(v)
}

/// Total number of cells of a grid, checked against overflow.
fn grid_cells(dims: &[usize]) -> Result<usize> {
    dims.iter()
        .try_fold(1_usize, |acc, &d| acc.checked_mul(d))
        .ok_or_else(|| BijliError::InvalidParameter {
            reason: format!("grid {dims:?} is too large"),
        })
}

/// Validate a grid spacing: finite and strictly positive.
#[inline]
fn validate_dx(dx: f64) -> Result<()> {
    if dx > 0.0 && dx.is_finite() {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!("dx must be finite and positive, got {dx}"),
        })
    }
}

/// Validate a relative permittivity: finite and strictly positive.
#[inline]
fn validate_eps_r(eps_r: f64) -> Result<()> {
    if eps_r > 0.0 && eps_r.is_finite() {
        Ok(())
    } else {
        Err(BijliError::InvalidPermittivity { value: eps_r })
    }
}

/// Validate a relative permeability: finite and strictly positive.
#[inline]
fn validate_mu_r(mu_r: f64) -> Result<()> {
    if mu_r > 0.0 && mu_r.is_finite() {
        Ok(())
    } else {
        Err(BijliError::InvalidPermeability { value: mu_r })
    }
}

/// Validate an electric conductivity: finite and non-negative (S/m).
#[inline]
fn validate_conductivity(sigma: f64) -> Result<()> {
    if sigma >= 0.0 && sigma.is_finite() {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!("conductivity must be finite and ≥ 0 S/m, got {sigma}"),
        })
    }
}

/// Validate the half-open index range `lo..hi` along an axis of `n` cells.
///
/// `lo` must lie inside the grid and `hi ≥ lo`; `hi` is clamped to `n`.
fn region_range(axis: &str, lo: usize, hi: usize, n: usize) -> Result<std::ops::Range<usize>> {
    if lo >= n || hi < lo {
        return Err(BijliError::InvalidParameter {
            reason: format!("region {axis} range {lo}..{hi} is invalid for a grid of {n} cells"),
        });
    }
    Ok(lo..hi.min(n))
}

/// E-update coefficients `(C_a, C_b)` for relative permittivity `eps_r` and
/// conductivity `sigma`:
/// C_a = (1 − σΔt/2ε)/(1 + σΔt/2ε), C_b = (Δt/(εΔx))/(1 + σΔt/2ε).
///
/// The lossless case returns exactly `(1, Δt/(ε₀Δx)/ε_r)`.
#[inline]
fn e_update_coefficients(dt: f64, dx: f64, eps_r: f64, sigma: f64) -> (f64, f64) {
    if sigma > 0.0 {
        let eps = eps_r * EPSILON_0;
        let half_loss = sigma * dt / (2.0 * eps);
        (
            (1.0 - half_loss) / (1.0 + half_loss),
            (dt / (eps * dx)) / (1.0 + half_loss),
        )
    } else {
        (1.0, dt / (EPSILON_0 * dx) / eps_r)
    }
}

/// `dst[i] ±= k[i]·(a[i] − b[i])` over the common length of the slices.
#[inline]
fn row_diff<const ADD: bool>(dst: &mut [f64], k: &[f64], a: &[f64], b: &[f64]) {
    for (((d, &k), &a), &b) in dst.iter_mut().zip(k).zip(a).zip(b) {
        let v = k * (a - b);
        if ADD {
            *d += v;
        } else {
            *d -= v;
        }
    }
}

/// Interior curl update of one row, x ∈ 1..n−1:
/// `dst[x] ±= k[x]·((a[x] − a[x−1]) − (b[x] − b_prev[x]))`.
#[inline]
fn row_curl<const ADD: bool>(dst: &mut [f64], k: &[f64], a: &[f64], b: &[f64], b_prev: &[f64]) {
    let n = dst.len();
    if n < 3 {
        return;
    }
    let (Some(d), Some(k), Some(a1), Some(b1), Some(b0)) = (
        dst.get_mut(1..n - 1),
        k.get(1..),
        a.get(1..),
        b.get(1..),
        b_prev.get(1..),
    ) else {
        return;
    };
    let cells = d
        .iter_mut()
        .zip(k)
        .zip(a1.iter().zip(a))
        .zip(b1.iter().zip(b0));
    for (((d, &k), (&a1, &a0)), (&b1, &b0)) in cells {
        let v = k * ((a1 - a0) - (b1 - b0));
        if ADD {
            *d += v;
        } else {
            *d -= v;
        }
    }
}

/// `f[i] *= factor[i]` (no-op when `factor` is empty, i.e. lossless).
#[inline]
fn scale_by(f: &mut [f64], factor: &[f64]) {
    for (v, &a) in f.iter_mut().zip(factor) {
        *v *= a;
    }
}

/// Zero the outermost ring of a row-major `nx`-wide 2D field.
#[inline]
fn zero_walls(f: &mut [f64], nx: usize) {
    if nx == 0 {
        return;
    }
    let len = f.len();
    if let Some(row) = f.get_mut(..nx) {
        row.fill(0.0);
    }
    if let Some(row) = f.get_mut(len.saturating_sub(nx)..) {
        row.fill(0.0);
    }
    for row in f.chunks_exact_mut(nx) {
        if let Some(v) = row.first_mut() {
            *v = 0.0;
        }
        if let Some(v) = row.last_mut() {
            *v = 0.0;
        }
    }
}

// ── 2D FDTD ──────────────────────────────────────────────────────

/// Polarization mode for 2D FDTD.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Mode2d {
    /// TM mode: E_z, H_x, H_y (E perpendicular to simulation plane).
    /// Default walls are PEC (E_z = 0).
    Tm,
    /// TE mode: H_z, E_x, E_y (H perpendicular to simulation plane).
    /// Default walls are PMC (H_z = 0).
    Te,
}

/// 2D FDTD simulation on a Yee grid.
///
/// - TM mode: E_z at integer nodes (x, y), H_x at (x, y + ½), H_y at (x + ½, y).
/// - TE mode: H_z at integer nodes (x, y), E_x at (x, y + ½), E_y at (x + ½, y).
///
/// Fields are stored in row-major order: index = y * nx + x.
///
/// # Boundaries
///
/// By default the outer walls are perfectly **reflecting**: the out-of-plane
/// field is held at zero on the outermost ring of nodes, which is a perfect
/// electric conductor (PEC) in TM mode (E_z = 0) and a perfect magnetic
/// conductor (PMC) in TE mode (H_z = 0). Call [`Fdtd2d::try_enable_cpml`] to
/// absorb outgoing waves instead: the CPML occupies the outermost `thickness`
/// cells of the grid and is terminated by the same wall.
///
/// # Materials
///
/// Per-cell ε_r, μ_r and electric conductivity σ (lossy E update with
/// C_a/C_b coefficients) through [`Fdtd2d::set_material`],
/// [`Fdtd2d::set_permittivity`], [`Fdtd2d::set_permeability`] and
/// [`Fdtd2d::set_conductivity`].
///
/// # Public fields
///
/// `nx`, `ny`, `dx`, `dt` and the field vectors are public for reading
/// (probes, visualisation) and for writing field *values* (initial
/// conditions). Treat the geometry and the vector lengths as read-only: the
/// update coefficients are precomputed from them. If they are changed so the
/// layout is inconsistent, stepping logs an error and leaves the state
/// untouched instead of panicking; changing `dx` or `dt` alone silently
/// invalidates the precomputed coefficients.
#[derive(Debug, Clone)]
pub struct Fdtd2d {
    /// Polarization mode.
    pub mode: Mode2d,
    /// Number of nodes along x (read-only in practice — see the type docs).
    pub nx: usize,
    /// Number of nodes along y (read-only in practice — see the type docs).
    pub ny: usize,
    /// Spatial step (uniform, m). Read-only in practice.
    pub dx: f64,
    /// Time step (s). Read-only in practice.
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

    /// E-field curl coefficient C_b per cell: Δt/(ε₀ε_r Δx), divided by
    /// (1 + σΔt/2ε) in lossy cells.
    ce: Vec<f64>,
    /// H-field update coefficient per cell: dt/(μ₀μ_r dx).
    ch: Vec<f64>,
    /// E-field decay coefficient C_a per cell; empty while all cells are lossless.
    ca: Vec<f64>,
    /// Relative permittivity per cell.
    permittivity: Vec<f64>,
    /// Relative permeability per cell.
    permeability: Vec<f64>,
    /// Electric conductivity σ per cell (S/m); empty while all cells are lossless.
    conductivity: Vec<f64>,
    /// Optional CPML absorbing boundaries.
    cpml: Option<Cpml>,
}

impl Fdtd2d {
    /// Create a new 2D FDTD simulation (vacuum, reflecting walls — see the
    /// type-level docs).
    ///
    /// Uses uniform grid spacing `dx` in both x and y. The time step is fixed
    /// at dt = Δx / (2√2 c), i.e. a Courant number S = c·dt/Δx = 1/(2√2) ≈ 0.354,
    /// half of the 2D stability limit 1/√2.
    ///
    /// # Errors
    ///
    /// - [`BijliError::InsufficientResolution`] if `nx < 3` or `ny < 3`.
    /// - [`BijliError::InvalidParameter`] if `dx` is not finite and positive,
    ///   or the grid is too large to allocate.
    pub fn new(nx: usize, ny: usize, dx: f64, mode: Mode2d) -> Result<Self> {
        if nx < 3 || ny < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx.min(ny),
                domain_size: nx.min(ny) as f64 * dx,
            });
        }
        validate_dx(dx)?;
        let n = grid_cells(&[nx, ny])?;

        // 2D Courant: dt ≤ dx/(c√2). Use half for safety (S = 1/(2√2)).
        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);

        let ce = try_filled(n, dt / (EPSILON_0 * dx))?;
        let ch = try_filled(n, dt / (MU_0 * dx))?;
        let fz = try_filled(n, 0.0)?;
        let fx = try_filled(n, 0.0)?;
        let fy = try_filled(n, 0.0)?;
        let permittivity = try_filled(n, 1.0)?;
        let permeability = try_filled(n, 1.0)?;

        tracing::debug!(nx, ny, dx, dt, ?mode, "2D FDTD simulation created");

        Ok(Self {
            mode,
            nx,
            ny,
            dx,
            dt,
            step: 0,
            fz,
            fx,
            fy,
            ce,
            ch,
            ca: Vec::new(),
            permittivity,
            permeability,
            conductivity: Vec::new(),
            cpml: None,
        })
    }

    /// Linear index from (x, y) coordinates.
    #[inline]
    #[must_use]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Check that the grid geometry and every per-cell vector agree.
    fn check_grid(&self) -> Result<usize> {
        let n = self
            .nx
            .checked_mul(self.ny)
            .filter(|_| self.nx >= 3 && self.ny >= 3);
        let consistent = n.is_some_and(|n| {
            [
                &self.fz,
                &self.fx,
                &self.fy,
                &self.ce,
                &self.ch,
                &self.permittivity,
                &self.permeability,
            ]
            .iter()
            .all(|v| v.len() == n)
                && ((self.ca.is_empty() && self.conductivity.is_empty())
                    || (self.ca.len() == n && self.conductivity.len() == n))
        });
        match n {
            Some(n) if consistent => Ok(n),
            _ => Err(BijliError::InvalidParameter {
                reason: format!(
                    "inconsistent 2D grid state: nx = {}, ny = {}, field lengths = ({}, {}, {})",
                    self.nx,
                    self.ny,
                    self.fz.len(),
                    self.fx.len(),
                    self.fy.len()
                ),
            }),
        }
    }

    /// [`Self::check_grid`] plus consistency of the attached CPML.
    fn check_layout(&self) -> Result<usize> {
        let n = self.check_grid()?;
        if self
            .cpml
            .as_ref()
            .is_some_and(|c| !c.is_attached_to(self.nx, self.ny))
        {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "CPML is not attached to the current {}×{} grid",
                    self.nx, self.ny
                ),
            });
        }
        Ok(n)
    }

    /// Write material values to an already validated region and refresh the
    /// update coefficients. `None` leaves that property unchanged.
    fn apply_region(
        &mut self,
        xs: std::ops::Range<usize>,
        ys: std::ops::Range<usize>,
        eps_r: Option<f64>,
        mu_r: Option<f64>,
        sigma: Option<f64>,
    ) -> Result<()> {
        let n = self.check_grid()?;
        if sigma.is_some_and(|s| s > 0.0) && self.ca.is_empty() {
            let ca = try_filled(n, 1.0)?;
            let conductivity = try_filled(n, 0.0)?;
            self.ca = ca;
            self.conductivity = conductivity;
        }
        tracing::debug!(
            ?xs,
            ?ys,
            ?eps_r,
            ?mu_r,
            ?sigma,
            "2D FDTD material region set"
        );
        let (dt, dx, nx) = (self.dt, self.dx, self.nx);
        let base_h = dt / (MU_0 * dx);
        for y in ys {
            for x in xs.clone() {
                let i = y * nx + x;
                if let Some(mu) = mu_r {
                    self.permeability[i] = mu;
                    self.ch[i] = base_h / mu;
                }
                if eps_r.is_none() && sigma.is_none() {
                    continue;
                }
                if let Some(eps) = eps_r {
                    self.permittivity[i] = eps;
                }
                if let (Some(s), Some(c)) = (sigma, self.conductivity.get_mut(i)) {
                    *c = s;
                }
                let s = self.conductivity.get(i).copied().unwrap_or(0.0);
                let (ca, cb) = e_update_coefficients(dt, dx, self.permittivity[i], s);
                self.ce[i] = cb;
                if let Some(a) = self.ca.get_mut(i) {
                    *a = ca;
                }
            }
        }
        Ok(())
    }

    /// Validate a rectangular region `[x0, x1) × [y0, y1)`.
    fn region(
        &self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
    ) -> Result<(std::ops::Range<usize>, std::ops::Range<usize>)> {
        Ok((
            region_range("x", x0, x1, self.nx)?,
            region_range("y", y0, y1, self.ny)?,
        ))
    }

    /// Set relative permittivity in the rectangular region `[x0, x1) × [y0, y1)`.
    ///
    /// The exclusive end coordinates are clamped to the grid. Any
    /// conductivity already set in the region is kept.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidPermittivity`] if `eps_r` is not finite and
    /// positive; [`BijliError::InvalidParameter`] if the region starts
    /// outside the grid or is inverted (`x1 < x0`, `y1 < y0`). Nothing is
    /// modified on error.
    pub fn set_permittivity(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        eps_r: f64,
    ) -> Result<()> {
        validate_eps_r(eps_r)?;
        let (xs, ys) = self.region(x0, y0, x1, y1)?;
        self.apply_region(xs, ys, Some(eps_r), None, None)
    }

    /// Set relative permeability in the rectangular region `[x0, x1) × [y0, y1)`.
    ///
    /// The exclusive end coordinates are clamped to the grid.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidPermeability`] if `mu_r` is not finite and
    /// positive; [`BijliError::InvalidParameter`] for a region that starts
    /// outside the grid or is inverted. Nothing is modified on error.
    pub fn set_permeability(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        mu_r: f64,
    ) -> Result<()> {
        validate_mu_r(mu_r)?;
        let (xs, ys) = self.region(x0, y0, x1, y1)?;
        self.apply_region(xs, ys, None, Some(mu_r), None)
    }

    /// Set electric conductivity σ (S/m) in the rectangular region
    /// `[x0, x1) × [y0, y1)`.
    ///
    /// Lossy cells use the standard semi-implicit E update
    /// E ← C_a E + C_b ∇×H with C_a = (1 − σΔt/2ε)/(1 + σΔt/2ε) and
    /// C_b = (Δt/εΔx)/(1 + σΔt/2ε). The exclusive end coordinates are clamped
    /// to the grid. The first lossy region allocates the per-cell C_a and σ
    /// arrays; fully lossless grids keep the cheaper lossless update.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `sigma` is negative or not finite,
    /// or the region starts outside the grid or is inverted. Nothing is
    /// modified on error.
    pub fn set_conductivity(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        sigma: f64,
    ) -> Result<()> {
        validate_conductivity(sigma)?;
        let (xs, ys) = self.region(x0, y0, x1, y1)?;
        self.apply_region(xs, ys, None, None, Some(sigma))
    }

    /// Set material properties in a rectangular region using a [`crate::material::Material`] struct.
    ///
    /// Sets permittivity, permeability and conductivity from the material
    /// (see [`Fdtd2d::set_conductivity`] for the lossy update). All inputs are
    /// validated before anything is modified.
    ///
    /// # Errors
    ///
    /// Same as [`Fdtd2d::set_permittivity`], [`Fdtd2d::set_permeability`] and
    /// [`Fdtd2d::set_conductivity`].
    pub fn set_material(
        &mut self,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        mat: &crate::material::Material,
    ) -> Result<()> {
        validate_eps_r(mat.eps_r)?;
        validate_mu_r(mat.mu_r)?;
        validate_conductivity(mat.conductivity)?;
        let (xs, ys) = self.region(x0, y0, x1, y1)?;
        self.apply_region(
            xs,
            ys,
            Some(mat.eps_r),
            Some(mat.mu_r),
            Some(mat.conductivity),
        )
    }

    /// Inject a soft source: adds `value` to the out-of-plane field (E_z in
    /// TM, H_z in TE) at grid point (x, y).
    ///
    /// Points outside the grid are ignored (no panic). Values injected on the
    /// outermost ring of nodes are cleared by the wall on the next step.
    #[inline]
    pub fn add_source(&mut self, x: usize, y: usize, value: f64) {
        if x < self.nx
            && y < self.ny
            && let Some(v) = y
                .checked_mul(self.nx)
                .and_then(|row| row.checked_add(x))
                .and_then(|i| self.fz.get_mut(i))
        {
            *v += value;
        }
    }

    /// Enable CPML absorbing boundaries with the classic grading
    /// ([`CpmlParams::new`]) and the given thickness (cells per side).
    ///
    /// Replaces the reflecting PEC (TM) / PMC (TE) walls by a convolutional
    /// PML occupying the outermost `thickness` cells. Typical thickness:
    /// 8–20 cells. Calling it again replaces the layer (Ψ restarts at zero).
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `thickness == 0` or
    /// `2·thickness + 3 > min(nx, ny)`, or the grid state is inconsistent.
    /// The current boundary configuration is unchanged on error.
    pub fn try_enable_cpml(&mut self, thickness: usize) -> Result<()> {
        self.try_enable_cpml_with(&CpmlParams::new(thickness))
    }

    /// Enable CPML absorbing boundaries with explicit grading parameters
    /// (κ and CFS-α profiles, grading order, σ_max scale).
    ///
    /// # Errors
    ///
    /// As [`Fdtd2d::try_enable_cpml`], plus any parameter rejected by
    /// [`Cpml::with_params`]. The current boundary configuration is unchanged
    /// on error.
    ///
    /// # Example
    ///
    /// ```
    /// use bijli::fdtd::{CpmlParams, Fdtd2d, Mode2d};
    ///
    /// let mut sim = Fdtd2d::new(100, 100, 1e-4, Mode2d::Tm)?;
    /// let mut params = CpmlParams::new(12);
    /// params.kappa_max = 5.0; // stronger absorption of evanescent fields
    /// sim.try_enable_cpml_with(&params)?;
    /// sim.run(100, Some((50, 50)), 1e10, 1.0);
    /// # Ok::<(), bijli::error::BijliError>(())
    /// ```
    pub fn try_enable_cpml_with(&mut self, params: &CpmlParams) -> Result<()> {
        self.check_grid()?;
        let (nx, ny) = (self.nx, self.ny);
        let t = params.thickness;
        if !cpml_fits(t, nx.min(ny)) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "CPML thickness {t} is invalid for a {nx}×{ny} grid \
                     (need thickness ≥ 1 and 2·thickness + {CPML_MIN_INTERIOR} ≤ min(nx, ny))"
                ),
            });
        }
        let mut cpml = Cpml::with_params(params, self.dx, self.dt)?;
        cpml.attach(nx, ny);
        if !cpml.is_attached_to(nx, ny) {
            return Err(BijliError::InvalidParameter {
                reason: format!("cannot allocate CPML auxiliary fields for a {nx}×{ny} grid"),
            });
        }
        tracing::debug!(
            thickness = t,
            sigma_max = cpml.sigma_max,
            kappa_max = cpml.kappa_max,
            alpha_max = cpml.alpha_max,
            mode = ?self.mode,
            "2D CPML enabled"
        );
        self.cpml = Some(cpml);
        Ok(())
    }

    /// Enable CPML absorbing boundaries with given thickness (number of cells).
    ///
    /// Convenience wrapper around [`Fdtd2d::try_enable_cpml`]: on invalid
    /// input (`thickness == 0` or `2·thickness + 3 > min(nx, ny)`) it logs a
    /// warning and leaves the current boundaries unchanged. Prefer
    /// `try_enable_cpml` to handle the error.
    /// Typical thickness: 8–20 cells.
    pub fn enable_cpml(&mut self, thickness: usize) {
        if let Err(e) = self.try_enable_cpml(thickness) {
            tracing::warn!(thickness, error = %e, "enable_cpml ignored; boundaries unchanged");
        }
    }

    /// Remove the CPML and return to reflecting PEC (TM) / PMC (TE) walls.
    pub fn disable_cpml(&mut self) {
        if self.cpml.take().is_some() {
            tracing::debug!("2D CPML disabled");
        }
    }

    /// The active CPML layer, if any.
    #[inline]
    #[must_use]
    pub fn cpml(&self) -> Option<&Cpml> {
        self.cpml.as_ref()
    }

    /// Advance the simulation by one time step.
    ///
    /// Applies the CPML if enabled; the outermost ring is always the PEC (TM)
    /// / PMC (TE) wall. If the public geometry fields or vector lengths were
    /// changed into an inconsistent state, logs an error and does nothing.
    pub fn step_once(&mut self) {
        if let Err(e) = self.check_layout() {
            tracing::error!(error = %e, step = self.step, "Fdtd2d step skipped");
            return;
        }
        match self.mode {
            Mode2d::Tm => self.step_tm(),
            Mode2d::Te => self.step_te(),
        }
        self.step += 1;
    }

    /// TM mode update: E_z, H_x, H_y. Requires a validated layout.
    #[inline]
    fn step_tm(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let Self {
            fz,
            fx,
            fy,
            ce,
            ch,
            ca,
            cpml,
            ..
        } = self;

        // H_x at (x, y+½): H_x −= C_h·(E_z[y+1] − E_z[y])
        let rows = fx
            .chunks_exact_mut(nx)
            .zip(ch.chunks_exact(nx))
            .zip(fz.chunks_exact(nx).zip(fz.chunks_exact(nx).skip(1)));
        for ((hx, k), (ez, ez_up)) in rows {
            row_diff::<false>(hx, k, ez_up, ez);
        }

        // H_y at (x+½, y): H_y += C_h·(E_z[x+1] − E_z[x])
        let rows = fy
            .chunks_exact_mut(nx)
            .zip(ch.chunks_exact(nx))
            .zip(fz.chunks_exact(nx));
        for ((hy, k), ez) in rows {
            row_diff::<true>(hy, k, ez.get(1..).unwrap_or(&[]), ez);
        }

        if let Some(p) = cpml.as_mut() {
            p.correct_fx(fx, fz, ch, (nx, ny), -1.0);
            p.correct_fy(fy, fz, ch, nx, 1.0);
        }

        // Lossy cells: E_z ← C_a·E_z before adding the curl term.
        scale_by(fz, ca);

        // E_z at (x, y): E_z += C_b·((H_y[x] − H_y[x−1]) − (H_x[y] − H_x[y−1]))
        let rows = fz
            .chunks_exact_mut(nx)
            .skip(1)
            .zip(ce.chunks_exact(nx).skip(1))
            .zip(fy.chunks_exact(nx).skip(1))
            .zip(fx.chunks_exact(nx).skip(1).zip(fx.chunks_exact(nx)))
            .take(ny - 2);
        for (((ez, k), hy), (hx, hx_down)) in rows {
            row_curl::<true>(ez, k, hy, hx, hx_down);
        }

        if let Some(p) = cpml.as_mut() {
            p.correct_fz(fz, (fx, fy), ce, (nx, ny), 1.0);
        }

        // PEC walls: tangential E_z = 0 on the outermost ring (fully reflecting).
        zero_walls(fz, nx);
    }

    /// TE mode update: H_z, E_x, E_y. Requires a validated layout.
    #[inline]
    fn step_te(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let Self {
            fz,
            fx,
            fy,
            ce,
            ch,
            ca,
            cpml,
            ..
        } = self;

        // Lossy cells: E ← C_a·E before adding the curl term.
        scale_by(fx, ca);
        scale_by(fy, ca);

        // E_x at (x, y+½): E_x += C_b·(H_z[y+1] − H_z[y])
        let rows = fx
            .chunks_exact_mut(nx)
            .zip(ce.chunks_exact(nx))
            .zip(fz.chunks_exact(nx).zip(fz.chunks_exact(nx).skip(1)));
        for ((ex, k), (hz, hz_up)) in rows {
            row_diff::<true>(ex, k, hz_up, hz);
        }

        // E_y at (x+½, y): E_y −= C_b·(H_z[x+1] − H_z[x])
        let rows = fy
            .chunks_exact_mut(nx)
            .zip(ce.chunks_exact(nx))
            .zip(fz.chunks_exact(nx));
        for ((ey, k), hz) in rows {
            row_diff::<false>(ey, k, hz.get(1..).unwrap_or(&[]), hz);
        }

        if let Some(p) = cpml.as_mut() {
            p.correct_fx(fx, fz, ce, (nx, ny), 1.0);
            p.correct_fy(fy, fz, ce, nx, -1.0);
        }

        // H_z at (x, y): H_z −= C_h·((E_y[x] − E_y[x−1]) − (E_x[y] − E_x[y−1]))
        let rows = fz
            .chunks_exact_mut(nx)
            .skip(1)
            .zip(ch.chunks_exact(nx).skip(1))
            .zip(fy.chunks_exact(nx).skip(1))
            .zip(fx.chunks_exact(nx).skip(1).zip(fx.chunks_exact(nx)))
            .take(ny - 2);
        for (((hz, k), ey), (ex, ex_down)) in rows {
            row_curl::<false>(hz, k, ey, ex, ex_down);
        }

        if let Some(p) = cpml.as_mut() {
            p.correct_fz(fz, (fx, fy), ch, (nx, ny), -1.0);
        }

        // PMC walls: tangential H_z = 0 on the outermost ring (fully reflecting).
        zero_walls(fz, nx);
    }

    /// Run the simulation for `steps` time steps with optional sinusoidal source.
    ///
    /// The source (soft, see [`Fdtd2d::add_source`]) is ignored if outside the
    /// grid. Does nothing (logging an error) if the grid state is inconsistent.
    pub fn run(
        &mut self,
        steps: usize,
        source: Option<(usize, usize)>,
        frequency: f64,
        amplitude: f64,
    ) {
        tracing::debug!(
            steps,
            ?source,
            frequency,
            amplitude,
            cpml = ?self.cpml.as_ref().map(|c| c.thickness),
            "2D FDTD run started"
        );
        if let Err(e) = self.check_layout() {
            tracing::error!(error = %e, "Fdtd2d run aborted");
            return;
        }
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

    /// Total electromagnetic energy in the grid (including the CPML region).
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        let area = self.dx * self.dx;
        let (z_mat, z_const, xy_mat, xy_const) = match self.mode {
            Mode2d::Tm => (&self.permittivity, EPSILON_0, &self.permeability, MU_0),
            Mode2d::Te => (&self.permeability, MU_0, &self.permittivity, EPSILON_0),
        };
        let mut energy = 0.0;
        let cells = self
            .fz
            .iter()
            .zip(&self.fx)
            .zip(&self.fy)
            .zip(z_mat.iter().zip(xy_mat));
        for (((&fz, &fx), &fy), (&mz, &mxy)) in cells {
            energy += 0.5 * mz * z_const * fz * fz;
            energy += 0.5 * mxy * xy_const * (fx * fx + fy * fy);
        }
        energy * area
    }
}

#[cfg(test)]
mod cpml_tests {
    use super::*;
    use crate::material::Material;

    const DX: f64 = 1e-4;
    /// CPML test grid (n × n) and layer thickness.
    const N: usize = 81;
    const T: usize = 10;
    /// Reference grid: large enough that no wall reflection reaches any probe
    /// within `STEPS` (S ≈ 0.354 cells/step).
    const N_REF: usize = 201;
    const STEPS: usize = 450;

    /// Probe offsets from the source (grid centre): just inside the CPML
    /// interface on the left and top, near a corner, and an oblique point.
    const PROBES: [(isize, isize); 4] = [(-28, 0), (0, 28), (-28, -28), (-15, -26)];

    /// Soft Gaussian point source at the grid centre; record f_z at the probes.
    fn record(mode: Mode2d, n: usize, cpml: Option<&CpmlParams>) -> Vec<Vec<f64>> {
        let mut sim = Fdtd2d::new(n, n, DX, mode).unwrap();
        if let Some(p) = cpml {
            sim.try_enable_cpml_with(p).unwrap();
        }
        let c = (n / 2) as isize;
        let probes: Vec<usize> = PROBES
            .iter()
            .map(|&(ox, oy)| sim.idx((c + ox) as usize, (c + oy) as usize))
            .collect();
        let tau = 25.0 * sim.dt;
        let t0 = 4.0 * tau;
        let mut out = vec![Vec::with_capacity(STEPS); probes.len()];
        for _ in 0..STEPS {
            sim.add_source(n / 2, n / 2, gaussian_pulse(sim.time(), t0, tau));
            sim.step_once();
            for (series, &i) in out.iter_mut().zip(&probes) {
                series.push(sim.fz[i]);
            }
        }
        out
    }

    /// Worst-case reflection over all probes: max|f − f_ref| / max|f_ref|.
    fn reflection(test: &[Vec<f64>], reference: &[Vec<f64>]) -> f64 {
        test.iter()
            .zip(reference)
            .map(|(a, r)| {
                let peak = r.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                let err = a
                    .iter()
                    .zip(r)
                    .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
                err / peak
            })
            .fold(0.0, f64::max)
    }

    fn db(r: f64) -> f64 {
        20.0 * r.log10()
    }

    fn check_reflection(mode: Mode2d) {
        let reference = record(mode, N_REF, None);
        let pec = reflection(&record(mode, N, None), &reference);
        let classic = reflection(&record(mode, N, Some(&CpmlParams::new(T))), &reference);
        let mut cfs = CpmlParams::new(T);
        cfs.kappa_max = 5.0;
        cfs.alpha_max = 2.0 * std::f64::consts::PI * EPSILON_0 * 1e9;
        let cfs_r = reflection(&record(mode, N, Some(&cfs)), &reference);
        eprintln!(
            "{mode:?}: wall {:.1} dB, CPML({T}) {:.1} dB, CFS-CPML({T}) {:.1} dB",
            db(pec),
            db(classic),
            db(cfs_r)
        );
        // Without CPML the wall reflects everything back into the probes.
        assert!(pec > 0.3, "{mode:?}: reference wall reflection only {pec}");
        // 10-cell CPML measures ≈ −88 dB here; require ≤ −60 dB.
        assert!(classic < 1e-3, "{mode:?}: CPML reflection {classic}");
        assert!(cfs_r < 1e-3, "{mode:?}: CFS-CPML reflection {cfs_r}");
    }

    #[test]
    fn cpml_reflection_tm() {
        check_reflection(Mode2d::Tm);
    }

    #[test]
    fn cpml_reflection_te() {
        check_reflection(Mode2d::Te);
    }

    /// Peak and final total energy for a zero-DC (Ricker) pulse.
    fn energy_history(mode: Mode2d, cpml: bool) -> (f64, f64) {
        let mut sim = Fdtd2d::new(N, N, DX, mode).unwrap();
        if cpml {
            sim.try_enable_cpml(T).unwrap();
        }
        let tau = 25.0 * sim.dt;
        let t0 = 4.0 * tau;
        let mut peak = 0.0_f64;
        for step in 0..1500 {
            sim.add_source(N / 2, N / 2, tau * ricker_wavelet(sim.time(), t0, tau));
            sim.step_once();
            if step % 5 == 0 {
                peak = peak.max(sim.total_energy());
            }
        }
        (peak, sim.total_energy())
    }

    #[test]
    fn cpml_energy_decays() {
        for mode in [Mode2d::Tm, Mode2d::Te] {
            let (peak_wall, end_wall) = energy_history(mode, false);
            let (peak_cpml, end_cpml) = energy_history(mode, true);
            eprintln!(
                "{mode:?}: wall end/peak {:.3e}, CPML end/peak {:.3e}",
                end_wall / peak_wall,
                end_cpml / peak_cpml
            );
            assert!(
                end_wall / peak_wall > 0.5,
                "{mode:?}: walls must trap energy"
            );
            assert!(
                end_cpml / peak_cpml < 1e-4,
                "{mode:?}: CPML must absorb energy"
            );
        }
    }

    #[test]
    fn cpml_stable_long_run_with_materials() {
        for mode in [Mode2d::Tm, Mode2d::Te] {
            let mut sim = Fdtd2d::new(60, 60, DX, mode).unwrap();
            sim.try_enable_cpml(8).unwrap();
            // Dielectric + lossy block that overlaps the PML strip.
            sim.set_material(0, 20, 30, 40, &Material::new(4.0, 1.0, 0.5))
                .unwrap();
            sim.run(3000, Some((40, 30)), 3e10, 1.0);
            assert!(
                sim.fz
                    .iter()
                    .chain(&sim.fx)
                    .chain(&sim.fy)
                    .all(|v| v.is_finite())
            );
        }
    }

    #[test]
    fn cpml_profile_is_graded() {
        let sim = {
            let mut s = Fdtd2d::new(60, 60, DX, Mode2d::Tm).unwrap();
            s.try_enable_cpml(10).unwrap();
            s
        };
        let cpml = sim.cpml().unwrap();
        // σ_opt = 0.8 (m + 1) / (η₀ Δx) for m = 3.
        let sigma_opt = 0.8 * 4.0 / (ETA_0 * DX);
        assert!((cpml.sigma_max() - sigma_opt).abs() < 1e-9 * sigma_opt);
        assert!((cpml.grading_order() - 3.0).abs() < 1e-15);
        assert!((cpml.kappa_max() - 1.0).abs() < 1e-15);
        assert!(cpml.alpha_max().abs() < 1e-15);
        // b increases (less loss) and c → 0 towards the interface.
        for w in cpml.b_node.windows(2).chain(cpml.b_half.windows(2)) {
            assert!(w[0] < w[1]);
        }
        for &b in cpml.b_node.iter().chain(&cpml.b_half) {
            assert!(b > 0.0 && b < 1.0);
        }
        for &c in cpml.c_node.iter().chain(&cpml.c_half) {
            assert!(c < 0.0 && c > -1.0);
        }
        // Half nodes lie between their integer neighbours.
        for l in 0..9 {
            assert!(cpml.b_node[l] < cpml.b_half[l] && cpml.b_half[l] < cpml.b_node[l + 1]);
        }
    }

    #[test]
    fn cpml_kappa_profile() {
        let mut p = CpmlParams::new(6);
        p.kappa_max = 7.0;
        let cpml = Cpml::with_params(&p, DX, 1e-13).unwrap();
        // At the wall: 1/κ − 1 = 1/7 − 1.
        assert!((cpml.km1_node[0] - (1.0 / 7.0 - 1.0)).abs() < 1e-12);
        assert!(cpml.km1_half.iter().all(|&k| k < 0.0 && k > -1.0));
    }

    #[test]
    fn cpml_params_validation() {
        let ok = CpmlParams::new(8);
        assert!(Cpml::with_params(&ok, DX, 1e-13).is_ok());
        assert!(Cpml::with_params(&ok, f64::NAN, 1e-13).is_err());
        assert!(Cpml::with_params(&ok, DX, 0.0).is_err());
        assert!(Cpml::with_params(&CpmlParams::new(0), DX, 1e-13).is_err());
        let mut bad = ok;
        bad.kappa_max = 0.5;
        assert!(Cpml::with_params(&bad, DX, 1e-13).is_err());
        let mut bad = ok;
        bad.alpha_max = -1.0;
        assert!(Cpml::with_params(&bad, DX, 1e-13).is_err());
        let mut bad = ok;
        bad.sigma_ratio = f64::NAN;
        assert!(Cpml::with_params(&bad, DX, 1e-13).is_err());
        let mut bad = ok;
        bad.grading_order = 0.0;
        assert!(Cpml::with_params(&bad, DX, 1e-13).is_err());
        // Infallible constructor degrades to an inert, unattachable layer.
        let mut inert = Cpml::new(8, -1.0, 1e-13);
        inert.attach(60, 60);
        assert!(!inert.is_attached_to(60, 60));
    }

    #[test]
    fn try_enable_cpml_validates_thickness() {
        let mut sim = Fdtd2d::new(30, 40, DX, Mode2d::Te).unwrap();
        assert!(sim.try_enable_cpml(0).is_err());
        // 2·14 + 3 = 31 > min(30, 40)
        assert!(sim.try_enable_cpml(14).is_err());
        assert!(sim.cpml().is_none());
        sim.try_enable_cpml(13).unwrap();
        assert_eq!(sim.cpml().map(|c| c.thickness), Some(13));
        // Invalid request keeps the active layer.
        assert!(sim.try_enable_cpml(100).is_err());
        sim.enable_cpml(100);
        assert_eq!(sim.cpml().map(|c| c.thickness), Some(13));
        sim.disable_cpml();
        assert!(sim.cpml().is_none());
        sim.enable_cpml(0);
        assert!(sim.cpml().is_none());
    }

    // ── Construction ──────────────────────────────────────────────

    #[test]
    fn new_rejects_non_finite_dx() {
        for dx in [f64::NAN, f64::INFINITY, -1e-3, 0.0] {
            assert!(Fdtd2d::new(10, 10, dx, Mode2d::Tm).is_err(), "2D dx = {dx}");
            assert!(Fdtd3d::new(5, 5, 5, dx).is_err(), "3D dx = {dx}");
        }
        assert!(Fdtd2d::new(usize::MAX, 3, 1e-3, Mode2d::Tm).is_err());
        assert!(Fdtd3d::new(usize::MAX, usize::MAX, 3, 1e-3).is_err());
    }

    #[test]
    fn courant_numbers_are_half_the_limit() {
        let s2 = Fdtd2d::new(10, 10, 1e-3, Mode2d::Tm).unwrap();
        let courant2 = SPEED_OF_LIGHT * s2.dt / s2.dx;
        assert!((courant2 - 0.5 / std::f64::consts::SQRT_2).abs() < 1e-12);
        let s3 = Fdtd3d::new(5, 5, 5, 1e-3).unwrap();
        let courant3 = SPEED_OF_LIGHT * s3.dt / s3.dx;
        assert!((courant3 - 0.5 / 3.0_f64.sqrt()).abs() < 1e-12);
    }

    // ── Materials ─────────────────────────────────────────────────

    #[test]
    fn set_material_validates_before_mutating() {
        let mut sim = Fdtd2d::new(20, 20, 1e-3, Mode2d::Tm).unwrap();
        let (eps, ce, ch) = (sim.permittivity.clone(), sim.ce.clone(), sim.ch.clone());
        let bad = [
            Material::new(4.0, 0.0, 0.0),
            Material::new(4.0, f64::NAN, 0.0),
            Material::new(f64::NAN, 1.0, 0.0),
            Material::new(4.0, 1.0, -1.0),
            Material::new(4.0, 1.0, f64::INFINITY),
        ];
        for mat in &bad {
            assert!(sim.set_material(0, 0, 20, 20, mat).is_err(), "{mat:?}");
        }
        assert!(
            sim.set_material(25, 0, 30, 20, &Material::dielectric(2.0))
                .is_err()
        );
        assert_eq!(sim.permittivity, eps);
        assert_eq!(sim.ce, ce);
        assert_eq!(sim.ch, ch);
        assert!(sim.ca.is_empty());

        let mut s3 = Fdtd3d::new(6, 6, 6, 1e-3).unwrap();
        let (eps3, ch3) = (s3.permittivity.clone(), s3.ch.clone());
        for mat in &bad {
            assert!(s3.set_material(0, 0, 0, 6, 6, 6, mat).is_err(), "{mat:?}");
        }
        assert_eq!(s3.permittivity, eps3);
        assert_eq!(s3.ch, ch3);
    }

    #[test]
    fn region_bounds() {
        let mut sim = Fdtd2d::new(10, 10, 1e-3, Mode2d::Tm).unwrap();
        // Start outside the grid or inverted → error.
        assert!(sim.set_permittivity(10, 0, 12, 5, 2.0).is_err());
        assert!(sim.set_permittivity(0, 10, 5, 12, 2.0).is_err());
        assert!(sim.set_permittivity(5, 0, 4, 5, 2.0).is_err());
        assert!(sim.set_conductivity(0, 0, 10, 10, -1.0).is_err());
        // End coordinates are clamped.
        sim.set_permittivity(8, 8, 100, 100, 2.0).unwrap();
        assert!((sim.permittivity[sim.idx(9, 9)] - 2.0).abs() < 1e-15);
        // Empty region is a no-op.
        sim.set_permittivity(3, 3, 3, 3, 5.0).unwrap();
        assert!(sim.permittivity.iter().all(|&e| e < 3.0));

        let mut s3 = Fdtd3d::new(5, 5, 5, 1e-3).unwrap();
        assert!(s3.set_permittivity(0, 0, 5, 5, 5, 6, 2.0).is_err());
        assert!(s3.set_permeability(0, 0, 0, 5, 5, 5, f64::NAN).is_err());
        s3.set_permeability(1, 1, 1, 9, 9, 9, 3.0).unwrap();
        assert!((s3.permeability[s3.idx(4, 4, 4)] - 3.0).abs() < 1e-15);
        assert!((s3.permeability[s3.idx(0, 0, 0)] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn conductivity_coefficients() {
        let mut sim = Fdtd2d::new(10, 10, 1e-3, Mode2d::Tm).unwrap();
        let lossless_ce = sim.ce[0];
        sim.set_material(2, 2, 5, 5, &Material::new(2.0, 1.0, 0.3))
            .unwrap();
        let (ca, cb) = LossyMaterial::new(2.0, 0.3).update_coefficients(sim.dt, sim.dx);
        let i = sim.idx(3, 3);
        assert!((sim.ca[i] - ca).abs() < 1e-15);
        assert!((sim.ce[i] - cb).abs() < 1e-12 * cb);
        assert!((sim.ca[0] - 1.0).abs() < 1e-15);
        assert!((sim.ce[0] - lossless_ce).abs() < 1e-30);
        // Changing ε keeps σ; setting σ = 0 restores the lossless coefficients.
        sim.set_permittivity(2, 2, 5, 5, 4.0).unwrap();
        let (ca4, _) = LossyMaterial::new(4.0, 0.3).update_coefficients(sim.dt, sim.dx);
        assert!((sim.ca[i] - ca4).abs() < 1e-15);
        sim.set_conductivity(2, 2, 5, 5, 0.0).unwrap();
        assert!((sim.ca[i] - 1.0).abs() < 1e-15);
        assert!((sim.ce[i] - lossless_ce / 4.0).abs() < 1e-12 * lossless_ce);
    }

    /// A spatially uniform field in a conductor obeys dE/dt = −σE/ε
    /// (no curl), so E(t) = E₀ exp(−σt/ε). Probe far enough from the walls
    /// that their influence (≤ 1 cell per step) cannot reach it.
    #[test]
    fn conductor_uniform_decay_matches_analytic_2d() {
        for mode in [Mode2d::Tm, Mode2d::Te] {
            let n = 41;
            let mut sim = Fdtd2d::new(n, n, 1e-3, mode).unwrap();
            let eps_r = 2.0;
            let sigma = 0.05 * eps_r * EPSILON_0 / sim.dt; // σΔt/ε = 0.05
            sim.set_material(0, 0, n, n, &Material::new(eps_r, 1.0, sigma))
                .unwrap();
            match mode {
                Mode2d::Tm => sim.fz.fill(1.0),
                Mode2d::Te => {
                    sim.fx.fill(1.0);
                    sim.fy.fill(1.0);
                }
            }
            let steps = 19;
            for _ in 0..steps {
                sim.step_once();
            }
            let c = sim.idx(n / 2, n / 2);
            let e = if mode == Mode2d::Tm {
                sim.fz[c]
            } else {
                sim.fx[c]
            };
            let analytic = (-sigma * steps as f64 * sim.dt / (eps_r * EPSILON_0)).exp();
            let (ca, _) = LossyMaterial::new(eps_r, sigma).update_coefficients(sim.dt, sim.dx);
            assert!((e - ca.powi(steps)).abs() < 1e-12, "{mode:?}: discrete");
            assert!(
                (e - analytic).abs() < 1e-3 * analytic,
                "{mode:?}: {e} vs {analytic}"
            );
        }
    }

    #[test]
    fn conductor_uniform_decay_matches_analytic_3d() {
        let n = 21;
        let mut sim = Fdtd3d::new(n, n, n, 1e-3).unwrap();
        let sigma = 0.1 * EPSILON_0 / sim.dt; // σΔt/ε₀ = 0.1
        sim.set_conductivity(0, 0, 0, n, n, n, sigma).unwrap();
        sim.ex.fill(1.0);
        sim.ez.fill(1.0);
        let steps = 9;
        for _ in 0..steps {
            sim.step_once();
        }
        let c = sim.idx(n / 2, n / 2, n / 2);
        let analytic = (-sigma * steps as f64 * sim.dt / EPSILON_0).exp();
        for e in [sim.ex[c], sim.ez[c]] {
            assert!((e - analytic).abs() < 2e-3 * analytic, "{e} vs {analytic}");
        }
    }

    #[test]
    fn conductor_slab_attenuates_wave() {
        // A good conductor between source and probe kills the transmitted field.
        let run = |sigma: f64| {
            let mut sim = Fdtd2d::new(80, 60, DX, Mode2d::Tm).unwrap();
            sim.try_enable_cpml(10).unwrap();
            sim.set_conductivity(40, 0, 46, 60, sigma).unwrap();
            let tau = 25.0 * sim.dt;
            let mut peak = 0.0_f64;
            for _ in 0..600 {
                sim.add_source(25, 30, gaussian_pulse(sim.time(), 4.0 * tau, tau));
                sim.step_once();
                peak = peak.max(sim.fz[sim.idx(60, 30)].abs());
            }
            peak
        };
        let open = run(0.0);
        let blocked = run(1e3);
        assert!(blocked < 1e-2 * open, "open {open}, blocked {blocked}");
    }

    // ── Robustness against desynchronised public fields ───────────

    #[test]
    fn desynced_fields_do_not_panic_2d() {
        let mut sim = Fdtd2d::new(20, 30, 1e-3, Mode2d::Tm).unwrap();
        sim.try_enable_cpml(4).unwrap();
        sim.nx = 1000;
        sim.step_once();
        sim.run(3, Some((500, 5)), 1e9, 1.0);
        sim.add_source(500, 5, 1.0);
        let _ = sim.total_energy();
        assert!(sim.set_permittivity(0, 0, 5, 5, 2.0).is_err());
        assert!(sim.try_enable_cpml(4).is_err());
        assert_eq!(sim.step, 0);
        // Same cell count, but the CPML no longer matches the geometry.
        sim.nx = 30;
        sim.ny = 20;
        sim.step_once();
        assert_eq!(sim.step, 0);
        sim.nx = 20;
        sim.ny = 30;
        sim.step_once();
        assert_eq!(sim.step, 1);
        sim.fz.truncate(10);
        sim.step_once();
        sim.add_source(15, 15, 1.0);
        let _ = sim.total_energy();
        assert_eq!(sim.step, 1);
        sim.nx = 0;
        sim.ny = 0;
        sim.step_once();
        assert_eq!(sim.step, 1);
    }

    #[test]
    fn desynced_fields_do_not_panic_3d() {
        let mut sim = Fdtd3d::new(5, 5, 5, 1e-3).unwrap();
        sim.nz = 0;
        sim.step_once();
        sim.run(2, Some((1, 1, 1)), 1e9, 1.0);
        sim.nz = 5;
        sim.ez.clear();
        sim.step_once();
        sim.add_source_ez(2, 2, 2, 1.0);
        let _ = sim.total_energy();
        assert!(sim.set_conductivity(0, 0, 0, 5, 5, 5, 1.0).is_err());
        sim.nx = usize::MAX;
        sim.add_source_ez(3, 3, 3, 1.0);
        sim.step_once();
        assert_eq!(sim.step, 0);
    }

    #[test]
    fn lossless_path_unchanged_by_material_api() {
        // A lossless set_material keeps the cheaper lossless update.
        let mut sim = Fdtd2d::new(10, 10, 1e-3, Mode2d::Te).unwrap();
        sim.set_material(0, 0, 10, 10, &Material::dielectric(2.0))
            .unwrap();
        assert!(sim.ca.is_empty() && sim.conductivity.is_empty());
        let mut s3 = Fdtd3d::new(5, 5, 5, 1e-3).unwrap();
        s3.set_material(0, 0, 0, 5, 5, 5, &Material::dielectric(2.0))
            .unwrap();
        assert!(s3.ca.is_empty() && s3.conductivity.is_empty());
    }
}

// ── 3D FDTD ──────────────────────────────────────────────────────

/// Geometry of one 3D CPML correction pass (see [`Cpml::correct_slabs`]).
#[derive(Debug, Clone)]
struct SlabPass {
    /// Grid nodes along x, y, z.
    dims: [usize; 3],
    /// Main-loop update ranges of the corrected field along x, y, z; the
    /// range along the derivative axis is replaced by the two slabs.
    ranges: [std::ops::Range<usize>; 3],
    /// Sign of this derivative term in the curl update.
    sign: f64,
}

/// A [`Cpml`] attached to an [`Fdtd3d`] grid.
///
/// Holds the 12 auxiliary Ψ fields of the Roden–Gedney CPML, named
/// `<component>_<derivative axis>` (e.g. `hx_z` is Ψ of ∂E_y/∂z in the H_x
/// update). Each is stored only in the two face slabs normal to its
/// derivative axis: the grid's shape with that axis cut to 2·thickness
/// layers. The slabs overlap in the edge and corner regions, where several
/// Ψ terms apply to the same component.
#[derive(Debug, Clone)]
struct Cpml3d {
    /// Grading profiles and parameters.
    layer: Cpml,
    /// Grid nodes along x, y, z the Ψ fields are allocated for.
    dims: [usize; 3],
    hx_y: Vec<f64>,
    hx_z: Vec<f64>,
    hy_x: Vec<f64>,
    hy_z: Vec<f64>,
    hz_x: Vec<f64>,
    hz_y: Vec<f64>,
    ex_y: Vec<f64>,
    ex_z: Vec<f64>,
    ey_x: Vec<f64>,
    ey_z: Vec<f64>,
    ez_x: Vec<f64>,
    ez_y: Vec<f64>,
}

impl Cpml3d {
    /// Ψ lengths of the slabs normal to x, y and z (`None` on overflow).
    fn slab_lens(thickness: usize, dims: [usize; 3]) -> Option<[usize; 3]> {
        let [nx, ny, nz] = dims;
        let t2 = thickness.checked_mul(2)?;
        let len = |a: usize, b: usize| t2.checked_mul(a)?.checked_mul(b);
        Some([len(ny, nz)?, len(nx, nz)?, len(nx, ny)?])
    }

    /// Allocate zeroed Ψ slabs for `layer` on a grid of `dims` nodes.
    fn attach(layer: Cpml, dims: [usize; 3]) -> Result<Self> {
        let [sx, sy, sz] =
            Self::slab_lens(layer.thickness, dims).ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("CPML slabs for a {dims:?} grid are too large"),
            })?;
        Ok(Self {
            layer,
            dims,
            hx_y: try_filled(sy, 0.0)?,
            hx_z: try_filled(sz, 0.0)?,
            hy_x: try_filled(sx, 0.0)?,
            hy_z: try_filled(sz, 0.0)?,
            hz_x: try_filled(sx, 0.0)?,
            hz_y: try_filled(sy, 0.0)?,
            ex_y: try_filled(sy, 0.0)?,
            ex_z: try_filled(sz, 0.0)?,
            ey_x: try_filled(sx, 0.0)?,
            ey_z: try_filled(sz, 0.0)?,
            ez_x: try_filled(sx, 0.0)?,
            ez_y: try_filled(sy, 0.0)?,
        })
    }

    /// Whether the profiles and every Ψ slab match a grid of `dims` nodes;
    /// the slab corrections rely on this (it bounds every index they use).
    fn is_attached_to(&self, dims: [usize; 3]) -> bool {
        let t = self.layer.thickness;
        let min_dim = dims[0].min(dims[1]).min(dims[2]);
        if !cpml_fits(t, min_dim) || self.dims != dims || !self.layer.has_profiles() {
            return false;
        }
        let Some([sx, sy, sz]) = Self::slab_lens(t, dims) else {
            return false;
        };
        [
            (&self.hx_y, sy),
            (&self.hx_z, sz),
            (&self.hy_x, sx),
            (&self.hy_z, sz),
            (&self.hz_x, sx),
            (&self.hz_y, sy),
            (&self.ex_y, sy),
            (&self.ex_z, sz),
            (&self.ey_x, sx),
            (&self.ey_z, sz),
            (&self.ez_x, sx),
            (&self.ez_y, sy),
        ]
        .iter()
        .all(|(v, len)| v.len() == *len)
    }

    /// CPML terms of the H update (half-cell profiles), applied after the
    /// main H loops with the same update ranges.
    ///
    /// Caller guarantees `is_attached_to(dims)` and slice lengths `nx·ny·nz`.
    fn correct_h(&mut self, h: [&mut [f64]; 3], e: [&[f64]; 3], ch: &[f64]) {
        let Self {
            layer,
            dims,
            hx_y,
            hx_z,
            hy_x,
            hy_z,
            hz_x,
            hz_y,
            ..
        } = self;
        let dims = *dims;
        let [nx, ny, nz] = dims;
        let [hx, hy, hz] = h;
        let [ex, ey, ez] = e;
        let pass = |ranges, sign| SlabPass { dims, ranges, sign };

        // H_x at (x, y+½, z+½) += C_h·(Δ_z E_y − Δ_y E_z)
        let r = [0..nx, 0..ny - 1, 0..nz - 1];
        layer.correct_slabs::<2, true>(hx_z, hx, ey, ch, &pass(r.clone(), 1.0));
        layer.correct_slabs::<1, true>(hx_y, hx, ez, ch, &pass(r, -1.0));
        // H_y at (x+½, y, z+½) += C_h·(Δ_x E_z − Δ_z E_x)
        let r = [0..nx - 1, 0..ny, 0..nz - 1];
        layer.correct_slabs::<0, true>(hy_x, hy, ez, ch, &pass(r.clone(), 1.0));
        layer.correct_slabs::<2, true>(hy_z, hy, ex, ch, &pass(r, -1.0));
        // H_z at (x+½, y+½, z) += C_h·(Δ_y E_x − Δ_x E_y)
        let r = [0..nx - 1, 0..ny - 1, 0..nz];
        layer.correct_slabs::<1, true>(hz_y, hz, ex, ch, &pass(r.clone(), 1.0));
        layer.correct_slabs::<0, true>(hz_x, hz, ey, ch, &pass(r, -1.0));
    }

    /// CPML terms of the E update (integer-node profiles), applied after the
    /// main E loops and before the PEC walls.
    ///
    /// Caller guarantees `is_attached_to(dims)` and slice lengths `nx·ny·nz`.
    fn correct_e(&mut self, e: [&mut [f64]; 3], h: [&[f64]; 3], ce: &[f64]) {
        let Self {
            layer,
            dims,
            ex_y,
            ex_z,
            ey_x,
            ey_z,
            ez_x,
            ez_y,
            ..
        } = self;
        let dims = *dims;
        let [nx, ny, nz] = dims;
        let [ex, ey, ez] = e;
        let [hx, hy, hz] = h;
        let pass = |ranges, sign| SlabPass { dims, ranges, sign };

        // E_x at (x+½, y, z) += C_b·(Δ_y H_z − Δ_z H_y)
        let r = [0..nx, 1..ny - 1, 1..nz - 1];
        layer.correct_slabs::<1, false>(ex_y, ex, hz, ce, &pass(r.clone(), 1.0));
        layer.correct_slabs::<2, false>(ex_z, ex, hy, ce, &pass(r, -1.0));
        // E_y at (x, y+½, z) += C_b·(Δ_z H_x − Δ_x H_z)
        let r = [1..nx - 1, 0..ny, 1..nz - 1];
        layer.correct_slabs::<2, false>(ey_z, ey, hx, ce, &pass(r.clone(), 1.0));
        layer.correct_slabs::<0, false>(ey_x, ey, hz, ce, &pass(r, -1.0));
        // E_z at (x, y, z+½) += C_b·(Δ_x H_y − Δ_y H_x)
        let r = [1..nx - 1, 1..ny - 1, 0..nz];
        layer.correct_slabs::<0, false>(ez_x, ez, hy, ce, &pass(r.clone(), 1.0));
        layer.correct_slabs::<1, false>(ez_y, ez, hx, ce, &pass(r, -1.0));
    }
}

/// 3D FDTD simulation on a full vector Yee grid.
///
/// Updates all six field components: E_x, E_y, E_z, H_x, H_y, H_z.
/// Fields are stored in row-major order: index = z * (nx * ny) + y * nx + x.
///
/// # Boundaries
///
/// By default the outer faces are perfect electric conductors (PEC):
/// tangential E is held at zero on all six faces, so outgoing waves are
/// fully reflected. Call [`Fdtd3d::try_enable_cpml`] to absorb outgoing waves
/// instead: the CPML occupies the outermost `thickness` cells on every face
/// (including the edge and corner regions) and is terminated by the same
/// PEC wall.
///
/// # Materials
///
/// Per-cell ε_r, μ_r and electric conductivity σ (lossy E update with
/// C_a/C_b coefficients) through [`Fdtd3d::set_material`],
/// [`Fdtd3d::set_permittivity`], [`Fdtd3d::set_permeability`] and
/// [`Fdtd3d::set_conductivity`].
///
/// # Public fields
///
/// `nx`, `ny`, `nz`, `dx`, `dt` and the field vectors are public for reading
/// and for writing field *values*. Treat the geometry and the vector lengths
/// as read-only: if they are changed so the layout is inconsistent, stepping
/// logs an error and leaves the state untouched instead of panicking;
/// changing `dx` or `dt` alone silently invalidates the precomputed
/// coefficients.
#[derive(Debug, Clone)]
pub struct Fdtd3d {
    /// Number of nodes along x (read-only in practice — see the type docs).
    pub nx: usize,
    /// Number of nodes along y (read-only in practice).
    pub ny: usize,
    /// Number of nodes along z (read-only in practice).
    pub nz: usize,
    /// Uniform spatial step (m). Read-only in practice.
    pub dx: f64,
    /// Time step (s). Read-only in practice.
    pub dt: f64,
    /// Current step index.
    pub step: usize,

    /// Electric field E_x.
    pub ex: Vec<f64>,
    /// Electric field E_y.
    pub ey: Vec<f64>,
    /// Electric field E_z.
    pub ez: Vec<f64>,
    /// Magnetic field H_x.
    pub hx: Vec<f64>,
    /// Magnetic field H_y.
    pub hy: Vec<f64>,
    /// Magnetic field H_z.
    pub hz: Vec<f64>,

    /// E-field curl coefficient C_b per cell: Δt/(ε₀ε_r Δx), divided by
    /// (1 + σΔt/2ε) in lossy cells.
    ce: Vec<f64>,
    /// H-field update coefficient per cell: dt/(μ₀μ_r dx).
    ch: Vec<f64>,
    /// E-field decay coefficient C_a per cell; empty while all cells are lossless.
    ca: Vec<f64>,
    /// Relative permittivity per cell.
    permittivity: Vec<f64>,
    /// Relative permeability per cell.
    permeability: Vec<f64>,
    /// Electric conductivity σ per cell (S/m); empty while all cells are lossless.
    conductivity: Vec<f64>,
    /// Optional CPML absorbing boundaries.
    cpml: Option<Cpml3d>,
}

/// Index ranges of a validated 3D region.
type Region3d = (
    std::ops::Range<usize>,
    std::ops::Range<usize>,
    std::ops::Range<usize>,
);

impl Fdtd3d {
    /// Create a new 3D FDTD simulation (vacuum, PEC walls).
    ///
    /// Uses uniform grid spacing `dx` in all three dimensions. The time step
    /// is fixed at dt = Δx / (2√3 c), i.e. a Courant number
    /// S = c·dt/Δx = 1/(2√3) ≈ 0.289, half of the 3D stability limit 1/√3.
    ///
    /// # Errors
    ///
    /// - [`BijliError::InsufficientResolution`] if any dimension is < 3.
    /// - [`BijliError::InvalidParameter`] if `dx` is not finite and positive,
    ///   or the grid is too large to allocate.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Result<Self> {
        if nx < 3 || ny < 3 || nz < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx.min(ny).min(nz),
                domain_size: nx.min(ny).min(nz) as f64 * dx,
            });
        }
        validate_dx(dx)?;
        let n = grid_cells(&[nx, ny, nz])?;

        // 3D Courant: dt ≤ dx/(c√3). Use half for safety (S = 1/(2√3)).
        let dt = dx / (2.0 * 3.0_f64.sqrt() * SPEED_OF_LIGHT);

        let ce = try_filled(n, dt / (EPSILON_0 * dx))?;
        let ch = try_filled(n, dt / (MU_0 * dx))?;
        let ex = try_filled(n, 0.0)?;
        let ey = try_filled(n, 0.0)?;
        let ez = try_filled(n, 0.0)?;
        let hx = try_filled(n, 0.0)?;
        let hy = try_filled(n, 0.0)?;
        let hz = try_filled(n, 0.0)?;
        let permittivity = try_filled(n, 1.0)?;
        let permeability = try_filled(n, 1.0)?;

        tracing::debug!(nx, ny, nz, dx, dt, n, "3D FDTD simulation created");

        Ok(Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            step: 0,
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
            ce,
            ch,
            ca: Vec::new(),
            permittivity,
            permeability,
            conductivity: Vec::new(),
            cpml: None,
        })
    }

    /// Linear index from (x, y, z) coordinates.
    #[inline]
    #[must_use]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * (self.nx * self.ny) + y * self.nx + x
    }

    /// [`Self::check_grid`] plus consistency of the attached CPML.
    fn check_layout(&self) -> Result<usize> {
        let n = self.check_grid()?;
        let dims = [self.nx, self.ny, self.nz];
        if self.cpml.as_ref().is_some_and(|c| !c.is_attached_to(dims)) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "CPML is not attached to the current {}×{}×{} grid",
                    self.nx, self.ny, self.nz
                ),
            });
        }
        Ok(n)
    }

    /// Check that the grid geometry and every per-cell vector agree.
    fn check_grid(&self) -> Result<usize> {
        let n = grid_cells(&[self.nx, self.ny, self.nz])
            .ok()
            .filter(|_| self.nx >= 3 && self.ny >= 3 && self.nz >= 3);
        let consistent = n.is_some_and(|n| {
            [
                &self.ex,
                &self.ey,
                &self.ez,
                &self.hx,
                &self.hy,
                &self.hz,
                &self.ce,
                &self.ch,
                &self.permittivity,
                &self.permeability,
            ]
            .iter()
            .all(|v| v.len() == n)
                && ((self.ca.is_empty() && self.conductivity.is_empty())
                    || (self.ca.len() == n && self.conductivity.len() == n))
        });
        match n {
            Some(n) if consistent => Ok(n),
            _ => Err(BijliError::InvalidParameter {
                reason: format!(
                    "inconsistent 3D grid state: {}×{}×{} grid, E lengths = ({}, {}, {})",
                    self.nx,
                    self.ny,
                    self.nz,
                    self.ex.len(),
                    self.ey.len(),
                    self.ez.len()
                ),
            }),
        }
    }

    /// Validate a box region `[x0, x1) × [y0, y1) × [z0, z1)`.
    #[allow(clippy::too_many_arguments)]
    fn region(
        &self,
        x0: usize,
        y0: usize,
        z0: usize,
        x1: usize,
        y1: usize,
        z1: usize,
    ) -> Result<Region3d> {
        Ok((
            region_range("x", x0, x1, self.nx)?,
            region_range("y", y0, y1, self.ny)?,
            region_range("z", z0, z1, self.nz)?,
        ))
    }

    /// Write material values to an already validated region and refresh the
    /// update coefficients. `None` leaves that property unchanged.
    fn apply_region(
        &mut self,
        region: Region3d,
        eps_r: Option<f64>,
        mu_r: Option<f64>,
        sigma: Option<f64>,
    ) -> Result<()> {
        let n = self.check_grid()?;
        if sigma.is_some_and(|s| s > 0.0) && self.ca.is_empty() {
            let ca = try_filled(n, 1.0)?;
            let conductivity = try_filled(n, 0.0)?;
            self.ca = ca;
            self.conductivity = conductivity;
        }
        tracing::debug!(
            ?region,
            ?eps_r,
            ?mu_r,
            ?sigma,
            "3D FDTD material region set"
        );
        let (xs, ys, zs) = region;
        let (dt, dx, nx, nxy) = (self.dt, self.dx, self.nx, self.nx * self.ny);
        let base_h = dt / (MU_0 * dx);
        for z in zs {
            for y in ys.clone() {
                for x in xs.clone() {
                    let i = z * nxy + y * nx + x;
                    if let Some(mu) = mu_r {
                        self.permeability[i] = mu;
                        self.ch[i] = base_h / mu;
                    }
                    if eps_r.is_none() && sigma.is_none() {
                        continue;
                    }
                    if let Some(eps) = eps_r {
                        self.permittivity[i] = eps;
                    }
                    if let (Some(s), Some(c)) = (sigma, self.conductivity.get_mut(i)) {
                        *c = s;
                    }
                    let s = self.conductivity.get(i).copied().unwrap_or(0.0);
                    let (ca, cb) = e_update_coefficients(dt, dx, self.permittivity[i], s);
                    self.ce[i] = cb;
                    if let Some(a) = self.ca.get_mut(i) {
                        *a = ca;
                    }
                }
            }
        }
        Ok(())
    }

    /// Set relative permittivity in the box `[x0, x1) × [y0, y1) × [z0, z1)`.
    ///
    /// The exclusive end coordinates are clamped to the grid. Any
    /// conductivity already set in the region is kept.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidPermittivity`] if `eps_r` is not finite and
    /// positive; [`BijliError::InvalidParameter`] if the region starts
    /// outside the grid or is inverted. Nothing is modified on error.
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
        validate_eps_r(eps_r)?;
        let region = self.region(x0, y0, z0, x1, y1, z1)?;
        self.apply_region(region, Some(eps_r), None, None)
    }

    /// Set relative permeability in the box `[x0, x1) × [y0, y1) × [z0, z1)`.
    ///
    /// The exclusive end coordinates are clamped to the grid.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidPermeability`] if `mu_r` is not finite and
    /// positive; [`BijliError::InvalidParameter`] if the region starts
    /// outside the grid or is inverted. Nothing is modified on error.
    #[allow(clippy::too_many_arguments)]
    pub fn set_permeability(
        &mut self,
        x0: usize,
        y0: usize,
        z0: usize,
        x1: usize,
        y1: usize,
        z1: usize,
        mu_r: f64,
    ) -> Result<()> {
        validate_mu_r(mu_r)?;
        let region = self.region(x0, y0, z0, x1, y1, z1)?;
        self.apply_region(region, None, Some(mu_r), None)
    }

    /// Set electric conductivity σ (S/m) in the box
    /// `[x0, x1) × [y0, y1) × [z0, z1)`.
    ///
    /// Lossy cells use E ← C_a E + C_b ∇×H with
    /// C_a = (1 − σΔt/2ε)/(1 + σΔt/2ε), C_b = (Δt/εΔx)/(1 + σΔt/2ε).
    /// The exclusive end coordinates are clamped to the grid.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `sigma` is negative or not finite,
    /// or the region starts outside the grid or is inverted. Nothing is
    /// modified on error.
    #[allow(clippy::too_many_arguments)]
    pub fn set_conductivity(
        &mut self,
        x0: usize,
        y0: usize,
        z0: usize,
        x1: usize,
        y1: usize,
        z1: usize,
        sigma: f64,
    ) -> Result<()> {
        validate_conductivity(sigma)?;
        let region = self.region(x0, y0, z0, x1, y1, z1)?;
        self.apply_region(region, None, None, Some(sigma))
    }

    /// Set material properties in a rectangular region using a [`crate::material::Material`] struct.
    ///
    /// Sets permittivity, permeability and conductivity. All inputs are
    /// validated before anything is modified.
    ///
    /// # Errors
    ///
    /// Same as [`Fdtd3d::set_permittivity`], [`Fdtd3d::set_permeability`] and
    /// [`Fdtd3d::set_conductivity`].
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
        validate_eps_r(mat.eps_r)?;
        validate_mu_r(mat.mu_r)?;
        validate_conductivity(mat.conductivity)?;
        let region = self.region(x0, y0, z0, x1, y1, z1)?;
        self.apply_region(
            region,
            Some(mat.eps_r),
            Some(mat.mu_r),
            Some(mat.conductivity),
        )
    }

    /// Inject a soft source at grid point (x, y, z) into E_z.
    ///
    /// Points outside the grid are ignored (no panic). Values injected where
    /// E_z is tangential to a wall (x or y faces) are cleared on the next step.
    #[inline]
    pub fn add_source_ez(&mut self, x: usize, y: usize, z: usize, value: f64) {
        if x < self.nx
            && y < self.ny
            && z < self.nz
            && let Some(v) = self
                .nx
                .checked_mul(self.ny)
                .and_then(|nxy| nxy.checked_mul(z))
                .and_then(|plane| plane.checked_add(y.checked_mul(self.nx)?))
                .and_then(|row| row.checked_add(x))
                .and_then(|i| self.ez.get_mut(i))
        {
            *v += value;
        }
    }

    /// Enable CPML absorbing boundaries on all six faces with the classic
    /// grading ([`CpmlParams::new`]) and the given thickness (cells per side).
    ///
    /// Replaces the reflecting PEC walls by a convolutional PML occupying the
    /// outermost `thickness` cells on every face. Typical thickness: 8–20
    /// cells. Calling it again replaces the layer (Ψ restarts at zero).
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `thickness == 0` or
    /// `2·thickness + 3 > min(nx, ny, nz)`, or the grid state is
    /// inconsistent. The current boundary configuration is unchanged on error.
    ///
    /// # Example
    ///
    /// ```
    /// use bijli::fdtd::Fdtd3d;
    ///
    /// let mut sim = Fdtd3d::new(32, 32, 32, 1e-3)?;
    /// sim.try_enable_cpml(8)?;
    /// sim.run(50, Some((16, 16, 16)), 3e10, 1.0);
    /// assert_eq!(sim.cpml().map(|c| c.thickness), Some(8));
    /// # Ok::<(), bijli::error::BijliError>(())
    /// ```
    pub fn try_enable_cpml(&mut self, thickness: usize) -> Result<()> {
        self.try_enable_cpml_with(&CpmlParams::new(thickness))
    }

    /// Enable CPML absorbing boundaries with explicit grading parameters
    /// (κ and CFS-α profiles, grading order, σ_max scale).
    ///
    /// # Errors
    ///
    /// As [`Fdtd3d::try_enable_cpml`], plus any parameter rejected by
    /// [`Cpml::with_params`], or if the Ψ slabs cannot be allocated. The
    /// current boundary configuration is unchanged on error.
    pub fn try_enable_cpml_with(&mut self, params: &CpmlParams) -> Result<()> {
        self.check_grid()?;
        let dims = [self.nx, self.ny, self.nz];
        let t = params.thickness;
        if !cpml_fits(t, self.nx.min(self.ny).min(self.nz)) {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "CPML thickness {t} is invalid for a {}×{}×{} grid (need thickness ≥ 1 \
                     and 2·thickness + {CPML_MIN_INTERIOR} ≤ min(nx, ny, nz))",
                    self.nx, self.ny, self.nz
                ),
            });
        }
        let cpml = Cpml3d::attach(Cpml::with_params(params, self.dx, self.dt)?, dims)?;
        if !cpml.is_attached_to(dims) {
            return Err(BijliError::InvalidParameter {
                reason: format!("cannot attach CPML to a {dims:?} grid"),
            });
        }
        tracing::debug!(
            thickness = t,
            sigma_max = cpml.layer.sigma_max,
            kappa_max = cpml.layer.kappa_max,
            alpha_max = cpml.layer.alpha_max,
            ?dims,
            "3D CPML enabled"
        );
        self.cpml = Some(cpml);
        Ok(())
    }

    /// Remove the CPML and return to reflecting PEC walls.
    pub fn disable_cpml(&mut self) {
        if self.cpml.take().is_some() {
            tracing::debug!("3D CPML disabled");
        }
    }

    /// The active CPML layer (thickness and grading), if any.
    #[inline]
    #[must_use]
    pub fn cpml(&self) -> Option<&Cpml> {
        self.cpml.as_ref().map(|c| &c.layer)
    }

    /// Advance the simulation by one time step using the 3D Yee algorithm.
    ///
    /// Updates all 6 field components and applies the CPML if enabled;
    /// tangential E is then zeroed on all faces (PEC). If the public geometry
    /// fields or vector lengths were changed into an inconsistent state, logs
    /// an error and does nothing.
    pub fn step_once(&mut self) {
        if let Err(e) = self.check_layout() {
            tracing::error!(error = %e, step = self.step, "Fdtd3d step skipped");
            return;
        }
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

        if let Some(p) = self.cpml.as_mut() {
            p.correct_h(
                [&mut self.hx, &mut self.hy, &mut self.hz],
                [&self.ex, &self.ey, &self.ez],
                &self.ch,
            );
        }

        // Lossy cells: E ← C_a·E before adding the curl term.
        scale_by(&mut self.ex, &self.ca);
        scale_by(&mut self.ey, &self.ca);
        scale_by(&mut self.ez, &self.ca);

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

        if let Some(p) = self.cpml.as_mut() {
            p.correct_e(
                [&mut self.ex, &mut self.ey, &mut self.ez],
                [&self.hx, &self.hy, &self.hz],
                &self.ce,
            );
        }

        // ── PEC walls: tangential E = 0 on all six faces ─────
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
    ///
    /// The source is ignored if outside the grid. Does nothing (logging an
    /// error) if the grid state is inconsistent.
    pub fn run(
        &mut self,
        steps: usize,
        source: Option<(usize, usize, usize)>,
        frequency: f64,
        amplitude: f64,
    ) {
        tracing::debug!(
            steps,
            ?source,
            frequency,
            amplitude,
            cpml = ?self.cpml.as_ref().map(|c| c.layer.thickness),
            "3D FDTD run started"
        );
        if let Err(e) = self.check_layout() {
            tracing::error!(error = %e, "Fdtd3d run aborted");
            return;
        }
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
        let vol = self.dx * self.dx * self.dx;
        let mut energy = 0.0;
        let e = self.ex.iter().zip(&self.ey).zip(&self.ez);
        let h = self.hx.iter().zip(&self.hy).zip(&self.hz);
        let mat = self.permittivity.iter().zip(&self.permeability);
        for ((((&ex, &ey), &ez), ((&hx, &hy), &hz)), (&eps, &mu)) in e.zip(h).zip(mat) {
            let e_sq = ex * ex + ey * ey + ez * ez;
            let h_sq = hx * hx + hy * hy + hz * hz;
            energy += 0.5 * eps * EPSILON_0 * e_sq + 0.5 * mu * MU_0 * h_sq;
        }
        energy * vol
    }
}

#[cfg(test)]
mod cpml3d_tests {
    use super::*;
    use crate::material::Material;

    const DX: f64 = 1e-3;
    /// CPML thickness of the reflection test.
    const T: usize = 10;
    /// Probe distance from the source along an axis.
    const P: isize = 4;
    /// Test grid (n³): the source at the centre lies P + 2 cells from the
    /// CPML interface, so axis probes sit 2 cells inside the interior.
    const N: usize = 2 * (P as usize + 2 + T) + 1;
    /// Gaussian width in time steps (≈ 4.3 cells at S ≈ 0.289).
    const TAU_STEPS: f64 = 15.0;
    /// Recording window: covers the pulse, the reflections from the whole
    /// graded layer and the PEC wall behind it, and ends before the wall
    /// reflection of the reference grid reaches any probe.
    const STEPS: usize = 170;
    /// Reference grid (no CPML): its walls are 26 cells from the source.
    /// Checked against a 97³ grid: the two agree to −85 dB over the window.
    const N_REF: usize = 53;

    /// Probe offsets from the source: all six faces, an edge and a corner.
    const PROBES: [(isize, isize, isize); 8] = [
        (-P, 0, 0),
        (P, 0, 0),
        (0, -P, 0),
        (0, P, 0),
        (0, 0, -P),
        (0, 0, P),
        (-P, P, 0),
        (P, -P, -P),
    ];

    /// Soft Gaussian E_z point source at the grid centre; record E_z at the
    /// probes.
    fn record(n: usize, cpml: Option<&CpmlParams>) -> Vec<Vec<f64>> {
        let mut sim = Fdtd3d::new(n, n, n, DX).unwrap();
        if let Some(p) = cpml {
            sim.try_enable_cpml_with(p).unwrap();
        }
        let c = (n / 2) as isize;
        let at = |o: isize| (c + o) as usize;
        let probes: Vec<usize> = PROBES
            .iter()
            .map(|&(ox, oy, oz)| sim.idx(at(ox), at(oy), at(oz)))
            .collect();
        let tau = TAU_STEPS * sim.dt;
        let t0 = 4.0 * tau;
        let mut out = vec![Vec::with_capacity(STEPS); probes.len()];
        for _ in 0..STEPS {
            sim.add_source_ez(n / 2, n / 2, n / 2, gaussian_pulse(sim.time(), t0, tau));
            sim.step_once();
            for (series, &i) in out.iter_mut().zip(&probes) {
                series.push(sim.ez[i]);
            }
        }
        out
    }

    /// Worst-case reflection over all probes: max|f − f_ref| / max|f_ref|.
    fn reflection(test: &[Vec<f64>], reference: &[Vec<f64>]) -> f64 {
        test.iter()
            .zip(reference)
            .map(|(a, r)| {
                let peak = r.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                let err = a
                    .iter()
                    .zip(r)
                    .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
                err / peak
            })
            .fold(0.0, f64::max)
    }

    fn db(r: f64) -> f64 {
        20.0 * r.log10()
    }

    /// CFS-CPML parameters that exercise the κ and α terms.
    fn cfs(thickness: usize, kappa_max: f64) -> CpmlParams {
        let mut p = CpmlParams::new(thickness);
        p.kappa_max = kappa_max;
        p.alpha_max = 2.0 * std::f64::consts::PI * EPSILON_0 * 1e9;
        p
    }

    #[test]
    fn cpml3d_reflection() {
        let classic = CpmlParams::new(T);
        // κ shortens the wavelength inside the layer; with this narrow pulse
        // (c·τ ≈ 4 cells) κ_max = 5 measures ≈ −64 dB (2D alike), κ_max = 3
        // ≈ −77 dB.
        let cfs = cfs(T, 3.0);
        // Independent runs; threads keep the debug-build test short.
        let [reference, wall, classic_rec, cfs_rec] = std::thread::scope(|s| {
            [
                s.spawn(|| record(N_REF, None)),
                s.spawn(|| record(N, None)),
                s.spawn(|| record(N, Some(&classic))),
                s.spawn(|| record(N, Some(&cfs))),
            ]
            .map(|h| h.join().unwrap())
        });
        let pec = reflection(&wall, &reference);
        let classic_r = reflection(&classic_rec, &reference);
        let cfs_r = reflection(&cfs_rec, &reference);
        eprintln!(
            "3D: wall {:.1} dB, CPML({T}) {:.1} dB, CFS-CPML({T}, κ=3) {:.1} dB",
            db(pec),
            db(classic_r),
            db(cfs_r)
        );
        // Without CPML the PEC walls reflect everything back into the probes.
        assert!(pec > 0.3, "reference wall reflection only {pec}");
        // 10-cell CPML measures ≈ −80 dB here (≈ −87 dB for a zero-DC
        // differentiated Gaussian, whose residual does not grow late in the
        // window like the plain Gaussian's static remnant); require ≤ −60 dB.
        assert!(classic_r < 1e-3, "CPML reflection {classic_r}");
        assert!(cfs_r < 1e-3, "CFS-CPML reflection {cfs_r}");
    }

    /// Peak and final total energy for a zero-DC (Ricker) pulse.
    fn energy_history(cpml: Option<usize>) -> (f64, f64) {
        let n = 25;
        let mut sim = Fdtd3d::new(n, n, n, DX).unwrap();
        if let Some(t) = cpml {
            sim.try_enable_cpml(t).unwrap();
        }
        let tau = TAU_STEPS * sim.dt;
        let t0 = 4.0 * tau;
        let mut peak = 0.0_f64;
        for step in 0..260 {
            sim.add_source_ez(n / 2, n / 2, n / 2, ricker_wavelet(sim.time(), t0, tau));
            sim.step_once();
            if step % 5 == 0 {
                peak = peak.max(sim.total_energy());
            }
        }
        (peak, sim.total_energy())
    }

    #[test]
    fn cpml3d_energy_decays() {
        let [(peak_wall, end_wall), (peak_cpml, end_cpml)] = std::thread::scope(|s| {
            [
                s.spawn(|| energy_history(None)),
                s.spawn(|| energy_history(Some(6))),
            ]
            .map(|h| h.join().unwrap())
        });
        eprintln!(
            "3D: wall end/peak {:.3e}, CPML end/peak {:.3e}",
            end_wall / peak_wall,
            end_cpml / peak_cpml
        );
        // The source cell stores most of the peak energy reactively, so a
        // closed PEC box keeps ≈ 10% of the peak once the pulse has left.
        assert!(end_wall / peak_wall > 0.03, "PEC walls must trap energy");
        assert!(
            end_cpml / peak_cpml < 1e-6,
            "CPML must absorb energy: {}",
            end_cpml / peak_cpml
        );
    }

    /// A z-directed source at the centre of a grid that is symmetric about
    /// the source must produce mirror-symmetric fields; any index error in
    /// the low/high slabs or between the x/y slab code paths breaks this.
    #[test]
    fn cpml3d_fields_are_mirror_symmetric() {
        // E_z sits at z + ½: an even nz puts the source E_z on the mid-plane.
        let (n, nz) = (21, 20);
        let (c, cz) = (n / 2, nz / 2 - 1);
        let mut sim = Fdtd3d::new(n, n, nz, DX).unwrap();
        sim.try_enable_cpml_with(&cfs(5, 3.0)).unwrap();
        let tau = 8.0 * sim.dt;
        for _ in 0..120 {
            sim.add_source_ez(c, c, cz, gaussian_pulse(sim.time(), 4.0 * tau, tau));
            sim.step_once();
        }
        let max_abs = |f: &[f64]| f.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let (tol, tol_h) = (1e-9 * max_abs(&sim.ez), 1e-9 * max_abs(&sim.hx));
        assert!(tol > 0.0 && tol_h > 0.0);
        let close = |a: f64, b: f64| (a - b).abs() <= tol;
        let i = |x, y, z| sim.idx(x, y, z);
        for z in 0..nz - 1 {
            for y in 0..n {
                for x in 0..n {
                    let ez = sim.ez[i(x, y, z)];
                    // E_z is even under x, y and z mirrors and under x ↔ y.
                    assert!(
                        close(ez, sim.ez[i(n - 1 - x, y, z)]),
                        "x mirror at {x},{y},{z}"
                    );
                    assert!(
                        close(ez, sim.ez[i(x, n - 1 - y, z)]),
                        "y mirror at {x},{y},{z}"
                    );
                    assert!(
                        close(ez, sim.ez[i(x, y, nz - 2 - z)]),
                        "z mirror at {x},{y},{z}"
                    );
                    assert!(close(ez, sim.ez[i(y, x, z)]), "x↔y at {x},{y},{z}");
                    // x ↔ y swaps E_x ↔ E_y and H_x ↔ −H_y.
                    assert!(close(sim.ex[i(x, y, z)], sim.ey[i(y, x, z)]), "E x↔y");
                    let (hx, hy) = (sim.hx[i(x, y, z)], sim.hy[i(y, x, z)]);
                    assert!((hx + hy).abs() <= tol_h, "H x↔y at {x},{y},{z}");
                }
            }
        }
    }

    #[test]
    fn cpml3d_stable_long_run_with_materials() {
        let mut sim = Fdtd3d::new(19, 19, 19, DX).unwrap();
        sim.try_enable_cpml_with(&cfs(4, 5.0)).unwrap();
        // Dielectric + lossy block that overlaps the x-low PML slab.
        sim.set_material(0, 6, 6, 9, 13, 13, &Material::new(4.0, 2.0, 0.5))
            .unwrap();
        sim.run(500, Some((12, 9, 9)), 3e10, 1.0);
        let fields = [&sim.ex, &sim.ey, &sim.ez, &sim.hx, &sim.hy, &sim.hz];
        assert!(fields.iter().all(|f| f.iter().all(|v| v.is_finite())));
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn cpml3d_enable_validates_and_disables() {
        let mut sim = Fdtd3d::new(20, 30, 40, DX).unwrap();
        assert!(sim.try_enable_cpml(0).is_err());
        // 2·9 + 3 = 21 > min(20, 30, 40)
        assert!(sim.try_enable_cpml(9).is_err());
        assert!(sim.cpml().is_none());
        sim.try_enable_cpml(8).unwrap();
        let cpml = sim.cpml().unwrap();
        assert_eq!(cpml.thickness, 8);
        let sigma_opt = 0.8 * 4.0 / (ETA_0 * DX);
        assert!((cpml.sigma_max() - sigma_opt).abs() < 1e-9 * sigma_opt);
        assert!((cpml.kappa_max() - 1.0).abs() < 1e-15);
        // Invalid requests keep the active layer.
        assert!(sim.try_enable_cpml(100).is_err());
        let mut bad = CpmlParams::new(4);
        bad.kappa_max = 0.5;
        assert!(sim.try_enable_cpml_with(&bad).is_err());
        bad.kappa_max = f64::NAN;
        assert!(sim.try_enable_cpml_with(&bad).is_err());
        assert_eq!(sim.cpml().map(|c| c.thickness), Some(8));
        // Re-enabling replaces the layer.
        sim.try_enable_cpml_with(&cfs(5, 7.0)).unwrap();
        let cpml = sim.cpml().unwrap();
        assert_eq!(cpml.thickness, 5);
        assert!((cpml.kappa_max() - 7.0).abs() < 1e-15);
        sim.disable_cpml();
        assert!(sim.cpml().is_none());
        sim.disable_cpml();
        assert!(sim.cpml().is_none());
    }

    #[test]
    fn cpml3d_disabled_matches_pec() {
        // Disabling the CPML restores the original PEC update exactly.
        let run = |toggle: bool| {
            let mut sim = Fdtd3d::new(11, 11, 11, DX).unwrap();
            if toggle {
                sim.try_enable_cpml(2).unwrap();
                sim.disable_cpml();
            }
            sim.run(40, Some((5, 5, 5)), 3e10, 1.0);
            sim.ez
        };
        assert_eq!(run(false), run(true));
    }

    #[test]
    fn cpml3d_desynced_fields_do_not_panic() {
        let mut sim = Fdtd3d::new(12, 12, 12, DX).unwrap();
        sim.try_enable_cpml(3).unwrap();
        // Same cell count, but the CPML no longer matches the geometry.
        sim.nx = 6;
        sim.ny = 24;
        sim.step_once();
        sim.run(3, Some((2, 2, 2)), 1e9, 1.0);
        assert_eq!(sim.step, 0);
        // The material API only needs a consistent grid.
        assert!(sim.set_permittivity(0, 0, 0, 2, 2, 2, 2.0).is_ok());
        sim.nx = 12;
        sim.ny = 12;
        sim.step_once();
        assert_eq!(sim.step, 1);
        sim.nz = 1000;
        assert!(sim.try_enable_cpml(3).is_err());
        sim.step_once();
        assert_eq!(sim.step, 1);
        sim.nz = 12;
        sim.ex.truncate(5);
        sim.step_once();
        assert!(sim.try_enable_cpml(3).is_err());
        assert_eq!(sim.step, 1);
    }
}

// ── Non-uniform grid ──────────────────────────────────────────────

/// Non-uniform grid specification for FDTD.
///
/// Stores per-cell spacing in each dimension, enabling mesh refinement
/// near material boundaries or objects of interest. A valid grid has at least
/// one cell per axis and only finite, positive spacings (see [`validate`](Self::validate)).
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
    ///
    /// Does not validate: an empty axis or a non-positive `spacing` produces an
    /// invalid grid (for which [`courant_dt`](Self::courant_dt) returns 0.0). Use
    /// [`try_uniform`](Self::try_uniform) for checked construction.
    #[must_use]
    pub fn uniform(nx: usize, ny: usize, nz: usize, spacing: f64) -> Self {
        Self {
            dx: vec![spacing; nx],
            dy: vec![spacing; ny],
            dz: vec![spacing; nz],
        }
    }

    /// Create a uniform grid, validating the parameters.
    ///
    /// # Errors
    /// [`BijliError::InsufficientResolution`] if any axis has zero cells;
    /// [`BijliError::InvalidParameter`] if `spacing` is not finite and positive.
    pub fn try_uniform(nx: usize, ny: usize, nz: usize, spacing: f64) -> Result<Self> {
        let min_cells = nx.min(ny).min(nz);
        if min_cells == 0 {
            return Err(BijliError::InsufficientResolution {
                cells: 0,
                domain_size: 0.0,
            });
        }
        if !is_positive_finite(spacing) {
            return Err(BijliError::InvalidParameter {
                reason: format!("grid spacing must be finite and positive, got {spacing}"),
            });
        }
        Ok(Self::uniform(nx, ny, nz, spacing))
    }

    /// Check that every axis has at least one cell and every spacing is finite
    /// and positive.
    ///
    /// # Errors
    /// [`BijliError::InsufficientResolution`] for an empty axis;
    /// [`BijliError::InvalidParameter`] for a non-finite or non-positive spacing.
    pub fn validate(&self) -> Result<()> {
        for (axis, d) in [("dx", &self.dx), ("dy", &self.dy), ("dz", &self.dz)] {
            if d.is_empty() {
                return Err(BijliError::InsufficientResolution {
                    cells: 0,
                    domain_size: 0.0,
                });
            }
            if let Some(&bad) = d.iter().find(|&&v| !is_positive_finite(v)) {
                return Err(BijliError::InvalidParameter {
                    reason: format!("{axis} contains invalid cell size {bad}"),
                });
            }
        }
        Ok(())
    }

    /// Create a grid with graded spacing (smooth transition from `d_min` to `d_max`).
    ///
    /// Uses geometric grading: d_i = d_min * ratio^i where ratio is computed
    /// to reach d_max at the end.
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] if `n < 2` or either size is not finite
    /// and positive.
    pub fn graded(n: usize, d_min: f64, d_max: f64) -> Result<Vec<f64>> {
        if n < 2 {
            return Err(BijliError::InvalidParameter {
                reason: format!("need at least 2 cells for grading, got {n}"),
            });
        }
        if !is_positive_finite(d_min) || !is_positive_finite(d_max) {
            return Err(BijliError::InvalidParameter {
                reason: format!("cell sizes must be finite and positive, got {d_min}, {d_max}"),
            });
        }
        let ratio = (d_max / d_min).powf(1.0 / (n - 1) as f64);
        Ok((0..n).map(|i| d_min * ratio.powi(i as i32)).collect())
    }

    /// Minimum cell size across all dimensions (`f64::INFINITY` if all axes are empty).
    #[must_use]
    pub fn min_spacing(&self) -> f64 {
        let min_dx = self.dx.iter().copied().fold(f64::INFINITY, f64::min);
        let min_dy = self.dy.iter().copied().fold(f64::INFINITY, f64::min);
        let min_dz = self.dz.iter().copied().fold(f64::INFINITY, f64::min);
        min_dx.min(min_dy).min(min_dz)
    }

    /// Time step at **half** the 3D Courant limit: dt = d_min/(2√3·c).
    ///
    /// (The limit itself is dt ≤ d_min/(√3·c), conservative for non-uniform
    /// grids.) Returns 0.0 — with a `tracing` warning — for an invalid grid
    /// (empty axis or non-positive/non-finite spacing); use
    /// [`try_courant_dt`](Self::try_courant_dt) to get an error instead.
    #[must_use]
    pub fn courant_dt(&self) -> f64 {
        match self.try_courant_dt() {
            Ok(dt) => dt,
            Err(e) => {
                tracing::warn!(%e, "NonUniformGrid::courant_dt on invalid grid, returning 0");
                0.0
            }
        }
    }

    /// Checked variant of [`courant_dt`](Self::courant_dt).
    ///
    /// # Errors
    /// Any error from [`validate`](Self::validate).
    pub fn try_courant_dt(&self) -> Result<f64> {
        self.validate()?;
        Ok(self.min_spacing() / (2.0 * 3.0_f64.sqrt() * SPEED_OF_LIGHT))
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
/// Uses the 4-point staggered stencil
/// f'(x_{i+½}) ≈ (c₁(f\[i+1\]-f\[i\]) + c₂(f\[i+2\]-f\[i-1\])) / dx
/// where c₁ = 9/8, c₂ = -1/24 (Yee fourth-order).
pub const FOURTH_ORDER_C1: f64 = 9.0 / 8.0;
/// Second coefficient for 4th-order stencil.
pub const FOURTH_ORDER_C2: f64 = -1.0 / 24.0;

/// Compute the 4th-order staggered spatial derivative of a 1D field array.
///
/// Output sample `i` approximates f'(x_{i+½}), the midpoint between samples
/// `i` and `i + 1` (the Yee staggering used to compute H from E):
/// - 1 ≤ i ≤ n − 3: (c₁(f\[i+1\]-f\[i\]) + c₂(f\[i+2\]-f\[i-1\])) / dx —
///   fourth-order accurate, exact for cubics;
/// - i = 0 and i = n − 2: (f\[i+1\] − f\[i\]) / dx — second-order accurate;
/// - i = n − 1: 0.0 (padding — x_{n−½} lies outside the data).
///
/// Returns a vector of length n. Fewer than 2 samples, or a `dx` that is not
/// finite and positive, gives all zeros (the latter with a `tracing` warning).
#[must_use]
pub fn fourth_order_derivative(field: &[f64], dx: f64) -> Vec<f64> {
    let n = field.len();
    let mut deriv = vec![0.0; n];
    if n < 2 {
        return deriv;
    }
    if !is_positive_finite(dx) {
        tracing::warn!(dx, "fourth_order_derivative: invalid dx, returning zeros");
        return deriv;
    }
    let inv_dx = 1.0 / dx;

    // Interior: 4th-order stencil (needs i-1 and i+2)
    if n >= 4 {
        for (d, w) in deriv[1..n - 2].iter_mut().zip(field.windows(4)) {
            // w = [f[i-1], f[i], f[i+1], f[i+2]]
            *d = inv_dx * (FOURTH_ORDER_C1 * (w[2] - w[1]) + FOURTH_ORDER_C2 * (w[3] - w[0]));
        }
    }

    // Edge midpoints: 2nd-order central difference about x_{i+½}
    deriv[0] = (field[1] - field[0]) * inv_dx;
    deriv[n - 2] = (field[n - 1] - field[n - 2]) * inv_dx;

    deriv
}

// ── Frequency-domain (CW) solver ──────────────────────────────────

/// Frequency-domain field accumulator for CW (continuous-wave) analysis.
///
/// Accumulates the running DFT F(ω) = Σₙ f(tₙ) e^{−jωtₙ} Δt of time-domain
/// fields at a single frequency, enabling extraction of steady-state amplitude
/// and phase after FDTD run (phases follow the e^{jωt} convention).
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
    /// Whether a field-length mismatch warning has already been emitted.
    mismatch_warned: bool,
}

impl CwAccumulator {
    /// Create a CW accumulator for a field of given size at target frequency.
    ///
    /// Does not validate (a non-finite frequency or non-positive `dt` is
    /// accepted with a `tracing` warning and yields meaningless results); use
    /// [`try_new`](Self::try_new) for checked construction.
    #[must_use]
    pub fn new(field_size: usize, frequency: f64, dt: f64) -> Self {
        if !frequency.is_finite() || !is_positive_finite(dt) {
            tracing::warn!(
                frequency,
                dt,
                "CwAccumulator created with invalid frequency/dt"
            );
        }
        tracing::debug!(field_size, frequency, dt, "CW accumulator created");
        Self {
            omega: 2.0 * std::f64::consts::PI * frequency,
            real: vec![0.0; field_size],
            imag: vec![0.0; field_size],
            n_samples: 0,
            dt,
            mismatch_warned: false,
        }
    }

    /// Create a CW accumulator, validating the parameters.
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] if `frequency` is negative or not finite,
    /// or `dt` is not finite and positive.
    pub fn try_new(field_size: usize, frequency: f64, dt: f64) -> Result<Self> {
        check_non_negative("CW frequency", frequency)?;
        check_positive("dt", dt)?;
        Ok(Self::new(field_size, frequency, dt))
    }

    /// Accumulate one time sample of the field into the running DFT.
    ///
    /// The sample is taken at t = `time_step`·dt. If `field` and the
    /// accumulator differ in length, only the common prefix is accumulated
    /// (a `tracing` warning is emitted once); use
    /// [`try_accumulate`](Self::try_accumulate) to reject mismatches.
    pub fn accumulate(&mut self, field: &[f64], time_step: usize) {
        if (field.len() != self.real.len() || field.len() != self.imag.len())
            && !self.mismatch_warned
        {
            tracing::warn!(
                field_len = field.len(),
                accumulator_len = self.real.len().min(self.imag.len()),
                "CwAccumulator::accumulate length mismatch, using common prefix"
            );
            self.mismatch_warned = true;
        }
        let t = time_step as f64 * self.dt;
        let (sin_wt, cos_wt) = (self.omega * t).sin_cos();
        let cos_dt = cos_wt * self.dt;
        let sin_dt = sin_wt * self.dt;

        for ((r, i), &f) in self.real.iter_mut().zip(self.imag.iter_mut()).zip(field) {
            *r += f * cos_dt;
            *i -= f * sin_dt;
        }
        self.n_samples += 1;
    }

    /// Checked variant of [`accumulate`](Self::accumulate).
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] (without modifying the accumulator) if
    /// `field.len()` differs from the accumulator size.
    pub fn try_accumulate(&mut self, field: &[f64], time_step: usize) -> Result<()> {
        if field.len() != self.real.len() || field.len() != self.imag.len() {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "CwAccumulator expects {} samples, got {}",
                    self.real.len().min(self.imag.len()),
                    field.len()
                ),
            });
        }
        self.accumulate(field, time_step);
        Ok(())
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
    /// Quality factor Q = f₀/Δf per mode, where Δf is the half-power
    /// (|F| = peak/√2) bandwidth around that peak; 0.0 when it cannot be
    /// determined (peak at the band edge or merged with a neighbouring mode).
    pub quality_factors: Vec<f64>,
    /// Spectral amplitudes at each frequency.
    pub amplitudes: Vec<f64>,
}

/// |Σᵢ sᵢ e^{−jω i dt}|·dt using a phasor recurrence (one complex multiply per
/// sample instead of two trig calls), reseeded exactly every block to bound
/// round-off drift.
#[inline]
#[must_use]
fn dft_magnitude(signal: &[f64], omega: f64, dt: f64) -> f64 {
    const RESEED: usize = 1024;
    let (s1, c1) = (omega * dt).sin_cos();
    let (wr, wi) = (c1, -s1); // e^{−jω dt}
    let mut re = 0.0;
    let mut im = 0.0;
    for (block, chunk) in signal.chunks(RESEED).enumerate() {
        let (s0, c0) = (omega * (block * RESEED) as f64 * dt).sin_cos();
        let (mut zr, mut zi) = (c0, -s0);
        for &s in chunk {
            re += s * zr;
            im += s * zi;
            (zr, zi) = (zr * wr - zi * wi, zr * wi + zi * wr);
        }
    }
    (re * re + im * im).sqrt() * dt
}

/// Interpolated frequency offset (in units of `df`, ≥ 0) from `peak` to the
/// half-power crossing, walking outward in direction `dir` (±1). `None` if the
/// spectrum edge or a rising flank (neighbouring mode) is reached first.
#[must_use]
fn half_power_offset(spectrum: &[f64], peak: usize, dir: isize) -> Option<f64> {
    let half_power = spectrum[peak] / std::f64::consts::SQRT_2;
    let mut i = peak;
    loop {
        let next = i.checked_add_signed(dir)?;
        let a_next = *spectrum.get(next)?;
        let a_cur = spectrum[i];
        if a_next < half_power {
            // Linear interpolation between i (≥ hp) and next (< hp).
            let frac = (a_cur - half_power) / (a_cur - a_next);
            return Some(i.abs_diff(peak) as f64 + frac);
        }
        if a_next > a_cur {
            return None; // rising again above half power: merged with another mode
        }
        i = next;
    }
}

/// Find resonant modes of a structure using FDTD + spectral analysis.
///
/// Excites the structure with a broadband pulse, records the time-domain
/// response at a probe point, then performs DFT to identify resonance peaks.
///
/// The spectrum is evaluated at `n_freq` points spanning
/// [`freq_min`, `freq_max`] inclusive. Peaks are strict local maxima above 10 %
/// of the strongest one; each peak's Q is f₀/Δf from its own half-power points,
/// found by walking outward from the peak and interpolating linearly.
///
/// # Parameters
/// - `time_signal`: recorded field values at probe point over time
/// - `dt`: time step between samples (s)
/// - `freq_min`, `freq_max`: frequency search range (Hz)
/// - `n_freq`: number of frequency points to evaluate (≥ 2)
///
/// # Errors
/// [`BijliError::InvalidParameter`] if the signal is empty or contains
/// non-finite values, `n_freq < 2`, `dt` is not finite and positive, or the
/// range is not 0 ≤ `freq_min` < `freq_max` < ∞.
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
    if n_freq < 2 {
        return Err(BijliError::InvalidParameter {
            reason: format!("need at least two frequency points, got {n_freq}"),
        });
    }
    if !is_positive_finite(dt) {
        return Err(BijliError::InvalidParameter {
            reason: format!("dt must be finite and positive, got {dt}"),
        });
    }
    if !(freq_min >= 0.0 && freq_min < freq_max && freq_max.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: format!(
                "need 0 <= freq_min < freq_max < inf, got freq_min = {freq_min}, freq_max = {freq_max}"
            ),
        });
    }
    if time_signal.iter().any(|s| !s.is_finite()) {
        return Err(BijliError::InvalidParameter {
            reason: "time signal contains non-finite values".into(),
        });
    }
    if freq_max > 0.5 / dt {
        tracing::warn!(
            freq_max,
            nyquist = 0.5 / dt,
            "find_eigenmodes: freq_max above Nyquist, spectrum is aliased"
        );
    }

    tracing::debug!(
        n_samples = time_signal.len(),
        dt,
        freq_min,
        freq_max,
        n_freq,
        "finding eigenmodes via DFT"
    );

    let df = (freq_max - freq_min) / (n_freq - 1) as f64;
    let freq_at = |i: usize| freq_min + i as f64 * df;

    // DFT magnitude at each frequency
    let spectrum: Vec<f64> = (0..n_freq)
        .map(|fi| dft_magnitude(time_signal, 2.0 * std::f64::consts::PI * freq_at(fi), dt))
        .collect();

    // Peak detection: find local maxima above threshold
    let max_amp = spectrum.iter().copied().fold(0.0, f64::max);
    let threshold = max_amp * 0.1; // 10% of peak

    let mut frequencies = Vec::new();
    let mut amplitudes = Vec::new();
    let mut quality_factors = Vec::new();

    for i in 1..n_freq - 1 {
        let amp = spectrum[i];
        if amp > threshold && amp > spectrum[i - 1] && amp > spectrum[i + 1] {
            let f0 = freq_at(i);
            let left = half_power_offset(&spectrum, i, -1);
            let right = half_power_offset(&spectrum, i, 1);
            let bw_bins = match (left, right) {
                (Some(l), Some(r)) => Some(l + r),
                (Some(one), None) | (None, Some(one)) => Some(2.0 * one),
                (None, None) => None,
            };
            let q = match bw_bins {
                Some(b) if b > 0.0 => f0 / (b * df),
                _ => 0.0,
            };
            frequencies.push(f0);
            amplitudes.push(amp);
            quality_factors.push(q);
        }
    }

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
        assert!((sim.abc_left - expected).abs() < 1e-15);
        assert!((sim.abc_right - expected).abs() < 1e-15);
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

    #[test]
    fn test_fdtd_invalid_dx_non_finite() {
        // NaN slipped through the old `dx <= 0.0` check.
        assert!(Fdtd1d::new(100, f64::NAN).is_err());
        assert!(Fdtd1d::new(100, f64::INFINITY).is_err());
    }

    #[test]
    fn test_fdtd_dt_is_half_courant_limit() {
        let sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!((SPEED_OF_LIGHT * sim.dt / sim.dx - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_set_material_rejects_non_finite() {
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        assert!(sim.set_permittivity(0, 100, f64::NAN).is_err());
        assert!(sim.set_permittivity(0, 100, f64::INFINITY).is_err());
        assert!(sim.set_permeability(0, 100, f64::NAN).is_err());
        assert!((sim.permittivity[50] - 1.0).abs() < 1e-15);
        assert!((sim.permeability[50] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_mur_coefficient_tracks_material() {
        // k = (v·dt − dx)/(v·dt + dx) with v = c/√(ε_r μ_r): S_eff = 0.5/√(ε_r μ_r)
        let mut sim = Fdtd1d::new(100, 1e-3).unwrap();
        sim.set_permittivity(0, 100, 4.0).unwrap();
        let k = |s: f64| (s - 1.0) / (s + 1.0);
        assert!((sim.abc_left - k(0.25)).abs() < 1e-12);
        assert!((sim.abc_right - k(0.25)).abs() < 1e-12);
        sim.set_permeability(50, 100, 4.0).unwrap();
        assert!((sim.abc_left - k(0.25)).abs() < 1e-12);
        assert!((sim.abc_right - k(0.125)).abs() < 1e-12);
    }

    /// Fraction of a Gaussian pulse's energy left in a homogeneous grid after
    /// both halves of the pulse have hit the Mur boundaries (and any reflection
    /// is still inside). A perfect ABC gives 0.
    fn fdtd1d_reflected_energy_fraction(eps_r: f64) -> f64 {
        let n = 600;
        let mut sim = Fdtd1d::new(n, 1e-4).unwrap();
        sim.set_permittivity(0, n, eps_r).unwrap();
        let cells_per_step = 0.5 / eps_r.sqrt();
        let tau = 40.0 * sim.dt;
        let t0 = 6.0 * tau;
        let src_steps = 480; // 12τ
        let travel = (300.0 / cells_per_step) as usize;
        let mut e_ref = 0.0;
        for step in 0..src_steps + travel + travel / 3 {
            if step < src_steps {
                let v = gaussian_pulse(sim.time(), t0, tau);
                sim.add_source(n / 2, v);
            }
            sim.step_once();
            if step == src_steps {
                e_ref = sim.total_energy();
            }
        }
        sim.total_energy() / e_ref
    }

    #[test]
    fn test_mur_abc_absorbs_in_vacuum() {
        let frac = fdtd1d_reflected_energy_fraction(1.0);
        assert!(frac < 1e-3, "reflected energy fraction {frac}");
    }

    #[test]
    fn test_mur_abc_absorbs_in_dielectric() {
        // Regression: the Mur coefficient was fixed at the vacuum value, giving
        // |R| ≈ 1/3 (≈ 11 % of the energy) in an ε_r = 4 medium.
        let frac = fdtd1d_reflected_energy_fraction(4.0);
        assert!(frac < 1e-3, "reflected energy fraction {frac}");
    }

    #[test]
    fn test_fdtd1d_inconsistent_public_fields_no_panic() {
        let mut sim = Fdtd1d::new(50, 1e-3).unwrap();
        sim.num_cells = 1000;
        sim.add_source(500, 1.0);
        sim.step_once();
        let _ = sim.total_energy();
        sim.e_field.truncate(1);
        sim.step_once();
        assert_eq!(sim.step, 2);
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
        assert_eq!(cpml.b_node.len(), 10);
    }

    #[test]
    fn test_cpml_coefficients_bounded() {
        let cpml = Cpml::new(10, 1e-4, 1e-13);
        for &b in &cpml.b_node {
            assert!((0.0..=1.0).contains(&b), "b_e out of range: {b}");
        }
        for &b in &cpml.b_half {
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
        assert!(!cpml.psi_fz_y.is_empty());
    }

    #[test]
    fn test_fdtd2d_cpml_runs_without_panic() {
        let mut sim = Fdtd2d::new(60, 60, 1e-4, Mode2d::Tm).unwrap();
        sim.enable_cpml(8);
        sim.run(30, Some((30, 30)), 1e9, 1.0);
        // Should complete without panic; the source is still radiating.
        assert!(sim.total_energy() > 0.0);
    }

    #[test]
    fn test_cpml_attach_sizes() {
        let mut cpml = Cpml::new(8, 1e-4, 1e-13);
        cpml.attach(100, 80);
        assert_eq!(cpml.psi_fz_y.len(), 2 * 8 * 100);
        assert_eq!(cpml.psi_fz_x.len(), 2 * 8 * 80);
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

    // These two tests previously asserted the odd first-derivative waveform
    // under the name `ricker_wavelet`; that waveform is `differentiated_gaussian`.
    #[test]
    fn test_differentiated_gaussian_zero_at_center() {
        assert!(differentiated_gaussian(1e-9, 1e-9, 1e-10).abs() < 1e-10);
    }

    #[test]
    fn test_differentiated_gaussian_antisymmetric() {
        let t0 = 1e-9;
        let tau = 1e-10;
        let v1 = differentiated_gaussian(t0 - 0.5e-10, t0, tau);
        let v2 = differentiated_gaussian(t0 + 0.5e-10, t0, tau);
        assert!((v1 + v2).abs() / v1.abs() < 1e-12); // antisymmetric
    }

    #[test]
    fn test_differentiated_gaussian_is_gaussian_derivative() {
        // Reference: central finite difference of gaussian_pulse.
        let (t0, tau) = (1e-9, 1e-10);
        let h = 1e-15;
        for i in -20..=20 {
            let t = t0 + f64::from(i) * 0.15 * tau;
            let fd = (gaussian_pulse(t + h, t0, tau) - gaussian_pulse(t - h, t0, tau)) / (2.0 * h);
            let v = differentiated_gaussian(t, t0, tau);
            assert!((v - fd).abs() < 1e-6 / tau, "t = {t}: {v} vs {fd}");
        }
        // Peak magnitude sqrt(2)·e^(−1/2)/τ at t = t0 − τ/sqrt(2)
        let peak = differentiated_gaussian(t0 - tau / 2f64.sqrt(), t0, tau);
        assert!((peak - 2f64.sqrt() * (-0.5f64).exp() / tau).abs() * tau < 1e-12);
    }

    #[test]
    fn test_ricker_wavelet_peak_and_symmetry() {
        // Ricker (Mexican hat): even, peak 1 at t0, zeros at t0 ± τ/√2,
        // minima −2e^{−3/2} at t0 ± τ√(3/2).
        let (t0, tau) = (1e-9, 1e-10);
        assert!((ricker_wavelet(t0, t0, tau) - 1.0).abs() < 1e-15);
        for i in 1..30 {
            let off = f64::from(i) * 0.1 * tau;
            let (a, b) = (
                ricker_wavelet(t0 - off, t0, tau),
                ricker_wavelet(t0 + off, t0, tau),
            );
            assert!((a - b).abs() < 1e-14, "not even at offset {off}");
            assert!(a <= 1.0);
        }
        assert!(ricker_wavelet(t0 + tau / 2f64.sqrt(), t0, tau).abs() < 1e-14);
        let min = ricker_wavelet(t0 + tau * 1.5f64.sqrt(), t0, tau);
        assert!((min + 2.0 * (-1.5f64).exp()).abs() < 1e-14);
    }

    #[test]
    fn test_ricker_wavelet_schneider_form_and_peak_frequency() {
        // Schneider's form (1 − 2π²f_p²t²)·exp(−π²f_p²t²) with f_p = 1/(πτ).
        let (t0, tau) = (2e-9, 1e-10);
        let fp = 1.0 / (std::f64::consts::PI * tau);
        let pi2fp2 = (std::f64::consts::PI * fp).powi(2);
        for i in -40..=40 {
            let t = t0 + f64::from(i) * 0.1 * tau;
            let u2 = pi2fp2 * (t - t0).powi(2);
            let reference = (1.0 - 2.0 * u2) * (-u2).exp();
            assert!((ricker_wavelet(t, t0, tau) - reference).abs() < 1e-12);
        }
        // Its amplitude spectrum ∝ f²·exp(−f²/f_p²) peaks at f_p: check by
        // direct numerical Fourier integration.
        let spectrum = |f: f64| {
            let dt = tau / 200.0;
            let (mut re, mut im) = (0.0, 0.0);
            for n in 0..8000 {
                let t = n as f64 * dt;
                let w = 2.0 * std::f64::consts::PI * f * t;
                re += ricker_wavelet(t, t0, tau) * w.cos() * dt;
                im -= ricker_wavelet(t, t0, tau) * w.sin() * dt;
            }
            re.hypot(im)
        };
        let at_peak = spectrum(fp);
        assert!(at_peak > spectrum(0.9 * fp));
        assert!(at_peak > spectrum(1.1 * fp));
        // Ratio vs analytic f²·exp(−f²/fp²) at f = fp/2: (1/4)·e^{3/4}
        let ratio = spectrum(0.5 * fp) / at_peak;
        assert!(
            (ratio - 0.25 * 0.75f64.exp()).abs() < 1e-6,
            "ratio = {ratio}"
        );
    }

    #[test]
    fn test_source_functions_zero_tau_no_nan() {
        for t in [0.0, 1e-9, 2e-9] {
            assert_eq!(gaussian_pulse(t, 1e-9, 0.0), 0.0);
            assert_eq!(ricker_wavelet(t, 1e-9, 0.0), 0.0);
            assert_eq!(differentiated_gaussian(t, 1e-9, 0.0), 0.0);
            assert_eq!(modulated_gaussian_pulse(t, 1e-9, 0.0, 1e9), 0.0);
            assert_eq!(gaussian_pulse(t, 1e-9, f64::NAN), 0.0);
            assert_eq!(ricker_wavelet(t, 1e-9, f64::NAN), 0.0);
        }
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

    #[test]
    fn test_tfsf_try_new_validation() {
        let dx = 1e-4;
        let dt = dx / (2.0 * SPEED_OF_LIGHT);
        assert!(TfsfSource::try_new(100, 50, dx, dt).is_ok());
        assert!(TfsfSource::try_new(100, 50, dx, dx / SPEED_OF_LIGHT).is_ok());
        assert!(TfsfSource::try_new(0, 0, dx, dt).is_err());
        assert!(TfsfSource::try_new(2, 1, dx, dt).is_err());
        assert!(TfsfSource::try_new(100, 0, dx, dt).is_err());
        assert!(TfsfSource::try_new(100, 100, dx, dt).is_err());
        assert!(TfsfSource::try_new(100, 50, 0.0, dt).is_err());
        assert!(TfsfSource::try_new(100, 50, f64::NAN, dt).is_err());
        assert!(TfsfSource::try_new(100, 50, dx, 0.0).is_err());
        assert!(TfsfSource::try_new(100, 50, dx, f64::NAN).is_err());
        assert!(TfsfSource::try_new(100, 50, dx, 1.01 * dx / SPEED_OF_LIGHT).is_err());
    }

    #[test]
    fn test_tfsf_new_sanitizes_without_panic() {
        // Each case used to panic (usize underflow / out-of-bounds) or diverge
        // (dt > dx/c).
        let cases = [
            (0, 0, 1e-4, 1e-13),
            (1, 0, 1e-4, 1e-13),
            (100, 500, 1e-4, 1e-13),
            (100, 50, f64::NAN, 1e-13),
            (100, 50, 1e-4, f64::NAN),
            (100, 50, 1e-4, 1e-11),
            (100, 50, -1.0, -1.0),
        ];
        for (nx, bx, dx, dt) in cases {
            let mut tfsf = TfsfSource::new(nx, bx, dx, dt);
            assert!(tfsf.inc_e.len() >= 3);
            assert!(tfsf.boundary_x >= 1 && tfsf.boundary_x < tfsf.inc_e.len());
            for step in 0..2000 {
                tfsf.step(gaussian_pulse(step as f64, 30.0, 10.0));
            }
            assert!(
                tfsf.inc_e.iter().chain(&tfsf.inc_h).all(|v| v.is_finite()),
                "diverged for {nx} {bx} {dx} {dt}"
            );
            assert!(tfsf.incident_e().is_finite() && tfsf.incident_h().is_finite());
        }
    }

    #[test]
    fn test_tfsf_mutated_public_fields_no_panic() {
        let mut tfsf = TfsfSource::try_new(50, 20, 1e-4, 1e-13).unwrap();
        tfsf.boundary_x = 1000;
        tfsf.source_cell = 1000;
        tfsf.step(1.0);
        assert_eq!(tfsf.incident_e(), 0.0);
        assert_eq!(tfsf.incident_h(), 0.0);
        tfsf.boundary_x = 0;
        assert_eq!(tfsf.incident_h(), 0.0);
        tfsf.inc_h.truncate(2);
        tfsf.step(1.0); // degenerate grid: no-op, no panic
    }

    #[test]
    fn test_tfsf_small_boundary_injects() {
        // Regression: boundary_x < 4 put the source on node 0, where the
        // boundary condition overwrote it and the incident field stayed 0.
        let dx = 1e-4;
        let dt = dx / (2.0 * SPEED_OF_LIGHT);
        for bx in 1..4 {
            let mut tfsf = TfsfSource::try_new(50, bx, dx, dt).unwrap();
            assert!(tfsf.source_cell >= 1);
            for step in 0..60 {
                tfsf.step(gaussian_pulse(step as f64 * dt, 20.0 * dt, 6.0 * dt));
            }
            assert!(tfsf.inc_e.iter().any(|v| v.abs() > 1e-3), "bx = {bx}");
        }
    }

    #[test]
    fn test_tfsf_aux_grid_energy_leaves() {
        // With Mur ABCs at both ends the auxiliary grid must empty once the
        // pulse has passed. Previously the right end was a PEC and the left a
        // PMC, so all energy stayed in the grid forever.
        let nx = 200;
        let dx = 1e-4;
        let dt = dx / (2.0 * SPEED_OF_LIGHT);
        let mut tfsf = TfsfSource::try_new(nx, 100, dx, dt).unwrap();
        let energy = |t: &TfsfSource| {
            t.inc_e.iter().map(|e| 0.5 * EPSILON_0 * e * e).sum::<f64>()
                + t.inc_h.iter().map(|h| 0.5 * MU_0 * h * h).sum::<f64>()
        };
        let tau = 20.0 * dt;
        let t0 = 6.0 * tau;
        let mut peak: f64 = 0.0;
        for step in 0..2000 {
            tfsf.step(gaussian_pulse(step as f64 * dt, t0, tau));
            peak = peak.max(energy(&tfsf));
        }
        let residual = energy(&tfsf) / peak;
        assert!(residual < 1e-6, "residual energy fraction {residual}");
    }

    #[test]
    fn test_tfsf_consistency_with_fdtd1d() {
        // TFSF consistency (Schneider §3.10): in a homogeneous main grid with
        // the left-edge corrections, the scattered-field region (x < bx) stays
        // empty and the total-field region reproduces the incident field. This
        // only holds when incident_h() is H^inc at bx − ½.
        let n = 300;
        let bx = 120;
        let dx = 1e-4;
        let mut sim = Fdtd1d::new(n, dx).unwrap();
        let mut tfsf = TfsfSource::try_new(n, bx, dx, sim.dt).unwrap();
        let ce = sim.dt / (EPSILON_0 * dx);
        let ch = sim.dt / (MU_0 * dx);
        let tau = 15.0 * sim.dt;
        let t0 = 6.0 * tau;
        let (mut max_inc, mut max_sf, mut max_tf_err) = (0.0_f64, 0.0_f64, 0.0_f64);
        for step in 0..900 {
            // H[bx−1] (scattered field) sees total E[bx]: subtract E^inc(bx, n).
            // Applied before step_once, which adds the regular curl term.
            sim.h_field[bx - 1] -= ch * tfsf.incident_e();
            tfsf.step(gaussian_pulse(step as f64 * sim.dt, t0, tau));
            sim.step_once();
            // E[bx] (total field) sees scattered H[bx−1]: subtract H^inc(bx−½, n+½).
            sim.e_field[bx] -= ce * tfsf.incident_h();

            max_inc = tfsf.inc_e[bx..].iter().fold(max_inc, |m, v| m.max(v.abs()));
            max_sf = sim.e_field[..bx].iter().fold(max_sf, |m, v| m.max(v.abs()));
            max_tf_err = sim.e_field[bx..]
                .iter()
                .zip(&tfsf.inc_e[bx..])
                .fold(max_tf_err, |m, (a, b)| m.max((a - b).abs()));
        }
        assert!(max_inc > 1e-2, "incident wave never reached the boundary");
        assert!(max_sf < 1e-9 * max_inc, "SF leakage {}", max_sf / max_inc);
        assert!(
            max_tf_err < 1e-9 * max_inc,
            "TF error {}",
            max_tf_err / max_inc
        );
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

    /// Complex amplitude X of a steady-state sequence x_n = Re[X e^{jωnΔt}],
    /// projected over the last `periods` whole periods (DFT orthogonality).
    fn steady_state_phasor(x: &[f64], samples_per_period: usize, periods: usize) -> (f64, f64) {
        let m = samples_per_period * periods;
        let start = x.len() - m;
        let w = 2.0 * std::f64::consts::PI / samples_per_period as f64;
        let (mut re, mut im) = (0.0, 0.0);
        for (k, &v) in x[start..].iter().enumerate() {
            let th = w * (start + k) as f64;
            re += v * th.cos();
            im -= v * th.sin();
        }
        (2.0 * re / m as f64, 2.0 * im / m as f64)
    }

    fn complex_rel_err(a: (f64, f64), b: (f64, f64)) -> f64 {
        (a.0 - b.0).hypot(a.1 - b.1) / b.0.hypot(b.1)
    }

    fn complex_div(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    }

    #[test]
    fn test_drude_ade_matches_analytic_response() {
        // Drive J^{n+1} = αJ^n + β(E^{n+1} + E^n) with E^n = cos(ωnΔt); the
        // steady state must match the analytic Drude current
        // J = ε₀ω_p²/(γ + jω)·E  (from dJ/dt + γJ = ε₀ω_p²E, e^{jωt}).
        // With β doubled (the formerly documented E^{n+1/2} reading) the
        // response would be off by a factor of 2.
        let mat = DrudeMaterial::try_new(1.0, 1.0e15, 1.0e14).unwrap();
        let dt = 1e-17;
        let ppp = 3000; // samples per period
        let omega = 2.0 * std::f64::consts::PI / (ppp as f64 * dt);
        let (alpha, beta, _) = mat.ade_coefficients(dt);
        let steps = 12 * ppp; // transient decays as e^{−γt}; γt = 36 at the end
        let e = |n: usize| (omega * n as f64 * dt).cos();
        let mut j = vec![0.0; steps];
        for n in 0..steps - 1 {
            j[n + 1] = alpha * j[n] + beta * (e(n + 1) + e(n));
        }
        let got = steady_state_phasor(&j, ppp, 2);
        let num = (EPSILON_0 * mat.omega_p * mat.omega_p, 0.0);
        let analytic = complex_div(num, (mat.gamma, omega));
        let err = complex_rel_err(got, analytic);
        assert!(err < 1e-4, "Drude ADE relative error {err}");
    }

    #[test]
    fn test_lorentz_ade_matches_analytic_response() {
        // P^{n+1} = αP^n + βP^{n−1} + γ_c E^n vs the analytic oscillator
        // response P = ε₀Δεω₀²/(ω₀² − ω² + jγω)·E.
        let mat = LorentzMaterial::try_new(2.0, 2.0, 1.0e15, 1.0e14).unwrap();
        let dt = 1e-17;
        let ppp = 785; // ω ≈ 0.8 ω₀
        let omega = 2.0 * std::f64::consts::PI / (ppp as f64 * dt);
        let (alpha, beta, gamma_c) = mat.ade_coefficients(dt);
        let steps = 110 * ppp; // transient decays as e^{−γt/2}; γt/2 ≈ 43 at the end
        let e = |n: usize| (omega * n as f64 * dt).cos();
        let mut p = vec![0.0; steps];
        for n in 1..steps - 1 {
            p[n + 1] = alpha * p[n] + beta * p[n - 1] + gamma_c * e(n);
        }
        let got = steady_state_phasor(&p, ppp, 4);
        let w0 = mat.omega_0;
        let num = (EPSILON_0 * mat.delta_eps * w0 * w0, 0.0);
        let analytic = complex_div(num, (w0 * w0 - omega * omega, mat.gamma * omega));
        let err = complex_rel_err(got, analytic);
        assert!(err < 1e-3, "Lorentz ADE relative error {err}");
    }

    #[test]
    fn test_debye_ade_matches_analytic_response() {
        // P^{n+1} = αP^n + β(E^{n+1} + E^n) vs P = ε₀(ε_s − ε_∞)/(1 + jωτ)·E.
        let mat = DebyeMaterial::try_new(4.9, 78.4, 1e-11).unwrap();
        let dt = 1e-14;
        let ppp = 6000; // ωτ ≈ 1.05
        let omega = 2.0 * std::f64::consts::PI / (ppp as f64 * dt);
        let (alpha, beta) = mat.ade_coefficients(dt);
        let steps = 8 * ppp; // transient decays as e^{−t/τ}; t/τ = 48 at the end
        let e = |n: usize| (omega * n as f64 * dt).cos();
        let mut p = vec![0.0; steps];
        for n in 0..steps - 1 {
            p[n + 1] = alpha * p[n] + beta * (e(n + 1) + e(n));
        }
        let got = steady_state_phasor(&p, ppp, 2);
        let num = (EPSILON_0 * (mat.eps_s - mat.eps_inf), 0.0);
        let analytic = complex_div(num, (1.0, omega * mat.tau));
        let err = complex_rel_err(got, analytic);
        assert!(err < 1e-4, "Debye ADE relative error {err}");
    }

    #[test]
    fn test_material_try_new_validation() {
        assert_eq!(
            DrudeMaterial::try_new(1.0, 1.37e16, 4.05e13).unwrap(),
            DrudeMaterial::gold()
        );
        assert!(DrudeMaterial::try_new(0.0, 1e16, 1e13).is_err()); // κ = inf
        assert!(DrudeMaterial::try_new(f64::NAN, 1e16, 1e13).is_err());
        assert!(DrudeMaterial::try_new(1.0, -1.0, 1e13).is_err());
        assert!(DrudeMaterial::try_new(1.0, 1e16, f64::INFINITY).is_err());

        assert!(LorentzMaterial::try_new(1.0, 2.0, 1e15, 1e13).is_ok());
        assert!(LorentzMaterial::try_new(0.0, 2.0, 1e15, 1e13).is_err());
        assert!(LorentzMaterial::try_new(1.0, -2.0, 1e15, 1e13).is_err());
        assert!(LorentzMaterial::try_new(1.0, 2.0, 0.0, 1e13).is_err());
        assert!(LorentzMaterial::try_new(1.0, 2.0, 1e15, f64::NAN).is_err());

        assert_eq!(
            DebyeMaterial::try_new(4.9, 78.4, 8.27e-12).unwrap(),
            DebyeMaterial::water()
        );
        assert!(DebyeMaterial::try_new(0.0, 80.0, 1e-11).is_err());
        assert!(DebyeMaterial::try_new(80.0, 2.0, 1e-11).is_err());
        assert!(DebyeMaterial::try_new(2.0, 80.0, 0.0).is_err());
        assert!(DebyeMaterial::try_new(2.0, f64::NAN, 1e-11).is_err());

        assert!(LossyMaterial::try_new(4.0, 0.5).is_ok());
        assert!(LossyMaterial::try_new(-1.0, 0.5).is_err());
        assert!(LossyMaterial::try_new(4.0, -0.5).is_err());
        assert!(LossyMaterial::try_new(4.0, f64::NAN).is_err());
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

    #[test]
    fn test_nf2ff_snapshot_validation() {
        let te = Fdtd2d::new(30, 30, 1e-4, Mode2d::Te).unwrap();
        assert!(nf2ff_2d(&te, 3, 8, 1e9).is_err()); // TE used to be treated as TM
        let tm = Fdtd2d::new(30, 30, 1e-4, Mode2d::Tm).unwrap();
        assert!(nf2ff_2d(&tm, 0, 8, 1e9).is_err()); // contour on the PEC edge
        assert!(nf2ff_2d(&tm, 3, 8, 0.0).is_err());
        assert!(nf2ff_2d(&tm, 3, 8, f64::NAN).is_err());
        let mut short = tm.clone();
        short.fx.truncate(10);
        assert!(nf2ff_2d(&short, 3, 8, 1e9).is_err());
    }

    #[test]
    fn test_nf2ff_snapshot_is_magnitude() {
        let n = 60;
        let mut sim = Fdtd2d::new(n, n, 1e-4, Mode2d::Tm).unwrap();
        sim.run(100, Some((n / 2, n / 2)), 3e9, 1.0);
        let result = nf2ff_2d(&sim, 5, 36, 3e9).unwrap();
        assert!(result.pattern.iter().all(|&v| v >= 0.0 && v.is_finite()));
        assert!(result.pattern.iter().any(|&v| v > 0.0));
    }

    /// Composite Simpson rule on [a, b] with `m` (even) intervals.
    fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, m: usize) -> f64 {
        let h = (b - a) / m as f64;
        let mut acc = f(a) + f(b);
        for i in 1..m {
            acc += f(a + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
        }
        acc * h / 3.0
    }

    /// Bessel J_n(x) from Bessel's integral (1/2π)∫₀^{2π} cos(nτ − x sin τ) dτ
    /// (periodic trapezoid rule, exponentially convergent).
    fn bessel_j(n: i32, x: f64) -> f64 {
        let m = 512;
        let h = 2.0 * std::f64::consts::PI / m as f64;
        (0..m)
            .map(|i| {
                let t = i as f64 * h;
                (f64::from(n) * t - x * t.sin()).cos()
            })
            .sum::<f64>()
            / m as f64
    }

    /// Bessel Y_n(x), x > 0, from the integral representation (DLMF §10.9)
    /// Y_n(x) = (1/π)∫₀^π sin(x sin τ − nτ)dτ
    ///          − (1/π)∫₀^∞ (e^{nt} + (−1)ⁿe^{−nt}) e^{−x sinh t} dt.
    fn bessel_y(n: i32, x: f64) -> f64 {
        let nf = f64::from(n);
        let pi = std::f64::consts::PI;
        let a = simpson(|t| (x * t.sin() - nf * t).sin(), 0.0, pi, 2000);
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        let upper = (60.0 / x).asinh() + 1.0; // e^{−x sinh t} < e^{−60} beyond
        let b = simpson(
            |t| ((nf * t).exp() + sign * (-nf * t).exp()) * (-x * t.sinh()).exp(),
            0.0,
            upper,
            2000,
        );
        (a - b) / pi
    }

    #[test]
    fn test_bessel_reference_helpers() {
        // Abramowitz & Stegun Table 9.1
        let table = [
            (
                1.0,
                0.765_197_686_6,
                0.088_256_964_2,
                0.440_050_585_7,
                -0.781_212_821_3,
            ),
            (
                10.0,
                -0.245_935_764_5,
                0.055_671_167_3,
                0.043_472_746_2,
                0.249_015_424_2,
            ),
        ];
        for (x, j0, y0, j1, y1) in table {
            assert!((bessel_j(0, x) - j0).abs() < 1e-9);
            assert!((bessel_y(0, x) - y0).abs() < 1e-9);
            assert!((bessel_j(1, x) - j1).abs() < 1e-9);
            assert!((bessel_y(1, x) - y1).abs() < 1e-9);
        }
    }

    /// Fill accumulators with exact TM phasors near the NF2FF contour of an
    /// nx×ny grid. `field(x, y, c)` returns the phasor of component c
    /// (0 = E_z, 1 = H_x, 2 = H_y) at the physical point (x, y); H is delayed
    /// by dt/2 like FDTD-accumulated H.
    fn nf2ff_fill(
        nx: usize,
        ny: usize,
        dx: f64,
        margin: usize,
        freq: f64,
        dt: f64,
        field: &dyn Fn(f64, f64, usize) -> (f64, f64),
    ) -> (CwAccumulator, CwAccumulator, CwAccumulator) {
        let n = nx * ny;
        let mut ez = CwAccumulator::new(n, freq, dt);
        let mut hx = CwAccumulator::new(n, freq, dt);
        let mut hy = CwAccumulator::new(n, freq, dt);
        let (x0, x1, y0, y1) = (margin, nx - 1 - margin, margin, ny - 1 - margin);
        let (sl, cl) = (-std::f64::consts::PI * freq * dt).sin_cos();
        let lag = |(r, i): (f64, f64)| (r * cl - i * sl, r * sl + i * cl);
        for y in 0..ny {
            for x in 0..nx {
                let in_box = x + 1 >= x0 && x <= x1 && y + 1 >= y0 && y <= y1;
                let near_edge = x <= x0 || x + 1 >= x1 || y <= y0 || y + 1 >= y1;
                if !(in_box && near_edge) {
                    continue; // never read by the transform
                }
                let i = y * nx + x;
                let (xf, yf) = (x as f64 * dx, y as f64 * dx);
                (ez.real[i], ez.imag[i]) = field(xf, yf, 0);
                (hx.real[i], hx.imag[i]) = lag(field(xf, yf + 0.5 * dx, 1));
                (hy.real[i], hy.imag[i]) = lag(field(xf + 0.5 * dx, yf, 2));
            }
        }
        (ez, hx, hy)
    }

    #[test]
    fn test_nf2ff_phasor_line_source() {
        // Exact fields of an off-centre line current, E_z = H₀⁽²⁾(kρ)
        // (e^{jωt}): H_x = −(j/η)H₁⁽²⁾ sin φ, H_y = (j/η)H₁⁽²⁾ cos φ. The far
        // field is isotropic with √ρ|E_z| → √(2/(πk)). Getting this needs both
        // J and M, complex accumulation and correct Yee positions.
        let freq = 1e9;
        let lambda = SPEED_OF_LIGHT / freq;
        let dx = lambda / 30.0;
        let (nx, ny, margin) = (72, 72, 6);
        let k = 2.0 * std::f64::consts::PI / lambda;
        let eta = (MU_0 / EPSILON_0).sqrt();
        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);
        let (sx, sy) = (35.5 * dx + 0.37 * lambda, 35.5 * dx - 0.21 * lambda);
        let field = |x: f64, y: f64, component: usize| {
            let (px, py) = (x - sx, y - sy);
            let r = px.hypot(py);
            if component == 0 {
                return (bessel_j(0, k * r), -bessel_y(0, k * r)); // H₀⁽²⁾ = J₀ − jY₀
            }
            let h1 = (bessel_j(1, k * r), -bessel_y(1, k * r));
            if component == 1 {
                let s = py / r;
                (h1.1 * s / eta, -h1.0 * s / eta) // −j·(a + jb) = b − ja
            } else {
                let c = px / r;
                (-h1.1 * c / eta, h1.0 * c / eta) // +j·(a + jb) = −b + ja
            }
        };
        let (ez, hx, hy) = nf2ff_fill(nx, ny, dx, margin, freq, dt, &field);
        let res = nf2ff_2d_phasor(&ez, &hx, &hy, nx, ny, dx, margin, 72).unwrap();
        let expected = (2.0 / (std::f64::consts::PI * k)).sqrt();
        for (&phi, &p) in res.angles.iter().zip(&res.pattern) {
            let err = (p - expected).abs() / expected;
            // Discretization leaves ≈ 0.4 %; J only or real part only ≈ 100 %,
            // no dt/2 H alignment ≈ 3 %.
            assert!(
                err < 0.01,
                "φ = {phi}: pattern {p}, expected {expected} (err {err})"
            );
        }
    }

    #[test]
    fn test_nf2ff_phasor_plane_wave_extinction() {
        // A plane wave crossing the contour has no source inside, so the
        // equivalent currents radiate nothing outside (extinction theorem).
        // J or M alone would not cancel.
        let freq = 1e9;
        let lambda = SPEED_OF_LIGHT / freq;
        let dx = lambda / 30.0;
        let (nx, ny, margin) = (72, 72, 6);
        let k = 2.0 * std::f64::consts::PI / lambda;
        let eta = (MU_0 / EPSILON_0).sqrt();
        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);
        let alpha: f64 = 0.3;
        let (sa, ca) = alpha.sin_cos();
        let field = |x: f64, y: f64, component: usize| {
            let (s, c) = (-k * (x * ca + y * sa)).sin_cos();
            match component {
                0 => (c, s),
                1 => (sa * c / eta, sa * s / eta),
                _ => (-ca * c / eta, -ca * s / eta),
            }
        };
        let (ez, hx, hy) = nf2ff_fill(nx, ny, dx, margin, freq, dt, &field);
        let res = nf2ff_2d_phasor(&ez, &hx, &hy, nx, ny, dx, margin, 72).unwrap();
        // Scale: one contour side radiating coherently, √(k/8π)·L·|E|.
        let side = (nx - 1 - 2 * margin) as f64 * dx;
        let scale = (k / (8.0 * std::f64::consts::PI)).sqrt() * side;
        let max = res.pattern.iter().copied().fold(0.0, f64::max);
        // Discretization leaves ≈ 0.7 %; without the dt/2 H alignment ≈ 5 %,
        // with J only ≈ 150 %.
        assert!(
            max < 0.02 * scale,
            "residual far field {} of scale",
            max / scale
        );
    }

    #[test]
    fn test_nf2ff_phasor_validation() {
        let (nx, ny, dx, dt) = (20, 20, 1e-3, 1e-12);
        let a = CwAccumulator::new(nx * ny, 1e9, dt);
        let b = CwAccumulator::new(nx * ny, 2e9, dt);
        let short = CwAccumulator::new(10, 1e9, dt);
        let zero = CwAccumulator::new(nx * ny, 0.0, dt);
        assert!(nf2ff_2d_phasor(&a, &a, &a, nx, ny, dx, 3, 8).is_ok());
        assert!(nf2ff_2d_phasor(&a, &b, &a, nx, ny, dx, 3, 8).is_err());
        assert!(nf2ff_2d_phasor(&a, &a, &short, nx, ny, dx, 3, 8).is_err());
        assert!(nf2ff_2d_phasor(&zero, &zero, &zero, nx, ny, dx, 3, 8).is_err());
        assert!(nf2ff_2d_phasor(&a, &a, &a, nx, ny, f64::NAN, 3, 8).is_err());
        assert!(nf2ff_2d_phasor(&a, &a, &a, nx, ny, dx, 10, 8).is_err());
        assert!(nf2ff_2d_phasor(&a, &a, &a, nx, ny, dx, 3, 0).is_err());
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

    // The two tests below previously asserted the swapped (wrong) assignment:
    // E ∥ interface sees the arithmetic mean, E ⊥ the harmonic mean
    // (Farjadpour et al. 2006).
    #[test]
    fn test_subpixel_perpendicular_less_than_parallel() {
        // Harmonic mean ≤ arithmetic mean (by AM-HM inequality)
        let par = subpixel_permittivity_parallel(1.0, 4.0, 0.5).unwrap();
        let perp = subpixel_permittivity_perpendicular(1.0, 4.0, 0.5);
        assert!(perp <= par);
    }

    #[test]
    fn test_subpixel_half_fill() {
        // At f=0.5: arithmetic = 2.5, harmonic = 8/5 = 1.6
        let par = subpixel_permittivity_parallel(1.0, 4.0, 0.5).unwrap();
        assert!((par - 2.5).abs() < 1e-10);
        let perp = subpixel_permittivity_perpendicular(1.0, 4.0, 0.5);
        assert!((perp - 1.6).abs() < 1e-10);
    }

    #[test]
    fn test_subpixel_matches_layered_capacitor() {
        // Independent reference: a parallel-plate capacitor filled with two
        // layers. Layers stacked across the gap (E ⊥ interfaces) are capacitors
        // in series; layers side by side (E ∥ interfaces) are in parallel.
        let (eps1, eps2, f) = (2.0, 9.0, 0.3);
        let (area, gap) = (2e-4, 1e-3);
        let c1 = EPSILON_0 * eps1 * area / ((1.0 - f) * gap);
        let c2 = EPSILON_0 * eps2 * area / (f * gap);
        let eps_perp_ref = (1.0 / (1.0 / c1 + 1.0 / c2)) * gap / (EPSILON_0 * area);
        let c_par = EPSILON_0 * eps1 * (1.0 - f) * area / gap + EPSILON_0 * eps2 * f * area / gap;
        let eps_par_ref = c_par * gap / (EPSILON_0 * area);
        let par = subpixel_permittivity_parallel(eps1, eps2, f).unwrap();
        let perp = subpixel_permittivity_perpendicular(eps1, eps2, f);
        assert!((par - eps_par_ref).abs() < 1e-12, "{par} vs {eps_par_ref}");
        assert!(
            (perp - eps_perp_ref).abs() < 1e-12,
            "{perp} vs {eps_perp_ref}"
        );
        let checked = subpixel_permittivity_perpendicular_checked(eps1, eps2, f).unwrap();
        assert!((checked - eps_perp_ref).abs() < 1e-12);
    }

    #[test]
    fn test_subpixel_fill_fraction_validation() {
        // f = 2 gave −2.0 and f = 4/3 gave inf before validation.
        for f in [2.0, 4.0 / 3.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(
                subpixel_permittivity_parallel(1.0, 4.0, f).is_err(),
                "f = {f}"
            );
            assert!(subpixel_permittivity_perpendicular_checked(1.0, 4.0, f).is_err());
        }
        assert!(subpixel_permittivity_parallel(f64::NAN, 4.0, 0.5).is_err());
        assert!(subpixel_permittivity_perpendicular_checked(1.0, -4.0, 0.5).is_err());
        // Unchecked variant: clamps f, never inf/negative for valid ε.
        assert!((subpixel_permittivity_perpendicular(1.0, 4.0, 4.0 / 3.0) - 4.0).abs() < 1e-12);
        assert!((subpixel_permittivity_perpendicular(1.0, 4.0, -3.0) - 1.0).abs() < 1e-12);
        assert!((subpixel_permittivity_perpendicular(1.0, 4.0, f64::NAN) - 1.0).abs() < 1e-12);
        assert!(subpixel_permittivity_perpendicular(0.0, 4.0, 0.5).is_nan());
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
        // boundary_x = 0 has no H^inc node at −½ and used to put the source on
        // node 0, where the boundary condition overwrote it (field stayed 0).
        // The checked constructor rejects it; the legacy one clamps to 1.
        assert!(TfsfSource::try_new(100, 0, 1e-4, 1e-13).is_err());
        let mut tfsf = TfsfSource::new(100, 0, 1e-4, 1e-13);
        assert_eq!(tfsf.boundary_x, 1);
        assert_eq!(tfsf.source_cell, 1);
        tfsf.step(1.0);
        assert!(tfsf.incident_e().abs() > 0.0);
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

    #[test]
    fn test_nonuniform_grid_validation() {
        assert!(NonUniformGrid::try_uniform(10, 10, 10, 1e-3).is_ok());
        assert!(NonUniformGrid::try_uniform(0, 10, 10, 1e-3).is_err());
        assert!(NonUniformGrid::try_uniform(10, 10, 10, 0.0).is_err());
        assert!(NonUniformGrid::try_uniform(10, 10, 10, f64::NAN).is_err());
        // Empty axis: min_spacing = ∞, so courant_dt used to be ∞.
        let empty = NonUniformGrid::uniform(0, 10, 10, 1e-3);
        assert!(empty.try_courant_dt().is_err());
        assert_eq!(empty.courant_dt(), 0.0);
        // Non-positive spacing: courant_dt used to be ≤ 0.
        let neg = NonUniformGrid::uniform(10, 10, 10, -1e-3);
        assert!(neg.validate().is_err());
        assert_eq!(neg.courant_dt(), 0.0);
        let mut g = NonUniformGrid::uniform(4, 4, 4, 1e-3);
        g.dy[2] = f64::NAN;
        assert!(g.try_courant_dt().is_err());
        assert!(NonUniformGrid::graded(10, f64::NAN, 1e-3).is_err());
        assert!(NonUniformGrid::graded(10, 1e-4, f64::INFINITY).is_err());
    }

    #[test]
    fn test_nonuniform_grid_courant_dt_value() {
        // Half the 3D limit d_min/(√3·c), d_min being the smallest cell anywhere.
        let mut g = NonUniformGrid::uniform(5, 5, 5, 1e-3);
        g.dz[3] = 2.5e-4;
        let dt = g.try_courant_dt().unwrap();
        let expected = 2.5e-4 / (2.0 * 3f64.sqrt() * SPEED_OF_LIGHT);
        assert!((dt - expected).abs() < 1e-12 * expected);
        assert!((g.courant_dt() - expected).abs() < 1e-12 * expected);
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
        // For f(x) = x², f'(x) = 2x. Output i is the derivative at the half
        // node x_{i+½} = (i + ½)dx, where the stencil is exact for quadratics
        // (the old test compared deriv[10] with f'(x₁₀) = 2.0 at tolerance 0.5).
        let dx = 0.1;
        let field: Vec<f64> = (0..20).map(|i| (i as f64 * dx).powi(2)).collect();
        let deriv = fourth_order_derivative(&field, dx);
        for (i, &d) in deriv.iter().enumerate().take(19) {
            let exact = 2.0 * (i as f64 + 0.5) * dx;
            assert!(
                (d - exact).abs() < 1e-12,
                "deriv[{i}] = {d}, expected {exact}"
            );
        }
        assert_eq!(deriv[19], 0.0); // padding
    }

    #[test]
    fn test_fourth_order_derivative_linear_all_points() {
        // Regression: the right-edge fallback was written to deriv[n−1]
        // instead of deriv[n−2] (f = x, n = 10 gave [1, …, 1, 0, 1]).
        let dx = 0.1;
        let field: Vec<f64> = (0..10).map(|i| i as f64 * dx).collect();
        let deriv = fourth_order_derivative(&field, dx);
        for &d in &deriv[..9] {
            assert!((d - 1.0).abs() < 1e-12, "{deriv:?}");
        }
        assert_eq!(deriv[9], 0.0);
    }

    #[test]
    fn test_fourth_order_derivative_cubic_interior_exact() {
        // The 4th-order staggered stencil is exact for cubics at x_{i+½}:
        // f = x³ → f'(x_{i+½}) = 3((i+½)dx)².
        let dx = 0.05;
        let n = 16;
        let field: Vec<f64> = (0..n).map(|i| (i as f64 * dx).powi(3)).collect();
        let deriv = fourth_order_derivative(&field, dx);
        for (i, &d) in deriv.iter().enumerate().take(n - 2).skip(1) {
            let exact = 3.0 * ((i as f64 + 0.5) * dx).powi(2);
            assert!((d - exact).abs() < 1e-12, "deriv[{i}] = {d} vs {exact}");
        }
        // Edge midpoints are 2nd order: (f[1] − f[0])/dx = dx² while the exact
        // value is 3dx²/4, an error of dx²/4.
        let edge_exact = 3.0 * (0.5 * dx).powi(2);
        assert!((deriv[0] - edge_exact - dx * dx / 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_fourth_order_derivative_small_and_invalid() {
        assert!(fourth_order_derivative(&[], 0.1).is_empty());
        assert_eq!(fourth_order_derivative(&[3.0], 0.1), vec![0.0]);
        let d2 = fourth_order_derivative(&[1.0, 2.0], 0.5);
        assert_eq!(d2, vec![2.0, 0.0]);
        let d3 = fourth_order_derivative(&[0.0, 1.0, 4.0], 1.0);
        assert_eq!(d3, vec![1.0, 3.0, 0.0]);
        assert!(
            fourth_order_derivative(&[0.0, 1.0, 2.0, 3.0], 0.0)
                .iter()
                .all(|v| *v == 0.0)
        );
        assert!(
            fourth_order_derivative(&[0.0, 1.0, 2.0, 3.0], f64::NAN)
                .iter()
                .all(|v| *v == 0.0)
        );
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

    #[test]
    fn test_cw_accumulator_exact_phasor() {
        // DFT orthogonality over whole periods: Σ cos(ωtₙ)e^{−jωtₙ}Δt = NΔt/2 and
        // Σ sin(ωtₙ)e^{−jωtₙ}Δt = −jNΔt/2 (sin has phasor −j in e^{jωt}).
        let freq = 1e9;
        let ppp = 100;
        let dt = 1.0 / (freq * ppp as f64);
        let n = 10 * ppp;
        let mut acc = CwAccumulator::try_new(2, freq, dt).unwrap();
        for step in 0..n {
            let th = 2.0 * std::f64::consts::PI * freq * step as f64 * dt;
            acc.accumulate(&[th.cos(), th.sin()], step);
        }
        let half = n as f64 * dt / 2.0;
        let tol = 1e-9 * half;
        assert!((acc.real[0] - half).abs() < tol && acc.imag[0].abs() < tol);
        assert!(acc.real[1].abs() < tol && (acc.imag[1] + half).abs() < tol);
        assert!((acc.amplitude()[0] - half).abs() < tol);
        assert!((acc.phase()[1] + std::f64::consts::FRAC_PI_2).abs() < 1e-9);
        assert_eq!(acc.n_samples, n);
    }

    #[test]
    fn test_cw_accumulator_length_mismatch() {
        let mut acc = CwAccumulator::new(3, 1e9, 1e-12);
        acc.accumulate(&[1.0; 5], 0); // longer field used to panic
        assert_eq!(acc.real, vec![1e-12; 3]);
        acc.accumulate(&[1.0], 0); // shorter: common prefix only
        assert!((acc.real[0] - 2e-12).abs() < 1e-24);
        assert!((acc.real[1] - 1e-12).abs() < 1e-24);
        let before = acc.clone();
        assert!(acc.try_accumulate(&[1.0; 4], 1).is_err());
        assert_eq!(acc.real, before.real);
        assert_eq!(acc.n_samples, before.n_samples);
        assert!(acc.try_accumulate(&[1.0; 3], 1).is_ok());
        assert_eq!(acc.n_samples, before.n_samples + 1);
    }

    #[test]
    fn test_cw_accumulator_try_new_validation() {
        assert!(CwAccumulator::try_new(4, 1e9, 1e-12).is_ok());
        assert!(CwAccumulator::try_new(4, f64::NAN, 1e-12).is_err());
        assert!(CwAccumulator::try_new(4, -1e9, 1e-12).is_err());
        assert!(CwAccumulator::try_new(4, 1e9, 0.0).is_err());
        assert!(CwAccumulator::try_new(4, 1e9, f64::INFINITY).is_err());
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

    #[test]
    fn test_eigenmode_two_mode_quality_factors() {
        // s = Σ sin(2πf_k t)·e^{−t/τ}: each mode is a Lorentzian line whose
        // half-power bandwidth is Δf = 1/(πτ), so Q = πf₀τ (47.1 and 110.0 for
        // 3 and 7 GHz at τ = 5 ns). The old whole-spectrum scan gave 0.37/0.87.
        let tau = 5e-9;
        let dt = 2e-11;
        let n = 5000; // 100 ns = 20τ
        let two_pi = 2.0 * std::f64::consts::PI;
        let signal: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                ((two_pi * 3e9 * t).sin() + (two_pi * 7e9 * t).sin()) * (-t / tau).exp()
            })
            .collect();
        let res = find_eigenmodes(&signal, dt, 1e9, 10e9, 1801).unwrap(); // df = 5 MHz
        assert_eq!(res.frequencies.len(), 2, "{:?}", res.frequencies);
        for (i, f0) in [3e9, 7e9].into_iter().enumerate() {
            let q_ref = std::f64::consts::PI * f0 * tau;
            let (f, q) = (res.frequencies[i], res.quality_factors[i]);
            assert!((f - f0).abs() < 1e6, "f = {f}");
            assert!(
                (q - q_ref).abs() / q_ref < 0.01,
                "f0 = {f0}: Q = {q}, expected {q_ref}"
            );
        }
    }

    #[test]
    fn test_eigenmode_frequency_grid_includes_endpoints() {
        // Regression: df = (fmax − fmin)/n_freq never evaluated fmax. Over
        // [1, 3] GHz with 3 points the grid is {1, 2, 3} GHz, so a 2 GHz mode is
        // reported at exactly 2 GHz (the old grid was 1, 1.67, 2.33 GHz).
        let dt = 1e-11;
        let tau = 3e-9;
        let signal: Vec<f64> = (0..3000)
            .map(|i| {
                let t = i as f64 * dt;
                (2.0 * std::f64::consts::PI * 2e9 * t).sin() * (-t / tau).exp()
            })
            .collect();
        let res = find_eigenmodes(&signal, dt, 1e9, 3e9, 3).unwrap();
        assert_eq!(res.frequencies.len(), 1);
        assert!((res.frequencies[0] - 2e9).abs() < 1.0);
    }

    #[test]
    fn test_eigenmode_validation() {
        let s = [1.0, 0.5, 0.25];
        assert!(find_eigenmodes(&s, 1e-12, 1e9, 10e9, 1).is_err());
        assert!(find_eigenmodes(&s, 0.0, 1e9, 10e9, 10).is_err());
        assert!(find_eigenmodes(&s, f64::NAN, 1e9, 10e9, 10).is_err());
        assert!(find_eigenmodes(&s, 1e-12, f64::NAN, 10e9, 10).is_err());
        assert!(find_eigenmodes(&s, 1e-12, 1e9, f64::NAN, 10).is_err());
        assert!(find_eigenmodes(&s, 1e-12, 1e9, f64::INFINITY, 10).is_err());
        assert!(find_eigenmodes(&s, 1e-12, -1e9, 10e9, 10).is_err());
        assert!(find_eigenmodes(&[1.0, f64::NAN], 1e-12, 1e9, 10e9, 10).is_err());
    }

    #[test]
    fn test_eigenmode_dft_recurrence_matches_direct() {
        // The phasor recurrence must match direct trig evaluation over a record
        // spanning several reseed blocks.
        let dt = 1e-12;
        let signal: Vec<f64> = (0..10_000)
            .map(|i| (i as f64 * 0.37).sin() + 0.1 * (i as f64 * 0.011).cos())
            .collect();
        for f in [1e9, 5.123e9, 47e9] {
            let omega = 2.0 * std::f64::consts::PI * f;
            let (mut re, mut im) = (0.0, 0.0);
            for (i, &s) in signal.iter().enumerate() {
                let th = omega * i as f64 * dt;
                re += s * th.cos() * dt;
                im -= s * th.sin() * dt;
            }
            let direct = re.hypot(im);
            let fast = dft_magnitude(&signal, omega, dt);
            assert!((fast - direct).abs() <= 1e-9 * direct, "{fast} vs {direct}");
        }
    }
}
