//! Finite-Difference Time-Domain (FDTD) solver for 1D EM wave propagation.
//!
//! Implements the Yee algorithm for solving Maxwell's equations on a
//! staggered grid. E and H fields are offset by half a cell in space
//! and half a time step in time.

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
}
