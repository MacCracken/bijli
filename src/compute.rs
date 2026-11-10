//! Generic compute backend for GPU-accelerated FDTD.
//!
//! Defines the [`ComputeBackend`] trait and provides:
//! - [`CpuBackend`] — reference CPU implementation (always available)
//! - `WgpuBackend` — **scaffold only** (requires the `compute-wgpu` feature): it
//!   does not yet implement [`ComputeBackend`] and performs no GPU work
//!
//! Use [`GpuFdtd2d`] for a backend-agnostic 2D FDTD solver.

use crate::error::{BijliError, Result};
use crate::field::{EPSILON_0, MU_0, SPEED_OF_LIGHT};

// ── ComputeBackend trait ──────────────────────────────────────────

/// Trait for compute backends (CPU or GPU) that can run FDTD kernels.
///
/// Implementations manage buffer storage and execute the stencil
/// update operations. The trait is generic over buffer type, allowing
/// CPU (`Vec<f64>`) and GPU (device buffers) to coexist.
///
/// Fields are 2D, row-major (index = y·nx + x), using the TM-mode Yee layout of
/// [`crate::fdtd::Fdtd2d`]: E_z at (x, y), H_x at (x, y + ½), H_y at (x + ½, y).
/// Implementations should return an error (not panic) when buffer lengths do
/// not match `nx·ny`.
pub trait ComputeBackend {
    /// Opaque buffer handle for field data.
    type Buffer;

    /// Create a buffer initialized with the given data.
    fn create_buffer(&self, data: &[f64]) -> Result<Self::Buffer>;

    /// Read buffer contents back to host memory (`dst` must match the buffer length).
    fn read_buffer(&self, buffer: &Self::Buffer, dst: &mut [f64]) -> Result<()>;

    /// Write host data into an existing buffer (`src` must match the buffer length).
    fn write_buffer(&self, buffer: &mut Self::Buffer, src: &[f64]) -> Result<()>;

    /// FDTD H-field update from E_z for one in-plane component.
    ///
    /// - `axis == 0` — H_x: `H[x, y] −= coeff[x, y] · (E[x, y+1] − E[x, y])`
    ///   for y in 0..ny−1 (i.e. −coeff·∂_y E_z·dy);
    /// - `axis == 1` — H_y: `H[x, y] += coeff[x, y] · (E[x+1, y] − E[x, y])`
    ///   for x in 0..nx−1 (i.e. +coeff·∂_x E_z·dx).
    ///
    /// Any other `axis` is an error.
    fn fdtd_update_h(
        &self,
        h: &mut Self::Buffer,
        e: &Self::Buffer,
        coeff: &Self::Buffer,
        nx: usize,
        ny: usize,
        axis: u8,
    ) -> Result<()>;

    /// FDTD E-field update from curl of H on interior nodes:
    /// `E[x, y] += coeff · ((Hy[x, y] − Hy[x−1, y]) − (Hx[x, y] − Hx[x, y−1]))`.
    fn fdtd_update_e(
        &self,
        e: &mut Self::Buffer,
        hx: &Self::Buffer,
        hy: &Self::Buffer,
        coeff: &Self::Buffer,
        nx: usize,
        ny: usize,
    ) -> Result<()>;

    /// Apply zero boundary conditions on a 2D field.
    fn apply_zero_boundary(&self, field: &mut Self::Buffer, nx: usize, ny: usize) -> Result<()>;

    /// Ensure all queued operations have completed.
    fn synchronize(&self) -> Result<()>;

    /// Add `value` to element `index` of `buffer` (point source injection).
    ///
    /// `scratch` is host memory with the buffer's length. The default
    /// implementation round-trips the whole buffer through it
    /// ([`read_buffer`](Self::read_buffer), modify, [`write_buffer`](Self::write_buffer));
    /// backends that can update one element in place should override this and
    /// may ignore `scratch`.
    ///
    /// # Errors
    /// [`BijliError::InvalidParameter`] if `index` is out of range, plus any
    /// error from the buffer transfers.
    fn add_to_element(
        &self,
        buffer: &mut Self::Buffer,
        index: usize,
        value: f64,
        scratch: &mut [f64],
    ) -> Result<()> {
        self.read_buffer(buffer, scratch)?;
        let len = scratch.len();
        let slot = scratch
            .get_mut(index)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("element index {index} out of range for buffer of {len}"),
            })?;
        *slot += value;
        self.write_buffer(buffer, scratch)
    }
}

// ── CpuBackend ────────────────────────────────────────────────────

/// Check that `nx·ny` is non-zero and each buffer holds exactly `nx·ny` values.
fn check_grid(nx: usize, ny: usize, lens: &[usize]) -> Result<usize> {
    let n = nx
        .checked_mul(ny)
        .filter(|&n| n > 0)
        .ok_or_else(|| BijliError::InvalidParameter {
            reason: format!("invalid grid dimensions {nx}x{ny}"),
        })?;
    if let Some(&bad) = lens.iter().find(|&&len| len != n) {
        return Err(BijliError::InvalidParameter {
            reason: format!("buffer length {bad} does not match grid {nx}x{ny} = {n}"),
        });
    }
    Ok(n)
}

/// Check that two buffers have equal length for a host transfer.
fn check_transfer(buffer_len: usize, host_len: usize) -> Result<()> {
    if buffer_len == host_len {
        Ok(())
    } else {
        Err(BijliError::InvalidParameter {
            reason: format!(
                "host slice length {host_len} does not match buffer length {buffer_len}"
            ),
        })
    }
}

/// CPU reference implementation of [`ComputeBackend`].
///
/// Uses `Vec<f64>` as buffers. Suitable for testing, small grids,
/// and systems without GPU hardware.
#[derive(Debug, Default)]
pub struct CpuBackend;

impl ComputeBackend for CpuBackend {
    type Buffer = Vec<f64>;

    fn create_buffer(&self, data: &[f64]) -> Result<Self::Buffer> {
        Ok(data.to_vec())
    }

    fn read_buffer(&self, buffer: &Self::Buffer, dst: &mut [f64]) -> Result<()> {
        check_transfer(buffer.len(), dst.len())?;
        dst.copy_from_slice(buffer);
        Ok(())
    }

    fn write_buffer(&self, buffer: &mut Self::Buffer, src: &[f64]) -> Result<()> {
        check_transfer(buffer.len(), src.len())?;
        buffer.copy_from_slice(src);
        Ok(())
    }

    fn fdtd_update_h(
        &self,
        h: &mut Self::Buffer,
        e: &Self::Buffer,
        coeff: &Self::Buffer,
        nx: usize,
        ny: usize,
        axis: u8,
    ) -> Result<()> {
        if axis > 1 {
            return Err(BijliError::InvalidParameter {
                reason: format!("axis must be 0 or 1, got {axis}"),
            });
        }
        check_grid(nx, ny, &[h.len(), e.len(), coeff.len()])?;
        if axis == 0 {
            // H_x update: H_x -= coeff * (E_z[y+1] - E_z[y])
            for y in 0..ny - 1 {
                for x in 0..nx {
                    let i = y * nx + x;
                    let i_yp = (y + 1) * nx + x;
                    h[i] -= coeff[i] * (e[i_yp] - e[i]);
                }
            }
        } else {
            // H_y update: H_y += coeff * (E_z[x+1] - E_z[x])
            for y in 0..ny {
                for x in 0..nx - 1 {
                    let i = y * nx + x;
                    let i_xp = y * nx + (x + 1);
                    h[i] += coeff[i] * (e[i_xp] - e[i]);
                }
            }
        }
        Ok(())
    }

    fn fdtd_update_e(
        &self,
        e: &mut Self::Buffer,
        hx: &Self::Buffer,
        hy: &Self::Buffer,
        coeff: &Self::Buffer,
        nx: usize,
        ny: usize,
    ) -> Result<()> {
        check_grid(nx, ny, &[e.len(), hx.len(), hy.len(), coeff.len()])?;
        for y in 1..ny.saturating_sub(1) {
            for x in 1..nx.saturating_sub(1) {
                let i = y * nx + x;
                let dhy_dx = hy[i] - hy[y * nx + (x - 1)];
                let dhx_dy = hx[i] - hx[(y - 1) * nx + x];
                e[i] += coeff[i] * (dhy_dx - dhx_dy);
            }
        }
        Ok(())
    }

    fn apply_zero_boundary(&self, field: &mut Self::Buffer, nx: usize, ny: usize) -> Result<()> {
        check_grid(nx, ny, &[field.len()])?;
        for x in 0..nx {
            field[x] = 0.0;
            field[(ny - 1) * nx + x] = 0.0;
        }
        for y in 0..ny {
            field[y * nx] = 0.0;
            field[y * nx + (nx - 1)] = 0.0;
        }
        Ok(())
    }

    fn synchronize(&self) -> Result<()> {
        Ok(()) // CPU is synchronous
    }

    fn add_to_element(
        &self,
        buffer: &mut Self::Buffer,
        index: usize,
        value: f64,
        _scratch: &mut [f64],
    ) -> Result<()> {
        let len = buffer.len();
        let slot = buffer
            .get_mut(index)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("element index {index} out of range for buffer of {len}"),
            })?;
        *slot += value;
        Ok(())
    }
}

// ── GpuFdtd2d ─────────────────────────────────────────────────────

/// Backend-agnostic 2D FDTD solver.
///
/// Same physics as [`crate::fdtd::Fdtd2d`] but field storage and computation
/// are delegated to a [`ComputeBackend`], enabling CPU or GPU execution.
///
/// Currently implements TM mode (E_z, H_x, H_y).
pub struct GpuFdtd2d<B: ComputeBackend> {
    /// The compute backend.
    pub backend: B,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Spatial step (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Current step index.
    pub step: usize,
    /// E_z field buffer.
    fz: B::Buffer,
    /// H_x field buffer.
    fx: B::Buffer,
    /// H_y field buffer.
    fy: B::Buffer,
    /// E-field coefficient buffer.
    ce: B::Buffer,
    /// H-field coefficient buffer.
    ch: B::Buffer,
    /// Host scratch for [`ComputeBackend::add_to_element`] (allocated on first use).
    scratch: Vec<f64>,
}

impl<B: ComputeBackend> GpuFdtd2d<B> {
    /// Create a new backend-accelerated 2D FDTD simulation.
    ///
    /// The time step is half the 2D Courant limit, dt = dx/(2√2·c).
    ///
    /// # Errors
    /// [`BijliError::InsufficientResolution`] if `nx < 3` or `ny < 3`;
    /// [`BijliError::InvalidParameter`] if `dx` is not finite and positive or
    /// `nx·ny` overflows; any backend buffer-creation error.
    pub fn new(nx: usize, ny: usize, dx: f64, backend: B) -> Result<Self> {
        if nx < 3 || ny < 3 {
            return Err(BijliError::InsufficientResolution {
                cells: nx.min(ny),
                domain_size: nx.min(ny) as f64 * dx,
            });
        }
        if !(dx > 0.0 && dx.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("dx must be finite and positive, got {dx}"),
            });
        }
        let n = nx
            .checked_mul(ny)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("grid {nx}x{ny} is too large"),
            })?;

        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);

        let zeros = vec![0.0; n];
        let ce_vals = vec![dt / (EPSILON_0 * dx); n];
        let ch_vals = vec![dt / (MU_0 * dx); n];

        let fz = backend.create_buffer(&zeros)?;
        let fx = backend.create_buffer(&zeros)?;
        let fy = backend.create_buffer(&zeros)?;
        let ce = backend.create_buffer(&ce_vals)?;
        let ch = backend.create_buffer(&ch_vals)?;

        tracing::debug!(nx, ny, dx, dt, "GpuFdtd2d created");

        Ok(Self {
            backend,
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
            scratch: Vec::new(),
        })
    }

    /// Advance one time step (TM mode).
    ///
    /// # Errors
    /// Any error reported by the backend kernels.
    pub fn step_once(&mut self) -> Result<()> {
        let (nx, ny) = (self.nx, self.ny);
        tracing::trace!(step = self.step, "GpuFdtd2d step");

        // Update H_x from E_z (y-direction differences)
        self.backend
            .fdtd_update_h(&mut self.fx, &self.fz, &self.ch, nx, ny, 0)?;
        // Update H_y from E_z (x-direction differences)
        self.backend
            .fdtd_update_h(&mut self.fy, &self.fz, &self.ch, nx, ny, 1)?;
        // Update E_z from curl(H)
        self.backend
            .fdtd_update_e(&mut self.fz, &self.fx, &self.fy, &self.ce, nx, ny)?;
        // Zero boundaries
        self.backend.apply_zero_boundary(&mut self.fz, nx, ny)?;

        self.step += 1;
        Ok(())
    }

    /// Inject a soft (additive) source into E_z at grid point (x, y).
    ///
    /// Points outside the grid are ignored, like [`crate::fdtd::Fdtd2d::add_source`].
    ///
    /// # Errors
    /// Any error reported by the backend while updating the buffer.
    pub fn add_source(&mut self, x: usize, y: usize, value: f64) -> Result<()> {
        if x >= self.nx || y >= self.ny {
            return Ok(());
        }
        let n = self
            .nx
            .checked_mul(self.ny)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("grid {}x{} is too large", self.nx, self.ny),
            })?;
        if self.scratch.len() != n {
            self.scratch = vec![0.0; n];
        }
        self.backend
            .add_to_element(&mut self.fz, y * self.nx + x, value, &mut self.scratch)
    }

    /// Run for multiple steps with optional sinusoidal source.
    ///
    /// A source outside the grid is ignored.
    ///
    /// # Errors
    /// Any error reported by the backend.
    pub fn run(
        &mut self,
        steps: usize,
        source: Option<(usize, usize)>,
        frequency: f64,
        amplitude: f64,
    ) -> Result<()> {
        tracing::debug!(
            steps,
            ?source,
            frequency,
            amplitude,
            "GpuFdtd2d run started"
        );
        for _ in 0..steps {
            if let Some((sx, sy)) = source {
                let t = self.step as f64 * self.dt;
                let value = amplitude * (2.0 * std::f64::consts::PI * frequency * t).sin();
                self.add_source(sx, sy, value)?;
            }
            self.step_once()?;
        }
        self.backend.synchronize()
    }

    /// Current simulation time.
    #[inline]
    #[must_use]
    pub fn time(&self) -> f64 {
        self.step as f64 * self.dt
    }

    /// Read the E_z field back to host memory.
    ///
    /// # Errors
    /// Backend error if `dst.len() != nx·ny` (for [`CpuBackend`]:
    /// [`BijliError::InvalidParameter`]).
    pub fn read_ez(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fz, dst)
    }

    /// Read the H_x field back to host memory.
    ///
    /// # Errors
    /// Backend error if `dst.len() != nx·ny`.
    pub fn read_hx(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fx, dst)
    }

    /// Read the H_y field back to host memory.
    ///
    /// # Errors
    /// Backend error if `dst.len() != nx·ny`.
    pub fn read_hy(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fy, dst)
    }
}

// ── wgpu backend (optional) ───────────────────────────────────────

#[cfg(feature = "compute-wgpu")]
mod wgpu_backend {
    use super::*;

    /// **Scaffold — not functional yet.** Placeholder for a wgpu GPU backend.
    ///
    /// It does not implement [`ComputeBackend`] and holds no device, queue or
    /// buffers, so it cannot run simulations; use [`CpuBackend`] for now. The
    /// type exists so the `compute-wgpu` feature and API surface are stable
    /// while the GPU implementation is developed.
    pub struct WgpuBackend {
        _placeholder: (), // Will hold wgpu::Device, Queue, etc.
    }

    impl WgpuBackend {
        /// Create the placeholder backend.
        ///
        /// Performs **no** GPU initialization (no adapter or device is
        /// requested) and currently always succeeds.
        ///
        /// # Errors
        /// None today; reserved for adapter/device failures once implemented.
        pub fn new() -> Result<Self> {
            tracing::debug!("WgpuBackend scaffold created (no GPU initialization performed)");
            Ok(Self { _placeholder: () })
        }
    }

    // TODO: Implement ComputeBackend for WgpuBackend. The implementation will:
    // 1. Create wgpu::Buffer storage
    // 2. Compile WGSL compute shaders for FDTD updates
    // 3. Dispatch compute passes and sync via queue.submit()
}

#[cfg(feature = "compute-wgpu")]
pub use wgpu_backend::WgpuBackend;

// ── Tests ─────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_backend_create_read() {
        let backend = CpuBackend;
        let data = vec![1.0, 2.0, 3.0];
        let buf = backend.create_buffer(&data).unwrap();
        let mut dst = vec![0.0; 3];
        backend.read_buffer(&buf, &mut dst).unwrap();
        assert_eq!(dst, data);
    }

    #[test]
    fn test_cpu_backend_write() {
        let backend = CpuBackend;
        let mut buf = backend.create_buffer(&[0.0; 3]).unwrap();
        backend.write_buffer(&mut buf, &[4.0, 5.0, 6.0]).unwrap();
        assert_eq!(buf, vec![4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_cpu_backend_synchronize() {
        let backend = CpuBackend;
        assert!(backend.synchronize().is_ok());
    }

    #[test]
    fn test_gpu_fdtd2d_creation() {
        let sim = GpuFdtd2d::new(50, 50, 1e-3, CpuBackend).unwrap();
        assert_eq!(sim.nx, 50);
        assert_eq!(sim.ny, 50);
        assert!(sim.dt > 0.0);
    }

    #[test]
    fn test_gpu_fdtd2d_too_small() {
        assert!(GpuFdtd2d::new(2, 50, 1e-3, CpuBackend).is_err());
    }

    #[test]
    fn test_gpu_fdtd2d_step() {
        let mut sim = GpuFdtd2d::new(30, 30, 1e-4, CpuBackend).unwrap();
        // Inject source directly into buffer
        let n = 30 * 30;
        let mut fz = vec![0.0; n];
        fz[15 * 30 + 15] = 1.0;
        sim.backend.write_buffer(&mut sim.fz, &fz).unwrap();

        sim.step_once().unwrap();

        let mut fz_out = vec![0.0; n];
        sim.read_ez(&mut fz_out).unwrap();
        // Energy should have spread
        let has_spread = fz_out
            .iter()
            .enumerate()
            .any(|(i, &v)| i != 15 * 30 + 15 && v.abs() > 1e-15);
        assert!(has_spread);
    }

    #[test]
    fn test_gpu_fdtd2d_run() {
        let mut sim = GpuFdtd2d::new(40, 40, 1e-4, CpuBackend).unwrap();
        sim.run(20, Some((20, 20)), 1e9, 1.0).unwrap();

        let mut fz = vec![0.0; 40 * 40];
        sim.read_ez(&mut fz).unwrap();
        let max: f64 = fz.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(max > 0.0);
    }

    #[test]
    fn test_gpu_fdtd2d_matches_fdtd2d() {
        // CpuBackend GpuFdtd2d should produce same results as Fdtd2d
        let nx = 40;
        let ny = 40;
        let dx = 1e-4;
        let steps = 30;
        let src = (20, 20);
        let freq = 1e9;

        // Standard Fdtd2d
        let mut std_sim = crate::fdtd::Fdtd2d::new(nx, ny, dx, crate::fdtd::Mode2d::Tm).unwrap();
        std_sim.run(steps, Some(src), freq, 1.0);

        // GpuFdtd2d with CpuBackend
        let mut gpu_sim = GpuFdtd2d::new(nx, ny, dx, CpuBackend).unwrap();
        gpu_sim.run(steps, Some(src), freq, 1.0).unwrap();

        // Compare E_z fields
        let mut gpu_fz = vec![0.0; nx * ny];
        gpu_sim.read_ez(&mut gpu_fz).unwrap();

        let max_diff: f64 = std_sim
            .fz
            .iter()
            .zip(gpu_fz.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        assert!(
            max_diff < 1e-10,
            "GpuFdtd2d(CpuBackend) doesn't match Fdtd2d: max_diff={max_diff}"
        );
    }

    #[test]
    fn test_cpu_backend_invalid_axis() {
        let backend = CpuBackend;
        let mut h = backend.create_buffer(&[0.0; 9]).unwrap();
        let e = backend.create_buffer(&[0.0; 9]).unwrap();
        let coeff = backend.create_buffer(&[0.0; 9]).unwrap();
        assert!(backend.fdtd_update_h(&mut h, &e, &coeff, 3, 3, 2).is_err());
    }

    #[test]
    fn test_gpu_fdtd2d_time() {
        let mut sim = GpuFdtd2d::new(10, 10, 1e-3, CpuBackend).unwrap();
        assert!(sim.time().abs() < 1e-30);
        sim.step_once().unwrap();
        assert!((sim.time() - sim.dt).abs() < 1e-30);
    }

    #[test]
    fn test_cpu_backend_transfer_length_mismatch() {
        // copy_from_slice used to panic on mismatched lengths.
        let backend = CpuBackend;
        let mut buf = backend.create_buffer(&[1.0, 2.0, 3.0]).unwrap();
        let mut short = vec![0.0; 2];
        assert!(backend.read_buffer(&buf, &mut short).is_err());
        assert!(backend.write_buffer(&mut buf, &[1.0; 4]).is_err());
        assert_eq!(buf, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_gpu_fdtd2d_read_wrong_length_errors() {
        let sim = GpuFdtd2d::new(10, 10, 1e-3, CpuBackend).unwrap();
        let mut dst = vec![0.0; 99];
        assert!(sim.read_ez(&mut dst).is_err());
        assert!(sim.read_hx(&mut dst).is_err());
        assert!(sim.read_hy(&mut dst).is_err());
        let mut ok = vec![0.0; 100];
        assert!(sim.read_ez(&mut ok).is_ok());
    }

    #[test]
    fn test_cpu_backend_kernels_validate_dimensions() {
        let backend = CpuBackend;
        let mut h = vec![0.0; 9];
        let e = vec![0.0; 9];
        let c = vec![0.0; 9];
        // ny = 0 used to underflow `ny - 1`
        assert!(backend.fdtd_update_h(&mut h, &e, &c, 3, 0, 0).is_err());
        assert!(backend.apply_zero_boundary(&mut h, 3, 0).is_err());
        assert!(backend.fdtd_update_e(&mut h, &e, &e, &c, 0, 3).is_err());
        // nx·ny larger than the buffers used to index out of bounds
        assert!(backend.fdtd_update_h(&mut h, &e, &c, 4, 4, 1).is_err());
        assert!(backend.fdtd_update_e(&mut h, &e, &e, &c, 4, 4).is_err());
        assert!(backend.apply_zero_boundary(&mut h, 4, 4).is_err());
        let short = vec![0.0; 8];
        assert!(backend.fdtd_update_h(&mut h, &short, &c, 3, 3, 0).is_err());
        assert!(
            backend
                .fdtd_update_h(&mut h, &e, &c, usize::MAX, 2, 0)
                .is_err()
        );
    }

    #[test]
    fn test_cpu_backend_update_h_matches_documented_stencil() {
        // 3×3 grid, E_z = x + 10y, coeff = 2: by the trait doc
        // axis 0: H_x[x,y] −= 2·(E[x,y+1] − E[x,y]) = −20 for y < 2, 0 on y = 2
        // axis 1: H_y[x,y] += 2·(E[x+1,y] − E[x,y]) = +2 for x < 2, 0 on x = 2
        let backend = CpuBackend;
        let e: Vec<f64> = (0..9)
            .map(|i| (i % 3) as f64 + 10.0 * (i / 3) as f64)
            .collect();
        let c = vec![2.0; 9];
        let mut hx = vec![0.0; 9];
        let mut hy = vec![0.0; 9];
        backend.fdtd_update_h(&mut hx, &e, &c, 3, 3, 0).unwrap();
        backend.fdtd_update_h(&mut hy, &e, &c, 3, 3, 1).unwrap();
        for i in 0..9 {
            let (x, y) = (i % 3, i / 3);
            assert_eq!(hx[i], if y < 2 { -20.0 } else { 0.0 }, "hx at ({x},{y})");
            assert_eq!(hy[i], if x < 2 { 2.0 } else { 0.0 }, "hy at ({x},{y})");
        }
    }

    #[test]
    fn test_gpu_fdtd2d_new_rejects_non_finite_dx() {
        // NaN passed the old `dx <= 0.0` check.
        assert!(GpuFdtd2d::new(10, 10, f64::NAN, CpuBackend).is_err());
        assert!(GpuFdtd2d::new(10, 10, f64::INFINITY, CpuBackend).is_err());
        assert!(GpuFdtd2d::new(10, 10, -1.0, CpuBackend).is_err());
        assert!(GpuFdtd2d::new(usize::MAX, 3, 1e-3, CpuBackend).is_err());
    }

    #[test]
    fn test_gpu_fdtd2d_source_outside_grid_ignored() {
        // (nx, 5) used to wrap to (0, 6); (0, ny) indexed out of bounds.
        let (nx, ny) = (20, 20);
        for src in [(nx, 5), (0, ny), (usize::MAX, 0)] {
            let mut sim = GpuFdtd2d::new(nx, ny, 1e-4, CpuBackend).unwrap();
            sim.run(10, Some(src), 1e9, 1.0).unwrap();
            let mut ez = vec![0.0; nx * ny];
            sim.read_ez(&mut ez).unwrap();
            assert!(
                ez.iter().all(|&v| v == 0.0),
                "source {src:?} leaked into grid"
            );
        }
    }

    #[test]
    fn test_gpu_fdtd2d_add_source() {
        let mut sim = GpuFdtd2d::new(10, 10, 1e-3, CpuBackend).unwrap();
        sim.add_source(3, 4, 2.5).unwrap();
        sim.add_source(3, 4, 0.5).unwrap();
        sim.add_source(10, 4, 1.0).unwrap(); // outside: ignored
        let mut ez = vec![0.0; 100];
        sim.read_ez(&mut ez).unwrap();
        assert_eq!(ez[4 * 10 + 3], 3.0);
        assert_eq!(ez.iter().filter(|&&v| v != 0.0).count(), 1);
    }

    /// Backend that only implements the required methods (delegating to
    /// `CpuBackend`), to exercise the provided `add_to_element`.
    struct RoundTripBackend;

    impl ComputeBackend for RoundTripBackend {
        type Buffer = Vec<f64>;
        fn create_buffer(&self, data: &[f64]) -> Result<Self::Buffer> {
            CpuBackend.create_buffer(data)
        }
        fn read_buffer(&self, buffer: &Self::Buffer, dst: &mut [f64]) -> Result<()> {
            CpuBackend.read_buffer(buffer, dst)
        }
        fn write_buffer(&self, buffer: &mut Self::Buffer, src: &[f64]) -> Result<()> {
            CpuBackend.write_buffer(buffer, src)
        }
        fn fdtd_update_h(
            &self,
            h: &mut Self::Buffer,
            e: &Self::Buffer,
            coeff: &Self::Buffer,
            nx: usize,
            ny: usize,
            axis: u8,
        ) -> Result<()> {
            CpuBackend.fdtd_update_h(h, e, coeff, nx, ny, axis)
        }
        fn fdtd_update_e(
            &self,
            e: &mut Self::Buffer,
            hx: &Self::Buffer,
            hy: &Self::Buffer,
            coeff: &Self::Buffer,
            nx: usize,
            ny: usize,
        ) -> Result<()> {
            CpuBackend.fdtd_update_e(e, hx, hy, coeff, nx, ny)
        }
        fn apply_zero_boundary(
            &self,
            field: &mut Self::Buffer,
            nx: usize,
            ny: usize,
        ) -> Result<()> {
            CpuBackend.apply_zero_boundary(field, nx, ny)
        }
        fn synchronize(&self) -> Result<()> {
            Ok(())
        }
    }

    #[test]
    fn test_default_add_to_element_matches_cpu() {
        let (nx, ny) = (30, 30);
        let mut a = GpuFdtd2d::new(nx, ny, 1e-4, CpuBackend).unwrap();
        let mut b = GpuFdtd2d::new(nx, ny, 1e-4, RoundTripBackend).unwrap();
        a.run(25, Some((12, 17)), 2e9, 1.0).unwrap();
        b.run(25, Some((12, 17)), 2e9, 1.0).unwrap();
        let (mut ea, mut eb) = (vec![0.0; nx * ny], vec![0.0; nx * ny]);
        a.read_ez(&mut ea).unwrap();
        b.read_ez(&mut eb).unwrap();
        assert_eq!(ea, eb);
        assert!(ea.iter().any(|&v| v != 0.0));
        let mut buf = vec![0.0; 4];
        let mut scratch = vec![0.0; 4];
        assert!(
            RoundTripBackend
                .add_to_element(&mut buf, 4, 1.0, &mut scratch)
                .is_err()
        );
        assert!(
            CpuBackend
                .add_to_element(&mut buf, 4, 1.0, &mut scratch)
                .is_err()
        );
    }
}
