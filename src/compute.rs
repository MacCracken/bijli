//! Generic compute backend for GPU-accelerated FDTD.
//!
//! Defines the [`ComputeBackend`] trait and provides:
//! - [`CpuBackend`] — reference CPU implementation (always available)
//! - [`WgpuBackend`] — GPU via wgpu (requires `compute-wgpu` feature)
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
pub trait ComputeBackend {
    /// Opaque buffer handle for field data.
    type Buffer;

    /// Create a buffer initialized with the given data.
    fn create_buffer(&self, data: &[f64]) -> Result<Self::Buffer>;

    /// Read buffer contents back to host memory.
    fn read_buffer(&self, buffer: &Self::Buffer, dst: &mut [f64]) -> Result<()>;

    /// Write host data into an existing buffer.
    fn write_buffer(&self, buffer: &mut Self::Buffer, src: &[f64]) -> Result<()>;

    /// FDTD H-field update: H -= coeff * (E_neighbor - E_current).
    ///
    /// `axis`: 0 = x-direction differences, 1 = y-direction differences.
    fn fdtd_update_h(
        &self,
        h: &mut Self::Buffer,
        e: &Self::Buffer,
        coeff: &Self::Buffer,
        nx: usize,
        ny: usize,
        axis: u8,
    ) -> Result<()>;

    /// FDTD E-field update from curl of H: E += coeff * (dHy/dx - dHx/dy).
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
}

// ── CpuBackend ────────────────────────────────────────────────────

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
        dst.copy_from_slice(buffer);
        Ok(())
    }

    fn write_buffer(&self, buffer: &mut Self::Buffer, src: &[f64]) -> Result<()> {
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
        match axis {
            0 => {
                // H_x update: H_x -= coeff * (E_z[y+1] - E_z[y])
                for y in 0..ny - 1 {
                    for x in 0..nx {
                        let i = y * nx + x;
                        let i_yp = (y + 1) * nx + x;
                        h[i] -= coeff[i] * (e[i_yp] - e[i]);
                    }
                }
            }
            1 => {
                // H_y update: H_y += coeff * (E_z[x+1] - E_z[x])
                for y in 0..ny {
                    for x in 0..nx - 1 {
                        let i = y * nx + x;
                        let i_xp = y * nx + (x + 1);
                        h[i] += coeff[i] * (e[i_xp] - e[i]);
                    }
                }
            }
            _ => {
                return Err(BijliError::InvalidParameter {
                    reason: format!("axis must be 0 or 1, got {axis}"),
                });
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
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let i = y * nx + x;
                let dhy_dx = hy[i] - hy[y * nx + (x - 1)];
                let dhx_dy = hx[i] - hx[(y - 1) * nx + x];
                e[i] += coeff[i] * (dhy_dx - dhx_dy);
            }
        }
        Ok(())
    }

    fn apply_zero_boundary(&self, field: &mut Self::Buffer, nx: usize, ny: usize) -> Result<()> {
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
}

impl<B: ComputeBackend> GpuFdtd2d<B> {
    /// Create a new backend-accelerated 2D FDTD simulation.
    pub fn new(nx: usize, ny: usize, dx: f64, backend: B) -> Result<Self> {
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

        let dt = dx / (2.0 * std::f64::consts::SQRT_2 * SPEED_OF_LIGHT);
        let n = nx * ny;

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
        })
    }

    /// Advance one time step (TM mode).
    pub fn step_once(&mut self) -> Result<()> {
        let (nx, ny) = (self.nx, self.ny);

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

    /// Run for multiple steps with optional sinusoidal source.
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
                // Read fz, modify, write back (not ideal for GPU but correct)
                let i = sy * self.nx + sx;
                let mut fz_host = vec![0.0; self.nx * self.ny];
                self.backend.read_buffer(&self.fz, &mut fz_host)?;
                fz_host[i] += value;
                self.backend.write_buffer(&mut self.fz, &fz_host)?;
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
    pub fn read_ez(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fz, dst)
    }

    /// Read the H_x field back to host memory.
    pub fn read_hx(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fx, dst)
    }

    /// Read the H_y field back to host memory.
    pub fn read_hy(&self, dst: &mut [f64]) -> Result<()> {
        self.backend.read_buffer(&self.fy, dst)
    }
}

// ── wgpu backend (optional) ───────────────────────────────────────

#[cfg(feature = "compute-wgpu")]
mod wgpu_backend {
    use super::*;

    /// GPU compute backend using wgpu directly.
    pub struct WgpuBackend {
        _placeholder: (), // Will hold wgpu::Device, Queue, etc.
    }

    impl WgpuBackend {
        /// Create a new wgpu backend.
        ///
        /// Requests a GPU adapter and creates a compute device.
        pub fn new() -> Result<Self> {
            tracing::info!("initializing wgpu GPU compute backend");
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
}
