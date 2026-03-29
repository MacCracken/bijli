//! Soorat integration — visualization data structures for electromagnetic analysis.
//!
//! Provides structured types that soorat can render: FDTD field slices,
//! field line traces, point charge positions, radiation patterns, and
//! vector field samples.

use serde::{Deserialize, Serialize};

// ── FDTD field slices ──────────────────────────────────────────────────────

/// A 2D slice of an FDTD field for heatmap rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldSlice2D {
    /// Field values at grid points (flattened row-major: `values[y * nx + x]`).
    pub values: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// World-space origin of the slice (min corner) `[x, y]` in metres.
    pub origin: [f64; 2],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Name of the field component (e.g. "Ez", "Hx").
    pub component: String,
    /// Time step index when this slice was captured.
    pub step: usize,
}

/// A 3D field volume for volumetric rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldSlice3D {
    /// Field values at grid points (flattened: `values[z * ny * nx + y * nx + x]`).
    pub values: Vec<f64>,
    /// Grid dimensions (nx, ny, nz).
    pub dimensions: [usize; 3],
    /// World-space origin `[x, y, z]` in metres.
    pub origin: [f64; 3],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Name of the field component.
    pub component: String,
    /// Time step index.
    pub step: usize,
}

#[cfg(feature = "fdtd")]
impl FieldSlice2D {
    /// Extract a field slice from a 2D FDTD simulation.
    ///
    /// `component`: which field to extract — "fz", "fx", or "fy".
    #[must_use]
    pub fn from_fdtd_2d(sim: &crate::fdtd::Fdtd2d, component: &str) -> Self {
        let data = match component {
            "fx" => &sim.fx,
            "fy" => &sim.fy,
            _ => &sim.fz,
        };
        let label = match component {
            "fx" => format!("{:?}_x", sim.mode),
            "fy" => format!("{:?}_y", sim.mode),
            _ => format!("{:?}_z", sim.mode),
        };
        Self {
            values: data.clone(),
            dimensions: [sim.nx, sim.ny],
            origin: [0.0, 0.0],
            spacing: sim.dx,
            component: label,
            step: sim.step,
        }
    }
}

#[cfg(feature = "fdtd")]
impl FieldSlice3D {
    /// Extract a field component from a 3D FDTD simulation.
    ///
    /// `component`: "ex", "ey", "ez", "hx", "hy", or "hz".
    #[must_use]
    pub fn from_fdtd_3d(sim: &crate::fdtd::Fdtd3d, component: &str) -> Self {
        let data = match component {
            "ex" => &sim.ex,
            "ey" => &sim.ey,
            "ez" => &sim.ez,
            "hx" => &sim.hx,
            "hy" => &sim.hy,
            "hz" => &sim.hz,
            _ => &sim.ez,
        };
        Self {
            values: data.clone(),
            dimensions: [sim.nx, sim.ny, sim.nz],
            origin: [0.0, 0.0, 0.0],
            spacing: sim.dx,
            component: component.to_string(),
            step: sim.step,
        }
    }
}

// ── Field line traces ──────────────────────────────────────────────────────

/// A collection of field line traces for polyline rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldLineVisualization {
    /// Each field line is a sequence of 3D points.
    pub lines: Vec<FieldLine>,
}

/// A single field line trace.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FieldLine {
    /// Points along the field line `[x, y, z]` in metres.
    pub points: Vec<[f64; 3]>,
    /// Field magnitude at each point (for color mapping).
    pub magnitudes: Vec<f64>,
}

// ── Point charge visualization ─────────────────────────────────────────────

/// Point charge data for particle rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChargeVisualization {
    /// Charges with position and properties.
    pub charges: Vec<ChargePoint>,
}

/// A single charge for rendering.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ChargePoint {
    /// World-space position `[x, y, z]` in metres.
    pub position: [f32; 3],
    /// Charge in coulombs (sign determines color: positive=red, negative=blue).
    pub charge: f32,
    /// Absolute charge magnitude for size scaling.
    pub magnitude: f32,
}

#[cfg(feature = "charge")]
impl ChargeVisualization {
    /// Create from a slice of bijli `PointCharge` values.
    #[must_use]
    pub fn from_point_charges(charges: &[crate::charge::PointCharge]) -> Self {
        Self {
            charges: charges
                .iter()
                .map(|pc| ChargePoint {
                    position: [
                        pc.position[0] as f32,
                        pc.position[1] as f32,
                        pc.position[2] as f32,
                    ],
                    charge: pc.charge as f32,
                    magnitude: pc.charge.abs() as f32,
                })
                .collect(),
        }
    }
}

// ── Radiation pattern ──────────────────────────────────────────────────────

/// Far-field radiation pattern for polar/balloon rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RadiationPattern {
    /// Observation angles in radians.
    pub angles: Vec<f64>,
    /// Pattern magnitude at each angle (linear scale).
    pub pattern: Vec<f64>,
    /// Maximum gain value (for normalization).
    pub max_gain: f64,
}

#[cfg(feature = "fdtd")]
impl RadiationPattern {
    /// Create from a bijli `Nf2ffResult`.
    #[must_use]
    pub fn from_nf2ff(result: &crate::fdtd::Nf2ffResult) -> Self {
        let max_gain = result.pattern.iter().cloned().fold(0.0_f64, f64::max);
        Self {
            angles: result.angles.clone(),
            pattern: result.pattern.clone(),
            max_gain,
        }
    }
}

// ── Vector field sampling ──────────────────────────────────────────────────

/// A regular grid of 3D vector field samples for arrow/streamline rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VectorFieldSample {
    /// Vector values `[vx, vy, vz]` at each grid point (flattened row-major).
    pub vectors: Vec<[f64; 3]>,
    /// Grid dimensions (nx, ny, nz). For 2D grids, nz = 1.
    pub dimensions: [usize; 3],
    /// World-space origin `[x, y, z]` in metres.
    pub origin: [f64; 3],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Maximum vector magnitude (for normalization).
    pub max_magnitude: f64,
}

impl VectorFieldSample {
    /// Create from parallel arrays of vector components on a regular grid.
    #[must_use]
    pub fn from_components(
        vx: &[f64],
        vy: &[f64],
        vz: &[f64],
        dimensions: [usize; 3],
        origin: [f64; 3],
        spacing: f64,
    ) -> Self {
        let count = vx.len().min(vy.len()).min(vz.len());
        let mut max_mag = 0.0_f64;
        let vectors: Vec<[f64; 3]> = (0..count)
            .map(|i| {
                let v = [vx[i], vy[i], vz[i]];
                let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if mag > max_mag {
                    max_mag = mag;
                }
                v
            })
            .collect();
        Self {
            vectors,
            dimensions,
            origin,
            spacing,
            max_magnitude: max_mag,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_slice_2d_serializes() {
        let slice = FieldSlice2D {
            values: vec![0.0; 4],
            dimensions: [2, 2],
            origin: [0.0, 0.0],
            spacing: 1.0,
            component: "Ez".to_string(),
            step: 0,
        };
        let json = serde_json::to_string(&slice);
        assert!(json.is_ok());
    }

    #[test]
    fn field_slice_3d_serializes() {
        let slice = FieldSlice3D {
            values: vec![0.0; 8],
            dimensions: [2, 2, 2],
            origin: [0.0; 3],
            spacing: 0.5,
            component: "ex".to_string(),
            step: 10,
        };
        let json = serde_json::to_string(&slice);
        assert!(json.is_ok());
    }

    #[test]
    fn field_line_visualization() {
        let viz = FieldLineVisualization {
            lines: vec![FieldLine {
                points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.5, 0.0]],
                magnitudes: vec![1.0, 0.8, 0.5],
            }],
        };
        assert_eq!(viz.lines.len(), 1);
        assert_eq!(viz.lines[0].points.len(), 3);
    }

    #[test]
    fn charge_visualization_manual() {
        let viz = ChargeVisualization {
            charges: vec![
                ChargePoint {
                    position: [0.0, 0.0, 0.0],
                    charge: 1.6e-19,
                    magnitude: 1.6e-19,
                },
                ChargePoint {
                    position: [1.0, 0.0, 0.0],
                    charge: -1.6e-19,
                    magnitude: 1.6e-19,
                },
            ],
        };
        assert_eq!(viz.charges.len(), 2);
        assert!(viz.charges[0].charge > 0.0);
        assert!(viz.charges[1].charge < 0.0);
    }

    #[cfg(feature = "charge")]
    #[test]
    fn charge_visualization_from_point_charges() {
        let charges = vec![
            crate::charge::PointCharge::proton([0.0, 0.0, 0.0]),
            crate::charge::PointCharge::electron([1e-10, 0.0, 0.0]),
        ];
        let viz = ChargeVisualization::from_point_charges(&charges);
        assert_eq!(viz.charges.len(), 2);
        assert!(viz.charges[0].charge > 0.0); // proton
        assert!(viz.charges[1].charge < 0.0); // electron
    }

    #[test]
    fn radiation_pattern_manual() {
        let pat = RadiationPattern {
            angles: vec![0.0, std::f64::consts::FRAC_PI_2, std::f64::consts::PI],
            pattern: vec![1.0, 0.5, 0.1],
            max_gain: 1.0,
        };
        assert_eq!(pat.angles.len(), 3);
        assert_eq!(pat.max_gain, 1.0);
    }

    #[test]
    fn vector_field_from_components() {
        let vx = vec![1.0, 0.0, 0.0, 1.0];
        let vy = vec![0.0, 1.0, 0.0, 0.0];
        let vz = vec![0.0, 0.0, 1.0, 0.0];
        let field = VectorFieldSample::from_components(&vx, &vy, &vz, [2, 2, 1], [0.0; 3], 1.0);
        assert_eq!(field.vectors.len(), 4);
        assert!((field.max_magnitude - 1.0).abs() < 0.001);
    }

    #[test]
    fn vector_field_empty() {
        let field = VectorFieldSample::from_components(&[], &[], &[], [0, 0, 0], [0.0; 3], 1.0);
        assert!(field.vectors.is_empty());
        assert_eq!(field.max_magnitude, 0.0);
    }

    #[cfg(feature = "fdtd")]
    #[test]
    fn field_slice_from_fdtd_2d() {
        let sim = crate::fdtd::Fdtd2d::new(10, 10, 0.01, crate::fdtd::Mode2d::Tm).unwrap();
        let slice = FieldSlice2D::from_fdtd_2d(&sim, "fz");
        assert_eq!(slice.dimensions, [10, 10]);
        assert_eq!(slice.values.len(), 100);
        assert!((slice.spacing - 0.01).abs() < 1e-10);
    }

    #[cfg(feature = "fdtd")]
    #[test]
    fn field_slice_from_fdtd_3d() {
        let sim = crate::fdtd::Fdtd3d::new(4, 4, 4, 0.01).unwrap();
        let slice = FieldSlice3D::from_fdtd_3d(&sim, "ez");
        assert_eq!(slice.dimensions, [4, 4, 4]);
        assert_eq!(slice.values.len(), 64);
    }
}
