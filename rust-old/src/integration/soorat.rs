//! Soorat integration — visualization data structures for electromagnetic analysis.
//!
//! Provides structured types that soorat can render: FDTD field slices,
//! field line traces, point charge positions, radiation patterns, and
//! vector field samples.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

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

/// Field axis of a 2D FDTD component.
#[cfg(feature = "fdtd")]
#[derive(Clone, Copy)]
enum Axis2d {
    X,
    Y,
    Z,
}

/// Resolve a 2D component name (case-insensitive) for the simulation mode.
///
/// Accepts the storage names `fx`/`fy`/`fz` and the physical names valid for
/// the mode (TM: `ez`, `hx`, `hy`; TE: `hz`, `ex`, `ey`).
#[cfg(feature = "fdtd")]
fn resolve_2d(mode: crate::fdtd::Mode2d, component: &str) -> Option<Axis2d> {
    use crate::fdtd::Mode2d;
    let name = component.trim().to_ascii_lowercase();
    match (mode, name.as_str()) {
        (_, "fx") | (Mode2d::Tm, "hx") | (Mode2d::Te, "ex") => Some(Axis2d::X),
        (_, "fy") | (Mode2d::Tm, "hy") | (Mode2d::Te, "ey") => Some(Axis2d::Y),
        (_, "fz") | (Mode2d::Tm, "ez") | (Mode2d::Te, "hz") => Some(Axis2d::Z),
        _ => None,
    }
}

#[cfg(feature = "fdtd")]
impl FieldSlice2D {
    fn build_2d(sim: &crate::fdtd::Fdtd2d, axis: Axis2d) -> Self {
        let (data, suffix) = match axis {
            Axis2d::X => (&sim.fx, "x"),
            Axis2d::Y => (&sim.fy, "y"),
            Axis2d::Z => (&sim.fz, "z"),
        };
        Self {
            values: data.clone(),
            dimensions: [sim.nx, sim.ny],
            origin: [0.0, 0.0],
            spacing: sim.dx,
            component: format!("{:?}_{suffix}", sim.mode),
            step: sim.step,
        }
    }

    /// Extract a field slice from a 2D FDTD simulation.
    ///
    /// `component` (case-insensitive): `"fz"`, `"fx"`, `"fy"`, or the
    /// physical name valid for the mode (TM: `"ez"`, `"hx"`, `"hy"`; TE:
    /// `"hz"`, `"ex"`, `"ey"`). Unrecognized names fall back to `fz`. The
    /// `component` label of the result always names the data actually
    /// returned, as `"<Mode>_<axis>"` (e.g. `"Tm_z"` for E_z in TM mode).
    /// Use [`FieldSlice2D::try_from_fdtd_2d`] to reject unknown names.
    #[must_use]
    pub fn from_fdtd_2d(sim: &crate::fdtd::Fdtd2d, component: &str) -> Self {
        let axis = resolve_2d(sim.mode, component).unwrap_or_else(|| {
            tracing::debug!(component, "unknown 2D FDTD component, falling back to fz");
            Axis2d::Z
        });
        Self::build_2d(sim, axis)
    }

    /// Extract a field slice from a 2D FDTD simulation, rejecting unknown
    /// component names.
    ///
    /// Accepts the same names as [`FieldSlice2D::from_fdtd_2d`].
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `component` does not name a field
    /// of the simulation's mode.
    pub fn try_from_fdtd_2d(sim: &crate::fdtd::Fdtd2d, component: &str) -> Result<Self> {
        let axis = resolve_2d(sim.mode, component).ok_or_else(|| BijliError::InvalidParameter {
            reason: format!(
                "unknown 2D FDTD component {component:?} for {:?} mode",
                sim.mode
            ),
        })?;
        tracing::trace!(component, "extracting 2D FDTD field slice");
        Ok(Self::build_2d(sim, axis))
    }
}

/// Canonical (lowercase) names of the 3D FDTD field components.
#[cfg(feature = "fdtd")]
const FDTD_3D_COMPONENTS: [&str; 6] = ["ex", "ey", "ez", "hx", "hy", "hz"];

#[cfg(feature = "fdtd")]
impl FieldSlice3D {
    /// Resolve a component name (case-insensitive) to its data and canonical
    /// lowercase label.
    fn select_3d<'a>(
        sim: &'a crate::fdtd::Fdtd3d,
        component: &str,
    ) -> Option<(&'a Vec<f64>, &'static str)> {
        let name = component.trim().to_ascii_lowercase();
        let label = FDTD_3D_COMPONENTS.iter().copied().find(|&c| c == name)?;
        let data = match label {
            "ex" => &sim.ex,
            "ey" => &sim.ey,
            "ez" => &sim.ez,
            "hx" => &sim.hx,
            "hy" => &sim.hy,
            _ => &sim.hz,
        };
        Some((data, label))
    }

    fn build_3d(sim: &crate::fdtd::Fdtd3d, data: &[f64], label: &str) -> Self {
        Self {
            values: data.to_vec(),
            dimensions: [sim.nx, sim.ny, sim.nz],
            origin: [0.0, 0.0, 0.0],
            spacing: sim.dx,
            component: label.to_string(),
            step: sim.step,
        }
    }

    /// Extract a field component from a 3D FDTD simulation.
    ///
    /// `component` (case-insensitive): `"ex"`, `"ey"`, `"ez"`, `"hx"`,
    /// `"hy"`, or `"hz"`. Unrecognized names fall back to `ez`. The
    /// `component` label of the result is the canonical lowercase name of
    /// the data actually returned (so `"Ex"` yields label `"ex"`, and an
    /// unknown name yields `"ez"`). Use [`FieldSlice3D::try_from_fdtd_3d`]
    /// to reject unknown names.
    #[must_use]
    pub fn from_fdtd_3d(sim: &crate::fdtd::Fdtd3d, component: &str) -> Self {
        let (data, label) = Self::select_3d(sim, component).unwrap_or_else(|| {
            tracing::debug!(component, "unknown 3D FDTD component, falling back to ez");
            (&sim.ez, "ez")
        });
        Self::build_3d(sim, data, label)
    }

    /// Extract a field component from a 3D FDTD simulation, rejecting
    /// unknown component names.
    ///
    /// Accepts the same (case-insensitive) names as
    /// [`FieldSlice3D::from_fdtd_3d`].
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `component` is not one of
    /// `ex`, `ey`, `ez`, `hx`, `hy`, `hz`.
    pub fn try_from_fdtd_3d(sim: &crate::fdtd::Fdtd3d, component: &str) -> Result<Self> {
        let (data, label) =
            Self::select_3d(sim, component).ok_or_else(|| BijliError::InvalidParameter {
                reason: format!(
                    "unknown 3D FDTD component {component:?} (expected one of {FDTD_3D_COMPONENTS:?})"
                ),
            })?;
        tracing::trace!(component = label, "extracting 3D FDTD field volume");
        Ok(Self::build_3d(sim, data, label))
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
    fn collect(vx: &[f64], vy: &[f64], vz: &[f64], count: usize) -> (Vec<[f64; 3]>, f64) {
        let mut max_mag = 0.0_f64;
        let vectors: Vec<[f64; 3]> = vx
            .iter()
            .zip(vy)
            .zip(vz)
            .take(count)
            .map(|((&x, &y), &z)| {
                let mag = (x * x + y * y + z * z).sqrt();
                if mag > max_mag {
                    max_mag = mag;
                }
                [x, y, z]
            })
            .collect();
        (vectors, max_mag)
    }

    /// Create from parallel arrays of vector components on a regular grid.
    ///
    /// Lenient: uses the first `min(vx.len(), vy.len(), vz.len())` samples
    /// and stores `dimensions` as given **without** checking that the sample
    /// count equals `nx·ny·nz`, so `vectors.len()` may disagree with
    /// `dimensions`. Prefer [`VectorFieldSample::try_from_components`], which
    /// rejects inconsistent input. Non-finite magnitudes are ignored when
    /// computing `max_magnitude`.
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
        let (vectors, max_magnitude) = Self::collect(vx, vy, vz, count);
        Self {
            vectors,
            dimensions,
            origin,
            spacing,
            max_magnitude,
        }
    }

    /// Create from parallel arrays of vector components on a regular grid,
    /// checking consistency.
    ///
    /// Samples are flattened row-major, `index = z·ny·nx + y·nx + x`.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if `vx`, `vy` and `vz` do not all
    /// have exactly `nx·ny·nz` elements (or that product overflows), if
    /// `spacing` is not finite and positive, or if `origin` is not finite.
    pub fn try_from_components(
        vx: &[f64],
        vy: &[f64],
        vz: &[f64],
        dimensions: [usize; 3],
        origin: [f64; 3],
        spacing: f64,
    ) -> Result<Self> {
        let [nx, ny, nz] = dimensions;
        let expected = nx
            .checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: format!("grid dimensions {dimensions:?} overflow usize"),
            })?;
        if vx.len() != expected || vy.len() != expected || vz.len() != expected {
            return Err(BijliError::InvalidParameter {
                reason: format!(
                    "component lengths ({}, {}, {}) must all equal nx·ny·nz = {expected}",
                    vx.len(),
                    vy.len(),
                    vz.len()
                ),
            });
        }
        if !(spacing.is_finite() && spacing > 0.0) {
            return Err(BijliError::InvalidParameter {
                reason: format!("grid spacing must be finite and positive, got {spacing}"),
            });
        }
        if !origin.iter().all(|v| v.is_finite()) {
            return Err(BijliError::InvalidParameter {
                reason: format!("grid origin must be finite, got {origin:?}"),
            });
        }
        tracing::trace!(?dimensions, "building vector field sample");
        let (vectors, max_magnitude) = Self::collect(vx, vy, vz, expected);
        Ok(Self {
            vectors,
            dimensions,
            origin,
            spacing,
            max_magnitude,
        })
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

    #[test]
    fn vector_field_try_from_components_checks_count() {
        let v4 = [1.0, 0.0, 0.0, 1.0];
        let v3 = [0.0, 1.0, 0.0];
        // Mismatched lengths: the lenient constructor truncated silently.
        assert!(
            VectorFieldSample::try_from_components(&v4, &v4, &v3, [2, 2, 1], [0.0; 3], 1.0)
                .is_err()
        );
        // Lengths agree with each other but not with nx·ny·nz.
        assert!(
            VectorFieldSample::try_from_components(&v4, &v4, &v4, [3, 2, 1], [0.0; 3], 1.0)
                .is_err()
        );
        assert!(
            VectorFieldSample::try_from_components(
                &v4,
                &v4,
                &v4,
                [usize::MAX, 2, 1],
                [0.0; 3],
                1.0
            )
            .is_err()
        );
        assert!(
            VectorFieldSample::try_from_components(&v4, &v4, &v4, [2, 2, 1], [0.0; 3], 0.0)
                .is_err()
        );
        let ok = VectorFieldSample::try_from_components(
            &[3.0, 0.0, 0.0, 0.0],
            &[4.0, 0.0, 0.0, 0.0],
            &[0.0, 0.0, 0.0, 2.0],
            [2, 2, 1],
            [0.0; 3],
            1.0,
        )
        .unwrap();
        assert_eq!(ok.vectors.len(), 4);
        assert_eq!(ok.vectors[0], [3.0, 4.0, 0.0]);
        assert!((ok.max_magnitude - 5.0).abs() < 1e-12);
    }

    #[test]
    fn vector_field_lenient_truncates_documented() {
        let field = VectorFieldSample::from_components(
            &[1.0, 2.0, 3.0],
            &[0.0, 0.0],
            &[0.0, 0.0, 0.0],
            [3, 1, 1],
            [0.0; 3],
            1.0,
        );
        assert_eq!(field.vectors.len(), 2);
        assert!((field.max_magnitude - 2.0).abs() < 1e-12);
    }

    #[cfg(feature = "fdtd")]
    #[test]
    fn field_slice_3d_component_is_case_insensitive_and_labelled() {
        let mut sim = crate::fdtd::Fdtd3d::new(4, 4, 4, 0.01).unwrap();
        sim.ex[0] = 1.0;
        sim.ez[0] = 2.0;
        sim.hy[0] = 3.0;
        // "Ex" used to return E_z data labelled "Ex".
        let slice = FieldSlice3D::from_fdtd_3d(&sim, "Ex");
        assert_eq!(slice.values[0], 1.0);
        assert_eq!(slice.component, "ex");
        let slice = FieldSlice3D::from_fdtd_3d(&sim, "HY");
        assert_eq!(slice.values[0], 3.0);
        assert_eq!(slice.component, "hy");
        // Unknown name: falls back to E_z and says so in the label.
        let slice = FieldSlice3D::from_fdtd_3d(&sim, "bogus");
        assert_eq!(slice.values[0], 2.0);
        assert_eq!(slice.component, "ez");
        assert!(FieldSlice3D::try_from_fdtd_3d(&sim, "bogus").is_err());
        let slice = FieldSlice3D::try_from_fdtd_3d(&sim, "Ex").unwrap();
        assert_eq!(slice.values[0], 1.0);
        assert_eq!(slice.component, "ex");
    }

    #[cfg(feature = "fdtd")]
    #[test]
    fn field_slice_2d_component_is_case_insensitive_and_mode_aware() {
        use crate::fdtd::{Fdtd2d, Mode2d};
        let mut tm = Fdtd2d::new(5, 5, 0.01, Mode2d::Tm).unwrap();
        tm.fx[0] = 1.0;
        tm.fz[0] = 3.0;
        let slice = FieldSlice2D::from_fdtd_2d(&tm, "FX");
        assert_eq!(slice.values[0], 1.0);
        assert_eq!(slice.component, "Tm_x");
        let slice = FieldSlice2D::try_from_fdtd_2d(&tm, "Ez").unwrap();
        assert_eq!(slice.values[0], 3.0);
        assert_eq!(slice.component, "Tm_z");
        let slice = FieldSlice2D::try_from_fdtd_2d(&tm, "Hx").unwrap();
        assert_eq!(slice.component, "Tm_x");
        // E_z is not a TE-mode field.
        let te = Fdtd2d::new(5, 5, 0.01, Mode2d::Te).unwrap();
        assert!(FieldSlice2D::try_from_fdtd_2d(&te, "Ez").is_err());
        assert_eq!(
            FieldSlice2D::try_from_fdtd_2d(&te, "ex").unwrap().component,
            "Te_x"
        );
        // Lenient fallback labels the data it returns.
        assert_eq!(FieldSlice2D::from_fdtd_2d(&te, "bogus").component, "Te_z");
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
