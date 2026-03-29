//! Bijli — Electromagnetism simulation for AGNOS
//!
//! Hindi: बिजली (bijli) — electricity, lightning
//!
//! Provides electric and magnetic field computation, Maxwell's equations,
//! charge dynamics, and electromagnetic wave propagation. Built on
//! [hisab](https://crates.io/crates/hisab) for numerical computing
//! (complex arithmetic, root-finding, transforms).
//!
//! # Modules
//!
//! - [`field`] — Electric and magnetic vector fields, field lines, potentials
//! - [`maxwell`] — Maxwell's equations, divergence, curl, wave equation
//! - [`charge`] — Point charges, Coulomb's law, Lorentz force, charge distributions
//! - [`wave`] — EM wave propagation, polarization, Poynting vector
//! - [`circuit`] — Ohm's law, capacitance, inductance, RC/RL/RLC circuits
//! - [`material`] — Dielectric and magnetic materials, polarization, magnetization
//! - [`fdtd`] — Finite-Difference Time-Domain solver for EM wave simulation
//! - [`relativity`] — Relativistic electrodynamics, EM tensor, Lorentz field transforms
//! - [`polarization`] — Jones vectors/matrices, Stokes parameters, Mueller matrices
//! - [`scattering`] — Mie scattering, Rayleigh scattering, cross-sections
//! - [`beam`] — Gaussian beams, ABCD matrices, resonator stability, HG/LG modes
//! - [`rf`] — Transmission lines, S-parameters, Smith chart, matching networks, Touchstone I/O
//! - [`error`] — Error types

pub mod error;

#[cfg(feature = "field")]
pub mod field;

#[cfg(feature = "maxwell")]
pub mod maxwell;

#[cfg(feature = "charge")]
pub mod charge;

#[cfg(feature = "wave")]
pub mod wave;

#[cfg(feature = "circuit")]
pub mod circuit;

#[cfg(feature = "material")]
pub mod material;

#[cfg(feature = "fdtd")]
pub mod fdtd;

#[cfg(feature = "relativity")]
pub mod relativity;

#[cfg(feature = "polarization")]
pub mod polarization;

#[cfg(feature = "scattering")]
pub mod scattering;

#[cfg(feature = "beam")]
pub mod beam;

#[cfg(feature = "rf")]
pub mod rf;

#[cfg(feature = "compute")]
pub mod compute;

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub use error::BijliError;
