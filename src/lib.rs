//! Bijli — Electromagnetism simulation for AGNOS
//!
//! Hindi: बिजली (bijli) — electricity, lightning
//!
//! Provides electric and magnetic field computation, Maxwell's equations,
//! charge dynamics, and electromagnetic wave propagation. Built on
//! [hisab](https://crates.io/crates/hisab) for PDE solvers and
//! [impetus](https://crates.io/crates/impetus) for particle interactions.
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

#[cfg(feature = "logging")]
pub mod logging;

#[cfg(feature = "ai")]
pub mod ai;

pub use error::BijliError;
