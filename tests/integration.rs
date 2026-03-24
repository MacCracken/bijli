//! Integration tests for bijli.

use bijli::field::{self, FieldVector, COULOMB_K, EPSILON_0, MU_0, SPEED_OF_LIGHT};
use bijli::maxwell;
use bijli::charge::{self, PointCharge, ELEMENTARY_CHARGE, ELECTRON_MASS};
use bijli::wave;

#[test]
fn test_maxwell_speed_of_light() {
    let v = maxwell::wave_speed(EPSILON_0, MU_0).unwrap();
    assert!((v - SPEED_OF_LIGHT).abs() / SPEED_OF_LIGHT < 1e-6);
}

#[test]
fn test_coulomb_inverse_square() {
    let q = 1e-6;
    let e1 = field::electric_field_point_charge(q, [0.0; 3], [1.0, 0.0, 0.0]).unwrap();
    let e2 = field::electric_field_point_charge(q, [0.0; 3], [2.0, 0.0, 0.0]).unwrap();
    let ratio = e1.magnitude() / e2.magnitude();
    assert!((ratio - 4.0).abs() < 1e-6);
}

#[test]
fn test_energy_conservation_em_wave() {
    let e0 = 1000.0;
    let b0 = wave::b_from_e(e0);
    let u_e = field::electric_energy_density(e0);
    let u_b = field::magnetic_energy_density(b0);
    // In an EM wave, electric and magnetic energy densities are equal
    assert!((u_e - u_b).abs() / u_e < 1e-6);
}

#[test]
fn test_lorentz_force_circular_motion() {
    // Electron in uniform B field should have circular motion
    let v = FieldVector::new(1e6, 0.0, 0.0);
    let b = FieldVector::new(0.0, 0.0, 1.0);
    let e = FieldVector::zero();
    let f = charge::lorentz_force(-ELEMENTARY_CHARGE, &v, &e, &b);
    // Force perpendicular to velocity → circular
    assert!(f.x.abs() < 1e-25);
    assert!(f.y.abs() > 0.0);
}

#[test]
fn test_hydrogen_atom_energy() {
    // Coulomb energy at Bohr radius
    let bohr_radius = 5.29177e-11;
    let u = charge::coulomb_potential_energy(
        ELEMENTARY_CHARGE,
        -ELEMENTARY_CHARGE,
        bohr_radius,
    ).unwrap();
    // Should be approximately -4.36e-18 J (-27.2 eV)
    let ev = u / ELEMENTARY_CHARGE;
    assert!((ev - (-27.2)).abs() < 0.5);
}

#[test]
fn test_superposition_principle() {
    let charges = vec![
        (1e-6, [0.0, 0.0, 0.0]),
        (-1e-6, [2.0, 0.0, 0.0]),
    ];
    let e = field::electric_field_superposition(&charges, [1.0, 0.0, 0.0]).unwrap();
    // At midpoint of dipole, both fields point in same direction (+x)
    assert!(e.x > 0.0);
}

#[test]
fn test_poynting_vector_direction() {
    // E in y, B in z → S in x (right-hand rule)
    let e = FieldVector::new(0.0, 100.0, 0.0);
    let b = FieldVector::new(0.0, 0.0, 1e-6);
    let s = wave::poynting_vector(&e, &b);
    assert!(s.x > 0.0);
    assert!(s.y.abs() < 1e-10);
    assert!(s.z.abs() < 1e-10);
}

#[test]
fn test_refractive_index_snells_law() {
    // Glass: ε_r ≈ 2.25 → n ≈ 1.5
    let n = maxwell::refractive_index(2.25, 1.0);
    assert!((n - 1.5).abs() < 1e-10);
}
