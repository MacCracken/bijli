//! Integration tests for bijli.

use bijli::charge::{self, ELEMENTARY_CHARGE};
use bijli::circuit;
use bijli::fdtd;
use bijli::field::{self, EPSILON_0, FieldVector, MU_0, SPEED_OF_LIGHT};
use bijli::material;
use bijli::maxwell;
use bijli::relativity;
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
    let u = charge::coulomb_potential_energy(ELEMENTARY_CHARGE, -ELEMENTARY_CHARGE, bohr_radius)
        .unwrap();
    // Should be approximately -4.36e-18 J (-27.2 eV)
    let ev = u / ELEMENTARY_CHARGE;
    assert!((ev - (-27.2)).abs() < 0.5);
}

#[test]
fn test_superposition_principle() {
    let charges = vec![(1e-6, [0.0, 0.0, 0.0]), (-1e-6, [2.0, 0.0, 0.0])];
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
    let n = maxwell::refractive_index(2.25, 1.0).unwrap();
    assert!((n - 1.5).abs() < 1e-10);
}

// ── Circuit integration tests ──────────────────────────────────────

#[test]
fn test_rc_full_charge_cycle() {
    // Charge to ~63%, then discharge from there
    let v0 = 12.0;
    let r = 1e3;
    let c = 1e-6;
    let tau = circuit::rc_time_constant(r, c);

    let v_charged = circuit::rc_charging_voltage(v0, r, c, 5.0 * tau);
    // After 5τ, should be ~99.3% charged
    assert!((v_charged / v0 - 0.9933).abs() < 0.001);

    let v_discharged = circuit::rc_discharging_voltage(v_charged, r, c, 5.0 * tau);
    // After 5τ discharge, should be ~0.7% of charged value
    assert!(v_discharged / v_charged < 0.01);
}

#[test]
fn test_rlc_resonance_matches_components() {
    // At resonance, impedance equals resistance
    let r = 50.0;
    let l = 10e-3;
    let c = 100e-9;
    let f0 = circuit::resonant_frequency(l, c).unwrap();
    let omega0 = 2.0 * std::f64::consts::PI * f0;
    let z = circuit::rlc_impedance(r, l, c, omega0).unwrap();
    assert!((z - r).abs() / r < 1e-6);
}

#[test]
fn test_series_parallel_identity() {
    // Two equal resistors: series = 4 × parallel
    let r = 100.0;
    let series = circuit::resistance_series(&[r, r]).unwrap();
    let parallel = circuit::resistance_parallel(&[r, r]).unwrap();
    assert!((series / parallel - 4.0).abs() < 1e-10);
}

#[test]
fn test_capacitor_energy_from_field() {
    // Energy stored in a parallel-plate capacitor should equal
    // the energy density × volume between the plates
    let eps_r = 1.0;
    let area = 0.01; // 100 cm²
    let sep = 1e-3; // 1mm
    let v = 100.0; // 100V

    let cap = circuit::parallel_plate_capacitance(EPSILON_0 * eps_r, area, sep).unwrap();
    let u_circuit = circuit::capacitor_energy(cap, v);

    // E = V/d, energy density = ε₀E²/2, total = u × A × d
    let e_field = v / sep;
    let u_field = material::dielectric_energy_density(eps_r, e_field) * area * sep;

    assert!((u_circuit - u_field).abs() / u_circuit < 1e-6);
}

// ── Material integration tests ─────────────────────────────────────

#[test]
fn test_d_field_equals_eps_e() {
    // D = ε₀ε_r E, verify consistency
    let eps_r = 4.7; // FR-4
    let e = FieldVector::new(500.0, 0.0, 0.0);
    let d = material::displacement_field(eps_r, &e);
    let p = material::polarization(eps_r, &e);
    // D = ε₀E + P
    let d_check = e.scale(EPSILON_0) + p;
    assert!((d.x - d_check.x).abs() < 1e-20);
}

#[test]
fn test_b_h_m_consistency() {
    // B = μ₀(H + M) = μ₀μ_r H
    let mu_r = 5000.0; // soft iron
    let h = FieldVector::new(100.0, 0.0, 0.0);
    let b = material::b_field_from_h(mu_r, &h);
    let m = material::magnetization(mu_r, &h);
    let b_check = (h + m).scale(MU_0);
    assert!((b.x - b_check.x).abs() / b.x < 1e-10);
}

#[test]
fn test_magnetic_energy_in_material() {
    // u = B²/(2μ) should equal u = μH²/2
    let mu_r = 1000.0;
    let h_mag = 100.0;
    let b_mag = mu_r * MU_0 * h_mag;
    let u1 = material::magnetic_energy_density_material(mu_r, b_mag).unwrap();
    let u2 = 0.5 * mu_r * MU_0 * h_mag * h_mag;
    assert!((u1 - u2).abs() / u1 < 1e-6);
}

// ── Wave propagation integration tests ─────────────────────────────

#[test]
fn test_snell_law_reciprocity() {
    // Snell's law is symmetric: air→glass→air should return to original angle
    let theta_i = 0.5; // ~28.6°
    let theta_t = wave::snell_refraction_angle(1.0, 1.5, theta_i)
        .unwrap()
        .unwrap();
    let theta_back = wave::snell_refraction_angle(1.5, 1.0, theta_t)
        .unwrap()
        .unwrap();
    assert!((theta_back - theta_i).abs() < 1e-10);
}

#[test]
fn test_fresnel_energy_conservation() {
    // R + T = 1 for normal incidence (no absorption)
    let r = wave::reflectance_normal(1.0, 1.5).unwrap();
    let t = wave::transmittance_normal(1.0, 1.5).unwrap();
    assert!((r + t - 1.0).abs() < 1e-10);
}

#[test]
fn test_waveguide_cutoff_vs_wavelength() {
    // f_c × λ_c = c for TE₁₀ mode
    let a = 22.86e-3;
    let b = 10.16e-3;
    let fc = wave::rectangular_waveguide_cutoff(a, b, 1, 0, SPEED_OF_LIGHT).unwrap();
    let lc = wave::rectangular_waveguide_cutoff_wavelength(a, b, 1, 0).unwrap();
    assert!((fc * lc - SPEED_OF_LIGHT).abs() / SPEED_OF_LIGHT < 1e-6);
}

#[test]
fn test_antenna_aperture_directivity_relation() {
    // A_e × 4π/λ² = D for any antenna
    let lambda = 0.1; // 3 GHz
    let d = wave::half_wave_dipole_directivity();
    let ae = wave::effective_aperture(lambda, d).unwrap();
    let d_back = ae * 4.0 * std::f64::consts::PI / (lambda * lambda);
    assert!((d_back - d).abs() < 1e-6);
}

// ── FDTD integration tests ─────────────────────────────────────────

#[test]
fn test_fdtd_wave_propagation() {
    // A pulse should propagate outward from the source
    let mut sim = fdtd::Fdtd1d::new(200, 1e-4).unwrap();
    sim.add_source(100, 1.0);

    for _ in 0..50 {
        sim.step_once();
    }

    // Field should have spread beyond the source cell
    // The field should be nonzero away from source
    let has_spread = sim
        .e_field
        .iter()
        .enumerate()
        .any(|(i, &v)| i != 100 && v.abs() > 1e-10);
    assert!(has_spread);
}

#[test]
fn test_fdtd_dielectric_slows_wave() {
    let n = 400;
    let dx = 1e-4;
    let steps = 100;

    // Vacuum simulation
    let mut vac = fdtd::Fdtd1d::new(n, dx).unwrap();
    // Dielectric simulation (ε_r = 4 → v = c/2)
    let mut die = fdtd::Fdtd1d::new(n, dx).unwrap();
    die.set_permittivity(0, n, 4.0).unwrap();

    for _ in 0..steps {
        let t_vac = vac.step as f64 * vac.dt;
        let t_die = die.step as f64 * die.dt;
        vac.add_source(50, (2.0 * std::f64::consts::PI * 1e9 * t_vac).sin());
        die.add_source(50, (2.0 * std::f64::consts::PI * 1e9 * t_die).sin());
        vac.step_once();
        die.step_once();
    }

    // Find furthest cell with significant field
    let threshold = 0.001;
    let vac_front = vac
        .e_field
        .iter()
        .rposition(|&v| v.abs() > threshold)
        .unwrap_or(50);
    let die_front = die
        .e_field
        .iter()
        .rposition(|&v| v.abs() > threshold)
        .unwrap_or(50);

    assert!(
        die_front <= vac_front,
        "dielectric wave should travel slower"
    );
}

// ── Relativity integration tests ───────────────────────────────────

#[test]
fn test_lorentz_transform_preserves_both_invariants() {
    let e = FieldVector::new(1e3, 2e3, -500.0);
    let b = FieldVector::new(1e-5, -3e-5, 2e-5);
    let f1 = relativity::EmTensor::from_fields(&e, &b);

    let (ep, bp) = relativity::lorentz_transform_x(&e, &b, 0.7 * SPEED_OF_LIGHT).unwrap();
    let f2 = relativity::EmTensor::from_fields(&ep, &bp);

    // Both Lorentz invariants must be preserved
    assert!(
        (f1.first_invariant() - f2.first_invariant()).abs() / (f1.first_invariant().abs() + 1e-30)
            < 1e-6
    );
    assert!(
        (f1.second_invariant() - f2.second_invariant()).abs()
            / (f1.second_invariant().abs() + 1e-30)
            < 1e-6
    );
}

#[test]
fn test_lw_velocity_field_matches_coulomb_at_rest() {
    // A static charge's LW field should match the Coulomb field exactly
    let q = 1e-6;
    let r_vec = FieldVector::new(0.0, 2.0, 0.0);
    let e_lw =
        relativity::lienard_wiechert_e(q, &r_vec, &FieldVector::zero(), &FieldVector::zero())
            .unwrap();
    let e_coulomb = field::electric_field_point_charge(q, [0.0; 3], [0.0, 2.0, 0.0]).unwrap();
    assert!((e_lw.x - e_coulomb.x).abs() < 1e-3);
    assert!((e_lw.y - e_coulomb.y).abs() / e_coulomb.y.abs() < 1e-6);
    assert!((e_lw.z - e_coulomb.z).abs() < 1e-3);
}

#[test]
fn test_four_vector_interval_lorentz_invariant() {
    // Spacetime interval must be preserved under boost
    let event = relativity::FourVector::from_time_position(1e-9, [0.3, 0.0, 0.0]);
    let s2_lab = event.interval_sq();
    let boosted = event.boost_x(0.9 * SPEED_OF_LIGHT).unwrap();
    let s2_boosted = boosted.interval_sq();
    assert!((s2_lab - s2_boosted).abs() / s2_lab.abs() < 1e-6);
}

#[test]
fn test_relativistic_larmor_matches_nonrelativistic_at_low_v() {
    let q = ELEMENTARY_CHARGE;
    let a = 1e15;
    let v = FieldVector::new(1e3, 0.0, 0.0); // v << c
    let acc = FieldVector::new(0.0, a, 0.0);
    let p_rel = relativity::relativistic_larmor_power(q, &v, &acc).unwrap();
    let p_nr = relativity::larmor_power(q, a);
    // At v << c, should be nearly equal
    assert!((p_rel - p_nr).abs() / p_nr < 0.01);
}

#[test]
fn test_pure_magnetic_field_transforms_to_electric() {
    // A pure B field, when boosted, gains an E field
    let e = FieldVector::zero();
    let b = FieldVector::new(0.0, 0.0, 1.0); // 1T in z
    let v = 0.5 * SPEED_OF_LIGHT;
    let (ep, _) = relativity::lorentz_transform_x(&e, &b, v).unwrap();
    // Should create E_y = −γvB_z
    let gamma = relativity::lorentz_factor(v).unwrap();
    let expected_ey = -gamma * v * b.z;
    assert!((ep.y - expected_ey).abs() / expected_ey.abs() < 1e-6);
}
