//! Basic electromagnetism with bijli.

fn main() {
    use bijli::field;
    use bijli::charge::{self, PointCharge};
    use bijli::maxwell;
    use bijli::wave;

    // Electric field from a point charge
    let e = field::electric_field_point_charge(
        1e-6,           // 1 μC
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ).unwrap();
    println!("E field at 1m from 1μC: {:.2} V/m", e.magnitude());

    // Coulomb force between two charges
    let q1 = PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]);
    let q2 = PointCharge::new(-1e-6, 1.0, [0.5, 0.0, 0.0], [0.0; 3]);
    let f = charge::coulomb_force(&q1, &q2).unwrap();
    println!("Coulomb force: {:.4} N", f.magnitude());

    // Speed of light from Maxwell's equations
    let c = maxwell::wave_speed(field::EPSILON_0, field::MU_0).unwrap();
    println!("Speed of light: {:.0} m/s", c);

    // EM wave intensity
    let intensity = wave::plane_wave_intensity(1000.0);
    println!("Plane wave intensity (E₀=1kV/m): {:.2} W/m²", intensity);
}
