use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};

fn field_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("field");

    group.bench_function("electric_field_point", |b| {
        b.iter(|| {
            black_box(bijli::field::electric_field_point_charge(
                1e-6,
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ))
        })
    });

    group.bench_function("electric_potential_point", |b| {
        b.iter(|| {
            black_box(bijli::field::electric_potential_point_charge(
                1e-6,
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ))
        })
    });

    group.bench_function("magnetic_field_moving", |b| {
        b.iter(|| {
            black_box(bijli::field::magnetic_field_moving_charge(
                1e-6,
                [0.0, 0.0, 0.0],
                [1e6, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ))
        })
    });

    let charges: Vec<(f64, [f64; 3])> = (0..10).map(|i| (1e-6, [i as f64, 0.0, 0.0])).collect();
    group.bench_function("superposition_10", |b| {
        b.iter(|| {
            black_box(bijli::field::electric_field_superposition(
                &charges,
                [5.0, 1.0, 0.0],
            ))
        })
    });

    let v = bijli::field::FieldVector::new(3.0, 4.0, 5.0);
    group.bench_function("vector_magnitude", |b| b.iter(|| black_box(v.magnitude())));

    let a = bijli::field::FieldVector::new(1.0, 2.0, 3.0);
    let bv = bijli::field::FieldVector::new(4.0, 5.0, 6.0);
    group.bench_function("vector_cross", |b| b.iter(|| black_box(a.cross(&bv))));

    let dipole = bijli::field::FieldVector::new(0.0, 0.0, 1e-30);
    group.bench_function("electric_dipole_field", |b| {
        b.iter(|| {
            black_box(bijli::field::electric_dipole_field(
                &dipole,
                [0.0, 0.0, 1.0],
            ))
        })
    });

    group.bench_function("magnetic_dipole_field", |b| {
        let m = bijli::field::FieldVector::new(0.0, 0.0, 1.0);
        b.iter(|| black_box(bijli::field::magnetic_dipole_field(&m, [0.0, 0.0, 1.0])))
    });

    group.bench_function("charged_sphere", |b| {
        b.iter(|| black_box(bijli::field::electric_field_charged_sphere(1e-6, 2.0)))
    });

    group.bench_function("ring_axis", |b| {
        b.iter(|| black_box(bijli::field::electric_field_ring_axis(1e-6, 1.0, 0.5)))
    });

    group.bench_function("disk_axis", |b| {
        b.iter(|| black_box(bijli::field::electric_field_disk_axis(1e-6, 1.0, 0.5)))
    });

    group.bench_function("trace_field_line_10", |b| {
        b.iter(|| {
            black_box(bijli::field::trace_field_line(
                [0.1, 0.0, 0.0],
                0.01,
                10,
                |pos| bijli::field::electric_field_point_charge(1e-6, [0.0; 3], pos),
            ))
        })
    });

    group.finish();
}

fn maxwell_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("maxwell");

    group.bench_function("wave_speed", |b| {
        b.iter(|| {
            black_box(bijli::maxwell::wave_speed(
                bijli::field::EPSILON_0,
                bijli::field::MU_0,
            ))
        })
    });

    group.bench_function("impedance", |b| {
        b.iter(|| {
            black_box(bijli::maxwell::impedance(
                bijli::field::EPSILON_0,
                bijli::field::MU_0,
            ))
        })
    });

    group.bench_function("skin_depth", |b| {
        b.iter(|| {
            black_box(bijli::maxwell::skin_depth(
                2.0 * std::f64::consts::PI * 60.0,
                bijli::field::MU_0,
                5.96e7,
            ))
        })
    });

    group.finish();
}

fn charge_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("charge");

    let q1 = bijli::charge::PointCharge::new(1e-6, 1.0, [0.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let q2 = bijli::charge::PointCharge::new(1e-6, 1.0, [1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    group.bench_function("coulomb_force", |b| {
        b.iter(|| black_box(bijli::charge::coulomb_force(&q1, &q2)))
    });

    let v = bijli::field::FieldVector::new(1e6, 0.0, 0.0);
    let e = bijli::field::FieldVector::new(1000.0, 0.0, 0.0);
    let bf = bijli::field::FieldVector::new(0.0, 0.0, 1.0);
    group.bench_function("lorentz_force", |b| {
        b.iter(|| {
            black_box(bijli::charge::lorentz_force(
                bijli::charge::ELEMENTARY_CHARGE,
                &v,
                &e,
                &bf,
            ))
        })
    });

    group.finish();
}

fn wave_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("wave");

    let e = bijli::field::FieldVector::new(0.0, 1000.0, 0.0);
    let b = bijli::field::FieldVector::new(0.0, 0.0, 1e-6);
    group.bench_function("poynting_vector", |b_bench| {
        b_bench.iter(|| black_box(bijli::wave::poynting_vector(&e, &b)))
    });

    group.bench_function("plane_wave_intensity", |b| {
        b.iter(|| black_box(bijli::wave::plane_wave_intensity(1000.0)))
    });

    group.bench_function("plane_wave_e", |b| {
        b.iter(|| black_box(bijli::wave::plane_wave_e(100.0, 1e7, 0.5, 3e15, 1e-15, 0.0)))
    });

    group.bench_function("fresnel_rs", |b| {
        b.iter(|| black_box(bijli::wave::fresnel_rs(1.0, 0.5, 1.5, 0.33)))
    });

    group.bench_function("snell_refraction", |b| {
        b.iter(|| black_box(bijli::wave::snell_refraction_angle(1.0, 1.5, 0.5)))
    });

    group.bench_function("half_wave_dipole_pattern", |b| {
        b.iter(|| black_box(bijli::wave::half_wave_dipole_pattern(1.0)))
    });

    group.bench_function("rect_waveguide_cutoff", |b| {
        b.iter(|| {
            black_box(bijli::wave::rectangular_waveguide_cutoff(
                22.86e-3,
                10.16e-3,
                1,
                0,
                bijli::field::SPEED_OF_LIGHT,
            ))
        })
    });

    group.finish();
}

fn fdtd_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("fdtd");

    group.bench_function("step_100_cells", |b| {
        let mut sim = bijli::fdtd::Fdtd1d::new(100, 1e-3).unwrap();
        sim.add_source(50, 1.0);
        b.iter(|| {
            sim.step_once();
            black_box(sim.e_field[50])
        })
    });

    group.bench_function("step_1000_cells", |b| {
        let mut sim = bijli::fdtd::Fdtd1d::new(1000, 1e-3).unwrap();
        sim.add_source(500, 1.0);
        b.iter(|| {
            sim.step_once();
            black_box(sim.e_field[500])
        })
    });

    group.bench_function("run_100_steps_500_cells", |b| {
        b.iter(|| {
            let mut sim = bijli::fdtd::Fdtd1d::new(500, 1e-4).unwrap();
            sim.run(100, Some(250), 1e9, 1.0);
            black_box(sim.total_energy())
        })
    });

    group.finish();
}

fn circuit_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("circuit");

    group.bench_function("rc_charging", |b| {
        b.iter(|| black_box(bijli::circuit::rc_charging_voltage(10.0, 1e3, 1e-6, 5e-4)))
    });

    group.bench_function("rl_current_rise", |b| {
        b.iter(|| black_box(bijli::circuit::rl_current_rise(10.0, 1.0, 1e-3, 5e-4)))
    });

    group.bench_function("rlc_impedance", |b| {
        b.iter(|| black_box(bijli::circuit::rlc_impedance(10.0, 1e-3, 1e-6, 31_416.0)))
    });

    group.bench_function("resonant_frequency", |b| {
        b.iter(|| black_box(bijli::circuit::resonant_frequency(1e-3, 1e-6)))
    });

    group.bench_function("resistance_parallel", |b| {
        let r = vec![10.0, 20.0, 30.0, 40.0, 50.0];
        b.iter(|| black_box(bijli::circuit::resistance_parallel(&r)))
    });

    group.finish();
}

fn material_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("material");

    let e = bijli::field::FieldVector::new(1000.0, 0.0, 0.0);
    group.bench_function("polarization", |b| {
        b.iter(|| black_box(bijli::material::polarization(4.0, &e)))
    });

    group.bench_function("displacement_field", |b| {
        b.iter(|| black_box(bijli::material::displacement_field(4.0, &e)))
    });

    let bf = bijli::field::FieldVector::new(0.0, 0.0, 1.0);
    group.bench_function("h_field_from_b", |b| {
        b.iter(|| black_box(bijli::material::h_field_from_b(1000.0, &bf)))
    });

    group.bench_function("clausius_mossotti", |b| {
        b.iter(|| black_box(bijli::material::clausius_mossotti(1e-40, 1e28)))
    });

    group.bench_function("curie_weiss", |b| {
        b.iter(|| black_box(bijli::material::curie_weiss(1.0, 1000.0, 770.0)))
    });

    group.finish();
}

fn relativity_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("relativity");

    group.bench_function("lorentz_factor", |b| {
        b.iter(|| {
            black_box(bijli::relativity::lorentz_factor(
                0.9 * bijli::field::SPEED_OF_LIGHT,
            ))
        })
    });

    let e = bijli::field::FieldVector::new(1e3, 2e3, 3e3);
    let bf = bijli::field::FieldVector::new(1e-6, 2e-6, 3e-6);
    group.bench_function("em_tensor_from_fields", |b| {
        b.iter(|| black_box(bijli::relativity::EmTensor::from_fields(&e, &bf)))
    });

    group.bench_function("lorentz_transform_x", |b| {
        b.iter(|| {
            black_box(bijli::relativity::lorentz_transform_x(
                &e,
                &bf,
                0.5 * bijli::field::SPEED_OF_LIGHT,
            ))
        })
    });

    let r = bijli::field::FieldVector::new(1.0, 0.0, 0.0);
    let beta_vec = bijli::field::FieldVector::new(0.1, 0.0, 0.0);
    let beta_dot = bijli::field::FieldVector::new(0.0, 1e-5, 0.0);
    group.bench_function("lienard_wiechert_e", |b| {
        b.iter(|| {
            black_box(bijli::relativity::lienard_wiechert_e(
                1e-6, &r, &beta_vec, &beta_dot,
            ))
        })
    });

    group.bench_function("four_vector_boost", |b| {
        let fv = bijli::relativity::FourVector::new(bijli::field::SPEED_OF_LIGHT, 1.0, 2.0, 3.0);
        b.iter(|| black_box(fv.boost_x(0.5 * bijli::field::SPEED_OF_LIGHT)))
    });

    group.finish();
}

criterion_group!(
    benches,
    field_benchmarks,
    maxwell_benchmarks,
    charge_benchmarks,
    wave_benchmarks,
    fdtd_benchmarks,
    circuit_benchmarks,
    material_benchmarks,
    relativity_benchmarks,
);
criterion_main!(benches);
