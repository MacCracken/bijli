use criterion::{Criterion, black_box, criterion_group, criterion_main};

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

    group.finish();
}

criterion_group!(
    benches,
    field_benchmarks,
    maxwell_benchmarks,
    charge_benchmarks,
    wave_benchmarks,
);
criterion_main!(benches);
