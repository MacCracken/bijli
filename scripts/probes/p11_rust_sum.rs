// Probe 5 (port plan §6 step 9): Rust's float `Sum` — its start value and fold order. Run with
//   rustc -O -o /tmp/p11 scripts/probes/p11_rust_sum.rs && /tmp/p11
fn main() {
    let show = |name: &str, x: f64| println!("{name:28} {:016x}", x.to_bits());
    let empty: [f64; 0] = [];
    show("sum([])", empty.iter().sum::<f64>());
    show("sum([-0.0])", [-0.0f64].iter().sum::<f64>());
    show("sum([-0.0, -0.0])", [-0.0f64, -0.0].iter().sum::<f64>());
    show("sum([0.0])", [0.0f64].iter().sum::<f64>());
    show("sum([-0.0, 0.0])", [-0.0f64, 0.0].iter().sum::<f64>());
    // Order: a left fold. (1e16 + 1) + -1e16 = 0 left-to-right; a pairwise sum would give 1.
    show("sum([1e16, 1, -1e16])", [1e16f64, 1.0, -1e16].iter().sum::<f64>());
    show("fold(-0.0) same input", [1e16f64, 1.0, -1e16].iter().fold(-0.0, |a, b| a + b));
    let v: Vec<f64> = (0..1000).map(|k| 1.0 / f64::from(k + 1)).collect();
    show("sum(1/k, k=1..1000)", v.iter().sum::<f64>());
    show("fold(-0.0) same", v.iter().fold(-0.0, |a, b| a + b));
    show("fold(+0.0) same", v.iter().fold(0.0, |a, b| a + b));
}
