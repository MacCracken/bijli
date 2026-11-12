// Probe 4 (port plan §6 step 9): the Rust reference for p10_fdtd_kernel.cyr. bijli 1.2.0's row_diff, verbatim, over
// the same seeded rows, run REPS times (dst accumulates), then the FNV-1a digest of dst's bits.
//   rustc -O -o /tmp/p10ref scripts/probes/p10_kernel_ref.rs && /tmp/p10ref
#[inline]
fn row_diff<const ADD: bool>(dst: &mut [f64], k: &[f64], a: &[f64], b: &[f64]) {
    for (((d, &k), &a), &b) in dst.iter_mut().zip(k).zip(a).zip(b) {
        let v = k * (a - b);
        if ADD { *d += v; } else { *d -= v; }
    }
}
fn lcg(s: &mut u64) -> f64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*s >> 11) as f64) * (1.0 / 9007199254740992.0) * 2.0 - 1.0
}
fn main() {
    const N: usize = 10_000;
    const REPS: usize = 100;
    let mut s = 12345u64;
    let k: Vec<f64> = (0..N).map(|_| lcg(&mut s) * 0.5).collect();
    let a: Vec<f64> = (0..N).map(|_| lcg(&mut s)).collect();
    let b: Vec<f64> = (0..N).map(|_| lcg(&mut s)).collect();
    let mut d: Vec<f64> = (0..N).map(|_| lcg(&mut s)).collect();
    let t = std::time::Instant::now();
    for r in 0..REPS {
        if r % 2 == 0 { row_diff::<true>(&mut d, &k, &a, &b); } else { row_diff::<false>(&mut d, &k, &b, &a); }
    }
    let ns = t.elapsed().as_nanos() as f64 / (REPS * N) as f64;
    let mut h: u64 = 0xcbf29ce484222325;
    for x in &d { h = (h ^ x.to_bits()).wrapping_mul(0x100000001b3); }
    println!("rust row_diff digest {h:016x}  {ns:.3} ns/cell");
}
