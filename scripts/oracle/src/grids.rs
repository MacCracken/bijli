//! Shared input grids.
//!
//! Every input here must be reproducible bit for bit on any machine, so the inputs of every
//! row are the same wherever the generator runs. Only correctly rounded or exact operations
//! are used: + − × ÷, `floor`, `str::parse::<f64>` (correctly rounded, in core) and bit
//! construction. **Never a transcendental std f64 method** — those call the platform's
//! glibc libm, which is what the `math` module *records*, not what a grid may depend on.
//! (tanmatra's grids used the pinned pure-Rust `libm` crate for `10^x`; bijli has no such
//! dependency, so log-uniform values are built from a binary exponent and a mantissa.)
//! Copied from impetus's oracle (itself kimiya's) unchanged in behaviour.

/// Special values every float function is probed with.
pub fn specials() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1), // smallest subnormal
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

/// 2^k exactly, for k in −1074 ..= 1023 (subnormal powers included).
pub fn pow2(k: i32) -> f64 {
    assert!((-1074..=1023).contains(&k), "pow2({k}) out of range");
    if k >= -1022 {
        let biased = u64::try_from(k + 1023).unwrap_or(0);
        f64::from_bits(biased << 52)
    } else {
        let shift = u32::try_from(k + 1074).unwrap_or(0);
        f64::from_bits(1_u64 << shift)
    }
}

/// Deterministic xorshift64* — fixed seeds keep the sweeps identical across runs.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    /// 64 random bits (a random double's bit pattern, say).
    #[allow(dead_code)] // first used by the M1 grids (random bit patterns)
    pub fn bits(&mut self) -> u64 {
        self.next_u64()
    }

    /// Uniform in [0, 1): the top 53 bits as a fraction.
    pub fn unit(&mut self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let top = (self.next_u64() >> 11) as f64;
        top * (1.0 / 9_007_199_254_740_992.0) // 2^-53
    }

    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform integer in [lo, hi].
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 {
        #[allow(clippy::cast_sign_loss)]
        let span = (hi - lo + 1) as u64;
        #[allow(clippy::cast_possible_wrap)]
        let off = (self.next_u64() % span) as i64;
        lo + off
    }

    /// Roughly 10^u, u uniform in [lo_exp10, hi_exp10): a binary exponent drawn uniformly
    /// over the matching binades, times a uniform mantissa in [1, 2). Log-uniform across
    /// binades, linear within one; built without a transcendental (module comment).
    pub fn log_uniform(&mut self, lo_exp10: f64, hi_exp10: f64) -> f64 {
        const LOG2_10: f64 = core::f64::consts::LOG2_10;
        #[allow(clippy::cast_possible_truncation)]
        let lo = (lo_exp10 * LOG2_10).floor() as i64;
        #[allow(clippy::cast_possible_truncation)]
        let hi = (hi_exp10 * LOG2_10).floor() as i64 - 1;
        let k = self.int(lo.max(-1074), hi.min(1023));
        let mant = 1.0 + self.unit();
        let k = i32::try_from(k).unwrap_or(0);
        // A subnormal 2^k times the mantissa rounds once, by IEEE rules: still deterministic.
        mant * pow2(k)
    }

    /// As `log_uniform`, with a random sign.
    pub fn signed_log_uniform(&mut self, lo_exp10: f64, hi_exp10: f64) -> f64 {
        let x = self.log_uniform(lo_exp10, hi_exp10);
        if self.next_u64() & 1 == 1 { -x } else { x }
    }
}

/// a·b exactly, as a double-double (a fused multiply-add is correctly rounded).
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

/// a + b exactly, as a double-double (Knuth's TwoSum).
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

/// 10^-30 to about 2^-104 relative: the double `1e-30` and its remainder, from 5^30 exactly (u128).
fn ten_m30() -> (f64, f64) {
    let hi = 1e-30_f64;
    let bits = hi.to_bits();
    let m = (bits & ((1_u64 << 52) - 1)) | (1_u64 << 52);
    let e = i32::try_from((bits >> 52) & 0x7ff).unwrap_or(0) - 1075; // hi = m·2^e
    let num = u128::from(m) * 5_u128.pow(30); // hi·10^30 = num·2^(e + 30)
    let sh = -(e + 30);
    let one = 1_u128 << sh;
    let d = i128::try_from(num).unwrap_or(0) - i128::try_from(one).unwrap_or(0); // hi·10^30 − 1 = d·2^-sh
    #[allow(clippy::cast_precision_loss)]
    let lo = -(d as f64) * pow2(-sh) * hi;
    (hi, lo)
}

/// |a − b|² as a double-double: each difference exactly (TwoSum), its square to about 2^-104.
fn dist_sq_dd(a: [f64; 3], b: [f64; 3]) -> (f64, f64) {
    let mut acc = (0.0, 0.0);
    for i in 0..3 {
        let (s, t) = two_sum(a[i], -b[i]);
        let (p, pe) = two_prod(s, s);
        let lo = pe + 2.0 * s * t + t * t;
        let (h, l) = two_sum(acc.0, p);
        acc = two_sum(h, l + acc.1 + lo);
    }
    acc
}

/// Knife-edge points for the singular cutoff (bijli ADR 0002 row 17, M1 part 5d-1), as (source, field point) pairs:
/// field points along four directions at |d| = 10^-15 (`p60`: 10^-30, for R), stepped by ulps, against three sources
/// (the origin, which leaves the difference exact, and two small offsets, which make it round), kept where the exact
/// |d − b|² lies within an ulp of 10^-30 (two of 10^-60). Every row on which Rust's rounded test and the exact one can
/// disagree is among them.
pub fn knife_edge(p60: bool) -> Vec<([f64; 3], [f64; 3])> {
    let (th, tl) = ten_m30();
    let (target, tol, scale) = if p60 {
        let (p, pe) = two_prod(th, th);
        (two_sum(p, pe + 2.0 * th * tl), pow2(-251), 1e-30)
    } else {
        ((th, tl), pow2(-152), 1e-15)
    };
    let step = |x: f64, k: i64| {
        if x == 0.0 {
            x
        } else {
            f64::from_bits(x.to_bits().wrapping_add_signed(k))
        }
    };
    let dirs = [
        [0.6, 0.8, 0.0],
        [0.28, 0.96, 0.0],
        [0.48, 0.64, 0.6],
        [1.0, 0.0, 0.0],
    ];
    let sources = [
        [0.0, 0.0, 0.0],
        [-3e-32 * scale * 1e15, 0.0, 0.0],
        [pow2(-1070), 0.0, 0.0],
    ];
    let mut out = Vec::new();
    for dir in dirs {
        for k in -24_i64..=24 {
            for j in -6_i64..=6 {
                let a = [
                    step(dir[0] * scale, k),
                    step(dir[1] * scale, j),
                    dir[2] * scale,
                ];
                for b in sources {
                    let (h, l) = dist_sq_dd(a, b);
                    let diff = (h - target.0) + (l - target.1);
                    if diff.abs() <= tol {
                        out.push((b, a));
                    }
                }
            }
        }
    }
    out
}
