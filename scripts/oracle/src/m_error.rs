//! `error.tsv` — every `BijliError` variant's Display text (`e.to_string()`), over its payload's grid.
//!
//! The port returns a `BJ_ERR_*` code and keeps the variant's fields and this text in the calling thread's
//! last-error record (ADR 0002 row 4), so each row's outputs are the variant, its fields and the text: class A on
//! the code and on the payload. The f64 fields run over the specials, a sweep of binades, the powers of ten and
//! seeded bit patterns, which also pins Rust's f64 Display inside a message; the strings include UTF-8, `%`,
//! spaces and the empty string. `InvalidParameter/f64` rows are a reason built as bijli's are,
//! `format!("{prefix}{x}")`, so the port's formatted reasons are pinned too.

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, err, f, i, s};
use bijli::BijliError;

fn display(e: &BijliError) -> String {
    e.to_string()
}

/// The f64 payload grid.
fn values() -> Vec<f64> {
    let mut xs = specials();
    for k in (-1074..=1023).step_by(37) {
        xs.push(pow2(k));
    }
    for e in [-300, -30, -7, -6, -5, -1, 0, 1, 2, 15, 16, 17, 22, 300] {
        let x: f64 = format!("1e{e}").parse().unwrap_or(0.0);
        xs.push(x);
    }
    xs.extend([
        0.1,
        -0.5,
        100.0,
        2.5e-8,
        1.0 / 3.0,
        299_792_458.0,
        8.854_187_812_8e-12,
    ]);
    let mut rng = Rng::new(0xe220);
    for _ in 0..200 {
        let x = f64::from_bits(rng.bits());
        if !x.is_nan() {
            xs.push(x);
        }
    }
    let negs: Vec<f64> = xs.iter().step_by(3).map(|x| -x).collect();
    xs.extend(negs);
    xs
}

/// The usize payload grid (every value an `i:` token can carry).
const COUNTS: [i64; 9] = [0, 1, 2, 3, 4, 10, 1000, 4_294_967_296, i64::MAX];

/// The string payload grid.
const TEXTS: [&str; 8] = [
    "",
    "negative frequency",
    "field magnitude",
    "refractive index n₂ cannot be zero",
    "100% — not a percentage",
    "  spaced  out  ",
    "quote \" and backslash \\",
    "daimon: registration request failed: connection refused",
];

pub fn emit() -> Module {
    let mut m = Module::new("error");
    for &a in &COUNTS {
        for &b in &COUNTS {
            let (ua, ub) = (
                usize::try_from(a).unwrap_or(0),
                usize::try_from(b).unwrap_or(0),
            );
            let e = BijliError::DimensionMismatch {
                expected: ua,
                got: ub,
            };
            m.row(
                "BijliError::DimensionMismatch",
                Class::A,
                &[i(a), i(b)],
                &[err(&e), i(a), i(b), s(display(&e))],
            );
        }
    }
    let e = BijliError::ZeroCharge;
    m.row(
        "BijliError::ZeroCharge",
        Class::A,
        &[],
        &[err(&e), s(display(&e))],
    );
    let e = BijliError::Singularity;
    m.row(
        "BijliError::Singularity",
        Class::A,
        &[],
        &[err(&e), s(display(&e))],
    );
    let xs = values();
    for (k, &x) in xs.iter().enumerate() {
        let c = COUNTS[k % COUNTS.len()];
        let e = BijliError::InsufficientResolution {
            cells: usize::try_from(c).unwrap_or(0),
            domain_size: x,
        };
        m.row(
            "BijliError::InsufficientResolution",
            Class::A,
            &[i(c), f(x)],
            &[err(&e), i(c), f(x), s(display(&e))],
        );
        let e = BijliError::InvalidPermittivity { value: x };
        m.row(
            "BijliError::InvalidPermittivity",
            Class::A,
            &[f(x)],
            &[err(&e), f(x), s(display(&e))],
        );
        let e = BijliError::InvalidPermeability { value: x };
        m.row(
            "BijliError::InvalidPermeability",
            Class::A,
            &[f(x)],
            &[err(&e), f(x), s(display(&e))],
        );
        let prefix = "frequency must be positive, got ";
        let e = BijliError::InvalidParameter {
            reason: format!("{prefix}{x}"),
        };
        m.row(
            "BijliError::InvalidParameter/f64",
            Class::A,
            &[s(prefix), f(x)],
            &[err(&e), s(format!("{prefix}{x}")), s(display(&e))],
        );
    }
    for t in TEXTS {
        let e = BijliError::DivisionByZero { context: t.into() };
        m.row(
            "BijliError::DivisionByZero",
            Class::A,
            &[s(t)],
            &[err(&e), s(t), s(display(&e))],
        );
        let e = BijliError::InvalidParameter { reason: t.into() };
        m.row(
            "BijliError::InvalidParameter",
            Class::A,
            &[s(t)],
            &[err(&e), s(t), s(display(&e))],
        );
    }
    m
}
