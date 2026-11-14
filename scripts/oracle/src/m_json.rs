//! `json.tsv` — serde_json 1.0.154 on the serializable types (port plan D5; M1 part 4, M2 part 1a, M3 part 1a-2):
//! Rust hisab 1.4.0's `Complex`, `field::FieldVector`, `charge::PointCharge`, `material::Material` and
//! polarization's `JonesVector`, `JonesMatrix`, `StokesVector` and `MuellerMatrix`, written with
//! `serde_json::to_string` and read with `serde_json::from_str`.
//!
//! - **Writer rows** (`json::<Type>::to_string`, class A): the value's fields in, serde_json's text out as *parsed
//!   tokens*: each number in the text (outside strings, a maximal run of `-0123456789.eE+`) is an `f:` token, the
//!   correctly rounded parse of its digits, and every run between numbers is an `s:` token. The port's text,
//!   tokenized the same way, must equal it: the same punctuation, keys and `null`s byte for byte, and the same bits
//!   for every number, whatever notation each printer chose (D5: bayan writes `10000000000000000.0` where ryu
//!   writes `1e16`). The note keeps serde_json's raw text.
//! - **Reader rows** (`json::<Type>::from_str`, class C): the document in (`s:`), then `i:<flags>`, what serde_json
//!   reads without checking (`unchecked`: 1 a repeated key, 2 a lone surrogate escape, 4 nesting past 128 levels; ADR
//!   0002 row 12: the port refuses those wherever they are, serde_json only in what it deserializes); the fields out,
//!   or `E:Deserialize` with serde_json's message in the note. serde_json's default float parse is best-effort, not correctly rounded: where it misrounds, the row
//!   is serde's and `json.cr.tsv` (validate.py, from Python's correctly rounded `float`) holds the correctly rounded
//!   value the port must give (ADR 0002 row 11).
//!
//! The documents: serde_json's own text for every writer value; the same values in Rust's `{}` (every digit, no
//! exponent), `{:e}` and bayan's notation (`js`); every number spelling worth a row (signs, `-0`, integers past
//! i64 and u64, overflow and underflow edges, half-way cases, long significands, refused spellings) in each field
//! position; seeded random decimals of up to 40 digits; and the shapes: field order, whitespace, unknown keys of
//! every JSON kind (deep, long, escaped, lone surrogates), missing and repeated fields, wrong types, the sequence
//! form, every prefix of a document, trailing data, and nesting at serde's recursion limit. The polarization types
//! nest, so their documents add the nested forms (`nested`): each inner struct and array in each of serde's forms,
//! short, long, of the wrong type, with its fields reordered, missing, repeated or unknown.

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, Tok, f, i, s};
use bijli::charge::{ELECTRON_MASS, ELEMENTARY_CHARGE, PROTON_MASS, PointCharge};
use bijli::field::{EPSILON_0, FieldVector, MU_0, SPEED_OF_LIGHT};
use bijli::material::Material;
use bijli::polarization::{Complex, JonesMatrix, JonesVector, MuellerMatrix, StokesVector};
use std::collections::BTreeSet;
use std::hint::black_box;

/// The types, by the name their rows carry. New ones go last: each type's rows draw on one seeded stream in turn.
#[derive(Clone, Copy)]
enum Kind {
    Complex,
    FieldVector,
    PointCharge,
    Material,
    JonesVector,
    JonesMatrix,
    StokesVector,
    MuellerMatrix,
}

/// `{"re":<a>,"im":<b>}`.
fn cx(a: &str, b: &str) -> String {
    format!("{{\"re\":{a},\"im\":{b}}}")
}

/// A Complex in serde's sequence form, `[<a>,<b>]`.
fn cxs(a: &str, b: &str) -> String {
    format!("[{a},{b}]")
}

/// `[<v0>,<v1>,..]`.
fn arr(v: &[String]) -> String {
    format!("[{}]", v.join(","))
}

/// A Jones matrix's `m` from eight spellings, its four Complex written by `c`.
fn jm(v: &[String], c: fn(&str, &str) -> String) -> String {
    format!("[[{},{}],[{},{}]]", c(&v[0], &v[1]), c(&v[2], &v[3]), c(&v[4], &v[5]), c(&v[6], &v[7]))
}

/// A Mueller matrix's `m` from sixteen spellings, row by row.
fn mm(v: &[String]) -> String {
    arr(&(0..4).map(|r| arr(&v[4 * r..4 * r + 4])).collect::<Vec<_>>())
}

impl Kind {
    const ALL: [Kind; 8] = [
        Kind::Complex,
        Kind::FieldVector,
        Kind::PointCharge,
        Kind::Material,
        Kind::JonesVector,
        Kind::JonesMatrix,
        Kind::StokesVector,
        Kind::MuellerMatrix,
    ];

    fn name(self) -> &'static str {
        match self {
            Kind::Complex => "Complex",
            Kind::FieldVector => "FieldVector",
            Kind::PointCharge => "PointCharge",
            Kind::Material => "Material",
            Kind::JonesVector => "JonesVector",
            Kind::JonesMatrix => "JonesMatrix",
            Kind::StokesVector => "StokesVector",
            Kind::MuellerMatrix => "MuellerMatrix",
        }
    }

    /// How many f64 the type holds.
    fn width(self) -> usize {
        match self {
            Kind::Complex => 2,
            Kind::FieldVector => 3,
            Kind::PointCharge => 8,
            Kind::Material => 3,
            Kind::JonesVector | Kind::StokesVector => 4,
            Kind::JonesMatrix => 8,
            Kind::MuellerMatrix => 16,
        }
    }

    /// How many values a row's set holds: eight, or the width when wider (the first `width()` used).
    fn draws(self) -> usize {
        self.width().max(8)
    }

    /// Its document with the number spellings `v` (the first `width()` used), as serde_json shapes it.
    fn doc(self, v: &[String]) -> String {
        match self {
            Kind::Complex => format!("{{\"re\":{},\"im\":{}}}", v[0], v[1]),
            Kind::FieldVector => format!("{{\"x\":{},\"y\":{},\"z\":{}}}", v[0], v[1], v[2]),
            Kind::PointCharge => format!(
                "{{\"charge\":{},\"mass\":{},\"position\":[{},{},{}],\"velocity\":[{},{},{}]}}",
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
            ),
            Kind::Material => format!("{{\"eps_r\":{},\"mu_r\":{},\"conductivity\":{}}}", v[0], v[1], v[2]),
            Kind::JonesVector => format!("{{\"x\":{},\"y\":{}}}", cx(&v[0], &v[1]), cx(&v[2], &v[3])),
            Kind::JonesMatrix => format!("{{\"m\":{}}}", jm(v, cx)),
            Kind::StokesVector => format!("{{\"s\":{}}}", arr(&v[..4])),
            Kind::MuellerMatrix => format!("{{\"m\":{}}}", mm(v)),
        }
    }

    /// The same values in serde's sequence form, at every level.
    fn seq(self, v: &[String]) -> String {
        match self {
            Kind::Complex => format!("[{},{}]", v[0], v[1]),
            Kind::FieldVector => format!("[{},{},{}]", v[0], v[1], v[2]),
            Kind::PointCharge => format!(
                "[{},{},[{},{},{}],[{},{},{}]]",
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
            ),
            Kind::Material => format!("[{},{},{}]", v[0], v[1], v[2]),
            Kind::JonesVector => format!("[{},{}]", cxs(&v[0], &v[1]), cxs(&v[2], &v[3])),
            Kind::JonesMatrix => format!("[{}]", jm(v, cxs)),
            Kind::StokesVector => format!("[{}]", arr(&v[..4])),
            Kind::MuellerMatrix => format!("[{}]", mm(v)),
        }
    }

    /// Its top-level fields with the spellings `v`: each name and its value's text, as `doc` writes them.
    fn members(self, v: &[String]) -> Vec<(&'static str, String)> {
        let a3 = |j: usize| format!("[{},{},{}]", v[j], v[j + 1], v[j + 2]);
        match self {
            Kind::Complex => vec![("re", v[0].clone()), ("im", v[1].clone())],
            Kind::FieldVector => vec![("x", v[0].clone()), ("y", v[1].clone()), ("z", v[2].clone())],
            Kind::PointCharge => {
                vec![("charge", v[0].clone()), ("mass", v[1].clone()), ("position", a3(2)), ("velocity", a3(5))]
            }
            Kind::Material => vec![("eps_r", v[0].clone()), ("mu_r", v[1].clone()), ("conductivity", v[2].clone())],
            Kind::JonesVector => vec![("x", cx(&v[0], &v[1])), ("y", cx(&v[2], &v[3]))],
            Kind::JonesMatrix => vec![("m", jm(v, cx))],
            Kind::StokesVector => vec![("s", arr(&v[..4]))],
            Kind::MuellerMatrix => vec![("m", mm(v))],
        }
    }

    /// The value, from up to eight f64.
    fn write(self, v: &[f64]) -> String {
        match self {
            Kind::Complex => serde_json::to_string(&Complex::new(black_box(v[0]), black_box(v[1]))),
            Kind::FieldVector => {
                serde_json::to_string(&FieldVector::new(black_box(v[0]), black_box(v[1]), black_box(v[2])))
            }
            Kind::PointCharge => serde_json::to_string(&PointCharge {
                charge: black_box(v[0]),
                mass: black_box(v[1]),
                position: [black_box(v[2]), black_box(v[3]), black_box(v[4])],
                velocity: [black_box(v[5]), black_box(v[6]), black_box(v[7])],
            }),
            Kind::Material => serde_json::to_string(&Material::new(black_box(v[0]), black_box(v[1]), black_box(v[2]))),
            Kind::JonesVector => serde_json::to_string(&JonesVector {
                x: Complex::new(black_box(v[0]), black_box(v[1])),
                y: Complex::new(black_box(v[2]), black_box(v[3])),
            }),
            Kind::JonesMatrix => {
                let c = |k: usize| Complex::new(black_box(v[2 * k]), black_box(v[2 * k + 1]));
                serde_json::to_string(&JonesMatrix { m: [[c(0), c(1)], [c(2), c(3)]] })
            }
            Kind::StokesVector => serde_json::to_string(&StokesVector {
                s: [black_box(v[0]), black_box(v[1]), black_box(v[2]), black_box(v[3])],
            }),
            Kind::MuellerMatrix => {
                let row = |r: usize| [0, 1, 2, 3].map(|c| black_box(v[4 * r + c]));
                serde_json::to_string(&MuellerMatrix { m: [row(0), row(1), row(2), row(3)] })
            }
        }
        .expect("serde_json writes every value")
    }

    /// serde_json's reading of `doc`: the f64 in field order, or its error's message.
    fn read(self, doc: &str) -> Result<Vec<f64>, String> {
        let doc = black_box(doc);
        match self {
            Kind::Complex => serde_json::from_str::<Complex>(doc).map(|c| vec![c.re, c.im]),
            Kind::FieldVector => serde_json::from_str::<FieldVector>(doc).map(|v| vec![v.x, v.y, v.z]),
            Kind::PointCharge => serde_json::from_str::<PointCharge>(doc).map(|p| {
                let mut out = vec![p.charge, p.mass];
                out.extend(p.position);
                out.extend(p.velocity);
                out
            }),
            Kind::Material => {
                serde_json::from_str::<Material>(doc).map(|m| vec![m.eps_r, m.mu_r, m.conductivity])
            }
            Kind::JonesVector => serde_json::from_str::<JonesVector>(doc).map(|j| vec![j.x.re, j.x.im, j.y.re, j.y.im]),
            Kind::JonesMatrix => serde_json::from_str::<JonesMatrix>(doc)
                .map(|j| j.m.iter().flatten().flat_map(|c| [c.re, c.im]).collect()),
            Kind::StokesVector => serde_json::from_str::<StokesVector>(doc).map(|s| s.s.to_vec()),
            Kind::MuellerMatrix => {
                serde_json::from_str::<MuellerMatrix>(doc).map(|m| m.m.iter().flatten().copied().collect())
            }
        }
        .map_err(|e| e.to_string())
    }
}

// ── the parsed-token form of a document ─────────────────────────────────────────────

/// The tokens of JSON text `t`: `f:` for each number outside strings (its correctly rounded value), `s:` for each run
/// between them.
fn tokens(t: &str) -> Vec<Tok> {
    let b = t.as_bytes();
    let mut out = Vec::new();
    let mut run = String::new();
    let mut k = 0;
    while k < b.len() {
        let c = b[k];
        if c == b'"' {
            let start = k;
            k += 1;
            while k < b.len() && b[k] != b'"' {
                k += if b[k] == b'\\' { 2 } else { 1 };
            }
            k = (k + 1).min(b.len());
            run.push_str(&t[start..k]);
        } else if c == b'-' || c.is_ascii_digit() {
            let start = k;
            while k < b.len() && (b[k].is_ascii_digit() || b".eE+-".contains(&b[k])) {
                k += 1;
            }
            if !run.is_empty() {
                out.push(s(std::mem::take(&mut run)));
            }
            out.push(f(t[start..k].parse::<f64>().expect("a JSON number parses")));
        } else {
            run.push(char::from(c));
            k += 1;
        }
    }
    if !run.is_empty() {
        out.push(s(run));
    }
    out
}

// ── number spellings ────────────────────────────────────────────────────────────────

/// x in bayan's notation (JavaScript's, plus `.0` on an integer): the shortest digits; positional while the decimal
/// point falls in (−6, 21], else `d.ddde±x`. `null` for NaN and ±∞.
fn js(x: f64) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0" } else { "0.0" }.to_string();
    }
    let e = format!("{:e}", x.abs());
    let (mant, exp) = e.split_once('e').expect("LowerExp has an exponent");
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let n = i32::try_from(digits.len()).expect("17 digits at most");
    let k = exp.parse::<i32>().expect("an integer exponent") + 1;
    let sign = if x < 0.0 { "-" } else { "" };
    let body = if k >= n && k <= 21 {
        format!("{digits}{}.0", "0".repeat(usize::try_from(k - n).expect("k ≥ n")))
    } else if k > 0 && k <= 21 {
        let at = usize::try_from(k).expect("k > 0");
        format!("{}.{}", &digits[..at], &digits[at..])
    } else if k > -6 && k <= 0 {
        format!("0.{}{digits}", "0".repeat(usize::try_from(-k).expect("k ≤ 0")))
    } else {
        let m = if n == 1 {
            digits.clone()
        } else {
            format!("{}.{}", &digits[..1], &digits[1..])
        };
        let ee = k - 1;
        format!("{m}e{}{}", if ee < 0 { '-' } else { '+' }, ee.abs())
    };
    format!("{sign}{body}")
}

/// The f64 values of the writer rows: the specials, the notation edges, SI and bijli constants, every 131st binade,
/// and seeded bit patterns of every class (NaNs of both signs and payloads among them).
fn values(rng: &mut Rng) -> Vec<f64> {
    let mut xs = specials();
    xs.extend([
        -f64::MAX,
        -f64::MIN_POSITIVE,
        -f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::from_bits(0xfff8_0000_0000_0001),
        f64::from_bits(0x7ff0_0000_0000_0001),
        0.1,
        0.2,
        0.3,
        1.0 / 3.0,
        2.0 / 3.0,
        1.5,
        -2.5,
        100.0,
        1e15,
        1e16,
        1e17,
        1e20,
        1e21,
        1e22,
        1e23,
        123_456_789_012_345_680.0,
        9_007_199_254_740_993.0,
        9_223_372_036_854_775_807.0,
        18_446_744_073_709_551_615.0,
        1e-5,
        1e-6,
        1e-7,
        1.5e-6,
        1.5e-7,
        5e-324,
        1e-300,
        1e300,
        1.797_693_134_862_315_7e308,
        2.225_073_858_507_201e-308,
        EPSILON_0,
        MU_0,
        SPEED_OF_LIGHT,
        ELEMENTARY_CHARGE,
        -ELEMENTARY_CHARGE,
        ELECTRON_MASS,
        PROTON_MASS,
    ]);
    for k in (-1074..=1023).step_by(131) {
        xs.push(pow2(k));
        xs.push(-pow2(k));
    }
    for _ in 0..600 {
        xs.push(f64::from_bits(rng.bits()));
    }
    for _ in 0..200 {
        xs.push(rng.signed_log_uniform(-30.0, 30.0));
    }
    xs
}

/// Number spellings for one field: every form worth a row, the refused ones among them.
fn spellings() -> Vec<String> {
    let mut out: Vec<String> = [
        "0", "-0", "0.0", "-0.0", "0e0", "-0e-0", "0E+5", "1", "-1", "1.5", "-1.5", "1e1", "1E1", "1e+1", "1e-1",
        "1E-001", "10", "100e-2", "12.5e-1", "9007199254740992", "9007199254740993", "9007199254740995",
        "9223372036854775807", "9223372036854775808", "-9223372036854775808", "-9223372036854775809",
        "18446744073709551615", "18446744073709551616", "-18446744073709551616", "123456789012345678901234567890",
        "1e308", "1e309", "-1e309", "1e400", "-1e400", "1.7976931348623157e308", "1.7976931348623158e308",
        "1.7976931348623159e308", "1.797693134862315807e308", "1.7976931348623158079e308",
        "1.797693134862315808e308", "179769313486231570000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        "4.9406564584124654e-324", "5e-324", "2.4703282292062327e-324", "2.4703282292062328e-324",
        "2.47032822920623272e-324", "1e-324", "3e-324", "2.2250738585072011e-308", "2.2250738585072012e-308",
        "2.2250738585072014e-308", "1e-400", "-1e-400", "1e-99999", "1e99999", "1e00000000000000000000001",
        "1e-00000000000000000000001", "0.1", "0.30000000000000004", "0.1000000000000000055511151231257827021181583404541015625",
        "3.14159265358979323846264338327950288419716939937510", "2.718281828459045235360287471352662497757",
        "0.000000000000000000000000000001e30", "100000000000000000000000000000000000000e-38",
        "7.2057594037927933e16", "8.98846567431158e307", "4.4501477170144023e-308", "9.5e-322",
        "123456789012345678e-400", "0.0000000000000000000000000000000000000000000001", "1.00000000000000011102230246251565404236316680908203125",
        "1.00000000000000011102230246251565404236316680908203124", "1.00000000000000011102230246251565404236316680908203126",
        "01", "-01", "00", "1.", ".1", "+1", "1e", "1e+", "1e-", "-", "--1", "0x10", "1_000", "1.0f", "NaN", "-NaN",
        "Infinity", "-Infinity", "nan", "inf", "null", "true", "false", "\"1\"", "\"\"", "[]", "[1]", "{}", "1 2",
        "\u{ff11}",
    ]
    .iter()
    .map(|x| (*x).to_string())
    .collect();
    out.push(format!("1{}", "0".repeat(400)));
    out.push(format!("0.{}1", "0".repeat(400)));
    out.push(format!("{}e-300", "9".repeat(320)));
    out.push(format!("1.{}e5", "0".repeat(800)));
    out
}

/// Seeded decimal spellings of up to 40 digits with exponents across the f64 range (correct rounding's hard cases).
fn random_decimals(rng: &mut Rng, n: usize) -> Vec<String> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let len = usize::try_from(rng.int(1, 40)).expect("1..=40");
        let mut d = String::with_capacity(len + 8);
        d.push(char::from(b'1' + u8::try_from(rng.int(0, 8)).expect("a digit")));
        for _ in 1..len {
            d.push(char::from(b'0' + u8::try_from(rng.int(0, 9)).expect("a digit")));
        }
        let point = rng.int(0, i64::try_from(len).expect("len"));
        let mut t = String::new();
        if rng.int(0, 1) == 1 {
            t.push('-');
        }
        let p = usize::try_from(point).expect("point ≥ 0");
        if p == 0 {
            t.push_str("0.");
            t.push_str(&d);
        } else if p >= len {
            t.push_str(&d);
        } else {
            t.push_str(&d[..p]);
            t.push('.');
            t.push_str(&d[p..]);
        }
        if rng.int(0, 3) != 0 {
            t.push_str(&format!("e{}", rng.int(-345, 310)));
        }
        out.push(t);
    }
    out
}

// ── what serde_json leaves unchecked ────────────────────────────────────────────────

/// The string starting after the quote at byte k of `t`: its code points (escapes decoded; a surrogate escape that
/// pairs joined, one that doesn't kept as its own unit), whether it holds a lone surrogate escape, and the index past
/// its closing quote.
fn string_at(t: &str, mut k: usize) -> (Vec<u32>, bool, usize) {
    let b = t.as_bytes();
    let mut units: Vec<u32> = Vec::new();
    let mut lone = false;
    let mut pending: Option<u32> = None; // a high surrogate escape waiting for its low half
    let hex = |at: usize| -> Option<u32> {
        let h = b.get(at..at + 4)?;
        u32::from_str_radix(std::str::from_utf8(h).ok()?, 16).ok()
    };
    while k < b.len() && b[k] != b'"' {
        let unit = if b[k] == b'\\' && k + 1 < b.len() && b[k + 1] == b'u' {
            let u = hex(k + 2);
            k += 6;
            match u {
                Some(u) => u,
                None => continue,
            }
        } else if b[k] == b'\\' && k + 1 < b.len() {
            let e = match b[k + 1] {
                b'b' => 8,
                b'f' => 12,
                b'n' => 10,
                b'r' => 13,
                b't' => 9,
                c => u32::from(c),
            };
            k += 2;
            e
        } else {
            let ch = t[k..].chars().next().expect("a char at a char boundary");
            k += ch.len_utf8();
            u32::from(ch)
        };
        if (0xdc00..0xe000).contains(&unit) {
            if let Some(h) = pending.take() {
                units.push(0x10000 + ((h - 0xd800) << 10) + (unit - 0xdc00));
            } else {
                lone = true;
                units.push(unit);
            }
            continue;
        }
        if let Some(h) = pending.take() {
            lone = true;
            units.push(h);
        }
        if (0xd800..0xdc00).contains(&unit) {
            pending = Some(unit);
        } else {
            units.push(unit);
        }
    }
    if let Some(h) = pending.take() {
        lone = true;
        units.push(h);
    }
    (units, lone, (k + 1).min(b.len()))
}

/// What serde_json reads in `doc` without checking, as flags (ADR 0002 row 12; the port refuses any of them): 1 a key
/// repeats in an object (compared decoded), 2 a string holds a lone surrogate escape, 4 containers nest past 128
/// levels. Asked only of a document serde_json reads (a row it refuses carries 0: the port refuses it too).
fn unchecked(doc: &str) -> i64 {
    let b = doc.as_bytes();
    let mut flags = 0;
    // The open containers: an object's keys so far and whether its next string is a key; None for an array.
    let mut stack: Vec<Option<(BTreeSet<Vec<u32>>, bool)>> = Vec::new();
    let mut k = 0;
    while k < b.len() {
        match b[k] {
            b'{' | b'[' => {
                stack.push((b[k] == b'{').then(|| (BTreeSet::new(), true)));
                if stack.len() > 128 {
                    flags |= 4;
                }
                k += 1;
            }
            b'}' | b']' => {
                stack.pop();
                k += 1;
            }
            b',' => {
                if let Some(Some((_, key_next))) = stack.last_mut() {
                    *key_next = true;
                }
                k += 1;
            }
            b'"' => {
                let (units, lone, end) = string_at(doc, k + 1);
                if lone {
                    flags |= 2;
                }
                if let Some(Some((keys, key_next))) = stack.last_mut()
                    && *key_next
                {
                    if !keys.insert(units) {
                        flags |= 1;
                    }
                    *key_next = false;
                }
                k = end;
            }
            _ => k += 1,
        }
    }
    flags
}

// ── rows ────────────────────────────────────────────────────────────────────────────

fn writer_row(m: &mut Module, kind: Kind, v: &[f64]) {
    let text = kind.write(v);
    let ins: Vec<Tok> = v[..kind.width()].iter().map(|&x| f(x)).collect();
    m.row_noted(
        &format!("json::{}::to_string", kind.name()),
        Class::A,
        &ins,
        &tokens(&text),
        &text,
    );
}

/// A document for a note: its first 160 characters, Debug-quoted, and its length when longer.
fn preview(doc: &str) -> String {
    let head: String = doc.chars().take(160).collect();
    if head.len() == doc.len() {
        format!("{head:?}")
    } else {
        format!("{head:?}… ({} bytes)", doc.len())
    }
}

fn reader_row(m: &mut Module, kind: Kind, doc: &str) {
    let func = format!("json::{}::from_str", kind.name());
    let read = kind.read(doc);
    let d = if read.is_ok() { unchecked(doc) } else { 0 };
    let ins = [s(doc), i(d)];
    let head = format!("{func}({}, unchecked {d})", preview(doc));
    match read {
        Ok(v) => {
            let outs: Vec<Tok> = v.iter().map(|&x| f(x)).collect();
            let vals: Vec<String> = v.iter().map(|x| format!("{x:?}")).collect();
            m.row_with_note(&func, Class::C, &ins, &outs, &format!("{head} -> Ok({})", vals.join(", ")));
        }
        Err(msg) => {
            let outs = [Tok::Err("Deserialize".to_string())];
            m.row_with_note(&func, Class::C, &ins, &outs, &format!("{head} -> Err({msg})"));
        }
    }
}

/// The plain spellings of a set of seeded values (`Kind::draws`) for a document's other fields.
fn fillers(kind: Kind, rng: &mut Rng) -> Vec<String> {
    (0..kind.draws()).map(|_| js(rng.signed_log_uniform(-20.0, 20.0))).collect()
}

/// The shapes: field order, whitespace, unknown keys, missing and repeated fields, wrong types, the sequence form,
/// prefixes, trailing data, nesting.
fn shapes(kind: Kind, rng: &mut Rng) -> Vec<String> {
    let v = fillers(kind, rng);
    let base = kind.doc(&v);
    let inner = &base[1..base.len() - 1];
    let fields = kind.members(&v);
    let mut docs = vec![base.clone(), kind.seq(&v)];

    // Field values as written in `base`, by name, for reordering and repeating.
    let members: Vec<String> = fields.iter().map(|(name, val)| format!("\"{name}\":{val}")).collect();
    let mut rev = members.clone();
    rev.reverse();
    docs.push(format!("{{{}}}", rev.join(",")));
    docs.push(format!(" \t\n\r{{ \t\n\r{} \t\n\r}} \t\n\r", members.join(" \t\n\r, \t\n\r").replace(':', " \t\n\r: \t\n\r")));
    for ws in ["\u{c}", "\u{b}", "\u{a0}", "\u{2028}", "\u{feff}"] {
        docs.push(format!("{ws}{base}"));
        docs.push(format!("{{{ws}{inner}}}"));
    }
    docs.push(format!("\u{feff}{base}"));

    // Unknown keys of every kind, before, between and after the fields.
    let unknown_values: Vec<String> = vec![
        "null".into(),
        "true".into(),
        "false".into(),
        "0".into(),
        "-0".into(),
        "1e400".into(),
        "-1e400".into(),
        "9223372036854775808".into(),
        format!("1{}", "0".repeat(400)),
        "\"str\"".into(),
        "\"esc\\\"\\\\\\/\\b\\f\\n\\r\\t\\u0041\\u00e9\"".into(),
        "\"\\ud800\"".into(),
        "\"\\udc00\"".into(),
        "\"\\ud83d\\ude00\"".into(),
        "\"\\ud800\\u0041\"".into(),
        "\"\u{e9}\u{1f600}\"".into(),
        "[]".into(),
        "{}".into(),
        "[1,[2,[3,{\"a\":[]}]]]".into(),
        "{\"a\":{\"b\":{\"x\":1}}}".into(),
        "{\"x\":\"not a number\"}".into(),
        "[1,2,3,]".into(),
        "{\"a\":1,}".into(),
        "\"a\u{1}b\"".into(),
        "\"\\x41\"".into(),
        "\"\\u004\"".into(),
        format!("\"{}\"", "a".repeat(20_000)),
        format!("[{}]", (0..2000).map(|k| k.to_string()).collect::<Vec<_>>().join(",")),
        format!("{{{}}}", (0..100).map(|k| format!("\"k{k}\":{k}")).collect::<Vec<_>>().join(",")),
    ];
    for u in &unknown_values {
        docs.push(format!("{{\"w\":{u},{inner}}}"));
        docs.push(format!("{{{inner},\"unknown\":{u}}}"));
    }
    let mid = members.len() / 2;
    docs.push(format!("{{{},\"w\":1,{}}}", members[..mid].join(","), members[mid..].join(",")));
    for key in ["\"X\"", "\"\"", "\"\\u0078\"", "\"\\ud800\"", "\"\\udc00x\"", "\"\\ud83d\\ude00\"", "\"a\\u0000b\""] {
        docs.push(format!("{{{key}:1,{inner}}}"));
    }
    // A field's name escaped: serde decodes keys before matching them.
    let first = fields[0].0;
    let esc: String = first.chars().map(|c| format!("\\u{:04x}", c as u32)).collect();
    docs.push(base.replacen(&format!("\"{first}\""), &format!("\"{esc}\""), 1));

    // Nesting about serde_json's recursion limit (128), in an unknown value: depth counts the document's own object.
    for depth in [120usize, 125, 126, 127, 128, 129, 130, 200] {
        let d = depth - 1;
        docs.push(format!("{{\"w\":{}{}, {inner}}}", "[".repeat(d), "]".repeat(d)));
        docs.push(format!("{{\"w\":{}0{}, {inner}}}", "[".repeat(d), "]".repeat(d)));
        docs.push(format!("{{{inner},\"w\":{}1{}}}", "{\"a\":".repeat(d), "}".repeat(d)));
    }

    // Missing fields, each; the empty object; not an object.
    for k in 0..members.len() {
        let mut ms = members.clone();
        ms.remove(k);
        docs.push(format!("{{{}}}", ms.join(",")));
    }
    for d in ["{}", "[]", "null", "true", "1", "\"x\"", "", " ", "{", "}", "[", "]"] {
        docs.push(d.to_string());
    }

    // Repeated fields: a known one (with the same and with another value), an unknown one, one in an unknown value.
    docs.push(format!("{{{inner},{}}}", members[0]));
    docs.push(format!("{{{inner},\"{}\":7}}", fields[0].0));
    docs.push(format!("{{\"w\":1,{inner},\"w\":2}}"));
    docs.push(format!("{{\"w\":{{\"a\":1,\"a\":2}},{inner}}}"));
    docs.push(format!("{{\"w\":[{{\"a\":1,\"b\":2,\"a\":3}}],{inner}}}"));

    // Wrong types in each field.
    for k in 0..members.len() {
        for bad in ["null", "true", "\"1\"", "{}", "[]", "[1]", "[1,2,3]", "[[1,2,3]]", "{\"0\":1}", "[1,2,3,4]", "[1,2]",
            "[1,null,3]", "[1,\"2\",3]"]
        {
            let mut ms = members.clone();
            ms[k] = format!("\"{}\":{bad}", fields[k].0);
            docs.push(format!("{{{}}}", ms.join(",")));
        }
    }

    // The sequence form: short, long, wrong element types, trailing comma, nested.
    let width_top = fields.len();
    let parts: Vec<String> = members.iter().map(|m| m.split_once(':').expect("a member").1.to_string()).collect();
    for take in 0..width_top {
        docs.push(format!("[{}]", parts[..take].join(",")));
    }
    docs.push(format!("[{},1]", parts.join(",")));
    docs.push(format!("[{},]", parts.join(",")));
    docs.push(format!("[{}]", base));
    docs.push(format!("[null,{}]", parts[1..].join(",")));
    docs.push(format!("[{}]", parts.iter().map(|p| format!("[{p}]")).collect::<Vec<_>>().join(",")));

    // Every prefix of the document, and trailing data.
    for k in 0..base.len() {
        if base.is_char_boundary(k) {
            docs.push(base[..k].to_string());
        }
    }
    for tail in [" ", "\n", "x", "}", "]", ",", "0", "{}", " {}", "null", "//", "/* */", "\u{0}"] {
        docs.push(format!("{base}{tail}"));
    }

    // Syntax: quotes, colons, commas, comments.
    docs.push(inner.to_string());
    docs.push(base.replace('"', "'"));
    docs.push(base.replacen(&format!("\"{first}\""), first, 1));
    docs.push(base.replacen(':', "", 1));
    docs.push(base.replacen(':', "::", 1));
    docs.push(base.replacen(',', "", 1));
    docs.push(base.replacen(',', ",,", 1));
    docs.push(format!("{{/*c*/{inner}}}"));
    docs.push(format!("{{{inner}}}//c"));
    docs.push(format!("[{base}"));
    docs
}

/// Wrong values for an inner struct or array: every JSON kind, and arrays of every length to six.
const INNER_BAD: [&str; 14] = [
    "null", "true", "1", "\"1\"", "{}", "[]", "[1]", "[1,2]", "[1,2,3]", "[1,2,3,4]", "[1,2,3,4,5]", "[1,2,3,4,5,6]",
    "[null,1]", "{\"0\":1,\"1\":2}",
];

/// An inner Complex holding (a, b) in every form serde reads or refuses: by position, reordered, a field missing,
/// repeated, unknown, escaped or of the wrong type, a key or string with a lone surrogate, and the INNER_BAD values.
fn complex_forms(a: &str, b: &str) -> Vec<String> {
    let mut out = vec![
        cxs(a, b),
        format!("{{\"im\":{b},\"re\":{a}}}"),
        format!("{{\"re\":{a}}}"),
        format!("{{\"im\":{b}}}"),
        format!("[{a}]"),
        format!("[{a},{b},{a}]"),
        format!("[[{a},{b}]]"),
        format!("{{\"re\":{a},\"im\":{b},\"re\":{a}}}"),
        format!("{{\"re\":{a},\"im\":{b},\"w\":1}}"),
        format!("{{\"w\":[1,{{\"x\":2}}],\"re\":{a},\"im\":{b}}}"),
        format!("{{\"re\":{a},\"im\":{b},\"w\":1,\"w\":2}}"),
        format!("{{\"\\u0072\\u0065\":{a},\"im\":{b}}}"),
        format!("{{\"\\ud800\":1,\"re\":{a},\"im\":{b}}}"),
        format!("{{\"re\":{a},\"im\":{b},\"w\":\"\\udc00\"}}"),
        format!("{{\"re\":\"{a}\",\"im\":{b}}}"),
        format!("{{\"re\":{a},\"im\":null}}"),
        format!("{{\"re\":{a},\"im\":[{b}]}}"),
        format!("{{\"Re\":{a},\"im\":{b}}}"),
    ];
    out.extend(INNER_BAD.iter().map(|x| (*x).to_string()));
    out
}

/// An inner `[f64; n]` holding v in every form serde reads or refuses: short, long, nested, an element of the wrong
/// type, an object, spaced, and the INNER_BAD values.
fn array_forms(v: &[String]) -> Vec<String> {
    let n = v.len();
    let head = v[..n - 1].join(",");
    let keyed: Vec<String> = v.iter().enumerate().map(|(k, x)| format!("\"{k}\":{x}")).collect();
    let mut out = vec![
        arr(&v[..n - 1]),
        format!("[{},{}]", v.join(","), v[0]),
        format!("[{}]", arr(v)),
        format!("[{head},null]"),
        format!("[{head},\"{}\"]", v[n - 1]),
        format!("[{head},[{}]]", v[n - 1]),
        format!("{{{}}}", keyed.join(",")),
        format!("[ {} ]", v.join(" , ")),
    ];
    out.extend(INNER_BAD.iter().map(|x| (*x).to_string()));
    out
}

/// The nested forms of a polarization document with the spellings v (none for the flat types): each inner Complex,
/// row and array in each of its forms and places, and the documents' mixed forms.
fn nested(kind: Kind, v: &[String]) -> Vec<String> {
    let mut docs = Vec::new();
    match kind {
        Kind::JonesVector => {
            let (x, y) = (cx(&v[0], &v[1]), cx(&v[2], &v[3]));
            docs.push(format!("{{\"x\":{},\"y\":{y}}}", cxs(&v[0], &v[1])));
            docs.push(format!("[{x},{y}]"));
            docs.push(format!("[{x},{}]", cxs(&v[2], &v[3])));
            docs.push(format!("{{\"y\":{y},\"x\":{x}}}"));
            for f in complex_forms(&v[0], &v[1]) {
                docs.push(format!("{{\"x\":{f},\"y\":{y}}}"));
                docs.push(format!("[{f},{y}]"));
            }
            for f in complex_forms(&v[2], &v[3]) {
                docs.push(format!("{{\"x\":{x},\"y\":{f}}}"));
            }
            docs.push(format!("{{\"x\":{x},\"y\":{y},\"x\":{x}}}"));
            docs.push(format!("{{\"x\":{x},\"y\":{y},\"z\":{x}}}"));
            docs.push(format!("[{x},{y},{x}]"));
        }
        Kind::JonesMatrix => {
            let c: Vec<String> = (0..4).map(|k| cx(&v[2 * k], &v[2 * k + 1])).collect();
            let (r0, r1) = (format!("[{},{}]", c[0], c[1]), format!("[{},{}]", c[2], c[3]));
            docs.push(format!("{{\"m\":{}}}", jm(v, cxs)));
            docs.push(format!("[{}]", jm(v, cx)));
            docs.push(format!("{{\"m\":[{r0},[{},{}]]}}", cxs(&v[4], &v[5]), c[3]));
            for f in complex_forms(&v[0], &v[1]) {
                docs.push(format!("{{\"m\":[[{f},{}],{r1}]}}", c[1]));
            }
            for f in complex_forms(&v[6], &v[7]) {
                docs.push(format!("{{\"m\":[{r0},[{},{f}]]}}", c[2]));
            }
            for f in INNER_BAD {
                docs.push(format!("{{\"m\":[{f},{r1}]}}"));
                docs.push(format!("{{\"m\":[{r0},{f}]}}"));
                docs.push(format!("{{\"m\":{f}}}"));
            }
            for bad in [
                format!("[{r0}]"),
                format!("[{r0},{r1},{r0}]"),
                format!("[[{}],{r1}]", c[0]),
                format!("[[{},{},{}],{r1}]", c[0], c[1], c[0]),
                format!("[{r0},[{}]]", c[2]),
                format!("[{},{},{},{}]", c[0], c[1], c[2], c[3]),
                format!("[[{r0},{r1}]]"),
                format!("{{\"0\":{r0},\"1\":{r1}}}"),
                format!("[{r0},{{\"0\":{},\"1\":{}}}]", c[2], c[3]),
                format!("[[{},null],{r1}]", c[0]),
            ] {
                docs.push(format!("{{\"m\":{bad}}}"));
            }
            docs.push(format!("{{\"m\":[{r0},{r1}],\"m\":[{r0},{r1}]}}"));
            docs.push(format!("{{\"M\":[{r0},{r1}]}}"));
            docs.push(format!("[[{r0},{r1}],1]"));
        }
        Kind::StokesVector => {
            let s4 = arr(&v[..4]);
            for f in array_forms(&v[..4]) {
                docs.push(format!("{{\"s\":{f}}}"));
                docs.push(format!("[{f}]"));
            }
            docs.push(format!("{{\"s\":{s4},\"s\":{s4}}}"));
            docs.push(format!("{{\"S\":{s4}}}"));
            docs.push(format!("[{s4},1]"));
            docs.push(s4);
        }
        Kind::MuellerMatrix => {
            let rows: Vec<String> = (0..4).map(|r| arr(&v[4 * r..4 * r + 4])).collect();
            for at in [0, 3] {
                for f in array_forms(&v[4 * at..4 * at + 4]) {
                    let mut rs = rows.clone();
                    rs[at] = f;
                    docs.push(format!("{{\"m\":{}}}", arr(&rs)));
                }
            }
            let keyed: Vec<String> = rows.iter().enumerate().map(|(k, r)| format!("\"{k}\":{r}")).collect();
            for bad in [
                arr(&rows[..3]),
                format!("[{},{}]", rows.join(","), rows[0]),
                arr(&v[..16]),
                format!("[{}]", mm(v)),
                format!("{{{}}}", keyed.join(",")),
                format!("[{},null]", rows[..3].join(",")),
            ] {
                docs.push(format!("{{\"m\":{bad}}}"));
                docs.push(format!("[{bad}]"));
            }
            for f in INNER_BAD {
                docs.push(format!("{{\"m\":{f}}}"));
            }
            docs.push(format!("{{\"m\":{},\"m\":{}}}", mm(v), mm(v)));
            docs.push(format!("{{\"M\":{}}}", mm(v)));
            docs.push(format!("[{},1]", mm(v)));
            docs.push(mm(v));
        }
        Kind::Complex | Kind::FieldVector | Kind::PointCharge | Kind::Material => {}
    }
    docs
}

/// The rows, each (function, inputs) once: a correction in `json.cr.tsv` keys on them.
struct Rows {
    m: Module,
    seen: BTreeSet<String>,
}

impl Rows {
    fn writer(&mut self, kind: Kind, v: &[f64]) {
        let key = format!("w{}{:?}", kind.name(), v[..kind.width()].iter().map(|x| x.to_bits()).collect::<Vec<_>>());
        if self.seen.insert(key) {
            writer_row(&mut self.m, kind, v);
        }
    }

    fn reader(&mut self, kind: Kind, doc: &str) {
        if self.seen.insert(format!("r{}{doc}", kind.name())) {
            reader_row(&mut self.m, kind, doc);
        }
    }
}

pub fn emit() -> Module {
    let mut rows = Rows { m: Module::new("json"), seen: BTreeSet::new() };
    let m = &mut rows;
    let mut rng = Rng::new(0x15_03d5);
    let xs = values(&mut rng);
    let n = xs.len();

    for kind in Kind::ALL {
        // Writer rows: each value, beside others from the list.
        let w = kind.width();
        let mut sets: Vec<Vec<f64>> = Vec::new();
        for j in 0..n {
            sets.push((0..kind.draws()).map(|k| xs[(j + k * 37) % n]).collect());
        }
        for _ in 0..300 {
            sets.push((0..kind.draws()).map(|_| f64::from_bits(rng.bits())).collect());
        }
        for v in &sets {
            m.writer(kind, v);
        }
        // Reader rows: serde's own text, and the same values in other notations (finite ones: `null` is refused, and
        // the shape rows cover it).
        for v in &sets {
            m.reader(kind, &kind.write(v));
            let finite = v[..w].iter().all(|x| x.is_finite());
            if finite {
                let d: Vec<String> = v.iter().map(|&x| js(x)).collect();
                m.reader(kind, &kind.doc(&d));
            }
        }
        for v in sets.iter().take(250) {
            if v[..w].iter().all(|x| x.is_finite()) {
                let plain: Vec<String> = v.iter().map(|x| format!("{x}")).collect();
                m.reader(kind, &kind.doc(&plain));
                let exp: Vec<String> = v.iter().map(|x| format!("{x:e}")).collect();
                m.reader(kind, &kind.doc(&exp));
            }
        }
        // Every spelling in each field position.
        let spell = spellings();
        let base = fillers(kind, &mut rng);
        for pos in 0..w {
            for sp in &spell {
                let mut v = base.clone();
                v[pos] = sp.clone();
                m.reader(kind, &kind.doc(&v));
            }
        }
        // Seeded decimals, a few in each position.
        for (k, d) in random_decimals(&mut rng, 1200).into_iter().enumerate() {
            let mut v = base.clone();
            v[k % w] = d;
            m.reader(kind, &kind.doc(&v));
        }
        for doc in shapes(kind, &mut rng) {
            m.reader(kind, &doc);
        }
        for doc in nested(kind, &base) {
            m.reader(kind, &doc);
        }
    }
    rows.m
}
