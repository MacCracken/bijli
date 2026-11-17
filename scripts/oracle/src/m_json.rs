//! `json.tsv` — serde_json 1.0.154 on the serializable types (port plan D5; M1 part 4, M2 part 1a, M3 parts 1a-2,
//! 1b-2 and 1c-2): Rust hisab 1.4.0's `Complex`, `field::FieldVector`, `charge::PointCharge`, `material::Material`,
//! polarization's `JonesVector`, `JonesMatrix`, `StokesVector` and `MuellerMatrix`, scattering's `MieResult`, and
//! beam's `GaussianBeam`, `AbcdMatrix` and `ResonatorStability`, written with `serde_json::to_string` and read with
//! `serde_json::from_str`.
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
//!
//! `MieResult` holds two `Vec<Complex>` of any length, so its rows have their own layout (`emit_mie`): a value is its
//! four f64, then `i:` the length of `a` and its parts (re, im of each), then `i:` the length of `b` and its parts. The
//! values: real Mie results, seeded and special parts with every pair of lengths to 7, and 500-element vectors; the
//! documents add arrays of every length and element form, a `Vec` given as an object, and a deep unknown value inside
//! an element.
//!
//! Beam's three structs (M3 part 1c-2) are `Kind`s emitted after MieResult (`Kind::BEAM`). `ResonatorStability`'s
//! `is_stable` is the first bool: in a value it is 0.0 or 1.0 (false, true), in a row's tokens `i:0` or `i:1`, and in
//! a document `false` or `true`. Its position takes its own spellings (`bool_spellings`: serde reads only the two
//! literals, never 0, 1 or a string), and the number spellings go to the f64 positions only.
//!
//! rf's three types (M4 part 1d, `emit_rf`, after beam's) have their own layouts, validate.py's JSON_FIELDS ones: a
//! `usize` is `i:` its word, a `Vec` is `i:` its length then its elements, so an `SMatrix` is `i:` ports, `i:` the
//! row count, each row as `i:` its length and its entries' parts, `i:` the impedance count and the impedances; a
//! `TouchstoneData` is `i:` ports, `f:` z0, the frequencies as a `Vec`, the matrices as a `Vec`; a `TouchstoneFormat`
//! is `i:` the port's enum value (0 MagnitudeAngle, 1 RealImag, 2 DbAngle), written as its name. A read may be ragged:
//! the documents include every shape serde reads and the port's fixed structs can't hold (ADR 0002 row 36).
//!
//! The FDTD types (M5 part 1f, after rf's): the four materials are `Kind`s, all-f64 structs (`LossyMaterial`,
//! `DrudeMaterial`, `LorentzMaterial`, `DebyeMaterial`), and `NonUniformGrid` and `EigenmodeResult` hold three
//! `Vec<f64>` each (`emit_fdtd_vecs`): a value is each vector as `i:` its length and its elements. Their documents add
//! every triple of lengths to 3; an eigenmode result whose vectors differ in length is read by serde and refused by the
//! port (ADR 0002 row 36), and its writer values are equal-length, the only ones the port can build.
//!
//! The 2D FDTD types (M6 part 1d, `emit_fdtd2d`, after the FDTD vectors): `Mode2d` is laid out as `TouchstoneFormat`,
//! `i:` the port's enum value (0 Tm, 1 Te), written as its name; `CpmlParams` is a `Kind` whose first value is a usize
//! (`Kind::is_usize`: `thickness`, carried in a value's f64 as its bits, `i:` its word in a row, an integer in a
//! document), then four f64. Its writer values take thicknesses at usize's edges once each (`THICKNESSES`) and small
//! ones after; a thickness past 2^63 − 1 is read by serde and refused by the port (ADR 0002 row 36).

use crate::grids::{Rng, pow2, specials};
use crate::tsv::{Class, Module, Tok, b, f, i, s};
use bijli::beam::{AbcdMatrix, GaussianBeam, ResonatorStability};
use bijli::charge::{ELECTRON_MASS, ELEMENTARY_CHARGE, PROTON_MASS, PointCharge};
use bijli::fdtd::{
    CpmlParams, DebyeMaterial, DrudeMaterial, EigenmodeResult, LorentzMaterial, LossyMaterial, Mode2d, NonUniformGrid,
};
use bijli::field::{EPSILON_0, FieldVector, MU_0, SPEED_OF_LIGHT};
use bijli::material::Material;
use bijli::polarization::{Complex, JonesMatrix, JonesVector, MuellerMatrix, StokesVector};
use bijli::rf::{SMatrix, TouchstoneData, TouchstoneFormat};
use bijli::scattering::{MieResult, mie};
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
    GaussianBeam,
    AbcdMatrix,
    ResonatorStability,
    LossyMaterial,
    DrudeMaterial,
    LorentzMaterial,
    DebyeMaterial,
    CpmlParams,
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

    /// Beam's three structs (M3 part 1c-2), emitted after MieResult's rows.
    const BEAM: [Kind; 3] = [Kind::GaussianBeam, Kind::AbcdMatrix, Kind::ResonatorStability];

    /// The FDTD materials (M5 part 1f), emitted after rf's rows.
    const FDTD: [Kind; 4] = [Kind::LossyMaterial, Kind::DrudeMaterial, Kind::LorentzMaterial, Kind::DebyeMaterial];

    /// The 2D FDTD's struct (M6 part 1d), emitted after the FDTD vectors.
    const FDTD2D: [Kind; 1] = [Kind::CpmlParams];

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
            Kind::GaussianBeam => "GaussianBeam",
            Kind::AbcdMatrix => "AbcdMatrix",
            Kind::ResonatorStability => "ResonatorStability",
            Kind::LossyMaterial => "LossyMaterial",
            Kind::DrudeMaterial => "DrudeMaterial",
            Kind::LorentzMaterial => "LorentzMaterial",
            Kind::DebyeMaterial => "DebyeMaterial",
            Kind::CpmlParams => "CpmlParams",
        }
    }

    /// How many values the type holds: its f64, and its bool or usize if any (`is_bool`, `is_usize`).
    fn width(self) -> usize {
        match self {
            Kind::Complex | Kind::LossyMaterial => 2,
            Kind::CpmlParams => 5,
            Kind::FieldVector | Kind::Material | Kind::GaussianBeam | Kind::DrudeMaterial | Kind::DebyeMaterial => 3,
            Kind::LorentzMaterial => 4,
            Kind::PointCharge => 8,
            Kind::JonesVector | Kind::StokesVector | Kind::AbcdMatrix | Kind::ResonatorStability => 4,
            Kind::JonesMatrix => 8,
            Kind::MuellerMatrix => 16,
        }
    }

    /// Whether value k is a bool (0.0 false, 1.0 true in a value; `i:0` or `i:1` in a row): ResonatorStability's
    /// `is_stable`.
    fn is_bool(self, k: usize) -> bool {
        matches!(self, Kind::ResonatorStability) && k == 3
    }

    /// Whether value k is a usize (M6 part 1d), carried in the value's f64 as its bits: CpmlParams' `thickness`.
    fn is_usize(self, k: usize) -> bool {
        matches!(self, Kind::CpmlParams) && k == 0
    }

    /// Whether every f64 among the values v is finite (a bool or a usize always is).
    fn finite(self, v: &[f64]) -> bool {
        v[..self.width()].iter().enumerate().all(|(k, x)| self.is_usize(k) || x.is_finite())
    }

    /// The row tokens of the values v: `f:` each f64, `i:` each bool and each usize (its word).
    fn toks(self, v: &[f64]) -> Vec<Tok> {
        v.iter()
            .enumerate()
            .map(|(k, &x)| {
                if self.is_bool(k) {
                    b(x != 0.0)
                } else if self.is_usize(k) {
                    usize_tok(word(x))
                } else {
                    f(x)
                }
            })
            .collect()
    }

    /// The spellings of the values v: each f64 by `g`, each bool as `true` or `false`, each usize in decimal.
    fn spell(self, v: &[f64], g: impl Fn(f64) -> String) -> Vec<String> {
        v.iter()
            .enumerate()
            .map(|(k, &x)| {
                if self.is_bool(k) {
                    (x != 0.0).to_string()
                } else if self.is_usize(k) {
                    word(x).to_string()
                } else {
                    g(x)
                }
            })
            .collect()
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
            Kind::GaussianBeam => {
                format!("{{\"wavelength\":{},\"waist\":{},\"rayleigh_range\":{}}}", v[0], v[1], v[2])
            }
            Kind::AbcdMatrix => format!("{{\"a\":{},\"b\":{},\"c\":{},\"d\":{}}}", v[0], v[1], v[2], v[3]),
            Kind::ResonatorStability => format!(
                "{{\"g1\":{},\"g2\":{},\"g_product\":{},\"is_stable\":{}}}",
                v[0], v[1], v[2], v[3]
            ),
            Kind::LossyMaterial => format!("{{\"eps_r\":{},\"conductivity\":{}}}", v[0], v[1]),
            Kind::DrudeMaterial => format!("{{\"eps_inf\":{},\"omega_p\":{},\"gamma\":{}}}", v[0], v[1], v[2]),
            Kind::LorentzMaterial => format!(
                "{{\"eps_inf\":{},\"delta_eps\":{},\"omega_0\":{},\"gamma\":{}}}",
                v[0], v[1], v[2], v[3]
            ),
            Kind::DebyeMaterial => format!("{{\"eps_inf\":{},\"eps_s\":{},\"tau\":{}}}", v[0], v[1], v[2]),
            Kind::CpmlParams => format!(
                "{{\"thickness\":{},\"grading_order\":{},\"sigma_ratio\":{},\"kappa_max\":{},\"alpha_max\":{}}}",
                v[0], v[1], v[2], v[3], v[4]
            ),
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
            Kind::GaussianBeam => format!("[{},{},{}]", v[0], v[1], v[2]),
            Kind::AbcdMatrix | Kind::ResonatorStability | Kind::LorentzMaterial => {
                format!("[{},{},{},{}]", v[0], v[1], v[2], v[3])
            }
            Kind::LossyMaterial => format!("[{},{}]", v[0], v[1]),
            Kind::DrudeMaterial | Kind::DebyeMaterial => format!("[{},{},{}]", v[0], v[1], v[2]),
            Kind::CpmlParams => format!("[{},{},{},{},{}]", v[0], v[1], v[2], v[3], v[4]),
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
            Kind::GaussianBeam => {
                vec![("wavelength", v[0].clone()), ("waist", v[1].clone()), ("rayleigh_range", v[2].clone())]
            }
            Kind::AbcdMatrix => {
                vec![("a", v[0].clone()), ("b", v[1].clone()), ("c", v[2].clone()), ("d", v[3].clone())]
            }
            Kind::ResonatorStability => vec![
                ("g1", v[0].clone()),
                ("g2", v[1].clone()),
                ("g_product", v[2].clone()),
                ("is_stable", v[3].clone()),
            ],
            Kind::LossyMaterial => vec![("eps_r", v[0].clone()), ("conductivity", v[1].clone())],
            Kind::DrudeMaterial => {
                vec![("eps_inf", v[0].clone()), ("omega_p", v[1].clone()), ("gamma", v[2].clone())]
            }
            Kind::LorentzMaterial => vec![
                ("eps_inf", v[0].clone()),
                ("delta_eps", v[1].clone()),
                ("omega_0", v[2].clone()),
                ("gamma", v[3].clone()),
            ],
            Kind::DebyeMaterial => vec![("eps_inf", v[0].clone()), ("eps_s", v[1].clone()), ("tau", v[2].clone())],
            Kind::CpmlParams => vec![
                ("thickness", v[0].clone()),
                ("grading_order", v[1].clone()),
                ("sigma_ratio", v[2].clone()),
                ("kappa_max", v[3].clone()),
                ("alpha_max", v[4].clone()),
            ],
        }
    }

    /// The value, from its `width()` values (a bool as 0.0 or 1.0). serde's derive never calls a constructor, so every
    /// field is set as given.
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
            Kind::GaussianBeam => serde_json::to_string(&GaussianBeam {
                wavelength: black_box(v[0]),
                waist: black_box(v[1]),
                rayleigh_range: black_box(v[2]),
            }),
            Kind::AbcdMatrix => serde_json::to_string(&AbcdMatrix::new(
                black_box(v[0]),
                black_box(v[1]),
                black_box(v[2]),
                black_box(v[3]),
            )),
            Kind::ResonatorStability => serde_json::to_string(&ResonatorStability {
                g1: black_box(v[0]),
                g2: black_box(v[1]),
                g_product: black_box(v[2]),
                is_stable: black_box(v[3] != 0.0),
            }),
            Kind::LossyMaterial => serde_json::to_string(&LossyMaterial::new(black_box(v[0]), black_box(v[1]))),
            Kind::DrudeMaterial => {
                serde_json::to_string(&DrudeMaterial::new(black_box(v[0]), black_box(v[1]), black_box(v[2])))
            }
            Kind::LorentzMaterial => serde_json::to_string(&LorentzMaterial::new(
                black_box(v[0]),
                black_box(v[1]),
                black_box(v[2]),
                black_box(v[3]),
            )),
            Kind::DebyeMaterial => {
                serde_json::to_string(&DebyeMaterial::new(black_box(v[0]), black_box(v[1]), black_box(v[2])))
            }
            Kind::CpmlParams => {
                let mut p = CpmlParams::new(black_box(word(v[0])));
                p.grading_order = black_box(v[1]);
                p.sigma_ratio = black_box(v[2]);
                p.kappa_max = black_box(v[3]);
                p.alpha_max = black_box(v[4]);
                serde_json::to_string(&p)
            }
        }
        .expect("serde_json writes every value")
    }

    /// serde_json's reading of `doc`: the values in field order (a bool as 0.0 or 1.0), or its error's message.
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
            Kind::GaussianBeam => serde_json::from_str::<GaussianBeam>(doc)
                .map(|g| vec![g.wavelength, g.waist, g.rayleigh_range]),
            Kind::AbcdMatrix => serde_json::from_str::<AbcdMatrix>(doc).map(|m| vec![m.a, m.b, m.c, m.d]),
            Kind::ResonatorStability => serde_json::from_str::<ResonatorStability>(doc)
                .map(|r| vec![r.g1, r.g2, r.g_product, f64::from(u8::from(r.is_stable))]),
            Kind::LossyMaterial => serde_json::from_str::<LossyMaterial>(doc).map(|m| vec![m.eps_r, m.conductivity]),
            Kind::DrudeMaterial => {
                serde_json::from_str::<DrudeMaterial>(doc).map(|m| vec![m.eps_inf, m.omega_p, m.gamma])
            }
            Kind::LorentzMaterial => serde_json::from_str::<LorentzMaterial>(doc)
                .map(|m| vec![m.eps_inf, m.delta_eps, m.omega_0, m.gamma]),
            Kind::DebyeMaterial => serde_json::from_str::<DebyeMaterial>(doc).map(|m| vec![m.eps_inf, m.eps_s, m.tau]),
            Kind::CpmlParams => serde_json::from_str::<CpmlParams>(doc)
                .map(|p| vec![carried(p.thickness), p.grading_order, p.sigma_ratio, p.kappa_max, p.alpha_max]),
        }
        .map_err(|e| e.to_string())
    }
}

/// A usize carried in a value's f64 (`Kind::is_usize`): its bits.
fn word(x: f64) -> usize {
    usize::try_from(x.to_bits()).expect("a 64-bit usize")
}

/// A usize as a value's f64: its word as the bits.
fn carried(n: usize) -> f64 {
    f64::from_bits(u64::try_from(n).expect("a 64-bit usize"))
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

/// Spellings for a bool field (M3 part 1c-2): serde reads only the literals `true` and `false`, and refuses every
/// other JSON value (numbers 0 and 1, strings, null, containers) and every near miss of the literals.
fn bool_spellings() -> Vec<String> {
    [
        "true", "false", " true ", "\t\nfalse\r ", "True", "TRUE", "False", "FALSE", "tru", "t", "f", "fals", "truee",
        "true1", "true0", "true.0", "falsee", "false0", "trUe", "\\u0074rue", "1", "0", "-0", "1.0", "0.0", "1e0", "-1",
        "2", "01", "1e400", "18446744073709551616", "-", "null", "nul", "\"true\"", "\"false\"", "\"\"",
        "\"\\u0074rue\"", "\"\\ud800\"", "[]", "[true]", "[false]", "{}", "{\"x\":true}", "yes", "on", "NaN",
        "Infinity",
    ]
    .iter()
    .map(|x| (*x).to_string())
    .collect()
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

// The notes stay short: a row's inputs and outputs already carry its document and values, and
// json.tsv must stay under GitHub's 50 MB file advisory (main.rs checks it).
fn writer_row(m: &mut Module, kind: Kind, v: &[f64]) {
    let text = kind.write(v);
    let ins = kind.toks(&v[..kind.width()]);
    m.row_with_note(&format!("json::{}::to_string", kind.name()), Class::A, &ins, &tokens(&text), &preview(&text));
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
    let head = format!("unchecked {d}");
    match read {
        Ok(v) => m.row_with_note(&func, Class::C, &ins, &kind.toks(&v), &format!("{head} -> Ok")),
        Err(msg) => {
            let outs = [Tok::Err("Deserialize".to_string())];
            m.row_with_note(&func, Class::C, &ins, &outs, &format!("{head} -> Err({msg})"));
        }
    }
}

/// The plain spellings of a set of seeded values (`Kind::draws`) for a document's other fields; a bool field is `true`,
/// a usize field `10`.
fn fillers(kind: Kind, rng: &mut Rng) -> Vec<String> {
    let mut v: Vec<String> = (0..kind.draws()).map(|_| js(rng.signed_log_uniform(-20.0, 20.0))).collect();
    for (k, x) in v.iter_mut().enumerate().take(kind.width()) {
        if kind.is_bool(k) {
            *x = "true".to_string();
        } else if kind.is_usize(k) {
            *x = "10".to_string();
        }
    }
    v
}

/// Unknown values of every JSON kind (deep, long, escaped, lone surrogates, refused spellings), for unknown keys.
fn unknown_values() -> Vec<String> {
    vec![
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
    ]
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
    let unknown_values = unknown_values();
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
        Kind::Complex
        | Kind::FieldVector
        | Kind::PointCharge
        | Kind::Material
        | Kind::GaussianBeam
        | Kind::AbcdMatrix
        | Kind::ResonatorStability
        | Kind::LossyMaterial
        | Kind::DrudeMaterial
        | Kind::LorentzMaterial
        | Kind::DebyeMaterial
        | Kind::CpmlParams => {}
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
    for kind in Kind::ALL {
        emit_kind(m, kind, &mut rng, &xs);
    }
    emit_mie(m, &mut rng, &xs);
    // M3 part 1c-2: beam's three structs, after MieResult, so the earlier rows and their seeded stream are unchanged.
    for kind in Kind::BEAM {
        emit_kind(m, kind, &mut rng, &xs);
    }
    // M4 part 1d: rf's three types, after beam's.
    emit_rf(m, &mut rng, &xs);
    // M5 part 1f: the FDTD materials, then the grid and the eigenmode result, after rf's.
    for kind in Kind::FDTD {
        emit_kind(m, kind, &mut rng, &xs);
    }
    emit_fdtd_vecs(m, &mut rng, &xs);
    // M6 part 1d: the 2D FDTD's Mode2d and CpmlParams, after the vectors.
    emit_fdtd2d(m, &mut rng, &xs);
    rows.m
}

/// One type's writer and reader rows: the value sets, serde's documents and the other notations, every spelling in
/// each position, seeded decimals, the shapes and the nested forms.
fn emit_kind(m: &mut Rows, kind: Kind, rng: &mut Rng, xs: &[f64]) {
    let n = xs.len();
    // Writer rows: each value, beside others from the list.
    let w = kind.width();
    let mut sets: Vec<Vec<f64>> = Vec::new();
    for j in 0..n {
        sets.push((0..kind.draws()).map(|k| xs[(j + k * 37) % n]).collect());
    }
    for _ in 0..300 {
        sets.push((0..kind.draws()).map(|_| f64::from_bits(rng.bits())).collect());
    }
    // A bool field alternates false, true over the sets; a usize field takes THICKNESSES once each, then 0 to 24.
    for (j, v) in sets.iter_mut().enumerate() {
        for (k, x) in v.iter_mut().enumerate().take(w) {
            if kind.is_bool(k) {
                *x = f64::from(u8::from(j % 2 == 1));
            } else if kind.is_usize(k) {
                *x = carried(THICKNESSES.get(j).copied().unwrap_or(j % 25));
            }
        }
    }
    for v in &sets {
        m.writer(kind, v);
    }
    // Reader rows: serde's own text, and the same values in other notations (finite ones: `null` is refused, and
    // the shape rows cover it).
    for v in &sets {
        m.reader(kind, &kind.write(v));
        if kind.finite(v) {
            m.reader(kind, &kind.doc(&kind.spell(v, js)));
        }
    }
    for v in sets.iter().take(250) {
        if kind.finite(v) {
            m.reader(kind, &kind.doc(&kind.spell(v, |x| format!("{x}"))));
            m.reader(kind, &kind.doc(&kind.spell(v, |x| format!("{x:e}"))));
        }
    }
    // Every spelling in each field position: a number's in an f64 position, a bool's in a bool position.
    let spell = spellings();
    let bool_spell = bool_spellings();
    let base = fillers(kind, rng);
    for pos in 0..w {
        for sp in if kind.is_bool(pos) { &bool_spell } else { &spell } {
            let mut v = base.clone();
            v[pos] = sp.clone();
            m.reader(kind, &kind.doc(&v));
        }
    }
    // Seeded decimals, a few in each f64 position.
    let nums: Vec<usize> = (0..w).filter(|&k| !kind.is_bool(k) && !kind.is_usize(k)).collect();
    for (k, d) in random_decimals(rng, 1200).into_iter().enumerate() {
        let mut v = base.clone();
        v[nums[k % nums.len()]] = d;
        m.reader(kind, &kind.doc(&v));
    }
    for doc in shapes(kind, rng) {
        m.reader(kind, &doc);
    }
    for doc in nested(kind, &base) {
        m.reader(kind, &doc);
    }
}

// ── MieResult (M3 part 1b-2) ──────────────────────────────────────────────────────────

/// A MieResult as its row tokens: the four f64, then each vector as `i:` its length and its parts.
fn mie_toks(r: &MieResult) -> Vec<Tok> {
    let mut out = vec![f(r.q_sca), f(r.q_ext), f(r.q_abs), f(r.g)];
    for v in [&r.a, &r.b] {
        out.push(i(i64::try_from(v.len()).expect("a length fits")));
        for z in v {
            out.extend([f(z.re), f(z.im)]);
        }
    }
    out
}

/// A MieResult from four f64 and the parts of its vectors.
fn mie_value(q: &[f64], a: &[f64], b: &[f64]) -> MieResult {
    let cv = |p: &[f64]| p.chunks(2).map(|c| Complex::new(black_box(c[0]), black_box(c[1]))).collect::<Vec<_>>();
    MieResult {
        q_sca: black_box(q[0]),
        q_ext: black_box(q[1]),
        q_abs: black_box(q[2]),
        g: black_box(q[3]),
        a: cv(a),
        b: cv(b),
    }
}

/// Its document with the spellings q (four) and the vectors' spellings a and b (re, im pairs), each element written
/// by `c`.
fn mie_doc(q: &[String], a: &[String], b: &[String], c: fn(&str, &str) -> String) -> String {
    let vec_of = |p: &[String]| arr(&p.chunks(2).map(|z| c(&z[0], &z[1])).collect::<Vec<_>>());
    format!(
        "{{\"q_sca\":{},\"q_ext\":{},\"q_abs\":{},\"g\":{},\"a\":{},\"b\":{}}}",
        q[0], q[1], q[2], q[3], vec_of(a), vec_of(b)
    )
}

impl Rows {
    fn mie_writer(&mut self, r: &MieResult) {
        let ins = mie_toks(r);
        if !self.seen.insert(format!("wMieResult{ins:?}")) {
            return;
        }
        let text = serde_json::to_string(r).expect("serde_json writes every MieResult");
        self.m.row_with_note("json::MieResult::to_string", Class::A, &ins, &tokens(&text), &preview(&text));
    }

    fn mie_reader(&mut self, doc: &str) {
        if !self.seen.insert(format!("rMieResult{doc}")) {
            return;
        }
        let func = "json::MieResult::from_str";
        let read = serde_json::from_str::<MieResult>(black_box(doc));
        let d = if read.is_ok() { unchecked(doc) } else { 0 };
        let ins = [s(doc), i(d)];
        let head = format!("unchecked {d}");
        match read {
            Ok(r) => self.m.row_with_note(func, Class::C, &ins, &mie_toks(&r), &format!("{head} -> Ok")),
            Err(e) => {
                let outs = [Tok::Err("Deserialize".to_string())];
                self.m.row_with_note(func, Class::C, &ins, &outs, &format!("{head} -> Err({e})"));
            }
        }
    }
}

/// MieResult's writer and reader rows: real results, seeded values with every pair of lengths to 7, the specials,
/// 500-element vectors; then serde's documents and the shapes.
fn emit_mie(m: &mut Rows, rng: &mut Rng, xs: &[f64]) {
    let n = xs.len();
    let mut values: Vec<MieResult> = Vec::new();
    for &x in &[1e-3, 0.1, 1.0, 5.0, 4.0 * std::f64::consts::PI, 50.0, 200.0] {
        for &(re, im) in &[(1.5, 0.0), (1.33, -0.1), (0.2, -3.0)] {
            values.push(mie(x, Complex::new(re, im)).expect("the sample spheres compute"));
        }
    }
    let lens: Vec<(usize, usize)> = (0..8).flat_map(|na| (0..8).map(move |nb| (na, nb))).collect();
    for j in 0..n {
        let pick = |k: usize| xs[(j + k * 37) % n];
        let (na, nb) = lens[j % lens.len()];
        let q: Vec<f64> = (0..4).map(pick).collect();
        let a: Vec<f64> = (0..2 * na).map(|k| pick(4 + k)).collect();
        let b: Vec<f64> = (0..2 * nb).map(|k| pick(4 + 2 * na + k)).collect();
        values.push(mie_value(&q, &a, &b));
    }
    for _ in 0..120 {
        let (na, nb) = lens[usize::try_from(rng.int(0, 63)).expect("an index")];
        let mut bits = |k: usize| (0..k).map(|_| f64::from_bits(rng.bits())).collect::<Vec<_>>();
        let q = bits(4);
        let a = bits(2 * na);
        let b = bits(2 * nb);
        values.push(mie_value(&q, &a, &b));
    }
    let long: Vec<f64> = (0..1000).map(|k| xs[(k * 7) % n]).collect();
    values.push(mie_value(&long[..4], &long[..1000], &[]));
    values.push(mie_value(&long[..4], &long[..2], &long[..1000]));
    for r in &values {
        m.mie_writer(r);
    }
    // Reader rows: serde's own text of each value, and the finite ones in bayan's notation, Rust's `{}` and `{:e}`.
    for (k, r) in values.iter().enumerate() {
        m.mie_reader(&serde_json::to_string(r).expect("serde_json writes every MieResult"));
        let toks = mie_toks(r);
        let fs: Vec<f64> = toks.iter().filter_map(|t| if let Tok::F(x) = t { Some(*x) } else { None }).collect();
        if !fs.iter().all(|x| x.is_finite()) {
            continue;
        }
        let na = r.a.len();
        let spell = |g: fn(f64) -> String| -> (Vec<String>, Vec<String>, Vec<String>) {
            let all: Vec<String> = fs.iter().map(|&x| g(x)).collect();
            (all[..4].to_vec(), all[4..4 + 2 * na].to_vec(), all[4 + 2 * na..].to_vec())
        };
        let (q, a, b) = spell(js);
        m.mie_reader(&mie_doc(&q, &a, &b, cx));
        if k < 200 {
            let (q, a, b) = spell(|x| format!("{x}"));
            m.mie_reader(&mie_doc(&q, &a, &b, cx));
            let (q, a, b) = spell(|x| format!("{x:e}"));
            m.mie_reader(&mie_doc(&q, &a, &b, cxs));
        }
    }
    // Every spelling in each of a two-and-two document's twelve number positions; seeded decimals.
    let base: Vec<String> = (0..12).map(|_| js(rng.signed_log_uniform(-20.0, 20.0))).collect();
    let at = |v: &[String]| mie_doc(&v[..4], &v[4..8], &v[8..12], cx);
    for pos in 0..12 {
        for sp in &spellings() {
            let mut v = base.clone();
            v[pos] = sp.clone();
            m.mie_reader(&at(&v));
        }
    }
    for (k, d) in random_decimals(rng, 600).into_iter().enumerate() {
        let mut v = base.clone();
        v[k % 12] = d;
        m.mie_reader(&at(&v));
    }
    for doc in mie_shapes(&base) {
        m.mie_reader(&doc);
    }
}

/// MieResult's shapes on the twelve spellings v: serde's two forms at every level, field order, whitespace, unknown
/// keys, missing, repeated and wrong-typed fields, vectors of every length and element form, a vector as an object,
/// deep unknown values inside an element, every prefix, trailing data.
fn mie_shapes(v: &[String]) -> Vec<String> {
    let (q, a, b) = (&v[..4], &v[4..8], &v[8..12]);
    let base = mie_doc(q, a, b, cx);
    let inner = base[1..base.len() - 1].to_string();
    let av = arr(&a.chunks(2).map(|z| cx(&z[0], &z[1])).collect::<Vec<_>>());
    let bv = arr(&b.chunks(2).map(|z| cx(&z[0], &z[1])).collect::<Vec<_>>());
    let as_ = arr(&a.chunks(2).map(|z| cxs(&z[0], &z[1])).collect::<Vec<_>>());
    let members = [
        format!("\"q_sca\":{}", q[0]),
        format!("\"q_ext\":{}", q[1]),
        format!("\"q_abs\":{}", q[2]),
        format!("\"g\":{}", q[3]),
        format!("\"a\":{av}"),
        format!("\"b\":{bv}"),
    ];
    let parts: Vec<String> = members.iter().map(|m| m.split_once(':').expect("a member").1.to_string()).collect();
    let with = |k: usize, val: &str| {
        let mut ms = members.to_vec();
        ms[k] = format!("{}:{val}", members[k].split_once(':').expect("a member").0);
        format!("{{{}}}", ms.join(","))
    };
    let mut docs = vec![
        base.clone(),
        mie_doc(q, a, b, cxs),
        format!("[{}]", parts.join(",")),
        format!("[{},{},{},{},{as_},{bv}]", q[0], q[1], q[2], q[3]),
        {
            let mut r = members.to_vec();
            r.reverse();
            format!("{{{}}}", r.join(","))
        },
        format!(" \t\n\r{{ {} }} \n", members.join(" ,\n ")),
        mie_doc(q, &[], &[], cx),
        mie_doc(q, &[], b, cx),
        mie_doc(q, a, &[], cx),
    ];
    // Unknown keys of every kind, before and after the fields.
    for u in &unknown_values() {
        docs.push(format!("{{\"w\":{u},{inner}}}"));
        docs.push(format!("{{{inner},\"unknown\":{u}}}"));
    }
    // Missing fields, each; repeated ones; the generic refusals.
    for k in 0..members.len() {
        let mut ms = members.to_vec();
        ms.remove(k);
        docs.push(format!("{{{}}}", ms.join(",")));
        docs.push(format!("{{{inner},{}}}", members[k]));
    }
    for d in ["{}", "[]", "null", "true", "1", "\"x\"", "", " ", "{", "}", "[", "]"] {
        docs.push(d.to_string());
    }
    // Wrong types: each f64 field, then each vector.
    for k in 0..4 {
        for bad in ["null", "true", "\"1\"", "{}", "[]", "[1]"] {
            docs.push(with(k, bad));
        }
    }
    let one = cx(&a[0], &a[1]);
    for k in 4..6 {
        for bad in [
            "null".to_string(),
            "true".to_string(),
            "1".to_string(),
            "\"1\"".to_string(),
            "{}".to_string(),
            format!("{{\"0\":{one}}}"),
            "[1]".to_string(),
            "[null]".to_string(),
            "[[1,2,3]]".to_string(),
            "[[1]]".to_string(),
            "[[]]".to_string(),
            "[{}]".to_string(),
            "[[1,2]]".to_string(),
            "[[1,2],{\"re\":3,\"im\":4}]".to_string(),
            format!("[{one},]"),
            format!("[{one},null]"),
            format!("[[{one}]]"),
        ] {
            docs.push(with(k, &bad));
        }
        for f in complex_forms(&a[0], &a[1]) {
            docs.push(with(k, &format!("[{f}]")));
            docs.push(with(k, &format!("[{one},{f}]")));
        }
    }
    // Vectors of many lengths, and long ones.
    for len in [3usize, 17, 128, 300] {
        let many = arr(&(0..len).map(|j| cx(&format!("{j}"), &format!("-{j}.5"))).collect::<Vec<_>>());
        docs.push(with(4, &many));
        docs.push(with(5, &many));
    }
    // A deep unknown value inside an element: the element is at depth 3 (the document, `a`, the element).
    for depth in [124usize, 125, 126, 127, 128, 129] {
        let d = depth - 3;
        let deep = format!("{{\"w\":{}0{},\"re\":1,\"im\":2}}", "[".repeat(d), "]".repeat(d));
        docs.push(with(4, &format!("[{deep}]")));
    }
    // Every prefix, and trailing data.
    for k in 0..base.len() {
        if base.is_char_boundary(k) {
            docs.push(base[..k].to_string());
        }
    }
    for tail in [" ", "x", "}", "]", ",", "{}", "//"] {
        docs.push(format!("{base}{tail}"));
    }
    // The sequence form short, long, with a trailing comma.
    docs.push(format!("[{}]", parts[..5].join(",")));
    docs.push(format!("[{},1]", parts.join(",")));
    docs.push(format!("[{},]", parts.join(",")));
    docs
}

// ── rf's three types (M4 part 1d) ──────────────────────────────────────────────────────

/// A `usize` as its row token (its 64-bit word).
fn usize_tok(n: usize) -> Tok {
    i(i64::from_ne_bytes(n.to_ne_bytes()))
}

/// A length as its row token.
fn count_tok(n: usize) -> Tok {
    i(i64::try_from(n).expect("a length fits"))
}

/// An SMatrix as its row tokens (validate.py's JSON_FIELDS layout): `i:` ports, `i:` the row count, each row as `i:`
/// its length then its entries' (re, im), `i:` the impedance count then the impedances. A read may be ragged.
fn smj_toks(m: &SMatrix) -> Vec<Tok> {
    let mut out = vec![usize_tok(m.ports), count_tok(m.s.len())];
    for row in &m.s {
        out.push(count_tok(row.len()));
        for z in row {
            out.extend([f(z.re), f(z.im)]);
        }
    }
    out.push(count_tok(m.z0.len()));
    out.extend(m.z0.iter().map(|&x| f(x)));
    out
}

/// A TouchstoneData as its row tokens: `i:` ports, `f:` z0, `i:` the frequency count and the frequencies, `i:` the
/// matrix count and each matrix (smj_toks).
fn tdj_toks(d: &TouchstoneData) -> Vec<Tok> {
    let mut out = vec![usize_tok(d.ports), f(d.z0), count_tok(d.frequencies.len())];
    out.extend(d.frequencies.iter().map(|&x| f(x)));
    out.push(count_tok(d.s_matrices.len()));
    for m in &d.s_matrices {
        out.extend(smj_toks(m));
    }
    out
}

const TS_FORMATS: [TouchstoneFormat; 3] =
    [TouchstoneFormat::MagnitudeAngle, TouchstoneFormat::RealImag, TouchstoneFormat::DbAngle];

/// A format as its row token: the port's enum value (0 MagnitudeAngle, 1 RealImag, 2 DbAngle).
fn fmt_tok(t: TouchstoneFormat) -> Tok {
    count_tok(TS_FORMATS.iter().position(|&x| x == t).expect("one of the three"))
}

/// A square ports-port matrix of the given entries (row-major (re, im)) and impedances.
fn smj(ports: usize, e: &[f64], z0: &[f64]) -> SMatrix {
    let s = (0..ports)
        .map(|r| (0..ports).map(|c| Complex::new(e[2 * (r * ports + c)], e[2 * (r * ports + c) + 1])).collect())
        .collect();
    SMatrix { ports, s, z0: z0.to_vec() }
}

impl Rows {
    fn rf_writer(&mut self, name: &str, ins: Vec<Tok>, text: &str) {
        if !self.seen.insert(format!("w{name}{ins:?}")) {
            return;
        }
        self.m.row_with_note(&format!("json::{name}::to_string"), Class::A, &ins, &tokens(text), &preview(text));
    }

    fn rf_reader<T>(&mut self, name: &str, doc: &str, read: Result<T, String>, toks: impl Fn(&T) -> Vec<Tok>) {
        if !self.seen.insert(format!("r{name}{doc}")) {
            return;
        }
        let func = format!("json::{name}::from_str");
        let d = if read.is_ok() { unchecked(doc) } else { 0 };
        let ins = [s(doc), i(d)];
        let head = format!("unchecked {d}");
        match read {
            Ok(v) => self.m.row_with_note(&func, Class::C, &ins, &toks(&v), &format!("{head} -> Ok")),
            Err(e) => {
                let outs = [Tok::Err("Deserialize".to_string())];
                self.m.row_with_note(&func, Class::C, &ins, &outs, &format!("{head} -> Err({e})"));
            }
        }
    }

    fn sm_json(&mut self, doc: &str) {
        let read = serde_json::from_str::<SMatrix>(black_box(doc)).map_err(|e| e.to_string());
        self.rf_reader("SMatrix", doc, read, smj_toks);
    }

    fn td_json(&mut self, doc: &str) {
        let read = serde_json::from_str::<TouchstoneData>(black_box(doc)).map_err(|e| e.to_string());
        self.rf_reader("TouchstoneData", doc, read, tdj_toks);
    }

    fn tf_json(&mut self, doc: &str) {
        let read = serde_json::from_str::<TouchstoneFormat>(black_box(doc)).map_err(|e| e.to_string());
        self.rf_reader("TouchstoneFormat", doc, read, |&t| vec![fmt_tok(t)]);
    }
}

/// A 2-port's document from its spellings: ports, the eight entry parts (row-major), the two impedances.
fn sm_doc(p: &str, e: &[String], z: &[String], c: fn(&str, &str) -> String) -> String {
    let row = |r: usize| arr(&[c(&e[4 * r], &e[4 * r + 1]), c(&e[4 * r + 2], &e[4 * r + 3])]);
    format!("{{\"ports\":{p},\"s\":{},\"z0\":{}}}", arr(&[row(0), row(1)]), arr(z))
}

/// rf's three types, after beam's: values from the shared list and seeded bits; serde's text of each, read back in
/// each notation; every spelling in each number position (`ports` among them, a usize); seeded decimals; and the shapes
/// a fixed matrix can't hold, which serde reads (ADR 0002 row 36): too few or too many rows, a ragged row, an impedance
/// list of another length, frequencies unlike the matrices in number, a matrix of another port count, a port count
/// past 2^63 − 1.
fn emit_rf(m: &mut Rows, rng: &mut Rng, xs: &[f64]) {
    let n = xs.len();
    // TouchstoneFormat: its three values, then every spelling of a variant and a near miss.
    for t in TS_FORMATS {
        let text = serde_json::to_string(&t).expect("serde_json writes a unit variant");
        m.rf_writer("TouchstoneFormat", vec![fmt_tok(t)], &text);
        m.tf_json(&text);
    }
    for d in [
        "{\"MagnitudeAngle\":null}", "{\"RealImag\":null}", "{\"DbAngle\":null}", " \"RealImag\" ",
        "{ \"DbAngle\" : null }", "\"\\u0052ealImag\"", "\"Real\\u0049mag\"", "{\"\\u0044bAngle\":null}",
        "{\"RealImag\":[]}", "{\"RealImag\":{}}", "{\"RealImag\":0}", "{\"RealImag\":false}", "{\"RealImag\":\"\"}",
        "{\"RealImag\":null,\"x\":1}", "{\"RealImag\":null,\"RealImag\":null}", "{\"x\":null}", "{}", "[\"RealImag\"]",
        "[]", "\"realimag\"", "\"REALIMAG\"", "\"Real Imag\"", "\"RealImag \"", "\" RealImag\"", "\"\"", "\"RI\"",
        "\"MA\"", "\"DB\"", "\"Magnitude\"", "\"\\ud800\"", "\"RealImag\\u0000\"", "0", "1", "2", "-0", "null",
        "true", "", " ", "\"RealImag", "RealImag", "\"RealImag\"x", "\"RealImag\",", "{\"RealImag\"}",
        "{\"RealImag\":}", "{\"RealImag\":null", "{\"RealImag\":null}}", "{\"RealImag\":nul}",
    ] {
        m.tf_json(d);
    }

    // SMatrix values: n-ports to 4 from the shared list, the Rust tests' 2-ports, seeded bits.
    let mut mats: Vec<SMatrix> = vec![SMatrix::new(0, 50.0), SMatrix::new(2, 50.0)];
    for j in 0..n {
        let ports = j % 5;
        let pick = |k: usize| xs[(j + k * 37) % n];
        let e: Vec<f64> = (0..2 * ports * ports).map(pick).collect();
        let z: Vec<f64> = (0..ports).map(|k| pick(2 * ports * ports + k)).collect();
        mats.push(smj(ports, &e, &z));
    }
    for _ in 0..60 {
        let ports = usize::try_from(rng.int(0, 4)).expect("0..=4");
        let e: Vec<f64> = (0..2 * ports * ports).map(|_| f64::from_bits(rng.bits())).collect();
        let z: Vec<f64> = (0..ports).map(|_| f64::from_bits(rng.bits())).collect();
        mats.push(smj(ports, &e, &z));
    }
    for (k, sm) in mats.iter().enumerate() {
        let text = serde_json::to_string(sm).expect("serde_json writes every SMatrix");
        m.rf_writer("SMatrix", smj_toks(sm), &text);
        m.sm_json(&text);
        let fs: Vec<f64> = sm.s.iter().flatten().flat_map(|z| [z.re, z.im]).chain(sm.z0.iter().copied()).collect();
        if k < 120 && fs.iter().all(|x| x.is_finite()) {
            let p = sm.ports;
            let spell = |g: &dyn Fn(f64) -> String, c: fn(&str, &str) -> String| {
                let rows: Vec<String> = sm
                    .s
                    .iter()
                    .map(|row| arr(&row.iter().map(|z| c(&g(z.re), &g(z.im))).collect::<Vec<_>>()))
                    .collect();
                let z: Vec<String> = sm.z0.iter().map(|&x| g(x)).collect();
                format!("{{\"ports\":{p},\"s\":{},\"z0\":{}}}", arr(&rows), arr(&z))
            };
            m.sm_json(&spell(&js, cx));
            m.sm_json(&spell(&|x| format!("{x}"), cx));
            m.sm_json(&spell(&|x| format!("{x:e}"), cxs));
        }
    }
    // Every spelling in each of a 2-port document's eleven number positions (ports, the eight parts, the two
    // impedances); seeded decimals in the f64 positions.
    let e: Vec<String> = (0..8).map(|_| js(rng.signed_log_uniform(-3.0, 1.0))).collect();
    let z: Vec<String> = vec!["50.0".to_string(), "50.0".to_string()];
    for sp in &spellings() {
        m.sm_json(&sm_doc(sp, &e, &z, cx));
        for k in 0..8 {
            let mut v = e.clone();
            v[k].clone_from(sp);
            m.sm_json(&sm_doc("2", &v, &z, cx));
        }
        for k in 0..2 {
            let mut v = z.clone();
            v[k].clone_from(sp);
            m.sm_json(&sm_doc("2", &e, &v, cx));
        }
    }
    for (k, d) in random_decimals(rng, 300).into_iter().enumerate() {
        let mut v = e.clone();
        v[k % 8].clone_from(&d);
        m.sm_json(&sm_doc("2", &v, &z, cx));
    }
    for doc in sm_shapes(&e) {
        m.sm_json(&doc);
    }

    // TouchstoneData values: 0- to 3-ports with 0 to 3 points from the shared list, impedances per matrix.
    let mut data: Vec<TouchstoneData> = Vec::new();
    for j in 0..n / 2 {
        let pick = |k: usize| xs[(j + k * 41) % n];
        let ports = j % 4;
        let points = (j / 4) % 4;
        let mut at = 0;
        let mut next = || {
            at += 1;
            pick(at)
        };
        let z0 = next();
        let frequencies: Vec<f64> = (0..points).map(|_| next()).collect();
        let s_matrices = (0..points)
            .map(|_| {
                let e: Vec<f64> = (0..2 * ports * ports).map(|_| next()).collect();
                let z: Vec<f64> = (0..ports).map(|_| next()).collect();
                smj(ports, &e, &z)
            })
            .collect();
        data.push(TouchstoneData { ports, z0, frequencies, s_matrices });
    }
    for (k, d) in data.iter().enumerate() {
        let text = serde_json::to_string(d).expect("serde_json writes every TouchstoneData");
        m.rf_writer("TouchstoneData", tdj_toks(d), &text);
        m.td_json(&text);
        let all = tdj_toks(d);
        let finite = all.iter().all(|t| if let Tok::F(x) = t { x.is_finite() } else { true });
        if k < 60 && finite {
            let doc = text.replace("\"z0\":", "\"z0\" : ").replace(',', " ,\n ");
            m.td_json(&doc);
        }
    }
    for doc in td_shapes(rng) {
        m.td_json(&doc);
    }
}

/// SMatrix's shapes on the 2-port spellings e: serde's two forms at every level, field order, unknown keys, missing,
/// repeated and wrong-typed fields, every row-36 shape serde reads, port counts at usize's edges, every prefix,
/// trailing data, deep unknown values.
fn sm_shapes(e: &[String]) -> Vec<String> {
    let z = ["50.0".to_string(), "75.0".to_string()];
    let base = sm_doc("2", e, &z, cx);
    let inner = base[1..base.len() - 1].to_string();
    let c0 = cx(&e[0], &e[1]);
    let row = |k: usize| format!("[{},{}]", cx(&e[4 * k], &e[4 * k + 1]), cx(&e[4 * k + 2], &e[4 * k + 3]));
    let s2 = format!("[{},{}]", row(0), row(1));
    let mut docs = vec![
        base.clone(),
        sm_doc("2", e, &z, cxs),
        format!("[2,{s2},[50.0,75.0]]"),
        format!("{{\"z0\":[50.0,75.0],\"s\":{s2},\"ports\":2}}"),
        format!(" \t\n{{ \"ports\" : 2 ,\n\"s\" : {s2} , \"z0\" : [ 50.0 , 75.0 ] }}\r\n"),
        "{\"ports\":0,\"s\":[],\"z0\":[]}".to_string(),
        "[0,[],[]]".to_string(),
        format!("{{\"ports\":1,\"s\":[[{c0}]],\"z0\":[50]}}"),
        format!("{{\"ports\":3,\"s\":[[{c0},{c0},{c0}],[{c0},{c0},{c0}],[{c0},{c0},{c0}]],\"z0\":[1,2,3]}}"),
    ];
    // Row 36: shapes serde reads and a fixed matrix can't hold.
    for d in [
        "{\"ports\":2,\"s\":[],\"z0\":[50,50]}".to_string(),
        format!("{{\"ports\":2,\"s\":[{}],\"z0\":[50,50]}}", row(0)),
        format!("{{\"ports\":2,\"s\":[{},{},{}],\"z0\":[50,50]}}", row(0), row(1), row(0)),
        format!("{{\"ports\":2,\"s\":[[{c0}],{}],\"z0\":[50,50]}}", row(1)),
        format!("{{\"ports\":2,\"s\":[{},[{c0},{c0},{c0}]],\"z0\":[50,50]}}", row(0)),
        "{\"ports\":2,\"s\":[[],[]],\"z0\":[50,50]}".to_string(),
        format!("{{\"ports\":2,\"s\":{s2},\"z0\":[]}}"),
        format!("{{\"ports\":2,\"s\":{s2},\"z0\":[50]}}"),
        format!("{{\"ports\":2,\"s\":{s2},\"z0\":[50,50,50]}}"),
        format!("{{\"ports\":1,\"s\":{s2},\"z0\":[50,50]}}"),
        format!("{{\"ports\":3,\"s\":{s2},\"z0\":[50,50]}}"),
        "{\"ports\":1,\"s\":[],\"z0\":[]}".to_string(),
        "{\"ports\":0,\"s\":[[]],\"z0\":[]}".to_string(),
        "{\"ports\":0,\"s\":[],\"z0\":[1]}".to_string(),
        "{\"ports\":4611686018427387904,\"s\":[],\"z0\":[]}".to_string(),
        "{\"ports\":9223372036854775807,\"s\":[],\"z0\":[]}".to_string(),
        "{\"ports\":9223372036854775808,\"s\":[],\"z0\":[]}".to_string(),
        "{\"ports\":18446744073709551615,\"s\":[],\"z0\":[]}".to_string(),
        "{\"ports\":18446744073709551616,\"s\":[],\"z0\":[]}".to_string(),
    ] {
        docs.push(d);
    }
    // Unknown keys of every kind, before and after the fields.
    for u in &unknown_values() {
        docs.push(format!("{{\"w\":{u},{inner}}}"));
        docs.push(format!("{{{inner},\"unknown\":{u}}}"));
    }
    // Missing fields, repeated ones, the generic refusals.
    let members = ["\"ports\":2".to_string(), format!("\"s\":{s2}"), "\"z0\":[50.0,75.0]".to_string()];
    for k in 0..3 {
        let mut ms = members.to_vec();
        ms.remove(k);
        docs.push(format!("{{{}}}", ms.join(",")));
        docs.push(format!("{{{inner},{}}}", members[k]));
    }
    for d in ["{}", "[]", "null", "true", "1", "\"x\"", "", " ", "{", "}", "[2,[],[]", "[2,[]]", "[2,[],[],1]"] {
        docs.push(d.to_string());
    }
    // Wrong types: ports (a usize), s at each level, z0.
    let with = |k: usize, val: &str| {
        let mut ms = members.to_vec();
        ms[k] = format!("{}:{val}", members[k].split_once(':').expect("a member").0);
        format!("{{{}}}", ms.join(","))
    };
    for bad in ["null", "true", "\"2\"", "{}", "[]", "[2]", "2.0", "2e0", "-2", "-0", "02", "2.5", "1e400"] {
        docs.push(with(0, bad));
    }
    for bad in [
        "null".to_string(), "1".to_string(), "{}".to_string(), "[null,null]".to_string(), "[1,2]".to_string(),
        format!("[{{}},{}]", row(1)), "[[1,2],[3,4]]".to_string(), "[[[1,2],[3,4]],[[5,6],[7,8]]]".to_string(),
        format!("[[{c0},null],{}]", row(1)), format!("[[{c0},[1]],{}]", row(1)), format!("[[{c0},[1,2,3]],{}]", row(1)),
        format!("[[{c0},{{\"re\":1}}],{}]", row(1)), format!("[[{c0},{{\"re\":1,\"im\":2,\"x\":3}}],{}]", row(1)),
        format!("[[{c0},{c0},],{}]", row(1)),
    ] {
        docs.push(with(1, &bad));
    }
    for bad in ["null", "50", "{}", "[null,1]", "[\"50\",50]", "[[50],50]", "[50,]", "[1e400,50]"] {
        docs.push(with(2, bad));
    }
    for f in complex_forms(&e[0], &e[1]) {
        docs.push(with(1, &format!("[[{f},{c0}],{}]", row(1))));
    }
    // A deep unknown value inside an element (at depth 4: the document, s, a row, the element).
    for depth in [124usize, 125, 126, 127, 128, 129] {
        let d = depth - 4;
        let deep = format!("{{\"w\":{}0{},\"re\":1,\"im\":2}}", "[".repeat(d), "]".repeat(d));
        docs.push(with(1, &format!("[[{deep},{c0}],{}]", row(1))));
    }
    // Every prefix, trailing data.
    for k in 0..base.len() {
        if base.is_char_boundary(k) {
            docs.push(base[..k].to_string());
        }
    }
    for tail in [" ", "x", "}", "]", ",", "{}"] {
        docs.push(format!("{base}{tail}"));
    }
    docs
}

/// TouchstoneData's documents beyond its writer's: a 1-port and a 2-port base, each field's wrong types, missing,
/// repeated and unknown fields, the sequence form, row 36's shapes (frequencies and matrices unlike in number, a
/// matrix of another port count or ragged, ports past 2^63 − 1 with no points), spellings in each number position,
/// every prefix of the 1-port base.
fn td_shapes(rng: &mut Rng) -> Vec<String> {
    let one = "{\"ports\":1,\"s\":[[{\"re\":0.5,\"im\":0.1}]],\"z0\":[50.0]}";
    let two = "{\"ports\":2,\"s\":[[{\"re\":0.1,\"im\":0.0},{\"re\":0.9,\"im\":0.0}],\
               [{\"re\":0.9,\"im\":0.0},{\"re\":0.1,\"im\":0.0}]],\"z0\":[50.0,50.0]}";
    let td = |p: &str, z: &str, f: &str, ms: &str| {
        format!("{{\"ports\":{p},\"z0\":{z},\"frequencies\":{f},\"s_matrices\":{ms}}}")
    };
    let base = td("1", "50.0", "[1000000000.0,2000000000.0]", &format!("[{one},{one}]"));
    let mut docs = vec![
        base.clone(),
        td("2", "50.0", "[1.0]", &format!("[{two}]")),
        td("0", "50.0", "[]", "[]"),
        td("3", "75.0", "[]", "[]"),
        format!("[1,50.0,[1e9,2e9],[{one},[1,[[[0.5,0.1]]],[50.0]]]]"),
        format!("{{\"s_matrices\":[{one}],\"frequencies\":[1],\"z0\":50,\"ports\":1}}"),
    ];
    // Row 36: shapes serde reads and fixed data can't hold.
    for d in [
        td("1", "50.0", "[1.0]", "[]"),
        td("1", "50.0", "[]", &format!("[{one}]")),
        td("1", "50.0", "[1.0,2.0]", &format!("[{one}]")),
        td("1", "50.0", "[1.0]", &format!("[{two}]")),
        td("2", "50.0", "[1.0,2.0]", &format!("[{two},{one}]")),
        td("1", "50.0", "[1.0]", "[{\"ports\":1,\"s\":[],\"z0\":[50.0]}]"),
        td("1", "50.0", "[1.0]", "[{\"ports\":1,\"s\":[[{\"re\":1,\"im\":0}]],\"z0\":[]}]"),
        td("0", "50.0", "[1.0]", "[{\"ports\":0,\"s\":[],\"z0\":[]}]"),
        td("4611686018427387904", "50.0", "[]", "[]"),
        td("9223372036854775807", "50.0", "[]", "[]"),
        td("9223372036854775808", "50.0", "[]", "[]"),
        td("18446744073709551615", "50.0", "[]", "[]"),
        td("18446744073709551616", "50.0", "[]", "[]"),
    ] {
        docs.push(d);
    }
    // Wrong types and missing, repeated, unknown fields.
    for bad in ["null", "\"1\"", "1.0", "-1", "-0", "[]", "{}"] {
        docs.push(td(bad, "50.0", "[1.0]", &format!("[{one}]")));
    }
    for bad in ["null", "\"50\"", "[]", "{}", "true"] {
        docs.push(td("1", bad, "[1.0]", &format!("[{one}]")));
    }
    for bad in ["null", "1.0", "{}", "[null]", "[\"1\"]", "[[1]]"] {
        docs.push(td("1", "50.0", bad, &format!("[{one}]")));
    }
    for bad in ["null", "{}", "[null]", "[1]", "[[]]", "[{}]", one] {
        docs.push(td("1", "50.0", "[1.0]", bad));
    }
    let inner = &base[1..base.len() - 1];
    for u in ["null", "[1,{\"a\":[]}]", "\"\\ud800\"", "1e400"] {
        docs.push(format!("{{\"w\":{u},{inner}}}"));
    }
    for k in ["ports", "z0", "frequencies", "s_matrices"] {
        docs.push(base.replacen(&format!("\"{k}\":"), "\"skip\":", 1));
        docs.push(format!("{{{inner},\"{k}\":null}}"));
    }
    // Spellings in each number position: ports, z0, a frequency, an entry.
    for sp in &spellings() {
        docs.push(td(sp, "50.0", "[]", "[]"));
        docs.push(td("1", sp, "[1.0]", &format!("[{one}]")));
        docs.push(td("1", "50.0", &format!("[{sp}]"), &format!("[{one}]")));
        docs.push(td("1", "50.0", "[1.0]", &format!("[{}]", one.replace("0.5", sp))));
    }
    for d in random_decimals(rng, 60) {
        docs.push(td("1", &d, "[1.0]", &format!("[{one}]")));
    }
    for k in 0..base.len() {
        if base.is_char_boundary(k) && k % 3 == 0 {
            docs.push(base[..k].to_string());
        }
    }
    docs
}

// ── the FDTD types with three Vec<f64> (M5 part 1f) ────────────────────────────────────────────────────────────────

/// NonUniformGrid (`dx`, `dy`, `dz`) or EigenmodeResult (`frequencies`, `quality_factors`, `amplitudes`): three
/// `Vec<f64>` fields. A value's row tokens are each vector as `i:` its length and its elements.
#[derive(Clone, Copy)]
enum Vecs {
    Grid,
    Eigen,
}

impl Vecs {
    fn name(self) -> &'static str {
        match self {
            Vecs::Grid => "NonUniformGrid",
            Vecs::Eigen => "EigenmodeResult",
        }
    }

    fn keys(self) -> [&'static str; 3] {
        match self {
            Vecs::Grid => ["dx", "dy", "dz"],
            Vecs::Eigen => ["frequencies", "quality_factors", "amplitudes"],
        }
    }

    fn write(self, v: &[Vec<f64>; 3]) -> String {
        let bx = |x: &Vec<f64>| x.iter().map(|&y| black_box(y)).collect::<Vec<f64>>();
        match self {
            Vecs::Grid => serde_json::to_string(&NonUniformGrid { dx: bx(&v[0]), dy: bx(&v[1]), dz: bx(&v[2]) }),
            Vecs::Eigen => serde_json::to_string(&EigenmodeResult {
                frequencies: bx(&v[0]),
                quality_factors: bx(&v[1]),
                amplitudes: bx(&v[2]),
            }),
        }
        .expect("serde_json writes every value")
    }

    fn read(self, doc: &str) -> Result<[Vec<f64>; 3], String> {
        let doc = black_box(doc);
        match self {
            Vecs::Grid => serde_json::from_str::<NonUniformGrid>(doc).map(|g| [g.dx, g.dy, g.dz]),
            Vecs::Eigen => serde_json::from_str::<EigenmodeResult>(doc)
                .map(|r| [r.frequencies, r.quality_factors, r.amplitudes]),
        }
        .map_err(|e| e.to_string())
    }

    /// Its document with the spellings of each vector, fields in declaration order.
    fn doc(self, v: &[Vec<String>; 3]) -> String {
        let k = self.keys();
        format!("{{\"{}\":{},\"{}\":{},\"{}\":{}}}", k[0], arr(&v[0]), k[1], arr(&v[1]), k[2], arr(&v[2]))
    }
}

fn vecs_toks(v: &[Vec<f64>; 3]) -> Vec<Tok> {
    let mut out = Vec::new();
    for x in v {
        out.push(i(i64::try_from(x.len()).expect("a length fits")));
        out.extend(x.iter().map(|&y| f(y)));
    }
    out
}

impl Rows {
    fn vecs_writer(&mut self, kind: Vecs, v: &[Vec<f64>; 3]) {
        let ins = vecs_toks(v);
        if !self.seen.insert(format!("w{}{ins:?}", kind.name())) {
            return;
        }
        let text = kind.write(v);
        let func = format!("json::{}::to_string", kind.name());
        self.m.row_with_note(&func, Class::A, &ins, &tokens(&text), &preview(&text));
    }

    fn vecs_reader(&mut self, kind: Vecs, doc: &str) {
        if !self.seen.insert(format!("r{}{doc}", kind.name())) {
            return;
        }
        let func = format!("json::{}::from_str", kind.name());
        let read = kind.read(doc);
        let d = if read.is_ok() { unchecked(doc) } else { 0 };
        let ins = [s(doc), i(d)];
        let head = format!("unchecked {d}");
        match read {
            Ok(v) => self.m.row_with_note(&func, Class::C, &ins, &vecs_toks(&v), &format!("{head} -> Ok")),
            Err(e) => {
                let outs = [Tok::Err("Deserialize".to_string())];
                self.m.row_with_note(&func, Class::C, &ins, &outs, &format!("{head} -> Err({e})"));
            }
        }
    }
}

/// The grid's and the eigenmode result's rows: values with every triple of lengths to 3 (the eigenmode result's equal
/// lengths and, read back, unequal ones: ADR 0002 row 36), seeded and special values, 300-element vectors; serde's
/// text and the other notations; every spelling in each position of a short document; seeded decimals; the shapes.
fn emit_fdtd_vecs(m: &mut Rows, rng: &mut Rng, xs: &[f64]) {
    let n = xs.len();
    for kind in [Vecs::Grid, Vecs::Eigen] {
        let mut values: Vec<[Vec<f64>; 3]> = Vec::new();
        let lens: Vec<[usize; 3]> = match kind {
            Vecs::Grid => (0..64).map(|j| [j % 4, (j / 4) % 4, j / 16]).collect(),
            Vecs::Eigen => (0..4).map(|j| [j, j, j]).collect(),
        };
        for j in 0..n {
            let ls = lens[j % lens.len()];
            let mut at = j;
            let mut pick = |k: usize| -> Vec<f64> {
                (0..k).map(|_| { at += 37; xs[at % n] }).collect()
            };
            values.push([pick(ls[0]), pick(ls[1]), pick(ls[2])]);
        }
        for _ in 0..120 {
            let ls = lens[usize::try_from(rng.int(0, i64::try_from(lens.len()).unwrap_or(1) - 1)).unwrap_or(0)];
            let mut bits = |k: usize| (0..k).map(|_| f64::from_bits(rng.bits())).collect::<Vec<_>>();
            let (a, b, c) = (bits(ls[0]), bits(ls[1]), bits(ls[2]));
            values.push([a, b, c]);
        }
        let long: Vec<f64> = (0..300).map(|k| xs[(k * 7) % n]).collect();
        // A value the port can build: the grid's axes any lengths, the eigenmode result's equal (ADR 0002 row 36).
        match kind {
            Vecs::Grid => values.push([long.clone(), long[..1].to_vec(), vec![]]),
            Vecs::Eigen => values.push([long[..2].to_vec(), long[2..4].to_vec(), long[4..6].to_vec()]),
        }
        values.push([long.clone(), long.clone(), long.clone()]);
        for v in &values {
            m.vecs_writer(kind, v);
        }
        for (k, v) in values.iter().enumerate() {
            m.vecs_reader(kind, &kind.write(v));
            if v.iter().flatten().all(|x| x.is_finite()) {
                let spell = |g: fn(f64) -> String| [0, 1, 2].map(|i2| v[i2].iter().map(|&x| g(x)).collect::<Vec<_>>());
                m.vecs_reader(kind, &kind.doc(&spell(js)));
                if k < 150 {
                    m.vecs_reader(kind, &kind.doc(&spell(|x| format!("{x}"))));
                    m.vecs_reader(kind, &kind.doc(&spell(|x| format!("{x:e}"))));
                }
            }
        }
        // Every spelling in each of a two-two-two document's positions; seeded decimals.
        let base: Vec<String> = (0..6).map(|_| js(rng.signed_log_uniform(-20.0, 20.0))).collect();
        let at = |v: &[String]| kind.doc(&[v[0..2].to_vec(), v[2..4].to_vec(), v[4..6].to_vec()]);
        for pos in 0..6 {
            for sp in &spellings() {
                let mut v = base.clone();
                v[pos] = sp.clone();
                m.vecs_reader(kind, &at(&v));
            }
        }
        for (k, d) in random_decimals(rng, 300).into_iter().enumerate() {
            let mut v = base.clone();
            v[k % 6] = d;
            m.vecs_reader(kind, &at(&v));
        }
        for doc in vecs_shapes(kind, &base) {
            m.vecs_reader(kind, &doc);
        }
    }
}

/// The shapes on six spellings: the sequence form at every level, field order, whitespace, unknown keys, missing,
/// repeated and wrong-typed fields, vectors of every length to 4 (unequal ones among them), a vector given as an
/// object, `null`, strings and nested arrays as elements, every prefix, trailing data.
fn vecs_shapes(kind: Vecs, v: &[String]) -> Vec<String> {
    let k = kind.keys();
    let (a, b, c) = (arr(&v[0..2]), arr(&v[2..4]), arr(&v[4..6]));
    let base = kind.doc(&[v[0..2].to_vec(), v[2..4].to_vec(), v[4..6].to_vec()]);
    let inner = base[1..base.len() - 1].to_string();
    let member = |j: usize, val: &str| format!("\"{}\":{val}", k[j]);
    let mut docs = vec![
        base.clone(),
        format!("[{a},{b},{c}]"),
        format!("{{{},{},{}}}", member(2, &c), member(0, &a), member(1, &b)),
        format!(" \t\n\r{{ {} , {} , {} }} \n", member(0, &a), member(1, &b), member(2, &c)),
        format!("{{{},{}}}", member(0, &a), member(1, &b)),
        format!("{{{},{}}}", member(1, &b), member(2, &c)),
        format!("{{{inner},{}}}", member(0, &a)),
        format!("{{{},{inner}}}", member(2, "[]")),
        format!("[{a},{b}]"),
        format!("[{a},{b},{c},[]]"),
        format!("[{a},{b},{c}] "),
        format!("{base}x"),
        format!("{base}{base}"),
        "{}".to_string(),
        "[]".to_string(),
        "null".to_string(),
        format!("{{\"{}\":{{}},\"{}\":{b},\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":[null],\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":[\"1\"],\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":[[1]],\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":1,\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":\"x\",\"{}\":{c}}}", k[0], k[1], k[2]),
        format!("{{\"{}\":{a},\"{}\":null,\"{}\":{c}}}", k[0], k[1], k[2]),
    ];
    for u in unknown_values().iter().take(12) {
        docs.push(format!("{{\"w\":{u},{inner}}}"));
        docs.push(format!("{{{inner},\"unknown\":{u}}}"));
    }
    // Every triple of lengths 0 to 3 (and one of 4): the eigenmode result's unequal ones serde reads, the port refuses.
    for l0 in 0..4 {
        for l1 in 0..4 {
            for l2 in 0..4 {
                let take = |l: usize, from: usize| arr(&(0..l).map(|q| v[(from + q) % 6].clone()).collect::<Vec<_>>());
                let (m0, m1, m2) = (member(0, &take(l0, 0)), member(1, &take(l1, 2)), member(2, &take(l2, 4)));
                docs.push(format!("{{{m0},{m1},{m2}}}"));
            }
        }
    }
    docs.push(format!("{{{},{},{}}}", member(0, &arr(&v[0..4])), member(1, &arr(&v[1..5])), member(2, &arr(&v[2..6]))));
    for cut in (1..base.len()).step_by(3) {
        if base.is_char_boundary(cut) {
            docs.push(base[..cut].to_string());
        }
    }
    docs
}


// ── the 2D FDTD types (M6 part 1d) ─────────────────────────────────────────────────────

/// CpmlParams' writer thicknesses at usize's edges, one set each before the small ones: zero (which `Cpml::with_params`
/// refuses), the typical layers, freelist's largest profile block (2^31/48 layers) and one past, 2^32, f64's integer
/// edge, and i64's and u64's (past 2^63 − 1 the port's reader refuses, ADR 0002 row 36).
const THICKNESSES: [usize; 16] = [
    0,
    1,
    8,
    10,
    20,
    44_739_242,
    44_739_243,
    1 << 32,
    (1 << 53) - 1,
    1 << 53,
    (1 << 53) + 1,
    (1 << 63) - 1,
    1 << 63,
    (1 << 63) + 1,
    usize::MAX - 1,
    usize::MAX,
];

const MODES: [Mode2d; 2] = [Mode2d::Tm, Mode2d::Te];

/// A mode as its row token: the port's enum value (0 Tm, 1 Te).
fn mode_tok(t: Mode2d) -> Tok {
    count_tok(MODES.iter().position(|&x| x == t).expect("one of the two"))
}

impl Rows {
    fn md_json(&mut self, doc: &str) {
        let read = serde_json::from_str::<Mode2d>(black_box(doc)).map_err(|e| e.to_string());
        self.rf_reader("Mode2d", doc, read, |&t| vec![mode_tok(t)]);
    }
}

/// The 2D FDTD's two types, after the FDTD vectors: Mode2d's two values, then every spelling of a variant and a near
/// miss (as TouchstoneFormat's); then CpmlParams, a `Kind` (`emit_kind`) whose thickness is a usize.
fn emit_fdtd2d(m: &mut Rows, rng: &mut Rng, xs: &[f64]) {
    for t in MODES {
        let text = serde_json::to_string(&t).expect("serde_json writes a unit variant");
        m.rf_writer("Mode2d", vec![mode_tok(t)], &text);
        m.md_json(&text);
    }
    for d in [
        "{\"Tm\":null}", "{\"Te\":null}", " \"Te\" ", "{ \"Tm\" : null }", "\"\\u0054m\"", "\"T\\u0065\"",
        "{\"\\u0054e\":null}", "{\"Te\":[]}", "{\"Te\":{}}", "{\"Te\":0}", "{\"Te\":false}", "{\"Te\":\"\"}",
        "{\"Te\":null,\"x\":1}", "{\"x\":1,\"Te\":null}", "{\"Te\":null,\"Te\":null}", "{\"Tm\":null,\"Te\":null}",
        "{\"x\":null}", "{}", "[\"Te\"]", "[]", "\"te\"", "\"TE\"", "\"tm\"", "\"TM\"", "\"T e\"", "\"Te \"",
        "\" Te\"", "\"\"", "\"T\"", "\"Tem\"", "\"Tx\"", "\"Tm\\u0000\"", "\"Mode2d::Te\"", "\"\\ud800\"",
        "{\"\\ud800\":null}", "0", "1", "2", "-0", "null", "true", "", " ", "\"Te", "Te", "\"Te\"x", "\"Te\",",
        "\"Tm\"\"Te\"", "{\"Te\"}", "{\"Te\":}", "{\"Te\":null", "{\"Te\":null}}", "{\"Te\":nul}",
    ] {
        m.md_json(d);
    }
    for kind in Kind::FDTD2D {
        emit_kind(m, kind, rng, xs);
    }
}
