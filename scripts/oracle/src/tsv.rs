//! Golden TSV, format 1. The spec is in `scripts/oracle/README.md`; keep the two in step.
//! The format is tanmatra's, as kimiya and impetus use it, unchanged: every port's harness reads
//! the same files.
//!
//! One row per call:
//!
//! ```text
//! fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
//! ```
//!
//! `inputs` and `outputs` are space-separated typed tokens (`-` when empty). `note` is a
//! human-readable rendering of the same call, for review only; the harness ignores it. A complex
//! value is two `f:` tokens, its real part first (port plan §6 step 5).

/// Parity tolerance class of a row (port plan §7).
#[allow(dead_code)] // B, D and T arrive with the modules that need them
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    /// Bit-identical (NaN matches any NaN), or equal integers / strings.
    A,
    /// Arithmetic: ≤ 2 ulp, expected 0.
    B,
    /// Transcendental: 0 ulp against the corrected oracle (kimiya ADR 0006).
    C,
    /// Iterative: relative ≤ 1e-12 (expected to stay empty in bijli).
    D,
    /// An FDTD trajectory frame: a digest of each field component (M5).
    T,
}

impl Class {
    fn letter(self) -> char {
        match self {
            Self::A => 'A',
            Self::B => 'B',
            Self::C => 'C',
            Self::D => 'D',
            Self::T => 'T',
        }
    }
}

/// One typed value.
#[allow(dead_code)] // E, S, None, Err and V arrive with the first module that returns one
#[derive(Clone, Debug)]
pub enum Tok {
    /// `f:<16 hex digits>` — the IEEE-754 bit pattern, lowercase.
    F(f64),
    /// `i:<decimal>` — any integer; booleans are 0 / 1.
    I(i64),
    /// `e:<Enum>.<Variant>` — a Rust enum variant by name.
    E(&'static str, String),
    /// `s:<percent-encoded UTF-8>`.
    S(String),
    /// `n:` — `Option::None`.
    None,
    /// `E:<Variant>` — `Err(BijliError::<Variant>)`. The payload is compared through its own
    /// tokens where ADR 0002 row 4 maps it to the last-error record.
    Err(String),
    /// `v:<len>` — a sequence; its elements' tokens follow, in order.
    V(usize),
}

pub fn f(x: f64) -> Tok {
    Tok::F(x)
}

pub fn i(x: impl Into<i64>) -> Tok {
    Tok::I(x.into())
}

pub fn b(x: bool) -> Tok {
    Tok::I(i64::from(x))
}

pub fn s(x: impl Into<String>) -> Tok {
    Tok::S(x.into())
}

/// `e:<Enum>.<Variant>`, the variant named by its `Debug` text.
#[allow(dead_code)] // first used by a module with a public enum (M2: MagneticType)
pub fn e(name: &'static str, v: impl std::fmt::Debug) -> Tok {
    Tok::E(name, format!("{v:?}"))
}

/// `E:<Variant>` for a `BijliError`: its `Debug` text up to the payload, so
/// `InvalidParameter { reason: ".." }` is `E:InvalidParameter`.
#[allow(dead_code)] // first used by the M1 modules
pub fn err(x: &bijli::BijliError) -> Tok {
    let text = format!("{x:?}");
    Tok::Err(
        text.split(['(', ' ', '{'])
            .next()
            .unwrap_or_default()
            .to_string(),
    )
}

/// Percent-encode every byte that is not printable ASCII, plus `%`.
fn pct(x: &str) -> String {
    let mut out = String::with_capacity(x.len());
    for &byte in x.as_bytes() {
        if (0x21..=0x7e).contains(&byte) && byte != b'%' {
            out.push(char::from(byte));
        } else {
            out.push_str(&format!("%{byte:02X}"));
        }
    }
    out
}

impl Tok {
    fn wire(&self) -> String {
        match self {
            Self::F(x) => format!("f:{:016x}", x.to_bits()),
            Self::I(x) => format!("i:{x}"),
            Self::E(name, v) => format!("e:{name}.{v}"),
            Self::S(x) => format!("s:{}", pct(x)),
            Self::None => "n:".to_string(),
            Self::Err(v) => format!("E:{v}"),
            Self::V(n) => format!("v:{n}"),
        }
    }

    fn human(&self) -> String {
        match self {
            Self::F(x) => format!("{x:?}"),
            Self::I(x) => format!("{x}"),
            Self::E(name, v) => format!("{name}::{v}"),
            Self::S(x) => format!("{x:?}"),
            Self::None => "None".to_string(),
            Self::Err(v) => format!("Err({v})"),
            Self::V(n) => format!("[{n}]"),
        }
    }
}

fn wire_list(toks: &[Tok]) -> String {
    if toks.is_empty() {
        return "-".to_string();
    }
    toks.iter().map(Tok::wire).collect::<Vec<_>>().join(" ")
}

fn human_list(toks: &[Tok]) -> String {
    toks.iter().map(Tok::human).collect::<Vec<_>>().join(", ")
}

/// The rows of one `tests/golden/<name>.tsv`.
pub struct Module {
    pub name: &'static str,
    rows: Vec<String>,
}

impl Module {
    pub fn new(name: &'static str) -> Self {
        Self {
            name,
            rows: Vec::new(),
        }
    }

    pub fn row(&mut self, func: &str, class: Class, ins: &[Tok], outs: &[Tok]) {
        let note = format!("{func}({}) -> {}", human_list(ins), human_list(outs));
        self.push(func, class, ins, outs, &note);
    }

    /// As [`Module::row`], with `extra` appended to the note in brackets: how the expected value
    /// was obtained (an ADR 0002 mapping, say). For review and `validate.py` only; the harness
    /// ignores the note.
    #[allow(dead_code)] // first used by a row whose expected value is mapped (ADR 0002)
    pub fn row_noted(&mut self, func: &str, class: Class, ins: &[Tok], outs: &[Tok], extra: &str) {
        let note = format!(
            "{func}({}) -> {} [{extra}]",
            human_list(ins),
            human_list(outs)
        );
        self.push(func, class, ins, outs, &note);
    }

    fn push(&mut self, func: &str, class: Class, ins: &[Tok], outs: &[Tok], note: &str) {
        self.rows.push(format!(
            "{func}\t{}\t{}\t{}\t{note}",
            class.letter(),
            wire_list(ins),
            wire_list(outs)
        ));
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// The file contents: a fixed header, then the rows in emission order. Deterministic —
    /// no timestamps or paths — so `--check` can compare byte for byte. `provenance` names
    /// the oracle and the glibc whose libm computed every transcendental (port plan §7).
    pub fn render(&self, provenance: &str) -> String {
        let mut out = String::new();
        out.push_str(&format!("# bijli golden corpus, format 1: {}\n", self.name));
        out.push_str(&format!(
            "# oracle: {provenance}. Generated by scripts/oracle; do not edit.\n"
        ));
        out.push_str(
            "# regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml\n",
        );
        out.push_str("# columns: fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note (tokens: scripts/oracle/README.md)\n");
        out.push_str(&format!("# rows: {}\n", self.rows.len()));
        for row in &self.rows {
            out.push_str(row);
            out.push('\n');
        }
        out
    }
}
