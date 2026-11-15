//! `src/touchstone_data.cyr` — the two Unicode sets the Touchstone reader's text follows, from the oracle's own core
//! (M4 part 1c-1):
//!
//! - `char::is_whitespace` (Unicode White_Space): `lines().trim()` and `split_whitespace()` split on it;
//! - the characters `{:?}` of a `str` writes as `\u{…}`: what isn't printable or is grapheme-extended. The error
//!   messages quote a token with `{:?}` (`snippet`). The probe for M4 part 1c-1 (`port-parity.md` §78) found the
//!   escape context-free: a character is escaped or not wherever it sits in the string.
//!
//! Both are computed over every scalar value by the compiler that builds the oracle, so they are that toolchain's
//! tables: rustc 1.99.1 and 1.100.0 differ by six escape ranges (a Unicode update). The header names the compiler
//! (`build.rs`), and `--check` regenerates and compares, so a toolchain that moves the tables fails it.

use crate::cyr::header;
use std::fmt::Write as _;

/// The scalar values Rust's `{:?}` of a one-character `str` writes as `\u{…}`, as inclusive ranges (the surrogate
/// gap bridged).
fn escape_ranges() -> Vec<(u32, u32)> {
    let mut ranges: Vec<(u32, u32)> = Vec::new();
    for u in 0..=0x0010_FFFF_u32 {
        let Some(c) = char::from_u32(u) else { continue };
        if !format!("{:?}", c.to_string()).starts_with("\"\\u{") {
            continue;
        }
        match ranges.last_mut() {
            Some(r) if r.1 + 1 == u || (r.1 == 0xD7FF && u == 0xE000) => r.1 = u,
            _ => ranges.push((u, u)),
        }
    }
    ranges
}

/// Every `char::is_whitespace` scalar value.
fn whitespace() -> Vec<u32> {
    (0..=0x0010_FFFF_u32)
        .filter(|&u| char::from_u32(u).is_some_and(char::is_whitespace))
        .collect()
}

pub fn emit(prov: &str) -> String {
    let esc = escape_ranges();
    let ws = whitespace();
    let mut out = header(
        "touchstone_data.cyr: the Unicode sets of src/touchstone.cyr's text (M4 part 1c-1)",
        prov,
        &[
            "_bj_ts_esc: the characters Rust's `{:?}` of a str writes as `\\u{…}` (not printable, or\n\
             grapheme-extended), as sorted inclusive ranges in one string literal: a newline, then each range as 12\n\
             uppercase hex digits (start, end; six each), nine to a line, so range r starts at byte\n\
             1 + 109·(r / 9) + 12·(r % 9). A string literal is not static data, where an i64 array of the 916 ranges\n\
             would cost every test binary 7.3 KB of it (the compiler's 128 KiB advisory). The ASCII escapes `\\0`,\n\
             `\\t`, `\\n`, `\\r`, `\\\\` and `\\\"` are src/touchstone.cyr's own.",
            "_bj_ts_ws: char::is_whitespace, Unicode's White_Space, ascending.",
            &format!("Computed by {}.", env!("ORACLE_RUSTC")),
        ],
    );
    let _ = writeln!(out, "const _BJ_TS_ESC_N = {};", esc.len());
    out.push_str("var _bj_ts_esc = \"\n");
    for (k, chunk) in esc.chunks(9).enumerate() {
        if k > 0 {
            out.push('\n');
        }
        for &(a, b) in chunk {
            let _ = write!(out, "{a:06X}{b:06X}");
        }
    }
    out.push_str("\";\n\n");
    let _ = writeln!(out, "const _BJ_TS_WS_N = {};", ws.len());
    let _ = writeln!(out, "var _bj_ts_ws: i64[{}] = {{", ws.len());
    for chunk in ws.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|u| format!("0x{u:X}")).collect();
        let _ = writeln!(out, "    {},", row.join(", "));
    }
    out.push_str("};\n");
    out
}
