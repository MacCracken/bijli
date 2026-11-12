//! `src/constants.cyr` — bijli's public constants, as `BJ_` consts.
//!
//! bijli 1.2.0 has twelve public scalar consts, all `f64`: the CODATA 2018 field constants
//! (`field`), the particle constants (`charge`), the fourth-order stencil weights (`fdtd`) and the
//! Mie limits (`scattering`). `all` lists every one, and `complete` proves the list against
//! `rust-old/src`. `tests/golden/constants.tsv` pins every value (`m_constants.rs`).
//!
//! Private constants (`ETA_0`, `CPML_MIN_INTERIOR`, `TRACE_INITIAL_CAPACITY`, …) are not public
//! API. They are extracted from the Rust source and proven against the oracle's output at the
//! milestone whose module uses them, as impetus did (port plan §6 step 6).
//!
//! Values come from the oracle's API; the comments are the Rust doc comments.

use crate::cyr::{doc_comment, f64_literal, header};
use std::fmt::Write as _;
use std::path::Path;

/// One public const of the oracle.
pub struct Constant {
    /// `field::EPSILON_0`, the golden row's input.
    pub path: &'static str,
    /// The Cyrius name.
    pub cyr: &'static str,
    /// The Rust file under `rust-old/src`, and the item line its doc comment sits on.
    file: &'static str,
    item: &'static str,
    pub value: f64,
}

/// Every public const of bijli 1.2.0, by file and then in source order.
pub fn all() -> Vec<Constant> {
    vec![
        Constant {
            path: "charge::ELEMENTARY_CHARGE",
            cyr: "BJ_ELEMENTARY_CHARGE",
            file: "charge.rs",
            item: "pub const ELEMENTARY_CHARGE",
            value: bijli::charge::ELEMENTARY_CHARGE,
        },
        Constant {
            path: "charge::ELECTRON_MASS",
            cyr: "BJ_ELECTRON_MASS",
            file: "charge.rs",
            item: "pub const ELECTRON_MASS",
            value: bijli::charge::ELECTRON_MASS,
        },
        Constant {
            path: "charge::PROTON_MASS",
            cyr: "BJ_PROTON_MASS",
            file: "charge.rs",
            item: "pub const PROTON_MASS",
            value: bijli::charge::PROTON_MASS,
        },
        Constant {
            path: "fdtd::FOURTH_ORDER_C1",
            cyr: "BJ_FOURTH_ORDER_C1",
            file: "fdtd.rs",
            item: "pub const FOURTH_ORDER_C1",
            value: bijli::fdtd::FOURTH_ORDER_C1,
        },
        Constant {
            path: "fdtd::FOURTH_ORDER_C2",
            cyr: "BJ_FOURTH_ORDER_C2",
            file: "fdtd.rs",
            item: "pub const FOURTH_ORDER_C2",
            value: bijli::fdtd::FOURTH_ORDER_C2,
        },
        Constant {
            path: "field::EPSILON_0",
            cyr: "BJ_EPSILON_0",
            file: "field.rs",
            item: "pub const EPSILON_0",
            value: bijli::field::EPSILON_0,
        },
        Constant {
            path: "field::MU_0",
            cyr: "BJ_MU_0",
            file: "field.rs",
            item: "pub const MU_0",
            value: bijli::field::MU_0,
        },
        Constant {
            path: "field::SPEED_OF_LIGHT",
            cyr: "BJ_SPEED_OF_LIGHT",
            file: "field.rs",
            item: "pub const SPEED_OF_LIGHT",
            value: bijli::field::SPEED_OF_LIGHT,
        },
        Constant {
            path: "field::COULOMB_K",
            cyr: "BJ_COULOMB_K",
            file: "field.rs",
            item: "pub const COULOMB_K",
            value: bijli::field::COULOMB_K,
        },
        Constant {
            path: "field::SPEED_OF_LIGHT_SQ",
            cyr: "BJ_SPEED_OF_LIGHT_SQ",
            file: "field.rs",
            item: "pub const SPEED_OF_LIGHT_SQ",
            value: bijli::field::SPEED_OF_LIGHT_SQ,
        },
        Constant {
            path: "scattering::MIE_MAX_SIZE_PARAMETER",
            cyr: "BJ_MIE_MAX_SIZE_PARAMETER",
            file: "scattering.rs",
            item: "pub const MIE_MAX_SIZE_PARAMETER",
            value: bijli::scattering::MIE_MAX_SIZE_PARAMETER,
        },
        Constant {
            path: "scattering::MIE_MAX_MX",
            cyr: "BJ_MIE_MAX_MX",
            file: "scattering.rs",
            item: "pub const MIE_MAX_MX",
            value: bijli::scattering::MIE_MAX_MX,
        },
    ]
}

/// The non-test lines of a Rust source file: every `#[cfg(test)] mod … { … }` block is skipped
/// by brace depth, wherever it sits (`fdtd.rs` has test modules between its sections).
fn library_lines(src: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut lines = src.lines().peekable();
    while let Some(line) = lines.next() {
        if line.trim() == "#[cfg(test)]"
            && lines
                .peek()
                .is_some_and(|l| l.trim_start().starts_with("mod ") && l.contains('{'))
        {
            let mut depth = 0_i64;
            for l in lines.by_ref() {
                depth += i64::try_from(l.matches('{').count()).unwrap_or(0);
                depth -= i64::try_from(l.matches('}').count()).unwrap_or(0);
                if depth <= 0 {
                    break;
                }
            }
            continue;
        }
        out.push(line);
    }
    out
}

/// Every `pub const NAME` item in `rust-old/src` outside test code must be in `all()`, and
/// every listed one must exist: a new public constant can't be missed.
fn complete(root: &Path) -> Result<(), String> {
    let listed: Vec<String> = all()
        .iter()
        .map(|c| c.item.trim_start_matches("pub const ").to_string())
        .collect();
    let mut found = Vec::new();
    let mut stack = vec![root.join("rust-old/src")];
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().is_none_or(|x| x != "rs") {
                continue;
            }
            let src =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            for line in library_lines(&src) {
                if let Some(rest) = line.trim_start().strip_prefix("pub const ")
                    && let Some((ident, _)) = rest.split_once(':')
                {
                    found.push(ident.trim().to_string());
                }
            }
        }
    }
    found.sort();
    let mut want = listed;
    want.sort();
    if found == want {
        Ok(())
    } else {
        Err(format!(
            "t_constants::all() lists {want:?}, but rust-old/src has public consts {found:?}"
        ))
    }
}

/// `src/constants.cyr`.
pub fn emit(root: &Path, provenance: &str) -> Result<String, String> {
    complete(root)?;
    let mut out = header(
        "public constants",
        provenance,
        &[
            "Every public constant of the Rust crate; t_constants.rs proves the list complete\n\
             against rust-old/src, and tests/golden/constants.tsv pins each value.",
            "Each is a positional decimal (the lexer has no exponent form) that the lexer rounds\n\
             back to the oracle's exact bits.",
            "Private constants are extracted and proven at the milestone whose module uses them.",
        ],
    );
    for c in &all() {
        if !c.value.is_finite() {
            return Err(format!("{} is not finite", c.path));
        }
        let file = root.join("rust-old/src").join(c.file);
        let src = std::fs::read_to_string(&file).map_err(|e| format!("{}: {e}", file.display()))?;
        out.push_str(&doc_comment(&src, c.item)?);
        let _ = writeln!(out, "# Rust: bijli::{}", c.path);
        let _ = writeln!(out, "const {} = {};", c.cyr, f64_literal(c.value));
        out.push('\n');
    }
    // One trailing newline, not a blank line, at the end of the file.
    while out.ends_with("\n\n") {
        out.pop();
    }
    Ok(out)
}
