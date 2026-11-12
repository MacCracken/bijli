"""Reading the frozen Rust oracle's source: shared by surface_parity.py and rust_test_ledger.py.

bijli's Rust files put `#[cfg(test)] mod … { … }` blocks between their sections (`fdtd.rs` has three), so a
"library part" is never "everything before the first test module". Every helper here skips or collects test
modules by brace depth, wherever they sit. Braces inside string or char literals are rare in bijli's source and
balanced where they occur; `split` checks that every block it skips closes.
"""
import re
from pathlib import Path


def rust_modules(root: Path) -> list:
    """Every library source file under rust-old/src, as `fdtd` or `integration/soorat`, sorted."""
    base = root / "rust-old/src"
    return sorted(str(p.relative_to(base).with_suffix("")) for p in base.rglob("*.rs"))


def _is_test_mod_start(lines, i) -> bool:
    if lines[i].strip() != "#[cfg(test)]" or i + 1 >= len(lines):
        return False
    nxt = lines[i + 1].strip()
    return nxt.startswith("mod ") and "{" in nxt


def split(text: str):
    """(library_lines, test_blocks): the file's lines outside its test modules as (line_no, text) pairs, and each
    test module as a list of (line_no, text). Line numbers are 1-based."""
    lines = text.split("\n")
    lib, blocks = [], []
    i = 0
    while i < len(lines):
        if _is_test_mod_start(lines, i):
            block = [(i + 1, lines[i])]
            depth = 0
            j = i + 1
            while j < len(lines):
                block.append((j + 1, lines[j]))
                depth += lines[j].count("{") - lines[j].count("}")
                if depth <= 0:
                    break
                j += 1
            if depth > 0:
                raise ValueError(f"unclosed #[cfg(test)] module at line {i + 1}")
            blocks.append(block)
            i = j + 1
            continue
        lib.append((i + 1, lines[i]))
        i += 1
    return lib, blocks


def library_text(path: Path) -> str:
    """The file with its test modules removed (their lines dropped)."""
    lib, _ = split(path.read_text())
    return "\n".join(t for _, t in lib)


def tests_in(path: Path, whole_file: bool = False) -> list:
    """The `#[test]` / `#[tokio::test]` fn names in a file's test modules (or anywhere in it, for an integration-test file), in source
    order, with their line numbers."""
    text = path.read_text()
    if whole_file:
        scan = list(enumerate(text.split("\n"), 1))
    else:
        _, blocks = split(text)
        scan = [pair for b in blocks for pair in b]
    out = []
    pending = False
    for no, line in scan:
        s = line.strip()
        if s == "#[test]" or re.fullmatch(r"#\[tokio::test(\(.*\))?\]", s):
            pending = True
            continue
        if pending:
            if s.startswith("#["):
                continue
            m = re.match(r"(?:pub )?(?:async )?fn (\w+)", s)
            if m:
                out.append((m.group(1), no))
            pending = False
    return out


def doc_tests(path: Path) -> list:
    """Runnable doc-test fences (```, ```rust, or ```no_run; not text / ignore / compile_fail) in `///` or `//!`
    comments, as the line numbers rustdoc reports (the opening fence)."""
    out = []
    inside = False
    for no, line in enumerate(path.read_text().split("\n"), 1):
        m = re.match(r"\s*(?:///|//!)\s?```\s*([\w,]*)\s*$", line)
        if not m:
            continue
        if inside:
            inside = False
            continue
        inside = True
        attrs = set(filter(None, m.group(1).split(",")))
        if not attrs or attrs <= {"rust", "no_run", "should_panic"}:
            out.append(no)
    return out
