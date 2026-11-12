"""The Rust-test ledger: every test of the frozen Rust oracle is ported, pending, or recorded with a reason.

Port plan §6 step 7 (brought forward from M8; the mela pattern, as impetus's scripts/rust_test_ledger.py). Usage,
from the repo root:
    python3 -I -B scripts/rust_test_ledger.py .           # write docs/development/rust-test-ledger.md
    python3 -I -B scripts/rust_test_ledger.py . --check   # exit 1 if the ledger is stale or a test is unaccounted

It lists each `#[test]` / `#[tokio::test]` in rust-old/src/**.rs (inside every `#[cfg(test)]` module, wherever it
sits: scripts/rustsrc.py), each test in rust-old/tests/*.rs, and each doc test (a runnable fence in a doc comment).
A test is **ported** when one of the Cyrius test files FILES maps its Rust file to has a `test_group("<name>")` of
the same name. A test that isn't is **pending** while its milestone is open, and must be ported or have a REASONS
entry once its milestone is in DONE. A REASONS entry for a test that is ported, or that no longer exists, is an
error, as is a Rust file with tests and no FILES entry.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rustsrc  # noqa: E402

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
check = "--check" in sys.argv[2:]
out_path = root / "docs/development/rust-test-ledger.md"

# Rust file -> (milestone, the Cyrius test files its tests are ported into). Port plan §4.1.
FILES = {
    "src/error.rs": ("M1", ["tests/error.tcyr"]),
    "src/logging.rs": ("M1", ["tests/log.tcyr"]),
    "src/field.rs": ("M1", ["tests/field.tcyr"]),
    "src/maxwell.rs": ("M1", ["tests/maxwell.tcyr"]),
    "src/charge.rs": ("M1", ["tests/charge.tcyr"]),
    "src/material.rs": ("M2", ["tests/material.tcyr"]),
    "src/wave.rs": ("M2", ["tests/wave.tcyr"]),
    "src/circuit.rs": ("M2", ["tests/circuit.tcyr"]),
    "src/relativity.rs": ("M2", ["tests/relativity.tcyr"]),
    "src/polarization.rs": ("M3", ["tests/polarization.tcyr"]),
    "src/scattering.rs": ("M3", ["tests/scattering.tcyr"]),
    "src/beam.rs": ("M3", ["tests/beam.tcyr"]),
    "src/rf.rs": ("M4", ["tests/rf.tcyr"]),
    # fdtd.rs's tests split by section: M5 (1D, sources, media, NF2FF, utilities) and M6 (CPML, 2D, 3D).
    "src/fdtd.rs": ("M6", ["tests/fdtd1d.tcyr", "tests/fdtd_media.tcyr", "tests/nf2ff.tcyr", "tests/fdtd_util.tcyr",
                           "tests/cpml.tcyr", "tests/fdtd2d.tcyr", "tests/fdtd3d.tcyr"]),
    "src/compute.rs": ("M6", ["tests/compute.tcyr"]),
    "src/integration/soorat.rs": ("M7", ["tests/soorat.tcyr"]),
    "src/bridge.rs": ("M7", ["tests/bridge.tcyr"]),
    "src/ai.rs": ("M7", ["tests/ai.tcyr"]),
    "tests/integration.rs": ("M8", ["tests/integration.tcyr"]),
}

# Milestones that are complete: a test of theirs that isn't ported needs a REASONS entry.
DONE = set()

# (Rust file, test name) -> why it has no Cyrius port. Added at the milestone that decides it (the ADR 0002 #1
# deprecations, #2 the per-step layout check, ...).
REASONS = {
    ("src/material.rs", "test_material_reflectance_same"):
        "Material::reflectance_at is deprecated in 1.2.0 and not carried (ADR 0002 row 1); tests/material.tcyr's "
        "cyrius_reflectance_via_try checks the same value through try_reflectance_at",
    ("src/material.rs", "test_material_reflectance_air_glass"):
        "deprecated reflectance_at (ADR 0002 row 1); its 0.04 is checked through try_reflectance_at",
    ("src/material.rs", "test_reflectance_at_deprecated_never_wrong"):
        "deprecated reflectance_at (ADR 0002 row 1): the port has no NaN-returning form; try_reflectance_at refuses "
        "the conductor (cyrius_reflectance_via_try)",
}


def cyrius_groups(files):
    groups = set()
    for rel in files:
        path = root / rel
        if path.exists():
            groups.update(re.findall(r'test_group\("([^"]+)"\)', path.read_text()))
    return groups


def collect():
    """Rust file -> [(test name, kind)], kind 'unit', 'integration' or 'doc'."""
    found = {}
    for module in rustsrc.rust_modules(root):
        path = root / "rust-old/src" / f"{module}.rs"
        tests = [(name, "unit") for name, _ in rustsrc.tests_in(path)]
        tests += [(f"doc test at line {no}", "doc") for no in rustsrc.doc_tests(path)]
        if tests:
            found[f"src/{module}.rs"] = tests
    for path in sorted((root / "rust-old/tests").glob("*.rs")):
        tests = [(name, "integration") for name, _ in rustsrc.tests_in(path, whole_file=True)]
        if tests:
            found[f"tests/{path.name}"] = tests
    return found


def main() -> int:
    found = collect()
    problems = []
    rows = []
    totals = {"tests": 0, "ported": 0, "na": 0, "pending": 0}
    for rel in sorted(found):
        if rel not in FILES:
            problems.append(f"{rel}: {len(found[rel])} Rust tests and no FILES entry")
            continue
        milestone, files = FILES[rel]
        groups = cyrius_groups(files)
        ported = na = pending = 0
        for name, _kind in found[rel]:
            key = (rel, name)
            if name in groups:
                ported += 1
                if key in REASONS:
                    problems.append(f"{rel}::{name}: ported, but REASONS still lists it")
            elif key in REASONS:
                na += 1
            elif milestone in DONE:
                problems.append(f"{rel}::{name}: {milestone} is done, but the test is neither ported nor in REASONS")
            else:
                pending += 1
        n = len(found[rel])
        totals["tests"] += n
        totals["ported"] += ported
        totals["na"] += na
        totals["pending"] += pending
        cy = ", ".join(f"`{f}`" for f in files)
        rows.append(f"| `{rel}` | {milestone} | {cy} | {n} | {ported} | {na} | {pending} |")
    for (rel, name) in REASONS:
        if rel not in found or name not in {t for t, _ in found[rel]}:
            problems.append(f"REASONS names {rel}::{name}, which is no Rust test")

    lines = [
        "# Rust-test ledger",
        "",
        "> **Generated** by `python3 -I -B scripts/rust_test_ledger.py .`; `--check` gates it (port plan §6 step 7, the",
        "> mela pattern). Every test of the frozen Rust oracle (`rust-old/`, bijli 1.2.0) is ported to a Cyrius test",
        "> group of the same name, recorded here with its reason, or pending until its milestone. The golden corpus",
        "> (`tests/golden/`) checks every ported function on its grid besides; this ledger is about the oracle's own",
        "> assertions.",
        "",
        f"**{totals['tests']} Rust tests:** {totals['ported']} ported, {totals['na']} not applicable, "
        f"{totals['pending']} pending.",
        "",
        "| Rust file | Milestone | Cyrius tests | Rust tests | Ported | N/A | Pending |",
        "|---|---|---|---:|---:|---:|---:|",
        *rows,
    ]
    if REASONS:
        lines += ["", "## Not applicable", "", "| Rust test | Reason |", "|---|---|"]
        lines += [f"| `{rel}::{name}` | {why} |" for (rel, name), why in sorted(REASONS.items())]
    text = "\n".join(lines) + "\n"

    if check:
        old = out_path.read_text() if out_path.exists() else ""
        if old != text:
            problems.append(f"{out_path.relative_to(root)} is stale: run scripts/rust_test_ledger.py . to rewrite it")
    else:
        out_path.write_text(text)
        print(f"wrote {out_path.relative_to(root)}")
    for p in problems:
        print("  " + p)
    print(f"ledger: {totals['tests']} Rust tests, {totals['ported']} ported, {totals['na']} n/a, "
          f"{totals['pending']} pending; {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
