"""The benchmark record: every Rust benchmark beside its Cyrius twin.

Port plan §6 step 7 (brought forward from M8; as impetus's scripts/bench_record.py). Usage, from the repo root:
    python3 -I -B scripts/bench_record.py .                 # write docs/benchmarks-rust-v-cyrius.md
    python3 -I -B scripts/bench_record.py . CYRIUS.log      # ... with the Cyrius column from a `cyrius bench` log
    python3 -I -B scripts/bench_record.py . --check         # exit 1 if a twin is missing or the record is stale

The Rust benches are read from rust-old/benches/benchmarks.rs (`benchmark_group("g")` then `bench_function("b")`,
criterion's id `g/b`). The Rust column is the run in rust-old/bench-history.csv stamped BASELINE (port plan §5 P4:
taken under load, so a ratio within that row's own noise is parity, not a finding). A twin is a
`bench_new("<name>")` in tests/*.bcyr; by default it has the Rust bench's own name. A Rust bench without a twin is
pending while its milestone is open, and an error once the milestone is in DONE unless REASONS explains it. A
Cyrius bench that is no twin must be in EXTRA.
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
args = sys.argv[2:]
check = "--check" in args
logs = [a for a in args if a != "--check"]
out_path = root / "docs/benchmarks-rust-v-cyrius.md"

BASELINE = "2026-11-12T03:15:21Z"
BASELINE_NOTE = (
    "`c95fd30` (the 1.2.0 tag's source), criterion medians. Taken under load (kimiya's gates ran beside it: "
    "load average about 2.6, median CI width 0.8 %, p90 4.2 %, max 13 %), so a ratio within the Rust row's own "
    "noise is read as parity (port plan §5 P4)."
)

# criterion group -> milestone (port plan §4.1); fdtd splits by bench name.
GROUP_MILESTONE = {
    "field": "M1", "maxwell": "M1", "charge": "M1",
    "wave": "M2", "circuit": "M2", "material": "M2", "relativity": "M2",
    "polarization": "M3", "scattering": "M3", "beam": "M3",
    "rf": "M4",
}
DONE = set()
# Rust bench -> its Cyrius twin's name, where it isn't the same.
TWINS = {}
# Rust bench -> why it has no twin.
REASONS = {}
# Cyrius benches that twin no Rust bench (the scaffold's harness check among them).
EXTRA = {
    "noop": "the scaffold's harness check (tests/bijli.bcyr)",
    "log_gate_off": "a disabled log site (src/log.cyr): Rust's tracing with no subscriber has no bench",
    "error_record": "making an error into the last-error record (src/error.cyr, ADR 0002 row 4)",
    "error_record_f64": "the same with an f64 field in Rust's Display form (src/fmtf.cyr)",
    "fmt_f64_display": "one f64 in Rust's Display form (src/fmtf.cyr)",
    "complex_mul_div": "a BjComplex multiply and divide, by value (src/complex.cyr)",
    "complex_exp": "BjComplex::exp: e^re and cos, sin of im, each correctly rounded (src/complex.cyr)",
    "sinh_cosh": "one sinh and one cosh over [0, 70), correctly rounded (src/crmath.cyr)",
    "complex_sin": "BjComplex::sin: sin, cos, sinh, cosh, each correctly rounded (src/complex.cyr)",
    "complex_exp_edge": "BjComplex::exp at e^1000, with an unbounded exponent (src/crmath_ux.cyr, ADR 0002 row 16)",
    "complex_powf_edge": "BjComplex::powf with |z|^2.5 past DBL_MAX, on the fixed point (src/crmath_ux.cyr, row 16)",
    "json/field_vector_to_json": "a FieldVector's JSON (src/json.cyr); rust-old benches no serde",
    "json/field_vector_from_json": "a FieldVector read from serde_json's document (src/json.cyr)",
    "json/point_charge_from_json": "a PointCharge read from serde_json's document (src/json.cyr)",
    "json/mueller_matrix_to_json": "a Mueller matrix's JSON, the largest document (src/json.cyr, M3 part 1a-2)",
    "json/mueller_matrix_from_json": "a Mueller matrix read back from its document (src/json.cyr)",
}


def milestone(bench: str) -> str:
    group, _, name = bench.partition("/")
    if group == "fdtd":
        return "M6" if name.startswith(("2d_", "3d_")) else "M5"
    return GROUP_MILESTONE[group]


def rust_benches() -> list:
    out, group = [], None
    for line in (root / "rust-old/benches/benchmarks.rs").read_text().split("\n"):
        m = re.search(r'benchmark_group\("([^"]+)"\)', line)
        if m:
            group = m.group(1)
        m = re.search(r'bench_function\("([^"]+)"', line)
        if m:
            out.append(f"{group}/{m.group(1)}")
    return out


def rust_ns() -> dict:
    ns = {}
    for line in (root / "rust-old/bench-history.csv").read_text().split("\n"):
        f = line.split(",")
        if len(f) == 5 and f[0] == BASELINE:
            ns[f[3]] = float(f[4])
    return ns


def cyrius_benches() -> set:
    names = set()
    for path in (root / "tests").glob("*.bcyr"):
        names.update(re.findall(r'bench_new\("([^"]+)"\)', path.read_text()))
    return names


def cyrius_ns(paths) -> dict:
    """`  name: 3ns avg (...)` lines from `cyrius bench` logs, in ns."""
    unit = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "µs": 1e3, "ms": 1e6, "s": 1e9}
    ns = {}
    for p in paths:
        for line in Path(p).read_text().split("\n"):
            m = re.match(r"\s*(\S+): ([\d.]+)(ps|ns|us|µs|ms|s) avg", line)
            if m:
                ns[m.group(1)] = float(m.group(2)) * unit[m.group(3)]
    return ns


def fmt_ns(x: float) -> str:
    if x >= 1e6:
        return f"{x / 1e6:.2f} ms"
    if x >= 1e3:
        return f"{x / 1e3:.2f} µs"
    return f"{x:.2f} ns"


def main() -> int:
    benches = rust_benches()
    rns = rust_ns()
    have = cyrius_benches()
    cns = cyrius_ns(logs)
    problems = []
    if len(benches) != len(set(benches)):
        problems.append("rust-old/benches/benchmarks.rs has duplicate bench ids")
    missing_rust = [b for b in benches if b not in rns]
    if missing_rust:
        problems.append(f"{len(missing_rust)} Rust benches have no {BASELINE} row: {missing_rust[:5]}")
    twins_used = set()
    rows, counts = [], {"twinned": 0, "pending": 0, "na": 0}
    for b in benches:
        twin = TWINS.get(b, b)
        ms = milestone(b)
        if twin in have:
            twins_used.add(twin)
            counts["twinned"] += 1
            c = fmt_ns(cns[twin]) if twin in cns else "—"
            ratio = f"{cns[twin] / rns[b]:.1f}×" if twin in cns and b in rns and rns[b] > 0 else "—"
            status = f"`{twin}`"
        elif b in REASONS:
            counts["na"] += 1
            c, ratio, status = "—", "—", f"n/a: {REASONS[b]}"
        elif ms in DONE:
            problems.append(f"{b}: {ms} is done, but the bench has no twin `{twin}` and no REASONS entry")
            c, ratio, status = "—", "—", "MISSING"
        else:
            counts["pending"] += 1
            c, ratio, status = "—", "—", f"pending ({ms})"
        r = fmt_ns(rns[b]) if b in rns else "—"
        rows.append(f"| `{b}` | {r} | {status} | {c} | {ratio} |")
    for name in sorted(have - twins_used):
        if name not in EXTRA:
            problems.append(f"Cyrius bench `{name}` twins no Rust bench and is not in EXTRA")
    for b in REASONS:
        if b not in benches:
            problems.append(f"REASONS names {b}, which is no Rust bench")

    lines = [
        "# Benchmarks: Rust 1.2.0 vs Cyrius",
        "",
        "> **Generated** by `python3 -I -B scripts/bench_record.py .` (pass a `cyrius bench` log for the Cyrius",
        "> column); `--check` gates the twin mapping (port plan §6 step 7). Speed is measured, not gated (port plan",
        "> §0.1): a slow row is a finding to file with its numbers.",
        "",
        f"**Rust column**: rust-old/bench-history.csv, run `{BASELINE}` at {BASELINE_NOTE}",
        "",
        f"**{len(benches)} Rust benches:** {counts['twinned']} twinned, {counts['na']} not applicable, "
        f"{counts['pending']} pending.",
        "",
        "| Rust bench | Rust | Cyrius twin | Cyrius | Cyrius / Rust |",
        "|---|---:|---|---:|---:|",
        *rows,
    ]
    text = "\n".join(lines) + "\n"
    if check:
        old = out_path.read_text() if out_path.exists() else ""
        # The Cyrius column depends on a log; --check compares the mapping only.
        strip = lambda t: re.sub(r"\| [^|]* \| [^|]* \|\n", "|\n", t)  # noqa: E731
        if strip(old) != strip(text):
            problems.append(f"{out_path.relative_to(root)} is stale: run scripts/bench_record.py . to rewrite it")
    else:
        out_path.write_text(text)
        print(f"wrote {out_path.relative_to(root)}")
    for p in problems:
        print("  " + p)
    print(f"bench record: {len(benches)} Rust benches, {counts['twinned']} twinned, {counts['na']} n/a, "
          f"{counts['pending']} pending; {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
