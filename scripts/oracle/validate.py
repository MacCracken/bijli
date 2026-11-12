#!/usr/bin/env python3
"""Independent check of bijli's golden corpus (port plan §6 step 5).

The generator (`scripts/oracle`, Rust) records what the frozen Rust bijli 1.2.0 returns. This
script recomputes every row from a separate transcription of the same semantics, so that a
defect in the generator (a wrong grid, a mislabelled function, a value written for the wrong
input) is caught before the Cyrius port is compared against it:

- **constants**: each public const's literal is read from `rust-old/src` and evaluated in
  Python with the Rust expression's operation order (`COULOMB_K = 1 / (4·π·ε₀)`);
- **math**: glibc's libm through `ctypes` for the transcendentals (the same glibc the oracle
  links; the header's platform is checked first), except `cbrt`: Rust std links compiler-builtins'
  own `cbrt` (its `libm` crate's port of CORE-MATH, correctly rounded), which MPFR reproduces
  (`mpfr_ref.py`). Python's IEEE doubles for the exact and
  arithmetic functions, with Rust's semantics transcribed (NaN-ignoring `min`/`max`, `clamp`,
  the folded `to_radians` / `to_degrees` constants, compiler-rt's `__powidf2` loop).

It also writes and checks two generated outputs (`--write` writes them; a plain run checks them):
- **the corrected oracle**, `tests/golden/math.cr.tsv` (kimiya ADR 0006, bijli ADR 0001): for every class-C row
  whose glibc result is not correctly rounded, the correctly rounded outputs (MPFR), which the harness compares
  against instead of glibc's;
- **src/crmath_data.cyr**, crtables.py's rendering, with its π and ln 2 and every table's hi checked against MPFR.

Run with `python3 -I -B scripts/oracle/validate.py . [--write]` from the repo root. Exit status 0 when every row
agrees and the generated files are current, 1 otherwise. Copied in spirit from impetus's and kimiya's validate.py.
"""

import ctypes
import math
import os
import platform
import re
import struct
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpfr_ref  # noqa: E402  (system MPFR through ctypes)
import crtables  # noqa: E402  (src/crmath_data.cyr's generator)
import crtables_k  # noqa: E402  (src/crmath_k_data.cyr's generator, kimiya's ln and atan tables)

WRITE = "--write" in sys.argv[2:]


def bits_to_f(h: str) -> float:
    return struct.unpack("<d", struct.pack("<Q", int(h, 16)))[0]


def f_to_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def same(a: float, b: float) -> bool:
    """Bit-identical, with any NaN matching any NaN (format 1's class-A rule)."""
    if math.isnan(a) and math.isnan(b):
        return True
    return f_to_bits(a) == f_to_bits(b)


# ── glibc ─────────────────────────────────────────────────────────────────────────────────────
LIBM = ctypes.CDLL("libm.so.6")
for _name in ("sin", "cos", "tan", "asin", "acos", "atan", "exp", "log", "log10", "cbrt", "sinh", "cosh"):
    getattr(LIBM, _name).restype = ctypes.c_double
    getattr(LIBM, _name).argtypes = [ctypes.c_double]
for _name in ("atan2", "pow", "hypot"):
    getattr(LIBM, _name).restype = ctypes.c_double
    getattr(LIBM, _name).argtypes = [ctypes.c_double, ctypes.c_double]
LIBC = ctypes.CDLL("libc.so.6")
LIBC.gnu_get_libc_version.restype = ctypes.c_char_p


def glibc_version() -> str:
    return LIBC.gnu_get_libc_version().decode()


# ── Rust semantics, transcribed ───────────────────────────────────────────────────────────────
PI = math.pi  # core::f64::consts::PI, the same double


def rust_max(a: float, b: float) -> float:
    """f64::max: NaN is ignored; between ±0 either may be returned. LLVM's `maxnum` on x86_64
    returns the second operand when the first is NaN, else the larger, preferring the first on a
    tie; the oracle's bits decide ±0 rows, so ±0 ties are accepted either way."""
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return a if a > b or (a == b) else b


def rust_min(a: float, b: float) -> float:
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return a if a < b or (a == b) else b


def rust_clamp(x: float, lo: float, hi: float) -> float:
    """f64::clamp: `if x < lo { lo } else if x > hi { hi } else { x }` (NaN passes through)."""
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


TO_RADIANS = PI / 180.0  # std: `self * (consts::PI / 180.0)`, folded
TO_DEGREES = 180.0 / PI  # std: `self * (180.0f64 / consts::PI)`, folded


def powidf2(a: float, b: int) -> float:
    """compiler-builtins `__powidf2`: square-and-multiply, then a reciprocal for b < 0."""
    recip = b < 0
    p = abs(b)
    mul = 1.0
    while True:
        if p & 1:
            mul *= a
        p >>= 1
        if p == 0:
            break
        a *= a
    return 1.0 / mul if recip else mul


def ieee_div(a: float, b: float) -> float:
    """IEEE division, with Python's ZeroDivisionError replaced by the IEEE results."""
    try:
        return a / b
    except ZeroDivisionError:
        if math.isnan(a) or a == 0.0:
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)


def powidf2_ieee(a: float, b: int) -> float:
    recip = b < 0
    p = abs(b)
    mul = 1.0
    while True:
        if p & 1:
            try:
                mul *= a
            except OverflowError:
                mul = math.copysign(math.inf, mul) * math.copysign(1.0, a)
        p >>= 1
        if p == 0:
            break
        try:
            a *= a
        except OverflowError:
            a = math.inf
    return ieee_div(1.0, mul) if recip else mul


def mul_ieee(x: float, k: float) -> float:
    try:
        return x * k
    except OverflowError:
        return math.copysign(math.inf, x) * math.copysign(1.0, k)


def sqrt_ieee(x: float) -> float:
    if math.isnan(x):
        return x
    if x < 0:
        return math.nan
    if x == 0:
        return x  # keeps -0.0
    if math.isinf(x):
        return x
    return math.sqrt(x)


def ceil_ieee(x: float) -> float:
    if math.isnan(x) or math.isinf(x) or x == 0.0:
        return x
    r = float(math.ceil(x))
    return -0.0 if r == 0.0 and x < 0 else r


MATH = {
    "f64::sqrt": lambda xs: [sqrt_ieee(xs[0])],
    "f64::abs": lambda xs: [abs(xs[0]) if not math.isnan(xs[0]) else bits_to_f("%016x" % (f_to_bits(xs[0]) & 0x7FFFFFFFFFFFFFFF))],
    "f64::ceil": lambda xs: [ceil_ieee(xs[0])],
    "f64::max": lambda xs: [rust_max(xs[0], xs[1])],
    "f64::min": lambda xs: [rust_min(xs[0], xs[1])],
    "f64::clamp": lambda xs: [rust_clamp(xs[0], xs[1], xs[2])],
    "f64::copysign": lambda xs: [bits_to_f("%016x" % ((f_to_bits(xs[0]) & 0x7FFFFFFFFFFFFFFF) | (f_to_bits(xs[1]) & 0x8000000000000000)))],
    "f64::to_radians": lambda xs: [mul_ieee(xs[0], TO_RADIANS)],
    "f64::to_degrees": lambda xs: [mul_ieee(xs[0], TO_DEGREES)],
    "f64::sin_cos": lambda xs: [LIBM.sin(xs[0]), LIBM.cos(xs[0])],
    "f64::sin": lambda xs: [LIBM.sin(xs[0])],
    "f64::cos": lambda xs: [LIBM.cos(xs[0])],
    "f64::tan": lambda xs: [LIBM.tan(xs[0])],
    "f64::sinh": lambda xs: [LIBM.sinh(xs[0])],
    "f64::cosh": lambda xs: [LIBM.cosh(xs[0])],
    "f64::exp": lambda xs: [LIBM.exp(xs[0])],
    "f64::atan2": lambda xs: [LIBM.atan2(xs[0], xs[1])],
    "f64::powf": lambda xs: [LIBM.pow(xs[0], xs[1])],
    "f64::hypot": lambda xs: [LIBM.hypot(xs[0], xs[1])],
    "f64::log10": lambda xs: [LIBM.log10(xs[0])],
    "f64::atan": lambda xs: [LIBM.atan(xs[0])],
    "f64::asin": lambda xs: [LIBM.asin(xs[0])],
    "f64::acos": lambda xs: [LIBM.acos(xs[0])],
    # compiler-builtins' correctly rounded cbrt, not glibc's (P3): MPFR gives its bits.
    "f64::cbrt": lambda xs: [mpfr_ref.cbrt(xs[0])],
}

# min / max between ±0: Rust's result is LLVM's choice, not specified; the row stands as recorded.
SIGNED_ZERO_TIE = {"f64::max", "f64::min"}


# ── the corpus ────────────────────────────────────────────────────────────────────────────────
def parse_tok(t: str):
    kind, _, body = t.partition(":")
    if kind == "f":
        return ("f", bits_to_f(body))
    if kind == "i":
        return ("i", int(body))
    return (kind, body)


def read_tsv(path: Path):
    header, rows = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            header.append(line)
            continue
        fn, cls, ins, outs, _note = line.split("\t")
        toks_in = [] if ins == "-" else [parse_tok(t) for t in ins.split(" ")]
        toks_out = [] if outs == "-" else [parse_tok(t) for t in outs.split(" ")]
        rows.append((fn, cls, toks_in, toks_out))
    return header, rows


def check_platform(header) -> list:
    errs = []
    want = f"glibc {glibc_version()} {platform.machine()}"
    line = next((h for h in header if h.startswith("# oracle: ")), "")
    if f"libm: {want}" not in line:
        errs.append(f"recorded on another platform ({line!r}); running {want}")
    return errs


def validate_math(path: Path) -> list:
    header, rows = read_tsv(path)
    errs = check_platform(header)
    if errs:
        return errs
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        if fn == "f64::is_finite" or fn == "f64::is_nan":
            x = ins[0][1]
            want = int(math.isfinite(x)) if fn == "f64::is_finite" else int(math.isnan(x))
            if outs != [("i", want)]:
                errs.append(f"math row {n}: {fn}({x!r}) recorded {outs}, transcription {want}")
            continue
        if fn == "f64::powi":
            x, k = ins[0][1], ins[1][1]
            got = [powidf2_ieee(x, k)]
        elif fn in MATH:
            got = MATH[fn]([t[1] for t in ins])
        else:
            errs.append(f"math row {n}: no transcription for {fn}")
            continue
        want = [t[1] for t in outs]
        if len(got) != len(want):
            errs.append(f"math row {n}: {fn} arity {len(want)} vs {len(got)}")
            continue
        for g, w in zip(got, want):
            if same(g, w):
                continue
            if fn in SIGNED_ZERO_TIE and g == 0.0 and w == 0.0:
                continue
            errs.append(
                f"math row {n}: {fn}{tuple(t[1] for t in ins)} recorded {w!r} "
                f"({f_to_bits(w):016x}), transcription {g!r} ({f_to_bits(g):016x})"
            )
    return errs


# ── the corrected oracle ──────────────────────────────────────────────────────────────────────
# A class-C math row's reference is the oracle's operation with its primitive correctly rounded. For a primitive
# row that is the MPFR value itself. Where glibc rounded correctly the oracle's row stands and nothing is written.
CR_UNARY = {"f64::sin": "sin", "f64::cos": "cos", "f64::tan": "tan", "f64::asin": "asin", "f64::acos": "acos",
            "f64::atan": "atan", "f64::exp": "exp", "f64::log10": "log10", "f64::cbrt": "cbrt",
            "f64::sinh": "sinh", "f64::cosh": "cosh"}
CR_BINARY = {"f64::atan2": "atan2", "f64::powf": "pow", "f64::hypot": "hypot"}


def corrections(rows) -> list:
    out = []
    for fn, cls, ins, outs in rows:
        if cls != "C":
            continue
        xs = [t[1] for t in ins]
        if fn == "f64::sin_cos":
            calls = [("sin", xs), ("cos", xs)]
        elif fn in CR_UNARY:
            calls = [(CR_UNARY[fn], xs)]
        elif fn in CR_BINARY:
            calls = [(CR_BINARY[fn], xs)]
        else:
            continue
        want, log = [], []
        for (name, args), got in zip(calls, (t[1] for t in outs)):
            w = mpfr_ref.cr(name, args[0]) if len(args) == 1 else mpfr_ref.cr2(name, args[0], args[1])
            want.append(w)
            if not same(w, got):
                arg = ", ".join(f"f:{f_to_bits(a):016x}" for a in args)
                log.append(f"{name}({arg}) glibc f:{f_to_bits(got):016x} correctly rounded f:{f_to_bits(w):016x}")
        if log:
            ins_s = " ".join(f"f:{f_to_bits(a):016x}" for a in xs)
            outs_old = " ".join(f"f:{f_to_bits(t[1]):016x}" for t in outs)
            outs_new = " ".join(f"f:{f_to_bits(w):016x}" for w in want)
            out.append((fn, ins_s, outs_new, f"oracle {outs_old}; " + "; ".join(log)))
    return out


def cr_text(module: str, header: list, rows, fixes_of=None) -> str:
    fixes = (fixes_of or corrections)(rows)
    if not fixes:
        return ""
    oracle = next(h for h in header if h.startswith("# oracle:")).replace(" Generated by scripts/oracle; do not edit.", "")
    head = [
        f"# bijli corrected oracle (ADR 0001), format 1: {module}",
        f"{oracle} Correctly rounded: MPFR {mpfr_ref.VERSION}. Generated by scripts/oracle/validate.py; do not edit.",
        "# regenerate: python3 -I -B scripts/oracle/validate.py . --write",
        "# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: the "
        "oracle's operation with its primitive correctly rounded; note: the oracle's outputs, then each misrounded call)",
        f"# rows: {len(fixes)}",
    ]
    return "\n".join(head + ["\t".join(f) for f in fixes]) + "\n"


def check_corrections(root: Path, module: str, path_tsv: Path, fixes_of=None) -> list:
    header, rows = read_tsv(path_tsv)
    text = cr_text(module, header, rows, fixes_of)
    path = path_tsv.with_name(f"{module}.cr.tsv")
    rel = f"tests/golden/{module}.cr.tsv"
    if WRITE:
        if text:
            path.write_text(text)
            print(f"wrote {rel} ({text.count(chr(10)) - 5} rows)")
        elif path.exists():
            path.unlink()
            print(f"removed {rel}")
        return []
    if not text:
        return [f"{rel}: no row needs a correction, but the file exists"] if path.exists() else []
    if not path.exists() or path.read_text() != text:
        return [f"{rel}: differs from the corrected oracle (regenerate with --write)"]
    print(f"ok   {module}.cr.tsv: {text.count(chr(10)) - 5} corrected rows")
    return []


def check_tables(root: Path) -> list:
    """src/crmath_data.cyr is crtables.py's rendering, and its π, ln 2 and every table's hi are MPFR's."""
    errs = []
    n = crtables.PREC
    for name, fn in (("pi", crtables.pi_bits), ("log2", crtables.ln2_bits)):
        if fn(n) != mpfr_ref.const_bits(name, n):
            errs.append(f"crtables: {name} to {n} bits differs from MPFR's")
    t = crtables.tables()
    for j, (hi, _lo) in enumerate(t["exp2"]):
        if hi != mpfr_ref.bits(mpfr_ref.cr("exp2", j / 64)):
            errs.append(f"crtables: 2^({j}/64) hi is not MPFR's correctly rounded value")
    for i, (sh, _sl, ch, _cl) in enumerate(t["sincos"]):
        if sh != mpfr_ref.bits(mpfr_ref.cr("sin", i / 64)) or ch != mpfr_ref.bits(mpfr_ref.cr("cos", i / 64)):
            errs.append(f"crtables: sin or cos of {i}/64 hi is not MPFR's correctly rounded value")
    for i, (hi, _lo) in enumerate(t["atan"]):
        if hi != mpfr_ref.bits(mpfr_ref.cr("atan", i / 64)):
            errs.append(f"crtables: atan({i}/64) hi is not MPFR's correctly rounded value")
    for name, want in (("PI_HI", mpfr_ref.atan2(0.0, -1.0)), ("PIO2_HI", mpfr_ref.atan2(1.0, 0.0)),
                       ("PIO4", mpfr_ref.atan2(1.0, 1.0)), ("PI3O4", mpfr_ref.atan2(1.0, -1.0))):
        if t[name] != mpfr_ref.bits(want):
            errs.append(f"crtables: {name} is not MPFR's correctly rounded value")
    path = root / "src/crmath_data.cyr"
    text = crtables.render()
    if WRITE:
        path.write_text(text)
        print(f"wrote src/crmath_data.cyr ({len(text.splitlines())} lines)")
    elif not path.exists() or path.read_text() != text:
        errs.append("src/crmath_data.cyr is not crtables.py's output (regenerate with --write)")
    elif not errs:
        print(f"ok   crmath_data: π and ln 2 to {n} bits, and every table's hi, equal MPFR {mpfr_ref.VERSION}'s")
    kpath = root / "src/crmath_k_data.cyr"
    ktext = crtables_k.render()
    if WRITE:
        kpath.write_text(ktext)
        print(f"wrote src/crmath_k_data.cyr ({len(ktext.splitlines())} lines)")
    elif not kpath.exists() or kpath.read_text() != ktext:
        errs.append("src/crmath_k_data.cyr is not crtables_k.py's output (regenerate with --write)")
    else:
        print(f"ok   crmath_k_data: crtables_k.py's output (MPFR {mpfr_ref.VERSION}, 400-bit values)")
    return errs


def rust_literal(text: str) -> float:
    return float(text.replace("_", ""))


def validate_constants(root: Path, path: Path) -> list:
    _header, rows = read_tsv(path)
    src = {}
    for f in ("field.rs", "charge.rs", "fdtd.rs", "scattering.rs"):
        for m in re.finditer(r"pub const ([A-Z0-9_]+): f64 = (.+);", (root / "rust-old/src" / f).read_text()):
            src[(f[:-3], m.group(1))] = m.group(2).strip()
    env = {}

    def value(mod: str, name: str) -> float:
        expr = src[(mod, name)]
        if re.fullmatch(r"-?[0-9_]+(\.[0-9_]*)?(e-?[0-9]+)?", expr):
            return rust_literal(expr)
        if expr == "1.0 / (4.0 * std::f64::consts::PI * EPSILON_0)":
            return 1.0 / (4.0 * PI * value("field", "EPSILON_0"))
        if expr == "SPEED_OF_LIGHT * SPEED_OF_LIGHT":
            c = value("field", "SPEED_OF_LIGHT")
            return c * c
        if expr == "9.0 / 8.0":
            return 9.0 / 8.0
        if expr == "-1.0 / 24.0":
            return -1.0 / 24.0
        raise ValueError(f"{mod}::{name}: no transcription for {expr!r}")

    errs = []
    for fn, _cls, ins, outs in rows:
        mod, _, name = ins[0][1].partition("::")
        try:
            want = value(mod, name)
        except (KeyError, ValueError) as e:
            errs.append(f"constants: {ins[0][1]}: {e}")
            continue
        env[name] = want
        if fn != "const" or not same(outs[0][1], want):
            errs.append(f"constants: {ins[0][1]} recorded {outs[0][1]!r}, transcription {want!r}")
    if len(rows) != len(src):
        errs.append(f"constants: {len(rows)} rows, rust-old/src has {len(src)} public f64 consts")
    return errs


# ── fmtf: Rust's Display and serde_json's f64 text ───────────────────────────────────────────
from fractions import Fraction  # noqa: E402
from urllib.parse import unquote  # noqa: E402


def tok_text(body: str) -> str:
    """An `s:` token's text: its %XX escapes decoded as UTF-8."""
    return unquote(body)


def shortest(a: float, rust_ties: bool):
    """(digits, k) with a = 0.d1…dn × 10^k: the shortest digits that read back to the positive finite a, the closest
    to it. Python's repr gives them (Gay's shortest mode, an exact tie between two closest candidates to the even
    one, as serde_json's zmij); Rust's core printer takes the upper one on such a tie (src/fmtf.cyr)."""
    mant, _, ex = repr(a).partition("e")
    ip, _, fp = mant.partition(".")
    alld = (ip + fp).lstrip("0")
    k = len(alld) + (int(ex) if ex else 0) - len(fp)
    ds = alld.rstrip("0")
    if rust_ties:
        n = len(ds)
        q = Fraction(a) * Fraction(10) ** (n - k)
        if (2 * q).denominator == 1 and (2 * q).numerator % 2 == 1:
            lo = q.numerator // q.denominator
            hi = lo + 1

            def back(m):
                return float(Fraction(m) * Fraction(10) ** (k - n)) == a
            if back(lo) and back(hi) and str(lo) == ds:
                h = str(hi)
                k += len(h) - n
                ds = h.rstrip("0")
    return ds, k


def positional(ds: str, k: int, int_suffix: bool) -> str:
    n = len(ds)
    if k <= 0:
        return "0." + "0" * (-k) + ds
    if k >= n:
        return ds + "0" * (k - n) + (".0" if int_suffix else "")
    return ds[:k] + "." + ds[k:]


def rust_display(x: float) -> str:
    """`format!("{x}")` for an f64."""
    if math.isnan(x):
        return "NaN"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if math.isinf(x):
        return sign + "inf"
    if x == 0:
        return sign + "0"
    ds, k = shortest(abs(x), True)
    return sign + positional(ds, k, False)


def serde_json_f64(x: float) -> str:
    """`serde_json::to_string(&x)` (1.0.154, zmij) for an f64."""
    if math.isnan(x) or math.isinf(x):
        return "null"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if x == 0:
        return sign + "0.0"
    ds, k = shortest(abs(x), False)
    ex = k - 1
    if -5 <= ex <= 15:
        return sign + positional(ds, k, True)
    return sign + ds[0] + ("." + ds[1:] if len(ds) > 1 else "") + "e" + ("+" if ex >= 0 else "-") + str(abs(ex))


def validate_fmtf(path: Path) -> list:
    """Both printers transcribed: every row's text, byte for byte."""
    _header, rows = read_tsv(path)
    errs = []
    ties = 0
    for fn, _cls, ins, outs in rows:
        x = ins[0][1]
        got = tok_text(outs[0][1])
        if fn == "f64::fmt_display":
            want = rust_display(x)
            if math.isfinite(x) and x != 0 and shortest(abs(x), True) != shortest(abs(x), False):
                ties += 1
        elif fn == "serde_json::f64":
            want = serde_json_f64(x)
        else:
            errs.append(f"fmtf: unknown function {fn}")
            continue
        if got != want:
            errs.append(f"fmtf: {fn}({x!r}): Rust {got!r}, transcription {want!r}")
    print(f"     fmtf: {ties} Display rows are exact ties (Rust rounds up, Python to even), each matched")
    return errs


# ── error: BijliError's Display text ──────────────────────────────────────────────────────────
ERROR_TEXT = {
    "BijliError::DimensionMismatch": lambda i: f"invalid field dimensions: expected {i[0]}D, got {i[1]}D",
    "BijliError::ZeroCharge": lambda i: "charge magnitude is zero — cannot compute field direction",
    "BijliError::Singularity": lambda i: "singularity: field evaluation at source point (distance = 0)",
    "BijliError::InsufficientResolution":
        lambda i: f"grid resolution too coarse: {i[0]} cells for domain size {rust_display(i[1])}",
    "BijliError::InvalidPermittivity": lambda i: f"invalid permittivity: {rust_display(i[0])} (must be > 0)",
    "BijliError::InvalidPermeability": lambda i: f"invalid permeability: {rust_display(i[0])} (must be > 0)",
    "BijliError::DivisionByZero": lambda i: f"division by zero in EM calculation: {i[0]}",
    "BijliError::InvalidParameter": lambda i: f"invalid parameter: {i[0]}",
    "BijliError::InvalidParameter/f64": lambda i: f"invalid parameter: {i[0]}{rust_display(i[1])}",
}


def validate_error(path: Path) -> list:
    """thiserror's messages transcribed from rust-old/src/error.rs, the f64 fields through rust_display; the
    variant token names the row's variant, and the payload tokens repeat the inputs (the reason, for the /f64
    rows, is the prefix and the value's Display)."""
    _header, rows = read_tsv(path)
    errs = []
    for fn, _cls, ins, outs in rows:
        vals = [tok_text(b) if k == "s" else b for k, b in ins]
        if fn not in ERROR_TEXT:
            errs.append(f"error: unknown function {fn}")
            continue
        variant = fn.split("::")[1].split("/")[0]
        if outs[0] != ("E", variant):
            errs.append(f"error: {fn}: variant token {outs[0]}")
        want_payload = vals
        if fn.endswith("/f64"):
            want_payload = [vals[0] + rust_display(vals[1])]
        got_payload = [tok_text(b) if k == "s" else b for k, b in outs[1:-1]]
        if len(got_payload) != len(want_payload) or any(
                not (same(g, w) if isinstance(w, float) else g == w) for g, w in zip(got_payload, want_payload)):
            errs.append(f"error: {fn}{tuple(vals)}: payload {got_payload!r}, inputs {want_payload!r}")
        text = tok_text(outs[-1][1])
        want = ERROR_TEXT[fn](vals)
        if text != want:
            errs.append(f"error: {fn}{tuple(vals)}: Rust {text!r}, transcription {want!r}")
    return errs


# ── complex: Rust hisab 1.4.0's Complex ───────────────────────────────────────────────────────
# Every operation transcribed from hisab-1.4.0/src/num/complex.rs in its order of operations, with IEEE double
# arithmetic (Python's floats; division and sqrt through the helpers below, since Python raises where IEEE gives
# ±∞ or NaN). Each one runs twice: with glibc's primitives it must reproduce the oracle's row exactly; with MPFR's
# it is the corrected oracle, written to complex.cr.tsv for the class-C rows where the two differ.
def fdiv(a: float, b: float) -> float:
    if b != 0.0 or math.isnan(b):
        return a / b if not math.isnan(b) else math.nan
    if a == 0.0 or math.isnan(a):
        return math.nan
    return math.copysign(math.inf, a) * math.copysign(1.0, b)


def fsqrt(x: float) -> float:
    return math.nan if x < 0.0 or math.isnan(x) else math.sqrt(x)


GLIBC_P = {"sin": LIBM.sin, "cos": LIBM.cos, "exp": LIBM.exp, "log": LIBM.log, "atan2": LIBM.atan2,
           "pow": LIBM.pow, "sinh": LIBM.sinh, "cosh": LIBM.cosh}
MPFR_P = {"sin": mpfr_ref.sin, "cos": mpfr_ref.cos, "exp": mpfr_ref.exp, "log": lambda x: mpfr_ref.cr("log", x),
          "atan2": mpfr_ref.atan2, "pow": mpfr_ref.pow, "sinh": lambda x: mpfr_ref.cr("sinh", x),
          "cosh": lambda x: mpfr_ref.cr("cosh", x)}
TINY_300 = bits_to_f("01a56e1fc2f8f359")  # 1e-300


def cx_abs(re, im):
    return fsqrt(re * re + im * im)


def cx_polar(P, r, t):
    return [r * P["cos"](t), r * P["sin"](t)]


def cx_div(ar, ai, br, bi):
    d = br * br + bi * bi
    return [fdiv(ar * br + ai * bi, d), fdiv(ai * br - ar * bi, d)]


def cx_powf(P, re, im, n):
    if re * re + im * im < TINY_300:
        return [0.0, 0.0]
    return cx_polar(P, P["pow"](cx_abs(re, im), n), P["atan2"](im, re) * n)


COMPLEX = {
    "Complex::from_real": lambda P, v: [v[0], 0.0],
    "polarization::complex_real": lambda P, v: [v[0], 0.0],
    "Complex::from_polar": lambda P, v: cx_polar(P, v[0], v[1]),
    "polarization::complex_from_polar": lambda P, v: cx_polar(P, v[0], v[1]),
    "Complex::default": lambda P, v: [0.0, 0.0],
    "polarization::complex_zero": lambda P, v: [0.0, 0.0],
    "Complex::new": lambda P, v: [v[0], v[1]],
    "Complex::abs": lambda P, v: [cx_abs(v[0], v[1])],
    "ComplexExt::norm": lambda P, v: [cx_abs(v[0], v[1])],
    "Complex::norm_sq": lambda P, v: [v[0] * v[0] + v[1] * v[1]],
    "Complex::conj": lambda P, v: [v[0], -v[1]],
    "Complex::neg": lambda P, v: [-v[0], -v[1]],
    "Complex::inv": lambda P, v: [fdiv(v[0], v[0] * v[0] + v[1] * v[1]), fdiv(-v[1], v[0] * v[0] + v[1] * v[1])],
    "Complex::arg": lambda P, v: [P["atan2"](v[1], v[0])],
    "Complex::exp": lambda P, v: [P["exp"](v[0]) * P["cos"](v[1]), P["exp"](v[0]) * P["sin"](v[1])],
    "Complex::ln": lambda P, v: [P["log"](cx_abs(v[0], v[1])), P["atan2"](v[1], v[0])],
    "Complex::sqrt": lambda P, v: cx_polar(P, fsqrt(cx_abs(v[0], v[1])), P["atan2"](v[1], v[0]) / 2.0),
    "Complex::powf": lambda P, v: cx_powf(P, v[0], v[1], v[2]),
    "Complex::sin": lambda P, v: [P["sin"](v[0]) * P["cosh"](v[1]), P["cos"](v[0]) * P["sinh"](v[1])],
    "Complex::cos": lambda P, v: [P["cos"](v[0]) * P["cosh"](v[1]), -P["sin"](v[0]) * P["sinh"](v[1])],
    "Complex::add": lambda P, v: [v[0] + v[2], v[1] + v[3]],
    "Complex::sub": lambda P, v: [v[0] - v[2], v[1] - v[3]],
    "Complex::mul": lambda P, v: [v[0] * v[2] - v[1] * v[3], v[0] * v[3] + v[1] * v[2]],
    "Complex::div": lambda P, v: cx_div(v[0], v[1], v[2], v[3]),
    "Complex::mul_f64": lambda P, v: [v[0] * v[2], v[1] * v[2]],
    "f64::mul_complex": lambda P, v: [v[0] * v[1], v[0] * v[2]],
    "Complex::div_f64": lambda P, v: [fdiv(v[0], v[2]), fdiv(v[1], v[2])],
    "Complex::is_zero": lambda P, v: [int(v[0] * v[0] + v[1] * v[1] < v[2] * v[2])],
}


def cx_fmt(re, im):
    return rust_display(re) + ("+" if im >= 0.0 else "") + rust_display(im) + "i"


def validate_complex(path: Path) -> list:
    header, rows = read_tsv(path)
    errs = check_platform(header)
    if errs:
        return errs
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [t[1] for t in ins]
        if fn == "Complex::fmt":
            if tok_text(outs[0][1]) != cx_fmt(v[0], v[1]):
                errs.append(f"complex row {n}: {fn}{tuple(v)}: Rust {tok_text(outs[0][1])!r}, "
                            f"transcription {cx_fmt(v[0], v[1])!r}")
            continue
        if fn not in COMPLEX:
            errs.append(f"complex row {n}: no transcription for {fn}")
            continue
        got = COMPLEX[fn](GLIBC_P, v)
        want = [t[1] for t in outs]
        if len(got) != len(want) or not all(
                (g == w if isinstance(w, int) else same(g, w)) for g, w in zip(got, want)):
            errs.append(f"complex row {n}: {fn}{tuple(v)} recorded {want}, transcription {got}")
    return errs


def complex_corrections(rows) -> list:
    out = []
    for fn, cls, ins, outs in rows:
        if cls != "C":
            continue
        v = [t[1] for t in ins]
        want = COMPLEX[fn](MPFR_P, v)
        old = [t[1] for t in outs]
        if all(same(a, b) for a, b in zip(want, old)):
            continue
        ins_s = " ".join(f"f:{f_to_bits(a):016x}" for a in v)
        outs_old = " ".join(f"f:{f_to_bits(a):016x}" for a in old)
        outs_new = " ".join(f"f:{f_to_bits(a):016x}" for a in want)
        out.append((fn, ins_s, outs_new, f"oracle {outs_old}; correctly rounded primitives give {outs_new}"))
    return out


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    golden = root / "tests/golden"
    checks = {
        "complex": validate_complex,
        "constants": lambda p: validate_constants(root, p),
        "error": validate_error,
        "fmtf": validate_fmtf,
        "math": validate_math,
    }
    listed = [l for l in (golden / "modules.txt").read_text().splitlines() if l and not l.startswith("#")]
    errs, total = [], 0
    for mod in listed:
        if mod.endswith(".cr"):
            continue
        path = golden / f"{mod}.tsv"
        if mod not in checks:
            errs.append(f"{mod}: no validator (add one when the module lands)")
            continue
        _h, rows = read_tsv(path)
        e = checks[mod](path)
        total += len(rows)
        print(f"{'ok ' if not e else 'BAD'}  {mod}: {len(rows)} rows" + (f", {len(e)} disagreements" if e else ""))
        errs.extend(e)
    errs += check_corrections(root, "math", golden / "math.tsv")
    errs += check_corrections(root, "complex", golden / "complex.tsv", complex_corrections)
    errs += check_tables(root)
    for e in errs[:40]:
        print("  " + e)
    if len(errs) > 40:
        print(f"  … {len(errs) - 40} more")
    print(f"validate: {total} rows, {len(errs)} problems")
    return 0 if not errs else 1


if __name__ == "__main__":
    sys.exit(main())
