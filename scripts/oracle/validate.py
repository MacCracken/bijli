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
import json
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


def cr_text(module: str, header: list, rows, fixes_of=None, source=None, columns=None) -> str:
    fixes = (fixes_of or corrections)(rows)
    if not fixes:
        return ""
    oracle = next(h for h in header if h.startswith("# oracle:")).replace(" Generated by scripts/oracle; do not edit.", "")
    source = source or f"MPFR {mpfr_ref.VERSION}"
    head = [
        f"# bijli corrected oracle (ADR 0001), format 1: {module}",
        f"{oracle} Correctly rounded: {source}. Generated by scripts/oracle/validate.py; do not edit.",
        "# regenerate: python3 -I -B scripts/oracle/validate.py . --write",
        "# columns: fn<TAB>inputs<TAB>outputs<TAB>note (fn, inputs: the golden row's, byte for byte; outputs: "
        + (columns or "the oracle's operation with its primitive correctly rounded; note: the oracle's outputs, then each "
           "misrounded call") + ")",
        f"# rows: {len(fixes)}",
    ]
    return "\n".join(head + ["\t".join(f) for f in fixes]) + "\n"


def check_corrections(root: Path, module: str, path_tsv: Path, fixes_of=None, source=None, columns=None) -> list:
    header, rows = read_tsv(path_tsv)
    text = cr_text(module, header, rows, fixes_of, source, columns)
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


def rust_debug(x: float) -> str:
    """`format!("{x:?}")` for an f64 (core's float_to_general_debug)."""
    if math.isnan(x):
        return "NaN"
    sign = "-" if math.copysign(1.0, x) < 0 else ""
    if math.isinf(x):
        return sign + "inf"
    if x == 0:
        return sign + "0.0"
    a = abs(x)
    ds, k = shortest(a, True)
    if a < 1e-4 or a >= 1e16:
        return sign + ds[0] + ("." + ds[1:] if len(ds) > 1 else "") + "e" + str(k - 1)
    return sign + positional(ds, k, True)


def debug3(p) -> str:
    """`format!("{p:?}")` for a [f64; 3]."""
    return "[" + ", ".join(rust_debug(c) for c in p) + "]"


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
        elif fn == "f64::fmt_debug":
            want = rust_debug(x)
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
        if cx_governed(fn, v):
            continue            # ADR 0002 rows 14 and 15 set this row's expectation (complex.adr.tsv)
        want = COMPLEX[fn](MPFR_P, v)
        old = [t[1] for t in outs]
        if all(same(a, b) for a, b in zip(want, old)):
            continue
        ins_s = " ".join(f"f:{f_to_bits(a):016x}" for a in v)
        outs_old = " ".join(f"f:{f_to_bits(a):016x}" for a in old)
        outs_new = " ".join(f"f:{f_to_bits(a):016x}" for a in want)
        out.append((fn, ins_s, outs_new, f"oracle {outs_old}; correctly rounded primitives give {outs_new}"))
    return out


# ── field (M1 part 3a): rust-old/src/field.rs ─────────────────────────────────────────────────
# Each function transcribed from field.rs in its order of operations, in IEEE double (Python's floats, with fdiv and
# fsqrt for the cases Python raises on). The constants are constants.tsv's, which validate_constants proves.
def golden_consts(root: Path) -> dict:
    _h, rows = read_tsv(root / "tests/golden/constants.tsv")
    return {ins[0][1]: outs[0][1] for _fn, _cls, ins, outs in rows}


def fv_text(x, y, z):
    return "(" + rust_display(x) + ", " + rust_display(y) + ", " + rust_display(z) + ")"


SINGULAR = [("E", "Singularity"), ("s", "singularity: field evaluation at source point (distance = 0)")]
R2_MIN = 1e-30
MIN_POSITIVE = 2.2250738585072014e-308


def ok3(v):
    return [("f", v[0]), ("f", v[1]), ("f", v[2])]


def point_e_k(K, q, a, b):
    """field::electric_field_point_charge: SINGULAR, or (its tokens, the vector)."""
    dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    r_sq = dx * dx + dy * dy + dz * dz
    if r_sq < R2_MIN:
        return SINGULAR
    inv_r = fdiv(1.0, fsqrt(r_sq))
    fct = K * q * inv_r * inv_r * inv_r
    return ok3((fct * dx, fct * dy, fct * dz)), (fct * dx, fct * dy, fct * dz)


def point_v_k(K, q, a, b):
    """field::electric_potential_point_charge's tokens."""
    dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    r_sq = dx * dx + dy * dy + dz * dz
    return SINGULAR if r_sq < R2_MIN else [("f", fdiv(K * q, fsqrt(r_sq)))]


def moving_b_k(MU0, q, a, vel, b):
    """field::magnetic_field_moving_charge's tokens."""
    d = [b[i] - a[i] for i in range(3)]
    r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    if r_sq < R2_MIN:
        return SINGULAR
    inv_r = fdiv(1.0, fsqrt(r_sq))
    r = [d[i] * inv_r for i in range(3)]
    cx = vel[1] * r[2] - vel[2] * r[1]
    cy = vel[2] * r[0] - vel[0] * r[2]
    cz = vel[0] * r[1] - vel[1] * r[0]
    fct = fdiv(MU0, 4.0 * PI) * q * inv_r * inv_r
    return ok3((fct * cx, fct * cy, fct * cz))


def ip_err(reason):
    return [("E", "InvalidParameter"), ("s", "invalid parameter: " + reason)]


def dz_err(context):
    return [("E", "DivisionByZero"), ("s", "division by zero in EM calculation: " + context)]


def radial_check(d):
    """field::check_radial_distance: its error's tokens, or None."""
    if math.isnan(d) or d < 0.0:
        return ip_err("distance must be non-negative, got " + rust_display(d))
    if d < R2_MIN:
        return SINGULAR
    return None


def rows_agree(fn, v, outs, want, label, n, errs):
    got = [(k, tok_text(b) if k == "s" else b) for k, b in outs]
    good = len(got) == len(want) and all(
        gk == wk and (same(gb, wb) if gk == "f" else gb == wb) for (gk, gb), (wk, wb) in zip(got, want))
    if not good:
        errs.append(f"{label} row {n}: {fn}{tuple(v)} recorded {got}, transcription {want}")


def validate_field(root: Path, path: Path) -> list:
    c = golden_consts(root)
    K, EPS0, MU0 = c["field::COULOMB_K"], c["field::EPSILON_0"], c["field::MU_0"]
    _header, rows = read_tsv(path)

    def point_e(q, a, b):
        return point_e_k(K, q, a, b)

    def rhat(p):
        r_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
        if r_sq < R2_MIN:
            return None
        inv_r = fdiv(1.0, fsqrt(r_sq))
        return inv_r, (p[0] * inv_r, p[1] * inv_r, p[2] * inv_r)

    def dipole(m, p, magnetic):
        h = rhat(p)
        if h is None:
            return SINGULAR
        inv_r, r = h
        inv_r3 = inv_r * inv_r * inv_r
        d = m[0] * r[0] + m[1] * r[1] + m[2] * r[2]
        fct = fdiv(MU0 * inv_r3, 4.0 * PI) if magnetic else fdiv(inv_r3, 4.0 * PI * EPS0)
        k = 3.0 * d * fct
        return ok3(tuple(r[i] * k - m[i] * fct for i in range(3)))

    def normalized(v):
        x, y, z = v
        ms = x * x + y * y + z * z
        if math.isfinite(ms) and ms >= MIN_POSITIVE:
            mag = math.sqrt(ms)
            return ok3((x / mag, y / mag, z / mag))
        if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
            return [("E", "InvalidParameter"), ("s", "invalid parameter: cannot normalize non-finite vector "
                                                     + fv_text(x, y, z))]
        s = rust_max(rust_max(abs(x), abs(y)), abs(z))
        if s == 0.0:
            return [("E", "DivisionByZero"), ("s", "division by zero in EM calculation: zero vector cannot be "
                                                   "normalized")]
        w = (x / s, y / s, z / s)
        mag = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
        return ok3((w[0] / mag, w[1] / mag, w[2] / mag))

    ZERO_CHARGE = [("E", "ZeroCharge"), ("s", "charge magnitude is zero — cannot compute field direction")]

    def superpose(cs, at):
        tot = (0.0, 0.0, 0.0)
        for q, c in cs:
            r = point_e(q, c, at)
            if r is SINGULAR:
                return SINGULAR
            tot = tuple(tot[i] + r[1][i] for i in range(3))
        return ok3(tot)

    def field_fn(sid, p, pos):
        """m_field.rs's field_fn scenarios: the field's tokens (ok3) or an error's."""
        if sid == 0:
            return ok3((p[0], p[1], p[2]))
        if sid == 1:
            r = point_e(p[0], p[1:4], pos)
            return r if r is SINGULAR else r[0]
        if sid == 2:
            return dipole(tuple(p[0:3]), pos, False)
        if sid == 3:
            return superpose([(p[0], (0.0, 0.0, p[1])), (p[0], (0.0, 0.0, -p[1])), (-2.0 * p[0], (0.0, 0.0, 0.0))],
                             pos)
        if sid == 4:
            return SINGULAR if pos[0] >= p[0] else ok3((1.0, 0.0, 0.0))
        if sid == 5:
            return ZERO_CHARGE
        if sid == 6:
            return ok3((math.nan, 0.0, 0.0)) if pos[0] >= p[0] else ok3((1.0, 0.0, 0.0))
        if sid == 7:
            return ok3((math.inf, 0.0, 0.0)) if pos[0] >= p[0] else ok3((1.0, 0.0, 0.0))
        return dipole(tuple(p[0:3]), pos, True)

    def trace(v, ins):
        """field.rs's trace_field_line over a scenario: the points' tokens, or the error's."""
        sid, start, step, maxw = v[0], tuple(v[1:4]), v[4], v[5] & 0xFFFFFFFFFFFFFFFF
        params = v[7:7 + int(ins[6][1])]
        if not (step > 0.0 and math.isfinite(step)):
            return [("E", "InvalidParameter"), ("s", "invalid parameter: step size must be positive and finite, got "
                                                     + rust_display(step))]
        if not all(math.isfinite(c) for c in start):
            return [("E", "InvalidParameter"), ("s", "invalid parameter: start position must be finite, got "
                                                     + debug3(start))]
        pos, pts = start, [start]
        for _ in range(maxw):
            r = field_fn(sid, params, pos)
            if r[0][0] == "E":
                if r is SINGULAR:
                    break
                return r
            fv = tuple(b for _k, b in r)
            if any(math.isnan(c) for c in fv):
                return [("E", "InvalidParameter"), ("s", "invalid parameter: field function returned NaN at "
                                                         + debug3(pos))]
            if any(math.isinf(c) for c in fv):
                break
            nd = normalized(fv)
            if nd[0][0] == "E":
                break
            d = tuple(b for _k, b in nd)
            pos = (pos[0] + step * d[0], pos[1] + step * d[1], pos[2] + step * d[2])
            pts.append(pos)
        return [("v", str(3 * len(pts)))] + [("f", c) for pt in pts for c in pt]

    errs = []
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [tok_text(b) if k == "s" else b for k, b in ins]
        V = tuple(v[0:3])
        W = tuple(v[3:6])
        w3b = field_3b(fn, v, ins, GLIBC_P, K, EPS0, MU0)
        if fn == "field::trace_field_line":
            want = trace(v, ins)
        elif w3b is not None:
            want = w3b
        elif fn == "FieldVector::zero" or fn == "FieldVector::default":
            want = ok3((0.0, 0.0, 0.0))
        elif fn in ("FieldVector::new", "FieldVector::from_array", "FieldVector::into_array"):
            want = ok3(V)
        elif fn == "FieldVector::magnitude":
            want = [("f", fsqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]))]
        elif fn == "FieldVector::magnitude_sq":
            want = [("f", V[0] * V[0] + V[1] * V[1] + V[2] * V[2])]
        elif fn == "FieldVector::normalized":
            want = normalized(V)
        elif fn == "FieldVector::dot":
            want = [("f", V[0] * W[0] + V[1] * W[1] + V[2] * W[2])]
        elif fn == "FieldVector::cross":
            want = ok3((V[1] * W[2] - V[2] * W[1], V[2] * W[0] - V[0] * W[2], V[0] * W[1] - V[1] * W[0]))
        elif fn == "FieldVector::add":
            want = ok3(tuple(V[i] + W[i] for i in range(3)))
        elif fn == "FieldVector::sub":
            want = ok3(tuple(V[i] - W[i] for i in range(3)))
        elif fn == "FieldVector::neg":
            want = ok3(tuple(-V[i] for i in range(3)))
        elif fn in ("FieldVector::scale", "FieldVector::mul"):
            want = ok3(tuple(V[i] * v[3] for i in range(3)))
        elif fn == "f64::mul_field_vector":
            want = ok3(tuple(v[0] * v[1 + i] for i in range(3)))
        elif fn == "FieldVector::fmt":
            want = [("s", fv_text(*V))]
        elif fn in ("field::electric_field_point_charge", "field::electric_potential_point_charge"):
            q, a, b = v[0], v[1:4], v[4:7]
            if fn.endswith("field_point_charge"):
                r = point_e(q, a, b)
                want = r if r is SINGULAR else r[0]
            else:
                want = point_v_k(K, q, a, b)
        elif fn == "field::magnetic_field_moving_charge":
            want = moving_b_k(MU0, v[0], v[1:4], v[4:7], v[7:10])
        elif fn == "field::electric_dipole_field":
            want = dipole(V, W, False)
        elif fn == "field::magnetic_dipole_field":
            want = dipole(V, W, True)
        elif fn == "field::electric_dipole_potential":
            h = rhat(W)
            if h is None:
                want = SINGULAR
            else:
                inv_r, r = h
                d = V[0] * r[0] + V[1] * r[1] + V[2] * r[2]
                want = [("f", fdiv(d * inv_r * inv_r, 4.0 * PI * EPS0))]
        elif fn == "field::electric_field_superposition":
            cnt = int(ins[0][1]) // 4
            cs = v[1:1 + 4 * cnt]
            at = v[1 + 4 * cnt:4 + 4 * cnt]
            tot = (0.0, 0.0, 0.0)
            want = None
            for j in range(cnt):
                r = point_e(cs[4 * j], cs[4 * j + 1:4 * j + 4], at)
                if r is SINGULAR:
                    want = SINGULAR
                    break
                tot = tuple(tot[i] + r[1][i] for i in range(3))
            if want is None:
                want = ok3(tot)
        elif fn == "field::electric_energy_density":
            want = [("f", 0.5 * EPS0 * v[0] * v[0])]
        elif fn == "field::magnetic_energy_density":
            want = [("f", fdiv(v[0] * v[0], 2.0 * MU0))]
        else:
            errs.append(f"field row {n}: no transcription for {fn}")
            continue
        rows_agree(fn, v, outs, want, "field", n, errs)
    return errs


def field_3b(fn, v, ins, P, K, EPS0, MU0):
    """Part 3b's functions transcribed from field.rs (None for another name), with the primitives P (glibc's or
    MPFR's) for the class-C ones."""
    ip = ip_err
    radial = radial_check

    def okf(x):
        return [("f", x)]

    def coords(p):
        rho_sq = p[0] * p[0] + p[1] * p[1]
        r = fsqrt(rho_sq + p[2] * p[2])
        if r == 0.0:
            return None
        return (r, rust_clamp(fdiv(p[2], r), -1.0, 1.0), rust_clamp(fdiv(fsqrt(rho_sq), r), 0.0, 1.0),
                P["atan2"](p[1], p[0]))

    def rpow(r, ll):
        return powidf2_ieee(r, ll) if ll <= 2147483647 else P["pow"](r, float(ll))

    def legendre(ll, m, x, s):
        if m > ll:
            return 0.0
        pmm, fact = 1.0, 1.0
        for _ in range(m):
            pmm *= -fact * s
            fact += 2.0
        if ll == m:
            return pmm
        mf = float(m)
        pmmp1 = x * (2.0 * mf + 1.0) * pmm
        if ll == m + 1:
            return pmmp1
        pll = 0.0
        for k in range(m + 2, ll + 1):
            kf = float(k)
            pll = fdiv(x * (2.0 * kf - 1.0) * pmmp1 - (kf + mf - 1.0) * pmm, kf - mf)
            pmm, pmmp1 = pmmp1, pll
        return pll

    def nlegendre(ll, m, x, s):
        if m > ll:
            return 0.0
        pmm = 1.0
        for k in range(1, m + 1):
            kf = float(k)
            pmm *= -s * fsqrt((2.0 * kf - 1.0) / (2.0 * kf))
        mf = float(m)
        pmm *= fsqrt((2.0 * mf + 1.0) / (4.0 * PI))
        if ll == m:
            return pmm
        pmmp1 = x * fsqrt(2.0 * mf + 3.0) * pmm
        if ll == m + 1:
            return pmmp1
        m_sq = mf * mf
        a_prev = fsqrt((4.0 * (mf + 1.0) * (mf + 1.0) - 1.0) / ((mf + 1.0) * (mf + 1.0) - m_sq))
        pll = 0.0
        for k in range(m + 2, ll + 1):
            kf = float(k)
            a = fsqrt((4.0 * kf * kf - 1.0) / (kf * kf - m_sq))
            pll = a * (x * pmmp1 - fdiv(pmm, a_prev))
            pmm, pmmp1, a_prev = pmmp1, pll, a
        return pll

    if fn == "field::electric_field_infinite_plane":
        return okf(fdiv(v[0], 2.0 * EPS0))
    if fn in ("field::electric_field_charged_sphere", "field::electric_field_infinite_line",
              "field::electric_field_outside_cylinder"):
        q, d = v[0], v[1]
        e = radial(d)
        if e:
            return e
        if fn.endswith("charged_sphere"):
            return okf(fdiv(K * q, d * d))
        return okf(fdiv(q, 2.0 * PI * EPS0 * d))
    if fn == "field::electric_field_inside_cylinder":
        q, d = v[0], v[1]
        if math.isnan(d) or d < 0.0:
            return ip("distance must be non-negative, got " + rust_display(d))
        return okf(fdiv(q * d, 2.0 * EPS0))
    if fn == "field::electric_field_inside_sphere":
        q, rr, d = v[0], v[1], v[2]
        if not (rr > 0.0 and math.isfinite(rr)):
            return ip("sphere radius must be positive and finite, got " + rust_display(rr))
        if not (d >= 0.0 and d <= rr):
            return ip(f"distance must lie in [0, R] = [0, {rust_display(rr)}] for the interior field, got "
                      + rust_display(d))
        return okf(fdiv(K * q * fdiv(d, rr), rr * rr))
    if fn == "field::electric_field_ring_axis":
        q, rr, z = v[0], v[1], v[2]
        if not (rr >= 0.0 and math.isfinite(rr)):
            return ip("ring radius must be non-negative and finite, got " + rust_display(rr))
        if not math.isfinite(z):
            return ip("axial distance must be finite, got " + rust_display(z))
        if z == 0.0:
            return SINGULAR if rr == 0.0 else okf(0.0)
        d_sq = z * z + rr * rr
        return okf(fdiv(K * q * z, d_sq * fsqrt(d_sq)))
    if fn == "field::electric_field_disk_axis":
        q, rr, z = v[0], v[1], v[2]
        if not (rr > 0.0 and math.isfinite(rr)):
            return ip("disk radius must be positive and finite, got " + rust_display(rr))
        if not math.isfinite(z):
            return ip("axial distance must be finite, got " + rust_display(z))
        r = fsqrt(z * z + rr * rr)
        geo = fdiv(rr, r) * fdiv(rr, r + abs(z))
        return okf(fdiv(math.copysign(1.0, z) * q, 2.0 * EPS0) * geo)
    if fn == "field::electric_quadrupole_potential":
        qt, p = v[0:9], v[9:12]
        r_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
        if r_sq < R2_MIN:
            return SINGULAR
        r = fsqrt(r_sq)
        r5 = r_sq * r_sq * r
        s = 0.0
        for i in range(3):
            for j in range(3):
                s += qt[3 * i + j] * p[i] * p[j]
        return okf(fdiv(s, 8.0 * PI * EPS0 * r5))
    if fn in ("field::multipole_moment", "field::spherical_multipole_moment"):
        ll, m = v[0], v[1]
        cnt = int(ins[2][1]) // 4
        cs = v[3:3 + 4 * cnt]
        ma = abs(m)
        if ma > ll:
            return ip(f"|m| must be ≤ l, got l={ll}, m={m}")
        mf = float(ma)
        if fn == "field::multipole_moment":
            mom = 0.0
            for j in range(cnt):
                q, pos = cs[4 * j], cs[4 * j + 1:4 * j + 4]
                c = coords(pos)
                if c is None:
                    if ll == 0:
                        mom += q
                    continue
                az = P["cos"](mf * c[3]) if m >= 0 else P["sin"](mf * c[3])
                mom += q * rpow(c[0], ll) * legendre(ll, ma, c[1], c[2]) * az
            return okf(mom)
        neg = 1.0 if ma % 2 == 0 else -1.0
        re = im = 0.0
        for j in range(cnt):
            q, pos = cs[4 * j], cs[4 * j + 1:4 * j + 4]
            c = coords(pos)
            if c is None:
                if ll == 0:
                    re += fdiv(q * 0.5, fsqrt(PI))
                continue
            w = q * rpow(c[0], ll) * nlegendre(ll, ma, c[1], c[2])
            sn, cn = P["sin"](mf * c[3]), P["cos"](mf * c[3])
            if m >= 0:
                re += w * cn
                im -= w * sn
            else:
                re += neg * w * cn
                im += neg * w * sn
        return [("f", re), ("f", im)]
    if fn == "field::image_charge_plane":
        q, p, pz = v[0], v[1:4], v[4]
        return [("f", -q), ("f", p[0]), ("f", p[1]), ("f", 2.0 * pz - p[2])]
    if fn == "field::image_charge_sphere":
        q, p, rr = v[0], v[1:4], v[4]
        if not (rr > 0.0 and math.isfinite(rr)):
            return ip("sphere radius must be positive and finite, got " + rust_display(rr))
        d = fsqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
        if math.isnan(d) or d <= rr:
            return ip(f"charge must be outside the sphere (|r| = {rust_display(d)} > R = {rust_display(rr)})")
        ratio = fdiv(rr, d)
        sc = ratio * ratio
        return [("f", -q * ratio), ("f", p[0] * sc), ("f", p[1] * sc), ("f", p[2] * sc)]
    if fn == "field::image_charge_dielectric":
        q, p, e1, e2 = v[0], v[1:4], v[4], v[5]
        if not (p[2] > 0.0 and math.isfinite(p[2])):
            return ip("charge must be in medium 1 (z > 0) for dielectric image method, got z = " + rust_display(p[2]))
        for eps in (e1, e2):
            if not (eps > 0.0 and math.isfinite(eps)):
                return [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(eps)} (must be > 0)")]
        return [("f", fdiv(-q * (e2 - e1), e1 + e2)), ("f", p[0]), ("f", p[1]), ("f", -p[2])]
    if fn in ("field::greens_function_scalar", "field::greens_function_static"):
        a, b = v[0:3], v[3:6]
        d = [a[i] - b[i] for i in range(3)]
        r = fsqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        if r < R2_MIN:
            return SINGULAR
        fct = fdiv(1.0, 4.0 * PI * r)
        if fn.endswith("static"):
            return okf(fct)
        kr = v[6] * r
        return [("f", fct * P["cos"](kr)), ("f", fct * P["sin"](kr))]
    if fn == "field::maxwell_stress_tensor":
        e, b = v[0:3], v[3:6]
        e_sq = e[0] * e[0] + e[1] * e[1] + e[2] * e[2]
        b_sq = b[0] * b[0] + b[1] * b[1] + b[2] * b[2]
        inv_mu0 = fdiv(1.0, MU0)
        out = []
        for i in range(3):
            for j in range(3):
                dl = 1.0 if i == j else 0.0
                out.append(("f", EPS0 * (e[i] * e[j] - 0.5 * dl * e_sq) + inv_mu0 * (b[i] * b[j] - 0.5 * dl * b_sq)))
        return out
    if fn == "field::stress_tensor_force":
        tt, n = v[0:9], v[9:12]
        out = []
        for i in range(3):
            acc = 0.0
            for j in range(3):
                acc += tt[3 * i + j] * n[j]
            out.append(("f", acc))
        return out
    if fn == "field::poynting_flux_surface":
        cnt = int(ins[0][1]) // 10
        sm = v[1:1 + 10 * cnt]
        inv_mu0 = fdiv(1.0, MU0)
        total = 0.0
        for j in range(cnt):
            e, b, n, da = sm[10 * j:10 * j + 3], sm[10 * j + 3:10 * j + 6], sm[10 * j + 6:10 * j + 9], sm[10 * j + 9]
            c = (e[1] * b[2] - e[2] * b[1], e[2] * b[0] - e[0] * b[2], e[0] * b[1] - e[1] * b[0])
            s = [c[i] * inv_mu0 for i in range(3)]
            total += (s[0] * n[0] + s[1] * n[1] + s[2] * n[2]) * da
        return okf(total)
    return None


# ── maxwell: rust-old/src/maxwell.rs ──────────────────────────────────────────────────────────
def validate_maxwell(root: Path, path: Path) -> list:
    """maxwell.rs transcribed in its order of operations (class A: + − × ÷ and sqrt), its checks and error texts."""
    c = golden_consts(root)
    EPS0, MU0 = c["field::EPSILON_0"], c["field::MU_0"]
    _header, rows = read_tsv(path)

    def material(e, m):
        if math.isnan(e) or e <= 0.0:
            return [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(e)} (must be > 0)")]
        if math.isnan(m) or m <= 0.0:
            return [("E", "InvalidPermeability"), ("s", f"invalid permeability: {rust_display(m)} (must be > 0)")]
        return None

    def positive(name, x, tail=""):
        if math.isnan(x) or x <= 0.0:
            return ip_err(f"{name} must be positive{tail}, got {rust_display(x)}")
        return None

    errs = []
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [b for _k, b in ins]
        if fn == "maxwell::gauss_electric_divergence":
            want = [("f", fdiv(v[0], EPS0))]
        elif fn == "maxwell::gauss_magnetic_divergence":
            want = [("f", 0.0)]
        elif fn == "maxwell::free_space_impedance":
            want = [("f", fsqrt(fdiv(MU0, EPS0)))]
        elif fn == "maxwell::displacement_current":
            want = ok3((v[0] * EPS0, v[1] * EPS0, v[2] * EPS0))
        elif fn in ("maxwell::wave_speed", "maxwell::refractive_index", "maxwell::impedance"):
            want = material(v[0], v[1])
            if want is None:
                if fn == "maxwell::wave_speed":
                    want = [("f", fdiv(1.0, fsqrt(v[0] * v[1])))]
                elif fn == "maxwell::refractive_index":
                    want = [("f", fsqrt(v[0] * v[1]))]
                else:
                    want = [("f", fsqrt(fdiv(v[1], v[0])))]
        elif fn in ("maxwell::wavelength", "maxwell::frequency"):
            first = "frequency" if fn == "maxwell::wavelength" else "wavelength"
            want = positive(first, v[0]) or positive("velocity", v[1]) or [("f", fdiv(v[1], v[0]))]
        elif fn == "maxwell::skin_depth":
            want = (positive("angular frequency", v[0], " for skin depth")
                    or positive("permeability", v[1], " for skin depth")
                    or positive("conductivity", v[2], " for skin depth"))
            if want is None:
                denom = v[0] * v[1] * v[2]
                want = [("f", fsqrt(fdiv(2.0, denom)))]
        else:
            errs.append(f"maxwell row {n}: no transcription for {fn}")
            continue
        rows_agree(fn, v, outs, want, "maxwell", n, errs)
    return errs


# ── charge: rust-old/src/charge.rs ────────────────────────────────────────────────────────────
def validate_charge(root: Path, path: Path) -> list:
    """charge.rs transcribed in its order of operations (class A), its checks and error texts; the PointCharge
    methods through field.rs's transcriptions above."""
    c = golden_consts(root)
    K, MU0 = c["field::COULOMB_K"], c["field::MU_0"]
    E, ME, MP = c["charge::ELEMENTARY_CHARGE"], c["charge::ELECTRON_MASS"], c["charge::PROTON_MASS"]
    _header, rows = read_tsv(path)

    def pc(q, mass, pos, vel):
        return [("f", x) for x in [q, mass, *pos, *vel]]

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

    def lorentz(q, vel, e, b):
        vb = cross(vel, b)
        return ok3(tuple((e[i] + vb[i]) * q for i in range(3)))

    def check_mass(mass, what):
        if mass == 0.0:
            return dz_err(f"mass cannot be zero for {what}")
        if math.isnan(mass) or mass <= 0.0:
            return ip_err(f"mass must be positive for {what}, got {rust_display(mass)}")
        return None

    errs = []
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [b for _k, b in ins]
        if fn == "PointCharge::new":
            mass = v[1]
            if math.isnan(mass) or mass <= 0.0:
                want = ip_err("mass must be positive, got " + rust_display(mass))
            else:
                want = pc(v[0], mass, v[2:5], v[5:8])
        elif fn == "PointCharge::electron":
            want = pc(-E, ME, v[0:3], (0.0, 0.0, 0.0))
        elif fn == "PointCharge::proton":
            want = pc(E, MP, v[0:3], (0.0, 0.0, 0.0))
        elif fn == "PointCharge::electric_field_at":
            r = point_e_k(K, v[0], v[2:5], v[8:11])
            want = r if r is SINGULAR else r[0]
        elif fn == "PointCharge::electric_potential_at":
            want = point_v_k(K, v[0], v[2:5], v[8:11])
        elif fn == "PointCharge::magnetic_field_at":
            want = moving_b_k(MU0, v[0], v[2:5], v[5:8], v[8:11])
        elif fn == "PointCharge::lorentz_force_in":
            want = lorentz(v[0], v[5:8], v[8:11], v[11:14])
        elif fn == "charge::coulomb_force":
            q1, a, q2, b = v[0], v[2:5], v[8], v[10:13]
            dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
            r_sq = dx * dx + dy * dy + dz * dz
            if r_sq < R2_MIN:
                want = SINGULAR
            else:
                inv_r = fdiv(1.0, fsqrt(r_sq))
                fct = K * q1 * q2 * inv_r * inv_r * inv_r
                want = ok3((fct * dx, fct * dy, fct * dz))
        elif fn == "charge::lorentz_force":
            want = lorentz(v[0], v[1:4], v[4:7], v[7:10])
        elif fn in ("charge::dipole_moment", "charge::potential_energy"):
            want = [("f", v[0] * v[1])]
        elif fn == "charge::coulomb_potential_energy":
            want = radial_check(v[2]) or [("f", fdiv(K * v[0] * v[1], v[2]))]
        elif fn == "charge::cyclotron_frequency":
            q, mass, bm = v
            want = check_mass(mass, "cyclotron frequency")
            if want is None and math.isnan(q):
                want = ip_err("charge must not be NaN")
            if want is None and (math.isnan(bm) or bm < 0.0):
                want = ip_err("|B| must be non-negative, got " + rust_display(bm))
            if want is None:
                want = [("f", fdiv(abs(q) * bm, mass))]
        elif fn == "charge::larmor_radius":
            mass, vp, q, bm = v
            want = check_mass(mass, "Larmor radius")
            if want is None and (math.isnan(vp) or math.isnan(q)):
                want = ip_err("v_perp and charge must not be NaN")
            if want is None and q == 0.0:
                want = dz_err("charge cannot be zero for Larmor radius")
            if want is None and bm == 0.0:
                want = dz_err("|B| cannot be zero for Larmor radius")
            if want is None and (math.isnan(bm) or bm <= 0.0):
                want = ip_err("|B| must be positive, got " + rust_display(bm))
            if want is None:
                want = [("f", fdiv(fdiv(mass * abs(vp), abs(q)), bm))]
        else:
            errs.append(f"charge row {n}: no transcription for {fn}")
            continue
        rows_agree(fn, v, outs, want, "charge", n, errs)
    return errs


def field_corrections(root: Path):
    """The field rows whose class-C outputs change when the primitives are correctly rounded (MPFR)."""
    c = golden_consts(root)
    K, EPS0, MU0 = c["field::COULOMB_K"], c["field::EPSILON_0"], c["field::MU_0"]

    def fixes(rows):
        out = []
        for fn, cls, ins, outs in rows:
            if cls != "C":
                continue
            v = [tok_text(b) if k == "s" else b for k, b in ins]
            want = field_3b(fn, v, ins, MPFR_P, K, EPS0, MU0)
            if want is None or any(k != "f" for k, _ in want):
                continue
            old = [b for k, b in outs]
            new = [b for _, b in want]
            if len(old) == len(new) and all(same(a, b) for a, b in zip(old, new)):
                continue
            raw_ins = []
            for k, b in ins:
                raw_ins.append(f"f:{f_to_bits(b):016x}" if k == "f" else f"{k}:{b}")
            outs_old = " ".join(f"f:{f_to_bits(a):016x}" for a in old)
            outs_new = " ".join(f"f:{f_to_bits(a):016x}" for a in new)
            out.append((fn, " ".join(raw_ins), outs_new, f"oracle {outs_old}; correctly rounded primitives give "
                        f"{outs_new}"))
        return out
    return fixes


# ── json (M1 part 4): serde_json 1.0.154 on Complex, FieldVector and PointCharge ─────────────────────────────────
# Recomputed from Python's json module and a transcription of serde's derive semantics, every number correctly
# rounded (Python's float): the writer rows' shape and values, the reader rows' outcome, and their `unchecked` flags
# (ADR 0002 row 12). Where serde_json's best-effort decimal parse misrounds a number, or judges an overflow on a
# misrounded value, the correctly rounded reading goes to tests/golden/json.cr.tsv (ADR 0002 row 11).

JSON_FIELDS = {
    "Complex": [("re", False), ("im", False)],
    "FieldVector": [("x", False), ("y", False), ("z", False)],
    "PointCharge": [("charge", False), ("mass", False), ("position", True), ("velocity", True)],
}


def json_writer_tokens(kind: str, xs: list) -> list:
    """The parsed-token form of serde_json's document for the fields xs: `s:` runs and `f:` numbers, with `null`
    for NaN and ±∞ inside a run."""
    parts, at = ["{"], 0
    for k, (name, arr) in enumerate(JSON_FIELDS[kind]):
        parts.append(("," if k else "") + f'"{name}":')
        if arr:
            parts += ["[", xs[at], ",", xs[at + 1], ",", xs[at + 2], "]"]
            at += 3
        else:
            parts.append(xs[at])
            at += 1
    parts.append("}")
    toks, run = [], ""
    for p in parts:
        if isinstance(p, str):
            run += p
        elif math.isfinite(p):
            if run:
                toks.append(("s", run))
                run = ""
            toks.append(("f", p))
        else:
            run += "null"
    if run:
        toks.append(("s", run))
    return toks


class JsonRefused(Exception):
    """A refused document; `overflow` when the cause is a number past f64's range."""

    def __init__(self, why: str, overflow: bool = False):
        super().__init__(why)
        self.overflow = overflow


class JsonObj(list):
    """A JSON object: its (key, value) pairs in document order. A plain list is a JSON array."""


def json_tree(doc: str):
    """The document as Python values (objects as JsonObj, numbers as ("num", text)), or JsonRefused. Python's json
    refuses what serde_json's grammar refuses (RFC 8259: JSON whitespace only, no BOM, no NaN or Infinity, no
    leading zero, no trailing comma, no raw control character in a string) and takes a lone surrogate escape, which
    serde_json refuses only where it reads the string (json_read)."""
    def constant(c):
        raise JsonRefused(f"no {c}")
    try:
        return json.loads(doc, object_pairs_hook=JsonObj, parse_int=lambda t: ("num", t),
                          parse_float=lambda t: ("num", t), parse_constant=constant)
    except (ValueError, RecursionError) as e:
        raise JsonRefused(str(e)) from None


def json_has_surrogate(s: str) -> bool:
    return any(0xD800 <= ord(c) <= 0xDFFF for c in s)


def json_f64(v) -> float:
    """An f64 field: a number, correctly rounded (`-0` is -0.0, as serde_json reads it), refused past f64's range."""
    if not (isinstance(v, tuple) and v[0] == "num"):
        raise JsonRefused("invalid type, expected f64")
    text = v[1]
    if text == "-0":
        return -0.0
    try:
        x = float(int(text)) if re.fullmatch(r"-?\d+", text) else float(text)
    except OverflowError:
        raise JsonRefused("number out of range", True) from None
    if math.isinf(x):
        raise JsonRefused("number out of range", True)
    return x


def json_read(kind: str, v) -> list:
    """serde's derive on a parsed document: an object by field name (unknown keys skipped unread, a repeated field
    refused, a key with a lone surrogate refused: serde_json parses every key it meets) or an array of exactly the
    fields; a `[f64; 3]` an array of exactly three numbers."""
    fields = JSON_FIELDS[kind]
    names = [n for n, _ in fields]
    if isinstance(v, JsonObj):
        found = {}
        for key, val in v:
            if json_has_surrogate(key):
                raise JsonRefused("lone surrogate in a key")
            if key in names:
                if key in found:
                    raise JsonRefused(f"duplicate field `{key}`")
                found[key] = val
        for n in names:
            if n not in found:
                raise JsonRefused(f"missing field `{n}`")
        vals = [found[n] for n in names]
    elif isinstance(v, list):
        if len(v) != len(fields):
            raise JsonRefused("invalid length")
        vals = v
    else:
        raise JsonRefused(f"invalid type, expected struct {kind}")
    out = []
    for (_name, arr), val in zip(fields, vals):
        if not arr:
            out.append(json_f64(val))
            continue
        if isinstance(val, JsonObj) or not isinstance(val, list):
            raise JsonRefused("invalid type, expected an array of length 3")
        if len(val) != 3:
            raise JsonRefused("invalid length, expected an array of length 3")
        out += [json_f64(e) for e in val]
    return out


def json_unchecked(tree) -> int:
    """ADR 0002 row 12's flags of a parsed document: 1 a key repeated in any object, 2 a lone surrogate in any
    string, 4 containers nested past 128 levels."""
    flags, stack = 0, [(tree, 1)]
    while stack:
        v, depth = stack.pop()
        if isinstance(v, list):
            if depth > 128:
                flags |= 4
            if isinstance(v, JsonObj):
                keys = [k for k, _ in v]
                if len(set(keys)) != len(keys):
                    flags |= 1
                if any(json_has_surrogate(k) for k in keys):
                    flags |= 2
                stack += [(x, depth + 1) for _k, x in v]
            else:
                stack += [(x, depth + 1) for x in v]
        elif isinstance(v, str) and json_has_surrogate(v):
            flags |= 2
    return flags


def json_model(kind: str, doc: str):
    """A reader row's expected reading: ("ok", the f64, flags) or ("err", the reason, overflow)."""
    try:
        tree = json_tree(doc)
        return "ok", json_read(kind, tree), json_unchecked(tree)
    except JsonRefused as e:
        return "err", str(e), e.overflow


def json_row(fn: str, ins):
    """(type, document, flags) of a reader row."""
    return fn.split("::")[1], tok_text(ins[0][1]), ins[1][1]


def validate_json(path: Path) -> list:
    _h, rows = read_tsv(path)
    notes = [l.split("\t")[4] for l in path.read_text().splitlines() if not l.startswith("#")]
    errs = []
    for n, ((fn, cls, ins, outs), note) in enumerate(zip(rows, notes), 1):
        if fn.endswith("::to_string") and cls == "A":
            want = json_writer_tokens(fn.split("::")[1], [t[1] for t in ins])
            got = [(k, tok_text(b) if k == "s" else b) for k, b in outs]
            if len(want) != len(got) or not all(w[0] == g[0] and (w[1] == g[1] if w[0] == "s" else same(w[1], g[1]))
                                                for w, g in zip(want, got)):
                errs.append(f"json row {n}: {fn} {want} != {got}")
            continue
        if not fn.endswith("::from_str") or cls != "C" or len(ins) != 2:
            errs.append(f"json row {n}: unexpected {fn} class {cls}")
            continue
        kind, doc, flags = json_row(fn, ins)
        outcome, val, extra = json_model(kind, doc)
        row_ok = not (outs and outs[0][0] == "E")
        if outcome == "ok" and row_ok:
            if flags != extra:
                errs.append(f"json row {n}: {fn} unchecked flags {flags}, Python's {extra}")
            if len(val) != len(outs):
                errs.append(f"json row {n}: {fn} {len(outs)} outputs, Python's {len(val)}")
        elif outcome == "err" and not row_ok:
            if flags != 0:
                errs.append(f"json row {n}: {fn} flags {flags} on a refused document")
        elif outcome == "ok" and "number out of range" in note:
            pass  # serde_json judged an overflow on a misrounded value: a correction
        elif outcome == "err" and extra:
            pass  # a number serde_json read overflows correctly rounded: a correction
        else:
            errs.append(f"json row {n}: {fn} oracle {'Ok' if row_ok else 'Err'}, Python {outcome} ({val}): {note[:160]}")
    return errs


def json_corrections(rows) -> list:
    """The reader rows serde_json misread, each with its correctly rounded reading: a number misrounded by its
    best-effort parse, or an overflow judged (or missed) on a misrounded value."""
    out = []
    for fn, cls, ins, outs in rows:
        if cls != "C" or not fn.endswith("::from_str"):
            continue
        kind, doc, flags = json_row(fn, ins)
        outcome, val, _extra = json_model(kind, doc)
        row_ok = not (outs and outs[0][0] == "E")
        old = " ".join(f"f:{f_to_bits(t[1]):016x}" if t[0] == "f" else f"{t[0]}:{t[1]}" for t in outs)
        new = " ".join(f"f:{f_to_bits(x):016x}" for x in val) if outcome == "ok" else "E:Deserialize"
        if outcome == "ok" and row_ok and all(same(a, b[1]) for a, b in zip(val, outs)):
            continue
        if outcome == "err" and not row_ok:
            continue
        why = {(True, True): "serde_json's best-effort parse misrounds",
               (True, False): "serde_json refuses it as out of range; every number is finite correctly rounded",
               (False, True): "a number serde_json read overflows correctly rounded"}[(outcome == "ok", row_ok)]
        out.append((fn, f"s:{ins[0][1]} i:{flags}", new, f"oracle {old}; {why}"))
    return out


# ── ADR 0002 mappings (M1 part 5): tests/golden/<module>.adr.tsv ─────────────────────────────────────────────────────
# A golden row whose port result an ADR 0002 row changes, with the port's expected outputs: fn<TAB>inputs<TAB>outputs
# <TAB>ADR 0002 row<TAB>note. Each expectation is computed exactly (scripts/oracle/exact.py), never from the port.
# Row 13 (the user's "round once" decision of 2026-11-12): one product, quotient or root of inputs and constants,
# correctly rounded where its operands are finite; an infinite or NaN operand, and every error, keeps Rust's result.

from exact import rn_mono, rn_sqrt  # noqa: E402


def _fin(*xs) -> bool:
    return all(math.isfinite(x) for x in xs)


def adr13(root: Path):
    """fn -> (inputs -> the port's f64, or None where the oracle's value stands)."""
    c = golden_consts(root)
    K, EPS0, MU0 = c["field::COULOMB_K"], c["field::EPSILON_0"], c["field::MU_0"]
    F = Fraction
    inf = float("inf")
    return {
        "maxwell::wave_speed": lambda x: 0.0 if inf in x else rn_sqrt(1 / (F(x[0]) * F(x[1]))),
        "maxwell::refractive_index": lambda x: inf if inf in x else rn_sqrt(F(x[0]) * F(x[1])),
        "maxwell::impedance": lambda x: None if inf in x else rn_sqrt(F(x[1]) / F(x[0])),
        "maxwell::free_space_impedance": lambda x: rn_sqrt(F(MU0) / F(EPS0)),
        "maxwell::skin_depth": lambda x: 0.0 if inf in x else rn_sqrt(2 / (F(x[0]) * F(x[1]) * F(x[2]))),
        "charge::coulomb_potential_energy": lambda x: rn_mono([K, x[0], x[1]], [x[2]]) if _fin(*x) else None,
        "charge::cyclotron_frequency": lambda x: rn_mono([abs(x[0]), x[2]], [x[1]]) if _fin(*x) else None,
        "charge::larmor_radius": lambda x: rn_mono([x[0], abs(x[1])], [abs(x[2]), x[3]]) if _fin(*x) else None,
        "field::electric_energy_density": lambda x: rn_mono([0.5, EPS0, x[0], x[0]], []) if _fin(*x) else None,
        "field::magnetic_energy_density": lambda x: rn_mono([x[0], x[0]], [2.0 * MU0]) if _fin(*x) else None,
        "field::electric_field_charged_sphere": lambda x: rn_mono([K, x[0]], [x[1], x[1]]) if _fin(*x) else None,
        "field::electric_field_inside_sphere":
            lambda x: rn_mono([K, x[0], x[2]], [x[1], x[1], x[1]]) if _fin(x[0]) else None,
        "field::electric_field_infinite_line":
            lambda x: rn_mono([x[0]], [2.0 * PI, EPS0, x[1]]) if _fin(*x) else None,
        "field::electric_field_outside_cylinder":
            lambda x: rn_mono([x[0]], [2.0 * PI, EPS0, x[1]]) if _fin(*x) else None,
        "field::electric_field_inside_cylinder": lambda x: rn_mono([x[0], x[1]], [2.0 * EPS0]) if _fin(*x) else None,
    }


def adr_rows(root: Path, module: str, path: Path) -> list:
    """(fn, raw inputs, outputs, ADR 0002 row, note) for each golden row of `module` an ADR 0002 row changes."""
    if module == "complex":
        return cx_rows(path)[0]
    maps = {13: adr13(root)}
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        for row, fns in maps.items():
            if fn not in fns:
                continue
            if outs.startswith("E:"):
                continue                                # an error keeps Rust's
            xs = [] if ins == "-" else [bits_to_f(t[2:]) for t in ins.split(" ")]
            want = fns[fn](xs)
            if want is None:
                continue
            got = outs.split(" ")[0]
            new = f"f:{f_to_bits(want):016x}"
            if got != new:
                out.append((fn, ins, new, str(row), f"oracle {outs}; ADR 0002 row {row}: rounded once"))
    return out


# Rows 14 and 15 (M1 part 5c): BjComplex's range edges. At an edge (a finite, nonzero component outside
# [2^-511, 2^511)) hisab's chain runs with an unbounded exponent, each operation rounded to 53 bits, the last one once
# into binary64 (IEEE's signed zeros throughout); in range the chain IS hisab's, which cx_model checks row by row.
# ln and powf take the unbounded |z| where it is a normal double. powf has no cutoff (row 15).
from exact import binade, rn_frac, rn_sqrt, rnd53, rnd53_sqrt  # noqa: E402

UX_LO, UX_HI = 2.0 ** -511, 2.0 ** 511


def cx_ok(x: float) -> bool:
    return x == 0 or UX_LO <= abs(x) < UX_HI


def cx_plain(x: float, y: float) -> bool:
    return not (math.isfinite(x) and math.isfinite(y)) or (cx_ok(x) and cx_ok(y))


def ux(x: float):
    """An exact value with its zero's sign: (Fraction, negative)."""
    return Fraction(x), math.copysign(1.0, x) < 0


def ux_mul(a, b, r=rnd53):
    v = a[0] * b[0]
    return r(v), (a[1] != b[1]) if v == 0 else v < 0


def ux_div(a, b, r=rnd53):
    if b[0] == 0:                                       # IEEE's x/±0: NaN for 0/0, else ±∞
        if a[0] == 0:
            return float("nan"), False
        return (-math.inf if (a[0] < 0) != b[1] else math.inf), False
    v = a[0] / b[0]
    return r(v), (a[1] != b[1]) if v == 0 else v < 0


def ux_add(a, b, r=rnd53):
    v = a[0] + b[0]
    if v == 0:
        both = a[0] == 0 and b[0] == 0
        return Fraction(0), (a[1] and b[1]) if both else False
    return r(v), v < 0


def ux_sub(a, b, r=rnd53):
    return ux_add(a, (-b[0], not b[1]), r)


def ux_sqrt(a):
    return rnd53_sqrt(a[0]), a[1] if a[0] == 0 else False


def ux_fin(a) -> float:
    """The value rounded once into binary64, a zero keeping its sign."""
    if isinstance(a[0], float):
        return a[0]
    if a[0] == 0:
        return -0.0 if a[1] else 0.0
    return rn_frac(a[0])


def ux_fin_sqrt(a) -> float:
    if a[0] == 0:
        return -0.0 if a[1] else 0.0
    return rn_sqrt(a[0])


def ux_sumsq(x, y):
    return ux_add(ux_mul(ux(x), ux(x)), ux_mul(ux(y), ux(y)))


def cx_abs_for(re, im):
    """|z| for ln and powf: hisab's where the components are in range or one isn't finite; else the unbounded |z|
    where it is a normal double, else hisab's."""
    if cx_plain(re, im):
        return cx_abs(re, im)
    a = ux_sqrt(ux_sumsq(re, im))
    if a[0] != 0 and -1022 <= binade(a[0]) <= 1023:
        return float(a[0])
    return cx_abs(re, im)


def cx_model(fn, v):
    """The port's outputs for a complex row with finite inputs, under rows 14 and 15 (MPFR's correctly rounded
    primitives for the transcendentals), or None for a function the rows don't touch."""
    P = MPFR_P
    if fn in ("Complex::abs", "ComplexExt::norm"):
        return [ux_fin_sqrt(ux_sumsq(v[0], v[1]))]
    if fn == "Complex::norm_sq":
        return [ux_fin(ux_add(ux_mul(ux(v[0]), ux(v[0])), ux_mul(ux(v[1]), ux(v[1])), lambda t: t))]
    if fn == "Complex::mul":
        a, b, c, d = (ux(x) for x in v)
        return [ux_fin(ux_sub(ux_mul(a, c), ux_mul(b, d), lambda t: t)),
                ux_fin(ux_add(ux_mul(a, d), ux_mul(b, c), lambda t: t))]
    if fn == "Complex::inv":
        d = ux_sumsq(v[0], v[1])
        return [ux_fin(ux_div(ux(v[0]), d, lambda t: t)), ux_fin(ux_div(ux(-v[1]), d, lambda t: t))]
    if fn == "Complex::div":
        a, b, c, dd = (ux(x) for x in v)
        den = ux_sumsq(v[2], v[3])
        nr = ux_add(ux_mul(a, c), ux_mul(b, dd))
        ni = ux_sub(ux_mul(b, c), ux_mul(a, dd))
        return [ux_fin(ux_div(nr, den, lambda t: t)), ux_fin(ux_div(ni, den, lambda t: t))]
    if fn == "Complex::is_zero":
        ns = ux_sumsq(v[0], v[1])[0]
        t2 = ux_mul(ux(v[2]), ux(v[2]))[0]
        return [int(ns < t2)]
    if fn == "Complex::sqrt":
        theta = P["atan2"](v[1], v[0]) / 2.0
        if cx_plain(v[0], v[1]):
            return cx_polar(P, fsqrt(cx_abs(v[0], v[1])), theta)
        r = ux_sqrt(ux_sqrt(ux_sumsq(v[0], v[1])))
        return [ux_fin(ux_mul(r, ux(P["cos"](theta)), lambda t: t)),
                ux_fin(ux_mul(r, ux(P["sin"](theta)), lambda t: t))]
    if fn == "Complex::ln":
        return [P["log"](cx_abs_for(v[0], v[1])), P["atan2"](v[1], v[0])]
    if fn == "Complex::powf":
        re, im, n = v
        if re == 0 and im == 0:
            if n != n:
                return [n, n]
            return [0.0 if n > 0 else (1.0 if n == 0 else float("inf")), 0.0]
        a = cx_abs_for(re, im)
        return cx_polar(P, P["pow"](a, n), P["atan2"](im, re) * n)
    return None


def cx_governed(fn, v) -> int:
    """The ADR 0002 row whose model sets this complex row's expectation (14 or 15), or 0 where hisab's chain (with
    the corrected oracle's primitives) still does."""
    if fn not in ("Complex::abs", "ComplexExt::norm", "Complex::norm_sq", "Complex::mul", "Complex::inv",
                  "Complex::div", "Complex::is_zero", "Complex::sqrt", "Complex::ln", "Complex::powf"):
        return 0
    if not all(math.isfinite(x) for x in v):
        return 0
    if fn == "Complex::powf":
        re, im, _n = v
        if (re == 0 and im == 0) or re * re + im * im < TINY_300:
            return 15
        return 14 if cx_abs_for(re, im) != cx_abs(re, im) else 0
    if fn in ("Complex::mul", "Complex::div"):
        return 0 if all(cx_ok(x) for x in v) else 14
    if fn == "Complex::is_zero":
        return 0 if all(cx_ok(x) for x in v) else 14
    if fn == "Complex::ln":
        return 14 if cx_abs_for(v[0], v[1]) != cx_abs(v[0], v[1]) else 0
    return 0 if cx_plain(v[0], v[1]) else 14


def cx_rows(path: Path):
    """(adr rows, problems) for complex.tsv: each governed row whose model differs from the oracle's outputs; and,
    for every in-scope row with finite inputs that hisab's chain still sets, the check that the unbounded model gives
    the same value (the guard's proof, row by row)."""
    out, errs = [], []
    for n, line in enumerate((l for l in path.read_text().splitlines() if not l.startswith("#")), 1):
        fn, _cls, ins, outs, _note = line.split("\t")
        v = [] if ins == "-" else [parse_tok(t)[1] for t in ins.split(" ")]
        model = cx_model(fn, v) if all(isinstance(x, float) and math.isfinite(x) for x in v) else None
        if model is None:
            continue
        row = cx_governed(fn, v)
        old = [parse_tok(t)[1] for t in outs.split(" ")]
        if row == 0:
            want = COMPLEX[fn](MPFR_P, v)
            if not all((a == b if isinstance(b, int) else same(a, b)) for a, b in zip(model, want)):
                errs.append(f"complex row {n}: {fn}{tuple(v)}: the unbounded model {model} differs from hisab's "
                            f"chain {want} where the port keeps the chain")
            continue
        if all((a == b if isinstance(b, int) else same(a, b)) for a, b in zip(model, old)):
            continue
        toks = " ".join(f"i:{x}" if isinstance(x, int) else f"f:{f_to_bits(x):016x}" for x in model)
        out.append((fn, ins, toks, str(row), f"oracle {outs}; ADR 0002 row {row}"))
    return out, errs


def check_adr(root: Path, module: str, path: Path) -> list:
    fixes = adr_rows(root, module, path)
    target = path.with_name(f"{module}.adr.tsv")
    rel = f"tests/golden/{module}.adr.tsv"
    text = ""
    if fixes:
        head = [
            f"# bijli ADR 0002 mappings, format 1: {module}",
            "# Generated by scripts/oracle/validate.py from exact values (scripts/oracle/exact.py); do not edit.",
            "# regenerate: python3 -I -B scripts/oracle/validate.py . --write",
            "# columns: fn<TAB>inputs<TAB>outputs<TAB>row<TAB>note (fn, inputs: the golden row's, byte for byte; outputs:"
            " the port's, as the ADR 0002 row gives them; row: that row's number; note: the oracle's outputs)",
            f"# rows: {len(fixes)}",
        ]
        text = "\n".join(head + ["\t".join(f) for f in fixes]) + "\n"
    if WRITE:
        if text:
            target.write_text(text)
            print(f"wrote {rel} ({len(fixes)} rows)")
        elif target.exists():
            target.unlink()
            print(f"removed {rel}")
        return []
    if not text:
        return [f"{rel}: no row needs a mapping, but the file exists"] if target.exists() else []
    if not target.exists() or target.read_text() != text:
        return [f"{rel}: differs from the ADR 0002 mappings (regenerate with --write)"]
    print(f"ok   {module}.adr.tsv: {len(fixes)} mapped rows")
    return []


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    golden = root / "tests/golden"
    checks = {
        "charge": lambda p: validate_charge(root, p),
        "complex": validate_complex,
        "constants": lambda p: validate_constants(root, p),
        "error": validate_error,
        "field": lambda p: validate_field(root, p),
        "fmtf": validate_fmtf,
        "json": validate_json,
        "math": validate_math,
        "maxwell": lambda p: validate_maxwell(root, p),
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
    errs += check_corrections(root, "field", golden / "field.tsv", field_corrections(root))
    errs += check_corrections(root, "json", golden / "json.tsv", json_corrections, "Python's float",
                              "the document read with every number correctly rounded; note: the oracle's outputs, "
                              "then why they differ")
    for mod in ("charge", "complex", "field", "maxwell"):
        errs += check_adr(root, mod, golden / f"{mod}.tsv")
    errs += cx_rows(golden / "complex.tsv")[1]
    errs += check_tables(root)
    for e in errs[:40]:
        print("  " + e)
    if len(errs) > 40:
        print(f"  … {len(errs) - 40} more")
    print(f"validate: {total} rows, {len(errs)} problems")
    return 0 if not errs else 1


if __name__ == "__main__":
    sys.exit(main())
