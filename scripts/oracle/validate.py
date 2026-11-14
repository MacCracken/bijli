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
- **src/crmath_data.cyr**, crtables.py's rendering, with its π and ln 2 and every table's hi checked against MPFR;
  and **src/crmath_ux_data.cyr**, its render_ux (bijli's longer 2/π), from π checked to 2,800 bits.

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
    upath = root / "src/crmath_ux_data.cyr"
    utext = crtables.render_ux()
    if crtables.pi_bits(crtables.UX_PREC) != mpfr_ref.const_bits("pi", crtables.UX_PREC):
        errs.append(f"crtables: π to {crtables.UX_PREC} bits (crmath_ux_data's 2/π) differs from MPFR's")
    if WRITE:
        upath.write_text(utext)
        print(f"wrote src/crmath_ux_data.cyr ({len(utext.splitlines())} lines)")
    elif not upath.exists() or upath.read_text() != utext:
        errs.append("src/crmath_ux_data.cyr is not crtables.py's render_ux output (regenerate with --write)")
    else:
        print(f"ok   crmath_ux_data: 2/π from π to {crtables.UX_PREC} bits, equal to MPFR's")
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


def rust_fixed(x: float, p: int) -> str:
    """`format!("{x:.p}")` for an f64: the exact value rounded half to even at p decimals, as Python's `f` format
    rounds it; Rust spells NaN and infinity its own way, and keeps the sign of a zero or of a value that rounds to
    one."""
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return f"{x:.{p}f}"


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
        elif fn == "f64::fmt_fixed":
            want = rust_fixed(x, ins[1][1])
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
           "pow": LIBM.pow, "sinh": LIBM.sinh, "cosh": LIBM.cosh, "asin": LIBM.asin, "atan": LIBM.atan}
MPFR_P = {"sin": mpfr_ref.sin, "cos": mpfr_ref.cos, "exp": mpfr_ref.exp, "log": lambda x: mpfr_ref.cr("log", x),
          "atan2": mpfr_ref.atan2, "pow": mpfr_ref.pow, "sinh": lambda x: mpfr_ref.cr("sinh", x),
          "cosh": lambda x: mpfr_ref.cr("cosh", x), "asin": lambda x: mpfr_ref.cr("asin", x),
          "atan": lambda x: mpfr_ref.cr("atan", x)}
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
            continue            # ADR 0002 rows 14 to 16 set this row's expectation (complex.adr.tsv)
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


# ADR 0002 row 17 (M1 part 5d-1): the Singularity tests on the exact squared distance of the inputs, against 10^-30
# (r² < 1e-30) or 10^-60 (R < 1e-30), the user's decision of 2026-11-12. The transcriptions below call cut_r2 and cut_r;
# with EXACT_THR off they are Rust's comparisons (validate_field and validate_charge check the oracle with it off), and
# adr17_rows turns it on to find the rows whose outcome the exact test changes. A non-finite point keeps Rust's test.
EXACT_THR = False
ORIGIN = (0.0, 0.0, 0.0)


def exact_sq(a, b):
    return sum((Fraction(x) - Fraction(y)) ** 2 for x, y in zip(a, b))


def cut_r2(r_sq, a, b) -> bool:
    if EXACT_THR and all(math.isfinite(x) for x in (*a, *b)):
        return exact_sq(a, b) < Fraction(1, 10 ** 30)
    return r_sq < R2_MIN


def cut_r(r, a, b) -> bool:
    if EXACT_THR and all(math.isfinite(x) for x in (*a, *b)):
        return exact_sq(a, b) < Fraction(1, 10 ** 60)
    return r < R2_MIN


def point_e_k(K, q, a, b):
    """field::electric_field_point_charge: SINGULAR, or (its tokens, the vector)."""
    dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    r_sq = dx * dx + dy * dy + dz * dz
    if cut_r2(r_sq, a, b):
        return SINGULAR
    inv_r = fdiv(1.0, fsqrt(r_sq))
    fct = K * q * inv_r * inv_r * inv_r
    return ok3((fct * dx, fct * dy, fct * dz)), (fct * dx, fct * dy, fct * dz)


def point_v_k(K, q, a, b):
    """field::electric_potential_point_charge's tokens."""
    dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    r_sq = dx * dx + dy * dy + dz * dz
    return SINGULAR if cut_r2(r_sq, a, b) else [("f", fdiv(K * q, fsqrt(r_sq)))]


def moving_b_k(MU0, q, a, vel, b):
    """field::magnetic_field_moving_charge's tokens."""
    d = [b[i] - a[i] for i in range(3)]
    r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    if cut_r2(r_sq, a, b):
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


def validate_field(root: Path, path: Path, sink=None, P=None) -> list:
    """field.rs transcribed (Rust's thresholds unless EXACT_THR). With a sink list, each row's expected tokens (None
    where nothing transcribes it) go to the sink instead of being checked, the class-C primitives P (glibc's by
    default)."""
    P = P or GLIBC_P
    c = golden_consts(root)
    K, EPS0, MU0 = c["field::COULOMB_K"], c["field::EPSILON_0"], c["field::MU_0"]
    _header, rows = read_tsv(path)

    def point_e(q, a, b):
        return point_e_k(K, q, a, b)

    def rhat(p):
        r_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
        if cut_r2(r_sq, p, ORIGIN):
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

    def trace(v, ins, port=False):
        """field.rs's trace_field_line over a scenario: the points' tokens, or the error's. With port, the port's
        (row 18): (its tokens, whether a field or a step left Rust's chain, the guard proof's problems)."""
        sid, start, step, maxw = v[0], tuple(v[1:4]), v[4], v[5] & 0xFFFFFFFFFFFFFFFF
        params = v[7:7 + int(ins[6][1])]
        if not (step > 0.0 and math.isfinite(step)):
            return [("E", "InvalidParameter"), ("s", "invalid parameter: step size must be positive and finite, got "
                                                     + rust_display(step))]
        if not all(math.isfinite(c) for c in start):
            return [("E", "InvalidParameter"), ("s", "invalid parameter: start position must be finite, got "
                                                     + debug3(start))]
        pos, pts = start, [start]
        gov, perr = False, []
        for _ in range(maxw):
            r = field_fn(sid, params, pos)
            if port:
                r, g = m5_trace_field(sid, params, pos, r)
                gov = gov or g
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
            rp = (pos[0] + step * d[0], pos[1] + step * d[1], pos[2] + step * d[2])
            if port:
                mp = m5_step(pos, fv, step)
                if not m5_step_plain(fv, step):
                    gov, rp = True, mp
                elif not toks_same(ok3(mp), ok3(rp)):
                    perr.append(f"the step from {pos} along {fv}: the unbounded step {mp} differs from Rust's {rp} "
                                f"where the port keeps Rust's chain")
            pos = rp
            pts.append(pos)
        toks = [("v", str(3 * len(pts)))] + [("f", c) for pt in pts for c in pt]
        return (toks, gov, perr) if port else toks

    global TRACE_FN
    TRACE_FN = trace

    errs = []
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [tok_text(b) if k == "s" else b for k, b in ins]
        V = tuple(v[0:3])
        W = tuple(v[3:6])
        w3b = field_3b(fn, v, ins, P, K, EPS0, MU0)
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
            if sink is not None:
                sink.append(None)
                continue
            errs.append(f"field row {n}: no transcription for {fn}")
            continue
        if sink is not None:
            sink.append(want)
            continue
        rows_agree(fn, v, outs, want, "field", n, errs)
    lo = Fraction(math.nextafter(R2_MIN, 0.0))
    if not (lo < Fraction(1, 10 ** 30) < Fraction(R2_MIN)):
        errs.append("field: 10^-30 is not between the double 1e-30 and its predecessor, so `distance < 1e-30` is not "
                    "the exact radial test (ADR 0002 row 17)")
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
        if cut_r2(r_sq, p, ORIGIN):
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
        if cut_r(r, a, b):
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


# ── material (M2 part 1a): rust-old/src/material.rs ───────────────────────────────────────────────────────────────
def material_rust(fn, v, EPS0, MU0):
    """material.rs transcribed in its order of operations (class A: + − × ÷ and sqrt), its checks and error texts;
    None for another name."""
    def perm_e(e):
        return [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(e)} (must be > 0)")]

    def perm_m(m):
        return [("E", "InvalidPermeability"), ("s", f"invalid permeability: {rust_display(m)} (must be > 0)")]

    def lossless(e, m, s):
        if not (math.isfinite(e) and e > 0.0):
            return perm_e(e), None
        if not (math.isfinite(m) and m > 0.0):
            return perm_m(m), None
        if s != 0.0:
            return ip_err(f"Fresnel coefficients from Material are lossless-only; conductivity σ = {rust_display(s)} "
                          "S/m needs a frequency-dependent complex refractive index"), None
        return None, fsqrt(fdiv(m, e))

    if fn == "material::electric_susceptibility" or fn == "material::magnetic_susceptibility":
        return [("f", v[0] - 1.0)]
    if fn == "material::relative_permittivity" or fn == "material::relative_permeability":
        return [("f", 1.0 + v[0])]
    if fn == "material::absolute_permittivity":
        return [("f", v[0] * EPS0)]
    if fn == "material::absolute_permeability":
        return [("f", v[0] * MU0)]
    if fn in ("material::polarization", "material::displacement_field", "material::magnetization",
              "material::b_field_from_h", "material::h_field_from_b"):
        x, w = v[0], v[1:4]
        if fn == "material::polarization":
            k = (x - 1.0) * EPS0
        elif fn == "material::displacement_field":
            k = x * EPS0
        elif fn == "material::magnetization":
            k = x - 1.0
        elif fn == "material::b_field_from_h":
            k = x * MU0
        else:
            if abs(x) < 1e-30:
                return dz_err("relative permeability cannot be zero")
            k = fdiv(1.0, x * MU0)
        return ok3(tuple(w[i] * k for i in range(3)))
    if fn == "material::bound_surface_charge":
        a, b = v[0:3], v[3:6]
        return [("f", a[0] * b[0] + a[1] * b[1] + a[2] * b[2])]
    if fn == "material::bound_surface_current":
        a, b = v[0:3], v[3:6]
        return ok3((a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]))
    if fn == "material::bound_volume_charge":
        return [("f", -(v[0] + v[1] + v[2]))]
    if fn == "material::dielectric_energy_density":
        return [("f", 0.5 * v[0] * EPS0 * v[1] * v[1])]
    if fn == "material::clausius_mossotti":
        a, n = v[0], v[1]
        if not math.isfinite(a) or not math.isfinite(n) or n < 0.0:
            return ip_err("polarizability and number density must be finite with N ≥ 0, got α = "
                          f"{rust_display(a)}, N = {rust_display(n)}")
        chi = fdiv(n * a, 3.0 * EPS0)
        if math.isnan(chi) or chi >= 1.0:
            return ip_err(f"Nα/(3ε₀) = {rust_display(chi)} ≥ 1: Clausius-Mossotti diverges (polarization "
                          "catastrophe), no physical ε_r")
        eps = fdiv(1.0 + 2.0 * chi, 1.0 - chi)
        if math.isnan(eps) or eps <= 0.0:
            return perm_e(eps)
        return [("f", eps)]
    if fn == "material::magnetic_energy_density_material":
        mu = v[0] * MU0
        if abs(mu) < 1e-30:
            return dz_err("permeability cannot be zero for energy density")
        return [("f", fdiv(v[1] * v[1], 2.0 * mu))]
    if fn == "material::curie_law":
        if v[1] <= 0.0:
            return ip_err(f"temperature must be positive, got {rust_display(v[1])} K")
        return [("f", fdiv(v[0], v[1]))]
    if fn == "material::curie_weiss":
        d = v[1] - v[2]
        if d <= 0.0:
            return ip_err(f"temperature ({rust_display(v[1])} K) must be above Curie temperature "
                          f"({rust_display(v[2])} K)")
        return [("f", fdiv(v[0], d))]
    if fn == "material::classify_magnetic":
        x = v[0]
        name = "Diamagnetic" if x < -1e-15 else "Ferromagnetic" if x > 1.0 else "Paramagnetic"
        return [("e", "MagneticType." + name)]
    if fn in ("Material::vacuum", "Material::default"):
        return [("f", 1.0), ("f", 1.0), ("f", 0.0)]
    if fn == "Material::dielectric":
        return [("f", v[0]), ("f", 1.0), ("f", 0.0)]
    if fn == "Material::conductor":
        return [("f", 1.0), ("f", 1.0), ("f", v[0])]
    if fn == "Material::new":
        return [("f", v[0]), ("f", v[1]), ("f", v[2])]
    e, m, s = v[0:3] if len(v) >= 3 else (None, None, None)
    if fn == "Material::refractive_index":
        return [("f", fsqrt(e * m))]
    if fn == "Material::permittivity":
        return [("f", e * EPS0)]
    if fn == "Material::permeability":
        return [("f", m * MU0)]
    if fn == "Material::impedance":
        if not (math.isfinite(e) and e > 0.0):
            return perm_e(e)
        if not (math.isfinite(m) and m > 0.0):
            return perm_m(m)
        return [("f", fsqrt(fdiv(m * MU0, e * EPS0)))]
    if fn == "Material::wave_speed":
        if not (e > 0.0 and m > 0.0):
            return ip_err("ε_r and μ_r must be positive")
        return [("f", fdiv(1.0, fsqrt(e * EPS0 * m * MU0)))]
    if fn == "Material::is_lossy":
        return [("i", 1 if s > 0.0 else 0)]
    if fn == "Material::try_reflectance_at":
        r1, eta1 = lossless(*v[0:3])
        if r1 is not None:
            return r1
        r2, eta2 = lossless(*v[3:6])
        if r2 is not None:
            return r2
        r = fdiv(eta2 - eta1, eta2 + eta1)
        return [("f", r * r)]
    return None


def validate_material(root: Path, path: Path) -> list:
    c = golden_consts(root)
    EPS0, MU0 = c["field::EPSILON_0"], c["field::MU_0"]
    _header, rows = read_tsv(path)
    errs = []
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [b for _k, b in ins]
        want = material_rust(fn, v, EPS0, MU0)
        if want is None:
            errs.append(f"material row {n}: no transcription for {fn}")
            continue
        rows_agree(fn, v, outs, want, "material", n, errs)
    return errs


# ADR 0002 rows 19–21 (M2 part 1a, the user's 2026-11-12 answers): the port's material.rs, where it differs from
# Rust's (and row 18 through FieldVector's dot and cross). Row 19: the Curie laws refuse a NaN temperature (`!(T > 0)`, `!(T − T_c > 0)`) with their own texts. Row 20:
# Material::wave_speed validates as Material::impedance does (ε_r, then μ_r, positive and finite; InvalidPermittivity /
# InvalidPermeability with the value). Row 21: classify_magnetic returns a Result, InvalidParameter for NaN.
# Row 13 (M2 part 2b-1, the user's "round once" decision of 2026-11-12): material.rs's products, quotients and roots of
# inputs and constants, correctly rounded where the operands are finite (exact values: Fraction, an integer square
# root); a non-finite operand keeps Rust's chain, and so does every error. D, B and H part by part.
def mat_cr_part(x, c, w, k, mode):
    """x·c·w (mode 0) or w/(x·c) (mode 1) rounded once for a finite w; Rust's w·k otherwise."""
    if not math.isfinite(w):
        return w * k
    return rn_mono([x, c, w], []) if mode == 0 else rn_mono([w], [x, c])


def mat_lossless_port(e, m, s):
    """The port's lossless optical constants: (error tokens, None) or (None, (n, η)), n and η correctly rounded."""
    err, _l = wave_lossless(e, m, s)
    if err:
        return err, None
    F = Fraction
    return None, (rn_sqrt(F(e) * F(m)), rn_sqrt(F(m) / F(e)))


def mat_reflectance(v, lossless, ux=False):
    """Material::try_reflectance_at on a lossless-constants function; ux: Γ as wave_gamma_port (row 18's rule, M2 part
    2d-2c)."""
    r1, l1 = lossless(*v[0:3])
    if r1:
        return r1
    r2, l2 = lossless(*v[3:6])
    if r2:
        return r2
    eta1, eta2 = l1[1], l2[1]
    r = wave_gamma_port(eta1, eta2, v[0:2], v[3:5]) if ux else fdiv(eta2 - eta1, eta2 + eta1)
    return [("f", r * r)]


def material_port13(fn, v, EPS0, MU0):
    """Row 13's port outputs for a material row, or None where Rust's result stands."""
    F = Fraction
    if fn in ("material::displacement_field", "material::b_field_from_h", "material::h_field_from_b"):
        x, w = v[0], v[1:4]
        if not math.isfinite(x):
            return None
        if fn == "material::h_field_from_b":
            if abs(x) < 1e-30:
                return None
            k, c, mode = fdiv(1.0, x * MU0), MU0, 1
        else:
            c = EPS0 if fn == "material::displacement_field" else MU0
            k, mode = x * c, 0
        return ok3(tuple(mat_cr_part(x, c, w[i], k, mode) for i in range(3)))
    if fn == "material::dielectric_energy_density":
        return okf(rn_mono([0.5, v[0], EPS0, v[1], v[1]], [])) if _fin(v[0], v[1]) else None
    if fn == "material::magnetic_energy_density_material":
        if not _fin(v[0], v[1]) or below30(Fraction(v[0]) * Fraction(MU0)):
            return None
        return okf(rn_mono([v[1], v[1]], [v[0], 2.0 * MU0]))
    if fn == "Material::refractive_index":
        e, m = v[0], v[1]
        if not _fin(e, m) or e == 0.0 or m == 0.0:
            return None
        if (e < 0.0) != (m < 0.0):
            return okf(math.nan)
        return okf(rn_sqrt(F(abs(e)) * F(abs(m))))
    if fn == "Material::impedance" or fn == "Material::wave_speed":
        e, m = v[0], v[1]
        if not (math.isfinite(e) and e > 0.0 and math.isfinite(m) and m > 0.0):
            return None
        if fn == "Material::impedance":
            return okf(rn_sqrt(F(m) * F(MU0) / (F(e) * F(EPS0))))
        return okf(rn_sqrt(1 / (F(e) * F(EPS0) * F(m) * F(MU0))))
    if fn == "Material::try_reflectance_at":
        return mat_reflectance(v, mat_lossless_port)
    return None


def material_cm_port(a, n, EPS0):
    """Clausius–Mossotti in the port (M2 part 2c-2, row 17): X ≥ 1 and ε_r ≤ 0 decided on the exact P = Nα, ε_r =
    (3ε₀ + 2P)/(3ε₀ − P) rounded once, and the catastrophe's X rounded once."""
    p, t = Fraction(n) * Fraction(a), 3 * Fraction(EPS0)
    if t - p <= 0:
        return ip_err(f"Nα/(3ε₀) = {rust_display(rn_frac(p / t))} ≥ 1: Clausius-Mossotti diverges (polarization "
                      "catastrophe), no physical ε_r")
    eps = rn_frac((t + 2 * p) / (t - p))
    if eps <= 0.0:
        return [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(eps)} (must be > 0)")]
    return okf(eps)


def material_port(fn, v, EPS0, MU0):
    """The port's outputs for a material row an ADR 0002 row changes, else None."""
    if fn == "material::classify_magnetic" and not math.isnan(v[0]):
        x = v[0]
        dia = Fraction(x) < -Fraction(1, 10 ** 15) if math.isfinite(x) else x < 0.0   # ±∞: Rust's test is exact
        exact = "Diamagnetic" if dia else "Ferromagnetic" if x > 1.0 else "Paramagnetic"
        return 17, [("e", "MagneticType." + exact)]
    if fn == "material::magnetic_energy_density_material" and _fin(v[0]):
        rust_below = abs(v[0] * MU0) < 1e-30
        exact_below = below30(Fraction(v[0]) * Fraction(MU0))
        if rust_below != exact_below:
            if exact_below:
                return 17, dz_err("permeability cannot be zero for energy density")
            m13 = material_port13(fn, v, EPS0, MU0)
            return 17, m13 if m13 is not None else okf(fdiv(v[1] * v[1], 2.0 * (v[0] * MU0)))
    if fn == "material::clausius_mossotti" and _fin(v[0], v[1]) and v[1] >= 0.0:
        return 17, material_cm_port(v[0], v[1], EPS0)
    if fn == "material::curie_law" and math.isnan(v[1]):
        return 19, ip_err(f"temperature must be positive, got {rust_display(v[1])} K")
    if fn == "material::curie_weiss" and math.isnan(v[1] - v[2]):
        return 19, ip_err(f"temperature ({rust_display(v[1])} K) must be above Curie temperature "
                          f"({rust_display(v[2])} K)")
    if fn == "Material::wave_speed":
        e, m = v[0], v[1]
        if not (math.isfinite(e) and e > 0.0):
            return 20, [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(e)} (must be > 0)")]
        if not (math.isfinite(m) and m > 0.0):
            return 20, [("E", "InvalidPermeability"), ("s", f"invalid permeability: {rust_display(m)} (must be > 0)")]
    if fn == "material::classify_magnetic" and math.isnan(v[0]):
        return 21, ip_err("magnetic susceptibility must not be NaN, got NaN")
    if fn in ("material::bound_surface_charge", "material::bound_surface_current"):
        # P·n̂ and M × n̂ are FieldVector::dot and ::cross, whose range edges row 18 already fixes (M1 part 5d-2).
        g = "FieldVector::dot" if fn.endswith("charge") else "FieldVector::cross"
        if all(math.isfinite(x) for x in v[0:6]) and fv_governed(g, v[0:6]) == 18:
            return 18, fv_model(g, v[0:6])
        return None
    if fn == "Material::try_reflectance_at":
        # Row 18 where the range edge's Γ (M2 part 2d-2c) moves the row, else row 13 (η rounded once).
        old13 = mat_reflectance(v, mat_lossless_port)
        new = mat_reflectance(v, mat_lossless_port, ux=True)
        return (18 if not toks_same(new, old13) else 13), new
    m13 = material_port13(fn, v, EPS0, MU0)
    return None if m13 is None else (13, m13)


def material_adr_rows(root: Path, path: Path) -> list:
    c = golden_consts(root)
    EPS0, MU0 = c["field::EPSILON_0"], c["field::MU_0"]
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        toks = [] if ins == "-" else [parse_tok(t) for t in ins.split(" ")]
        r = material_port(fn, [b for _k, b in toks], EPS0, MU0)
        if r is None:
            continue
        row, model = r
        old = [(k, tok_text(b) if k == "s" else b) for k, b in (parse_tok(t) for t in outs.split(" "))]
        if not toks_same(model, old):
            out.append((fn, ins, wire(model), str(row), f"oracle {outs}; ADR 0002 row {row}"))
    return out

# ── wave (M2 part 1b): rust-old/src/wave.rs ───────────────────────────────────────────────────────────────────────
FRAC_PI_2 = PI / 2.0          # core::f64::consts::FRAC_PI_2, the same double (a halving is exact)
COS_TOL = 1e-6                # wave.rs's COS_THETA_TOLERANCE
WAVE_C_FNS = ("wave::plane_wave_e", "wave::snell_refraction_angle", "wave::critical_angle", "wave::brewster_angle",
              "wave::fresnel_rs", "wave::fresnel_rp", "wave::reflectance_s", "wave::reflectance_p",
              "wave::hertzian_dipole_pattern", "wave::half_wave_dipole_pattern")


def wave_consts(root: Path) -> dict:
    c = golden_consts(root)
    return {"EPS0": c["field::EPSILON_0"], "MU0": c["field::MU_0"], "C": c["field::SPEED_OF_LIGHT"]}


def okf(x):
    return [("f", x)]


def wave_lossless(e, m, s):
    """Material::lossless_optical_constants: (error tokens, None) or (None, (n, η))."""
    if not (math.isfinite(e) and e > 0.0):
        return [("E", "InvalidPermittivity"), ("s", f"invalid permittivity: {rust_display(e)} (must be > 0)")], None
    if not (math.isfinite(m) and m > 0.0):
        return [("E", "InvalidPermeability"), ("s", f"invalid permeability: {rust_display(m)} (must be > 0)")], None
    if s != 0.0:
        return ip_err(f"Fresnel coefficients from Material are lossless-only; conductivity σ = {rust_display(s)} "
                      "S/m needs a frequency-dependent complex refractive index"), None
    return None, (fsqrt(e * m), fsqrt(fdiv(m, e)))


def wave_scaled_indices(v):
    """src/wave.cyr's indices at an edge (M2 part 2d-2c): n = √(ε_r μ_r) to 53 bits with an unbounded exponent, both
    scaled by the larger's binade into [1, 2), the smaller clamped at 2^-200 below it; and whether they are equal."""
    F = Fraction
    a1, a2 = rnd53_sqrt(F(v[0]) * F(v[1])), rnd53_sqrt(F(v[3]) * F(v[4]))
    b1, b2 = binade(a1), binade(a2)
    top = max(b1, b2)
    n1 = float(a1 / F(2) ** b1 * F(2) ** max(b1 - top, -200))
    n2 = float(a2 / F(2) ** b2 * F(2) ** max(b2 - top, -200))
    return n1, n2, a1 == a2


def wave_gamma_port(eta1, eta2, m1, m2):
    """src/crmath_vx.cyr's _bj_vx_gamma (M2 part 2d-2c, row 18's rule): Rust's Γ where both η's lie in (2^-1022,
    DBL_MAX] and their sum is finite; else on η = √(μ_r/ε_r) to 53 bits (m1, m2: (ε_r, μ_r)), the difference and
    sum to 53 bits, the quotient rounded once."""
    den = eta2 + eta1
    if VX_MN < eta1 <= DBL_MAX and VX_MN < eta2 <= DBL_MAX and den <= DBL_MAX:
        return fdiv(eta2 - eta1, den)
    F = Fraction
    h1, h2 = rnd53_sqrt(F(m1[1]) / F(m1[0])), rnd53_sqrt(F(m2[1]) / F(m2[0]))
    num = rnd53(h2 - h1)
    return 0.0 if num == 0 else rn_frac(num / rnd53(h2 + h1))


def wave_interface(v, lossless, exact=False, ux=False):
    """reflectance_at_interface on a lossless-constants function (Rust's, or the port's rounded once), the cosines
    Rust's or (exact) the port's; ux: the port's range edges (M2 part 2d-2c, row 18's rule): indices at or below
    2^-1022 as wave_scaled_indices, Γ as wave_gamma_port, the ratios on 53-bit η's where an η is outside (2^-1022,
    DBL_MAX]."""
    e1, l1 = lossless(*v[0:3])
    if e1:
        return e1
    e2, l2 = lossless(*v[3:6])
    if e2:
        return e2
    (n1, eta1), (n2, eta2) = l1, l2
    same = n1 == n2
    if ux and not (n1 > VX_MN and n2 > VX_MN):
        n1, n2, same = wave_scaled_indices(v)
    e, ci, ct = wave_cosines(n1, n2, v[6], exact)
    if e:
        return e
    if ct is None:
        return okf(1.0)
    if same:
        g = wave_gamma_port(eta1, eta2, v[0:2], v[3:5]) if ux else fdiv(eta2 - eta1, eta2 + eta1)
        return okf(g * g)
    if ux and not (VX_MN < eta1 <= DBL_MAX and VX_MN < eta2 <= DBL_MAX):
        F = Fraction
        h1, h2 = rnd53_sqrt(F(v[1]) / F(v[0])), rnd53_sqrt(F(v[4]) / F(v[3]))
        e, rs = wave_ratio_ux(h2, ci, h1, ct, False)
        if e:
            return e
        e, rp = wave_ratio_ux(h1, ci, h2, ct, False)
        if e:
            return e
        return okf(0.5 * (rs * rs + rp * rp))
    return wave_unpolarized(eta2, eta1, ci, ct, ux)


def wave_cosines(n1, n2, c, exact=False):
    """fresnel_cosines: (error tokens, None, None) or (None, cos θ_i, cos θ_t or None under TIR); with exact, the
    port's (M2 part 2c-3, row 17): TIR when n₂² − n₁²(1 − c²) < 0 exactly, cos θ_t = √(that)/n₂ rounded once."""
    if n2 == 0.0:
        return dz_err("refractive index n₂ cannot be zero"), None, None
    if not (math.isfinite(n1) and n1 > 0.0 and math.isfinite(n2) and n2 > 0.0):
        return ip_err(f"refractive indices must be positive and finite, got n₁={rust_display(n1)}, "
                      f"n₂={rust_display(n2)}"), None, None
    if 0.0 <= c <= 1.0:
        ci = c
    elif -COS_TOL <= c <= 1.0 + COS_TOL:
        ci = 0.0 if c < 0.0 else 1.0
    else:
        return ip_err("cos θ_i must lie in [0, 1] (measured against the incident-side normal), got "
                      + rust_display(c)), None, None
    if exact:
        d = Fraction(n2) ** 2 - Fraction(n1) ** 2 * (1 - Fraction(ci) ** 2)
        if d < 0:
            return None, ci, None
        return None, ci, 0.0 if d == 0 else rn_sqrt(d / Fraction(n2) ** 2)
    ratio = n1 / n2
    sin2_t = ratio * ratio * (1.0 - ci * ci)
    if sin2_t > 1.0:
        return None, ci, None
    return None, ci, fsqrt(1.0 - sin2_t)


def wave_ratio(num, den):
    if den == 0.0:
        return dz_err("Fresnel denominator is zero"), None
    return None, fdiv(num, den)


VX_MN = 2.0 ** -1022
DBL_MAX = sys.float_info.max


def vx_exact_p(a, b, p) -> bool:
    """src/crmath_vx.cyr's _bj_vx_exact: a zero factor, or |p| in (2^-1022, DBL_MAX]."""
    return a == 0.0 or b == 0.0 or VX_MN < abs(p) <= DBL_MAX


def wave_ratio_of(a, x, b, y, strict, ux):
    """(error tokens, None) or (None, r) for (a·x − b·y)/(a·x + b·y): Rust's chain, or with ux (ADR 0002 row 18's rule,
    M2 part 2d-2a) the unbounded one where Rust's is not (a product not exact for the chain, a sum past DBL_MAX): each
    step to 53 bits (exact.py's rnd53), the quotient rounded once. The denominator test: |den| < 1e-30 (strict, the angle
    forms; on the chain's den) or den = 0."""
    p, q = a * x, b * y
    num, den = p - q, p + q
    plain = (not _fin(a, x, b, y) or (vx_exact_p(a, x, p) and vx_exact_p(b, y, q)
                                      and abs(num) <= DBL_MAX and abs(den) <= DBL_MAX))
    if not ux or plain:
        if (strict and abs(den) < 1e-30) or den == 0.0:
            return dz_err("Fresnel denominator is zero"), None
        return None, fdiv(num, den)
    return wave_ratio_ux(Fraction(a), x, Fraction(b), y, strict)


def wave_ratio_ux(a, x, b, y, strict):
    """wave_ratio_of's unbounded chain for Fractions a and b (53-bit values) and doubles x and y."""
    P, Q = rnd53(a * Fraction(x)), rnd53(b * Fraction(y))
    N, D = rnd53(P - Q), rnd53(P + Q)
    if D == 0 or (strict and abs(D) < Fraction(1e-30)):
        return dz_err("Fresnel denominator is zero"), None
    if N == 0:
        return None, (0.0 if D > 0 else -0.0)
    return None, rn_frac(N / D)


def wave_flux_part(a, b, c, d, ci, k, k2, mode):
    """src/wave.cyr's _bj_wave_flux_part (ADR 0002 row 18's rule, M2 part 2d-2b): one part of S = cross·k (mode 0) or
    g = S·k2 (mode 1) from the port's cross component ci: Rust's chain where ci is normal and (for g) S finite, where
    ci is an exact 0 or past DBL_MAX for S, or an operand is not finite; else the chain from the products, each step
    to 53 bits (rnd53), the last product rounded once."""
    sv = ci * k
    r = sv if mode == 0 else sv * k2
    ac = abs(ci)
    if VX_MN <= ac <= DBL_MAX and (mode == 0 or abs(sv) <= DBL_MAX):
        return r
    if not _fin(a, b, c, d) or (mode == 0 and ac > DBL_MAX):
        return r
    if ci == 0.0 and vx_exact_p(a, b, a * b) and vx_exact_p(c, d, c * d):
        return r
    cc = rnd53(rnd53(Fraction(a) * Fraction(b)) - rnd53(Fraction(c) * Fraction(d)))
    if mode == 0:
        return rn_frac(cc * Fraction(k)) if cc != 0 else 0.0
    sw = rnd53(cc * Fraction(k))
    return rn_frac(sw * Fraction(k2)) if sw != 0 else 0.0


def wave_unpolarized(a1, a2, ci, ct, ux=False):
    """0.5·(r_s² + r_p²) with r_s = (a₁ cos θ_i − a₂ cos θ_t)/(a₁ cos θ_i + a₂ cos θ_t) and r_p the same with a₁, a₂
    swapped: the index form takes (n₁, n₂), the material form (η₂, η₁); ux as wave_ratio_of."""
    e, rs = wave_ratio_of(a1, ci, a2, ct, False, ux)
    if e:
        return e
    e, rp = wave_ratio_of(a2, ci, a1, ct, False, ux)
    if e:
        return e
    return okf(0.5 * (rs * rs + rp * rp))


def wave_direct(fn, n1, n2, c, exact=False, ux=False):
    """The trig-free family, transcribed (exact: on the port's cosines; ux: its ratios as wave_ratio_of)."""
    e, ci, ct = wave_cosines(n1, n2, c, exact)
    if e:
        return e
    if fn == "wave::snell_cos_theta_t":
        return [("n", "")] if ct is None else okf(ct)
    if fn in ("wave::fresnel_rs_direct", "wave::fresnel_rp_direct", "wave::reflectance_s_direct",
              "wave::reflectance_p_direct"):
        if ct is None:
            return [("n", "")] if "fresnel" in fn else okf(1.0)
        if n1 == n2:
            r = 0.0
        else:
            a, b = (n1, n2) if "_s_" in fn or "rs" in fn else (n2, n1)
            err, r = wave_ratio_of(a, ci, b, ct, False, ux)
            if err:
                return err
        return okf(r) if "fresnel" in fn else okf(r * r)
    if fn in ("wave::reflectance_unpolarized", "wave::transmittance_unpolarized"):
        if ct is None:
            out = okf(1.0)
        elif n1 == n2:
            out = okf(0.0)
        else:
            out = wave_unpolarized(n1, n2, ci, ct, ux)
        if fn.startswith("wave::trans") and out[0][0] == "f":
            return okf(1.0 - out[0][1])
        return out
    return None


def wave_reflectance_normal(n1, n2, exact=False, ux=False):
    den = n1 + n2
    if (below30(Fraction(n1) + Fraction(n2)) if exact and _fin(n1, n2) else abs(den) < 1e-30):
        return dz_err("sum of refractive indices cannot be zero")
    _e, r = wave_ratio_of(n1, 1.0, n2, 1.0, False, ux)
    return okf(r * r)


def wave_schlick(n1, n2, c, port=False, exact=False, ux=False):
    """schlick_reflectance; with port, R₀ without reflectance_normal's |n₁ + n₂| < 1e-30 guard (ADR 0002 row 24);
    with exact, on the port's cosines (row 17)."""
    e, ci, ct = wave_cosines(n1, n2, c, exact)
    if e:
        return e
    if n1 > n2:
        if ct is None:
            return okf(1.0)
        cos = ct
    else:
        cos = ci
    if port:
        _e, r = wave_ratio_of(n1, 1.0, n2, 1.0, False, ux)
        r0 = r * r
    else:
        r0t = wave_reflectance_normal(n1, n2)
        if r0t[0][0] != "f":
            return r0t
        r0 = r0t[0][1]
    omc = 1.0 - cos
    omc2 = omc * omc
    omc5 = omc2 * omc2 * omc
    return okf(r0 + (1.0 - r0) * omc5)


def wave_rect_check(a, b, m, n):
    if not (math.isfinite(a) and a > 0.0 and math.isfinite(b) and b > 0.0):
        return ip_err(f"waveguide dimensions must be positive and finite, got a={rust_display(a)}, b={rust_display(b)}")
    if m == 0 and n == 0:
        return ip_err("both m and n cannot be zero (no TE₀₀/TM₀₀ mode)")
    return None


def wave_positive(what, x):
    if not (math.isfinite(x) and x > 0.0):
        return ip_err(f"{what} must be positive and finite, got {rust_display(x)}")
    return None


def wave_rust(fn, v, P, K, ux=False):
    """wave.rs transcribed in its order of operations, its checks and error texts, on the primitives P (glibc's for
    the oracle, MPFR's for the corrected oracle); None for another name."""
    EPS0, MU0, C = K["EPS0"], K["MU0"], K["C"]
    if fn in ("wave::poynting_vector", "wave::momentum_density"):
        a, b = v[0:3], v[3:6]
        cr = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        s = [x * fdiv(1.0, MU0) for x in cr]
        if fn == "wave::momentum_density":
            s = [x * fdiv(1.0, C * C) for x in s]
        return ok3(s)
    if fn == "wave::plane_wave_intensity":
        return okf(0.5 * C * EPS0 * v[0] * v[0])
    if fn == "wave::radiation_pressure_absorbed":
        return okf(fdiv(v[0], C))
    if fn == "wave::radiation_pressure_reflected":
        return okf(fdiv(2.0 * v[0], C))
    if fn == "wave::e_from_b":
        return okf(C * v[0])
    if fn == "wave::b_from_e":
        return okf(fdiv(v[0], C))
    if fn == "wave::plane_wave_e":
        e0, k, x, om, tt, ph = v
        kx, wt = k * x, om * tt
        d = kx - wt
        arg = d + ph
        plain = ((VX_MN < abs(kx) and VX_MN < abs(wt) and abs(d) <= DBL_MAX and abs(arg) <= DBL_MAX)
                 or not _fin(k, x, om, tt, ph, e0)
                 or (vx_exact_p(k, x, kx) and vx_exact_p(om, tt, wt) and abs(d) <= DBL_MAX and abs(arg) <= DBL_MAX))
        if not ux or plain:
            return okf(e0 * P["sin"](arg))
        # Row 18's rule (M2 part 2d-2b): the argument to 53 bits with an unbounded exponent, its sine to 53 bits
        # (MPFR's widest range past DBL_MAX), the product by E₀ rounded once.
        A = rnd53(rnd53(rnd53(Fraction(k) * Fraction(x)) - rnd53(Fraction(om) * Fraction(tt))) + Fraction(ph))
        if A == 0:
            return okf(e0 * 0.0)
        if abs(A) < Fraction(VX_MN):
            sn = A
        else:
            sn = mpfr_ref.wide("sin", A)
        return okf(rn_frac(Fraction(e0) * sn))
    if fn == "wave::wave_number":
        if abs(v[1]) < 1e-30:
            return dz_err("velocity cannot be zero for wave number")
        return okf(fdiv(2.0 * PI * v[0], v[1]))
    if fn == "wave::angular_frequency":
        return okf(2.0 * PI * v[0])
    if fn == "wave::snell_refraction_angle":
        n1, n2, th = v
        if abs(n2) < 1e-30:
            return dz_err("refractive index n₂ cannot be zero")
        sn = P["sin"](th)
        pp = n1 * sn
        st = fdiv(pp, n2)
        if ux and not abs(pp) > VX_MN and _fin(n1, sn, n2) and not vx_exact_p(n1, sn, pp):
            st = rn_frac(rnd53(Fraction(n1) * Fraction(sn)) / Fraction(n2))   # row 18's rule (M2 part 2d-2b)
        return [("n", "")] if abs(st) > 1.0 else okf(P["asin"](st))
    if fn == "wave::critical_angle":
        n1, n2 = v
        if abs(n1) < 1e-30:
            return dz_err("refractive index n₁ cannot be zero")
        r = fdiv(n2, n1)
        return [("n", "")] if r >= 1.0 else okf(P["asin"](r))
    if fn == "wave::brewster_angle":
        n1, n2 = v
        if abs(n1) < 1e-30:
            return dz_err("refractive index n₁ cannot be zero")
        return okf(P["atan"](fdiv(n2, n1)))
    if fn in ("wave::fresnel_rs", "wave::fresnel_rp", "wave::reflectance_s", "wave::reflectance_p"):
        n1, ti, n2, tt = v
        ci, ct = P["cos"](ti), P["cos"](tt)
        a, b = (n1, n2) if fn.endswith("s") else (n2, n1)
        e, r = wave_ratio_of(a, ci, b, ct, True, ux)
        if e:
            return e
        return okf(r) if "fresnel" in fn else okf(r * r)
    if fn == "wave::reflectance_normal":
        return wave_reflectance_normal(v[0], v[1])
    if fn == "wave::transmittance_normal":
        r = wave_reflectance_normal(v[0], v[1])
        return okf(1.0 - r[0][1]) if r[0][0] == "f" else r
    if fn == "wave::schlick_reflectance":
        return wave_schlick(v[0], v[1], v[2])
    direct = wave_direct(fn, *v) if fn.endswith(("_direct", "_unpolarized", "cos_theta_t")) else None
    if direct is not None:
        return direct
    if fn == "wave::reflectance_at_interface":
        return wave_interface(v, wave_lossless)
    if fn == "wave::reflectance_normal_materials":
        return material_rust("Material::try_reflectance_at", v, EPS0, MU0)
    if fn == "wave::material_reflection_coefficient":
        for part in (v[0:3], v[3:6]):
            e, _l = wave_lossless(*part)
            if e:
                return e
        eta1 = material_rust("Material::impedance", v[0:3], EPS0, MU0)[0][1]
        eta2 = material_rust("Material::impedance", v[3:6], EPS0, MU0)[0][1]
        den = eta2 + eta1
        if abs(den) < 1e-30:
            return dz_err("η₁ + η₂ = 0")
        return okf(fdiv(eta2 - eta1, den))
    if fn in ("wave::rectangular_waveguide_cutoff", "wave::rectangular_waveguide_cutoff_wavelength"):
        a, b, m, n = v[0:4]
        e = wave_rect_check(a, b, m, n)
        if e:
            return e
        if fn.endswith("cutoff"):
            e = wave_positive("wave velocity", v[4])
            if e:
                return e
        ma, nb = fdiv(float(m), a), fdiv(float(n), b)
        m2, n2 = ma * ma, nb * nb
        root = fsqrt(m2 + n2)
        # Row 18's rule (M2 part 2d-2c): where a square is at or below 2^-1022 (and not an exact 0) or their sum past
        # DBL_MAX, or v/2 at or below 2^-1022, the root runs to 53 bits with an unbounded exponent, f_c and λ_c
        # rounded once.
        plain = m2 + n2 <= DBL_MAX and (m == 0 or m2 > VX_MN) and (n == 0 or n2 > VX_MN)
        F = Fraction
        R = None
        if ux and not plain:
            q, p = rnd53(F(m) / F(a)), rnd53(F(n) / F(b))
            R = rnd53_sqrt(rnd53(rnd53(q * q) + rnd53(p * p)))
        if fn.endswith("cutoff"):
            hv = 0.5 * v[4]
            if ux and not (plain and hv > VX_MN):
                if R is None:
                    q, p = rnd53(F(m) / F(a)), rnd53(F(n) / F(b))
                    R = rnd53_sqrt(rnd53(rnd53(q * q) + rnd53(p * p)))
                return okf(rn_frac(F(v[4]) / 2 * R))
            return okf(hv * root)
        if R is not None:
            lam = rn_frac(2 / R)
            if lam <= DBL_MAX:
                return okf(lam)
            root = 0.0
        if math.isnan(root) or root <= 0.0:
            return ip_err(f"cutoff wavelength not representable for a={rust_display(a)}, b={rust_display(b)}")
        return okf(fdiv(2.0, root))
    if fn == "wave::guide_wavelength":
        lam, lc = v
        if lc == 0.0:
            return dz_err("cutoff wavelength cannot be zero")
        if math.isnan(lc) or lc < 0.0:
            return ip_err(f"cutoff wavelength must be positive, got {rust_display(lc)}")
        if not (math.isfinite(lam) and lam > 0.0):
            return ip_err(f"wavelength must be positive and finite, got {rust_display(lam)}")
        ratio = fdiv(lam, lc)
        rad = 1.0 - ratio * ratio
        if math.isnan(rad) or rad <= 0.0:
            return ip_err(f"frequency below cutoff: λ/λ_c = {rust_fixed(ratio, 4)} (must be < 1)")
        return okf(fdiv(lam, fsqrt(rad)))
    if fn == "wave::cylindrical_waveguide_cutoff":
        r, x, vel = v
        e = (wave_positive("waveguide radius", r) or wave_positive("Bessel root", x)
             or wave_positive("wave velocity", vel))
        return e or okf(fdiv(x * vel, 2.0 * PI * r))
    if fn == "wave::hertzian_dipole_pattern":
        s = P["sin"](v[0])
        return okf(s * s)
    if fn == "wave::half_wave_dipole_pattern":
        s = P["sin"](v[0])
        if abs(s) < 1e-15:
            return okf(0.0)
        f = fdiv(P["cos"](FRAC_PI_2 * P["cos"](v[0])), s)
        return okf(f * f)
    if fn == "wave::hertzian_dipole_directivity":
        return okf(1.5)
    if fn == "wave::half_wave_dipole_directivity":
        return okf(1.6409)
    if fn == "wave::half_wave_dipole_radiation_resistance":
        return okf(73.1)
    if fn == "wave::hertzian_dipole_radiation_resistance":
        if abs(v[1]) < 1e-30:
            return dz_err("wavelength cannot be zero")
        ratio = fdiv(v[0], v[1])
        return okf(80.0 * PI * PI * ratio * ratio)
    if fn == "wave::effective_aperture":
        lam, d = v
        if lam < 0.0:
            return ip_err(f"wavelength must be non-negative, got {rust_display(lam)}")
        return okf(fdiv(lam * lam * d, 4.0 * PI))
    if fn == "wave::friis_transmission":
        gt, gr, lam, d = v
        if abs(d) < 1e-30:
            return SINGULAR
        ratio = fdiv(lam, 4.0 * PI * d)
        return okf(gt * gr * ratio * ratio)
    return None


def wave_toks(outs):
    return [(k, tok_text(b) if k == "s" else b) for k, b in outs]


def validate_wave(root: Path, path: Path) -> list:
    header, rows = read_tsv(path)
    errs = check_platform(header)
    if errs:
        return errs
    K = wave_consts(root)
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [b for _k, b in ins]
        want = wave_rust(fn, v, GLIBC_P, K)
        if want is None:
            errs.append(f"wave row {n}: no transcription for {fn}")
            continue
        if not toks_same(want, wave_toks(outs)):
            errs.append(f"wave row {n}: {fn}{tuple(v)} recorded {wave_toks(outs)}, transcription {want}")
    return errs


def wave_corrections(root: Path, path: Path):
    """The class-C wave rows whose outputs change when sin, cos, asin and atan are correctly rounded (MPFR)."""
    K = wave_consts(root)

    def fixes(_rows):
        out = []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                continue
            fn, cls, ins, outs, _note = line.split("\t")
            if cls != "C" or fn == "wave::half_wave_dipole_pattern":
                continue                    # (the pattern's rows are ADR 0002 row 28's mappings)
            v = [b for _k, b in (parse_tok(t) for t in ins.split(" "))]
            want = wave_rust(fn, v, MPFR_P, K)
            if not toks_same(want, wave_rust(fn, v, MPFR_P, K, ux=True)):
                continue                    # row 18's (M2 part 2d-2a): the range edge's chain moves it
            old = wave_toks([parse_tok(t) for t in outs.split(" ")])
            if toks_same(want, old):
                continue
            out.append((fn, ins, wire(want), f"oracle {outs}; correctly rounded primitives give {wire(want)}"))
        return out
    return fixes


# ADR 0002 rows 22–25 (M2 part 1b, the user's 2026-11-12 answers), and row 18 through FieldVector::cross. Row 22: the
# dipoles on bijli's η₀ = √(μ₀/ε₀) (maxwell's free_space_impedance) with Cin(2π) exact, where Rust has 120π's textbook
# figures: D = 4/Cin(2π), R = η₀·Cin(2π)/(4π), and the Hertzian (2π/3)·η₀·(dl/λ)² for 80π²(dl/λ)², each constant
# rounded once (src/wave_data.cyr). Row 23: effective_aperture refuses a NaN wavelength (`!(λ ≥ 0)`), as row 19.
# Row 24: Schlick's R₀ without reflectance_normal's |n₁ + n₂| < 1e-30 guard. Row 25: material_reflection_coefficient's
# Γ on the relative impedances √(μ_r/ε_r) its validation computes, so Γ² is reflectance_normal_materials' R.
def cin_2pi_frac(bits=700):
    """Cin(2π) = Σ_{k ≥ 1} (−1)^(k+1) (2π)^(2k) / (2k·(2k)!) as a Fraction, with π = crtables.pi_bits(bits)/2^bits
    (MPFR's const_bits must agree); its error is below 2^-(bits − 8)."""
    pb = crtables.pi_bits(bits)
    if pb != mpfr_ref.const_bits("pi", bits):
        raise SystemExit("π: crtables and MPFR disagree")
    x = 2 * Fraction(pb, 1 << bits)
    x2, term, total = x * x, Fraction(1), Fraction(0)
    for k in range(1, 120):
        term = term * x2 / ((2 * k - 1) * (2 * k))
        total += (term if k % 2 else -term) / (2 * k)
    return total, Fraction(pb, 1 << bits)


def frac_sqrt(q: Fraction, bits=700) -> Fraction:
    return Fraction(math.isqrt(q.numerator * (1 << (2 * bits)) // q.denominator), 1 << bits)


def frac_rn(q: Fraction, rel=Fraction(1, 1 << 600)) -> float:
    """q rounded to the nearest double, refused if q's error bound (relative `rel`) straddles a rounding boundary."""
    lo, hi = float(q * (1 - rel)), float(q * (1 + rel))
    if lo != hi:
        raise SystemExit(f"dipole constant {float(q)!r}: not decided at this precision")
    return lo


def dipole_consts(root: Path):
    """(D, R, K) of row 22, each the double nearest the exact value."""
    c = golden_consts(root)
    cin, pi = cin_2pi_frac()
    eta = frac_sqrt(Fraction(c["field::MU_0"]) / Fraction(c["field::EPSILON_0"]))
    return frac_rn(4 / cin), frac_rn(eta * cin / (4 * pi)), frac_rn(2 * pi / 3 * eta)


def render_wave_data(root: Path) -> str:
    d, r, k = dipole_consts(root)
    return "\n".join([
        "# wave_data.cyr — the dipole constants of src/wave.cyr (ADR 0002 row 22). Generated by scripts/oracle/validate.py",
        "# (python3 -I -B scripts/oracle/validate.py . --write) from π (Machin, checked against MPFR), Cin(2π) by its",
        "# series and η₀ = √(μ₀/ε₀) of BJ_MU_0 and BJ_EPSILON_0, all in exact rational arithmetic; do not edit.",
        "",
        "# The doubles, as bit patterns in typed f64 globals.",
        f"var _BJ_WAVE_HALF_WAVE_D: f64 = 0x{f_to_bits(d):016X};  # 4/Cin(2π) = {d!r}, nearest",
        f"var _BJ_WAVE_HALF_WAVE_R: f64 = 0x{f_to_bits(r):016X};  # η₀·Cin(2π)/(4π) = {r!r} Ω, nearest",
        f"var _BJ_WAVE_HERTZ_K: f64 = 0x{f_to_bits(k):016X};  # (2π/3)·η₀ = {k!r} Ω, nearest",
        "",
    ])


def check_wave_data(root: Path) -> list:
    path = root / "src/wave_data.cyr"
    text = render_wave_data(root)
    if WRITE:
        if not path.exists() or path.read_text() != text:
            path.write_text(text)
            print("wrote src/wave_data.cyr")
        return []
    if not path.exists() or path.read_text() != text:
        return ["src/wave_data.cyr: differs from its rendering (regenerate with --write)"]
    print("ok   src/wave_data.cyr: the dipole constants")
    return []


def wave_half_wave_port(theta):
    """ADR 0002 row 28 (M2's range-edge Part, the user's answer of 2026-11-12): the pattern's numerator
    cos((π/2)cos θ) as sin(π·sin²(θ/2)), or sin(π·cos²(θ/2)) where cos θ < 0, on correctly rounded sin and cos."""
    P = MPFR_P
    s = P["sin"](theta)
    if abs(s) < 1e-15:
        return okf(0.0)
    c = P["cos"](theta)
    h = 0.5 * theta
    if c >= 0.0:
        sh = P["sin"](h)
        m = sh * sh
    else:
        ch = P["cos"](h)
        m = ch * ch
    f = fdiv(P["sin"](PI * m), s)
    return okf(f * f)


# Row 13 (M2 part 2b-2, the user's "round once" decision of 2026-11-12): wave.rs's products and quotients of inputs and
# constants, correctly rounded where the operands are finite (exact values); a constant enters as its double, 2π and
# 4π as Rust's π times 2 and 4 (exact), the Hertzian factor as row 22's double; every error keeps Rust's result.
def wave_port13(fn, v, K, kh):
    """Row 13's port outputs for a wave row, or None where Rust's result stands."""
    C, EPS0 = K["C"], K["EPS0"]
    if fn == "wave::plane_wave_intensity":
        return okf(rn_mono([0.5, C, EPS0, v[0], v[0]], [])) if _fin(v[0]) else None
    if fn == "wave::radiation_pressure_reflected":
        return okf(rn_mono([2.0, v[0]], [C])) if _fin(v[0]) else None
    if fn == "wave::wave_number":
        if abs(v[1]) < 1e-30 or not _fin(v[0], v[1]):
            return None
        return okf(rn_mono([2.0 * PI, v[0]], [v[1]]))
    if fn == "wave::cylindrical_waveguide_cutoff":
        r, x, vel = v
        if (wave_positive("waveguide radius", r) or wave_positive("Bessel root", x)
                or wave_positive("wave velocity", vel)):
            return None
        return okf(rn_mono([x, vel], [2.0 * PI, r]))
    if fn == "wave::hertzian_dipole_radiation_resistance":
        if abs(v[1]) < 1e-30 or not _fin(v[0], v[1]):
            return None
        return okf(rn_mono([kh, v[0], v[0]], [v[1], v[1]]))
    if fn == "wave::effective_aperture":
        lam, d = v
        if not lam >= 0.0 or not _fin(lam, d):
            return None
        return okf(rn_mono([lam, lam, d], [4.0 * PI]))
    if fn == "wave::friis_transmission":
        gt, gr, lam, d = v
        if abs(d) < 1e-30 or not _fin(gt, gr, lam, d):
            return None
        return okf(rn_mono([gt, gr, lam, lam], [4.0 * PI, 4.0 * PI, d, d]))
    return None


WAVE_TRIG_FREE = ("wave::snell_cos_theta_t", "wave::fresnel_rs_direct", "wave::fresnel_rp_direct",
                  "wave::reflectance_s_direct", "wave::reflectance_p_direct", "wave::reflectance_unpolarized",
                  "wave::transmittance_unpolarized")


def wave_port(fn, v, K, consts):
    """The port's outputs for a wave row an ADR 0002 row changes, else None."""
    d, r, kh = consts
    if fn == "wave::half_wave_dipole_pattern":
        return 28, wave_half_wave_port(v[0])
    if fn in ("wave::poynting_vector", "wave::momentum_density"):
        # Row 18 through FieldVector::cross, and (M2 part 2d-2b) the flux parts' own edges.
        e, b = v[0:3], v[3:6]
        if all(math.isfinite(x) for x in v[0:6]) and fv_governed("FieldVector::cross", v[0:6]) == 18:
            cr = [x for _k, x in fv_model("FieldVector::cross", v[0:6])]
        else:
            cr = [e[1] * b[2] - e[2] * b[1], e[2] * b[0] - e[0] * b[2], e[0] * b[1] - e[1] * b[0]]
        k, k2 = fdiv(1.0, K["MU0"]), fdiv(1.0, K["C"] * K["C"])
        mode = 1 if fn == "wave::momentum_density" else 0
        parts = [(e[1], b[2], e[2], b[1]), (e[2], b[0], e[0], b[2]), (e[0], b[1], e[1], b[0])]
        out = ok3([wave_flux_part(*parts[i], cr[i], k, k2, mode) for i in range(3)])
        return None if toks_same(out, wave_rust(fn, v, None, K)) else (18, out)
    if fn == "wave::half_wave_dipole_directivity":
        return 22, okf(d)
    if fn == "wave::half_wave_dipole_radiation_resistance":
        return 22, okf(r)
    if fn == "wave::hertzian_dipole_radiation_resistance" and not abs(v[1]) < 1e-30:
        ratio = fdiv(v[0], v[1])
        old22 = okf(kh * ratio * ratio)
        new = wave_port13(fn, v, K, kh) or old22
        # Row 22 where bijli's factor alone moves the row, else row 13 (the product rounded once).
        return (22 if not toks_same(old22, wave_rust(fn, v, None, K)) else 13), new
    if fn == "wave::effective_aperture" and math.isnan(v[0]):
        return 23, ip_err("wavelength must be non-negative, got NaN")
    if fn in ("wave::fresnel_rs", "wave::fresnel_rp", "wave::reflectance_s", "wave::reflectance_p",
              "wave::plane_wave_e", "wave::snell_refraction_angle", "wave::rectangular_waveguide_cutoff",
              "wave::rectangular_waveguide_cutoff_wavelength"):
        # Row 18 (M2 parts 2d-2a to 2d-2c) where the range edge's chain moves the row, on correctly rounded primitives.
        new = wave_rust(fn, v, MPFR_P, K, ux=True)
        return (18, new) if not toks_same(new, wave_rust(fn, v, MPFR_P, K)) else None
    if fn == "wave::schlick_reflectance":
        # Row 18 where the range edge's chain alone moves the row, else row 24 where R₀'s guard does, else row 17 (the
        # exact TIR and cos θ_t).
        old24 = wave_schlick(v[0], v[1], v[2], port=True)
        old18 = wave_schlick(v[0], v[1], v[2], port=True, exact=True)
        new = wave_schlick(v[0], v[1], v[2], port=True, exact=True, ux=True)
        if not toks_same(new, old18):
            return 18, new
        return (24 if not toks_same(old24, wave_rust(fn, v, None, K)) else 17), new
    if fn in WAVE_TRIG_FREE:
        old18 = wave_direct(fn, v[0], v[1], v[2], exact=True)
        new = wave_direct(fn, v[0], v[1], v[2], exact=True, ux=True)
        return (18 if not toks_same(new, old18) else 17), new
    if fn == "wave::material_reflection_coefficient":
        for part in (v[0:3], v[3:6]):
            e, _l = wave_lossless(*part)
            if e:
                return None

        def gamma(eta1, eta2):
            den = eta2 + eta1
            if abs(den) < 1e-30:
                return dz_err("η₁ + η₂ = 0")
            return okf(fdiv(eta2 - eta1, den))
        rust = wave_rust(fn, v, None, K)
        old25 = gamma(wave_lossless(*v[0:3])[1][1], wave_lossless(*v[3:6])[1][1])
        p1, p2 = mat_lossless_port(*v[0:3])[1][1], mat_lossless_port(*v[3:6])[1][1]
        guarded = gamma(p1, p2)
        old = okf(fdiv(p2 - p1, p2 + p1))               # row 25 amended (M2 part 2d-1): the guard is gone
        new = okf(wave_gamma_port(p1, p2, v[0:2], v[3:5]))
        if not toks_same(new, old):
            return 18, new                              # the range edge's Γ (M2 part 2d-2c)
        # Row 25 where the relative η or the guard's removal alone moves the row, else row 13 (η rounded once).
        return (25 if not toks_same(old25, rust) or not toks_same(guarded, new) else 13), new
    if fn == "wave::reflectance_at_interface":
        # Row 18 where the range edge's chain alone moves the row, else row 17 where the exact TIR and cos θ_t do,
        # else row 13 (n and η rounded once).
        old17 = wave_interface(v, wave_lossless, exact=True)
        old18 = wave_interface(v, mat_lossless_port, exact=True)
        new = wave_interface(v, mat_lossless_port, exact=True, ux=True)
        if not toks_same(new, old18):
            return 18, new
        return (17 if not toks_same(old17, wave_rust(fn, v, None, K)) else 13), new
    if fn == "wave::reflectance_normal_materials":
        old13 = mat_reflectance(v, mat_lossless_port)
        new = mat_reflectance(v, mat_lossless_port, ux=True)
        return (18 if not toks_same(new, old13) else 13), new
    if fn in ("wave::reflectance_normal", "wave::transmittance_normal"):
        old18 = wave_reflectance_normal(v[0], v[1], exact=True)
        new = wave_reflectance_normal(v[0], v[1], exact=True, ux=True)
        if fn == "wave::transmittance_normal":
            old18 = okf(1.0 - old18[0][1]) if old18[0][0] == "f" else old18
            new = okf(1.0 - new[0][1]) if new[0][0] == "f" else new
        # Row 18 where the range edge's chain moves the row, else row 17 (the exact |n₁ + n₂| test).
        return (18 if not toks_same(new, old18) else 17), new
    m13 = wave_port13(fn, v, K, kh)
    return None if m13 is None else (13, m13)


def wave_adr_rows(root: Path, path: Path) -> list:
    K = wave_consts(root)
    consts = dipole_consts(root)
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        toks = [] if ins == "-" else [parse_tok(t) for t in ins.split(" ")]
        r = wave_port(fn, [b for _k, b in toks], K, consts)
        if r is None:
            continue
        row, model = r
        old = wave_toks([parse_tok(t) for t in outs.split(" ")])
        if not toks_same(model, old):
            out.append((fn, ins, wire(model), str(row), f"oracle {outs}; ADR 0002 row {row}"))
    return out

# ── circuit (M2 part 1c): rust-old/src/circuit.rs ─────────────────────────────────────────────────────────────────
def circuit_lc(l, c):
    if not (math.isfinite(l) and l > 0.0 and math.isfinite(c) and c > 0.0):
        return ip_err(f"L and C must be finite and positive, got L={rust_display(l)}, C={rust_display(c)}")
    return None


def circuit_r(r):
    if not (math.isfinite(r) and r >= 0.0):
        return ip_err(f"resistance must be finite and non-negative, got {rust_display(r)}")
    return None


def circuit_combo(xs, what, parallel, which, ux=False):
    """The series/parallel combinations: a sum (Rust's Sum, a left fold from −0.0) or 1/Σ(1/x) with its two cuts; ux:
    a fold of finite values past DBL_MAX runs with an unbounded exponent, each step to 53 bits (exact.py's rnd53), the
    last rounded once (ADR 0002 row 18's rule, M2 part 2d-3a)."""
    if not xs:
        return ip_err(f"at least one {what} value is required")
    if not parallel:
        acc = -0.0
        for x in xs:
            acc = acc + x
        if ux and not math.isfinite(acc) and all(math.isfinite(x) for x in xs):
            a = Fraction(xs[0])
            for x in xs[1:-1]:
                a = rnd53(a + Fraction(x))
            return okf(rn_frac(a + Fraction(xs[-1])))
        return okf(acc)
    s = 0.0
    for x in xs:
        if abs(x) < 1e-30:
            return dz_err(f"{what} cannot be zero in {which} combination")
        s += fdiv(1.0, x)
    if abs(s) < 1e-30:
        return dz_err(f"sum of inverse {what}s is zero")
    return okf(fdiv(1.0, s))


def circuit_plain(tau, x, w, t, mode) -> bool:
    """src/circuit.cyr's _bj_circuit_plain: whether a transient's Rust chain is the unbounded one (M2 part 2d-3b)."""
    if not abs(tau) <= DBL_MAX:
        return False
    if not VX_MN < abs(x) <= DBL_MAX and t != 0.0:
        return False
    if not abs(w) <= DBL_MAX:
        return False
    return mode == 0 or abs(w) > VX_MN


def circuit_edge(av, asg, tau, t, mode):
    """A transient at a range edge (ADR 0002 row 18's rule, M2 part 2d-3b): asg·av·w for x = −t/τ to 53 bits,
    w = −expm1(x) (mode 0) or e^x (mode 1) at 53 bits in MPFR's widest range (mpfr_ref.wide), the product rounded once.
    av ≥ 0 and tau are 53-bit Fractions, asg the amplitude's sign; a zero amplitude or t gives IEEE's signed zero."""
    xsg = -math.copysign(1.0, t) * (1.0 if tau > 0 else -1.0)
    if t == 0:
        if mode == 1:
            return asg * 0.0 if av == 0 else math.copysign(rn_frac(av), asg)
        return asg * -xsg * 0.0
    if av == 0:
        return asg * (-xsg if mode == 0 else 1.0) * 0.0
    x = rnd53(Fraction(-t) / tau)
    if abs(x) > 100000:
        # e^±10^5 lies past 2^±144000: a charge's w is 1 for x < 0; otherwise the product overflows or rounds to zero
        # whatever the amplitude (between 2^-2098 and 2^1124). MPFR's widest range would hold e^x, but not as a Fraction.
        if mode == 0 and x < 0:
            return rn_frac(Fraction(asg) * av)
        if mode == 1 and x < 0:
            return asg * 0.0
        return asg * (-math.inf if mode == 0 else math.inf)
    w = mpfr_ref.wide("expm1" if mode == 0 else "exp", x)
    if mode == 0:
        w = -w
    if isinstance(w, float):                    # past even MPFR's widest range: ±∞, or e^x = 0
        return asg * w if w != 0 else asg * 0.0
    return rn_frac(Fraction(asg) * av * w)


def circuit_rust(fn, v, P, expm1=False, exact=False, ux=False):
    """circuit.rs transcribed in its order of operations, its checks and error texts, on the primitives P; None for
    another name. ux: the port's range edges (ADR 0002 row 18's rule, M2 part 2d-3a), the series folds and the
    impedance."""
    name = fn.split("::")[1]
    if name in ("resistance_series", "resistance_parallel", "capacitance_series", "capacitance_parallel",
                "inductance_series", "inductance_parallel"):
        what, kind = name.split("_")
        sums = {"resistance": "parallel", "capacitance": "series", "inductance": "parallel"}
        parallel = kind == sums[what]
        return circuit_combo(v, what, parallel, kind, ux)
    if name == "voltage" or name == "power" or name == "capacitor_charge" or name == "rc_time_constant":
        return okf(v[0] * v[1])
    if name == "current":
        return dz_err("resistance cannot be zero for Ohm's law") if abs(v[1]) < 1e-30 else okf(fdiv(v[0], v[1]))
    if name == "resistance":
        if abs(v[1]) < 1e-30:
            return dz_err("current cannot be zero for resistance calculation")
        return okf(fdiv(v[0], v[1]))
    if name == "resistance_from_geometry":
        if v[2] <= 0.0:
            return ip_err(f"cross-section area must be positive, got {rust_display(v[2])}")
        return okf(fdiv(v[0] * v[1], v[2]))
    if name == "parallel_plate_capacitance":
        if v[2] <= 0.0:
            return ip_err(f"plate separation must be positive, got {rust_display(v[2])}")
        return okf(fdiv(v[0] * v[1], v[2]))
    if name in ("capacitor_energy", "inductor_energy"):
        return okf(0.5 * v[0] * v[1] * v[1])
    if name == "solenoid_inductance":
        return okf(v[0] * v[1] * v[1] * v[2] * v[3])
    if name in ("rc_charging_voltage", "rc_discharging_voltage"):
        v0, r, c, tt = v
        tau = r * c
        instant = below30(Fraction(r) * Fraction(c)) if exact and _fin(r, c) else abs(tau) < 1e-30
        if instant:
            return okf(v0 if name.startswith("rc_charging") else 0.0)
        x = fdiv(-tt, tau)
        charge = name.startswith("rc_charging")
        if expm1 and charge:
            w = -mpfr_ref.cr("expm1", x)                # ADR 0002 row 29: −expm1(x), correctly rounded
        else:
            e = P["exp"](x)
            w = 1.0 - e if charge else e
        mode = 0 if charge else 1
        if ux and _fin(v0, r, c, tt) and not circuit_plain(tau, x, w, tt, mode):
            return okf(circuit_edge(abs(Fraction(v0)), math.copysign(1.0, v0), rnd53(Fraction(r) * Fraction(c)),
                                    tt, mode))
        return okf(v0 * w)
    if name == "rl_time_constant":
        if abs(v[1]) < 1e-30:
            return dz_err("resistance cannot be zero for RL time constant")
        return okf(fdiv(v[0], v[1]))
    if name in ("rl_current_rise", "rl_current_decay"):
        a, r, l, tt = v
        if abs(r) < 1e-30:
            return dz_err("resistance cannot be zero for RL current" if name.endswith("rise")
                          else "resistance cannot be zero for RL decay")
        tau = fdiv(l, r)
        instant = below30(Fraction(l) / Fraction(r)) if exact and _fin(l, r) else abs(tau) < 1e-30
        if instant:
            return okf(fdiv(a, r) if name.endswith("rise") else 0.0)
        x = fdiv(-tt, tau)
        rise = name.endswith("rise")
        amp = fdiv(a, r) if rise else a
        if expm1 and rise:
            w = -mpfr_ref.cr("expm1", x)
        else:
            e = P["exp"](x)
            w = 1.0 - e if rise else e
        mode = 0 if rise else 1
        amp_ok = not rise or a == 0.0 or VX_MN < abs(amp) <= DBL_MAX
        if ux and _fin(a, r, l, tt) and not (amp_ok and circuit_plain(tau, x, w, tt, mode)):
            F = Fraction
            av = abs(rnd53(F(a) / F(r))) if rise else abs(F(a))
            asg = math.copysign(1.0, a) * (math.copysign(1.0, r) if rise else 1.0)
            return okf(circuit_edge(av, asg, rnd53(F(l) / F(r)), tt, mode))
        return okf(amp * w)
    if name in ("resonant_frequency", "resonant_angular_frequency"):
        e = circuit_lc(v[0], v[1])
        if e:
            return e
        s = fsqrt(v[0] * v[1])
        return okf(fdiv(1.0, 2.0 * PI * s) if name == "resonant_frequency" else fdiv(1.0, s))
    if name == "quality_factor":
        r, l, c = v
        e = circuit_r(r)
        if e:
            return e
        if r < 1e-30:
            return dz_err("resistance cannot be zero for quality factor")
        e = circuit_lc(l, c)
        return e or okf(fdiv(fsqrt(fdiv(l, c)), r))
    if name == "damping_ratio":
        r, l, c = v
        e = circuit_r(r) or circuit_lc(l, c)
        return e or okf(r * 0.5 * fsqrt(fdiv(c, l)))
    if name == "rlc_impedance":
        r, l, c, w = v
        if abs(w) < 1e-30:
            return dz_err("angular frequency cannot be zero for impedance")
        if abs(c) < 1e-30:
            return dz_err("capacitance cannot be zero for impedance")
        x = w * l - fdiv(1.0, w * c)
        sq = r * r + x * x
        if ux and not 2.0 ** -900 <= sq <= DBL_MAX and _fin(r, l, c, w):
            # The unbounded chain, each step to 53 bits, √ rounded once (src/circuit.cyr's band guard on Rust's sum).
            F = Fraction
            xx = rnd53(rnd53(F(w) * F(l)) - rnd53(1 / rnd53(F(w) * F(c))))
            sv = rnd53(rnd53(F(r) * F(r)) + rnd53(xx * xx))
            return okf(0.0 if sv == 0 else rn_sqrt(sv))
        return okf(fsqrt(sq))
    return None


def circuit_vals(ins):
    """A row's inputs as values; a `v:n` slice becomes its list."""
    if ins and ins[0][0] == "v":
        return [b for _k, b in ins[1:]]
    return [b for _k, b in ins]


def validate_circuit(path: Path) -> list:
    header, rows = read_tsv(path)
    errs = check_platform(header)
    if errs:
        return errs
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = circuit_vals(ins)
        want = circuit_rust(fn, v, GLIBC_P)
        if want is None:
            errs.append(f"circuit row {n}: no transcription for {fn}")
            continue
        if not toks_same(want, wave_toks(outs)):
            errs.append(f"circuit row {n}: {fn}{tuple(v)} recorded {wave_toks(outs)}, transcription {want}")
    return errs


def circuit_corrections(path: Path):
    """The class-C circuit rows (the transients) whose outputs change when exp is correctly rounded (MPFR)."""
    def fixes(_rows):
        out = []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                continue
            fn, cls, ins, outs, _note = line.split("\t")
            if cls != "C" or fn in ("circuit::rc_charging_voltage", "circuit::rl_current_rise"):
                continue                    # (their rows are ADR 0002 row 29's mappings)
            v = circuit_vals([parse_tok(x) for x in ins.split(" ")])
            if fn in CIRCUIT_TRANSIENTS and circuit_flip17(fn, v):
                continue                    # (row 17's mappings: the exact τ test)
            want = circuit_rust(fn, v, MPFR_P)
            old = wave_toks([parse_tok(x) for x in outs.split(" ")])
            if toks_same(want, old):
                continue
            out.append((fn, ins, wire(want), f"oracle {outs}; correctly rounded primitives give {wire(want)}"))
        return out
    return fixes


# ADR 0002 row 26 (M2 part 1c, rows 19 and 23's rule): resistance_from_geometry and parallel_plate_capacitance test
# `A <= 0` and `d <= 0`, so a NaN area or separation passed as Ok(NaN); the port refuses unless it is positive.
# Row 13 (M2 part 2b-3, the user's "round once" decision of 2026-11-12): circuit.rs's products, quotients and roots of
# inputs and constants, correctly rounded where the operands are finite (exact values); 2π is Rust's π doubled; every
# error, and a zero R in the damping ratio, keeps Rust's result.
def circuit_port13(fn, v):
    """Row 13's port outputs for a circuit row, or None where Rust's result stands."""
    F = Fraction
    if fn in ("circuit::resistance_from_geometry", "circuit::parallel_plate_capacitance"):
        a, b, c = v
        if not c > 0.0 or not _fin(a, b, c):
            return None
        return okf(rn_mono([a, b], [c]))
    if fn in ("circuit::capacitor_energy", "circuit::inductor_energy"):
        return okf(rn_mono([0.5, v[0], v[1], v[1]], [])) if _fin(v[0], v[1]) else None
    if fn == "circuit::solenoid_inductance":
        return okf(rn_mono([v[0], v[1], v[1], v[2], v[3]], [])) if _fin(*v) else None
    if fn in ("circuit::resonant_frequency", "circuit::resonant_angular_frequency"):
        l, c = v
        if circuit_lc(l, c):
            return None
        if fn == "circuit::resonant_frequency":
            return okf(rn_sqrt(1 / (F(2.0 * PI) ** 2 * F(l) * F(c))))
        return okf(rn_sqrt(1 / (F(l) * F(c))))
    if fn in ("circuit::quality_factor", "circuit::damping_ratio"):
        r, l, c = v
        if circuit_r(r) or circuit_lc(l, c):
            return None
        if fn == "circuit::quality_factor":
            return None if r < 1e-30 else okf(rn_sqrt(F(l) / (F(c) * F(r) ** 2)))
        return None if r == 0.0 else okf(rn_sqrt(F(r) ** 2 * F(c) / (4 * F(l))))
    return None


CIRCUIT_TRANSIENTS = ("circuit::rc_charging_voltage", "circuit::rc_discharging_voltage", "circuit::rl_current_rise",
                      "circuit::rl_current_decay")


def circuit_flip17(fn, v) -> bool:
    """Whether the exact τ test (row 17's rule) decides a transient row differently from Rust's rounded one."""
    return not toks_same(circuit_rust(fn, v, MPFR_P, expm1=True), circuit_rust(fn, v, MPFR_P, expm1=True, exact=True))


# Row 17 (M2 part 2c-5, the user's "exact tests, accurate values" decision of 2026-11-12): the inverse combinations
# take 1/Σ 1/xᵢ from the exact sum (an ±∞ adds 0) rounded once, and refuse |Σ 1/xᵢ| < 10^-30 exactly; an empty slice,
# a value refused by the element test and a NaN value keep Rust's result (src/crmath_edge.cyr's _bj_cr_parallel).
CIRCUIT_INVERSE = ("circuit::resistance_parallel", "circuit::capacitance_series", "circuit::inductance_parallel")


def circuit_inverse17(fn, xs):
    """The port's result for an inverse combination row, or None where Rust's stands."""
    if not xs or any(abs(x) < 1e-30 for x in xs) or any(math.isnan(x) for x in xs):
        return None
    what = fn.split("::")[1].split("_")[0]
    s = sum(Fraction(1) / Fraction(x) for x in xs if math.isfinite(x))
    if abs(s) < E30:
        return dz_err(f"sum of inverse {what}s is zero")
    return okf(rn_frac(1 / s))


CIRCUIT_SERIES = ("circuit::resistance_series", "circuit::capacitance_parallel", "circuit::inductance_series")


def circuit_port(fn, v):
    if fn in CIRCUIT_SERIES or fn == "circuit::rlc_impedance":
        # Row 18 (M2 part 2d-3a) where the range edge's chain moves the row.
        new = circuit_rust(fn, v, MPFR_P, ux=True)
        return None if toks_same(new, circuit_rust(fn, v, MPFR_P)) else (18, new)
    if fn in CIRCUIT_INVERSE:
        m17 = circuit_inverse17(fn, v)
        return None if m17 is None else (17, m17)
    if fn in CIRCUIT_TRANSIENTS:
        model = circuit_rust(fn, v, MPFR_P, expm1=True, exact=True)
        new = circuit_rust(fn, v, MPFR_P, expm1=True, exact=True, ux=True)
        if not toks_same(new, model):
            return 18, new                              # the range edge's chain (M2 part 2d-3b)
        if circuit_flip17(fn, v):
            return 17, model
        if fn in ("circuit::rc_charging_voltage", "circuit::rl_current_rise"):
            return 29, model
        return None
    if fn == "circuit::resistance_from_geometry" and math.isnan(v[2]):
        return 26, ip_err("cross-section area must be positive, got NaN")
    if fn == "circuit::parallel_plate_capacitance" and math.isnan(v[2]):
        return 26, ip_err("plate separation must be positive, got NaN")
    m13 = circuit_port13(fn, v)
    return None if m13 is None else (13, m13)


def circuit_adr_rows(path: Path) -> list:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        toks = [] if ins == "-" else [parse_tok(x) for x in ins.split(" ")]
        r = circuit_port(fn, circuit_vals(toks))
        if r is None:
            continue
        row, model = r
        old = wave_toks([parse_tok(x) for x in outs.split(" ")])
        if not toks_same(model, old):
            out.append((fn, ins, wire(model), str(row), f"oracle {outs}; ADR 0002 row {row}"))
    return out

# ── relativity (M2 part 1d-1): rust-old/src/relativity.rs ─────────────────────────────────────────────────────────
# Each function transcribed in Rust's order on a set of vector operations: Rust's plain FieldVector chains (which must
# reproduce the oracle), or the port's, where FieldVector's magnitude, magnitude_sq, dot and cross take ADR 0002 row
# 18's unbounded chain at a range edge (fv_governed / fv_model, as material's and wave's). The port's model adds row 27.
def rel_exact_gamma(vs, c, mode):
    """ADR 0002 row 30: RN(c/√(c² − Σv²)) (mode 0) or RN(c²/(c² − Σv²)) (mode 1) from the exact difference, or None
    unless |v| < c exactly (src/crmath_gamma.cyr's _bj_cr_lorentz)."""
    if any(x != x or abs(x) >= c for x in vs):
        return None
    d = Fraction(c) ** 2 - sum(Fraction(x) ** 2 for x in vs)
    if d <= 0:
        return None
    q = Fraction(c) ** 2 / d
    return rn_sqrt(q) if mode == 0 else rn_frac(q)


# ── an unbounded-exponent chain with IEEE's signed zeros (M2 part 2d-3d) ──────────────────────────────────────────
# A chain value is a Fraction (nonzero, 53 bits) or the float ±0, so that zeros keep IEEE's signs as kimiya's
# _bj_k_ux_* do: a product's or quotient's zero takes the xor of the signs, x + ±0 = x, two zeros add as IEEE's,
# and an exact cancellation is +0. final=True rounds the operation once into binary64 instead of to 53 bits.
def uz(x):
    return x if x == 0 else Fraction(x)


def _uz_sign(a) -> float:
    return math.copysign(1.0, a) if isinstance(a, float) else (1.0 if a > 0 else -1.0)


def _uz_round(v, final):
    return rn_frac(v) if final else rnd53(v)


def uz_mul(a, b, final=False):
    if isinstance(a, float) or isinstance(b, float):
        return _uz_sign(a) * _uz_sign(b) * 0.0
    return _uz_round(a * b, final)


def uz_div(a, b, final=False):
    if isinstance(a, float):
        return _uz_sign(a) * _uz_sign(b) * 0.0
    return _uz_round(a / b, final)


def uz_add(a, b, final=False):
    if isinstance(a, float) and isinstance(b, float):
        return a + b
    if isinstance(b, float):
        return _uz_round(a, final) if final else a
    if isinstance(a, float):
        return _uz_round(b, final) if final else b
    v = a + b
    return 0.0 if v == 0 else _uz_round(v, final)


def uz_sub(a, b, final=False):
    return uz_add(a, -b, final)


def uz_out(a) -> float:
    """A chain value already rounded into binary64 (a final operation's result) as a float."""
    return a if isinstance(a, float) else float(a)


def uz_cross(a, b, i):
    j, k = (i + 1) % 3, (i + 2) % 3
    return uz_sub(uz_mul(uz(a[j]), uz(b[k])), uz_mul(uz(a[k]), uz(b[j])))


def uz_dot3(a, b):
    return uz_add(uz_add(uz_mul(uz(a[0]), uz(b[0])), uz_mul(uz(a[1]), uz(b[1]))), uz_mul(uz(a[2]), uz(b[2])))


def uz_sumsq(a):
    return uz_dot3(a, a)


def rel_box(xs, lo, hi) -> bool:
    """src/relativity.cyr's _bj_rel_box_n: each value ±0 or of biased exponent in [lo, hi] (NaN and ±∞ fail)."""
    for x in xs:
        if x == 0:
            continue
        if not math.isfinite(x):
            return False
        be = (struct.unpack("<Q", struct.pack("<d", abs(x)))[0] >> 52) & 0x7FF
        if be < lo or be > hi:
            return False
    return True


def rel_band(s) -> bool:
    """src/crmath_vx.cyr's band: |s| in [2^-900, DBL_MAX]."""
    return 2.0 ** -900 <= abs(s) <= DBL_MAX


def rel_dot4(a, b, d):
    """src/relativity.cyr's −a₀b₀ + a₁b₁ + a₂b₂ + a₃b₃ (M2 part 2d-3c): Rust's d in the band, for an exact 0 of exact
    products or a non-finite operand; else each product and sum to 53 bits (rnd53), the last sum rounded once."""
    fs = [-a[0]] + list(a[1:])
    if rel_band(d) or not _fin(*a, *b):
        return d
    if d == 0 and all(vx_exact_p(x, y, x * y) for x, y in zip(fs, b)):
        return d
    F = Fraction
    ps = [rnd53(F(x) * F(y)) for x, y in zip(fs, b)]
    return rn_frac(rnd53(rnd53(ps[0] + ps[1]) + ps[2]) + ps[3])


def rel_ff_row(r, u, q, sv):
    """One row of the port's four-force, q·(r₀u₀ + r₁u₁ + r₂u₂ + r₃u₃) (M2 part 2d-3c), given Rust's sum sv."""
    k = q * sv
    if rel_band(sv) or not _fin(*r, *u, q):
        return k
    if sv == 0 and all(vx_exact_p(x, y, x * y) for x, y in zip(r, u)):
        return k
    F = Fraction
    ps = [rnd53(F(x) * F(y)) for x, y in zip(r, u)]
    S = rnd53(rnd53(rnd53(ps[0] + ps[1]) + ps[2]) + ps[3])
    if S == 0:
        return q * 0.0
    if q == 0:
        return q * (1.0 if S > 0 else -1.0)
    return rn_frac(F(q) * S)


def rel_boost_edge(g, sp, pv, qv, d, C):
    """src/relativity.cyr's _bj_rel_boost_edge: γ(p − β·q) with β = v/c (M2 part 2d-3c)."""
    r = g * d
    if not _fin(pv, qv):
        return r
    b = fdiv(sp, C)
    if (abs(b) > VX_MN or sp == 0) and (vx_exact_p(b, qv, b * qv) or rel_band(d)):
        return r
    F = Fraction
    D = rnd53(F(pv) - rnd53(rnd53(F(sp) / F(C)) * F(qv)))
    return g * 0.0 if D == 0 else rn_frac(F(g) * D)


def rel_ops(port: bool, gamma=None, r13=None, r17=None):
    def gov(fn, args, plain):
        if port and all(math.isfinite(x) for x in args) and fv_governed(fn, args) == 18:
            return [b for _k, b in fv_model(fn, args)]
        return plain

    def mag(a):
        return gov("FieldVector::magnitude", list(a), [fsqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])])[0]

    def mag_sq(a):
        return gov("FieldVector::magnitude_sq", list(a), [a[0] * a[0] + a[1] * a[1] + a[2] * a[2]])[0]

    def dot(a, b):
        return gov("FieldVector::dot", list(a) + list(b), [a[0] * b[0] + a[1] * b[1] + a[2] * b[2]])[0]

    def cross(a, b):
        plain = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        return gov("FieldVector::cross", list(a) + list(b), plain)

    return {"mag": mag, "mag_sq": mag_sq, "dot": dot, "cross": cross, "port": port, "ux": port,
            "gamma": port if gamma is None else gamma, "r13": port if r13 is None else r13,
            "r17": port if r17 is None else r17}


def rel_r_below(a, b, r_mag, O) -> bool:
    """The Singularity test |a − b| < 10^-30: exact on the inputs (Σ (aᵢ − bᵢ)² < 10^-60, row 17's rule) for the port,
    Rust's r_mag < 1e-30 otherwise or for a non-finite input."""
    if O["r17"] and _fin(*a, *b):
        return sum((Fraction(a[i]) - Fraction(b[i])) ** 2 for i in range(3)) < E30 * E30
    return r_mag < 1e-30


# Row 13 (M2 part 2b-4, the user's "round once" decision of 2026-11-12): relativity.rs's products and quotients of
# inputs and constants, correctly rounded where the operands are finite (exact values; c³ of c's double, 4π and 6·π on
# Rust's π); the dual tensor's components are F's own (a non-finite one keeps Rust's chain).
REL13_FNS = ("EmTensor::dual", "relativity::retarded_scalar_potential", "relativity::retarded_vector_potential",
             "relativity::larmor_power")


def rel_dual13(c, rust):
    """The port's dual: Rust's tensor with each component whose source in F is finite replaced by ±that source."""
    out = [b for _k, b in rust]
    for k, src, s in ((1, 14, -1), (2, 7, -1), (3, 9, -1), (4, 14, 1), (8, 7, 1), (12, 9, 1), (6, 12, 1),
                      (7, 8, -1), (9, 12, -1), (11, 4, 1), (13, 8, 1), (14, 4, -1)):
        if math.isfinite(c[src]):
            out[k] = -c[src] if s < 0 else c[src]
    return [("f", x) for x in out]


def rel_consts(root: Path) -> dict:
    c = golden_consts(root)
    C, C2 = c["field::SPEED_OF_LIGHT"], c["field::SPEED_OF_LIGHT_SQ"]
    return {"C": C, "C2": C2, "C3": C * C2, "EPS0": c["field::EPSILON_0"], "MU0": c["field::MU_0"]}


def rel_from_fields(e, b, K):
    C = K["C"]
    ex, ey, ez = fdiv(e[0], C), fdiv(e[1], C), fdiv(e[2], C)
    return [0.0, -ex, -ey, -ez, ex, 0.0, -b[2], b[1], ey, b[2], 0.0, -b[0], ez, -b[1], b[0], 0.0]


def rel_gamma(speed, K, exact=False):
    if exact:
        g = rel_exact_gamma([speed], K["C"], 0)
        if g is None:
            return ip_err(f"speed ({rust_display(speed)} m/s) must be less than c"), None
        return None, g
    beta_sq = fdiv(speed * speed, K["C2"])
    if math.isnan(beta_sq) or beta_sq >= 1.0:
        return ip_err(f"speed ({rust_display(speed)} m/s) must be less than c"), None
    return None, fdiv(1.0, fsqrt(1.0 - beta_sq))


def rel_rust(fn, v, O, K):
    """relativity.rs transcribed on the vector operations O; None for another name."""
    C, C2, EPS0, MU0 = K["C"], K["C2"], K["EPS0"], K["MU0"]
    if fn == "relativity::lorentz_factor":
        e, g = rel_gamma(v[0], K, O["gamma"])
        return e or okf(g)
    if fn == "relativity::beta":
        return okf(fdiv(v[0], C))
    if fn == "EmTensor::from_fields":
        return [("f", x) for x in rel_from_fields(v[0:3], v[3:6], K)]
    if fn.startswith("EmTensor::"):
        c = v[0:16]
        ef = [c[4] * C, c[8] * C, c[12] * C]
        bf = [c[14], c[7], c[9]]
        if fn == "EmTensor::electric_field":
            return ok3(ef)
        if fn == "EmTensor::magnetic_field":
            return ok3(bf)
        F = Fraction
        fin = _fin(c[4], c[8], c[12], *bf)
        if fn == "EmTensor::first_invariant":
            d = O["mag_sq"](bf) - fdiv(O["mag_sq"](ef), C2)
            if O["ux"] and fin and not rel_band(d):
                # Row 18's rule (M2 part 2d-3c): E = c·F, E², B², E²/c² and the difference to 53 bits, 2D rounded once.
                E = [rnd53(F(x) * F(C)) for x in (c[4], c[8], c[12])]
                e2 = rnd53(rnd53(rnd53(E[0] ** 2) + rnd53(E[1] ** 2)) + rnd53(E[2] ** 2))
                b2 = rnd53(rnd53(rnd53(F(bf[0]) ** 2) + rnd53(F(bf[1]) ** 2)) + rnd53(F(bf[2]) ** 2))
                return okf(rn_frac(2 * rnd53(b2 - rnd53(e2 / F(C2)))))
            return okf(2.0 * d)
        if fn == "EmTensor::second_invariant":
            d = O["dot"](ef, bf)
            exact = all(vx_exact_p(x, C, x * C) for x in (c[4], c[8], c[12]))
            if O["ux"] and fin and not (VX_MN < abs(d) <= DBL_MAX and exact):
                # Row 18's rule: E = c·F, E⋅B to 53 bits, −E⋅B/c rounded once.
                E = [rnd53(F(x) * F(C)) for x in (c[4], c[8], c[12])]
                ps = [rnd53(E[k] * F(bf[k])) for k in range(3)]
                S = rnd53(rnd53(ps[0] + ps[1]) + ps[2])
                if S == 0:
                    # an exact 0: IEEE's signs where every product is 0, else a cancellation's +0
                    return okf(fdiv(-d, C)) if all(x == 0 for x in ps) else okf(-0.0)
                return okf(rn_frac(-S / F(C)))
            return okf(fdiv(-d, C))
        if fn == "EmTensor::dual":
            inv = -fdiv(1.0, C)
            rust = [("f", x) for x in rel_from_fields([x * C for x in bf], [x * inv for x in ef], K)]
            return rel_dual13(c, rust) if O["r13"] else rust
        if fn == "EmTensor::four_force":
            q, u = v[16], v[17:21]
            ul = [u[0], -u[1], -u[2], -u[3]]
            out = []
            for i in range(4):
                r = c[4 * i:4 * i + 4]
                sv = r[0] * ul[0] + r[1] * ul[1] + r[2] * ul[2] + r[3] * ul[3]
                out.append(("f", rel_ff_row(r, ul, q, sv) if O["ux"] else q * sv))
            return out
        return None
    if fn == "relativity::lorentz_transform_fields":
        e, b, w = v[0:3], v[3:6], v[6:9]
        if O["gamma"]:
            g = rel_exact_gamma(w, C, 0)
            if g is None:
                return ip_err("boost velocity must be less than c")
        else:
            v_sq = O["mag_sq"](w)
            if math.isnan(v_sq) or v_sq >= C2:
                return ip_err("boost velocity must be less than c")
            g = fdiv(1.0, fsqrt(1.0 - fdiv(v_sq, C2)))
        k = fdiv(g * g, (g + 1.0) * C2)
        vxb = O["cross"](w, b)
        ke = k * O["dot"](w, e)
        ep = [(e[i] + vxb[i]) * g - w[i] * ke for i in range(3)]
        vxe = O["cross"](w, e)
        inv = fdiv(1.0, C2)
        kb = k * O["dot"](w, b)
        bp = [(b[i] - vxe[i] * inv) * g - w[i] * kb for i in range(3)]
        if O["ux"] and not (rel_box(e + b + w, 723, 1923) and g < 2.0 ** 27) and _fin(*e, *b, *w):
            # Row 18's rule (M2 part 2d-3d): outside the box every step to 53 bits, each output's difference once.
            G = uz(g)
            kk = uz_div(uz_mul(G, G), uz_mul(uz_add(G, Fraction(1)), uz(C2)))
            ke2, kb2, iv = uz_mul(kk, uz_dot3(w, e)), uz_mul(kk, uz_dot3(w, b)), uz(inv)
            ep = [uz_out(uz_sub(uz_mul(uz_add(uz(e[i]), uz_cross(w, b, i)), G), uz_mul(uz(w[i]), ke2), True))
                  for i in range(3)]
            bp = [uz_out(uz_sub(uz_mul(uz_sub(uz(b[i]), uz_mul(uz_cross(w, e, i), iv)), G), uz_mul(uz(w[i]), kb2),
                                True)) for i in range(3)]
        return ok3(ep) + ok3(bp)
    if fn == "relativity::lorentz_transform_x":
        e, b, s = v[0:3], v[3:6], v[6]
        err, g = rel_gamma(s, K, O["gamma"])
        if err:
            return err
        ep = [e[0], g * (e[1] - s * b[2]), g * (e[2] + s * b[1])]
        bp = [b[0], g * (b[1] + fdiv(s * e[2], C2)), g * (b[2] - fdiv(s * e[1], C2))]
        if O["ux"] and not rel_box(e + b + [s], 723, 1923) and _fin(*e, *b):
            G, S, c2 = uz(g), uz(s), uz(C2)
            ep = [e[0], uz_out(uz_mul(G, uz_sub(uz(e[1]), uz_mul(S, uz(b[2]))), True)),
                  uz_out(uz_mul(G, uz_add(uz(e[2]), uz_mul(S, uz(b[1]))), True))]
            bp = [b[0], uz_out(uz_mul(G, uz_add(uz(b[1]), uz_div(uz_mul(S, uz(e[2])), c2)), True)),
                  uz_out(uz_mul(G, uz_sub(uz(b[2]), uz_div(uz_mul(S, uz(e[1])), c2)), True))]
        return ok3(ep) + ok3(bp)
    if fn == "FourVector::new":
        return [("f", x) for x in v[0:4]]
    if fn == "FourVector::from_time_position":
        return [("f", C * v[0]), ("f", v[1]), ("f", v[2]), ("f", v[3])]
    if fn == "FourVector::dot":
        a, b = v[0:4], v[4:8]
        d = -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
        return okf(rel_dot4(a, b, d) if O["ux"] else d)
    if fn == "FourVector::interval_sq":
        a = v[0:4]
        d = -a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3]
        return okf(rel_dot4(a, a, d) if O["ux"] else d)
    if fn == "FourVector::spatial_magnitude":
        a = v[0:4]
        return okf(O["mag"](a[1:4]) if O["ux"] else fsqrt(a[1] * a[1] + a[2] * a[2] + a[3] * a[3]))
    if fn == "FourVector::boost_x":
        a, s = v[0:4], v[4]
        err, g = rel_gamma(s, K, O["gamma"])
        if err:
            return err
        bv = fdiv(s, C)
        d1, d2 = a[0] - bv * a[1], a[1] - bv * a[0]
        t2, x2 = g * d1, g * d2
        if O["ux"] and not (abs(bv) > VX_MN and rel_band(d1) and rel_band(d2)):
            t2, x2 = rel_boost_edge(g, s, a[0], a[1], d1, C), rel_boost_edge(g, s, a[1], a[0], d2, C)
        return [("f", t2), ("f", x2), ("f", a[2]), ("f", a[3])]
    if fn == "relativity::retarded_scalar_potential":
        q, d = v
        e = radial_check(d)
        if e:
            return e
        if O["r13"] and _fin(q, d):
            return okf(rn_mono([q], [4.0 * PI, EPS0, d]))
        return okf(fdiv(q, 4.0 * PI * EPS0 * d))
    if fn == "relativity::retarded_vector_potential":
        q, w, d = v[0], v[1:4], v[4]
        e = radial_check(d)
        if e:
            return e
        k = fdiv(MU0 * q, 4.0 * PI * d)
        if O["r13"] and _fin(q, d):
            return ok3([rn_mono([MU0, q, x], [4.0 * PI, d]) if math.isfinite(x) else x * k for x in w])
        return ok3([x * k for x in w])
    if fn in ("relativity::lienard_wiechert_scalar_potential", "relativity::lienard_wiechert_vector_potential"):
        q, r, bv = v[0], v[1:4], v[4:7]
        r_mag = O["mag"](r)
        if math.isnan(r_mag):
            return ip_err("R vector must not contain NaN")
        if O["port"] and not _fin(*r):
            return ip_err("R vector must be finite")    # row 27 amended (M2 part 2d-3e, the user's answer)
        if rel_r_below(r, [0.0, 0.0, 0.0], r_mag, O):
            return SINGULAR
        if O["gamma"]:
            if rel_exact_gamma(bv, 1.0, 1) is None:
                return ip_err("β must be < 1")
        else:
            b_sq = O["mag_sq"](bv)
            if math.isnan(b_sq) or b_sq >= 1.0:
                return ip_err("β must be < 1")
        if O["r17"] and _fin(*r):
            r_eff = rn_lw(r, bv, 0)          # row 17 (M2 part 2c-4): R_eff from the exact S and P, rounded once
        else:
            r_eff = r_mag - O["dot"](r, bv)
        if fn.endswith("scalar_potential"):
            if O["ux"] and r_eff > DBL_MAX and math.isfinite(q):
                # Row 18's rule (M2 part 2d-3e): R_eff past DBL_MAX to 53 bits, 4πε₀R_eff too, the quotient once.
                d = uz_mul(uz(4.0 * PI * EPS0), rnd53_lw(r, bv))
                return okf(uz_out(uz_div(uz(q), d, True)))
            return okf(fdiv(q, 4.0 * PI * EPS0 * r_eff))
        mq = MU0 * q
        num = mq * C
        k = fdiv(num, 4.0 * PI * r_eff)
        plain = r_eff <= DBL_MAX and (q == 0 or (abs(mq) > VX_MN and VX_MN < abs(num) <= DBL_MAX
                                                 and VX_MN < abs(k) <= DBL_MAX))
        if O["ux"] and not plain and math.isfinite(q) and r_eff > 0:   # (R_eff ≥ 2^-367 in the port)
            reff = uz(r_eff) if r_eff <= DBL_MAX else rnd53_lw(r, bv)
            kk = uz_div(uz_mul(uz_mul(uz(MU0), uz(q)), uz(C)), uz_mul(uz(4.0 * PI), reff))
            return ok3([uz_out(uz_mul(uz(x), kk, True)) for x in bv])
        return ok3([x * k for x in bv])
    if fn == "relativity::lienard_wiechert_e":
        q, r, bv, bd = v[0], v[1:4], v[4:7], v[7:10]
        r_mag = O["mag"](r)
        if O["port"] and math.isnan(r_mag):
            return ip_err("R vector must not contain NaN")
        if O["port"] and not _fin(*r):
            return ip_err("R vector must be finite")    # row 27 amended (M2 part 2d-3e)
        if rel_r_below(r, [0.0, 0.0, 0.0], r_mag, O):
            return SINGULAR
        if O["gamma"]:
            g2 = rel_exact_gamma(bv, 1.0, 1)
            if g2 is None:
                return ip_err("β must be < 1")
        else:
            b_sq = O["mag_sq"](bv)
            if math.isnan(b_sq) or b_sq >= 1.0:
                return ip_err("β must be < 1")
            g2 = fdiv(1.0, 1.0 - b_sq)
        inv = fdiv(1.0, r_mag)
        n = [x * inv for x in r]
        if O["r17"] and _fin(*r):
            if lw_below(r, bv):              # row 17 (M2 part 2c-4): κ < 10^-30 exactly, κ rounded once
                return dz_err("Liénard-Wiechert denominator (1 − n̂⋅β) is zero")
            kappa = rn_lw(r, bv, 1)
        else:
            kappa = 1.0 - O["dot"](n, bv)
            if abs(kappa) < 1e-30:
                return dz_err("Liénard-Wiechert denominator (1 − n̂⋅β) is zero")
        k3 = kappa * kappa * kappa
        pre = fdiv(q, 4.0 * PI * EPS0)
        nmb = [n[i] - bv[i] for i in range(3)]
        if O["r17"] and _fin(*r):                # row 17 (M2 part 2d-1): n̂ − β rounded once, part by part
            nmb = [rn_nb(r, bv, i) if r[i] != 0.0 else nmb[i] for i in range(3)]
        fv = fdiv(pre, g2 * k3 * r_mag * r_mag)
        vt = [x * fv for x in nmb]
        co = O["cross"](n, O["cross"](nmb, bd))
        fa = fdiv(pre, k3 * r_mag * C)
        at = [x * fa for x in co]
        box = ((q == 0 or 2.0 ** -100 <= abs(q) < 2.0 ** 201) and g2 < 2.0 ** 100 and rel_box(r + bd, 923, 1223)
               and rel_box(nmb, 923, 2046))
        if O["ux"] and not box and _fin(q, *bd):
            # Row 18's rule (M2 part 2d-3f): outside the box the whole tail to 53 bits, each part's sum rounded once.
            M = rnd53_sqrt(uz_sumsq(r))
            iv = uz_div(Fraction(1), M)
            nn = [uz_mul(uz(x), iv) for x in r]
            K = uz(kappa)
            k3u = uz_mul(uz_mul(K, K), K)
            pr = uz_div(uz(q), uz(4.0 * PI * EPS0))
            fvu = uz_div(pr, uz_mul(uz_mul(uz_mul(uz(g2), k3u), M), M))
            fau = uz_div(pr, uz_mul(uz_mul(k3u, M), uz(C)))
            inner = [uz_cross(nmb, bd, i) for i in range(3)]
            out = []
            for i in range(3):
                j, k = (i + 1) % 3, (i + 2) % 3
                o = uz_sub(uz_mul(nn[j], inner[k]), uz_mul(nn[k], inner[j]))
                out.append(uz_out(uz_add(uz_mul(uz(nmb[i]), fvu), uz_mul(o, fau), True)))
            return ok3(out)
        return ok3([vt[i] + at[i] for i in range(3)])
    if fn == "relativity::lienard_wiechert_b":
        cr = O["cross"](v[0:3], v[3:6])
        inv = fdiv(1.0, C)
        return ok3([x * inv for x in cr])
    if fn == "relativity::larmor_power":
        q, a = v
        if O["r13"] and _fin(q, a):
            return okf(rn_mono([q, q, a, a], [6.0, PI, EPS0, C, C, C]))
        return okf(fdiv(q * q * a * a, 6.0 * PI * EPS0 * K["C3"]))
    if fn == "relativity::relativistic_larmor_power":
        q, w, a = v[0], v[1:4], v[4:7]
        if O["gamma"]:
            g = rel_exact_gamma(w, C, 0)
            if g is None:
                return ip_err("velocity must be less than c")
        else:
            v_sq = O["mag_sq"](w)
            if math.isnan(v_sq) or v_sq >= C2:
                return ip_err("velocity must be less than c")
            g = fdiv(1.0, fsqrt(1.0 - fdiv(v_sq, C2)))
        g2 = g * g
        g6 = g2 * g2 * g2
        a_sq = O["mag_sq"](a)
        vxa_sq = O["mag_sq"](O["cross"](w, a))
        den = 6.0 * PI * EPS0 * K["C3"]
        pre = fdiv(q * q, den)
        if (O["ux"] and not (rel_box([q], 823, 1423) and rel_box(a, 823, 1503) and rel_box(w, 823, 2046)
                             and g < 2.0 ** 27) and _fin(q, *a)):
            G = uz(g)
            g2u = uz_mul(G, G)
            g6u = uz_mul(uz_mul(g2u, g2u), g2u)
            x = [uz_cross(w, a, i) for i in range(3)]
            x2 = uz_add(uz_add(uz_mul(x[0], x[0]), uz_mul(x[1], x[1])), uz_mul(x[2], x[2]))
            pr = uz_div(uz_mul(uz(q), uz(q)), uz(den))
            dd = uz_sub(uz_sumsq(a), uz_div(x2, uz(C2)))
            return okf(uz_out(uz_mul(uz_mul(pr, g6u), dd, True)))
        return okf(pre * g6 * (a_sq - fdiv(vxa_sq, C2)))
    return None


def validate_relativity(root: Path, path: Path) -> list:
    _header, rows = read_tsv(path)
    K, O = rel_consts(root), rel_ops(False)
    errs = check_platform(_header)
    if errs:
        return errs
    for n, (fn, _cls, ins, outs) in enumerate(rows, 1):
        v = [b for _k, b in ins]
        if fn in ("relativity::solve_retarded_time", "relativity::lienard_wiechert_fields"):
            want = rel_solver_rust(fn, ins, GLIBC_P, O, K)
        else:
            want = rel_rust(fn, v, O, K)
        if want is None:
            errs.append(f"relativity row {n}: no transcription for {fn}")
            continue
        if not toks_same(want, wave_toks(outs)):
            errs.append(f"relativity row {n}: {fn}{tuple(v)} recorded {wave_toks(outs)}, transcription {want}")
    return errs


# ── relativity's solver (M2 part 1d-2): solve_retarded_time and lienard_wiechert_fields ──────────────────────────
# The trajectories of m_relativity.rs's `trajectory`, by number, on the primitives P (sin and cos for 2 and 4).
def rel_trajectory(tid, p, t, P, C):
    if tid == 0:
        return [p[0], p[1], p[2]], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
    if tid == 1:
        dt = t - p[6]
        return [p[0] + p[3] * dt, p[1] + p[4] * dt, p[2] + p[5] * dt], [p[3], p[4], p[5]], [0.0, 0.0, 0.0]
    if tid == 2:
        r, w = p[0], p[1]
        s, c = P["sin"](w * t), P["cos"](w * t)
        a2 = r * w * w
        return [r * c, r * s, 0.0], [-r * w * s, r * w * c, 0.0], [-a2 * c, -a2 * s, 0.0]
    if tid == 3:
        a = p[0]
        u = fdiv(a * t, C)
        g = fsqrt(1.0 + u * u)
        return [fdiv(C * C, a) * (g - 1.0), 0.0, 0.0], [fdiv(a * t, g), 0.0, 0.0], [fdiv(a, g * g * g), 0.0, 0.0]
    if tid == 4:
        amp, w = p[0], p[1]
        s, c = P["sin"](w * t), P["cos"](w * t)
        return [0.0, 0.0, amp * s], [0.0, 0.0, amp * w * c], [0.0, 0.0, -(amp * w * w) * s]
    return [p[0], p[1], p[2]], [0.0, 0.0, 0.0], [p[3], p[4], p[5]]


EPS52 = 2.0 ** -52
MIN_NORMAL = 2.0 ** -1022


def rel_sample_plain(dx, dy, dz, d2, light, dt, vel, mvt, m4, at) -> bool:
    """src/relativity.cyr's _bj_rel_sample_edge test (_bj_rel_sample's inline test is a case of it): whether the
    light-cone sample's Rust chain is the unbounded one."""
    if not _fin(dx, dy, dz):
        return False
    if not 2.0 ** -900 <= d2 <= DBL_MAX and not (dx == 0 and dy == 0 and dz == 0):
        return False
    if not (light <= DBL_MAX and (light > VX_MN or dt == 0)):
        return False
    if d2 > 0:
        ps = [dx * vel[0], dy * vel[1], dz * vel[2]]
        rv = ps[0] + ps[1] + ps[2]
        if not 2.0 ** -900 <= abs(rv) <= DBL_MAX:
            if rv != 0 or not all(vx_exact_p(a, b, a * b) for a, b in zip((dx, dy, dz), vel)):
                return False
    if not mvt <= DBL_MAX:
        return False
    if m4 >= 2.0 ** -900:
        return True
    return m4 == 0 and (at == 0 or max(abs(x) for x in vel) == 0)


def rel_sample_ux(obs, t, tau, pos, vel, C):
    """src/relativity.cyr's _bj_rel_sample_ux (M2 part 2d-3g): the light-cone sample with an unbounded exponent."""
    t_r = t - tau
    dt = t - t_r
    d = [uz_sub(uz(o), uz(p)) for o, p in zip(obs, pos)]
    s2 = uz_add(uz_add(uz_mul(d[0], d[0]), uz_mul(d[1], d[1])), uz_mul(d[2], d[2]))
    m = 0.0 if isinstance(s2, float) else rnd53_sqrt(s2)
    lt = uz_mul(uz(C), uz(dt))
    gg = uz_sub(m, lt)
    g = uz_out(uz_sub(m, lt, True))
    if isinstance(m, float):
        dg = -C
    else:
        rv = uz_add(uz_add(uz_mul(d[0], uz(vel[0])), uz_mul(d[1], uz(vel[1]))), uz_mul(d[2], uz(vel[2])))
        dg = uz_out(uz_sub(uz_div(rv, m), uz(C), True))

    def mabs(xs):
        return max(abs(x) for x in xs)
    cands = [m, lt, uz(mabs(obs)), uz(mabs(pos)), uz_mul(uz(mabs(vel)), uz(abs(t_r)))]
    mx = max(cands, key=lambda z: Fraction(0) if isinstance(z, float) else z)
    noise = uz_out(uz_mul(Fraction(1, 2 ** 51), mx, True)) if not isinstance(mx, float) else 0.0
    return {"tau": tau, "g": g, "G": gg, "dg": dg, "noise": noise}


def _gval(z):
    return Fraction(0) if isinstance(z, float) else z


BRACKET_ERR = ip_err("retarded time solver could not bracket the light-cone root (is the trajectory faster than light?)")


def rel_solve(obs, t, traj, max_it, C, port=False):
    """solve_retarded_time transcribed: ("ok", t_r) or ("err", tokens). traj(t) -> (pos, vel). port: a non-finite state
    met while the bracket expands reports the bracket error (ADR 0002 row 31, the user's answer of 2026-11-13), and
    (row 18's rule, M2 part 2d-3g) the light-cone sample runs with an unbounded exponent outside its guard, every
    decision on g taken on its 53-bit value G (τ₀ = G/c, the comparisons, the Newton point and step)."""
    if not (math.isfinite(t) and all(math.isfinite(c) for c in obs)):
        return "err", ip_err("observer position and time must be finite, got [" + ", ".join(rust_debug(c) for c in obs)
                             + f"] at t = {rust_display(t)}")

    def ev(tau):
        t_r = t - tau
        pos, vel = traj(t_r)
        if not all(math.isfinite(c) for c in pos + vel):
            return None, ip_err(f"trajectory returned a non-finite state at t = {rust_display(t_r)}")
        dx, dy, dz = obs[0] - pos[0], obs[1] - pos[1], obs[2] - pos[2]
        d2 = dx * dx + dy * dy + dz * dz
        dist = fsqrt(d2)
        dt = t - t_r
        light = C * dt
        dg = fdiv(dx * vel[0] + dy * vel[1] + dz * vel[2], dist) - C if dist > 0.0 else -C

        def mabs(v):
            return rust_max(rust_max(abs(v[0]), abs(v[1])), abs(v[2]))
        mvt = mabs(vel) * abs(t_r)
        m4 = rust_max(rust_max(rust_max(rust_max(dist, light), mabs(obs)), mabs(pos)), mvt)
        noise = 2.0 * EPS52 * m4
        g = dist - light
        if port and not rel_sample_plain(dx, dy, dz, d2, light, dt, vel, mvt, m4, abs(t_r)):
            return rel_sample_ux(obs, t, tau, pos, vel, C), None
        return {"tau": tau, "g": g, "G": uz(g) if port else None, "dg": dg, "noise": noise}, None

    def tol(tau):
        return 4.0 * EPS52 * rust_max(abs(t), tau)

    def sign(s):
        z = _gval(s["G"]) if port else s["g"]
        return (z > 0) - (z < 0)

    def abs_lt(a, b):
        return abs(_gval(a["G"])) < abs(_gval(b["G"])) if port else abs(a["g"]) < abs(b["g"])

    def le_noise(s):
        return abs(_gval(s["G"])) <= Fraction(s["noise"]) if port else abs(s["g"]) <= s["noise"]

    def newton_pt(s):
        """(τ − g/dg, |g/dg|): Rust's doubles, or (port) g/dg to 53 bits and each rounded once."""
        if not port:
            gd = fdiv(s["g"], s["dg"])
            return s["tau"] - gd, abs(gd)
        q = uz_div(s["G"], uz(s["dg"]))
        return uz_out(uz_sub(uz(s["tau"]), q, True)), abs(uz_out(uz_mul(q, Fraction(1), True)))

    def g2_le(s, last_step):
        if not port:
            return abs(2.0 * s["g"]) <= abs(last_step * s["dg"])
        return abs(2 * _gval(s["G"])) <= abs(_gval(uz_mul(uz(last_step), uz(s["dg"]))))

    s0, e = ev(0.0)
    if e:
        return "err", e
    if sign(s0) == 0:
        return "ok", t
    lo = s0
    tau0 = uz_out(uz_div(s0["G"], uz(C), True)) if port else fdiv(s0["g"], C)
    hi, e = ev(rust_max(tau0, MIN_NORMAL))
    if e:
        return "err", e
    expansions = 0
    while sign(hi) > 0:
        nxt = 2.0 * hi["tau"]
        if expansions >= 2100 or not math.isfinite(nxt):
            return "err", BRACKET_ERR
        lo = hi
        hi, e = ev(nxt)
        if e:
            return "err", BRACKET_ERR if port else e
        expansions += 1
    cur = lo if abs_lt(lo, hi) else hi
    last_step = hi["tau"] - lo["tau"]
    it = 0
    while True:
        if sign(cur) == 0:
            conv = cur["tau"]
            break
        if le_noise(cur):
            if cur["dg"] < 0.0:
                conv = rust_min(rust_max(newton_pt(cur)[0], lo["tau"]), hi["tau"])
            else:
                conv = cur["tau"]
            break
        if it == max_it:
            return "err", ip_err(f"retarded time solver did not converge in {max_it} iterations")
        it += 1
        newton, nstep = newton_pt(cur)
        use_newton = (cur["dg"] < 0.0 and newton >= lo["tau"] and newton <= hi["tau"] and g2_le(cur, last_step))
        if use_newton:
            nxt, step = newton, nstep
        else:
            half = 0.5 * (hi["tau"] - lo["tau"])
            nxt, step = lo["tau"] + half, half
        if step <= tol(nxt) or nxt == cur["tau"]:
            conv = nxt
            break
        last_step = step
        cur, e = ev(nxt)
        if e:
            return "err", e
        if sign(cur) > 0:
            lo = cur
        else:
            hi = cur
    return "ok", t - conv


def rel_solver_rust(fn, ins, P, O, K):
    """A solver or fields row's outputs, on the primitives P and vector operations O."""
    C = K["C"]
    if fn == "relativity::solve_retarded_time":
        tid, obs, t, mx = ins[0][1], [ins[1][1], ins[2][1], ins[3][1]], ins[4][1], ins[5][1]
        ps = [b for _k, b in ins[7:]]
        max_it = mx if mx >= 0 else mx + (1 << 64)
        kind, r = rel_solve(obs, t, lambda tt: rel_trajectory(tid, ps, tt, P, C)[0:2], max_it, C, O["port"])
        return okf(r) if kind == "ok" else r
    q, tid, obs, t = ins[0][1], ins[1][1], [ins[2][1], ins[3][1], ins[4][1]], ins[5][1]
    ps = [b for _k, b in ins[7:]]
    kind, r = rel_solve(obs, t, lambda tt: rel_trajectory(tid, ps, tt, P, C)[0:2], 200, C, O["port"])
    if kind == "err":
        return r
    t_r = r
    pos, vel, acc = rel_trajectory(tid, ps, t_r, P, C)
    if not all(math.isfinite(a) for a in acc):
        return ip_err(f"trajectory returned a non-finite acceleration at t = {rust_display(t_r)}")
    rv = [obs[i] - pos[i] for i in range(3)]
    r_mag = O["mag"](rv)
    if rel_r_below(obs, pos, r_mag, O):
        return SINGULAR
    inv = fdiv(1.0, r_mag)
    n = [x * inv for x in rv]
    bv = [fdiv(x, C) for x in vel]
    bd = [fdiv(x, C) for x in acc]
    e = rel_rust("relativity::lienard_wiechert_e", [q] + rv + bv + bd, O, K)
    if e[0][0] == "E":
        return e
    ev3 = [b for _k, b in e]
    b = rel_rust("relativity::lienard_wiechert_b", n + ev3, O, K)
    return e + b


# ADR 0002 row 27 (M2 part 1d-1, row 20's rule): lienard_wiechert_e refuses a NaN R as its sibling potentials do
# ("R vector must not contain NaN"), where Rust's went on to Ok(NaN). Row 18 through FieldVector's operations.
def relativity_corrections(root: Path, path: Path):
    """The class-C solver rows (the circular orbit and the oscillation) whose outputs change when sin and cos are
    correctly rounded (MPFR), on the port's vector operations (ADR 0002 row 18). A row whose γ² (row 30) changes it is
    an ADR 0002 mapping instead (relativity_adr_rows), and so is one the exact |R| test (row 17) changes."""
    K, O, O18 = rel_consts(root), rel_ops(True), rel_ops(True, gamma=False, r17=False)
    O17 = rel_ops(True, r17=False)

    def fixes(_rows):
        out = []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                continue
            fn, cls, ins, outs, _note = line.split("\t")
            if cls != "C":
                continue
            toks = [parse_tok(x) for x in ins.split(" ")]
            want = rel_solver_rust(fn, toks, MPFR_P, O, K)
            if not toks_same(want, rel_solver_rust(fn, toks, MPFR_P, O17, K)):
                continue                    # row 17's
            if not toks_same(want, rel_solver_rust(fn, toks, MPFR_P, O18, K)):
                continue                    # row 30's
            old = wave_toks([parse_tok(x) for x in outs.split(" ")])
            if toks_same(want, old):
                continue
            out.append((fn, ins, wire(want), f"oracle {outs}; correctly rounded primitives give {wire(want)}"))
        return out
    return fixes


def relativity_adr_rows(root: Path, path: Path) -> list:
    K, R, P = rel_consts(root), rel_ops(False), rel_ops(True)
    P18 = rel_ops(True, gamma=False, r17=False)
    P17 = rel_ops(True, r17=False)                   # the port without row 17's exact |R| test
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        toks = [] if ins == "-" else [parse_tok(x) for x in ins.split(" ")]
        v = [b for _k, b in toks]
        if fn in ("relativity::solve_retarded_time", "relativity::lienard_wiechert_fields"):
            if _cls == "C":
                # its trajectory's sin and cos: on correctly rounded primitives, a mapping only where γ² (row 30) or
                # the exact |R| test (row 17) changes it; a change from the primitives alone is the corrected
                # oracle's row
                full = rel_solver_rust(fn, toks, MPFR_P, P, K)
                if not toks_same(full, rel_solver_rust(fn, toks, MPFR_P, P17, K)):
                    out.append((fn, ins, wire(full), "17", f"oracle {outs}; ADR 0002 row 17"))
                elif not toks_same(full, rel_solver_rust(fn, toks, MPFR_P, P18, K)):
                    out.append((fn, ins, wire(full), "30", f"oracle {outs}; ADR 0002 row 30"))
                continue
            rust, port = rel_solver_rust(fn, toks, GLIBC_P, R, K), rel_solver_rust(fn, toks, GLIBC_P, P, K)
            r17 = rel_solver_rust(fn, toks, GLIBC_P, P17, K)
            r18 = rel_solver_rust(fn, toks, GLIBC_P, P18, K)
        else:
            rust, port = rel_rust(fn, v, R, K), rel_rust(fn, v, P, K)
            r17, r18 = rel_rust(fn, v, P17, K), rel_rust(fn, v, P18, K)
        if toks_same(rust, port):
            continue
        # 27: the NaN R; 17: the exact |R| test; 13: the round-once functions; 18: the vector operations alone; 30: γ
        # (alone or with the vector operations).
        if (fn == "relativity::lienard_wiechert_e" and port[0] == ("E", "InvalidParameter")
                and port[1][1].endswith("must not contain NaN") and rust[0][0] != "E"):
            row = 27
        elif port[0] == ("E", "InvalidParameter") and port[1][1].endswith("R vector must be finite"):
            row = 27                                     # amended (M2 part 2d-3e): an ∞ R refused
        elif (toks_same(port, BRACKET_ERR) and rust[0] == ("E", "InvalidParameter")
              and "trajectory returned a non-finite state" in str(rust[1][1])):
            row = 31                                     # M2 part 2d-3e, the user's answer of 2026-11-13
        elif not toks_same(port, r17):
            row = 17
        elif fn in REL13_FNS:
            row = 13
        elif toks_same(port, r18):
            row = 18
        else:
            row = 30
        out.append((fn, ins, wire(port), str(row), f"oracle {outs}; ADR 0002 row {row}"))
    return out


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
def validate_charge(root: Path, path: Path, sink=None, P=None) -> list:
    """charge.rs transcribed in its order of operations (class A), its checks and error texts; the PointCharge
    methods through field.rs's transcriptions above. A sink as validate_field's (charge has no class-C row)."""
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
            if cut_r2(r_sq, a, b):
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
            if sink is not None:
                sink.append(None)
                continue
            errs.append(f"charge row {n}: no transcription for {fn}")
            continue
        if sink is not None:
            sink.append(want)
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
            if want is None or any(k != "f" for k, _ in want) or any(k != "f" for k, _ in outs):
                continue
            global EXACT_THR
            EXACT_THR = True
            w17 = field_3b(fn, v, ins, MPFR_P, K, EPS0, MU0)
            EXACT_THR = False
            if any(k != "f" for k, _ in w17):
                continue                    # the exact threshold makes it Singularity: row 17's mapping
            if fn in FV_FNS and all(isinstance(x, float) and math.isfinite(x) for x in v) and fv_governed(fn, v):
                continue                    # row 18's unbounded chain sets it
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
    "Material": [("eps_r", False), ("mu_r", False), ("conductivity", False)],
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

from exact import lw_below, rn_lw, rn_mono, rn_nb, rn_sqrt, rnd53_lw  # noqa: E402


def _fin(*xs) -> bool:
    return all(math.isfinite(x) for x in xs)


# Row 17's rule over M2 (part 2c-1, the user's "exact thresholds" decision of 2026-11-12): a test against 10^-30 of one
# rounded product, quotient or sum of two inputs, or of |R|, decided on the exact value (Fraction); a non-finite input
# keeps Rust's test.
E30 = Fraction(1, 10 ** 30)


def below30(x) -> bool:
    """|x| < 10^-30 exactly, for an exact rational x."""
    return abs(x) < E30


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


def pct(text: str) -> str:
    """tsv.rs's `s:` body: bytes 0x21–0x7e as they are (but %), the rest as %XX."""
    return "".join(chr(b) if 0x21 <= b <= 0x7E and b != 0x25 else f"%{b:02X}" for b in text.encode())


def wire(toks) -> str:
    out = []
    for k, b in toks:
        if k == "f":
            out.append(f"f:{f_to_bits(b):016x}")
        elif k == "s":
            out.append("s:" + pct(b))
        else:
            out.append(f"{k}:{b}")
    return " ".join(out)


def toks_same(x, y) -> bool:
    if x is None or y is None:
        return x is y
    return len(x) == len(y) and all(a[0] == b[0] and (same(a[1], b[1]) if a[0] == "f" else a[1] == b[1])
                                    for a, b in zip(x, y))


def adr17_rows(root: Path, module: str, path: Path) -> list:
    """Row 17's mappings: the golden rows whose outcome the exact Singularity threshold changes (Rust's chain with the
    exact test against Rust's, both on glibc's primitives as the oracle ran them), with the port's outputs (the exact
    test, MPFR's correctly rounded primitives)."""
    global EXACT_THR
    check = validate_field if module == "field" else validate_charge
    rust, exact, port = [], [], []
    try:
        EXACT_THR = False
        check(root, path, rust)
        EXACT_THR = True
        check(root, path, exact)
        check(root, path, port, MPFR_P)
    finally:
        EXACT_THR = False
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    out = []
    for line, w0, w1, wp in zip(lines, rust, exact, port):
        if toks_same(w0, w1):
            continue
        fn, _cls, ins, outs, _note = line.split("\t")
        out.append((fn, ins, wire(wp), "17", f"oracle {outs}; ADR 0002 row 17: the exact threshold"))
    return out


# Row 18 (M1 part 5d-2): field.rs's and charge.rs's vector chains at a range edge, as row 14's complex ones. The model
# is the chain with an unbounded exponent, exact values rounded to 53 bits (rnd53), the last operation once into
# binary64; normalized is its first branch, |v| then each quotient, wherever the port leaves Rust's code.
PC_FNS = ("field::electric_field_point_charge", "field::electric_potential_point_charge",
          "field::magnetic_field_moving_charge", "field::electric_field_superposition", "charge::coulomb_force",
          "PointCharge::electric_field_at", "PointCharge::electric_potential_at", "PointCharge::magnetic_field_at",
          "field::greens_function_static", "field::greens_function_scalar")
D4_FNS = ("field::electric_dipole_field", "field::magnetic_dipole_field", "field::electric_dipole_potential",
          "field::electric_quadrupole_potential", "field::electric_field_ring_axis", "field::electric_field_disk_axis",
          "field::image_charge_sphere")
FV_FNS = ("FieldVector::magnitude", "FieldVector::magnitude_sq", "FieldVector::dot", "FieldVector::cross",
          "FieldVector::normalized", "field::image_charge_plane", "field::image_charge_dielectric",
          "charge::lorentz_force", "PointCharge::lorentz_force_in", "field::trace_field_line") + PC_FNS + D4_FNS + (
    "field::maxwell_stress_tensor", "field::stress_tensor_force", "field::poynting_flux_surface",
    "field::multipole_moment", "field::spherical_multipole_moment")


VX_LO, VX_HI = 2.0 ** -900, sys.float_info.max


def vx_band(x) -> bool:
    """src/crmath_vx.cyr's band: |x| in [2^-900, DBL_MAX] (NaN: False)."""
    return VX_LO <= abs(x) <= VX_HI


def vx_exact(a, b) -> bool:
    p = a * b
    return a == 0.0 or b == 0.0 or MIN_NORMAL < abs(p) <= VX_HI


def vx_nonfinite(*xs) -> bool:
    return not all(math.isfinite(x) for x in xs)


def sumsq_plain(x, y, z) -> bool:
    s = x * x + y * y + z * z
    return vx_band(s) or (x == 0.0 and y == 0.0 and z == 0.0) or vx_nonfinite(x, y, z)


def cross1_plain(a, b, c, d) -> bool:
    r = a * b - c * d
    return vx_band(r) or vx_nonfinite(a, b, c, d) or (r == 0.0 and vx_exact(a, b) and vx_exact(c, d))


def lorentz1_plain(e, a, b, c, d, q) -> bool:
    s = e + (a * b - c * d)
    return (vx_band(s) or vx_nonfinite(q, e) or vx_nonfinite(a, b, c, d)
            or (s == 0.0 and vx_exact(a, b) and vx_exact(c, d)))


def fv_sumsq(x, y, z):
    return ux_add(ux_add(ux_mul(ux(x), ux(x)), ux_mul(ux(y), ux(y))), ux_mul(ux(z), ux(z)))


def fv_cross1(a, b, c, d, r=None):
    return ux_sub(ux_mul(ux(a), ux(b)), ux_mul(ux(c), ux(d)), r or rnd53)


def fv_lorentz(q, vel, e, b):
    ident = lambda t: t                                  # noqa: E731
    vb = [fv_cross1(vel[1], b[2], vel[2], b[1]), fv_cross1(vel[2], b[0], vel[0], b[2]),
          fv_cross1(vel[0], b[1], vel[1], b[0])]
    return ok3(tuple(ux_fin(ux_mul(ux_add(ux(e[i]), vb[i]), ux(q), ident)) for i in range(3)))


# M1 part 5d-3's chains: the point charges, coulomb_force and the Green's functions. pc_parts gives the inputs in the
# port's terms; pc_plain mirrors src/field.cyr's and src/charge.cyr's guards on Rust's binary64 values; pc_model is the
# exact chain rounded to 53 bits, the exact threshold first (row 17).
def pc_parts(fn, v):
    """(kind, charges [(q, pos)], second charge, field point, velocity, k) of a 5d-3 row."""
    if fn in ("field::electric_field_point_charge", "field::electric_potential_point_charge"):
        return v[0], v[1:4], v[4:7], None, None
    if fn == "field::magnetic_field_moving_charge":
        return v[0], v[1:4], v[7:10], v[4:7], None
    if fn in ("PointCharge::electric_field_at", "PointCharge::electric_potential_at", "PointCharge::magnetic_field_at"):
        return v[0], v[2:5], v[8:11], v[5:8], None
    if fn in ("field::greens_function_static", "field::greens_function_scalar"):
        return None, v[3:6], v[0:3], None, (v[6] if fn.endswith("scalar") else None)
    return None


def nrm(x) -> bool:
    """A biased exponent in [2, 2046]: 2^-1021 ≤ |x| ≤ DBL_MAX."""
    return 2.0 ** -1021 <= abs(x) <= VX_HI


def unit_ok(h, d) -> bool:
    return d == 0.0 or nrm(h)


def pc_rust(c, f):
    """Rust's offset, r², 1/r in binary64."""
    d = [f[i] - c[i] for i in range(3)]
    r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    rt = fsqrt(r_sq)
    return d, r_sq, rt, fdiv(1.0, rt)


def pc_plain(fn, v) -> bool:
    """The port takes Rust's chain for this (finite, non-singular) row."""
    c0 = golden_consts_cached()
    K, MU0 = c0["field::COULOMB_K"], c0["field::MU_0"]
    if fn == "charge::coulomb_force":
        q1, a, q2, b = v[0], v[2:5], v[8], v[10:13]
        d, r_sq, rt, inv = pc_rust(b, a)
        kq = K * q1
        kq2 = kq * q2
        fct = kq2 * inv * inv * inv
        return r_sq <= VX_HI and ((nrm(kq) and nrm(kq2) and nrm(fct)) or q1 == 0.0 or q2 == 0.0)
    if fn == "field::electric_field_superposition":
        cnt = int(len(v) - 3) // 4
        cs, at = v[0:4 * cnt], v[4 * cnt:4 * cnt + 3]
        tot = [0.0, 0.0, 0.0]
        ok = True
        for j in range(cnt):
            q, c = cs[4 * j], cs[4 * j + 1:4 * j + 4]
            d, r_sq, rt, inv = pc_rust(c, at)
            kq = K * q
            fct = kq * inv * inv * inv
            e = [fct * d[i] for i in range(3)]
            tot = [tot[i] + e[i] for i in range(3)]
            if not (r_sq <= VX_HI and (q == 0.0 or (nrm(kq) and nrm(fct) and all(unit_ok(e[i], d[i]) for i in range(3))))):
                ok = False
        return ok and all(math.isfinite(t) for t in tot)
    q, c, f, vel, k = pc_parts(fn, v)
    d, r_sq, rt, inv = pc_rust(c, f)
    if fn.startswith("field::greens_function"):
        den = 4.0 * PI * rt
        if fn.endswith("static"):
            return den <= VX_HI
        fct = fdiv(1.0, den)
        kr = k * rt
        return den <= VX_HI and nrm(fct) and (nrm(kr) or k == 0.0)
    if fn.endswith("magnetic_field_moving_charge") or fn.endswith("magnetic_field_at"):
        h = [d[i] * inv for i in range(3)]
        cr = [vel[1] * h[2] - vel[2] * h[1], vel[2] * h[0] - vel[0] * h[2], vel[0] * h[1] - vel[1] * h[0]]
        mq = fdiv(MU0, 4.0 * PI) * q
        fct = mq * inv * inv
        pairs = [(vel[1], h[2], vel[2], h[1]), (vel[2], h[0], vel[0], h[2]), (vel[0], h[1], vel[1], h[0])]
        crok = all(vx_band(cr[i]) or (cr[i] == 0.0 and (a == 0.0 or b == 0.0) and (x == 0.0 or y == 0.0))
                   for i, (a, b, x, y) in enumerate(pairs))
        return r_sq <= VX_HI and (q == 0.0 or (nrm(mq) and nrm(fct))) and all(unit_ok(h[i], d[i]) for i in range(3)) and crok
    kq = K * q
    if fn.endswith("potential_point_charge") or fn.endswith("potential_at"):
        return r_sq <= VX_HI and (q == 0.0 or nrm(kq))
    fct = kq * inv * inv * inv
    return r_sq <= VX_HI and (q == 0.0 or (nrm(kq) and nrm(fct)))


_CONSTS = {}


def golden_consts_cached():
    if not _CONSTS:
        _CONSTS.update(golden_consts(Path(__file__).resolve().parents[2]))
    return _CONSTS


ONE = (Fraction(1), False)


def pc_offset_ux(c, f):
    d = [ux_sub(ux(f[i]), ux(c[i])) for i in range(3)]
    r2 = ux_add(ux_add(ux_mul(d[0], d[0]), ux_mul(d[1], d[1])), ux_mul(d[2], d[2]))
    rt = ux_sqrt(r2)
    return d, rt, ux_div(ONE, rt)


def pc_factor_ux(K, q1, q2, inv):
    m = ux_mul(ux_mul(ux(K), ux(q1)), ux(q2))
    for _ in range(3):
        m = ux_mul(m, inv)
    return m


def pc_model(fn, v):
    """The port's outputs for a 5d-3 row with finite inputs: the exact threshold, then the unbounded chain."""
    ident = lambda t: t                                  # noqa: E731
    c0 = golden_consts_cached()
    K, MU0 = c0["field::COULOMB_K"], c0["field::MU_0"]
    if fn == "charge::coulomb_force":
        q1, a, q2, b = v[0], v[2:5], v[8], v[10:13]
        if exact_sq(a, b) < Fraction(1, 10 ** 30):
            return SINGULAR
        d, rt, inv = pc_offset_ux(b, a)
        fct = pc_factor_ux(K, q1, q2, inv)
        return ok3(tuple(ux_fin(ux_mul(fct, d[i], ident)) for i in range(3)))
    if fn == "field::electric_field_superposition":
        cnt = int(len(v) - 3) // 4
        cs, at = v[0:4 * cnt], v[4 * cnt:4 * cnt + 3]
        acc = [(Fraction(0), False)] * 3
        fin = [0.0, 0.0, 0.0]
        for j in range(cnt):
            q, c = cs[4 * j], cs[4 * j + 1:4 * j + 4]
            if exact_sq(c, at) < Fraction(1, 10 ** 30):
                return SINGULAR
            d, rt, inv = pc_offset_ux(c, at)
            fct = pc_factor_ux(K, q, 1.0, inv)
            for i in range(3):
                g = ux_mul(fct, d[i])
                if j == cnt - 1:
                    fin[i] = ux_fin(ux_add(acc[i], g, ident))
                else:
                    acc[i] = ux_add(acc[i], g)
        return ok3(tuple(fin))
    q, c, f, vel, k = pc_parts(fn, v)
    if fn.startswith("field::greens_function"):
        if exact_sq(c, f) < Fraction(1, 10 ** 60):
            return SINGULAR
        d, rt, inv = pc_offset_ux(c, f)
        den = ux_mul(ux(4.0 * PI), rt)
        if fn.endswith("static"):
            return [("f", ux_fin(ux_div(ONE, den, ident)))]
        fct = ux_div(ONE, den)
        kr = ux_mul(ux(k), rt)
        return [("f", ux_fin(ux_mul(fct, w_trig("cos", kr), ident))),
                ("f", ux_fin(ux_mul(fct, w_trig("sin", kr), ident)))]
    if exact_sq(c, f) < Fraction(1, 10 ** 30):
        return SINGULAR
    d, rt, inv = pc_offset_ux(c, f)
    if fn.endswith("magnetic_field_moving_charge") or fn.endswith("magnetic_field_at"):
        h = [ux_mul(d[i], inv) for i in range(3)]
        cr = [ux_sub(ux_mul(ux(vel[1]), h[2]), ux_mul(ux(vel[2]), h[1])),
              ux_sub(ux_mul(ux(vel[2]), h[0]), ux_mul(ux(vel[0]), h[2])),
              ux_sub(ux_mul(ux(vel[0]), h[1]), ux_mul(ux(vel[1]), h[0]))]
        fct = ux_mul(ux_mul(ux_mul(ux(fdiv(MU0, 4.0 * PI)), ux(q)), inv), inv)
        return ok3(tuple(ux_fin(ux_mul(fct, cr[i], ident)) for i in range(3)))
    kq = ux_mul(ux(K), ux(q))
    if fn.endswith("potential_point_charge") or fn.endswith("potential_at"):
        return [("f", ux_fin(ux_div(kq, rt, ident)))]
    fct = pc_factor_ux(K, q, 1.0, inv)
    return ok3(tuple(ux_fin(ux_mul(fct, d[i], ident)) for i in range(3)))


# M1 part 5d-4's chains: the dipoles, the quadrupole, the ring and disk axes, the sphere image. d4_plain mirrors
# src/field.cyr's guards on Rust's binary64 values (rng: a biased exponent in [lo, 2046]); d4_model is the exact chain.
def rng(x, lo) -> bool:
    return (2.0 ** -1021 if lo == 2 else 2.0 ** -900) <= abs(x) <= VX_HI


def d4_rust_dipole(M, W, magnetic, EPS0, MU0):
    r_sq = W[0] * W[0] + W[1] * W[1] + W[2] * W[2]
    inv = fdiv(1.0, fsqrt(r_sq))
    h = [W[i] * inv for i in range(3)]
    pd = M[0] * h[0] + M[1] * h[1] + M[2] * h[2]
    inv3 = inv * inv * inv
    fct = fdiv(MU0 * inv3, 4.0 * PI) if magnetic else fdiv(inv3, 4.0 * PI * EPS0)
    return inv, h, pd, inv3, fct


def d4_plain(fn, v) -> bool:
    c0 = golden_consts_cached()
    K, EPS0, MU0 = c0["field::COULOMB_K"], c0["field::EPSILON_0"], c0["field::MU_0"]
    if fn in ("field::electric_dipole_field", "field::magnetic_dipole_field", "field::electric_dipole_potential"):
        M, W = v[0:3], v[3:6]
        inv, h, pd, inv3, fct = d4_rust_dipole(M, W, fn == "field::magnetic_dipole_field", EPS0, MU0)
        zd = all(M[i] == 0.0 or W[i] == 0.0 for i in range(3))
        units = all(rng(h[i], 2) or W[i] == 0.0 for i in range(3)) and (rng(pd, 123) or zd)
        if fn == "field::electric_dipole_potential":
            t2 = pd * inv * inv
            return units and (rng(t2, 2) or zd)
        t3 = 3.0 * pd
        k = t3 * fct
        e = [h[i] * k - M[i] * fct for i in range(3)]
        return (units and rng(inv3, 2) and rng(fct, 2) and (rng(t3, 2) or zd) and (rng(k, 2) or zd)
                and all(rng(e[i], 123) or ((W[i] == 0.0 or zd) and M[i] == 0.0) for i in range(3)))
    if fn == "field::electric_quadrupole_potential":
        qt, p = v[0:9], v[9:12]
        r_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
        r = fsqrt(r_sq)
        r5 = r_sq * r_sq * r
        d = 8.0 * PI * EPS0 * r5
        tot = 0.0
        ok = rng(r_sq, 2) and rng(r5, 2) and rng(d, 2)
        for i in range(3):
            for j in range(3):
                a = qt[3 * i + j] * p[i]
                t = a * p[j]
                tot += t
                if qt[3 * i + j] != 0.0 and p[i] != 0.0 and not (rng(a, 2) and (p[j] == 0.0 or rng(t, 2))):
                    ok = False
        return ok and math.isfinite(tot)
    q, rr, z = v[0], v[1], v[2]
    if fn == "field::electric_field_ring_axis":
        if not (rr >= 0.0 and math.isfinite(rr)) or not math.isfinite(z) or z == 0.0:
            return True
        d_sq = z * z + rr * rr
        denom = d_sq * fsqrt(d_sq)
        kq = K * q
        return rng(d_sq, 123) and rng(denom, 2) and (q == 0.0 or (rng(kq, 2) and rng(kq * z, 2)))
    if fn == "field::electric_field_disk_axis":
        if not (rr > 0.0 and math.isfinite(rr)) or not math.isfinite(z):
            return True
        s2 = z * z + rr * rr
        r = fsqrt(s2)
        g1, rz = fdiv(rr, r), r + abs(z)
        g2 = fdiv(rr, rz)
        sd = fdiv(math.copysign(1.0, z) * q, 2.0 * EPS0)
        return (rng(s2, 123) and rng(g1, 2) and rz <= VX_HI and rng(g2, 2) and rng(g1 * g2, 2)
                and (q == 0.0 or rng(sd, 2)))
    q, p, rr = v[0], v[1:4], v[4]
    if not (rr > 0.0 and math.isfinite(rr)):
        return True
    d_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
    if not rng(d_sq, 123):
        return False
    d = fsqrt(d_sq)
    if d <= rr:
        return True
    ratio = fdiv(rr, d)
    return rng(ratio, 2) and rng(ratio * ratio, 2)


def d4_model(fn, v):
    """The port's outputs for a 5d-4 row with finite inputs."""
    ident = lambda t: t                                  # noqa: E731
    c0 = golden_consts_cached()
    K, EPS0, MU0 = c0["field::COULOMB_K"], c0["field::EPSILON_0"], c0["field::MU_0"]
    if fn in ("field::electric_dipole_field", "field::magnetic_dipole_field", "field::electric_dipole_potential"):
        M, W = v[0:3], v[3:6]
        if exact_sq(W, ORIGIN) < Fraction(1, 10 ** 30):
            return SINGULAR
        d, rt, inv = pc_offset_ux(ORIGIN, W)
        h = [ux_mul(d[i], inv) for i in range(3)]
        pd = ux_add(ux_add(ux_mul(ux(M[0]), h[0]), ux_mul(ux(M[1]), h[1])), ux_mul(ux(M[2]), h[2]))
        if fn == "field::electric_dipole_potential":
            t2 = ux_mul(ux_mul(pd, inv), inv)
            return [("f", ux_fin(ux_div(t2, ux(4.0 * PI * EPS0), ident)))]
        inv3 = ux_mul(ux_mul(inv, inv), inv)
        if fn == "field::magnetic_dipole_field":
            fct = ux_div(ux_mul(ux(MU0), inv3), ux(4.0 * PI))
        else:
            fct = ux_div(inv3, ux(4.0 * PI * EPS0))
        k = ux_mul(ux_mul(ux(3.0), pd), fct)
        return ok3(tuple(ux_fin(ux_sub(ux_mul(h[i], k), ux_mul(ux(M[i]), fct), ident)) for i in range(3)))
    if fn == "field::electric_quadrupole_potential":
        qt, p = v[0:9], v[9:12]
        if exact_sq(p, ORIGIN) < Fraction(1, 10 ** 30):
            return SINGULAR
        r2 = fv_sumsq(*p)
        r = ux_sqrt(r2)
        r5 = ux_mul(ux_mul(r2, r2), r)
        d = ux_mul(ux(8.0 * PI * EPS0), r5)
        acc = (Fraction(0), False)
        for i in range(3):
            for j in range(3):
                acc = ux_add(acc, ux_mul(ux_mul(ux(qt[3 * i + j]), ux(p[i])), ux(p[j])))
        return [("f", ux_fin(ux_div(acc, d, ident)))]
    if fn in ("field::electric_field_ring_axis", "field::electric_field_disk_axis"):
        q, rr, z = v[0], v[1], v[2]
        w = field_3b(fn, v, None, MPFR_P, K, EPS0, MU0)
        if w[0][0] == "E" or (fn.endswith("ring_axis") and z == 0.0):
            return w
        r2 = ux_add(ux_mul(ux(z), ux(z)), ux_mul(ux(rr), ux(rr)))
        r = ux_sqrt(r2)
        if fn.endswith("ring_axis"):
            num = ux_mul(ux_mul(ux(K), ux(q)), ux(z))
            return [("f", ux_fin(ux_div(num, ux_mul(r2, r), ident)))]
        geo = ux_mul(ux_div(ux(rr), r), ux_div(ux(rr), ux_add(r, ux(abs(z)))))
        sd = ux_div(ux(math.copysign(1.0, z) * q), ux(2.0 * EPS0))
        return [("f", ux_fin(ux_mul(sd, geo, ident)))]
    q, p, rr = v[0], v[1:4], v[4]
    w = field_3b(fn, v, None, MPFR_P, K, EPS0, MU0)
    if not (rr > 0.0 and math.isfinite(rr)):
        return w
    d = ux_sqrt(fv_sumsq(*p))
    if d[0] <= Fraction(rr):
        dd = ux_fin(d)
        return ip_err(f"charge must be outside the sphere (|r| = {rust_display(dd)} > R = {rust_display(rr)})")
    ratio = ux_div(ux(rr), d)
    sc = ux_mul(ratio, ratio)
    return ([("f", ux_fin(ux_mul(ux(-q), ratio, ident)))]
            + [("f", ux_fin(ux_mul(ux(p[i]), sc, ident))) for i in range(3)])


# M1 part 5d-5's chains: the stress tensor, the force from it, the Poynting flux, the multipoles and a field line's
# steps. m5_plain mirrors src/field.cyr's guards on Rust's binary64 values; m5_model is the whole unbounded chain.
M5_FNS = ("field::maxwell_stress_tensor", "field::stress_tensor_force", "field::poynting_flux_surface",
          "field::multipole_moment", "field::spherical_multipole_moment")
HALF = (Fraction(1, 2), False)
TRACE_FN = None                                         # validate_field's trace, registered by each run


def chk(r, a, b) -> bool:
    """src/field.cyr's _bj_fv_chk: the product or quotient r = a ⊗ b is normal above 2^-1022, or 0 with a zero
    operand."""
    if r == 0.0:
        return a == 0.0 or b == 0.0
    return nrm(r)


def ux_neg(a):
    return -a[0], not a[1]


def ux_sum(terms):
    """Σ terms from 0 in order, each addition to 53 bits but the last, the result rounded once (+0 for none)."""
    if not terms:
        return 0.0
    acc = (Fraction(0), False)
    for t in terms[:-1]:
        acc = ux_add(acc, t)
    return ux_fin(ux_add(acc, terms[-1], lambda t: t))


def m5_box(xs, lo, hi) -> bool:
    """_bj_fv_box: each part ±0 or with its biased exponent in [lo, hi)."""
    for x in xs:
        a = f_to_bits(x) & 0x7FFFFFFFFFFFFFFF
        if a != 0 and not (lo <= a >> 52 < hi):
            return False
    return True


def m5_consts():
    c0 = golden_consts_cached()
    return c0["field::EPSILON_0"], c0["field::MU_0"]


def m5_coords_k(p, ok):
    """_bj_spherical_coords on Rust's values: (r, cos θ, sin θ, φ) or None at the origin; clears ok[0]."""
    x, y, z = p
    xx, yy, zz = x * x, y * y, z * z
    ok[0] &= chk(xx, x, x) and chk(yy, y, y) and chk(zz, z, z)
    rho_sq = xx + yy
    sq = rho_sq + zz
    r = fsqrt(sq)
    if not math.isfinite(sq):
        ok[0] = False
    if r == 0.0:
        if x != 0.0 or y != 0.0 or z != 0.0:
            ok[0] = False
        return None
    c0 = fdiv(z, r)
    rho = fsqrt(rho_sq)
    s0 = fdiv(rho, r)
    ok[0] &= chk(c0, z, r) and chk(s0, rho, r)
    g = MPFR_P["atan2"](y, x)
    if y != 0.0 and abs(g) < 2.0 ** -1021:
        ok[0] = False
    return r, rust_clamp(c0, -1.0, 1.0), rust_clamp(s0, 0.0, 1.0), g


def m5_rpow(r, ll):
    return powidf2_ieee(r, ll) if ll <= 2147483647 else MPFR_P["pow"](r, float(ll))


def m5_legendre_k(ll, m, x, s, ok):
    if m > ll:
        return 0.0
    pmm, fact = 1.0, 1.0
    for _ in range(m):
        t = -fact * s
        p2 = pmm * t
        ok[0] &= chk(t, -fact, s) and chk(p2, pmm, t)
        pmm = p2
        fact += 2.0
    if ll == m:
        return pmm
    mf = float(m)
    c1 = 2.0 * mf + 1.0
    u = x * c1
    pmmp1 = u * pmm
    ok[0] &= chk(u, x, c1) and chk(pmmp1, u, pmm)
    if ll == m + 1:
        return pmmp1
    pll = 0.0
    for k in range(m + 2, ll + 1):
        kf = float(k)
        c2 = 2.0 * kf - 1.0
        u2 = x * c2
        a = u2 * pmmp1
        c3 = kf + mf - 1.0
        b = c3 * pmm
        c4 = kf - mf
        dd = a - b
        pll = fdiv(dd, c4)
        ok[0] &= chk(u2, x, c2) and chk(a, u2, pmmp1) and chk(b, c3, pmm) and chk(pll, dd, c4)
        pmm, pmmp1 = pmmp1, pll
    return pll


def m5_nlegendre_k(ll, m, x, s, ok):
    if m > ll:
        return 0.0
    pmm = 1.0
    for k in range(1, m + 1):
        kf = float(k)
        sr = fsqrt(fdiv(2.0 * kf - 1.0, 2.0 * kf))
        t = -s * sr
        p2 = pmm * t
        ok[0] &= chk(t, -s, sr) and chk(p2, pmm, t)
        pmm = p2
    mf = float(m)
    c = fsqrt(fdiv(2.0 * mf + 1.0, 4.0 * PI))
    p3 = pmm * c
    ok[0] &= chk(p3, pmm, c)
    pmm = p3
    if ll == m:
        return pmm
    s3 = fsqrt(2.0 * mf + 3.0)
    u = x * s3
    pmmp1 = u * pmm
    ok[0] &= chk(u, x, s3) and chk(pmmp1, u, pmm)
    if ll == m + 1:
        return pmmp1
    m_sq = mf * mf
    a_prev = fsqrt(fdiv(4.0 * (mf + 1.0) * (mf + 1.0) - 1.0, (mf + 1.0) * (mf + 1.0) - m_sq))
    pll = 0.0
    for k in range(m + 2, ll + 1):
        kf = float(k)
        a = fsqrt(fdiv(4.0 * kf * kf - 1.0, kf * kf - m_sq))
        xa = x * pmmp1
        pb = fdiv(pmm, a_prev)
        inner = xa - pb
        pll = a * inner
        ok[0] &= chk(xa, x, pmmp1) and chk(pb, pmm, a_prev) and chk(pll, a, inner)
        pmm, pmmp1, a_prev = pmmp1, pll, a
    return pll


def m5_charges(v):
    cnt = int(v[2]) // 4
    return v[0], v[1], v[3:3 + 4 * cnt]


def m5_plain(fn, v) -> bool:
    """The port keeps Rust's chain for this row (finite parts, no error)."""
    EPS0, MU0 = m5_consts()
    inv_mu0 = fdiv(1.0, MU0)
    if fn == "field::maxwell_stress_tensor":
        return m5_box(v[0:6], 773, 1273)
    if fn == "field::stress_tensor_force":
        tt, n = v[0:9], v[9:12]
        for i in range(3):
            acc = 0.0
            for j in range(3):
                acc += tt[3 * i + j] * n[j]
            if not (vx_band(acc) or (acc == 0.0 and all(vx_exact(tt[3 * i + j], n[j]) for j in range(3)))):
                return False
        return True
    if fn == "field::poynting_flux_surface":
        cnt = int(v[0]) // 10
        sm = v[1:1 + 10 * cnt]
        total, ok, tz = 0.0, True, True
        for j in range(cnt):
            p = sm[10 * j:10 * j + 10]
            dt = 0.0
            for k in range(3):
                j1, j2 = (k + 1) % 3, (k + 2) % 3
                ck = p[j1] * p[3 + j2] - p[j2] * p[3 + j1]
                if not (vx_band(ck) or (ck == 0.0 and vx_exact(p[j1], p[3 + j2]) and vx_exact(p[j2], p[3 + j1]))):
                    ok = False
                sk = ck * inv_mu0
                dt = sk * p[6 + k] if k == 0 else dt + sk * p[6 + k]
            if not vx_band(dt):
                if dt != 0.0:
                    ok = False
                else:
                    for k in range(3):
                        j1, j2 = (k + 1) % 3, (k + 2) % 3
                        sk = (p[j1] * p[3 + j2] - p[j2] * p[3 + j1]) * inv_mu0
                        if not vx_exact(sk, p[6 + k]):
                            ok = False
            tm = dt * p[9]
            tz = tz and vx_exact(dt, p[9])
            total += tm
        return ok and (vx_band(total) or (total == 0.0 and tz))
    ll, m, cs = m5_charges(v)
    ma = abs(m)
    mf = float(ma)
    ok = [True]
    if fn == "field::multipole_moment":
        mom = 0.0
        for j in range(len(cs) // 4):
            q = cs[4 * j]
            c = m5_coords_k(cs[4 * j + 1:4 * j + 4], ok)
            if c is None:
                if ll == 0:
                    mom += q
                continue
            ang = mf * c[3]
            ok[0] &= chk(ang, mf, c[3])
            az = MPFR_P["cos"](ang) if m >= 0 else MPFR_P["sin"](ang)
            rp = m5_rpow(c[0], ll)
            if ll != 0 and not nrm(rp):
                ok[0] = False
            pl = m5_legendre_k(ll, ma, c[1], c[2], ok)
            t1 = q * rp
            t2 = t1 * pl
            term = t2 * az
            ok[0] &= chk(t1, q, rp) and chk(t2, t1, pl) and chk(term, t2, az)
            mom += term
        return ok[0] and math.isfinite(mom)
    neg = 1.0 if ma % 2 == 0 else -1.0
    sp = fsqrt(PI)
    re = im = 0.0
    for j in range(len(cs) // 4):
        q = cs[4 * j]
        c = m5_coords_k(cs[4 * j + 1:4 * j + 4], ok)
        if c is None:
            if ll == 0:
                qh = q * 0.5
                y00 = fdiv(qh, sp)
                ok[0] &= chk(qh, q, 0.5) and chk(y00, qh, sp)
                re += y00
            continue
        rp = m5_rpow(c[0], ll)
        if ll != 0 and not nrm(rp):
            ok[0] = False
        pl = m5_nlegendre_k(ll, ma, c[1], c[2], ok)
        t1 = q * rp
        w = t1 * pl
        ang = mf * c[3]
        ok[0] &= chk(t1, q, rp) and chk(w, t1, pl) and chk(ang, mf, c[3])
        sv, cv = MPFR_P["sin"](ang), MPFR_P["cos"](ang)
        if m >= 0:
            wc, ws = w * cv, w * sv
            ok[0] &= chk(wc, w, cv) and chk(ws, w, sv)
            re += wc
            im -= ws
        else:
            g = neg * w
            nc, ns = g * cv, g * sv
            ok[0] &= chk(nc, g, cv) and chk(ns, g, sv)
            re += nc
            im += ns
    return ok[0] and math.isfinite(re) and math.isfinite(im)


def m5_clamp(a, lo, hi):
    if a[0] < lo[0]:
        return lo
    if a[0] > hi[0]:
        return hi
    return a


def m5_coords_ux(p):
    """_bj_ux_sph: (r, cos θ, sin θ, φ) as model pairs, None at the origin."""
    x, y, z = p
    if x == 0.0 and y == 0.0 and z == 0.0:
        return None
    r2 = ux_add(ux_mul(ux(x), ux(x)), ux_mul(ux(y), ux(y)))
    r = ux_sqrt(ux_add(r2, ux_mul(ux(z), ux(z))))
    c = m5_clamp(ux_div(ux(z), r), (Fraction(-1), True), (Fraction(1), False))
    sn = ux_div(ux_sqrt(r2), r)
    if sn[0] > 1:
        sn = (Fraction(1), False)
    g = MPFR_P["atan2"](y, x)
    phi = ux(g) if (y == 0.0 or abs(g) > 2.0 ** -1022) else ux_div(ux(y), ux(x))
    return r, c, sn, phi


def m5_rpow_ux(r, ll):
    if ll > 2147483647:
        return wv_pow(r, float(ll))
    a, acc, p = r, (Fraction(1), False), ll
    while True:
        if p & 1:
            acc = ux_mul(acc, a)
        p >>= 1
        if p == 0:
            break
        a = ux_mul(a, a)
    return acc


def m5_legendre_ux(ll, m, x, s):
    if m > ll:
        return Fraction(0), False
    pmm, fact = (Fraction(1), False), 1.0
    for _ in range(m):
        pmm = ux_mul(pmm, ux_mul(ux(-fact), s))
        fact += 2.0
    if ll == m:
        return pmm
    mf = float(m)
    pmmp1 = ux_mul(ux_mul(x, ux(2.0 * mf + 1.0)), pmm)
    if ll == m + 1:
        return pmmp1
    pll = (Fraction(0), False)
    for k in range(m + 2, ll + 1):
        kf = float(k)
        a = ux_mul(ux_mul(x, ux(2.0 * kf - 1.0)), pmmp1)
        b = ux_mul(ux(kf + mf - 1.0), pmm)
        pll = ux_div(ux_sub(a, b), ux(kf - mf))
        pmm, pmmp1 = pmmp1, pll
    return pll


def m5_nlegendre_ux(ll, m, x, s):
    if m > ll:
        return Fraction(0), False
    pmm = (Fraction(1), False)
    for k in range(1, m + 1):
        kf = float(k)
        pmm = ux_mul(pmm, ux_mul(ux_neg(s), ux(fsqrt(fdiv(2.0 * kf - 1.0, 2.0 * kf)))))
    mf = float(m)
    pmm = ux_mul(pmm, ux(fsqrt(fdiv(2.0 * mf + 1.0, 4.0 * PI))))
    if ll == m:
        return pmm
    pmmp1 = ux_mul(ux_mul(x, ux(fsqrt(2.0 * mf + 3.0))), pmm)
    if ll == m + 1:
        return pmmp1
    m_sq = mf * mf
    a_prev = fsqrt(fdiv(4.0 * (mf + 1.0) * (mf + 1.0) - 1.0, (mf + 1.0) * (mf + 1.0) - m_sq))
    pll = (Fraction(0), False)
    for k in range(m + 2, ll + 1):
        kf = float(k)
        a = fsqrt(fdiv(4.0 * kf * kf - 1.0, kf * kf - m_sq))
        inner = ux_sub(ux_mul(x, pmmp1), ux_div(pmm, ux(a_prev)))
        pll = ux_mul(ux(a), inner)
        pmm, pmmp1, a_prev = pmmp1, pll, a
    return pll


def m5_model(fn, v):
    """The port's outputs under row 18 (finite parts, no error)."""
    ident = lambda t: t                                  # noqa: E731
    EPS0, MU0 = m5_consts()
    inv_mu0 = fdiv(1.0, MU0)
    if fn == "field::maxwell_stress_tensor":
        e, b = v[0:3], v[3:6]
        he, hb = ux_mul(HALF, fv_sumsq(*e)), ux_mul(HALF, fv_sumsq(*b))
        out = []
        for i in range(3):
            for j in range(3):
                p, qq = ux_mul(ux(e[i]), ux(e[j])), ux_mul(ux(b[i]), ux(b[j]))
                t = ux_sub(p, he if i == j else ux(0.0))
                u = ux_sub(qq, hb if i == j else ux(0.0))
                out.append(("f", ux_fin(ux_add(ux_mul(ux(EPS0), t), ux_mul(ux(inv_mu0), u), ident))))
        return out
    if fn == "field::stress_tensor_force":
        tt, n = v[0:9], v[9:12]
        return ok3(tuple(ux_sum([ux_mul(ux(tt[3 * i + j]), ux(n[j])) for j in range(3)]) for i in range(3)))
    if fn == "field::poynting_flux_surface":
        cnt = int(v[0]) // 10
        sm = v[1:1 + 10 * cnt]
        terms = []
        for j in range(cnt):
            p = sm[10 * j:10 * j + 10]
            dt = None
            for k in range(3):
                j1, j2 = (k + 1) % 3, (k + 2) % 3
                ck = fv_cross1(p[j1], p[3 + j2], p[j2], p[3 + j1])
                dk = ux_mul(ux_mul(ck, ux(inv_mu0)), ux(p[6 + k]))
                dt = dk if k == 0 else ux_add(dt, dk)
            terms.append(ux_mul(dt, ux(p[9])))
        return [("f", ux_sum(terms))]
    ll, m, cs = m5_charges(v)
    ma = abs(m)
    mf = ux(float(ma))
    if fn == "field::multipole_moment":
        terms = []
        for j in range(len(cs) // 4):
            q = ux(cs[4 * j])
            c = m5_coords_ux(cs[4 * j + 1:4 * j + 4])
            if c is None:
                if ll == 0:
                    terms.append(q)
                continue
            az = w_trig("cos" if m >= 0 else "sin", ux_mul(mf, c[3]))
            pl = m5_legendre_ux(ll, ma, c[1], c[2])
            terms.append(ux_mul(ux_mul(ux_mul(q, m5_rpow_ux(c[0], ll)), pl), az))
        return [("f", ux_sum(terms))]
    neg = ux(1.0 if ma % 2 == 0 else -1.0)
    sp = ux(fsqrt(PI))
    tr, ti = [], []
    for j in range(len(cs) // 4):
        q = ux(cs[4 * j])
        c = m5_coords_ux(cs[4 * j + 1:4 * j + 4])
        if c is None:
            if ll == 0:
                tr.append(ux_div(ux_mul(q, HALF), sp))
            continue
        w = ux_mul(ux_mul(q, m5_rpow_ux(c[0], ll)), m5_nlegendre_ux(ll, ma, c[1], c[2]))
        ang = ux_mul(mf, c[3])
        sv, cv = w_trig("sin", ang), w_trig("cos", ang)
        if m >= 0:
            tr.append(ux_mul(w, cv))
            ti.append(ux_neg(ux_mul(w, sv)))
        else:
            g = ux_mul(neg, w)
            tr.append(ux_mul(g, cv))
            ti.append(ux_mul(g, sv))
    return [("f", ux_sum(tr)), ("f", ux_sum(ti))]


def m5_step(pos, f, h):
    """One field-line step with an unbounded exponent: pos + h·F/|F|, each new position rounded once."""
    mg = ux_sqrt(fv_sumsq(*f))
    return tuple(ux_fin(ux_add(ux(pos[i]), ux_mul(ux(h), ux_div(ux(f[i]), mg)), lambda t: t)) for i in range(3))


def m5_step_plain(f, h) -> bool:
    """The trace's step guard on Rust's values: |F|² in the band, h < 2^1023, and each part F_i = 0 or v_i normal above
    2^-1022, v_i = h·F_i/|F| for h ≤ 1 (then |F_i/|F|| ≥ |v_i|) and F_i/|F| for h > 1 (then |h·F_i/|F|| ≥ it)."""
    ms = f[0] * f[0] + f[1] * f[1] + f[2] * f[2]
    if not vx_band(ms) or h >= 2.0 ** 1023:
        return False
    mag = math.sqrt(ms)
    return all(f[i] == 0.0 or 2.0 ** -1021 <= abs(h * (f[i] / mag) if h <= 1.0 else f[i] / mag) for i in range(3))


def m5_trace_field(sid, p, pos, rust):
    """The port's field for a trace scenario: the field function's row-18 model where its guard leaves Rust's chain,
    else Rust's tokens (`rust`); and whether it was governed."""
    if sid == 1:
        fn, vv = "field::electric_field_point_charge", [p[0], p[1], p[2], p[3], *pos]
    elif sid in (2, 8):
        fn = "field::magnetic_dipole_field" if sid == 8 else "field::electric_dipole_field"
        vv = [p[0], p[1], p[2], *pos]
    elif sid == 3:
        fn = "field::electric_field_superposition"
        vv = [p[0], 0.0, 0.0, p[1], p[0], 0.0, 0.0, -p[1], -2.0 * p[0], 0.0, 0.0, 0.0, *pos]
    else:
        return rust, False
    if not all(math.isfinite(x) for x in vv) or fv_governed(fn, vv) == 0:
        return rust, False
    return fv_model(fn, vv), True


def fv_args(fn, v):
    """The parts a row-18 function's chain reads (None for the images)."""
    if fn in D4_FNS:
        return list(v[0:12] if fn.endswith("quadrupole_potential") else v[0:6] if "dipole" in fn else
                    v[0:5] if fn.endswith("sphere") else v[0:3])
    if fn == "charge::coulomb_force":
        return [v[0], v[8]] + v[2:5] + v[10:13]
    if fn == "field::electric_field_superposition":
        return list(v)
    if fn in PC_FNS:
        q, c, f, vel, k = pc_parts(fn, v)
        return [x for x in [q, k] if x is not None] + list(c) + list(f) + (list(vel) if vel is not None else [])
    if fn in ("FieldVector::magnitude", "FieldVector::magnitude_sq", "FieldVector::normalized"):
        return v[0:3]
    if fn in ("FieldVector::dot", "FieldVector::cross"):
        return v[0:6]
    if fn == "charge::lorentz_force":
        return [v[0]] + v[1:10]
    if fn == "PointCharge::lorentz_force_in":
        return [v[0]] + v[5:8] + v[8:14]
    if fn == "field::maxwell_stress_tensor":
        return v[0:6]
    if fn == "field::stress_tensor_force":
        return v[0:12]
    if fn == "field::poynting_flux_surface":
        return v[1:1 + int(v[0])]
    if fn in ("field::multipole_moment", "field::spherical_multipole_moment"):
        return m5_charges(v)[2]
    return None


def fv_governed(fn, v) -> int:
    """18 where the port takes the unbounded chain for this row (any component), else 0 (Rust's chain): the guards of
    src/field.cyr and src/charge.cyr, on Rust's binary64 values."""
    if fn in PC_FNS:
        return 0 if (pc_model(fn, v) is SINGULAR or pc_plain(fn, v)) else 18
    if fn in D4_FNS:
        return 0 if (d4_model(fn, v) is SINGULAR or d4_plain(fn, v)) else 18
    if fn in M5_FNS:
        return 0 if m5_plain(fn, v) else 18
    if fn == "field::image_charge_plane":
        pz, z = v[4], v[3]
        return 18 if math.isfinite(pz) and abs(pz) >= 2.0 ** 1023 and math.isfinite(z) else 0
    if fn == "field::image_charge_dielectric":
        q, p, e1, e2 = v[0], v[1:4], v[4], v[5]
        if not (p[2] > 0.0 and math.isfinite(p[2]) and all(x > 0.0 and math.isfinite(x) for x in (e1, e2))):
            return 0
        if not math.isfinite(q):
            return 0
        nq, diff = -q, e2 - e1
        num, den = nq * diff, e1 + e2
        if not math.isfinite(den):
            return 18
        if num == 0.0:
            return 0 if (nq == 0.0 or diff == 0.0) else 18
        return 0 if MIN_NORMAL < abs(num) < math.inf else 18
    if fn in ("FieldVector::magnitude", "FieldVector::magnitude_sq", "FieldVector::normalized"):
        return 0 if sumsq_plain(*v[0:3]) else 18
    if fn == "FieldVector::dot":
        a, b = v[0:3], v[3:6]
        s = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        ok = vx_band(s) or vx_nonfinite(*a, *b) or (s == 0.0 and all(vx_exact(a[i], b[i]) for i in range(3)))
        return 0 if ok else 18
    if fn == "FieldVector::cross":
        a, b = v[0:3], v[3:6]
        ok = (cross1_plain(a[1], b[2], a[2], b[1]) and cross1_plain(a[2], b[0], a[0], b[2])
              and cross1_plain(a[0], b[1], a[1], b[0]))
        return 0 if ok else 18
    q, vel, e, b = (v[0], v[1:4], v[4:7], v[7:10]) if fn == "charge::lorentz_force" else (v[0], v[5:8], v[8:11],
                                                                                         v[11:14])
    ok = (lorentz1_plain(e[0], vel[1], b[2], vel[2], b[1], q) and lorentz1_plain(e[1], vel[2], b[0], vel[0], b[2], q)
          and lorentz1_plain(e[2], vel[0], b[1], vel[1], b[0], q))
    return 0 if ok else 18


def fv_model(fn, v):
    """The port's outputs under row 18 for a governed row (finite parts)."""
    ident = lambda t: t                                  # noqa: E731
    if fn in PC_FNS:
        return pc_model(fn, v)
    if fn in D4_FNS:
        return d4_model(fn, v)
    if fn in M5_FNS:
        return m5_model(fn, v)
    if fn == "FieldVector::magnitude":
        return [("f", ux_fin_sqrt(fv_sumsq(*v[0:3])))]
    if fn == "FieldVector::magnitude_sq":
        x, y, z = v[0:3]
        return [("f", ux_fin(ux_add(ux_add(ux_mul(ux(x), ux(x)), ux_mul(ux(y), ux(y))), ux_mul(ux(z), ux(z)),
                                    ident)))]
    if fn == "FieldVector::dot":
        a, b = v[0:3], v[3:6]
        return [("f", ux_fin(ux_add(ux_add(ux_mul(ux(a[0]), ux(b[0])), ux_mul(ux(a[1]), ux(b[1]))),
                                    ux_mul(ux(a[2]), ux(b[2])), ident)))]
    if fn == "FieldVector::cross":
        a, b = v[0:3], v[3:6]
        return ok3((ux_fin(fv_cross1(a[1], b[2], a[2], b[1], ident)), ux_fin(fv_cross1(a[2], b[0], a[0], b[2], ident)),
                    ux_fin(fv_cross1(a[0], b[1], a[1], b[0], ident))))
    if fn == "FieldVector::normalized":
        m = ux_sqrt(fv_sumsq(*v[0:3]))
        return ok3(tuple(ux_fin(ux_div(ux(c), m, ident)) for c in v[0:3]))
    if fn == "field::image_charge_plane":
        q, p, pz = v[0], v[1:4], v[4]
        two = ux(pz)
        iz = ux_fin(ux_sub((2 * two[0], two[1]), ux(p[2]), ident))
        return [("f", -q), ("f", p[0]), ("f", p[1]), ("f", iz)]
    if fn == "field::image_charge_dielectric":
        q, p, e1, e2 = v[0], v[1:4], v[4], v[5]
        diff = ux_sub(ux(e2), ux(e1))
        num = ux_mul(ux(-q), diff)
        den = ux_add(ux(e1), ux(e2))
        return [("f", ux_fin(ux_div(num, den, ident))), ("f", p[0]), ("f", p[1]), ("f", -p[2])]
    if fn == "charge::lorentz_force":
        return fv_lorentz(v[0], v[1:4], v[4:7], v[7:10])
    return fv_lorentz(v[0], v[5:8], v[8:11], v[11:14])


def fv_rows(root: Path, module: str, path: Path):
    """(row-18 mappings, problems) for field.tsv or charge.tsv: each governed row whose model differs from the oracle's
    outputs is mapped; on every other row of these functions, the port keeps Rust's chain, and the check that the
    unbounded chain gives the same value there is the guard's proof, row by row (Rust's chain from validate_field's
    or validate_charge's transcription)."""
    global EXACT_THR
    check = validate_field if module == "field" else validate_charge
    rust, rounded = [], []
    try:
        EXACT_THR = True
        check(root, path, rust, MPFR_P)
        EXACT_THR = False
        check(root, path, rounded, MPFR_P)
    finally:
        EXACT_THR = False
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    out, errs = [], []
    for n, (line, want, w_round) in enumerate(zip(lines, rust, rounded), 1):
        fn, _cls, ins, outs, _note = line.split("\t")
        if fn not in FV_FNS or ins == "-":
            continue
        if not toks_same(want, w_round):
            continue                                    # row 17's (the exact threshold changes it)
        toks = [parse_tok(t) for t in ins.split(" ")]
        v = [b for _k, b in toks]
        if fn == "field::electric_field_superposition":
            v = v[1:]                                   # the charges' count token
        if fn == "field::trace_field_line":
            try:
                EXACT_THR = True
                res = TRACE_FN(v, toks, port=True)
            finally:
                EXACT_THR = False
            model, gov, perr = res if isinstance(res, tuple) else (res, False, [])     # an error: no step ran
            errs += [f"{module} row {n}: {e}" for e in perr]
            if not gov:
                if not toks_same(model, want):
                    errs.append(f"{module} row {n}: {fn}: the port's trace {model} differs from Rust's {want} where "
                                f"every step keeps Rust's chain")
                continue
            old = [(k, tok_text(b) if k == "s" else b) for k, b in (parse_tok(t) for t in outs.split(" "))]
            if not toks_same(model, old):
                out.append((fn, ins, wire(model), "18", f"oracle {outs}; ADR 0002 row 18"))
            continue
        if fn in M5_FNS and want[0][0] == "E":
            continue                                    # an order error: no chain
        args = fv_args(fn, v)
        if args is not None and not all(isinstance(x, float) and math.isfinite(x) for x in args):
            continue                                    # a non-finite part keeps Rust's chain
        row = fv_governed(fn, v)
        if row == 0:
            if fn == "field::image_charge_plane":
                proof = all(math.isfinite(x) for x in v)
            elif fn == "field::image_charge_dielectric":
                proof = all(math.isfinite(x) for x in v) and v[3] > 0.0 and v[4] > 0.0 and v[5] > 0.0
            elif fn == "FieldVector::normalized":
                proof = all(math.isfinite(x) for x in args) and any(x != 0.0 for x in args)
            else:
                proof = all(math.isfinite(x) for x in args)
            if proof:
                if not toks_same(fv_model(fn, v), want):
                    errs.append(f"{module} row {n}: {fn}{tuple(v)}: the unbounded model {fv_model(fn, v)} differs "
                                f"from Rust's chain {want} where the port keeps the chain")
            continue
        model = fv_model(fn, v)
        old = [(k, tok_text(b) if k == "s" else b) for k, b in (parse_tok(t) for t in outs.split(" "))]
        if toks_same(model, old):
            continue
        out.append((fn, ins, wire(model), "18", f"oracle {outs}; ADR 0002 row 18"))
    return out, errs


def adr_rows(root: Path, module: str, path: Path) -> list:
    """(fn, raw inputs, outputs, ADR 0002 row, note) for each golden row of `module` an ADR 0002 row changes."""
    if module == "complex":
        return cx_rows(path)[0]
    if module == "material":
        return material_adr_rows(root, path)
    if module == "wave":
        return wave_adr_rows(root, path)
    if module == "circuit":
        return circuit_adr_rows(path)
    if module == "relativity":
        return relativity_adr_rows(root, path)
    maps = {13: adr13(root)}
    out = adr17_rows(root, module, path) if module in ("field", "charge") else []
    if module in ("field", "charge"):
        out += fv_rows(root, module, path)[0]
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


# Rows 14, 15 and 16 (M1 part 5c): BjComplex's range edges. At an edge (a finite, nonzero component outside
# [2^-511, 2^511)) hisab's chain runs with an unbounded exponent, each operation rounded to 53 bits, the last one once
# into binary64 (IEEE's signed zeros throughout); in range the chain IS hisab's, which cx_model checks row by row.
# powf has no cutoff (row 15). Row 16 extends the unbounded chain to its transcendental operations: e^re, cosh and
# sinh, ln and pow of an |z| or a power past either end, and the angles arg z / 2 and n·arg z.
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


# Row 16 (M1 part 5c-2): the transcendental operations of the chain with an unbounded exponent too. Each is MPFR's
# value at 53 bits with its widest exponent range (mpfr_ref.wide), whatever the port's guard: the model is the whole
# unbounded chain, so on the rows the port keeps hisab's chain for, cx_rows proves the guard as it does row 14's.
# A value past 2^5000 (or below 2^-5000) stands as 2^±6000: any product of it with a nonzero double overflows (rounds
# to zero) exactly as the true value's does, and an exact one would be billions of bits.
W_FAR = 5000


def wv(name, *xs):
    """mpfr_ref.wide as a model pair (Fraction, negative): float arguments keep their zeros' signs."""
    r = mpfr_ref.wide(name, *xs)
    if isinstance(r, float):
        if r != r or r in (math.inf, -math.inf):
            raise ValueError(f"{name}{xs}: {r} in the widest exponent range")
        return Fraction(0), math.copysign(1.0, r) < 0
    return r, r < 0


def w_log2(a) -> float:
    """log2 of a positive Fraction, to a double's precision (it may be far outside the binary64 range)."""
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if a < Fraction(2) ** e:
        e -= 1
    return e + math.log2(float(a / Fraction(2) ** e))


def wv_exp(x: float):
    if x > W_FAR:
        return Fraction(2) ** 6000, False
    if x < -W_FAR:
        return Fraction(2) ** -6000, False
    return wv("exp", x)


def wv_hyp(name, x: float):
    if abs(x) > W_FAR:
        return (Fraction(2) ** 6000) * (-1 if name == "sinh" and x < 0 else 1), name == "sinh" and x < 0
    return wv(name, x)


def wv_pow(a, n: float):
    """|z|^n for the model pair a (positive) and a finite n."""
    if n == 0:
        return Fraction(1), False
    t = n * w_log2(a[0])
    if t > 4 * W_FAR:                                   # past 2^20000: its products overflow
        return Fraction(2) ** 30000, False
    if t < -4 * W_FAR:
        return Fraction(2) ** -30000, False
    return wv("pow", a[0], n)


def w_trig(name, t):
    """cos or sin of the model pair t, to 53 bits: of ±0, 1 and ±0."""
    if t[0] == 0:
        return (Fraction(1), False) if name == "cos" else t
    return wv(name, t[0])


def ux_abs(re, im):
    return ux_sqrt(ux_sumsq(re, im))


def ux_is_normal(a) -> bool:
    return a[0] != 0 and -1022 <= binade(a[0]) <= 1023


MIN_NORMAL = 2.0 ** -1022


def powf_plain(r, g, th, n, im) -> bool:
    """src/complex.cyr's _bj_cx_powf_plain: every intermediate of powf's binary64 chain a 53-bit value."""
    def nrm(x):
        return MIN_NORMAL < abs(x) < math.inf
    if not nrm(r):
        return False
    if im != 0 and not abs(g) > MIN_NORMAL:
        return False
    if th == 0:
        return g == 0 or n == 0
    return nrm(th)


def cx_powf_path(re, im, n) -> str:
    """The path src/complex.cyr's powf takes for a finite z ≠ 0 and a finite n: "plain" (hisab's chain), "abs" (the
    unbounded |z|, rounded, then hisab's chain: row 14) or "edge" (the unbounded chain: row 16)."""
    P = MPFR_P
    g = P["atan2"](im, re)
    th = g * n
    if cx_plain(re, im):
        return "plain" if powf_plain(P["pow"](cx_abs(re, im), n), g, th, n, im) else "edge"
    a = ux_abs(re, im)
    if ux_is_normal(a) and powf_plain(P["pow"](float(a[0]), n), g, th, n, im):
        return "abs"
    return "edge"


def cx_sqrt_tiny(re, im) -> bool:
    """src/complex.cyr's sqrt: im ≠ 0 and |arg z| < 2^-1020 take the unbounded angle (row 16)."""
    return im != 0 and abs(MPFR_P["atan2"](im, re)) < 2.0 ** -1020


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
    ident = lambda t: t                                  # noqa: E731
    if fn == "Complex::sqrt":
        r = ux_sqrt(ux_abs(v[0], v[1]))
        g = wv("atan2", v[1], v[0])
        th = (g[0] / 2, g[1])
        return [ux_fin(ux_mul(r, w_trig("cos", th), ident)), ux_fin(ux_mul(r, w_trig("sin", th), ident))]
    if fn == "Complex::ln":
        a = ux_abs(v[0], v[1])
        lv = float(mpfr_ref.wide("log", a[0])) if a[0] != 0 else -math.inf
        return [lv, P["atan2"](v[1], v[0])]
    if fn == "Complex::powf":
        re, im, n = v
        if re == 0 and im == 0:
            if n != n:
                return [n, n]
            return [0.0 if n > 0 else (1.0 if n == 0 else float("inf")), 0.0]
        pw = wv_pow(ux_abs(re, im), n)
        g = wv("atan2", im, re)
        th = ux_mul(g, ux(n))
        return [ux_fin(ux_mul(pw, w_trig("cos", th), ident)), ux_fin(ux_mul(pw, w_trig("sin", th), ident))]
    if fn == "Complex::exp":
        e = wv_exp(v[0])
        return [ux_fin(ux_mul(e, ux(P["cos"](v[1])), ident)), ux_fin(ux_mul(e, ux(P["sin"](v[1])), ident))]
    if fn in ("Complex::sin", "Complex::cos"):
        ch, sh = wv_hyp("cosh", v[1]), wv_hyp("sinh", v[1])
        s, c = P["sin"](v[0]), P["cos"](v[0])
        if fn == "Complex::sin":
            return [ux_fin(ux_mul(ux(s), ch, ident)), ux_fin(ux_mul(ux(c), sh, ident))]
        return [ux_fin(ux_mul(ux(c), ch, ident)), ux_fin(ux_mul(ux(-s), sh, ident))]
    return None


def cx_governed(fn, v) -> int:
    """The ADR 0002 row whose model sets this complex row's expectation (14, 15 or 16), or 0 where hisab's chain (with
    the corrected oracle's primitives) still does: the path src/complex.cyr takes."""
    if fn not in ("Complex::abs", "ComplexExt::norm", "Complex::norm_sq", "Complex::mul", "Complex::inv",
                  "Complex::div", "Complex::is_zero", "Complex::sqrt", "Complex::ln", "Complex::powf",
                  "Complex::exp", "Complex::sin", "Complex::cos"):
        return 0
    if not all(math.isfinite(x) for x in v):
        return 0
    if fn == "Complex::exp":
        return 0 if -708.0 <= v[0] <= 709.0 else 16
    if fn in ("Complex::sin", "Complex::cos"):
        return 0 if abs(v[1]) <= 710.0 else 16
    if fn == "Complex::powf":
        re, im, n = v
        if (re == 0 and im == 0) or re * re + im * im < TINY_300:
            return 15
        return {"plain": 0, "abs": 14, "edge": 16}[cx_powf_path(re, im, n)]
    if fn in ("Complex::mul", "Complex::div"):
        return 0 if all(cx_ok(x) for x in v) else 14
    if fn == "Complex::is_zero":
        return 0 if all(cx_ok(x) for x in v) else 14
    if fn == "Complex::ln":
        if cx_plain(v[0], v[1]):
            return 0
        return 14 if ux_is_normal(ux_abs(v[0], v[1])) else 16
    if fn == "Complex::sqrt" and cx_sqrt_tiny(v[0], v[1]):
        return 16
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
        "circuit": validate_circuit,
        "complex": validate_complex,
        "constants": lambda p: validate_constants(root, p),
        "error": validate_error,
        "field": lambda p: validate_field(root, p),
        "fmtf": validate_fmtf,
        "json": validate_json,
        "math": validate_math,
        "material": lambda p: validate_material(root, p),
        "maxwell": lambda p: validate_maxwell(root, p),
        "relativity": lambda p: validate_relativity(root, p),
        "wave": lambda p: validate_wave(root, p),
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
    errs += check_corrections(root, "wave", golden / "wave.tsv", wave_corrections(root, golden / "wave.tsv"))
    errs += check_corrections(root, "circuit", golden / "circuit.tsv", circuit_corrections(golden / "circuit.tsv"))
    errs += check_corrections(root, "relativity", golden / "relativity.tsv",
                              relativity_corrections(root, golden / "relativity.tsv"))
    for mod in ("charge", "circuit", "complex", "field", "material", "maxwell", "relativity", "wave"):
        errs += check_adr(root, mod, golden / f"{mod}.tsv")
    errs += cx_rows(golden / "complex.tsv")[1]
    errs += fv_rows(root, "field", golden / "field.tsv")[1] + fv_rows(root, "charge", golden / "charge.tsv")[1]
    errs += check_tables(root)
    errs += check_wave_data(root)
    for e in errs[:40]:
        print("  " + e)
    if len(errs) > 40:
        print(f"  … {len(errs) - 40} more")
    print(f"validate: {total} rows, {len(errs)} problems")
    return 0 if not errs else 1


if __name__ == "__main__":
    sys.exit(main())
