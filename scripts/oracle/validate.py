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

Run with `python3 -I -B scripts/oracle/validate.py .` from the repo root. Exit status 0 when
every row agrees, 1 otherwise. Copied in spirit from impetus's and kimiya's validate.py.
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
for _name in ("sin", "cos", "tan", "asin", "acos", "atan", "exp", "log10", "cbrt"):
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


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    golden = root / "tests/golden"
    checks = {
        "constants": lambda p: validate_constants(root, p),
        "math": validate_math,
    }
    listed = [l for l in (golden / "modules.txt").read_text().splitlines() if l and not l.startswith("#")]
    errs, total = [], 0
    for mod in listed:
        path = golden / f"{mod}.tsv"
        if mod not in checks:
            errs.append(f"{mod}: no validator (add one when the module lands)")
            continue
        _h, rows = read_tsv(path)
        e = checks[mod](path)
        total += len(rows)
        print(f"{'ok ' if not e else 'BAD'}  {mod}: {len(rows)} rows" + (f", {len(e)} disagreements" if e else ""))
        errs.extend(e)
    for e in errs[:40]:
        print("  " + e)
    if len(errs) > 40:
        print(f"  … {len(errs) - 40} more")
    print(f"validate: {total} rows, {len(errs)} problems")
    return 0 if not errs else 1


if __name__ == "__main__":
    sys.exit(main())
