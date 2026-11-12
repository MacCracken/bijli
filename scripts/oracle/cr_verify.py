#!/usr/bin/env python3
"""Verifies bijli's correctly rounded primitives against MPFR, beyond the golden corpus (bijli ADR 0001).

Usage, from the repo root:
    python3 -I -B scripts/oracle/cr_verify.py . [--n N] [--aarch64] [fn ...]

For each function it draws N arguments (default 100,000) from seeded sets that cover the whole domain, every binade
and the implementation's thresholds, writes them to build/cr_in.txt, builds scripts/oracle/cr_driver.cyr (against
the math sources) natively, and with --aarch64 also for aarch64 under qemu-aarch64, runs it, and compares every result
with MPFR's correctly rounded value (mpfr_ref.py). Any difference is printed and fails the run (exit 1). This is how
a new primitive earns its 0-ulp gate (port plan §6 step 8: ≥ 10⁶ arguments for asin and tan).
"""
import math
import os
from fractions import Fraction
import random
import struct
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpfr_ref  # noqa: E402

root = Path(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else ".").resolve()
args = [a for a in sys.argv[2:]]
n = 100_000
if "--n" in args:
    n = int(args[args.index("--n") + 1])
aarch64 = "--aarch64" in args
wanted = [a for a in args if not a.startswith("--") and not a.isdigit()]


def f2b(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def b2f(b: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", b & 0xFFFFFFFFFFFFFFFF))[0]


def rand_double(rng) -> float:
    while True:
        x = b2f(rng.getrandbits(64))
        if x == x and abs(x) != float("inf"):
            return x


def around(x):
    b = f2b(x)
    return [b2f(b + d) for d in (-2, -1, 0, 1, 2)]


SPECIALS = [0.0, -0.0, 1.0, -1.0, 2.0, 0.5, float("inf"), float("-inf"), float("nan"), 5e-324, -5e-324,
            2.2250738585072014e-308, 1.7976931348623157e308, -1.7976931348623157e308]


def unary_args(name, rng, count):
    xs = list(SPECIALS)
    share = count // 4
    if name == "exp":
        xs += [rng.uniform(-745.2, 709.8) for _ in range(share)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(0, 60) for _ in range(share)]
        for e in (2.0 ** -55, 708.0, 709.7, 710.0, -746.0, -708.4, -745.13):
            xs += around(e)
    elif name in ("sin", "cos"):
        xs += [rng.uniform(-7.0, 7.0) for _ in range(share)]
        xs += [rng.uniform(-1e6, 1e6) for _ in range(share // 2)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** rng.randint(-60, 1023) for _ in range(share // 2)]
        pio2 = b2f(0x3FF921FB54442D18)
        for k in range(1, 9):
            xs += around(k * pio2)
        for e in (2.0 ** -27, 0.78, 2.0 ** 20):
            xs += around(e) + around(-e)
    elif name == "atan":
        xs += [rng.uniform(-10.0, 10.0) for _ in range(share)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** rng.randint(-1074, 1023) for _ in range(share)]
    elif name in ("asin", "acos"):
        xs += [rng.uniform(-1.0, 1.0) for _ in range(share)]
        xs += [math.copysign(1.0 - rng.uniform(0, 1) * 2.0 ** -rng.randint(1, 53), rng.uniform(-1, 1))
               for _ in range(share // 2)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(0, 80) for _ in range(share // 2)]
        for e in (2.0 ** -26, 2.0 ** -56, 0.5, 2.0 ** -0.5, 1.0):
            xs += around(e) + around(-e)
        xs += [b2f(f2b(1.0) - k) for k in range(1, 2000)] + [-b2f(f2b(1.0) - k) for k in range(1, 2000)]
        xs += [rng.uniform(-3.0, 3.0) for _ in range(1000)]
    elif name in ("log10", "ln"):
        xs += [rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1074, 1023) for _ in range(share)]
        xs += [1.0 + rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(1, 60) for _ in range(share // 2)]
        xs += [b2f(rng.getrandbits(52)) for _ in range(share // 4)]          # subnormals
        for k in range(-30, 31):
            xs += around(10.0 ** k)
        xs += [-rng.uniform(0.0, 10.0) for _ in range(500)]
    elif name in ("tan", "tan_slow"):
        xs += [rng.uniform(-7.0, 7.0) for _ in range(share)]
        xs += [rng.uniform(-1e6, 1e6) for _ in range(share // 2)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** rng.randint(-60, 1023) for _ in range(share // 2)]
        pio2 = b2f(0x3FF921FB54442D18)
        # The zeros and poles: the doubles around k·π/2 (k·fl(π/2) is within a few ulps), to |x| = 2^20 and past.
        for k in range(1, 65):
            xs += around(k * pio2) + around(-k * pio2)
        xs += [b2f(f2b(rng.randint(1, 700000) * pio2) + rng.randint(-3, 3)) for _ in range(share // 4)]
        xs += around(6381956970095103.0 * 2.0 ** 797)    # the double closest to a multiple of π/2
        for e in (2.0 ** -27, 2.0 ** -24, 0.78, 2.0 ** 20):
            xs += around(e) + around(-e)
    elif name in ("sinh", "cosh", "sinh_slow", "cosh_slow"):
        xs += [rng.uniform(-1.0, 1.0) for _ in range(share)]
        xs += [rng.uniform(-30.0, 30.0) for _ in range(share // 2)]
        xs += [rng.uniform(-711.0, 711.0) for _ in range(share // 2)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(0, 40) for _ in range(share // 4)]
        xs += [math.copysign(rng.uniform(700.0, 712.0), rng.uniform(-1, 1)) for _ in range(share // 8)]
        # The thresholds: the tiny shortcut, the series' end, e^-ax's last term, exp's phases' range, overflow.
        for e in (2.0 ** -26, 0.5, 708.0, 709.7, 709.78, 710.4758600739439, 710.4758600739440, 711.0):
            xs += around(e) + around(-e)
    xs += [rand_double(rng) for _ in range(count - len(xs))]
    if name in ("sinh_slow", "cosh_slow"):    # the fixed-point path's domain: finite, 2^-26 ≤ |x| < 711
        xs = [x for x in xs if x == x and 2.0 ** -26 <= abs(x) < 711.0]
    if name == "tan_slow":    # the fixed-point path's domain: finite, |x| ≥ 2^-27
        xs = [x for x in xs if x == x and 2.0 ** -27 <= abs(x) < float("inf")]
    return [(x,) for x in xs[:count]]


def binary_args(name, rng, count):
    pairs = [(a, b) for a in SPECIALS for b in SPECIALS]
    share = count // 3
    if name == "atan2":
        pairs += [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(share)]
    elif name == "hypot":
        pairs += [(rng.uniform(-1e3, 1e3), rng.uniform(-1e3, 1e3)) for _ in range(share)]
        pairs += [(a * 3.0, a * 4.0) for a in (rng.uniform(1, 2) * 2.0 ** rng.randint(-1000, 1000) for _ in range(500))]
    elif name == "pow":
        pairs += [(rng.uniform(0, 4), rng.uniform(-10, 10)) for _ in range(share)]
        pairs += [(rng.uniform(0.5, 2.0), rng.uniform(-1000, 1000)) for _ in range(share // 2)]
        pairs += [(-rng.uniform(0.1, 10), float(rng.randint(-20, 20))) for _ in range(share // 4)]
    while len(pairs) < count:
        a, b = rand_double(rng), rand_double(rng)
        if name == "pow":
            a = abs(a)
        pairs.append((a, b))
    return pairs[:count]


FUNCS = {
    "exp": (1, "exp"), "sin": (1, "sin"), "cos": (1, "cos"), "atan": (1, "atan"),
    "atan2": (2, "atan2"), "hypot": (2, "hypot"), "pow": (2, "pow"),
    "acos": (1, "acos"), "asin": (1, "asin"), "log10": (1, "log10"),
    "tan": (1, "tan"), "tan_slow": (1, "tan"), "ln": (1, "log"),
    "sinh": (1, "sinh"), "cosh": (1, "cosh"), "sinh_slow": (1, "sinh"), "cosh_slow": (1, "cosh"),
}


# ── kimiya's exact rationals and unbounded exponent (src/crmath_kx.cyr, M1 part 5): exact references ──────────────
# Each family's reference is exact: Python's Fraction for the monomials and the unbounded operations (CPython's int
# division rounds once, to nearest even, onto the subnormal grid), an integer square root for the roots. A call line
# carries words: an f64's bits, or an integer in two's complement.

from exact import INF, rn_frac, rn_sqrt  # noqa: E402  (shared with validate.py)


def w64(x) -> str:
    """A call word: an f64's bits, or an int in two's complement."""
    return f"{(x if isinstance(x, int) else f2b(x)) & 0xFFFFFFFFFFFFFFFF:016x}"


def any_finite(rng):
    """A finite double from any binade (subnormals included), either sign, sometimes a small integer."""
    c = rng.random()
    if c < 0.1:
        return float(rng.choice([1, 2, 3, 5, 7, 9, 10, 100, 1000])) * rng.choice([1.0, -1.0])
    if c < 0.2:
        return b2f(rng.getrandbits(52)) * rng.choice([1.0, -1.0])
    return rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1022, 1023) * rng.choice([1.0, -1.0])


def near_one(rng):
    return rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 60) * rng.choice([1.0, -1.0])


def mono_ref(num, den):
    v = Fraction(1)
    for x in num:
        v *= Fraction(x)
    for y in den:
        v /= Fraction(y)
    if v == 0:          # a zero factor: the signed zero of the factors' signs
        neg = sum(math.copysign(1.0, x) < 0 for x in num + den) % 2
        return -0.0 if neg else 0.0
    return rn_frac(v)


SHAPES = {  # name: (numerator count, denominator count)
    "mul_div": (2, 1), "mul3": (3, 0), "div2": (1, 2), "mul3_div": (3, 1), "mul5": (5, 0),
    "mul2_div2": (2, 2), "mul3_div2": (3, 2), "mul4": (4, 0), "div3": (1, 3),
}


def factor_sets(rng, nn, nd, count):
    """Factor lists: any binade, near one (ties and the fast path's rounding test), products that leave the range."""
    out = []
    for i in range(count):
        mode = i % 4
        pick = any_finite if mode in (0, 1) else near_one
        num = [pick(rng) for _ in range(nn)]
        den = [pick(rng) for _ in range(nd)]
        den = [d if d != 0 else 1.0 for d in den]
        if nn == 5:     # the exact path's capacity: four doubles and a 41-bit integer per side (kimiya's domain)
            num[4] = float(rng.randint(1, (1 << 41) - 1)) * rng.choice([1.0, -1.0])
        if mode == 3 and nn >= 1:
            # a midpoint: 3·y with 3y a 54-bit integer over a power-of-two denominator
            num[0] = float(rng.getrandbits(52) | (1 << 52)) * 2.0 ** rng.randint(-30, 30)
            if nn >= 2:
                num[1] = 3.0
        out.append((num, den))
    return out


def exact_families(rng, count):
    """(family, call words, reference) for every exact family."""
    calls = []
    for name, (nn, nd) in SHAPES.items():
        for num, den in factor_sets(rng, nn, nd, count):
            calls.append((name, num + den, mono_ref(num, den)))
    for i in range(count):
        nn, nd = rng.randint(0, 4), rng.randint(0, 4)
        num, den = factor_sets(rng, nn, nd, 1)[0]
        calls.append(("mono", [nn, nd] + num + den, mono_ref(num, den)))
    for i in range(count):
        if i % 3 == 0:
            x, y = abs(any_finite(rng)), abs(any_finite(rng))
        elif i % 3 == 1:
            x, y = abs(near_one(rng)), abs(near_one(rng))
        else:
            m = rng.getrandbits(26) | 1
            y = abs(near_one(rng))
            x = float(m * m) * y if float(m * m) * y != INF else 4.0   # exact squares
        x, y = x or 1.0, y or 1.0
        calls.append(("sqrt_div", [x, y], rn_sqrt(Fraction(x) / Fraction(y))))
    for i in range(count):
        w, x, y = (abs(any_finite(rng)) or 1.0, abs(any_finite(rng)) or 1.0, abs(any_finite(rng)) or 1.0)
        if i % 2:
            w, x, y = abs(near_one(rng)), abs(near_one(rng)), abs(near_one(rng))
        if i % 7 == 0:
            x, y = 9.0, 1.0                     # w·3: a midpoint when 3w has 54 bits
            w = float(rng.getrandbits(52) | (1 << 52)) * 2.0 ** rng.randint(-40, 40)
        e = rng.randint(-2200, 2200) if i % 3 == 0 else rng.randint(-8, 8)
        ref = rn_sqrt(Fraction(w) ** 2 * Fraction(x) * Fraction(2) ** e / Fraction(y))
        calls.append(("wsqrt", [w, x, y, e], ref))
    # bijli's RN(√(Π num / Π den)) (src/crmath_mono.cyr): the full function and the exact path alone, on factor
    # sets of every size in the domain (≤ 3 over ≤ 4), any binade (subnormal factors, roots past DBL_MAX and on the
    # subnormal grid), near one, exact squares, and maxwell's two shapes at SI scales and at the range ends.
    for i in range(2 * count):
        mode = i % 6
        if mode == 4:      # wave_speed: 1/√(ε·μ)
            num, den = [], [abs(any_finite(rng)) or 1.0, abs(any_finite(rng)) or 1.0]
        elif mode == 5:    # skin_depth: √(2/(ω·μ·σ))
            num, den = [2.0], [abs(any_finite(rng)) or 1.0 for _ in range(3)]
        else:
            nn, nd = rng.randint(0, 3), rng.randint(0, 4)
            if nn + nd == 0:
                nd = 1
            pick = any_finite if mode < 2 else near_one
            num = [abs(pick(rng)) or 1.0 for _ in range(nn)]
            den = [abs(pick(rng)) or 1.0 for _ in range(nd)]
            if mode == 3 and nn:
                m = rng.getrandbits(26) | (1 << 25)
                num[0] = float(m * m) * 2.0 ** (2 * rng.randint(-200, 200))    # an exact square times …
                den = [d * d for d in den[:1]] if den and den[0] * den[0] not in (0.0, INF) else []
                num = num[:1]
        v = Fraction(1)
        for x in num:
            v *= Fraction(x)
        for y in den:
            v /= Fraction(y)
        fam = "sqrt_mono" if i % 2 == 0 else "sqrt_mono_acc"
        calls.append((fam, [len(num), len(den)] + num + den, rn_sqrt(v)))
    for op in ("ux_mul", "ux_div", "ux_add", "ux_sub"):
        for i in range(count):
            s1 = rng.uniform(1.0, 2.0) * rng.choice([1.0, -1.0])
            s2 = rng.uniform(1.0, 2.0) * rng.choice([1.0, -1.0])
            if i % 5 == 0:
                s2 = s1 * rng.choice([1.0, -1.0])   # cancellation and exact doubling
            e1 = rng.randint(-2600, 2600)
            e2 = e1 + rng.randint(-70, 70) if i % 2 else rng.randint(-2600, 2600)
            if op in ("ux_mul", "ux_div") and i % 3 == 0:
                e2 = rng.choice([-1, 1]) * rng.randint(900, 1200) - e1 * (1 if op == "ux_mul" else -1)
            a, b = Fraction(s1) * Fraction(2) ** e1, Fraction(s2) * Fraction(2) ** e2
            v = {"ux_mul": a * b, "ux_div": a / b, "ux_add": a + b, "ux_sub": a - b}[op]
            calls.append((op, [s1, e1, s2, e2], rn_frac(v)))
    for i in range(count):
        s1 = rng.uniform(1.0, 2.0)
        e1 = rng.randint(-4300, 2100)
        calls.append(("ux_sqrt", [s1, e1], rn_sqrt(Fraction(s1) * Fraction(2) ** e1)))
    return calls


def mpfr(name, a):
    return mpfr_ref.cr(name, a[0]) if len(a) == 1 else mpfr_ref.cr2(name, a[0], a[1])


def run(target, lines):
    out_bin = root / f"build/cr_driver_{target}"
    build = ["cyrius", "build"] + (["--aarch64"] if target == "aarch64" else []) + \
            ["scripts/oracle/cr_driver.cyr", str(out_bin.relative_to(root))]
    b = subprocess.run(build, cwd=root, capture_output=True, text=True)
    if b.returncode != 0 or "warning" in (b.stdout + b.stderr).lower():
        sys.exit("cr_verify: build failed or warned:\n" + b.stdout + b.stderr)
    cmd = (["qemu-aarch64", str(out_bin)] if target == "aarch64" else [str(out_bin)])
    r = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"cr_verify: driver exit {r.returncode}: {r.stdout[-300:]}{r.stderr[-300:]}")
    got = r.stdout.split("\n")[:-1]
    if len(got) != len(lines):
        sys.exit(f"cr_verify: {target}: {len(got)} results for {len(lines)} calls (truncated output?)")
    return [int(g, 16) for g in got]


EXACT = ["mono"] + list(SHAPES) + ["sqrt_div", "wsqrt", "sqrt_mono", "sqrt_mono_acc", "ux_mul", "ux_div", "ux_add",
         "ux_sub", "ux_sqrt"]


def main() -> int:
    (root / "build").mkdir(exist_ok=True)
    names = wanted or list(FUNCS) + ["exact"]
    rng = random.Random(0xB1_1E_C0DE)
    calls = []
    for name in names:
        if name == "exact" or name in EXACT:
            continue
        arity, mname = FUNCS[name]
        argsets = unary_args(name, rng, n) if arity == 1 else binary_args(name, rng, n)
        calls += [(name, (lambda m, a: f2b(mpfr(m, a)))(mname, a), a) for a in argsets]
    if "exact" in names or any(x in EXACT for x in names):
        erng = random.Random(0xC0DE_5A)
        fams = exact_families(erng, max(1, n // 5))
        calls += [(f, f2b(ref), ws) for f, ws, ref in fams if "exact" in names or f in names]
    lines = [" ".join([name] + [w64(x) for x in a]) for name, _, a in calls]
    (root / "build/cr_in.txt").write_text("\n".join(lines) + "\n")
    targets = ["x86_64"] + (["aarch64"] if aarch64 else [])
    results = {t: run(t, lines) for t in targets}
    bad = 0
    per = {}
    for i, (name, want, a) in enumerate(calls):
        for t in targets:
            g = results[t][i]
            ok = g == want or (b2f(g) != b2f(g) and b2f(want) != b2f(want))
            per.setdefault((name, t), [0, 0])[0 if ok else 1] += 1
            if not ok:
                bad += 1
                if bad <= 20:
                    print(f"MISMATCH {t} {name}{tuple(a)}: got {g:016x}, want {want:016x}")
    for (name, t), (ok, ko) in sorted(per.items()):
        ref = "the exact value" if name in EXACT else "MPFR"
        print(f"  {name:9} {t:7} {ok + ko:8} args: {ok} equal to {ref}, {ko} differ")
    print(f"cr_verify: {len(calls)} calls x {len(targets)} targets, {bad} mismatches (MPFR {mpfr_ref.VERSION}; "
          "Fraction and isqrt for the exact families)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
