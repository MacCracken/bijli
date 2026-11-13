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
    elif name == "expm1":
        xs += [rng.uniform(-40.0, 709.8) for _ in range(share)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(0, 60) for _ in range(share)]
        xs += [rng.uniform(-1.0, 1.0) * 2.0 ** -rng.randint(60, 1074) for _ in range(share // 4)]
        # The thresholds: the tiny shortcut, the small-|x| fast path's end, overflow, the −1 saturation.
        for e in (2.0 ** -54, 2.0 ** -27, 2.0 ** -9, 0.5, 1.0, 709.78, 710.0, 36.04, 37.43, 38.0):
            xs += around(e) + around(-e)
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
    "expm1": (1, "expm1"),
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
    # src/crmath_cmp.cyr's exact Σ (aᵢ − bᵢ)² < 10^-p (p = 30, 60): points crowded around the cutoff a few ulps either
    # side, differences that round (an offset with a tiny or subnormal partner), and far ones (the result as 0 or 1).
    for i in range(count):
        p = 30 if i % 2 == 0 else 60
        scale = 1e-15 if p == 30 else 1e-30
        u = [rng.uniform(-1.0, 1.0) for _ in range(3)]
        nrm = math.sqrt(sum(x * x for x in u)) or 1.0
        mode = i % 7
        a = [x / nrm * scale for x in u]
        if mode == 1:
            a = [b2f(f2b(x) + rng.randint(-40, 40)) if x else x for x in a]
        b = [0.0, 0.0, 0.0]
        if mode in (2, 3):
            tiny = [5e-324, 3e-321, 2.2250738585072014e-308, 1e-300, 3e-32, 1e-47][rng.randint(0, 5)]
            b = [rng.choice([-1.0, 1.0]) * tiny * rng.randint(1, 9) if rng.random() < 0.7 else 0.0 for _ in range(3)]
        if mode == 4:                                         # big, nearly equal points: exact differences
            off = [rng.uniform(-1e3, 1e3) for _ in range(3)]
            b = [o for o in off]
            a = [b2f(f2b(o) + rng.randint(-3, 3)) for o in off]
        if mode == 5:
            a = [x * rng.uniform(0.5, 2.0) for x in a]
        if mode == 6:
            a = [x * 2.0 ** rng.randint(-60, 60) for x in a]
        ex = sum((Fraction(x) - Fraction(y)) ** 2 for x, y in zip(a, b)) < Fraction(1, 10 ** p)
        calls.append(("r2_below", a + b + [p], 1.0 if ex else 0.0))
    # src/crmath_cmp.cyr's exact |x·y|, |x/y| and |x + y| < 10^-p (M2 part 2c-1): pairs whose exact value lies within a
    # few ulps of 10^-p (p = 30 mostly, and every p up to 60), sums that round (a tiny or subnormal partner) or cancel,
    # and far ones; any binade, either sign.
    for i in range(3 * count):
        op = ("below_prod", "below_quot", "below_sum")[i % 3]
        p = 30 if i % 4 else rng.randint(1, 60)
        target = Fraction(1, 10 ** p)
        mode = (i // 3) % 5
        x = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-200, 200) * rng.choice([1.0, -1.0])
        if op == "below_prod":
            y = float(target / Fraction(abs(x)))
        elif op == "below_quot":
            y = float(Fraction(abs(x)) / target)
            x, y = float(target * Fraction(abs(y))) * rng.choice([1.0, -1.0]), y
        else:
            x = float(target) * rng.uniform(0.0, 1.0) * rng.choice([1.0, -1.0])
            y = float(target - Fraction(abs(x))) if x > 0 else -float(target - Fraction(abs(x)))
            if mode == 2:
                x, y = float(target), rng.choice([1.0, -1.0]) * [5e-324, 3e-321, 1e-300, 1e-47][rng.randint(0, 3)]
            if mode == 3:
                big = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 40)
                x, y = big, -big + float(target) * rng.uniform(0.5, 1.5)
        y = b2f(f2b(y) + rng.randint(-3, 3)) if y and mode != 4 else y
        if mode == 4:
            y = y * 2.0 ** rng.randint(-100, 100) if y else y
        if op != "below_sum" and (x == 0.0 or y == 0.0 or not math.isfinite(x * 1.0) or not math.isfinite(y)):
            continue
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        fx, fy = Fraction(x), Fraction(y)
        v = fx * fy if op == "below_prod" else fx / fy if op == "below_quot" else fx + fy
        calls.append((op, [x, y, p], 1.0 if abs(v) < target else 0.0))
    return calls


UX_HUGE = 1 << 20      # src/crmath_ux.cyr's stand-in exponent for a value beyond every product's range


def split53(v):
    """(s, e) with v = s·2^e and 1 ≤ |s| < 2, for a nonzero Fraction of at most 53 bits."""
    a = abs(v)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if a < Fraction(2) ** e:
        e -= 1
    return float(v / Fraction(2) ** e), e


def ux_families(rng, count):
    """(family, call words, wanted result word) for src/crmath_ux.cyr's functions, against MPFR at 53 bits with its
    widest exponent range (mpfr_ref.wide): e^x, cosh and sinh as (s, e) pairs (two calls each, part 0 and 1), ln of
    s·2^e, (s·2^e)^y as a pair, and sin and cos of a 53-bit angle past DBL_MAX."""
    W = mpfr_ref.wide
    M = 0xFFFFFFFFFFFFFFFF
    calls = []

    def pair(fam, ws, se):
        calls.append((fam, ws + [0], f2b(se[0])))
        calls.append((fam, ws + [1], se[1] & M))

    sign = lambda: rng.choice([1.0, -1.0])
    for i in range(count):
        c = i % 6
        if c == 0:
            x = rng.uniform(-1500.0, 1500.0)
        elif c == 1:
            x = rng.uniform(700.0, 760.0) * sign()           # binary64's ends
        elif c == 2:
            x = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-300, 10) * sign()
        elif c == 3:
            x = float(rng.randint(-1500, 1500))
        elif c == 4:
            x = rng.uniform(1400.0, 1500.0) * sign()
        else:
            x = rng.uniform(1500.001, 1e6) * sign()          # the stand-in
        pair("ux_exp", [x], (1.0, UX_HUGE if x > 0 else -UX_HUGE) if abs(x) > 1500 else split53(W("exp", x)))
        y = rng.uniform(710.000001, 1500.0) * sign() if i % 5 else rng.uniform(1500.001, 1e5) * sign()
        if abs(y) > 1500:
            ch, sh = (1.0, UX_HUGE), (math.copysign(1.0, y), UX_HUGE)
        else:
            ch, sh = split53(W("cosh", y)), split53(W("sinh", y))
        pair("ux_cosh", [y], ch)
        pair("ux_sinh", [y], sh)
    for i in range(count):
        c = i % 7
        s1 = rng.uniform(1.0, 2.0)
        e1 = rng.randint(-2200, 2200)
        if c == 0:
            s1 = 1.0                                          # powers of two
        elif c == 1:
            s1, e1 = 2.0 - rng.uniform(0.0, 1e-9), rng.choice([-1, -1, 0, 2200, -2200])   # near 1 from below
        elif c == 2:
            s1, e1 = 1.0 + rng.uniform(0.0, 1e-12), rng.choice([0, 0, 1, -1, 1500])       # near 1 from above
        elif c == 3:
            s1, e1 = abs(any_finite(rng)) or 1.0, rng.randint(-1000, 1000)               # any normal s
            if s1 < 2.2250738585072014e-308:
                s1 = 1.5
        v = Fraction(s1) * Fraction(2) ** e1
        lv = W("log", v)
        calls.append(("ux_ln", [s1, e1], f2b(float(lv))))
    exact_y = [1.0, 2.0, -1.0, -2.0, 0.5, 1.5, 0.25, -0.5, 3.0, 0.75, -1.5, 1.25]
    i = 0
    while i < count:
        c = i % 8
        s1, e1 = rng.uniform(1.0, 2.0), rng.randint(-2200, 2200)
        if c == 3:                                            # exact powers: squares and fourth powers
            m = rng.getrandbits(13) | 1
            s1, e0 = split53(Fraction(m ** 4))
            e1 = 4 * rng.randint(-500, 500) + e0
        elif c == 5:                                          # 53-bit midpoints: w³ and a² of 54 bits
            if i % 16 < 8:
                w = rng.randrange(185365, 262144, 2) | 1     # w³ in [2^53, 2^54)
                s1, e0 = split53(Fraction(w * w))
                e1 = e0 + 2 * rng.randint(-1000, 1000)
            else:
                a = rng.randrange(94906267, 134217728, 2) | 1   # a² in [2^53, 2^54)
                s1, e0 = split53(Fraction(a))
                e1 = e0 + rng.randint(-1000, 1000)
        lnx = math.log(s1) + e1 * math.log(2.0)
        if lnx == 0.0:
            continue
        if c in (0, 1, 2):
            y = rng.uniform(-1480.0, 1480.0) / lnx if c < 2 else rng.uniform(-9900.0, 9900.0) / lnx
        elif c == 3:
            y = rng.choice(exact_y)
        elif c == 4:
            y = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-75, -50) * sign() / max(1.0, abs(lnx))   # near 1
        elif c == 5:
            y = 1.5 if i % 16 < 8 else 2.0
        elif c == 6:
            y = rng.uniform(10050.0, 1e6) * sign() / lnx      # the stand-in
        else:
            y = rng.uniform(-3.0, 3.0)
        t = y * lnx
        if 9999.0 < abs(t) < 10001.0:
            continue
        i += 1
        if t > 10000:
            pe = (1.0, UX_HUGE)
        elif t < -10000:
            pe = (1.0, -UX_HUGE)
        else:
            pe = split53(W("pow", Fraction(s1) * Fraction(2) ** e1, y))
        pair("ux_pow", [s1, e1, y], pe)
    for i in range(count):
        if i % 4 == 1:
            ts, te = rng.uniform(1.0, 2.0) * sign(), rng.choice([1024, 1025])
        elif i % 4 == 3:                                      # the Green's function's k·R, up to 2^2052
            ts, te = rng.uniform(1.0, 2.0) * sign(), rng.randint(1024, 2052)
        elif i % 4 == 2:
            kk = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(500, 1023)
            rr = Fraction(rng.uniform(1.0, 2.0)) * Fraction(2) ** rng.randint(524, 1025)
            ts, te = split53(rnd53_frac(Fraction(kk) * rr))
            if te < 1024:
                continue
        else:                                                 # powf's angle: RN53(arg·n) past DBL_MAX
            g = rng.choice([math.pi, -math.pi, 3.0 * math.pi / 4.0, math.pi / 2.0, 2.5, 3.1])
            n = rng.uniform(5.8e307, 1.7976931348623157e308) * sign()
            ts, te = split53(rnd53_frac(Fraction(g) * Fraction(n)))
            if te < 1024:
                continue
        th = Fraction(ts) * Fraction(2) ** te
        calls.append(("ux_sin", [ts, te], f2b(float(W("sin", th)))))
        calls.append(("ux_cos", [ts, te], f2b(float(W("cos", th)))))
    return calls


def rnd53_frac(v):
    """v rounded to 53 bits, ties to even, with an unbounded exponent (a nonzero Fraction)."""
    a = abs(v)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    if a < Fraction(2) ** e:
        e -= 1
    sc = a / Fraction(2) ** (e - 52)
    q, r = divmod(sc.numerator, sc.denominator)
    if 2 * r > sc.denominator or (2 * r == sc.denominator and q & 1):
        q += 1
    out = Fraction(q) * Fraction(2) ** (e - 52)
    return out if v > 0 else -out


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
         "ux_sub", "ux_sqrt", "r2_below", "gamma", "gamma_slow", "gamma2", "gamma2_slow", "gamma_ok", "prod5",
         "prod5_acc", "mono_wide", "mono_wide_acc", "below_prod", "below_quot", "below_sum", "cm", "cm_slow"]
UX = ["ux_exp", "ux_cosh", "ux_sinh", "ux_ln", "ux_pow", "ux_sin", "ux_cos"]


# ── the correctly rounded Lorentz factor (src/crmath_gamma.cyr, M2's range-edge Part, ADR 0002 row 30) ───────────────
C_LIGHT = 299792458.0
GAMMA_INVALID = (1 << 64) - 1
GAMMA_C_LO = 2.0 ** -30                     # the range of c that _bj_cr_lorentz accepts
GAMMA_C_HI = 2.0 ** 30


def gamma_ref(vs, c, mode):
    """RN(c/√(c² − Σv²)) (mode 0) or RN(c²/(c² − Σv²)) (mode 1) exactly, or the invalid marker (c outside
    [2^-30, 2^30] too)."""
    if not GAMMA_C_LO <= c <= GAMMA_C_HI:
        return GAMMA_INVALID
    if any(x != x or abs(x) >= c for x in vs):
        return GAMMA_INVALID
    d = Fraction(c) ** 2 - sum(Fraction(x) ** 2 for x in vs)
    if d <= 0:
        return GAMMA_INVALID
    q = Fraction(c) ** 2 / d
    return f2b(rn_sqrt(q) if mode == 0 else rn_frac(q))


def gamma_cases(rng, count, c):
    """Speeds and vectors: uniform in |v| < c, β = 1 − 2^-k, components of mixed binades, subnormal and tiny parts, the
    c boundary ulp by ulp, and vectors whose squares nearly cancel c² (|v| within an ulp or two of c)."""
    out = []
    for i in range(count):
        kind = i % 8
        n = rng.randint(1, 3)
        if kind == 0:
            vs = [rng.uniform(-c, c) for _ in range(n)]
            vs = [x / math.sqrt(n) for x in vs]
        elif kind == 1:
            b = 1.0 - rng.uniform(0.5, 1.0) * 2.0 ** -rng.randint(1, 52)
            vs = [math.copysign(b * c, rng.uniform(-1, 1))]
        elif kind == 2:
            u = [rng.gauss(0, 1) for _ in range(3)]
            s = math.sqrt(sum(x * x for x in u)) or 1.0
            b = 1.0 - 2.0 ** -rng.randint(1, 50)
            vs = [x / s * b * c for x in u]
        elif kind == 3:
            vs = [rng.uniform(-c, c) * 0.5, rng.uniform(-1, 1) * 2.0 ** -rng.randint(400, 1074), 0.0]
        elif kind == 4:
            x = b2f(f2b(c) - rng.randint(1, 8))
            vs = [x] + [rng.uniform(-1, 1) * 2.0 ** -rng.randint(0, 40) for _ in range(n - 1)]
        elif kind == 5:
            vs = [rng.uniform(-1, 1) * c * 2.0 ** -rng.randint(20, 40) for _ in range(n)]
        elif kind == 6:
            vs = [b2f(f2b(c) + rng.randint(-2, 2)) * rng.choice([1.0, -1.0])]
        else:
            vx = rng.uniform(0.3, 0.9) * c
            vy = math.sqrt(c * c - vx * vx) * (1.0 - 2.0 ** -rng.randint(30, 53))
            vs = [vx, vy]
        out.append(vs)
    return out


# ── bijli's monomial of up to six factors a side (src/crmath_mono.cyr, M2 part 2b-4) ───────────────────────────────
def odd_bits(x):
    """The bit length of x's integer significand with its trailing zeros dropped (kimiya's capacity measure)."""
    if x == 0 or not math.isfinite(x):
        return 0
    m = Fraction(abs(x))
    p, q = m.numerator, m.denominator
    while p % 2 == 0:
        p //= 2
    return p.bit_length()


def mono_wide_families(rng, count):
    """mono_wide (the full function) and mono_wide_acc (the exact path alone): up to six factors a side, at most four of
    them full-precision doubles (any binade, near one, subnormal) and the rest small integers or powers of two, within
    the exact path's 253 bits a side; the Larmor power's shape q·q·a·a/(6·π·ε₀·c·c·c) at SI scales and at the range
    ends; exact ties (a 54-bit odd 3y). The reference is the exact quotient, rounded once (Fraction)."""
    C, PI, EPS0 = 299792458.0, math.pi, 8.8541878128e-12
    out = []
    for i in range(count):
        kind = i % 5
        if kind == 4:
            q = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-100, 100) * rng.choice([1.0, -1.0])
            a = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-200, 200)
            if i % 3 == 0:
                q, a = 1.602176634e-19 * rng.uniform(0.5, 2.0), rng.uniform(1.0, 2.0) * 10.0 ** rng.randint(0, 25)
            num, den = [q, q, a, a], [6.0, PI, EPS0, C, C, C]
        else:
            nn, nd = rng.randint(0, 6), rng.randint(0, 6)
            pick = any_finite if kind in (0, 1) else near_one
            num = [pick(rng) for _ in range(min(nn, 4))] + [float(rng.randint(1, 1 << 30)) * 2.0 ** rng.randint(-8, 8)
                                                            for _ in range(nn - min(nn, 4))]
            den = [pick(rng) or 1.0 for _ in range(min(nd, 4))] + [float(rng.randint(1, 1 << 30))
                                                                   for _ in range(nd - min(nd, 4))]
            if kind == 3 and nn >= 2:
                num[0] = float(rng.getrandbits(52) | (1 << 52)) * 2.0 ** rng.randint(-30, 30)
                num[1] = 3.0
            rng.shuffle(num)
            rng.shuffle(den)
        if sum(odd_bits(x) for x in num) > 253 or sum(odd_bits(x) for x in den) > 253:
            continue
        fam = "mono_wide" if i % 2 == 0 else "mono_wide_acc"
        out.append((fam, [len(num), len(den)] + num + den, mono_ref(num, den)))
    return out


# ── bijli's product of five full-precision doubles (src/crmath_mono.cyr, M2 part 2b-3) ────────────────────────────────
def prod5_families(rng, count):
    """prod5 (the full function: kimiya's fast path, then the exact one) and prod5_acc (the exact path alone): five
    doubles of any binade (subnormal factors, products past DBL_MAX and on the subnormal grid), near one, products
    steered into the subnormal range and its boundary, exact ties (a 54-bit odd integer over a power of two), and
    signed zeros. The reference is the exact product, rounded once (Fraction)."""
    out = []
    for i in range(count):
        kind = i % 6
        if kind == 0:
            xs = [any_finite(rng) for _ in range(5)]
        elif kind == 1:
            xs = [near_one(rng) for _ in range(5)]
        elif kind == 2:
            # the product's exponent T steered to the subnormal range and its boundary, spread over the five factors
            t = rng.randint(-1080, -1015)
            es = [t // 5 + rng.randint(-100, 100) for _ in range(4)]
            es.append(t - sum(es))
            xs = [rng.uniform(1.0, 2.0) * 2.0 ** e * rng.choice([1.0, -1.0]) for e in es]
        elif kind == 3:
            # an exact tie: q odd with 2^52 ≤ q < 2^54/3, so 3q is a 54-bit odd integer, times powers of two
            q = rng.randrange((1 << 52) | 1, (1 << 54) // 3, 2)
            xs = [float(q), 3.0] + [2.0 ** rng.randint(-200, 200) for _ in range(3)]
        elif kind == 4:
            xs = [any_finite(rng) for _ in range(5)]
            xs[rng.randrange(5)] = rng.choice([0.0, -0.0])
        else:
            xs = [b2f(rng.getrandbits(52)) * rng.choice([1.0, -1.0])] + [near_one(rng) for _ in range(4)]
        fam = "prod5" if i % 2 == 0 else "prod5_acc"
        out.append((fam, xs, mono_ref(xs, [])))
    return out


def gamma_families(rng, count):
    """Each family over c = the speed of light, 1, the range's two edges and four random c in it (count/8 cases each),
    plus c just outside the range, ∞ and NaN, which are refused. gamma_ok is the validity test alone
    (_bj_cr_lorentz_ok): 1 when the value exists, else 0."""
    fams = []
    for name, mode in (("gamma", 0), ("gamma_slow", 0), ("gamma2", 1), ("gamma2_slow", 1), ("gamma_ok", 0)):
        cs = [C_LIGHT, 1.0, GAMMA_C_HI, GAMMA_C_LO] + [2.0 ** rng.uniform(-30.0, 30.0) for _ in range(4)]
        ref = gamma_ref
        if name == "gamma_ok":
            ref = lambda vs, c, mode: 0 if gamma_ref(vs, c, mode) == GAMMA_INVALID else 1
        for c in cs:
            for vs in gamma_cases(rng, max(1, count // len(cs)), c):
                fams.append((name, [len(vs), c] + vs, ref(vs, c, mode)))
        for c in (2.0 * GAMMA_C_HI, 0.5 * GAMMA_C_LO, math.inf, math.nan):
            fams.append((name, [1, c, 0.0], ref([0.0], c, mode)))
    return fams



# ── Clausius–Mossotti's ε_r (src/crmath_edge.cyr, M2 part 2c-2, ADR 0002 row 17) ──────────────────────────────────────
CM_INVALID = (1 << 64) - 1                  # X ≥ 1
CM_EPS0 = float(next(ln for ln in (root / "src/constants.cyr").read_text().splitlines()
                     if ln.startswith("public const BJ_EPSILON_0 = ")).split("= ")[1].rstrip(";"))


def cm_ref(a, nd):
    """RN((3ε₀ + 2P)/(3ε₀ − P)) for P = Nα exactly, ε₀ its double, or the invalid marker when 3ε₀ − P ≤ 0."""
    p = Fraction(nd) * Fraction(a)
    t = 3 * Fraction(CM_EPS0)
    if t - p <= 0:
        return CM_INVALID
    return f2b(rn_frac((t + 2 * p) / (t - p)))


def cm_ulps(x, k):
    """The double k steps from the positive double x."""
    return b2f(f2b(x) + k)


def cm_cases(rng, count):
    """(α, N) pairs: physical ones; X uniform in (−4, 1) over N of every binade; X within a few ulps of 1 and of −1/2
    (ε_r's pole and zero); any finite α and N ≥ 0; the exact ties (ε_r = 9e·2^-j − 2 near the pole and −2 + d·2^-53
    below the zero, e ε₀'s odd significand and d = 3f for f dividing it), at scales 2^±s; the fast region's edges
    (χ near −0.24 and 0.49) and the shortcuts' (|χ| near 2^-60 and 2^136); and the exact zero and pole (P = −1.5ε₀,
    P = 3ε₀), zeros, Nα overflowing and underflowing."""
    b = f2b(CM_EPS0)
    m = (b & ((1 << 52) - 1)) | (1 << 52)
    ex = ((b >> 52) & 0x7FF) - 1075
    while m % 2 == 0:
        m //= 2
        ex += 1
    divs = [f for f in range(1, 10000) if m % f == 0] + [m // f for f in range(1, 10000) if m % f == 0]
    divs = sorted(set(f for f in divs if 3 * f <= 1 << 53))
    out = []
    for i in range(count):
        kind = i % 8
        if kind == 0:
            out.append((rng.choice([1.0, -1.0]) * 10.0 ** rng.uniform(-42.0, -38.0), 10.0 ** rng.uniform(24.0, 29.0)))
        elif kind == 1:
            nd = 2.0 ** rng.uniform(-1000.0, 1000.0)
            x = rng.uniform(-4.0, 1.0)
            out.append((x * 3.0 * CM_EPS0 / nd, nd))
        elif kind == 2:
            nd = rng.choice([1.0, 2.0 ** rng.randint(-900, 900), 10.0 ** rng.uniform(-300.0, 300.0)])
            out.append((cm_ulps(3.0 * CM_EPS0 / nd, rng.randint(-6, 6)), nd))
        elif kind == 3:
            nd = rng.choice([1.0, 2.0 ** rng.randint(-900, 900), 10.0 ** rng.uniform(-300.0, 300.0)])
            out.append((-cm_ulps(1.5 * CM_EPS0 / nd, rng.randint(-6, 6)), nd))
        elif kind == 4:
            out.append((any_finite(rng), abs(any_finite(rng))))
        elif kind == 5:
            sc = 2.0 ** rng.randint(-400, 400)
            if rng.random() < 0.5:
                j = rng.randint(0, 50)
                out.append((float((3 * m - (1 << j)) * Fraction(2) ** ex) / sc, sc))
            else:
                f = rng.choice(divs)
                out.append((-float(3 * (m // f) * Fraction(2) ** ex) / sc, float((1 << 53) - f) * sc))
        elif kind == 6:
            x = rng.choice([-0.24, 0.49, 2.0 ** -60, -(2.0 ** -60), 2.0 ** 136, -(2.0 ** 136)])
            nd = 2.0 ** rng.uniform(-200.0, 200.0)
            out.append((cm_ulps(abs(x) * 3.0 * CM_EPS0 / nd, rng.randint(-4, 4)) * math.copysign(1.0, x), nd))
        else:
            sc = 2.0 ** rng.randint(-500, 500)
            out.append(rng.choice([(-CM_EPS0 / sc, 1.5 * sc), (CM_EPS0 / sc, 3.0 * sc), (0.0, 1e28), (-0.0, 1e28),
                                   (1e-40, 0.0), (-1e10, 1e300), (1e10, 1e300), (-1e-300, 1e-300),
                                   (5e-324, 1e300), (-5e-324, 1e300), (-1.7976931348623157e308, 1e10)]))
    return out


def cm_families(rng, count):
    """cm (the full function) and cm_slow (the exact check alone) on each of cm_cases' pairs, against cm_ref."""
    out = []
    for a, nd in cm_cases(rng, count):
        want = cm_ref(a, nd)
        out.append(("cm", [a, nd], want))
        out.append(("cm_slow", [a, nd], want))
    return out


def main() -> int:
    (root / "build").mkdir(exist_ok=True)
    names = wanted or list(FUNCS) + ["exact", "ux"]
    rng = random.Random(0xB1_1E_C0DE)
    calls = []
    for name in names:
        if name in ("exact", "ux") or name in EXACT or name in UX:
            continue
        arity, mname = FUNCS[name]
        argsets = unary_args(name, rng, n) if arity == 1 else binary_args(name, rng, n)
        calls += [(name, (lambda m, a: f2b(mpfr(m, a)))(mname, a), a) for a in argsets]
    if "exact" in names or any(x in EXACT for x in names):
        erng = random.Random(0xC0DE_5A)
        fams = exact_families(erng, max(1, n // 5))
        calls += [(f, f2b(ref), ws) for f, ws, ref in fams if "exact" in names or f in names]
        wrng = random.Random(0x3_1D_E6)
        calls += [(f, f2b(ref), ws) for f, ws, ref in mono_wide_families(wrng, 2 * max(1, n // 5))
                  if "exact" in names or f in names]
        prng = random.Random(0x9_5D_57)
        calls += [(f, f2b(ref), ws) for f, ws, ref in prod5_families(prng, 2 * max(1, n // 5))
                  if "exact" in names or f in names]
        grng = random.Random(0x6A_77A)
        calls += [(f, ref, ws) for f, ws, ref in gamma_families(grng, max(1, n // 5))
                  if "exact" in names or f in names]
        crng = random.Random(0xC1_A05)
        calls += [(f, ref, ws) for f, ws, ref in cm_families(crng, max(1, n // 5))
                  if "exact" in names or f in names]
    if "ux" in names or any(x in UX for x in names):
        urng = random.Random(0x5C_2E_D6E)
        fams = ux_families(urng, max(1, n // 20))
        calls += [(f, want, ws) for f, ws, want in fams if "ux" in names or f in names]
    too_long = [name for name, _, a in calls if len(a) > 16]
    if too_long:                                # cr_driver.cyr reads at most sixteen words a call
        print(f"cr_verify: {len(too_long)} calls pass more than sixteen words ({too_long[0]})")
        return 1
    lines = [" ".join([name] + [w64(x) for x in a]) for name, _, a in calls]
    (root / "build/cr_in.txt").write_text("\n".join(lines) + "\n")
    targets = ["x86_64"] + (["aarch64"] if aarch64 else [])
    results = {t: run(t, lines) for t in targets}
    bad = 0
    per = {}
    for i, (name, want, a) in enumerate(calls):
        for t in targets:
            g = results[t][i]
            any_nan = name not in UX and not name.startswith(("gamma", "cm"))  # their "invalid" is a word
            ok = g == want or (any_nan and b2f(g) != b2f(g) and b2f(want) != b2f(want))
            per.setdefault((name, t), [0, 0])[0 if ok else 1] += 1
            if not ok:
                bad += 1
                if bad <= 20:
                    print(f"MISMATCH {t} {name}{tuple(a)}: got {g:016x}, want {want:016x}")
    for (name, t), (ok, ko) in sorted(per.items()):
        ref = "the exact value" if name in EXACT else ("MPFR (widest range)" if name in UX else "MPFR")
        print(f"  {name:9} {t:7} {ok + ko:8} args: {ok} equal to {ref}, {ko} differ")
    print(f"cr_verify: {len(calls)} calls x {len(targets)} targets, {bad} mismatches (MPFR {mpfr_ref.VERSION}; "
          "Fraction and isqrt for the exact families; MPFR's widest exponent range for the ux families)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
