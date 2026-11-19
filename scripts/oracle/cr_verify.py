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

from exact import (INF, binade, lw_below, lw_sign, rn_frac, rn_lw, rn_nb, rn_sqrt, rnd53,  # noqa: E402
                   rnd53_lw, rnd53_sqrt)


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
UX_FAR = 1 << 62       # its far stand-in's, beyond |x| = 2^47 (_bj_ux_exp_wide, M3 part 2d-3b-1)
UX_BIG = 2.0 ** 47


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
        # ux_exp_w (M2 part 2d-3b): the same with the stand-in beyond 2300, mostly between 1500 and 2300.
        xw = rng.uniform(1500.0, 2300.0) * sign() if i % 3 else rng.choice([x, rng.uniform(2300.001, 1e5) * sign()])
        pair("ux_exp_w", [xw], (1.0, UX_HUGE if xw > 0 else -UX_HUGE) if abs(xw) > 2300 else split53(W("exp", xw)))
        # ux_exp_big (M3 part 2d-3b-1): to |x| = 2^47 by its own reduction, every binade from 2^11 up, integers and
        # halves (k·ln 2 nowhere near them), around 2300 and 2^47, and the far stand-in beyond.
        cb = i % 5
        if cb == 0:
            xb = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(11, 46) * sign()
        elif cb == 1:
            xb = float(rng.randint(2300, 1 << 46)) * sign() + rng.choice([0.0, 0.5])
        elif cb == 2:
            xb = rng.choice([rng.uniform(2200.0, 2400.0), UX_BIG - rng.randint(0, 4096) * 2.0 ** -5]) * sign()
        elif cb == 3:
            xb = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(47, 60) * sign()
        else:
            xb = rng.uniform(-2300.0, 2300.0)
        if abs(xb) > UX_BIG:
            want = (1.0, UX_FAR if xb > 0 else -UX_FAR)
        else:
            want = mpfr_ref.wide_split("exp", xb)
        pair("ux_exp_big", [xb], want)
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
        # ux_log10 (M4 part 2d-4b): the same arguments, and the Touchstone writer's |S| just past DBL_MAX.
        calls.append(("ux_log10", [s1, e1], f2b(float(W("log10", v)))))
        s2 = rng.uniform(1.0, 1.4142135623730951)
        calls.append(("ux_log10", [s2, 1024], f2b(float(W("log10", Fraction(s2) * Fraction(2) ** 1024)))))
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
        if i % 5 == 1:
            ts, te = rng.uniform(1.0, 2.0) * sign(), rng.choice([1024, 1025])
        elif i % 5 == 3:                                      # the Green's function's k·R, up to 2^2052
            ts, te = rng.uniform(1.0, 2.0) * sign(), rng.randint(1024, 2052)
        elif i % 5 == 4:                                      # the HG field's phase, up to 2^4452 (M3 part 2d-3c)
            ts, te = rng.uniform(1.0, 2.0) * sign(), rng.choice([rng.randint(2053, 4452), rng.randint(4400, 4452)])
        elif i % 5 == 2:
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
         "prod5_acc", "mono_wide", "mono_wide_acc", "below_prod", "below_quot", "below_sum", "cm", "cm_slow", "cos_t",
         "cos_t_slow", "lw_reff", "lw_reff_slow", "lw_kappa", "lw_kappa_slow", "par", "par_x", "par_slow",
         "par_fx", "par_rep", "par_rep_x", "par_rep_slow", "par_rep_fx", "lw_nb", "lw_nb_slow", "lw_reff_ux",
         "vx_sqrt", "vx_sqrt_e",
         "vx_gamma", "lwo_reff", "lwo_reff_slow", "lwo_kappa", "lwo_kappa_slow", "lwo_nb", "lwo_nb_slow", "lwo_ux",
         "gammao", "gammao_slow", "gammao2", "gammao2_slow", "sumsq", "sumsq_exact", "abcd", "abcd_slow", "cascade",
         "cascade_slow", "s2z", "s2z_slow", "rat", "rat_x", "rat_slow", "dot4", "jcdot", "jint", "jnorm",
         "jover", "sjones", "sdop", "sellip", "sorient", "rsig", "rq", "rphase", "rang", "bspot", "bpeak", "broc",
         "bgouy", "bci", "bdot2", "bw0", "bgp", "bfsr", "bhg", "blg", "bhgf"]
UX = ["ux_exp", "ux_exp_w", "ux_exp_big", "ux_cosh", "ux_sinh", "ux_ln", "ux_log10", "ux_pow", "ux_sin", "ux_cos"]


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



# ── the trig-free Fresnel cosine (src/crmath_edge.cyr, M2 part 2c-3, ADR 0002 row 17) ─────────────────────────────────
COS_T_TIR = (1 << 64) - 1


def cos_t_ref(n1, n2, c):
    """RN(√(n₂² − n₁²(1 − c²))/n₂) exactly, or the TIR marker when that difference is negative."""
    d = Fraction(n2) ** 2 - Fraction(n1) ** 2 * (1 - Fraction(c) ** 2)
    if d < 0:
        return COS_T_TIR
    return f2b(0.0 if d == 0 else rn_sqrt(d / Fraction(n2) ** 2))


def cos_t_cases(rng, count):
    """(n₁, n₂, c): indices and angles of every kind; c within ulps of the critical cosine √(1 − (n₂/n₁)²), at scales
    2^±s; n₁ within ulps of n₂ with small and tiny c (Rust's cancellation; c below 2^-300 among them); any positive
    finite indices with any c in [0, 1]; ratios at the shortcuts' exponent differences (27, 28, −61, −62); the fast
    region's edge a²(1 − c²) = b²/2; and c = 0, c = 1, n₁ = n₂."""
    out = []
    for i in range(count):
        kind = i % 8
        if kind == 0:
            out.append((rng.uniform(1.0, 3.0), rng.uniform(1.0, 3.0), rng.random()))
        elif kind == 1:
            n2 = rng.uniform(1.0, 3.0)
            n1 = n2 * rng.choice([rng.uniform(1.0, 3.0), 1.0 + 2.0 ** -rng.uniform(1.0, 50.0),
                                  2.0 ** rng.uniform(0.0, 26.0)])
            r = n2 / n1
            cc = math.sqrt(1.0 - r * r)
            sc = 2.0 ** rng.randint(-990, 990) if rng.random() < 0.3 else 1.0
            if cc > 0.0:
                out.append((n1 * sc, n2 * sc, min(1.0, cm_ulps(cc, rng.randint(-8, 8)))))
        elif kind == 2:
            n2 = rng.uniform(0.5, 4.0)
            n1 = cm_ulps(n2, rng.randint(-6, 6))
            c = rng.choice([2.0 ** -rng.uniform(0.0, 60.0), 2.0 ** -rng.uniform(60.0, 1074.0), rng.random() * 1e-8])
            out.append((n1, n2, c))
        elif kind == 3:
            x1, x2 = abs(any_finite(rng)), abs(any_finite(rng))
            if x1 > 0.0 and x2 > 0.0 and math.isfinite(x1) and math.isfinite(x2):
                out.append((x1, x2, rng.choice([rng.random(), abs(any_finite(rng)) % 1.0, 1.0 - 2.0 ** -53])))
        elif kind == 4:
            d = rng.choice([27, 28, -61, -62, 26, -60])
            n2 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-500, 500)
            n1 = rng.uniform(1.0, 2.0) * n2 * 2.0 ** d
            out.append((n1, n2, rng.choice([1.0 - 2.0 ** -53, 1.0 - 2.0 ** -40, rng.random(), 0.0])))
        elif kind == 5:
            n2 = rng.uniform(1.0, 2.0)
            n1 = n2 * rng.uniform(0.75, 3.0)
            q = 1.0 - (n2 * n2) / (2.0 * n1 * n1)
            if q > 0.0:
                out.append((n1, n2, min(1.0, cm_ulps(math.sqrt(q), rng.randint(-8, 8)))))
        elif kind == 6:
            n = rng.uniform(0.1, 10.0)
            out.append(rng.choice([(n, n, rng.random()), (n, rng.uniform(0.1, 10.0), 0.0),
                                   (n, rng.uniform(0.1, 10.0), 1.0), (1e300, 1e-300, 1.0), (5e-324, 1e308, 0.5)]))
        else:
            n2 = abs(any_finite(rng))
            if n2 > 0.0 and math.isfinite(n2):
                n1 = n2 * rng.uniform(0.5, 2.0)
                if 0.0 < n1 < math.inf:
                    out.append((n1, n2, rng.random()))
    return out


def cos_t_families(rng, count):
    """cos_t (the full function) and cos_t_slow (the exact check alone) on each of cos_t_cases' triples."""
    out = []
    for n1, n2, c in cos_t_cases(rng, count):
        want = cos_t_ref(n1, n2, c)
        out.append(("cos_t", [n1, n2, c], want))
        out.append(("cos_t_slow", [n1, n2, c], want))
    return out


# ── the Liénard–Wiechert R_eff and κ (src/crmath_edge.cyr, M2 part 2c-4, ADR 0002 row 17) ────────────────────────────
LW_BELOW = (1 << 64) - 1                    # κ < 10^-30


def lw_valid(r, b) -> bool:
    """The callers' preconditions: finite, |R| ≥ 10^-30 and |β| < 1 exactly."""
    if not all(math.isfinite(x) for x in list(r) + list(b)):
        return False
    return (sum(Fraction(x) ** 2 for x in r) >= Fraction(1, 10 ** 60)
            and sum(Fraction(x) ** 2 for x in b) < 1)


def lw_unit(rng):
    """A random unit direction (doubles, not exactly unit)."""
    while True:
        v = [rng.gauss(0.0, 1.0) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if n > 1e-3:
            return [x / n for x in v]


def lw_knife(rng):
    """(R, β) with κ on either side of 10^-30, within a lattice step: β₁ = 1 − m·2^-53, β₂ leaving
    T₂ = 1 − β₁² − β₂² a little above D* = 2·10^-30, β₃ walked over its ulps (each moves κ by about 2^-54 of it) to
    the last points on either side of κ = 10^-30, R = β·2^s along β (κ = 1 − |β|), the parts permuted."""
    b1 = 1.0 - 2.0 ** -53 * rng.randint(1, 4)
    t1 = 1 - Fraction(b1) ** 2
    dstar = Fraction(2, 10 ** 30)
    b2 = rn_sqrt(t1 - dstar - Fraction(rng.randint(10, 90), 10 ** 32))
    t2 = t1 - Fraction(b2) ** 2
    if t2 <= dstar:
        return None
    base = f2b(rn_sqrt(t2 - dstar))
    s = 2.0 ** rng.randint(-60, 60)

    def side(k):
        bb = [b1, b2, b2f(base + k)]
        return lw_sign([x * s for x in bb], bb, 1, Fraction(1, 10 ** 30))
    lo, hi = -512, 512                          # κ falls as β₃ grows: side(lo) > 0 ≥ side(hi)
    if side(lo) <= 0 or side(hi) > 0:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if side(mid) > 0:
            lo = mid
        else:
            hi = mid
    k = rng.choice([lo, hi, lo - 1, hi + 1])
    bb = [b1, b2, b2f(base + k)]
    perm = rng.choice([(0, 1, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)])
    bb = [bb[j] for j in perm]
    return [x * s for x in bb], bb


def lw_tie(rng, mode):
    """(R, β) with the value exactly on a midpoint: R along a Pythagorean direction (1, b, c) of norm d, β₂ and β₃ on a
    coarse grid, and β₁ solving for the midpoint (None when β₁ would not be a double or |β| ≥ 1)."""
    _a, bq, cq, d = rng.choice([(1, 2, 2, 3), (1, 4, 8, 9), (1, 12, 12, 17), (1, 6, 18, 19), (1, 8, 4, 9)])
    sc = 2.0 ** rng.randint(-20, 20)
    b2 = rng.randint(-64, 64) / 1024.0
    b3 = rng.randint(-64, 64) / 1024.0
    if mode == 1:
        y = rng.uniform(0.05, 1.5)
    else:
        y = rng.uniform(0.05, 1.5) * d
    ulp = math.ulp(y)
    mu = Fraction(y) + Fraction(ulp) / 2
    target = d * (1 - mu) if mode == 1 else d - mu      # P/|R| = 1 − μ, or P = d − μ (R before scaling)
    b1 = target - bq * Fraction(b2) - cq * Fraction(b3)
    if float(b1) != b1 or abs(b1) >= 1:
        return None
    r = [sc, bq * sc, cq * sc]                           # R_eff scales with R: its midpoint too
    bb = [float(b1), b2, b3]
    if not lw_valid(r, bb):
        return None
    return r, bb


def lw_cases(rng, count):
    """(R, β): generic directions and speeds; near the cone's edge (β along R within 2^-k, |β| = 1 − 2^-e); an
    ultra-relativistic charge seen near its axis; κ's threshold knife edge; exact midpoints (ties); the near/far
    region's edge P = |R|/2; extreme scales and spreads (R past DBL_MAX's root, R near 10^-30, subnormal and zero
    parts); any finite parts."""
    out = []
    i = 0
    while len(out) < count:
        kind = i % 9
        i += 1
        r = b = None
        if kind == 0:
            r = [rng.uniform(-10.0, 10.0) for _ in range(3)]
            u = lw_unit(rng)
            sp = rng.random() ** 0.5
            b = [x * sp for x in u]
        elif kind == 1:
            n = lw_unit(rng)
            sc = 2.0 ** rng.randint(-90, 90)
            r = [x * sc for x in n]
            sp = 1.0 - 2.0 ** -rng.uniform(1.0, 52.0)
            e = lw_unit(rng)
            th = 2.0 ** -rng.uniform(0.0, 60.0)
            b = [(n[j] + th * e[j]) * sp for j in range(3)]
            m = math.sqrt(sum(x * x for x in b))
            b = [x * sp / m for x in b]
        elif kind == 2:
            x = 2.0 ** -rng.uniform(20.0, 40.0)
            y = 2.0 ** -rng.uniform(24.0, 30.0)
            r = [rng.choice([1.0, 3.0, 1e-20, 1e20]), x * rng.choice([1.0, -1.0]), 0.0]
            r[1] *= abs(r[0])
            b = [1.0 - 2.0 ** -53 * rng.randint(1, 8), y, 0.0]
            if rng.random() < 0.5:
                r = [r[1], r[2], r[0]]
                b = [b[1], b[2], b[0]]
        elif kind == 3:
            t = lw_knife(rng)
            if t is None:
                continue
            r, b = t
        elif kind == 4:
            t = lw_tie(rng, rng.choice([0, 1]))
            if t is None:
                continue
            r, b = t
        elif kind == 5:
            r = [rng.uniform(-3.0, 3.0) for _ in range(3)]
            n = math.sqrt(sum(x * x for x in r))
            f = 0.5 * (1.0 + rng.choice([1.0, -1.0]) * 2.0 ** -rng.uniform(1.0, 60.0))
            b = [x / n * f for x in r]
        elif kind == 6:
            r = rng.choice([
                [1e308, 1e308, -1e308], [1.7e308, 1.7e308, 1.7e308], [1e-30, 0.0, 0.0], [6e-31, 6e-31, 6e-31],
                [1e300, 1e-300, 5e-324], [5e-324, 1.0, 0.0], [1.0, 5e-324, 5e-324], [2.0 ** -100, 2.0 ** 1000, 1.0]])
            r = [x * rng.choice([1.0, -1.0]) for x in r]
            b = rng.choice([[0.0, 0.0, 0.0], [5e-324, 0.0, -5e-324], [0.999, 1e-200, 0.0],
                            [1.0 - 2.0 ** -53, 0.0, 0.0], [-0.6, 0.7, 0.3], [0.0, 1.0 - 2.0 ** -53, 2.0 ** -27]])
            if rng.random() < 0.5:
                b = [r[j] / math.sqrt(3.0) * 0.9 if 0 < abs(r[j]) < 1e300 else b[j] for j in range(3)]
        elif kind == 7:
            r = [any_finite(rng) for _ in range(3)]
            u = lw_unit(rng)
            sp = rng.random()
            b = [x * sp for x in u]
            if rng.random() < 0.3:
                b[rng.randrange(3)] = any_finite(rng) % 0.5
        else:
            r = [any_finite(rng) for _ in range(3)]
            n = math.sqrt(sum((x / 2.0 ** 600) ** 2 for x in r)) if max(abs(x) for x in r) > 1e200 else \
                math.sqrt(sum(x * x for x in r))
            if n == 0.0 or not math.isfinite(n):
                continue
            sp = 1.0 - 2.0 ** -rng.uniform(1.0, 50.0)
            big = max(abs(x) for x in r) > 1e200
            b = [(x / 2.0 ** 600 if big else x) / n * sp for x in r]
        if r is not None and lw_valid(r, b):
            out.append((r, b))
    return out


def lw_families(rng, count):
    """lw_reff and lw_kappa (the full function) and their *_slow (the exact checks alone) on each of lw_cases' pairs."""
    out = []
    for r, b in lw_cases(rng, count):
        reff = f2b(rn_lw(r, b, 0))
        kap = LW_BELOW if lw_below(r, b) else f2b(rn_lw(r, b, 1))
        out.append(("lw_reff", r + b, reff))
        out.append(("lw_reff_slow", r + b, reff))
        out.append(("lw_kappa", r + b, kap))
        out.append(("lw_kappa_slow", r + b, kap))
    return out


def lw_ux_families(rng, count):
    """lw_reff_ux (M2 part 2d-3e): _bj_cr_lw's mode 3, R_eff to 53 bits with an unbounded exponent as its significand
    (part 0) and exponent (part 1), against exact.py's rnd53_lw: on lw_cases' pairs and on them scaled by 2^k toward
    DBL_MAX (R_eff past it), where R stays finite."""
    out = []
    M = (1 << 64) - 1
    for i, (r, b) in enumerate(lw_cases(rng, count)):
        rs = [r]
        mx = max(abs(x) for x in r)
        if mx > 0 and i % 2 == 0:
            sh = 1023 - binade(Fraction(mx)) - rng.randint(0, 3)
            if sh > 0:
                try:
                    rr = [math.ldexp(x, sh) for x in r]
                except OverflowError:
                    rr = []
                if rr and all(Fraction(x) == Fraction(y) * Fraction(2) ** sh for x, y in zip(rr, r)):
                    rs.append(rr)
        for rv in rs:
            if lw_below(rv, b) or all(x == 0 for x in rv):
                continue
            a = rnd53_lw(rv, b)
            e = binade(a)
            out.append(("lw_reff_ux", rv + b + [0], f2b(float(a / Fraction(2) ** e))))
            out.append(("lw_reff_ux", rv + b + [1], e & M))
    return out


# ── the parallel combinations (src/crmath_edge.cyr, M2 part 2c-5, ADR 0002 row 17) ───────────────────────────────────
PAR_BELOW = (1 << 64) - 1                   # |Σ 1/xᵢ| < 10^-30
E30F = Fraction(1, 10 ** 30)


def par_ref(xs, rep=1):
    """RN(1/Σ 1/xᵢ) exactly (an ±∞ adds 0), the values repeated rep times, or the marker when |Σ| < 10^-30."""
    s = rep * sum(Fraction(1) / Fraction(x) for x in xs if math.isfinite(x))
    if abs(s) < E30F:
        return PAR_BELOW
    return f2b(rn_frac(1 / s))


def par_big(rng):
    """A finite double with |x| ≥ 10^-30 (the callers' cut), any binade above it, either sign."""
    while True:
        x = any_finite(rng)
        if math.isfinite(x) and abs(x) >= 1e-30:
            return x


def par_triple(rng):
    """Three values whose reciprocals cancel exactly: a, −b, −ab/2^j for odd a near 2^25.5 and b = a + 2^j (1/a − 1/b =
    2^j/(ab)), scaled by 2^t so that the reciprocals reach near 1 (2.0.1, H18: Σ|1/x| past the double-double bounds)."""
    a = 2 * rng.randint(2 ** 24, 2 ** 25 - 1) + 1
    j = rng.randint(1, 20)
    b = a + 2 ** j
    t = rng.randint(-28, -20)
    return [math.ldexp(a, t), -math.ldexp(b, t), -math.ldexp(a * b, t - j)]


def par_cases(rng, count):
    """Value lists (n ≤ 15): physical positive values; mixed signs; pairs and triples that cancel to ulps; |Σ 1/x|
    within ulps of 10^-30 (one value near 10^30, k equal values near k·10^30, a pair a, −(a + j ulps) near 3e14); exact
    midpoints of 1/S (M₁ and 2^28 − M₁ for odd M₁ near 2^27, at scales 2^±s, either sign order); extremes (±DBL_MAX,
    1e-30 itself, ±∞ among them, all ∞); repeated values with alternating signs; any finite values; up to four exactly
    cancelling triples (par_triple), alone or with one more value v, whose combination is then v itself (2.0.1)."""
    out = []
    i = 0
    while len(out) < count:
        kind = i % 9
        i += 1
        xs = None
        if kind == 0:
            xs = [rng.uniform(0.1, 1000.0) * 10.0 ** rng.randint(-6, 9) for _ in range(rng.randint(1, 15))]
        elif kind == 1:
            xs = [rng.uniform(0.1, 10.0) * rng.choice([1.0, -1.0]) for _ in range(rng.randint(2, 15))]
        elif kind == 2:
            a = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 60)
            c = rng.choice([0, 1, 2])
            if c == 0:
                xs = [a, -cm_ulps(a, rng.randint(-4, 4))]
            elif c == 1:
                b = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 60)
                t = -1.0 / (1.0 / a + 1.0 / b)
                xs = [a, b, cm_ulps(abs(t), rng.randint(-4, 4)) * -1.0]
            else:
                xs = [a, -a, rng.uniform(1e29, 1e31)]
        elif kind == 3:
            c = rng.choice([0, 1, 2])
            if c == 0:
                xs = [cm_ulps(1e30, rng.randint(-8, 8)) * rng.choice([1.0, -1.0])]
            elif c == 1:
                k = rng.randint(2, 15)
                xs = [cm_ulps(k * 1e30, rng.randint(-8, 8))] * k
            else:
                a = rng.uniform(1.0, 9.0) * 1e14
                xs = [a, -cm_ulps(a, rng.randint(1, 40))]
        elif kind == 4:
            m1 = 2 * rng.randint(2 ** 25 + 2 ** 24, 2 ** 26 + 2 ** 25 - 1) + 1
            sc = 2.0 ** rng.randint(-60, 60)
            xs = [m1 * sc, (2 ** 28 - m1) * sc]
            if rng.random() < 0.5:
                xs = xs[::-1]
            if rng.random() < 0.3:
                xs = [-x for x in xs]
        elif kind == 5:
            pool = [1.7976931348623157e308, -1.7976931348623157e308, 1e-30, -1e-30, 1e300, math.inf, -math.inf,
                    1.0, 3.0, -3.0]
            xs = [rng.choice(pool) for _ in range(rng.randint(1, 6))]
            if rng.random() < 0.1:
                xs = [math.inf] * rng.randint(1, 4)
        elif kind == 6:
            v = rng.choice([3.0, 0.1, 7.0, 1e20, 5e-30])
            xs = [v if j % 2 == 0 else -v for j in range(rng.randint(2, 15))]
        elif kind == 7:
            xs = [par_big(rng) for _ in range(rng.randint(1, 15))]
        else:
            xs = [x for _ in range(rng.randint(1, 4)) for x in par_triple(rng)]
            if rng.random() < 0.5:
                xs.insert(rng.randrange(len(xs) + 1), rng.uniform(0.5, 2.0) * 2.0 ** rng.randint(-40, 40))
            if rng.random() < 0.5:
                xs = xs[::-1]
        if xs and all(math.isinf(x) or abs(x) >= 1e-30 for x in xs):
            out.append(xs)
    return out


def par_families(rng, count):
    """par (every tier), par_x (the expansion, then the exact path), par_slow (the exact path alone) and par_fx (the
    fixed-point tier, then the exact path; 2.0.1) on par_cases' lists, and par_rep / par_rep_x / par_rep_slow /
    par_rep_fx on short patterns repeated up to 20,000 times (sums that cancel to 0 or to ulps, knife edges at 10^-30,
    pairs of exactly cancelling triples)."""
    out = []
    for xs in par_cases(rng, count):
        want = par_ref(xs)
        ws = [len(xs)] + xs
        out.append(("par", ws, want))
        out.append(("par_x", ws, want))
        out.append(("par_slow", ws, want))
        out.append(("par_fx", ws, want))
    for _ in range(max(1, count // 200)):
        m = rng.choice([2, 10, 1000, 20000])
        c = rng.randrange(6)
        if c == 5:
            pat = par_triple(rng) + par_triple(rng)
        elif c == 0:
            pat = [3.0, -3.0]
        elif c == 1:
            pat = [10.0, 20.0, 30.0, 40.0, 50.0]
        elif c == 2:
            pat = [3.0, -cm_ulps(3.0, rng.randint(1, 3))]
        elif c == 3:
            pat = [cm_ulps(m * 1e30, rng.randint(-4, 4))]
        else:
            pat = [rng.uniform(1.0, 2.0) * rng.choice([1.0, -1.0]) for _ in range(rng.randint(1, 6))]
        want = par_ref(pat, m)
        ws = [m, len(pat)] + pat
        out.append(("par_rep", ws, want))
        out.append(("par_rep_x", ws, want))
        out.append(("par_rep_slow", ws, want))
        out.append(("par_rep_fx", ws, want))
    return out


# ── M3's exact sum-of-squares threshold (src/crmath_cmp.cyr, M3 part 2c, ADR 0002 row 17 amended) ──────────────────
def sumsq_ref(a, b, c, e, f, p):
    """The sign of (a·b + c)² + (e·f)² − 10^p as the word of −1.0, 0.0 or 1.0, exactly."""
    v = (Fraction(a) * Fraction(b) + Fraction(c)) ** 2 + (Fraction(e) * Fraction(f)) ** 2 - Fraction(10) ** p
    return f2b(-1.0 if v < 0 else 1.0 if v > 0 else 0.0)


def sumsq_elim(p):
    """The least elim with 2^elim ≥ 2·√(10^p): the exact test's bound on |X| and |Y|."""
    k = math.ceil(1 + p * math.log2(10) / 2) - 2
    while Fraction(2) ** (2 * k) < 4 * Fraction(10) ** p:
        k += 1
    return k


def sumsq_ulps(x, k):
    return b2f(f2b(x) + k) if x != 0.0 else x


def sumsq_families(rng, count):
    """sumsq (the band on the rounded values, then the exact test) and sumsq_exact (the exact test alone): operands
    whose sum of squares lies within a few ulps of 10^p (p = −60 and 14, propagate_beam's and Mie's, and any p in
    [−60, 30]), X = a·b + c cancelling to its product's rounding error, subnormal and tiny partners, exact ties
    (Pythagorean X and Y, 10^p = X² + Y²) and their one-ulp neighbours, Mie-shaped (a·x, e·x) and propagate-shaped
    (C·q₀ + D, C·q₁) pairs, and far operands from any binade."""
    out = []
    for i in range(count):
        mode = i % 7
        p = (-60, 14)[i % 2] if i % 3 else rng.randint(-60, 30)
        root = math.sqrt(10.0 ** p)
        th = rng.uniform(0.0, math.pi / 2)
        xt, yt = root * math.cos(th), root * math.sin(th)
        sg = rng.choice([1.0, -1.0])
        if mode == 0:                                   # near the threshold, no cancellation
            a = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-400, 400) * sg
            b = sumsq_ulps(xt / a, rng.randint(-3, 3))
            c = 0.0 if i % 4 else xt * 2.0 ** -rng.randint(1, 60) * rng.choice([1.0, -1.0])
            e = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-400, 400)
            f = sumsq_ulps(yt / e, rng.randint(-3, 3)) if rng.random() < 0.8 else 0.0
        elif mode == 1:                                 # X = a·b + c is the product's rounding error, near √(10^p)
            j = rng.randint(0, 40)
            a = rng.uniform(1.0, 2.0) * sg
            b = rng.uniform(1.0, 2.0) * xt * 2.0 ** (52 + j) / a
            h = a * b
            c = sumsq_ulps(-h, rng.choice([0, 0, 0, 1, -1]))
            lo = Fraction(a) * Fraction(b) + Fraction(c)
            rest = Fraction(10) ** p - lo * lo
            yv = math.sqrt(float(rest)) if rest > 0 else 0.0
            e = rng.uniform(1.0, 2.0)
            f = sumsq_ulps(yv / e, rng.randint(-2, 2)) if yv else 0.0
        elif mode == 2:                                 # subnormal and tiny partners
            tiny = [5e-324, 3e-321, 2.2250738585072014e-308, 1e-300, 1e-200]
            c = sumsq_ulps(xt * sg, rng.randint(-3, 3))
            a = rng.choice(tiny) * rng.randint(1, 9) * rng.choice([1.0, -1.0])
            b = rng.choice([5e-324, 1e-310, 1e-160, 1.0, 1e100])
            if rng.random() < 0.5:
                e, f = rng.choice(tiny) * rng.randint(1, 9), rng.choice(tiny) * rng.randint(1, 9)
                c = sumsq_ulps(root * sg, rng.randint(-3, 3))
            else:
                e = rng.uniform(1.0, 2.0) * 2.0 ** -1050
                f = sumsq_ulps(yt / e, rng.randint(-2, 2)) if math.isfinite(yt / e) else 1.0
        elif mode == 3:                                 # far: any binade, either sign
            a, b, c, e, f = (any_finite(rng) for _ in range(5))
        elif mode == 4:                                 # exact ties for p ≥ 0, and neighbours
            p = rng.randint(0, 30)
            h = p // 2
            if p % 2:
                xs, ys = 10.0 ** h, 3.0 * 10.0 ** h
            elif p >= 2 and rng.random() < 0.6:
                xs, ys = 6.0 * 10.0 ** (h - 1), 8.0 * 10.0 ** (h - 1)
            else:
                xs, ys = 10.0 ** h, 0.0
            if rng.random() < 0.5:
                xs, ys = ys, xs
            b = 2.0 ** rng.randint(-20, 20)
            a, c = xs / b * sg, 0.0
            if rng.random() < 0.3:
                a, c = xs / 2.0 / b * sg, xs / 2.0 * sg
            f = 2.0 ** rng.randint(-20, 20)
            e = ys / f * rng.choice([1.0, -1.0])
            k = rng.choice([0, 0, 0, 1, -1])
            if k:
                which = rng.randrange(3)
                if which == 0 and a:
                    a = sumsq_ulps(a, k)
                elif which == 1 and e:
                    e = sumsq_ulps(e, k)
                else:
                    b = sumsq_ulps(b, k)
        elif mode == 5:                                 # Mie: (Re m·x)² + (Im m·x)² against 10^14
            p = 14
            x = rng.choice([1e5, 1e5, rng.uniform(1.0, 1e5)])
            mr, mi = 1e7 / x * math.cos(th), -1e7 / x * math.sin(th)
            if rng.random() < 0.3:
                mr, mi = 1e7 / x * 0.6, -1e7 / x * 0.8
            a, b, c = sumsq_ulps(mr, rng.randint(-2, 2)), x, 0.0
            e, f = sumsq_ulps(mi, rng.randint(-2, 2)), x
        else:                                           # propagate_beam: X = C·q₀ + D, Y = C·q₁ against 10^-60
            p = -60
            cm = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 60) * sg
            q0 = rng.uniform(-1.0, 1.0) * 2.0 ** rng.randint(-60, 60)
            d = sumsq_ulps(-cm * q0 + xt * rng.choice([0.0, 1.0, -1.0]), rng.randint(-2, 2))
            yv = 1e-30 if rng.random() < 0.5 else yt
            a, b, c, e, f = cm, q0, d, cm, sumsq_ulps(yv / cm, rng.randint(-2, 2))
        if not all(math.isfinite(v) for v in (a, b, c, e, f)):
            continue
        elim = sumsq_elim(p)
        t = float(Fraction(10) ** p)
        want = sumsq_ref(a, b, c, e, f, p)
        out.append(("sumsq", [a, b, c, e, f, t, p, elim], want))
        out.append(("sumsq_exact", [a, b, c, e, f, p, elim], want))
    return out


# ── M3's ABCD propagation (src/crmath_edge.cyr, M3 part 2c-2, ADR 0002 row 17 amended) ──────────────────────────────
ABCD_REFUSED = (1 << 64) - 1                # parts 0 and 1 when |C·q + D|² < 10^-60


def abcd_ref(v):
    """(RN(Re q_out), RN(Im q_out), flags) as words for q_out = (A·q + B)/(C·q + D), exactly; the refused marker for
    |C·q + D|² < 10^-60."""
    a, b, c, d, q0, q1 = (Fraction(x) for x in v)
    x, y = c * q0 + d, c * q1
    den = x * x + y * y
    if den < Fraction(1, 10 ** 60):
        return ABCD_REFUSED, ABCD_REFUSED, 0
    nr = (a * q0 + b) * x + (a * q1) * y
    ni = q1 * (a * d - b * c)
    re = 0.0 if nr == 0 else rn_frac(nr / den)
    im = 0.0 if ni == 0 else rn_frac(ni / den)
    return f2b(re), f2b(im), 1 + (2 if nr == 0 else 0) + (4 if ni == 0 else 0)


def abcd_compose(m, n):
    """m·n as Rust's compose rounds it (n first)."""
    return [m[0] * n[0] + m[1] * n[2], m[0] * n[1] + m[1] * n[3], m[2] * n[0] + m[3] * n[2],
            m[2] * n[1] + m[3] * n[3]]


def abcd_cases(rng, count):
    """(A, B, C, D, q₀, q₁): composed optical systems and Gaussian q's, free space with sums that tie, cancellations
    at the pole C·q₀ + D → 0 at every depth, the 10^-60 threshold, constructed midpoints of Im q_out = det/q₁ and their
    neighbours, exact zeros, entries at the fast path's box edges, scaled systems, any finite values, a lens followed
    by the distance to its focus (Re q_out cancelling in U·X + V·Y), near-singular matrices (AD − BC cancelling) and
    C = 0 (free space and interfaces: any D, A = ±1, q₀ = 0, ties)."""
    out = []
    for i in range(count):
        mode = i % 12
        sg = rng.choice([1.0, -1.0])
        if mode == 0:                                   # optical systems
            m = [1.0, 0.0, 0.0, 1.0]
            for _ in range(rng.randint(1, 4)):
                k = rng.randrange(3)
                if k == 0:
                    e = [1.0, rng.uniform(1e-3, 2.0), 0.0, 1.0]
                elif k == 1:
                    e = [1.0, 0.0, -1.0 / (rng.uniform(0.01, 1.0) * rng.choice([1.0, -1.0])), 1.0]
                else:
                    e = [1.0, 0.0, 0.0, rng.uniform(0.5, 2.0)]
                m = abcd_compose(e, m)
            q = [rng.uniform(-2.0, 2.0) if rng.random() < 0.7 else 0.0, rng.uniform(1e-4, 10.0)]
            v = m + q
        elif mode == 1:                                 # free space: Re q_out = (A·q₀ + B)/D, ties included
            j = rng.randint(-60, 60)
            q0 = rng.choice([1.0, 3.0, rng.uniform(1.0, 2.0)]) * 2.0 ** j
            b = rng.choice([2.0 ** (j - 53), 3 * 2.0 ** (j - 54),
                            rng.uniform(-1.0, 1.0) * 2.0 ** (j - rng.randint(0, 60))])
            a = rng.choice([1.0, 2.0, 0.5, rng.uniform(0.5, 2.0)])
            d = rng.choice([1.0, 2.0, 0.25, rng.uniform(0.5, 2.0)])
            v = [a, b * sg, 0.0, d, q0, rng.uniform(1e-3, 10.0)]
        elif mode == 2:                                 # C·q₀ + D cancelling to its rounding error or a few ulps
            c = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-30, 30) * sg
            q0 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-30, 30) * rng.choice([1.0, -1.0])
            d = b2f(f2b(-(c * q0)) + rng.randint(-4, 4)) if c * q0 != 0 else 1.0
            q1 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-120, 0)
            v = [rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), c, d, q0, q1]
        elif mode == 3:                                 # |C·q + D|² near 10^-60
            th = rng.uniform(0.0, math.pi / 2)
            xt, yt = 1e-30 * math.cos(th), 1e-30 * math.sin(th)
            c = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-20, 20)
            q0 = rng.uniform(-1.0, 1.0)
            d = b2f(f2b(xt - c * q0) + rng.randint(-2, 2)) if xt - c * q0 else xt
            q1 = b2f(f2b(yt / c) + rng.randint(-2, 2))
            v = [rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), c, d, q0, q1]
        elif mode == 4:                                 # Im q_out = det/q₁ at or near a midpoint (X = 0, C = 1)
            j = rng.randint(-40, 40)
            k = rng.choice([0, 0, 1, -1, 2])
            bb = -2.0 ** -53 * rng.choice([1.0, 3.0, 5.0])
            v = [1.0, b2f(f2b(bb) + k), 1.0, 1.0, -1.0, 2.0 ** j * sg]
        elif mode == 5:                                 # any finite values
            v = [any_finite(rng) for _ in range(6)]
        elif mode == 6:                                 # a system and q scaled past the box, toward the range ends
            m = [rng.uniform(0.5, 2.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(0.5, 2.0)]
            s1, s2 = rng.randint(-1100, 1000), rng.randint(-1100, 1000)
            v = [m[0], m[1] * 2.0 ** s1 if s1 > -1074 else m[1], m[2] * 2.0 ** -s1 if -s1 < 1000 else m[2], m[3],
                 rng.uniform(-1.0, 1.0) * 2.0 ** max(s2, -1074), rng.uniform(0.0, 1.0) * 2.0 ** max(s2, -1074)]
            if rng.random() < 0.3:
                v[rng.randrange(6)] = rng.choice([2.0 ** 200, 2.0 ** -200, 2.0 ** 201, 2.0 ** -201, 5e-324,
                                                  1.7976931348623157e308])
        elif mode == 7:                                 # exact zeros: U = 0, AD = BC, q₁ = 0
            a = float(rng.randint(1, 9)) * sg
            q0 = float(rng.randint(-9, 9)) * 2.0 ** rng.randint(-4, 4)
            c = float(rng.randint(-9, 9))
            d = float(rng.randint(1, 9))
            b = -a * q0
            z = rng.randrange(3)
            if z == 1:
                b = a * d / c if c else 1.0
            q1 = 0.0 if z == 2 else float(rng.randint(1, 9)) * rng.choice([1.0, -1.0, 0.0])
            v = [a, b, c, d, q0, q1]
        elif mode == 9:                                 # Re q_out at a waist: a lens, then the distance to the focus
            f = rng.uniform(0.01, 1.0) * rng.choice([1.0, -1.0])
            q = [rng.uniform(-0.5, 0.5) if rng.random() < 0.5 else 0.0, rng.uniform(1e-4, 1.0)]
            lens = [1.0, 0.0, -1.0 / f, 1.0]
            mq = abcd_ref(lens + q)
            dist = -b2f(mq[0]) if mq[2] else 0.0
            dist = b2f(f2b(dist) + rng.randint(-3, 3)) if dist else dist
            v = abcd_compose([1.0, dist, 0.0, 1.0], lens) + q
        elif mode == 11:                                # C = 0: free space and interfaces, any D, A = ±1, ties
            j = rng.randint(-50, 50)
            q0 = rng.choice([0.0, 1.0, rng.uniform(1.0, 2.0), 3.0]) * 2.0 ** j * rng.choice([1.0, -1.0])
            a = rng.choice([1.0, -1.0, rng.uniform(0.25, 4.0) * sg])
            b = rng.choice([2.0 ** (j - 53), -(2.0 ** (j - 54)) * 3, rng.uniform(-1.0, 1.0) * 2.0 ** j,
                            -a * q0, 0.0])
            d = rng.choice([1.0, -1.0, 2.0, 1.5, rng.uniform(-3.0, 3.0), 1e-30, 9.999999999999999e-31])
            q1 = rng.choice([0.0, rng.uniform(-2.0, 2.0) * 2.0 ** rng.randint(-60, 60)])
            v = [a, b, 0.0, d, q0, q1]
        elif mode == 10:                                # near-singular: AD − BC cancelling to a few ulps of AD
            a = rng.uniform(0.5, 2.0) * 2.0 ** rng.randint(-20, 20) * sg
            b = rng.uniform(-2.0, 2.0) * 2.0 ** rng.randint(-20, 20)
            c = rng.uniform(-2.0, 2.0) * 2.0 ** rng.randint(-20, 20)
            d = b * c / a
            d = b2f(f2b(d) + rng.randint(-3, 3)) if d else d
            v = [a, b, c, d, rng.uniform(-1.0, 1.0), rng.uniform(1e-3, 2.0)]
        else:                                           # a cancellation factor of 2^k in X, around the fast budget
            c = rng.uniform(1.0, 2.0)
            q0 = rng.uniform(1.0, 2.0)
            kk = rng.randint(0, 30)
            d = -(c * q0) * (1.0 - 2.0 ** -kk * rng.uniform(0.5, 1.0))
            v = [rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), c, d, q0, rng.uniform(1e-6, 1.0) * 2.0 ** -kk]
        if all(math.isfinite(x) for x in v):
            out.append(v)
    return out


def abcd_families(rng, count):
    """abcd (the fast path, then the exact one) and abcd_slow (the exact path alone): each case's three parts."""
    out = []
    for v in abcd_cases(rng, count):
        ref = abcd_ref(v)
        for part in range(3):
            out.append(("abcd", v + [part], ref[part]))
            out.append(("abcd_slow", v + [part], ref[part]))
    return out


# ── M4's S-parameter cascade (src/crmath_edge.cyr, M4 part 2a-2, ADR 0002 row 17 amended) ──────────────────────────
CASC_REFUSED = (1 << 64) - 1                # parts 0 to 7 when |T22|² < 10^-60


def _cmul(x, y):
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def cascade_ref(v):
    """The eight parts of Redheffer's star product of a = v[0:8] then b = v[8:16] (S11, S12, S21, S22 as (re, im)),
    each RN of the exact value, then the return (1 + 4·m, bit k of m set where part k is exactly 0); the refused
    marker and 0 when 10^60·|D|² < |A21|²·|B21|²."""
    f = [Fraction(x) for x in v]
    a = [(f[2 * k], f[2 * k + 1]) for k in range(4)]
    b = [(f[8 + 2 * k], f[9 + 2 * k]) for k in range(4)]
    p = _cmul(a[3], b[0])
    d = (1 - p[0], -p[1])
    den = d[0] * d[0] + d[1] * d[1]
    if 10 ** 60 * den < (a[2][0] ** 2 + a[2][1] ** 2) * (b[2][0] ** 2 + b[2][1] ** 2):
        return [CASC_REFUSED] * 8 + [0]

    zeros = []

    def quot(pp, a0=(0, 0)):
        nr = a0[0] * den + pp[0] * d[0] + pp[1] * d[1]
        ni = a0[1] * den + pp[1] * d[0] - pp[0] * d[1]
        zeros.extend([nr == 0, ni == 0])
        return [f2b(0.0 if nr == 0 else rn_frac(nr / den)), f2b(0.0 if ni == 0 else rn_frac(ni / den))]
    out = quot(_cmul(_cmul(a[1], a[2]), b[0]), a[0]) + quot(_cmul(a[1], b[1])) + quot(_cmul(a[2], b[2]))
    out += quot(_cmul(_cmul(a[3], b[1]), b[2]), b[3])
    return out + [1 + 4 * sum(1 << k for k, z in enumerate(zeros) if z)]


def cascade_cases(rng, count):
    """Two 2-ports' 16 parts: passive networks; resonances (A22·B11 → 1 to 2^-60 and past); the 10^-60 threshold
    straddled ulp by ulp; exact zeros and thrus; S11's numerator cancelling; the fast path's box edges at 2^±120;
    any finite values (10^±300) and subnormals."""
    def polar(r, th):
        return [r * math.cos(th), r * math.sin(th)]

    def port(mag=1.0):
        return sum((polar(rng.uniform(0.0, mag), rng.uniform(-math.pi, math.pi)) for _ in range(4)), [])
    out = []
    for i in range(count):
        mode = i % 8
        a, b = port(), port()
        if mode == 1:                                       # a resonance: A22·B11 = 1 − δ
            r, th = rng.uniform(0.5, 1.0), rng.uniform(-math.pi, math.pi)
            a[6:8] = polar(r, th)
            k = rng.randint(0, 70)
            z = complex(*a[6:8])
            w = (1.0 - complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) * 2.0 ** -k) / z
            b[0:2] = [w.real, w.imag]
            if rng.random() < 0.5:
                a[4:6] = polar(math.sqrt(1 - r * r), rng.uniform(-math.pi, math.pi))
        elif mode == 2:                                     # the threshold: D = d exactly against 10^-30·|S21·S21'|
            kd = rng.randint(1, 80)
            dv = 2.0 ** -kd * rng.choice([1.0, 1.5, 1.25])
            a[6:8], b[0:2] = [1.0, 0.0], [1.0 - dv, 0.0]
            x = float(Fraction(dv) * 10 ** 30)
            for _ in range(rng.randint(-3, 3) + 3):
                x = math.nextafter(x, math.inf if rng.random() < 0.5 else 0.0)
            a[4:6], b[4:6] = [1.0, 0.0], [x, 0.0] if rng.random() < 0.5 else [0.0, x]
        elif mode == 3:                                     # zeros and thrus
            a = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0] if rng.random() < 0.5 else a
            for j in rng.sample(range(16), rng.randint(1, 8)):
                (a if j < 8 else b)[j % 8] = 0.0
            if a[4] == a[5] == 0.0:
                a[4] = 0.5
            if b[4] == b[5] == 0.0:
                b[4] = 0.5
        elif mode == 4:                                     # S11's numerator cancelling: A11 ≈ −A12·A21·B11/D
            z = complex(*a[2:4]) * complex(*a[4:6]) * complex(*b[0:2]) / (1 - complex(*a[6:8]) * complex(*b[0:2]))
            a[0:2] = [-z.real, -z.imag]
        elif mode == 5:                                     # the box edges
            for j in rng.sample(range(16), rng.randint(1, 4)):
                e = rng.choice([120, -120, 121, -121, 119, -119])
                s = rng.choice([1.0, -1.0, 1.5, -1.75])
                (a if j < 8 else b)[j % 8] = s * 2.0 ** e
        elif mode == 6:                                     # any finite values
            for j in range(16):
                (a if j < 8 else b)[j % 8] = rng.choice([1.0, -1.0]) * 10.0 ** rng.uniform(-300, 300)
        elif mode == 7:                                     # subnormal and tiny parts
            for j in rng.sample(range(16), rng.randint(1, 6)):
                (a if j < 8 else b)[j % 8] = rng.choice([1.0, -1.0]) * 2.0 ** rng.randint(-1074, -1000)
        v = a + b
        assert len(v) == 16
        if all(math.isfinite(x) for x in v):
            out.append(v)
    return out


def cascade_families(rng, count):
    """cascade (the fast path, then the exact one) and cascade_slow (the exact path alone): each case's nine parts."""
    out = []
    for v in cascade_cases(rng, count):
        ref = cascade_ref(v)
        for part in range(9):
            out.append(("cascade", v + [part], ref[part]))
            out.append(("cascade_slow", v + [part], ref[part]))
    return out


# ── M4's S → Z conversion (src/crmath_edge.cyr, M4 part 2a-3, ADR 0002 row 17 amended) ────────────────────────────
def s2z_ref(v):
    """The eight parts of Z = Z₀(I + S)(I − S)⁻¹ for S = v[0:8] (S11, S12, S21, S22 as (re, im)) and Z₀ = v[8], each RN
    of the exact value, then the return (1 + 4·m, bit k of m set where part k is exactly 0); the refused marker and 0
    when |det(I − S)|² < 10^-60."""
    f = [Fraction(x) for x in v]
    a, b, c, d = [(f[2 * k], f[2 * k + 1]) for k in range(4)]
    z0 = f[8]
    ad, bc = _cmul(a, d), _cmul(b, c)
    det = (1 - a[0] - d[0] + ad[0] - bc[0], -a[1] - d[1] + ad[1] - bc[1])
    den = det[0] ** 2 + det[1] ** 2
    if 10 ** 60 * den < 1:
        return [CASC_REFUSED] * 8 + [0]
    zeros = []

    def quot(n):
        nr = z0 * (n[0] * det[0] + n[1] * det[1])
        ni = z0 * (n[1] * det[0] - n[0] * det[1])
        zeros.extend([nr == 0, ni == 0])
        return [f2b(0.0 if nr == 0 else rn_frac(nr / den)), f2b(0.0 if ni == 0 else rn_frac(ni / den))]
    n11 = (1 + a[0] - d[0] - ad[0] + bc[0], a[1] - d[1] - ad[1] + bc[1])
    n22 = (1 - a[0] + d[0] - ad[0] + bc[0], -a[1] + d[1] - ad[1] + bc[1])
    out = quot(n11) + quot((2 * b[0], 2 * b[1])) + quot((2 * c[0], 2 * c[1])) + quot(n22)
    return out + [1 + 4 * sum(1 << k for k, z in enumerate(zeros) if z)]


def s2z_cases(rng, count):
    """(S11, S12, S21, S22, Z₀): passive networks; near-thrus (det → 0 to 2^-70); |det| straddling 10^-30 ulp by ulp
    (det = (1 − S11)(1 − S22) − S21 with S11 = S22 = 1 + 2^-52·k and S12 = 1); lossless networks (Re Z ≈ 0); exact
    zeros; the box edges at 2^±120; any finite values; subnormals."""
    def polar(r, th):
        return [r * math.cos(th), r * math.sin(th)]
    out = []
    for i in range(count):
        mode = i % 8
        s = sum((polar(rng.uniform(0.0, 1.0), rng.uniform(-math.pi, math.pi)) for _ in range(4)), [])
        z0 = rng.choice([50.0, 75.0, 1.0, rng.uniform(0.1, 1000.0)])
        if mode == 1:                                       # near-thru: S12·S21 → (1 − S11)(1 − S22)
            k = rng.randint(0, 70)
            e = [rng.uniform(-1, 1) * 2.0 ** -k for _ in range(4)]
            s = [e[0], e[1], 1.0, 0.0, 1.0 - abs(e[2]), e[3], -e[0], e[1]]
        elif mode == 2:                                     # |det| at 10^-30: P − c with P = (2^-52·k)²
            kk = rng.randint(1, 3000)
            pv = Fraction(kk * kk, 2 ** 104)
            c = float(pv - Fraction(1, 10 ** 30) * rng.choice([1, -1]))
            for _ in range(rng.randint(0, 6)):
                c = math.nextafter(c, math.inf if rng.random() < 0.5 else -math.inf)
            g = 1.0 + kk * 2.0 ** -52
            s = [g, 0.0, 1.0, 0.0, c, 0.0, g, 0.0]
        elif mode == 3:                                     # lossless: Z = jX, S from it in floats
            x11, x12, x22 = (rng.uniform(-200, 200) for _ in range(3))
            zz = [[complex(0, x11), complex(0, x12)], [complex(0, x12), complex(0, x22)]]
            a = [[zz[r][q] - (z0 if r == q else 0) for q in range(2)] for r in range(2)]
            bm = [[zz[r][q] + (z0 if r == q else 0) for q in range(2)] for r in range(2)]
            dt = bm[0][0] * bm[1][1] - bm[0][1] * bm[1][0]
            inv = [[bm[1][1] / dt, -bm[0][1] / dt], [-bm[1][0] / dt, bm[0][0] / dt]]
            sm = [[sum(a[r][m] * inv[m][q] for m in range(2)) for q in range(2)] for r in range(2)]
            s = [sm[0][0].real, sm[0][0].imag, sm[0][1].real, sm[0][1].imag, sm[1][0].real, sm[1][0].imag,
                 sm[1][1].real, sm[1][1].imag]
        elif mode == 4:                                     # exact zeros
            for j in rng.sample(range(8), rng.randint(1, 6)):
                s[j] = 0.0
        elif mode == 5:                                     # the box edges
            for j in rng.sample(range(9), rng.randint(1, 3)):
                x = rng.choice([1.0, -1.0, 1.5]) * 2.0 ** rng.choice([120, -120, 121, -121, 119, -119])
                if j == 8:
                    z0 = abs(x)
                else:
                    s[j] = x
        elif mode == 6:                                     # any finite values
            s = [rng.choice([1.0, -1.0]) * 10.0 ** rng.uniform(-300, 300) for _ in range(8)]
            z0 = 10.0 ** rng.uniform(-300, 300)
        elif mode == 7:                                     # subnormal parts
            for j in rng.sample(range(8), rng.randint(1, 4)):
                s[j] = rng.choice([1.0, -1.0]) * 2.0 ** rng.randint(-1074, -1000)
        v = s + [z0]
        assert len(v) == 9
        if all(math.isfinite(x) for x in v):
            out.append(v)
    return out


def s2z_families(rng, count):
    """s2z (the fast path, then the exact one) and s2z_slow (the exact path alone): each case's nine parts."""
    out = []
    for v in s2z_cases(rng, count):
        ref = s2z_ref(v)
        for part in range(9):
            out.append(("s2z", v + [part], ref[part]))
            out.append(("s2z_slow", v + [part], ref[part]))
    return out


# ── M5's rationals of signed monomials (src/crmath_edge.cyr, M5 part 2a-1a, ADR 0002 row 40) ─────────────────────────
RAT_EPS0 = 0.0000000000088541878128
RAT_C = 299792458.0
RAT_KINDS = 14
RAT_K = list(range(RAT_KINDS)) + [17]   # 14–16 (the step time, two thresholds) are checked by their own rows


def rat_ref(k, v):
    """RN of rational k over v (src/crmath_edge.cyr's _bj_cr_rat_spec), from the coefficient's formula: the media's C_a,
    C_b; Drude α, β, κ; Lorentz α, β, γ_c; Debye α, β; the TF/SF Mur coefficient and the 1D one's numerator (M5 part
    2a-1b); the parallel and perpendicular subpixel means (M5 part 2a-2); the CPML's 1/κ − 1 (M6 part 4b). All ones when
    the denominator is exactly 0. An exact zero takes IEEE's signs: a one-product numerator's factors' (and
    coefficient's), a cancelled sum +0, then D's."""
    F = [Fraction(x) for x in v]
    sg = [math.copysign(1.0, x) < 0 for x in v]
    if k in (0, 1):
        a, b, c, d, e = F[:5]
        num = 2 * e * a - b * c if k == 0 else 2 * c
        den = 2 * e * a + b * c if k == 0 else 2 * e * a * d + b * c * d
        zs = None if k == 0 else sg[2]
    elif k in (2, 3, 4):
        a, b, c, d, e = F[:5]
        num, den, zs = [(2 - c * d, 2 + c * d, None), (e * b * b * d, 2 + c * d, sg[4] ^ sg[3]),
                        (d, e * a, sg[3])][k - 2]
    elif k in (5, 6, 7):
        _, b, c, d, e, f = F[:6]
        num, den, zs = [(4 - 2 * c * c * e * e, 2 + d * e, None), (d * e - 2, 2 + d * e, None),
                        (2 * f * b * c * c * e * e, 2 + d * e, sg[5] ^ sg[1])][k - 5]
    elif k in (8, 9):
        a, b, c, d, e = F[:5]
        num = 2 * c - d if k == 8 else e * b * d - e * a * d
        den, zs = 2 * c + d, None
    elif k == 10:
        dt, dx, c = F[:3]
        num, den, zs = c * dt - dx, c * dt + dx, None
    elif k == 17:                                   # −(κ_max − 1)·g/(1 + (κ_max − 1)·g), g = b·c·d
        a, g = F[0], F[1] * F[2] * F[3]
        num, den, zs = -(a - 1) * g, 1 + (a - 1) * g, None
    elif k in (12, 13):
        e1, e2, fr = F[:3]
        num, den = ((1 - fr) * e1 + fr * e2, Fraction(1)) if k == 12 else (e1 * e2, (1 - fr) * e2 + fr * e1)
        zs = None if k == 12 else sg[0] ^ sg[1]
    else:
        c, dt, e, m, dx = F[:5]
        num, den, zs = c * c * dt * dt - e * m * dx * dx, Fraction(1), None
    if den == 0:
        return (1 << 64) - 1
    if num == 0:
        neg = (zs or False) ^ (den < 0)
        return f2b(-0.0 if neg else 0.0)
    return f2b(rn_frac(num / den))


def rat_cases(rng, count):
    """(k, v): each rational on typical media and FDTD steps; at its cancellation (s → 1, γΔt → 2, ω₀Δt → √2,
    2τ → Δt, ε_s → ε_∞) by relative offsets and ulp by ulp, with exact zeros; the fast path's box edges at 2^±120;
    signed zeros; any finite values and subnormals."""
    out = []
    for i in range(count):
        k = RAT_K[i % len(RAT_K)]
        mode = (i // len(RAT_K)) % 7
        dt = 10.0 ** rng.uniform(-18, -10)
        e0 = RAT_EPS0
        dxx = 10.0 ** rng.uniform(-6, -1)
        if k in (12, 13):                           # the subpixel means: equal, adjacent or any media, any fill
            e1 = 10.0 ** rng.uniform(-1, 2)
            e2 = rng.choice([e1, math.nextafter(e1, math.inf), 10.0 ** rng.uniform(-1, 2), e1 * 1.000001])
            fr = rng.choice([rng.random(), rng.randint(0, 10) / 10, rng.randint(0, 3) / 3, 0.0, 1.0,
                             math.nextafter(1.0, 0.0), 5e-324, rng.random() * 1e-20])
            v = [e1, e2, fr, 0.0, 0.0, 0.0]
            if mode >= 3:                           # wide and edge media, the fill kept in [0, 1]
                v = [10.0 ** rng.uniform(-300, 300), 10.0 ** rng.uniform(-300, 300), rng.random(), 0.0, 0.0, 0.0]
                if mode == 6:
                    v[rng.randint(0, 1)] = rng.choice([5e-324, 2.2250738585072014e-308, 1.7976931348623157e308])
            out.append((k, v))
            continue
        if k == 17:                                 # the CPML's 1/κ − 1: κ_max near 1, CFS gradings, any grade
            km = rng.choice([1.0, math.nextafter(1.0, 2.0), 1 + 10.0 ** -rng.uniform(1, 15), rng.uniform(1, 20),
                             10.0 ** rng.uniform(0, 300), 1.7976931348623157e308])
            g = rng.choice([rng.random(), 10.0 ** -rng.uniform(0, 30), 1.0, 2.0 ** -rng.randint(1, 1022)])
            v = [km, g, 1.0, 1.0, 0.0, 0.0]
            if mode in (3, 6):                      # an edge grade s·2^e, e ≤ −1023, its power split as the port's
                e = -rng.randint(1023, 2100)
                p = max(e, -1022)
                q = max(e - p, -1022)
                v = [km, (1 + rng.random()) * 2.0 ** p, 2.0 ** q, 2.0 ** (e - p - q), 0.0, 0.0]
            elif mode == 4:                         # the cancelled sum: κ_max = 1, or a zero grade
                v[rng.randint(0, 1)] = rng.choice([1.0, 0.0, -0.0]) if v[0] != 1.0 else 0.0
            elif mode == 5:                         # any finite values, either sign
                v = [rng.choice([1.0, -1.0]) * 10.0 ** rng.uniform(-300, 300) for _ in range(4)] + [0.0, 0.0]
            out.append((k, v))
            continue
        if k == 10:                                 # the TF/SF Mur coefficient near S = 1
            v = [dxx / RAT_C * (1 - 10.0 ** -rng.uniform(0, 16)), dxx, RAT_C, 0.0, 0.0, 0.0]
        elif k == 11:                               # the 1D one's numerator near ε_rμ_r = 1/4 at dt = dx/2c
            q = 0.25 * (1 + rng.choice([1.0, -1.0]) * 10.0 ** -rng.uniform(0, 16))
            er = rng.uniform(0.3, 3)
            v = [RAT_C, dxx / (2 * RAT_C), er, q / er, dxx, 0.0]
        elif k in (0, 1):
            v = [rng.uniform(1, 20), 10.0 ** rng.uniform(-6, 8), dt, 10.0 ** rng.uniform(-6, -2), e0, 0.0]
        elif k in (2, 3, 4):
            v = [rng.uniform(1, 10), 10.0 ** rng.uniform(14, 17), 10.0 ** rng.uniform(11, 15), dt, e0, 0.0]
        elif k in (5, 6, 7):
            v = [0.0, rng.uniform(0, 10), 10.0 ** rng.uniform(13, 16), 10.0 ** rng.uniform(11, 15), dt, e0]
        else:
            ei = rng.uniform(1, 10)
            v = [ei, rng.uniform(ei, 100), 10.0 ** rng.uniform(-13, -9), dt, e0, 0.0]
        if k >= 10 and mode in (1, 2):              # already at their cancellation; walk an operand ulp by ulp
            j = 0 if k == 10 else 3
            for _ in range(rng.randint(0, 6)):
                v[j] = math.nextafter(v[j], math.inf if rng.random() < 0.5 else 0.0)
        elif mode in (1, 2):                        # the cancellation, by a relative offset or ulp by ulp
            rho = rng.choice([1.0, -1.0]) * 10.0 ** -rng.uniform(1, 16)
            if k in (0, 1):
                x = 2 * e0 * v[0] / v[2]
                j = 1
            elif k in (2, 3, 4):
                x, j = 2 / v[3], 2
            elif k in (5, 7):
                x, j = math.sqrt(2) / v[4], 2
            elif k == 6:
                x, j = 2 / v[4], 3
            elif k == 8:
                x, j = v[3] / 2, 2
            else:
                x, j = v[0], 1
            if mode == 1:
                x *= 1 + rho
            else:
                for _ in range(rng.randint(0, 6)):
                    x = math.nextafter(x, math.inf if rng.random() < 0.5 else 0.0)
            v[j] = x
            if rng.random() < 0.2 and k in (2, 6, 8):   # an exact zero: γΔt = 2, 2τ = Δt
                dt2 = 2.0 ** -rng.randint(30, 60)
                if k == 8:
                    v[3], v[2] = dt2, dt2 / 2
                else:
                    v[3 if k == 2 else 4] = dt2
                    v[2 if k == 2 else 3] = 2 / dt2
        elif mode == 3:                             # the box edges
            for j in rng.sample(range(6), rng.randint(1, 3)):
                v[j] = rng.choice([1.0, -1.0, 1.5]) * 2.0 ** rng.choice([120, -120, 121, -121, 119, -119])
        elif mode == 4:                             # zeros, signed
            for j in rng.sample(range(6), rng.randint(1, 2)):
                v[j] = rng.choice([0.0, -0.0])
        elif mode == 5:                             # any finite values, either sign
            v = [rng.choice([1.0, -1.0]) * 10.0 ** rng.uniform(-300, 300) for _ in range(6)]
        elif mode == 6:                             # subnormal and huge parts
            for j in rng.sample(range(6), rng.randint(1, 3)):
                v[j] = rng.choice([1.0, -1.0]) * (2.0 ** rng.randint(-1074, -1000) if rng.random() < 0.6
                                                  else 2.0 ** rng.randint(1000, 1023))
        out.append((k, v))
    return out


def rat_families(rng, count):
    """rat (the fast path's two tiers, then the exact path), rat_x (its second tier, then the exact path) and rat_slow
    (the exact path alone)."""
    out = []
    for k, v in rat_cases(rng, count):
        ref = rat_ref(k, v)
        out.append(("rat", [k] + v, ref))
        out.append(("rat_x", [k] + v, ref))
        out.append(("rat_slow", [k] + v, ref))
    return out


def nb_families(rng, count):
    """lw_nb and lw_nb_slow (M2 part 2d-1): each component of n̂ − β on lw_cases' pairs, and on pairs built to cancel
    (β along n̂ to 2^-k, a Pythagorean R with βᵢ its exact n̂ᵢ's double, zero parts)."""
    out = []
    cases = lw_cases(rng, count)
    for j in range(max(1, count // 4)):
        c = j % 4
        if c == 3:
            # the widest spread on doubles (M2 part 2d-3i's regression): R'ᵢ near 2^-2096 beside R near DBL_MAX, βᵢ
            # subnormal or tiny of Rᵢ's sign, where the exact path's candidate overflowed on R'ᵢ's scale
            sg = rng.choice([1.0, -1.0])
            r = [rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(1015, 1023) * rng.choice([1.0, -1.0]),
                 sg * b2f(rng.getrandbits(rng.randint(1, 20)) | 1), rng.choice([0.0, rng.uniform(-1.0, 1.0)])]
            b = [rng.uniform(-0.9, 0.9), sg * b2f(rng.getrandbits(rng.randint(1, 52)) | 1), rng.uniform(-0.3, 0.3)]
            if rng.random() < 0.3:
                b[1] = sg * rng.uniform(1.0, 2.0) * 2.0 ** -rng.randint(1000, 1022)
        elif c == 0:
            n = lw_unit(rng)
            sc = 2.0 ** rng.randint(-60, 60)
            sp = 1.0 - 2.0 ** -rng.uniform(1.0, 52.0)
            r = [x * sc for x in n]
            b = [cm_ulps(abs(x * sp), rng.randint(-3, 3)) * math.copysign(1.0, x) if x != 0.0 else 0.0 for x in n]
        elif c == 1:
            a, bq, cq, d = rng.choice([(3, 4, 0, 5), (1, 2, 2, 3), (2, 3, 6, 7), (1, 4, 8, 9)])
            sc = 2.0 ** rng.randint(-30, 30)
            r = [a * sc, bq * sc, cq * sc]
            b = [a / d * 0.9, bq / d * 0.9, cq / d * 0.9] if rng.random() < 0.5 else [a / d, bq / d * 0.5, 0.0]
        else:
            r = [rng.uniform(-3.0, 3.0), 0.0, rng.choice([0.0, -0.0, rng.uniform(-3.0, 3.0)])]
            b = [rng.uniform(-0.5, 0.5), rng.choice([0.0, -0.0, 0.3]), rng.uniform(-0.5, 0.5)]
        if lw_valid(r, b):
            cases.append((r, b))
    for r, b in cases:
        below = lw_below(r, b)
        for i in range(3):
            if below:
                want = LW_BELOW
            elif r[i] == 0.0:
                want = f2b(r[i] - b[i])
            else:
                want = f2b(rn_nb(r, b, i))
            out.append(("lw_nb", [i] + r + b, want))
            out.append(("lw_nb_slow", [i] + r + b, want))
    return out


# ── the offset forms of _bj_cr_lw and _bj_cr_lorentz (M2 part 2d-3i, for lienard_wiechert_fields' chain) ──────────
OFF_BIG = 2.0 - 2.0 ** -52                  # 2·DBL_MAX = OFF_BIG·2^1024


def off_big(rng):
    """A part past DBL_MAX in the offset form, |R| in (DBL_MAX, 2·DBL_MAX]: (significand in [1, 2), 1024)."""
    m = rng.choice([1.0, OFF_BIG, rng.uniform(1.0, OFF_BIG), 1.0 + 2.0 ** -52])
    return m * rng.choice([1.0, -1.0]), 1024


def off_tiny(rng):
    """A β part on the 53-bit grid below 2^-1022 (a quotient v/c), in the offset form: (significand, exponent) with
    the value in [2^-1110, 2^-1022)."""
    m = rng.choice([1.0, rng.uniform(1.0, 2.0), 2.0 - 2.0 ** -52, b2f(f2b(1.0) | rng.getrandbits(52) | 1)])
    return m * rng.choice([1.0, -1.0]), rng.randint(-1110, -1023)


def off_val(x, o):
    return Fraction(x) * Fraction(2) ** o


def lwo_cases(rng, count):
    """(r, b, offsets) in _bj_cr_lw's offset form, with R and β exact: R with one to three parts past DBL_MAX beside
    any finite parts (odd subnormals among them); lw_cases' pairs scaled by 2^k so their largest part passes DBL_MAX
    (the knife edges and ties carried along); ordinary pairs written with offsets; β with parts below 2^-1022 on the
    53-bit grid, near the cone's edge too; the widest exponent spreads (R past DBL_MAX beside 2^-1074, β at 2^-1110)."""
    out = []
    base = lw_cases(rng, count)
    i = 0
    while len(out) < count:
        kind = i % 6
        i += 1
        r, b, o = [0.0] * 3, [0.0] * 3, [0] * 6
        if kind == 0:
            r = [any_finite(rng) for _ in range(3)]
            j = rng.randrange(3)
            r[j], o[j] = off_big(rng)
            u = lw_unit(rng)
            sp = rng.uniform(0.0, 0.999)
            b = [x * sp for x in u]
        elif kind == 1:
            u = lw_unit(rng)
            for j in range(3):
                if rng.random() < 0.6:
                    r[j], o[j] = abs(u[j]) * 1.5 * rng.choice([1.0, 1.0 + 2.0 ** -52]) or 1.0, 1024
                    r[j] = math.copysign(min(r[j], OFF_BIG), u[j])
                else:
                    r[j] = u[j] * 2.0 ** 1000
            if all(x == 0 for x in o[:3]):
                r[0], o[0] = off_big(rng)
            sp = 1.0 - 2.0 ** -rng.uniform(1.0, 52.0)
            b = [x * sp for x in u]
        elif kind == 2:
            rv, bv = base[i % len(base)]
            mx = max(abs(x) for x in rv)
            if mx == 0:
                continue
            sh = 1024 - binade(Fraction(mx))       # the largest part in [2^1024, 2^1025)
            if max(abs(off_val(x, sh)) for x in rv) > Fraction(sys.float_info.max) * 2:
                sh -= 1
            r, b, o = list(rv), list(bv), [sh, sh, sh, 0, 0, 0]
        elif kind == 3:
            rv, bv = base[i % len(base)]
            sh = rng.randint(-40, 40)
            try:
                r = [math.ldexp(x, -sh) for x in rv]
            except OverflowError:
                continue
            if any(Fraction(y) * Fraction(2) ** sh != Fraction(x) for x, y in zip(rv, r)):
                continue
            b, o = list(bv), [sh, sh, sh, 0, 0, 0]
        elif kind == 4:
            r = [any_finite(rng), rng.choice([0.0, -0.0, b2f(rng.getrandbits(52) | 1)]), any_finite(rng)]
            r[0], o[0] = off_big(rng)
            b = [rng.uniform(-0.99, 0.99), 0.0, rng.choice([0.0, rng.uniform(-0.1, 0.1)])]
        else:
            r = [0.0, b2f(1), 0.0]
            r[0], o[0] = off_big(rng)
            r[2] = rng.choice([0.0, b2f(rng.getrandbits(52) | 1), -1e-300])
            b = [math.copysign(1.0 - 2.0 ** -rng.randint(1, 53), r[0]), 0.0, 0.0]
        for j in range(3):
            if kind in (0, 1, 4, 5) and ((b[j] == 0.0 and rng.random() < 0.8) or rng.random() < 0.15):
                b[j], o[3 + j] = off_tiny(rng)
        rt = [off_val(x, e) for x, e in zip(r, o[:3])]
        bt = [off_val(x, e) for x, e in zip(b, o[3:])]
        if not all(math.isfinite(x) for x in r + b):
            continue
        if any(abs(x) > Fraction(sys.float_info.max) * 2 for x in rt):
            continue
        if sum(x * x for x in rt) < Fraction(1, 10 ** 60) or sum(x * x for x in bt) >= 1:
            continue
        out.append((r, b, o, rt, bt))
    return out


def lwo_families(rng, count):
    """lwo_reff, lwo_kappa and lwo_nb (each with its *_slow) and lwo_ux on lwo_cases, against exact.py's rn_lw, rn_nb
    and rnd53_lw on the exact R and β: n̂ − β of a zero Rᵢ is RN(±0 − βᵢ), signed zeros included."""
    out = []
    M = (1 << 64) - 1
    for r, b, o, rt, bt in lwo_cases(rng, count):
        below = lw_below(rt, bt)
        reff = f2b(rn_lw(rt, bt, 0))
        kap = LW_BELOW if below else f2b(rn_lw(rt, bt, 1))
        ws = r + b + o
        out.append(("lwo_reff", [0] + ws, reff))
        out.append(("lwo_reff_slow", [0] + ws, reff))
        out.append(("lwo_kappa", [0] + ws, kap))
        out.append(("lwo_kappa_slow", [0] + ws, kap))
        for i in range(3):
            if below:
                want = LW_BELOW
            elif rt[i] == 0:
                want = f2b(-rn_frac(bt[i])) if bt[i] != 0 else f2b(r[i] - b[i])
            else:
                want = f2b(rn_nb(rt, bt, i))
            out.append(("lwo_nb", [i] + ws, want))
            out.append(("lwo_nb_slow", [i] + ws, want))
        if not below:
            a = rnd53_lw(rt, bt)
            e = binade(a)
            out.append(("lwo_ux", [0] + ws, f2b(float(a / Fraction(2) ** e))))
            out.append(("lwo_ux", [1] + ws, e & M))
    return out


def gammao_families(rng, count):
    """gammao and gammao2 (γ and γ², each with its *_slow) on _bj_cr_lorentz's offset form, c = 1 (a β) or the speed
    of light: one or two normal parts (near c, nearly cancelling c², random) beside parts below 2^-1022 on the 53-bit
    grid, or those parts alone (γ = 1), against gamma_ref on the exact values."""
    out = []
    for i in range(count):
        kind = i % 4
        c = 1.0 if i % 2 == 0 else C_LIGHT
        v, o = [0.0, 0.0, 0.0], [0, 0, 0]
        if kind == 0:
            v[0] = math.copysign(c * (1.0 - rng.uniform(0.5, 1.0) * 2.0 ** -rng.randint(1, 52)), rng.uniform(-1, 1))
        elif kind == 1:
            vx = rng.uniform(0.3, 0.9) * c
            v[0], v[1] = vx, math.sqrt(c * c - vx * vx) * (1.0 - 2.0 ** -rng.randint(30, 53))
        elif kind == 2:
            v[0], v[1] = rng.uniform(-c, c) * 0.6, rng.uniform(-c, c) * 0.6
        for j in range(3):
            if v[j] == 0.0 and (kind == 3 or rng.random() < 0.7):
                v[j], o[j] = off_tiny(rng)
        vt = [off_val(x, e) for x, e in zip(v, o)]
        for name, mode in (("gammao", 0), ("gammao_slow", 0), ("gammao2", 1), ("gammao2_slow", 1)):
            out.append((name, [3, c] + v + o, gamma_ref(vt, c, mode)))
    return out


VX_MN = 2.0 ** -1022
DBL_MAX = sys.float_info.max


def vx_pos(rng):
    """A positive finite double: subnormal, a special, at either end of the exponent range, or anywhere."""
    c = rng.random()
    if c < 0.2:
        return b2f(rng.randint(1, (1 << 52) - 1))
    if c < 0.3:
        return rng.choice([2.0 ** -1074, 2.0 ** -1073, VX_MN, DBL_MAX, DBL_MAX / 2.0, 1.0, 2.0 ** 1023])
    if c < 0.6:
        return b2f((rng.choice([rng.randint(1, 40), rng.randint(2006, 2046)]) << 52) | rng.getrandbits(52))
    return b2f((rng.randint(1, 2046) << 52) | rng.getrandbits(52))


def vx_near(x, k):
    """x moved k ulps, kept positive and finite."""
    b = min(max(f2b(x) + k, 1), f2b(DBL_MAX))
    return b2f(b)


def vx_families(rng, count):
    """vx_sqrt, vx_sqrt_e and vx_gamma (M2 part 2d-2c): √(x/y) and √(x·y) to 53 bits with an unbounded exponent (the
    exact root by isqrt, exact.py's rnd53_sqrt), and Γ against the unbounded chain alone: η = √(μ_r/ε_r) to 53 bits,
    their difference and sum to 53 bits, the quotient rounded once. That Γ is the expected value everywhere, so the rows
    where _bj_vx_gamma takes Rust's quotient check its guard's claim too. Pairs of materials independent, sharing ε_r
    with μ_r a few ulps apart (η's within an ulp: cancellation), and scaled by a power of 4 (equal η's)."""
    F = Fraction
    out = []
    for _ in range(count):
        x, y = vx_pos(rng), vx_pos(rng)
        for mode in (0, 1):
            a = rnd53_sqrt(F(x) / F(y)) if mode == 0 else rnd53_sqrt(F(x) * F(y))
            e = binade(a)
            out.append(("vx_sqrt", [x, y, mode], f2b(float(a / F(2) ** e))))
            out.append(("vx_sqrt_e", [x, y, mode], e & ((1 << 64) - 1)))
    for j in range(count):
        e1, m1 = vx_pos(rng), vx_pos(rng)
        c = j % 3
        if c == 0:
            e2, m2 = vx_pos(rng), vx_pos(rng)
        elif c == 1:
            e2, m2 = e1, vx_near(m1, rng.randint(-3, 3))
        else:
            k = rng.randint(-8, 8)
            e2, m2 = e1 * 4.0 ** k, m1 * 4.0 ** k
            if not (0.0 < e2 <= DBL_MAX and 0.0 < m2 <= DBL_MAX):
                e2, m2 = e1, m1
        h1, h2 = rnd53_sqrt(F(m1) / F(e1)), rnd53_sqrt(F(m2) / F(e2))
        num = rnd53(h2 - h1)
        want = 0.0 if num == 0 else rn_frac(num / rnd53(h2 + h1))
        eta1, eta2 = rn_sqrt(F(m1) / F(e1)), rn_sqrt(F(m2) / F(e2))
        out.append(("vx_gamma", [eta1, eta2, e1, m1, e2, m2], f2b(want)))
    return out


def dot4_ref(xs, ys) -> int:
    """The bits of Σ xs[k]·ys[k] from −0.0 with an unbounded exponent (M3 part 2d-1a): each product and partial sum
    rounded to 53 bits, the last sum once into binary64, IEEE's signed zeros throughout (finite inputs)."""
    F = Fraction

    def neg(x):
        return math.copysign(1.0, x) < 0

    def add(a, b, r=rnd53):
        v = a[0] + b[0]
        if v == 0:
            return F(0), (a[1] and b[1]) if a[0] == 0 and b[0] == 0 else False
        return r(v), v < 0

    acc = (F(0), True)
    for k in range(4):
        v = F(xs[k]) * F(ys[k])
        acc = add(acc, (rnd53(v), (neg(xs[k]) != neg(ys[k])) if v == 0 else v < 0), (lambda t: t) if k == 3 else rnd53)
    if acc[0] == 0:
        return f2b(-0.0 if acc[1] else 0.0)
    return f2b(rn_frac(acc[0]))


def dot4_cases(rng, count):
    """Four products around the range edges: anywhere in the range (zeros of both signs among them); near 2^-1022,
    some below it and some normal; past DBL_MAX, cancelling or not; cancelling to small sums that subnormal products
    then meet; and products just below 2^-1022 that round up to it."""
    def sig():
        return 1.0 + rng.getrandbits(52) * 2.0 ** -52

    def pm(x):
        return -x if rng.random() < 0.5 else x

    def at(e):                                  # a random double in [2^e, 2^(e + 1)), e within the normal range
        return pm(sig() * 2.0 ** e)

    def split(e):                               # a factor pair whose product lies near 2^e
        ex = rng.randint(max(-1022, e - 1023), min(1023, e + 1022))
        return at(ex), at(e - ex)

    out = []
    for j in range(count):
        c = j % 5
        xs, ys = [], []
        for k in range(4):
            if c == 0:
                if rng.random() < 0.1:
                    x, y = pm(0.0), at(rng.randint(-1022, 1023))
                else:
                    x, y = pm(vx_pos(rng)), pm(vx_pos(rng))
            elif c == 1:
                x, y = split(rng.randint(-1110, -1000) if rng.random() < 0.7 else rng.randint(-1000, -780))
            elif c == 2:
                x, y = split(rng.randint(990, 1030))
                if k > 0 and rng.random() < 0.6:
                    x, y = xs[0], -ys[0] * (1.0 + rng.randint(-4, 4) * 2.0 ** -52)
            elif c == 3:
                if k == 1:
                    x, y = xs[0], -ys[0] * (1.0 + rng.randint(-3, 3) * 2.0 ** -52)
                else:
                    x, y = split(rng.randint(-1100, -900) if k > 1 else rng.randint(-1000, -850))
            else:
                if k % 2 == 0:
                    x, y = pm(1.0 - rng.randint(1, 3) * 2.0 ** -53), VX_MN
                else:
                    x, y = pm(1.0), pm(VX_MN - rng.randint(0, 2) * 2.0 ** -1074)
            if not (math.isfinite(x) and math.isfinite(y)):
                x, y = 1.0, 1.0
            xs.append(x)
            ys.append(y)
        out.append((xs, ys))
    return out


def dot4_families(rng, count):
    """dot4 (M3 part 2d-1a): the Mueller sum and its guard (cr_driver.cyr's _drv_dot4) against the unbounded chain
    alone, with b contiguous (mode 0) and at a stride of 4 (mode 1): the cases the guard keeps Rust's sum on check its
    claim."""
    out = []
    for xs, ys in dot4_cases(rng, count):
        ref = dot4_ref(xs, ys)
        for mode in (0, 1):
            out.append(("dot4", xs + ys + [mode], ref))
    return out


# ── the Jones chains (src/crmath_vx.cyr, M3 part 2d-1b): exact references ─────────────────────────────────────────
# Each value is a pair (Fraction, negative zero), every operation rounded to 53 bits with an unbounded exponent (IEEE's
# signed zeros throughout), the last one once into binary64: the unbounded chain alone, so the calls the kernels'
# boxes keep on Rust's chain check that claim.

def _jx(x):
    return Fraction(x), math.copysign(1.0, x) < 0


def _jmul(a, b, r=rnd53):
    v = a[0] * b[0]
    return r(v), (a[1] != b[1]) if v == 0 else v < 0


def _jadd(a, b, r=rnd53):
    v = a[0] + b[0]
    if v == 0:
        return Fraction(0), (a[1] and b[1]) if a[0] == 0 and b[0] == 0 else False
    return r(v), v < 0


def _jsub(a, b, r=rnd53):
    return _jadd(a, (-b[0], not b[1]), r)


def _jfin(a) -> float:
    if a[0] == 0:
        return -0.0 if a[1] else 0.0
    return rn_frac(a[0])


def _jcprod(a, x):
    """The parts of the complex product a·x, each operation to 53 bits."""
    re = _jsub(_jmul(_jx(a[0]), _jx(x[0])), _jmul(_jx(a[1]), _jx(x[1])))
    im = _jadd(_jmul(_jx(a[0]), _jx(x[1])), _jmul(_jx(a[1]), _jx(x[0])))
    return re, im


def _jinten(p):
    def sq2(u, v):
        return _jadd(_jmul(_jx(u), _jx(u)), _jmul(_jx(v), _jx(v)))
    return _jadd(sq2(p[0], p[1]), sq2(p[2], p[3]))


def jcdot_ref(a, x, b, y, mode):
    """A·X ± B·Y (mode & 2: −; mode & 1: A and B conjugated), each part's last sum rounded once."""
    if mode & 1:
        a, b = (a[0], -a[1]), (b[0], -b[1])
    p, q = _jcprod(a, x), _jcprod(b, y)
    op = _jsub if mode & 2 else _jadd
    return [_jfin(op(p[0], q[0], lambda t: t)), _jfin(op(p[1], q[1], lambda t: t))]


def jint_ref(p):
    def sq2(u, v):
        return _jadd(_jmul(_jx(u), _jx(u)), _jmul(_jx(v), _jx(v)))
    return _jfin(_jadd(sq2(p[0], p[1]), sq2(p[2], p[3]), lambda t: t))


def jnorm_ref(p):
    k = rnd53(1 / rnd53_sqrt(_jinten(p)[0]))
    return [_jfin(_jmul(_jx(t), (k, False), lambda u: u)) for t in p]


def jover_ref(p, q):
    den = _jmul(_jinten(p), _jinten(q))
    u, v = _jcprod((p[0], -p[1]), (q[0], q[1])), _jcprod((p[2], -p[3]), (q[2], q[3]))
    wr, wi = _jadd(u[0], v[0]), _jadd(u[1], v[1])
    n = _jadd(_jmul(wr, wr), _jmul(wi, wi))
    return rn_frac(n[0] / den[0])


def jpart(rng):
    """A part around the kernels' box: ±0, inside [2^-511, 2^511), at its edges, tiny or subnormal, huge."""
    c = rng.random()
    sg = -1.0 if rng.random() < 0.5 else 1.0
    if c < 0.12:
        return sg * 0.0
    if c < 0.2:
        return sg * b2f(rng.randint(1, (1 << 52) - 1))
    if c < 0.45:
        e = rng.randint(-511, 510)
    elif c < 0.65:
        e = rng.choice([rng.randint(-530, -500), rng.randint(500, 520)])
    elif c < 0.8:
        e = rng.randint(-1022, -530)
    elif c < 0.9:
        e = rng.randint(520, 1023)
    else:
        e = rng.randint(-30, 30)
    return sg * (1.0 + rng.getrandbits(52) * 2.0 ** -52) * 2.0 ** e


def jbox(rng):
    """A part inside the box: ±0, or 2^-511 ≤ |x| < 2^511, often of moderate size."""
    sg = -1.0 if rng.random() < 0.5 else 1.0
    if rng.random() < 0.12:
        return sg * 0.0
    e = rng.randint(-60, 60) if rng.random() < 0.5 else rng.randint(-511, 510)
    return sg * (1.0 + rng.getrandbits(52) * 2.0 ** -52) * 2.0 ** e


def jones_families(rng, count):
    """jcdot, jint, jnorm and jover (M3 part 2d-1b): src/crmath_vx.cyr's Jones kernels against the unbounded chain:
    parts around the box, products cancelling (B·Y near −A·X), intensities around 1e-290 and past DBL_MAX, and nearly
    orthogonal pairs whose inner product lies below 2^-511."""
    out = []
    for j in range(count):
        part = jbox if j % 6 in (1, 4) else jpart
        a, x = (part(rng), part(rng)), (part(rng), part(rng))
        if j % 3 == 0:
            b = a
            k = 1.0 + rng.randint(-4, 4) * 2.0 ** -52
            y = (-x[0] * k, -x[1] * k)
            if not all(math.isfinite(t) for t in y):
                y = (jpart(rng), jpart(rng))
        else:
            b, y = (part(rng), part(rng)), (part(rng), part(rng))
        mode = j % 4
        ref = jcdot_ref(a, x, b, y, mode)
        for part in (0, 1):
            out.append(("jcdot", [a[0], a[1], x[0], x[1], b[0], b[1], y[0], y[1], mode, part], f2b(ref[part])))
    for j in range(count):
        if j % 3 == 0:                          # an intensity near 1e-290, or near DBL_MAX
            base = 1e-145 if j % 2 == 0 else 1.0e154
            pv = [base * rng.uniform(0.5, 1.2), base * rng.uniform(-0.3, 0.3), base * rng.uniform(-0.3, 0.3), 0.0]
        else:
            pv = [jpart(rng) for _ in range(4)]
        if all(t == 0 for t in pv):
            pv[0] = 1.0
        out.append(("jint", pv, f2b(jint_ref(pv))))
        ref = jnorm_ref(pv)
        for k in range(4):
            out.append(("jnorm", pv + [k], f2b(ref[k])))
    for j in range(count):
        part = jbox if j % 2 == 1 else jpart
        pv = [part(rng) for _ in range(4)]
        if j % 3 == 0:                          # q nearly orthogonal to p: (−p_y*, p_x*) scaled, one part nudged
            sc = 2.0 ** (rng.randint(-40, 40) if part is jbox else rng.randint(-400, 400))
            qv = [-pv[2] * sc, pv[3] * sc, pv[0] * sc, -pv[1] * sc]
            i = rng.randrange(4)
            qv[i] = qv[i] * (1.0 + rng.randint(-3, 3) * 2.0 ** -52) + part(rng) * 2.0 ** -60
        elif j % 6 == 1:                        # boxed pairs with a small inner product, I_p·I_q near 1e-290
            pv = [rng.uniform(1.0, 2.0) * 2.0 ** -19, 0.0, rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(20, 60), 0.0]
            qv = [rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-511, -480), rng.uniform(-1.0, 1.0) * 2.0 ** -500, 0.0,
                  0.0]
        else:
            qv = [part(rng) for _ in range(4)]
        if not all(math.isfinite(t) for t in qv):
            qv = [jpart(rng) for _ in range(4)]
        if all(t == 0 for t in pv):
            pv[0] = 1.0
        if all(t == 0 for t in qv):
            qv[2] = 1.0
        out.append(("jover", pv + qv, f2b(jover_ref(pv, qv))))
    return out


# ── the Stokes kernels (src/crmath_vx.cyr, M3 part 2d-1c): exact references ───────────────────────────────────────
# from_jones on exact.py's chain, the degree of polarization and the angles on MPFR at 53 bits with its widest exponent
# range (mpfr_ref.wide), each last operation rounded once: the unbounded chain alone.

def _whalf(t) -> float:
    return 0.5 * t if isinstance(t, float) else rn_frac(t / 2)


def sjones_ref(p):
    def sq2(u, v):
        return _jadd(_jmul(_jx(u), _jx(u)), _jmul(_jx(v), _jx(v)))
    ix, iy = sq2(p[0], p[1]), sq2(p[2], p[3])
    re, im = _jcprod((p[0], p[1]), (p[2], -p[3]))
    return [_jfin(_jadd(ix, iy, lambda t: t)), _jfin(_jsub(ix, iy, lambda t: t)), _jfin((2 * re[0], re[1])),
            _jfin((2 * im[0], im[1]))]


def sdop_ref(s):
    h = mpfr_ref.wide("hypot", mpfr_ref.wide("hypot", s[1], s[2]), s[3])
    return 0.0 if h == 0 else rn_frac(Fraction(h) / Fraction(s[0]))


def sellip_ref(s):
    def sq(u):
        return _jmul(_jx(u), _jx(u))
    q = _jadd(sq(s[1]), sq(s[2]))
    lin = rnd53_sqrt(q[0]) if q[0] != 0 else 0.0
    return _whalf(mpfr_ref.wide("atan2", s[3], lin))


def sorient_ref(s):
    return _whalf(mpfr_ref.wide("atan2", s[2], s[1]))


def spart(rng):
    """A Stokes part: jpart's, or a double near DBL_MAX, or one whose quotient by a neighbour is tiny."""
    c = rng.random()
    if c < 0.15:
        return (-1.0 if rng.random() < 0.5 else 1.0) * DBL_MAX * rng.uniform(0.5, 1.0)
    return jpart(rng)


def stokes_families(rng, count):
    """sjones, sdop, sellip and sorient (M3 part 2d-1c): from_jones (each part), the degree of polarization and the two
    angles against the unbounded chain: parts around the box, subnormal, near DBL_MAX, hypot operands 2^10 to 2^40
    apart, and angles below 2^-1022 (S₂ or S₃ tiny against S₁)."""
    out = []
    for j in range(count):
        part = jbox if j % 3 == 1 else jpart
        pv = [part(rng) for _ in range(4)]
        ref = sjones_ref(pv)
        for k in range(4):
            out.append(("sjones", pv + [k], f2b(ref[k])))
    for j in range(count):
        sv = [abs(spart(rng)) or 1.0, spart(rng), spart(rng), spart(rng)]
        if j % 4 == 1:                          # the second hypot's operands 2^10 to 2^40 apart
            sv[3] = sv[1] * 2.0 ** -rng.randint(10, 40) * rng.uniform(1.0, 2.0)
        elif j % 4 == 2:                        # an angle below 2^-1022
            sv[1] = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-20, 300)
            sv[2] = rng.uniform(-1.0, 1.0) * sv[1] * 2.0 ** -rng.randint(1023, 1100)
            sv[3] = rng.uniform(-1.0, 1.0) * sv[1] * 2.0 ** -rng.randint(1023, 1100)
        if not all(math.isfinite(t) for t in sv) or sv[0] <= 0:
            sv = [1.0, jpart(rng), jpart(rng), jpart(rng)]
        out.append(("sdop", sv, f2b(sdop_ref(sv))))
        out.append(("sellip", sv, f2b(sellip_ref(sv))))
        out.append(("sorient", sv, f2b(sorient_ref(sv))))
    return out


# ── the Rayleigh kernels (src/crmath_vx.cyr, M3 part 2d-2): exact references ───────────────────────────────────────
# σ, Q, the phase function and the angular intensity on exact.py's chain: every operation (Smith's division on its
# branch decided exactly) to 53 bits with an unbounded exponent, the last one rounded once. A zero divisor gives the
# driver's all-ones word.
RAY_DZ = (1 << 64) - 1
RAY_PI = struct.unpack("<d", bytes.fromhex("182d4454fb210940"))[0]


def _rk2(mre, mim):
    """|K|² as a pair (Fraction, negative zero), or None for m² + 2 = 0."""
    re = _jsub(_jmul(_jx(mre), _jx(mre)), _jmul(_jx(mim), _jx(mim)))
    im = _jadd(_jmul(_jx(mre), _jx(mim)), _jmul(_jx(mim), _jx(mre)))
    a = (_jsub(re, _jx(1.0)), _jsub(im, _jx(0.0)))
    b = (_jadd(re, _jx(2.0)), _jadd(im, _jx(0.0)))
    if b[0][0] == 0 and b[1][0] == 0:
        return None

    def div(u, v):
        q = u[0] / v[0]
        return rnd53(q), (u[1] != v[1]) if q == 0 else q < 0
    if abs(b[0][0]) >= abs(b[1][0]):
        r = div(b[1], b[0])
        den = _jadd(b[0], _jmul(b[1], r))
        kr = div(_jadd(a[0], _jmul(a[1], r)), den)
        ki = div(_jsub(a[1], _jmul(a[0], r)), den)
    else:
        r = div(b[0], b[1])
        den = _jadd(_jmul(b[0], r), b[1])
        kr = div(_jadd(_jmul(a[0], r), a[1]), den)
        ki = div(_jsub(_jmul(a[1], r), a[0]), den)
    return _jadd(_jmul(kr, kr), _jmul(ki, ki))


def rsig_ref(r, lam, mre, mim) -> int:
    k = _rk2(mre, mim)
    if k is None:
        return RAY_DZ
    x = _jx(0.0)
    q = _jmul(_jx(2.0 * RAY_PI), _jx(r))[0] / Fraction(lam)
    x = (rnd53(q), False)
    x4 = _jmul(_jmul(_jmul(x, x), x), x)
    c = _jmul(_jmul(_jmul(_jmul(_jx(8.0 / 3.0), _jx(RAY_PI)), _jx(r)), _jx(r)), x4)
    return f2b(_jfin(_jmul(c, k, lambda t: t)))


def rq_ref(x, mre, mim) -> int:
    k = _rk2(mre, mim)
    if k is None:
        return RAY_DZ
    c = _jmul(_jmul(_jmul(_jmul(_jx(8.0 / 3.0), _jx(x)), _jx(x)), _jx(x)), _jx(x))
    return f2b(_jfin(_jmul(c, k, lambda t: t)))


def _rphase(c):
    return _jmul(_jx(0.75), _jadd(_jx(1.0), _jmul(_jx(c), _jx(c))))


def rphase_ref(c) -> int:
    return f2b(_jfin(_jmul(_jx(0.75), _jadd(_jx(1.0), _jmul(_jx(c), _jx(c))), lambda t: t)))


def rang_ref(cs, d, c) -> int:
    d2 = _jmul(_jx(d), _jx(d))
    q = Fraction(cs) / d2[0]
    qq = (rnd53(q), math.copysign(1.0, cs) < 0 if q == 0 else q < 0)
    t = _jmul(qq, _rphase(c))
    v = t[0] / Fraction(4.0 * RAY_PI)
    return f2b(-0.0 if v == 0 and t[1] else (0.0 if v == 0 else rn_frac(v)))


def rpos(rng):
    """A positive double around σ's box (2^-100, 2^100), far past it, or ordinary."""
    c = rng.random()
    if c < 0.3:
        e = rng.choice([rng.randint(-104, -96), rng.randint(96, 104)])
    elif c < 0.6:
        e = rng.randint(-330, 330)
    else:
        e = rng.randint(-40, 10)
    return rng.uniform(1.0, 2.0) * 2.0 ** e


def rindex(rng):
    """An index m around the box (parts ±0 or in [2^-40, 2^40]): ordinary, at the box's edges, huge, tiny, near 1."""
    c = rng.random()
    if c < 0.3:
        return rng.uniform(1.0, 3.0), -rng.choice([0.0, rng.uniform(0.0, 2.0)])
    if c < 0.5:
        e = rng.choice([rng.randint(-44, -36), rng.randint(36, 44)])
        return rng.uniform(1.0, 2.0) * 2.0 ** e, -rng.choice([0.0, rng.uniform(0.0, 1.0) * 2.0 ** e])
    if c < 0.7:
        e = rng.randint(100, 520)
        return rng.uniform(1.0, 2.0) * 2.0 ** e, -rng.uniform(0.0, 2.0) * 2.0 ** (e - rng.randint(0, 60))
    if c < 0.85:
        return rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-520, -100), -rng.choice([0.0, 2.0 ** rng.randint(-600, -100)])
    return 1.0 + rng.randint(-8, 8) * 2.0 ** -52, -rng.choice([0.0, 2.0 ** rng.randint(-200, -20)])


def rayleigh_families(rng, count):
    """rsig, rq, rphase and rang (M3 part 2d-2): σ, Q, the phase function and the angular intensity against the
    unbounded chain alone, so the calls the boxes keep on Rust's chain check that claim: radii, wavelengths and size
    parameters around the boxes and far past them, indices around theirs, cosines around 2^511, and distances and
    cross-sections whose σ/d² and its product by P leave the range."""
    out = []
    for _ in range(count):
        r, lam = rpos(rng), rpos(rng)
        mre, mim = rindex(rng)
        out.append(("rsig", [r, lam, mre, mim], rsig_ref(r, lam, mre, mim)))
        x = rpos(rng) * rng.choice([1.0, 2.0 ** -100, 2.0 ** 100])
        out.append(("rq", [x, mre, mim], rq_ref(x, mre, mim)))
        c = rng.choice([rng.uniform(-1.0, 1.0), rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(505, 515),
                        -rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(505, 515), 2.0 ** rng.randint(-600, 600)])
        if not math.isfinite(c):
            c = 0.5
        out.append(("rphase", [c], rphase_ref(c)))
        cs = rng.choice([rng.uniform(0.0, 1.0) * 2.0 ** rng.randint(-1074, 1023), 0.0,
                         -rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1000, 0)])
        d = rng.choice([1.0, 1e-30, rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-99, 515)]) * rng.choice([1.0, -1.0])
        if not (math.isfinite(cs) and math.isfinite(d)):
            cs, d = 1e-20, 1.0
        out.append(("rang", [cs, d, c], rang_ref(cs, d, c)))
    return out


# ── the beam's chains at their range edges (src/crmath_vx.cyr's _bj_vx_bm_*, M3 part 2d-3a) ─────────────────────────
BM_PI = RAY_PI
BM_MIN = Fraction(2) ** -1022
BM_MAX = Fraction(sys.float_info.max)
BM_FAMS = ("bspot", "bpeak", "broc", "bgouy", "bci", "bdot2", "bw0", "bgp", "bfsr", "bhg", "blg", "bhgf")


def _jdiv(a, b, r=rnd53):
    v = a[0] / b[0]
    return r(v), (a[1] != b[1]) if v == 0 else v < 0


def _bone_sq(a, b):
    """1 + (a/b)², each operation to 53 bits."""
    q = _jdiv(_jx(a), _jx(b))
    return _jadd(_jx(1.0), _jmul(q, q))


def _bw53(w0, z, zr):
    """w = w₀·√(1 + (z/z_R)²) to 53 bits."""
    o = _bone_sq(z, zr)
    return _jmul(_jx(w0), (rnd53_sqrt(o[0]), False))


def bspot_ref(w0, z, zr) -> int:
    o = _bone_sq(z, zr)
    return f2b(_jfin(_jmul(_jx(w0), (rnd53_sqrt(o[0]), False), lambda t: t)))


def bpeak_ref(w0, z, zr) -> int:
    w = _bw53(w0, z, zr)
    return f2b(_jfin(_jdiv(_jx(2.0), _jmul(_jmul(_jx(BM_PI), w), w), lambda t: t)))


def broc_ref(z, zr) -> int:
    return f2b(_jfin(_jmul(_jx(z), _bone_sq(zr, z), lambda t: t)))


def bgouy_ref(z, zr) -> int:
    """arctan of the 53-bit z/z_R rounded once: below 2^-1022 it lies within |t|³ inside t."""
    v = _jdiv(_jx(z), _jx(zr))
    if v[0] == 0:
        return f2b(-0.0 if v[1] else 0.0)
    if abs(v[0]) < BM_MIN:
        return f2b(rn_frac(v[0] * (1 - Fraction(2) ** -2000)))
    if abs(v[0]) > BM_MAX:
        return f2b(mpfr_ref.cr("atan", -INF if v[0] < 0 else INF))
    return f2b(mpfr_ref.cr("atan", float(v[0])))


def bci_ref(n1, n2, rad) -> int:
    return f2b(_jfin(_jdiv(_jsub(_jx(n1), _jx(n2)), _jmul(_jx(n2), _jx(rad)), lambda t: t)))


def bdot2_ref(a, x, b, y) -> int:
    return f2b(_jfin(_jadd(_jmul(_jx(a), _jx(x)), _jmul(_jx(b), _jx(y)), lambda t: t)))


def bw0_ref(zr2, lam) -> int:
    return f2b(rn_sqrt(rnd53(rnd53(Fraction(zr2) * Fraction(lam)) / Fraction(BM_PI))))


def _bg53(length, r):
    if math.isinf(r):
        return Fraction(1), False
    return _jsub(_jx(1.0), _jdiv(_jx(length), _jx(r)))


def bgp_ref(length, r1, r2) -> int:
    return f2b(_jfin(_jmul(_bg53(length, r1), _bg53(length, r2), lambda t: t)))


def bfsr_ref(fsr, g1, g2) -> int:
    """(Δν_FSR/π)·arccos(±√(g₁g₂)) with each operation to 53 bits: a root below 2^-511 has arccos RN(π/2), as 0."""
    r53 = rnd53_sqrt(rnd53(Fraction(g1) * Fraction(g2)))
    arg = 0.0 if r53 < BM_MIN else float(r53)
    if g1 < 0.0 or g2 < 0.0:
        arg = -arg
    ac = mpfr_ref.cr("acos", arg)
    q = rnd53(Fraction(fsr) / Fraction(BM_PI))
    return f2b(rn_frac(q * Fraction(ac)) if ac != 0.0 else 0.0)


def _bsc(x, e):
    """x·2^e, ±∞ past DBL_MAX (Python's float power raises there)."""
    try:
        return math.ldexp(x, e)
    except OverflowError:
        return INF if x > 0 else -INF


def _bpos(rng, lo, hi):
    """A positive double 2^e·[1, 2) with e uniform in [lo, hi] (at least 2^-1074)."""
    e = rng.randint(lo, hi)
    return max(rng.uniform(1.0, 2.0) * 2.0 ** e, 5e-324) if e >= -1074 else 5e-324


def beam_families(rng, count):
    """bspot, bpeak, broc, bgouy, bci, bdot2, bw0, bgp and bfsr (M3 part 2d-3a): beam.cyr's tests copied into
    cr_driver.cyr, then the helpers, against the unbounded chain alone, so the calls the tests keep on Rust's chain
    check that claim: ratios around 2^512 and far past it, waists around DBL_MAX's root, quotients below 2^-1022 (exact
    midpoints of the subnormal grid among them), products past DBL_MAX cancelling and below 2^-1022, cavities whose
    L/R passes DBL_MAX, and free spectral ranges below π·2^-1022."""
    out = []
    sg = lambda: rng.choice([1.0, -1.0])                    # noqa: E731
    for i in range(count):
        zr = _bpos(rng, -1074, 1023) * sg()
        w0 = rng.choice([0.0, _bpos(rng, -1074, 1023), _bpos(rng, -560, -500), _bpos(rng, 500, 545)])
        e = rng.choice([rng.randint(505, 520), rng.randint(-40, 40), rng.randint(500, 1100)])
        z = _bsc(rng.uniform(1.0, 2.0) * abs(zr), e) * sg()
        if not math.isfinite(z) or z == 0.0:
            z = rng.uniform(-3.0, 3.0)
        out.append(("bspot", [w0, z, zr], bspot_ref(w0, z, zr)))
        if w0 != 0.0:
            out.append(("bpeak", [w0, z, zr], bpeak_ref(w0, z, zr)))
        zc = rng.choice([_bpos(rng, -99, 40), 1e-30, _bpos(rng, -99, -60)]) * sg()
        zrc = _bsc(rng.uniform(1.0, 2.0) * abs(zc), rng.randint(500, 600)) * sg()
        if not math.isfinite(zrc):
            zrc = 1.0
        out.append(("broc", [zc, zrc], broc_ref(zc, zrc)))
        if i % 3 == 0:
            j = rng.randint(0, 2000)
            u = rng.randint(1 << 40, 1 << 41) / float(1 << 40)
            zg, zrg = (2 * j + 1) * 2.0 ** -75 * u * sg(), 2.0 ** 1000 * u * sg()
        else:
            zrg = _bpos(rng, 900, 1023) * sg()
            zg = _bsc(rng.uniform(1.0, 2.0) * abs(zrg), rng.randint(-1150, -1000)) * sg()
            if not math.isfinite(zg) or zg == 0.0:
                zg = 1e-300
        out.append(("bgouy", [zg, zrg], bgouy_ref(zg, zrg)))
        n2 = rng.choice([_bpos(rng, 900, 1023), _bpos(rng, -100, 100), 1.5]) * sg()
        rad = rng.choice([_bpos(rng, 900, 1023), _bpos(rng, -99, 100)]) * sg()
        n1 = rng.choice([_bpos(rng, 1000, 1023) * sg(), -n2 * rng.uniform(0.5, 1.0), _bpos(rng, -100, 100)])
        out.append(("bci", [n1, n2, rad], bci_ref(n1, n2, rad)))
        k = i % 4
        if k == 0:
            a, x = _bpos(rng, 500, 600) * sg(), _bpos(rng, 500, 600) * sg()
            b, y = (-a, x) if rng.random() < 0.5 else (a * sg(), _bpos(rng, 500, 600))
        elif k == 1:
            a, x, b, y = (_bpos(rng, -560, -500) * sg() for _ in range(4))
        elif k == 2:
            a, x = _bpos(rng, -560, -480) * sg(), _bpos(rng, -560, -480) * sg()
            b, y = rng.uniform(-2.0, 2.0), _bpos(rng, -1100, -900) * sg()
        else:
            a, x, b, y = (rng.uniform(-3.0, 3.0) for _ in range(4))
        if not all(math.isfinite(t) for t in (a, x, b, y)):
            a, x, b, y = 1.0, 2.0, 3.0, 4.0
        out.append(("bdot2", [a, x, b, y], bdot2_ref(a, x, b, y)))
        zr2 = rng.choice([_bpos(rng, 500, 1023), _bpos(rng, -1074, -500), _bpos(rng, -100, 100)])
        lam = rng.choice([_bpos(rng, 500, 1023), _bpos(rng, -1074, -500), _bpos(rng, -100, 100)])
        out.append(("bw0", [zr2, lam], bw0_ref(zr2, lam)))
        ln = _bpos(rng, 950, 1023)
        r1 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-99, -40) * sg()
        r2 = rng.choice([ln, b2f(f2b(ln) + rng.randint(-3, 3)), 2.0 * ln, INF, -ln, r1, _bpos(rng, -99, 1023) * sg()])
        if not math.isfinite(r2) and r2 != INF:
            r2 = ln
        out.append(("bgp", [ln, r1, r2], bgp_ref(ln, r1, r2)))
        fsr = rng.choice([_bpos(rng, -1074, -1020), _bpos(rng, -1030, -1015), _bpos(rng, -100, 100)])
        if rng.random() < 0.3:
            g1 = rng.choice([1.0, -1.0]) * _bpos(rng, -1074, -100)
            g2 = math.copysign(rng.uniform(0.0, 1.0), g1)
        else:
            g1 = rng.uniform(-1.0, 1.0)
            g2 = math.copysign(rng.uniform(0.0, 1.0), g1)
        out.append(("bfsr", [fsr, g1, g2], bfsr_ref(fsr, g1, g2)))
    return out


BM_SQRT2 = struct.unpack("<d", bytes.fromhex("cd3b7f669ea0f63f"))[0]


def _bherm(n, s):
    """H_n(s) by its recurrence, each operation to 53 bits (s a pair; 2s and 2(k − 1) exact)."""
    if n == 0:
        return Fraction(1), False
    t = (2 * s[0], s[1])
    if n == 1:
        return t
    hp, hc = (Fraction(1), False), t
    for k in range(2, n + 1):
        hp, hc = hc, _jsub(_jmul(t, hc), _jmul(_jx(2.0 * (k - 1)), hp))
    return hc


def bhg_ref(w0, zr, z, m, n, x, y) -> int:
    """The Hermite–Gaussian amplitude's chain with an unbounded exponent (src/crmath_vx.cyr's _bj_vx_bm_hg): every
    operation to 53 bits, e^{−q} from MPFR at 53 bits with its widest exponent range (1 below 2^-60, the far stand-in
    from 2^47: a product whose binade lies more than 1200 below e^q's is ±0), the last product rounded once."""
    o = _bone_sq(z, zr)
    w = _jmul(_jx(w0), (rnd53_sqrt(o[0]), False))
    hm = _bherm(m, _jdiv(_jmul(_jx(BM_SQRT2), _jx(x)), w))
    hn = _bherm(n, _jdiv(_jmul(_jx(BM_SQRT2), _jx(y)), w))
    q = _jdiv(_jadd(_jmul(_jx(x), _jx(x)), _jmul(_jx(y), _jx(y))), _jmul(w, w))
    p2 = _jmul(_jmul(_jdiv(_jx(w0), w), hm), hn)
    if p2[0] == 0:
        return f2b(-0.0 if p2[1] else 0.0)
    if q[0] < Fraction(2) ** -60:
        return f2b(rn_frac(p2[0]))
    if q[0] >= Fraction(2) ** 47 or binade(p2[0]) - 1.4426950408889634 * float(q[0]) < -1200:
        return f2b(-0.0 if p2[0] < 0 else 0.0)
    es, ee = mpfr_ref.wide_split("exp", -q[0])
    return f2b(rn_frac(p2[0] * Fraction(es) * Fraction(2) ** ee))


def _blag(p, l, x):
    """L_p^l(x) by its recurrence, each operation to 53 bits (x a pair; the integer parts exact)."""
    if p == 0:
        return Fraction(1), False
    lf = float(l)
    lc = _jsub(_jx(1.0 + lf), x)
    if p == 1:
        return lc
    lp = (Fraction(1), False)
    for k in range(2, p + 1):
        kf = float(k)
        a = _jmul(_jsub(_jx(2.0 * kf + lf - 1.0), x), lc)
        b = _jmul(_jx(kf + lf - 1.0), lp)
        lp, lc = lc, _jdiv(_jsub(a, b), _jx(kf))
    return lc


def _bme(a):
    """A pair as (m, e, negative) with 1 ≤ |m| < 2 (m = 0 for a zero), for products past what a Fraction holds."""
    if a[0] == 0:
        return Fraction(0), 0, a[1]
    e = binade(a[0])
    return a[0] / Fraction(2) ** e, e, a[0] < 0


def _bme_mul(a, b):
    m = a[0] * b[0]
    if m == 0:
        return Fraction(0), 0, a[2] != b[2]
    m, e = rnd53(m), a[1] + b[1]
    if abs(m) >= 2:
        m, e = m / 2, e + 1
    return m, e, m < 0


def blg_ref(w0, zr, z, p, la, r) -> int:
    """The Laguerre–Gaussian intensity's chain with an unbounded exponent (src/crmath_vx.cyr's _bj_vx_bm_lg): every
    operation to 53 bits, the power in __powidf2's order and the products on (m, e) triples, e^{−q} from MPFR at 53 bits
    with its widest exponent range (1 below 2^-60, the far stand-in from 2^47), the last product rounded once."""
    o = _bone_sq(z, zr)
    w = _jmul(_jx(w0), (rnd53_sqrt(o[0]), False))
    rho = _jdiv(_jmul(_jx(BM_SQRT2), _jx(r)), w)
    x = _jmul(rho, rho)
    lpv = _bme(_blag(p, la, x))
    q = _jdiv(_jmul((2 * Fraction(r), r < 0), _jx(r)), _jmul(w, w))
    sc = _bme(_jdiv(_jmul(_jx(w0), _jx(w0)), _jmul(w, w)))
    pw, a, n = (Fraction(1), 0, False), _bme(x), la
    while n:
        if n & 1:
            pw = _bme_mul(pw, a)
        n >>= 1
        if n:
            a = _bme_mul(a, a)
    p3 = _bme_mul(_bme_mul(_bme_mul(sc, pw), lpv), lpv)
    zero = f2b(-0.0 if p3[2] else 0.0)
    if p3[0] == 0:
        return zero
    if q[0] < Fraction(2) ** -60:
        es, ee = 1.0, 0
    elif q[0] >= Fraction(2) ** 47:
        return zero
    else:
        es, ee = mpfr_ref.wide_split("exp", -q[0])
    e = p3[1] + ee
    if e > 1100:
        return f2b(-INF if p3[2] else INF)
    if e < -1200:
        return zero
    return f2b(rn_frac(p3[0] * Fraction(es) * Fraction(2) ** e))


def bhgf_ref(lam, w0, zr, x, y, z, m, n):
    """The Hermite–Gaussian field's chain with an unbounded exponent (src/crmath_vx.cyr's _bj_vx_bm_hgf): the amplitude
    to 53 bits (an exponent past ∓1300 decides ±0 or ±∞ alone, |cos φ| and |sin φ| exceeding 2^-80 but at φ = 0), k,
    each term of φ and ζ to 53 bits, cos φ and sin φ from MPFR at 53 bits with its widest exponent range, the two last
    products rounded once: [re, −im] as bits."""
    o = _bone_sq(z, zr)
    w = _jmul(_jx(w0), (rnd53_sqrt(o[0]), False))
    hm = _bherm(m, _jdiv(_jmul(_jx(BM_SQRT2), _jx(x)), w))
    hn = _bherm(n, _jdiv(_jmul(_jx(BM_SQRT2), _jx(y)), w))
    q = _jdiv(_jadd(_jmul(_jx(x), _jx(x)), _jmul(_jx(y), _jx(y))), _jmul(w, w))
    p2 = _jmul(_jmul(_jdiv(_jx(w0), w), hm), hn)
    kind, amp = "v", p2
    if p2[0] != 0 and q[0] >= Fraction(2) ** -60:
        if q[0] >= Fraction(2) ** 47:
            kind, amp = "tiny", p2[0] < 0
        else:
            lb = binade(p2[0]) - 1.4426950408889634 * float(q[0])
            if lb < -1300:
                kind, amp = "tiny", p2[0] < 0
            elif lb > 1300:
                kind, amp = "huge", p2[0] < 0
            else:
                es, ee = mpfr_ref.wide_split("exp", -q[0])
                amp = _jmul(p2, (Fraction(es) * Fraction(2) ** ee, False))
    k = _jdiv(_jx(2.0 * BM_PI), _jx(lam))
    zp = _jx(z)
    t1 = _jmul(k, zp)
    inv = _jdiv(zp, _jadd(_jmul(zp, zp), _jmul(_jx(zr), _jx(zr))))
    r2 = _jadd(_jmul(_jx(x), _jx(x)), _jmul(_jx(y), _jx(y)))
    t2 = _jmul(_jmul((k[0] / 2, k[1]), r2), inv)
    qz = _jdiv(zp, _jx(zr))
    g = qz if qz[0] == 0 else (mpfr_ref.wide("atan", qz[0]), qz[0] < 0)
    t3 = _jmul(_jx(float(m) + float(n) + 1.0), g)
    ph = _jsub(_jadd(t1, t2), t3)
    if ph[0] == 0:
        c, sn = (Fraction(1), False), ph
    else:
        cv, sv = mpfr_ref.wide("cos", ph[0]), mpfr_ref.wide("sin", ph[0])
        c, sn = (cv, cv < 0), (sv, sv < 0)
    out = []
    for f in (c, sn):
        if kind == "v":
            v = _jfin(_jmul(amp, f, lambda t: t))
        else:
            neg = amp != f[1]
            if kind == "huge" and f[0] != 0:
                v = -INF if neg else INF
            else:
                v = -0.0 if neg else 0.0
        out.append(v)
    return [f2b(out[0]), f2b(-out[1])]


def hgf_families(rng, count):
    """bhgf (M3 part 2d-3c): the HG field's unbounded chain on beams from 2^±600 with wavelengths from 2^-1074 to
    2^1023 (k past DBL_MAX), z from z_R·2^-40 to 2^1000 and to z_R·2^1000, points to √2x/w = 2^6, orders to 12: phases
    past DBL_MAX and to 2^4100, envelopes to the stand-in, each part asked separately."""
    out = []
    for i in range(count):
        lam = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1074, 1023)
        if lam == 0.0 or not math.isfinite(lam):
            lam = 1e-6
        w0 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-600, 600)
        zr = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1074, 600)
        if zr == 0.0:
            zr = 1.0
        z = rng.choice([0.0, zr * rng.uniform(-3.0, 3.0), zr * 2.0 ** rng.randint(-40, 40),
                        rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-1074, 1000)]) * rng.choice([1.0, -1.0])
        if not math.isfinite(z):
            z = zr
        m, n = rng.randint(0, 12), rng.randint(0, 12)
        wv = abs(w0) * math.sqrt(1.0 + min((z / zr) ** 2, 1e300)) if abs(z / zr) < 1e150 else abs(w0)

        def pt():
            if rng.random() < 0.2:
                return 0.0
            v = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-30, 6) * wv / BM_SQRT2 * rng.choice([1.0, -1.0])
            return v if math.isfinite(v) and v != 0.0 else 1.0
        x, y = pt(), pt()
        ref = bhgf_ref(lam, w0, zr, x, y, z, m, n)
        for part in (0, 1):
            out.append(("bhgf", [lam, w0, zr, x, y, z, m, n, part], ref[part]))
    return out


def lg_families(rng, count):
    """blg (M3 part 2d-3b-2): the LG intensity's unbounded chain on beams from 2^±600 (z_R independent of w₀), the
    waist, near and far field, radial orders to 60 and some to 300, azimuthal ones to 60 and some to 2³¹ (the power's
    squares), and radii from ρ² = 2^-120 to 2^80, envelopes past e^-2300 and past the stand-in among them."""
    out = []
    for i in range(count):
        w0 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-600, 600)
        zr = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-600, 600)
        z = rng.choice([0.0, zr * rng.uniform(-3.0, 3.0), zr * 2.0 ** rng.randint(0, 900)]) * rng.choice([1.0, -1.0])
        if not math.isfinite(z):
            z = zr
        p = rng.choice([rng.randint(0, 12), rng.randint(0, 60), rng.randint(100, 300) if i % 9 == 0 else 2])
        la = rng.choice([rng.randint(0, 12), rng.randint(0, 60), rng.randint(1, 1 << 31) if i % 11 == 0 else 1])
        wv = abs(w0) * math.sqrt(1.0 + min((z / zr) ** 2, 1e300)) if abs(z / zr) < 1e150 else abs(w0)
        c = rng.random()
        if c < 0.08:
            r = 0.0
        else:
            x = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-120, 80 if c < 0.3 else 12)
            r = math.sqrt(x / 2.0) * wv * rng.choice([1.0, -1.0])
            if not math.isfinite(r) or r == 0.0:
                r = 1.0
        out.append(("blg", [w0, zr, z, p, la, r], blg_ref(w0, zr, z, p, la, r)))
    return out


def hg_families(rng, count):
    """bhg (M3 part 2d-3b-1): the HG amplitude's unbounded chain on beams from 2^±600 (z_R independent of w₀, as a
    hand-built beam's may be), the waist, the near and the far field (w past DBL_MAX), orders to 60 and some to 400,
    and transverse points from √2x/w = 2^-60 to 2^40, envelopes past e^-2300 and past the stand-in among them."""
    out = []
    for i in range(count):
        w0 = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-600, 600)
        zr = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-600, 600)
        z = rng.choice([0.0, zr * rng.uniform(-3.0, 3.0), zr * 2.0 ** rng.randint(0, 900)]) * rng.choice([1.0, -1.0])
        if not math.isfinite(z):
            z = zr
        m = rng.choice([rng.randint(0, 12), rng.randint(0, 60), rng.randint(100, 400) if i % 7 == 0 else 3])
        n = rng.choice([rng.randint(0, 12), rng.randint(0, 60), 0, 1])
        wv = abs(w0) * math.sqrt(1.0 + min((z / zr) ** 2, 1e300)) if abs(z / zr) < 1e150 else abs(w0)

        def pt():
            c = rng.random()
            if c < 0.1:
                return 0.0
            s = rng.uniform(1.0, 2.0) * 2.0 ** rng.randint(-60, 40 if c < 0.3 else 6)
            v = s * wv / BM_SQRT2 * rng.choice([1.0, -1.0])
            return v if math.isfinite(v) and v != 0.0 else 1.0
        x, y = pt(), pt()
        out.append(("bhg", [w0, zr, z, m, n, x, y], bhg_ref(w0, zr, z, m, n, x, y)))
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
        trng = random.Random(0xC0_57)
        calls += [(f, ref, ws) for f, ws, ref in cos_t_families(trng, max(1, n // 5))
                  if "exact" in names or f in names]
        lrng = random.Random(0x1E_4A_2C4)
        calls += [(f, ref, ws) for f, ws, ref in lw_families(lrng, max(1, n // 10))
                  if "exact" in names or f in names]
        xrng = random.Random(0x2D_3E)
        calls += [(f, ref, ws) for f, ws, ref in lw_ux_families(xrng, max(1, n // 10))
                  if "exact" in names or f in names]
        nrng = random.Random(0x2D_1B)
        calls += [(f, ref, ws) for f, ws, ref in nb_families(nrng, max(1, n // 30))
                  if "exact" in names or f in names]
        vrng = random.Random(0x2D_2C)
        calls += [(f, ref, ws) for f, ws, ref in vx_families(vrng, max(1, n // 10))
                  if "exact" in names or f in names]
        orng = random.Random(0x2D_3F)
        calls += [(f, ref, ws) for f, ws, ref in lwo_families(orng, max(1, n // 20))
                  if "exact" in names or f in names]
        grng2 = random.Random(0x2D_3F6)
        calls += [(f, ref, ws) for f, ws, ref in gammao_families(grng2, max(1, n // 10))
                  if "exact" in names or f in names]
        qrng = random.Random(0x3_2C)
        calls += [(f, ref, ws) for f, ws, ref in sumsq_families(qrng, max(1, n // 5))
                  if "exact" in names or f in names]
        arng = random.Random(0x3_2C2)
        calls += [(f, ref, ws) for f, ws, ref in abcd_families(arng, max(1, n // 20))
                  if "exact" in names or f in names]
        crng = random.Random(0x4_2A2)
        calls += [(f, ref, ws) for f, ws, ref in cascade_families(crng, max(1, n // 50))
                  if "exact" in names or f in names]
        zrng = random.Random(0x4_2A3)
        calls += [(f, ref, ws) for f, ws, ref in s2z_families(zrng, max(1, n // 50))
                  if "exact" in names or f in names]
        mrng = random.Random(0x5_2A1A)
        calls += [(f, ref, ws) for f, ws, ref in rat_families(mrng, max(1, n // 5))
                  if "exact" in names or f in names]
        jrng = random.Random(0x3_2D1B)
        calls += [(f, ref, ws) for f, ws, ref in jones_families(jrng, max(1, n // 10))
                  if "exact" in names or f in names]
        rrng = random.Random(0x3_2D2A)
        calls += [(f, ref, ws) for f, ws, ref in rayleigh_families(rrng, max(1, n // 10))
                  if "exact" in names or f in names]
        frng = random.Random(0x3_2D3C)
        calls += [(f, ref, ws) for f, ws, ref in hgf_families(frng, max(1, n // 50))
                  if "exact" in names or f in names]
        lrng = random.Random(0x3_2D3B2)
        calls += [(f, ref, ws) for f, ws, ref in lg_families(lrng, max(1, n // 50))
                  if "exact" in names or f in names]
        hrng = random.Random(0x3_2D3B1)
        calls += [(f, ref, ws) for f, ws, ref in hg_families(hrng, max(1, n // 50))
                  if "exact" in names or f in names]
        brng = random.Random(0x3_2D3A)
        calls += [(f, ref, ws) for f, ws, ref in beam_families(brng, max(1, n // 10))
                  if "exact" in names or f in names]
        srng = random.Random(0x3_2D1C)
        calls += [(f, ref, ws) for f, ws, ref in stokes_families(srng, max(1, n // 10))
                  if "exact" in names or f in names]
        drng = random.Random(0x3_2D1A)
        calls += [(f, ref, ws) for f, ws, ref in dot4_families(drng, max(1, n // 2))
                  if "exact" in names or f in names]
        prng2 = random.Random(0x9A_2C5)
        calls += [(f, ref, ws) for f, ws, ref in par_families(prng2, max(1, n // 10))
                  if "exact" in names or f in names]
    if "ux" in names or any(x in UX for x in names):
        urng = random.Random(0x5C_2E_D6E)
        fams = ux_families(urng, max(1, n // 20))
        calls += [(f, want, ws) for f, ws, want in fams if "ux" in names or f in names]
    too_long = [name for name, _, a in calls if len(a) > 24]
    if too_long:                                # cr_driver.cyr reads at most twenty-four words a call
        print(f"cr_verify: {len(too_long)} calls pass more than twenty-four words ({too_long[0]})")
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
            # "invalid" is a word for these
            any_nan = name not in UX and not name.startswith(("gamma", "cm", "cos_t", "lw_", "lwo_", "par", "vx_",
                                                              "abcd", "cascade", "dot4", "j", "s", "r") + BM_FAMS)
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
