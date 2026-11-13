"""Exact binary64 rounding of rationals and their square roots, shared by cr_verify.py and validate.py.

rn_frac(v): RN(v) for a Fraction v, rounded once to nearest even onto binary64 (CPython's int division does exactly
that, the subnormal grid included), ±∞ past DBL_MAX's rounding boundary. rn_sqrt(v): RN(√v) for a Fraction v > 0,
exactly (an integer square root and an exact midpoint comparison). rn_mono(num, den): RN(Π num / Π den) for finite
doubles, a zero factor giving IEEE's signed zero (the xor of the factors' signs). rn_lw(r, b, mode): the
Liénard–Wiechert RN(|R| − R⋅β) (mode 0) or RN(1 − R⋅β/|R|) (mode 1) exactly, and lw_sign the exact signs it rests on;
rn_nb(r, b, i): RN(Rᵢ/|R| − βᵢ), a component of n̂ − β, with nb_sign.
"""
import math
import struct
from fractions import Fraction

INF = float("inf")


def rn_frac(v: Fraction) -> float:
    if v == 0:
        return 0.0
    try:
        return v.numerator / v.denominator
    except OverflowError:
        return INF if v > 0 else -INF


def rn_sqrt(v: Fraction) -> float:
    p, q = v.numerator, v.denominator
    lf = p.bit_length() - q.bit_length()
    if (p << max(0, -lf)) < (q << max(0, lf)):
        lf -= 1                                 # 2^lf ≤ v < 2^(lf + 1)
    e = lf // 2                                 # 2^e ≤ √v < 2^(e + 1)
    k = max(e - 52, -1074)
    num, den = (p, q << (2 * k)) if k >= 0 else (p << (-2 * k), q)
    m = math.isqrt(num // den)
    lhs, rhs = 4 * num, (2 * m + 1) ** 2 * den
    if lhs > rhs or (lhs == rhs and m % 2 == 1):
        m += 1
    r = Fraction(m) * Fraction(2) ** k
    return INF if r >= Fraction(2) ** 1024 else float(r)


def rn_mono(num, den) -> float:
    v = Fraction(1)
    for x in num:
        v *= Fraction(x)
    for y in den:
        v /= Fraction(y)
    if v == 0:
        neg = sum(math.copysign(1.0, x) < 0 for x in list(num) + list(den)) % 2
        return -0.0 if neg else 0.0
    return rn_frac(v)


def _log2_floor(p: int, q: int) -> int:
    """⌊log₂(p/q)⌋ for positive integers p and q."""
    lf = p.bit_length() - q.bit_length()
    if (p << max(0, -lf)) < (q << max(0, lf)):
        lf -= 1
    return lf


def rnd53(v: Fraction) -> Fraction:
    """v rounded to 53 significant bits, to nearest even, with an unbounded exponent (no subnormal, no overflow)."""
    if v == 0:
        return v
    p, q = abs(v.numerator), v.denominator
    k = _log2_floor(p, q) - 52
    num, den = (p, q << k) if k >= 0 else (p << -k, q)
    m, r = divmod(num, den)
    if 2 * r > den or (2 * r == den and m % 2 == 1):
        m += 1
    out = Fraction(m) * Fraction(2) ** k
    return out if v > 0 else -out


def rnd53_sqrt(v: Fraction) -> Fraction:
    """√v rounded to 53 significant bits, to nearest even, with an unbounded exponent, for v ≥ 0."""
    if v == 0:
        return v
    p, q = v.numerator, v.denominator
    k = _log2_floor(p, q) // 2 - 52
    num, den = (p, q << (2 * k)) if k >= 0 else (p << (-2 * k), q)
    m = math.isqrt(num // den)
    lhs, rhs = 4 * num, (2 * m + 1) ** 2 * den
    if lhs > rhs or (lhs == rhs and m % 2 == 1):
        m += 1
    return Fraction(m) * Fraction(2) ** k


def binade(v: Fraction) -> int:
    """⌊log₂|v|⌋ for v ≠ 0."""
    return _log2_floor(abs(v.numerator), v.denominator)


def lw_sign(r, b, mode, mu) -> int:
    """The sign (−1, 0 or 1), exactly, of |R| − R⋅β − μ (mode 0) or of 1 − R⋅β/|R| − μ (mode 1), for a finite R ≠ 0,
    finite β and a Fraction μ (by the signs of the parts, or of the squared comparison)."""
    s = sum(Fraction(x) ** 2 for x in r)
    p = sum(Fraction(x) * Fraction(y) for x, y in zip(r, b))
    if mode == 0:
        bb = p + mu                             # sign(√S − (P + μ))
        if bb <= 0:
            return 1
        d = s - bb * bb
    else:
        a = 1 - mu                              # sign(a√S − P)
        if a >= 0 and p <= 0:
            return 0 if a == 0 and p == 0 else 1
        if a <= 0 and p >= 0:
            return -1
        d = a * a * s - p * p if a > 0 else p * p - a * a * s
    return (d > 0) - (d < 0)


def rn_lw(r, b, mode) -> float:
    """RN(|R| − R⋅β) (mode 0, +∞ past DBL_MAX's rounding boundary) or RN(1 − R⋅β/|R|) (mode 1), exactly, for a finite
    R ≠ 0 and |β| < 1: an approximation within 2^-433 of the value (relative; it lies above 2^-367), rounded, then
    moved by the exact signs at its neighbouring midpoints, ties to even."""
    s = sum(Fraction(x) ** 2 for x in r)
    p = sum(Fraction(x) * Fraction(y) for x, y in zip(r, b))
    k = max(0, (1400 - binade(s)) // 2 + 1)
    sq = Fraction(math.isqrt(math.floor(s * 4 ** k)), 2 ** k)
    v = sq - p if mode == 0 else 1 - p / sq
    y = rn_frac(v)
    top = Fraction(2) ** 1024 - Fraction(2) ** 970     # the boundary past DBL_MAX
    for _ in range(4):
        if y == INF:
            return INF if lw_sign(r, b, mode, top) >= 0 else float.fromhex("0x1.fffffffffffffp+1023")
        up = math.nextafter(y, INF)
        mu = top if up == INF else (Fraction(y) + Fraction(up)) / 2
        c = lw_sign(r, b, mode, mu)
        if c > 0:
            y = up
            continue
        if c == 0:
            return up if (struct.unpack("<Q", struct.pack("<d", y))[0] & 1) else y
        dn = math.nextafter(y, -INF)
        c = lw_sign(r, b, mode, (Fraction(y) + Fraction(dn)) / 2)
        if c < 0:
            y = dn
            continue
        if c == 0:
            return dn if (struct.unpack("<Q", struct.pack("<d", dn))[0] & 1) == 0 else y
        return y
    raise AssertionError("rn_lw: the candidate did not settle")


def lw_below(r, b) -> bool:
    """κ = 1 − R⋅β/|R| < 10^-30 exactly (lienard_wiechert_e's DivisionByZero)."""
    return lw_sign(r, b, 1, Fraction(1, 10 ** 30)) < 0


def nb_sign(r, b, i, sg, mu) -> int:
    """The sign (−1, 0 or 1), exactly, of sg·(Rᵢ/|R| − βᵢ) − μ for Rᵢ ≠ 0 and μ ≥ 0: that of R̃ − B·|R|, R̃ = sg·Rᵢ and
    B = sg·βᵢ + μ, by their signs or by R̃² − B²·|R|²."""
    s = sum(Fraction(x) ** 2 for x in r)
    rt = sg * Fraction(r[i])
    bb = sg * Fraction(b[i]) + mu
    if rt > 0 and bb <= 0:
        return 1
    if rt < 0 and bb >= 0:
        return -1
    d = rt * rt - bb * bb * s
    c = (d > 0) - (d < 0)
    return c if rt > 0 else -c


def rn_nb(r, b, i) -> float:
    """RN(Rᵢ/|R| − βᵢ) exactly, for a finite R with Rᵢ ≠ 0 and |β| < 1: an approximation within 2^-699 of the value
    (relative; through (Rᵢ² − βᵢ²|R|²)/(|R|(Rᵢ + βᵢ|R|)) when Rᵢ and βᵢ share a sign), rounded, then moved by the exact
    signs at the neighbouring midpoints of its magnitude, ties to even."""
    s = sum(Fraction(x) ** 2 for x in r)
    x, y = Fraction(r[i]), Fraction(b[i])
    k = max(0, (1400 - binade(s)) // 2 + 1)
    sq = Fraction(math.isqrt(math.floor(s * 4 ** k)), 2 ** k)
    if y == 0 or (x < 0) != (y < 0):
        v = x / sq - y
    else:
        v = (x * x - y * y * s) / (sq * (x + y * sq))
    sg = 1 if v > 0 else -1
    a = abs(rn_frac(abs(v)))
    for _ in range(12):
        up = math.nextafter(a, INF)
        c = nb_sign(r, b, i, sg, (Fraction(a) + Fraction(up)) / 2)
        if c > 0:
            a = up
            continue
        if c == 0:
            a = up if (struct.unpack("<Q", struct.pack("<d", a))[0] & 1) else a
            break
        if a == 0.0:
            break
        dn = math.nextafter(a, -INF)
        c = nb_sign(r, b, i, sg, (Fraction(a) + Fraction(dn)) / 2)
        if c < 0:
            a = dn
            continue
        if c == 0:
            a = dn if (struct.unpack("<Q", struct.pack("<d", dn))[0] & 1) == 0 else a
        break
    return a if sg > 0 else -a
