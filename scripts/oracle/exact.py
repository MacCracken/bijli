"""Exact binary64 rounding of rationals and their square roots, shared by cr_verify.py and validate.py.

rn_frac(v): RN(v) for a Fraction v, rounded once to nearest even onto binary64 (CPython's int division does exactly
that, the subnormal grid included), ±∞ past DBL_MAX's rounding boundary. rn_sqrt(v): RN(√v) for a Fraction v > 0,
exactly (an integer square root and an exact midpoint comparison). rn_mono(num, den): RN(Π num / Π den) for finite
doubles, a zero factor giving IEEE's signed zero (the xor of the factors' signs).
"""
import math
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
