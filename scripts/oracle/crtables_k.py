"""Copied from kimiya's scripts/oracle/mathx_tables.py (its ADR 0006); bijli renders only the ln and atan
tables that src/crmath_k.cyr reads (render() below), under the `_BJ_K_CR_` names.

Generate src/mathx_data.cyr: the tables and constants of the correctly rounded exp, expm1, exp10, ln
and log10 (ADR 0006, M3c-2a), sin, cos, atan2 and acos (M4a), x·2^(1/6) (M4d) and 1/(100·x) (M4e) in
src/mathx.cyr.

Every value is computed from MPFR (through scripts/oracle/mpfr_ref.py's ctypes binding, at 400 bits)
or from exact rational arithmetic (Python's Fraction), then split into doubles, each the round-to-
nearest of what is left. So a k-word entry is the value to about 2^(-53k) relative. The output is
deterministic: no randomness, no time stamps, no dict ordering.

Usage, from the repo root:
    python3 -I scripts/oracle/mathx_tables.py            # write src/mathx_data.cyr
    python3 -I scripts/oracle/mathx_tables.py --check    # exit 1 if src/mathx_data.cyr is stale

Besides writing the values, the generator checks the properties the code's error analysis relies on
(asserting, so a table that would break a proof is never written):
- exp: A_j has at most 26 significant bits and A_j + B_j + C_j is 2^(j/512) to 2^-128 relative;
  L1 has at most 33 bits (K·L1 is exact for |K| < 2^20); the hi part of ln 10 splits 26 + 27 bits.
- ln: for each of the 512 intervals of m, r = m·c_j − 1 is exactly representable for every double m
  in the interval and |r| < 2^-9 (a bit-count argument on the interval's end points); c = 1 and
  L = 0 on the interval holding 1; |L_j| ≥ max |r| elsewhere (the Fast2Sum precondition).
- sin and cos: the 26-bit heads of sin(i/256) and cos(i/256) (A·t1 exact), |B| ≤ 2^-26·|A|, A + B to 2^-78
  and A + B + C + D to 2^-180; 201 rows cover |r| ≤ π/4 + 2^-30; π/2's Cody–Waite parts P1, P2 have at most
  33 bits, |P2| < 2^-32, |P3| < 2^-65 and the tail below 2^-225; the bits of 2/π agree with Machin's
  formula computed in integers and cover every exponent; no double is closer than 4.687·10^-19 to a nonzero
  multiple of π/2 (a continued-fraction bound in each binade, attained by Lefèvre and Muller's worst case
  6381956970095103·2^797: the reduction's precision rests on it).
- atan: atan(j/256) as triple-doubles to 2^-158.
- 2^(1/6) (M4d, ADR 0002 #16): H + M + L to 2^-158 from an exact integer root, H Rust's RN(2^(1/6)), and the
  continued-fraction bound that makes the triple-double product x·(H + M + L) round correctly once.
- 1/100 (M4e, ADR 0002 #19): H + L to 2^-116 from the exact rational, H Rust's 0.01 (the wavenumber conversions'
  correctly rounded 1/(100·x), _km_cr_hundredth_div).
"""
from fractions import Fraction
import ctypes
import struct
import sys
from pathlib import Path

sys.dont_write_bytecode = True   # no __pycache__ in scripts/oracle/
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpfr_ref  # noqa: E402  (the binding: emin/emax already set for binary64)

_m = mpfr_ref._m
_P = ctypes.c_void_p
_RNDN = 0
PREC = 400

_m.mpfr_exp2.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_log.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_const_log2.argtypes = [_P, ctypes.c_int]
_m.mpfr_sub.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_add.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_mul.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_div.argtypes = [_P, _P, _P, ctypes.c_int]
_m.mpfr_set.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_neg.argtypes = [_P, _P, ctypes.c_int]
_m.mpfr_clear.argtypes = [_P]


class _V:
    """An MPFR number at `prec` bits (400 by default)."""

    def __init__(self, prec=PREC):
        self.t = (ctypes.c_byte * 32)()
        _m.mpfr_init2(ctypes.byref(self.t), prec)

    @property
    def p(self):
        return ctypes.byref(self.t)

    def __del__(self):
        _m.mpfr_clear(ctypes.byref(self.t))


def _d(x):
    v = _V()
    _m.mpfr_set_d(v.p, x, _RNDN)
    return v


def _op(fn, a, b=None):
    r = _V()
    if b is None:
        getattr(_m, fn)(r.p, a.p, _RNDN)
    else:
        getattr(_m, fn)(r.p, a.p, b.p, _RNDN)
    return r


def _ln2():
    r = _V()
    _m.mpfr_const_log2(r.p, _RNDN)
    return r


def _get(v):
    return _m.mpfr_get_d(v.p, _RNDN)


def _round_bits(v, bits):
    """v rounded to nearest at `bits` significant bits, as a double (exact for bits <= 53)."""
    r = _V(bits)
    _m.mpfr_set(r.p, v.p, _RNDN)
    return _m.mpfr_get_d(r.p, _RNDN)


def _split(v, n, first_bits=53):
    """n doubles: the first rounded to `first_bits` bits, each next the nearest double to what is left."""
    out = []
    rest = _op("mpfr_set", v)
    for i in range(n):
        h = _round_bits(rest, first_bits) if i == 0 else _get(rest)
        out.append(h)
        rest = _op("mpfr_sub", rest, _d(h))
    return out


def _fsplit(q, n, first_bits=53):
    """The same split of an exact Fraction."""
    out = []
    rest = q
    for i in range(n):
        if i == 0 and first_bits < 53:
            h = _fround_bits(rest, first_bits)
        else:
            h = float(rest)       # Fraction -> float is correctly rounded
        out.append(h)
        rest = rest - Fraction(h)
    return out


def _fround_bits(q, bits):
    """A Fraction rounded to nearest-even at `bits` significant bits."""
    if q == 0:
        return 0.0
    s = -1 if q < 0 else 1
    a = abs(q)
    e = a.numerator.bit_length() - a.denominator.bit_length()
    while Fraction(2) ** e > a:
        e -= 1
    while Fraction(2) ** (e + 1) <= a:
        e += 1
    scale = Fraction(2) ** (bits - 1 - e)
    n = a * scale
    fl = n.numerator // n.denominator
    rem = n - fl
    if rem > Fraction(1, 2) or (rem == Fraction(1, 2) and fl % 2 == 1):
        fl += 1
    return s * float(Fraction(fl) / scale)


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def _sigbits(x):
    """The number of significant bits of a double (0 for 0)."""
    if x == 0.0:
        return 0
    f = Fraction(abs(x))
    n, d = f.numerator, f.denominator
    while n % 2 == 0:
        n //= 2
    return n.bit_length()


def _exact(xs):
    return sum((Fraction(x) for x in xs), Fraction(0))


def _rel(approx, v):
    """|approx - v| / |v| as a float, for an exact Fraction approx and an MPFR v."""
    d = _op("mpfr_sub", _d_frac(approx), v)
    return abs(_get(d) / _get(v))


def _d_frac(q):
    """An MPFR number holding the Fraction q (a sum of a few doubles) exactly, at 400 bits."""
    acc = _d(0.0)
    rest = q
    for _ in range(12):
        if rest == 0:
            break
        h = float(rest)
        acc = _op("mpfr_add", acc, _d(h))
        rest -= Fraction(h)
    assert rest == 0, "Fraction is not a short sum of doubles"
    return acc


# ── exp: 2^(j/512) = A + B + C, the reduction and the polynomials ────────────────────────────────

N_EXP = 512


def exp_table():
    rows = []
    for j in range(N_EXP):
        v = _op("mpfr_exp2", _d(j / N_EXP))
        a, b, c, d = _split(v, 4, first_bits=26)
        assert _sigbits(a) <= 26
        assert _rel(_exact([a, b, c]), v) < 2.0 ** -128, j
        assert _rel(_exact([a, b, c, d]), v) < 2.0 ** -180, j
        if j == 0:
            assert (a, b, c, d) == (1.0, 0.0, 0.0, 0.0)
        rows.append((a, b, c, d))
    return rows


def exp_consts():
    ln2 = _ln2()
    ln2n = _op("mpfr_div", ln2, _d(float(N_EXP)))
    inv = _get(_op("mpfr_div", _d(float(N_EXP)), ln2))
    l1, l2, l3, l4 = _split(ln2n, 4, first_bits=33)
    assert _sigbits(l1) <= 33
    ln10 = _op("mpfr_log", _d(10.0))
    t10 = _split(ln10, 3)
    h1, h2 = _fsplit(Fraction(t10[0]), 2, first_bits=26)
    assert _sigbits(h1) <= 26 and _sigbits(h2) <= 27 and h1 + h2 == t10[0]
    return {
        "INVLN2N": inv, "L1": l1, "L2": l2, "L3": l3, "L4": l4,
        "P3": float(Fraction(1, 6)), "P4": float(Fraction(1, 24)), "P5": float(Fraction(1, 120)),
        "LN10H": t10[0], "LN10M": t10[1], "LN10L": t10[2], "LN10H1": h1, "LN10H2": h2,
    }


def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def expm1_poly():
    """1/3!, 1/4!, ..., 1/9!: the coefficients of expm1's small-|x| fast path, x³·S(x)."""
    return [float(Fraction(1, _fact(k))) for k in range(3, 10)]


def exp_acc():
    """1/k! for k = 2 ... 17 as triple-doubles: the accurate paths' Taylor series."""
    return [_fsplit(Fraction(1, _fact(k)), 3) for k in range(2, 18)]


# ── ln: m·c_j ≈ 1 over 512 intervals of m in [0.7495, 1.4990), −ln c_j as a triple-double ───────

N_LN = 512
LN_OFF = 0x3FE7FC0000000000   # bits(1) − 256.5·2^43: interval 256 is centred (in bits) on 1


def _m_of_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def ln_table():
    rows = []
    for j in range(N_LN):
        lo_b = LN_OFF + j * (1 << 43)
        hi_b = lo_b + (1 << 43) - 1                 # the last double of the interval
        m_lo, m_hi = _m_of_bits(lo_b), _m_of_bits(hi_b)
        if j == 256:
            c = 1.0
        else:
            mid = (Fraction(m_lo) + Fraction(m_hi)) / 2
            inv = 1 / mid
            # c on the grid 2^-9 (c >= 1, m < 1) or 2^-10 (c < 1, m >= 1): r = m·c − 1 then fits 53 bits.
            g = 9 if m_hi < 1.0 else 10
            c = float(Fraction(round(inv * 2 ** g), 2 ** g))
        cf = Fraction(c)
        r_lo = Fraction(m_lo) * cf - 1
        r_hi = Fraction(m_hi) * cf - 1
        rmax = max(abs(r_lo), abs(r_hi))
        assert rmax < Fraction(1, 2 ** 9), (j, float(rmax))
        # Exactness of r = m·c − 1 for every double m of the interval: r is a multiple of
        # ulp(m)·ulp(c), and |r| < 2^-9, so it has at most log2(2^-9 / (ulp(m)·ulp(c))) bits.
        for mm in (m_lo, m_hi):
            if j == 256:
                continue
            um = Fraction(2) ** (_exp2(mm) - 52)
            uc = _ulp_grid(c)
            nbits = (Fraction(1, 2 ** 9) / (um * uc)).numerator.bit_length() - 1
            assert nbits <= 53, (j, nbits)
            assert _exp2(mm) == _exp2(m_lo), j    # the interval stays in one binade
        if j == 256:
            lh, ll, lt = 0.0, 0.0, 0.0
            assert m_lo < 1.0 < m_hi
        else:
            lv = _op("mpfr_neg", _op("mpfr_log", _d(c)))
            lh, ll, lt = _split(lv, 3)
            assert abs(lh) >= float(rmax), j      # Fast2Sum(L, r) precondition
        rows.append((c, lh, ll, lt))
    return rows


def _exp2(x):
    """The binary exponent of a positive double."""
    return ((bits(x) >> 52) & 0x7FF) - 1023


def _ulp_grid(c):
    """The coarsest power of two c is a multiple of."""
    f = Fraction(c)
    g = Fraction(1)
    while (f / g).denominator != 1:
        g /= 2
    while (f / (g * 2)).denominator == 1:
        g *= 2
    return g


def ln_consts():
    ln2 = _ln2()
    h, mid, lo = _split(ln2, 3, first_bits=42)
    assert _sigbits(h) <= 42
    ln10 = _op("mpfr_log", _d(10.0))
    iv = _op("mpfr_div", _d(1.0), ln10)
    ivh, ivm, ivl = _split(iv, 3)
    i1, i2 = _fsplit(Fraction(ivh), 2, first_bits=26)
    assert _sigbits(i1) <= 26 and _sigbits(i2) <= 27 and i1 + i2 == ivh
    return {
        "LN2H": h, "LN2M": mid, "LN2L": lo,
        "C3": float(Fraction(1, 3)), "C4": -0.25, "C5": float(Fraction(1, 5)),
        "C6": float(Fraction(-1, 6)), "C7": float(Fraction(1, 7)), "C8": -0.125,
        "IVH": ivh, "IVM": ivm, "IVL": ivl, "IVH1": i1, "IVH2": i2,
    }


def ln_acc():
    """(−1)^(k+1)/k for k = 2 ... 16 as triple-doubles: the accurate log1p series."""
    return [_fsplit(Fraction((-1) ** (k + 1), k), 3) for k in range(2, 17)]


# ── sin and cos (M4a): π/2 reduction, sin and cos of i/256, the bits of 2/π ────────────────────

for _f in ("mpfr_sin", "mpfr_cos", "mpfr_atan", "mpfr_floor"):
    getattr(_m, _f).argtypes = [_P, _P, ctypes.c_int] if _f != "mpfr_floor" else [_P, _P]
_m.mpfr_const_pi.argtypes = [_P, ctypes.c_int]
_m.mpfr_mul_2ui.argtypes = [_P, _P, ctypes.c_ulong, ctypes.c_int]
_m.mpfr_div_2ui.argtypes = [_P, _P, ctypes.c_ulong, ctypes.c_int]
_m.mpfr_ui_div.argtypes = [_P, ctypes.c_ulong, _P, ctypes.c_int]
_m.mpfr_mul_ui.argtypes = [_P, _P, ctypes.c_ulong, ctypes.c_int]

N_SC = 256         # the fast and accurate paths' table step: a = i/256
SC_ROWS = 202      # i = 0 … 201
N_TWOPI = 44       # 32-bit words of 2/π (after 64 zero bits): bits 1 … 1344
TWOPI_PREC = 1600


def _pi(prec=PREC):
    r = _V(prec)
    _m.mpfr_const_pi(r.p, _RNDN)
    return r


def _pio2():
    p = _pi()
    r = _V()
    _m.mpfr_div_2ui(r.p, p.p, 1, _RNDN)
    return r


def sc_table():
    """sin(i/256) and cos(i/256) for i = 0 … 201, each split A (26 bits) + B + C + D."""
    pio4 = _get(_pio2()) / 2
    # |r| ≤ π/4 + 2^-30 (the reduction's slack) gives |round(r·256)| ≤ 201 < SC_ROWS − 1/2.
    assert (pio4 + 2.0 ** -30) * N_SC < SC_ROWS - 0.5
    rows = []
    for i in range(SC_ROWS):
        a = _d(i / N_SC)                              # i/256 is exact
        sv = _op("mpfr_sin", a)
        cv = _op("mpfr_cos", a)
        s = _split(sv, 4, first_bits=26)
        c = _split(cv, 4, first_bits=26)
        assert _sigbits(s[0]) <= 26 and _sigbits(c[0]) <= 26, i
        assert abs(s[1]) <= 2.0 ** -26 * abs(s[0]) and abs(c[1]) <= 2.0 ** -26 * abs(c[0]), i
        if i == 0:
            assert s == [0.0, 0.0, 0.0, 0.0] and c == [1.0, 0.0, 0.0, 0.0]
        else:
            assert _rel(_exact(s[:2]), sv) < 2.0 ** -78 and _rel(_exact(c[:2]), cv) < 2.0 ** -78, i
            assert _rel(_exact(s), sv) < 2.0 ** -180 and _rel(_exact(c), cv) < 2.0 ** -180, i
        rows.append(s + c)
    return rows


def sc_consts():
    """π/2 = P1 + P2 + P3 + P4 + P5 (Cody–Waite: 33, 33, 53, 53, 53 bits), 2/π, π/2 and π as triple-doubles,
    the fast paths' coefficients and the atan2 results at ±∞."""
    pio2 = _pio2()
    parts = []
    rest = _op("mpfr_set", pio2)
    for nb in (33, 33, 53, 53, 53):
        h = _round_bits(rest, nb)
        parts.append(h)
        rest = _op("mpfr_sub", rest, _d(h))
    p1, p2, p3, p4, p5 = parts
    assert _sigbits(p1) <= 33 and _sigbits(p2) <= 33
    assert abs(p2) < 2.0 ** -32 and abs(p3) < 2.0 ** -65       # the fast path's bounds on k·P2, k·P3
    assert abs(_get(rest)) < 2.0 ** -225                       # the dropped tail, times |k| < 2^20
    pi = _pi()                                   # kept alive: a temporary's __del__ would free its limbs
    two = _V()
    _m.mpfr_ui_div(two.p, 2, pi.p, _RNDN)
    pih, pim, pil = _split(pi, 3)
    three = _V()
    _m.mpfr_mul_ui(three.p, pi.p, 3, _RNDN)
    _m.mpfr_div_2ui(three.p, three.p, 2, _RNDN)
    return {
        "P1": p1, "P2": p2, "P3": p3, "P4": p4, "P5": p5, "INVPIO2": _get(two),
        "PIH": pih, "PIM": pim, "PIL": pil, "PIO2H": pih / 2, "PIO2M": pim / 2, "PIO2L": pil / 2,
        "PIO4": pih / 4, "PI3O4": _get(three),
        "S3": float(Fraction(-1, 6)), "S5": float(Fraction(1, 120)), "S7": float(Fraction(-1, 5040)),
        "C4": float(Fraction(1, 24)), "C6": float(Fraction(-1, 720)),
        "A3": float(Fraction(-1, 3)), "A5": float(Fraction(1, 5)), "A7": float(Fraction(-1, 7)),
        "A9": float(Fraction(1, 9)),
    }


def _machin_pi(bits):
    """floor(π·2^bits), from Machin's formula π = 16·atan(1/5) − 4·atan(1/239) in integers (a check of MPFR's)."""
    guard = 64
    one = 1 << (bits + guard)

    def atan_inv(n):
        total, term, k, n2 = 0, one // n, 1, n * n
        while term:
            total += term // k if (k // 2) % 2 == 0 else -(term // k)
            term //= n2
            k += 2
        return total
    return (16 * atan_inv(5) - 4 * atan_inv(239)) >> guard


def twopi_words():
    """2/π after 64 zero bits, in 32-bit words: word w holds bits 32w − 63 … 32w − 32 of 2/π (bit j has
    weight 2^-j), most significant first. From MPFR at 1600 bits, checked against Machin's formula."""
    pi = _pi(TWOPI_PREC)
    v = _V(TWOPI_PREC)
    _m.mpfr_ui_div(v.p, 2, pi.p, _RNDN)
    words = [0, 0]
    fl = _V(TWOPI_PREC)
    for _ in range(N_TWOPI - 2):
        _m.mpfr_mul_2ui(v.p, v.p, 32, _RNDN)
        _m.mpfr_floor(fl.p, v.p)
        w = int(_get(fl))
        assert 0 <= w < 2 ** 32 and float(w) == _get(fl)
        words.append(w)
        _m.mpfr_sub(v.p, v.p, fl.p, _RNDN)
    assert ((971 + 62) >> 5) + 9 < N_TWOPI    # _km_ph_reduce reads words w … w + 9, w ≤ 32 (|x| < 2^1024)
    nbits = 32 * (N_TWOPI - 2)
    acc = 0
    for w in words[2:]:
        acc = (acc << 32) | w
    pis = _machin_pi(nbits + 80)
    want = (1 << (2 * nbits + 81)) // pis      # floor(2/π·2^nbits), within one unit of the exact value
    assert abs(acc - want) <= 1, "2/π: MPFR and Machin disagree"
    return words


def atan_table():
    """atan(j/256) for j = 0 … 256 as triple-doubles."""
    rows = []
    for j in range(N_SC + 1):
        v = _op("mpfr_atan", _d(j / N_SC))
        t = _split(v, 3)
        if j == 0:
            assert t == [0.0, 0.0, 0.0]
        else:
            assert _rel(_exact(t), v) < 2.0 ** -158, j
        rows.append(t)
    return rows


def atan_acc():
    """(−1)^k/(2k + 1) for k = 1 … 8 as triple-doubles: the accurate atan series."""
    return [_fsplit(Fraction((-1) ** k, 2 * k + 1), 3) for k in range(1, 9)]


def reduction_worst_case():
    """Lefèvre and Muller's (2003, §2.3) double closest to a nonzero multiple of π/2, 6381956970095103·2^797:
    the bound |x − kπ/2| > 2^-61 that the accurate reduction's precision rests on, recomputed here."""
    x = _V(2200)
    _m.mpfr_set_d(x.p, 6381956970095103.0, _RNDN)
    _m.mpfr_mul_2ui(x.p, x.p, 797, _RNDN)
    pi = _pi(2200)
    q = _V(2200)
    _m.mpfr_div(q.p, x.p, pi.p, _RNDN)           # x/π
    _m.mpfr_mul_2ui(q.p, q.p, 1, _RNDN)          # x/(π/2)
    fl = _V(2200)
    _m.mpfr_floor(fl.p, q.p)
    _m.mpfr_sub(q.p, q.p, fl.p, _RNDN)           # frac(x/(π/2)) in [0, 1)
    if _get(q) > 0.5:
        one = _V(2200)
        _m.mpfr_set_d(one.p, 1.0, _RNDN)
        _m.mpfr_sub(q.p, one.p, q.p, _RNDN)
    d = _get(q) * _get(_pio2())
    assert 4.6e-19 < d < 4.8e-19 and d > 2.0 ** -61, d
    assert reduction_bound() > 0.99 * d
    return d


def reduction_bound():
    """A lower bound on |x − kπ/2| over every double x ≥ 1 (k ≥ 1; below 1 the distance is at least π/2 − 1),
    binade by binade: x = m·2^(e − 52) with 2^52 ≤ m < 2^53 and |x − kπ/2| = (π/2)·|m·α − k|, α = 2^(e−52)·2/π.
    If q_n ≤ 2^53 − 1 < q_(n+1) are consecutive continued-fraction denominators of α, every m < q_(n+1) has
    |m·α − k| ≥ |q_n·α − p_n| (best approximations of the second kind), so (π/2)·|q_n·α − p_n| bounds the binade.
    α is taken from the 1,344 bits of 2/π (an error below 2^-320 after multiplying by m)."""
    t = 0
    for w in twopi_words()[2:]:
        t = (t << 32) | w
    nb = 32 * (N_TWOPI - 2)
    best = None
    for e in range(0, 1024):
        num, den = t, 1 << (nb - (e - 52))        # α = num/den
        a, b = num, den
        p0, q0, p1, q1 = 0, 1, 1, 0
        while b:
            c = a // b
            a, b = b, a - c * b
            p2, q2 = c * p1 + p0, c * q1 + q0
            if q2 >= 1 << 53:
                break
            p0, q0, p1, q1 = p1, q1, p2, q2
        dist = abs(Fraction(q1 * num, den) - p1)    # |q_n·α − p_n|, q_n the last denominator below 2^53
        bound = float(dist) * 1.5707963267948966
        if best is None or bound < best[0]:
            best = (bound, e)
    assert best[0] > 2.0 ** -61 and best[1] == 849, best   # the binade of 6381956970095103·2^797
    return best[0]


# ── 2^(1/6) (M4d): ADR 0002 #16's lj_equilibrium_distance ───────────────────────────────────────

ROOT6_BITS = 450


def _iroot(n, k):
    """⌊n^(1/k)⌋ for an integer n ≥ 0 (Newton's method in integers)."""
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def root6_2_bounds():
    """2^(1/6) between two exact rationals 2^-450 apart: ⌊2^(1/6)·2^450⌋ = ⌊(2^2701)^(1/6)⌋, an integer root."""
    c = _iroot(1 << (6 * ROOT6_BITS + 1), 6)
    assert c ** 6 <= 1 << (6 * ROOT6_BITS + 1) < (c + 1) ** 6
    return Fraction(c, 1 << ROOT6_BITS), Fraction(c + 1, 1 << ROOT6_BITS)


def root6_2_min_distance():
    """min over 1 ≤ q ≤ 2^54 of |q·2^(1/6) − p| (p the nearest integer): |q_N·c − p_N| for the last continued-fraction
    denominator q_N ≤ 2^54 of c = 2^(1/6) (best approximations of the second kind), on the two rational bounds of c,
    whose expansions agree far past q_N. Every midpoint test of x·2^(1/6) (port: _km_cr_sixth_root2_mul) reduces to it:
    with x = S·2^e (S < 2^53) the result in ulps is S·c (or S·c/2 in the binade [2, 4)), and its distance δ to a
    midpoint k + 1/2 gives ||q·c|| ≤ 2δ for q = 2S (or S) ≤ 2^54, so δ ≥ (this)/2."""
    lo, hi = root6_2_bounds()
    out = []
    for v in (lo, hi):
        a, b = v.numerator, v.denominator
        p0, q0, p1, q1 = 0, 1, 1, 0
        cf = []
        while b:
            t = a // b
            a, b = b, a - t * b
            p2, q2 = t * p1 + p0, t * q1 + q0
            if q2 > 1 << 54:
                break
            cf.append(t)
            p0, q0, p1, q1 = p1, q1, p2, q2
        out.append((cf, p1, q1))
    assert out[0] == out[1], "2^(1/6)'s continued fraction is not settled at 450 bits"
    _, p, q = out[0]
    d = min(abs(q * lo - p), abs(q * hi - p))
    assert d > 0
    return d, q


def root6_2():
    """2^(1/6) = H + M + L (each the nearest double to what is left), with the properties _km_cr_sixth_root2_mul's proof
    uses: H is Rust's RN(2^(1/6)) = 2.0.powf(1.0 / 6.0); |H + M + L − 2^(1/6)| < 2^-158; and no x·2^(1/6) is within
    2^-100 ulp of a rounding boundary (the continued-fraction bound above), so the triple-double product, accurate to
    2^-102.6 ulp, rounds correctly by _km_cr_finish with no rounding test."""
    lo, hi = root6_2_bounds()
    h, m, l_ = _fsplit(lo, 3)
    assert bits(h) == 0x3FF1F59AC3C7D6C0, hex(bits(h))
    assert abs(_exact([h, m, l_]) - lo) < Fraction(1, 1 << 158)
    d, q = root6_2_min_distance()
    assert d / 2 > Fraction(1, 1 << 100), float(d)
    return {"H": h, "M": m, "L": l_}, d, q


# ── 1/100 (M4e): ADR 0002 #19's wavenumber conversions ──────────────────────────────────────────

def hundredth():
    """1/100 = H + L + t as a double-double (each the nearest double to what is left), with the properties
    _km_cr_hundredth_div's fast path uses: H is Rust's 0.01 literal, |t| < 2^-116 (so the residual 1/100 − y0·x,
    formed from H, L and an exact product, is within 2^-110.2 of the exact one), and 1/(100·x) is never a midpoint
    (100·x = 2^2·25·m·2^e has the odd factor 25, so no dyadic value equals 1/(100·x)), which is why the integer
    path's sticky bit is never 0."""
    q = Fraction(1, 100)
    h, l_ = _fsplit(q, 2)
    assert h == 0.01 and bits(h) == 0x3F847AE147AE147B, hex(bits(h))
    assert abs(q - _exact([h, l_])) < Fraction(1, 1 << 116)
    assert 0 < abs(l_) < 2.0 ** -61
    return {"H": h, "L": l_}


# ── thresholds, from the correctly rounded MPFR functions ───────────────────────────────────────

def _last_true(pred, lo_b, hi_b):
    """The largest bit pattern b in [lo_b, hi_b] with pred(b), for pred monotone (true then false)."""
    assert pred(lo_b) and not pred(hi_b)
    while hi_b - lo_b > 1:
        mid = (lo_b + hi_b) // 2
        if pred(mid):
            lo_b = mid
        else:
            hi_b = mid
    return lo_b


def limits():
    fb = _m_of_bits
    inf = float("inf")
    # Positive arguments: bit patterns increase with x.
    exp_ovf = _last_true(lambda b: mpfr_ref.exp(fb(b)) < inf, bits(700.0), bits(720.0))
    exp10_ovf = _last_true(lambda b: mpfr_ref.exp10(fb(b)) < inf, bits(300.0), bits(320.0))
    expm1_ovf = _last_true(lambda b: mpfr_ref.expm1(fb(b)) < inf, bits(700.0), bits(720.0))
    # Negative arguments: bit patterns increase with |x|.
    exp_unf = _last_true(lambda b: mpfr_ref.exp(fb(b)) > 0.0, bits(-740.0), bits(-750.0))
    exp10_unf = _last_true(lambda b: mpfr_ref.exp10(fb(b)) > 0.0, bits(-320.0), bits(-330.0))
    expm1_m1 = _last_true(lambda b: mpfr_ref.expm1(fb(b)) > -1.0, bits(-30.0), bits(-40.0))
    return {
        "EXP_OVF": exp_ovf, "EXP_UNF": exp_unf, "EXP10_OVF": exp10_ovf, "EXP10_UNF": exp10_unf,
        "EXPM1_OVF": expm1_ovf, "EXPM1_M1": expm1_m1,
    }


# ── output ───────────────────────────────────────────────────────────────────────────────────────

def _hex(x):
    return f"0x{bits(x):016X}"


def _table(name, words, comment, per_line=4):
    out = [f"# {line}" if line else "#" for line in comment]
    out.append(f"var {name}: i64[{len(words)}] = {{")
    for i in range(0, len(words), per_line):
        chunk = words[i:i + per_line]
        sep = "," if i + per_line < len(words) else ""
        out.append("    " + ", ".join(chunk) + sep)
    out.append("};")
    out.append("")
    return out


def _scalars(prefix, d, order, comment):
    """One i64 global per value (`var _KM_CR_<prefix>_<NAME> = 0x…;`): a scalar global reads in one load,
    a table element in an address computation and a load, which the fast paths feel."""
    lines = [f"# {line}" if line else "#" for line in comment]
    for k in order:
        lines.append(f"var _BJ_K_CR_{prefix}{k} = {_hex(d[k])};   # {d[k]!r}")
    lines.append("")
    return lines


EXP_K_ORDER = ["INVLN2N", "L1", "L2", "L3", "L4", "P3", "P4", "P5",
               "LN10H", "LN10M", "LN10L", "LN10H1", "LN10H2"]
EXPM1_ORDER = ["Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"]
LN_K_ORDER = ["LN2H", "LN2M", "LN2L", "C3", "C4", "C5", "C6", "C7", "C8",
              "IVH", "IVM", "IVL", "IVH1", "IVH2"]
LIM_ORDER = ["EXP_OVF", "EXP_UNF", "EXP10_OVF", "EXP10_UNF", "EXPM1_OVF", "EXPM1_M1"]
SC_K_ORDER = ["P1", "P2", "P3", "P4", "P5", "INVPIO2", "PIH", "PIM", "PIL", "PIO2H", "PIO2M", "PIO2L",
              "PIO4", "PI3O4", "S3", "S5", "S7", "C4", "C6", "A3", "A5", "A7", "A9"]


def render():
    """src/crmath_k_data.cyr: the subset of kimiya's tables that bijli's acos, asin, log10 and expm1 read (the
    exp, ln and atan tables, their constants and the limits), under bijli's `_BJ_K_CR_` names."""
    et = exp_table()
    ek = exp_consts()
    lt = ln_table()
    lk = ln_consts()
    lim = limits()
    L = []
    L += [
        "# bijli — tables of the correctly rounded acos, asin, log10 and expm1 (src/crmath_k.cyr, src/crmath_kx.cyr)",
        "#",
        f"# GENERATED by scripts/oracle/crtables_k.py (kimiya's scripts/oracle/mathx_tables.py, its exp, ln and",
        f"# atan tables) from MPFR {mpfr_ref.VERSION} (400-bit values) and exact rationals. Do not edit.",
        "# Regenerate: python3 -I -B scripts/oracle/crtables_k.py   (--check verifies this file is current)",
        "#",
        "# Every value is an IEEE-754 binary64 bit pattern, in an i64 table (read with load64) or an i64 global,",
        "# assigned to an f64 local (a hex literal in an f64 initialiser would be an integer). A multi-word entry",
        "# is a value split into doubles, each the nearest double to what the previous ones leave: a + b + c is",
        "# the value to about 2^(-53·k) relative. kimiya ADR 0006; the error analysis is in src/crmath_k.cyr.",
        "",
    ]
    L += _table("_BJ_K_CR_EXP_T", [_hex(v) for row in et for v in row], [
        "2^(j/512) for j = 0 … 511, four words per j at index 4j: A_j (rounded to 26 significant bits,",
        "so A_j·r1 is exact for a 27-bit r1), then B_j, C_j and D_j (each the nearest double to the rest).",
        "The fast path reads A_j + B_j (2^(j/512) to 2^-79), the accurate path all four (to 2^-180).",
        "Row 0 is (1, 0, 0, 0).",
    ], per_line=4)
    L += _scalars("EXP_", ek, EXP_K_ORDER, [
        "exp / exp10 reduction: x = (512·E + j)·ln2/512 + r. INVLN2N = 512/ln 2 rounded; ln2/512 =",
        "L1 + L2 + L3 + L4 (L1 at 33 bits, so K·L1 is exact for |K| < 2^20). P3 … P5: 1/6, 1/24, 1/120,",
        "the fast path's Taylor coefficients (1/2 is exact). ln 10 = LN10H + LN10M + LN10L, and",
        "LN10H = LN10H1 + LN10H2 (26 + 27 bits, Dekker's product of x·LN10H without a split of LN10H).",
    ])
    L += _scalars("EXPM1_", dict(zip(EXPM1_ORDER, expm1_poly())), EXPM1_ORDER, [
        "Q3 … Q9 = 1/3!, …, 1/9!: expm1's small-|x| fast path, x + x²/2 + x³·(Q3 + x·Q4 + … + x⁶·Q9).",
    ])
    L += _table("_BJ_K_CR_EXP_ACC", [_hex(v) for row in exp_acc() for v in row], [
        "1/k! for k = 2 … 17 as triple-doubles (index 3·(k − 2)): the accurate paths' Taylor series of",
        "exp(r) − 1 and expm1(x).",
    ], per_line=3)
    L += _table("_BJ_K_CR_LN_T", [_hex(v) for row in lt for v in row], [
        "ln: x = 2^E·m, m in [0.74951171875, 1.4990234375), j = the 9 bits of bits(x) − 0x3FE7FC0000000000",
        "below the exponent field. Four words per j at index 4j: c_j ≈ 1/m on interval j (a multiple of",
        "2^-9 if c_j > 1, of 2^-10 if c_j < 1, so r = m·c_j − 1 is exact and |r| < 2^-9), then",
        "−ln(c_j) as a triple-double. Interval 256 holds m = 1: c = 1, −ln c = 0.",
    ], per_line=4)
    L += _scalars("LN_", lk, LN_K_ORDER, [
        "ln 2 = LN2H + LN2M + LN2L (LN2H at 42 bits, so E·LN2H is exact for |E| < 2^11). C3 … C8: the",
        "Taylor coefficients 1/3, −1/4, 1/5, −1/6, 1/7, −1/8 of log1p(r) − r + r²/2. 1/ln 10 = IVH + IVM +",
        "IVL, and IVH = IVH1 + IVH2 (26 + 27 bits).",
    ])
    L += _table("_BJ_K_CR_LN_ACC", [_hex(v) for row in ln_acc() for v in row], [
        "(−1)^(k+1)/k for k = 2 … 16 as triple-doubles (index 3·(k − 2)): the accurate log1p series.",
    ], per_line=3)
    lim_d = {k: _m_of_bits(v) for k, v in lim.items()}
    L += _scalars("", lim_d, LIM_ORDER, [
        "Thresholds, each the last double before the result changes (bisection over bit patterns with",
        "MPFR's correctly rounded function): EXP_OVF the largest x with exp(x) finite, EXP_UNF the",
        "most negative x with exp(x) > 0 (below it exp(x) < 2^-1075 rounds to +0), the same for exp10,",
        "EXPM1_OVF the largest x with expm1(x) finite, EXPM1_M1 the most negative x with expm1(x) > −1.",
    ])
    L += _scalars("SC_", sc_consts(), SC_K_ORDER, [
        "π = PIH + PIM + PIL and π/2 = PIO2H + PIO2M + PIO2L (triple-doubles); A3 … A9 (−1/3, 1/5, −1/7, 1/9):",
        "the fast atan polynomial. (The other SC_ values are kimiya's sin and cos constants, kept so the",
        "table set stays kimiya's; bijli's sin and cos are src/crmath.cyr's.)",
    ])
    L += _table("_BJ_K_CR_ATAN_T", [_hex(v) for row in atan_table() for v in row], [
        "atan(j/256) for j = 0 … 256 as triple-doubles (index 3j): the atan table of acos and asin.",
    ], per_line=3)
    L += _table("_BJ_K_CR_ATAN_ACC", [_hex(v) for row in atan_acc() for v in row], [
        "(−1)^k/(2k + 1) for k = 1 … 8 as triple-doubles (index 3·(k − 1)): the accurate atan series.",
    ], per_line=3)
    return "\n".join(L).rstrip("\n") + "\n"


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "src/crmath_k_data.cyr"
    text = render()
    if "--check" in sys.argv[1:]:
        old = out.read_text() if out.exists() else ""
        if old != text:
            print("STALE    src/crmath_k_data.cyr (run python3 -I -B scripts/oracle/crtables_k.py)")
            sys.exit(1)
        print(f"crtables_k: src/crmath_k_data.cyr current (MPFR {mpfr_ref.VERSION})")
        return
    out.write_text(text)
    print(f"crtables_k: wrote src/crmath_k_data.cyr ({text.count(chr(10))} lines, MPFR {mpfr_ref.VERSION})")


if __name__ == "__main__":
    main()
