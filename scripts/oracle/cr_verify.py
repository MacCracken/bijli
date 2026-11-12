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
    elif name == "log10":
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
    xs += [rand_double(rng) for _ in range(count - len(xs))]
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
    "tan": (1, "tan"), "tan_slow": (1, "tan"),
}


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


def main() -> int:
    (root / "build").mkdir(exist_ok=True)
    names = wanted or list(FUNCS)
    rng = random.Random(0xB1_1E_C0DE)
    calls = []
    for name in names:
        arity, mname = FUNCS[name]
        argsets = unary_args(name, rng, n) if arity == 1 else binary_args(name, rng, n)
        calls += [(name, mname, a) for a in argsets]
    lines = [" ".join([name] + [f"{f2b(x):016x}" for x in a]) for name, _, a in calls]
    (root / "build/cr_in.txt").write_text("\n".join(lines) + "\n")
    targets = ["x86_64"] + (["aarch64"] if aarch64 else [])
    results = {t: run(t, lines) for t in targets}
    bad = 0
    per = {}
    for i, (name, mname, a) in enumerate(calls):
        want = f2b(mpfr(mname, a))
        for t in targets:
            g = results[t][i]
            ok = g == want or (b2f(g) != b2f(g) and b2f(want) != b2f(want))
            per.setdefault((name, t), [0, 0])[0 if ok else 1] += 1
            if not ok:
                bad += 1
                if bad <= 20:
                    print(f"MISMATCH {t} {name}{tuple(a)}: got {g:016x}, MPFR {want:016x}")
    for (name, t), (ok, ko) in sorted(per.items()):
        print(f"  {name:6} {t:7} {ok + ko:8} args: {ok} equal to MPFR, {ko} differ")
    print(f"cr_verify: {len(calls)} calls x {len(targets)} targets, {bad} mismatches (MPFR {mpfr_ref.VERSION})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
