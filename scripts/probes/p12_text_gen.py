"""Probe 6 inputs (port plan §6 step 9): Touchstone-style number text. Writes build/p12_in.txt, one decimal per
line, in every shape Touchstone files and Rust's formatting produce: Python repr (shortest round-trip, the same
digits Rust's {} / {:e} print), %.6E / %.9E / %.15E / %g, integers, signed zeros, exponent-only forms. Python's
float() is correctly rounded, as Rust's str::parse::<f64> is, so its bits are the reference (p12_check.py)."""
import random, struct, sys
rng = random.Random(20261111)
out = []
def rnd_double():
    while True:
        x = struct.unpack("<d", struct.pack("<Q", rng.getrandbits(64)))[0]
        if x == x and abs(x) != float("inf"):
            return x
for _ in range(20000):
    x = rnd_double()
    out += [repr(x), "%.6E" % x, "%.15e" % x]
for _ in range(20000):
    x = rng.uniform(-1e3, 1e3) * 10 ** rng.randint(-12, 12)
    out += [repr(x), "%.9E" % x, "%g" % x]
out += ["0", "-0", "0.0", "-0.0", "1", "-1", "1e9", "1E+09", "1.0E9", "  2.5", "5e-324", "1e-400", "1.7976931348623157e308",
        "123456789012345678901234567890", "0.000000000000000000000000000001", "-1.000000000E+00", "4.9e-324", "2.2250738585072014e-308"]
open(sys.argv[1], "w").write("\n".join(out) + "\n")
print(len(out), "lines")
