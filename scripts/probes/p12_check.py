"""Probe 6 check: compare p12_text.cyr's output with Python's float() (parse) and repr() (digits)."""
import re, struct, sys
ins = open(sys.argv[1]).read().split("\n")[:-1]
outs = open(sys.argv[2]).read().split("\n")[:-1]
assert len(ins) == len(outs), (len(ins), len(outs))
def bits(x): return struct.unpack("<Q", struct.pack("<d", x))[0]
def sig(s):
    """(significant digits, decimal exponent) of a finite decimal string."""
    s = s.strip().lstrip("+-").lower()
    m, _, e = s.partition("e")
    e = int(e) if e else 0
    ip, _, fp = m.partition(".")
    digits = (ip + fp).lstrip("0")
    lead = len(ip.lstrip("0")) if ip.strip("0") else -(len(fp) - len(fp.lstrip("0")))
    return digits.rstrip("0"), e + lead
bad_parse = bad_round = longer = shorter = 0
for s, o in zip(ins, outs):
    hexpart, _, js = o.partition(" ")
    got = int(hexpart, 16)
    want = bits(float(s))
    if got != want:
        bad_parse += 1
        if bad_parse <= 5: print("PARSE", repr(s), hex(got), hex(want))
    x = struct.unpack("<d", struct.pack("<Q", got))[0]
    if js == "null":
        continue
    if bits(float(js)) != got:
        bad_round += 1
        if bad_round <= 5: print("ROUNDTRIP", repr(s), js)
    d_js, _ = sig(js)
    d_py, _ = sig(repr(x))
    if len(d_js) > len(d_py): longer += 1
    if len(d_js) < len(d_py): shorter += 1
print(f"{len(ins)} lines | parse mismatches vs float(): {bad_parse} | bayan round-trip failures: {bad_round} | "
      f"bayan longer than repr: {longer}, shorter: {shorter}")
