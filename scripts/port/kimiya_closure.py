"""Copies the call closure of kimiya's correctly rounded rationals and unbounded-exponent helpers into bijli.

Usage, from the repo root:
    python3 -I -B scripts/port/kimiya_closure.py ~/Repos/kimiya [--check]

Reads kimiya's src/mathx.cyr, takes every function, const and var that ROOTS reach (transitively), renames the
prefixes (_km_ → _bj_k_, _KM_ → _BJ_K_) as src/crmath_k.cyr's copy did, leaves out what src/crmath_k.cyr already
defines (checked identical after the rename), and writes src/crmath_kx.cyr in kimiya's order, each item with the
comment lines above it (section headers kept). Nothing is edited by hand: re-run it after kimiya changes. With
--check it exits 1 when the file differs.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOTS = [
    "_km_cr_mono", "_km_cr_mul_div", "_km_cr_mul3", "_km_cr_div2", "_km_cr_mul3_div", "_km_cr_mul5",
    "_km_cr_mul2_div2", "_km_cr_mul3_div2", "_km_cr_mul4", "_km_cr_div3", "_km_cr_sqrt_div",
    "_km_cr_mul_sqrt_ratio", "_km_cr_ratio",
    "_km_ux_from", "_km_ux_neg", "_km_ux_add", "_km_ux_sub", "_km_ux_mul", "_km_ux_div", "_km_ux_sqrt",
    "_km_ux_cmp", "_km_ux_le", "_km_ux_max", "_km_ux_is_positive", "_km_ux_normal", "_km_ux_round",
    "_km_ux_add_f64", "_km_ux_sub_f64", "_km_ux_mul_f64", "_km_ux_div_f64", "_km_ux_sqrt_f64",
    # Not reached from the above but tested by kimiya's groups that tests/crmath_kx.tcyr copies.
    "_km_bn_add", "_km_bn_sub", "_km_bn_mul_small", "_km_dd_div_cr", "_km_dd_mul", "_km_max",
]

# kimiya's tests/mathx.tcyr from this header to its main(): the ux, multi-word, ratio, monomial, B2 and B3 groups,
# copied to tests/crmath_kx.tcyr with the same renames. Its ux atan2 checks are dropped (bijli doesn't copy kimiya's
# atan2: src/crmath.cyr's correctly rounded atan2 is scale-free and no bijli formula needs an unbounded one); each
# drop is an exact text replacement below, refused when kimiya's text changes.
TESTS_FROM = "# ── binary64 with an unbounded exponent (M4e-2): ADR 0002 #20"
TESTS_TO = "\nfn main() {"
TEST_EDITS = [
    ("    var f: f64 = _km_ux_atan2_f64(s1, e1, s2, e2);\n    return f;",
     "    var f: f64 = 0.0 / 0.0;                        # bijli: op 4 (atan2) is not copied\n    return f;"),
    ("""        var r: f64 = uxf(load64(at), load64(at + 8), load64(at + 16) - 8192, load64(at + 24), load64(at + 32) - 8192);
        assert_bits_eq(r, load64(at + 40), "an unbounded-exponent final operation against MPFR");
        if (load64(at) == 4) {
            var s1: f64 = from_bits(load64(at + 8));
            var s2: f64 = from_bits(load64(at + 24));
            var a: f64 = _km_ux_atan2_acc(s1, load64(at + 16) - 8192, s2, load64(at + 32) - 8192);
            assert_bits_eq(a, load64(at + 40), "ux atan2, accurate path");
        }""",
     """        if (load64(at) != 4) {                          # bijli: the atan2 rows (op 4) are not copied
            var r: f64 = uxf(load64(at), load64(at + 8), load64(at + 16) - 8192, load64(at + 24), load64(at + 32) - 8192);
            assert_bits_eq(r, load64(at + 40), "an unbounded-exponent final operation against MPFR");
        }"""),
    ("""    var f8: f64 = _km_ux_atan2_f64(inf, 0, inf, 0);
    assert_bits_eq(f8, 0x3FE921FB54442D18, "atan2(+∞, +∞) = π/4");
    var f9: f64 = _km_ux_atan2_f64(nz, 0, m11, 9000);
    assert_bits_eq(f9, 0xC00921FB54442D18, "atan2(−0, −1.1·2^9000) = −π");
""", ""),
    ("    var bad_a = 0;\n", ""),
    ("""        # atan2: the correctly rounded atan2 of the doubles, and of the pair shifted by any common amount
        var t0: f64 = _km_cr_atan2(y, x);
        var t1: f64 = _km_ux_atan2_f64(ys, ye, xs, xe);
        var t2: f64 = _km_ux_atan2_f64(ys, ye + sh, xs, xe + sh);
        var t3: f64 = _km_ux_atan2_acc(ys, ye + sh, xs, xe + sh);
        if (ulp_bits_match(_km_bits(t1), _km_bits(t0)) != 1) { bad_a = bad_a + 1; }
        if (ulp_bits_match(_km_bits(t2), _km_bits(t0)) != 1) { bad_a = bad_a + 1; }
        if (ulp_bits_match(_km_bits(t3), _km_bits(t0)) != 1) { bad_a = bad_a + 1; }
""", ""),
    ("""    assert_eq(bad_a, 0, "ux atan2 = _km_cr_atan2 at any common shift (fast and accurate paths)");
""", ""),
]

NAME = re.compile(r"\b(_km_\w+|_KM_\w+)\b")


def items(text: str):
    """[(name, kind, first_line, last_line)] of the top-level fn / const / var definitions."""
    lines = text.split("\n")
    out, k = [], 0
    while k < len(lines):
        line = lines[k]
        m = re.match(r"(?:public )?(fn|const|var) (\w+)", line)
        if m:
            start = k
            if m.group(1) == "fn":
                while not lines[k].startswith("}"):
                    k += 1
            else:
                depth = line.count("[") - line.count("]") + line.count("(") - line.count(")")
                while depth > 0 or not lines[k].rstrip().endswith(";"):
                    k += 1
                    depth += lines[k].count("[") - lines[k].count("]") + lines[k].count("(") - lines[k].count(")")
            out.append((m.group(2), m.group(1), start, k))
        k += 1
    return out, lines


def rename(s: str) -> str:
    return s.replace("_km_", "_bj_k_").replace("_KM_", "_BJ_K_")


WIDTH = 120


def rewrap(text: str) -> str:
    """Lines the rename pushed past WIDTH characters, rewrapped: a comment's overflow words move to the start of the
    next line of its paragraph (a new line when the paragraph ends there), cascading; a code line breaks at its last
    ", " that fits, the rest indented two more (cyrius fmt's continuation)."""
    lines = text.split("\n")
    k = 0
    while k < len(lines):
        line = lines[k]
        if len(line) <= WIDTH:
            k += 1
            continue
        stripped = line.lstrip()
        ind = line[:len(line) - len(stripped)]
        if stripped.startswith("#"):
            cut = line.rfind(" ", 0, WIDTH + 1)
            if cut <= len(ind) + 1:
                k += 1
                continue
            head, rest = line[:cut], line[cut + 1:]
            lines[k] = head
            nxt = lines[k + 1] if k + 1 < len(lines) else ""
            cont = nxt.startswith(ind + "# ") and not nxt.startswith(ind + "#  ") and nxt[len(ind) + 2:][:2] not in ("- ", "──") \
                and nxt.strip() != "#"
            if cont:
                lines[k + 1] = ind + "# " + rest + " " + nxt[len(ind) + 2:]
            else:
                lines.insert(k + 1, ind + "# " + rest)
            k += 1
            continue
        cut, sep = line.rfind(", ", 0, WIDTH - 1), ", "
        if cut <= len(ind):
            cut, sep = line.rfind(" + ", 0, WIDTH - 2), " + "
        if cut <= len(ind):
            k += 1                                       # no break point: left long (lint reports it)
            continue
        depth = max(1, line[:cut].count("(") - line[:cut].count(")"))
        lines[k] = line[:cut + 1] if sep == ", " else line[:cut + 2]
        lines.insert(k + 1, ind + "  " * depth + line[cut + len(sep):])
        k += 1
    return "\n".join(lines)


def fmt(root: Path, text: str, ext: str) -> str:
    """text as `cyrius fmt` writes it (its continuation indents), through a scratch file under build/."""
    tmp = root / "build" / f"kimiya_closure_fmt{ext}"
    tmp.parent.mkdir(exist_ok=True)
    tmp.write_text(text)
    subprocess.run(["cyrius", "fmt", str(tmp)], cwd=root, capture_output=True, check=True)
    out = tmp.read_text()
    tmp.unlink()
    return out


def main() -> int:
    kimiya = Path(sys.argv[1]).expanduser()
    root = Path(".").resolve()
    text = (kimiya / "src/mathx.cyr").read_text()
    its, lines = items(text)
    by_name = {n: i for i, (n, *_r) in enumerate(its)}
    body = {n: "\n".join(lines[a:b + 1]) for n, _k, a, b in its}
    have_text = (root / "src/crmath_k.cyr").read_text()
    have_items, have_lines = items(have_text)
    have = {n: "\n".join(have_lines[a:b + 1]) for n, _k, a, b in have_items}

    want, todo = set(), list(ROOTS)
    while todo:
        n = todo.pop()
        if n in want:
            continue
        if n not in by_name:
            print(f"not defined in kimiya: {n}")
            return 1
        want.add(n)
        for ref in NAME.findall(body[n]):
            if ref != n and ref in by_name:
                todo.append(ref)

    problems, chunks, prev_end, carried = [], [], -1, []
    for n, kind, a, b in its:
        gap = lines[prev_end + 1:a]
        prev_end = b
        # The comment lines directly above the item (back to a blank line), and any section headers in the gap.
        own = []
        j = len(gap) - 1
        while j >= 0 and (gap[j].startswith("#") and not gap[j].startswith("#include")):
            own.insert(0, gap[j])
            j -= 1
        headers = [g for g in gap[:j + 1] if g.startswith("# ──")]
        if n not in want:
            carried += headers
            continue
        renamed = rename(n)
        if renamed in have:
            if have[renamed] != rename(body[n]):
                problems.append(f"{renamed}: src/crmath_k.cyr's copy differs from kimiya's")
            carried += headers
            continue
        block = carried + headers
        carried = []
        chunks.append("\n".join(([""] + block if block else []) + [""] + own + lines[a:b + 1]))

    head = (
        "# crmath_kx.cyr — kimiya's correctly rounded monomials, square roots of quotients and binary64 with an\n"
        "# unbounded exponent (kimiya ADR 0002 #20, #27–#32), for bijli's range-edge Part (M1 part 5).\n"
        "#\n"
        "# GENERATED by scripts/port/kimiya_closure.py from kimiya's src/mathx.cyr: the call closure of its _km_cr_mono\n"
        "# family, _km_cr_sqrt_div, _km_cr_mul_sqrt_ratio and _km_ux_*, function for function, with the prefixes renamed\n"
        "# (_km_ → _bj_k_, _KM_ → _BJ_K_) as src/crmath_k.cyr's copy was; what src/crmath_k.cyr already holds is used from\n"
        "# there. Do not edit: re-run the script (port plan §10 G1: a shared module would replace both copies).\n"
    )
    out = head + rewrap(rename("\n".join(chunks))) + "\n"
    tests = (kimiya / "tests/mathx.tcyr").read_text()
    a = tests.index(TESTS_FROM)
    b = tests.index(TESTS_TO, a)
    section = tests[a:b].rstrip()
    while section.split("\n")[-1].startswith("#"):        # main()'s own comment
        section = section.rsplit("\n", 1)[0].rstrip()
    section += "\n\n"
    for old, new in TEST_EDITS:
        if section.count(old) != 1:
            problems.append(f"kimiya's tests changed: a dropped atan2 check is not where it was: {old[:60]!r}")
            continue
        section = section.replace(old, new)
    groups = re.findall(r"^fn (_t_\w+)\(\)", section, re.M)
    tout = (
        "# crmath_kx.tcyr — kimiya's tests of src/crmath_kx.cyr (bijli M1 part 5).\n"
        "#\n"
        "# GENERATED by scripts/port/kimiya_closure.py from kimiya's tests/mathx.tcyr: its ux, multi-word integer, RN(p/q·2^e),\n"
        "# monomial, B2 and B3 groups with the same renames, their ux atan2 checks dropped (bijli copies no unbounded atan2).\n"
        "# Every expected bit pattern is kimiya's (MPFR 4.2.2, or Python's Fraction). Do not edit: re-run the script.\n"
        "\n"
        "include \"tests/lib.cyr\"\n"
        "\n"
        "include \"tests/ulp.cyr\"\n"
        "\n"
        "# The f64 with bit pattern b.\n"
        "#must_use\n"
        "fn from_bits(b): f64 {\n"
        "    var x: f64 = 0.0;\n"
        "    store64(&x, b);\n"
        "    return x;\n"
        "}\n"
        "\n" + rewrap(rename(section)) +
        "# Runs every copied group; returns the failure count from assert_summary().\n"
        "fn main() {\n" + "".join(f"    {g}();\n" for g in groups) +
        "    return assert_summary();\n"
        "}\n"
        "\n"
        "var exit_code = main();\n"
        "if (exit_code > 0) { exit_code = 1; }\n"
        "sys_exit_group(exit_code);\n"
    )
    outs = {root / "src/crmath_kx.cyr": fmt(root, out, ".cyr"), root / "tests/crmath_kx.tcyr": fmt(root, tout, ".tcyr")}
    for p in problems:
        print("  " + p)
    if "--check" in sys.argv[2:]:
        stale = [str(t.relative_to(root)) for t, text in outs.items() if not t.exists() or t.read_text() != text]
        print(f"kimiya closure: {len(chunks)} items, {len(groups)} test groups; {'stale: ' + ', '.join(stale) if stale else 'current'}")
        return 0 if not stale and not problems else 1
    for t, text in outs.items():
        t.write_text(text)
    print(f"wrote src/crmath_kx.cyr ({len(chunks)} items, {out.count(chr(10))} lines) and tests/crmath_kx.tcyr "
          f"({len(groups)} groups); {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
