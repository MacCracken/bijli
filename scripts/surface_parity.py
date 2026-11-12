"""Surface parity: every public item of the ported Rust modules has a Cyrius counterpart.

Port plan §8 (a per-milestone gate) and M8 (the release gate), brought forward to M0 (§6 step 7). Adapted from
impetus's scripts/surface_parity.py (itself kimiya's and tanmatra's). Usage, from the repo root:
    python3 -I -B scripts/surface_parity.py .

For each module in MODULES it lists rust-old/src/<module>.rs's public items outside every `#[cfg(test)]` module
(scripts/rustsrc.py skips them by brace depth, wherever they sit) and checks that src/ defines the Cyrius name the
naming rules give (port plan §9.4):
- free `pub fn f`                         -> `fn bj_f(`
- `impl Struct { pub fn f(...) }`, no self -> `fn BjStruct_f(`
- `impl Struct { pub fn m(&self ...) }`    -> `fn m(` inside `impl BjStruct {`
- `impl Default for X`                     -> `fn BjX_default(`
- `impl Add / Sub / Mul / Neg / Div for X` -> the operator fn `BjX_add` / `_sub` / `_mul` / `_neg` / `_div`
- `pub const X`                            -> `const BJ_X =` (src/constants.cyr, generated)
- `pub enum E` and each variant V          -> `enum <ENUM_NAMES[E]> {` and `<ENUM_PREFIX[E]><SNAKE(V)> =`
- `pub struct S` and each `pub` field f    -> `struct BjS {` with a field f, or FIELD_RENAMES' fields
A function the rules would name otherwise is in FN_RENAMES; anything else needs an EXCEPTIONS entry with its reason
(an ADR 0002 change, say). Every file under rust-old/src must be in MODULES, PENDING_MODULES (a later milestone) or
NOT_APPLICABLE_MODULES.

It also checks two things the port relies on:
- tests/lib.cyr includes `[lib] modules` in cyrius.cyml, in order, and there is no `[test] modules`: tests see
  exactly the library, each src/ file as a file of its own (port-parity.md §8; impetus's tests/lib_*.cyr pattern);
- no `s.f = g(.., t.e, ..)` where g returns a struct and t.e is struct-typed (impetus's scan of the open cyrius
  field-store defect).
Exits 1 on a missing item, an EXCEPTIONS entry that names no Rust item, a file unaccounted for, or a failed check.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rustsrc  # noqa: E402

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()

# Ported modules (their milestone done or in progress, item by item).
MODULES = ["error", "logging"]
# Rust module -> the milestone that ports it (port plan §4.1).
PENDING_MODULES = {
    "maxwell": "M1", "charge": "M1",
    "material": "M2", "wave": "M2", "circuit": "M2", "relativity": "M2",
    "polarization": "M3", "scattering": "M3", "beam": "M3",
    "rf": "M4",
    "fdtd": "M5–M6", "compute": "M6",
    "integration/soorat": "M7", "bridge": "M7", "ai": "M7",
}
# Modules being ported Part by Part: module -> (the Parts that finish it, the Rust items still to come, by name). Their
# other items are checked as a ported module's now; the listed ones count as pending, and one found ported is a
# problem (drop it from the list).
IN_PROGRESS = {
    "field": ("M1 part 3c", {"trace_field_line"}),
}
NOT_APPLICABLE_MODULES = {
    "lib": "crate root: module list and re-exports",
    "integration/mod": "module root of integration/ (no items)",
}

ENUM_NAMES = {"BijliError": "BjError"}
ENUM_PREFIX = {"BijliError": "BJ_ERR_"}
OPERATORS = {"Add": "add", "Sub": "sub", "Mul": "mul", "Div": "div", "Neg": "neg"}
# Cyrius has one namespace per program, so logging's `init` / `init_with_level` keep their module in the name.
FN_RENAMES = {"init": "bj_log_init", "init_with_level": "bj_log_init_with_level",
              # pub(crate): bijli-internal in the port
              "check_radial_distance": "_bj_check_radial_distance"}
FIELD_RENAMES = {}
EXCEPTIONS = {
    "field::impl From for FieldVector": "BjFieldVector_from_array(p); the other direction, From<FieldVector> for "
                                        "[f64; 3], is the method to_array(p)",
    "field::impl AddAssign for FieldVector": "`a = a + b`: Cyrius refuses `+=` on a struct value (guide, Compound "
                                             "assignment); the same operations",
    "field::impl SubAssign for FieldVector": "`a = a - b`, as AddAssign",
    "field::impl MulAssign for FieldVector": "`a = a * s`, as AddAssign",
    "field::impl Mul for f64": "BjFieldVector_scale_left(s, v): an operator fn belongs to its left operand's type",
    "field::impl Display for FieldVector": "the method display(t), into a BjText (src/fmtf.cyr)",
}


def snake(camel):
    out = []
    for k, ch in enumerate(camel):
        if ch.isupper() and k > 0:
            prev = camel[k - 1]
            nxt = camel[k + 1] if k + 1 < len(camel) else ""
            if prev.islower() or prev.isdigit() or (prev.isupper() and nxt.islower()):
                out.append("_")
        out.append(ch.upper())
    return "".join(out)


def cyrius_index():
    """Every fn, impl method (as Type::m), const, f64 global, enum member, enum and struct field in src/."""
    fns, consts, members, enums, structs = set(), set(), set(), set(), {}
    for path in sorted((root / "src").glob("*.cyr")):
        impl, depth = None, 0
        text = path.read_text()
        for line in text.splitlines():
            s = line.strip()
            m = re.match(r"(?:public |pub )?impl (?:\w+ for )?(\w+) \{", s)
            if m and impl is None:
                impl, depth = m.group(1), 0
            m = re.match(r"(?:public |pub )?fn (\w+)", s)
            if m:
                fns.add(f"{impl}::{m.group(1)}" if impl else m.group(1))
            m = re.match(r"(?:public |pub )?const (\w+) =", s)
            if m:
                consts.add(m.group(1))
            m = re.match(r"(?:public |pub )?var (\w+): f64 =", s)
            if m:
                consts.add(m.group(1))
            m = re.match(r"(BJ_\w+) = -?\d+;", s)
            if m:
                members.add(m.group(1))
            m = re.match(r"(?:public |pub )?enum (\w+)", s)
            if m:
                enums.add(m.group(1))
            if impl is not None:
                depth += line.count("{") - line.count("}")
                if depth <= 0 and s == "}":
                    impl = None
        for m in re.finditer(r"(?:public |pub )?struct (\w+) \{(.*?)\}", text, re.S):
            structs[m.group(1)] = {f.strip().split(":")[0].strip() for f in m.group(2).replace("\n", ";").split(";")
                                   if f.strip() and not f.strip().startswith("#")}
    return fns, consts, members, enums, structs


VIS = r"pub(?:\((?:crate|super)\))?"


def rust_items(module):
    """(kind, owner, name, has_self) for each public item of rust-old/src/<module>.rs outside test modules."""
    text = rustsrc.library_text(root / "rust-old/src" / f"{module}.rs")
    items = []
    owner = enum_name = struct_name = None
    lines = text.splitlines()
    for li, line in enumerate(lines):
        s = line.strip()
        if (m := re.match(r"impl(?:<[^>]*>)? (?:std::ops::|std::fmt::)?(\w+)(?:<[^>]*>)? for (\w+)", s)):
            owner = m.group(2)
            items.append(("trait", owner, m.group(1), None))
            continue
        if (m := re.match(r"impl(?:<[^>]*>)? (\w+)(?:<[^>]*>)? \{", s)):
            owner = m.group(1)
            continue
        if (m := re.match(VIS + r" enum (\w+)", s)):
            enum_name, struct_name = m.group(1), None
            items.append(("enum", None, enum_name, None))
            continue
        if (m := re.match(VIS + r" struct (\w+)", s)):
            struct_name, enum_name = m.group(1), None
            items.append(("struct", None, struct_name, None))
            continue
        if s.startswith("}") and not line.startswith(" "):
            owner = enum_name = struct_name = None
            continue
        if enum_name and (m := re.match(r"([A-Z]\w*)(?:\(|,|\s*\{|$)", s)) and not s.startswith("#"):
            items.append(("variant", enum_name, m.group(1), None))
        if struct_name and (m := re.match(r"pub (\w+):", s)):
            items.append(("field", struct_name, m.group(1), None))
        if (m := re.match(VIS + r" (?:const )?fn (\w+)(?:<[^>]*>)?\((&?(?:mut )?self)?", s)):
            has_self = bool(m.group(2))
            if not has_self and s.endswith("(") and li + 1 < len(lines):
                has_self = bool(re.match(r"&?(?:mut )?self\b", lines[li + 1].strip()))
            items.append(("fn", owner, m.group(1), has_self))
        if (m := re.match(r"pub const (\w+):", s)):
            items.append(("const", owner, m.group(1), None))
    return items


def lib_sync():
    """tests/lib.cyr includes [lib] modules, in order, and cyrius.cyml has no [test] modules (port-parity.md §8)."""
    sections, current = {}, None
    for line in (root / "cyrius.cyml").read_text().split("\n"):
        if re.match(r"\[[\w.]+\]\s*$", line):
            current = line.strip()[1:-1]
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    lists = {}
    for table in ("test", "lib"):
        body = "\n".join(l.split("#")[0] for l in sections.get(table, []))
        m = re.search(r"^modules\s*=\s*\[(.*?)\]", body, re.S | re.M)
        lists[table] = re.findall(r'"([^"]+)"', m.group(1)) if m else None
    if lists["lib"] is None:
        return ["cyrius.cyml has no [lib] modules list"]
    if lists["test"] is not None:
        return ["cyrius.cyml [test] modules is back: test units include tests/lib.cyr instead"]
    mirror = re.findall(r'^include "([^"]+)"', (root / "tests/lib.cyr").read_text(), re.M)
    if mirror != lists["lib"]:
        return [f"tests/lib.cyr includes {mirror}, not [lib] modules {lists['lib']}"]
    return []


def field_store_scan():
    """No `s.f = g(.., t.e, ..)` where g returns a struct and t.e is struct-typed (cyrius
    docs/development/issues/2026-11-09-field-store-of-struct-call-copies-a-field-argument.md): bind the call to a
    local first. Types come from struct fields, fn parameters, `var x: T` and `self` in an impl."""
    files = sorted(list((root / "src").glob("*.cyr")) + list((root / "tests").glob("*.*cyr")))
    texts = {f: f.read_text() for f in files}
    fields = {}
    for text in texts.values():
        for m in re.finditer(r"^struct (\w+)\s*\{(.*?)\}", text, re.M | re.S):
            fields[m.group(1)] = dict(re.findall(r"(\w+)\s*:\s*\*?(\w+)\s*;", m.group(2)))
    sret = set()
    for text in texts.values():
        for m in re.finditer(r"^\s*(?:public |pub )?fn (\w+)\([^)]*\)\s*:\s*(\w+)\s*\{", text, re.M):
            if m.group(2) in fields:
                sret.add(m.group(1))
    store = re.compile(r"([a-z_]\w*(?:\.\w+)+)\s*=\s*(\w+)\(([^;]*)\);")
    ref = re.compile(r"\b([a-z_]\w*)((?:\.[a-z_]\w*)+)\b(?!\s*\()")
    out = []
    for f, text in texts.items():
        types, impl = {}, None
        for i, line in enumerate(text.split("\n"), 1):
            m = re.match(r"^impl (\w+)", line)
            if m:
                impl = m.group(1)
            elif line.startswith("}"):
                impl = None
            m = re.match(r"^\s*(?:public |pub )?fn \w+\(([^)]*)\)", line)
            if m:
                types = dict(re.findall(r"(\w+)\s*:\s*\*?(\w+)", m.group(1)))
                if impl:
                    types["self"] = impl
            for name, ty in re.findall(r"\bvar (\w+)\s*:\s*\*?(\w+)", line):
                types[name] = ty
            for sm in store.finditer(line):
                if sm.group(2) not in sret:
                    continue
                for base, chain in ref.findall(sm.group(3)):
                    ty = types.get(base)
                    for name in chain.split(".")[1:]:
                        ty = fields.get(ty, {}).get(name)
                    if ty in fields:
                        out.append(f"{f.relative_to(root)}:{i}: a struct call with a struct-field argument stored into "
                                   f"a field (cyrius 6.7.18 misstores it; bind to a local): {line.strip()[:90]}")
                        break
    return out


def expected(module, kind, owner, name, has_self):
    """The Cyrius name the rules give, and which index holds it."""
    if kind == "fn":
        key = f"{owner}::{name}" if owner else name
        if key in FN_RENAMES:
            return "fn", FN_RENAMES[key]
        if owner is None:
            return "fn", f"bj_{name}"
        return "fn", (f"Bj{owner}::{name}" if has_self else f"Bj{owner}_{name}")
    if kind == "const":
        return "const", f"BJ_{name}"
    if kind == "enum":
        return "enum", ENUM_NAMES.get(name, f"Bj{name}")
    if kind == "variant":
        return "member", f"{ENUM_PREFIX.get(owner, 'BJ_' + snake(owner) + '_')}{snake(name)}"
    if kind == "struct":
        return "struct", f"Bj{name}"
    if kind == "field":
        return "field", (f"Bj{owner}", FIELD_RENAMES.get((owner, name), [name]))
    if kind == "trait":
        if name == "Default":
            return "fn", f"Bj{owner}_default"
        if name in OPERATORS:
            return "fn", f"Bj{owner}_{OPERATORS[name]}"
        return "exception", None
    return "exception", None


def main() -> int:
    fns, consts, members, enums, structs = cyrius_index()
    problems, ok = [], 0
    used = set()
    modules = rustsrc.rust_modules(root)
    for module in modules:
        if module not in MODULES and module not in PENDING_MODULES and module not in NOT_APPLICABLE_MODULES \
                and module not in IN_PROGRESS:
            problems.append(f"rust-old/src/{module}.rs is in no module list")
    pending_counts = {}
    for module in modules:
        items = rust_items(module)
        if module in PENDING_MODULES:
            pending_counts[module] = len(items)
            continue
        if module not in MODULES and module not in IN_PROGRESS:
            continue
        still = IN_PROGRESS.get(module, ("", set()))[1]
        for kind, owner, name, has_self in items:
            if name in still and kind == "fn":
                pending_counts[module] = pending_counts.get(module, 0) + 1
                where, want = expected(module, kind, owner, name, has_self)
                if where == "fn" and want in fns:
                    problems.append(f"{module}::{name} is ported ({want}): drop it from IN_PROGRESS")
                continue
            label = f"{module}::{owner + '::' if owner else ''}{name}"
            key = label if kind != "trait" else f"{module}::impl {name} for {owner}"
            if key in EXCEPTIONS:
                used.add(key)
                ok += 1
                continue
            where, want = expected(module, kind, owner, name, has_self)
            hit = {
                "fn": lambda w: w in fns,
                "const": lambda w: w in consts,
                "enum": lambda w: w in enums,
                "member": lambda w: w in members,
                "struct": lambda w: w in structs,
                "field": lambda w: w[0] in structs and all(f in structs[w[0]] for f in w[1]),
            }.get(where, lambda w: False)(want)
            if hit:
                ok += 1
            else:
                problems.append(f"{label} ({kind}) -> {want}")
    for key in EXCEPTIONS:
        if key not in used:
            problems.append(f"EXCEPTIONS entry {key!r} names no Rust item of a ported module")
    problems += lib_sync()
    problems += field_store_scan()
    for p in problems:
        print("  " + p)
    pend = sum(pending_counts.values())
    print(f"surface parity: {len(MODULES)} modules ported, {len(IN_PROGRESS)} in progress, {ok} items matched; "
          f"{len(pending_counts)} modules pending ({pend} public items); {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
