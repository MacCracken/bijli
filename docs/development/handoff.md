# Handoff: the Cyrius port, M6's range-edge review (Part 4) under way (2026-11-17)

For the next agent picking up the port. **Read first:** `CLAUDE.md` (the rules), [`state.md`](state.md) (live state),
[`bijli-port-plan.md`](bijli-port-plan.md) (§0.1 holds the user's M5 Part 2 answers; the tracking table is near the
end), [`port-parity.md`](port-parity.md) §94–§115 (M5) and §116 on (M6), and [ADR 0002](../adr/0002-port-api-divergences.md) (rows 1–42;
M5 touched 2, 3, 9, 13, 17, 18, 36, 40, 41 and 42). This file is the working knowledge those don't hold: how a Part 2
bite is done, the machinery, the scripts and the traps. Delete it once it has served (the M3 one lives at
`git show 01bb43e:docs/development/handoff.md`; its bite workflow for a module's API still applies).

## Where it stands

- **Committed through `53b06f0`** ("M6 pt3", Part 3: compute.rs); part 4a is in the working tree for the user to
  commit.
- **M6 Part 1** (the 2D CPML and `Fdtd2d`, §116–§119) runs in bites. Done: 1a the CPML layer (`src/cpml.cyr`); 1b-1
  `Fdtd2d` without the CPML (`src/fdtd2d.cyr`; its scripts `Fdtd2d::script` in `fdtd_2d3d.tsv`, replayed by
  validate.py's `Fd2Scripts` in pure Python, no NumPy on the host); 1b-2 the CPML attached (the y-strips packed, the
  x-strips scalar: the CPML step 5.7–6.2× Rust's); 1c `nf2ff_2d` (on the phasor form's core with the rotation (1, 0),
  `src/nf2ff.cyr` after `src/fdtd2d.cyr` in the bundle; its edge rows mapped in `fdtd_2d3d.adr.tsv`); 1d the JSON
  of `CpmlParams` and `Mode2d` (src/json.cyr; `m_json.rs`'s `Kind::is_usize` carries the thickness in an f64's bits;
  a thickness past 2^63 − 1 is refused, ADR 0002 row 36, §120). **Part 1 is done.** Part 2 runs in bites: 2a
  `Fdtd3d` without the CPML (`src/fdtd3d.cyr`, §121: `Fdtd3d::script` rows replayed by validate.py's `Fd3Scripts` and
  the harness's `_par_fd3_script`; the 2D file's helpers public as `_bj_fd2_*`); 2b the 3D CPML (§123: `_bj_fd3_slabs`,
  the twelve Ψ slabs in one mapping; the first rows from `scripts/oracle/emit.py`, whose `emit_2b` is the template for
  the next part's emitter). **Part 2 is done**; fdtd.rs is ported whole. Part 3, compute.rs, is done (§124:
  `src/compute.cyr`; probe p20 for the trait and generic shape; `emit_3`'s GpuFdtd2d scripts and CPU kernel rows).
  **Part 4, M6's range-edge review** (plan §0.1's M6 paragraph: the standing answers applied, nothing new asked), runs
  in bites: 4a round once (done, §125: `_bj_e_fd_coeff`, `_bj_e_fd_time`, `_bj_e_fd_courant` in
  `src/crmath_edge.cyr`, the lossy pair through `BjLossyMaterial`; validate.py's `FD2_PORT` and `fd23_runs`, the
  mappings named 40, then 13, else 18; `emit_4a`), then 4b the CPML profile (σ_max and each layer's b, c, 1/κ − 1 at
  their range edges; c's b − 1 as `_bj_expm1`, 1/κ − 1 as (1 − κ)/κ rounded once; every CPML trajectory will change
  with c), 4c the per-cell rule for the 2D updates and the CPU kernels (`_bj_e_cell_row` needs a subtract form and a
  curl form: d ± k·((a1 − a0) − (b1 − b0)), the marker on the product with w = Δa − Δb ≠ 0 as the factor test), 4d
  the 3D updates (the same curl form), 4e the CPML corrections (two statements a cell: Ψ ← bΨ + cΔ, then f +=
  (s·k)·(km1·Δ + Ψ); a product's marker tests it against the other addend, as `_bj_e_cell_edge` against d; Ψ's
  geometric decay makes b·Ψ subnormal in quiet layers), 4f `total_energy` (the 1D's `_bj_fd1_energy_ux` pattern) and
  `run`'s source (`_bj_fd1_source_ux`, to share) in the three grids. Then M6 is done.
- **No Rust (ADR 0003, 2026-11-16).** The Rust tree and the Rust generator are retired; the corpus is frozen data and
  `scripts/oracle/validate.py` the reference model. New rows come from `scripts/oracle/emit.py` (2b's `emit_2b` is
  the first part): add a part function, register it in `main()`'s module list, run `emit.py .`, then `validate.py .
  --write` for the corrections; `emit.py . --check` must pass. A module's frozen rows (FROZEN) are never rewritten. The tools read `scripts/rust_1_2_0.json`;
  the Rust source is `git show 1.2.0:src/fdtd.rs`. Steps 4 and 9 below and the final job have no `cargo` step now.
- **M6's corpus** is `tests/golden/fdtd_2d3d.tsv` (its frozen rows from the retired `m_fdtd_2d3d.rs`, which imported `m_fdtd`'s
  `u`, `seq`, `res_toks`, `digest` and `dbg_*`), checked by validate.py's `validate_fdtd_2d3d` (Rust's chain on
  glibc's primitives) and corrected by `fdtd_2d3d_corrections` (MPFR's), replayed by the harness's `_par_fdtd_2d3d`.
  A class-C function needs its line in `_par_c_bound` (tests/parity.tcyr), or its rows stay "measured".
  `cpml_fits` and `is_attached_to` are public as `_bj_cpml_fits` and `_bj_cpml_is_attached_to` for the grid files.
- **M5 Part 1 is done** (§94–§99): the 1D grid, the pulses and TF/SF, the media, the utilities, NF2FF's phasor form,
  their JSON.
- **M5 Part 2, the range-edge review**, runs on the user's answers of 2026-11-15 (plan §0.1). Done:
  - 2a-1a, the media coefficients rounded once from exact (row 40, §100);
  - 2a-1b, the Mur coefficients: TF/SF rounded once, 1D cancellation-free in its band (row 40, §101);
  - 2a-2, the subpixel means rounded once (row 41, §102);
  - 2a-3, `graded`'s last cell d_max, no cell past it (asked again in the bite), its index past 2³¹ (row 42, §103);
  - 2b, round once (row 13: the update coefficients, `courant_dt`, `time()` past 2⁵³; §104);
  - 2c, exact thresholds (row 17: the TF/SF Courant test, NF2FF's 10⁻¹² agreement, the Nyquist warning; §105);
  - 2d-1, the pulses' range edges (row 18, §106);
  - 2d-2a, the 1D grid's per-cell updates at their range edges (the user's answer; row 18, §107);
  - 2d-2b, the 1D Mur coefficient outside its band and its band at every scale, `total_energy`, `run`'s source
    (rows 18 and 40, §108);
  - 2d-3, the TF/SF grid's updates on the per-cell rule, which moved into `src/crmath_edge.cyr` and was completed for
    products of nonzero factors that round to ±0 (−0 on a +0 field; a ±2^-1075 tie against an odd field below
    2^-1021): 10 1D rows changed with it (row 18, §109);
  - 2d-4a, the CW accumulator: ω kept to 53 bits beside Rust's rounded `omega` (two private fields), the scalars with
    an unbounded exponent at their edges, the sums on the per-cell rule as d + k·a (b = 0), `amplitude` as BjComplex's
    abs; the rule's marker now tests a zero product's factors (row 18, §110);
  - 2d-4b-1, `dft_magnitude`: `_bj_fu_dft` takes ω as a 53-bit pair; Rust's recurrence where it is the unbounded chain
    (ω·dt and the block seeds doubles, |w_r|, |w_i| ≥ 2^-959 once a call, the products s·z checked once a sample,
    re/im finite, the magnitude in range), else `_bj_fu_dft_ux` (row 18, §111);
  - 2d-4b-2, `find_eigenmodes` around it: df, f_i and ω_i to 53 bits (`_bj_fu_freq`), the frequencies rounded once,
    the 10 % floor's 53-bit value compared exactly (`_bj_fu_floor`), Q by Rust's chain or `_bj_fu_quality_ux` (row 18,
    §112). 2d-4 is done;
  - 2d-5, the utilities: `graded` by cell (`_bj_fu_graded_ux`, `_bj_ux_pow` for the ratio), `extents`' sum rerun when
    it ends non-finite, `fourth_order_derivative` by cell (`_bj_fu_fourth_ux`), each check one exponent-field test
    (row 18, §113);
  - 2d-6, NF2FF: ω to 53 bits for the agreement (`_bj_nf_same53`), the frequency, k and the rotation; the core on
    Rust's chain where every operation is its 53-bit value (mode 0, checked per sample by one exponent test and per
    angle), else the call or the angle unbounded (`_bj_nf_angle_ux`); a non-finite field keeps Rust's (row 18, §114).
    **Part 2 is done.**
- **Numbers after M6 part 4a**: 19 modules, 644,480 ok, 0 failed, 122,697 accepted (ADR 01: 43,096; ADR 0002:
  79,601) over 767,177 rows; 488 digests equal on x86_64 ELF, x86_64 PE and aarch64. Ledger 778 ported, 12 n/a, 68
  pending of 858; API snapshot 2,139; surface parity 653 items (15 modules); benches 80 of 80 twinned; execution
  coverage 394 of 394. `json.tsv` is 40.7 MB of the 48 MiB cap: the rest of M6 adds no JSON (`Fdtd3d` and compute.rs
  have no serde), but M7's eight serde types (integration/) will want their own module.
- **M5 is done** (§115). The benchmark record's Cyrius column is each row's least of four full runs: rows vary up to
  1.5× between runs under other repos' loads, so take the least of several when regenerating it.
- **Next**: part 4b, the CPML profile (above).
  `nf2ff_2d` (the snapshot form) calls `_bj_nf_core` with the rotation (1, 0) and k from `_bj_nf_validate` on its
  double frequency (`_bj_k_ux_from`): the core's range edges come with it.
- **Corpus size**: `fdtd.tsv` is 38.8 MB (107,341 rows) and `json.tsv` 40.7 MB, against the 48 MiB cap that
  validate.py's `too_big` guards (§81); M6's rows are in `fdtd_2d3d.tsv` (2.9 MB after 4a).
- **A recorded cost**: the 1D step runs 3.5–3.8× HEAD's with 2d-2a's scan (§107), the TF/SF step about 2× with
  2d-3's (§109) and `accumulate` about 3× with 2d-4a's (§110). The user chose that fix; I had estimated 1.5–2×.
  Cyrius 6.7.18 has no vector verb that classifies a lane (plan §10, listed as a gap). M6's grids will pay the same.
  Don't argue the cost (memory: no cost argument); measure and record it. In a scan, every local and compare on the
  common path costs about a nanosecond a cell, and a rerun about 100 ns: keep the common path to the fewest tests,
  move the rare work into a helper (§109's zero path), and never mark a cell whose result can't change (an exact zero
  product: §110's −0 marker made `accumulate` 12× HEAD's).

## How a Part 2 bite is done

1. **The finding**: read the port's chain (it mirrors Rust's) and list where an intermediate leaves the normal range
   (±∞, NaN, a nonzero subnormal, a product of nonzero factors that underflowed to 0), or where a sum cancels. Write
   down one concrete input per case with Rust's value and the true one.
2. **The port**: keep Rust's chain and compute it first; where the edge is detected for finite inputs, rerun with an
   unbounded exponent (kimiya's `_bj_k_ux_*`) and round the last step once. In-range inputs must keep Rust's bits
   exactly. Non-finite inputs keep Rust's result. Narrow the detection to what actually changes a result (the Gaussian
   pulse needs only t − t₀ = ±∞; see §106), and confirm the narrowing leaves `fdtd.adr.tsv` byte for byte the same.
3. **validate.py**: an independent model on Fractions (`ux`, `ux_mul`, `ux_add`, `ux_div`, `ux_sqrt`, `rnd53`,
   `ux_fin`; `mpfr_ref.wide_split` for e^x and sin past the range). Wire it into the port set: `FD_PORT =
   {"40", "13", "17", "18"}`. `Fd1`, `Tf` and `FdScripts` take a port set, and `fdtd_adr_rows` names a trajectory by
   the first subset that accounts for its whole change: `fs40`, then `fs13`, then `fs17`, else 18. Free functions
   get a `*_port` beside the Rust transcription, mapped in `fdtd_adr_rows`. **Every new mapping must make
   `fdtd_corrections` skip the row**: a golden row never has both a `.cr.tsv` and an `.adr.tsv` entry.
4. **The grid**: a Python emitter in `scripts/oracle/` (ADR 0003): inputs from a seeded grid, each output from
   validate.py's model, rows appended to the module's TSV with `wire()`; skip a repeated (function, inputs). Check
   that the old body is an exact prefix of the new file.
5. **Map and check**: `python3 -I -B scripts/oracle/validate.py . --write` (about 10–15 min now; run it in the
   background), then `cyrius test tests/parity.tcyr`: every `fdtd.adr.tsv` row must be accepted, 0 warnings.
6. **Tests**: `cyrius_*` groups pinning values computed from validate.py's model. Load the model with the exec trick
   (below), not by retyping it. Give long doubles by their bits (`_bj_from_bits(0x…)`); negate a negative's
   positive bits.
7. **Benchmarks**: the interleaved A/B against HEAD (two rounds; a third if a round was disturbed), noting the load,
   plus a focused scratch `.bcyr` for the changed row where useful. Never wait for a quiet machine; check `ps` for
   other repos' jobs and say so.
8. **Docs**: ADR 0002 row (new or `*Amended …*` in both the Cyrius and golden columns), a CHANGELOG Changed entry and
   a Fixed entry (with the cost), `port-parity.md` §N with placeholders until the numbers are in, `state.md`'s M5 row
   and counts, the plan's tracking row, `scripts/oracle/README.md`'s fdtd row. Keep lines ≤ 120 (tables exempt).
9. **The final job** (background, one command, below): the gate, validate.py, then
   wine, then aarch64, one at a time, each with its `=== parity summary` line and the 475 digests compared with ELF's.
   Fill §N's result line, then wrap up. The user commits each bite; never commit or push.

## The machinery Part 2 built

- **`_bj_cr_rat(v, n, spec, slow, out)`** (`src/crmath_edge.cyr`): RN(N/D) for N and D sums of signed monomials in the
  doubles v, by spec string ("+2ea-1bc/+2ea+1bc": sign, coefficient digit 1/2/4/8, a letter per factor, 'a' = v[0]).
  It runs in three paths: a double-double tier, an expansion tier, then exact bignums. `slow` is 0 for all three, 2
  for the expansion tier then exact, 1 for exact only. It returns 0 when D is exactly 0. An exact zero takes IEEE's
  signs: a one-term N's factors', a cancelled sum +0, then D's. `_bj_cr_rat_spec(k)` numbers the specs:
  - 0–9: the media coefficients;
  - 10: the TF/SF Mur coefficient;
  - 11: the 1D Mur numerator;
  - 12 and 13: the subpixel means;
  - 14: step·dt split;
  - 15: the Courant sign;
  - 16: the Nyquist sign.

  `cr_verify.py`'s rat, rat_x and rat_slow families check every spec against Fractions. Add a spec there too, with an
  independent `rat_ref` formula.
- **`_bj_e_within12(a, b, m)`** (`src/crmath_edge.cyr`): |a − b| ≤ 10⁻¹²·m exactly, shared by rf and NF2FF.
- **The per-cell rule** (`src/crmath_edge.cyr`, since 2d-3): `_bj_e_cell_row(d, k, ks, a, b, w, n)` runs d += k·(a − b)
  on the packed verbs with a coefficient row (ks = 8) or one coefficient at k (ks = 0); `_bj_e_cell_fix` scans the
  product row between the multiply and the add, `_bj_e_cell_at` reruns a marked cell (`_bj_e_cell`) and sets its
  product slot to −0 so the packed add leaves the rerun value; `_bj_e_cell_one` is one update (Mur ends, the tail).
  b = 0 means d += k·a (the CW sums). The marker: the product ±∞, NaN or a nonzero subnormal, or ±0 that could move the
  sum (−0 on a +0 field, or against a nonzero field of at most 2^-1021) from nonzero factors (k ≠ 0, a ≠ b). The scan
  takes untyped addresses on purpose (byte arithmetic). validate.py's `fd_cell18` is the rule itself (Fractions
  wherever Rust's product isn't normal or an exact zero), not the port's marker.
- **kimiya's unbounded exponent**:
  - `_bj_k_ux_from` comes **first**: `_bj_k_ux_mul` multiplies significands as doubles, so raw doubles overflow, or
    underflow into a "special".
  - The operations: `_bj_k_ux_mul`, `_div`, `_add`, `_sub`, `_sqrt`, `_neg`, `_cmp` (−1, 0, 1, or 2 for NaN).
  - `_bj_k_ux_round(s, e, 0.0)` rounds into binary64, `_bj_k_ux_mul_f64`, `_div_f64` and `_add_f64` round the last
    step once, `_bj_ux_exp_wide(x)` gives e^x (its far stand-in past 2⁴⁷ is (1, ±2⁶²)), and `_bj_vx_sin53(s, e)` the
    sine (Payne–Hanek past DBL_MAX).
- **kimiya's `_bj_k_cr_mono(&num, nn, &den, nd)`**: RN(Π num / Π den), row 13's workhorse. It returns NaN for a zero
  denominator factor, so guard it.

## Scripts (recreate them in your scratchpad)

The gate (`gate.sh`; run with `export S=<scratch>`, in the background, about 10 minutes):

```bash
#!/bin/bash
cd /home/macro/Repos/bijli
S=${S:?set S to a scratch directory}
F="src/*.cyr tests/*.tcyr tests/*.bcyr tests/*.fcyr tests/*.cyr scripts/oracle/cr_driver.cyr"
run() {
  echo "### $*"; eval "$@" > $S/g.out 2>&1; rc=$?; tail -n 4 $S/g.out; echo "rc=$rc"
  grep -c "warning:" $S/g.out
  grep -E "warnings total|warn line|lint warnings" $S/g.out | sort | uniq -c | head -5
}
run cyrius fmt --check $F
run cyrius lint $F
run cyrius vet src/main.cyr
run cyrius build src/main.cyr build/bijli
run "cyrius test > $S/test.out 2>&1; cat $S/test.out"
run cyrius test --bounds
run cyrius coverage --exec --min 80
run "cyrius distlib --all && cyrius distlib --check"
run cyrius deps --locked
run cyrius api-surface
run python3 -I -B scripts/surface_parity.py .
run python3 -I -B scripts/rust_test_ledger.py . --check
run python3 -I -B scripts/doc_examples.py .
run "cyrius audit > $S/audit.txt 2>&1"
run python3 -I -B scripts/bench_record.py . --check
```

The final job (after the gate): validate.py, then the cross-targets one at a time.

```bash
cd /home/macro/Repos/bijli
python3 -I -B scripts/oracle/validate.py . > $S/validate.txt 2>&1; echo "validate rc=$?"; tail -1 $S/validate.txt
grep -a '^  digest' $S/test.out > $S/dig.elf
cyrius test --win > $S/win.out 2>&1; grep -a '=== parity summary' $S/win.out
grep -a '^  digest' $S/win.out | tr -d '\r' | diff -q $S/dig.elf - && echo "win digests equal"
cyrius test --aarch64 > $S/arm.out 2>&1; grep -a '=== parity summary' $S/arm.out
grep -a '^  digest' $S/arm.out | tr -d '\r' | diff -q $S/dig.elf - && echo "arm digests equal"
```

The A/B (`ab.sh`), with HEAD's tree from `rm -rf $S/head && mkdir -p $S/head && git archive HEAD | tar -x -C $S/head
&& cp -a lib $S/head/`:

```bash
#!/bin/bash
S=<scratch>
for r in 1 2; do
  cat /proc/loadavg > $S/ab_load_h$r.txt
  (cd $S/head && cyrius bench > $S/ab_head_$r.txt 2>&1)
  cat /proc/loadavg > $S/ab_load_n$r.txt
  (cd /home/macro/Repos/bijli && cyrius bench > $S/ab_new_$r.txt 2>&1)
done
cat /proc/loadavg > $S/ab_load_end.txt
```

Its comparison (`abcmp.py head.txt new.txt`: the median ratio and the extremes):

```python
import re, sys
unit = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "µs": 1e3, "ms": 1e6, "s": 1e9}
def load(p):
    d = {}
    for l in open(p, encoding="utf-8", errors="replace"):
        m = re.match(r"\s*(\S+): ([\d.]+)(ps|ns|us|µs|ms|s) avg", l)
        if m: d[m.group(1)] = float(m.group(2)) * unit[m.group(3)]
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
rs = sorted((b[k] / a[k], k) for k in a if k in b)
print(len(rs), "common; median ratio", round(rs[len(rs)//2][0], 3))
for r, k in rs[:3] + rs[-3:]: print(f"{k:45s} {a[k]:10.1f} {b[k]:10.1f} {r:.3f}")
```

The prefix check after a regeneration (copy `fdtd.tsv` to the scratchpad first):

```python
def body(p): return [l for l in open(p, 'rb').read().split(b'\n') if l and not l.startswith(b'#')]
o, n = body(OLD), body('tests/golden/fdtd.tsv')
print("added", len(n) - len(o), "prefix", n[:len(o)] == o)
```

validate.py's model from a script (the exec trick; run from `scripts/oracle`):

```python
import sys; sys.argv = ["validate.py", "/home/macro/Repos/bijli"]
src = open("validate.py").read().replace('if __name__ == "__main__":', 'if False:')
g = {"__name__": "v", "__file__": "/home/macro/Repos/bijli/scripts/oracle/validate.py"}
exec(compile(src, "validate.py", "exec"), g)
# e.g. g["fd_pulse_port"](name, v), g["Fd1"](n, dx, K, frozenset({"40", "13", "17", "18"}))
```

fdtd_2d3d's mappings and corrections alone (M6 part 4a: about a minute, against the full `validate.py --write`'s 15),
then `cyrius test tests/parity.tcyr` (13 s) to check the port against them before the full run. Run from
`scripts/oracle`; drop `"--write"` from argv to print the counts without writing:

```python
import sys, collections
sys.argv = ["validate.py", "/home/macro/Repos/bijli", "--write"]
src = open("validate.py").read().replace('if __name__ == "__main__":', 'if False:')
g = {"__name__": "v", "__file__": "/home/macro/Repos/bijli/scripts/oracle/validate.py"}
exec(compile(src, "validate.py", "exec"), g)
from pathlib import Path
root = Path("/home/macro/Repos/bijli"); path = root / "tests/golden/fdtd_2d3d.tsv"
rows = g["fdtd_2d3d_adr_rows"](root, path)
print("adr rows", len(rows), collections.Counter((r[0], r[3]) for r in rows))
print(g["check_corrections"](root, "fdtd_2d3d", path, g["fdtd_2d3d_corrections"](root, path)))
print(g["check_adr"](root, "fdtd_2d3d", path))
```

The final job as part 4a ran it (`S` set; after the gate, which writes `$S/test.out`):

```bash
cd /home/macro/Repos/bijli
python3 -I -B scripts/oracle/validate.py . --write > $S/validate_w.txt 2>&1; echo "validate --write rc=$?"
tail -3 $S/validate_w.txt; git status --short tests/golden
python3 -I -B scripts/oracle/validate.py . > $S/validate.txt 2>&1; echo "validate rc=$?"; tail -1 $S/validate.txt
grep -a '^  digest' $S/test.out > $S/dig.elf; wc -l < $S/dig.elf
grep -a '=== parity summary' $S/test.out
cyrius test --win > $S/win.out 2>&1; grep -a '=== parity summary' $S/win.out
grep -a '^  digest' $S/win.out | tr -d '\r' | diff -q $S/dig.elf - && echo "win digests equal"
cyrius test --aarch64 > $S/arm.out 2>&1; grep -a '=== parity summary' $S/arm.out
grep -a '^  digest' $S/arm.out | tr -d '\r' | diff -q $S/dig.elf - && echo "arm digests equal"
echo FINAL-DONE
```

A focused bench outside the tree (a stray `tests/*.bcyr` breaks `bench_record.py --check`): put it in
`$S/zzb/` and run `cyrius bench $S/zzb/<name>.bcyr` from each tree in turn (HEAD's in `$S/head`). Part 4a's, the
setup paths (the body's pattern for any one-off row):

```cyrius
include "lib/alloc.cyr"
include "tests/lib.cyr"

fn main() {
    var g2: BjFdtd2d;
    var t2, v2 = BjFdtd2d_new(100, 100, 0.001, BJ_MODE2D_TM, &g2);
    var b2 = bench_new("setup/fdtd2d_set_conductivity_100x100");
    bench_batch_start(b2);
    for (var k2 = 0; k2 < 100; k2 = k2 + 1) {
        var s2: f64 = f64_from(k2 + 1);
        var r2, w2 = g2.set_conductivity(0, 0, 100, 100, s2);
    }
    bench_batch_stop(b2, 100);
    bench_report(b2);
    BjFdtd2d_free(&g2);
    return 0;
}

var rc = main();
sys_exit_group(rc);
```

Rust 1.2.0's sources for reading: `git show 1.2.0:src/fdtd.rs > $S/fdtd.rs` (and `compute.rs`); the CPML profile is
`cpml_coefficients` (line 2040), the strips `correct_fx`/`fy`/`fz` (1869–1987), the 3D slabs `correct_slabs` (1772).

## Traps found in M5 (Cyrius 6.7.18)

- **Waiting on a background job**: never wait on `pgrep -f <pattern>` where the pattern is in your own command line,
  because it matches itself (one waiter sat 2.5 hours on itself). Run the job itself with `run_in_background`, or
  wait with `until grep -q '<its last line>' <its output>; do sleep 15; done`.
- **validate.py and huge powers of two**: a `Fraction(2) ** e` with e near 2⁴⁷ takes trillions of bits; the run was
  killed for memory with **no output at all**. Keep e^x split as (s, e) and decide ±0/±∞ from the binade first
  (`fs_exp53`, `fs_times`; M3's `bm_hg_ux` did the same).
- **Typed pointer arithmetic scales**: `&v + 16` on an `f64` array or a `*f64` adds 16 elements, not bytes. Pass
  untyped words where byte offsets are meant, or index.
- **File-private helpers**: a `private` file's helpers can't be called from another file. Share them through an
  internals file with no marker (`src/crmath_edge.cyr`).
- **Literals**: a small hex value stored in a local `f64` warns ("an integer stored into an f64 slot"); write
  `_bj_from_bits(0x…)`. There are no exponent literals: give long values by their bits.
- **A narrowing must cover zero products**: 2d-2a's marker took ±∞, NaN and nonzero subnormal products and missed a
  product of nonzero factors that rounds to ±0, which still moves d + p (−0 on a +0 field; a ±2^-1075 tie). Model the
  rule itself in validate.py (Fractions wherever Rust's chain might differ), not the port's narrowing, or the two agree
  on the same gap.
- **A finiteness test that overflows**: `dt * dx * 0.0` is NaN at dx = DBL_MAX; test each operand on its own (2b's
  bug, caught by 11 failing mappings).
- **Variables are function-wide**: a repeated `var` name in a test is a compile error.
- **An emitter**: a repeated (function, inputs) row breaks the corrections' keys; skip one already written.
- **Cost model**: a scalar loop over a row costs about 2.6–2.9 ns an element; the `f64v_*` verbs about 0.35 ns. There
  is no vector compare or classify (plan §10).

## Standing preferences (also in the memory directory)

- Full rigor: correctness first; speed is measured and recorded, never a reason to skip a fix (memory: no cost
  argument against edge fixes).
- Ask the user only what is genuinely debatable, with the recommended option first. When a recorded answer turns
  out to have a side effect the question didn't show (2a-3's monotonicity), ask again rather than decide.
- The port removes Rust's faults: each becomes an ADR 0002 row with exact oracle mappings. Never keep one "as Rust",
  never widen a tolerance.
- Review the Cyrius language (guide and a probe on the pinned compiler) before deferring anything.
- Never propose side tasks or separate sessions: findings go in the docs.
- Interleaved A/B, never wait for a quiet machine, note the load and other repos' jobs.
- Never `gh`, never commit or push: the user commits each bite (titles like "m5 2d-2b done").
