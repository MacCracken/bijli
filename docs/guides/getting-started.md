# Getting started with bijli

## Build

```sh
cyrius deps                              # resolve dependencies
cyrius build src/main.cyr build/bijli    # compile
cyrius test                              # run tests/*.tcyr
```

## Layout

- `src/main.cyr` — entry point. Top-level `var r = main(); sys_exit_group(r);`.
  Use `sys_exit_group`, not `syscall(SYS_EXIT, r)` — the latter is exit(2) and ends
  only the calling thread, so it hangs once the program spawns its first worker.
- `tests/` — test suite (`.tcyr` files, auto-discovered by `cyrius test`).
- `tests/golden/` — the golden corpus the parity harness compares against; `scripts/oracle/validate.py` is its
  reference model. Rust bijli 1.2.0, the specification, is read from the `1.2.0` tag (ADR 0003).

## Adding a feature

1. Edit `src/main.cyr` (or add a new module and `include` it).
2. Cross-check parity: golden rows written from validate.py's model of `git show 1.2.0:…`.
3. Add a test case to `tests/bijli.tcyr`.
4. Run `cyrius test`.
5. Bump `VERSION` and add a CHANGELOG entry before tagging.

See [`../adr/template.md`](../adr/template.md) when a non-trivial design choice deserves an ADR.
