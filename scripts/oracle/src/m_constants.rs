//! `constants.tsv` — every public const of bijli 1.2.0, class A.
//!
//! One row per constant, `const A s:<path> f:<bits>`. The list is shared with the
//! `src/constants.cyr` generator (`t_constants::all`), so the two can't drift apart.

use crate::t_constants;
use crate::tsv::{Class, Module, f, s};

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    for c in t_constants::all() {
        m.row("const", Class::A, &[s(c.path)], &[f(c.value)]);
    }
    m
}
