//! `src/fdtd_media_data.cyr` — the dispersive material presets of fdtd.rs (M5 part 1c): `DrudeMaterial::gold()` and
//! `silver()`, and `DebyeMaterial::water()`, each field as the oracle's constructor returns it. The values come from
//! the oracle's API, written as positional decimals that read back to the same bits (`cyr::f64_literal`);
//! `tests/golden/fdtd.tsv` pins the three presets besides (`m_fdtd.rs`).

use crate::cyr::{f64_literal, header};
use bijli::fdtd::{DebyeMaterial, DrudeMaterial};
use std::fmt::Write as _;

pub fn emit(prov: &str) -> String {
    let mut out = header(
        "fdtd_media_data.cyr: the dispersive presets of src/fdtd_media.cyr (M5 part 1c)",
        prov,
        &["Each field of DrudeMaterial::gold(), DrudeMaterial::silver() and DebyeMaterial::water()."],
    );
    let gold = DrudeMaterial::gold();
    let silver = DrudeMaterial::silver();
    let water = DebyeMaterial::water();
    let rows: [(&str, f64); 9] = [
        ("_BJ_GOLD_EPS_INF", gold.eps_inf),
        ("_BJ_GOLD_OMEGA_P", gold.omega_p),
        ("_BJ_GOLD_GAMMA", gold.gamma),
        ("_BJ_SILVER_EPS_INF", silver.eps_inf),
        ("_BJ_SILVER_OMEGA_P", silver.omega_p),
        ("_BJ_SILVER_GAMMA", silver.gamma),
        ("_BJ_WATER_EPS_INF", water.eps_inf),
        ("_BJ_WATER_EPS_S", water.eps_s),
        ("_BJ_WATER_TAU", water.tau),
    ];
    for (name, value) in rows {
        assert!(value.is_finite(), "{name} is finite");
        let _ = writeln!(out, "var {name}: f64 = {};", f64_literal(value));
    }
    out
}
