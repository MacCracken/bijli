//! Records the compiler that builds the oracle (`ORACLE_RUSTC`): Rust's Unicode tables are the compiler's own, so
//! `src/touchstone_data.cyr` names it, and a toolchain that moves the tables shows up in `--check`.

use std::process::Command;

fn main() {
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".to_string());
    let version = Command::new(rustc)
        .arg("--version")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .unwrap_or_default();
    println!("cargo:rustc-env=ORACLE_RUSTC={}", version.trim());
    println!("cargo:rerun-if-changed=build.rs");
}
